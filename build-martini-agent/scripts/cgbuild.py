#!/usr/bin/env python3
"""Deterministic bookkeeping for build-martini-agent generate mode (batch).

This file contains NO chemical decisions. Every bead, type, analogue choice,
force constant, evidence class and source is read from the molecule's own
decisions/design.json (authored per molecule by the agent). Bonded terms are not
agent decisions: `bonded` derives them with derive_bonded.py from the SMILES and
the mapping, without MD, under the environment's python_command. This code
otherwise parses SMILES, embeds conformers, runs index queries, writes ITP/GRO/SVG/JSON
records, runs GROMACS and parses its output.

The environment file is $MARTINI_ENV, or environment.json next to the skill folder.

Subcommands (run with the environment's python_command):
  input    <moldir> --smiles S --id ID --name NAME [--class C] [--stereo-note N]
  search   <moldir>
  bonded   <moldir>
  build    <moldir>
  test     <moldir>
  finalize <moldir>
  summary  <batchdir>
"""
import argparse, datetime, hashlib, json, math, os, shutil, subprocess, sys, time
from pathlib import Path

import numpy as np
from rdkit import Chem, RDLogger, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Draw import rdMolDraw2D

RDLogger.DisableLog('rdApp.*')
_TAUT = rdMolStandardize.TautomerEnumerator()

SCRIPT_DIR = Path(__file__).resolve().parent
ENV_PATH = Path(os.environ.get('MARTINI_ENV', SCRIPT_DIR.parents[1] / 'environment.json')).resolve()
PROJECT = ENV_PATH.parent
ENV = json.loads(ENV_PATH.read_text())
GMX = ENV['gromacs_command']
SEED = ENV['random_seed']
CORE_ITP = PROJECT / ENV['force_field']['core_itp']
SOLV_ITP = PROJECT / ENV['force_field']['solvent_itp']
IONS_ITP = PROJECT / ENV['force_field']['ions_itp']
WATER_GRO = PROJECT / ENV['force_field']['solvent_coordinates']
PROD_MDP = PROJECT / ENV['force_field']['production_mdp']
TUT_DIR = PROJECT / 'references/original/tutorial/parametrize_new_small_molecule/ENAP-in-water/6_CG-takeCURRENT'
INDEX = PROJECT / ENV['reference_index']
BIBLE_XLSX = PROJECT / 'references/original/bible/src/MartiniCGSmilesDB_v6.xlsx'
sys.path.insert(0, str(SCRIPT_DIR))
from d5bonded.terms import SIZE_MASS, bead_masses  # noqa: E402


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def wjson(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def rjson(path):
    return json.loads(Path(path).read_text())


TERM_RULES_SOURCE = [{'path': 'build-martini-agent/scripts/d5bonded/terms.py', 'locator': 'module rules'},
                     {'path': 'references/original/small-molecules/models/itps/opt-mono/PBEN.itp', 'locator': '[exclusions] 1 4'},
                     {'path': 'references/derived/tutorial-page.txt', 'locator': 'line 278: exclusions between all beads of a stiff molecule'}]


def bead_size(t):
    return t[0] if t[0] in 'ST' else 'R'


def core_types():
    types, sec = {}, None
    for line in open(CORE_ITP):
        s = line.split(';')[0].strip()
        if not s:
            continue
        if s.startswith('['):
            sec = s.strip('[] ').lower()
            continue
        if sec == 'atomtypes':
            p = s.split()
            types[p[0]] = float(p[1])
    return types


def run_cmd(moldir, args, cwd, tag, stdin=None):
    """Run a command, keep stdout/stderr files and append to logs/commands.jsonl."""
    logdir = Path(moldir) / 'logs'
    logdir.mkdir(exist_ok=True)
    out = logdir / f'{tag}.stdout'
    err = logdir / f'{tag}.stderr'
    t0 = time.time()
    start = now()
    with open(out, 'w') as fo, open(err, 'w') as fe:
        p = subprocess.run(args, cwd=cwd, stdout=fo, stderr=fe, input=stdin, text=True)
    rec = {'tag': tag, 'args': [str(a) for a in args], 'cwd': str(cwd), 'stdin': stdin,
           'exit_code': p.returncode, 'start_utc': start, 'end_utc': now(),
           'wall_s': round(time.time() - t0, 2),
           'stdout': str(out.relative_to(moldir)), 'stderr': str(err.relative_to(moldir))}
    with open(logdir / 'commands.jsonl', 'a') as f:
        f.write(json.dumps(rec) + '\n')
    return p.returncode, out, err


# ---------------------------------------------------------------- geometry
def dist(a, b):
    return float(np.linalg.norm(a - b))


def angle(a, b, c):
    v1, v2 = a - b, c - b
    cosv = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return float(np.degrees(np.arccos(np.clip(cosv, -1, 1))))


def dihedral(a, b, c, d):
    """IUPAC/GROMACS convention, degrees in (-180, 180]. Checked against a GROMACS improper energy."""
    b1, b2, b3 = b - a, c - b, d - c
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return float(np.degrees(np.arctan2(np.linalg.norm(b2) * np.dot(b1, n2), np.dot(n1, n2))))


def circ_mean(vals):
    r = np.radians(vals)
    m = math.degrees(math.atan2(np.mean(np.sin(r)), np.mean(np.cos(r))))
    R = math.hypot(np.mean(np.sin(r)), np.mean(np.cos(r)))
    sd = math.degrees(math.sqrt(max(0.0, -2 * math.log(max(R, 1e-12)))))
    return m, sd


def measure(section, funct, pos):
    n = len(pos)
    if section in ('bonds', 'constraints'):
        return dist(*pos)
    if section == 'angles':
        return angle(*pos)
    if section == 'dihedrals':
        return dihedral(*pos)
    raise ValueError(section)


# ---------------------------------------------------------------- input
def cmd_input(a):
    md = Path(a.moldir)
    (md / 'input').mkdir(parents=True, exist_ok=True)
    mol0 = Chem.MolFromSmiles(a.smiles)
    if mol0 is None:
        raise SystemExit('SMILES parse failed')
    Chem.SanitizeMol(mol0)
    can = Chem.MolToSmiles(mol0, isomericSmiles=True)
    centres = Chem.FindMolChiralCenters(mol0, includeUnassigned=True, useLegacyImplementation=False)
    dbonds = []
    for b in mol0.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE:
            dbonds.append({'atoms_1based': [b.GetBeginAtomIdx() + 1, b.GetEndAtomIdx() + 1], 'stereo': str(b.GetStereo())})
    frags = Chem.GetMolFrags(mol0)
    molH = Chem.AddHs(mol0)
    params = AllChem.ETKDGv3()
    params.randomSeed = SEED
    cids = list(AllChem.EmbedMultipleConfs(molH, numConfs=ENV['max_conformers'], params=params))
    if not cids:
        wjson(md / 'input/embedding_failure.json', {'smiles': a.smiles, 'error': 'embedding failed'})
        raise SystemExit('embedding failed')
    res = AllChem.MMFFOptimizeMoleculeConfs(molH, maxIters=5000)
    confs = []
    for cid, (nc, e) in zip(cids, res):
        m2 = Chem.Mol(molH, confId=cid)
        mc = Chem.Mol(molH)
        mc.RemoveAllConformers()
        mc.AddConformer(Chem.Conformer(molH.GetConformer(cid)), assignId=True)
        Chem.AssignStereochemistryFrom3D(mc)
        cip = Chem.FindMolChiralCenters(mc, includeUnassigned=True, useLegacyImplementation=False)
        confs.append({'conf_id': int(cid), 'mmff_energy_kcal_mol': round(float(e), 4), 'converged': nc == 0,
                      'cip_from_3d': [[i + 1, c] for i, c in cip]})
    ok = [c for c in confs if c['converged']] or confs
    sel = min(ok, key=lambda c: c['mmff_energy_kcal_mol'])
    sel_id = sel['conf_id']
    # SDF: selected and all conformers (atom order = RDKit order, H appended after heavy atoms)
    w = Chem.SDWriter(str(md / 'input/atomistic.sdf'))
    molH.SetProp('_Name', a.id)
    w.write(molH, confId=sel_id)
    w.close()
    w = Chem.SDWriter(str(md / 'input/conformers.sdf'))
    for c in confs:
        molH.SetProp('_Name', f"{a.id}_conf{c['conf_id']}")
        molH.SetProp('mmff_energy_kcal_mol', str(c['mmff_energy_kcal_mol']))
        w.write(molH, confId=c['conf_id'])
    w.close()
    conf = molH.GetConformer(sel_id)
    atoms = []
    for at in molH.GetAtoms():
        i = at.GetIdx()
        p = conf.GetAtomPosition(i)
        atoms.append({'atom_id': i + 1, 'rdkit_index': i, 'element': at.GetSymbol(),
                      'formal_charge': at.GetFormalCharge(), 'aromatic': at.GetIsAromatic(),
                      'hybridization': str(at.GetHybridization()), 'num_h': at.GetTotalNumHs(),
                      'neighbors': [n.GetIdx() + 1 for n in at.GetNeighbors()],
                      'chiral_tag': str(at.GetChiralTag()),
                      'xyz_nm_selected': [round(p.x / 10, 5), round(p.y / 10, 5), round(p.z / 10, 5)]})
    wjson(md / 'input/atoms.json', {
        'atom_index_base': 1, 'rdkit_index_base': 0, 'conversion': 'rdkit_index = atom_id - 1',
        'coordinate_unit': 'nm (converted from RDKit Angstrom by /10)',
        'selected_conformer_id': sel_id, 'atoms': atoms})
    req = {'molecule_id': a.id, 'name': a.name, 'class': a.cls, 'original_smiles': a.smiles,
           'canonical_isomeric_smiles': can, 'formal_charge': Chem.GetFormalCharge(mol0),
           'formula': rdMolDescriptors.CalcMolFormula(mol0),
           'heavy_atoms': mol0.GetNumHeavyAtoms(), 'disconnected_components': len(frags),
           'stereocentres': [[i + 1, c] for i, c in centres],
           'unspecified_stereocentres': [i + 1 for i, c in centres if c == '?'],
           'double_bond_stereo': dbonds, 'test_set_stereo_note': a.stereo_note,
           'overrides': {}, 'environment_json': str(ENV_PATH), 'environment_sha256': sha256(ENV_PATH),
           'config_snapshot': ENV}
    wjson(md / 'input/request.json', req)
    wjson(md / 'input/conformers.json', {
        'method': 'RDKit ETKDGv3 EmbedMultipleConfs + MMFF94 MMFFOptimizeMoleculeConfs(maxIters=5000)',
        'rdkit_version': rdBase.rdkitVersion, 'random_seed': SEED, 'requested': ENV['max_conformers'],
        'embedded': len(cids), 'selected_conformer_id': sel_id,
        'selection_rule': 'lowest MMFF94 energy among converged conformers',
        'note': 'Minimised conformers are not a thermodynamic ensemble.', 'conformers': confs})
    print('input ok', a.id, can, 'confs', len(cids), 'selected', sel_id)


# ---------------------------------------------------------------- search
def bible_rows():
    import openpyxl
    wb = openpyxl.load_workbook(BIBLE_XLSX, read_only=True)
    rows = list(wb.worksheets[0].iter_rows(values_only=True))
    out = []
    for i, r in enumerate(rows[1:], start=2):
        if not r[2]:
            continue
        out.append({'reference_id': f'BIBLE:row{i}', 'name': r[0], 'moltag': r[1], 'smiles': r[2], 'cgsmiles': r[3],
                    'source_path': str(BIBLE_XLSX.relative_to(PROJECT)), 'locator': f'sheet All, row {i}'})
    return out


def tautomer_key(mol):
    """Canonical-tautomer SMILES without stereo, so tautomers of one molecule compare equal; None if it fails."""
    try:
        m = Chem.Mol(mol)
        Chem.RemoveStereochemistry(m)
        return Chem.MolToSmiles(_TAUT.Canonicalize(m), isomericSmiles=False)
    except Exception:  # noqa: BLE001  (recorded as no tautomer match, exact comparison still applies)
        return None


def cmd_search(a):
    md = Path(a.moldir)
    des = rjson(md / 'decisions/design.json')
    req = rjson(md / 'input/request.json')
    target = Chem.MolFromSmiles(req['original_smiles'])
    tcan = Chem.MolToSmiles(target, isomericSmiles=False)
    entries = [json.loads(l) for l in open(INDEX)]
    bib = bible_rows()
    hits, exposure = [], []
    ann = des['analogue_annotations']
    for q in des['analogue_queries']:
        patt = Chem.MolFromSmarts(q['smarts']) if q.get('smarts') else None
        tags = set(q.get('motif_tags', []))
        for e in entries:
            s = e.get('canonical_isomeric_smiles')
            m = Chem.MolFromSmiles(s) if s else None
            matched = []
            if patt is not None and m is not None:
                matched = [[i + 1 for i in t] for t in m.GetSubstructMatches(patt)]
            tag_hit = sorted(tags & set(e.get('motifs', [])))
            if matched or tag_hit:
                hits.append({'query': q['name'], 'reference_id': e['reference_id'], 'name': e['name'],
                             'source_path': e['source_path'], 'canonical_smiles': s,
                             'matched_atoms_1based_in_reference_smiles': matched[:6], 'motif_tags_hit': tag_hit})
        if patt is not None and q.get('search_bible', True):
            for r in bib:
                m = Chem.MolFromSmiles(r['smiles'])
                if m is None:
                    continue
                matched = [[i + 1 for i in t] for t in m.GetSubstructMatches(patt)]
                if matched:
                    hits.append({'query': q['name'], 'reference_id': r['reference_id'], 'name': r['name'],
                                 'source_path': r['source_path'], 'locator': r['locator'], 'canonical_smiles': r['smiles'],
                                 'cgsmiles': r['cgsmiles'], 'matched_atoms_1based_in_reference_smiles': matched[:6],
                                 'motif_tags_hit': []})
    # exposure: same constitution (exact), or the same molecule in another tautomer (RDKit canonical tautomer)
    ttau = tautomer_key(target)
    tautomer_exposure = []
    refs = [(e['reference_id'], e.get('canonical_isomeric_smiles')) for e in entries] + \
           [(r['reference_id'], r['smiles']) for r in bib]
    for rid, s in refs:
        m = Chem.MolFromSmiles(s) if s else None
        if m is None:
            continue
        if Chem.MolToSmiles(m, isomericSmiles=False) == tcan:
            exposure.append(rid)
        elif ttau is not None and tautomer_key(m) == ttau:
            tautomer_exposure.append(rid)
    # merge by reference id (one record per reference, listing the queries it matched)
    merged = {}
    for h in hits:
        k = h['reference_id']
        if k not in merged:
            merged[k] = dict(h)
            merged[k]['queries'] = [h['query']]
            del merged[k]['query']
        elif h['query'] not in merged[k]['queries']:
            merged[k]['queries'].append(h['query'])
    out_hits = []
    for k, h in merged.items():
        an = ann.get(k)
        if an:
            h.update({'used': an['used'], 'used_for': an.get('used_for', []), 'reason': an['reason'],
                      'opened': an.get('opened', [])})
        else:
            h.update({'used': False, 'used_for': [], 'reason': des['default_rejection_reason'], 'opened': []})
        out_hits.append(h)
    missing = [k for k in ann if k not in merged]
    wjson(md / 'decisions/analogue_search.json', {
        'schema_version': 1, 'molecule_id': req['molecule_id'], 'run_utc': now(),
        'code_path': 'scripts/cgbuild.py search', 'index': str(INDEX.relative_to(PROJECT)), 'index_sha256': sha256(INDEX),
        'bible_db': str(BIBLE_XLSX.relative_to(PROJECT)), 'bible_db_sha256': sha256(BIBLE_XLSX),
        'queries': des['analogue_queries'], 'n_hits': len(out_hits),
        'n_used': sum(1 for h in out_hits if h['used']),
        'exact_target_matches': exposure, 'tautomer_target_matches': tautomer_exposure,
        'tautomer_method': f'rdMolStandardize.TautomerEnumerator canonical tautomer, no stereo (RDKit {rdBase.rdkitVersion})',
        'reference_exposed': bool(exposure or tautomer_exposure),
        'annotated_but_not_hit_by_queries': [{'reference_id': k, **ann[k]} for k in missing],
        'hits': out_hits})
    print('search ok', req['molecule_id'], 'hits', len(out_hits), 'exposed', exposure, 'tautomer', tautomer_exposure,
          'annot-not-hit', missing)


# ---------------------------------------------------------------- bonded
def cmd_bonded(a):
    md = Path(a.moldir).resolve()
    args = ENV['python_command'] + [str(SCRIPT_DIR / 'derive_bonded.py'), str(md), '--env', str(ENV_PATH)]
    rc, out, err = run_cmd(md, args, md, 'derive_bonded')
    print(Path(out).read_text().strip())
    if rc != 0:
        wjson(md / 'decisions/bonded_failure.json', {'exit_code': rc, 'stderr_tail': Path(err).read_text()[-1500:]})
        raise SystemExit(f'bonded failed for {md.name}: see logs/derive_bonded.stderr')


# ---------------------------------------------------------------- build
def load_conformers(md):
    sup = Chem.SDMolSupplier(str(md / 'input/conformers.sdf'), removeHs=False)
    mols = [m for m in sup]
    sel = Chem.SDMolSupplier(str(md / 'input/atomistic.sdf'), removeHs=False)[0]
    return mols, sel


def bead_weights(design, mol):
    """Return list of dict atom_idx(0-based)->weight for each bead, H inheriting parent allocation."""
    heavy = {}
    for bi, b in enumerate(design['beads']):
        for aid, wgt in b['members'].items():
            heavy.setdefault(int(aid), {})[bi] = float(wgt)
    out = [dict() for _ in design['beads']]
    excl = set(design.get('hydrogen_excluded_parents', {}).get('atom_ids', []))
    for at in mol.GetAtoms():
        i = at.GetIdx()
        if at.GetAtomicNum() == 1:
            parent = at.GetNeighbors()[0].GetIdx() + 1
            alloc = {} if parent in excl else heavy.get(parent, {})
        else:
            alloc = heavy.get(i + 1, {})
        for bi, wgt in alloc.items():
            out[bi][i] = wgt
    return out, heavy


def cog(weights, xyz):
    idx = np.array(list(weights.keys()))
    w = np.array(list(weights.values()))
    return (xyz[idx] * w[:, None]).sum(0) / w.sum()


def bead_positions(design, weights, xyz):
    pos = []
    for bi, b in enumerate(design['beads']):
        pos.append(cog(weights[bi], xyz))
    pos = np.array(pos)
    mapped = pos.copy()
    for vs in design.get('virtual_sites', []):
        site = vs['site'] - 1
        cons = [c - 1 for c in vs['constructing']]
        if vs['construction'] == 'virtual_sitesn_cog':
            pos[site] = mapped[cons].mean(0)
        elif vs['construction'] == 'virtual_sites3':
            j, k, l = cons
            pos[site] = mapped[j] + vs['a'] * (mapped[k] - mapped[j]) + vs['b'] * (mapped[l] - mapped[j])
        else:
            raise ValueError(vs['construction'])
    return pos, mapped


def cmd_build(a):
    md = Path(a.moldir)
    des = rjson(md / 'decisions/design.json')
    req = rjson(md / 'input/request.json')
    bonded = rjson(md / 'decisions/bonded.json')
    for vs in des.get('virtual_sites', []):
        vs.update(next((x for x in bonded.get('virtual_site_parameters', []) if x['site'] == vs['site']), {}))
    mols, sel = load_conformers(md)
    types = core_types()
    weights, heavy = bead_weights(des, sel)
    nb = len(des['beads'])
    static = []

    def chk(name, ok, detail):
        static.append({'check': name, 'result': 'pass' if ok else 'fail', 'detail': detail})

    # --- static mapping checks
    heavy_ids = [at.GetIdx() + 1 for at in sel.GetAtoms() if at.GetAtomicNum() > 1]
    cover = {h: round(sum(heavy.get(h, {}).values()), 6) for h in heavy_ids}
    chk('heavy_atom_coverage_and_weight_sum', all(abs(v - 1) < 1e-6 for v in cover.values()),
        {str(k): v for k, v in cover.items()})
    conn = {}
    for bi, b in enumerate(des['beads']):
        ids = [int(k) for k in b['members']]
        sub = set(ids)
        seen, stack = {ids[0]}, [ids[0]]
        while stack:
            x = stack.pop()
            for n in sel.GetAtomWithIdx(x - 1).GetNeighbors():
                if n.GetIdx() + 1 in sub and n.GetIdx() + 1 not in seen:
                    seen.add(n.GetIdx() + 1)
                    stack.append(n.GetIdx() + 1)
        conn[b['name']] = len(seen) == len(sub)
    chk('bead_members_connected', all(conn.values()), conn)
    bad_types = [b['type'] for b in des['beads'] if b['type'] not in types]
    chk('bead_types_exist_in_core_itp', not bad_types, {'core_itp': str(CORE_ITP.relative_to(PROJECT)), 'unknown': bad_types})
    netq = sum(b.get('charge_e', 0.0) for b in des['beads'])
    chk('net_charge_matches_formal_charge', abs(netq - req['formal_charge']) < 1e-9,
        {'model_charge_e': netq, 'formal_charge_e': req['formal_charge']})
    # --- geometry over conformers
    all_pos, all_mapped = [], []
    for m in mols:
        xyz = m.GetConformer().GetPositions() / 10.0
        p, mp = bead_positions(des, weights, xyz)
        all_pos.append(p)
        all_mapped.append(mp)
    selxyz = sel.GetConformer().GetPositions() / 10.0
    selpos, selmapped = bead_positions(des, weights, selxyz)
    vs_disc = []
    for vs in des.get('virtual_sites', []):
        s = vs['site'] - 1
        d = [float(np.linalg.norm(p[s] - mp[s])) for p, mp in zip(all_pos, all_mapped)]
        vs_disc.append({'site': vs['site'], 'max_nm': round(max(d), 4), 'mean_nm': round(float(np.mean(d)), 4),
                        'selected_nm': round(float(np.linalg.norm(selpos[s] - selmapped[s])), 4)})
    # --- masses
    # computed, never taken from design.json: the stage-4 record, or the same rule if it predates it
    if 'masses' in bonded:
        masses = [float(x) for x in bonded['masses']['bead_amu']]
        mass_redistribution = bonded['masses']['virtual_site_redistribution']
    else:
        m_arr, mass_redistribution = bead_masses([b['type'] for b in des['beads']], des.get('virtual_sites', []),
                                                 [[x - 1 for x in g] for g in bonded['rigid_groups']])
        masses = [float(x) for x in m_arr]
    for b, m in zip(des['beads'], masses):
        b['mass_amu'] = m
    ignored = [{'bead': b['name'], 'design_amu': b0['mass_amu'], 'used_amu': m}
               for b, b0, m in zip(des['beads'], rjson(md / 'decisions/design.json')['beads'], masses)
               if b0.get('mass_amu') is not None and abs(float(b0['mass_amu']) - m) > 1e-6]
    std_total = sum(SIZE_MASS[bead_size(b['type'])] for b in des['beads'])
    chk('total_mass_preserved_under_vs_redistribution', abs(sum(masses) - std_total) < 1e-6,
        {'model_total_amu': sum(masses), 'standard_size_total_amu': std_total,
         'virtual_site_redistribution': mass_redistribution, 'design_mass_amu_ignored': ignored})
    # --- terms: all values come from decisions/bonded.json
    params = []
    for t in bonded['terms']:
        rec = dict(t)
        rec['bead_ids'] = rec.pop('beads')
        rec['bead_names'] = [des['beads'][i - 1]['name'] for i in rec['bead_ids']]
        params.append(rec)
    exclusions = bonded['exclusions']
    close = [e['beads'] for e in bonded.get('close_pair_exclusions', [])]
    for i, (x, y) in enumerate(exclusions):
        params.append({'term_id': f'x{i + 1}', 'section': 'exclusions', 'bead_ids': [x, y],
                       'evidence_class': 'ensemble_derived', 'source': TERM_RULES_SOURCE,
                       'reason': 'close pair: mapped distance below the LJ sigma (Souza SI C3 p.44)' if [x, y] in close else 'pair within a rigid group, or a substituent and a rigid-group bead it is not bonded to'})
    params.append({'term_id': 'nrexcl', 'section': 'moleculetype', 'value': bonded['nrexcl'],
                   'evidence_class': 'documented_default',
                   'source': [{'path': 'references/derived/tutorial-page.txt', 'locator': 'line 233: standard exclusion number 1'}]})
    for vs in des.get('virtual_sites', []):
        three = vs['construction'] == 'virtual_sites3'
        params.append({'term_id': f"vs{vs['site']}", 'section': 'virtual_sites3' if three else 'virtual_sitesn', 'site': vs['site'],
                       'function': 1, 'function_meaning': 'in-plane site, coefficients a and b from the ensemble' if three else 'centre of geometry of constructing particles',
                       'a': vs.get('a'), 'b': vs.get('b'),
                       'constructing': vs['constructing'], 'evidence_class': vs['evidence_class'], 'source': vs['source'],
                       'discrepancy_vs_mapped_cog': [d for d in vs_disc if d['site'] == vs['site']]})
    for b, m in zip(des['beads'], masses):
        params.append({'term_id': f"mass_{b['name']}", 'section': 'atoms.mass', 'bead': b['name'], 'mass_amu': m,
                       'evidence_class': 'computed',
                       'rule': 'standard size mass; virtual sites 0, their mass split equally over the real beads of the rigid group (d5bonded.terms.bead_masses)',
                       'source': des['mass_policy']['source']})
    # --- static term checks
    keyset, dup = set(), []
    for t in bonded['terms']:
        k = (t['section'], tuple(t['beads']) if t['beads'][0] < t['beads'][-1] else tuple(reversed(t['beads'])), t['function'])
        if k in keyset:
            dup.append(t['term_id'])
        keyset.add(k)
    chk('no_duplicate_terms', not dup, dup)
    # SI C4 p.46: dihedrals i-j-k-l are unstable when angle i-j-k or j-k-l reaches 0/180 deg
    near180 = []
    for d in [t for t in bonded['terms'] if t['section'] == 'dihedrals']:
        q = d['beads']
        for p in [t for t in bonded['terms'] if t['section'] == 'angles']:
            ab = p['beads']
            if (ab in (q[:3], q[1:], q[:3][::-1], q[1:][::-1])) and p['parameters'][0]['value'] > 160:
                near180.append([d['term_id'], p['term_id']])
    chk('no_dihedral_spanning_an_angle_with_theta0_above_160deg (SI C4 p.46)', not near180, near180)
    k_low = [t['term_id'] for t in bonded['terms'] if t['section'] == 'bonds' and t['parameters'][1]['value'] < 700]
    chk('bond_k_not_below_700 (Souza SI C3 p.44)', not k_low, k_low)
    vs_ok = all(d['max_nm'] < 0.05 for d in vs_disc)
    chk('virtual_site_vs_mapped_discrepancy_below_0.05nm', vs_ok, vs_disc)
    # graph connectivity of the CG model (bonds+constraints+vs)
    adj = {i: set() for i in range(1, nb + 1)}
    for t in bonded['terms']:
        if t['section'] in ('bonds', 'constraints'):
            i, j = t['beads']
            adj[i].add(j)
            adj[j].add(i)
    for vs in des.get('virtual_sites', []):
        for c in vs['constructing']:
            adj[vs['site']].add(c)
            adj[c].add(vs['site'])
    seen, stack = {1}, [1]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    chk('cg_graph_connected', len(seen) == nb, {'reached': sorted(seen), 'n_beads': nb})
    big = [f['beads'] for f in bonded['flags'] if f['flag'] == 'rigid_group_above_four_beads']
    chk('rigid_groups_have_at_most_four_real_beads (Souza SI C4 p.45: use a virtual-site hinge)', not big, big)
    # --- mapping.json
    resname = des['resname']
    beads_out = []
    for bi, b in enumerate(des['beads']):
        vsdef = None
        for vs in des.get('virtual_sites', []):
            if vs['site'] == bi + 1:
                vsdef = {'function': 'virtual_sites3 funct 1' if vs['construction'] == 'virtual_sites3' else 'virtual_sitesn funct 1 (COG)',
                         'constructing': vs['constructing'], 'a': vs.get('a'), 'b': vs.get('b')}
        hmem = sorted(i + 1 for i in weights[bi] if sel.GetAtomWithIdx(i).GetAtomicNum() == 1)
        beads_out.append({'bead_id': bi + 1, 'name': b['name'], 'type': b['type'], 'size': bead_size(b['type']),
                          'charge_e': b.get('charge_e', 0.0), 'mass_amu': masses[bi],
                          'members': [{'atom_id': int(k), 'element': sel.GetAtomWithIdx(int(k) - 1).GetSymbol(),
                                       'allocation_weight': float(v)} for k, v in b['members'].items()],
                          'hydrogen_members': {'atom_ids': hmem, 'convention': 'each H inherits the allocation weight of its bonded heavy atom, except H on atoms listed in hydrogen_excluded_parents'},
                          'virtual_site': vsdef, 'sources': b['sources'], 'rationale': b['rationale'],
                          'evidence_status': b['evidence_status'], 'typing_ladder': b['typing_ladder'],
                          'provisional_type': b['provisional_type'],
                          'cog_selected_conformer_nm': [round(float(x), 4) for x in selmapped[bi]]})
    wjson(md / 'decisions/mapping.json', {
        'schema_version': 1, 'molecule_id': req['molecule_id'], 'atom_index_base': 1,
        'atomistic_structure': 'input/atomistic.sdf', 'atomistic_structure_sha256': sha256(md / 'input/atomistic.sdf'),
        'coordinate_unit': 'nm',
        'coordinate_convention': 'bead position = centre of geometry of member atoms (heavy + H), allocation weights normalised within the bead (Souza SI C3 p.41-42 COG rule); RDKit Angstrom converted to nm by /10',
        'hydrogen_excluded_parents': des.get('hydrogen_excluded_parents'),
        'excluded_hydrogens': sorted(at.GetIdx() + 1 for at in sel.GetAtoms() if at.GetAtomicNum() == 1 and
                                     at.GetNeighbors()[0].GetIdx() + 1 in set(des.get('hydrogen_excluded_parents', {}).get('atom_ids', []))),
        'chosen_candidate': des['chosen_candidate'], 'n_beads': nb, 'net_charge_e': netq,
        'resname': resname, 'beads': beads_out})
    wjson(md / 'decisions/mapping_candidates.json', {
        'schema_version': 1, 'molecule_id': req['molecule_id'], 'candidates': des['mapping_candidates'],
        'chosen': des['chosen_candidate'], 'selection_rationale': des['selection_rationale']})
    wjson(md / 'decisions/parameters.json', {
        'schema_version': 1, 'molecule_id': req['molecule_id'], 'units': {'length': 'nm', 'angle': 'deg', 'energy': 'kJ mol-1'},
        'nrexcl': bonded['nrexcl'], 'terms': params, 'bonded_record': 'decisions/bonded.json'})
    # --- ITP
    L = []
    L.append(f"; {req['name']} ({req['molecule_id']}) Martini 3.0.0 candidate model")
    L.append(f"; SMILES: {req['original_smiles']}  canonical: {req['canonical_isomeric_smiles']}")
    L.append('; Generated by build-martini-agent generate mode: agent mapping and types; bonded terms from an')
    L.append(f"; energy-weighted conformer ensemble, without MD ({bonded['energy_model']['name']}).")
    L.append('; Provenance of every term: decisions/parameters.json. Status: see validation.json. Not physically validated.')
    L.append('')
    L.append('[ moleculetype ]')
    L.append('; molname  nrexcl')
    L.append(f'  {resname}  {bonded["nrexcl"]}')
    L.append('')
    L.append('[ atoms ]')
    L.append('; id type resnr resname atom cgnr charge mass')
    for bi, b in enumerate(des['beads']):
        L.append(f"  {bi+1:<3d} {b['type']:<6s} 1  {resname:<5s} {b['name']:<4s} {bi+1:<3d} {b.get('charge_e',0.0):6.3f} {masses[bi]:6.1f}")
    by = {sec: [t for t in bonded['terms'] if t['section'] == sec] for sec in ('bonds', 'constraints', 'angles', 'dihedrals')}
    val = lambda t, n: next(x['value'] for x in t['parameters'] if x['name'] == n)
    if by['bonds']:
        L += ['', '[ bonds ]', '; i j funct b0(nm) k(kJ/mol/nm2)']
        for t in by['bonds']:
            L.append(f"  {t['beads'][0]} {t['beads'][1]}  1  {val(t, 'b0'):.3f}  {val(t, 'k'):.0f}  ; {t['term_id']} {t['role']}")
    if by['constraints']:
        L += ['', '#ifdef FLEXIBLE', '; flexible variant used only for energy minimisation (define = -DFLEXIBLE)', '[ bonds ]']
        for t in by['constraints']:
            L.append(f"  {t['beads'][0]} {t['beads'][1]}  1  {val(t, 'b0'):.3f}  {t['flexible_variant']['k']}  ; {t['term_id']}")
        L += ['#else', '[ constraints ]', '; i j funct b0(nm)']
        for t in by['constraints']:
            L.append(f"  {t['beads'][0]} {t['beads'][1]}  1  {val(t, 'b0'):.3f}  ; {t['term_id']} {t['role']}")
        L.append('#endif')
    if by['angles']:
        L += ['', '[ angles ]', '; i j k funct theta0(deg) k   (funct 2 = cosine-harmonic)']
        for t in by['angles']:
            L.append(f"  {' '.join(map(str, t['beads']))}  2  {val(t, 'theta0'):.1f}  {val(t, 'k'):.2f}  ; {t['term_id']}")
    if by['dihedrals']:
        L += ['', '[ dihedrals ]', '; i j k l funct phi0(deg) k   (funct 2 = improper); i j k l 9 phi_s k n (Fourier)']
        for t in by['dihedrals']:
            q = ' '.join(map(str, t['beads']))
            if t['function'] == 2:
                L.append(f"  {q}  2  {val(t, 'phi0'):.1f}  {val(t, 'k'):.1f}  ; {t['term_id']} {t['role']}")
            else:
                for x in t['parameters'][::2]:
                    n = x['multiplicity']
                    L.append(f"  {q}  9  {x['value']:.2f}  {val(t, f'k_{n}'):.3f}  {n}  ; {t['term_id']}")
    sites3 = [v for v in des.get('virtual_sites', []) if v['construction'] == 'virtual_sites3']
    sitesn = [v for v in des.get('virtual_sites', []) if v['construction'] != 'virtual_sites3']
    if sites3:
        L += ['', '[ virtual_sites3 ]', '; site j k l funct a b   (r_i = r_j + a r_jk + b r_jl)']
        for vs in sites3:
            L.append(f"  {vs['site']}  {' '.join(map(str, vs['constructing']))}  1  {vs['a']:.4f}  {vs['b']:.4f}")
    if sitesn:
        L += ['', '[ virtual_sitesn ]', '; site funct(1=COG) constructing']
        for vs in sitesn:
            L.append(f"  {vs['site']}  1  {' '.join(map(str, vs['constructing']))}")
    if exclusions:
        L += ['', '[ exclusions ]']
        for x, y in exclusions:
            L.append(f'  {x} {y}')
    (md / 'model').mkdir(exist_ok=True)
    (md / 'model/molecule.itp').write_text('\n'.join(L) + '\n')
    # --- GRO (selected conformer, centred in 4.5 nm box as in tutorial run_1mol_CG_system.sh)
    box = max(4.5, float(np.ptp(selpos, axis=0).max()) + 2 * 1.5)
    c = selpos.mean(0)
    G = [f"{req['name']} CG from selected conformer COG", f'{nb:5d}']
    for bi, b in enumerate(des['beads']):
        x, y, z = selpos[bi] - c + box / 2
        G.append(f"{1:5d}{resname:<5s}{b['name']:>5s}{bi+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}")
    G.append(f'{box:10.5f}{box:10.5f}{box:10.5f}')
    (md / 'model/molecule.gro').write_text('\n'.join(G) + '\n')
    # --- drawing
    draw_mapping(md, des, req)
    wjson(md / 'test/static_checks.json', {'checks': static, 'virtual_site_discrepancy': vs_disc})
    nfail = sum(1 for s in static if s['result'] == 'fail')
    print('build ok', req['molecule_id'], 'beads', nb, 'static fails', nfail,
          [s['check'] for s in static if s['result'] == 'fail'])


def draw_mapping(md, des, req):
    mol = Chem.MolFromSmiles(req['original_smiles'])
    AllChem.Compute2DCoords(mol)
    palette = [(0.90, 0.38, 0.33), (0.30, 0.60, 0.85), (0.40, 0.75, 0.40), (0.95, 0.70, 0.25),
               (0.65, 0.45, 0.80), (0.35, 0.80, 0.80), (0.85, 0.50, 0.70), (0.60, 0.60, 0.35)]
    hl, cols, radii = [], {}, {}
    for bi, b in enumerate(des['beads']):
        for aid, wgt in b['members'].items():
            i = int(aid) - 1
            prev = mol.GetAtomWithIdx(i).GetPropsAsDict().get('molAtomMapNumber')
            note = mol.GetAtomWithIdx(i).GetProp('atomNote') if mol.GetAtomWithIdx(i).HasProp('atomNote') else ''
            lab = f"{b['name']}" + ('' if float(wgt) == 1 else f"({float(wgt):.2g})")
            mol.GetAtomWithIdx(i).SetProp('atomNote', (note + '/' if note else '') + lab)
            if i not in cols:
                hl.append(i)
                cols[i] = palette[bi % len(palette)]
                radii[i] = 0.45
    for at in mol.GetAtoms():
        at.SetProp('atomLabel', f"{at.GetSymbol()}{at.GetIdx()+1}")
    d = rdMolDraw2D.MolDraw2DSVG(560, 420)
    opts = d.drawOptions()
    opts.annotationFontScale = 0.6
    opts.legendFontSize = 14
    legend = f"{req['molecule_id']}: " + ', '.join(f"{b['name']}={b['type']}" for b in des['beads'])
    d.DrawMolecule(mol, highlightAtoms=hl, highlightAtomColors=cols, highlightAtomRadii=radii, legend=legend)
    d.FinishDrawing()
    (md / 'model/mapping.svg').write_text(d.GetDrawingText())


# ---------------------------------------------------------------- test
def read_mdp(path):
    out = []
    for line in open(path):
        out.append(line.rstrip('\n'))
    return out


def mdp_set(lines, key, value, note):
    norm = lambda s: s.split('=')[0].strip().replace('_', '-').lower()
    k = key.replace('_', '-').lower()
    done = False
    res = []
    for l in lines:
        s = l.split(';')[0]
        if '=' in s and norm(s) == k:
            if value is None:
                res.append(f'; REMOVED ({note}): {l.strip()}')
            else:
                res.append(f'{key:<24s} = {value} ; ADAPTED: {note}')
            done = True
        else:
            res.append(l)
    if not done and value is not None:
        res.append(f'{key:<24s} = {value} ; ADDED: {note}')
    return res


def make_mdps(td, nsteps_smoke, uses_flexible):
    adapt = []
    em = read_mdp(TUT_DIR / 'martini_em.mdp')
    emad = [('ns_type', None, 'obsolete option removed in GROMACS 2026 (unknown-option warning)'),
            ('energygrps', None, 'not supported with Verlet/GPU-era GROMACS; no energy groups needed')]
    if not uses_flexible:
        emad.append(('define', None, 'model has no constraints, so no FLEXIBLE macro; GROMACS 2026 grompp warns on unused -D macros'))
    for k, v, n in emad:
        em = mdp_set(em, k, v, n)
        adapt.append({'mdp': 'em.mdp', 'key': k, 'value': v, 'reason': n})
    (td / 'em.mdp').write_text('\n'.join(em) + '\n')
    eq = read_mdp(TUT_DIR / 'martini_eq.mdp')
    eqa = [('ns_type', None, 'obsolete option removed in GROMACS 2026'),
           ('energygrps', None, 'not needed'),
           ('Pcoupl', 'c-rescale', 'Berendsen barostat replaced by stochastic cell rescaling (same as adapted production MDP)'),
           ('ref_t', '298.15', 'environment default_temperature_K'),
           ('gen_vel', 'yes', 'start dynamics from minimised coordinates with Maxwell velocities'),
           ('gen_temp', '298.15', 'environment default_temperature_K'),
           ('gen_seed', str(SEED), 'environment random_seed'),
           ('dt', '0.01', 'stage 1: half the target timestep after minimisation (validation.md: smaller equilibration steps where needed)'),
           ('nsteps', '10000', 'stage 1: 100 ps at 10 fs')]
    eq1 = eq
    for k, v, n in eqa:
        eq1 = mdp_set(eq1, k, v, n)
        adapt.append({'mdp': 'eq1.mdp', 'key': k, 'value': v, 'reason': n})
    (td / 'eq1.mdp').write_text('\n'.join(eq1) + '\n')
    eq2 = eq
    for k, v, n in eqa[:4] + [('gen_vel', 'no', 'continue from stage 1 velocities'),
                               ('continuation', 'yes', 'continue from stage 1 (constraints not re-solved at start)'),
                               ('dt', '0.02', 'stage 2 at target timestep; tutorial value'),
                               ('nsteps', '25000', 'tutorial value (500 ps)')]:
        eq2 = mdp_set(eq2, k, v, n)
        adapt.append({'mdp': 'eq2.mdp', 'key': k, 'value': v, 'reason': n})
    (td / 'eq2.mdp').write_text('\n'.join(eq2) + '\n')
    pr = read_mdp(PROD_MDP)
    for k, v, n in [('nsteps', str(nsteps_smoke), f'smoke test {nsteps_smoke*0.02:.0f} ps at 20 fs'),
                    ('nstxout-compressed', '50', 'write a frame every 1 ps for geometry checks'),
                    ('compressed-x-precision', '1000', 'xtc precision 0.001 nm so constraint deviations are measurable (0.01 nm default rounding exceeds the check tolerance)'),
                    ('nstlog', '500', 'more frequent logging for a short run'),
                    ('nstenergy', '50', 'energy every 1 ps for averages'),
                    ('continuation', 'yes', 'continue from equilibration')]:
        pr = mdp_set(pr, k, v, n)
        adapt.append({'mdp': 'smoke.mdp', 'key': k, 'value': v, 'reason': n})
    (td / 'smoke.mdp').write_text('\n'.join(pr) + '\n')
    return adapt


def count_in(paths, needle):
    n = 0
    for p in paths:
        if Path(p).exists():
            n += Path(p).read_text(errors='replace').count(needle)
    return n


def grompp_notes(errpath):
    txt = Path(errpath).read_text(errors='replace')
    out = []
    for kind in ('NOTE', 'WARNING', 'ERROR'):
        i = 0
        while True:
            j = txt.find(kind + ' ', i)
            if j < 0:
                break
            seg = txt[j:j + 400].split('\n\n')[0]
            if seg.startswith(kind + ' ') and ('\n' in seg or ':' in seg[:20]):
                out.append({'kind': kind, 'text': ' '.join(seg.split())[:400]})
            i = j + 1
    return out


def parse_em(logpath):
    t = Path(logpath).read_text(errors='replace')
    res = {'converged': None, 'fmax': None, 'epot': None, 'steps': None, 'message': None}
    for line in t.splitlines():
        if 'converged to Fmax' in line or 'did not converge' in line or 'converged to machine precision' in line:
            res['message'] = line.strip()
            res['converged'] = ('converged to Fmax' in line)
            try:
                res['steps'] = int(line.split(' in ')[-1].split()[0])
            except Exception:
                pass
        if line.strip().startswith('Potential Energy'):
            res['epot'] = float(line.split('=')[1])
        if line.strip().startswith('Maximum force'):
            res['fmax'] = float(line.split('=')[1].split()[0])
    return res


def read_gro_frames(path, nat):
    frames, box = [], []
    with open(path) as f:
        while True:
            title = f.readline()
            if not title:
                break
            n = int(f.readline())
            xyz = []
            for i in range(n):
                l = f.readline()
                if i < nat:
                    xyz.append([float(l[20:28]), float(l[28:36]), float(l[36:44])])
            box.append([float(x) for x in f.readline().split()[:3]])
            frames.append(np.array(xyz))
    return frames, box


def xvg_avg(path):
    vals = []
    for l in open(path):
        if l[0] in '#@':
            continue
        vals.append([float(x) for x in l.split()])
    a = np.array(vals)
    return a


def cmd_test(a):
    md = Path(a.moldir)
    des = rjson(md / 'decisions/design.json')
    par = rjson(md / 'decisions/parameters.json')
    td = md / 'test'
    if td.exists() and (td / 'execution.json').exists():
        k = 1
        while (md / 'test_attempts' / f'attempt{k}').exists():
            k += 1
        (md / 'test_attempts').mkdir(exist_ok=True)
        shutil.move(str(td), str(md / 'test_attempts' / f'attempt{k}'))
        td.mkdir()
        shutil.copy(md / 'test_attempts' / f'attempt{k}' / 'static_checks.json', td / 'static_checks.json')
    td.mkdir(exist_ok=True)
    nb = len(des['beads'])
    resname = des['resname']
    nsm = int(ENV['smoke_test_duration_ps'] / ENV['target_test_timestep_ps'])
    adapt = make_mdps(td, nsm, 'FLEXIBLE' in (md / 'model/molecule.itp').read_text())
    wjson(td / 'mdp_adaptations.json', {'em_source': str((TUT_DIR / 'martini_em.mdp').relative_to(PROJECT)),
                                         'eq_source': str((TUT_DIR / 'martini_eq.mdp').relative_to(PROJECT)),
                                         'smoke_source': str(PROD_MDP.relative_to(PROJECT)), 'adaptations': adapt})
    shutil.copy(md / 'model/molecule.gro', td / 'solute.gro')
    shutil.copy(md / 'model/molecule.itp', td / 'molecule.itp')
    exe = {'steps': []}
    box = float((md / 'model/molecule.gro').read_text().splitlines()[-1].split()[0])
    rc, so, se = run_cmd(md, GMX + ['solvate', '-cp', 'solute.gro', '-cs', str(WATER_GRO), '-o', 'solvated.gro',
                                    '-box', str(box), str(box), str(box)], td, 'solvate')
    nW = sum(1 for l in open(td / 'solvated.gro').readlines()[2:-1] if l[5:10].strip() == 'W') if rc == 0 else 0
    top = ['; minimal test system: one solute in Martini 3.0.0 water',
           f'#include "{CORE_ITP}"', f'#include "{SOLV_ITP}"', f'#include "{IONS_ITP}"', '#include "molecule.itp"', '',
           '[ system ]', f'{resname} in Martini water', '', '[ molecules ]', f'{resname}  1', f'W  {nW}']
    (td / 'system.top').write_text('\n'.join(top) + '\n')
    exe['solvate'] = {'exit_code': rc, 'n_water_beads': nW, 'box_nm': box,
                      'solvent_coordinates': str(WATER_GRO.relative_to(PROJECT)), 'solvent_sha256': sha256(WATER_GRO)}
    # counterions: the temporary tpr is only an input to genion and is never simulated
    netq = sum(b.get('charge_e', 0.0) for b in des['beads'])
    start_coords = 'solvated.gro'
    if rc == 0 and abs(netq) > 1e-6:
        r1, _, e1 = run_cmd(md, GMX + ['grompp', '-f', 'em.mdp', '-c', 'solvated.gro', '-p', 'system.top', '-o', 'ions.tpr',
                                       '-po', 'ions_out.mdp', '-maxwarn', '1'], td, 'grompp_ions')
        r2 = 1
        if r1 == 0:
            r2, _, _ = run_cmd(md, GMX + ['genion', '-s', 'ions.tpr', '-o', 'solvated_ions.gro', '-p', 'system.top',
                                          '-pname', 'NA', '-nname', 'CL', '-neutral', '-seed', str(SEED)],
                               td, 'genion', stdin='W\n')
        rc = r1 or r2
        start_coords = 'solvated_ions.gro'
        exe['counterions'] = {'net_charge_e': netq, 'ion_names': ['NA', 'CL'], 'grompp_exit': r1, 'genion_exit': r2,
                              'grompp_maxwarn_note': 'maxwarn 1 only for the non-neutral ions.tpr, which is not run',
                              'ions_itp': str(IONS_ITP.relative_to(PROJECT)), 'ions_itp_sha256': sha256(IONS_ITP)}
    stages = [('em', 'em.mdp', start_coords, None), ('eq1', 'eq1.mdp', 'em.gro', None),
              ('eq2', 'eq2.mdp', 'eq1.gro', 'eq1.cpt'), ('smoke', 'smoke.mdp', 'eq2.gro', 'eq2.cpt')]
    status_ok = rc == 0
    for name, mdp, conf, cpt in stages:
        st = {'stage': name}
        if not status_ok:
            st['result'] = 'not_run'
            exe['steps'].append(st)
            continue
        args = GMX + ['grompp', '-f', mdp, '-c', conf, '-p', 'system.top', '-o', f'{name}.tpr', '-po', f'{name}_out.mdp',
                      '-pp', f'{name}_processed.top', '-maxwarn', '0']
        if cpt:
            args += ['-t', cpt]
        rc1, so1, se1 = run_cmd(md, args, td, f'grompp_{name}')
        st['grompp_exit'] = rc1
        st['grompp_messages'] = grompp_notes(se1)
        if rc1 != 0:
            st['result'] = 'fail'
            status_ok = False
            exe['steps'].append(st)
            continue
        rc2, so2, se2 = run_cmd(md, GMX + ['mdrun', '-deffnm', name, '-ntmpi', '1', '-ntomp', str(os.environ.get('MARTINI_GMX_THREADS', 4)), '-nb', 'cpu', '-pin', 'off'],
                                td, f'mdrun_{name}')
        st['mdrun_exit'] = rc2
        logf = td / f'{name}.log'
        st['lincs_warnings'] = count_in([logf, se2], 'LINCS WARNING')
        st['finished'] = logf.exists() and 'Finished mdrun' in logf.read_text(errors='replace')
        if name == 'em':
            st['em'] = parse_em(logf)
            st['em']['emtol_kJ_mol_nm'] = 10.0
            st['em']['nsteps_limit'] = 2000
        if rc2 != 0 or not st['finished'] or not (td / f'{name}.gro').exists():
            st['result'] = 'fail'
            status_ok = False
        else:
            st['result'] = 'pass'
        exe['steps'].append(st)
    # --- analysis of smoke trajectory
    geo = None
    thermo = None
    if status_ok:
        rc, _, _ = run_cmd(md, GMX + ['trjconv', '-f', 'smoke.xtc', '-s', 'smoke.tpr', '-pbc', 'mol', '-o', 'smoke_frames.gro'],
                           td, 'trjconv_smoke', stdin='0\n')
        frames, _ = read_gro_frames(td / 'smoke_frames.gro', nb)
        geo = []
        for p in par['terms']:
            if p.get('section') not in ('bonds', 'constraints', 'angles', 'dihedrals'):
                continue
            ids = [i - 1 for i in p['bead_ids']]
            vals = [measure(p['section'], p['function'], [f[i] for i in ids]) for f in frames]
            ref = None if p['function'] == 9 else p['parameters'][0]['value']
            if ref is None:
                m, sd = circ_mean(vals)
                dev = None
            elif p['section'] == 'dihedrals':
                m, sd = circ_mean(vals)
                dev = max(abs((v - ref + 180) % 360 - 180) for v in vals)
            else:
                m, sd = float(np.mean(vals)), float(np.std(vals))
                dev = max(abs(v - ref) for v in vals)
            geo.append({'term_id': p['term_id'], 'section': p['section'], 'function': p['function'], 'reference': ref,
                        'mean': round(m, 4), 'sd': round(sd, 4), 'min': round(min(vals), 4), 'max': round(max(vals), 4),
                        'max_abs_dev_from_reference': None if dev is None else round(dev, 4), 'n_frames': len(vals),
                        'ensemble': p['ensemble']})
        finite = all(np.isfinite(f).all() for f in frames)
        thermo = {'finite_coordinates': finite, 'n_frames': len(frames)}
        rc, _, _ = run_cmd(md, GMX + ['energy', '-f', 'smoke.edr', '-o', 'smoke_energy.xvg'], td, 'energy_smoke',
                           stdin='Potential\nTemperature\nPressure\nDensity\n\n')
        if rc == 0:
            arr = xvg_avg(td / 'smoke_energy.xvg')
            names = ['potential_kJ_mol', 'temperature_K', 'pressure_bar', 'density_kg_m3']
            for k, nm in enumerate(names):
                col = arr[:, k + 1]
                thermo[nm] = {'mean': round(float(col.mean()), 3), 'sd': round(float(col.std()), 3),
                              'finite': bool(np.isfinite(col).all())}
    exe['geometry_smoke'] = geo
    exe['thermo_smoke'] = thermo
    wjson(td / 'execution.json', exe)
    print('test', des['resname'], 'ok' if status_ok else 'FAILED', [(s['stage'], s.get('result')) for s in exe['steps']])


# ---------------------------------------------------------------- finalize
def file_hashes(root, sub):
    out = {}
    for p in sorted((root / sub).rglob('*')):
        if p.is_file():
            out[str(p.relative_to(root))] = sha256(p)
    return out


def skill_hash():
    """Tree hash of the skill this script belongs to, over source files only (no __pycache__), so it can be
    checked against config/skill.*.lock.json wherever the skill is installed."""
    sd = SCRIPT_DIR.parent
    h = hashlib.sha256()
    for p in sorted(sd.rglob('*')):
        if p.is_file() and '__pycache__' not in p.parts:
            h.update(str(p.relative_to(sd)).encode())
            h.update(sha256(p).encode())
    return h.hexdigest(), sha256(sd / 'SKILL.md')


def verify_reference_lock(lock):
    p = subprocess.run([sys.executable, str(SCRIPT_DIR / 'reference_lock.py'), 'verify', str(PROJECT / ENV['reference_root']), str(lock)],
                       capture_output=True, text=True)
    return {'result': 'pass' if p.returncode == 0 else 'fail', 'detail': (p.stdout + p.stderr).strip()[-300:]}


def gmx_version():
    p = subprocess.run(GMX + ['--version'], capture_output=True, text=True)
    for l in (p.stdout + p.stderr).splitlines():
        if l.startswith('GROMACS version'):
            return l.split(':', 1)[1].strip()
    return 'unknown'


def classify(md, des, par, exe, static):
    checks = []
    for s in static['checks']:
        checks.append({'level': 'static', 'check': s['check'], 'result': s['result'], 'evidence': 'test/static_checks.json',
                       'value': s['detail']})
    steps = {s['stage']: s for s in exe.get('steps', [])}
    for st in ('em', 'eq1', 'eq2', 'smoke'):
        s = steps.get(st, {})
        res = s.get('result', 'not_run')
        checks.append({'level': 'execution', 'check': f'{st}_completed', 'result': res if res in ('pass', 'fail') else 'not_run',
                       'evidence': f'test/{st}.log', 'value': {k: s.get(k) for k in ('grompp_exit', 'mdrun_exit', 'finished', 'lincs_warnings')}})
        if s.get('grompp_messages') is not None:
            warn = [m for m in s['grompp_messages'] if m['kind'] in ('WARNING', 'ERROR')]
            checks.append({'level': 'execution', 'check': f'{st}_grompp_no_warnings', 'result': 'pass' if not warn else 'fail',
                           'evidence': f'logs/grompp_{st}.stderr', 'value': s['grompp_messages']})
        if s.get('lincs_warnings') is not None:
            checks.append({'level': 'execution', 'check': f'{st}_no_lincs_warnings', 'result': 'pass' if s['lincs_warnings'] == 0 else 'fail',
                           'evidence': f'test/{st}.log', 'value': s['lincs_warnings']})
    em = steps.get('em', {}).get('em')
    if em:
        checks.append({'level': 'execution', 'check': 'em_converged_to_emtol', 'result': 'pass' if em['converged'] else 'fail',
                       'evidence': 'test/em.log', 'value': em})
    geo = exe.get('geometry_smoke')
    if geo:
        cons = [g for g in geo if g['section'] == 'constraints']
        cdev = max([g['max_abs_dev_from_reference'] for g in cons], default=0.0)
        checks.append({'level': 'execution', 'check': 'constraint_max_abs_deviation_below_0.005nm_in_smoke',
                       'result': 'pass' if cdev < 0.005 else 'fail', 'evidence': 'test/execution.json', 'value': cdev})
        bonds = [g for g in geo if g['section'] == 'bonds']
        bdev = max([abs(g['mean'] - g['reference']) for g in bonds], default=0.0)
        checks.append({'level': 'execution', 'check': 'bond_mean_within_0.02nm_of_b0_in_smoke',
                       'result': 'pass' if bdev < 0.02 else 'fail', 'evidence': 'test/execution.json', 'value': round(bdev, 4)})
        dih = [g for g in geo if g['section'] == 'dihedrals' and g['function'] == 2]
        ddev = max([abs((g['mean'] - g['reference'] + 180) % 360 - 180) for g in dih], default=0.0)
        checks.append({'level': 'execution', 'check': 'improper_mean_within_20deg_of_phi0_in_smoke (planarity kept)',
                       'result': 'pass' if ddev < 20 else 'fail', 'evidence': 'test/execution.json', 'value': round(ddev, 2)})
        dih_trip = set()
        for t in par['terms']:
            if t.get('section') == 'dihedrals':
                q = t['bead_ids']
                dih_trip |= {tuple(q[:3]), tuple(q[1:]), tuple(q[:3][::-1]), tuple(q[1:][::-1])}
        ang = [g for g in geo if g['section'] == 'angles']
        risky = [g for g in ang if tuple(next(t['bead_ids'] for t in par['terms'] if t['term_id'] == g['term_id'])) in dih_trip]
        amax = max([g['max'] for g in risky], default=0.0)
        checks.append({'level': 'execution', 'check': 'dihedral_spanning_angles_never_reached_178deg_in_smoke (SI C4 p.46)', 'result': 'pass' if amax < 178 else 'fail',
                       'evidence': 'test/execution.json', 'value': {'max_deg': round(amax, 2), 'n_angles_checked': len(risky),
                                                                    'all_angles_max_deg_info': {g['term_id']: g['max'] for g in ang}}})
    th = exe.get('thermo_smoke')
    if th:
        fin = th['finite_coordinates'] and all(v.get('finite', True) for k, v in th.items() if isinstance(v, dict))
        checks.append({'level': 'execution', 'check': 'finite_energies_and_coordinates_in_smoke', 'result': 'pass' if fin else 'fail',
                       'evidence': 'test/smoke_energy.xvg', 'value': th})
    unresolved = []
    static_fail = [c for c in checks if c['level'] == 'static' and c['result'] == 'fail']
    run_ok = all(steps.get(s, {}).get('result') == 'pass' for s in ('em', 'eq1', 'eq2', 'smoke'))
    pre_ok = steps.get('em', {}).get('grompp_exit') == 0
    hard_exec = [c for c in checks if c['level'] == 'execution' and c['result'] == 'fail'
                 and c['check'] not in ('em_converged_to_emtol',)]
    if unresolved or static_fail:
        status = 'blocked'
    elif run_ok and not hard_exec:
        status = 'execution_checked'
    elif pre_ok:
        status = 'preprocessed'
    else:
        status = 'constructed_unexecuted'
    return checks, status, unresolved


def cmd_finalize(a):
    md = Path(a.moldir)
    des = rjson(md / 'decisions/design.json')
    req = rjson(md / 'input/request.json')
    par = rjson(md / 'decisions/parameters.json')
    mp = rjson(md / 'decisions/mapping.json')
    ana = rjson(md / 'decisions/analogue_search.json')
    bonded = rjson(md / 'decisions/bonded.json')
    static = rjson(md / 'test/static_checks.json')
    exe = rjson(md / 'test/execution.json') if (md / 'test/execution.json').exists() else {}
    checks, status, unresolved = classify(md, des, par, exe, static)
    steps = {s['stage']: s for s in exe.get('steps', [])}
    tested = []
    for st, dt, ns in (('eq1', 0.01, 10000), ('eq2', 0.02, 25000), ('smoke', 0.02, int(ENV['smoke_test_duration_ps'] / 0.02))):
        if steps.get(st, {}).get('result') == 'pass':
            tested.append({'stage': st, 'timestep_ps': dt, 'nsteps': ns, 'duration_ps': round(dt * ns, 1)})
    provisional_types = [b['name'] + ':' + b['type'] for b in mp['beads'] if b['provisional_type']]
    ev = {}
    for p in par['terms']:
        if p.get('section') in ('bonds', 'constraints', 'angles', 'dihedrals'):
            ev.setdefault(p['evidence_class'], []).append(p['term_id'])
    phys = {'observables': [], 'note': 'No physical evaluation (partitioning, density, SASA vs AA) was run in this generation pass. Smoke-test thermodynamic averages are numerical diagnostics only.'}
    if exe.get('thermo_smoke'):
        phys['smoke_diagnostics_not_validation'] = exe['thermo_smoke']
    val = {'schema_version': 1, 'molecule_id': req['molecule_id'], 'status': status,
           'status_definition': 'validation.md: blocked | constructed_unexecuted | preprocessed | execution_checked',
           'checks': checks, 'unresolved_terms': unresolved, 'provisional_types': provisional_types,
           'tested_dynamics': tested, 'outstanding_warnings': [m for s in exe.get('steps', []) for m in s.get('grompp_messages', []) if m['kind'] != 'ERROR'],
           'physical_evaluation': phys}
    wjson(md / 'validation.json', val)
    sk, skmd = skill_hash()
    lock = PROJECT / ENV['reference_lock']
    man = {
        'schema_version': 1, 'molecule_id': req['molecule_id'], 'created_utc': now(),
        'reference_lock': str(lock.relative_to(PROJECT)), 'reference_lock_sha256': sha256(lock),
        'reference_lock_verify': verify_reference_lock(lock),
        'dependencies': {str(p.relative_to(PROJECT)): sha256(p) for p in (CORE_ITP, SOLV_ITP, WATER_GRO, PROD_MDP, INDEX, BIBLE_XLSX,
                                                                          TUT_DIR / 'martini_em.mdp', TUT_DIR / 'martini_eq.mdp')},
        'skill': {'path': SCRIPT_DIR.parent.name, 'tree_sha256': sk, 'SKILL.md_sha256': skmd},
        'software': {'python': sys.executable, 'python_version': sys.version.split()[0], 'rdkit': rdBase.rdkitVersion,
                     'numpy': np.__version__, 'gromacs': gmx_version(), 'gromacs_command': GMX,
                     'bonded_stage': bonded['software']},
        'model_identifier': des.get('authored_by', 'unknown'),
        'seeds': {'conformer_embedding': SEED, 'gen_seed_eq1': SEED},
        'hashes': {'request_json': sha256(md / 'input/request.json'), 'environment_json': sha256(ENV_PATH),
                   'design_json': sha256(md / 'decisions/design.json'), 'bonded_json': sha256(md / 'decisions/bonded.json'),
                   'cgbuild_py': sha256(SCRIPT_DIR / 'cgbuild.py'), 'derive_bonded_py': sha256(SCRIPT_DIR / 'derive_bonded.py')},
        'revisions': des.get('revisions', []),
        'runtime': {'commands_wall_s': round(sum(json.loads(l)['wall_s'] for l in open(md / 'logs/commands.jsonl')), 1)
                    if (md / 'logs/commands.jsonl').exists() else None},
        'reproduce': [f"{ENV['python_command'][0]} scripts/cgbuild.py input . --smiles '{req['original_smiles']}' --id {req['molecule_id']} --name '{req['name']}'",
                      f"{ENV['python_command'][0]} scripts/cgbuild.py search .", f"{ENV['python_command'][0]} scripts/cgbuild.py bonded .",
                      f"{ENV['python_command'][0]} scripts/cgbuild.py build .",
                      f"{ENV['python_command'][0]} scripts/cgbuild.py test .", f"{ENV['python_command'][0]} scripts/cgbuild.py finalize .",
                      'decisions/design.json is the agent-authored decision record and is an input, not regenerated'],
        'flags': {'reference_exposed': ana['reference_exposed'], 'reference_exposed_evidence': ana['exact_target_matches'] + ana.get('tautomer_target_matches', []),
                  'assisted_refinement': False, 'comparison_assisted': False, 'automatic_cg_generators_used': False,
                  'pretraining_contamination': 'unknown'},
        'manual_interventions': des.get('manual_interventions', []),
        'unverified_statements': des.get('unverified_statements', []),
    }
    wjson(md / 'manifest.json', man)
    write_report(md, req, des, mp, par, ana, val, man, exe, bonded)
    # the code that produced this model, kept with the run
    if (md / 'scripts').exists():
        shutil.rmtree(md / 'scripts')
    shutil.copytree(SCRIPT_DIR, md / 'scripts', ignore=shutil.ignore_patterns('__pycache__'))
    # baseline freeze
    bl = md / 'baseline'
    if bl.exists():
        shutil.rmtree(bl)
    for sub in ('input', 'decisions', 'model'):
        shutil.copytree(md / sub, bl / sub)
    for f in ('validation.json', 'manifest.json', 'report.md'):
        shutil.copy(md / f, bl / f)
    hashes = file_hashes(md, 'baseline')
    wjson(md / 'baseline_hashes.json', {'frozen_utc': now(), 'files': hashes})
    print('finalize', req['molecule_id'], status, 'provisional', provisional_types)


def fmt_params(p):
    out = []
    for x in p['parameters']:
        out.append(f"{x['name']}={x['value']} {x['unit']}")
    return ', '.join(out)


def write_report(md, req, des, mp, par, ana, val, man, exe, bonded):
    R = []
    R.append(f"# {req['name']} ({req['molecule_id']}) — Martini 3 candidate model")
    R.append('')
    R.append(f"Status: **{val['status']}** (from validation.json). Reference exposure: **{man['flags']['reference_exposed']}**.")
    R.append('')
    R.append(f"- SMILES: `{req['original_smiles']}` (canonical `{req['canonical_isomeric_smiles']}`), formula {req['formula']}, {req['heavy_atoms']} heavy atoms")
    R.append(f"- Net charge: {mp['net_charge_e']} e (formal {req['formal_charge']})")
    R.append(f"- Unspecified stereocentres (atom ids): {req['unspecified_stereocentres'] or 'none'}; double-bond stereo: {req['double_bond_stereo'] or 'none'}")
    R.append(f"- Beads: {mp['n_beads']} ({', '.join(b['name']+'='+b['type'] for b in mp['beads'])}); "
             f"virtual sites: {sum(1 for b in mp['beads'] if b['virtual_site'])}")
    R.append(f"- Provisional types: {', '.join(val['provisional_types']) or 'none'}")
    R.append(f"- Unresolved terms: {', '.join(val['unresolved_terms']) or 'none'}")
    used = [h for h in ana['hits'] if h['used']]
    R.append(f"- Analogue search: {ana['n_hits']} hits, {len(used)} used: " + ', '.join(h['reference_id'] for h in used))
    R.append('')
    R.append('## Mapping')
    R.append('')
    R.append(f"Chosen candidate {mp['chosen_candidate']}. Drawing: `model/mapping.svg`.")
    R.append('')
    R.append('| Bead | Type | Atoms (id:weight) | Mass | Ladder rung | Provisional | Rationale |')
    R.append('|---|---|---|---|---|---|---|')
    for b in mp['beads']:
        atoms = ' '.join(f"{m['element']}{m['atom_id']}" + ('' if m['allocation_weight'] == 1 else f":{m['allocation_weight']}") for m in b['members'])
        R.append(f"| {b['name']} | {b['type']} | {atoms} | {b['mass_amu']} | {b['typing_ladder']['rung_used']} | {b['provisional_type']} | {b['rationale']} |")
    R.append('')
    R.append('## Bonded terms')
    R.append('')
    en, em = bonded['ensemble'], bonded['energy_model']
    R.append(f"Derived without MD and without agent decisions. Energy model: {em['name']}. "
             f"{en['n_conf']} conformers, {en['n_used']} minima used, effective minima {en['n_eff_minima']:.1f}, "
             f"top weight {en['top_weight']:.2f}, {bonded['temperature_K']} K. Engines tried: "
             + ', '.join(f"{x['engine']} ({x['result']})" for x in em['engines_tried']) + '.')
    R.append(f"Rigid groups: {bonded['rigid_groups'] or 'none'}. Terms considered and dropped: {len(bonded['dropped'])} (reasons in `decisions/bonded.json`).")
    limited = [t for t in bonded['terms'] if t.get('fitted_before_limit') and t['section'] != 'constraints']
    if limited:
        R.append('Terms changed from their fitted value: ' + '; '.join(f"{t['term_id']} ({t['limitations']}, fitted {t['fitted_before_limit']})" for t in limited) + '.')
    for f in bonded['flags']:
        R.append(f"**Flag `{f['flag']}`:** {f.get('note', '')} {({k: v for k, v in f.items() if k not in ('flag', 'note')})}")
    R.append('')
    R.append('| Term | Section | Beads | Funct | Parameters | Evidence | Source |')
    R.append('|---|---|---|---|---|---|---|')
    for p in par['terms']:
        if p.get('section') in ('bonds', 'constraints', 'angles', 'dihedrals'):
            src = p['source'] if isinstance(p['source'], str) else '; '.join(f"{s.get('path','')} {s.get('locator','')}".strip() for s in p['source'])
            R.append(f"| {p['term_id']} | {p['section']} | {'-'.join(p['bead_names'])} | {p['function']} | {fmt_params(p)} | {p['evidence_class']} | {src} |")
    R.append('')
    ex = [p for p in par['terms'] if p.get('section') == 'exclusions']
    if ex:
        R.append('Exclusions: ' + '; '.join(f"{'-'.join(map(str, e['bead_ids']))} ({e['evidence_class']})" for e in ex) + f". nrexcl = {par['nrexcl']}.")
    else:
        R.append(f'Exclusions: none beyond nrexcl = {par["nrexcl"]}.')
    R.append('')
    R.append('## Checks')
    R.append('')
    for lvl in ('static', 'execution'):
        cs = [c for c in val['checks'] if c['level'] == lvl]
        R.append(f"**{lvl.capitalize()}** ({sum(1 for c in cs if c['result']=='pass')}/{len(cs)} pass)")
        R.append('')
        for c in cs:
            v = c['value']
            vs = json.dumps(v) if not isinstance(v, (int, float)) else str(v)
            if len(vs) > 160:
                vs = vs[:157] + '...'
            R.append(f"- {c['check']}: {c['result']} — `{vs}`")
        R.append('')
    R.append('Tested dynamics: ' + ('; '.join(f"{t['stage']} {t['duration_ps']} ps at {t['timestep_ps']*1000:.0f} fs" for t in val['tested_dynamics']) or 'none'))
    R.append('')
    R.append('Physical evaluation: ' + val['physical_evaluation']['note'])
    R.append('')
    R.append('## Unverified statements ("agent background knowledge, unverified")')
    R.append('')
    if man['unverified_statements']:
        for u in man['unverified_statements']:
            R.append(f"- {u['text']} — at {u['location']} — {u['later_status']}")
    else:
        R.append('- none')
    R.append('')
    R.append('## Deviations and limitations')
    R.append('')
    for d in des.get('limitations', []):
        R.append(f'- {d}')
    R.append('')
    R.append('Files: `model/molecule.itp`, `model/molecule.gro`, `model/mapping.svg`, `decisions/*.json`, `test/` (system.top, MDPs, logs, trajectories), `manifest.json`, `validation.json`, `baseline/` (frozen, hashes in `baseline_hashes.json`).')
    (md / 'report.md').write_text('\n'.join(R) + '\n')


# ---------------------------------------------------------------- summary
def cmd_summary(a):
    bd = Path(a.batchdir)
    rows = []
    for md in sorted(p for p in bd.iterdir() if p.is_dir() and (p / 'validation.json').exists()):
        val = rjson(md / 'validation.json')
        mp = rjson(md / 'decisions/mapping.json')
        man = rjson(md / 'manifest.json')
        par = rjson(md / 'decisions/parameters.json')
        ana = rjson(md / 'decisions/analogue_search.json')
        req = rjson(md / 'input/request.json')
        ev = {}
        for p in par['terms']:
            if p.get('section') in ('bonds', 'constraints', 'angles', 'dihedrals'):
                ev[p['evidence_class']] = ev.get(p['evidence_class'], 0) + 1
        sm = {c['check']: c for c in val['checks']}
        rows.append({'molecule_id': req['molecule_id'], 'name': req['name'], 'smiles': req['original_smiles'],
                     'status': val['status'], 'n_beads': mp['n_beads'], 'bead_types': [b['type'] for b in mp['beads']],
                     'net_charge_e': mp['net_charge_e'], 'provisional_types': val['provisional_types'],
                     'unresolved_terms': val['unresolved_terms'], 'term_evidence_counts': ev,
                     'analogues_used': [h['reference_id'] for h in ana['hits'] if h['used']],
                     'reference_exposed': man['flags']['reference_exposed'],
                     'em_converged': sm.get('em_converged_to_emtol', {}).get('result'),
                     'tested_dynamics': val['tested_dynamics'],
                     'n_unverified_statements': len(man['unverified_statements']),
                     'failed_checks': [c['check'] for c in val['checks'] if c['result'] == 'fail']})
    wjson(bd / 'summary.json', {'schema_version': 1, 'batch': bd.name, 'created_utc': now(), 'molecules': rows,
                                'status_counts': {s: sum(1 for r in rows if r['status'] == s) for s in sorted({r['status'] for r in rows})}})
    R = [f'# Batch {bd.name}', '', 'Built only from each molecule\'s records (validation.json, mapping.json, parameters.json, analogue_search.json, manifest.json).', '',
         '| ID | Name | Status | Beads | Types | Charge | Provisional types | Unresolved | Term evidence | EM converged | Failed checks | Exposed |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        R.append(f"| {r['molecule_id']} | {r['name']} | {r['status']} | {r['n_beads']} | {' '.join(r['bead_types'])} | {r['net_charge_e']} | "
                 f"{', '.join(r['provisional_types']) or '-'} | {', '.join(r['unresolved_terms']) or '-'} | "
                 f"{', '.join(f'{k}:{v}' for k, v in sorted(r['term_evidence_counts'].items()))} | {r['em_converged']} | {', '.join(r['failed_checks']) or '-'} | {r['reference_exposed']} |")
    R.append('')
    R.append('Per-molecule details: `<ID>/report.md`.')
    (bd / 'report.md').write_text('\n'.join(R) + '\n')
    print('summary', {s: sum(1 for r in rows if r['status'] == s) for s in {r['status'] for r in rows}})


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('input')
    p.add_argument('moldir')
    p.add_argument('--smiles', required=True)
    p.add_argument('--id', required=True)
    p.add_argument('--name', required=True)
    p.add_argument('--class', dest='cls', default='')
    p.add_argument('--stereo-note', default='')
    for c in ('search', 'bonded', 'build', 'test', 'finalize'):
        q = sp.add_parser(c)
        q.add_argument('moldir')
    q = sp.add_parser('summary')
    q.add_argument('batchdir')
    a = ap.parse_args()
    {'input': cmd_input, 'search': cmd_search, 'bonded': cmd_bonded, 'build': cmd_build, 'test': cmd_test,
     'finalize': cmd_finalize, 'summary': cmd_summary}[a.cmd](a)


if __name__ == '__main__':
    main()
