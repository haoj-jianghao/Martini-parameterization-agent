#!/usr/bin/env python3
"""Write decisions/design.json for ONE molecule, one decision at a time (SKILL.md stages 2-3).

The agent makes every chemical decision; this tool records it, supplies the evidence the decision needs
and checks that the record says only what the evidence supports:
- `init` derives the index queries from the molecule itself (its ring systems, ring-substituent junctions,
  functional groups with their neighbouring carbons, aliphatic chains) and runs the search;
- `open` shows a reference model bead by bead (type, member atoms, fragment signature) and records that it
  was opened;
- `match` lists every indexed bead that holds exactly the same fragment as a proposed bead;
- `bead --rung 1 --from REF:BEAD` copies the type from that source bead after verifying that it holds the
  same atoms in the same bonding context; rung 2/3 record the named group, rule or table, and the rung-1
  precedents the agent decided not to use, with the reason.
Every command is appended to decisions/design_log.jsonl with the hash of the design it leaves, so
`cgbuild.py build` can tell whether design.json was written by this tool.

Run with the environment's python_command. One molecule folder per agent session: when MARTINI_MOLDIR is
set (scripts/batch.py sets it), every command must target that folder.

  design.py init      <moldir> [--agent A] [--model M]
  design.py atoms     <moldir>
  design.py hits      <moldir> [--query NAME] [--limit N]
  design.py query     <moldir> --name N --smarts S [--tags a,b] [--no-bible]
  design.py open      <moldir> REF_ID
  design.py annotate  <moldir> REF_ID --reason TEXT [--used-for TEXT]
  design.py match     <moldir> --atoms 1,2,3[:w]
  design.py candidate <moldir> --id ID --beads "1,2,3|4,5:0.5,6" --summary S --pros "a;b" --cons "c;d"
  design.py choose    <moldir> --candidate ID --rationale TEXT
  design.py bead      <moldir> --name B1 --atoms 1,2,3 --rung 1 --from REF_ID:BEAD --rationale TEXT
  design.py bead      <moldir> --name B1 --atoms 1,2,3 --rung 2 --type T --group G --rule R --source PATH:LOCATOR
                      [--provisional] [--rung1-rejection TEXT] --rationale TEXT
  design.py bead      <moldir> --name B1 --atoms ... --rung 3 --type T --table TABLE --values V --rung2-finding TEXT
                      --source PATH:LOCATOR [--provisional] [--rung1-rejection TEXT] --rationale TEXT
                      (all rungs: [--charge Q] [--size-exception TEXT] [--replace])
  design.py vsite     <moldir> --site N --construction virtual_sitesn_cog|virtual_sites3 --constructing 1,2,3
                      --source PATH:LOCATOR --evidence-class C
  design.py hydrogens <moldir> --excluded-parents 3,4 --reason TEXT
  design.py note      <moldir> (--unverified TEXT | --limitation TEXT | --manual TEXT)
  design.py show      <moldir>
  design.py check     <moldir>
"""
import argparse, datetime, hashlib, json, os, sys
from pathlib import Path
from types import SimpleNamespace

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import cgbuild as C  # noqa: E402  (shared paths, index search, core types)
import refbeads as R  # noqa: E402

DESIGN = 'decisions/design.json'
LOG = 'decisions/design_log.jsonl'
REFCACHE = 'decisions/reference_beads.json'
_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

# Functional groups searched in every target; each matched group becomes a query made of the group and its
# neighbouring heavy atoms in the target (so a hydroxyl on a CH2 queries C-C-O, not every O).
GROUPS = [
    ('carboxylic_acid', '[CX3](=O)[OX2H1,OX1-]'), ('ester', '[#6][CX3](=O)[OX2][#6]'), ('amide', '[CX3](=O)[NX3]'),
    ('aldehyde', '[CX3H1](=O)'), ('ketone', '[#6][CX3](=O)[#6]'), ('hydroxyl', '[OX2H1][#6]'),
    ('ether', '[#6][OX2][#6]'), ('amine', '[NX3;!$(N-C=O);!$(N-a);!$([N]=*)]'), ('aniline_n', '[NX3;$(N-a)]'),
    ('ammonium', '[NX4+]'), ('nitrile', 'C#N'), ('nitro', '[N+](=O)[O-]'), ('halide', '[#6][F,Cl,Br,I]'),
    ('thiol', '[SX2H1]'), ('thioether', '[#6][SX2][#6]'), ('sulfonyl', '[SX4](=O)(=O)'), ('phosphate', 'P(=O)(O)O'),
    ('carbonate', 'O=C(O)O'), ('imine', '[CX3]=[NX2]'), ('alkene', '[CX3;!a]=[CX3;!a]'), ('alkyne', 'C#C'),
]


# ------------------------------------------------------------------ helpers

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def die(msg):
    print('design.py: ' + msg, file=sys.stderr)
    sys.exit(2)


def session():
    return {k: os.environ.get(k) for k in ('MARTINI_BATCH', 'MARTINI_MOLDIR', 'MARTINI_AGENT', 'MARTINI_MODEL')
            if os.environ.get(k)}


def bound_check(md):
    b = os.environ.get('MARTINI_MOLDIR')
    if b and Path(b).resolve() != md.resolve():
        die(f'this agent session is bound to {b}; it may not write decisions for {md}. '
            'One molecule per session (SKILL.md "Rules"): other molecules get their own session via scripts/batch.py.')


def load(md):
    p = md / DESIGN
    if not p.exists():
        die(f'{p} does not exist; run `design.py init {md}` first')
    return json.loads(p.read_text())


def save(md, des, command, args):
    C.wjson(md / DESIGN, des)
    args = {k: v for k, v in args.items() if k != 'fn'}
    rec = {'utc': now(), 'command': command, 'args': args, 'session': session(), 'design_sha256': sha(md / DESIGN)}
    with open(md / LOG, 'a') as f:
        f.write(json.dumps(rec) + '\n')


def target(md):
    mol = Chem.MolFromMolFile(str(md / 'input/atomistic.sdf'), removeHs=False)
    if mol is None:
        die('cannot read input/atomistic.sdf; run cgbuild.py input first')
    return mol


def heavy_ids(mol):
    return [a.GetIdx() + 1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1]


def parse_atoms(spec, mol):
    """'1,2,3:0.5' -> {1: 1.0, 2: 1.0, 3: 0.5} (one-based heavy-atom IDs)."""
    out = {}
    hv = set(heavy_ids(mol))
    for tok in [t for t in spec.replace(' ', '').split(',') if t]:
        i, _, w = tok.partition(':')
        i = int(i)
        if i not in hv:
            die(f'atom {i} is not a heavy atom of this molecule (see `design.py atoms`)')
        out[i] = float(w) if w else 1.0
    if not out:
        die('no atoms given')
    return out


def frag_smiles(mol, ids):
    return Chem.MolFragmentToSmiles(mol, [i - 1 for i in ids], canonical=True)


def signature(mol, members):
    return R.fragment_signature(mol, {i - 1: w for i, w in members.items()})


def ref_beads(md):
    """Reference beads with fragment signatures, cached per molecule folder (keyed by index and Bible hashes)."""
    key = {'index_sha256': C.sha256(C.INDEX), 'bible_sha256': C.sha256(C.BIBLE_XLSX),
           'refbeads_sha256': C.sha256(SCRIPT_DIR / 'refbeads.py')}
    p = md / REFCACHE
    if p.exists():
        d = json.loads(p.read_text())
        if d.get('key') == key:
            return d['references']
    refs = R.reference_beads(C.INDEX, C.bible_rows(), C.core_types(), C.PROJECT)
    C.wjson(p, {'key': key, 'note': 'written by design.py: bead membership and fragment signatures of indexed models '
                                    '(scripts/refbeads.py); entries with `error` cannot support rung 1', 'references': refs})
    return refs


def precedents(md, sig):
    out = []
    for rid, r in ref_beads(md).items():
        for b in r.get('beads', []):
            if b['signature'] == sig:
                out.append({'reference_id': rid, 'bead': b['name'], 'type': b['type'], 'name': r['name']})
    return out


def parse_source(s):
    path, _, loc = s.partition(':')
    if not (C.PROJECT / path).exists():
        die(f'source path {path} does not exist under {C.PROJECT}')
    if not loc.strip():
        die(f'source {s!r} needs a locator after the colon (section, page, table, moleculetype or line)')
    return {'path': path, 'locator': loc.strip()}


# ------------------------------------------------------------------ query derivation

def derive_queries(mol):
    m = Chem.RemoveHs(mol)
    qs, seen = [], set()

    def add(name, atoms, tags, origin):
        s = Chem.MolFragmentToSmiles(m, sorted(atoms), canonical=True)
        if s in seen or Chem.MolFromSmarts(s) is None:
            return
        seen.add(s)
        qs.append({'name': name, 'smarts': s, 'motif_tags': tags, 'search_bible': True, 'origin': origin,
                   'target_atoms': sorted(a + 1 for a in atoms)})

    ri = m.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    systems = []
    for r in rings:
        cur = set(r)
        for s in [s for s in systems if s & cur]:
            cur |= s
            systems.remove(s)
        systems.append(cur)
    for k, sys_atoms in enumerate(systems):
        add(f'ring_system_{k + 1}', sys_atoms, [], 'derived: ring system')
        if sum(1 for r in rings if r <= sys_atoms) > 1:
            for j, r in enumerate([r for r in rings if r <= sys_atoms]):
                add(f'ring_system_{k + 1}_ring_{j + 1}', r, [], 'derived: single ring of a fused system')
        for a in sys_atoms:
            for n in m.GetAtomWithIdx(a).GetNeighbors():
                if n.GetIdx() in sys_atoms:
                    continue
                sub = {n.GetIdx()} | {x.GetIdx() for x in n.GetNeighbors() if x.GetIdx() not in sys_atoms}
                add(f'ring_{k + 1}_substituent_at_{a + 1}', sys_atoms | sub, [], 'derived: ring-substituent junction')
    for name, smarts in GROUPS:
        patt = Chem.MolFromSmarts(smarts)
        matches = m.GetSubstructMatches(patt)
        if matches and smarts not in seen:
            seen.add(smarts)
            qs.append({'name': name, 'smarts': smarts, 'motif_tags': [], 'search_bible': True,
                       'origin': 'derived: functional group (library SMARTS)',
                       'target_atoms': sorted({a + 1 for t in matches for a in t})})
        for t in matches:
            env = set(t)
            for a in t:
                env |= {n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors()}
            add(f'{name}_at_{"-".join(str(a + 1) for a in t)}', env, [], f'derived: functional group {name} with neighbouring atoms')
    # aliphatic chains: maximal runs of acyclic sp3 carbons, 3 atoms or more (capped at 6 atoms)
    chain = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == 'C' and not a.IsInRing()
             and a.GetHybridization() == Chem.HybridizationType.SP3]
    cs, left = [], set(chain)
    while left:
        st = [left.pop()]
        comp = set(st)
        while st:
            x = st.pop()
            for n in m.GetAtomWithIdx(x).GetNeighbors():
                if n.GetIdx() in left:
                    left.discard(n.GetIdx())
                    comp.add(n.GetIdx())
                    st.append(n.GetIdx())
        cs.append(comp)
    for k, comp in enumerate(cs):
        if len(comp) >= 3:
            add(f'aliphatic_chain_{k + 1}', set(sorted(comp)[:6]), [], 'derived: aliphatic carbon chain')
    return qs


# ------------------------------------------------------------------ search and display

def run_search(md):
    import contextlib, io
    with contextlib.redirect_stdout(io.StringIO()):
        C.cmd_search(SimpleNamespace(moldir=str(md)))
    return json.loads((md / 'decisions/analogue_search.json').read_text())


def ranked_hits(md, ana, query=None):
    tgt = Chem.MolFromSmiles(json.loads((md / 'input/request.json').read_text())['original_smiles'])
    tfp = _FP.GetFingerprint(tgt)
    rows = []
    for h in ana['hits']:
        if query and query not in h['queries']:
            continue
        m = Chem.MolFromSmiles(h['canonical_smiles']) if h.get('canonical_smiles') else None
        sim = DataStructs.TanimotoSimilarity(tfp, _FP.GetFingerprint(m)) if m else 0.0
        rows.append((sim, h))
    rows.sort(key=lambda x: -x[0])
    return rows


def print_hits(md, ana, query=None, limit=15):
    rows = ranked_hits(md, ana, query)
    refs = ref_beads(md)
    print(f'{len(rows)} hits' + (f' for query {query}' if query else '') + ' (ranked by Morgan similarity to the target):')
    for sim, h in rows[:limit]:
        r = refs.get(h['reference_id'], {})
        types = ' '.join(b['type'] or '?' for b in r.get('beads', [])) or f"(no bead view: {r.get('error', '?')})"
        print(f"  {sim:.2f}  {h['reference_id']:<28} {h['name'][:30]:<30} {h.get('canonical_smiles') or '':<28} {types}")
    if len(rows) > limit:
        print(f'  ... {len(rows) - limit} more (design.py hits <moldir> --limit N)')


def atom_table(mol):
    m = Chem.Mol(mol)
    for a in m.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1 if a.GetAtomicNum() > 1 else 0)
    print('heavy atoms (one-based IDs as atom-map numbers):', Chem.MolToSmiles(Chem.RemoveHs(m), canonical=False))
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            nb = [n.GetIdx() + 1 for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
            print(f'  {a.GetIdx() + 1:>3} {a.GetSymbol():<2} {"arom" if a.GetIsAromatic() else "    "} '
                  f'H{a.GetTotalNumHs(includeNeighbors=True)} ring={"y" if a.IsInRing() else "n"} neighbours={nb}')


# ------------------------------------------------------------------ commands

def cmd_init(a):
    md = Path(a.moldir)
    bound_check(md)
    if (md / DESIGN).exists():
        die(f'{md / DESIGN} exists; continue with the other commands (`design.py show`)')
    req = json.loads((md / 'input/request.json').read_text())
    mol = target(md)
    agent = a.agent or os.environ.get('MARTINI_AGENT') or 'unknown'
    model = a.model or os.environ.get('MARTINI_MODEL') or 'unknown'
    des = {'schema_version': 1, 'molecule_id': req['molecule_id'], 'resname': req['molecule_id'][:5],
           'authored_by': {'agent': agent, 'model': model, 'tool': 'scripts/design.py', 'session': session()},
           'analogue_queries': derive_queries(mol), 'analogue_annotations': {},
           'default_rejection_reason': 'not opened with design.py open (every opened reference has its own reason)',
           'mapping_candidates': [], 'chosen_candidate': None, 'selection_rationale': None, 'beads': [],
           'virtual_sites': [], 'mass_policy': {'source': 'computed by stage 4 (d5bonded.terms.bead_masses); none authored'},
           'hydrogen_excluded_parents': {'atom_ids': []}, 'unverified_statements': [], 'limitations': [],
           'manual_interventions': [], 'revisions': []}
    (md / 'decisions').mkdir(exist_ok=True)
    save(md, des, 'init', vars(a))
    ana = run_search(md)
    atom_table(mol)
    print(f"\n{len(des['analogue_queries'])} queries derived from the molecule:")
    for q in des['analogue_queries']:
        n = sum(1 for h in ana['hits'] if q['name'] in h['queries'])
        print(f"  {q['name']:<34} {q['smarts']:<30} {n} hits")
    print(f"exposure: exact {ana['exact_target_matches']} tautomer {ana['tautomer_target_matches']}")
    print()
    print_hits(md, ana)
    print('\nNext: `design.py open` the relevant hits, `design.py match` proposed beads, then candidate/choose/bead.')


def cmd_atoms(a):
    atom_table(target(Path(a.moldir)))


def cmd_hits(a):
    md = Path(a.moldir)
    print_hits(md, json.loads((md / 'decisions/analogue_search.json').read_text()), a.query, a.limit)


def cmd_query(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    if Chem.MolFromSmarts(a.smarts) is None:
        die(f'invalid SMARTS {a.smarts!r}')
    if any(q['name'] == a.name for q in des['analogue_queries']):
        die(f'query {a.name} exists')
    des['analogue_queries'].append({'name': a.name, 'smarts': a.smarts,
                                    'motif_tags': [t for t in (a.tags or '').split(',') if t],
                                    'search_bible': not a.no_bible, 'origin': 'agent'})
    save(md, des, 'query', vars(a))
    ana = run_search(md)
    print_hits(md, ana, a.name)


def cmd_open(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    refs = ref_beads(md)
    if a.ref not in refs:
        die(f'{a.ref} is not in the index or the Bible')
    r = refs[a.ref]
    print(f"{a.ref}: {r['name']}\n  source: {r['source_path']}  {r.get('locator') or ''}")
    if r.get('mapping_path'):
        print(f"  mapping: {r['mapping_path']}\n  figure: {r.get('mapping_figure')}")
    opened = [r['source_path']] + [p for p in (r.get('mapping_path'), r.get('mapping_figure')) if p]
    if 'error' in r:
        print(f"  no bead view: {r['error']}. Read the source file itself; this model cannot support rung 1.")
    else:
        print('  beads (fragment signature: heavy atoms, * = bond out of the bead, [50c] = atom shared at weight 0.5):')
        for b in r['beads']:
            extra = f"  atoms in reference SMILES {b['heavy_atoms_in_reference_smiles']}" if 'heavy_atoms_in_reference_smiles' in b else \
                    (f"  CGSmiles fragment {b['fragment_cgsmiles']}" if 'fragment_cgsmiles' in b else '')
            print(f"    {b['name']:<8} {b['type'] or '?':<6} {b['signature']:<24}{extra}")
    ann = des['analogue_annotations'].setdefault(a.ref, {'used': False, 'used_for': [], 'opened': [], 'reason': None})
    ann['opened'] = sorted(set(ann.get('opened', [])) | set(opened))
    save(md, des, 'open', vars(a))
    print(f'\nRecorded as opened. Give it a reason with `design.py annotate {md} {a.ref} --reason ...` '
          '(used beads are annotated automatically by `design.py bead`).')


def cmd_annotate(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    ann = des['analogue_annotations'].get(a.ref)
    if ann is None:
        die(f'{a.ref} was not opened; `design.py open` it first')
    ann['reason'] = a.reason
    if a.used_for:
        ann['used'] = True
        ann['used_for'] = sorted(set(ann.get('used_for', [])) | {a.used_for})
    save(md, des, 'annotate', vars(a))
    print('ok')


def cmd_match(a):
    md = Path(a.moldir)
    mol = target(md)
    mem = parse_atoms(a.atoms, mol)
    sig = signature(mol, mem)
    pre = precedents(md, sig)
    print(f'target fragment {frag_smiles(mol, list(mem))}  signature {sig}  ({sum(mem.values()):g} heavy atoms)')
    if not pre:
        n = sum(len(r.get('beads', [])) for r in ref_beads(md).values())
        print(f'no indexed bead holds this fragment ({n} reference beads checked): rung 1 is not available for it')
        return
    by = {}
    for p in pre:
        by.setdefault(p['type'], []).append(f"{p['reference_id']}:{p['bead']}")
    for t, lst in sorted(by.items(), key=lambda x: -len(x[1])):
        print(f'  {t:<6} {len(lst):>3}  {", ".join(lst[:8])}{" ..." if len(lst) > 8 else ""}')
    print('Rung 1 needs one of these, opened with `design.py open`, cited as --from REF_ID:BEAD.')


def partition(spec, mol):
    beads = [parse_atoms(p, mol) for p in spec.split('|') if p.strip()]
    tot = {}
    for b in beads:
        for i, w in b.items():
            tot[i] = tot.get(i, 0.0) + w
    miss = [i for i in heavy_ids(mol) if abs(tot.get(i, 0.0) - 1.0) > 1e-6]
    if miss:
        die(f'allocation weights do not sum to one for atoms {miss}')
    return beads


def cmd_candidate(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    mol = target(md)
    beads = partition(a.beads, mol)
    cands = [c for c in des['mapping_candidates'] if c['id'] != a.id]
    if len(cands) >= 3:
        die('at most three mapping candidates (SKILL.md stage 2.4)')
    rec = {'id': a.id, 'beads': [{'name': f'B{k + 1}', 'members': {str(i): w for i, w in b.items()}} for k, b in enumerate(beads)],
           'summary': a.summary, 'pros': [x.strip() for x in a.pros.split(';') if x.strip()],
           'cons': [x.strip() for x in a.cons.split(';') if x.strip()]}
    des['mapping_candidates'] = cands + [rec]
    save(md, des, 'candidate', vars(a))
    for k, b in enumerate(beads):
        sig = signature(mol, b)
        pre = precedents(md, sig)
        types = sorted({p['type'] for p in pre})
        print(f"  B{k + 1}: {frag_smiles(mol, list(b)):<14} {sum(b.values()):g} heavy atoms  signature {sig:<22} "
              f"rung-1 precedents: {len(pre)} {types if types else ''}")


def cmd_choose(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    if a.candidate not in [c['id'] for c in des['mapping_candidates']]:
        die(f'no candidate {a.candidate}')
    des['chosen_candidate'] = a.candidate
    des['selection_rationale'] = a.rationale
    save(md, des, 'choose', vars(a))
    print('ok')


def cmd_bead(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    mol = target(md)
    mem = parse_atoms(a.atoms, mol)
    sig = signature(mol, mem)
    pre = precedents(md, sig)
    types = C.core_types()
    others = [b for b in des['beads'] if b['name'] != a.name]
    exists = len(others) != len(des['beads'])
    if exists and not a.replace:
        die(f'bead {a.name} exists; pass --replace to redefine it')
    if len(a.rationale.strip()) < 25:
        die('the rationale must say, for this bead, why these atoms and this type')
    if any(b['rationale'].strip() == a.rationale.strip() for b in others):
        die('another bead of this molecule has the same rationale; write the reason for this bead')
    rec = {'name': a.name, 'charge_e': a.charge, 'members': {str(i): w for i, w in mem.items()},
           'rationale': a.rationale, 'fragment_signature': sig}
    pre_txt = ', '.join(f"{p['reference_id']}:{p['bead']} {p['type']}" for p in pre[:12]) + (' ...' if len(pre) > 12 else '')
    refs = ref_beads(md)
    used = []
    if a.rung == 1:
        if not a.source_bead or ':' not in a.source_bead:
            die('rung 1 needs --from REF_ID:BEAD (see `design.py match`)')
        rid, _, bname = a.source_bead.rpartition(':')
        r = refs.get(rid)
        if r is None or 'beads' not in r:
            die(f'{rid} has no bead view ({(r or {}).get("error", "unknown reference")}); it cannot support rung 1. Use rung 2 and cite it.')
        if rid not in des['analogue_annotations']:
            die(f'{rid} was not opened; `design.py open {md} {rid}` first')
        sb = next((b for b in r['beads'] if b['name'] == bname), None)
        if sb is None:
            die(f'{rid} has no bead {bname}; its beads: {[b["name"] for b in r["beads"]]}')
        if sb['signature'] != sig:
            die(f'not the same fragment: target {sig} vs {rid}:{bname} {sb["signature"]}. Rung 1 needs the same heavy atoms '
                f'in the same bonding context. Precedents that do match: {pre_txt or "none"}; otherwise use rung 2.')
        if a.type and a.type != sb['type']:
            die(f'rung 1 copies the source type {sb["type"]}; --type {a.type} contradicts it')
        rec.update({'type': sb['type'], 'evidence_status': 'transferred', 'provisional_type': False,
                    'typing_ladder': {'rung_used': 1, 'rung_results': {
                        'rung_1': f'same fragment {sig} as {rid} bead {bname} ({sb["type"]}); verified by design.py '
                                  f'signature match. All indexed beads holding it: {pre_txt}'}},
                    'sources': [{'path': r['source_path'], 'locator': f"{r.get('locator') or ''}; bead {bname}".strip('; ')}]
                    + ([{'path': r['mapping_path'], 'locator': f'group of bead {bname}'}] if r.get('mapping_path') else [])})
        used.append(rid)
    elif a.rung in (2, 3):
        if not a.type or a.type not in types:
            die(f'--type must be a bead type of the core ITP ({a.type!r} is not)')
        if pre and not a.rung1_rejection:
            die(f'rung 1 is available for this fragment ({pre_txt}). Use it, or say why not with --rung1-rejection.')
        if not a.sources:
            die('--source PATH:LOCATOR is required (repeatable)')
        srcs = [parse_source(s) for s in a.sources]
        r1 = (f'indexed beads holding {sig}: {pre_txt}; not used: {a.rung1_rejection}' if pre else
              f'no indexed bead holds fragment {sig} ({sum(len(x.get("beads", [])) for x in refs.values())} reference beads checked)')
        if a.rung == 2:
            if not (a.group and a.rule):
                die('rung 2 needs --group (the documented chemical group, Souza SI Table 24/25) and --rule (the SI C2 proximity rule applied)')
            res = {'rung_1': r1, 'rung_2': f'group: {a.group}; rule: {a.rule}'}
        else:
            if not (a.table and a.values and a.rung2_finding):
                die('rung 3 needs --table, --values and --rung2-finding (why no documented group fits)')
            res = {'rung_1': r1, 'rung_2': a.rung2_finding, 'rung_3': f'{a.table}: {a.values}'}
        rec.update({'type': a.type, 'evidence_status': 'adapted', 'provisional_type': bool(a.provisional),
                    'typing_ladder': {'rung_used': a.rung, 'rung_results': res}, 'sources': srcs})
        index_paths = {x['source_path']: k for k, x in refs.items()}
        used += [index_paths[s['path']] for s in srcs if s['path'] in index_paths and index_paths[s['path']] in des['analogue_annotations']]
    else:
        die('--rung must be 1, 2 or 3 (a rung-4 bead is unresolved: record it with `design.py note --limitation` and stop; the model is blocked)')
    if a.size_exception:
        rec['size_exception'] = a.size_exception
    for rid in used:
        ann = des['analogue_annotations'][rid]
        ann['used'] = True
        ann['used_for'] = sorted(set(ann.get('used_for', [])) | {f'{a.name} type {rec["type"]} (rung {a.rung})'})
        ann['reason'] = ann.get('reason') or f'source of bead {a.name}'
    des['beads'] = [rec if b['name'] == a.name else b for b in des['beads']] if exists else des['beads'] + [rec]
    save(md, des, 'bead', vars(a))
    print(f"{a.name}: {rec['type']} rung {a.rung}  fragment {frag_smiles(mol, list(mem))}  signature {sig}")


def cmd_vsite(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    src = parse_source(a.source)
    des['virtual_sites'] = [v for v in des['virtual_sites'] if v['site'] != a.site] + [{
        'site': a.site, 'construction': a.construction, 'constructing': [int(x) for x in a.constructing.split(',')],
        'evidence_class': a.evidence_class, 'source': src}]
    save(md, des, 'vsite', vars(a))
    print('ok')


def cmd_hydrogens(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    des['hydrogen_excluded_parents'] = {'atom_ids': [int(x) for x in a.excluded_parents.split(',') if x], 'reason': a.reason}
    save(md, des, 'hydrogens', vars(a))
    print('ok')


def cmd_note(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    if a.unverified:
        des['unverified_statements'].append({'statement': a.unverified, 'label': 'agent background knowledge, unverified'})
    if a.limitation:
        des['limitations'].append(a.limitation)
    if a.manual:
        des['manual_interventions'].append({'utc': now(), 'what': a.manual})
    save(md, des, 'note', vars(a))
    print('ok')


def problems(md, des, mol):
    p = []
    if not des['mapping_candidates']:
        p.append('no mapping candidate (`design.py candidate`)')
    if not des['chosen_candidate']:
        p.append('no chosen candidate (`design.py choose`)')
    else:
        cand = next(c for c in des['mapping_candidates'] if c['id'] == des['chosen_candidate'])
        cs = sorted(tuple(sorted(b['members'].items())) for b in cand['beads'])
        bs = sorted(tuple(sorted(b['members'].items())) for b in des['beads'])
        if cs != bs:
            p.append(f"the beads do not partition the molecule as candidate {des['chosen_candidate']} does")
    tot = {}
    for b in des['beads']:
        for i, w in b['members'].items():
            tot[int(i)] = tot.get(int(i), 0.0) + w
    miss = [i for i in heavy_ids(mol) if abs(tot.get(i, 0.0) - 1.0) > 1e-6]
    if miss:
        p.append(f'heavy atoms not fully allocated to beads: {miss}')
    for rid, ann in des['analogue_annotations'].items():
        if not ann.get('reason'):
            p.append(f'opened reference {rid} has no reason (`design.py annotate`)')
    req = json.loads((md / 'input/request.json').read_text())
    q = sum(b.get('charge_e', 0.0) for b in des['beads'])
    if abs(q - req['formal_charge']) > 1e-9:
        p.append(f"bead charges sum to {q}, formal charge is {req['formal_charge']}")
    return p


def cmd_show(a):
    md = Path(a.moldir)
    des = load(md)
    mol = target(md)
    print(json.dumps({k: des[k] for k in ('molecule_id', 'authored_by', 'chosen_candidate', 'selection_rationale')}, indent=1))
    print(f"queries {len(des['analogue_queries'])}; opened {list(des['analogue_annotations'])}")
    for c in des['mapping_candidates']:
        print(f"candidate {c['id']}: {' | '.join(frag_smiles(mol, [int(i) for i in b['members']]) for b in c['beads'])}  {c['summary']}")
    for b in des['beads']:
        print(f"  {b['name']} {b['type']} rung {b['typing_ladder']['rung_used']} {frag_smiles(mol, [int(i) for i in b['members']])}: {b['rationale']}")
    for x in problems(md, des, mol):
        print('  TODO:', x)


def cmd_check(a):
    md = Path(a.moldir)
    bound_check(md)
    des = load(md)
    mol = target(md)
    p = problems(md, des, mol)
    run_search(md)
    save(md, des, 'check', vars(a))
    if p:
        print('design incomplete:')
        for x in p:
            print('  -', x)
        sys.exit(1)
    print(f"design complete: {len(des['beads'])} beads; analogue_search.json refreshed. Next: cgbuild.py bonded {md}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)

    def p(name, fn):
        q = sp.add_parser(name)
        q.add_argument('moldir')
        q.set_defaults(fn=fn)
        return q
    q = p('init', cmd_init); q.add_argument('--agent'); q.add_argument('--model')
    p('atoms', cmd_atoms)
    q = p('hits', cmd_hits); q.add_argument('--query'); q.add_argument('--limit', type=int, default=25)
    q = p('query', cmd_query); q.add_argument('--name', required=True); q.add_argument('--smarts', required=True)
    q.add_argument('--tags'); q.add_argument('--no-bible', action='store_true')
    q = p('open', cmd_open); q.add_argument('ref')
    q = p('annotate', cmd_annotate); q.add_argument('ref'); q.add_argument('--reason', required=True); q.add_argument('--used-for')
    q = p('match', cmd_match); q.add_argument('--atoms', required=True)
    q = p('candidate', cmd_candidate); q.add_argument('--id', required=True); q.add_argument('--beads', required=True)
    q.add_argument('--summary', required=True); q.add_argument('--pros', required=True); q.add_argument('--cons', required=True)
    q = p('choose', cmd_choose); q.add_argument('--candidate', required=True); q.add_argument('--rationale', required=True)
    q = p('bead', cmd_bead); q.add_argument('--name', required=True); q.add_argument('--atoms', required=True)
    q.add_argument('--rung', type=int, required=True); q.add_argument('--rationale', required=True)
    q.add_argument('--from', dest='source_bead'); q.add_argument('--type'); q.add_argument('--group'); q.add_argument('--rule')
    q.add_argument('--table'); q.add_argument('--values'); q.add_argument('--rung2-finding')
    q.add_argument('--source', dest='sources', action='append'); q.add_argument('--provisional', action='store_true')
    q.add_argument('--rung1-rejection'); q.add_argument('--charge', type=float, default=0.0)
    q.add_argument('--size-exception'); q.add_argument('--replace', action='store_true')
    q = p('vsite', cmd_vsite); q.add_argument('--site', type=int, required=True)
    q.add_argument('--construction', required=True, choices=['virtual_sitesn_cog', 'virtual_sites3'])
    q.add_argument('--constructing', required=True); q.add_argument('--source', required=True)
    q.add_argument('--evidence-class', required=True)
    q = p('hydrogens', cmd_hydrogens); q.add_argument('--excluded-parents', required=True); q.add_argument('--reason', required=True)
    q = p('note', cmd_note); q.add_argument('--unverified'); q.add_argument('--limitation'); q.add_argument('--manual')
    p('show', cmd_show)
    p('check', cmd_check)
    a = ap.parse_args()
    a.fn(a)


if __name__ == '__main__':
    main()
