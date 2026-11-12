"""Bead-level view of the indexed reference models, for rung-1 typing (SKILL.md stage 3).

Rung 1 transfers a type only from a source bead that holds the same heavy atoms in the same bonding
context as the target bead. This module turns each indexed model into beads with member atoms and a
*fragment signature*: the canonical SMILES of the bead's heavy atoms, with one dummy atom (*) per bond to
an atom outside the bead and a mass label on atoms the bead only shares (allocation weight below one,
label = round(100 * weight)). Two beads hold the same fragment exactly when their signatures are equal.

Sources of atom membership, by index entry:
- `mapping_path` (.ndx, CG-builder) over `mapping_atom_reference` (all-atom ITP): atom i of the ndx group
  of bead b has weight count(i in b) / count(i in all groups) (repeated indices carry weight);
- `mapping_cgsmiles` and Bible rows: the CGSmiles fragment strings (bonding descriptors become dummies,
  `[C;0.5]` weights become labels, weight-0 atoms are context only);
- single-bead models without a mapping file: the whole molecule.
A multi-bead model without a mapping file cannot establish exact transfer (data-contract.md) and yields
no beads. Deterministic; contains no chemical decisions.
"""
import json, re
from pathlib import Path

from rdkit import Chem, RDLogger

RDLogger.DisableLog('rdApp.*')
PT = Chem.GetPeriodicTable()


# ------------------------------------------------------------------ signatures

def _normalise(m):
    """Strip H bookkeeping and make dummy bonds generic, so both sides compare on heavy-atom graph,
    element, aromaticity, formal charge, bond order and weight label only."""
    m = Chem.RWMol(m)
    for at in m.GetAtoms():
        at.SetNoImplicit(True)
        at.SetNumExplicitHs(0)
        at.SetAtomMapNum(0)
        at.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
        if at.GetAtomicNum() == 0:
            at.SetIsotope(0)
            at.SetIsAromatic(False)
    for b in m.GetBonds():
        b.SetStereo(Chem.BondStereo.STEREONONE)
        b.SetBondDir(Chem.BondDir.NONE)
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if a1.GetAtomicNum() == 0 or a2.GetAtomicNum() == 0:
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)


def fragment_signature(mol, weights):
    """Signature of the bead holding `weights` {atom index (0-based): allocation weight} of `mol`.
    Hydrogens are ignored; each bond from a member heavy atom to a non-member heavy atom is a dummy."""
    members = {i: w for i, w in weights.items() if w > 1e-9 and mol.GetAtomWithIdx(i).GetAtomicNum() > 1}
    rw = Chem.RWMol()
    idx = {}
    for i, w in members.items():
        a = mol.GetAtomWithIdx(i)
        na = Chem.Atom(a.GetAtomicNum())
        na.SetFormalCharge(a.GetFormalCharge())
        na.SetIsAromatic(a.GetIsAromatic())
        if w < 0.995:
            na.SetIsotope(int(round(100 * w)))
        idx[i] = rw.AddAtom(na)
    for i in members:
        for b in mol.GetAtomWithIdx(i).GetBonds():
            j = b.GetOtherAtomIdx(i)
            if mol.GetAtomWithIdx(j).GetAtomicNum() <= 1:
                continue
            if j in members:
                if i < j:
                    rw.AddBond(idx[i], idx[j], b.GetBondType())
                    if b.GetIsAromatic():
                        rw.GetBondWithIdx(rw.GetNumBonds() - 1).SetIsAromatic(True)
            else:
                d = rw.AddAtom(Chem.Atom(0))
                rw.AddBond(idx[i], d, Chem.BondType.SINGLE)
    return _normalise(rw.GetMol())


# ------------------------------------------------------------------ CGSmiles fragments

_TOK = re.compile(r'(\[[^\]]*\]|Br|Cl|[BCNOPSFI]|[bcnops]|%\d\d|\d|[-=#:/\\.$]|\(|\))')


def _cgsmiles_fragment_to_smiles(s):
    """One CGSmiles fragment (e.g. `[$]COC[$]`, `O[C;0.5][$a]`) -> SMILES with dummies for bonding
    descriptors, isotope labels for partial weights (weight 0 -> dummy). Raises ValueError if unparsable."""
    toks = _TOK.findall(s)
    if ''.join(toks) != s:
        raise ValueError(f'unparsed characters in {s!r}')
    out, leading, pending_bond, seen_atom = [], 0, '', False
    for t in toks:
        if t.startswith('[') and t[1:2] in '$<>!':
            if not seen_atom:
                leading += 1
            else:
                out.append(f'({pending_bond}*)')
                pending_bond = ''
            continue
        if t in '-=#:/\\' and len(t) == 1:
            pending_bond = t
            continue
        if t == '.':
            raise ValueError('disconnected fragment')
        if t.startswith('['):
            body = t[1:-1]
            if ';' in body:
                f = body.split(';')
                w = float(f[1]) if re.match(r'^[0-9.]+$', f[1]) else 1.0
                body = f[0]
                if w < 0.005:
                    t = '*'
                elif w < 0.995:
                    t = f'[{int(round(100 * w))}{body}]'
                else:
                    t = f'[{body}]'
        if pending_bond:
            out.append(pending_bond)
            pending_bond = ''
        out.append(t)
        if (t[0].isalpha() or t.startswith('[') or t == '*') and not seen_atom:
            seen_atom = True
            out.extend(['(*)'] * leading)
            leading = 0
    return ''.join(out)


def cgsmiles_beads(cgsmiles, valid_types):
    """Beads of a CGSmiles string: [{name, type, signature}] in CG-graph order, or raises ValueError."""
    parts = re.findall(r'\{[^{}]*\}', cgsmiles)
    if len(parts) < 2:
        raise ValueError('no fragment block')
    frags = {}
    for frag in parts[1].strip('{}').split(','):
        m = re.match(r'#([^=]+)=(.*)$', frag.strip())
        if not m:
            raise ValueError(f'bad fragment {frag!r}')
        frags[m.group(1)] = m.group(2)
    cg_names = re.findall(r'\[#([^\]]+)\]', parts[0])
    out = []
    for k, name in enumerate(cg_names):
        if name not in frags:
            raise ValueError(f'fragment {name} not defined')
        typ = name if name in valid_types else re.sub(r'[A-Z0-9]$', '', name)
        if typ not in valid_types:
            typ = re.sub(r'[A-Z]+$', '', name)
        if typ not in valid_types:
            raise ValueError(f'fragment name {name} is not a bead type')
        smi = _cgsmiles_fragment_to_smiles(frags[name])
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is None:
            raise ValueError(f'RDKit cannot parse {smi!r}')
        # drop context beyond dummies (weight-0 atoms written as * may carry further atoms)
        rw = Chem.RWMol(m)
        for at in list(rw.GetAtoms())[::-1]:
            if at.GetAtomicNum() == 0 and all(n.GetAtomicNum() == 0 for n in at.GetNeighbors()):
                rw.RemoveAtom(at.GetIdx())
        rw2 = Chem.RWMol(rw.GetMol())
        rm = []
        for at in list(rw2.GetAtoms()):
            if at.GetAtomicNum() > 0:
                continue
            heavy_n = [n for n in at.GetNeighbors() if n.GetAtomicNum() > 0]
            if len(heavy_n) > 1:  # a weight-0 atom bonded to several members: one dummy per bond
                for n in heavy_n[1:]:
                    d = rw2.AddAtom(Chem.Atom(0))
                    rw2.AddBond(n.GetIdx(), d, Chem.BondType.SINGLE)
                    rm.append((at.GetIdx(), n.GetIdx()))
        for a, b in rm:
            rw2.RemoveBond(a, b)
        fm = rw2.GetMol()
        for at in fm.GetAtoms():
            if at.GetAtomicNum() > 0 and at.GetSymbol().islower():
                at.SetIsAromatic(True)
        out.append({'name': f'{name}#{k + 1}', 'type': typ, 'signature': _normalise(_aromatic_bonds(fm)),
                    'fragment_cgsmiles': frags[name]})
    return out


def _aromatic_bonds(m):
    rw = Chem.RWMol(m)
    for b in rw.GetBonds():
        if b.GetBeginAtom().GetIsAromatic() and b.GetEndAtom().GetIsAromatic() and b.GetBondType() == Chem.BondType.SINGLE:
            b.SetBondType(Chem.BondType.AROMATIC)
            b.SetIsAromatic(True)
    return rw.GetMol()


# ------------------------------------------------------------------ ITP / ndx models

def _itp_atoms_bonds(path):
    sec, atoms, bonds = None, [], []
    for line in open(path):
        s = line.split(';')[0].strip()
        if not s or s.startswith('#'):
            continue
        m = re.match(r'\[\s*(\w+)\s*\]', s)
        if m:
            sec = m.group(1).lower()
            continue
        f = s.split()
        if sec == 'atoms':
            atoms.append({'nr': int(f[0]), 'type': f[1], 'name': f[4],
                          'mass': float(f[7]) if len(f) > 7 else None})
        elif sec in ('bonds', 'constraints'):
            bonds.append((int(f[0]), int(f[1])))
    return atoms, bonds


def _element_from_mass(mass):
    best = min(range(1, 55), key=lambda z: abs(PT.GetAtomicWeight(z) - mass))
    return best


def _parse_ndx(path):
    groups, cur = [], None
    for line in open(path):
        s = line.strip()
        if s.startswith('['):
            cur = []
            groups.append((s.strip('[] '), cur))
        elif s and cur is not None:
            cur.extend(int(x) for x in s.split())
    return groups


def _aa_to_smiles_perm(smiles, aa_itp):
    """Map all-atom ITP atom numbers (1-based) to atom indices of Chem.AddHs(MolFromSmiles(smiles))."""
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    atoms, bonds = _itp_atoms_bonds(aa_itp)
    if len(atoms) != mol.GetNumAtoms():
        raise ValueError(f'{len(atoms)} AA atoms vs {mol.GetNumAtoms()} in SMILES')
    q = Chem.RWMol()
    for a in atoms:
        q.AddAtom(Chem.Atom(_element_from_mass(a['mass'])))
    for i, j in bonds:
        q.AddBond(i - 1, j - 1, Chem.BondType.SINGLE)
    g = Chem.RWMol(mol)
    Chem.Kekulize(g, clearAromaticFlags=True)
    for b in g.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
    for at in g.GetAtoms():
        at.SetFormalCharge(0)
        at.SetNoImplicit(True)
    for at in q.GetAtoms():
        at.SetNoImplicit(True)
    gm, qm = g.GetMol(), q.GetMol()
    gm.UpdatePropertyCache(strict=False)
    qm.UpdatePropertyCache(strict=False)
    match = gm.GetSubstructMatch(qm)
    if not match:
        raise ValueError('AA topology graph does not match the SMILES')
    return mol, {k + 1: match[k] for k in range(len(atoms))}


def ndx_beads(entry, project):
    mol, perm = _aa_to_smiles_perm(entry['canonical_isomeric_smiles'], project / entry['mapping_atom_reference'])
    cg_atoms, _ = _itp_atoms_bonds(project / entry['source_path'])
    groups = _parse_ndx(project / entry['mapping_path'])
    if len(groups) != len(cg_atoms):
        raise ValueError(f'{len(groups)} ndx groups vs {len(cg_atoms)} CG atoms')
    total = {}
    for _, g in groups:
        for a in g:
            total[a] = total.get(a, 0) + 1
    out = []
    for (gname, g), cg in zip(groups, cg_atoms):
        cnt = {}
        for a in g:
            cnt[a] = cnt.get(a, 0) + 1
        w = {perm[a]: c / total[a] for a, c in cnt.items()}
        heavy = sorted(i for i in w if mol.GetAtomWithIdx(i).GetAtomicNum() > 1)
        if not heavy:
            continue
        out.append({'name': cg['name'], 'type': cg['type'], 'signature': fragment_signature(mol, w),
                    'heavy_atoms_in_reference_smiles': {str(i + 1): round(w[i], 3) for i in heavy}})
    return out


def single_bead(entry, project):
    cg_atoms, _ = _itp_atoms_bonds(project / entry['source_path']) if str(entry['source_path']).endswith('.itp') else ([], [])
    mol = Chem.MolFromSmiles(entry['canonical_isomeric_smiles'])
    return [{'name': cg_atoms[0]['name'] if cg_atoms else 'B1', 'type': cg_atoms[0]['type'] if cg_atoms else None,
             'signature': fragment_signature(mol, {a.GetIdx(): 1.0 for a in mol.GetAtoms()})}]


def _one_moleculetype_atoms(entry, project):
    """CG atoms of the entry's moleculetype (some source files hold several)."""
    path = project / entry['source_path']
    if not str(path).endswith('.itp'):
        return None
    want = entry['source_locator'].split()[1] if entry.get('source_locator', '').startswith('moleculetype') else None
    cur, sec, atoms = None, None, {}
    for line in open(path):
        s = line.split(';')[0].strip()
        if not s or s.startswith('#'):
            continue
        m = re.match(r'\[\s*(\w+)\s*\]', s)
        if m:
            sec = m.group(1).lower()
            continue
        if sec == 'moleculetype':
            cur = s.split()[0]
            atoms[cur] = []
        elif sec == 'atoms' and cur:
            f = s.split()
            atoms[cur].append({'name': f[4], 'type': f[1]})
    if want and want in atoms:
        return atoms[want]
    return next(iter(atoms.values())) if len(atoms) == 1 else None


def reference_beads(index_path, bible_rows, valid_types, project):
    """All verifiable reference beads: {reference_id: {'beads': [...]} or {'error': why}}."""
    out = {}
    for line in open(index_path):
        e = json.loads(line)
        rid = e['reference_id']
        try:
            if not e.get('canonical_isomeric_smiles'):
                raise ValueError('no verified SMILES')
            if e.get('mapping_path') and e.get('mapping_atom_reference'):
                beads, how = ndx_beads(e, project), 'ndx'
            elif e.get('mapping_cgsmiles'):
                beads, how = cgsmiles_beads(e['mapping_cgsmiles'], valid_types), 'cgsmiles'
            else:
                atoms = _one_moleculetype_atoms(e, project)
                if atoms is not None and len(atoms) == 1:
                    mol = Chem.MolFromSmiles(e['canonical_isomeric_smiles'])
                    beads = [{'name': atoms[0]['name'], 'type': atoms[0]['type'],
                              'signature': fragment_signature(mol, {a.GetIdx(): 1.0 for a in mol.GetAtoms()})}]
                    how = 'single_bead'
                else:
                    raise ValueError('multi-bead model without an atom-level mapping: exact transfer cannot be established')
            out[rid] = {'how': how, 'source_path': e['source_path'], 'locator': e.get('source_locator'),
                        'mapping_path': e.get('mapping_path'), 'mapping_figure': e.get('mapping_figure'),
                        'name': e['name'], 'beads': beads}
        except Exception as ex:  # noqa: BLE001  (recorded: such an entry cannot support rung 1)
            out[rid] = {'error': str(ex), 'source_path': e['source_path'], 'name': e['name']}
    for r in bible_rows:
        try:
            if not r.get('cgsmiles'):
                raise ValueError('no CGSmiles')
            out[r['reference_id']] = {'how': 'cgsmiles', 'source_path': r['source_path'], 'locator': r['locator'],
                                      'name': r['name'], 'beads': cgsmiles_beads(r['cgsmiles'], valid_types)}
        except Exception as ex:  # noqa: BLE001
            out[r['reference_id']] = {'error': str(ex), 'source_path': r['source_path'], 'name': r['name']}
    return out
