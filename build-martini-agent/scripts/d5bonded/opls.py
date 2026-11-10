"""OPLS-AA topology for stage 4 from opls-ff-gen (OPLS/2020), reduced to one interface (no MD).

opls-ff-gen gives 1.14*CM1A-LBCC charges for uncharged species. Its default GFN2-xTB geometry step runs in vacuum
and moves protons in zwitterions and ammonium ions (NH3+ -> COO-), which it then rejects; charged species
therefore use an MMFF94 geometry (bonding fixed) and CM1A charges. opls-ff-gen imports OpenMM only for its --qm
torsion fit, which is not used here.
Atom order is engine-defined. It is mapped onto the RDKit AddHs order by graph isomorphism, so a changed bond graph
is an error rather than a silent mislabelling. Equivalent atoms may swap; they carry equal parameters, and
chirality comes from the conformers, not from the force field.
The energy model built on this topology is d5bonded/energy.py.
"""
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from rdkit import Chem

from .energy import TopologyError

ENGINE = "opls-ff-gen"
MAX_ABS_CHARGE = 3.4   # e; P(V) charges near +3.5 from opls-ff-gen blew up in AA simulation
# LJ given to sigma = 0 polar hydrogens in the energy model only (engine files unchanged).
# epsilon is the CHARMM36 polar-H value (0.046 kcal/mol). sigma = 0.1 nm is the smallest tested value that stops
# H-O-P=O / H-N-P=O 1-4 collapse (UDHX, PCRE, ATPH, ADPH); the CHARMM sigma 0.04 nm does not (H stops at 0.8 A).
# Cost against 50 AA references: dihedral W1 sum +0.9%, angle +0.4% (runs/polar_h_benchmark-20261107).
POLAR_H_LJ = (0.1, 0.192)
MASS_TO_Z = {1: 1, 12: 6, 14: 7, 16: 8, 19: 9, 28: 14, 31: 15, 32: 16, 35: 17, 80: 35, 127: 53}


@dataclass
class OplsTopology:
    engine: str
    charge_model: str
    top: Path
    itp: Path
    perm: np.ndarray            # perm[i] = RDKit index of topology atom i
    charges: list
    log: Path = None


def parse_itp(itp):
    """(atomic numbers, 0-based bonds, charges) from a GROMACS itp, using masses for elements."""
    z, q, bonds, sec = [], [], [], None
    for line in open(itp):
        m = re.match(r"^\s*\[\s*(\w+)\s*\]", line)
        if m:
            sec = m.group(1).lower()
            continue
        s = line.split(";")[0].strip()
        if not s or s.startswith("#"):
            continue
        f = s.split()
        if sec == "atoms":
            q.append(float(f[6]))
            z.append(MASS_TO_Z[round(float(f[7]))])
        elif sec in ("bonds", "constraints"):
            bonds.append(tuple(sorted((int(f[0]) - 1, int(f[1]) - 1))))
    return z, sorted(set(bonds)), q


def _bare_graph(z, bonds):
    m = Chem.RWMol()
    for a in z:
        at = Chem.Atom(int(a))
        at.SetNoImplicit(True)
        m.AddAtom(at)
    for i, j in bonds:
        m.AddBond(int(i), int(j), Chem.BondType.SINGLE)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def atom_permutation(mol_h, z, bonds):
    """perm[i] = RDKit index of topology atom i; raises if the bond graphs differ."""
    rd = _bare_graph([a.GetAtomicNum() for a in mol_h.GetAtoms()],
                     [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol_h.GetBonds()])
    if len(z) != rd.GetNumAtoms() or len(bonds) != rd.GetNumBonds():
        raise TopologyError(f"topology has {len(z)} atoms/{len(bonds)} bonds, "
                            f"SMILES has {rd.GetNumAtoms()}/{rd.GetNumBonds()}")
    match = rd.GetSubstructMatch(_bare_graph(z, bonds), useChirality=False)
    if not match:
        raise TopologyError("topology bond graph is not the SMILES bond graph")
    return np.array(match, int)


def has_charged_group(mol):
    """True if any atom outside a nitro group carries a formal charge (ions and zwitterions).
    LBCC is fitted to uncharged species; a nitro group has formal charges but is an uncharged group."""
    nitro = {x for m in mol.GetSubstructMatches(Chem.MolFromSmarts("[N+](=O)[O-]")) for x in m}
    return any(a.GetFormalCharge() and a.GetIdx() not in nitro for a in mol.GetAtoms())


def run_opls_ff_gen(smiles, workdir, command):
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    charged = has_charged_group(Chem.MolFromSmiles(smiles))
    # charged: MMFF94 geometry keeps the protonation state; unscaled CM1A
    extra = ["--optimize", "mmff", "--charge-model", "cm1a"] if charged else []
    log = workdir / "opls_ff_gen.log"
    with open(log, "w") as fh:
        p = subprocess.run([*command, "--smiles", smiles, "--name", "MOL", "--output-dir", ".", *extra],
                           cwd=workdir, stdout=fh, stderr=subprocess.STDOUT)
    itp, top = workdir / "mol.itp", workdir / "mol.top"
    if p.returncode != 0 or not itp.exists():
        raise TopologyError(f"opls-ff-gen exit {p.returncode}: {log.read_text()[-400:]}")
    cm = "CM1A, MMFF94 geometry (OPLS/2020)" if charged else "1.14*CM1A-LBCC (OPLS/2020)"
    return top, itp, cm, log


def build_topology(smiles, workdir, command):
    """opls-ff-gen topology for smiles in workdir/opls-ff-gen; command: opls-ff-gen argv list."""
    mol_h = Chem.AddHs(Chem.MolFromSmiles(smiles))
    top, itp, cm, log = run_opls_ff_gen(smiles, Path(workdir) / ENGINE, command)
    z, bonds, q = parse_itp(itp)
    if max(abs(x) for x in q) > MAX_ABS_CHARGE:
        raise TopologyError(f"atomic charge {max(abs(x) for x in q):.2f} e exceeds {MAX_ABS_CHARGE} e")
    return OplsTopology(ENGINE, cm, top, itp, atom_permutation(mol_h, z, bonds), q, log)
