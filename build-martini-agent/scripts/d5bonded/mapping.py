"""Atom-to-bead mapping as a weight matrix M (n_beads x n_atoms): R_bead = M @ x_atoms.

Mapping JSON schema (written by scripts/extract_mapping.py; any mapping method
should write the same schema):
  mol_id, smiles, method, n_atoms, atoms [element per atom, AA order],
  beads: [{index, name, type, virtual_site, atoms [0-based], weights [sum 1]}]
Weights are the center-of-geometry weights including hydrogens (Martini 3);
an atom shared between beads carries a fractional weight in each.
"""
import hashlib
import json

import numpy as np


class Mapping:
    def __init__(self, M, atoms=None, bead_names=None, bead_types=None, virtual_site=None,
                 mol_id=None, method=None, source=None, sha256=None):
        M = np.asarray(M, float)
        if M.ndim != 2:
            raise ValueError("M must be 2D")
        if np.any(M < 0):
            raise ValueError("negative mapping weight")
        if not np.allclose(M.sum(1), 1.0, atol=1e-9):
            raise ValueError("each bead's weights must sum to 1")
        self.M = M
        self.atoms = atoms
        nb = M.shape[0]
        self.bead_names = bead_names or [f"B{i + 1}" for i in range(nb)]
        self.bead_types = bead_types or [""] * nb
        self.virtual_site = np.asarray(virtual_site if virtual_site is not None else [False] * nb, bool)
        self.mol_id, self.method, self.source, self.sha256 = mol_id, method, source, sha256

    @property
    def n_beads(self):
        return self.M.shape[0]

    @property
    def n_atoms(self):
        return self.M.shape[1]

    def bead_atoms(self, b):
        return np.nonzero(self.M[b])[0]

    @classmethod
    def from_groups(cls, groups, n_atoms, weights=None, **kw):
        """groups: list of atom-index lists; weights: optional matching lists (normalized here)."""
        M = np.zeros((len(groups), n_atoms))
        for b, g in enumerate(groups):
            w = np.ones(len(g)) if weights is None else np.asarray(weights[b], float)
            np.add.at(M[b], np.asarray(g, int), w)
            M[b] /= M[b].sum()
        return cls(M, **kw)

    @classmethod
    def from_json(cls, path):
        raw = open(path, "rb").read()
        d = json.loads(raw)
        beads = sorted(d["beads"], key=lambda b: b["index"])
        mp = cls.from_groups([b["atoms"] for b in beads], d["n_atoms"],
                             weights=[b["weights"] for b in beads],
                             atoms=d.get("atoms"), bead_names=[b["name"] for b in beads],
                             bead_types=[b.get("type", "") for b in beads],
                             virtual_site=[b.get("virtual_site", False) for b in beads],
                             mol_id=d.get("mol_id"), method=d.get("method"), source=str(path),
                             sha256=hashlib.sha256(raw).hexdigest())
        mp.raw = d
        return mp

    def check_atoms(self, atoms):
        """Elements of the trajectory must match the mapping's AA atom list."""
        if len(atoms) != self.n_atoms:
            raise ValueError(f"mapping has {self.n_atoms} atoms, trajectory {len(atoms)}")
        if self.atoms is not None and list(self.atoms) != list(atoms):
            bad = [i for i, (a, b) in enumerate(zip(self.atoms, atoms)) if a != b]
            raise ValueError(f"element mismatch at atoms {bad[:10]}")

    def apply(self, xyz):
        """xyz [..., n_atoms, 3] -> bead positions [..., n_beads, 3] (float64 accumulation)."""
        xyz = np.asarray(xyz)
        if xyz.shape[-2] != self.n_atoms:
            raise ValueError(f"expected {self.n_atoms} atoms, got {xyz.shape[-2]}")
        return np.einsum("ba,...ai->...bi", self.M, xyz.astype(np.float64))
