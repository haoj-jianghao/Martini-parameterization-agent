"""Vectorized bonded geometry on bead positions R [n_frames, n_beads, 3] (nm).

distances: |r_j - r_i| (nm)
angles:    angle i-j-k at j, degrees in [0, 180]
dihedrals: IUPAC/GROMACS signed dihedral i-j-k-l, degrees in (-180, 180]
           (angle between planes i-j-k and j-k-l; 0 = cis, 180 = trans; positive
           = clockwise rotation of l about j->k looking from j). The same formula
           gives the GROMACS improper (funct 2/4) angle for an i-j-k-l quadruplet.
           Also returns the two bending angles (i-j-k, j-k-l): the dihedral is
           ill-defined when either approaches 180 degrees.
"""
import numpy as np


def _take(R, idx, col):
    return R[:, idx[:, col], :]


def distances(R, pairs):
    pairs = np.asarray(pairs, int).reshape(-1, 2)
    return np.linalg.norm(_take(R, pairs, 1) - _take(R, pairs, 0), axis=-1)


def angles(R, triples):
    t = np.asarray(triples, int).reshape(-1, 3)
    a = _take(R, t, 0) - _take(R, t, 1)
    b = _take(R, t, 2) - _take(R, t, 1)
    # atan2(|a x b|, a.b) is accurate near 0 and 180, unlike arccos
    return np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b), axis=-1), np.sum(a * b, -1)))


def dihedrals(R, quads, return_bends=False):
    q = np.asarray(quads, int).reshape(-1, 4)
    b1 = _take(R, q, 1) - _take(R, q, 0)
    b2 = _take(R, q, 2) - _take(R, q, 1)
    b3 = _take(R, q, 3) - _take(R, q, 2)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    y = np.linalg.norm(b2, axis=-1) * np.sum(b1 * n2, -1)
    x = np.sum(n1 * n2, -1)
    phi = np.degrees(np.arctan2(y, x))
    phi = np.where(phi <= -180.0, phi + 360.0, phi)  # (-180, 180]
    if not return_bends:
        return phi
    bends = np.stack([angles(R, q[:, [0, 1, 2]]), angles(R, q[:, [1, 2, 3]])], -1)
    return phi, bends
