"""Rule-based Martini 3 bonded topology from a mapping, the atomistic graph and an ensemble.

Which terms exist is decided here, from ring structure and ensemble statistics. The values come from
fit.py. Conventions follow the published small-molecule library (models/itps/opt-mono) and Souza 2021
SI C3-C4:

  Rigid rings   aromatic rings, or rings whose atoms are all sp2. Rings sharing a bead merge into one
                rigid group, which is triangulated with constraints.
                  3 beads: triangle (BENZ, TOLU)
                  4 beads: rhombus, 5 constraints, 1 improper, 1 excluded diagonal (OMA, EBEN)
                  more:    not allowed. Published polycyclics keep four real beads as the hinge and make
                           the rest virtual (NAPH, ANTH, CAFF, TECE; SI C4), a mapping decision. A larger
                           group is still triangulated, but flagged, and the builder blocks the model.
                Virtual sites stay out of the constraint network (NAPH). All pairs within a rigid
                group are excluded.
  Substituent   a bead bonded to a rigid group bridges two adjacent group beads with bonds, gets one
                position term (improper if narrow, Fourier if it has several out-of-plane states), and
                is excluded from the other beads (PBEN, CUME, CLPR).
  Other edges   bonds, never constraints. Stiff bonds stay bonds, with k capped per bead type.
  Angles        every i-j-k path, unless i-k is itself connected (the triangle fixes it) or the angle
                is nearly unrestrained (k below kT).
  Exclusions    rigid groups and substituents as above; also 2nd and 3rd neighbours whose ensemble mean
                distance is below their Lennard-Jones sigma (Souza SI C3 p.44).
  Dihedrals     at most one per central edge, the one with the largest potential range; dropped when
                the range is below kT, or when a flanking angle reaches BEND_LIMIT_DEG in more than
                BEND_FRACTION of the ensemble (SI C4, angles around 180 deg).
"""
import math
from collections import defaultdict
from itertools import combinations

import numpy as np

from . import fit as F
from .geometry import angles, dihedrals, distances

BEND_LIMIT_DEG = 160.0
BEND_LOW_DEG = 20.0           # SI C4 warns about angles around 0 and 180 deg alike
NEAR_LINEAR_TRIANGLE_DEG = 150.0
BEND_FRACTION = 0.05
TIMESTEP_PS = 0.02
MIN_PERIOD_STEPS = 9.0        # with T beads, this gives kb = 2e4 kJ/mol/nm2, the Bartender constraint rule
BOND_K_FLOOR = 700.0          # Souza SI C3 p.44: weaker bonds should be avoided
DIAGONAL_TIE_NM = 0.01
MAX_HULL_EDGE_NM = 0.6        # a longer Delaunay hull edge is a chord across a concave outline
SIZE_NOTE = "Souza SI C4 p.45: use four real hinge beads and make the other ring beads virtual sites"
BOND_K_CAP = {"T": 20000.0, "S": 30000.0, "R": 40000.0}   # kJ/mol/nm2, same-size pairs; proportional to mass


def bead_size(bead_type):
    return bead_type[0] if bead_type[0] in "ST" else "R"


def bond_k_cap(type_i, type_j):
    """Harmonic mean of the two beads' caps, which equals scaling by the reduced mass (T-S: 24000)."""
    a, b = BOND_K_CAP[bead_size(type_i)], BOND_K_CAP[bead_size(type_j)]
    return 2 * a * b / (a + b)


def _omega_limit_sq():
    return (2 * math.pi / (MIN_PERIOD_STEPS * TIMESTEP_PS)) ** 2       # ps^-2


def angle_k_limit(pos, trip, masses, theta0_deg):
    """Largest G96 k (kJ/mol) whose local vibration keeps MIN_PERIOD_STEPS at the timestep.
    Wilson G element of the bend, with k_harmonic = k sin^2(theta0). Masses in amu and nm give ps^-2."""
    i, j, k = trip
    r1, r2 = np.linalg.norm(pos[i] - pos[j]), np.linalg.norm(pos[k] - pos[j])
    d2 = np.sum((pos[i] - pos[k]) ** 2)
    g = 1 / (masses[i] * r1 ** 2) + 1 / (masses[k] * r2 ** 2) + d2 / (masses[j] * r1 ** 2 * r2 ** 2)
    s2 = max(math.sin(math.radians(theta0_deg)) ** 2, 0.05)
    return _omega_limit_sq() / g / s2


def dihedral_k_limit(pos, quad, masses):
    """Largest harmonic dihedral k (kJ/mol/rad2) for the out-of-plane motion of the end beads."""
    i, j, k, l = quad
    u = pos[k] - pos[j]
    u = u / np.linalg.norm(u)

    def perp2(a):
        v = pos[a] - pos[j]
        v = v - (v @ u) * u
        return max(v @ v, 1e-4)
    g = 1 / (masses[i] * perp2(i)) + 1 / (masses[l] * perp2(l))
    return _omega_limit_sq() / g


def read_sigma(core_itp, types):
    """sigma(i, j) in nm for beads of the given types, from [ nonbond_params ] of the core itp."""
    wanted, table, sec = set(types), {}, None
    for line in open(core_itp):
        s = line.split(";")[0].strip()
        if s.startswith("["):
            sec = s.strip("[] ").lower()
            continue
        if sec == "nonbond_params" and s:
            f = s.split()
            if f[0] in wanted and f[1] in wanted:
                table[frozenset((f[0], f[1]))] = float(f[3])

    def sigma(i, j):
        return table[frozenset((types[i], types[j]))]
    return sigma


def bead_graph(M, atom_bonds):
    """Sorted bead pairs joined by an atom bond or a shared atom."""
    owners = [set(np.nonzero(M[:, a])[0].tolist()) for a in range(M.shape[1])]
    edges = set()
    for a, b in atom_bonds:
        for x in owners[a]:
            for y in owners[b]:
                if x != y:
                    edges.add((min(x, y), max(x, y)))
    for o in owners:
        edges.update(combinations(sorted(o), 2))
    return sorted(edges)


SIZE_MASS = {"R": 72.0, "S": 54.0, "T": 36.0}   # martini_v3.0.0.itp [atomtypes]; Souza SI C4 p.47


def bead_masses(types, virtual_sites, groups):
    """Bead masses with virtual-site mass redistributed (Souza SI C4: hinge beads carry the mass).

    types: bead types in bead order. virtual_sites: design.json entries (1-based `site`, `constructing`).
    groups: rigid groups as 0-based bead lists (real and virtual members).
    Each virtual site gets 0; its standard mass is split equally over the real beads of its rigid
    group, or over its constructing beads if it is in no rigid group (CAFF 4 x 63, BZIM 4 x 45).
    Returns (masses, record) where record lists each site and the beads that received its mass.
    """
    std = np.array([SIZE_MASS[t[0] if t[0] in "ST" else "R"] for t in types])
    virtual = {v["site"] - 1 for v in virtual_sites}
    m = std.copy()
    for s in virtual:
        m[s] = 0.0
    record = []
    for v in virtual_sites:
        s = v["site"] - 1
        grp = next((g for g in groups if s in g), None)
        to = [b for b in grp if b not in virtual] if grp else []
        if not to:
            to = [c - 1 for c in v["constructing"] if c - 1 not in virtual]
        for b in to:
            m[b] += std[s] / len(to)
        record.append({"site": s + 1, "mass_amu": float(std[s]), "to_beads": [b + 1 for b in to],
                       "share_amu": round(float(std[s]) / len(to), 6)})
    return m, record


def rigid_groups(mol_h, M, is_virtual):
    """Lists of real beads in one rigid group, with the real+virtual members, merged over shared beads."""
    from rdkit import Chem
    SP2 = Chem.HybridizationType.SP2
    rings = []
    for r in mol_h.GetRingInfo().AtomRings():
        atoms = [mol_h.GetAtomWithIdx(a) for a in r]
        if all(a.GetIsAromatic() for a in atoms) or all(a.GetHybridization() == SP2 for a in atoms):
            rings.append(set(np.nonzero(M[:, list(r)].sum(1))[0].tolist()))
    merged = []
    for s in rings:
        s = set(s)
        rest = []
        for m in merged:
            if m & s:
                s |= m
            else:
                rest.append(m)
        merged = rest + [s]
    return [sorted(s) for s in merged if len(s) >= 2]


def triangulate(pos, beads):
    """Constraint pairs and triangles of a rigid group, in the group's best-fit plane."""
    if len(beads) == 2:
        return [tuple(beads)], []
    if len(beads) == 3:
        return sorted(combinations(beads, 2)), [tuple(beads)]
    if len(beads) == 4:
        # the diagonals are the pairing with the largest summed length; keep the shorter one, and on a tie
        # (a rectangle, as in NAPH) the first in bead order, so the choice is reproducible
        d = lambda e: float(np.linalg.norm(pos[e[0]] - pos[e[1]]))
        a, b, c, e = sorted(beads)
        pairings = [((a, b), (c, e)), ((a, c), (b, e)), ((a, e), (b, c))]
        d1, d2 = max(pairings, key=lambda p: d(p[0]) + d(p[1]))
        keep, other = (d1, d2) if d(d1) <= d(d2) + DIAGONAL_TIE_NM else (d2, d1)
        tris = [tuple(sorted(keep + (other[0],))), tuple(sorted(keep + (other[1],)))]
        return sorted(set(combinations(tris[0], 2)) | set(combinations(tris[1], 2))), tris
    from scipy.spatial import Delaunay
    x = pos[beads]
    x = x - x.mean(0)
    _, _, vt = np.linalg.svd(x, full_matrices=False)
    tri = Delaunay(x @ vt[:2].T, qhull_options="QJ").simplices
    tris = [tuple(sorted(beads[i] for i in t)) for t in tri]

    def edges_of(ts):
        out = defaultdict(int)
        for t in ts:
            for e in combinations(t, 2):
                out[e] += 1
        return out

    # peel hull triangles that only span a concave outline
    changed = True
    while changed and len(tris) > 1:
        changed = False
        cnt = edges_of(tris)
        for t in list(tris):
            long_hull = [e for e in combinations(t, 2) if cnt[e] == 1
                         and np.linalg.norm(pos[e[0]] - pos[e[1]]) > MAX_HULL_EDGE_NM]
            if long_hull:
                tris.remove(t)
                changed = True
                break
    return sorted(edges_of(tris)), tris


def build_terms(mol_h, M, xyz_nm, weights, masses, is_virtual, temperature, bead_types, atom_bonds=None, sigma=None):
    """Returns dict(terms, dropped, exclusions, close_pairs, groups, flags). Bead indices are 0-based.

    M: [n_beads, n_atoms]; xyz_nm: [n_frames, n_atoms, 3]; weights: [n_frames]; bead_types: Martini type names;
    sigma(i, j): Lennard-Jones sigma of a bead pair in nm, or None to skip the close-pair rule.
    """
    close_pairs = []
    is_virtual = np.asarray(is_virtual, bool)
    if atom_bonds is None:
        atom_bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol_h.GetBonds()]
    R = np.einsum("ba,fai->fbi", M, xyz_nm)
    w = np.asarray(weights, float)
    w = w / w.sum()
    # Frames come from different minima in arbitrary orientations, some mirrored, so positions cannot be
    # averaged. Distances can: mean_d chooses neighbours, and one representative frame supplies 3D geometry.
    pick = np.random.default_rng(0).choice(len(w), size=min(2000, len(w)), p=w)
    D = np.linalg.norm(R[pick, :, None] - R[pick, None, :], axis=-1)
    mean_d = D.mean(0)
    rep = R[pick[np.argmin(np.sum((D - mean_d) ** 2, axis=(1, 2)))]]
    kT = F.KB * temperature
    nb = M.shape[0]
    real = [b for b in range(nb) if not is_virtual[b]]

    terms, dropped, flags = [], [], []
    groups = rigid_groups(mol_h, M, is_virtual)
    group_of = {b: gi for gi, g in enumerate(groups) for b in g}
    all_edges = bead_graph(M, atom_bonds)
    graph = [e for e in all_edges if not (is_virtual[e[0]] or is_virtual[e[1]])]
    for x, y in all_edges:           # a bead on a virtual site attaches to the nearest real bead of its group
        if is_virtual[x] == is_virtual[y]:
            continue
        v, j = (x, y) if is_virtual[x] else (y, x)
        if v in group_of and group_of.get(j) != group_of[v]:
            near = min((b for b in groups[group_of[v]] if not is_virtual[b]),
                       key=lambda b: mean_d[b, j])
            graph.append((min(near, j), max(near, j)))
    graph = sorted(set(graph))

    def dist_vals(pair):
        return distances(R, [pair])[:, 0]

    def bond_term(pair, role, force_constraint=False):
        pair = tuple(pair)
        p = F.fit_bond(dist_vals(pair), w, temperature)
        t = dict(section="bonds", beads=list(pair), funct=1, role=role, b0=p["b0"], k=p["kb"],
                 std=float(np.sqrt(kT / p["kb"])))
        cap = bond_k_cap(bead_types[pair[0]], bead_types[pair[1]])
        if force_constraint:
            t.update(section="constraints", k=None, kb_fit=p["kb"], rule="rigid ring")
        elif p["kb"] > cap:
            t.update(k=cap, kb_fit=p["kb"], capped=f"bead-type cap {cap:.0f} kJ/mol/nm2 ({bead_types[pair[0]]}-{bead_types[pair[1]]})")
        elif p["kb"] < BOND_K_FLOOR:
            t.update(k=BOND_K_FLOOR, kb_fit=p["kb"], floored="Souza SI C3 p.44: bonds below 700 kJ/mol avoided")
        return t

    # --- rigid groups
    connected = set()
    tri_of = {}
    for gi, g in enumerate(groups):
        reals = [b for b in g if not is_virtual[b]]
        if len(reals) < 2:
            continue
        if len(reals) > 4:
            flags.append(dict(flag="rigid_group_above_four_beads", beads=reals, note=SIZE_NOTE))
        pairs, tris = triangulate(rep, reals)
        tri_of[gi] = tris
        for tr in tris:
            sides = sorted(float(mean_d[a, b]) for a, b in combinations(tr, 2))
            cosine = (sides[0] ** 2 + sides[1] ** 2 - sides[2] ** 2) / (2 * sides[0] * sides[1])
            widest = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
            if widest > NEAR_LINEAR_TRIANGLE_DEG:
                flags.append(dict(flag="near_linear_rigid_triangle", beads=list(tr), widest_angle_deg=round(widest, 1),
                                  note="a nearly straight constraint triangle is numerically fragile; prefer a mapping with a wider one"))
        for e in pairs:
            terms.append(bond_term(e, "rigid_group", force_constraint=True))
            connected.add(e)
        # one planarity improper across every interior edge
        by_edge = defaultdict(list)
        for t in tris:
            for e in combinations(t, 2):
                by_edge[e].append(t)
        for e, ts in sorted(by_edge.items()):
            if len(ts) == 2:
                a = next(x for x in ts[0] if x not in e)
                d = next(x for x in ts[1] if x not in e)
                terms.append(_improper(R, w, (a, e[0], e[1], d), rep, masses, temperature, "ring_planarity"))

    # --- connections outside rigid groups
    substituents = defaultdict(set)
    for i, j in graph:
        gi, gj = group_of.get(i), group_of.get(j)
        if gi is not None and gi == gj:
            continue
        if gi is not None and gj is None:
            substituents[j].add((gi, i))
        elif gj is not None and gi is None:
            substituents[i].add((gj, j))
        elif (i, j) not in connected:
            terms.append(bond_term((i, j), "chain"))
            connected.add((i, j))

    excl = set()
    for g in groups:
        excl.update(combinations(sorted(g), 2))
    for j, anchors in sorted(substituents.items()):
        by_group = defaultdict(set)
        for gi, a in anchors:
            by_group[gi].add(a)
        for gi, attach in sorted(by_group.items()):
            reals = [b for b in groups[gi] if not is_virtual[b]]
            if gi not in tri_of or len(reals) < 3:
                for a in attach:
                    terms.append(bond_term(sorted((a, j)), "ring_substituent"))
                    connected.add(tuple(sorted((a, j))))
                continue
            attach = sorted(attach, key=lambda a: mean_d[a, j])[:2]
            if len(attach) == 1:          # bridge to the closest group neighbour of the attachment bead
                a = attach[0]
                cons = {e for e in connected if a in e and set(e) <= set(reals)}
                nbrs = [x for e in cons for x in e if x != a]
                attach.append(min(nbrs, key=lambda n: mean_d[n, j]))
            for a in attach:
                terms.append(bond_term(sorted((a, j)), "ring_substituent"))
                connected.add(tuple(sorted((a, j))))
            # planarity improper against the third bead of the triangle holding the bridged edge
            e = tuple(sorted(attach))
            third = [x for t in tri_of[gi] if set(e) <= set(t) for x in t if x not in e]
            if third:
                t = _position_term(R, w, (third[0], e[0], e[1], j), rep, masses, temperature, kT, dropped)
                if t:
                    terms.append(t)
            for x in groups[gi]:
                if x not in attach and x != j:
                    excl.add(tuple(sorted((x, j))))

    # --- angles
    nbr = defaultdict(set)
    for i, j in connected:
        nbr[i].add(j)
        nbr[j].add(i)
    angle_theta0 = {}
    for j in sorted(nbr):
        for i, k in combinations(sorted(nbr[j]), 2):
            if tuple(sorted((i, k))) in connected:
                dropped.append(dict(kind="angle", beads=[i, j, k], reason="i-k connected: triangle fixes it"))
                continue
            if i in group_of and group_of[i] == group_of.get(j) == group_of.get(k):
                dropped.append(dict(kind="angle", beads=[i, j, k], reason="inside a rigid group: constraints and improper fix it"))
                continue
            gs = [group_of.get(b) for b in (i, j, k)]
            big = [g for g in set(gs) if g is not None and gs.count(g) == 2 and g in tri_of and tri_of[g]]
            if big:
                dropped.append(dict(kind="angle", beads=[i, j, k], reason="substituent on a rigid group: bridge bonds and position term fix it"))
                continue
            th = angles(R, [(i, j, k)])[:, 0]
            p = F.fit_angle(th, w, temperature)
            if p["ktheta"] < kT:
                dropped.append(dict(kind="angle", beads=[i, j, k], reason=f"unrestrained (k {p['ktheta']:.2f} below kT)",
                                    theta0=p["theta0"], k=p["ktheta"]))
                continue
            t = dict(section="angles", beads=[i, j, k], funct=2, role="chain", theta0=p["theta0"], k=p["ktheta"])
            angle_theta0[(i, j, k)] = p["theta0"]
            lim = angle_k_limit(rep, (i, j, k), masses, p["theta0"])
            if p["ktheta"] > lim:
                t.update(k=lim, k_fit=p["ktheta"], capped=f"timestep limit: period at least {MIN_PERIOD_STEPS:.0f} steps of {TIMESTEP_PS*1000:.0f} fs")
            terms.append(t)

    # --- dihedrals: best per central edge
    in_group = lambda a, b: a in group_of and group_of.get(a) == group_of.get(b)
    for j, k in sorted(connected):
        if in_group(j, k):
            continue
        cands = []
        for i in sorted(nbr[j] - {k}):
            for l in sorted(nbr[k] - {j, i}):
                q = (i, j, k, l)
                if _within_group_plus_one(q, group_of):
                    dropped.append(dict(kind="dihedral", beads=list(q), reason="rigid group plus one bead: bridge bonds and position term fix it"))
                    continue
                phi, bends = dihedrals(R, [q], return_bends=True)
                phi, bends = phi[:, 0], bends[:, 0, :]
                bad = (bends.max(-1) >= BEND_LIMIT_DEG) | (bends.min(-1) <= BEND_LOW_DEG)
                frac = float(np.sum(w * bad))
                flanks = [tuple(sorted((i, k))) for i, k in ((q[0], q[2]), (q[1], q[3]))]
                th0 = [angle_theta0.get((a, c, b)) for (a, b), c in zip(flanks, (q[1], q[2]))]
                if any(x is not None and x > BEND_LIMIT_DEG for x in th0):
                    dropped.append(dict(kind="dihedral", beads=list(q), reason=f"flanking angle theta0 above {BEND_LIMIT_DEG:.0f} deg (SI C4)"))
                    continue
                if frac > BEND_FRACTION:
                    dropped.append(dict(kind="dihedral", beads=list(q),
                                        reason=f"flanking angle near 0 or 180 deg in {frac:.0%} of ensemble (SI C4)"))
                    continue
                ww = w * ~bad
                p = F.fit_dihedral(phi, ww, temperature)
                g = np.linspace(-180, 180, 361)
                V = F.potential(p, g)
                cands.append((float(V.max() - V.min()), q, p))
        if not cands:
            continue
        cands.sort(key=lambda c: -c[0])
        rng, q, p = cands[0]
        for c in cands[1:]:
            dropped.append(dict(kind="dihedral", beads=list(c[1]), reason="another dihedral on this central edge has a larger range"))
        if rng < kT:
            dropped.append(dict(kind="dihedral", beads=list(q), reason=f"flat: potential range {rng:.2f} kJ/mol below kT"))
            continue
        t = dict(section="dihedrals", beads=list(q), role="torsion", range_kJ_mol=rng)
        if p["form"] == "harmonic":
            lim = dihedral_k_limit(rep, q, masses)
            t.update(funct=2, phi0=p["phi0"], k=min(p["k"], lim))
            if p["k"] > lim:
                t.update(k_fit=p["k"], capped="timestep limit")
        else:
            t.update(funct=9, fourier=[x for x in p["terms"] if x["k"] >= F.FOURIER_MIN_K])
        terms.append(t)

    # 2nd and 3rd neighbours whose mapped distance is inside their repulsive wall (Souza SI C3 p.44):
    # at nrexcl 1 the CG non-bonded interaction would hold them apart and distort the fitted terms.
    if sigma is not None:
        hop = {}
        for s in range(nb):
            seen, frontier, level = {s}, {s}, 0
            while frontier and level < 3:
                level += 1
                frontier = {y for x in frontier for y in nbr[x]} - seen
                seen |= frontier
                for y in frontier:
                    hop[(s, y)] = level
        for (x, y), level in sorted(hop.items()):
            if x < y and level in (2, 3) and (x, y) not in excl and mean_d[x, y] < sigma(x, y):
                excl.add((x, y))
                close_pairs.append(dict(beads=[x, y], neighbour_order=level, mean_distance_nm=round(float(mean_d[x, y]), 4),
                                        sigma_nm=round(float(sigma(x, y)), 4)))
    excl -= connected
    # two beads on one ring: COG distances are short, and published models refine them against SASA/volume
    for r in mol_h.GetRingInfo().AtomRings():
        owners = sorted(set(np.nonzero(M[:, list(r)].sum(1))[0].tolist()))
        if len(owners) == 2 and tuple(owners) in connected:
            flags.append(dict(flag="two_bead_ring_cog_length", beads=owners, length_nm=round(float(mean_d[owners[0], owners[1]]), 3),
                              note="COG length is not volume-refined (Souza SI C3 p.42); check SASA or density before production"))
    return dict(terms=terms, dropped=dropped, exclusions=sorted(excl), close_pairs=close_pairs, groups=groups, flags=flags)


def _within_group_plus_one(q, group_of):
    groups = [group_of.get(b) for b in q]
    return any(g is not None and groups.count(g) >= len(q) - 1 for g in set(groups))


def _position_term(R, w, quad, pos, masses, temperature, kT, dropped):
    """Out-of-plane position of a substituent: harmonic if narrow; Fourier if it has several states."""
    phi, bends = dihedrals(R, [quad], return_bends=True)
    phi, bends = phi[:, 0], bends[:, 0, :]
    bad = (bends.max(-1) >= BEND_LIMIT_DEG) | (bends.min(-1) <= BEND_LOW_DEG)
    frac = float(np.sum(w * bad))
    if frac > BEND_FRACTION:
        dropped.append(dict(kind="dihedral", beads=list(quad), reason=f"flanking angle near 0 or 180 deg in {frac:.0%} of ensemble (SI C4)"))
        return None
    ww = w * ~bad
    p = F.fit_dihedral(phi, ww, temperature)
    if p["form"] == "harmonic":
        return _improper(R, ww / ww.sum(), quad, pos, masses, temperature, "substituent_position")
    V = F.potential(p, np.linspace(-180, 180, 361))
    rng = float(V.max() - V.min())
    if rng < kT:
        dropped.append(dict(kind="dihedral", beads=list(quad), reason=f"flat: potential range {rng:.2f} kJ/mol below kT"))
        return None
    return dict(section="dihedrals", beads=list(quad), funct=9, role="substituent_position", range_kJ_mol=rng,
                fourier=[x for x in p["terms"] if x["k"] >= F.FOURIER_MIN_K])


def _improper(R, w, quad, pos, masses, temperature, role):
    phi = dihedrals(R, [quad])[:, 0]
    mu, sd = F._circ_stats(np.radians(phi), w)
    k = F.KB * temperature / max(sd, 1e-3) ** 2
    lim = dihedral_k_limit(pos, quad, masses)
    t = dict(section="dihedrals", beads=list(quad), funct=2, role=role, phi0=math.degrees(mu), k=min(k, lim),
             std_deg=math.degrees(sd))
    if k > lim:
        t.update(k_fit=k, capped="timestep limit")
    return t
