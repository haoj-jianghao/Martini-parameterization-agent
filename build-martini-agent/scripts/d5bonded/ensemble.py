"""D5 ensemble: harmonic superposition approximation (HSA) around MMFF94 minima.

Classical configurational distribution, so masses do not enter:
  p(x) = sum_i  W_i N(x; x_i, kT H_i^+)          (internal coordinates; rigid-body modes projected out)
  W_i  ∝ exp(-E_i/kT) * sqrt(det I_i) * prod_k (2 pi kT / lambda_ik)^(1/2)
x_i: unique MMFF minima, E_i their energies, lambda_ik the 3N-6 nonzero eigenvalues of
the Hessian H_i projected onto internal motions, and I_i the unit-mass geometric
inertia tensor (the rotational part of the configurational integral). Permutation
symmetry numbers are ignored: equivalent minima are merged, and all minima of
a molecule share its symmetry.

Samples: n_samples per minimum, frame weight W_i / n_samples. Straight-line
(Cartesian) displacement along soft torsional modes stretches bonds badly (a
C=O bond came out 0.134 +/- 0.019 nm instead of ~0.121 +/- 0.003), so the
internal space is split into
  T: generators of rotation about every rotatable bond (acyclic single bonds with
     substituents on both ends; methyl/OH/NH2 rotors included), each rotating the
     smaller side; and
  S: the orthonormal complement (stiff motions).
In the basis C = [T, S] the harmonic covariance is kT (C^T H C)^-1 (torsion-stiff
coupling kept). Each sample applies the S part as a small linear displacement,
then the torsion angles as exact finite rotations, so bonds and angles are not
distorted. Soft ring modes (e.g. cyclopentane pseudorotation) are not rotors and
would stretch ring bonds when applied linearly, so each eigenmode of the stiff block
is capped at max_stiff_std_A (default 0.3 A; real stiff vibrations are 0.03-0.1 A),
and any sample with a bond off by more than 0.3 A from the minimum is rejected
(counted in info["rejected_frac"]). Torsion std is capped at max_torsion_std_deg (default 60 deg) so the
harmonic well of one minimum does not spill into its neighbours (they are separate
minima). For achiral molecules every sample is also added reflected (z -> -z) with
the same weight, which mirror-symmetrizes the dihedrals exactly.

energy_model (diagnostic): an object with energy/grad/minimize in A and kcal/mol (e.g.
d5bonded.energy.NumpyEnergy for the AA OPLS + OBC2 Hamiltonian). ETKDG + MMFF stay the
conformer search; every unique MMFF minimum is re-minimized with the energy model, duplicates
are merged again, and E_i, H_i (so the weights and the sampling covariance) come from it.
If the model has .hessian(x) (e.g. hsa.xtb_energy.XTBEnergy: GFN2-xTB / GFN-FF + ALPB via the
xtb program, no MD engine) that is used; otherwise finite differences of .grad.
refine_saddles (external models only): a minimum whose Hessian has a negative mode is displaced
by refine_step_A along it both ways, re-minimized (lower end kept) and rechecked, up to this many
times; refined minima are merged again with the others. MMFF minima are never refined (baseline).
hess_e_window (kcal/mol): minima more than this above the lowest one get no Hessian and no weight
(saves cost for expensive models; at 300 K, 8 kcal/mol is a Boltzmann factor of ~1e-6).

Units: Angstrom and kcal/mol internally (RDKit); the trajectory is returned in nm.
"""
import time
from dataclasses import dataclass, field

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from rdkit.Chem import rdForceFieldHelpers as FF

KB_KCAL = 0.0019872043  # kcal mol-1 K-1
COLLAPSE_DISTANCE_A = 1.0  # non-bonded atoms closer than this are a force-field collapse, not a contact


@dataclass
class D5Ensemble:
    xyz_nm: np.ndarray            # [n_frames, n_atoms, 3]
    weights: np.ndarray           # [n_frames], sum 1
    atoms: list
    minima: list                  # per unique minimum: dict(energy, weight, n_neg, ...)
    info: dict = field(default_factory=dict)

    def as_trajectory(self):
        """Duck-typed cglabels Trajectory (no box)."""
        class T:
            pass
        t = T()
        n = self.xyz_nm.shape[0]
        t.xyz, t.weights, t.atoms, t.box = self.xyz_nm, self.weights, self.atoms, None
        t.steps, t.time_ps, t.run_index = np.arange(n), np.zeros(n), np.zeros(n, int)
        t.meta, t.n_frames = [dict(self.info, solvent="none (MMFF94)", run_dir=None)], n
        return t


def _rigid_basis(x):
    """Orthonormal basis [3N, 6] of rigid translations and rotations at geometry x [N, 3]."""
    n = len(x)
    c = x - x.mean(0)
    B = []
    for k in range(3):
        t = np.zeros((n, 3)); t[:, k] = 1; B.append(t.ravel())
    for k in range(3):
        e = np.zeros(3); e[k] = 1
        B.append(np.cross(e, c).ravel())
    Q, _ = np.linalg.qr(np.array(B).T)
    return Q


def _grad(ff, x):
    """Analytic MMFF gradient at x. RDKit (2025.09) ForceField.CalcGrad(pos) partly reuses
    state from the previous CalcEnergy call, so it is wrong unless CalcEnergy(pos) is
    called on the same positions first (verified against finite differences of the energy)."""
    p = x.tolist()
    ff.CalcEnergy(p)
    return np.array(ff.CalcGrad(p))


def _hessian(ff, x, h=1e-4):
    """Central finite differences of the analytic gradient (kcal/mol/A^2). ff: RDKit MMFF
    force field, or an energy model with .grad(x) (see energy.py)."""
    x = np.asarray(x, float).ravel()
    n = len(x)
    g = ff.grad if hasattr(ff, "grad") else (lambda v: _grad(ff, v))
    H = np.empty((n, n))
    for k in range(n):
        xp = x.copy(); xp[k] += h
        xm = x.copy(); xm[k] -= h
        H[k] = (g(xp) - g(xm)) / (2 * h)
    return 0.5 * (H + H.T)


def _rotors(mol):
    """[(axis_from, axis_to, moving_atoms)] for rotatable bonds; the smaller side moves."""
    out = []
    n = mol.GetNumAtoms()
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.SINGLE or b.IsInRing():
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if mol.GetAtomWithIdx(i).GetDegree() < 2 or mol.GetAtomWithIdx(j).GetDegree() < 2:
            continue
        side, stack = {j}, [j]                       # atoms reachable from j without crossing i-j
        while stack:
            a = stack.pop()
            for nb in mol.GetAtomWithIdx(a).GetNeighbors():
                k = nb.GetIdx()
                if k not in side and not (a == j and k == i):
                    side.add(k); stack.append(k)
        if i in side:
            continue                                  # ring (should not happen for acyclic bonds)
        if len(side) <= n - len(side):
            out.append((i, j, np.array(sorted(side))))
        else:
            out.append((j, i, np.array(sorted(set(range(n)) - side))))
    return out


def _rotate(x, a, b, moving, angle):
    """Rotate atoms `moving` of x [N,3] about the axis a->b by angle (rad), in place."""
    u = x[b] - x[a]; u = u / np.linalg.norm(u)
    p = x[moving] - x[b]
    c, s_ = np.cos(angle), np.sin(angle)
    x[moving] = x[b] + p * c + np.cross(u, p) * s_ + np.outer(p @ u, u) * (1 - c)


def _torsion_generators(x, rotors):
    """[3N, n_rot]: displacement per radian of rotation about each rotor axis."""
    T = np.zeros((x.size, len(rotors)))
    for r, (a, b, mov) in enumerate(rotors):
        u = x[b] - x[a]; u = u / np.linalg.norm(u)
        d = np.zeros_like(x); d[mov] = np.cross(u, x[mov] - x[b])
        T[:, r] = d.ravel()
    return T


def _is_achiral(mol):
    m = Chem.RemoveHs(mol)
    mirror = Chem.Mol(m)
    for a in mirror.GetAtoms():
        t = a.GetChiralTag()
        if t == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif t == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            a.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    return Chem.MolToSmiles(m) == Chem.MolToSmiles(mirror)


def d5_ensemble(smiles, n_conf=300, temperature=300.0, n_samples=200, seed=20261103,
                dielectric="distance", eps=4.0, e_tol=1e-3, rms_tol=0.1, max_torsion_std_deg=60.0,
                max_stiff_std_A=0.3, mirror=True, return_minima_only=False, energy_model=None,
                hess_e_window=None, refine_saddles=2, refine_step_A=0.3):
    """dielectric: "constant" (eps) or "distance" (eps*r, MMFF distance-dependent screening)."""
    t0 = time.perf_counter()
    kT = KB_KCAL * temperature
    rng = np.random.default_rng(seed)
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3()
    p.randomSeed, p.pruneRmsThresh, p.numThreads = seed, -1.0, 1
    cids = list(AllChem.EmbedMultipleConfs(mol, n_conf, p))
    if not cids:
        raise RuntimeError(f"ETKDG failed for {smiles}")
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise RuntimeError(f"no MMFF94 parameters for {smiles}")
    props.SetMMFFDielectricModel(2 if dielectric == "distance" else 1)  # RDKit: 1 constant, 2 distance-dependent
    props.SetMMFFDielectricConstant(eps)
    energies = []
    for c in cids:
        ff = FF.MMFFGetMoleculeForceField(mol, props, confId=c)
        for _ in range(20):          # minimize to convergence (Minimize returns 0 when converged)
            if ff.Minimize(maxIts=2000, forceTol=1e-6, energyTol=1e-12) == 0:
                break
        energies.append(ff.CalcEnergy())
    energies = np.array(energies)
    # unique minima: energy within e_tol kcal/mol and heavy-atom best RMSD (symmetry-aware, no reflection) < rms_tol
    heavy = Chem.RemoveHs(mol)

    def dedup(cands):
        out = []
        for c in sorted(cands, key=lambda c: energies[c]):
            dup = False
            for u in out:
                if abs(energies[c] - energies[u]) < e_tol and \
                        rdMolAlign.GetBestRMS(heavy, heavy, prbId=cids[u], refId=cids[c]) < rms_tol:
                    dup = True; break
            if not dup:
                out.append(c)
        return out
    uniq = dedup(range(len(cids)))
    n_unique_mmff = len(uniq)
    n_collapsed = 0
    if energy_model is not None:     # re-minimize the MMFF minima with the external model, merge again
        from rdkit.Geometry import Point3D
        far = Chem.GetDistanceMatrix(mol) >= 3
        kept = []
        for c in uniq:
            conf = mol.GetConformer(cids[c])
            x = energy_model.minimize(conf.GetPositions())
            d = np.linalg.norm(x[:, None] - x[None], axis=-1)
            if d[far].min() < COLLAPSE_DISTANCE_A or not np.isfinite(x).all():
                n_collapsed += 1      # a sigma-0 hydrogen fell onto a charged atom: unphysical, energy diverges
                continue
            kept.append(c)
            for a, (px, py, pz) in enumerate(x):
                conf.SetAtomPosition(a, Point3D(float(px), float(py), float(pz)))
            energies[c] = energy_model.energy(x)
        if not kept:
            raise RuntimeError(f"all {n_collapsed} minima collapsed under {getattr(energy_model, 'name', 'the energy model')}")
        heavy = Chem.RemoveHs(mol)   # refresh heavy-atom coordinates from the updated conformers
        uniq = dedup(kept)
    minima = []
    e_min = min(energies[c] for c in uniq)
    n_refine_max = refine_saddles if (energy_model is not None and refine_saddles is not None) else 0
    n_refined = n_rescued = 0

    def modes(x, H):
        Q = _rigid_basis(x)
        P = np.eye(len(H)) - Q @ Q.T
        lam, V = np.linalg.eigh(P @ H @ P)
        inside = np.linalg.norm(Q.T @ V, axis=0) > 0.5     # drop the 6 rigid-body directions
        return Q, lam[~inside], V[:, ~inside]

    for c in uniq:
        x = mol.GetConformer(cids[c]).GetPositions()
        if hess_e_window is not None and energies[c] - e_min > hess_e_window:
            # far above the lowest minimum: negligible weight, skip its (expensive) Hessian
            minima.append(dict(conf=int(cids[c]), energy=float(energies[c]), n_neg=0, used=False,
                               skipped="e_window"))
            continue
        ff = energy_model or FF.MMFFGetMoleculeForceField(mol, props, confId=cids[c])
        H = ff.hessian(x) if hasattr(ff, "hessian") else _hessian(ff, x)
        Q, lam, V = modes(x, H)
        tries = 0
        while (lam <= 1e-6).any() and tries < n_refine_max:
            # stopped on a saddle / flat shoulder: step along the most negative mode both ways,
            # re-minimize, keep the lower end point, recheck
            v = V[:, np.argmin(lam)].reshape(-1, 3)
            v = refine_step_A * v / np.linalg.norm(v)
            ends = [energy_model.minimize(x + s * v) for s in (1.0, -1.0)]
            es = [energy_model.energy(y) for y in ends]
            x = ends[int(np.argmin(es))]; energies[c] = min(es)
            H = energy_model.hessian(x) if hasattr(energy_model, "hessian") else _hessian(energy_model, x)
            Q, lam, V = modes(x, H)
            tries += 1
        if tries:
            n_refined += 1
            from rdkit.Geometry import Point3D
            conf = mol.GetConformer(cids[c])
            for a, (px, py, pz) in enumerate(x):
                conf.SetAtomPosition(a, Point3D(float(px), float(py), float(pz)))
        n_neg = int((lam <= 1e-6).sum())
        if n_neg:                    # not a true minimum (saddle or flat): skip it
            minima.append(dict(conf=int(cids[c]), energy=float(energies[c]), n_neg=n_neg, used=False,
                               lam_min=float(lam.min()), refine_tries=tries))
            continue
        n_rescued += bool(tries)
        cen = x - x.mean(0)
        I = np.eye(3) * (cen ** 2).sum() - cen.T @ cen
        lw = -energies[c] / kT + 0.5 * np.log(np.linalg.det(I)) + 0.5 * np.sum(np.log(2 * np.pi * kT / lam))
        minima.append(dict(conf=int(cids[c]), energy=float(energies[c]), n_neg=0, used=True,
                           x=x, H=H, Q=Q, logw=lw, refine_tries=tries))
    if n_refined:                    # refined minima may have fallen into an existing one: merge again
        heavy = Chem.RemoveHs(mol)
        idx = {int(cid): c for c, cid in enumerate(cids)}
        keep = {int(cids[c]) for c in dedup([idx[m["conf"]] for m in minima if m["used"]])}
        for m in minima:
            if m["used"] and m["conf"] not in keep:
                m["used"] = False; m["skipped"] = "dup_after_refine"
    used = [m for m in minima if m["used"]]
    if not used:
        low = sorted((m for m in minima if m["n_neg"]), key=lambda m: m["energy"])[:5]
        raise RuntimeError(f"no true minima for {smiles}: {len(minima)} candidates, lowest with negative modes: "
                           + str([(round(m["energy"], 2), m["n_neg"], round(m["lam_min"], 5)) for m in low]))
    logW = [m["logw"] for m in used]
    logW = np.array(logW)
    W = np.exp(logW - logW.max()); W /= W.sum()
    for m, w in zip(used, W):
        m["weight"] = float(w)
    if return_minima_only:
        return minima
    frames, wts = [], []
    n_atoms = mol.GetNumAtoms()
    rotors = _rotors(mol)
    max_t = np.radians(max_torsion_std_deg)
    bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()])
    n_rej = n_tot = 0
    for m in used:
        x0, H, Q = m["x"], m["H"], m["Q"]
        T = _torsion_generators(x0, rotors)
        T = T - Q @ (Q.T @ T)                             # remove rigid-body parts
        keep = np.linalg.norm(T, axis=0) > 1e-6           # rotors about linear groups do nothing
        T, rot = T[:, keep], [r for r, k in zip(rotors, keep) if k]
        # stiff complement S: orthonormal basis of the internal space orthogonal to span(T)
        U, sv, _ = np.linalg.svd(np.hstack([Q, T]), full_matrices=True)
        rank = int((sv > 1e-8 * sv.max()).sum())
        S = U[:, rank:]
        C = np.hstack([T, S])
        Hc = C.T @ H @ C
        cov = kT * np.linalg.inv(0.5 * (Hc + Hc.T))
        nt = T.shape[1]
        sd_t = np.sqrt(np.diag(cov)[:nt])
        scale = np.minimum(1.0, max_t / np.maximum(sd_t, 1e-12))   # cap torsion std, keep correlations
        Lc = np.linalg.cholesky(cov)
        z = rng.standard_normal((n_samples, len(cov))) @ Lc.T
        z[:, :nt] *= scale
        # cap soft modes of the stiff block along the eigenvectors of its marginal covariance
        ls, Vs = np.linalg.eigh(cov[nt:, nt:])
        cs = np.minimum(1.0, max_stiff_std_A / np.sqrt(np.maximum(ls, 1e-30)))
        z[:, nt:] = ((z[:, nt:] @ Vs) * cs) @ Vs.T
        out = np.empty((n_samples, n_atoms, 3))
        for k in range(n_samples):
            y = x0 + (S @ z[k, nt:]).reshape(n_atoms, 3)
            for (a, b, mov), ang in zip(rot, z[k, :nt]):
                _rotate(y, a, b, mov, ang)
            out[k] = y
        L0 = np.linalg.norm(x0[bonds[:, 1]] - x0[bonds[:, 0]], axis=-1)
        L = np.linalg.norm(out[:, bonds[:, 1]] - out[:, bonds[:, 0]], axis=-1)
        ok = np.abs(L - L0).max(-1) < 0.3
        n_rej += int((~ok).sum()); n_tot += n_samples
        if ok.any():
            frames.append(out[ok])
            wts.append(np.full(int(ok.sum()), m["weight"] / ok.sum()))
        m["n_stiff_capped"] = int((cs < 1).sum())
        m["torsion_std_deg"] = np.degrees(sd_t).round(1).tolist()
    xyz = np.concatenate(frames); w = np.concatenate(wts); w = w / w.sum()
    achiral = _is_achiral(mol)
    if mirror and achiral:
        xyz = np.concatenate([xyz, xyz * np.array([1.0, 1.0, -1.0])])
        w = np.concatenate([w, w]) / 2
    info = dict(method="D5 harmonic superposition", smiles=smiles, n_conf=len(cids), n_unique=len(uniq),
                n_used=len(used), temperature_K=temperature, dielectric=f"{dielectric} eps={eps}",
                n_samples_per_min=n_samples, max_torsion_std_deg=max_torsion_std_deg,
                n_rotors=len(rotors), max_stiff_std_A=max_stiff_std_A, n_unique_mmff=n_unique_mmff,
                energy_model=getattr(energy_model, "name", "MMFF94"), hess_e_window=hess_e_window,
                n_e_window_skipped=sum(m.get("skipped") == "e_window" for m in minima),
                n_refined=n_refined, n_rescued=n_rescued, n_collapsed=n_collapsed,
                n_dup_after_refine=sum(m.get("skipped") == "dup_after_refine" for m in minima),
                n_saddle_rejected=sum(m["n_neg"] > 0 for m in minima),
                rejected_frac=n_rej / max(n_tot, 1), achiral=achiral,
                mirrored=bool(mirror and achiral), top_weight=float(W.max()),
                n_eff_minima=float(1 / np.sum(W ** 2)), wall_s=time.perf_counter() - t0)
    for m in minima:
        for k in ("x", "H", "Q"):
            m.pop(k, None)
    return D5Ensemble(xyz_nm=xyz / 10.0, weights=w, atoms=[a.GetSymbol() for a in mol.GetAtoms()],
                      minima=minima, info=info)
