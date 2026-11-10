"""OPLS-AA + GBSA-OBC2 single-molecule energy model in NumPy (energies, analytic gradients, minimization; no MD).

The stage-4 energy model, read directly from the opls-ff-gen GROMACS topology:
  bonded       bonds (funct 1), angles (funct 1), proper/improper dihedrals: Ryckaert-Bellemans (funct 3) and
               periodic (funct 1, 4, 9)
  nonbonded    Coulomb and Lennard-Jones over all pairs beyond nrexcl bonds, no cutoff; comb-rule 3 (geometric
               sigma and epsilon); [ pairs ] (1-4) scaled by fudgeQQ and fudgeLJ (gen-pairs)
  polar H      sigma = epsilon = 0 hydrogens bonded to N, O or S get POLAR_H_LJ, in the pair LJ and the 1-4 LJ
               (opls.POLAR_H_LJ); engine files are not changed
  GBSA-OBC2    Onufriev-Bashford-Case II Born radii (alpha 1, beta 0.8, gamma 4.85, offset 0.009 nm), mbondi2
               radii and OBC screening factors, solute dielectric 1, solvent 78.5, all pairs without exclusions,
               plus the ACE nonpolar term (2.25936 kJ/mol/nm^2, probe 0.14 nm)
These match the OpenMM 8.x defaults (GBSAOBCForce, customgbforces.GBSAOBC2Force) that the workflow used before
2026-11-07; energies and forces agree with OpenMM to round-off on 426 engine topologies
(runs/openmm_parity-20261107/report.md). OpenMM is not needed.

Internal units are nm and kJ/mol. The energy-model interface (energy, grad, minimize) takes and returns Angstrom,
kcal/mol and coordinates in RDKit AddHs order.
"""
import re
from pathlib import Path

import numpy as np

KJ_TO_KCAL = 1 / 4.184
ONE_4PI_EPS0 = 138.93545764438198  # kJ mol^-1 nm e^-2: 1/(4 pi eps0), CODATA 2018
OBC = (1.0, 0.8, 4.85)             # OBC type II alpha, beta, gamma
GB_OFFSET = 0.009                  # nm, dielectric offset
SOLUTE_DIELECTRIC, SOLVENT_DIELECTRIC = 1.0, 78.5
ACE_ENERGY = 2.25936               # kJ mol^-1 nm^-2
ACE_PROBE = 0.14                   # nm
# mbondi2 radii (nm) and OBC screening factors by atomic number; H on N is 0.13 nm. Unlisted: 0.15 nm, 0.8.
MBONDI2 = {1: 0.12, 6: 0.17, 7: 0.155, 8: 0.15, 9: 0.15, 14: 0.21, 15: 0.185, 16: 0.18, 17: 0.17}
SCREEN = {1: 0.85, 6: 0.72, 7: 0.79, 8: 0.85, 9: 0.88, 15: 0.86, 16: 0.96}
MASS_TO_Z = {1: 1, 12: 6, 14: 7, 16: 8, 19: 9, 28: 14, 31: 15, 32: 16, 35: 17, 80: 35, 127: 53}
POLAR_H_PARTNERS = (7, 8, 16)
PRECOND_EPS = 0.1                  # minimizer preconditioner regularization, fraction of the mean diagonal
WARMUP_ITER = 30                   # plain L-BFGS-B iterations before preconditioning, to settle in the start's basin
POLISH_GMAX = 1e-4                 # kcal/mol/A; above this after the preconditioned run, polish in plain coordinates

SUPPORTED = {"defaults", "atomtypes", "moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals",
             "exclusions", "system", "molecules"}


class TopologyError(RuntimeError):
    pass


def _read_lines(path, seen=None):
    """Lines of a GROMACS topology with #include expanded (relative to the including file); no #define/#ifdef."""
    path = Path(path).resolve()
    seen = seen or set()
    if path in seen:
        raise TopologyError(f"recursive #include of {path}")
    seen = seen | {path}
    out = []
    for line in path.read_text().splitlines():
        s = line.split(";")[0].strip()
        if s.startswith("#include"):
            inc = s.split(None, 1)[1].strip().strip('"<>')
            p = path.parent / inc
            if not p.exists():
                raise TopologyError(f"#include {inc} not found next to {path.name}")
            out += _read_lines(p, seen)
        elif s.startswith("#"):
            raise TopologyError(f"preprocessor directive not supported: {s}")
        else:
            out.append(s)
    return out


class GmxTopology:
    """One molecule (count 1) from a GROMACS topology, 0-based indices."""

    def __init__(self, top_path):
        self.path = Path(top_path)
        types, sec = {}, None
        self.defaults = None
        self.nrexcl = None
        self.charge, self.mass, self.z, self.sigma, self.eps = [], [], [], [], []
        self.bonds, self.pairs, self.angles, self.rb, self.periodic, self.exclusions = [], [], [], [], [], []
        molecules = []
        n_moltypes = 0
        for s in _read_lines(top_path):
            if not s:
                continue
            m = re.match(r"^\[\s*(\w+)\s*\]$", s)
            if m:
                sec = m.group(1).lower()
                if sec not in SUPPORTED:
                    raise TopologyError(f"section [ {sec} ] not supported")
                n_moltypes += sec == "moleculetype"
                continue
            f = s.split()
            if sec == "defaults":
                self.defaults = dict(nbfunc=int(f[0]), comb_rule=int(f[1]), gen_pairs=f[2].lower() == "yes",
                                     fudge_lj=float(f[3]), fudge_qq=float(f[4]))
            elif sec == "atomtypes":
                k = f.index("A") if "A" in f else len(f) - 3      # ptype column; sigma and epsilon follow it
                types[f[0]] = (float(f[k + 1]), float(f[k + 2]), float(f[k - 2]))   # sigma, eps, mass
            elif sec == "moleculetype":
                self.nrexcl = int(f[1])
            elif sec == "atoms":
                t = types[f[1]]
                mass = float(f[7]) if len(f) > 7 else t[2]
                self.charge.append(float(f[6]))
                self.mass.append(mass)
                self.z.append(MASS_TO_Z[round(mass)])
                self.sigma.append(t[0])
                self.eps.append(t[1])
            elif sec == "bonds":
                self._need(f, 2, (1,), "bonds")
                self.bonds.append((int(f[0]) - 1, int(f[1]) - 1, float(f[3]), float(f[4])))
            elif sec == "pairs":
                self._need(f, 2, (1,), "pairs")
                if len(f) > 3:
                    raise TopologyError("explicit [ pairs ] parameters not supported (gen-pairs expected)")
                self.pairs.append((int(f[0]) - 1, int(f[1]) - 1))
            elif sec == "angles":
                self._need(f, 3, (1,), "angles")
                self.angles.append((int(f[0]) - 1, int(f[1]) - 1, int(f[2]) - 1, np.radians(float(f[4])), float(f[5])))
            elif sec == "dihedrals":
                fn = int(f[4])
                ijkl = tuple(int(x) - 1 for x in f[:4])
                if fn == 3:
                    self.rb.append((*ijkl, *[float(x) for x in f[5:11]]))
                elif fn in (1, 4, 9):
                    self.periodic.append((*ijkl, np.radians(float(f[5])), float(f[6]), int(f[7])))
                else:
                    raise TopologyError(f"dihedral funct {fn} not supported")
            elif sec == "exclusions":
                a = int(f[0]) - 1
                self.exclusions += [(a, int(x) - 1) for x in f[1:]]
            elif sec == "molecules":
                molecules.append((f[0], int(f[1])))
        d = self.defaults
        if d is None or d["nbfunc"] != 1 or d["comb_rule"] != 3 or not d["gen_pairs"]:
            raise TopologyError(f"expected [ defaults ] nbfunc 1, comb-rule 3, gen-pairs yes; got {d}")
        if n_moltypes != 1 or (molecules and (len(molecules) != 1 or molecules[0][1] != 1)):
            raise TopologyError("expected one moleculetype, one molecule")
        self.n = len(self.charge)
        self.charge, self.mass, self.z = np.array(self.charge), np.array(self.mass), np.array(self.z, int)
        self.sigma, self.eps = np.array(self.sigma), np.array(self.eps)

    @staticmethod
    def _need(f, natoms, functs, sec):
        if int(f[natoms]) not in functs:
            raise TopologyError(f"[ {sec} ] funct {f[natoms]} not supported")

    def neighbours(self):
        adj = [set() for _ in range(self.n)]
        for i, j, *_ in self.bonds:
            adj[i].add(j)
            adj[j].add(i)
        return adj

    def excluded(self):
        """Pairs within nrexcl bonds, plus [ exclusions ], as a boolean matrix (diagonal True)."""
        adj = self.neighbours()
        ex = np.eye(self.n, dtype=bool)
        for i in range(self.n):
            seen, front = {i}, {i}
            for _ in range(self.nrexcl):
                front = {y for x in front for y in adj[x]} - seen
                seen |= front
            ex[i, list(seen)] = True
        for i, j in self.exclusions:
            ex[i, j] = ex[j, i] = True
        return ex


# ------------------------------------------------------------------ geometry with gradients
def _angles(x, i, j, k):
    """theta and d(theta)/d(x_i, x_j, x_k)."""
    u, v = x[i] - x[j], x[k] - x[j]
    lu, lv = np.linalg.norm(u, axis=1), np.linalg.norm(v, axis=1)
    c = np.clip(np.einsum("ij,ij->i", u, v) / (lu * lv), -1.0, 1.0)
    th = np.arccos(c)
    s = np.sqrt(np.maximum(1 - c * c, 1e-24))
    dci = (v / lv[:, None] - c[:, None] * u / lu[:, None]) / lu[:, None]
    dck = (u / lu[:, None] - c[:, None] * v / lv[:, None]) / lv[:, None]
    di, dk = -dci / s[:, None], -dck / s[:, None]
    return th, di, -(di + dk), dk


def _dihedrals(x, i, j, k, l):
    """IUPAC phi (trans = pi) and d(phi)/d(x_i, x_j, x_k, x_l)."""
    b1, b2, b3 = x[j] - x[i], x[k] - x[j], x[l] - x[k]
    m, n = np.cross(b1, b2), np.cross(b2, b3)
    lb2 = np.linalg.norm(b2, axis=1)
    phi = np.arctan2(lb2 * np.einsum("ij,ij->i", b1, n), np.einsum("ij,ij->i", m, n))
    mm, nn = np.einsum("ij,ij->i", m, m), np.einsum("ij,ij->i", n, n)
    di = -(lb2 / mm)[:, None] * m
    dl = (lb2 / nn)[:, None] * n
    p = (np.einsum("ij,ij->i", b1, b2) / lb2 ** 2)[:, None]
    q = (np.einsum("ij,ij->i", b3, b2) / lb2 ** 2)[:, None]
    dj = (-p - 1) * di + q * dl
    dk = (-q - 1) * dl + p * di
    return phi, di, dj, dk, dl


class OplsObc2:
    """Energy and analytic gradient (kJ/mol, nm) of one molecule from a GmxTopology."""

    def __init__(self, top, polar_h_lj=None, solvent="obc2"):
        if solvent not in ("obc2", "vacuum"):
            raise ValueError(solvent)
        self.top, self.solvent = top, solvent
        n = top.n
        sigma, eps = top.sigma.copy(), top.eps.copy()
        adj = top.neighbours()
        self.polar_h = sorted(h for h in range(n) if top.z[h] == 1 and any(top.z[x] in POLAR_H_PARTNERS for x in adj[h]))
        self.polar_h_changed = []
        if polar_h_lj:
            for h in self.polar_h:
                if sigma[h] == 0 and eps[h] == 0:
                    sigma[h], eps[h] = polar_h_lj
                    self.polar_h_changed.append(h)
        q = top.charge
        se = np.sqrt(eps)
        # all non-excluded pairs i < j
        ex = top.excluded()
        iu, ju = np.triu_indices(n, 1)
        keep = ~ex[iu, ju]
        self.nb_i, self.nb_j = iu[keep], ju[keep]
        self.nb_qq = ONE_4PI_EPS0 * q[self.nb_i] * q[self.nb_j]
        s_ij = np.sqrt(sigma[self.nb_i] * sigma[self.nb_j])
        e_ij = se[self.nb_i] * se[self.nb_j]
        self.nb_a, self.nb_c = 4 * e_ij * s_ij ** 12, 4 * e_ij * s_ij ** 6
        # 1-4 pairs, generated from the atom types
        d = top.defaults
        pr = np.array(top.pairs, int).reshape(-1, 2)
        self.p_i, self.p_j = pr[:, 0], pr[:, 1]
        self.p_qq = d["fudge_qq"] * ONE_4PI_EPS0 * q[self.p_i] * q[self.p_j]
        s_ij = np.sqrt(sigma[self.p_i] * sigma[self.p_j])
        e_ij = se[self.p_i] * se[self.p_j]
        self.p_a, self.p_c = d["fudge_lj"] * 4 * e_ij * s_ij ** 12, d["fudge_lj"] * 4 * e_ij * s_ij ** 6
        # bonded
        b = np.array(top.bonds, float).reshape(-1, 4)
        self.b_i, self.b_j, self.b_0, self.b_k = b[:, 0].astype(int), b[:, 1].astype(int), b[:, 2], b[:, 3]
        a = np.array(top.angles, float).reshape(-1, 5)
        self.a_idx, self.a_0, self.a_k = a[:, :3].astype(int).T, a[:, 3], a[:, 4]
        r = np.array(top.rb, float).reshape(-1, 10)
        self.rb_idx, self.rb_c = r[:, :4].astype(int).T, r[:, 4:]
        p = np.array(top.periodic, float).reshape(-1, 7)
        self.pd_idx, self.pd_phase, self.pd_k, self.pd_n = p[:, :4].astype(int).T, p[:, 4], p[:, 5], p[:, 6]
        # GBSA-OBC2
        rad = np.array([MBONDI2.get(int(zz), 0.15) for zz in top.z])
        for h in range(n):
            if top.z[h] == 1 and adj[h] and top.z[min(adj[h])] == 7:
                rad[h] = 0.13
        self.gb_radius = rad
        self.gb_or = rad - GB_OFFSET
        self.gb_sr = self.gb_or * np.array([SCREEN.get(int(zz), 0.8) for zz in top.z])
        self.gb_pre = -ONE_4PI_EPS0 * (1 / SOLUTE_DIELECTRIC - 1 / SOLVENT_DIELECTRIC)
        self.q = q

    # -------------------------------------------------------------- terms
    def _pairs(self, x, i, j, qq, A, C, g):
        d = x[i] - x[j]
        r2 = np.einsum("ij,ij->i", d, d)
        ir2 = 1 / r2
        ir = np.sqrt(ir2)
        ir6 = ir2 ** 3
        e_c, e_lj = qq * ir, (A * ir6 - C) * ir6
        dEdr_over_r = -(qq * ir + 12 * A * ir6 * ir6 - 6 * C * ir6) * ir2
        f = dEdr_over_r[:, None] * d
        np.add.at(g, i, f)
        np.add.at(g, j, -f)
        return e_c.sum(), e_lj.sum()

    def _bonded(self, x, g):
        e = {}
        d = x[self.b_i] - x[self.b_j]
        r = np.linalg.norm(d, axis=1)
        dr = r - self.b_0
        e["bond"] = 0.5 * (self.b_k * dr * dr).sum()
        f = (self.b_k * dr / r)[:, None] * d
        np.add.at(g, self.b_i, f)
        np.add.at(g, self.b_j, -f)
        i, j, k = self.a_idx
        th, di, dj, dk = _angles(x, i, j, k)
        dth = th - self.a_0
        e["angle"] = 0.5 * (self.a_k * dth * dth).sum()
        w = (self.a_k * dth)[:, None]
        for idx, dd in ((i, di), (j, dj), (k, dk)):
            np.add.at(g, idx, w * dd)
        e["rb"] = 0.0
        if len(self.rb_c):
            i, j, k, l = self.rb_idx
            phi, *dphi = _dihedrals(x, i, j, k, l)
            cpsi = -np.cos(phi)                       # psi = phi - 180
            powers = cpsi[:, None] ** np.arange(6)
            e["rb"] = (self.rb_c * powers).sum()
            dV = (self.rb_c[:, 1:] * np.arange(1, 6) * powers[:, :5]).sum(1) * np.sin(phi)
            for idx, dd in zip((i, j, k, l), dphi):
                np.add.at(g, idx, dV[:, None] * dd)
        e["periodic"] = 0.0
        if len(self.pd_k):
            i, j, k, l = self.pd_idx
            phi, *dphi = _dihedrals(x, i, j, k, l)
            arg = self.pd_n * phi - self.pd_phase
            e["periodic"] = (self.pd_k * (1 + np.cos(arg))).sum()
            dV = -self.pd_k * self.pd_n * np.sin(arg)
            for idx, dd in zip((i, j, k, l), dphi):
                np.add.at(g, idx, dV[:, None] * dd)
        return e

    def born_radii(self, x):
        """Born radii (nm), with dB_i/dr_ij pieces for the gradient."""
        n = self.top.n
        d = x[:, None, :] - x[None, :, :]
        r = np.sqrt((d * d).sum(-1))
        np.fill_diagonal(r, 1.0)
        orr, sr = self.gb_or[:, None], self.gb_sr[None, :]
        U = r + sr
        D = np.abs(r - sr)
        L = np.maximum(orr, D)
        on = (r + sr > orr) & ~np.eye(n, dtype=bool)
        I = 0.5 * (1 / L - 1 / U + 0.25 * (r - sr ** 2 / r) * (1 / U ** 2 - 1 / L ** 2) + 0.5 * np.log(L / U) / r)
        # d I / d r (L depends on r only where |r - sr| > or)
        dL = np.where(D > orr, np.sign(r - sr), 0.0)
        dI = 0.5 * (-dL / L ** 2 + 1 / U ** 2
                    + 0.25 * (1 + sr ** 2 / r ** 2) * (1 / U ** 2 - 1 / L ** 2)
                    + 0.25 * (r - sr ** 2 / r) * (-2 / U ** 3 + 2 * dL / L ** 3)
                    + 0.5 * (dL / L - 1 / U) / r - 0.5 * np.log(L / U) / r ** 2)
        I, dI = np.where(on, I, 0.0), np.where(on, dI, 0.0)
        a, b, c = OBC
        psi = I.sum(1) * self.gb_or
        t = np.tanh(a * psi - b * psi ** 2 + c * psi ** 3)
        B = 1 / (1 / self.gb_or - t / self.gb_radius)
        dB_dpsi = B ** 2 * (1 - t ** 2) * (a - 2 * b * psi + 3 * c * psi ** 2) / self.gb_radius
        # dB_i/dr_ij = dB_dpsi_i * or_i * dI_ij
        return B, (dB_dpsi * self.gb_or)[:, None] * dI, d, r

    def _gb(self, x, g):
        B, dBdr, d, r = self.born_radii(x)
        n = self.top.n
        q, pre = self.q, self.gb_pre
        np.fill_diagonal(r, 0.0)
        r2 = r * r
        BB = B[:, None] * B[None, :]
        ex = np.exp(-r2 / (4 * BB))
        f2 = r2 + BB * ex
        f = np.sqrt(f2)
        qq = q[:, None] * q[None, :]
        off = ~np.eye(n, dtype=bool)
        e_pair = 0.5 * pre * (qq / f)[off].sum()            # i < j counted once
        e_self = 0.5 * pre * (q * q / B).sum()
        rad = self.gb_radius
        sa = 4 * np.pi * ACE_ENERGY * (rad + ACE_PROBE) ** 2 * (rad / B) ** 6
        e_sa = sa.sum()
        # dE/dB_k
        Dm = r2 / (4 * BB)
        dEdB = np.where(off, -0.5 * pre * qq * f ** -3 * (B[None, :] * ex * (1 + Dm)), 0.0).sum(1) \
            - 0.5 * pre * q * q / B ** 2 - 6 * sa / B
        # explicit r dependence of the pair term at fixed B: d(1/f)/d vec r = -f^-3 (1 - ex/4) vec r
        w = np.where(off, -pre * qq * f ** -3 * (1 - ex / 4), 0.0)      # full sum over j (both orders) / 2 * 2
        g += (w[:, :, None] * d).sum(1)
        # through Born radii: dE = sum_i dEdB_i * sum_j dBdr_ij * d r_ij
        np.fill_diagonal(r, 1.0)
        u = d / r[:, :, None]
        c_ij = dEdB[:, None] * dBdr
        g += ((c_ij + c_ij.T)[:, :, None] * u).sum(1)
        return {"gb": e_pair + e_self, "ace": e_sa}

    def bonded_gauss_newton(self, x_nm):
        """Positive semi-definite model Hessian (3N x 3N, kJ/mol/nm^2) from the bonded terms: sum_t k_t J_t^T J_t,
        J_t the gradient of the bond length, angle or dihedral. Used only to precondition minimization."""
        x = np.asarray(x_nm, float).reshape(self.top.n, 3)
        n3 = 3 * self.top.n
        rows, ks = [], []

        def add(idx_list, grads, k):
            J = np.zeros((len(k), n3))
            for idx, gr in zip(idx_list, grads):
                for c in range(3):
                    np.add.at(J, (np.arange(len(k)), 3 * idx + c), gr[:, c])
            rows.append(J)
            ks.append(k)
        d = x[self.b_i] - x[self.b_j]
        u = d / np.linalg.norm(d, axis=1)[:, None]
        add((self.b_i, self.b_j), (u, -u), self.b_k)
        i, j, k = self.a_idx
        _, di, dj, dk = _angles(x, i, j, k)
        add((i, j, k), (di, dj, dk), self.a_k)
        if len(self.rb_c):
            i, j, k, l = self.rb_idx
            _, *dphi = _dihedrals(x, i, j, k, l)
            add((i, j, k, l), dphi, (np.abs(self.rb_c) * np.arange(6) ** 2).sum(1))
        if len(self.pd_k):
            i, j, k, l = self.pd_idx
            _, *dphi = _dihedrals(x, i, j, k, l)
            add((i, j, k, l), dphi, np.abs(self.pd_k) * self.pd_n ** 2)
        J, kk = np.vstack(rows), np.concatenate(ks)
        return J.T @ (kk[:, None] * J)

    def energy_terms(self, x_nm):
        """(dict of term energies in kJ/mol, gradient (n, 3) in kJ/mol/nm), topology atom order."""
        x = np.asarray(x_nm, float).reshape(self.top.n, 3)
        g = np.zeros_like(x)
        e = self._bonded(x, g)
        e["coulomb"], e["lj"] = self._pairs(x, self.nb_i, self.nb_j, self.nb_qq, self.nb_a, self.nb_c, g)
        e["coulomb14"], e["lj14"] = self._pairs(x, self.p_i, self.p_j, self.p_qq, self.p_a, self.p_c, g)
        if self.solvent == "obc2":
            e.update(self._gb(x, g))
        return e, g


class NumpyEnergy:
    """Energy model for d5_ensemble: Angstrom and kcal/mol, coordinates in RDKit atom order."""

    backend = "numpy"

    def __init__(self, model, perm, name, gtol=1e-6):
        self.model = model
        self.n = model.top.n
        self.perm = np.asarray(perm, int)       # perm[i] = RDKit index of topology atom i
        self.name = name
        self.gtol = gtol

    def energy_grad(self, x_a):
        x = np.asarray(x_a, float).reshape(self.n, 3)[self.perm] / 10.0
        e, g_t = self.model.energy_terms(x)
        g = np.empty_like(g_t)
        g[self.perm] = g_t * KJ_TO_KCAL / 10.0
        return sum(e.values()) * KJ_TO_KCAL, g.ravel()

    def energy(self, x_a):
        return self.energy_grad(x_a)[0]

    def grad(self, x_a):
        """dE/dx in kcal/mol/A, flattened [3N], RDKit order."""
        return self.energy_grad(x_a)[1]

    def preconditioner(self, x_a, eps_frac=PRECOND_EPS):
        """Upper Cholesky factor U of the bonded Gauss-Newton Hessian (RDKit order, kcal/mol/A^2), with eps_frac
        times its mean diagonal added so rigid-body and soft modes stay finite: H + eps I = U^T U."""
        from scipy.linalg import cholesky
        xt = np.asarray(x_a, float).reshape(self.n, 3)[self.perm] / 10.0
        Ht = self.model.bonded_gauss_newton(xt) * KJ_TO_KCAL / 100.0
        p3 = (3 * self.perm[:, None] + np.arange(3)).ravel()
        H = np.empty_like(Ht)
        H[np.ix_(p3, p3)] = Ht
        return cholesky(H + eps_frac * np.trace(H) / len(H) * np.eye(len(H)), lower=False)

    def minimize(self, x_a):
        """WARMUP_ITER plain L-BFGS-B iterations, then L-BFGS-B in coordinates whitened by the bonded Gauss-Newton
        Hessian (x = x1 + U^-1 y), which removes the bond/torsion stiffness gap; plain polish only if the gradient
        is still above POLISH_GMAX. The warm-up keeps large early preconditioned steps from leaving the start's
        basin (runs/minimizer_benchmark-20261107)."""
        from scipy.linalg import solve_triangular
        from scipy.optimize import minimize
        x0 = minimize(self.energy_grad, np.asarray(x_a, float).ravel(), jac=True, method="L-BFGS-B",
                      options=dict(maxiter=WARMUP_ITER, gtol=self.gtol, ftol=1e-14)).x
        U = self.preconditioner(x0)

        def fy(y):
            e, g = self.energy_grad(x0 + solve_triangular(U, y, lower=False))
            return e, solve_triangular(U, g, trans="T", lower=False)
        r = minimize(fy, np.zeros_like(x0), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=20000, gtol=self.gtol * 1e-2, ftol=1e-14, maxcor=30))
        x = x0 + solve_triangular(U, r.x, lower=False)
        if np.abs(self.grad(x)).max() > POLISH_GMAX:
            x = minimize(self.energy_grad, x, jac=True, method="L-BFGS-B",
                         options=dict(maxiter=20000, gtol=self.gtol, ftol=1e-14, maxcor=30)).x
        return x.reshape(self.n, 3)


def build_model(top_path, solvent="obc2", polar_h_lj=None):
    return OplsObc2(GmxTopology(top_path), polar_h_lj=polar_h_lj, solvent=solvent)
