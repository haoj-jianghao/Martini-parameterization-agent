"""Fit GROMACS CG bonded parameters to a (weighted) distribution by Boltzmann inversion. numpy only.

Forms (standard Martini 3 / Bartender choices):
  bond      funct 1   V = 1/2 kb (r - b0)^2               b0 = <r>, kb = kT / var(r)
            -> written as a constraint if kb > KB_CONSTRAINT (2e4 kJ mol-1 nm-2, Bartender's rule)
  angle     funct 2   V = 1/2 ktheta (cos t - cos t0)^2    (G96 cosine-harmonic)
            P(t)/sin t ∝ exp(-V/kT) is Gaussian in c = cos t, truncated to c in [-1, 1]:
            maximum-likelihood fit of that truncated Gaussian (moments are biased when the
            well reaches t = 180 deg: soft angles near linear, common in Martini)
  dihedral  funct 2   V = 1/2 k (phi - phi0)^2             if the distribution is single-peaked and
                                                          narrow (circular std < NARROW_DEG)
            funct 9   V = sum_n k_n [1 + cos(n phi - phi_n)], n = 1..4, by least squares on
                      U(phi) = -kT ln P(phi) (72-bin histogram). Empty bins are held at a ceiling
                      U_CEILING_KT above the minimum with a small weight FLOOR_WEIGHT: fitting only
                      the populated bins leaves the series free elsewhere, and narrow distributions
                      then gave cancelling amplitudes of order 1e5 kJ/mol.
The dihedral form can be forced (form=...), so two distributions can be fitted with the same form.
Units: nm, degrees in and out (internally radians), kJ/mol.
"""
import math

import numpy as np

KB = 0.0083144626          # kJ mol-1 K-1
KB_CONSTRAINT = 2.0e4      # kJ mol-1 nm-2
NARROW_DEG = 20.0
N_FOURIER = 4
U_CEILING_KT = 6.0         # unpopulated dihedral region held this far above the minimum
FLOOR_WEIGHT = 0.01        # relative to the highest histogram bin
FOURIER_MIN_K = 0.05       # kJ/mol: smaller Fourier components are not written


def _norm(w, n):
    w = np.ones(n) if w is None else np.asarray(w, float)
    return w / w.sum()


def fit_bond(r, w=None, T=300.0):
    r = np.asarray(r, float); w = _norm(w, len(r))
    m = float(np.sum(w * r)); var = float(np.sum(w * (r - m) ** 2))
    kb = KB * T / var
    return dict(kind="bond", funct=1, b0=m, kb=kb, constraint=bool(kb > KB_CONSTRAINT))


def fit_angle(theta_deg, w=None, T=300.0):
    from scipy.optimize import minimize
    from scipy.special import log_ndtr, ndtr
    c = np.cos(np.radians(np.asarray(theta_deg, float))); w = _norm(w, len(c))
    m = float(np.sum(w * c)); sd = math.sqrt(float(np.sum(w * (c - m) ** 2)))

    def nll(p):
        mu, ls = p
        s = math.exp(ls)
        a, b = (-1 - mu) / s, (1 - mu) / s
        # log(Phi(b) - Phi(a)), stable: use the upper tail when both bounds are above the mean
        if a > 0:
            Z = ndtr(-a) - ndtr(-b)
        else:
            Z = ndtr(b) - ndtr(a)
        logZ = math.log(Z) if Z > 1e-300 else log_ndtr(b if b < 0 else -a)
        return float(np.sum(w * 0.5 * ((c - mu) / s) ** 2) + ls + logZ)

    res = minimize(nll, [m, math.log(max(sd, 1e-4))], method="Nelder-Mead",
                   options=dict(xatol=1e-7, fatol=1e-10, maxiter=4000))
    mu, s = res.x[0], math.exp(res.x[1])
    return dict(kind="angle", funct=2, theta0=math.degrees(math.acos(max(-1.0, min(1.0, mu)))),
                ktheta=KB * T / s ** 2, cos_theta0=float(mu))


def _circ_stats(phi, w):
    C, S = np.sum(w * np.cos(phi)), np.sum(w * np.sin(phi))
    R = min(math.hypot(C, S), 1.0)
    return math.atan2(S, C), math.sqrt(-2 * math.log(max(R, 1e-12)))


def _hist(phi, w, bins=72):
    edges = np.linspace(-math.pi, math.pi, bins + 1)
    h, _ = np.histogram(phi, edges, weights=w)
    h = h / h.sum() / (edges[1] - edges[0])
    # light circular smoothing (3-bin kernel) to tame empty bins
    hs = 0.25 * np.roll(h, 1) + 0.5 * h + 0.25 * np.roll(h, -1)
    return 0.5 * (edges[1:] + edges[:-1]), hs


def n_peaks_circular(phi_deg, w=None, bins=72, min_frac=0.05):
    phi = np.radians(np.asarray(phi_deg, float)); w = _norm(w, len(phi))
    _, h = _hist(phi, w, bins)
    peak = (h > np.roll(h, 1)) & (h >= np.roll(h, -1)) & (h > min_frac * h.max())
    return int(peak.sum())


def fit_dihedral(phi_deg, w=None, T=300.0, form=None):
    phi = np.radians(np.asarray(phi_deg, float)); w = _norm(w, len(phi))
    kT = KB * T
    mu, sd = _circ_stats(phi, w)
    if form is None:
        form = "harmonic" if (n_peaks_circular(phi_deg, w) == 1 and math.degrees(sd) < NARROW_DEG) else "fourier"
    if form == "harmonic":
        return dict(kind="dihedral", funct=2, phi0=math.degrees(mu), k=kT / sd ** 2, form="harmonic")
    x, P = _hist(phi, w)
    U = -kT * np.log(np.maximum(P, 1e-300))
    U = np.minimum(U, U.min() + U_CEILING_KT * kT)
    A = [np.ones_like(x)]
    for n in range(1, N_FOURIER + 1):
        A += [np.cos(n * x), np.sin(n * x)]
    A = np.array(A).T
    sw = np.sqrt(P / P.max() + FLOOR_WEIGHT)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], U * sw, rcond=None)
    terms = []
    for n in range(1, N_FOURIER + 1):
        a, b = coef[2 * n - 1], coef[2 * n]
        k = math.hypot(a, b)                      # a cos + b sin = k cos(n phi - phi_n)
        terms.append(dict(n=n, k=k, phi_n=math.degrees(math.atan2(b, a))))
    return dict(kind="dihedral", funct=9, form="fourier", terms=terms)


def fit_term(kind, values, w=None, T=300.0, form=None):
    if kind == "bond":
        return fit_bond(values, w, T)
    if kind == "angle":
        return fit_angle(values, w, T)
    return fit_dihedral(values, w, T, form)


# ------------------------------------------------------ fitted potentials --
def potential(p, x_deg_or_nm):
    """V(x) in kJ/mol for fitted parameters p (x in nm for bonds, degrees otherwise)."""
    x = np.asarray(x_deg_or_nm, float)
    if p["kind"] == "bond":
        return 0.5 * p["kb"] * (x - p["b0"]) ** 2
    if p["kind"] == "angle":
        return 0.5 * p["ktheta"] * (np.cos(np.radians(x)) - math.cos(math.radians(p["theta0"]))) ** 2
    if p["form"] == "harmonic":
        d = (np.radians(x) - math.radians(p["phi0"]) + math.pi) % (2 * math.pi) - math.pi
        return 0.5 * p["k"] * d ** 2
    xr = np.radians(x)
    return sum(t["k"] * (1 + np.cos(t["n"] * xr - math.radians(t["phi_n"]))) for t in p["terms"])


def implied_distribution(p, T=300.0, n=3601, bond_range=None):
    """(grid, density) of the coordinate under the fitted CG potential alone (Jacobian included:
    sin(theta) for angles; for bonds the r^2 factor is omitted, as in the harmonic fit)."""
    kT = KB * T
    if p["kind"] == "bond":
        s = math.sqrt(kT / p["kb"])
        lo, hi = bond_range or (p["b0"] - 8 * s, p["b0"] + 8 * s)
        x = np.linspace(max(lo, 1e-4), hi, n); J = np.ones_like(x)
    elif p["kind"] == "angle":
        x = np.linspace(0, 180, n); J = np.sin(np.radians(x))
    else:
        x = np.linspace(-180, 180, n); J = np.ones_like(x)
    V = potential(p, x)
    rho = J * np.exp(-(V - V.min()) / kT)
    rho /= np.trapezoid(rho, x)
    return x, rho
