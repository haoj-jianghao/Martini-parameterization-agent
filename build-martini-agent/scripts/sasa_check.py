#!/usr/bin/env python3
"""SASA of a built model against its D5 atomistic ensemble. Report-only: it changes no model or status.

A: atoms of the D5 ensemble, Rowland-Taylor radii from the tutorial's vdwradii_AA.dat.
B: the same frames mapped to beads, with R/S/T radii 0.264/0.230/0.191 nm (Souza SI C3 p.42).
Bands on B against A: pass within 5%, review to 10%, flag above 10% (tutorial section 7).
Runs under python_command. No dynamics.

Frames come from bonded/ensemble_atoms.npz, or the ensemble is regenerated into OUT from the same AA topology
(the supplied one, or opls-ff-gen) and the NumPy energy model (only for records made with them). Regenerated frames must match bonded/ensemble_beads.npz (see frames_vs_saved), or the
molecule fails, because that would mean comparing a different ensemble from the one the model came from.

usage: sasa_check.py OUT MOLDIR [MOLDIR ...] [--frames 200] [--jobs 4]
"""
import argparse
import json
import os
import re
import subprocess
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from derive_bonded import mapping_matrix  # noqa: E402
from d5bonded.ensemble import d5_ensemble  # noqa: E402
from derive_bonded import get_topology  # noqa: E402
from rdkit import Chem  # noqa: E402

ENV_PATH = Path(os.environ.get("MARTINI_ENV", Path(__file__).resolve().parents[2] / "environment.json")).resolve()
ENV = json.loads(ENV_PATH.read_text())
GMX = ENV["gromacs_command"]
PROBE, NDOTS = 0.191, 4800
BEAD_RADII = {"R": 0.264, "S": 0.230, "T": 0.191}
AA_RADII_FILE = ENV_PATH.parent / "references/original/tutorial/parametrize_new_small_molecule/ENAP-in-water/7_SASA/vdwradii_AA.dat"
CROSS_CHECK_FRAMES = 3
CROSS_CHECK_TOL = 0.02
MAX_DISTANCE_DRIFT_NM = 2e-3    # measured: 0.0012 nm or less over reruns of the most sensitive molecule, RIBF


def aa_radii():
    out = {}
    for line in AA_RADII_FILE.read_text().splitlines():
        f = line.split(";")[0].split()
        if len(f) == 3:
            out[f[1]] = float(f[2])
    return out


def bead_radius(bead_type):
    return BEAD_RADII[bead_type[0] if bead_type[0] in "ST" else "R"]


def band(pct):
    return "pass" if abs(pct) <= 5 else "review" if abs(pct) <= 10 else "flag"


def sphere_points(n):
    i = np.arange(n) + 0.5
    phi, th = np.arccos(1 - 2 * i / n), np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def shrake_rupley(pos, radii):
    """Independent of gmx sasa: same radii and probe, same number of dots."""
    pts = sphere_points(NDOTS)
    r = np.asarray(radii) + PROBE
    total = 0.0
    for i in range(len(pos)):
        p = pos[i] + r[i] * pts
        free = np.ones(len(p), bool)
        for j in range(len(pos)):
            if j != i:
                free &= np.sum((p - pos[j]) ** 2, axis=1) >= r[j] ** 2
        total += free.mean() * 4 * np.pi * r[i] ** 2
    return float(total)


def gmx_sasa(workdir, frames, radii):
    """Per-frame SASA (nm2) from gmx sasa on plain .gro files with a local vdwradii.dat."""
    wd = Path(workdir)
    wd.mkdir(parents=True, exist_ok=True)
    n = len(radii)
    (wd / "vdwradii.dat").write_text("".join(f"MOL  A{i + 1}  {r:.4f}\n" for i, r in enumerate(radii)))

    def block(x, title):
        c = x - x.mean(0) + 15.0
        body = "".join(f"{1:5d}{'MOL':<5s}{'A' + str(i + 1):>5s}{(i + 1) % 100000:5d}{a:8.3f}{b:8.3f}{d:8.3f}\n"
                       for i, (a, b, d) in enumerate(c))
        return f"{title}\n{n:5d}\n{body}  30.0 30.0 30.0\n"
    (wd / "ref.gro").write_text(block(frames[0], "ref"))
    (wd / "frames.gro").write_text("".join(block(x, f"frame t= {k}.0") for k, x in enumerate(frames)))
    p = subprocess.run([*GMX, "sasa", "-s", "ref.gro", "-f", "frames.gro", "-o", "sasa.xvg", "-probe", str(PROBE),
                        "-ndots", str(NDOTS), "-surface", "all", "-nopbc"], cwd=wd, capture_output=True, text=True)
    (wd / "gmx_sasa.log").write_text(p.stdout + p.stderr)
    if p.returncode:
        raise RuntimeError(f"gmx sasa failed in {wd}: {p.stderr[-300:]}")
    v = np.array([float(l.split()[1]) for l in open(wd / "sasa.xvg") if l[0] not in "#@"])
    if len(v) != len(frames):
        raise RuntimeError(f"gmx sasa returned {len(v)} values for {len(frames)} frames in {wd}")
    return v


def checked_sasa(workdir, frames, radii):
    """gmx sasa, cross-checked on a few frames by the independent calculation."""
    v = gmx_sasa(workdir, frames, radii)
    checks = []
    for k in np.linspace(0, len(frames) - 1, min(CROSS_CHECK_FRAMES, len(frames))).astype(int):
        sr = shrake_rupley(np.round(frames[k] - frames[k].mean(0) + 15.0, 3), radii)
        rel = abs(v[k] - sr) / sr
        checks.append(dict(frame=int(k), gmx_nm2=round(float(v[k]), 4), independent_nm2=round(sr, 4), relative=round(float(rel), 5)))
        if rel > CROSS_CHECK_TOL:
            raise RuntimeError(f"cross-check failed in {workdir}: {checks[-1]}")
    return v, checks


def self_test(workdir):
    """Isolated beads must give 4 pi (r + probe)^2: gmx sasa falls back to guessed radii silently."""
    out = {}
    for size, r in BEAD_RADII.items():
        got = float(gmx_sasa(Path(workdir) / size, np.zeros((1, 1, 3)), [r])[0])
        want = 4 * np.pi * (r + PROBE) ** 2
        out[size] = dict(gmx_nm2=round(got, 4), exact_nm2=round(want, 4))
        if abs(got - want) / want > 0.005:
            raise RuntimeError(f"SASA self-test failed for {size}: {got} vs {want}")
    return out


def regenerate_atoms(md, out, req, bonded):
    """The 4000 weighted frames behind bonded/ensemble_beads.npz, by the same seed and settings."""
    em_rec = bonded["energy_model"]
    if ((em_rec.get("topology_source", "engine") == "engine" and em_rec["engine"] != "opls-ff-gen")
            or "openmm" in json.dumps(bonded.get("software", {})).lower()):
        raise RuntimeError("bonded.json predates the opls-ff-gen + NumPy energy model and has no ensemble_atoms.npz; "
                           "re-run stage 4 (cgbuild.py bonded) before the SASA check")
    t = get_topology(req["original_smiles"], req, md if req.get("aa_topology") else out / "topology", ENV)
    from derive_bonded import energy_model
    em = energy_model(t)
    b = ENV["bonded"]
    ens = d5_ensemble(req["original_smiles"], n_conf=b["n_conformers"], temperature=bonded["temperature_K"],
                      n_samples=b["n_samples_per_minimum"], seed=bonded["seed"], energy_model=em,
                      hess_e_window=b["hess_e_window_kcal_mol"])
    w = ens.weights / ens.weights.sum()
    keep = np.random.default_rng(bonded["seed"]).choice(len(w), size=min(4000, len(w)), p=w)
    return ens.xyz_nm[keep]


def frames_vs_saved(mapped, saved):
    """The ensemble is not bit-reproducible between runs: small numerical noise changes minimum weights and some
    frames. SASA is an ensemble mean, so the gate is on mean bead distances; changed frames are only recorded."""
    if mapped.shape != saved.shape:
        return 1.0, float("inf")
    changed = float((np.abs(mapped - saved).max((1, 2)) > 1e-3).mean())
    dm = lambda r: np.linalg.norm(r[:, :, None] - r[:, None, :], axis=-1).mean(0)
    return changed, float(np.abs(dm(mapped) - dm(saved)).max())


def check_one(job):
    md, base, n_frames = job
    md = Path(md)
    label = f"{md.parent.name}/{md.name}"
    out = Path(base) / md.parent.name / md.name
    try:
        req = json.loads((md / "input/request.json").read_text())
        des = json.loads((md / "decisions/design.json").read_text())
        bonded = json.loads((md / "decisions/bonded.json").read_text())
        mol = Chem.AddHs(Chem.MolFromSmiles(req["original_smiles"]))
        M = mapping_matrix(des, mol)
        saved = md / "bonded/ensemble_atoms.npz"
        source = "saved"
        if saved.exists():
            atoms = np.load(saved)["atom_xyz_nm"].astype(float)
        else:
            atoms, source = regenerate_atoms(md, out, req, bonded), "regenerated"
        beads = np.load(md / "bonded/ensemble_beads.npz")["bead_xyz_nm"]
        mapped = np.einsum("ba,fai->fbi", M, atoms)
        changed, drift = frames_vs_saved(mapped, beads)
        if drift > MAX_DISTANCE_DRIFT_NM:
            raise RuntimeError(f"frames do not match ensemble_beads.npz: {changed:.1%} of frames differ, "
                               f"mean bead distances differ by up to {drift:.4f} nm")
        pick = np.random.default_rng(bonded["seed"]).choice(len(atoms), size=min(n_frames, len(atoms)), replace=False)
        table = aa_radii()
        ra = [table[a.GetSymbol()] for a in mol.GetAtoms()]
        rb = [bead_radius(b["type"]) for b in des["beads"]]
        a, ca = checked_sasa(out / "atoms", atoms[pick], ra)
        b, cb = checked_sasa(out / "beads", mapped[pick], rb)
        pct = 100 * (b.mean() - a.mean()) / a.mean()
        res = dict(id=label, molecule_id=req["molecule_id"], smiles=req["original_smiles"], n_frames=len(pick), frames=source,
                   n_beads=len(rb), heavy_atoms=req["heavy_atoms"], engine=bonded["energy_model"]["engine"],
                   aa_nm2=round(float(a.mean()), 4), aa_sd=round(float(a.std()), 4),
                   mapped_nm2=round(float(b.mean()), 4), mapped_sd=round(float(b.std()), 4),
                   difference_pct=round(float(pct), 2), band=band(pct),
                   difference_sd_pct=round(float(100 * np.std(b - a) / a.mean()), 2),
                   flags=sorted({f["flag"] for f in bonded["flags"]}), n_used_minima=bonded["ensemble"]["n_used"],
                   frames_changed_fraction=round(changed, 4), mean_distance_drift_nm=round(drift, 5), cross_check=dict(atoms=ca, beads=cb),
                   settings=dict(probe_nm=PROBE, ndots=NDOTS, bead_radii_nm=BEAD_RADII, atom_radii=str(AA_RADII_FILE.relative_to(ENV_PATH.parent)),
                                 source="Souza SI C3 p.42; tutorial section 7"))
        (out / "sasa.json").write_text(json.dumps(res, indent=1))
        print(label, res["aa_nm2"], res["mapped_nm2"], f"{pct:+.1f}%", res["band"], source, flush=True)
        return res
    except Exception as e:  # noqa: BLE001
        print(label, "ERROR", type(e).__name__, str(e)[:300], flush=True)
        return dict(id=label, error=f"{type(e).__name__}: {e}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("moldirs", nargs="+", type=Path)
    ap.add_argument("--frames", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    test = self_test(a.out / "self_test")
    print("self-test passed", test, flush=True)
    with Pool(a.jobs, maxtasksperchild=1) as pool:
        res = pool.map(check_one, [(str(m.resolve()), str(a.out.resolve()), a.frames) for m in a.moldirs], chunksize=1)
    (a.out / "summary.json").write_text(json.dumps(dict(self_test=test, molecules=res), indent=1))


if __name__ == "__main__":
    main()
