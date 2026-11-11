#!/usr/bin/env python3
"""Agent-free bonded terms for one molecule. Runs under the environment's python_command (RDKit, NumPy, SciPy).

Reads   input/request.json, decisions/design.json (beads, members, virtual sites); masses are computed (terms.bead_masses)
Writes  decisions/bonded.json        terms, exclusions, dropped terms, ensemble record
        bonded/ensemble_beads.npz    resampled bead coordinates (nm), for later analysis
        bonded/ensemble_atoms.npz    the same frames as atoms (nm), for the SASA check
        bonded/opls-ff-gen/          the engine's OPLS topology and its log (only when no topology was supplied)

The OPLS-AA topology is the one supplied at `cgbuild.py input --aa-topology` (input/aa_topology.top), used as given;
without one, the engine in opls_ff_gen_command is run (d5bonded/opls.py).
No molecular dynamics is run. Every value is a reproducible calculation from the SMILES, the topology and the mapping.

usage: derive_bonded.py MOLDIR --env ENVIRONMENT_JSON
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from d5bonded import fit as F  # noqa: E402
from d5bonded.ensemble import d5_ensemble  # noqa: E402
from d5bonded.geometry import angles, dihedrals, distances  # noqa: E402
from d5bonded.energy import NumpyEnergy, build_model  # noqa: E402
from d5bonded.opls import POLAR_H_LJ, TopologyError, build_topology, load_supplied_topology  # noqa: E402
from d5bonded.terms import bead_masses, build_terms, read_sigma, rigid_groups  # noqa: E402
from rdkit import Chem, rdBase  # noqa: E402

SOURCE = [{"path": "decisions/bonded.json", "locator": "ensemble record; D5 harmonic superposition + Boltzmann inversion"}]
FLEX_K = 1000000
FLEX_SOURCE = "references/original/small-molecules/models/itps/opt-poly/NAPH.itp, [constraints] k column"


def mapping_matrix(design, mol_h):
    """[n_beads, n_atoms], with H inheriting the allocation of its heavy atom, normalized per bead."""
    excl = set(design.get("hydrogen_excluded_parents", {}).get("atom_ids", []) or [])
    heavy = {}
    for bi, b in enumerate(design["beads"]):
        for aid, w in b["members"].items():
            heavy.setdefault(int(aid), {})[bi] = float(w)
    M = np.zeros((len(design["beads"]), mol_h.GetNumAtoms()))
    for a in mol_h.GetAtoms():
        i = a.GetIdx()
        if a.GetAtomicNum() == 1:
            parent = a.GetNeighbors()[0].GetIdx() + 1
            alloc = {} if parent in excl else heavy.get(parent, {})
        else:
            alloc = heavy.get(i + 1, {})
        for bi, w in alloc.items():
            M[bi, i] = w
    return M / M.sum(1, keepdims=True)


def energy_model(t):
    """The OPLS-AA + GBSA-OBC2 + polar-H LJ energy model of topology t (d5bonded/energy.py)."""
    name = (f"{t.engine} OPLS-AA ({t.charge_model}) + GBSA-OBC2"
            f" + polar-H LJ sigma {POLAR_H_LJ[0]} nm eps {POLAR_H_LJ[1]} kJ/mol")
    return NumpyEnergy(build_model(t.top, polar_h_lj=POLAR_H_LJ), t.perm, name)


def get_topology(smiles, req, md, env):
    """The supplied topology if the request has one, else the engine's (opls_ff_gen_command)."""
    if req.get("aa_topology"):
        return load_supplied_topology(smiles, req["aa_topology"], md)
    if not env.get("opls_ff_gen_command"):
        raise TopologyError("no AA topology supplied (cgbuild.py input --aa-topology) and no opls_ff_gen_command set")
    return build_topology(smiles, Path(md) / "bonded", env["opls_ff_gen_command"])


def derive_ensemble(smiles, req, md, env, T):
    """Topology (supplied or engine), energy model, ensemble; the attempt is recorded."""
    b = env["bonded"]
    try:
        t = get_topology(smiles, req, md, env)
        em = energy_model(t)
        ens = d5_ensemble(smiles, n_conf=b["n_conformers"], temperature=T, n_samples=b["n_samples_per_minimum"],
                          seed=env["random_seed"], energy_model=em, hess_e_window=b["hess_e_window_kcal_mol"])
    except (TopologyError, RuntimeError, OSError, KeyError, ValueError, TypeError) as e:
        src = ("supplied topology" if req.get("aa_topology") else "opls-ff-gen" if env.get("opls_ff_gen_command")
               else "stage 4")
        raise TopologyError(f"{src}: {type(e).__name__}: {str(e)[:300]}") from e
    return em, t, ens, [dict(engine=t.engine, source=t.source, result="used")]


def software_versions():
    import scipy
    return dict(python=sys.executable, rdkit=rdBase.rdkitVersion, numpy=np.__version__, scipy=scipy.__version__,
                energy_model="d5bonded/energy.py (NumPy)")


def stats(R, w, kind, beads):
    if kind == "dihedrals":
        phi = np.radians(dihedrals(R, [beads])[:, 0])
        C, S = np.sum(w * np.cos(phi)), np.sum(w * np.sin(phi))
        return dict(mean=round(float(np.degrees(np.arctan2(S, C))), 2),
                    std=round(float(np.degrees(np.sqrt(-2 * np.log(max(np.hypot(C, S), 1e-12))))), 2), unit="deg")
    v = distances(R, [beads])[:, 0] if kind in ("bonds", "constraints") else angles(R, [beads])[:, 0]
    m = float(np.sum(w * v))
    return dict(mean=round(m, 4), std=round(float(np.sqrt(np.sum(w * (v - m) ** 2))), 4),
                unit="nm" if kind in ("bonds", "constraints") else "deg")


def to_record(t, tid, R, w):
    sec, beads = t["section"], [int(b) for b in t["beads"]]
    if sec == "constraints":
        params = [dict(name="b0", value=round(t["b0"], 4), unit="nm")]
    elif sec == "bonds":
        params = [dict(name="b0", value=round(t["b0"], 4), unit="nm"),
                  dict(name="k", value=round(t["k"], 1), unit="kJ mol-1 nm-2")]
    elif sec == "angles":
        params = [dict(name="theta0", value=round(t["theta0"], 2), unit="deg"),
                  dict(name="k", value=round(t["k"], 2), unit="kJ mol-1")]
    elif t["funct"] == 2:
        params = [dict(name="phi0", value=round(t["phi0"], 2), unit="deg"),
                  dict(name="k", value=round(t["k"], 2), unit="kJ mol-1 rad-2")]
    else:
        params = []
        for x in t["fourier"]:
            params += [dict(name=f"phi_s_{x['n']}", value=round(x["phi_n"], 2), unit="deg", multiplicity=x["n"]),
                       dict(name=f"k_{x['n']}", value=round(x["k"], 3), unit="kJ mol-1", multiplicity=x["n"])]
    notes = [t[k] for k in ("capped", "floored", "rule") if t.get(k)]
    rec = dict(term_id=tid, section=sec, beads=[b + 1 for b in beads], function=t["funct"], role=t["role"],
               parameters=params, evidence_class="ensemble_derived", source=SOURCE,
               limitations="; ".join(notes) or None, ensemble=stats(R, w, sec, beads),
               fitted_before_limit={k: round(t[k], 2) for k in ("k_fit", "kb_fit") if t.get(k) is not None} or None)
    if sec == "constraints":
        rec["flexible_variant"] = dict(k=FLEX_K, source=FLEX_SOURCE)
    if sec == "dihedrals" and t["funct"] == 9:
        rec["potential_range_kJ_mol"] = round(t["range_kJ_mol"], 2)
    return rec


def virtual_site_parameters(design, R, w):
    """virtual_sites3 funct 1, r_i = r_j + a r_jk + b r_jl: (a, b) by weighted least squares on the ensemble.
    Each frame is fitted in its own coordinates, so frame orientations do not matter."""
    out = []
    for vs in design.get("virtual_sites", []):
        if vs.get("construction") != "virtual_sites3":
            continue
        i, j, k, l = vs["site"] - 1, *[c - 1 for c in vs["constructing"]]
        A = np.stack([R[:, k] - R[:, j], R[:, l] - R[:, j]], -1).reshape(-1, 2)
        y = (R[:, i] - R[:, j]).reshape(-1)
        sw = np.sqrt(np.repeat(w, 3))
        (a, b), *_ = np.linalg.lstsq(A * sw[:, None], y * sw, rcond=None)
        rms = float(np.sqrt(np.sum(np.repeat(w, 3) * (A @ np.array([a, b]) - y) ** 2) * 3))
        out.append(dict(site=vs["site"], function=1, a=round(float(a), 4), b=round(float(b), 4), rms_nm=round(rms, 4)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("moldir", type=Path)
    ap.add_argument("--env", type=Path, required=True)
    a = ap.parse_args()
    env = json.loads(a.env.read_text())
    md = a.moldir
    req = json.loads((md / "input/request.json").read_text())
    des = json.loads((md / "decisions/design.json").read_text())
    smiles = req["original_smiles"]
    if len(Chem.GetMolFrags(Chem.MolFromSmiles(smiles))) != 1:
        raise SystemExit("derive_bonded: one connected molecule expected")
    mol_h = Chem.AddHs(Chem.MolFromSmiles(smiles))
    b = env["bonded"]
    T = float(env["default_temperature_K"])
    t0 = time.time()
    em, top, ens, tried = derive_ensemble(smiles, req, md, env, T)
    M = mapping_matrix(des, mol_h)
    virtual = [i for i, bd in enumerate(des["beads"]) if any(v["site"] == i + 1 for v in des.get("virtual_sites", []))]
    is_virtual = np.isin(np.arange(M.shape[0]), virtual)
    masses, mass_record = bead_masses([bd["type"] for bd in des["beads"]], des.get("virtual_sites", []),
                                      rigid_groups(mol_h, M, is_virtual))
    sigma = read_sigma(a.env.parent / env["force_field"]["core_itp"], [bd["type"] for bd in des["beads"]])
    types = [bd["type"] for bd in des["beads"]]
    built = build_terms(mol_h, M, ens.xyz_nm, ens.weights, masses, is_virtual, T, types, sigma=sigma)

    R = np.einsum("ba,fai->fbi", M, ens.xyz_nm)
    w = ens.weights / ens.weights.sum()
    site_params = virtual_site_parameters(des, R, w)
    counters, records = {}, []
    prefix = {"bonds": "b", "constraints": "c", "angles": "a", "dihedrals": "d"}
    order = {"constraints": 0, "bonds": 1, "angles": 2, "dihedrals": 3}
    for t in sorted(built["terms"], key=lambda x: (order[x["section"]], x["beads"])):
        counters[t["section"]] = counters.get(t["section"], 0) + 1
        records.append(to_record(t, f"{prefix[t['section']]}{counters[t['section']]}", R, w))

    rng = np.random.default_rng(env["random_seed"])
    keep = rng.choice(len(w), size=min(4000, len(w)), p=w)
    (md / "bonded").mkdir(exist_ok=True)
    np.savez_compressed(md / "bonded/ensemble_beads.npz", bead_xyz_nm=R[keep].astype(np.float32))
    np.savez_compressed(md / "bonded/ensemble_atoms.npz", atom_xyz_nm=ens.xyz_nm[keep].astype(np.float32))
    q_max = float(np.abs(top.charges).max())
    used = [m for m in ens.minima if m["used"]]
    flags = list(built["flags"])
    if q_max > 2.0:
        flags.append(dict(flag="large_atomic_charge", max_abs_charge_e=round(q_max, 3),
                          note="hypervalent P or S: conformer weights are less certain"))
    odd = [n for n in top.notes or [] if not n.startswith("no [ defaults ]")]
    if odd:
        flags.append(dict(flag="supplied_topology_notes", notes=odd,
                          note="the supplied AA topology was used as given; see energy_model.topology_notes"))
    if ens.info["n_collapsed"]:
        flags.append(dict(flag="collapsed_minima_discarded", n_collapsed=ens.info["n_collapsed"], n_used_minima=len(used),
                          note="minima where a sigma-0 hydrogen fell onto a charged atom were removed; conformer weights are less certain"))
    if len(used) >= 0.9 * b["n_conformers"]:
        flags.append(dict(flag="conformer_search_not_saturated", n_used_minima=len(used),
                          note="nearly every conformer is a distinct minimum"))
    out = dict(
        schema_version=1, molecule_id=req["molecule_id"], method="D5 harmonic superposition + Boltzmann inversion; no MD",
        temperature_K=T,
        ensemble={k: ens.info[k] for k in ("n_conf", "n_unique", "n_used", "n_e_window_skipped", "n_saddle_rejected",
                                           "n_collapsed", "n_eff_minima", "top_weight", "rejected_frac", "achiral", "mirrored",
                                           "n_samples_per_min")},
        energy_model=dict(name=em.name, engine=top.engine, charge_model=top.charge_model, topology_source=top.source,
                          topology_notes=top.notes or [], engines_tried=tried,
                          max_abs_atomic_charge_e=round(q_max, 3), topology=str(top.itp.relative_to(md)),
                          conformer_search="ETKDGv3 + MMFF94, re-minimized under the energy model"),
        software=software_versions(),
        seed=env["random_seed"], wall_s=round(time.time() - t0, 1),
        rigid_groups=[[int(x) + 1 for x in g] for g in built["groups"]],
        terms=records,
        exclusions=[[int(i) + 1, int(j) + 1] for i, j in built["exclusions"]],
        virtual_site_parameters=site_params,
        masses=dict(bead_amu=[float(x) for x in masses], virtual_site_redistribution=mass_record,
                    rule="virtual-site mass split equally over the real beads of its rigid group (d5bonded.terms.bead_masses)"),
        close_pair_exclusions=[dict(c, beads=[int(x) + 1 for x in c["beads"]]) for c in built["close_pairs"]],
        nrexcl=1,
        dropped=[dict(d, beads=[int(x) + 1 for x in d["beads"]]) for d in built["dropped"]],
        flags=flags,
    )
    (md / "decisions").mkdir(exist_ok=True)
    (md / "decisions/bonded.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"bonded {req['molecule_id']}: {len(records)} terms, {len(out['exclusions'])} exclusions, "
          f"{len(out['dropped'])} dropped, engine {top.engine}, {len(used)} minima, "
          f"neff {ens.info['n_eff_minima']:.1f}, flags {[f['flag'] for f in flags]}")


if __name__ == "__main__":
    main()
