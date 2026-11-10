#!/usr/bin/env python3
"""Run EM, eq1 and eq2 of built models with several velocity seeds. One passing run can be chance.

usage: multi_seed_check.py OUT_DIR MOLDIR [MOLDIR ...] [--seeds 11 22 33 44] [--jobs 4] [--threads 2]
The environment file is $MARTINI_ENV, or environment.json next to the skill folder.
Writes OUT_DIR/results.json. A run fails on a non-zero mdrun exit. LINCS warnings are counted separately.
"""
import argparse, json, os, re, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ENV_PATH = Path(os.environ.get("MARTINI_ENV", Path(__file__).resolve().parents[2] / "environment.json"))
GMX = json.loads(ENV_PATH.read_text())["gromacs_command"]


def run(job):
    md, seed, out, threads = job
    src, d = md / "test", out / md.parent.name / md.name / f"seed{seed}"
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    for f in ("em.mdp", "eq2.mdp", "system.top", "molecule.itp"):
        shutil.copy(src / f, d / f)
    (d / "eq1.mdp").write_text(re.sub(r"(?mi)^gen[_-]seed.*$", f"gen_seed = {seed}", (src / "eq1.mdp").read_text()))
    start = "solvated_ions.gro" if (src / "solvated_ions.gro").exists() else "solvated.gro"
    shutil.copy(src / start, d / "start.gro")
    lincs = 0
    for st, mdp, c, cpt in (("em", "em.mdp", "start.gro", None), ("eq1", "eq1.mdp", "em.gro", None), ("eq2", "eq2.mdp", "eq1.gro", "eq1.cpt")):
        args = [*GMX, "grompp", "-f", mdp, "-c", c, "-p", "system.top", "-o", f"{st}.tpr", "-maxwarn", "0"] + (["-t", cpt] if cpt else [])
        if subprocess.run(args, cwd=d, capture_output=True).returncode:
            return f"{md.parent.name}/{md.name}", seed, f"grompp_{st}", lincs
        r = subprocess.run([*GMX, "mdrun", "-deffnm", st, "-ntmpi", "1", "-ntomp", str(threads), "-nb", "cpu", "-pin", "off"], cwd=d, capture_output=True, text=True)
        log = d / f"{st}.log"
        lincs += log.read_text(errors="replace").count("LINCS WARNING") if log.exists() else 0
        if r.returncode:
            return f"{md.parent.name}/{md.name}", seed, f"fail_{st}", lincs
    return f"{md.parent.name}/{md.name}", seed, "ok", lincs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("moldirs", nargs="+", type=Path)
    ap.add_argument("--seeds", nargs="+", type=int, default=[11, 22, 33, 44])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    jobs = [(m.resolve(), s, a.out.resolve(), a.threads) for m in a.moldirs for s in a.seeds]
    res = {}
    with ThreadPoolExecutor(a.jobs) as ex:
        for name, seed, result, lincs in ex.map(run, jobs):
            res.setdefault(name, []).append(dict(seed=seed, result=result, lincs_warnings=lincs))
            print(name, seed, result, lincs, flush=True)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "results.json").write_text(json.dumps(res, indent=1))
    print("molecules passing every seed:", sum(all(x["result"] == "ok" for x in r) for r in res.values()), "of", len(res))
    for n, r in res.items():
        k = sum(x["result"] == "ok" for x in r)
        if k < len(r):
            print(" ", n, f"{k}/{len(r)}", [x["result"] for x in r])


if __name__ == "__main__":
    main()
