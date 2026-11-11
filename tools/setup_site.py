#!/usr/bin/env python3
"""Set up a site for the Martini parameterization agent (no LLM involved).

Run from a clone of this repository, with the references bundle downloaded:

    python3 tools/setup_site.py --bundle /path/to/martini-agent-references-v7.tar.gz \
        [--root ~/martini-agent-lab] [--gmx /path/to/gmx]

Each step is skipped when its result is already in place, so the script can be re-run
after you supply what it reports as missing. Steps:

  1. Python environment: create the conda env from environment.yml (or use --python),
     check RDKit, NumPy, SciPy, openpyxl.
  2. Skill: copy build-martini-agent/ into the skills directory and verify it against
     locks/skill.lock.json (an older copy is kept as build-martini-agent.bak-<time>).
  3. References: check the bundle against SHA256SUMS and unpack it under --root.
  4. GROMACS: find it and read its version.
  5. environment.json: create it from the template and fill python_command,
     gromacs_command, the GROMACS version and the force-field paths.
  6. Site downloads: the GROMACS manual pages and PDF for that version and the polyply SI
     (skip with --no-fetch); derive their text with pdftotext.
  7. Papers: check the publisher downloads you placed (Souza 2021, Alessandri 2022) and
     derive their text.
  8. MDP check: grompp -maxwarn 0 and a 10 ps run of the adapted production MDP on the
     bundled water box.
  9. Index and manifest: build_index.py check, reference_manifest.py verify.
 10. Lock (only with --lock, and only when nothing above failed, the papers are present
     and references/index/readiness.json exists).

It never edits bundled reference files and never runs downloaded code. The scientific
review (paper identity, sources.json, readiness.json) is prepare mode's job:
references/environment.md, section 0.

No OPLS engine is installed: stage 4 uses an OPLS-AA topology that the user makes
elsewhere for each molecule and passes to `cgbuild.py input --aa-topology`.
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SKILL_SRC = REPO / "build-martini-agent"
MANIFEST = REPO / "reference_manifest.json"
SKILL_LOCK = REPO / "locks" / "skill.lock.json"
TESTED_RDKIT = "2025.09.6"
TESTED_GMX = "2026.3"
FF_FILES = {"core_itp": "martini_v3.0.0.itp", "solvent_itp": "martini_v3.0.0_solvents_v1.itp",
            "ions_itp": "martini_v3.0.0_ions_v1.itp"}
GMX_PAGES = {"topology-file-formats.html": "topologies/topology-file-formats.html",
             "bonded-interactions.html": "functions/bonded-interactions.html",
             "interaction-methods.html": "functions/interaction-methods.html"}
PAPERS = {"Souza2021": ("article.pdf", "supporting-information.pdf", "souza2021"),
          "alessandri2022": ("article.pdf", "supporting-information.pdf", "alessandri2022")}

RESULTS = []          # (step, status, message); status: ok, skip, todo, warn, fail


def note(step, status, msg):
    RESULTS.append((step, status, msg))
    print(f"[{status:4s}] {step}: {msg}", flush=True)


def run(cmd, **kw):
    env = dict(os.environ)
    env.pop("OMP_NUM_THREADS", None)          # GROMACS rejects a mismatch with -ntomp
    env.update(kw.pop("env", {}))
    return subprocess.run([str(c) for c in cmd], capture_output=True, text=True, env=env, **kw)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check_sums(sums_file, folder):
    """True if every file listed in a sha256sum file matches; missing files fail."""
    bad = []
    for line in Path(sums_file).read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        p = Path(folder) / name.strip().lstrip("*")
        if not p.is_file() or sha256(p) != digest:
            bad.append(p.name)
    return bad


def tool(py, name):
    """An executable next to the interpreter (conda env or venv bin), else on PATH."""
    p = Path(py).parent / name
    return p if p.exists() else (Path(shutil.which(name)) if shutil.which(name) else None)


# 1 -------------------------------------------------------------------------------------
def python_env(a):
    if a.python:
        py = Path(a.python).expanduser().absolute()   # not resolve(): a venv python is a symlink
        note("python", "skip", f"using the given interpreter {py}")
    else:
        conda = shutil.which("mamba") or shutil.which("conda")
        if not conda:
            note("python", "fail", "conda (or mamba) not found; install it, or pass --python")
            return None
        envs = json.loads(run([conda, "env", "list", "--json"]).stdout or "{}").get("envs", [])
        prefix = next((Path(e) for e in envs if Path(e).name == a.env_name), None)
        if prefix is None:
            print(f"creating conda env {a.env_name} from environment.yml (several minutes) ...", flush=True)
            r = subprocess.run([conda, "env", "create", "-n", a.env_name, "-f", str(REPO / "environment.yml")])
            if r.returncode:
                note("python", "fail", f"conda env create failed (exit {r.returncode})")
                return None
            envs = json.loads(run([conda, "env", "list", "--json"]).stdout).get("envs", [])
            prefix = next(Path(e) for e in envs if Path(e).name == a.env_name)
            note("python", "ok", f"created conda env {a.env_name} at {prefix}")
        else:
            note("python", "skip", f"conda env {a.env_name} exists at {prefix}")
        py = prefix / "bin" / "python"
    r = run([py, "-c", "import rdkit, numpy, scipy, openpyxl, sys; "
                       "print(rdkit.__version__, numpy.__version__, scipy.__version__, sys.version.split()[0])"])
    if r.returncode:
        note("python", "fail", f"{py} cannot import rdkit/numpy/scipy/openpyxl: {r.stderr.strip()[-200:]}")
        return None
    rd, npv, spv, pyv = r.stdout.split()
    note("python", "ok" if rd == TESTED_RDKIT else "warn",
         f"Python {pyv}, RDKit {rd}, NumPy {npv}, SciPy {spv}"
         + ("" if rd == TESTED_RDKIT else f" (the index was built with RDKit {TESTED_RDKIT}; the index check may differ)"))
    return py


# 2 -------------------------------------------------------------------------------------
def skill(a, py):
    dest = Path(a.skills_dir).expanduser() / SKILL_SRC.name
    lock_tool = SKILL_SRC / "scripts" / "reference_lock.py"
    if dest.exists():
        if run([py, lock_tool, "verify", dest, SKILL_LOCK]).returncode == 0:
            note("skill", "skip", f"{dest} matches locks/skill.lock.json")
            return dest
        old = dest.with_name(f"{dest.name}.bak-{dt.datetime.now():%Y%m%d%H%M%S}")
        dest.rename(old)
        note("skill", "warn", f"an older copy was moved to {old}")
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(SKILL_SRC, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    ok = run([py, lock_tool, "verify", dest, SKILL_LOCK]).returncode == 0
    note("skill", "ok" if ok else "fail", f"installed to {dest}" + ("" if ok else "; it does not match the skill lock"))
    return dest


# 3 -------------------------------------------------------------------------------------
def references(a, root):
    refs = root / "references"
    for d in ("config", "runs"):
        (root / d).mkdir(parents=True, exist_ok=True)
    if refs.exists():
        note("references", "skip", f"{refs} exists (checked against the manifest in step 10)")
        return refs
    if not a.bundle:
        note("references", "todo", "download martini-agent-references-*.tar.gz from the release and pass --bundle")
        return None
    bundle = Path(a.bundle).expanduser().resolve()
    sums = Path(a.sums).expanduser() if a.sums else bundle.parent / "SHA256SUMS"
    if not sums.is_file():
        note("references", "fail", f"no SHA256SUMS beside the bundle ({sums}); pass --sums")
        return None
    listed = {l.split(None, 1)[1].strip().lstrip("*"): l.split()[0] for l in sums.read_text().splitlines() if l.strip()}
    if listed.get(bundle.name) != sha256(bundle):
        note("references", "fail", f"{bundle.name} is not listed in, or does not match, {sums}")
        return None
    with tarfile.open(bundle) as tar:
        for m in tar.getmembers():
            p = Path(m.name)
            if p.is_absolute() or ".." in p.parts or not (m.isfile() or m.isdir()) or p.parts[0] != "references":
                note("references", "fail", f"unsafe or unexpected member in the bundle: {m.name}")
                return None
        if sys.version_info >= (3, 12):
            tar.extractall(root, filter="data")
        else:
            tar.extractall(root)
    note("references", "ok", f"unpacked {bundle.name} into {root}")
    return refs


# 4 -------------------------------------------------------------------------------------
def gromacs(a):
    cands = [Path(a.gmx).expanduser()] if a.gmx else [Path(p) for p in map(shutil.which, ("gmx", "gmx_mpi")) if p]
    for g in cands:
        r = run([g, "--version"])
        m = re.search(r"GROMACS version:\s+(\S+)", r.stdout + r.stderr)
        if r.returncode == 0 and m:
            full = m.group(1)
            doc = re.match(r"\d{4}(\.\d+)?", full).group(0)
            note("gromacs", "ok" if doc == TESTED_GMX else "warn",
                 f"{g.resolve()}: GROMACS {full}" + ("" if doc == TESTED_GMX else
                 f" (the bundled MDP and index locators are for {TESTED_GMX}; step 9 re-checks the MDP)"))
            return g.resolve(), doc
    note("gromacs", "fail", "GROMACS not found: install it and pass --gmx /path/to/gmx")
    return None, None


# 5 -------------------------------------------------------------------------------------
def configure(a, root, py, gmx, gmx_doc):
    path = root / "environment.json"
    if not path.exists():
        shutil.copy(SKILL_SRC / "assets" / "environment.example.json", path)
        note("environment.json", "ok", f"created {path} from the template")
    cfg = json.loads(path.read_text())
    changed = []

    def put(key, value, placeholders, sub=None):
        d = cfg[sub] if sub else cfg
        if value is not None and (a.reconfigure or d.get(key) in placeholders):
            if d.get(key) != value:
                d[key] = value
                changed.append(f"{sub + '.' if sub else ''}{key}")
    put("python_command", [str(py)] if py else None, (None, [], ["python3"], ["python"]))
    put("gromacs_command", [str(gmx)] if gmx else None, (None, [], ["gmx"]))
    if gmx_doc and cfg.get("gromacs_documentation_version") != gmx_doc:
        cfg["gromacs_documentation_version"] = gmx_doc
        changed.append("gromacs_documentation_version")
    ff_root = root / "references" / "original" / "forcefield"
    for key, name in FF_FILES.items():
        hits = sorted(ff_root.rglob(name)) if ff_root.exists() else []
        put(key, str(hits[0].relative_to(root)) if hits else None, (None, ""), sub="force_field")
    path.write_text(json.dumps(cfg, indent=2) + "\n")
    note("environment.json", "ok" if changed else "skip",
         f"set {', '.join(changed)}" if changed else "nothing to change (use --reconfigure to overwrite)")
    missing = [f"{k}: {v}" for k, v in list(cfg["force_field"].items()) + [("reference_index", cfg["reference_index"])]
               if k != "version" and (not v or not (root / v).exists())]
    if missing:
        note("environment.json", "fail", "paths not found: " + "; ".join(missing))
    return cfg


# 6, 7 ----------------------------------------------------------------------------------
def pdf_text(py, pdf, txt):
    if txt.exists() or not pdf.exists():
        return None
    pt = tool(py, "pdftotext")
    if pt is None:
        return "pdftotext not found (install poppler)"
    r = run([pt, "-layout", pdf, txt])
    return None if r.returncode == 0 else f"pdftotext failed on {pdf.name}"


def downloads(a, root, py, gmx_doc):
    refs = root / "references"
    fetch = SKILL_SRC / "scripts" / "fetch_reference.py"
    man = json.loads(MANIFEST.read_text())
    si_url = next(s["acquire"]["url"] for s in man["sources"] if s["id"] == "polyply-si")
    jobs = []
    if gmx_doc:
        base = f"https://manual.gromacs.org/documentation/{gmx_doc}/reference-manual/"
        jobs += [(base + rel, refs / "original" / "gromacs" / name, "html") for name, rel in GMX_PAGES.items()]
        jobs.append((f"https://ftp.gromacs.org/manual/manual-{gmx_doc}.pdf",
                     refs / "original" / "gromacs" / f"manual-{gmx_doc}.pdf", "pdf"))
    jobs.append((si_url, refs / "original" / "papers" / "grunewald2022-polyply" / "supplementary-information.pdf", "pdf"))
    for url, dest, kind in jobs:
        if dest.exists():
            note("downloads", "skip", f"{dest.relative_to(root)} present")
            continue
        if a.no_fetch:
            note("downloads", "todo", f"{dest.relative_to(root)} missing (fetch skipped)")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        r = run([py, fetch, url, dest, "--kind", kind])
        note("downloads", "ok" if r.returncode == 0 else "fail",
             f"{dest.relative_to(root)}" + ("" if r.returncode == 0 else f": {(r.stdout + r.stderr).strip()[-200:]}"))
    for pdf, txt in ([(refs / "original/gromacs" / f"manual-{gmx_doc}.pdf", refs / "derived" / f"gromacs-manual-{gmx_doc}.txt")]
                     if gmx_doc else []) + [(refs / "original/papers/grunewald2022-polyply/supplementary-information.pdf",
                                             refs / "derived/grunewald2022-polyply-si.txt")]:
        err = pdf_text(py, pdf, txt)
        if err:
            note("downloads", "warn", err)


def papers(root, py):
    p = root / "references" / "original" / "papers"
    for folder, (art, si, stem) in PAPERS.items():
        need = [p / folder / art, p / folder / si]
        missing = [str(f.relative_to(root)) for f in need if not f.is_file()]
        if missing:
            note("papers", "todo", "download from the publisher and place: " + ", ".join(missing))
            continue
        for f, tag in ((need[0], "article"), (need[1], "si")):
            err = pdf_text(py, f, root / "references" / "derived" / f"{stem}-{tag}.txt")
            if err:
                note("papers", "warn", err)
        note("papers", "ok" if (p / folder / "sources.json").is_file() else "todo",
             f"{folder}: PDFs present" + ("" if (p / folder / "sources.json").is_file()
                                         else "; sources.json still to write (prepare-mode review)"))


# 8 -------------------------------------------------------------------------------------
def mdp_check(a, root, cfg, gmx):
    if gmx is None:
        return
    ff = cfg["force_field"]
    paths = {k: root / ff[k] for k in ("core_itp", "solvent_itp", "solvent_coordinates", "production_mdp") if ff.get(k)}
    if len(paths) < 4 or not all(p.exists() for p in paths.values()):
        note("mdp check", "fail", "force-field, water-box or MDP path missing in environment.json")
        return
    wd = root / "prepare" / f"setup-{dt.date.today():%Y%m%d}" / "mdp_check"
    wd.mkdir(parents=True, exist_ok=True)
    gro = paths["solvent_coordinates"].read_text().splitlines()
    n = int(gro[1].split()[0])
    (wd / "water.top").write_text(f'#include "{paths["core_itp"]}"\n#include "{paths["solvent_itp"]}"\n\n'
                                  f"[ system ]\nwater\n\n[ molecules ]\nW {n}\n")
    lines = [l for l in paths["production_mdp"].read_text().splitlines()
             if l.split(";")[0].split("=")[0].strip().replace("_", "-").lower() != "nsteps"]
    (wd / "check.mdp").write_text("\n".join(lines + ["nsteps = 500  ; setup check: 10 ps at 20 fs"]) + "\n")
    r = run([gmx, "grompp", "-f", "check.mdp", "-c", paths["solvent_coordinates"], "-p", "water.top",
             "-o", "check.tpr", "-maxwarn", "0"], cwd=wd)
    (wd / "grompp.log").write_text(r.stdout + r.stderr)
    notes = len(re.findall(r"^NOTE \d+", r.stdout + r.stderr, re.M))
    if r.returncode:
        note("mdp check", "fail", f"grompp -maxwarn 0 failed with {paths['production_mdp'].name}; see {wd / 'grompp.log'}")
        return
    r = run([gmx, "mdrun", "-deffnm", "check", "-nt", str(a.threads)], cwd=wd)
    (wd / "mdrun.log").write_text(r.stdout + r.stderr)
    note("mdp check", "ok" if r.returncode == 0 else "fail",
         f"grompp passed ({notes} note(s)), 10 ps run " + ("completed" if r.returncode == 0 else "failed")
         + f"; logs in {wd.relative_to(root)}")


# 9, 10 --------------------------------------------------------------------------------
def checks(a, root, py):
    refs = root / "references"
    r = run([py, SKILL_SRC / "scripts" / "build_index.py", "check", refs])
    note("index", "ok" if r.returncode == 0 else "warn",
         "build_index.py check: identical" if r.returncode == 0 else
         "build_index.py check reports differences (expected with another RDKit version); review them in prepare mode")
    r = run([py, SKILL_SRC / "scripts" / "reference_manifest.py", "verify", MANIFEST, refs])
    rows = [l.split(None, 2) for l in r.stdout.splitlines() if not l.startswith((" ", "files"))]
    bundled_bad = [x[0] for x in rows if len(x) == 3 and x[1] == "bundled" and x[2] != "identical"]
    site = {x[0]: x[2] for x in rows if len(x) == 3 and x[1] == "site"}
    note("manifest", "fail" if bundled_bad else "ok",
         ("bundled sources changed: " + ", ".join(bundled_bad)) if bundled_bad else
         "bundled sources identical; site sources: " + ", ".join(f"{k} {v}" for k, v in site.items()))
    if not a.lock:
        return
    lock = root / "config" / "references.lock.json"
    blockers = [s for s, st, _ in RESULTS if st in ("fail", "todo")]
    if not (refs / "index" / "readiness.json").is_file():
        blockers.append("references/index/readiness.json")
    if blockers:
        note("lock", "todo", "not created; first resolve: " + ", ".join(dict.fromkeys(blockers)))
        return
    if lock.exists():
        note("lock", "skip", f"{lock.relative_to(root)} exists; make a new snapshot rather than relocking")
        return
    lt = SKILL_SRC / "scripts" / "reference_lock.py"
    ok = run([py, lt, "create", refs, lock]).returncode == 0 and run([py, lt, "verify", refs, lock]).returncode == 0
    note("lock", "ok" if ok else "fail", f"created and verified {lock.relative_to(root)}" if ok else "lock failed")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.environ.get("MARTINI_AGENT_ROOT", "~/martini-agent-lab"),
                    help="workspace (default $MARTINI_AGENT_ROOT or ~/martini-agent-lab)")
    ap.add_argument("--bundle", help="martini-agent-references-*.tar.gz from the release")
    ap.add_argument("--sums", help="SHA256SUMS for the bundle (default: beside the bundle)")
    ap.add_argument("--env-name", default="martini-agent", help="conda env to create or use")
    ap.add_argument("--python", help="use this interpreter instead of a conda env")
    ap.add_argument("--gmx", help="GROMACS executable (default: gmx or gmx_mpi on PATH)")
    ap.add_argument("--skills-dir", default="~/.claude/skills", help="agent skills directory")
    ap.add_argument("--threads", type=int, default=4, help="mdrun threads for the MDP check")
    ap.add_argument("--no-fetch", action="store_true", help="skip network downloads")
    ap.add_argument("--reconfigure", action="store_true", help="overwrite values already set in environment.json")
    ap.add_argument("--lock", action="store_true", help="create the site's reference lock when everything is in place")
    a = ap.parse_args()
    root = Path(a.root).expanduser().resolve()
    root.mkdir(parents=True, exist_ok=True)

    py = python_env(a)
    if py is None:
        return summary(root)
    skill(a, py)
    refs = references(a, root)
    gmx, gmx_doc = gromacs(a)
    if refs is None:
        return summary(root)
    cfg = configure(a, root, py, gmx, gmx_doc)
    downloads(a, root, py, gmx_doc)
    papers(root, py)
    mdp_check(a, root, cfg, gmx)
    checks(a, root, py)
    return summary(root)


def summary(root):
    fails = [r for r in RESULTS if r[1] == "fail"]
    todos = [r for r in RESULTS if r[1] == "todo"]
    print("\n== Summary ==")
    print(f"{sum(r[1] in ('ok', 'skip') for r in RESULTS)} done, {len(todos)} to do, {len(fails)} failed")
    for step, st, msg in fails + todos:
        print(f"  {st.upper():4s} {step}: {msg}")
    print(f"\nNext: export MARTINI_ENV={root / 'environment.json'}; cd {root}; claude\n"
          "then ask for a prepare-mode review (README, 'Generating a topology', step 4), write "
          "references/index/readiness.json, and re-run this script with --lock.")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
