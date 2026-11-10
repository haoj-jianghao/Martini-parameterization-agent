#!/usr/bin/env python3
"""Reference manifest for releasing the Martini agent workflow.

  make_manifest.py build  REFERENCES LOCK OUT.json   write the manifest from a verified snapshot
  make_manifest.py verify MANIFEST.json REFERENCES   same as the skill's scripts/reference_manifest.py verify
  make_manifest.py bundle MANIFEST.json REFERENCES SKILL_DIR OUT_DIR
                                                     copy skill + manifest + bundled references

The manifest groups every locked reference file into one source (a download, a git
clone, a manual download or a file we authored) and records how to acquire it, its
licence status, whether it is bundled, and the SHA-256 and size of each file.
Standard library only. Never runs downloaded code.
"""
import hashlib
import json
import shutil
import sys
from pathlib import Path

SKILL_SCRIPTS = Path(__file__).resolve().parent.parent / "build-martini-agent" / "scripts"
sys.path.insert(0, str(SKILL_SCRIPTS))
from reference_manifest import verify  # noqa: E402  (site-side check lives in the skill)

# Source table. Files are assigned to the first source whose prefix matches.
# Licence wording: "stated by the owner" = the project owner's statement, not a
# licence file in the download.
SOURCES = [
    {"id": "martini-3.0.0-forcefield",
     "title": "Martini 3.0.0 force-field parameters (martini_v300.zip)",
     "prefixes": ["archives/martini_v300.zip", "original/forcefield/"],
     "acquire": {"method": "fetch-zip", "source_record": "archives/martini_v300.zip.source.json",
                 "extract_to": "original/forcefield"},
     "licence": "MIT or similar (stated by the owner); no licence file in the archive",
     "bundle": True},
    {"id": "tutorial-zip",
     "title": "Martini 3 small-molecule parametrization tutorial files",
     "prefixes": ["archives/parametrize_new_small_molecule.zip",
                  "original/tutorial/parametrize_new_small_molecule/"],
     "acquire": {"method": "fetch-zip",
                 "source_record": "archives/parametrize_new_small_molecule.zip.source.json",
                 "extract_to": "original/tutorial"},
     "licence": "MIT or similar (stated by the owner); no licence file in the archive",
     "bundle": True,
     "note": "Includes ENAP-in-water/6_CG-takeCURRENT/box_CG_W_eq.gro, the water box used for "
             "execution checks (environment.json solvent_coordinates)."},
    {"id": "tutorial-page",
     "title": "Tutorial web page (HTML snapshot)",
     "prefixes": ["original/tutorial/page.html"],
     "acquire": {"method": "fetch", "source_record": "original/tutorial/page.html.source.json"},
     "licence": "MIT or similar (stated by the owner)",
     "bundle": True,
     "note": "Web page; a fresh fetch will normally hash differently."},
    {"id": "tutorial-pdf",
     "title": "Tutorial page printed to PDF from a browser (image-only, 17 pages)",
     "prefixes": ["original/tutorial/tutorial.pdf"],
     "acquire": {"method": "manual", "source_record": "original/tutorial/tutorial.pdf.source.json",
                 "instructions": "Print the tutorial page to PDF from a browser; check that the "
                                 "Section 5.1 hinge diagram (Figure 3) is legible."},
     "licence": "MIT or similar (stated by the owner)",
     "bundle": True,
     "note": "A site-made print will hash differently; the bundled copy avoids that."},
    {"id": "production-mdp",
     "title": "Official Martini 3 production MDP",
     "prefixes": ["original/mdp/"],
     "acquire": {"method": "fetch", "source_record": "original/mdp/martini_v3.0_prod.mdp.source.json"},
     "licence": "MIT or similar (stated by the owner)",
     "bundle": True,
     "note": "Upstream header shows an update on 20/12/2024; the file can change without notice."},
    {"id": "gromacs-manual-html",
     "title": "GROMACS reference-manual pages (topology formats, bonded interactions, interaction methods)",
     "prefixes": ["original/gromacs/topology-file-formats.html", "original/gromacs/bonded-interactions.html",
                  "original/gromacs/interaction-methods.html"],
     "acquire": {"method": "fetch", "version_specific": True,
                 "url_template": "https://manual.gromacs.org/documentation/{GMX_VERSION}/reference-manual/..."},
     "licence": "not stated in the download",
     "bundle": False,
     "note": "Fetch at the site for the installed GROMACS version (environment.md section 6)."},
    {"id": "gromacs-manual-pdf",
     "title": "GROMACS full PDF manual",
     "prefixes": ["original/gromacs/manual-"],
     "acquire": {"method": "fetch", "version_specific": True,
                 "url_template": "https://ftp.gromacs.org/manual/manual-{GMX_VERSION}.pdf",
                 "source_record": "original/gromacs/manual-2026.3.pdf.source.json"},
     "licence": "not stated in the download",
     "bundle": False,
     "note": "Fetch at the site for the installed GROMACS version. Indexed in sections.json "
             "(Table 5.14, topology directives)."},
    {"id": "m3-bible",
     "title": "Martini 3 Bible (building-block table and sources)",
     "prefixes": ["original/bible.commit.txt", "original/bible/"],
     "acquire": {"method": "git", "url": "https://github.com/Martini-Force-Field-Initiative/M3-Bible.git",
                 "commit_file": "original/bible.commit.txt", "dest": "original/bible"},
     "licence": "MIT or similar (stated by the owner); no licence file in the repository",
     "bundle": True},
    {"id": "m3-small-molecules",
     "title": "Martini 3 small-molecule library",
     "prefixes": ["original/small-molecules.commit.txt", "original/small-molecules/"],
     "acquire": {"method": "git", "url": "https://github.com/Martini-Force-Field-Initiative/M3-Small-Molecules.git",
                 "commit_file": "original/small-molecules.commit.txt", "dest": "original/small-molecules"},
     "licence": "Apache-2.0 (original/small-molecules/LICENSE)",
     "bundle": True},
    {"id": "polyply-peo",
     "title": "Polyply Martini 3 PEO data file, with licence, README and citations",
     "prefixes": ["original/polymers-m3/"],
     "acquire": {"method": "fetch-pinned", "commit_file": "original/polymers-m3/polyply_1.0.commit.txt",
                 "url_template": "https://raw.githubusercontent.com/marrink-lab/polyply_1.0/{COMMIT}/..."},
     "licence": "Apache-2.0 (original/polymers-m3/polyply_1.0/LICENSE)",
     "bundle": True,
     "note": "Data file only; polyply itself is never run (automatic generator)."},
    {"id": "polyply-paper",
     "title": "Grunewald et al. 2022, Nat. Commun. (polyply) article",
     "prefixes": ["original/papers/grunewald2022-polyply/article.pdf",
                  "derived/grunewald2022-polyply-article.txt"],
     "acquire": {"method": "fetch",
                 "source_record": "original/papers/grunewald2022-polyply/article.pdf.source.json",
                 "doi": "10.1038/s41467-021-27627-4"},
     "licence": "CC BY 4.0 (stated in the article PDF)",
     "bundle": True},
    {"id": "polyply-si",
     "title": "Grunewald et al. 2022, Nat. Commun. (polyply) supplementary information",
     "prefixes": ["original/papers/grunewald2022-polyply/", "derived/grunewald2022-polyply-si.txt"],
     "acquire": {"method": "fetch",
                 "source_record": "original/papers/grunewald2022-polyply/supplementary-information.pdf.source.json",
                 "doi": "10.1038/s41467-021-27627-4"},
     "licence": "open access; no licence statement found in the SI file",
     "bundle": False,
     "note": "Fetch at the site (open URL, automated). Also holds the folder's sources.json."},
    {"id": "souza2021",
     "title": "Souza et al. 2021, Nat. Methods (Martini 3) article and SI",
     "prefixes": ["original/papers/Souza2021/", "derived/souza2021-"],
     "acquire": {"method": "manual", "doi": "10.1038/s41592-021-01098-3",
                 "dest": "original/papers/Souza2021/{article.pdf,supporting-information.pdf,sources.json}",
                 "instructions": "Download from the publisher; derive text with pdftotext -layout."},
     "licence": "publisher copyright",
     "bundle": False},
    {"id": "alessandri2022",
     "title": "Alessandri et al. 2022, Adv. Theory Simul. (Martini 3 small molecules) article and SI",
     "prefixes": ["original/papers/alessandri2022/", "derived/alessandri2022-"],
     "acquire": {"method": "manual", "doi": "10.1002/adts.202100391",
                 "dest": "original/papers/alessandri2022/{article.pdf,supporting-information.pdf,sources.json}",
                 "instructions": "Download from the publisher; derive text with pdftotext -layout."},
     "licence": "publisher copyright",
     "bundle": False},
    {"id": "derived-gromacs-text",
     "title": "Searchable text of the GROMACS PDF manual",
     "prefixes": ["derived/gromacs-manual-"],
     "acquire": {"method": "derive", "command": "pdftotext -layout original/gromacs/manual-{GMX_VERSION}.pdf "
                                                 "derived/gromacs-manual-{GMX_VERSION}.txt"},
     "licence": "follows the GROMACS manual",
     "bundle": False},
    {"id": "derived-authored",
     "title": "Derived text of bundled sources, adapted MDP and derived README",
     "prefixes": ["derived/"],
     "acquire": {"method": "authored",
                 "note": "bible.txt and tutorial-page.txt are derived from bundled sources; the adapted MDP "
                         "and README are ours. Each change in the MDP is marked '; ADAPTED'."},
     "licence": "ours (derived text follows its source)",
     "bundle": True,
     "note": "The adapted MDP was validated on GROMACS 2026.3; re-run grompp -maxwarn 0 at a site "
             "with another version."},
    {"id": "readiness",
     "title": "Readiness review (index/readiness.json)",
     "prefixes": ["index/readiness.json"],
     "acquire": {"method": "site-review",
                 "instructions": "Each site writes its own after prepare mode: every required resource "
                                 "with pass/missing/content-review-needed and evidence (environment.md section 7)."},
     "licence": "ours",
     "bundle": False,
     "note": "Site-specific: records this site's interpreter, GROMACS build and manual downloads."},
    {"id": "index",
     "title": "Reference index: models.jsonl, sections.json",
     "prefixes": ["index/"],
     "acquire": {"method": "authored",
                 "command": "python $SKILL_DIR/scripts/build_index.py build|check references"},
     "licence": "ours",
     "bundle": True,
     "note": "Valid only where the indexed sources match this manifest. build_index.py check reports "
             "identical with the same bytes and RDKit version."},
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def tree_digest(files):
    lines = "".join(f"{p}\t{v['sha256']}\n" for p, v in sorted(files.items()))
    return hashlib.sha256(lines.encode()).hexdigest()


def assign(path):
    for src in SOURCES:
        if any(path.startswith(p) for p in src["prefixes"]):
            return src["id"]
    return None


def read_json(path):
    with open(path) as f:
        return json.load(f)


def build(refs, lock_path, out):
    refs = Path(refs)
    lock = read_json(lock_path)["files"]
    by_src = {s["id"]: {} for s in SOURCES}
    unassigned = []
    for rel, digest in sorted(lock.items()):
        disk = sha256(refs / rel)
        if disk != digest:
            sys.exit(f"lock mismatch on disk: {rel}; verify the lock first")
        sid = assign(rel)
        if sid is None:
            unassigned.append(rel)
            continue
        by_src[sid][rel] = {"sha256": digest, "bytes": (refs / rel).stat().st_size}
    if unassigned:
        sys.exit("files not assigned to a source:\n  " + "\n  ".join(unassigned))

    sources = []
    for s in SOURCES:
        files = by_src[s["id"]]
        if not files:
            sys.exit(f"source {s['id']} has no files")
        acq = dict(s["acquire"])
        rec = acq.get("source_record")
        if rec:
            r = read_json(refs / rec)
            acq["url"] = r.get("url")
            acq["retrieved_utc"] = r.get("retrieved_utc") or r.get("retrieval_date")
        if "commit_file" in acq:
            acq["commit"] = (refs / acq["commit_file"]).read_text().strip()
        entry = {k: s[k] for k in ("id", "title", "licence", "bundle")}
        entry["acquire"] = acq
        if "note" in s:
            entry["note"] = s["note"]
        entry["file_count"] = len(files)
        entry["bytes"] = sum(v["bytes"] for v in files.values())
        entry["tree_sha256"] = tree_digest(files)
        entry["files"] = files
        sources.append(entry)

    manifest = {
        "schema_version": 1,
        "description": "Reference sources of the Martini agent workflow: how each is acquired, its "
                       "licence status, whether it is bundled, and the hash of every file. Paths are "
                       "relative to the references root.",
        "snapshot": {"lock": Path(lock_path).name, "lock_sha256": sha256(lock_path),
                     "file_count": len(lock),
                     "gromacs_documentation_version": "2026.3"},
        "excluded_from_lock": [".git metadata", "__pycache__"],
        "sources": sources,
    }
    Path(out).write_text(json.dumps(manifest, indent=1) + "\n")
    nb = sum(s["bytes"] for s in sources if s["bundle"])
    print(f"wrote {out}: {len(sources)} sources, {len(lock)} files; bundled {nb / 1e6:.1f} MB")


def bundle(manifest_path, refs, skill_dir, out_dir):
    refs, skill_dir, out_dir = Path(refs), Path(skill_dir), Path(out_dir)
    if out_dir.exists():
        sys.exit(f"{out_dir} exists; choose a new directory")
    m = read_json(manifest_path)
    shutil.copytree(skill_dir, out_dir / skill_dir.name,
                    ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    shutil.copy2(manifest_path, out_dir / "reference_manifest.json")
    n = 0
    for s in m["sources"]:
        if not s["bundle"]:
            continue
        for rel, v in s["files"].items():
            dst = out_dir / "references" / rel
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(refs / rel, dst)
            if sha256(dst) != v["sha256"]:
                sys.exit(f"copy mismatch: {rel}")
            n += 1
    print(f"bundle in {out_dir}: skill, manifest, {n} reference files")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["build"] and len(a) == 4:
        build(*a[1:])
    elif a[:1] == ["verify"] and len(a) == 3:
        sys.exit(verify(*a[1:]))
    elif a[:1] == ["bundle"] and len(a) == 5:
        bundle(*a[1:])
    else:
        sys.exit(__doc__)
