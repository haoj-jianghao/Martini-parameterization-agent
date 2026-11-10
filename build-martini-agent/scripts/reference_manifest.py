#!/usr/bin/env python3
"""Compare a site's references with the release reference manifest.

  reference_manifest.py verify MANIFEST.json REFERENCES

For each manifest source, reports identical, absent (no file present) or differs
(some files changed or missing), then lists files not in the manifest. Exit status
0 unless a source is partly present or changed. Reads files only; standard library.
Differences in publisher downloads, web pages or another GROMACS version are
expected; record them in readiness.json. A match is not a readiness review.
"""
import hashlib
import json
import sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def verify(manifest_path, refs):
    refs = Path(refs)
    with open(manifest_path) as f:
        m = json.load(f)
    known = set()
    status = 0
    for s in m["sources"]:
        missing, changed = [], []
        for rel, v in s["files"].items():
            known.add(rel)
            p = refs / rel
            if not p.is_file():
                missing.append(rel)
            elif sha256(p) != v["sha256"]:
                changed.append(rel)
        if not missing and not changed:
            state = "identical"
        elif len(missing) == len(s["files"]):
            state = "absent"
        else:
            state = f"differs ({len(changed)} changed, {len(missing)} missing)"
            status = 1
        flag = "bundled" if s.get("bundle") else "site"
        print(f"{s['id']:28s} {flag:8s} {state}")
        for rel in (changed + missing)[:5]:
            print(f"    {rel}")
    extra = sorted(str(p.relative_to(refs)) for p in refs.rglob("*")
                   if p.is_file() and ".git" not in p.parts and "__pycache__" not in p.parts
                   and str(p.relative_to(refs)) not in known)
    if extra:
        print(f"files not in the manifest: {len(extra)}")
        for rel in extra[:10]:
            print(f"    {rel}")
    return status


if __name__ == "__main__":
    if len(sys.argv) == 4 and sys.argv[1] == "verify":
        sys.exit(verify(sys.argv[2], sys.argv[3]))
    sys.exit(__doc__)
