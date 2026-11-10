#!/usr/bin/env python3
"""Build the reference index (prepare mode, references/environment.md section 7).

Writes REFERENCE_ROOT/index/models.jsonl and REFERENCE_ROOT/index/sections.json from the fetched references and
the reviewed curation decisions in assets/index_curation.json. Deterministic: the same reference bytes, curation
file and RDKit version give the same index. Uses only the standard library and RDKit; never runs downloaded code.

Sources indexed (paths under REFERENCE_ROOT):
  original/small-molecules/LIBRARY.md and models/   one entry per published model variant (SM3:<uID>:opt|cog)
  original/forcefield/martini_v300/*solvents*, *ions*, *nucleobases* one entry per moleculetype (FF300:<kind>:<name>)
  original/bible/src/MartiniCGSmilesDB_v6.xlsx       SMILES and CGSmiles for solvent and ion identities
  curation extra_entries                             hand-reviewed entries (for example polyply PEO)

Each SMILES is checked: library SMILES against the element formula of the model's AA topology; solvent and ion
SMILES against the ITP charge sum. Failures are kept in the entry and listed in the summary; they are not hidden.
Paths in the index are written relative to the parent of REFERENCE_ROOT, as environment.json expects.

usage:
  build_index.py build REFERENCE_ROOT [--curation FILE] [--force]   write index/models.jsonl and index/sections.json
  build_index.py check REFERENCE_ROOT [--curation FILE]             rebuild in memory and diff against the files on disk
Exit status of check: 0 identical, 1 different. After build, review the summary and write index/readiness.json.
"""
import argparse
import hashlib
import json
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase

RDLogger.DisableLog("rdApp.*")  # parse failures are recorded in smiles_verification

SKILL = Path(__file__).resolve().parents[1]
MT_RE = re.compile(r"\[\s*moleculetype\s*\][^\n]*\n(?:\s*;[^\n]*\n)*\s*(\S+)")
MASS2EL = [(1.008, "H"), (12.011, "C"), (14.007, "N"), (15.999, "O"), (18.998, "F"), (32.06, "S"), (35.45, "Cl"),
           (79.904, "Br"), (126.904, "I"), (30.974, "P")]
XNS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def aa_formula(itp):
    """Element counts of an AA topology's [ atoms ], elements assigned from the mass column."""
    counts, on = {}, False
    for line in open(itp):
        s = line.split(";")[0].strip()
        if s.startswith("["):
            on = s.replace(" ", "") == "[atoms]"
            continue
        if on and s:
            m = float(s.split()[7])
            el = min(MASS2EL, key=lambda x: abs(x[0] - m))[1]
            counts[el] = counts.get(el, 0) + 1
    return counts


def smi_formula(mol):
    c = {}
    for a in Chem.AddHs(mol).GetAtoms():
        c[a.GetSymbol()] = c.get(a.GetSymbol(), 0) + 1
    return c


def itp_charge(txt, start):
    """Charge sum of the [ atoms ] of the moleculetype that starts at txt[start]."""
    q, on = 0.0, False
    for line in txt[start:].splitlines()[1:]:
        s = line.split(";")[0].strip()
        if s.startswith("[") and "moleculetype" in s:
            break
        if s.startswith("["):
            on = s.replace(" ", "") == "[atoms]"
            continue
        if on and s:
            q += float(s.split()[6])
    return round(q, 3)


def header_dois(itp):
    return sorted(set(re.findall(r"DOI:\s*([0-9.]+/[^\s\]\),;]+)", open(itp, errors="replace").read())))


def read_bible(xlsx):
    """First row per MOLTAG of the Bible spreadsheet: [name, MOLTAG, SMILES, CGSmiles, ...]."""
    z = zipfile.ZipFile(xlsx)
    ss = ["".join(t.text or "" for t in si.iter(XNS + "t"))
          for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(XNS + "si")]
    out = {}
    for r in ET.fromstring(z.read("xl/worksheets/sheet1.xml")).iter(XNS + "row"):
        vals = []
        for c in r.findall(XNS + "c"):
            v = c.find(XNS + "v")
            vals.append(ss[int(v.text)] if (v is not None and c.get("t") == "s") else (v.text if v is not None else ""))
        if len(vals) >= 4 and vals[1] != "MOLTAG":
            out.setdefault(vals[1], vals)
    return out


class Builder:
    def __init__(self, ref, cur):
        self.ref, self.cur = ref, cur
        self.base = ref.parent
        self.motif_pat = {k: Chem.MolFromSmarts(v) for k, v in cur["motif_smarts"].items()}
        bad = [k for k, p in self.motif_pat.items() if p is None]
        if bad:
            raise SystemExit(f"invalid motif SMARTS in curation: {bad}")

    def rel(self, p):
        return str(Path(p).relative_to(self.base))

    def motifs(self, mol):
        return [k for k, p in self.motif_pat.items() if mol.HasSubstructMatch(p)]

    # ---------------------------------------------------------------- small-molecule library
    def small_molecules(self):
        sm = self.ref / "original/small-molecules"
        cfg = self.cur["small_molecules"]
        proposed = cfg["proposed_smiles"]
        rows = []
        for line in open(sm / "LIBRARY.md"):
            if not line.startswith("|") or "uID" in line or "---" in line:
                continue
            cells = [c.strip() for c in line.strip().strip("|").split("|")]
            rows.append((cells[0], cells[1], cells[2].strip("`"), re.findall(r"\[([^\]]+)\]\(\./([^)]+)\)", cells[3])))

        def proposed_check(uid, aa):
            if uid not in proposed:
                return None
            m = Chem.MolFromSmiles(proposed[uid])
            ok = m is not None and aa.exists() and smi_formula(m) == aa_formula(aa)
            return "formula matches AA itp; adopted as canonical (user decision)" if ok else \
                "formula does NOT match AA itp; not adopted"

        out = []
        for uid, name, smi, links in rows:
            mol = Chem.MolFromSmiles(smi)
            can = Chem.MolToSmiles(mol, isomericSmiles=True) if mol else None
            aa = sm / "models/AA/oplsaa-LigParGen" / f"{uid}_LigParGen.itp"
            mapdir = sm / "models/mapping-cgbuilder/oplsaa-LigParGen" / uid
            if mol and aa.exists():
                f_aa, f_smi = aa_formula(aa), smi_formula(mol)
                ver = "formula_match_AA_itp" if f_aa == f_smi else f"FORMULA_MISMATCH smiles={f_smi} aa={f_aa}"
            else:
                ver = "unverified (no AA reference itp)" if mol else "LIBRARY_SMILES_INVALID (RDKit parse/kekulize failure)"
            pc = proposed_check(uid, aa)
            if ver.startswith("formula_match"):
                canon = can
            elif uid in proposed and "adopted as canonical" in str(pc):
                canon = Chem.MolToSmiles(Chem.MolFromSmiles(proposed[uid]))
            else:
                canon = None
            nd = sorted(f.name for f in mapdir.iterdir() if f.suffix == ".ndx") if mapdir.is_dir() else []
            png = sorted(f.name for f in mapdir.iterdir() if f.suffix == ".png") if mapdir.is_dir() else []
            for label, path in links:
                p = sm / path
                if not p.exists():
                    out.append({"reference_id": f"SM3:{uid}:{label}", "name": name, "error": f"missing file {path}"})
                    continue
                mt = MT_RE.search(p.read_text())
                dois = header_dois(p)
                out.append({
                    "reference_id": f"SM3:{uid}:{'opt' if 'optimized' in label else 'cog'}", "name": name,
                    "canonical_isomeric_smiles": canon,
                    "library_smiles": smi, "smiles_verification": ver,
                    "proposed_smiles": proposed.get(uid), "proposed_smiles_check": pc,
                    "formal_charge": Chem.GetFormalCharge(mol) if mol else None,
                    "source_path": self.rel(p), "sha256": sha(p), "doi_url": dois or [cfg["default_doi"]],
                    "doi_note": None if dois else "no DOI in ITP header; library paper (Alessandri 2022) assumed as collection source",
                    "source_locator": f"moleculetype {mt.group(1) if mt else '?'} in {p.name}; LIBRARY.md row {uid}",
                    "mapping_path": self.rel(mapdir / nd[0]) if nd else None,
                    "mapping_figure": self.rel(mapdir / png[0]) if png else None,
                    "mapping_atom_reference": self.rel(aa) if aa.exists() else None,
                    "forcefield_version": "Martini 3.0.0 (M3-Small-Molecules repo, commit in original/small-molecules.commit.txt)",
                    "motifs": self.motifs(mol) if mol else [], "model_variant": label,
                    "notes": "mapping .ndx indices are 1-based atoms of mapping_atom_reference; repeated indices = weights",
                    "role": "generation"})
        return out

    # ---------------------------------------------------------------- solvents and ions of the force-field release
    def forcefield(self):
        cfg = self.cur["forcefield"]
        bible = read_bible(self.ref / "original/bible/src/MartiniCGSmilesDB_v6.xlsx")
        ff = self.ref / "original/forcefield/martini_v300"
        out = []
        for fname, kind in [("martini_v3.0.0_solvents_v1.itp", "solvent"), ("martini_v3.0.0_ions_v1.itp", "ion"),
                            ("martini_v3.0.0_nucleobases_v1.itp", "nucleobase")]:
            p = ff / fname
            txt, h = p.read_text(), sha(p)
            for m in MT_RE.finditer(txt):
                mt = m.group(1)
                prev = [l for l in txt[:m.start()].splitlines() if l.startswith(";;;;;;")]
                cname = cfg["user_names"].get(mt, prev[-1].strip(";").strip() if prev else mt)
                qitp = itp_charge(txt, m.start())
                b = bible.get(mt)
                smi = can = None
                ver = "no SMILES source (name from ITP comment only; identity not inferred)"
                if not b and mt in cfg["agent_smiles"]:
                    mol = Chem.MolFromSmiles(cfg["agent_smiles"][mt])
                    ok = mol is not None and abs(Chem.GetFormalCharge(mol) - qitp) < 1e-6
                    smi, can = cfg["agent_smiles"][mt], (Chem.MolToSmiles(mol, isomericSmiles=True) if ok else None)
                    ver = (f"agent-assigned from ITP comment name '{cname}'; RDKit parse OK; formal charge "
                           f"{Chem.GetFormalCharge(mol):+d} equals ITP charge sum {qitp:+g}"
                           if ok else f"AGENT_SMILES_CHECK_FAILED (ITP charge {qitp})")
                    if mt in cfg.get("unspecified_double_bond_stereo", []):
                        ver += "; double-bond stereo unspecified in source name"
                    if mt in cfg.get("agent_smiles_notes", {}):
                        ver += "; " + cfg["agent_smiles_notes"][mt]
                elif b:
                    mol = Chem.MolFromSmiles(b[2])
                    if mol:
                        smi, can = b[2], Chem.MolToSmiles(mol, isomericSmiles=True)
                        ver = (f"SMILES from Bible xlsx row MOLTAG={mt} ('{b[0]}'); name checked by agent against ITP "
                               f"comment '{cname}'; formal charge {Chem.GetFormalCharge(mol):+d} vs ITP {qitp:+g}")
                        if mt in cfg["user_names"]:
                            ver += "; identity per user decision"
                out.append({
                    "reference_id": f"FF300:{kind}:{mt}", "name": cname, "canonical_isomeric_smiles": can,
                    "library_smiles": smi, "smiles_verification": ver,
                    "formal_charge": Chem.GetFormalCharge(Chem.MolFromSmiles(smi)) if smi else None,
                    "source_path": self.rel(p), "sha256": h, "doi_url": [cfg["doi"]],
                    "source_locator": f"moleculetype {mt} in {fname}", "itp_charge_sum": qitp,
                    "mapping_path": None, "mapping_figure": None,
                    "mapping_cgsmiles": b[3] if b else None,
                    "forcefield_version": "Martini 3.0.0 (martini_v300.zip)",
                    "motifs": self.motifs(Chem.MolFromSmiles(smi)) if smi else [],
                    "notes": "topology without an atom-level mapping file; exact motif transfer only via Bible CGSmiles if present",
                    "role": "generation"})
        return out

    def extra(self):
        out = []
        for e in self.cur.get("extra_entries", []):
            p = self.ref / e["source_path"]
            if not p.exists():
                print(f"extra entry {e['reference_id']}: {e['source_path']} absent, skipped", file=sys.stderr)
                continue
            d = {}
            for k, v in e.items():
                d[k] = self.rel(p) if k == "source_path" else v
                if k == "source_path":
                    d["sha256"] = sha(p)
            out.append(d)
        return out

    def models(self):
        return self.small_molecules() + self.forcefield() + self.extra()


def models_text(entries):
    return "".join(json.dumps(e, ensure_ascii=False) + "\n" for e in entries)


def sections_text(cur):
    return json.dumps(cur["sections"], indent=2, ensure_ascii=False) + "\n"


def summary(entries):
    v = lambda e: str(e.get("smiles_verification"))
    return {"n_entries": len(entries), "rdkit": rdBase.rdkitVersion,
            "small_molecule_entries": sum(e["reference_id"].startswith("SM3:") for e in entries),
            "solvent_ion_entries": sum(e["reference_id"].startswith("FF300:") for e in entries),
            "solvent_ion_with_smiles": sum(1 for e in entries if e["reference_id"].startswith("FF300:") and e.get("canonical_isomeric_smiles")),
            "extra_entries": sum(1 for e in entries if not e["reference_id"].startswith(("SM3:", "FF300:"))),
            "formula_mismatches": [e["reference_id"] for e in entries if "MISMATCH" in v(e)],
            "failed_checks": [e["reference_id"] for e in entries
                              if ("FAILED" in v(e) or "INVALID" in v(e)) and not e.get("canonical_isomeric_smiles")],
            "replaced_by_curation": [e["reference_id"] for e in entries if e.get("proposed_smiles") and e.get("canonical_isomeric_smiles")],
            "missing_files": [e["reference_id"] for e in entries if "error" in e],
            "null_smiles": [e["reference_id"] for e in entries if not e.get("canonical_isomeric_smiles")]}


def diff_models(new, path):
    old = {}
    for l in open(path):
        e = json.loads(l)
        old[e["reference_id"]] = e
    newd = {e["reference_id"]: e for e in new}
    rep = {"only_on_disk": sorted(set(old) - set(newd)), "only_rebuilt": sorted(set(newd) - set(old)), "changed": {}}
    for k in sorted(set(old) & set(newd)):
        ch = sorted(f for f in set(old[k]) | set(newd[k]) if old[k].get(f) != newd[k].get(f))
        if ch:
            rep["changed"][k] = ch
    return rep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["build", "check"])
    ap.add_argument("reference_root", type=Path)
    ap.add_argument("--curation", type=Path, default=SKILL / "assets/index_curation.json")
    ap.add_argument("--force", action="store_true", help="overwrite an existing index (build)")
    a = ap.parse_args()
    ref = a.reference_root.resolve()
    cur = json.loads(a.curation.read_text())
    entries = Builder(ref, cur).models()
    mtext, stext = models_text(entries), sections_text(cur)
    idx = ref / "index"
    mp, sp = idx / "models.jsonl", idx / "sections.json"
    info = summary(entries)
    info["curation"] = {"path": str(a.curation), "sha256": sha(a.curation)}
    if a.mode == "build":
        for p in (mp, sp):
            if p.exists() and not a.force:
                raise SystemExit(f"{p} exists; use check to compare, or --force to overwrite")
        idx.mkdir(parents=True, exist_ok=True)
        mp.write_text(mtext, encoding="utf-8")
        sp.write_text(stext, encoding="utf-8")
        info["written"] = {str(mp): sha(mp), str(sp): sha(sp)}
        print(json.dumps(info, indent=1))
        return 0
    same_m = mp.exists() and mp.read_text(encoding="utf-8") == mtext
    same_s = sp.exists() and sp.read_text(encoding="utf-8") == stext
    info["models.jsonl"] = "identical" if same_m else ("absent" if not mp.exists() else diff_models(entries, mp))
    info["sections.json"] = "identical" if same_s else ("absent" if not sp.exists() else "different")
    print(json.dumps(info, indent=1))
    return 0 if same_m and same_s else 1


if __name__ == "__main__":
    sys.exit(main())
