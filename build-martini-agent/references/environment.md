# Prepare a Martini agent environment

This is the acquisition and placement procedure, not a claim that references are bundled. Download once, review, lock and then generate offline. Do not execute scripts from downloaded repositories without inspecting them. Preserve original files, citations and licenses. A checksum records the acquired bytes; it is not a publisher signature or proof of scientific accuracy.

## 0. Starting from the release bundle

The release bundle holds this skill, `reference_manifest.json` and a `references/` tree with the openly licensed sources: the force field, tutorial (archive, page and PDF), production MDP, Bible, small-molecule library, polyply PEO data file, polyply article, our derived text and adapted MDP, and the index (`models.jsonl`, `sections.json`). Each manifest source records how it was acquired, its licence status, whether it is bundled, and the SHA-256 of every file.

With the bundle, do section 1 (software), place `references/` under `$MARTINI_AGENT_ROOT`, then do only the site steps:

1. Fetch the GROMACS manual pages and PDF for the installed version and derive the PDF text (section 6).
2. Fetch the polyply SI and download the Souza 2021 and Alessandri 2022 articles and SI; derive text (section 5).
3. If GROMACS is not the version named in the adapted MDP's file name, re-check it (section 6).
4. Configure, check the index, write readiness and lock (section 7).

The water box for execution checks is the tutorial's, already in the bundle (section 2).

Before locking, compare the site's references with the manifest:

```bash
python "$SKILL_DIR/scripts/reference_manifest.py" verify \
  /path/to/reference_manifest.json "$MARTINI_AGENT_ROOT/references"
```

It reports each source as `identical`, `absent` or `differs`. Bundled sources must be identical. Pinned git commits and the polyply files should be identical. Publisher downloads, web pages and a manual for another GROMACS version may differ, which is expected; record which sources differ in `readiness.json`. Do not edit bundled files to match. The site still creates its own lock (section 7). Without the bundle, follow sections 1 to 7 in full.

## 1. Workspace and software

Suggested Linux/WSL location: `$HOME/martini-agent-lab`. Keep it separate from the installed skill. On Windows, use the WSL Linux filesystem for simulation work. Python 3.10+ is required by the provided helpers.

Use an existing environment with RDKit if available. Otherwise one possible setup is:

```bash
conda create -n martini-agent -c conda-forge python=3.11 rdkit numpy scipy matplotlib pyyaml
conda activate martini-agent
python -c 'from rdkit import Chem; print(Chem.MolFromSmiles("CCO") is not None)'
gmx --version
```

Record the absolute path of the interpreter that imports RDKit (for example a conda env's `bin/python`) as `python_command` in environment.json; the system `python3` often lacks RDKit. Use that interpreter for every helper command below.

Install/locate GROMACS separately if `gmx` is unavailable. Poppler's `pdftotext`, `pdfinfo` and `pdftoppm` are useful for paper inspection. Record actual versions; do not infer GROMACS version from this guide's example documentation URLs.

`scripts/cgbuild.py search` reads the Bible spreadsheet, so `python_command` also needs openpyxl.

Bonded terms (stage 4) run with `python_command`, which also needs NumPy and SciPy. They need an OPLS-AA GROMACS topology of each molecule, which the user makes elsewhere with any OPLS engine (for example opls-ff-gen, opls-ff-lite or the LigParGen server) and passes to `cgbuild.py input --aa-topology FILE`. A `.top` (with `#include` files next to it) or a self-contained `.itp` with `[ atomtypes ]` works; a bare `.itp` without `[ defaults ]` gets the OPLS-AA defaults. The topology must describe the SMILES species, protonation state included; any charge method is accepted (`references/bonded.md`). No OPLS engine is needed at the site.

Batches run one agent session per molecule through `scripts/batch.py`, which has built-in command lines for `claude` (Claude Code, headless `-p`) and `codex` (`codex exec`). A site may add or override them in environment.json as `agent_commands`, a map from a name to an argv template with the fields `{prompt}`, `{model}`, `{moldir}`, `{project}` and `{budget}`. The agent CLI must be installed and signed in by the user.

Optionally, a site with an OPLS/2020 engine installed may record it in environment.json as `opls_ff_gen_command` (an argv list taking `--smiles`, `--name`, `--output-dir`, `--optimize` and `--charge-model`, as opls-ff-gen and opls-ff-lite do). It is used only for molecules given without `--aa-topology`; leave it `null` otherwise. The energy model is NumPy code in the skill (`scripts/d5bonded/energy.py`); OpenMM is not used, and no engine runs here unless `opls_ff_gen_command` is set. The energy calculations run on CPU, and no dynamics is run.

Set `SKILL_DIR` to the directory containing this installed skill's SKILL.md (or its complete copy in a local agent's skills directory):

```bash
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
export SKILL_DIR="/absolute/path/to/build-martini-agent"
mkdir -p "$MARTINI_AGENT_ROOT/config" "$MARTINI_AGENT_ROOT/runs"
mkdir -p "$MARTINI_AGENT_ROOT/references/original/papers" "$MARTINI_AGENT_ROOT/references/original/mdp"
mkdir -p "$MARTINI_AGENT_ROOT/references/original/gromacs" "$MARTINI_AGENT_ROOT/references/archives"
mkdir -p "$MARTINI_AGENT_ROOT/references/derived" "$MARTINI_AGENT_ROOT/references/index"
cp "$SKILL_DIR/assets/environment.example.json" "$MARTINI_AGENT_ROOT/environment.json"
```

Do not leave the illustrative SKILL_DIR path unchanged. If porting to another agent harness, copy the complete skill folder, including references/assets/scripts; a single SKILL.md is not self-contained. Refer to that harness's skill-directory convention. The skill's three modes do not depend on ChatGPT-specific APIs.

## 2. Obtain the official numerical release

Authoritative landing page: https://cgmartini.nl/docs/downloads/force-field-parameters/martini3/particle-definitions.html

```bash
python "$SKILL_DIR/scripts/fetch_reference.py" \
  https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/ff_parameters/martini3/martini_v300.zip \
  "$MARTINI_AGENT_ROOT/references/archives/martini_v300.zip" --kind zip \
  --extract "$MARTINI_AGENT_ROOT/references/original/forcefield"
```

Inspect actual archive structure rather than assuming internal paths:

```bash
rg --files "$MARTINI_AGENT_ROOT/references/original/forcefield" -g '*.itp' -g '*.gro' -g '*.top'
```

Select `martini_v3.0.0.itp` and compatible solvent/ion ITPs from the actual listing. Populate environment.json with exact paths relative to the project root. Identify solvent molecule names from their moleculetype sections. Preserve any required includes. Do not import all molecular ITP files into every test system.

For solvent coordinates, use the equilibrated water box in the tutorial archive (step 4): `original/tutorial/parametrize_new_small_molecule/ENAP-in-water/6_CG-takeCURRENT/box_CG_W_eq.gro`, 2704 W beads in a 6.858 nm cubic box (about 1002 kg/m3). The tutorial's `run_1mol_CG_system.sh` uses it with `gmx solvate -cs`, as `cgbuild.py test` does. Set `solvent_coordinates` to that path (the template already does), and check that its residue and atom names match the `W` moleculetype of the selected solvent ITP. Search the extracted tree if the archive layout differs. If the box is ever missing, build one with deterministic `gmx` commands (insert W beads, minimise, equilibrate at the default temperature and pressure, check the density), record the commands with it, and report the preparation. Do not substitute atomistic `spc216.gro` as if it were Martini water.

Alternative source: https://github.com/marrink-lab/martini-forcefields . Choose one numerical baseline; do not silently replace the original archive with changing repository files.

## 3. Acquire Bible and expert molecular examples

```bash
git clone https://github.com/Martini-Force-Field-Initiative/M3-Bible.git \
  "$MARTINI_AGENT_ROOT/references/original/bible"
git clone https://github.com/Martini-Force-Field-Initiative/M3-Small-Molecules.git \
  "$MARTINI_AGENT_ROOT/references/original/small-molecules"
git -C "$MARTINI_AGENT_ROOT/references/original/bible" rev-parse HEAD > "$MARTINI_AGENT_ROOT/references/original/bible.commit.txt"
git -C "$MARTINI_AGENT_ROOT/references/original/small-molecules" rev-parse HEAD > "$MARTINI_AGENT_ROOT/references/original/small-molecules.commit.txt"
```

After acquisition, do not pull updates in an experiment snapshot. To reproduce elsewhere, check out the recorded commits. Keep repositories intact including license and README files. The helpers do not execute any repository programs.

Key paths:
- `original/bible/building_block_table.pdf` and `original/bible/src/`.
- `original/small-molecules/LIBRARY.md`.
- `original/small-molecules/models/martini_v3.0.0_small_molecules_v1.itp` and accompanying model/mapping files.

The current Bible includes later publications; verify all chosen types against the numerical baseline. For holdout testing, prepare a separate curated reference root omitting exact target models (including their entries in combined ITPs, indexes and figures), with a record of the filtering. Keep the full upstream snapshot outside the generation filesystem view. A full-repository setup is convenient development mode, not automatically a blind benchmark.

### Martini 3 PEO data file

Ether and PEG/PEO chains need Martini 3 chain terms. The pinned source is the `PEO.martini3.ff` data file from polyply_1.0 (Apache-2.0), indexed as `POLYPLY-M3:PEO`. It ships in the release bundle under `original/polymers-m3/polyply_1.0/` with its LICENSE, README, citations.bib and `polyply_1.0.commit.txt`. If building without the bundle, fetch those four files with `fetch_reference.py --kind text` from `https://raw.githubusercontent.com/marrink-lab/polyply_1.0/<commit>/...` at the commit and paths recorded in the manifest, and write the commit to `original/polymers-m3/polyply_1.0.commit.txt`. Use the data file only; do not install or run polyply (an automatic generator). Add `references/original/polymers-m3` to `generation_reference_paths`.

## 4. Obtain the worked bonded-construction example

Official tutorial: https://cgmartini.nl/docs/tutorials/Martini3/Small_Molecule_Parametrization/

```bash
python "$SKILL_DIR/scripts/fetch_reference.py" \
  https://cgmartini-library.s3.ca-central-1.amazonaws.com/0_Tutorials/m3_tutorials/Parametrize_Small_Molecule/parametrize_new_small_molecule.zip \
  "$MARTINI_AGENT_ROOT/references/archives/parametrize_new_small_molecule.zip" --kind zip \
  --extract "$MARTINI_AGENT_ROOT/references/original/tutorial"
python "$SKILL_DIR/scripts/fetch_reference.py" \
  https://cgmartini.nl/docs/tutorials/Martini3/Small_Molecule_Parametrization/ \
  "$MARTINI_AGENT_ROOT/references/original/tutorial/page.html" --kind html
```

Save a browser-rendered PDF of the full tutorial as `original/tutorial/tutorial.pdf`; verify mapping/bonded diagrams are visible. Download HTML alone does not save linked images. Alternatively save a complete webpage with its assets. Retain `ENAP-worked` and all supplied mapping, topology, geometry, radius and simulation files. Search the extracted tree rather than assuming its outer folder name. Section 5.1 is the central bonded-construction reference. This tutorial contains legacy/open-beta discussion: use its construction principles while verifying numerical files against the chosen released force field. Its illustrative initial parameters are not automatically final fitted values.

## 5. Obtain papers and SI manually when needed

Use authorized publisher downloads; do not bypass access controls. Save the actual PDFs, not login HTML. Preserve a `sources.json` alongside them with title, DOI, publisher URL, retrieval date, original downloaded filename and local filename for each item.

| Destination under original/papers/ | Source |
|---|---|
| `souza2021/article.pdf` | https://doi.org/10.1038/s41592-021-01098-3 |
| `souza2021/supplementary-information.pdf` | Same article's Supplementary Information attachment |
| `alessandri2022/article.pdf` | https://doi.org/10.1002/adts.202100391 |
| `alessandri2022/supporting-information.pdf` | Same article's Supporting Information attachment |

The polyply paper (Grünewald et al. 2022, Nat. Commun., doi 10.1038/s41467-021-27627-4) is open access. Its article ships in the bundle; fetch its SI with `fetch_reference.py --kind pdf` from the URL in the manifest (source `polyply-si`) to `original/papers/grunewald2022-polyply/supplementary-information.pdf`, and keep the folder's `sources.json`. SI Table 1 lists PEO in the Martini 3 library.

These local filenames are our convention, not asserted publisher attachment names. Retain additional supporting data/code archives under the corresponding paper directory if supplied. Read attachment titles to distinguish SI from reporting summaries. If a separate SI is unavailable, record it as missing; do not substitute an unrelated document or mark preparation complete.

Create paper subdirectories before copying downloads. Inspect title/DOI and open every PDF. Use `pdfinfo` to confirm parseability. For searchable copies:

```bash
pdftotext -layout "$MARTINI_AGENT_ROOT/references/original/papers/souza2021/article.pdf" \
  "$MARTINI_AGENT_ROOT/references/derived/souza2021-article.txt"
pdftotext -layout "$MARTINI_AGENT_ROOT/references/original/papers/souza2021/supplementary-information.pdf" \
  "$MARTINI_AGENT_ROOT/references/derived/souza2021-si.txt"
pdftotext -layout "$MARTINI_AGENT_ROOT/references/original/bible/building_block_table.pdf" \
  "$MARTINI_AGENT_ROOT/references/derived/bible.txt"
```

Repeat for the second paper/SI and tutorial. Text does not preserve chemical diagrams: keep PDFs available for rendering at generation time. Index original SI sections A1-A4 and C1-C3, tutorial Section 5.1, and the small-molecule mapping figures. Do not assume extracted page numbers equal printed numbering; record both where needed.

## 6. Acquire execution references

```bash
python "$SKILL_DIR/scripts/fetch_reference.py" \
  https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/example_input_files/mdps/martini_v3.0_prod.mdp \
  "$MARTINI_AGENT_ROOT/references/original/mdp/martini_v3.0_prod.mdp" --kind text
```

Download GROMACS manual pages matching `gmx --version`. Example for 2026.0:

```bash
GMX_DOC_VERSION=2026.0
python "$SKILL_DIR/scripts/fetch_reference.py" \
  "https://manual.gromacs.org/documentation/$GMX_DOC_VERSION/reference-manual/topologies/topology-file-formats.html" \
  "$MARTINI_AGENT_ROOT/references/original/gromacs/topology-file-formats.html" --kind html
python "$SKILL_DIR/scripts/fetch_reference.py" \
  "https://manual.gromacs.org/documentation/$GMX_DOC_VERSION/reference-manual/functions/bonded-interactions.html" \
  "$MARTINI_AGENT_ROOT/references/original/gromacs/bonded-interactions.html" --kind html
python "$SKILL_DIR/scripts/fetch_reference.py" \
  "https://manual.gromacs.org/documentation/$GMX_DOC_VERSION/reference-manual/functions/interaction-methods.html" \
  "$MARTINI_AGENT_ROOT/references/original/gromacs/interaction-methods.html" --kind html
```

HTML may require linked assets for equations/figures, so also fetch the complete PDF manual for the same version and derive its text:

```bash
python "$SKILL_DIR/scripts/fetch_reference.py" \
  "https://ftp.gromacs.org/manual/manual-$GMX_DOC_VERSION.pdf" \
  "$MARTINI_AGENT_ROOT/references/original/gromacs/manual-$GMX_DOC_VERSION.pdf" --kind pdf
pdftotext -layout "$MARTINI_AGENT_ROOT/references/original/gromacs/manual-$GMX_DOC_VERSION.pdf" \
  "$MARTINI_AGENT_ROOT/references/derived/gromacs-manual-$GMX_DOC_VERSION.txt"
```

Confirm equations and topology parameter tables are readable. `sections.json` points to the topology-directive table (Table 5.14 in the 2026.3 manual); check the table number and pages in the site's version and record any difference in `readiness.json`.

### Adapt the production MDP

The official MDP is written for membranes at 310 K and gives grompp warnings with current GROMACS. Keep it unchanged under `original/mdp/` and write an adapted copy to `derived/mdp/martini_v3.0_prod_single-solute_gmx<VERSION>.mdp`, marking every changed line `; ADAPTED` and listing the changes in a header:

1. Remove `ns_type` (obsolete).
2. `tc-grps = System` (a single small solute is too small for its own coupling group); `ref_t` and `gen_temp` 298.15.
3. `Pcoupltype = isotropic` (no membrane), with single `ref_p` (1.0) and `compressibility` (3e-4) values.
4. `nstpcouple = 5`, so that C-rescale's `tau_p >= 25*nstpcouple*dt` holds.
5. `comm-grps = System`.

`nsteps` is set per run. The release bundle includes the copy adapted for GROMACS 2026.3. On any version, check it on a small solute-in-water system: `grompp -maxwarn 0` must pass (notes are acceptable; record them), and a 10 ps run must complete. If the version differs and changes are needed, write a new adapted file named for that version rather than editing the bundled one. Set `production_mdp` in environment.json to the file that passed. Record the check in `derived/README.md` and `readiness.json`. A successful HTTP response may still be an error page: inspect each title/content. Do not use atomistic bonded parameter files as Martini parameter defaults.

## 7. Configure, index and lock

Edit root `environment.json`. Null FF paths are intentional placeholders, not guessed filenames. Set core, solvent, ions and solvent box after inspection; set actual GROMACS command/version and `python_command`. Keep benchmark mode off and generators/AA fitting disabled for the first experiment. JSON does not expand environment variables or '~'; use absolute or project-relative paths.

Build the index with the skill's script, using `python_command`:

```bash
python "$SKILL_DIR/scripts/build_index.py" build "$MARTINI_AGENT_ROOT/references"
```

It writes `index/models.jsonl` (data-contract.md) from the small-molecule library, the force-field solvent and ion ITPs and the Bible spreadsheet, and `index/sections.json`. Reviewed curation decisions (corrected SMILES, names, agent-assigned solvent and ion SMILES, motif SMARTS, section locators, hand-written entries) are in `assets/index_curation.json`; change that file, not the script, and record why. The script re-checks each SMILES against the AA topology formula or the ITP charge sum and lists failures in its summary. Review every listed failure, null SMILES and missing file. `build_index.py check` rebuilds in memory and diffs against the index on disk; with the same reference bytes and RDKit version as the authors' snapshot it reports `identical`. Spot-check entries against the original mapping and topology. Create `references/index/readiness.json` listing every required resource with pass/missing/content-review-needed and evidence. Do not confuse file integrity with readiness. Run `python -m json.tool environment.json` from the project root to validate JSON syntax; then inspect all referenced paths.

Once original and derived references/indexes are reviewed:

```bash
python "$SKILL_DIR/scripts/reference_lock.py" create \
  "$MARTINI_AGENT_ROOT/references" "$MARTINI_AGENT_ROOT/config/references.lock.json"
python "$SKILL_DIR/scripts/reference_lock.py" verify \
  "$MARTINI_AGENT_ROOT/references" "$MARTINI_AGENT_ROOT/config/references.lock.json"
```

The lock is outside references to avoid self-hashing. It rejects missing, new and modified files and deliberately refuses to overwrite an existing lock. Excludes .git metadata and __pycache__, rejects symlinks. If references need updating, create a new snapshot and lock, not a silent relock of an old experiment. The current environment JSON and lock hashes are copied into each run manifest. Keep run-specific derived work under runs, so it does not mutate locked references.

## 8. Invoke

Prepare: `Use $build-martini-agent in prepare mode for /home/me/martini-agent-lab/environment.json. Inspect existing files, acquire missing public references, and list any papers/SI requiring my download.`

Generate: `Use $build-martini-agent with environment /home/me/martini-agent-lab/environment.json and SMILES CCOCCO. Generate independently, using no automatic CG generator, and run the available execution checks.`

Compare: `Use $build-martini-agent in compare mode for the frozen run at PATH. Compare with my published topology at BENCHMARK_PATH; preserve the independent result.`

The skill is an agent workflow, not a standalone force-field prediction executable. The bundled scripts acquire and lock references only. The agent must implement/preserve deterministic per-run molecular geometry and file-writing code and execute real GROMACS tests. If RDKit/GROMACS or references are unavailable, report the corresponding limitation; do not claim a tested model.
