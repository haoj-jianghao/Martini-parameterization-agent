# Martini parameterization agent

An agent skill (`build-martini-agent`) for agentic coding tools such as Claude Code. It turns a SMILES string into a Martini 3 GROMACS model (mapping, bead types, charges, bonded terms, topology, coordinates) and checks that the model runs in water. It is distributed with a reference manifest and a bundle of the openly licensed references.

The agent makes the chemical decisions: mapping, bead types, charges. Deterministic code in the skill handles geometry, bonded terms, file writing and GROMACS. Bonded terms come from an MD-free method: OPLS + GBSA-OBC2 conformer ensembles and Hessians, evaluated on an all-atom OPLS-AA topology of the molecule that you supply (`build-martini-agent/references/bonded.md`). Make it with any OPLS engine (for example opls-ff-gen, opls-ff-lite or the LigParGen server), with any charge method; no OPLS engine is needed at the site. The only dynamics are coarse-grained checks of the finished model.

Modes: **prepare** (a pinned, locked reference library, once per site), **generate** (one model per SMILES), **compare** (a frozen result against a published model).


## What is in this release

| Item | Where |
| --- | --- |
| The skill | `build-martini-agent/` (copy the whole folder; `SKILL.md` alone is not self-contained) |
| Reference manifest | `reference_manifest.json`: every reference source, how to acquire it, licence status, bundled or not, SHA-256 of every file |
| Locks | `locks/skill.lock.json` (skill v14), `locks/references.lock.json` (reference snapshot v7, our site) |
| Maintainer tool | `tools/make_manifest.py` (build the manifest, assemble the bundle) |
| Python environment | `environment.yml` |
| References bundle | Release asset `martini-agent-references-v7.tar.gz` with `SHA256SUMS`, not in git |

The references bundle holds the openly licensed sources:
- the Martini 3.0.0 force field;
- the small-molecule parametrization tutorial, including the water box used for the execution checks;
- the official production MDP and our copy adapted for GROMACS 2026.3;
- the Martini 3 Bible and the small-molecule library at pinned commits;
- the polyply Martini 3 PEO data file and the polyply article;
- the reference index.

Each site fetches or downloads the rest itself: the GROMACS manual for its version, the polyply SI, and the Souza 2021 and Alessandri 2022 papers with their SI (publisher copyright).

## Requirements

| Need | Notes |
| --- | --- |
| Linux or WSL | On Windows, keep the workspace in the WSL filesystem |
| Python env with RDKit, NumPy, SciPy, openpyxl | `conda env create -f environment.yml`. The reference index was built with RDKit 2025.09.6; other versions can change conformer searches and the index check |
| An OPLS-AA topology of each molecule | GROMACS format, made by you with any OPLS engine; see "Generating a topology", step 2. Not needed to install |
| GROMACS | Tested with 2026.3. A GPU is optional |
| Poppler (`pdftotext`, `pdfinfo`, `pdftoppm`) | PDF search and rendering |
| An agentic coding tool with a shell | Claude Code is tested. Codex is likely to work but untested |

## Install

### With the setup script (recommended)

Clone this repository and download the release assets `martini-agent-references-v7.tar.gz` and `SHA256SUMS` into one folder. Install GROMACS yourself. Then, from the repository folder:

```bash
python3 tools/setup_site.py \
    --bundle /path/to/martini-agent-references-v7.tar.gz \
    --gmx /path/to/gmx                 # optional if gmx is on PATH
```

The script needs only Python 3.8+ and conda (or mamba), and installs into `~/martini-agent-lab` unless you pass `--root`. It creates the conda env `martini-agent` from `environment.yml` (or uses `--python`), copies the skill into the workspace (`<root>/build-martini-agent/`) and checks it against the skill lock, writes `AGENTS.md` and `CLAUDE.md` in the workspace so that a coding agent started there finds the skill, unpacks and checks the references bundle, finds GROMACS, writes `environment.json`, downloads the GROMACS manual for your version and the polyply SI, re-checks the adapted MDP (`grompp -maxwarn 0` and a 10 ps run on the bundled water box), runs the index check and compares your references with the manifest. Steps already done are skipped, so re-run it after you supply what it reports as missing:

- the Souza 2021 and Alessandri 2022 articles and SI, placed under `references/original/papers/Souza2021/` and `references/original/papers/alessandri2022/` as `article.pdf` and `supporting-information.pdf` (the script derives their text);
- the prepare-mode review ("Generating a topology", step 4), which writes the papers' `sources.json` and `references/index/readiness.json`.

Everything goes under `--root`; nothing is written to an agent's own configuration (such as `~/.claude`), so the workspace works with any coding agent that has a shell:

```text
<root>/
  AGENTS.md, CLAUDE.md     point the agent to build-martini-agent/SKILL.md
  build-martini-agent/     the skill
  environment.json         site configuration (the skill's scripts find it here)
  references/  config/     references and their lock
  runs/                    one folder per model
```

Then run it once more with `--lock` to create the site's reference lock. It refuses while anything is failed or missing. `python3 tools/setup_site.py --help` lists every option.

### By hand

The same steps, one by one:

```bash
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
export SKILL_DIR="$MARTINI_AGENT_ROOT/build-martini-agent"

# 1. Python environment
conda env create -f environment.yml && conda activate martini-agent

# 2. Skill, inside the workspace, and the notes that point agents to it
mkdir -p "$MARTINI_AGENT_ROOT" && cp -r build-martini-agent "$SKILL_DIR"
cp tools/workspace-AGENTS.md "$MARTINI_AGENT_ROOT/AGENTS.md"
echo "@AGENTS.md" > "$MARTINI_AGENT_ROOT/CLAUDE.md"

# 3. References bundle (download the release asset first)
sha256sum -c SHA256SUMS
mkdir -p "$MARTINI_AGENT_ROOT/config" "$MARTINI_AGENT_ROOT/runs"
tar -xzf martini-agent-references-v7.tar.gz -C "$MARTINI_AGENT_ROOT"

# 4. Configuration
cp "$SKILL_DIR/assets/environment.example.json" "$MARTINI_AGENT_ROOT/environment.json"
```

Then follow `build-martini-agent/references/environment.md`, section 0:

1. In `environment.json`, set the force-field ITP paths, `python_command` (absolute path of the env's `python`), `gromacs_command` and `gromacs_documentation_version`. Leave `opls_ff_gen_command` as `null`: topologies are supplied per molecule. (A site with its own OPLS/2020 engine may set it as a fallback for molecules given without a topology; `references/environment.md`.) The template already points to the water box, the adapted MDP and the PEO data.
2. Fetch the GROMACS manual pages and PDF for your version and the polyply SI. Download the two Martini papers and SI from the publishers. Derive text with `pdftotext -layout` (environment.md sections 5 and 6).
3. If GROMACS is not 2026.3, re-check the adapted MDP: `grompp -maxwarn 0` and a 10 ps run (section 6).
4. Run `python "$SKILL_DIR/scripts/build_index.py" check "$MARTINI_AGENT_ROOT/references"`. It reports `identical` with the same bytes and RDKit version. Then write `references/index/readiness.json` (section 7).
5. Compare your references with ours, then lock them:

   ```bash
   python "$SKILL_DIR/scripts/reference_manifest.py" verify reference_manifest.json "$MARTINI_AGENT_ROOT/references"
   python "$SKILL_DIR/scripts/reference_lock.py" create "$MARTINI_AGENT_ROOT/references" "$MARTINI_AGENT_ROOT/config/references.lock.json"
   python "$SKILL_DIR/scripts/reference_lock.py" verify "$MARTINI_AGENT_ROOT/references" "$MARTINI_AGENT_ROOT/config/references.lock.json"
   ```

   Bundled sources must be `identical`. Site downloads may differ; record which ones in `readiness.json`.

6. Ask the agent to review the setup once in prepare mode ("Generating a topology", step 4), then generate models.

## Generating a topology

The agent turns one SMILES into a tested Martini 3 topology. It makes the chemical decisions itself (mapping, bead types, charges, virtual sites) and runs the skill's scripts for everything else. A run takes a few minutes per molecule. For now, only claude is tested, but codex is likely to work as well. 

### 1. Start a session

Start the agent in the workspace. Every shell it uses needs the Python environment:

```bash
conda activate martini-agent
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
# export MARTINI_ENV=/other/environment.json               # optional: default is environment.json beside the skill folder
export MARTINI_GMX_THREADS=4                                 # optional: mdrun threads (default 4)
unset OMP_NUM_THREADS                                        # GROMACS rejects a mismatch with -ntomp
cd "$MARTINI_AGENT_ROOT"
claude        # or codex, or another coding agent with a shell
```

Starting in the workspace keeps the skill, references and `runs/` inside the agent's working folder. The agent learns about the skill from `AGENTS.md` (Codex and most other agents) or `CLAUDE.md` (Claude Code), which tell it to read `build-martini-agent/SKILL.md`; naming the skill in the prompt also works. Claude Code users who want the `/build-martini-agent` command can link the skill: `ln -s "$MARTINI_AGENT_ROOT/build-martini-agent" ~/.claude/skills/`. The agent runs many shell commands (Python, GROMACS); approve them as they come, or allow shell and file edits for this workspace.

### 2. Make the all-atom topology

Stage 4 derives the CG bonded terms from an OPLS-AA model of the molecule, which you make beforehand with any OPLS engine: for example opls-ff-gen or opls-ff-lite (`--smiles ... --output-dir ...` writes `mol.top` and `mol.itp`), or the LigParGen web server (`MOL.gmx.itp`). Any charge method is accepted. The workflow uses the file as given:

- A `.top` (its `#include` files next to it) or a self-contained `.itp` with `[ atomtypes ]`. A bare `.itp` without `[ defaults ]` gets the OPLS-AA defaults (`1 3 yes 0.5 0.5`).
- It must be the same species as the SMILES, protonation state included, with the hydrogens explicit; the atom order may differ (atoms are matched by bond graph). The net charge must equal the SMILES formal charge within 0.02 e.
- comb-rule 3 and gen-pairs are required, so non-OPLS force fields (GAFF, CHARMM) are rejected.

`cgbuild.py input --aa-topology FILE` checks the file, rejects a mismatch with a clear message, and copies it into the run (`input/aa_topology.top`, with its SHA-256 in `input/request.json`). The bonded terms were benchmarked with OPLS/2020 topologies (1.14×CM1A-LBCC for neutral molecules, CM1A for ions); other charge methods work but were not benchmarked.

### 3. Ask for a model

Give the molecule's SMILES, the path of its topology, and an ID and a name if you want them in the file names:

```text
Use the build-martini-agent skill in Generate mode with environment $MARTINI_ENV.
Molecule: ID DEET, name N,N-diethyl-meta-toluamide, SMILES CCN(CC)C(=O)c1cccc(C)c1,
OPLS-AA topology /path/to/DEET/mol.top (made with opls-ff-gen, 1.14*CM1A-LBCC charges).
Make the mapping, bead-type and charge decisions yourself, use no automatic CG generator,
run the full sequence (input, design, search, bonded, build, test, finalize), then the
multi-seed check and the SASA check, and summarise the result.
```

Write the SMILES in the charge and protonation state you want simulated; the agent keeps it as given and asks only when the species is ambiguous. For several molecules, list them all in one request, each with its topology; each gets its own decisions and folder.

### 4. What the agent does

| Step | Who | What happens | Main output |
| --- | --- | --- | --- |
| Preflight | agent | Checks the Python env, GROMACS, the reference lock and that each molecule has a topology | |
| `cgbuild.py input` | script | Parses the SMILES, freezes an atom-indexed structure and conformers; checks and copies the OPLS-AA topology | `input/` |
| Mapping, types, charges | **agent** | Searches the reference index for analogues, reads the papers, Bible and analogous models, proposes up to three mappings, picks one and types each bead, citing a pinned file for every decision | `decisions/design.json` |
| `cgbuild.py search` | script | Re-runs the agent's index queries; records every hit and whether the target itself is in the references | `decisions/analogue_search.json` |
| `cgbuild.py bonded` | script | Derives all bonds, angles, dihedrals and constraints without MD: the supplied OPLS-AA topology, conformer minima, Hessians, Boltzmann-weighted samples | `decisions/bonded.json` |
| `cgbuild.py build` | script | Writes the topology, coordinates, drawing and static checks | `model/` |
| `cgbuild.py test` | script | Solvates in Martini water, adds counterions, runs minimisation, two equilibration stages and a 100 ps smoke test | `test/` |
| `cgbuild.py finalize` | script | Writes the validation record, manifest and report; freezes `baseline/` | `report.md` |
| `multi_seed_check.py` | script | Repeats minimisation and equilibration with four other velocity seeds (11, 22, 33, 44) | `multi_seed/` |
| `sasa_check.py` | script | Compares the CG SASA with the atomistic ensemble (report only) | `sasa/` |

The agent never writes bonded terms; if a check fails, it repairs the mapping or virtual sites (at most three cycles) and records each repair.

To prepare a new site, ask instead: `Use the build-martini-agent skill in prepare mode for $MARTINI_ENV. Inspect the existing files and list anything missing.`

### 5. Run without interaction (scripts, batches)

The benchmark ran one headless Claude Code session per molecule, from inside an empty molecule folder (other agents have similar non-interactive modes, e.g. `codex exec`):

```bash
mkdir -p "$MARTINI_AGENT_ROOT/runs/DEET" && cd "$MARTINI_AGENT_ROOT/runs/DEET"
claude -p "$PROMPT" --model claude-haiku-5-5 --output-format json \
  --allowedTools "Bash Read Write Edit Glob Grep" --disallowedTools "Agent Task WebFetch WebSearch" \
  --add-dir "$MARTINI_AGENT_ROOT" --max-budget-usd 10 > DEET.agent.json
```

`$PROMPT` is the request from step 3, plus: the path of `SKILL.md` to read first, the topology path, the molecule folder to work in, the `python_command` to use, and an instruction not to modify anything outside that folder and to report script failures instead of patching code. Each finished folder is self-contained; `python "$SKILL_DIR/scripts/cgbuild.py" summary <batch folder>` writes a batch report.

### 6. The result

Each molecule gets a folder under `runs/` (`runs/<name>-<timestamp>/`, or the folder you gave):

| File | What it is |
| --- | --- |
| `model/molecule.itp` | **The Martini 3 topology** (moleculetype, atoms, bonds, constraints, angles, dihedrals, virtual sites, exclusions) |
| `model/molecule.gro` | CG coordinates of one molecule |
| `model/mapping.svg` | Drawing of the mapping over the molecule |
| `test/system.top` | A working system topology: includes the core, solvent and ion ITPs, then `molecule.itp` |
| `test/smoke.gro` | The solvated system at the end of the smoke test |
| `report.md` | Status, beads with types and rationale, analogues used, bonded terms, flags, checks |
| `decisions/` | The agent's design record, the bonded record and the parameter records with their sources |
| `validation.json`, `manifest.json` | Every check with its value; software, hashes and seeds needed to reproduce the run |

To use the model in your own system, include it after the force field, as `test/system.top` does:

```text
#include "martini_v3.0.0.itp"
#include "martini_v3.0.0_solvents_v1.itp"
#include "martini_v3.0.0_ions_v1.itp"
#include "molecule.itp"
```

### 7. Check before you use it

Open `report.md` first.

- **Status:** `execution_checked` (minimisation, equilibration and smoke test passed), `preprocessed` (`grompp` passed, dynamics not fully checked), `constructed_unexecuted` (files and static checks only) or `blocked` (the report says why). Only `execution_checked` models have run.
- **Seeds:** how many of the four seeds passed. Treat a model that passes only some seeds with care.
- **Provisional types:** beads typed without a direct precedent are marked; check them against the cited rule.
- **Bonded flags:** read the `flags` and `dropped` lists in `decisions/bonded.json`, for example `conformer_search_not_saturated` or `large_atomic_charge`.
- **Reference exposure:** `True` means the target, or a tautomer of it, already exists in the references, so the run is not an independent prediction.
- **SASA band:** pass (within 5% of the atomistic value), review (5–10%) or flag (above 10%). It is reported only and changes nothing.
- **Unverified statements:** anything the agent could not cite to a pinned file is labelled "agent background knowledge, unverified" and listed at the end.

`execution_checked` means the model is numerically stable in water for these short runs. It is not a physical validation; partitioning, density or other observables need their own tests.

### Troubleshooting

| Symptom | Cause and fix |
| --- | --- |
| Every GROMACS test fails at once | `OMP_NUM_THREADS` is set; unset it |
| Scripts cannot find `environment.json` | The skill folder is not beside `environment.json` (`<root>/build-martini-agent/`); move it there, or export `MARTINI_ENV` in the agent's shell |
| `cgbuild.py input` rejects the topology | The message says why: different atoms or bonds than the SMILES (often a different protonation state or missing hydrogens), a net charge that differs from the SMILES, or an unsupported section or function type. Regenerate the topology for the exact SMILES |
| Stage 4 (`bonded`) fails with "no AA topology supplied" | The molecule was started without `--aa-topology`; re-run `cgbuild.py input` with it |
| `reference_lock.py verify` fails | A reference file changed after locking; restore it or make a new lock, never relock silently |
| Status `blocked` | Read the reason in `report.md`; most are mapping problems the agent should repair (for example a rigid group above four beads) |

## Maintainers: building a release

From the development workspace, after verifying the skill and reference locks:

```bash
python3 tools/make_manifest.py build references config/references.vN.lock.json reference_manifest.json
python3 tools/make_manifest.py bundle reference_manifest.json references build-martini-agent OUT_DIR
tar --sort=name --mtime='2026-11-09 00:00Z' --owner=0 --group=0 --numeric-owner \
    -C OUT_DIR -cf - references | gzip -n > martini-agent-references-vN.tar.gz
sha256sum martini-agent-references-vN.tar.gz > SHA256SUMS
```

`tools/make_manifest.py` expects the skill folder beside `tools/`.
