# Martini parameterization agent

An agent skill (`build-martini-agent`) for agentic coding tools such as Claude Code. It turns a SMILES string into a Martini 3 GROMACS model (mapping, bead types, charges, bonded terms, topology, coordinates) and checks that the model runs in water. It is distributed together with its OPLS engine, `opls-ff-lite` 0.1.0 (`vendor/`), a reference manifest and a bundle of the openly licensed references.

The agent makes the chemical decisions: mapping, bead types, charges. Deterministic code in the skill handles geometry, bonded terms, file writing and GROMACS. Bonded terms come from an MD-free method: OPLS + GBSA-OBC2 conformer ensembles and Hessians, with topologies from `opls-ff-lite`, the reduced copy of opls-ff-gen that needs only RDKit, NumPy and MOPAC (`build-martini-agent/references/bonded.md`). The only dynamics are coarse-grained checks of the finished model.

Modes: **prepare** (a pinned, locked reference library, once per site), **generate** (one model per SMILES), **compare** (a frozen result against a published model).

> Private pre-release for test installs. Licences of the bundled references and of the skill are not settled yet.

## What is in this release

| Item | Where |
| --- | --- |
| The skill | `build-martini-agent/` (copy the whole folder; `SKILL.md` alone is not self-contained) |
| Reference manifest | `reference_manifest.json`: every reference source, how to acquire it, licence status, bundled or not, SHA-256 of every file |
| Locks | `locks/skill.lock.json` (skill v13), `locks/references.lock.json` (reference snapshot v7, our site) |
| OPLS engine | `vendor/opls_ff_lite-0.1.0-py3-none-any.whl` and its source package, with `vendor/SHA256SUMS` |
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
| Python env with RDKit, NumPy, SciPy, openpyxl, MOPAC | `conda env create -f environment.yml`. The reference index was built with RDKit 2025.09.6; other versions can change conformer searches and the index check |
| `opls-ff-lite` 0.1.0 (in `vendor/`) | OPLS/2020 topologies for bonded terms; needs only RDKit, NumPy and MOPAC, so it installs into the same environment. MIT; two AmberTools data files are GPL-3.0 (see the package's `NOTICE`). It is the reduced copy of opls-ff-gen and writes the same topologies with MMFF94 geometry |
| OPLS/2020 parameter files | Supporting Information of Jorgensen et al., *J. Phys. Chem. B* 2024 (https://doi.org/10.1021/acs.jpcb.3c06602): `jp3c06602_si_002.txt` and `jp3c06602_si_003.txt`. They may not be redistributed, so download them yourself |
| GROMACS | Tested with 2026.3. A GPU is optional |
| Poppler (`pdftotext`, `pdfinfo`, `pdftoppm`) | PDF search and rendering |
| An agentic coding tool with a shell | Claude Code is tested. Codex is likely to work but untested |

## Install

```bash
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
export SKILL_DIR="$HOME/.claude/skills/build-martini-agent"

# 1. Python environment and the OPLS engine
conda env create -f environment.yml && conda activate martini-agent
(cd vendor && sha256sum -c SHA256SUMS)
python -m pip install --no-deps vendor/opls_ff_lite-0.1.0-py3-none-any.whl
opls-ff-lite install-si /path/to/jp3c06602_si_002.txt /path/to/jp3c06602_si_003.txt
opls-ff-lite doctor          # every item must be ok

# 2. Skill
mkdir -p "$(dirname "$SKILL_DIR")" && cp -r build-martini-agent "$SKILL_DIR"

# 3. References bundle (download the release asset first)
sha256sum -c SHA256SUMS
mkdir -p "$MARTINI_AGENT_ROOT/config" "$MARTINI_AGENT_ROOT/runs"
tar -xzf martini-agent-references-v7.tar.gz -C "$MARTINI_AGENT_ROOT"

# 4. Configuration
cp "$SKILL_DIR/assets/environment.example.json" "$MARTINI_AGENT_ROOT/environment.json"
```

Then follow `build-martini-agent/references/environment.md`, section 0:

1. In `environment.json`, set the force-field ITP paths, `python_command` (absolute path of the env's `python`), `gromacs_command` and `gromacs_documentation_version`. `opls_ff_gen_command` is already `["opls-ff-lite"]`; give its absolute path if the agent's shell does not activate the environment. The template already points to the water box, the adapted MDP and the PEO data.
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

6. Ask the agent to review the setup once in prepare mode (next section, step 3), then generate models.

## Generating a topology

The agent turns one SMILES into a tested Martini 3 topology. It makes the chemical decisions itself (mapping, bead types, charges, virtual sites) and runs the skill's scripts for everything else. A run takes a few minutes per molecule: in our 26-molecule headless benchmark with Claude Haiku 5.5, the median was about 5 minutes (80–460 s) and the total cost $6.99.

### 1. Start a session

Every shell the agent uses needs the environment file and the Python environment:

```bash
conda activate martini-agent
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
export MARTINI_ENV="$MARTINI_AGENT_ROOT/environment.json"   # required: the skill's scripts do not find it otherwise
export MARTINI_GMX_THREADS=4                                 # optional: mdrun threads (default 4)
unset OMP_NUM_THREADS                                        # GROMACS rejects a mismatch with -ntomp
cd "$MARTINI_AGENT_ROOT"
claude
```

Start the agent in `$MARTINI_AGENT_ROOT`, so the references and `runs/` are inside its workspace. Claude Code finds the skill in `~/.claude/skills/build-martini-agent/`; you can name it in the prompt or invoke it as `/build-martini-agent`. The agent runs many shell commands (Python, GROMACS); approve them as they come, or allow `Bash`, `Read`, `Write` and `Edit` for this workspace.

### 2. Ask for a model

Give the molecule's SMILES, plus an ID and a name if you want them in the file names:

```text
Use the build-martini-agent skill in Generate mode with environment $MARTINI_ENV.
Molecule: ID DEET, name N,N-diethyl-meta-toluamide, SMILES CCN(CC)C(=O)c1cccc(C)c1.
Make the mapping, bead-type and charge decisions yourself, use no automatic CG generator,
run the full sequence (input, design, search, bonded, build, test, finalize), then the
multi-seed check and the SASA check, and summarise the result.
```

Write the SMILES in the charge and protonation state you want simulated; the agent keeps it as given and asks only when the species is ambiguous. For several molecules, list them all in one request; each gets its own decisions and folder.

### 3. What the agent does

| Step | Who | What happens | Main output |
| --- | --- | --- | --- |
| Preflight | agent | Checks the Python env, `opls-ff-lite`, GROMACS and the reference lock | |
| `cgbuild.py input` | script | Parses the SMILES, freezes an atom-indexed structure and conformers | `input/` |
| Mapping, types, charges | **agent** | Searches the reference index for analogues, reads the papers, Bible and analogous models, proposes up to three mappings, picks one and types each bead, citing a pinned file for every decision | `decisions/design.json` |
| `cgbuild.py search` | script | Re-runs the agent's index queries; records every hit and whether the target itself is in the references | `decisions/analogue_search.json` |
| `cgbuild.py bonded` | script | Derives all bonds, angles, dihedrals and constraints without MD: `opls-ff-lite` topology, conformer minima, Hessians, Boltzmann-weighted samples | `decisions/bonded.json` |
| `cgbuild.py build` | script | Writes the topology, coordinates, drawing and static checks | `model/` |
| `cgbuild.py test` | script | Solvates in Martini water, adds counterions, runs minimisation, two equilibration stages and a 100 ps smoke test | `test/` |
| `cgbuild.py finalize` | script | Writes the validation record, manifest and report; freezes `baseline/` | `report.md` |
| `multi_seed_check.py` | script | Repeats minimisation and equilibration with four other velocity seeds (11, 22, 33, 44) | `multi_seed/` |
| `sasa_check.py` | script | Compares the CG SASA with the atomistic ensemble (report only) | `sasa/` |

The agent never writes bonded terms; if a check fails, it repairs the mapping or virtual sites (at most three cycles) and records each repair.

To prepare a new site, ask instead: `Use the build-martini-agent skill in prepare mode for $MARTINI_ENV. Inspect the existing files and list anything missing.`

### 4. Run without interaction (scripts, batches)

The benchmark ran one headless session per molecule, from inside an empty molecule folder:

```bash
mkdir -p "$MARTINI_AGENT_ROOT/runs/DEET" && cd "$MARTINI_AGENT_ROOT/runs/DEET"
claude -p "$PROMPT" --model claude-haiku-5-5 --output-format json \
  --allowedTools "Bash Read Write Edit Glob Grep" --disallowedTools "Agent Task WebFetch WebSearch" \
  --add-dir "$MARTINI_AGENT_ROOT" --max-budget-usd 10 > DEET.agent.json
```

`$PROMPT` is the request from step 2, plus: the path of `SKILL.md` to read first, the molecule folder to work in, the `python_command` to use, and an instruction not to modify anything outside that folder and to report script failures instead of patching code. Each finished folder is self-contained; `python "$SKILL_DIR/scripts/cgbuild.py" summary <batch folder>` writes a batch report.

### 5. The result

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

### 6. Check before you use it

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
| Scripts cannot find `environment.json` | `MARTINI_ENV` is not exported in the agent's shell |
| Stage 4 (`bonded`) fails | Run `opls-ff-lite doctor`: MOPAC and the OPLS/2020 files must both be ok |
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
