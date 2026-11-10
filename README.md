# Martini parameterization agent

An agent skill (`build-martini-agent`) for agentic coding tools such as Claude Code. It turns a SMILES string into a Martini 3 GROMACS model (mapping, bead types, charges, bonded terms, topology, coordinates) and checks that the model runs in water.

The agent makes the chemical decisions: mapping, bead types, charges. Deterministic code in the skill handles geometry, bonded terms, file writing and GROMACS. Bonded terms come from an MD-free method: OPLS + GBSA-OBC2 conformer ensembles and Hessians, with topologies from `opls-ff-gen` (`build-martini-agent/references/bonded.md`). The only dynamics are coarse-grained checks of the finished model.

Modes: **prepare** (a pinned, locked reference library, once per site), **generate** (one model per SMILES), **compare** (a frozen result against a published model).

> Private pre-release for test installs. Licences of the bundled references and of the skill are not settled yet.

## What is in this release

| Item | Where |
| --- | --- |
| The skill | `build-martini-agent/` (copy the whole folder; `SKILL.md` alone is not self-contained) |
| Reference manifest | `reference_manifest.json`: every reference source, how to acquire it, licence status, bundled or not, SHA-256 of every file |
| Locks | `locks/skill.lock.json` (skill v12), `locks/references.lock.json` (reference snapshot v7, our site) |
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
| `opls-ff-gen` with MOPAC (and xTB) | OPLS topologies for bonded terms; its own environment is fine. Install from the OPLS-ff-gen repository with its `environment.yml`; `pip install --no-deps` avoids OpenMM, which the workflow does not use. **Version tested: commit `0b925c5` (0.3.0 plus later changes), not yet released; a pinned release is pending.** |
| GROMACS | Tested with 2026.3. A GPU is optional |
| Poppler (`pdftotext`, `pdfinfo`, `pdftoppm`) | PDF search and rendering |
| An agentic coding tool with a shell | Claude Code is tested. Codex is likely to work but untested |

## Install

```bash
export MARTINI_AGENT_ROOT="$HOME/martini-agent-lab"
export SKILL_DIR="$HOME/.claude/skills/build-martini-agent"

# 1. Skill
mkdir -p "$(dirname "$SKILL_DIR")" && cp -r build-martini-agent "$SKILL_DIR"

# 2. References bundle (download the release asset first)
sha256sum -c SHA256SUMS
mkdir -p "$MARTINI_AGENT_ROOT/config" "$MARTINI_AGENT_ROOT/runs"
tar -xzf martini-agent-references-v7.tar.gz -C "$MARTINI_AGENT_ROOT"

# 3. Configuration
cp "$SKILL_DIR/assets/environment.example.json" "$MARTINI_AGENT_ROOT/environment.json"
```

Then follow `build-martini-agent/references/environment.md`, section 0:

1. In `environment.json`, set the force-field ITP paths, `python_command` (absolute path of the env's `python`), `gromacs_command`, `opls_ff_gen_command` and `gromacs_documentation_version`. The template already points to the water box, the adapted MDP and the PEO data.
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

6. Ask the agent to review the setup in prepare mode, then generate, for example:

   `Use $build-martini-agent with environment $HOME/martini-agent-lab/environment.json and SMILES CCOCCO. Generate independently, using no automatic CG generator, and run the available execution checks.`

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
