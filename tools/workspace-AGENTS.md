# Martini parameterization workspace

This folder is a workspace for the `build-martini-agent` skill (Martini 3 models from SMILES).

- Before any Martini parameterization task, read `build-martini-agent/SKILL.md` and follow it.
  It is the skill's entry point; its `references/` and `scripts/` sit beside it.
- One molecule per agent session. For several molecules, start `build-martini-agent/scripts/batch.py run`
  (one session per molecule); do not build them yourself in one session. Record decisions only with
  `build-martini-agent/scripts/design.py`. The reasons are in SKILL.md "Rules".
- Configuration: `environment.json` in this folder (the skill's scripts find it here; `MARTINI_ENV` may override it).
  Run Python with its `python_command`.
- References: `references/` (locked by `config/references.lock.json`). Do not edit them.
- Put each run under `runs/`. Leave `OMP_NUM_THREADS` unset when running GROMACS.
