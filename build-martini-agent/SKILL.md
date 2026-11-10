---
name: build-martini-agent
description: Prepare a pinned Martini 3 reference environment and generate GROMACS coarse-grained molecular models from SMILES through agent-led mapping and bead typing, deterministic MD-free bonded terms and execution checks, without automated CG generators. Optionally compare a frozen result with withheld published or Auto-Martini models.
---

# Build Martini models with an agent

## Modes and contract

Use the LLM for mapping, bead types and charges. Use deterministic code for bonded terms, graph inspection, geometry, file writing and numerical checks. The agent does not choose bonded terms or their values. Deliver a candidate model or an explicit partial result that says what is unresolved. Never claim physical validation from generated files or short dynamics alone.

- **Prepare:** follow `references/environment.md` (section 0 when starting from the release bundle; check it with `scripts/reference_manifest.py verify`); use `assets/environment.example.json`. Acquire fixed references, build the index with `scripts/build_index.py` and lock them before generation.
- **Generate:** accept SMILES and environment JSON path. After configuration, SMILES is the only required molecular input. Follow the stages below.
- **Compare:** only after freezing an independent result, follow `references/comparison.md`. It matches beads through their atoms, then reports bead counts, bead-vs-bead size/type/polarity differences, bonded-term differences, and radius of gyration and SASA without MD. A reference model (for example Auto-Martini's) is a comparator, not ground truth, and is never simulated.

Generate mode runs these commands in order, each with `python_command`, from one molecule folder `<moldir>`:

| Command | Does | Needs |
| --- | --- | --- |
| `scripts/cgbuild.py input <moldir> --smiles S --id ID --name NAME` | Freezes the indexed structure and conformers | SMILES |
| You write `decisions/design.json` | Mapping, types, charges, virtual sites (stages 2 and 3) | Reference search |
| `scripts/cgbuild.py search <moldir>` | Runs your index queries and records every hit | `design.json` |
| `scripts/cgbuild.py bonded <moldir>` | Derives all bonded terms (stage 4) | OPLS engine (`opls-ff-lite`) |
| `scripts/cgbuild.py build <moldir>` | Writes ITP, GRO, drawing, records and static checks | `bonded.json` |
| `scripts/cgbuild.py test <moldir>` | Solvates, adds counterions, runs EM, equilibration and the smoke test | GROMACS |
| `scripts/cgbuild.py finalize <moldir>` | Writes validation, manifest and report; freezes `baseline/` | |
| `scripts/multi_seed_check.py OUT <moldir>` | Repeats EM and equilibration with other seeds (stage 5) | GROMACS |
| `scripts/sasa_check.py OUT <moldir>` | Report-only SASA of the mapped model against its D5 atomistic ensemble (stage 6) | GROMACS |

The environment file is `$MARTINI_ENV`, or `environment.json` next to the skill folder. `$MARTINI_GMX_THREADS` sets the `mdrun` threads (default 4). Leave `OMP_NUM_THREADS` unset: GROMACS rejects a mismatch with `-ntomp` and every test fails. `scripts/cgbuild.py summary <batchdir>` writes a batch report.

Read `references/data-contract.md` before writing run files and `references/validation.md` before simulation. Resolve relative paths against the environment JSON directory. Do not assume the skill directory is the molecular workspace.

## Experiment boundary

Allow RDKit, Python/NumPy/SciPy, generic graph/geometry libraries, drawing, reference search, PDF rendering, GROMACS, and for stage 4 the OPLS engine named by `opls_ff_gen_command` (`opls-ff-lite`, shipped with the release; the full `opls-ff-gen` also works). Do not invoke, import, copy decision algorithms from, or use generated output from Auto-Martini, SMartini, SwarmCG, CGCompiler, Martini Mapper, martinize or other automatic CG parameter generators during generation. Generic parsers and deterministic bookkeeping are allowed.

Allow documented reference fragments and analogous published models. For strict holdout tests use a precurated reference subset; do not search a full topology file containing the target. If the target topology has been seen, mark `reference_exposed`; do not call the run blind independent generation. Do not silently clear that flag. General guidelines and bead definitions remain allowed.

Cite every factual claim, chemical decision and numerical value to a pinned file under the reference root (path plus locator: section, page, table, moleculetype or line), or to a reproducible calculation in the run. Label anything not cited that way, including recalled literature, citations, model names, experimental values and rules of thumb, exactly "agent background knowledge, unverified" wherever it appears (decisions, parameter records, comments, report). Unverified knowledge may suggest what to look for, but it cannot by itself justify a bead type, a mapping or a numerical parameter: such a term is `unresolved`. If a pinned file later confirms or contradicts the statement, cite it; correct visibly and do not delete the original statement. List every unverified statement in the report.

Parameterization is MD-free: stage 4 evaluates conformer energies, minima and Hessians with an OPLS force field and runs no dynamics. Do not run AA MD, Bartender or iterative bonded fitting by default. If explicitly enabled, mark `assisted_refinement` and preserve the baseline. The only dynamics in the workflow are the CG execution checks of the finished model. Never alter the nonbonded interaction matrix or introduce new bead types in this workflow.

## 1. Preflight and identity

1. Read environment configuration. Run all Python (helpers, RDKit, run scripts) with the interpreter in `python_command`; confirm it imports RDKit before generation and record its path and version. Confirm that it also imports NumPy and SciPy, and that `opls_ff_gen_command` runs. Run `scripts/reference_lock.py verify REFERENCE_ROOT LOCK_FILE`. Explain missing/changed references before numerical assignment; do not silently relock. The script checks integrity, not scientific completeness.
2. Confirm selected core ITP, solvent/ion dependencies, relevant paper/SI, Bible, examples, tutorial and matching GROMACS documentation exist. Inspect actual content. Use only approved paths. Record source hashes, software versions, skill hash, seeds, agent/model identifier if available and tool actions.
3. Parse/sanitize SMILES using RDKit. Preserve encoded charge, protonation, tautomer and stereochemistry; save original and canonical isomeric SMILES. Document disconnected components and unspecified stereochemistry. Ask only where ambiguity changes the intended species; never silently neutralize or remove salts.
4. Freeze an explicit-H atom-indexed structure and atom table. Use stable one-based IDs, preserving conversion to RDKit zero-based indices. Generate at most 10 conformers with seed 20261031 by default. Record geometry optimization method/convergence. Preserve SDF atom ordering and chirality. A minimized conformer is not a thermodynamic ensemble. Report failed embedding rather than inventing coordinates.

## 2. Choose mapping

1. Search the reference index (`reference_index`, `models.jsonl`) before choosing a mapping. Query it with code, not by memory: substructure matches of the target's functional groups, ring systems and junctions against each entry's canonical SMILES, plus its motif tags. Open the topology, mapping file and figure of every relevant hit. Record the queries, every hit, and which were used or rejected and why in `decisions/analogue_search.json`. A target found in the index or the Bible, as the same constitution or as another tautomer of it (RDKit canonical tautomer, checked by `cgbuild.py search`), sets `reference_exposed`.
2. Read relevant Bible entries, original SI and the analogous model mappings found above. Cite exact file/page/table/molecule section for transferred decisions. Render figures when text cannot convey chemical mapping.
3. Propose at most three plausible mappings without exhaustive global enumeration. Consider functional groups, symmetry, charge localization, branching, ring shape, resolution, overall volume and motif junctions.
4. Work in full-molecule context. Local decomposition is allowed; inspect each junction and conjugated group. Do not blindly concatenate fragment topologies.
   - **Choosing among precedents.** When several documented models cover the same fragment, start from the closest context, in this order: (a) the same ring system with the same fusion, as a published model (for example the force-field nucleobases for purines and pyrimidines, CAFF for xanthines); (b) the same fused ring pair with other substituents (BZIM, INDO, QUIN); (c) the isolated ring (IMID, PYMI, PYRI). Take its bead split and hinge layout first. Moving down the order, or departing from that precedent's split, is allowed when chemistry or geometry justifies it (a substituent, tautomer or protonation change, or a virtual site that misses the 0.05 nm limit); state the reason in the bead rationale and in the rejection reason of the closer precedent. Different defensible mappings are acceptable, as between chemists; an unexplained departure from a closer precedent is not.
   - **Junctions between a rigid group and a flexible part** (for example the glycosidic N9–C1′ of a nucleoside): keep each atom in a bead of its own part. Do not put an atom of a flexible ring or chain into a rigid-group bead. Otherwise several flexible beads connect to the rigid group, stage 4 bridges each to it, and the flexible part is over-restrained.
5. Write candidates, chosen mapping and concise selection rationale in `decisions/design.json`; its fields are in `references/data-contract.md`. Declare all atom contributions and intentional shared-atom allocation weights. Draw the mapping over the indexed molecule.
6. Every bead must contain at least one atom; `cgbuild` cannot hold a bead with no atoms (for example the central TC4 of the Bible glucose). When a documented model has one, transfer the other beads and record the omission as a deviation.
7. Check connected groups, atom coverage and intentional overlap. Allocation weights sum to one across beads for each AA atom; normalize within each bead to compute its center of geometry. Include hydrogens per the source convention. Mapping weights do not define bead masses.

## 3. Assign types and charge

Parse allowed types from the pinned core ITP. Use documented chemical/size/polarity/label precedents, including substituent effects. Verify every chosen name exists. Do not borrow a later extension silently. Assign net charge from intended chemical state and compatible documented models, not arbitrary sums of atomistic partial charges. Record a source-backed assignment for every bead.

Type each bead by this ladder and record which rung was used:
1. The same fragment typed in an indexed model or Bible entry → `transferred`.
2. The nearest documented chemical group (Souza SI Table 24/25), refined with the SI proximity rules (section C2: inductive/conjugation, local dipole, internal H-bonds, resolution) → `adapted`, naming the group and the rule applied.
3. No matching group: select the type from the partitioning data of the isolated fragment, as the SI prescribes (Souza SI section C introduction, step 2, p. 35; Supplementary Tables 17–19) → `adapted`, citing the table and values.
4. Only if all three fail: `unresolved`, with the record of what each rung found.

A type from rung 2 or 3 for chemistry with no direct precedent (for example a carbonate or a conjugated lactone) is a provisional type: it does not block execution, but it must be flagged `provisional_type` and listed in the report.

## 4. Derive bonded model

Run `scripts/cgbuild.py bonded <moldir>`. It takes the SMILES and the frozen mapping, and writes `decisions/bonded.json`. Do not write bonded terms, exclusions or `nrexcl` in `design.json`, and do not edit `bonded.json`. The method is in `references/bonded.md`.

Your remaining decisions are in the mapping:
1. **Fused rings.** Follow the SI hinge (Souza SI section C4, p. 45-46), as the published polycyclics do (NAPH, QUIN, ANTH, CAFF, TECE): four real beads forming the hinge, and the other ring beads as virtual sites. Declare each in `virtual_sites` with `construction` set to `virtual_sitesn_cog` (centre of the constructing beads) or `virtual_sites3` (an in-plane position). You choose the site and its constructing beads; stage 4 computes the coefficients. A rigid group above four real beads is `blocked`, so change the mapping. If a site is over 0.05 nm from its mapped position, change the construction or move side-chain atoms out of that bead. Avoid a near-straight rigid triangle (`near_linear_rigid_triangle`).
2. **Masses** are computed, not chosen. Each bead gets its standard size mass (R 72, S 54, T 36), each virtual site 0, and each site's standard mass is split equally over the real beads of its rigid group, keeping the total (CAFF 4 × 63, BZIM 4 × 45). A site outside any rigid group passes its mass to its constructing beads. Stage 4 and `build` use the same rule (`d5bonded/terms.py`, `bead_masses`) and record it in `bonded.json` `masses`. Do not write `mass_amu` in `design.json`; `build` ignores it and lists any difference under the mass check.

Read `bonded.json` before building: its `flags`, its `dropped` list, and any term whose `limitations` says a value was capped or floored. If stage 4 fails (the OPLS engine or the energy model), the model is `blocked`; report it and do not substitute hand-chosen values.

## Batch mode

When several SMILES are requested, run stages 2–3 as a separate decision pass for each molecule: its own index search, mapping candidates and typing ladder, each citing that molecule's own hits. Shared code must not hard-code analogue choices or bead types across molecules; any such value must be read from that molecule's decision records. Stage 4 is shared deterministic code by design.

## 5. Emit and check

`finalize` copies the deterministic generation code into the run's `scripts/`. Emit ITP/GRO, mapping and parameter records, drawing, minimal system TOP, adapted EM/equilibration/test MDPs and dependency manifest. Ensure include ordering, one defaults definition, no duplicate molecule definitions and identical topology/coordinate ordering. Bundle compatible dependencies if permitted or record exact paths/hashes.

Run `scripts/cgbuild.py build`, `test` and `finalize`. Follow `references/validation.md`. Default to one molecule in compatible Martini water; `test` neutralizes charged molecules with NA or CL. Report missing solvent or GROMACS; never imply execution. Use at most three repair cycles with retained versions and reasons. For stiff or fused molecules, also run `scripts/multi_seed_check.py`, because one passing run can be chance, and report how many seeds passed. Do not achieve stability through indiscriminate stiffening, freezing flexible groups, suppressing warnings or deleting interactions. If a model fails, report the failing term and bonded rule. A repair changes the mapping or virtual sites, or is recorded as a manual intervention.

Freeze `baseline/` with hashes before opening benchmarks. Report chemistry checks, numerical execution and physical evaluation separately. Isolate failures in batch mode. Never promote generated provisional models to the reference library automatically.

## 6. Deliver

Provide model files, test inputs, mapping, provenance, dependencies and report. State bead count, charge, analogues, provisional types, the OPLS engine used, bonded flags, seeds tested, tested timestep/duration, actual checks, deviations and reference exposure. Compute every status, count and list in the report from the run's records (mapping, parameters, analogue search, validation, manifest); never write them from a template. Reports carry no claim about a molecule that its records do not support. List all statements labelled "agent background knowledge, unverified". Report statuses `blocked`, `constructed_unexecuted`, `preprocessed` or `execution_checked`; store physical observables separately. Report the SASA band (pass within 5%, review to 10%, flag above 10%) if `sasa_check.py` was run. It changes no model or status, and do not change bond lengths to meet it without a recorded decision. If dependencies prevent completion, list concrete pending actions and deliver available work without claiming success. Do not invent calibrated confidence scores.
