# Checks and bounded repair

## Static

Check input graph; mapping coverage/connectivity and shared weights; intended charge; valid bead names; finite geometry; units; atom order; indices; functional-form arguments; duplicates; missing terms; mass redistribution; virtual-site geometry; redundant constraints; ring planarity; exclusions; motif junctions. `cgbuild.py build` records these in `test/static_checks.json`; the OPLS topology's bond graph must equal the SMILES graph, or stage 4 stops. Inspect the annotated mapping. A custom linter does not replace GROMACS preprocessing.

## Execution

Use approved compatible Martini solvent topology and coordinates, never atomistic water. Inspect include order/residue names. Prepare a periodic solvated box with adequate separation and counterions where needed. Preserve commands/counts. If an approved solvent box is absent, acquire it during preparation or construct/equilibrate one explicitly and report the preparation. Generation uses fixed references offline.

Use tutorial EM/equilibration inputs and official production MDP as starting points. Inspect compatibility with installed GROMACS, chosen model and ensemble. Record adaptations. Do not fabricate universal settings; use a valid thermostat/barostat for each stage and avoid inappropriate tiny temperature-coupling groups.

Run in sequence, retaining all command output:
1. `gmx grompp -f em.mdp -c solvated.gro -p system.top -o em.tpr -pp processed.top`
2. `gmx mdrun -deffnm em`
3. Preprocess/run equilibration using actual minimized coordinates.
4. Preprocess/run a default 100 ps smoke test, targeting 0.02 ps timestep after smaller equilibration steps where needed.

This duration is only a numerical smoke test. Compare each term's smoke-test mean and width with its `ensemble` values in `test/execution.json`. A large shift means CG non-bonded interactions are competing with the bonded term, which is a limitation of the model to report, not a reason to change its values by hand. Record actual timesteps, durations, completion times, warnings and geometry checks. Reduced-timestep success is not success at the requested timestep. Inspect finite energies/coordinates, LINCS/constraint failures, extreme deformation and ring/virtual-site geometry. Record minimization convergence criterion and actual result; step-limit termination is not convergence. Investigate warnings individually instead of blanket `-maxwarn`.

A stiff or fused model can pass or fail depending on the velocity seed, so one pass is weak evidence for these. Run `scripts/multi_seed_check.py` and record the number of seeds passed in the report. The check covers EM and both equilibration stages, not the smoke test.

Preserve failures and permit at most the configured repair count. Explain changes and rerun affected checks. Do not freeze flexible terms or add arbitrary restraints just to obtain a stable trajectory.

## Status

`blocked`: unresolved chemistry/parameter/dependency or failure preventing the requested output, including any failed static check, such as a rigid group above four real beads or a virtual site over 0.05 nm from its mapped position. `constructed_unexecuted`: candidate and static checks only. `preprocessed`: GROMACS preprocessing succeeded, dynamics not fully checked; explain failures. `execution_checked`: documented minimization, equilibration and smoke test succeeded and geometry/logs inspected. Failure to reach `emtol` is reported but does not block this status. State outstanding warnings and provisional parameters at every level.

## Optional physical evaluation

`scripts/sasa_check.py` compares the mapped model (bead radii R 0.264, S 0.230, T 0.191 nm) with the D5 atomistic ensemble (Rowland-Taylor radii from the tutorial's `vdwradii_AA.dat`) on the same frames, with probe 0.191 nm and 4800 dots. It first checks that `gmx sasa` used the given radii, and cross-checks values against an independent calculation within 2%. Bands: pass within 5%, review to 10%, flag above 10%. It measures mapping and bead size only: no CG dynamics enters, and the target comes from D5, not MD. Results are report-only, in `sasa.json`.

Declare reference conditions, targets and tolerances before comparison. Compare bonded distributions, conformational populations, SASA/shape, density, transfer free energies or application observables as available; report uncertainty and sampling. Minimized-conformer geometry/SASA is not ensemble validation. Use the tutorial's explicit AA/CG radii and probe convention for SASA comparisons. Do not treat RDKit conformer frequency as Boltzmann population.
