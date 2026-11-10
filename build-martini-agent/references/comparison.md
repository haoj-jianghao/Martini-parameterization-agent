# Compare mode

Compare a frozen generated model with a reference model of the same molecule (for example an Auto-Martini topology or a published one). The reference is a comparator, not ground truth: report differences, never "errors", and never revise the frozen model from the comparison without marking the revision.

## 1. Preconditions and integrity

1. The generated run must have a frozen `baseline/` with hashes. Verify them before reading any reference file, and again after the comparison. A mismatch stops the comparison.
2. Record each reference file's path, SHA-256, producing tool, version and arguments (from its own README or log, when supplied) in `comparison/inputs.json`. Auto-Martini or another CG generator may be used only here, and only when explicitly requested.
3. Do not run MD on reference topologies: no grompp, minimization or dynamics. All reference-side quantities come from the topology text, the reference coordinate file, and the generated run's atomistic structure.
4. Blindness: a strong blind test uses a fresh generation session whose filesystem holds only curated approved references, without benchmark data or network access. An instruction not to read a folder is not access isolation. Record pretraining contamination as unknown.

## 2. Parse both topologies

- Parse the `[ moleculetype ]` of each ITP with `FLEXIBLE` undefined (the production variant): `#ifdef/#ifndef/#else/#endif` are honoured, so constraint blocks written as `#ifndef FLEXIBLE [constraints] #endif` are read as constraints, and extra columns after a constraint length are ignored.
- Collect beads (index, type, charge, mass), bonds, constraints, angles, dihedrals (proper and improper by function number), virtual sites and exclusions.
- Split each bead type into size (R regular, S small, T tiny), chemical class (P, N, C, X, Q, D, W), polarity level (the digit) and labels (a, d, r, h, v, e, q, p, n), using the pinned core ITP type list to reject unknown names.

## 3. Map beads to atoms and match beads

1. Generated beads: atom members come from `decisions/mapping.json` (one-based atom_ids of `input/atomistic.sdf`).
2. Reference beads: if the reference records a bead-to-atom correspondence (Auto-Martini writes `; atoms: C0, O1, C2,` after each bead: element plus zero-based heavy-atom index), parse it. Establish the reference atom order by parsing the SMILES the reference used (its header, e.g. `; SMILES code :`) with RDKit, then align that molecule to the generated structure with an RDKit substructure match on heavy atoms. Check that every recorded element matches the aligned atom; a mismatch stops the comparison for that molecule. When symmetry gives several matches, keep the one with the highest total bead overlap (below) and record that the choice was symmetric.
3. If the reference has no bead-to-atom record, mark atom-level matching `not_available`, skip §3 matching and bead-level bonded comparison, and report only bead counts, type composition and whole-molecule observables.
4. Overlap between a generated bead g and a reference bead r = Jaccard index of their heavy-atom sets (weighting shared atoms by their allocation weights).
5. Match beads one-to-one by maximising total overlap (Hungarian assignment over the overlap matrix). A pair with overlap ≥ 0.5 is `matched`; a pair below 0.5 is `partial`; beads left without a partner are `unmatched` and listed with the atoms they cover. Also report, for each unmatched bead, the bead on the other side with which it shares the most atoms.

## 4. Metrics

For each molecule, compute:

a. **Bead count**: generated, reference, difference; heavy atoms per bead on each side.

b. **Bead vs bead** (matched and partial pairs): overlap; atoms of each bead; size difference (R/S/T, plus the change in step: T→S = +1, T→R = +2); chemical class same or different; polarity-level difference within a class (signed, generated minus reference); labels added or removed; charge difference. Summaries: fraction of matched pairs with identical type, identical size, identical class.

c. **Bonded terms**: align terms through matched bead pairs (a reference term on beads r_i…r_k corresponds to the generated term on their matched beads g_i…g_k). For each aligned term report: kind on each side (bond, constraint, angle function, dihedral function), equilibrium value difference (nm or degrees), force-constant difference where both have one, and whether the functional forms agree. List terms present only in the generated model and only in the reference. Compare parameters only for terms whose beads are all `matched`; for other terms report presence only. Summaries: counts per section on each side, aligned-term coverage, median absolute difference of bond lengths and of angle values.

d. **Whole-molecule observables**, computed without MD, on two bases:
   1. Common basis (primary): place both models' beads at the centres of geometry of their atoms (hydrogens included as in the generated mapping) on the same atomistic structures: the generated run's selected conformer, and each of its regenerated conformers. For the conformer set, compute Rg and SASA on every conformer and report mean and standard deviation; never average coordinates across conformers. Regenerate conformers with the run's own recorded procedure and seed, and check that the atom order and the selected conformer reproduce `input/atomistic.sdf`.
   2. Own coordinates (secondary): each model's own coordinate file (the generated `model/molecule.gro` and the reference `.gro`).
   For each basis report the radius of gyration (unweighted and mass-weighted, nm) and SASA with probe 0.191 nm, ndots 4800 and bead radii R 0.264, S 0.230, T 0.191 nm (Souza SI section C3, p. 42). Virtual-site beads are included in SASA with their bead radius; with zero mass they drop out of mass-weighted Rg.
   SASA method (GROMACS 2026.3 has no radii flag for `gmx sasa`):
   - Run `gmx sasa` on a plain `.gro` file, with no topology and no `grompp`. This also avoids `-maxwarn` and virtual-site mass errors.
   - Put a `vdwradii.dat` in the working directory. GROMACS matches it by residue and atom name, so use a short fixed residue name (≤ 5 characters, e.g. `MOL`) and unique atom names (`B1`…`Bn`), written identically in both files.
   - GROMACS prints "radii will be guessed" even when the file is used, and silently falls back to guessed radii when names do not match. Therefore, before any comparison, run a self-test: single isolated R, S and T beads must give 4π(r + 0.191)² within 0.5%. Abort if not.
   - Cross-check every `gmx sasa` value against an independent Shrake–Rupley calculation with the same radii and probe, and fail any value that differs by more than 2%.
   Optionally add the atomistic SASA of the same conformer for context, stating the atom radii used.

When the two mappings differ greatly, the whole-molecule observables are the meaningful comparison; do not over-interpret element-by-element parameter differences.

## 5. Outputs

Per molecule, in the generated run's `comparison/` folder:
- `inputs.json`: reference files with SHA-256 and producing tool details; baseline hash verification before and after.
- `atom_correspondence.json`: the reference-to-generated atom alignment and any symmetric alternatives.
- `bead_matching.json`: overlap matrix, assignment, matched/partial/unmatched lists.
- `comparison.json`: all metrics of §4 in the schema of `data-contract.md`.
- `report.md`: generated from `comparison.json`; a bead-matching table, a bonded-term table and the observable table, followed by stated caveats.

Per batch, in `runs/<batch>/comparison/`: `summary.json` and `report.md` tabulating, per molecule, bead counts, matched fraction, identical-type fraction, aligned bonded coverage, Rg and SASA on both bases, and any molecule where matching was not possible.

## 6. Afterwards

Save any post-comparison revision of a model as a new revision with `comparison_assisted=true`; never overwrite the baseline. Record runtime, manual interventions, failures and repeated-run variability. Do not add the reference models to the generation references.
