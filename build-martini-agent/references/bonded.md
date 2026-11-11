# Bonded terms without MD

`scripts/derive_bonded.py` turns a SMILES and a mapping into Martini 3 bonded terms. It makes no agent decisions and runs no dynamics. The rules are in `scripts/d5bonded/terms.py`.

## Ensemble

1. Up to 300 ETKDGv3 conformers, minimized with MMFF94 (`d5bonded/ensemble.py`). This is only the conformer search.
2. Each unique minimum is re-minimized under the energy model, OPLS-AA + GBSA-OBC2 (OBC II, mbondi2 radii, ACE nonpolar term), computed in NumPy from the molecule's OPLS-AA GROMACS topology with analytic gradients (`d5bonded/energy.py`). Minimization is L-BFGS-B preconditioned by the bonded terms' Gauss-Newton Hessian, after 30 plain iterations (`runs/minimizer_benchmark-20261107/report.md`). It replaced an OpenMM implementation of the same model on 2026-11-07; the two agree to round-off on 426 engine topologies (`runs/openmm_parity-20261107/report.md`), and 20 re-derived models are compared in `runs/backend_regression_v2-20261107/report.md`. The topology is supplied by the user (`cgbuild.py input --aa-topology`) and used as given, whatever OPLS engine and charge method made it; it is flattened into `input/aa_topology.top`, and a bare engine `.itp` without `[ defaults ]` gets the OPLS-AA defaults (1 3 yes 0.5 0.5). Without a supplied topology, the engine in `opls_ff_gen_command` is run (`d5bonded/opls.py`). The topology's bond graph, protonation state included, must equal the SMILES graph, and its net charge must equal the SMILES formal charge within 0.02 e. Atomic charges above 3.4 e are rejected; above 2 e the molecule is flagged. fudgeLJ/fudgeQQ other than 0.5/0.5 are flagged (`supplied_topology_notes`). The parser needs comb-rule 3 and gen-pairs, so non-OPLS force fields are rejected.
   The benchmarks below used OPLS/2020 topologies from opls-ff-gen: 1.14×CM1A-LBCC charges for neutral species, CM1A on an MMFF94 geometry for charged species. That choice belongs to the engine fallback only: the full tool's vacuum xTB step moved NH3+ protons, so the fallback runs it with `--optimize mmff --charge-model cm1a` for charged species. A supplied topology is used with whatever charges it has; other charge methods are accepted but were not benchmarked. Supplying an engine's own topology gives the same terms as running the engine here (GLYC, CYTS, LPHE, 2026-11-10).
   Polar hydrogens (on N, O, S) have σ = ε = 0 in OPLS-AA. In the energy model only, they get σ = 0.1 nm and ε = 0.192 kJ/mol, in the pair LJ and in the 1–4 LJ (fudgeLJ 0.5); engine files are not changed. Without this, the hydrogen of P–OH or N(H)–P falls onto a phosphate oxygen three bonds away. Evidence and cost: `runs/polar_h_benchmark-20261107/report.md`.
   Minima with non-bonded atoms closer than 1.0 Å are still discarded (`collapsed_minima_discarded`). If the engine or the energy model fails, stage 4 fails.
3. Hessians come from finite differences of the gradient. Minima more than 8 kcal/mol above the lowest get no weight.
4. The weight of minimum i is W_i ∝ exp(−E_i/kT) · sqrt(det I_i) · ∏_k (2πkT/λ_ik)^½. This is the harmonic superposition approximation, with λ the internal Hessian eigenvalues.
5. 200 thermal samples per minimum: exact rotations about rotatable bonds, linear displacements for stiff modes. Achiral molecules also get mirror copies.
6. The mapping weights, including hydrogens, turn samples into bead centres of geometry.

Frames come from different minima in arbitrary orientations, so positions must never be averaged across frames. Use distances, or one representative frame.

An OPLS engine may run MOPAC or xTB internally to prepare a topology. They produce force-field parameters only, so the `allow_aa_md_or_qm` setting, which concerns AA MD and QM sampling, does not apply.

## Which terms exist

| Structure | Rule | Published pattern |
| --- | --- | --- |
| Aromatic or all-sp2 ring | Rigid group; rings sharing a bead merge | |
| Rigid group, 3 beads | Constraint triangle | BENZ, TOLU |
| Rigid group, 4 beads | Rhombus: 5 constraints, 1 improper, 1 excluded diagonal. Shorter diagonal; on a tie, the first in bead order | OMA, EBEN |
| Rigid group, more real beads | Not allowed: the model is `blocked`. Keep four real hinge beads and make the rest virtual sites (SI C4) | NAPH, ANTH, CAFF, TECE |
| Virtual sites | Left out of the constraint network; all group pairs excluded | NAPH |
| Bead attached to a rigid group | Bonds to two adjacent group beads, chosen by mean distance; one position term; excluded from other group beads | PBEN, CUME, CLPR |
| Other connected beads | Bond, never a constraint | |
| Angle i-j-k | Kept unless i-k is connected, it lies in a rigid group or its substituent triangle, or k is below kT | |
| Dihedral | Largest-range candidate per central edge. Dropped if flat (range below kT), or if a flanking angle is below 20° or above 160° in over 5% of the ensemble, or has θ0 above 160° (SI C4) | |
| 2nd or 3rd neighbours | Excluded if their mean distance is below their Lennard-Jones σ (SI C3 p. 44) | |

## Values

| Term | Form | Value |
| --- | --- | --- |
| Bond | funct 1 | b0 = mean; k = kT/variance, capped by bead type: T 20,000, S 30,000, R 40,000 kJ mol⁻¹ nm⁻². Mixed pairs use the harmonic mean (T-S 24,000). Floor 700 (SI C3 p. 44) |
| Virtual site | `virtual_sites3` funct 1 | Coefficients a and b by least squares on the ensemble, per frame |
| Constraint in a rigid group | | Mean length, regardless of stiffness |
| Angle | funct 2 (cosine-harmonic) | Maximum-likelihood fit of a Gaussian in cos θ, truncated to [−1, 1] |
| Improper or narrow dihedral | funct 2 | Circular mean; k = kT/variance |
| Other dihedral | funct 9, n = 1–4 | Least squares on −kT ln P. Unpopulated regions are held 6 kT above the minimum, with low weight |

**Timestep limit:** each mode keeps a period of at least 9 steps at 20 fs. This is about the T-bead bond cap, and it caps angles and impropers from their masses and geometry. A capped bond is softer than fitted, which is accepted. The fitted value is kept in `fitted_before_limit`.

## Evidence and limits

Evaluated against 20 ns AA references on 50 neutral molecules of up to 25 heavy atoms (flexible median W1, bond/angle/dihedral): `opls-ff-gen` 0.0026 nm / 3.1° / 6.9°. The references were run with `opls-ff-gen` topologies (xTB geometry for the CM1A charges of these neutral molecules), and with OpenMM as the energy model, which the NumPy model reproduces. `opls-ff-lite` was not re-evaluated against them; for neutral molecules its MMFF94 geometry changes CM1A charges by a median 0.03 e (largest 0.2 e) across 36 acceptance molecules.

- Distances are COG values with no SASA or volume refinement, so saturated two-bead rings come out much shorter than published `opt` models (`two_bead_ring_cog_length`).
- Angle force constants were softer than OPLS MD on flexible molecules, on a small sample.
- Charged molecules have no AA reference here: they are checked numerically, not physically. Phosphate conformer weights are less certain.
- Molecules above about 25 heavy atoms are not validated.
- Stability can depend on the random seed. One passing run is weak evidence for a stiff, fused molecule, so test several seeds (`scripts/multi_seed_check.py`).
- Bonded terms come from the atomistic ensemble, and CG non-bonded interactions can shift them. Compare smoke-test and ensemble values.
