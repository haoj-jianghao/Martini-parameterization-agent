# Data contract

Use UTF-8 JSON/JSONL, named units and stable IDs. Paths below are prescribed harness outputs, not upstream filenames. Resolve config paths relative to the environment JSON directory.

## Reference index

Store verified entries in `references/index/models.jsonl`: `reference_id`, name, verified canonical isomeric SMILES (null if unknown), formal charge, `source_path`, SHA-256, DOI/URL, source locator (molecule section/page/figure), mapping path/figure, force-field version, motifs, notes and role (`generation` or `withheld_benchmark`). Do not invent identity from an abbreviation. Index relevant models incrementally. A topology without chemical mapping cannot establish exact motif transfer. Never index withheld answers into the generation view.

## Runs

Create `runs/<safe-name>-<timestamp>/` with:
- `input/request.json`: original SMILES, overrides, config snapshot/hash, and `aa_topology` (null, or `file`, `sha256`, `source_path`, `source_sha256`, `engine`, `charge_model`, `n_atoms`, `net_charge_e`, `notes`).
- `input/aa_topology.top`: the supplied OPLS-AA topology, flattened (`#include` expanded), when one was given.
- `input/atomistic.sdf`, `input/atoms.json`: explicit-H indexed reference structure; atoms have stable one-based ID, RDKit index, element, charge, stereochemistry and neighbors.
- `decisions/analogue_search.json`: the index queries run (substructure SMARTS or SMILES, motif tags, code path), every hit with `reference_id`, `source_path` and the matched atoms, and for each hit whether it was used, for which beads or terms, or why it was rejected. It also lists `exact_target_matches` (same constitution) and `tautomer_target_matches` (same RDKit canonical tautomer, stereo ignored); either sets `reference_exposed`.
- `decisions/design.json`: the agent's record, holding beads, members, types, charges, virtual sites, candidates, analogue queries and annotations, written only by `scripts/design.py`. It holds no bonded terms.
- `decisions/design_log.jsonl`: one line per `design.py` command (UTC time, command, arguments, agent session variables, SHA-256 of the design it left). `build` compares the last hash with `design.json` (`design_record_written_by_design_tool`).
- `decisions/reference_beads.json`: written by `design.py`; each indexed model bead by bead (type, member atoms, fragment signature), or why it has no bead view (`scripts/refbeads.py`). CGSmiles beads (Bible rows, `mapping_cgsmiles`) are placed on the source SMILES and take the signature they have there; a bead that cannot be placed unambiguously has `signature_error` and no signature, and cannot support rung 1.
- `agent_prompt.txt`: in a batch, the prompt that started this molecule's session.
- `decisions/bonded.json`: written by `cgbuild.py bonded`. Agent-free bonded record (see "Bonded record").
- `decisions/mapping_candidates.json`, `decisions/mapping.json`, `decisions/parameters.json`: written by `cgbuild.py build`.
- `bonded/`: the engine's OPLS topology and its log (`bonded/opls-ff-gen/`, only without a supplied topology), `ensemble_beads.npz` (resampled bead coordinates in nm) and `ensemble_atoms.npz` (the same frames as atoms).
- `model/molecule.itp`, `model/molecule.gro`, `model/mapping.svg` (or PNG).
- `test/`: system TOP, MDPs, preprocessed topology, construction records, logs and trajectories.
- `scripts/`: actual deterministic code used to construct outputs.
- `logs/commands.jsonl`: argument arrays, working directories, exit codes, timestamps and stdout/stderr paths, without credentials.
- `manifest.json`, `validation.json`, `report.md`.
- `baseline/`: immutable pre-comparison deliverables with hashes.
- `comparison/`: compare-mode outputs (`inputs.json`, `atom_correspondence.json`, `bead_matching.json`, `comparison.json`, `report.md`; see `comparison.md` §5).

A batch uses `runs/<batch-name>-<timestamp>/<molecule-id>/`, each molecule folder holding the full layout above, plus `batch.json` (agent, model, command template, per-molecule exit code and wall time), `logs/<ID>.out`/`.err` of each agent session, and a batch `report.md` and `summary.json` built only from the molecule records. `manifest.json` records `agent_session` (the folder the session was bound to by `MARTINI_MOLDIR`) and `design_provenance`.

## Design record

`decisions/design.json` records the agent's stage 2 and 3 decisions; `scripts/design.py` writes it, one decision per command. Top-level keys: `schema_version`, `molecule_id`, `resname`, `authored_by` (agent, model, tool, session), `analogue_queries` (each with `name`, `smarts`, `motif_tags`, and `origin`: derived from the molecule by `design.py init`, or added by the agent), `analogue_annotations` (by `reference_id`: `used`, `used_for`, `opened` (files shown by `design.py open`), `reason`), `default_rejection_reason`, `mapping_candidates` (each with `id`, `beads`, `summary`, `pros`, `cons`), `chosen_candidate`, `selection_rationale`, `beads`, `virtual_sites`, `mass_policy` (`source`), `unverified_statements`, `limitations`, `manual_interventions`, `revisions`.

Each bead has `name`, `type`, `charge_e`, `members` (heavy-atom one-based ID to allocation weight), `rationale`, `evidence_status`, `typing_ladder` (`rung_used`, `rung_results`), `provisional_type`, `sources` and `fragment_signature`. The signature is the canonical SMILES of the bead's heavy atoms with one `*` per bond leaving the bead and a mass label (100 × weight) on shared atoms; two beads hold the same fragment exactly when their signatures are equal (`scripts/refbeads.py`). A rung-1 bead cites the source bead (`REF_ID`, bead name and mapping file), whose signature equals its own. Optional `size_exception`: why a bead's size does not fit its heavy atoms, citing a published model with the same atoms in a bead of that size (SKILL.md stage 2); it clears the `bead_size_matches_heavy_atoms` flag. Do not give `mass_amu`: masses are computed (SKILL.md stage 4). Hydrogens inherit the allocation of their heavy atom; list exceptions under `hydrogen_excluded_parents` (`atom_ids`).

Each virtual site has `site` (bead number), `construction` (`virtual_sitesn_cog` or `virtual_sites3`), `constructing` (real bead numbers, three for `virtual_sites3`), `evidence_class` and `source`. Its mass is 0 and is redistributed by rule (SKILL.md stage 4). Write no bonded terms, exclusions or `nrexcl` here.

## Mapping schema

Use mapping keys `schema_version`, `atom_index_base: 1`, `atomistic_structure_sha256`, `coordinate_unit: "nm"`, `coordinate_convention`, `beads`.
Each bead has `bead_id`, name, type, `charge_e`, `mass_amu`, `members` (atom_id and positive allocation_weight), virtual-site definition or null, sources, rationale and evidence status, plus `typing_ladder` (the rung used, 1–4, and what each earlier rung found) and `provisional_type` (true when the type comes from rung 2 or 3 for chemistry with no direct precedent). Unshared atoms have allocation 1; shared atoms sum to one across beads. Normalize the allocations within each bead to compute COG. Convert Angstrom to nm explicitly. Virtual-site dynamic coordinates follow the declared construction; quantify discrepancy from mapped target geometry. Virtual sites may represent mapped atoms, but never double-count mass or charge.

## Parameter schema

`parameters.json` is generated. Entries have `term_id`, topology section, bead IDs, function number, named parameter value/unit pairs, evidence class, source, `limitations`, and the `ensemble` mean and standard deviation of the term. Bonded terms have evidence class `ensemble_derived`, with `decisions/bonded.json` as their source. Also include exclusions, nrexcl, masses and virtual sites. Constraints carry a length and no force constant; their `flexible_variant` is used only for energy minimization. Fourier dihedrals list `phi_s_n` and `k_n` for each multiplicity n. Each virtual site still needs a cited source in `design.json`.

## Bonded record

`decisions/bonded.json` keys: `schema_version`, `molecule_id`, `method`, `temperature_K`, `ensemble` (conformers, minima used, effective minima, top weight, rejected fraction, achiral, mirrored), `energy_model` (name, engine, charge model, `engines_tried` with each result, largest atomic charge, topology path), `software`, `seed`, `wall_s`, `rigid_groups`, `terms`, `exclusions`, `close_pair_exclusions` (each with neighbour order, mean distance, sigma), `nrexcl`, `dropped` (each considered term and why), `flags`. A term changed from its fitted value keeps that value in `fitted_before_limit` and the reason in `limitations`. `virtual_site_parameters` holds the computed `virtual_sites3` coefficients with their fit error. `masses` holds `bead_amu` (bead order) and `virtual_site_redistribution` (each site, its mass and the beads that received it). Flags: `rigid_group_above_four_beads` (blocks the model), `near_linear_rigid_triangle`, `two_bead_ring_cog_length`, `large_atomic_charge`, `collapsed_minima_discarded`, `conformer_search_not_saturated`.

## Manifest

Record reference lock hash, selected dependency hashes, skill hash/version, software versions, model identifier if exposed (otherwise unknown), seeds, request/config hashes, revisions, runtime and how to reproduce the test. Track `reference_exposed`, `assisted_refinement`, `comparison_assisted` and manual interventions independently. Record `unverified_statements`: a list of {text, location, later_status (unverified, confirmed by path, contradicted by path)} for every item labelled "agent background knowledge, unverified"; record an empty list when there are none. Validation checks use pass/fail/not_run with evidence paths and actual measured values. Physical evaluation lists observables/conditions/results, not a global validated boolean.

## Comparison schema

`comparison/comparison.json` keys: `schema_version`, `molecule_id`, `generated` and `reference` (paths, SHA-256, producing tool and version), `atom_matching` (`available` | `not_available`, alignment method, symmetric_choice), and:
- `bead_count`: {generated, reference, difference, heavy_atoms_per_bead_generated, heavy_atoms_per_bead_reference}.
- `bead_pairs`: list of {generated_bead, reference_bead, status (`matched` | `partial`), overlap, generated_type, reference_type, size_step_difference, class_same, polarity_difference (null across classes), labels_added, labels_removed, charge_difference_e}; plus `unmatched_generated` and `unmatched_reference` lists with atoms and best-overlapping partner; plus `bead_summary` {matched, partial, unmatched_generated, unmatched_reference, identical_type_fraction, identical_size_fraction, identical_class_fraction}.
- `bonded_terms`: list of aligned terms {section, generated_beads, reference_beads, generated_kind, reference_kind, equilibrium_difference with unit, force_constant_difference with unit or null, form_agrees}; `only_generated`; `only_reference`; `bonded_summary` {counts per section per side, aligned_coverage, median_abs_bond_length_difference_nm, median_abs_angle_difference_deg}.
- `observables`: for each basis (`common_selected_conformer`, `common_conformer_mean`, `own_coordinates`) and each side: {rg_nm, rg_mass_weighted_nm, sasa_nm2}, with `sasa_settings` {probe_nm: 0.191, ndots: 4800, bead_radii_nm: {R: 0.264, S: 0.230, T: 0.191}, source: Souza SI section C3 p. 42}; optional `atomistic_sasa_nm2` with its radii.
- `caveats`: list of strings (reference is a comparator; no MD run on the reference; mappings differ; etc.).
Values that cannot be computed are null with a reason, never omitted.

## Report

`report.md` is generated from the records above. Status, bead count, charge, analogues used, provisional types, provisional and unresolved terms, and unverified statements are each read from `validation.json`, `mapping.json`, `parameters.json`, `analogue_search.json` or `manifest.json`, never typed as fixed text. Free-text explanation is allowed but must not contradict those records.
