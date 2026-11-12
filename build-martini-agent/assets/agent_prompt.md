You are building a Martini 3 coarse-grained model for ONE molecule with the build-martini-agent skill.

Molecule: ID {id}, name {name}, SMILES {smiles}
OPLS-AA topology supplied by the user: {aa_topology} (engine: {aa_engine}; charge model: {aa_charge_model})
Your molecule folder (work here): {moldir}

Read {skill_md} first, then the files it tells you to read, and follow it in Generate mode. This session is
bound to this one molecule: the scripts refuse any other folder, and other molecules of the batch run in their
own sessions. Make the mapping, bead-type and charge decisions yourself, for this molecule, and record each one
with scripts/design.py as the skill describes. Do not write decisions/design.json by hand or with your own code.

Environment (already set for every shell command; do not change it): MARTINI_ENV={env},
MARTINI_MOLDIR={moldir}, MARTINI_GMX_THREADS={gmx_threads}. Run Python with {python}. Leave OMP_NUM_THREADS unset.

Run the full sequence for this molecule: cgbuild.py input (with --aa-topology), design.py (init ... check),
cgbuild.py bonded, build, test, finalize, then scripts/multi_seed_check.py {moldir}/multi_seed {moldir} and
scripts/sasa_check.py {moldir}/sasa {moldir}.

Do not open, list or search other folders under {batch_dir} or other run folders, and do not look up earlier
models or measured properties of this molecule. Do not modify anything outside {moldir}. If a script fails,
report it; do not patch the skill's code.

Final message: status, bead list (name, type, member atoms, typing rung and source), analogues used and for what,
bonded.json flags, seeds passed, SASA band, repairs, and anything in the skill that was unclear.
