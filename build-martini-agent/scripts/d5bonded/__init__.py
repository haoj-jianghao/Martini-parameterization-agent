"""Deterministic, agent-free Martini 3 bonded terms from an energy-weighted conformer ensemble.

No molecular dynamics: ETKDG conformers -> minima and Hessians under an OPLS + OBC2 energy model
-> harmonic superposition samples (D5) -> bead-mapped distributions -> Boltzmann inversion.
"""
