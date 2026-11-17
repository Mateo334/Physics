# Calculations

- `code/`: Python scripts; shared imports remain alongside their callers.
- `code/four_model_setup.py` holds the four-model observational-entropy preparations; `code/chernoff_utils.py` evaluates classical outcome Chernoff entropy.
- `results/`: generated PNG/PDF figures, NPZ data, JSON summaries and CSV tables. Existing filenames are retained.
- [Calculation notes](../notes/INDEX.md#numerical-calculations): current summaries, interpretation and historical details.

## Entry points
Run from the workspace root with NumPy, SciPy and Matplotlib installed:

- `py -3.13 calculations/code/compare_obsen_four_models_100_states.py` — 100-state observational entropy for SYK, kicked tops, Bose--Hubbard and GOE.
- `py -3.13 calculations/code/kicked_top_chernoff_observational_povm_ranks.py` — j=30, 100-state Chernoff observational entropy for several fixed PVMs and a full-rank unsharp Jz POVM. Saves `results/kicked_top_chernoff_observational_povm_ranks_j30_100_states.png` and `.npz`; the plot does not display a difference statistic.
- `py -3.13 calculations/code/kicked_top_fixed_povm_backaction_comparison.py` — fixed-POVM outcome KL with/without discarded-outcome backaction.
- `py -3.13 calculations/code/kicked_top_husimi_parameter_animation.py` — five-kick Husimi MP4 animation; requires ffmpeg. Options: `--resolution 540p --k-start 3 --k-stop 6 --fps 30`.
- `python calculations/code/kicked_top_relative_entropy.py` — fixed-Jz final-only outcome KL for 200 kicks.

Scripts contain numerical checks; saved plots and finite-time fits do not establish universal chaos or quantum Pesin relations. Read the linked notes before interpreting slopes.
