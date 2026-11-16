# Calculations

- `code/`: Python scripts; shared imports remain alongside their callers.
- `results/`: generated PNG/PDF figures, NPZ data, JSON summaries and CSV tables. Existing filenames are retained.
- [Calculation notes](../notes/INDEX.md#numerical-calculations): current summaries, interpretation and historical details.

## Main commands
Run from the workspace root with NumPy, SciPy and Matplotlib installed:

```powershell
python calculations/code/kicked_top_simple.py
python calculations/code/kicked_top_random_pvms.py
```

`kicked_top_simple.py` computes the four-measurement, 20-state kicked-top comparison and classical exponents. Results use prefix `results/kicked_top_simple`.

`kicked_top_random_pvms.py` computes the 20-random-spin-axis average for k=4,…,20. It requires these existing files in `results/`:
- `kicked_top_simple_data.npz`
- `kicked_top_simple_summary.json`
- `kicked_top_random_pvms_data.npz` — saved prior k=4..7 regression fixture.

The simple command supplies its own data and summary. The random-PVM command checks against the prior fixture; it does **not** independently regenerate that fixture from scratch. Extended outputs use prefix `results/kicked_top_random_pvms_k4_20`; older four-k results remain available.

## Other entry points
- `python calculations/code/kicked_top_relative_entropy.py` — minimal direct evolution of reference and perturbed states; fixed-Jz outcome relative-entropy plot, 200 kicks.
- `python calculations/code/kicked_top_frechet.py` — original derivative check, 30 kicks.
- `python calculations/code/kicked_top_growth.py` — 200-kick extension and exponential fits.
- `python calculations/code/kicked_top_observational_entropy.py` — j=200 single-state entropy.
- `python calculations/code/kicked_top_entropy_state_average.py` — its 20-state average.
- `python calculations/code/plot_kicked_top_four_k.py` — replot saved single-state CSV.
- `python calculations/code/syk_observational_entropy.py` — N=16 SYK fixed-POVM entropy.

Scripts contain numerical checks; saved plots and finite-time fits do not establish universal chaos or quantum Pesin relations. Read the linked notes before interpreting slopes.
