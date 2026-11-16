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
- `py -3.13 calculations/code/kicked_top_husimi_k3_k6.py` — two spin-Husimi phase-space heatmaps in q=azimuth, p=cos(polar angle) at k=3 and k=6 after five kicks; writes `results/kicked_top_husimi_k3_t5.png` and `results/kicked_top_husimi_k6_t5.png`.
- `py -3.13 calculations/code/kicked_top_husimi_parameter_animation.py` — Matplotlib MP4 sweeping k=2 to 10 over 200 frames at 50 fps (4 s), now with a 512x361 Husimi grid and 1920x1080 video; requires ffmpeg. Older 540p movies are preserved; `--resolution 540p --k-start 3 --k-stop 6 --fps 30` reproduces the previous k=3..6 movie.
- `py -3.13 calculations/code/compare_obsen_four_models_100_states.py` — one figure comparing 100-state mean normalized observational entropy for Majorana SYK, kicked tops at k=2,4,5,6, three-site Bose--Hubbard and GOE, at t=0..5; writes `results/obsen_four_models_k2_4_5_6_5kicks_100_states.png` (older 20-step figure retained).
- `py -3.13 calculations/code/kicked_top_fixed_povm_backaction_comparison.py` — side-by-side comparison of the two fixed-POVM single-time KL protocols, with and without discarded-outcome Lüders backaction; same 100 initial pairs and same fixed rank-7 POVM, 20 kicks; saves `results/kicked_top_fixed_povm_backaction_comparison_100_states.png`.
- `py -3.13 calculations/code/coupled_kicked_tops_obsen_lyapunov.py` — two coupled kicked tops, average over 100 product coherent initial states with a fixed product phase-space PGM, versus the mean classical positive Lyapunov sum; writes `results/coupled_kicked_tops_obsen_lyapunov_pgm_100_states.png` (earlier eight-state figure preserved).
- `python calculations/code/kicked_top_relative_entropy.py` — minimal direct evolution of reference and perturbed states; fixed-Jz outcome relative-entropy plot, 200 kicks.
- `python calculations/code/kicked_top_repeated_measurement.py` — a new seeded random projective basis at each time, shared by both preparations and independent of the dynamics; saves `results/kicked_top_random_repeated_measurement.png` (marginal outcome divergence, not full-record divergence).
- `py -3.13 calculations/code/kicked_top_rank_m_povm.py` — evolve then apply the same fixed rank-7 POVM at each of five kicks; keeps all outcome branches, computes exact full-history KL for 100 state pairs, then averages. Saves `results/kicked_top_fixed_rank_m_povm_evolve_then_measure_100_states.png`; previous-order figure retained.
- `python calculations/code/kicked_top_frechet.py` — original derivative check, 30 kicks.
- `python calculations/code/kicked_top_growth.py` — 200-kick extension and exponential fits.
- `python calculations/code/kicked_top_observational_entropy.py` — j=200 single-state entropy.
- `python calculations/code/kicked_top_entropy_state_average.py` — its 20-state average.
- `python calculations/code/plot_kicked_top_four_k.py` — replot saved single-state CSV.
- `python calculations/code/syk_observational_entropy.py` — N=16 SYK fixed-POVM entropy.

Scripts contain numerical checks; saved plots and finite-time fits do not establish universal chaos or quantum Pesin relations. Read the linked notes before interpreting slopes.
