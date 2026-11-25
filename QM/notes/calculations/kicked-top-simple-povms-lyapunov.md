# Kicked-top entropy and classical exponents: current summary

## Run and find results
- [Unresolved] Reproducibility blocker (2026-11-25): both scripts below are absent from this checkout. Commands are historical, not runnable here; recover verified originals before rerunning. Recorded scientific results below were not reassessed by the workflow audit.
- From the workspace: `python calculations/code/kicked_top_simple.py`, then `python calculations/code/kicked_top_random_pvms.py`.
- Dependencies: NumPy, SciPy, Matplotlib. Outputs: `calculations/results/`.
- Random-PVM execution requires `kicked_top_simple_data.npz`, `kicked_top_simple_summary.json` and the saved prior `kicked_top_random_pvms_data.npz` regression fixture. It does not independently regenerate that fixture from scratch.
- Simple outputs: `kicked_top_simple*`; extended random-PVM outputs: `kicked_top_random_pvms_k4_20*`. NPZ holds trajectories/sampling data; CSV holds means/fits; JSON records summaries; PNG/PDF are figures.

## Model and completed computations [Numerically tested]
- j=199.5, d=400; rotation π/2 followed by the Jz² twist; 20 sphere-uniform coherent states, seed 20261115; 0–20 kicks. Final-only fixed measurements, no backaction. Average individual entropies, not the entropy of a mixed state; units nats/kick.
- Simple run: k=4,5,6,7; paper-style Jz rank-2/rank-4 partition, two random spin-axis rank-4 PVMs, and an unsharp version of axis A with noise weight 0.25.
- Random-PVM run: 20 fixed sphere-uniform spin axes, seed 20261116, each with rank-2/rank-4 bins; average over axes and the same states. These are physical spin axes, not Haar-random Hilbert-space bases. Extension covers every integer k=4,…,20 and reproduces the saved k=4..7 means within 1e-10.
- Short-time affine fits use kicks 1–4, with windows 1–3, 2–4, 1–5 for sensitivity. Paper-PVM slopes at k=4..7: 0.5153, 0.6452, 0.6527, 0.7339. Random-PVM slopes: 0.41442, 0.47028, 0.45356, 0.46355. Curvature is substantial; high-k slopes include saturation.
- Classical exact-Jacobian tangent evolution: same 20 sphere points, 1000 burn-in and 10000 measured kicks. Mean exponents at k=4..7: 0.661727, 0.875858, 0.961774, 1.099673. Additional k=8,12,16,20 means: 1.246870, 1.620066, 1.839376, 2.057074.
- `_ehrenfest_classical` figures use ln(400)/mean classical exponent; `_ehrenfest` figures use the paper estimate ln(d)/ln(k/2). Both are reference estimates, not measured break times. Older fit CSV nulls above k=7 predate the separate `_classical_markers` files.

## Interpretation and checks
No stable linear-regime quantum rate, quantum Pesin equality, j-scaling or full phase-space convergence established. State/axis SEMs are conditional sampling diagnostics, not 400 independent samples or total confidence intervals. Classical estimates are finite-time, with noticeable half/full-run shift at k=6.

Recorded checks cover coherent states, unitarity, measurement partitions, probabilities/entropy bounds, affine-fit toys, classical Jacobian and sphere-map checks, averaging and regression agreement. Figures were visually inspected. The shortened POVM LaTeX was not compile-tested.

Read [the detailed historical record](../archive/kicked-top-simple-povms-lyapunov-history.md) for tables, tolerances, source locations, environment and update chronology; [source note](../papers/observational-entropy-kicked-top.md) for the paper conventions.
