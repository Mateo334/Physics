# Historical record: kicked-top observational entropy and classical exponents

Archived during workspace reorganization, 2026-11-16. Original detailed content follows; later updates supersede earlier status statements. Read the [current summary](../calculations/kicked-top-simple-povms-lyapunov.md) first.

Historical paths below: `calculations/*.py` now means `calculations/code/*.py`; all generated files and output prefixes now live in `calculations/results/`. Filenames are unchanged. The random-PVM runner requires the saved simple data/summary AND the prior k=4..7 random-PVM data regression fixture; this record is not a fresh-fixture regeneration recipe.

LaTeX presentation update: at user request, `latex/kicked_top_random_projective_povms.tex` was greatly shortened to only the POVM construction, ranks, brief motivation/source, and averaging convention. Removed dynamics, fit discussion and unrelated caveats; retained standalone preamble. Other LaTeX files unchanged. Compiler availability checked again: neither pdflatex nor tectonic available, so NOT compile-tested.

## Random-projective-measurement average

Classical-exponent marker update [Numerically tested]: `plot_results` now uses t_E=ln(400)/mean(lambda_cl), where the mean is over the same 20 initial sphere points, with 1000 burn-in and 10000 measured kicks. It reuses k=4 data and newly computes k=8,12,16,20 using the existing tangent-map routine. Mean exponents: 0.661727, 1.246870, 1.620066, 1.839376, 2.057074 per kick. Corresponding markers: 9.05429, 4.80520, 3.69828, 3.25734, 2.91262 kicks. Half-run means: 0.659934, 1.247687, 1.622624, 1.837028, 2.057868; no full phase-space convergence claim. k=0 rotation check, sphere-norm monitoring, and positive/finite rate/time checks passed. Plot visually inspected. No quantum evolution recomputed.

New outputs with prefix `calculations/kicked_top_random_pvms_k4_20`: `_ehrenfest_classical.png`, `_ehrenfest_classical.pdf`, `_classical_markers.npz` (all running per-orbit estimates and initial points), and `_classical_markers.json`. Earlier paper-estimate figures and fit CSVs remain unchanged; old null CSV entries above k=7 do not incorporate these newer marker computations. The marker is based on the inverse mean exponent, not an average of inverse exponents. The ln(d) convention is retained as requested; these remain semiclassical reference estimates, not measured break times. This update supersedes earlier statements that no exponents beyond k=7 had been computed, only for the four newly evaluated k values.

Ehrenfest annotation: `plot_results` in `kicked_top_random_pvms.py` replots saved extended data with vertical dashed markers and a labeled top axis. Uses the reference estimate t_E=ln(d)/ln(k/2) stated in `papers/2212.01585v3.pdf`, p. 7, not the measured classical exponents or a saturation fit. Representative markers k=4,8,12,16,20 are at 8.64386, 4.32193, 3.34390, 2.88129, 2.60206 kicks. This formula is a convention-dependent semiclassical estimate, not a measured exact break time; particularly at k=4 the underlying large-k Lyapunov approximation should not be treated as exact. New files: `calculations/kicked_top_random_pvms_k4_20_ehrenfest.png` and `.pdf`. Marker positivity/monotonicity checks passed; plot visually inspected. No new evolution computed and no plateau line added.

Extension [Numerically tested]: `kicked_top_random_pvms.py` now evaluates every integer k=4,...,20, preserving all earlier output files. New prefix: `calculations/kicked_top_random_pvms_k4_20`, with PNG/PDF, data NPZ, mean CSV, fit CSV and summary JSON. Same 20 states, directions, volumes and 0–20 kicks. A colorbar identifies k on the single graph. Overlap k=4,5,6,7 reproduces prior means within absolute tolerance 1e-10; directions match at 1e-14. Script now additionally requires the saved earlier random-PVM NPZ for this regression check. All other checks passed; generated PNG visually inspected.

At high k the curves reach near-plateau values within very few kicks, so the fixed 1–4 slope includes substantial saturation curvature: at k=20, slope=0.26154 nats/kick, R2=0.82211. This is not evidence of a smaller classical Lyapunov exponent. No new classical exponents were computed for k>7; their CSV cells are empty and JSON fields null, not extrapolated. Existing values for k=4..7 are reused.

POVM explanation: standalone `latex/kicked_top_random_projective_povms.tex`, including explicit projectors, volumes, axis sampling, entropy average, source attribution, and limitations. Neither pdflatex nor tectonic was found on PATH; this file is NOT compile-tested.

[Numerically tested] `python calculations/kicked_top_random_pvms.py` (requires the prior `kicked_top_simple_data.npz` and summary) samples 20 independent sphere-uniform spin directions by normalizing Gaussian vectors, axis seed 20261116. Each direction defines a fixed projective measurement: ordered eigenstates of n dot J, lower half grouped into rank-2 cells and upper half into rank-4 cells. These are random physical spin axes, not Haar-random bases in dimension 400. Uses exactly the prior 20 coherent states, j=199.5, k=4,5,6,7, kicks 0–20; the same axes/states are reused across k and time. Entropies are evaluated separately and then averaged over axes and states, not on averaged effects or a mixed state. Quantum dynamics are recomputed; classical estimates reused unchanged.

Outputs share prefix `calculations/kicked_top_random_pvms`: one-graph PNG/PDF, `_means.csv`, `_fits.csv`, `_data.npz`, `_summary.json`. The NPZ preserves all axis/state/time entropies, all directions and state centers. The single plot has four k curves; no uncertainty bands. Slope fits and window sensitivity use the shared `affine_fit` function. Primary kicks 1–4 slopes (nats/kick) are 0.41442, 0.47028, 0.45356, 0.46355, with R2 0.98541, 0.97986, 0.98366, 0.92386. These do not increase monotonically with k. Two conditional SEMs are saved: state sampling with axes held fixed, and axis sampling with states held fixed. They must not be interpreted as independent errors of 400 iid samples or as total confidence intervals. Fits remain short-time proxies, not established linear-regime rates.

Passed eigenbasis orthonormality (absolute tolerance 1e-11), spin-component spectrum (1e-10), disjoint complete partition volumes, maximally mixed entropy, nonnegative/normalized outcome probabilities and entropy bounds (1e-10), common initial values across k, and consistency of the double average (1e-12). Same float64/complex128 environment as below. Plot visually inspected. A half-versus-full-axis mean-difference diagnostic is stored, but no multi-seed convergence study was performed.

## Short-time affine-fit update

[Numerically tested] `kicked_top_simple.py` now includes `FIT_WINDOW=(1,4)`, `affine_fit`, and `growth_fits`, called by main. These fit S(n)=intercept+slope*n to the ensemble mean. Fits were run on the saved NPZ without redoing dynamics; the integrated main was not rerun for this update. Files: `calculations/kicked_top_simple_fits.csv`, `.png`, `.pdf`. Plot visually inspected. Tests recover known affine slopes/intercepts, R2, and across-state slope SEM at absolute tolerance 1e-12, plus zero slope for constant trajectories.

Window 1–4 is a shared short-time proxy, not an automatically detected linear regime or a best-R2 selection. Also report windows 1–3, 2–4, 1–5. R2 is computed for the ensemble mean; uncertainty is sample standard deviation of individual-state OLS slopes divided by sqrt(20), not regression error from treating time points as independent. It excludes window bias and finite-j effects.

| k | Paper Jz | Axis A | Axis B | Unsharp A | Classical mean |
|---|---:|---:|---:|---:|---:|
| 4 | 0.5153 | 0.4050 | 0.4570 | 0.2278 | 0.6617 |
| 5 | 0.6452 | 0.4795 | 0.3904 | 0.2663 | 0.8759 |
| 6 | 0.6527 | 0.5062 | 0.2877 | 0.2811 | 0.9618 |
| 7 | 0.7339 | 0.4864 | 0.3623 | 0.2676 | 1.0997 |

Entropy slopes in nats/kick. Paper-PVM R2 values: 0.9797, 0.9770, 0.9924, 0.9531; state SEMs: 0.0532, 0.0350, 0.0313, 0.0278. At k=7 its slopes are 0.9322 (1–3), 0.7339 (1–4), 0.5870 (2–4), 0.5665 (1–5), so curvature/window sensitivity is substantial despite high R2. Other POVMs need not show monotone growth of slope with k. No stable linear-regime rate or Pesin equality established. Earlier statements below that no fits were obtained describe the original run and are superseded by this update; no scaling or asymptotic conclusion has changed.

[Numerically tested] Self-contained script: `calculations/kicked_top_simple.py`; run `python calculations/kicked_top_simple.py`. Dependencies: NumPy, SciPy, Matplotlib. Constants at the top control parameters. Outputs: `calculations/kicked_top_simple.png`, `.pdf`, `_data.npz` (individual and mean entropies, SEMs, running per-trajectory classical exponents, initial angles and measurement axes), `_summary.json`. Existing simulations were not overwritten.

## Dynamics, states and POVMs

j=199.5, d=400, k=4,5,6,7, rotation pi/2; U=exp(-i k Jz^2/(2j))exp(-i pi Jy/2), rotation first. Twenty pure spin-coherent states with iid sphere-uniform centers, seed 20261115, same states for every k and POVM, 0–20 quantum kicks. Entropies are averaged after evaluation, not on the mixed state. Units: nats and kicks. No backaction. Compared with previous plots, d and seed have changed.

1. Paper PVM: consecutive Jz cells, rank 2 over the lower-m half of the Hilbert space, rank 4 over the upper-m half. This implements the rank allocation of `papers/2212.01585v3.pdf`, pp. 5–6. The source does not specify in the checked prose which half carries each rank; our choice is explicit. Initial-state sampling/count is not an exact replication of the source.
2. Spin component along random unit axis A; eigenstates sorted by eigenvalue, grouped into rank-4 projectors.
3. Same construction for independently chosen axis B.
4. Unsharp version of A: E_a=0.75 P_a+0.25 Tr(P_a)I/d.

Random axes are normalized independent Gaussian vectors, hence sphere-uniform and fixed across all k and times. A=(0.860055,-0.045449,0.508174), B=(-0.439166,0.871870,0.216738). These are physically natural spin-component measurements, not Haar-random Hilbert-space projectors. The rank-4 bins compare orientation at fixed resolution; A versus unsharp A compares noise. Source rationale for the paper PVM is recorded in `../papers/observational-entropy-kicked-top.md`.

## Classical exponents

Uses the sphere map in source Eq. (11), p. 3: X=z cos(kx)+y sin(kx), Y=-z sin(kx)+y cos(kx), Z=-x. Propagates an exact-Jacobian tangent vector, projects out radial numerical error, normalizes each step and accumulates log stretching. Uses the same 20 initial sphere points as the quantum coherent states. Burns in 1000 steps for alignment, then measures 10000 steps. No stochastic noise or finite-separation trajectory approximation. Classical trajectories are evaluated in double precision; individual orbits are not claimed accurate indefinitely.

| k | Mean exponent per kick, 10000 measured steps | Mean after first 5000 measured steps | Between-state SEM, 10000 steps |
|---|---:|---:|---:|
| 4 | 0.661727 | 0.659934 | 0.001363 |
| 5 | 0.875858 | 0.874593 | 0.001415 |
| 6 | 0.961774 | 0.976184 | 0.009377 |
| 7 | 1.099673 | 1.096823 | 0.002943 |

These are finite-time ensemble estimates, not proven asymptotic exponents or rigorous error bars. k=6 has a noticeable half/full-run shift (~0.0144). A 20-point ensemble can miss regular islands; the SEM omits finite-time bias and inadequate phase-space sampling. No claim of monotonic exponents for every orbit.

For the smooth area-preserving sphere map, the classical Pesin comparison is the invariant-area average of positive Lyapunov exponents, not a single selected orbit. Its applicability does not imply equality with a fixed-POVM quantum entropy slope. No OE linear-window fit, quantum Pesin theorem or j-scaling result has been obtained.

## Checks

Verified coherent-state formula against a matrix-exponential state (absolute tolerance 1e-11), rotation/basis unitarity (1e-11), volumes and complete disjoint partitions, maximally mixed entropy, binary entropy, underflow-safe entropy evaluation, probability normalization and entropy bounds (1e-10), equal initial entropy across k, and entropy increase under added noise. POVM positivity/completeness follow from orthonormal eigenbases and disjoint positive-rank bins, with convex white-noise mixing. Maximum evolved norm error 3.81e-13. Classical Jacobian checked against central differences at h=1e-6 (absolute tolerance 2e-8); sphere tangency, area preservation, and k=0 zero-exponent rotation checked at 1e-12. Sphere norm monitored at 1e-10 throughout evolution. Plot visually inspected.

Environment: Python 3.11.9; NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0. float64/complex128. No new dependencies. Plot shows means without uncertainty bands; SEMs are saved. No classical phase-space convergence sweep or entropy-slope fitting performed.
