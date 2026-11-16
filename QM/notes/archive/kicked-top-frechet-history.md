# Historical record: kicked-top initial-state Frechet derivatives

Archived during workspace reorganization, 2026-11-16. Original detailed content follows, including the original 30-kick experiment and its 200-kick extension. Read the [current summary](../calculations/kicked-top-frechet.md) first.

Historical paths below: `calculations/*.py` now means `calculations/code/*.py`; generated CSV, JSON and PNG files now live in `calculations/results/`. Filenames are unchanged.

[Numerically tested] The fixed-POVM first differential and diagonal Hessian agree with central finite differences for the 62 tested (kappa,time) pairs. This is a finite-dimensional verification of the formulas, not a proof or a chaos diagnostic.

## Reproduction and outputs
- Run `python calculations/kicked_top_frechet.py` from the workspace.
- Script: `calculations/kicked_top_frechet.py`.
- Complete data: `calculations/kicked_top_frechet_results.csv`.
- Parameters, environment, convergence and selected results: `calculations/kicked_top_frechet_summary.json`.
- Figure: `calculations/kicked_top_frechet.png` (visually checked).
- Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0; double precision.

## Model and conventions
- hbar=1; spin j=10, dimension 21, Jz basis ordered m=-j,...,j.
- Floquet unitary U=exp[-i kappa Jz^2/(2j)] exp[-i (pi/2) Jy]; rotation acts first. U_t=U^t, t=0,...,30; kappa=0.5,6.0.
- Fixed projective POVM E_m=|m><m|, not measurement histories and not optimized over measurements.
- Perturbation V=exp[-i epsilon Jy], epsilon=0.07. The second state is U_t V rho V† U_t†.
- Initial state rho=0.7|psi><psi|+0.3 I/21, with |psi>=exp[-i phi Jz] exp[-i theta Jy]|j,j>, theta=1.1, phi=0.7. Minimum eigenvalue 0.0142857142857142.
- A is a seeded random Hermitian traceless matrix, normalized to operator norm 1; NumPy default_rng seed 20261112. The script specifies its construction exactly.
- Vary rho -> rho+h A in BOTH occurrences, fixing epsilon,t,U,V and E_m. These are NOT epsilon derivatives.
- Natural logarithms. First derivative sum a_i log(p_i/q_i)-(p_i/q_i)b_i; Hessian sum p_i(a_i/p_i-b_i/q_i)^2, where a_i,b_i are the two pushed-forward directional probabilities.

## Checks and convergence
Central differences evaluate the divergence independently on rho+hA and rho-hA. All tested states are checked for Hermiticity, trace one, strict positivity; probabilities are positive. Also checked spin commutators, unitary time evolution, traceless A, nonnegative entropy/Hessian, and the exact epsilon=0 null case. Matrix checks use absolute tolerances 1e-12 to 1e-11 (np.allclose also retains its default relative tolerance); trace and imaginary-part checks are absolute. At h=1e-4 derivative assertions use atol=2e-6, rtol=2e-4.

Maximum absolute finite-difference errors across both kappa values and all 31 times:

| h | First derivative | Second derivative |
|---|---:|---:|
| 1e-3 | 2.09060e-5 | 9.83687e-4 |
| 3e-4 | 1.88086e-6 | 8.84995e-5 |
| 1e-4 | 2.08978e-7 | 9.83008e-6 |
| 3e-5 | 1.88111e-8 | 5.85546e-7 |

The decrease is consistent with central-difference truncation error, with cancellation/roundoff affecting the smallest-step second differences.

Selected analytic values:

| kappa | t | Divergence | First derivative | Hessian |
|---|---:|---:|---:|---:|
| either | 0 | 0.013766319 | 0.037410075 | 1.764739604 |
| 0.5 | 10 | 0.005295633 | -0.013967440 | 3.094173842 |
| 6.0 | 10 | 0.007922623 | -0.026957714 | 1.608026978 |
| 0.5 | 30 | 0.004209923 | 0.035994582 | 3.842562854 |
| 6.0 | 30 | 0.008219380 | 0.052873019 | 4.770144850 |

## Limitations
One initial state, one direction, one measurement, one spin size, and 30 kicks only. No semiclassical scaling, ensemble statistics, optimized-measurement calculation, or Lyapunov fit. The observed responses fluctuate; these data do not establish exponential growth or a universal ordering by kappa. No new analytic derivation is needed beyond `latex/fixed_povm_relative_entropy_derivative.tex`.

## Extension to 200 kicks and exponential fits

### Bound for a future repeated-measurement experiment
[Proved] For any common instrument protocol applied to the two initial states, the full classical history distributions satisfy D(P_n||Q_n) <= D(rho||V rho V†), by quantum relative-entropy data processing. This bounds history relative entropy, not automatically its initial-state Hessian, and is independent of kappa, time, and measurement strength. No repeated weak-measurement simulation has yet been performed.

[Numerically evaluated] For the parameters above the quantum bound is 0.04460888859994434 nats. Evaluated with SciPy expm/logm in dimension 21; the computed imaginary residual was 1.3e-16. As an independent reproducible scalar expression, for rho=0.7|psi><psi|+0.3I/21 the eigenvalue ratio is 50, and the bound is 0.7 log(50)(1-|<psi|V|psi>|²). The squared overlap is 0.9837099670360365. Equivalently it is [1-(1-sin(1.1)^2 sin(0.7)^2)sin(0.07/2)^2]^20. The common measurement channel may yield a strictly smaller attainable value; this bound is not a saturation prediction.

[Numerically tested] Extended the same parameters and direction to t=0,...,200 (402 cases), preserving the original 30-kick files. The preceding 30-kick limitation describes the original experiment; this extension changes only the time horizon. Run `python calculations/kicked_top_growth.py`. It reuses `kicked_top_frechet.main(max_time=200)`. New raw outputs are `calculations/kicked_top_frechet_200_results.csv`, `calculations/kicked_top_frechet_200_summary.json` and `calculations/kicked_top_frechet_200.png`.

Plots: `calculations/kicked_top_growth_semilog.png` and `calculations/kicked_top_growth_early.png`. The latter was visually inspected. They show divergence, absolute first derivative (sign discarded only for plotting), and Hessian. Lines connect actual integer-time samples as guides; no fractional-kick dynamics or additional interpolated data are asserted.

Fit output: `calculations/kicked_top_growth_fits.csv`. Fit log(y)=log(C)+lambda*t by ordinary least squares, for divergence and Hessian separately, in the prespecified windows [0,10], [0,20], [0,30], [10,50], [0,200]. No additive offset, smoothing, extrapolation, or selected rising subsequence. R-squared is calculated in log space, not the original scale. An exact exponential toy case checks rate, intercept and R-squared to 1e-12. Time samples are deterministic and correlated; no independent-sample statistical uncertainty is claimed.

Hessian fit results:

| kappa | window | rate per kick | log-space R-squared |
|---|---|---:|---:|
| 0.5 | 0–20 | -0.0384115 | 0.202018 |
| 6.0 | 0–20 | 0.00769781 | 0.0144230 |
| 0.5 | 0–200 | 0.000327730 | 0.00212379 |
| 6.0 | 0–200 | -0.00000976914 | 0.00000271579 |

Conclusion: no convincing sustained exponential growth in these tested windows. The small positive early-window slope at kappa=6 explains only about 1.4% of log-Hessian variation; its sign is not stable under changing the window. This does not exclude a transient in other states, directions, measurements, spin sizes, or perturbation scales, and does not refute classical kicked-top chaos. The quantity is an initial-state directional Hessian, not a classical trajectory-separation derivative.

All finite-difference assertions and state checks still pass. At h=3e-5, maximum absolute first/Hessian errors across the 402 cases are 1.88111e-8 and 8.09769e-7. Environment and all other tolerances are unchanged.
