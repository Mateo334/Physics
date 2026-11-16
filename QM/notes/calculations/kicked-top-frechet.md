# Kicked-top initial-state Frechet check: current summary

## Reproduce
- `python calculations/code/kicked_top_frechet.py`: original 30-kick experiment.
- `python calculations/code/kicked_top_growth.py`: extension through 200 kicks and exponential fits.
- Results in `calculations/results/`: `kicked_top_frechet*` CSV/JSON/PNG and `kicked_top_growth*` figures/fit CSV. Original 30-kick outputs were preserved.

## Setup and findings
- j=10, d=21; kappa=0.5,6; rotation π/2 followed by Jz² twist. Fixed final Jz projective measurement; no histories or backaction.
- Initial state: 70% coherent projector (theta=1.1, phi=0.7), 30% I/21. Preparation perturbation is a Jy rotation by 0.07. One traceless Hermitian initial-state direction, operator norm 1, seed 20261112.
- Vary the initial rho in both arguments, holding time, unitary and measurement fixed. These are not time or epsilon derivatives.
- [Numerically tested] First differential and Hessian agree with central differences across 402 cases through 200 kicks. At h=3e-5, maximum absolute errors: 1.88111e-8 and 8.09769e-7. State, unitarity and epsilon=0 checks pass.
- No convincing sustained exponential growth in tested windows. For kappa=6, the 0–20 Hessian log-fit rate is 0.00769781 per kick, with R²=0.0144230: not a Lyapunov exponent. This does not refute classical chaos or exclude other parameter regimes.
- [Proved bound; numerically evaluated value] Common-instrument full-history relative entropy is bounded by the initial quantum relative entropy, here 0.04460888859994434 nats. This does not automatically bound its initial-state Hessian or predict saturation. Repeated weak measurements were not simulated.

See [full historical details](../archive/kicked-top-frechet-history.md) for formulas, convergence tables, fit windows, tolerances and limitations; [theorem note](../theorems/fixed-povm-relative-entropy-derivative.md) for the analytic statements.
