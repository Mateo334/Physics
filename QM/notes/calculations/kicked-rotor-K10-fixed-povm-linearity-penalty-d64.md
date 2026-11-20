# K=10, d=64: fixed four-outcome POVM with derivative-linearity penalty

[Numerically tested] Code: `calculations/code/kicked_rotor_K10_linearity_penalty_d64.py`. Saved outputs: `calculations/results/kicked_rotor_K10_d64_fixed_povm_linearity_mu1_50.{json,npz,png}`. Run: `OPENBLAS_NUM_THREADS=1; py -3.13 calculations/code/kicked_rotor_K10_linearity_penalty_d64.py --mu 1 --restarts 4 --steps 400`. Script refuses to overwrite these results. The saved plot was visually inspected; its derivative has a sharp initial peak and a sharp drop close to the second kick, even though the entropy curve appears nearly straight.

The previous Gaussian state, K=10, d=64, and hbar_eff=2pi/64 were retained. Preparation at t=0, first delta kick at t=1, second at t=2; free propagation between kicks. One **fixed**, not time-dependent, four-effect POVM is used for the complete unmeasured state evolution; the objective probes only the interval t=1+ through t=2-. Entropy is `S=-sum_i p_i log(p_i/Tr E_i)` (natural log); derivative is the analytic free-Hamiltonian derivative. Objective `J=R-mu*P`, `R=S(2-)-S(1+)`, `P=int_0^1 (dot S(1+s)-R)^2 ds`, with **user-unspecified mu chosen as 1**. The discrete objective uses 50 equally spaced samples including both one-sided endpoints and trapezoidal quadrature. No physical measurement/backaction is performed during the evolution.

Exact reduction for this interval objective: the post-first-kick free orbit lies in a 33-dimensional subspace `F` (distinct squared momentum eigenvalues). Twirling full 64-dimensional effects on its orthogonal complement preserves all free-orbit outcome probabilities and the effect volumes, so an optimizer can be taken in the form `E_i=F A_i F†+alpha_i(I-FF†)` with `sum_i A_i=I33`, `sum_i alpha_i=1`. This reduction does not constrain the four-effect optimum for the stated interval objective, though the effects still define the same fixed measurement before kick 1. Analytical matrix and complement-weight gradients were checked by a random-direction centered finite difference (relative discrepancy 2.0e-10 on the run), and the 33x33 effects were projected with Dykstra; complement weights with simplex projection. Four seeded random starts (seed 20261128; the separate gradient-check draws from the same generator), 400 accepted steps each. All hit the iteration limit.

| Quantity | Best attained value |
|---|---:|
| `R=S(2-)-S(1+)` | 2.28173594 nats/free-time |
| `S(1+)`, `S(2-)` | 1.86678543, 4.14852138 nats |
| 50-point `P` | 0.04214406 |
| 50-point RMS derivative deviation `sqrt(P)` | 0.20529019 nats/free-time |
| 50-point discrete objective `J=R-P` | 2.23959188 |
| 501-point check of `P` (same fixed effects) | 0.03320819 |
| 501-point RMS deviation | 0.18223114 nats/free-time |
| maximum absolute 501-point derivative deviation | 1.20460588 nats/free-time |
| 50-point derivative trapezoid minus exact secant | -0.00746355 |
| best projected-gradient residual (step 0.03) | 1.14246 |

**Limitations:** [Numerically tested] The 50-point penalty differs from the 501-point check by about 27% of the fine-grid value (or the fine-grid estimate is about 21% below the coarse result). Neither quadrature certifies the continuous integral. [Unresolved] Because the best start was still improving at 400 steps, with a substantial projected-gradient residual, `R`, `J`, and `P` are **attained candidate values, not a proven or converged optimum**; 10-start/600-step attempt timed out before saving and is not counted. Compared to the *user-provided* classical lambda=1.619, the attained secant R is ~0.663 higher (about 41%), but this is not a justified chaos/Lyapunov identification: a nonstationary one-interval measured-entropy slope need not equal a classical exponent. [Numerically tested] After first kick, 75.67% of the momentum probability lies in the outer 20% of the d=64 momentum range, so wraparound makes this finite-torus model unreliable as an approximation of an unbounded K=10 rotor.

## The same best attained POVM plotted from preparation through kick 3

[Numerically tested] `calculations/code/kicked_rotor_K10_linearity_best_povm_full_curve.py` loads exactly the saved best-effects NPZ without reoptimizing and plots `calculations/results/kicked_rotor_K10_d64_fixed_povm_linearity_mu1_full_evolution.png` (visually inspected). The reconstructed compact-effect formula uses the **33-dimensional** complement volume `Tr A_i+31 alpha_i`; it agrees with stored first-to-second-kick endpoint entropy to 1e-10. Values: `S(0)=4.12350561`, `S(1-)=4.11912859`, `S(1+)=1.86678543`, `S(2-)=4.14852138`, `S(2+)=4.15631836`, `S(3-)=4.14511469`, `S(3+)=4.15168776` nats. The steep drop at kick 1 and near-flat sections elsewhere demonstrate that only the highlighted first-to-second free interval has the targeted growth; the POVM is nevertheless unchanged at all times. No backaction from observations is applied.
