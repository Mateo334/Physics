# Fixed-POVM entropy curves with best attained initial rate

[Numerically tested] Requested plots of observational entropy for models realizing the optimized rate. Script `calculations/code/plot_obsen_attained_rate_models_m4.py`; outputs `calculations/results/obsen_attained_rate_models_d256_m4.{png,pdf,npz,json}`. Run with `py -3.11` and single-threaded OPENBLAS/OMP. Python 3.11.9, NumPy 2.2.3, Matplotlib 3.10.0. Figure visually inspected.

Uses d=256 and four outcomes, A_d256/alpha_d256 from `obsen_constant_m4_dimension.npz`, transported with the saved state-tangent frames from the SYK and kicked-model computations. No new optimization. Each model has its own POVM optimized at its initial plotted instant, then held fixed; final-only measurement with no backaction. Trace volumes include the full 256-dimensional Hilbert space. Common attained coefficient R/(DeltaH/hbar)=5.8783229964001205. Initial observational entropy 2.518546192812 nats for all three.

| Model | Initial physical time | Attained initial rate | Plotted elapsed-time endpoint |
|---|---:|---:|---:|
| N16 quartic SYK, seed 20261113 | 0 | 2.43781186004 | 4.82262236291 |
| Kicked top j=127.5, K6 | 1.1 | 475.045261464 | 0.024748475454 |
| Kicked rotor d256, K10 | 1+ | 364.500395139 | 0.032254137854 |

Each panel spans elapsed time times generator uncertainty from 0 to 2, with 801 samples, but axes display actual elapsed time. Top remains entirely in rotation phase H=pi Jy; rotor remains entirely in free phase H_phys=p^2/2, hbar=2pi/64. SYK uses its full saved Hamiltonian. Orange initial tangents displayed only through dimensionless elapsed time .4. Dotted line is ln(256). Curves bend away from the initial tangent; no sustained maximal rate is claimed.

Checks: state normalization, frame orthonormality, reduced effect positivity/completeness and scalar simplex, probability normalization/positivity, entropy between 0 and ln256. Central derivative step=1e-5/generator uncertainty, relative agreement below 1.6e-10. For rotor the negative-time check is a mathematical continuation of the free branch, not a derivative crossing the kick. Rotor cutoff was checked in the prior comparison; no new cutoff enlargement here.

The rates remain numerical lower bounds on the global instantaneous optimum, not certified maxima. The plot shows fixed measurements reaching those rates at the marked initial instants, not a time-dependent measurement that maximizes continuously.
