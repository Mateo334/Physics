# Matched d=256, ten-outcome instantaneous OE comparison

[Numerically tested] `calculations/code/kicked_models_d256_m10_compare.py`; outputs `calculations/results/kicked_models_d256_m10_compare.{json,npz,csv}`. Run: `py -3.11 calculations/code/kicked_models_d256_m10_compare.py`, with OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1. NumPy 2.2.3, SciPy 1.15.1, complex128. Reuses saved SYK results without overwriting them.

All models use pure states, full d=256 trace volumes and ten nominal POVM outcomes, unrestricted effects. A fresh canonical reduced PGA search used four starts (one binary seed with effect splitting, three random), 300-iteration cap, seed 20261124. Attained coefficient 5.878322996400126, best residual 1.66e-8. All starts stopped with no improving step, coefficients in [5.878322995567,5.878322996400]. Reconstructed full effects separately for each state and generator. This shares one canonical optimization across models by exact unitary equivalence; it is not four independent full-matrix searches. The binary-seeded best candidate matches the earlier SYK candidate.

Top: j=127.5, initial coherent state theta=1.1, phi=0.7, K=3,6. Half-period rotation H=pi Jy then half-period twist H=K Jz^2/j, hbar=1. Evaluate t=1.1, in the second rotation phase.

Rotor: K=6,10, Gaussian q0=pi/2, p0=0, sigma_q=sqrt(hbar/2), hbar fixed at 2pi/64. This is a d=256 momentum cutoff, not the hbar=2pi/256 torus convention. Free evolution first, delta kick at t=1; evaluate free right-hand derivative at t=1+, H_phys=p^2/2. Central finite-difference checks extend the free branch backward and forward locally; they do not differentiate across the kick.

| Model | Physical DeltaH | DeltaH/hbar | Physical-time QFI | Attained rate (nats/time) |
|---|---:|---:|---:|---:|
| Top K3 | 12.4845111778 | 12.4845111778 | 623.452077392 | 73.3879891551 |
| Top K6 | 80.8130587168 | 80.8130587168 | 26123.0018366 | 475.045261464 |
| Rotor K6 | 2.39743888903 | 24.4201119968 | 2385.36747974 | 143.549305925 |
| Rotor K10 | 6.08757678583 | 62.0075479626 | 15379.7440173 | 364.500395139 |
| SYK seed 20261113 | 0.414712131595 | 0.414712131595 | 0.687944608368 | 2.43781186004 |
| SYK seed 20261123 | 0.387457031901 | 0.387457031901 | 0.600491806277 | 2.27759758074 |

Every candidate has rate/(DeltaH/hbar)=2 rate/sqrt(QFI)=5.8783229964, agreeing across all cases within 7e-14. Physical energy and time scales are model-dependent; rate/DeltaH alone is not comparable when hbar differs.

Checks: full-effect completeness residuals below 3.42e-15, minimum effect eigenvalues above -1.06e-15 (roundoff), normalized pure states, Hermitian generators, positive probabilities. Entropy central differences with dt=1e-5/max(1,DeltaH/hbar) agree to relative errors below 1.7e-10. Rotor d512 check at the SAME hbar: relative generator-spread change 4.44e-16 (K6), 4.01e-12 (K10); probability outside d256 momentum indices below 5.24e-12. The larger cutoff checks the dynamics at the selected time only, not later kicks or the dimension-dependent coefficient.

[Proved, existing result] Universality of the exact pure-state unrestricted instantaneous supremum coefficient follows from unitary equivalence of state-tangent pairs; see `../theorems/observational-entropy-qfi-pesin.md`. Numerical agreement here verifies transported candidates and normalization, not an independent proof of global optimality. No global optimum certificate, new Lyapunov calculation, or sustained fixed-POVM entropy-growth claim. Ten effects may contain proportional outcomes. No new derivation or LaTeX document required.
