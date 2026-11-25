# Resolved rotor: 100-state mean energy uncertainty times ln2

[Numerically tested] User selected a resolved d512 rotor rather than d2 when requesting a 100-initial-coherent-state average of the preceding DeltaH ln2 curve. Agreed ensemble: uniformly random initial position centers, initial mean momentum zero. Same hbar=2pi/64, evaluate immediately after the first kick, K in [4,12].

Script `calculations/code/rotor_energy_100_coherent.py`; outputs `calculations/results/rotor_d512_100coherent_energy_ln2_K4to12.{png,pdf,npz,csv,json}`. Run with `py -3.11` and single-threaded OPENBLAS/OMP. Python 3.11.9, NumPy 2.2.3, Matplotlib 3.10.0. Figure visually checked.

100 independent q0 uniform in [0,2pi), seed 20261123; periodized Gaussian sigma_q=sqrt(hbar/2)=0.221556731363, p0=0. The same packets are reused at every K; this is not a full phase-space ensemble. Free evolution for one unit followed by exp(-i K cos(q)/hbar). Post-kick energy operator H_phys=p^2/2. K grid has 401 points, spacing .02. Quantity plotted is mean of individual pure-state DeltaH times ln2, NOT energy uncertainty of their mixed ensemble. Shading is sample standard deviation divided by sqrt(100), one standard error for the position-center sampling, not a numerical-error bound.

| K | Mean DeltaH ln2 | Standard error |
|---|---:|---:|
| 4 | 1.17072369886 | 0.03519111880 |
| 6 | 2.51634450867 | 0.07215856707 |
| 8 | 4.38901493797 | 0.12337922397 |
| 10 | 6.78953997454 | 0.18889207339 |
| 12 | 9.71816517853 | 0.26870839712 |

The sampled mean increases at every grid step. All 100 states at all 401 K values checked against d1024 at fixed hbar: maximum relative DeltaH difference 3.04e-15, roundoff level. Largest outer-20%-momentum probability at d512: 1.12e-28. Initial energy uncertainty checked against hbar/(2sqrt(2)) to 1e-11, all state normalizations to 1e-11. Full per-state arrays for both cutoffs stored in NPZ.

This replaces the severely aliased d2 model for this task; both dimension AND preparation ensemble changed, so do not attribute the qualitative difference solely to averaging. Monotonicity is observed on the chosen grid, not a general theorem or a chaos diagnosis. DeltaH ln2 is not the optimized observational-entropy rate; no POVM or rate optimization performed here.
