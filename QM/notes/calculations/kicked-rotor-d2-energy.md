# Two-state rotor energy uncertainty, K=8 through 13

[Numerically tested] User explicitly selected keeping the previous hbar=2pi/64 rather than using a d2 torus hbar=pi. Computed immediately after the first kick, t=1+, using the previous periodized Gaussian (q0=pi/2, p0=0, sigma_q=sqrt(hbar/2)) sampled on two position points q=0,pi. FFT momentum indices 0,-1, physical kinetic Hamiltonian H=p^2/2. Free evolution for one unit precedes the kick exp(-iK cos(q)/hbar).

Script `calculations/code/kicked_rotor_d2_energy_uncertainty.py`; outputs `calculations/results/kicked_rotor_d2_K8to13_energy.{json,csv}`. Reuses `rotor_state` in `kicked_models_d256_m10_compare.py`, which uses the established rotor preparation. Python 3.11.9, NumPy 2.2.3. Run `py -3.11 calculations/code/kicked_rotor_d2_energy_uncertainty.py`.

| K | Physical DeltaH | DeltaH/hbar |
|---|---:|---:|
| 8 | 0.000911980550870 | 0.00928935761118 |
| 9 | 0.00218346340689 | 0.0222405756330 |
| 10 | 0.00112377867523 | 0.0114467155907 |
| 11 | 0.00207445575653 | 0.0211302328242 |
| 12 | 0.00132500294898 | 0.0134963692123 |
| 13 | 0.00194592915129 | 0.0198210715734 |

The sampled Gaussian is constant on the two position points, so the initial momentum state is zero momentum. The post-kick negative-momentum probability is sin^2(K/hbar), giving exact two-state check DeltaH=(hbar^2/4)|sin(2K/hbar)|. Numerical differences below 9e-19, normalization within 1e-12. The omitted physical kick dynamics is severely aliased: this is a finite two-state grid toy, NOT a converged kicked rotor. Nonmonotonicity here must not be interpreted as a chaos signature. No delta-impulse energy uncertainty is defined or computed; all values refer to the following free phase.

Follow-up: user requested DeltaH times natural log(2), then a noninteger K sweep from 4 to 12. Evaluated the verified analytic two-state expression at 4001 equally spaced K values (spacing .002). Plot/data: `calculations/results/kicked_rotor_d2_DeltaH_ln2_K4to12.{png,csv}`. Reproduction script: `calculations/code/plot_rotor_d2_energy_ln2.py`, refuses overwrite. Same fixed hbar, preparation and post-first-kick convention; no continuum-rotor inference.
