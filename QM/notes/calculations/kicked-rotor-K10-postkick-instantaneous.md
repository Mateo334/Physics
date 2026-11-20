# K=10 kicked rotor: instantaneous POVM optimization after kicks 1 and 2

[Numerically tested] `calculations/code/kicked_rotor_K10_postkick_instantaneous.py`; run with `OPENBLAS_NUM_THREADS=1; py -3.13 calculations/code/kicked_rotor_K10_postkick_instantaneous.py`. Output `calculations/results/kicked_rotor_K10_postkick_d64_d512_m4.{json,npz,png}`. Loads earlier optimized four-outcome coefficients `canonical_obsen_coefficient_d64_m4.{json,npz}` and freshly run `canonical_obsen_coefficient_d512_m4.{json,npz}` (12-start, 800-step reduced canonical optimization; d512 best C=6.539866994104168 and projected residual 6.85e-9; d64 C=4.573808360217963). Compact NPZ stores the canonical 2x2 effects, complement weights and 64/512 component state-tangent frames for each time; reconstruct each full effect by `E_i=F A_i F†+alpha_i(I-FF†)`. No full 512x512 four-effect arrays were required for the plot. No data overwritten.

The same initial Gaussian as earlier: q0=pi/2, p0=0, sigma_q=sqrt(hbar_eff/2), hbar_eff=2pi/64 fixed **for both dimensions**, initial preparation at t=0. Standard kicked-rotor model `H_phys=p²/2+K cos(q) sum_n delta(t-n)` with K=10 and **first** delta kick shifted to t=1, then t=2 and t=3; unitary kick in the q grid `exp(-i K cos q/hbar_eff)`, free momentum unitary `exp(-i p_n² Delta t/(2 hbar_eff))`. The q grid has d points and momentum p_n=hbar_eff*n with n in FFT ordering. Changing d from 64 to 512 while holding hbar fixed increases the momentum cutoff (it is *not* the torus convention hbar=2pi/d for d512). At d512 this is a larger momentum cutoff for the same physical packet and kick strength. Each instant t=1+ and t=2+ is evaluated using the **right-hand free-interval derivative**, with H_generator=p²/(2hbar_eff). There is no ordinary instantaneous derivative during ideal delta kicks.

| d | after kick | measured Var(p) | measured Delta H_phys | attained instant R=C_d,4 Delta H_phys/hbar_eff (nats/free time) | edge mass outer 20% of momentum cutoff |
|---:|---:|---:|---:|---:|---:|
| 64 | 1 | 4.38600116 | 0.97226771 | 45.29642578 | 0.7566912 |
| 64 | 2 | 1.85490572 | 1.32726741 | 61.83530392 | 0.1298156 |
| 512 | 1 | 0.48628907 | 6.08757679 | 405.52111631 | 8.5e-30 |
| 512 | 2 | 9.83595040 | 20.72110617 | 1380.32691846 | 2.0e-8 |

**Severe cutoff warning**: the original d64 torus calculation has 75.7% of first-post-kick probability in the outer 20% of its momentum range: its post-kick moments are badly affected by wraparound and should **not** represent the unbounded-momentum kicked rotor at K=10. Increasing the cutoff to d512 with the **same hbar** nearly eliminates edge mass for both kicks and dramatically changes the uncertainty. The d512 coefficient C_d differs because the available POVM dimension differs. Neither huge instantaneous rate is a classical Lyapunov or KS entropy rate: OE is bounded by log d and the displayed free-interval entropy very quickly bends/saturates and can fluctuate.

Plot is not one continuously reoptimized POVM and not one common POVM optimized over (1,3): solid blue uses the POVM individually optimized at 1+ held fixed only through 1<t<2; solid orange uses the different POVM optimized at 2+ held fixed through 2<t<3; dashed blue shows what happens when the 1+ POVM is **retained** after the second kick. The apparent solid-curve discontinuity at t=2 includes both the physical kick and the intentional measurement switch; do not interpret it as an entropy jump of one fixed measurement. A finite-window rate optimized over a *single shared* POVM is a separate problem. Figure inspected visually.

[Proved] During any one free interval, Delta(p²/2) is constant because the free Hamiltonian commutes with its own propagator. Thus for a pure state and POVMs allowed to depend on the evaluation instant, the **supremum of the instantaneous derivative** has one constant value throughout `(1,2)` and another throughout `(2,3)`. The POVM achieving that instantaneous rate generally **must rotate with the state in time**; the fixed-POVM entropy traces in the plot bend and their derivatives are not constant. This is not a rate integrated across the whole two-kick window.

Checks: kicks/free transforms unitary via norm 1 within 1e-11; each state-tangent frame orthonormal within 1e-10; momentum-fourth-moment energy uncertainty agrees with the direct generator norm; direct compressed-effect outcome derivative matches C*DeltaH/hbar within 1e-9, centered fixed-effect entropy finite difference at h=1e-7 agrees within ~3e-6 absolute even at rate ~1380; all plotted probabilities positive and sum to 1 within 1e-10. Canonical optimizer outcomes numerically stable but *not globally certified*.

## Corrected interpretation: fixed-POVM finite-interval slopes

[Numerically tested] `py -3.13 calculations/code/kicked_rotor_fixed_povm_interval_slopes.py` reads the existing compact POVM/state files, evaluates exact endpoint entropies **with the same measurement held fixed** across free evolution and across the second kick, and saves `calculations/results/kicked_rotor_K10_postkick_d64_d512_m4_fixed_interval_slopes.json` (no overwrite). Unit-time free-interval **secant** = `S(2^-)-S(1^+)` or `S(3^-)-S(2^+)`; separately ordinary least-squares fit of 401 free-time points per interval. A second-kick state jump is excluded from each free-interval slope and recorded separately. The first measurement E1 was optimized only for the *instantaneous derivative at 1+*; E2 was optimized only at 2+. Neither was optimized for these interval-average slopes. Both were evaluated before and after the second kick with **no measurement changes** along each full history.

| d | fixed POVM chosen at | secant 1+ to 2- | secant 2+ to 3- | OLS 1..2 | OLS 2..3 |
|---:|---|---:|---:|---:|---:|
| 64 | 1+ | 2.292849 | -0.001006 | 0.481640 | 0.010322 |
| 64 | 2+ | 0.005424 | 2.246561 | -0.006745 | 0.655026 |
| 512 | 1+ | **3.381816** | **0.000471** | **0.001232** | **0.000265** |
| 512 | 2+ | 0.001280 | 3.237073 | 0.000815 | 1.023393 |

At d512 with E1 retained, `S(1+)=2.856035`, `S(2-)=6.237850`, `S(2+)=6.236370`, `S(3-)=6.236841` nats; the entropy is already near ln(512)=6.238325 by the second kick. Strong curvature and fluctuations make OLS slopes differ dramatically from endpoint secants (notably d512 E1, first free interval): **no robust linear-production rate**. Treat both slope definitions as descriptive, not maxima over fixed POVMs for a finite window. The d64 results retain the severe K10 momentum-wraparound artifact.
