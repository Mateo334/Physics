# d=64 free-rotor Gaussian and four-outcome instantaneous OE coefficient

[Numerically tested] The user specified a much lower dimension and one initial coherent state. We **choose** a d=64 torus quantization, which supplies an extra convention `hbar_eff=2pi/d=0.09817477042468103`; d alone would not fix hbar on the infinite-cylinder kicked rotor. Choose q0=pi/2 and p0=0, with a periodized Gaussian position wavefunction of standard deviation parameter `sigma_q=sqrt(hbar_eff/2)` and momentum labels `n=-32,...,31`, `p_n=hbar_eff*n`; use a unitary FFT (`norm='ortho'`) to get the momentum distribution. Evaluate **before the first delta kick**, during a kick-free interval `H_phys=p²/2`; the time generator is `H_phys/hbar_eff` in `dot rho=-i[H_phys/hbar_eff,rho]`. Hence the number is independent of K and is not a post-kick entropy-production fit.

Code: `calculations/code/kicked_rotor_gaussian_free_obsen_d64.py`, after `calculations/code/canonical_instantaneous_obsen_coefficient.py --dimension 64 --outcomes 4 --restarts 12 --iterations 800`. Saved `calculations/results/canonical_obsen_coefficient_d64_m4.{json,npz}`, `calculations/results/kicked_rotor_gaussian_free_d64_m4.{json,npz}`. Both refuse to overwrite existing outputs.

| Quantity | d=64 value |
|---|---:|
| measured position variance | 0.049087385212340545 |
| measured momentum variance | 0.04908738521234053 |
| inferred `2 Delta q Delta p` | 0.09817477042468109 |
| physical kinetic-energy uncertainty `sqrt(<p^4>-<p²>²)/2` | 0.034710022954362235 |
| rate-generator uncertainty `Delta H_phys/hbar_eff` | 0.3535533905932738 |
| attained four-outcome coefficient `C_64,4` lower bound | 4.573808360217963 |
| attained instantaneous maximum candidate `C Delta H_phys/hbar_eff` | **1.6170854536789228 nats/free-time-unit** |

The discrete distribution is Gaussian to reported precision; the ideal minimum-uncertainty prediction `Delta(p²/2)=hbar_eff/(2sqrt(2))` matches the measured fourth-moment result. `Delta q Delta p=hbar_eff/2` is a *narrow-packet check* with wrapped coordinate variance, not a global angle-operator uncertainty theorem on a circle. Twelve independent/structured reduced-POVM starts yield the same coefficient to ~7e-11, best projected-gradient residual 1.58e-8. The binary-projective analytical **lower bound** coincides numerically with the optimized value, but global optimality is not proved. Compact canonical effects were embedded into an explicit 64x64 four-outcome POVM: smallest eigenvalues >=-3e-16, Hilbert--Schmidt completeness residual 7.59e-16; direct outcome-time derivative gives 1.6170854536789228, central finite-difference h=1e-5 gives 1.617085453364808. State norm and Fourier transform checked. For a later kick-dependent result, specify K and cycle phase and recompute the evolved momentum fourth moment; the ideal delta kick has no instantaneous differentiable generator at the kick itself.

Plot update [Numerically tested]: `py -3.13 calculations/code/plot_kicked_rotor_pre_first_kick_obsen.py` uses the same **fixed** four-outcome POVM for free evolution, preparing the wavepacket at t=0 and placing the first ideal kick at t=1. Saves `calculations/results/kicked_rotor_gaussian_free_d64_m4_pre_first_kick.png` and `.npz` (501 time samples, `0<=t<1`); refuses to overwrite them. S(0)=1.8521981408667316, S(t approximately 0.998)=3.0765771666871866; local tangent at preparation has slope 1.6170854536789228. Figure inspected visually: entropy bends over the free interval; zoom 0..0.06 nearly overlays the initial tangent. This is a free-rotor entropy curve independent of K, **not** post-kick chaos-driven growth or a sustained maximum rate; at t>0 the same POVM need not remain optimal.
