# Fixed-POVM gradient-optimized kicked-top observational-entropy slopes

[Numerically tested] Script: `calculations/code/optimize_kicked_top_obsen_rate.py`; outputs: `calculations/results/kicked_top_obsen_optimized_j20_m4.json`, `kicked_top_obsen_optimized_k{3,6}_j20_m4.npz`, and separate common-window `kicked_top_obsen_optimized_k3_j20_m4_T3.json/.npz` (the latter from importing `optimize(3.,3,np.random.default_rng(20261121))`). Original `grad_obsen.py` optimizes a continuous-time Hamiltonian derivative, so is not directly applicable to discrete kicks.

Setup: spin j=20 (d=41), pure coherent initial state theta=1.1, phi=0.7; U=exp(-ik Jz²/(2j)) exp(-i pi Jy/2). Four-outcome time-independent POVM, final-only probabilities and no measurement backaction. E_i=M^{-1/2}A_i†A_i M^{-1/2}, M=sum A_i†A_i; S(t)=-sum_i p_i(t) log[p_i(t)/Tr(E_i)]. Optimize the ordinary least-squares slope over integer t=0,...,T using complex128 PyTorch automatic differentiation and Adam, 100 iterations, learning rate 0.05, two random starts, seed 20261120. Not an optimization over time-dependent POVMs or a globally certified maximum. Heuristic Ehrenfest reference t_E=ln(d)/ln(k/2), *not* a measured break time; especially unreliable as a Lyapunov estimate at k=3. T=floor(t_E) gives T=9 and 3, respectively.

| k | Fit t | best slope (nats/kick) | optimized S(t), abbreviated |
|---|---|---:|---|
| 3 | 0–9 | 0.47673749 | S(0)=0.01591; S(5)=3.67734; S(9)=3.65105 |
| 6 | 0–3 | 1.47436542 | S(0)=0.00782; S(1)=0.02421; S(2)=3.69030; S(3)=3.70035 |
| 3 | 0–3 (matched-window auxiliary run) | 1.47915182 | S(0)=0.00129; S(1)=0.02759; S(2)=3.68692; S(3)=3.71203 |

For comparison, fixed contiguous Jz four-bin PVM slopes for k=3/T=3, k=3/T=9, k=6/T=3 are respectively 0.14618, 0.11622, 0.30522 nats/kick. Random-restart optimized slopes at k=3/T=9: 0.47673749 and 0.41440014; k=6/T=3: 1.47436542 and 1.47436285; k=3/T=3: 1.47915182 and 1.20075825. Software: NumPy 2.3.5, SciPy 1.16.3, PyTorch 2.9.1+cpu; double complex, one CPU thread. POVM completeness/positivity, unitary and state norms, entropy bounds and fit formula asserted. Autograd directional derivative tested against centered differences (h=1e-5, independent seed 12): relative error ~2.6e-7.

Interpretation: This unrestricted Hilbert-space POVM can distinguish the early orbit states and make entropy jump sharply from near zero to near ln(41), even at k=3. Thus fitted slopes are **not evidence of a linear growth regime, a Lyapunov exponent, or a chaos contrast**; unequal fit windows cannot be compared as intrinsic rates. A physically restricted phase-space POVM or a fixed POVM shared across k and an initial-state ensemble would answer a different and more discriminating question. Finite random starts and iterations give only achieved lower bounds on the maxima.
