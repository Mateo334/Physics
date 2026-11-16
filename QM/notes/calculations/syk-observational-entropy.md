# Fixed-POVM observational entropy in SYK

[Numerically tested] Script: `calculations/syk_observational_entropy.py`. Run with `python calculations/syk_observational_entropy.py`. Outputs with the same stem: four-panel PNG, PDF, CSV and `_summary.json`.

## Setup

q=4 Majorana SYK, N=16, full Fock dimension d=256, hbar=J=1. Convention: {chi_i,chi_j}=delta_ij; H=sum over i<j<k<l of J_ijkl chi_i chi_j chi_k chi_l. Independent real Gaussian couplings of variance 6/N^3, NumPy default_rng seed 20261113. Initial state is the occupation vacuum. Exact dense diagonalization evolves a pure state at 801 times Jt in [0,40], with no measurement backaction.

Four fixed POVMs: rank-one occupation projectors; total-number projectors (rank binomial(8,k)); and unsharp versions E_a=0.75 P_a+0.25 Tr(P_a)I/d. Effects are represented by their diagonal entries. Observational entropy is -sum_a p_a ln(p_a/Tr(E_a)), in nats, with zero terms omitted. Volumes use the full Fock space, not the even-parity sector. The state stays even, so ln(d) is only a universal bound, not necessarily reachable for a given measurement; particle number itself is not conserved by this Majorana SYK Hamiltonian.

## Checks and observations

Passed checks of the Majorana anticommutators and Hermiticity, explicit four-Majorana Hamiltonian product, POVM positivity/completeness and volumes, projective versus nonprojective effects, entropy bounds, maximally mixed state, stationary-state toy case and initial entropy. Typical matrix checks use 1e-12 absolute tolerance (some toy checks use NumPy default relative tolerance); evolved normalization error was 1.89e-15, versus a 1e-11 threshold. Odd-sector probability per basis state was below 1e-20. Pure-state density matrices are positive by construction. No time-integration error (spectral evolution); no sampling-convergence or large-N check. The generated PNG was visually inspected.

Curves rise initially and then fluctuate rather than growing monotonically. At Jt=40: occupation 4.3970, particle number 4.8119, unsharp occupation 4.9423, unsharp particle number 5.2023 nats. Both sharp entropies start at zero; both unsharp ones start at 1.9412 nats because measurement noise already obscures the vacuum. These are absolute entropies, not entropy increases with their initial baselines removed. No disorder averaging or chaos inference.

Environment: NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0; float64/complex128. No external literature used.
