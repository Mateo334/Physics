# Gaussian stretching versus state distance; ordered derivative check

[Numerically tested] 2026-11-24. Code: `../../calculations/code/quantum_tangent_geometry.py`.

Run from workspace: `py -3.11 calculations/code/quantum_tangent_geometry.py`. Refuses to overwrite any of its three outputs.

## Parameters and outputs
- Exact Gaussian formulas for H=(P²−Q²)/2, ħ=1, initial covariance I/2, initial displacement direction (1,1)/sqrt(2), 121 times from 0 to 3. Squared overlap uses initial displacement 0.1. No Hilbert-space truncation, numerical ODE solver or random Gaussian sampling.
- Independent 5×5 random Hermitian Q,G; NumPy default_rng seed 20261124; V=[Q,G]/i. Central-difference step 1e−5 tests the derivative of Q³. This finite matrix test is not a representation of canonical commutation.
- Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0; float64/complex128.
- Results: `../../calculations/results/quantum_tangent_geometry.json`, `.csv`, `.png`.

## Checks
- Centroid tangent agrees with exp(t) to relative tolerance 1e−12.
- Displacement QFI remains 2; maximum absolute error 8.8818e−16.
- Maximum symplectic residual 1.2057e−13.
- Maximum covariance determinant deviation from 1/4: 2.4075e−12; all sampled covariances positive.
- Ordered cubic derivative agrees with commutator Leibniz identity to relative error 1.7505e−16; central difference error 1.8175e−9.
- Incorrect commuting derivative 3Q²V has relative error 1.98709 against the correct Q²V+QVQ+VQ².

These are numerical sanity checks of identities proved in `../../latex/quantum_hamiltonian_tangent_geometry.tex`, not evidence of chaotic dynamics. The integrable inverted oscillator deliberately separates instability from chaos. No nonlinear-model convergence or semiclassical-time-window study was attempted.
