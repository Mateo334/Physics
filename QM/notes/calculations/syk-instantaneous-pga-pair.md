# Two SYK realizations: instantaneous OE, energy uncertainty and time QFI

[Numerically tested] User selected two independent disorder realizations of quartic SYK, N=16, and physical-time QFI. Script: `calculations/code/syk_instantaneous_pga_pair.py`. Run with `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1`, `py -3.11 calculations/code/syk_instantaneous_pga_pair.py`. Outputs: `calculations/results/syk_N16_instantaneous_pga_pair.{json,npz,csv}`. Python 3.11, NumPy 2.2.3, SciPy 1.15.1; complex128.

Setup: occupation vacuum, J=hbar=1; H=sum over increasing quartets J_ijkl chi_i chi_j chi_k chi_l, {chi_i,chi_j}=delta_ij, independent real Gaussian couplings with variance 6/N^3. No per-realization spectral rescaling. Full Fock dimension d=256, ten effects, volumes defined using full-space trace. Optimized effects commute with parity, but replacing the full-space trace by the even-sector trace would change this optimization problem.

PGA reuses `canonical_instantaneous_obsen_coefficient.run`, with the exact two-dimensional pure-state tangent reduction and scalar complement, four starts, 300-iteration cap, optimizer seed 20261123. One start is the binary-seed POVM with proportional effect splitting; three are random. Dykstra PSD/completeness projection plus simplex projection, with the full Hilbert–Schmidt metric. All starts stopped with no improving step; attained unit-energy-spread rates ranged from 5.878322995596 to 5.878322996400. Best projected residual at step .03: 1.66e-8. Final roundoff PSD repair followed by common normalization. Best feasible coefficient 5.878322996400126. This is an attained lower bound, not a globally certified optimum; ten nominal outcomes need not give ten independent measurement directions.

| Disorder seed | Mean energy | DeltaH | Time QFI | Attained OE rate (nats/time) |
|---|---:|---:|---:|---:|
| 20261113 | -0.05268938294 | 0.414712131595 | 0.687944608368 | 2.437811860041 |
| 20261123 | -0.11472958931 | 0.387457031901 | 0.600491806277 | 2.277597580739 |

Full 256-dimensional effects were reconstructed separately for each Hamiltonian and tested at t=0. Completeness errors below 3.22e-15, minimum effect eigenvalues above -2.28e-16 (roundoff). Central entropy differences at h=1e-5 agree with analytic rates to 1.30e-10 absolute. State normalization, Hamiltonian Hermiticity, parity conservation, effect parity commutation and positive probabilities passed asserts. Time QFI computed from evolved state tangents at 81 times in [0,40] agrees with 4 DeltaH^2 within 1.34e-15.

[Proved, using the existing pure-state reduction] Both rates have R/DeltaH=5.878322996400124. The common coefficient follows from unitary equivalence of pure-state tangent pairs, not from a comparison of chaoticity. DeltaH and time QFI are constant under these autonomous Hamiltonians. The supremum of instantaneous rates is also constant if the measurement can be reoptimized at each time; this does not imply constant growth for one fixed POVM.

Finite-size diagnostic only: mean adjacent-gap ratios in the central 96 eigenvalues of the even-parity sector are 0.498117 and 0.514016. These two finite samples are compatible with level repulsion but do not establish chaos or an asymptotic exponent. No disorder-ensemble study or Lyapunov computation performed. No new substantive derivation or LaTeX deliverable was needed.
