# SYK single-Majorana final-probe QFI

[Numerically tested] `calculations/code/syk_majorana_probe_qfi.py`; run `py -3.11 calculations/code/syk_majorana_probe_qfi.py`. Outputs `calculations/results/syk_N16_majorana_probe_qfi.{png,csv,json}`; refuses overwrite. Reuses the Jordan-Wigner action from `four_model_setup.majorana`, but not its energy rescaling or restricted even-sector unitary.

N=16, q=4, full dimension 256, hbar=J=1. H=sum_{i<j<k<l} J_ijkl chi_i chi_j chi_k chi_l, real independent Gaussian couplings with variance 6/N^3, default_rng seed 20261113. Initial occupation vacuum. Spectral evolution at 801 equally spaced times Jt in [0,40]. Hypothetical final family exp(-i theta chi_1)exp(-i Ht)|vacuum>; no probe actually applied and no time-parameter QFI.

[Proved] With {chi_i,chi_j}=delta_ij, chi_1^2=I/2. Quartic H preserves parity, while chi_1 is parity odd, so its mean vanishes in the evolved vacuum. Pure-state QFI=4 Var(chi_1)=2 identically. This holds for any parity-definite pure initial state and any even Hamiltonian, not just chaotic SYK. Computing the variance requires the full-space action: projecting chi into the even sector would incorrectly give zero.

[Numerically tested] Hermiticity, chi square, parity commutator/anticommutator, state parity, norm and exact QFI checked with tolerances 1e-12 to 1e-10. Maximum norm error 1.33e-15; maximum deviation from QFI=2 was 2.66e-15. Python 3.11.9, NumPy 2.2.3, Matplotlib 3.10.0, complex128; plot visually inspected. No disorder averaging or exponent fitting.

Physical caveat: the odd probe mixes parity sectors. Its unrestricted mathematical QFI presumes access to this encoding and suitable measurements; it is not automatically operational under strict fermion-parity superselection without a parity-exchanging reference/ancilla.
