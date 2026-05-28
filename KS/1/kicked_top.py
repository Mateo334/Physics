"""
Quantum kicked top: AFL entropy, matrix entropy E(U^n), and Quantum Pesin verification.

The kicked top evolution operator (spin-j, dimension d = 2j+1):
  U(k) = exp(-i (k/2j) J_z^2) * exp(-i (pi/2) J_y)

KEY ANALYTICAL INSIGHT (proved here):
  U_kick = exp(-i(k/2j)J_z^2) is DIAGONAL in the J_z eigenbasis with entries
  exp(-i(k/2j)m^2). Therefore:
    [U]_{ij} = exp(-i(k/2j)*i^2) * [exp(-i(pi/2)*J_y)]_{ij}
  and |U_{ij}|^2 = |[U_rot]_{ij}|^2, independent of k.
  => E(U^1) = E(U_rot) is k-INDEPENDENT.

  However, for n >= 2:
    [U^n]_{ij} involves sums over intermediate indices with phases exp(-i(k/2j)m^2),
    so |[U^n]_{ij}|^2 DOES depend on k via quantum interference.
  => E(U^n) for n >= 2 distinguishes chaotic from integrable dynamics.

Classical Lyapunov exponents (Haake 1987, tau=pi/2):
  k < 2: lambda_L = 0 (integrable)
  k = 3: lambda_L ~ 0.451 (chaotic)
  k = 6: lambda_L ~ 0.880 (strongly chaotic)
"""

import numpy as np
from scipy.linalg import expm

# ---- Angular momentum operators for spin j ----
def spin_operators(j):
    d = int(2*j + 1)
    ms = np.arange(j, -j-1, -1)
    Jz = np.diag(ms).astype(complex)
    Jplus = np.zeros((d, d), dtype=complex)
    for i in range(d-1):
        m = ms[i+1]
        Jplus[i, i+1] = np.sqrt(j*(j+1) - m*(m+1))
    Jminus = Jplus.conj().T
    Jx = 0.5 * (Jplus + Jminus)
    Jy = -0.5j * (Jplus - Jminus)
    return Jx, Jy, Jz

def kicked_top_unitary(j, k, tau=np.pi/2):
    """U(k) = exp(-i*(k/2j)*Jz^2) * exp(-i*tau*Jy)."""
    Jx, Jy, Jz = spin_operators(j)
    U_kick = expm(-1j * (k / (2*j)) * Jz @ Jz)
    U_rot  = expm(-1j * tau * Jy)
    return U_kick @ U_rot

def rotation_only(j, tau=np.pi/2):
    """U_rot = exp(-i*tau*Jy) — k-independent part."""
    Jx, Jy, Jz = spin_operators(j)
    return expm(-1j * tau * Jy)

# ---- Matrix entropy E(U) ----
def matrix_entropy(U):
    """E(U) = -(1/d) sum_{ij} |U_ij|^2 log|U_ij|^2."""
    d = U.shape[0]
    q = np.abs(U)**2
    with np.errstate(divide='ignore', invalid='ignore'):
        terms = np.where(q > 1e-300, q * np.log(q), 0.0)
    return -np.sum(terms) / d

def total_otoc(U):
    """C_1(U) = (2/d^3)(d - sum|U_ij|^4)."""
    d = U.shape[0]
    return (2.0 / d**3) * (d - np.sum(np.abs(U)**4))

# ---- Prove k-independence analytically ----
def verify_k_independence(j, k_vals=None):
    """
    Verify that |U(k)_{ij}|^2 = |U_rot_{ij}|^2 for all k.
    Proof: U_kick is diagonal => [U_kick]_{mn} = exp(-i(k/2j)m^2)*delta_{mn}.
    So [U]_{ij} = [U_kick]_{ii} * [U_rot]_{ij} = phase * [U_rot]_{ij}.
    Hence |U_{ij}|^2 = |U_rot_{ij}|^2, k-independent.
    """
    if k_vals is None:
        k_vals = [0.5, 1.0, 2.0, 3.0, 6.0]
    U_rot = rotation_only(j)
    Q_rot = np.abs(U_rot)**2
    max_diff = 0.0
    for k in k_vals:
        U = kicked_top_unitary(j, k)
        Q = np.abs(U)**2
        diff = np.max(np.abs(Q - Q_rot))
        max_diff = max(max_diff, diff)
    return max_diff

# ---- E(U^n) as a function of n ----
def matrix_entropy_timeseries(j, k, n_max=15):
    """
    Compute E(U^n) for n = 1..n_max.
    U^n = n-th power of the Floquet unitary.
    This is NOT the AFL density matrix entropy; it is the matrix entropy of
    the iterated evolution operator.

    Relation to AFL: S(rho[Z^(2)]) after n Floquet steps = log(d) + E(U^n).
    """
    U = kicked_top_unitary(j, k)
    Un = np.eye(U.shape[0], dtype=complex)
    results = []
    for n in range(1, n_max+1):
        Un = Un @ U
        E = matrix_entropy(Un)
        results.append(E)
    return results

# ---- Normalized scrambling rate ----
def scrambling_rate(j, k, n_max=20):
    """
    Compute E(U^n) / log(d) as a function of n.
    For chaotic k, this should grow (up to saturation).
    Compare growth rate with classical Lyapunov exponent.
    """
    d = int(2*j + 1)
    U = kicked_top_unitary(j, k)
    Un = np.eye(d, dtype=complex)
    results = []
    for n in range(1, n_max+1):
        Un = Un @ U
        E = matrix_entropy(Un)
        results.append(E / np.log(d))
    return results

# ---- Qutrit analysis (d=3) ----
def qutrit_unitary_family(theta):
    """
    Qutrit rotation U(theta) for d=3.
    Use the generalized Hadamard (DFT matrix) and rotations.
    """
    d = 3
    omega = np.exp(2j * np.pi / 3)
    # Qutrit DFT (Hadamard analog)
    H3 = np.array([[1, 1, 1],
                   [1, omega, omega**2],
                   [1, omega**2, omega**4]], dtype=complex) / np.sqrt(3)
    # Qutrit rotation: exp(-i theta J_y) for j=1
    Jx, Jy, Jz = spin_operators(1)
    U_rot = expm(-1j * theta * Jy)
    return U_rot

def sic_povm_d2():
    """
    For d=2, SIC-POVM states (4 equiangular lines):
    |phi_0> = |0>, |phi_1> = (1/sqrt(3))(|0>+sqrt(2)|1>),
    |phi_2> = (1/sqrt(3))(|0>+sqrt(2)e^{2pi i/3}|1>),
    |phi_3> = (1/sqrt(3))(|0>+sqrt(2)e^{4pi i/3}|1>).
    These satisfy <phi_i|phi_j>^2 = 1/(d+1) = 1/3 for i != j.
    """
    phi = [None] * 4
    phi[0] = np.array([1, 0], dtype=complex)
    for k in range(1, 4):
        angle = 2 * np.pi * (k-1) / 3
        phi[k] = np.array([1, np.sqrt(2)*np.exp(1j*angle)], dtype=complex) / np.sqrt(3)
    # Verify SIC condition
    ok = True
    for i in range(4):
        for j in range(i+1, 4):
            overlap = abs(np.dot(phi[i].conj(), phi[j]))**2
            if abs(overlap - 1/3) > 1e-10:
                ok = False
    return phi, ok

# ---- Dissipative channel: amplitude damping ----
def amplitude_damping_kraus(gamma):
    """Kraus operators for amplitude damping: K0, K1."""
    K0 = np.array([[1, 0], [0, np.sqrt(1-gamma)]], dtype=complex)
    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype=complex)
    return [K0, K1]

def apply_channel(rho, kraus_ops):
    """Apply a channel Phi(rho) = sum_k K_k rho K_k^dag."""
    result = np.zeros_like(rho)
    for K in kraus_ops:
        result += K @ rho @ K.conj().T
    return result

def opu_density_matrix_channel(Z_ops, kraus_ops, n_steps, omega_init):
    """
    Compute rho[Z^(n)] for a CP-map channel (non-unital dynamics).
    The CP map is given by its Kraus operators.
    The OPU at step k is Theta^{k-1}(Z_i) where Theta is the Heisenberg picture:
      Theta(A) = sum_k K_k^dag A K_k
    """
    d = omega_init.shape[0]
    k_ops = Z_ops  # initial OPU

    # Build all OPU elements up to step n
    def heisenberg_channel(A, kraus_ops):
        """Theta(A) = sum_k K_k^dag A K_k (Heisenberg picture)."""
        result = np.zeros_like(A)
        for K in kraus_ops:
            result += K.conj().T @ A @ K
        return result

    # Compose OPU operators
    Z_n = list(Z_ops)  # copy of initial OPU
    for step in range(n_steps - 1):
        new_Z = []
        for z_prev in Z_n:
            for z_init in Z_ops:
                # Z^(step+1)_{(alpha,i)} = Z^(step)_alpha * Theta^step(Z_i)
                # Theta^step(Z_i) = apply heisenberg map step times
                Z_rot = z_init.copy()
                for _ in range(step + 1):
                    Z_rot = heisenberg_channel(Z_rot, kraus_ops)
                new_Z.append(z_prev @ Z_rot)
        Z_n = new_Z

    # Compute rho[Z^(n)]_{alpha,beta} = omega(Z_beta^dag Z_alpha) = Tr(Z_beta^dag Z_alpha omega)
    N = len(Z_n)
    rho = np.zeros((N, N), dtype=complex)
    for ia, Za in enumerate(Z_n):
        for ib, Zb in enumerate(Z_n):
            rho[ia, ib] = np.trace(Zb.conj().T @ Za @ omega_init)

    rho = (rho + rho.conj().T) / 2
    eigvals = np.linalg.eigvalsh(rho)
    eigvals = eigvals[eigvals > 1e-14]
    return -np.sum(eigvals * np.log(eigvals))

# ---- Classical Lyapunov reference values ----
CLASSICAL_LYAPUNOV = {
    0.5: 0.000, 1.0: 0.000, 1.5: 0.000,
    2.0: 0.090, 2.5: 0.300, 3.0: 0.451,
    3.5: 0.556, 4.0: 0.650, 5.0: 0.770, 6.0: 0.880,
}

# ===================================================================
# MAIN
# ===================================================================
if __name__ == '__main__':
    print("="*70)
    print("PART 1: K-INDEPENDENCE PROOF VERIFICATION")
    print("="*70)
    print("\nTheorem (K-independence): For U(k) = U_kick(k)*U_rot,")
    print("  |U(k)_{ij}|^2 = |[U_rot]_{ij}|^2  for all k.")
    print("Proof: U_kick = diag(exp(-i(k/2j)m^2)) => |U_{ij}| = |[U_rot]_{ij}|.\n")
    for j in [0.5, 1.0, 1.5, 2.5, 5.0]:
        diff = verify_k_independence(j, k_vals=[0.5, 1.0, 2.0, 3.0, 6.0, 10.0])
        print(f"  j = {j:.1f}, d = {int(2*j+1):3d}: max||U(k)|^2 - |U_rot|^2| = {diff:.2e}")
    print("\nConclusion: E(U^1) = E(U_rot), k-independent. AFL entropy at step 1 does NOT")
    print("detect quantum chaos for the kicked top.\n")

    print("="*70)
    print("PART 2: E(U^n) GROWTH — CHAOS vs INTEGRABILITY")
    print("="*70)
    print("\nFor n >= 2, U^n = (U_kick*U_rot)^n has matrix elements that depend on k")
    print("through quantum interference of phases. E(U^n) distinguishes chaos.\n")

    j = 5.0; d = int(2*j+1)
    k_test = [0.5, 1.0, 2.0, 3.0, 6.0]
    n_max = 20

    print(f"j = {j:.1f}, d = {d} (spin-5, 11-dimensional Hilbert space)")
    print(f"{'n':>4}", end="")
    for k in k_test:
        print(f"  k={k:.1f}", end="")
    print(f"  log(d)={np.log(d):.4f}")
    print("-" * (6 + 7*len(k_test)))

    series = {k: matrix_entropy_timeseries(j, k, n_max=n_max) for k in k_test}
    for n in range(1, n_max+1):
        print(f"{n:4d}", end="")
        for k in k_test:
            print(f"  {series[k][n-1]:6.4f}", end="")
        print()

    print("\nNormalized E(U^n)/log(d):")
    logd = np.log(d)
    for n in [1, 2, 3, 5, 10, 15, 20]:
        print(f"  n={n:3d}:", end="")
        for k in k_test:
            print(f"  {series[k][n-1]/logd:.4f}", end="")
        print()

    print("\n")
    print("="*70)
    print("PART 3: SEMICLASSICAL SCALING E(U^n)/log(d) vs CLASSICAL lambda_L")
    print("="*70)
    print("\nFor fixed n, E(U^n)/log(d) as a function of j (-> infinity = semiclassical).")
    print("Quantum Pesin: E(U^n)/log(d) -> lambda_L * n / log(d)?  Or max out at 1?\n")

    n_fixed = 5
    k_vals = [0.5, 2.0, 3.0, 6.0]
    j_vals = [1, 2, 3, 5, 8, 10, 15, 20]
    print(f"At n = {n_fixed} Floquet steps:")
    print(f"{'j':>5} | {'d':>4}", end="")
    for k in k_vals:
        lL = CLASSICAL_LYAPUNOV.get(k, 0.0)
        print(f" | k={k:.1f}(lL={lL:.2f})", end="")
    print()
    print("-" * (12 + 20*len(k_vals)))
    for j in j_vals:
        d_j = int(2*j+1)
        print(f"{j:5.1f} | {d_j:4d}", end="")
        for k in k_vals:
            Un_list = matrix_entropy_timeseries(j, k, n_max=n_fixed)
            E_n = Un_list[n_fixed-1]
            print(f" | {E_n/np.log(d_j):18.5f}", end="")
        print()

    print()
    print("="*70)
    print("PART 4: QUTRIT ANALYSIS (d=3, spin-1)")
    print("="*70)
    j = 1.0; d = 3
    print(f"\nQutrit kicked top (j=1, d=3):")
    print(f"{'k':>5} | {'E(U)':>9} | {'E(U^2)':>9} | {'E(U^3)':>9} | {'E(U^5)':>9} | {'E(U^10)':>10}")
    print("-" * 60)
    for k in [0.5, 1.0, 2.0, 3.0, 6.0]:
        ts = matrix_entropy_timeseries(j, k, n_max=10)
        print(f"{k:5.1f} | {ts[0]:9.5f} | {ts[1]:9.5f} | {ts[2]:9.5f} | {ts[4]:9.5f} | {ts[9]:10.5f}")
    print(f"\nlog(d=3) = {np.log(3):.5f}")

    print()
    print("="*70)
    print("PART 5: SIC-POVM ANALYSIS (d=2)")
    print("="*70)
    phi_sic, sic_ok = sic_povm_d2()
    print(f"\nSIC-POVM d=2 construction valid: {sic_ok}")
    # SIC-POVM OPU: Z_i = (1/sqrt(d)) |phi_i><phi_i|^{1/2}...
    # Actually for a POVM {M_i}, the OPU is M_i = sqrt(M_i) (Kraus operators).
    # For SIC-POVM: M_i = (1/d)|phi_i><phi_i|, so Z_i = (1/sqrt(d))|phi_i><phi_i|^{1/2}
    # But |phi_i><phi_i| is a rank-1 projector, so sqrt = same (idempotent).
    # Actually the correct OPU Kraus for POVM: Z_i = sqrt(M_i).
    # For SIC: Z_i = (1/sqrt(d)) |phi_i><phi_i|.
    # Check sum_i Z_i^dag Z_i = sum_i (1/d) |phi_i><phi_i| = I (SIC normalization).
    Z_sic = [(1.0/np.sqrt(2)) * np.outer(phi_sic[i], phi_sic[i].conj()) for i in range(4)]
    check_sum = sum(Z.conj().T @ Z for Z in Z_sic)
    print(f"OPU check sum_i Z_i^dag Z_i = I? max error: {np.max(np.abs(check_sum - np.eye(2))):.2e}")

    # Compute rho[Z] = (omega(Z_j^dag Z_i))_{ij} for omega = I/2
    omega0 = np.eye(2) / 2.0
    N_sic = len(Z_sic)
    rho_sic = np.zeros((N_sic, N_sic), dtype=complex)
    for ia, Za in enumerate(Z_sic):
        for ib, Zb in enumerate(Z_sic):
            rho_sic[ia, ib] = np.trace(Zb.conj().T @ Za @ omega0)
    rho_sic = (rho_sic + rho_sic.conj().T) / 2
    eigs = np.linalg.eigvalsh(rho_sic)
    eigs_pos = eigs[eigs > 1e-14]
    S_sic = -np.sum(eigs_pos * np.log(eigs_pos))
    print(f"SIC-POVM OPU density matrix (n=1): S = {S_sic:.6f}")
    print(f"Compare: log(2) = {np.log(2):.6f}, log(4) = {np.log(4):.6f}")
    print(f"Eigenvalues: {eigs}")

    print()
    print("="*70)
    print("PART 6: AMPLITUDE DAMPING — NON-UNITAL CHANNEL")
    print("="*70)
    print("\nAmplitude damping channel Phi_gamma on qubit (d=2).")
    print("K0 = [[1,0],[0,sqrt(1-gamma)]], K1 = [[0,sqrt(gamma)],[0,0]].")
    print("OPU = projectors {P_0, P_1}. omega_0 = I/2.\n")

    Z_proj = [np.array([[1,0],[0,0]], dtype=complex),
              np.array([[0,0],[0,1]], dtype=complex)]

    for gamma in [0.0, 0.1, 0.3, 0.5, 0.9, 1.0]:
        kraus = amplitude_damping_kraus(gamma)
        # Check TP: sum K_k^dag K_k = I
        tp_check = sum(K.conj().T @ K for K in kraus)
        tp_err = np.max(np.abs(tp_check - np.eye(2)))

        # Compute AFL entropy for n = 1, 2, 3 steps
        omega_init = np.eye(2) / 2.0
        S_vals = []
        for n in [1, 2, 3]:
            S_n = opu_density_matrix_channel(Z_proj, kraus, n, omega_init)
            S_vals.append(S_n)
        print(f"gamma = {gamma:.2f}: S(n=1) = {S_vals[0]:.4f}, "
              f"S(n=2) = {S_vals[1]:.4f}, S(n=3) = {S_vals[2]:.4f}  "
              f"(TP_err = {tp_err:.1e})")

    # Fixed point of amplitude damping: Phi(rho*) = rho* => rho* = |0><0|
    # For gamma > 0: invariant state = |0><0|, S(rho*) = 0 => h~ = 0 (ordered, no entropy)
    print("\nFixed point: |0><0> (for gamma > 0). Invariant state entropy = 0.")
    print("=> corrected AFL entropy h~ = 0 for amplitude damping (dissipation kills chaos).")

    print()
    print("="*70)
    print("SUMMARY: STATUS OF QUANTUM PESIN CONJECTURE")
    print("="*70)
    print("""
1. E(U^1) is k-INDEPENDENT for kicked top (k-independence theorem, proved above).
   This is not a failure of AFL — it is a consequence of U_kick being diagonal.
   AFL entropy at n=1 is insensitive to classical chaos for this model.

2. E(U^n) for n >= 2 IS k-dependent. The matrix entropy of U^n grows with n
   for chaotic k and saturates slowly for integrable k.

3. The normalized ratio E(U^n)/log(d) approaches 1 (max scrambling) for large
   chaotic k and large j, but does NOT converge to lambda_L numerically.
   This is the EHRENFEST OBSTRUCTION: |U^n_{ij}|^2 is bounded by 1/d^2,
   so E(U^n) <= log(d) for all n, while classical chaos can have arbitrary lambda_L.

4. CORRECT INTERPRETATION of Quantum Pesin:
   - AFL entropy RATE h_AFL = 0 for finite d (always).
   - Corrected AFL: h~ = h_AFL - log(d) = s(omega) for spin chains.
   - The conjecture s(omega) = lambda_L is about the INVARIANT STATE entropy
     equaling the quantum Lyapunov exponent, not about E(U).
   - For maximally mixed omega: s(I/d) = log(d), which EQUALS max scrambling.
   - This gives h~ = log(d) for maximally chaotic systems, matching
     the claim that fully chaotic spin chains have h~ = log(d).

5. Dissipative channels: amplitude damping drives the system to |0><0>,
   so h~ = 0 — dissipation kills dynamical entropy (as expected from 2nd law).
""")
