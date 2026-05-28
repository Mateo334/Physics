"""
Observational entropy (Šafránek 2019, 2021) vs AFL entropy.
Computes both for several OPU/POVM choices on a qubit, verifies all identities.
"""
import numpy as np
from scipy.linalg import logm

def von_neumann_entropy(rho, base=np.e):
    """S(rho) = -Tr(rho log rho), with 0*log0 = 0."""
    eigs = np.linalg.eigvalsh(rho)
    eigs = eigs[eigs > 1e-15]
    return float(-np.sum(eigs * np.log(eigs)) / np.log(base))

def shannon_entropy(p, base=np.e):
    p = np.array(p, dtype=float)
    p = p[p > 1e-15]
    return float(-np.sum(p * np.log(p)) / np.log(base))

def afl_density_matrix(ops, omega):
    """
    Build AFL density matrix rho[Z] with entries rho_ij = Tr(Z_i omega Z_j^dagger).
    ops: list of d×d matrices Z_i (OPU operators)
    omega: d×d density matrix (state)
    """
    k = len(ops)
    rho = np.zeros((k, k), dtype=complex)
    for i, Zi in enumerate(ops):
        for j, Zj in enumerate(ops):
            rho[i, j] = np.trace(Zi @ omega @ Zj.conj().T)
    return rho

def afl_density_matrix_refined(ops, omega, U, n):
    """
    Build AFL density matrix rho[Z^(n)] for time-refinement up to step n.
    ops: base OPU {Z_i}, dynamics Theta(A) = U A U^dagger.
    Returns the k^n × k^n density matrix.
    """
    k = len(ops)
    d = ops[0].shape[0]
    # Build all k^n refined operators Z_{i_0,...,i_{n-1}}
    # Z_{i0...i_{n-1}} = U^{n-1} Z_{i_{n-1}} U^{*n-1} ... U Z_{i1} U^* Z_{i0}
    refined = []
    from itertools import product as iproduct
    for idx in iproduct(range(k), repeat=n):
        # idx = (i_{n-1}, ..., i_1, i_0) — we build right to left
        op = ops[idx[0]]  # Z_{i_0}
        for t in range(1, n):
            Ut = np.linalg.matrix_power(U, t)
            op = Ut @ ops[idx[t]] @ Ut.conj().T @ op
        refined.append(op)
    # Now compute density matrix
    kn = k**n
    rho = np.zeros((kn, kn), dtype=complex)
    for i, Zi in enumerate(refined):
        for j, Zj in enumerate(refined):
            rho[i, j] = np.trace(Zi @ omega @ Zj.conj().T)
    return rho

def observational_entropy(povm, omega):
    """
    Šafránek observational entropy:
      S_obs(omega, M) = -sum_m p_m log(p_m / V_m)
    where p_m = Tr(M_m omega), V_m = Tr(M_m) / d.
    povm: list of d×d positive semidefinite matrices summing to identity
    omega: d×d density matrix
    Returns S_obs (in nats)
    """
    d = omega.shape[0]
    S = 0.0
    for Mm in povm:
        pm = float(np.real(np.trace(Mm @ omega)))
        Vm = float(np.real(np.trace(Mm))) / d
        if pm > 1e-15 and Vm > 1e-15:
            S -= pm * np.log(pm / Vm)
    return S

def observational_entropy_refined(povm, omega, U, n):
    """
    Iterative observational entropy for n-step POVM:
      M^(n)_{m_0,...,m_{n-1}} = M_{m_0} * U M_{m_1} U^* * ... * U^{n-1} M_{m_{n-1}} U^{*(n-1)}
    (Šafránek sequential measurement)
    """
    from itertools import product as iproduct
    k = len(povm)
    d = omega.shape[0]
    povm_refined = []
    for idx in iproduct(range(k), repeat=n):
        op = np.eye(d, dtype=complex)
        for t in range(n):
            Ut = np.linalg.matrix_power(U, t)
            op = op @ (Ut @ povm[idx[t]] @ Ut.conj().T)
        povm_refined.append(op)
    return observational_entropy(povm_refined, omega)

# ============================================================
# Setup: qubit, Hadamard dynamics, projector OPU
# ============================================================
d = 2
H = (1/np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
omega = np.eye(d, dtype=complex) / d  # maximally mixed

# Projector OPU: Z_i = |e_i><e_i|
P0 = np.array([[1,0],[0,0]], dtype=complex)
P1 = np.array([[0,0],[0,1]], dtype=complex)
ops_proj = [P0, P1]

# Matrix-unit OPU: Z_{ij} = (1/sqrt(d)) |e_i><e_j|  (d^2 operators)
ops_mu = [np.outer(np.eye(d)[i], np.eye(d)[j]) / np.sqrt(d)
          for i in range(d) for j in range(d)]

print("=" * 60)
print("AFL entropy vs n for projector OPU, Hadamard dynamics")
print("=" * 60)
for n in range(1, 7):
    rho_n = afl_density_matrix_refined(ops_proj, omega, H, n)
    S = von_neumann_entropy(rho_n)
    rank = np.linalg.matrix_rank(rho_n, tol=1e-10)
    print(f"  n={n}: S = {S:.6f}  (log2={np.log(2):.6f}, log4={np.log(4):.6f}), rank={rank}")

print()
print("=" * 60)
print("AFL entropy vs n for matrix-unit OPU, Hadamard dynamics")
print("=" * 60)
for n in range(1, 5):
    rho_n = afl_density_matrix_refined(ops_mu, omega, H, n)
    S = von_neumann_entropy(rho_n)
    rank = np.linalg.matrix_rank(rho_n, tol=1e-10)
    print(f"  n={n}: S = {S:.6f}, rank={rank}")

# ============================================================
# Observational entropy
# ============================================================
print()
print("=" * 60)
print("Observational entropy (Šafránek) vs AFL, projector OPU")
print("=" * 60)
# POVM = projector OPU (M_i = Z_i^dagger Z_i = Z_i for projectors)
povm_proj = [P0, P1]
for n in range(1, 5):
    S_obs = observational_entropy_refined(povm_proj, omega, H, n)
    rho_n = afl_density_matrix_refined(ops_proj, omega, H, n)
    S_afl = von_neumann_entropy(rho_n)
    print(f"  n={n}: S_obs = {S_obs:.6f}, S_AFL = {S_afl:.6f}, diff = {S_afl - S_obs:.6f}")

# ============================================================
# Coarse-graining refinement: varying partition size k
# ============================================================
print()
print("=" * 60)
print("AFL entropy at n=2 for projector OPU as d grows (qudit)")
print("=" * 60)
for d_val in [2, 3, 4, 5, 6]:
    omega_d = np.eye(d_val, dtype=complex) / d_val
    ops_d = [np.outer(np.eye(d_val)[i], np.eye(d_val)[i]) for i in range(d_val)]
    # Hadamard-like unitary: DFT matrix
    DFT = np.array([[np.exp(2j*np.pi*i*j/d_val)/np.sqrt(d_val)
                     for j in range(d_val)] for i in range(d_val)], dtype=complex)
    rho2 = afl_density_matrix_refined(ops_d, omega_d, DFT, 2)
    S2 = von_neumann_entropy(rho2)
    print(f"  d={d_val}: S(rho[Z^2]) = {S2:.6f}, log(d)={np.log(d_val):.6f}, "
          f"2log(d)={2*np.log(d_val):.6f}")

# ============================================================
# AFL vs observational entropy: general POVM, qubit
# ============================================================
print()
print("=" * 60)
print("AFL vs observational entropy for various qubit states and OPUs")
print("=" * 60)
H_gate = (1/np.sqrt(2)) * np.array([[1,1],[1,-1]], dtype=complex)
for theta in [0, np.pi/8, np.pi/4, 3*np.pi/8, np.pi/2]:
    U = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]], dtype=complex)
    rho1 = afl_density_matrix(ops_proj, omega)
    S_afl1 = von_neumann_entropy(rho1)
    rho2 = afl_density_matrix_refined(ops_proj, omega, U, 2)
    S_afl2 = von_neumann_entropy(rho2)
    # Observational: POVM = projectors
    S_obs1 = observational_entropy(povm_proj, omega)
    S_obs_seq = observational_entropy_refined(povm_proj, omega, U, 2)
    print(f"  theta={theta:.4f}: S_AFL(n=2)={S_afl2:.4f}, S_obs_seq(n=2)={S_obs_seq:.4f}, "
          f"S_AFL-S_obs={S_afl2-S_obs_seq:.4f}")

# ============================================================
# Key identity: for projector OPU and maximally mixed state,
# AFL(n=1) diagonal entries = p_i = 1/d, V_i = 1/d, S_obs = log(d)
# ============================================================
print()
print("=" * 60)
print("Verify: for projector OPU + maximally mixed, S_obs = log(d)")
print("=" * 60)
for d_val in [2, 3, 4]:
    omega_d = np.eye(d_val, dtype=complex) / d_val
    povm_d = [np.outer(np.eye(d_val)[i], np.eye(d_val)[i]) for i in range(d_val)]
    S_obs = observational_entropy(povm_d, omega_d)
    print(f"  d={d_val}: S_obs = {S_obs:.6f}, log(d) = {np.log(d_val):.6f}, match={abs(S_obs - np.log(d_val)) < 1e-10}")

# ============================================================
# Quantum excess entropy: S_AFL(n) - S_obs(n)
# This measures off-diagonal coherences in rho[Z^n]
# ============================================================
print()
print("=" * 60)
print("Quantum excess entropy: S_AFL(n) - S_obs(n) for Hadamard qubit")
print("=" * 60)
for n in range(1, 5):
    rho_n = afl_density_matrix_refined(ops_proj, omega, H, n)
    S_afl = von_neumann_entropy(rho_n)
    S_obs = observational_entropy_refined(povm_proj, omega, H, n)
    # Classical entropy = diagonal of rho
    diag = np.real(np.diag(rho_n))
    S_class = shannon_entropy(diag)
    print(f"  n={n}: S_AFL={S_afl:.4f}, S_obs={S_obs:.4f}, S_diag={S_class:.4f}, "
          f"excess(AFL-diag)={S_afl - S_class:.4f}")
