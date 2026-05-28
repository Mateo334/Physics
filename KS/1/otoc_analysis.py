"""
OTOC analysis: verify Theorem (OTOC for projector OPU at n=1)
and compare AFL entropy (Rényi-1) with OTOC-based Rényi-2 measure.

OTOC(P_i, P_j; 1) = (2/d)|U_{ji}|^2 (1 - |U_{ji}|^2)   [Theorem 21.1]
Total OTOC C1(U)  = (2/d^2)(d - sum_{ij}|U_{ij}|^4)       [Corollary 21.2]
S(rho[Z^2])       = log d + E(U)   where E(U) = H1(q_{ij}) [Theorem 18.2]
"""

import numpy as np
from itertools import product as iproduct
from scipy.linalg import eigvalsh

# ── helpers ──────────────────────────────────────────────────────────────────

def vn_entropy(rho):
    ev = np.real(eigvalsh(rho))
    ev = ev[ev > 1e-14]
    return float(-np.sum(ev * np.log(ev)))

def opu_density_matrix(elements, omega):
    k = len(elements)
    rho = np.zeros((k, k), dtype=complex)
    for i, xi in enumerate(elements):
        for j, xj in enumerate(elements):
            rho[i, j] = np.trace(xi @ omega @ xj.conj().T)
    rho = (rho + rho.conj().T) / 2
    return rho

def time_refined_opu(base_ops, U, n):
    k = len(base_ops)
    d = base_ops[0].shape[0]
    Upow = [np.eye(d, dtype=complex)]
    for _ in range(n - 1):
        Upow.append(U @ Upow[-1])
    theta = [[Upow[t] @ base_ops[i] @ Upow[t].conj().T
              for i in range(k)] for t in range(n)]
    elements = []
    for idx in iproduct(range(k), repeat=n):
        Z = theta[n - 1][idx[n - 1]].copy()
        for t in range(n - 2, -1, -1):
            Z = Z @ theta[t][idx[t]]
        elements.append(Z)
    return elements

def otoc_projectors(U, d, omega):
    """
    Compute OTOC(P_i, P_j; 1) = (1/d) Tr([U P_i U*, P_j]† [U P_i U*, P_j])
    for all i,j and return the full d×d matrix plus the total.
    """
    otoc_mat = np.zeros((d, d))
    for i in range(d):
        Pi = np.zeros((d, d), dtype=complex)
        Pi[i, i] = 1.0
        UPiU = U @ Pi @ U.conj().T  # rank-1 projector
        for j in range(d):
            Pj = np.zeros((d, d), dtype=complex)
            Pj[j, j] = 1.0
            comm = UPiU @ Pj - Pj @ UPiU
            C = comm.conj().T @ comm
            otoc_mat[i, j] = np.real(np.trace(C)) / d
    return otoc_mat

def matrix_entropy_EU(U, d):
    """E(U) = -(1/d) sum_{ij} |U_{ij}|^2 log|U_{ij}|^2"""
    Umod2 = np.abs(U) ** 2
    mask = Umod2 > 1e-15
    return -float(np.sum(Umod2[mask] * np.log(Umod2[mask]))) / d

def renyi2_entropy(U, d):
    """H_2(q) = -log( sum_{ij} |U_{ij}|^4 / d^2 )"""
    s = np.sum(np.abs(U) ** 4)
    return -np.log(s / d**2)

def total_otoc_formula(U, d):
    """C1(U) = (1/d^2) sum_{i,j} OTOC(P_i,P_j;1) = (2/d^3)(d - sum|U_{ij}|^4)"""
    s4 = np.sum(np.abs(U) ** 4)
    return 2.0 * (d - s4) / d**3

def Utheta(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=complex)

# ── Gate library ─────────────────────────────────────────────────────────────

H_gate = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
T_gate = np.array([[1, 0], [0, np.exp(1j*np.pi/4)]], dtype=complex)
S_gate = np.array([[1, 0], [0, 1j]], dtype=complex)
X_gate = np.array([[0, 1], [1, 0]], dtype=complex)
I_gate = np.eye(2, dtype=complex)

gates_d2 = [
    ("Identity",  I_gate),
    ("Pauli X",   X_gate),
    ("T gate",    T_gate),
    ("S gate",    S_gate),
    ("Hadamard",  H_gate),
    ("U(pi/4)",   Utheta(np.pi/4)),
    ("U(pi/8)",   Utheta(np.pi/8)),
    ("U(pi/3)",   Utheta(np.pi/3)),
]

d = 2
omega2 = np.eye(d, dtype=complex) / d
proj_opu = [np.outer(np.eye(d)[i], np.eye(d)[i]).astype(complex) for i in range(d)]

# ── Table 5: OTOC formula verification ───────────────────────────────────────
print("=" * 75)
print("Table 5: OTOC(P_i, P_j; 1) formula verification")
print("Formula: OTOC(P_i,P_j;1) = (2/d)|U_{ji}|^2 (1-|U_{ji}|^2)")
print()
for name, U in gates_d2:
    otoc_mat_num = otoc_projectors(U, d, omega2)
    # Formula value
    Umod2 = np.abs(U)**2
    otoc_formula = (2.0/d) * Umod2.T * (1 - Umod2.T)
    max_err = np.max(np.abs(otoc_mat_num - otoc_formula))
    print(f"  {name:<14}: max|OTOC_num - OTOC_formula| = {max_err:.2e}")

# ── Table 6: AFL vs OTOC comparison ──────────────────────────────────────────
print()
print("=" * 90)
print("Table 6: AFL entropy, H1(q), H2(q), total OTOC C1(U), MUB check for d=2")
print("  q_{ij} = |U_{ij}|^2/d,  H1(q)=E(U)+log(d),  H2(q)=-log(sum q_ij^2)")
print(f"{'Gate':<14}  {'S(Z^2)':>8}  {'H1(q)':>8}  {'H2(q)':>8}  {'C1_num':>8}  {'C1_form':>8}  {'H1-H2':>7}  {'MUB?':>5}")
print("-" * 90)
for name, U in gates_d2:
    EU = matrix_entropy_EU(U, d)
    H1 = EU + np.log(d)   # H1(q) = E(U) + log d
    H2 = renyi2_entropy(U, d)
    otoc_num = np.sum(otoc_projectors(U, d, omega2)) / d**2
    C1_formula = total_otoc_formula(U, d)
    is_mub = abs(EU - np.log(d)) < 1e-6
    S2 = np.log(d) + EU
    print(f"{name:<14}  {S2:>8.5f}  {H1:>8.5f}  {H2:>8.5f}  {otoc_num:>8.6f}  {C1_formula:>8.6f}  {H1-H2:>7.5f}  {'YES' if is_mub else 'no':>5}")

# ── Table 7: Rényi comparison vs rotation angle ───────────────────────────────
print()
print("=" * 75)
print("Table 7: H1(q)=S(Z^2), H2(q), C1(U) for U(theta), d=2")
print(f"{'theta/pi':<10}  {'S(Z^2)=H1':>12}  {'H2(q)':>10}  {'C1(U)':>10}  {'H1-H2':>10}")
print("-" * 75)
for frac in [0, 1/8, 1/6, 1/4, 1/3, 3/8, 1/2]:
    theta = frac * np.pi
    U = Utheta(theta)
    EU = matrix_entropy_EU(U, d)
    H1 = EU + np.log(d)
    H2 = renyi2_entropy(U, d)
    C1 = total_otoc_formula(U, d)
    print(f"{frac:<10.4f}  {H1:>12.6f}  {H2:>10.6f}  {C1:>10.6f}  {H1-H2:>10.6f}")

# ── Check: maximum of total OTOC and E(U) for d=2 ────────────────────────────
print()
print("Maximum values (d=2):")
print(f"  max S(Z^2)   = 2*log(d)   = {2*np.log(d):.6f}  (Hadamard/MUB)")
print(f"  max C1(U)   = 2(d-1)/d^3 = {2*(d-1)/d**3:.6f}  (Hadamard/MUB)")
print(f"  max H2(q)   = log(d^2)   = {np.log(d**2):.6f}  (Hadamard/MUB)")
print(f"  min S(Z^2)  = log(d)     = {np.log(d):.6f}  (permutation)")
print(f"  min C1(U)   = 0          (permutation/monomial)")

# ── Verify: OTOC for Hadamard, all (i,j) pairs ───────────────────────────────
print()
print("OTOC matrix for Hadamard gate, d=2:")
otoc_had = otoc_projectors(H_gate, d, omega2)
print("  Numerical:")
for i in range(d):
    for j in range(d):
        print(f"    OTOC(P_{i},P_{j};1) = {otoc_had[i,j]:.6f}", end="")
    print()
print(f"  Formula (2/d)|U_ji|^2(1-|U_ji|^2): all = {2/d * 0.5 * 0.5:.6f} (Hadamard |U_ij|^2=1/2)")

# ── Summary table for LaTeX ───────────────────────────────────────────────────
print()
print("=" * 75)
print("Summary for LaTeX Table 8: S(rho[Z^2]), E(U), H2, C1(U) for d=2")
print(f"{'Gate':<14}  {'S(Z^2)':>8}  {'E(U)':>8}  {'H2':>8}  {'C1':>8}")
print("-" * 75)
for name, U in gates_d2:
    elems2 = time_refined_opu(proj_opu, U, 2)
    rho2 = opu_density_matrix(elems2, omega2)
    S2 = vn_entropy(rho2)
    EU = matrix_entropy_EU(U, d)
    H2 = renyi2_entropy(U, d)
    C1 = total_otoc_formula(U, d)
    print(f"{name:<14}  {S2:>8.5f}  {EU:>8.5f}  {H2:>8.5f}  {C1:>8.6f}")
