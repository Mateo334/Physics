"""
AFL entropy computation for finite-dimensional quantum systems.
Mateo's framework: rho[X]_{ij} = Tr(x_i omega x_j*)
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
    """rho[X]_{ij} = Tr(x_i omega x_j*)"""
    k = len(elements)
    rho = np.zeros((k, k), dtype=complex)
    for i, xi in enumerate(elements):
        for j, xj in enumerate(elements):
            rho[i, j] = np.trace(xi @ omega @ xj.conj().T)
    # symmetrise tiny numerical errors
    rho = (rho + rho.conj().T) / 2
    return rho

def time_refined_opu(base_ops, U, n):
    """
    Z^(n)_{(i0,i1,...,i_{n-1})} = Theta^{n-1}(Z_{i_{n-1}}) @ ... @ Z_{i0}
    iproduct index tuple idx = (i0, i1, ..., i_{n-1})
    """
    k = len(base_ops)
    d = base_ops[0].shape[0]
    Upow = [np.eye(d, dtype=complex)]
    for _ in range(n - 1):
        Upow.append(U @ Upow[-1])
    theta = [[Upow[t] @ base_ops[i] @ Upow[t].conj().T
              for i in range(k)] for t in range(n)]
    elements = []
    for idx in iproduct(range(k), repeat=n):
        # idx[t] = i_t; product from left = Theta^{n-1}(Z_{i_{n-1}}) @ ... @ Z_{i0}
        Z = theta[n - 1][idx[n - 1]].copy()
        for t in range(n - 2, -1, -1):
            Z = Z @ theta[t][idx[t]]
        elements.append(Z)
    return elements

def entropy_sequence(base_ops, U, omega, max_n):
    results = []
    for n in range(1, max_n + 1):
        elems = time_refined_opu(base_ops, U, n)
        rho = opu_density_matrix(elems, omega)
        S = vn_entropy(rho)
        rank = int(np.sum(np.real(eigvalsh(rho)) > 1e-10))
        results.append((n, S, S / n, rank))
    return results

# ── systems ───────────────────────────────────────────────────────────────────

# Hadamard gate
H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)

# T gate
T = np.array([[1, 0], [0, np.exp(1j * np.pi / 4)]], dtype=complex)

# S gate (90-degree phase)
S_gate = np.array([[1, 0], [0, 1j]], dtype=complex)

# Rotation by angle theta in xz plane: U(theta)
def Utheta(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=complex)

# Projector OPU for qubit
P0 = np.array([[1, 0], [0, 0]], dtype=complex)
P1 = np.array([[0, 0], [0, 1]], dtype=complex)
proj_opu = [P0, P1]

# Matrix-unit OPU for qubit: e_{ij}/sqrt(2)
mat_units = [np.outer([1, 0], [1, 0]) / np.sqrt(2),
             np.outer([1, 0], [0, 1]) / np.sqrt(2),
             np.outer([0, 1], [1, 0]) / np.sqrt(2),
             np.outer([0, 1], [0, 1]) / np.sqrt(2)]

omega2 = np.eye(2, dtype=complex) / 2  # maximally mixed qubit state

# ── Table 1: Hadamard + projector OPU ────────────────────────────────────────
print("=" * 65)
print("Table 1: Hadamard qubit, projector OPU, omega = I/2")
print(f"{'n':>3}  {'k^n':>5}  {'S(rho[Z^n])':>14}  {'S/n':>10}  {'rank':>6}")
print("-" * 65)
for n, S, Sn, rank in entropy_sequence(proj_opu, H, omega2, 6):
    kn = 2 ** n
    print(f"{n:>3}  {kn:>5}  {S:>14.6f}  {Sn:>10.6f}  {rank:>6}")

# ── Table 2: Hadamard + matrix-unit OPU ──────────────────────────────────────
print()
print("=" * 65)
print("Table 2: Hadamard qubit, matrix-unit OPU, omega = I/2")
print(f"{'n':>3}  {'k^n':>5}  {'S(rho[Z^n])':>14}  {'S/n':>10}  {'rank':>6}")
print("-" * 65)
for n, S, Sn, rank in entropy_sequence(mat_units, H, omega2, 4):
    kn = 4 ** n
    print(f"{n:>3}  {kn:>5}  {S:>14.6f}  {Sn:>10.6f}  {rank:>6}")

# ── Table 3: Various gates, projector OPU, n=2 ───────────────────────────────
print()
print("=" * 65)
print("Table 3: S(rho[Z^2]) for various unitaries (projector OPU)")
print(f"{'Gate':>14}  {'S(n=2)':>10}  {'S(n=2)/log(d^2)':>18}  {'MUB?':>6}")
print("-" * 65)
log4 = np.log(4)
gates = [("Identity", np.eye(2, dtype=complex)),
         ("Pauli X",  np.array([[0,1],[1,0]], dtype=complex)),
         ("T gate",   T),
         ("S gate",   S_gate),
         ("Hadamard", H),
         ("U(pi/4)",  Utheta(np.pi / 4)),
         ("U(pi/8)",  Utheta(np.pi / 8)),
         ("U(pi/3)",  Utheta(np.pi / 3))]
for name, U in gates:
    elems = time_refined_opu(proj_opu, U, 2)
    rho = opu_density_matrix(elems, omega2)
    S = vn_entropy(rho)
    ratio = S / log4
    is_mub = abs(ratio - 1.0) < 1e-6
    print(f"{name:>14}  {S:>10.6f}  {ratio:>18.6f}  {'YES' if is_mub else 'no':>6}")

# ── Table 4: S(rho[Z^2]) vs rotation angle theta ─────────────────────────────
print()
print("=" * 65)
print("Table 4: S(rho[Z^2]) vs rotation angle for U(theta)")
print(f"{'theta/pi':>10}  {'S(n=2)':>10}  {'E(U)':>10}")
print("-" * 65)
thetas = [0, 1/8, 1/6, 1/4, 1/3, 3/8, 1/2]
for frac in thetas:
    theta = frac * np.pi
    U = Utheta(theta)
    elems = time_refined_opu(proj_opu, U, 2)
    rho = opu_density_matrix(elems, omega2)
    S = vn_entropy(rho)
    # matrix entropy E(U) = -(1/d) sum_{ij} |U_{ij}|^2 log|U_{ij}|^2
    Umod2 = np.abs(U) ** 2
    E_U = -float(np.sum(Umod2[Umod2 > 1e-14] * np.log(Umod2[Umod2 > 1e-14]))) / 2
    print(f"{frac:>10.4f}  {S:>10.6f}  {E_U:>10.6f}")

# ── Verify OPU property ───────────────────────────────────────────────────────
print()
print("OPU property check for proj OPU:")
print("  sum Z_i* Z_i =", np.round(sum(z.conj().T @ z for z in proj_opu), 6))

print()
print("OPU property check for matrix-unit OPU:")
print("  sum Z_i* Z_i =", np.round(sum(z.conj().T @ z for z in mat_units), 6))

# ── Verify rank formula ───────────────────────────────────────────────────────
print()
print("Rank of rho[Z^(n)] for Hadamard (max = d^2 = 4):")
for n, S, _, rank in entropy_sequence(proj_opu, H, omega2, 6):
    print(f"  n={n}: rank={rank}, S={S:.4f}, S/log(d^2)={S/np.log(4):.4f}")

# ── Verify matrix entropy formula S = log d + E(U) ───────────────────────────
print()
print("Verify: S(rho[Z^2]) = log(d) + E(U) for various U")
for name, U in gates:
    elems = time_refined_opu(proj_opu, U, 2)
    rho = opu_density_matrix(elems, omega2)
    S_num = vn_entropy(rho)
    Umod2 = np.abs(U) ** 2
    E_U = -float(np.sum(Umod2[Umod2 > 1e-14] * np.log(Umod2[Umod2 > 1e-14]))) / 2
    S_formula = np.log(2) + E_U
    print(f"  {name:>14}: S_num={S_num:.6f}, log2+E(U)={S_formula:.6f}, diff={abs(S_num-S_formula):.2e}")
