"""
opu_comparison.py — Compare h_AFL^time for three OPU types on KI and XX chains.

OPU types:
  1. Projector OPU: {P_0, P_1} = {|0><0|, |1><1|} on site 0
  2. SIC-POVM OPU: 4-element rank-1 POVM on site 0 (d=2 SIC)
  3. Matrix-unit OPU: {|i><j|} on site 0 (4 elements for d=2)

Task 10, Step 32 — numerical study of OPU-choice dependence.
"""

import numpy as np
from scipy.linalg import expm
from itertools import product

# ─────────────────────────────────────────────────
# Pauli matrices
# ─────────────────────────────────────────────────
I2 = np.eye(2, dtype=complex)
X = np.array([[0,1],[1,0]], dtype=complex)
Y = np.array([[0,-1j],[1j,0]], dtype=complex)
Z = np.array([[1,0],[0,-1]], dtype=complex)

def kron_n(ops):
    result = ops[0]
    for op in ops[1:]:
        result = np.kron(result, op)
    return result

# ─────────────────────────────────────────────────
# Kicked Ising unitary (L sites, periodic)
# ─────────────────────────────────────────────────
def kicked_ising(L, J, g):
    D = 2**L
    # ZZ coupling: H_ZZ = J * sum_{i} Z_i Z_{i+1}
    ZZ = np.zeros((D,D), dtype=complex)
    for i in range(L):
        j = (i+1)%L
        ops = [I2]*L
        ops[i] = Z; ops[j] = Z
        ZZ += kron_n(ops)
    U_ZZ = expm(-1j * J * ZZ)
    # Transverse field: H_X = g * sum_i X_i
    HX = np.zeros((D,D), dtype=complex)
    for i in range(L):
        ops = [I2]*L; ops[i] = X
        HX += kron_n(ops)
    U_X = expm(-1j * g * HX)
    return U_X @ U_ZZ

# ─────────────────────────────────────────────────
# XX chain unitary (free fermion Floquet)
# ─────────────────────────────────────────────────
def xx_chain(L, J, dt):
    D = 2**L
    H = np.zeros((D,D), dtype=complex)
    for i in range(L-1):
        for s in [X, Y]:
            ops = [I2]*L; ops[i] = s; ops[i+1] = s
            H += 0.5*J*kron_n(ops)
    return expm(-1j * dt * H)

# ─────────────────────────────────────────────────
# Build OPUs on site 0 of L-site chain
# ─────────────────────────────────────────────────

def projector_opu(L):
    """2-element projector OPU {P_0, P_1} on site 0."""
    D = 2**L
    ops = []
    for k in range(2):
        P = np.zeros((2,2), dtype=complex)
        P[k,k] = 1.0
        rest = np.eye(2**(L-1), dtype=complex)
        ops.append(np.kron(P, rest))
    return ops

def sic_povm_opu(L):
    """4-element SIC-POVM OPU on site 0 (d=2).
    SIC vectors (Hoggar): tetrahedron on Bloch sphere.
    Z_i = (1/sqrt(2)) |phi_i><phi_i|  (rank-1, POVM condition verified)
    """
    # Tetrahedron SIC-POVM vectors for d=2
    phi = [
        np.array([1, 0], dtype=complex),
        np.array([1/np.sqrt(3), np.sqrt(2/3)], dtype=complex),
        np.array([1/np.sqrt(3), np.sqrt(2/3)*np.exp(2j*np.pi/3)], dtype=complex),
        np.array([1/np.sqrt(3), np.sqrt(2/3)*np.exp(4j*np.pi/3)], dtype=complex),
    ]
    # SIC operator: Z_i = (1/sqrt(2)) |phi_i><phi_i|
    # Check: sum Z_i^dag Z_i = (1/2) sum |phi_i><phi_i| = I  (SIC resolution of identity)
    D = 2**L
    rest = np.eye(2**(L-1), dtype=complex)
    ops = []
    for v in phi:
        P = (1/np.sqrt(2)) * np.outer(v, v.conj())
        ops.append(np.kron(P, rest))
    # Verify OPU condition
    check = sum(op.conj().T @ op for op in ops)
    err = np.max(np.abs(check - np.eye(D)))
    if err > 1e-10:
        raise ValueError(f"SIC-POVM OPU condition failed: max error {err}")
    return ops

def matrix_unit_opu(L):
    """4-element matrix-unit OPU {|i><j|} on site 0."""
    D = 2**L
    rest = np.eye(2**(L-1), dtype=complex)
    ops = []
    for i in range(2):
        for j in range(2):
            eij = np.zeros((2,2), dtype=complex)
            eij[i,j] = 1.0
            ops.append(np.kron(eij, rest))
    # Verify OPU condition: sum e_{ij}^dag e_{ij} = sum e_{ji} e_{ij} = sum delta_{ii} I = 2*I ≠ I
    # Wait — matrix-unit OPU requires normalization: Z_{ij} = (1/sqrt(d)) |i><j|
    ops = [(1/np.sqrt(2)) * op for op in ops]
    check = sum(op.conj().T @ op for op in ops)
    err = np.max(np.abs(check - np.eye(D)))
    if err > 1e-10:
        raise ValueError(f"Matrix-unit OPU condition failed: max error {err}")
    return ops

# ─────────────────────────────────────────────────
# Build time-refined OPU density matrix
# ─────────────────────────────────────────────────

def opu_density_matrix(opu_ops, U, omega, n):
    """
    Compute rho[Z^(n)] for time-AFL:
      Z^(n)_{i_1,...,i_n} = Z_{i_n}(t_{n-1}) ... Z_{i_1}(t_0)
                           = U^{n-1,dag} Z_{i_n} U^{n-1} ... Z_{i_1}
    rho_{I,J} = Tr(Z_J^dag Z_I  omega)

    Returns (k^n x k^n) density matrix where k = len(opu_ops).
    """
    k = len(opu_ops)
    kn = k**n
    D = opu_ops[0].shape[0]

    # Build all composed operators Z^(n)_I for all multi-indices I
    # I = (i_1, ..., i_n) with i_m in {0,...,k-1}
    # Z^(n)_I = U^{-(n-1)} Z_{i_n} U^{n-1} ... U^{-1} Z_{i_2} U  Z_{i_1}
    # = U^dag ... U^dag Z_{i_n} U ... U  Z_{i_2}  U  Z_{i_1}
    # Heisenberg evolved: Z_{i_m}(t_{m-1}) = (U^dag)^{m-1} Z_{i_m} U^{m-1}

    # Build Heisenberg-evolved operators at each time step
    Udag = U.conj().T
    evolved = []
    Ut = np.eye(D, dtype=complex)  # U^t
    for t in range(n):
        Utdag = Ut.conj().T  # (U^t)^dag
        evolved.append([Utdag @ Z @ Ut for Z in opu_ops])
        Ut = U @ Ut

    # Build all Z^(n)_I operators
    def build_op(multi_index):
        # Z^(n)_I = Z_{i_n}(n-1) * Z_{i_{n-1}}(n-2) * ... * Z_{i_1}(0)
        op = np.eye(D, dtype=complex)
        for t in range(n):
            op = evolved[t][multi_index[t]] @ op
        return op

    # Enumerate multi-indices
    all_indices = list(product(range(k), repeat=n))
    assert len(all_indices) == kn

    # Build density matrix
    rho = np.zeros((kn, kn), dtype=complex)
    ops_list = [build_op(idx) for idx in all_indices]

    for a, Ia in enumerate(all_indices):
        for b, Ib in enumerate(all_indices):
            rho[a, b] = np.trace(ops_list[b].conj().T @ ops_list[a] @ omega)

    # Symmetrize to enforce Hermitian (numerical noise)
    rho = (rho + rho.conj().T) / 2
    return rho

def von_neumann_entropy(rho):
    evals = np.linalg.eigvalsh(rho)
    evals = evals[evals > 1e-15]
    return float(-np.sum(evals * np.log(evals)))

def renyi2_entropy(rho):
    return float(-np.log(np.trace(rho @ rho).real))

def afl_time_entropy(opu_ops, U, omega, n):
    rho = opu_density_matrix(opu_ops, U, omega, n)
    return von_neumann_entropy(rho) / n

def renyi2_afl_time_entropy(opu_ops, U, omega, n):
    rho = opu_density_matrix(opu_ops, U, omega, n)
    return renyi2_entropy(rho) / n

# ─────────────────────────────────────────────────
# Main computation
# ─────────────────────────────────────────────────

def main():
    L = 4  # chain length (small for tractability with SIC and matrix-unit)
    D = 2**L
    omega = np.eye(D, dtype=complex) / D  # maximally mixed state

    # Dynamics
    J_ki = np.pi/4; g_ki = np.pi/4  # dual-unitary KI
    J_xx = 1.0; dt_xx = 0.5

    U_ki = kicked_ising(L, J_ki, g_ki)
    U_xx = xx_chain(L, J_xx, dt_xx)

    # OPUs
    opu_proj = projector_opu(L)
    opu_sic  = sic_povm_opu(L)
    opu_mu   = matrix_unit_opu(L)

    print(f"L={L}, D={D}")
    print(f"OPU sizes: proj={len(opu_proj)}, SIC={len(opu_sic)}, mu={len(opu_mu)}")
    print()

    # ─── h_AFL^time (von Neumann) vs n for KI and XX, three OPUs ───
    print("=" * 68)
    print("h_AFL^time / log(2)  [von Neumann AFL]")
    print(f"{'n':>3}  {'KI-proj':>10} {'KI-SIC':>10} {'KI-mu':>10}  "
          f"{'XX-proj':>10} {'XX-SIC':>10} {'XX-mu':>10}")
    print("-" * 68)

    log2 = np.log(2)
    results_ki = {name: [] for name in ['proj','sic','mu']}
    results_xx = {name: [] for name in ['proj','sic','mu']}

    for n in range(1, 6):
        h_ki_p = afl_time_entropy(opu_proj, U_ki, omega, n)
        h_ki_s = afl_time_entropy(opu_sic,  U_ki, omega, n)
        h_ki_m = afl_time_entropy(opu_mu,   U_ki, omega, n)
        h_xx_p = afl_time_entropy(opu_proj, U_xx, omega, n)
        h_xx_s = afl_time_entropy(opu_sic,  U_xx, omega, n)
        h_xx_m = afl_time_entropy(opu_mu,   U_xx, omega, n)
        print(f"{n:>3}  {h_ki_p/log2:>10.4f} {h_ki_s/log2:>10.4f} {h_ki_m/log2:>10.4f}  "
              f"{h_xx_p/log2:>10.4f} {h_xx_s/log2:>10.4f} {h_xx_m/log2:>10.4f}")
        for lst, val in zip([results_ki['proj'],results_ki['sic'],results_ki['mu']],
                            [h_ki_p, h_ki_s, h_ki_m]):
            lst.append(val/log2)
        for lst, val in zip([results_xx['proj'],results_xx['sic'],results_xx['mu']],
                            [h_xx_p, h_xx_s, h_xx_m]):
            lst.append(val/log2)

    print()
    print("OPU spread Δh = max - min over {proj, SIC, mu} at n=4:")
    n_idx = 3  # n=4 is index 3
    ki_vals = [results_ki[k][n_idx] for k in ['proj','sic','mu']]
    xx_vals = [results_xx[k][n_idx] for k in ['proj','sic','mu']]
    print(f"  KI (dual-unitary): Δh = {max(ki_vals)-min(ki_vals):.4f}")
    print(f"  XX (integrable):   Δh = {max(xx_vals)-min(xx_vals):.4f}")

    # ─── Rényi-2 AFL (GK entropy) comparison ───
    print()
    print("=" * 68)
    print("h_AFL^{(2),time} / log(2)  [Rényi-2 AFL = GK entropy]")
    print(f"{'n':>3}  {'KI-proj':>10} {'KI-SIC':>10} {'KI-mu':>10}  "
          f"{'XX-proj':>10} {'XX-SIC':>10} {'XX-mu':>10}")
    print("-" * 68)
    for n in range(1, 6):
        h_ki_p = renyi2_afl_time_entropy(opu_proj, U_ki, omega, n)
        h_ki_s = renyi2_afl_time_entropy(opu_sic,  U_ki, omega, n)
        h_ki_m = renyi2_afl_time_entropy(opu_mu,   U_ki, omega, n)
        h_xx_p = renyi2_afl_time_entropy(opu_proj, U_xx, omega, n)
        h_xx_s = renyi2_afl_time_entropy(opu_sic,  U_xx, omega, n)
        h_xx_m = renyi2_afl_time_entropy(opu_mu,   U_xx, omega, n)
        print(f"{n:>3}  {h_ki_p/log2:>10.4f} {h_ki_s/log2:>10.4f} {h_ki_m/log2:>10.4f}  "
              f"{h_xx_p/log2:>10.4f} {h_xx_s/log2:>10.4f} {h_xx_m/log2:>10.4f}")

    # ─── Verify S1 >= H2 (AFL >= GK) ───
    print()
    print("Verify S1(rho) >= H2(rho) for rho[Z^(2)] on KI with each OPU:")
    for name, opu in [('proj', opu_proj), ('SIC', opu_sic), ('mu', opu_mu)]:
        rho = opu_density_matrix(opu, U_ki, omega, 2)
        s1 = von_neumann_entropy(rho)
        h2 = renyi2_entropy(rho)
        print(f"  {name}: S1={s1/log2:.4f}, H2={h2/log2:.4f}, "
              f"S1-H2={( s1-h2)/log2:.4f} {'OK' if s1 >= h2-1e-10 else 'FAIL'}")

    # ─── Dual-unitary: all OPUs should give h = log(d) ───
    print()
    print("Dual-unitary check (KI): h_AFL / log(2) should = 1.0 for all OPUs")
    for name, opu in [('proj',opu_proj),('SIC',opu_sic),('mu',opu_mu)]:
        vals = [afl_time_entropy(opu, U_ki, omega, n)/log2 for n in range(1,5)]
        print(f"  {name}: n=1..4 -> {[f'{v:.4f}' for v in vals]}")

    # ─── Integrable chain: OPU spread shows h depends on OPU ───
    print()
    print("XX chain (integrable): OPU comparison at n=3")
    for name, opu in [('proj',opu_proj),('SIC',opu_sic),('mu',opu_mu)]:
        h = afl_time_entropy(opu, U_xx, omega, 3)/log2
        print(f"  {name}: h(n=3)={h:.4f}")

if __name__ == "__main__":
    main()
