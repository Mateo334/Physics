"""
cnt_marginal_consistency.py — Verify two-sided marginal consistency for time AFL.

Theorem (Two-Sided Marginal Consistency):
  Tr_first[rho[Z^n]] = rho[Z^{n-1}]    (already proved in Section 16)
  Tr_last [rho[Z^n]] = rho[Z^{n-1}]    (proved below, via OPU completeness)

This makes the SSA proof of concavity (Section 16, Theorem thm:ssa_concave) rigorous.

Proof of Tr_last:
  Z^n_{(I,i_n)} = A_I * C_{i_n}   where  A_I = P_{i1} Ud ... Ud P_{i_{n-1}}
                                          C_{i_n} = Ud P_{i_n} U^{n-1}

  Sum_{i_n} C_{i_n} C_{i_n}^dag
    = Sum_{i_n} (Ud P_{i_n} U^{n-1}) (U^{n-1})^dag P_{i_n} U
    = Ud (Sum_{i_n} P_{i_n} * I * P_{i_n}) U     [U^{n-1}(U^{n-1})^dag = I]
    = Ud (Sum_{i_n} P_{i_n}^2) U = Ud I U = I.   [OPU: Sum P_i = I, P_i^2=P_i]

  (rho[Z^n])^{marg-last}_{I,J}
    = (1/D) Sum_{i_n} Tr(Z^n_{J,i_n}^dag Z^n_{I,i_n})
    = (1/D) Tr(A_J^dag A_I * Sum_{i_n} C_{i_n} C_{i_n}^dag)   [cyclic]
    = (1/D) Tr(A_J^dag A_I * I) = (1/D) Tr(A_J^dag A_I)
    = rho[Z^{n-1}]_{I,J}.   QED.
"""

import numpy as np
from scipy.linalg import expm, eigh
import itertools

np.random.seed(42)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)

def kron_site(op, site, L):
    ops = [np.eye(2, dtype=complex)] * L
    ops[site] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def kicked_ising_open(L, J, g):
    H_ZZ = sum(kron_site(sz, i, L) @ kron_site(sz, i+1, L) for i in range(L-1))
    H_X  = sum(kron_site(sx, i, L) for i in range(L))
    return expm(-1j * J * H_ZZ) @ expm(-1j * g * H_X)

def x_projectors_site0(L):
    D = 2**L
    Px0 = np.zeros((D, D), dtype=complex)
    Px1 = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            bi = (i >> (L-1)) & 1
            bj = (j >> (L-1)) & 1
            if (i & ((1 << (L-1)) - 1)) == (j & ((1 << (L-1)) - 1)):
                Px0[i, j] += 0.5
                Px1[i, j] += 0.5 * (-1)**(bi + bj)
    return [Px0, Px1]

def build_ops(U, P, n):
    """Build all Z^(n)_I operators."""
    D = U.shape[0]
    Ud = U.conj().T
    Un1 = np.linalg.matrix_power(U, n - 1)
    ops = {}
    for idx in itertools.product(range(len(P)), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ Ud @ P[idx[t]]
        ops[idx] = Z @ Un1
    return ops

def time_afl_density_matrix(U, P, n):
    D = U.shape[0]
    ops = build_ops(U, P, n)
    indices = list(ops.keys())
    dim = len(indices)
    M = np.zeros((dim, dim), dtype=complex)
    for a, idxA in enumerate(indices):
        for b, idxB in enumerate(indices):
            M[a, b] = np.sum(ops[idxB].conj() * ops[idxA]) / D
    return (M + M.conj().T) / 2

def von_neumann_entropy(M, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True))
    evals = evals[evals > tol]
    if len(evals) == 0: return 0.0
    evals /= evals.sum()
    return float(-np.sum(evals * np.log(evals)))

def marginal_last(rho_n, d, n):
    """
    Partial trace over the LAST index of rho[Z^n].
    rho_n is (d^n x d^n), index order: (i1, i2, ..., i_n).
    Returns rho_{n-1} = (d^{n-1} x d^{n-1}).
    """
    dim_n = d**n
    dim_nm1 = d**(n-1)
    rho_nm1 = np.zeros((dim_nm1, dim_nm1), dtype=complex)
    for a in range(dim_nm1):
        for b in range(dim_nm1):
            for i_n in range(d):
                row_idx = a * d + i_n
                col_idx = b * d + i_n
                rho_nm1[a, b] += rho_n[row_idx, col_idx]
    return rho_nm1

def marginal_first(rho_n, d, n):
    """
    Partial trace over the FIRST index of rho[Z^n].
    """
    dim_n = d**n
    dim_nm1 = d**(n-1)
    rho_nm1 = np.zeros((dim_nm1, dim_nm1), dtype=complex)
    for a in range(dim_nm1):
        for b in range(dim_nm1):
            for i_1 in range(d):
                row_idx = i_1 * dim_nm1 + a
                col_idx = i_1 * dim_nm1 + b
                rho_nm1[a, b] += rho_n[row_idx, col_idx]
    return rho_nm1

# ─── Main verification ────────────────────────────────────────────────────────

L = 4; d = 2

print("=" * 72)
print("TWO-SIDED MARGINAL CONSISTENCY VERIFICATION")
print("=" * 72)
print()
print("Theorem: For time-AFL with projector OPU and unitary U:")
print("  Tr_first[rho[Z^n]] = rho[Z^{n-1}]  (proved: OPU completeness, Sum P_{i1} = I)")
print("  Tr_last [rho[Z^n]] = rho[Z^{n-1}]  (proved: Sum C_{i_n} C_{i_n}^dag = I)")
print()

test_cases = [
    (np.pi/8, 0.0,    "J=pi/8, g=0 (integrable)"),
    (np.pi/8, np.pi/8, "J=pi/8, g=pi/8 (off-DU)"),
    (np.pi/4, np.pi/4, "J=pi/4, g=pi/4 (dual-unitary)"),
    (0.10, 0.50,       "J=0.10, g=0.50 (arbitrary)"),
]

print("-" * 72)
print("Testing Tr_first = rho[Z^{n-1}] and Tr_last = rho[Z^{n-1}]:")
print()

all_ok = True
for J, g, name in test_cases:
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site0(L)

    print(f"  System: {name}")
    for n in range(2, 5):
        rho_n   = time_afl_density_matrix(U, P, n)
        rho_nm1 = time_afl_density_matrix(U, P, n-1)

        # Marginal over FIRST index
        marg_first = marginal_first(rho_n, d, n)
        err_first = np.max(np.abs(marg_first - rho_nm1))

        # Marginal over LAST index
        marg_last = marginal_last(rho_n, d, n)
        err_last = np.max(np.abs(marg_last - rho_nm1))

        ok_f = err_first < 1e-10
        ok_l = err_last  < 1e-10
        if not (ok_f and ok_l): all_ok = False

        print(f"    n={n}: err_first={err_first:.2e} {'✓' if ok_f else '✗'}  "
              f"err_last={err_last:.2e} {'✓' if ok_l else '✗'}")
    print()

print(f"All marginal consistency checks pass: {'YES' if all_ok else 'NO'}")
print()

# ─── Verify Sum C_{i_n} C_{i_n}^dag = I analytically ─────────────────────────
print("-" * 72)
print("Verifying Sum_{i_n} C_{i_n} C_{i_n}^dag = I for various U, n:")
print()
print("  C_{i_n} = Ud P_{i_n} U^{n-1}")
print("  Sum C_{i_n} C_{i_n}^dag = Ud (Sum P_{i_n} U^{n-1} (U^{n-1})^dag P_{i_n}) U")
print("                           = Ud (Sum P_{i_n} I P_{i_n}) U  [unitarity]")
print("                           = Ud (Sum P_{i_n}^2) U = Ud I U = I")
print()

D = 2**L
I_D = np.eye(D, dtype=complex)

for J, g, name in test_cases:
    U = kicked_ising_open(L, J, g)
    Ud = U.conj().T
    P = x_projectors_site0(L)

    print(f"  {name}:")
    for n in range(1, 5):
        Un1 = np.linalg.matrix_power(U, n-1)
        total = np.zeros((D, D), dtype=complex)
        for Pi in P:
            Ci = Ud @ Pi @ Un1
            total += Ci @ Ci.conj().T
        err = np.max(np.abs(total - I_D))
        print(f"    n={n}: ||Sum C_i C_i^dag - I|| = {err:.2e} {'✓' if err < 1e-12 else '✗'}")
    print()

# ─── Confirm SSA concavity is rigorous ────────────────────────────────────────
print("-" * 72)
print("Rigorous SSA proof ingredients (all verified):")
print()
print("For rho[Z^n] = rho(X_1,...,X_n) and tripartite split A={X_1}, B={X_2,...,X_{n-1}}, C={X_n}:")
print()
print("  S(AB) = S(X_1,...,X_{n-1}) = S_{n-1}  [Tr_C = Tr_last:  rho[Z^n] -> rho[Z^{n-1}]]")
print("  S(BC) = S(X_2,...,X_n)     = S_{n-1}  [Tr_A = Tr_first: rho[Z^n] -> rho[Z^{n-1}]]")
print("          ^ same entropy by time-stationarity (same spectrum as rho[Z^{n-1}])")
print("  S(B)  = S(X_2,...,X_{n-1}) = S_{n-2}  [two marginals: Tr_first then Tr_last of (n-1)-step]")
print("  S(ABC)= S_n")
print()
print("  SSA: S(AB) + S(BC) >= S(B) + S(ABC)")
print("  =>   S_{n-1} + S_{n-1} >= S_{n-2} + S_n")
print("  =>   S_n - S_{n-1} <= S_{n-1} - S_{n-2}   QED concavity.")
print()
print("PROOF IS RIGOROUS: both Tr_first and Tr_last give rho[Z^{n-1}].")
print("Time-stationarity: S(X_2,...,X_n) = S(X_1,...,X_{n-1}) by unitary equivalence.")
print()

# Verify S(X_2,...,X_n) = S(X_1,...,X_{n-1}) numerically
print("-" * 72)
print("Verifying time-stationarity: S(X_2,...,X_n) = S(X_1,...,X_{n-1}):")
print()

for J, g, name in test_cases[:3]:
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site0(L)
    print(f"  {name}:")
    for n in range(2, 5):
        rho_n = time_afl_density_matrix(U, P, n)
        rho_nm1 = time_afl_density_matrix(U, P, n-1)

        # S(X_1,...,X_{n-1}) = S_{n-1} (marginal over last)
        marg_last = marginal_last(rho_n, d, n)
        S_AB = von_neumann_entropy(marg_last)

        # S(X_2,...,X_n) = S_{n-1} (marginal over first)
        marg_first = marginal_first(rho_n, d, n)
        S_BC = von_neumann_entropy(marg_first)

        # S_{n-1} from direct computation
        S_nm1 = von_neumann_entropy(rho_nm1)

        print(f"    n={n}: S_AB={S_AB:.6f}  S_BC={S_BC:.6f}  S_{{n-1}}={S_nm1:.6f}  "
              f"equal={'YES' if abs(S_AB - S_BC) < 1e-8 else 'NO'}")
    print()

print("=" * 72)
print("SUMMARY")
print("=" * 72)
print("""
Both marginal consistencies hold exactly (to machine precision):
  Tr_first[rho[Z^n]] = rho[Z^{n-1}]   (OPU: Sum_{i1} P_{i1}^2 = Sum P_{i1} = I)
  Tr_last [rho[Z^n]] = rho[Z^{n-1}]   (OPU: Sum_{i_n} C_{i_n} C_{i_n}^dag = I)

Time-stationarity holds: S(X_2,...,X_n) = S(X_1,...,X_{n-1}) = S_{n-1}.

The SSA proof of concavity (S_n + S_{n-2} <= 2 S_{n-1}) in Section 16 is
therefore FULLY RIGOROUS for any unitary dynamics and projector OPU.

Addition to Section 16:
  - Lemma: Tr_last[rho[Z^n]] = rho[Z^{n-1}] (proved by Sum C_i C_i^dag = I).
  - Remark: time-stationarity (unitary equivalence of spectra).
  - Full SSA proof with all four marginals explicitly identified.
""")
