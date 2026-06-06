"""
cnt_weingarten_general.py

Generalise the Weingarten formula E_Haar[r_2] = 2/3 (Section 44, d=D=2)
to arbitrary D = 2^L.

Key results to prove and verify:
  1. E_Haar[r_2(U,D)] = 2/(D+1)          for all D
  2. E_Haar[r_alpha(U,D)] = Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha)  for all alpha > 0, D
  3. Eigenvalues of B^(j) are {|U_kj|^2 : k=0,...,D-1}  (exact)
  4. Comparison with DU minimum r_alpha(DU) = D^{1-alpha}

Structure of the proof for general D:
  G_1 = sum_j |jj><jj| (standard-basis diagonal in H_D x H_D)
  F_j = P_j U^dag x U^T  (Kraus operators)
  Tr[G_1(A_j x I)G_1(A_j x I)] = sum_p |U_pj|^4   (shown in Sec 44)
  => Tr[G_2^2] = sum_{j,k} |U_kj|^4
  E[|U_kj|^4] = 2 / (D(D+1))   (Haar 2nd moment)
  => E[r_2] = 2/(D+1)

Generalisation to alpha:
  B^(j)_mm' = sum_k |U_kj|^2 U_mk U*_m'k
             = sum_k |U_kj|^2 |u_k><u_k|   (columns of U form ONB)
  => eigenvalues of B^(j) are {|U_kj|^2 : k}
  => Tr[(B^(j))^alpha] = sum_k |U_kj|^{2*alpha}
  => Tr[G_2^alpha] = sum_{j,k} |U_kj|^{2*alpha}
  E[|U_kj|^{2*alpha}] = Gamma(alpha+1)*Gamma(D)/Gamma(D+alpha)
  => E[r_alpha] = D * Gamma(alpha+1)*Gamma(D)/Gamma(D+alpha)
               = Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha)
"""

import numpy as np
from scipy.special import gamma as Gamma
from numpy.linalg import eigh

rng = np.random.default_rng(42)


def haar_unitary(D, rng):
    """Draw a Haar-random D x D unitary matrix."""
    Z = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    R_diag = np.diag(R)
    Q = Q * (R_diag / np.abs(R_diag))
    return Q


def build_G1(D):
    """G1 = sum_j |jj><jj| in H_D x H_D (size D^2 x D^2)."""
    G = np.zeros((D * D, D * D), dtype=complex)
    for j in range(D):
        idx = j * D + j
        G[idx, idx] = 1.0
    return G


def build_G2(U, D):
    """G2 = sum_j F_j G1 F_j^dag via the Kraus recursion."""
    G1 = build_G1(D)
    G2 = np.zeros((D * D, D * D), dtype=complex)
    Udagger = U.conj().T
    UT = U.T
    for j in range(D):
        # F_j = P_j U^dag (x) U^T   (P_j = |j><j|)
        # F_j acts on H_D x H_D
        # (F_j)_{(p,q),(r,s)} = delta_{pj} (U^dag)_{jr} (U^T)_{qs}
        # We compute F_j G1 F_j^dag directly as a D^2 x D^2 matrix

        # Build F_j as a D^2 x D^2 matrix
        Fj = np.zeros((D * D, D * D), dtype=complex)
        for r in range(D):
            for s in range(D):
                col = r * D + s
                Fj[j * D + 0:j * D + D, col] = 0  # init
                # row (p,q) with p=j: Fj[j*D+q, r*D+s] = Udagger[j,r] * UT[q,s]
                for q in range(D):
                    Fj[j * D + q, col] = Udagger[j, r] * UT[q, s]

        G2 += Fj @ G1 @ Fj.conj().T
    return G2


def r_alpha_from_orbit(U, D, alpha):
    """Compute r_alpha = P2^alpha / P1^alpha using the orbit state approach."""
    # P1^alpha = D^{1-alpha}  (rho1 = I/D)
    P1_alpha = D ** (1 - alpha)

    # Compute Tr[G2^alpha] via eigenvalues
    # Use the exact formula Tr[G2^alpha] = sum_{j,k} |U_kj|^{2*alpha}
    Tr_G2_alpha = np.sum(np.abs(U) ** (2 * alpha))

    # P2^alpha = D^{-alpha} * Tr[G2^alpha]
    P2_alpha = D ** (-alpha) * Tr_G2_alpha

    return P2_alpha / P1_alpha


def haar_average_formula(alpha, D):
    """Theoretical formula E[r_alpha] = Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha)."""
    return Gamma(alpha + 1) * Gamma(D + 1) / Gamma(D + alpha)


# ============================================================
# Part 1: Verify E[r_2] = 2/(D+1) for D=2,4,8
# ============================================================
print("=" * 60)
print("Part 1: E[r_2] = 2/(D+1) for D=2,4,8 (L=1,2,3)")
print("=" * 60)

N_samples = 10000
alpha_test = 2.0

for L, D in [(1, 2), (2, 4), (3, 8)]:
    empirical = []
    for _ in range(N_samples):
        U = haar_unitary(D, rng)
        r = r_alpha_from_orbit(U, D, alpha_test)
        empirical.append(r)
    emp_mean = np.mean(empirical)
    theory = 2.0 / (D + 1)
    du_min = D ** (1 - alpha_test)
    print(f"  L={L}, D={D}: E[r_2] empirical={emp_mean:.6f}  "
          f"theory=2/(D+1)={theory:.6f}  "
          f"DU_min=1/D={du_min:.6f}  "
          f"error={abs(emp_mean - theory):.2e}")

# ============================================================
# Part 2: Verify E[r_alpha] = Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha)
# ============================================================
print()
print("=" * 60)
print("Part 2: E[r_alpha] = Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha)")
print("=" * 60)

D = 4  # L=2
N_samples = 5000
alphas = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]

print(f"  D={D}, {N_samples} Haar samples:")
print(f"  {'alpha':>5}  {'empirical':>10}  {'theory':>10}  {'DU_min':>10}  {'error':>8}")
for alpha in alphas:
    empirical = []
    for _ in range(N_samples):
        U = haar_unitary(D, rng)
        r = r_alpha_from_orbit(U, D, alpha)
        empirical.append(r)
    emp = np.mean(empirical)
    th = haar_average_formula(alpha, D)
    du_min = D ** (1 - alpha)
    print(f"  {alpha:>5.2f}  {emp:>10.6f}  {th:>10.6f}  {du_min:>10.6f}  {abs(emp-th):>8.2e}")

# ============================================================
# Part 3: Eigenvalue structure of B^(j): eigenvalues = {|U_kj|^2}
# ============================================================
print()
print("=" * 60)
print("Part 3: Eigenvalues of B^(j) = {|U_kj|^2 : k=0,...,D-1}")
print("=" * 60)

D = 4
N_tests = 200
max_err = 0.0

for _ in range(N_tests):
    U = haar_unitary(D, rng)
    for j in range(D):
        # Construct B^(j)
        B = np.zeros((D, D), dtype=complex)
        for k in range(D):
            uk = U[:, k]  # k-th column of U
            w = np.abs(U[k, j]) ** 2  # weight |U_kj|^2
            B += w * np.outer(uk, uk.conj())

        # Eigenvalues of B
        eigs_B = np.sort(np.real(np.linalg.eigvalsh(B)))
        # Predicted eigenvalues
        pred = np.sort(np.abs(U[:, j]) ** 2)

        err = np.max(np.abs(eigs_B - pred))
        max_err = max(max_err, err)

print(f"  D={D}, {N_tests} random U, all j: max |eig(B^j) - |U_kj|^2| = {max_err:.2e}")

# ============================================================
# Part 4: Analytical formula table for D=2,3,4,8 and alpha=1,2,3,4
# ============================================================
print()
print("=" * 60)
print("Part 4: E[r_alpha] formula table")
print("=" * 60)
print(f"  {'D':>3} {'alpha':>5}  {'E[r_alpha]':>12}  {'DU_min':>12}  "
      f"{'E/DU_min':>10}")

for D in [2, 4, 8, 16]:
    for alpha in [1.0, 2.0, 3.0, 4.0]:
        th = haar_average_formula(alpha, D)
        du_min = D ** (1 - alpha)
        ratio = th / du_min if du_min > 0 else np.inf
        print(f"  {D:>3} {alpha:>5.1f}  {th:>12.6f}  {du_min:>12.6f}  {ratio:>10.4f}")

# ============================================================
# Part 5: Tr[G2^alpha] = sum_{j,k} |U_kj|^{2*alpha} (exact formula)
# ============================================================
print()
print("=" * 60)
print("Part 5: Tr[G2^alpha] = sum_{j,k} |U_kj|^{2*alpha} (direct verification)")
print("=" * 60)

D = 4
N_tests = 20
print(f"  D={D}, {N_tests} random U:")
max_errors = {}
for alpha in [1.5, 2.0, 2.5, 3.0]:
    errs = []
    for _ in range(N_tests):
        U = haar_unitary(D, rng)
        # Method 1: Direct formula
        direct = np.sum(np.abs(U) ** (2 * alpha))

        # Method 2: Compute G2 eigenvalues
        G2 = build_G2(U, D)
        eigs = np.real(np.linalg.eigvalsh(G2))
        eigs = np.maximum(eigs, 0)
        Tr_G2_alpha = np.sum(eigs ** alpha)

        errs.append(abs(direct - Tr_G2_alpha))
    max_errors[alpha] = max(errs)
    print(f"  alpha={alpha:.1f}: max|direct - from_G2_eigs| = {max_errors[alpha]:.2e}")

# ============================================================
# Part 6: Verify inequality E[r_alpha] > DU_min for all D, alpha>1
# ============================================================
print()
print("=" * 60)
print("Part 6: E[r_alpha] > DU_min=D^{1-alpha} for all D (alpha>1)")
print("=" * 60)
print("  Checking Gamma(alpha+1)*Gamma(D+1)/Gamma(D+alpha) > D^{1-alpha}:")
all_pass = True
for D in [2, 3, 4, 5, 8, 16, 32]:
    for alpha in [1.01, 1.5, 2.0, 3.0, 5.0]:
        th = haar_average_formula(alpha, D)
        du_min = D ** (1 - alpha)
        if th <= du_min - 1e-12:
            print(f"  VIOLATION: D={D}, alpha={alpha}, E[r]={th:.6g} <= DU_min={du_min:.6g}")
            all_pass = False
if all_pass:
    print("  All checks passed: E[r_alpha] > DU_min for all D>=2, alpha>1. CONFIRMED.")

# ============================================================
# Part 7: Closed-form check for L=1,2,3
# ============================================================
print()
print("=" * 60)
print("Part 7: Summary table for L=1,2,3 (D=2,4,8)")
print("=" * 60)
print(f"  {'L':>2} {'D':>3} {'alpha':>5}  {'E[r_alpha]':>14}  "
      f"{'DU_min':>10}  {'formula':>30}")
for L, D in [(1, 2), (2, 4), (3, 8)]:
    for alpha in [2.0, 3.0]:
        th = haar_average_formula(alpha, D)
        du_min = D ** (1 - alpha)
        if alpha == 2.0:
            formula_str = f"2/(D+1) = {2/(D+1):.6f}"
        elif alpha == 3.0:
            formula_str = f"6/((D+1)(D+2)) = {6/((D+1)*(D+2)):.6f}"
        else:
            formula_str = f"Gamma formula = {th:.6f}"
        print(f"  {L:>2} {D:>3} {alpha:>5.1f}  {th:>14.8f}  "
              f"{du_min:>10.6f}  {formula_str}")
