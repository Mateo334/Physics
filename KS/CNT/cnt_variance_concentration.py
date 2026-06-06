"""
cnt_variance_concentration.py

Variance and concentration of r_alpha under Haar measure.

Key results:
  E[r_2] = 2/(D+1)   (Section 45)
  Var[r_2] = E[r_2^2] - (2/(D+1))^2

  For D=2: Var[r_2] = 1/45  (exact, using F=A for 2x2 unitaries)
  Concentration: std[r_2]/E[r_2] ~ C/D -> 0 as D->inf

Variance decomposition:
  r_2 = (1/D) * sum_{j,k} |U_kj|^4

  E[r_2^2] = A + 2*(D-1)*B + (D-1)^2*F   (using B=C by symmetry)

  where:
    A = E[|U_11|^8]               = Gamma(5)*Gamma(D)/Gamma(D+4)  = 24/poch(D,4)
    B = E[|U_11|^4 |U_21|^4]      = Gamma(3)^2*Gamma(D)/Gamma(D+4) = 4/poch(D,4)
      (same column j, different row k)
    C = E[|U_11|^4 |U_12|^4] = B  (same row k, different col j, by symmetry)
    F = E[|U_11|^4 |U_22|^4]      (different row AND different col)

  For D=2: |U_11|^2 = |U_22|^2 for all 2x2 unitaries => F = A = 1/5.
    => E[r_2^2] = A + 2*B + B... wait:
    For D=2: E[r_2^2] = A + 2*(1)*B + (1)^2*F = A + 2B + F = 2A+2B (since F=A)
    = 2*(24 + 4)/120 = 56/120 = 7/15
    Var[r_2] = 7/15 - (2/3)^2 = 7/15 - 4/9 = 21/45 - 20/45 = 1/45 (EXACT).

General formula for Var[r_2]:
  Use the F = A - correction identity to be derived.
"""

import numpy as np
from scipy.special import gamma as Gamma

rng = np.random.default_rng(42)


def haar_unitary(D, rng):
    Z = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    R_diag = np.diag(R)
    return Q * (R_diag / np.abs(R_diag))


def r2(U, D):
    return np.sum(np.abs(U) ** 4) / D


def poch(a, k):
    return float(np.prod([a + i for i in range(int(k))]))


def haar_mean_formula(alpha, D):
    return Gamma(alpha + 1) * Gamma(D + 1) / Gamma(D + alpha)


# ============================================================
# Part 1: Joint moments A, B, C, F for D=4 (using SINGLE matrix)
# ============================================================
print("=" * 60)
print("Part 1: Joint moments A,B,C,F — all from single U matrix")
print("=" * 60)

D = 4
N = 500000
A_vals, B_vals, C_vals, F_vals = [], [], [], []
for _ in range(N):
    U = haar_unitary(D, rng)
    A_vals.append(np.abs(U[0, 0]) ** 8)
    B_vals.append(np.abs(U[0, 0]) ** 4 * np.abs(U[1, 0]) ** 4)  # same col
    C_vals.append(np.abs(U[0, 0]) ** 4 * np.abs(U[0, 1]) ** 4)  # same row
    F_vals.append(np.abs(U[0, 0]) ** 4 * np.abs(U[1, 1]) ** 4)  # diff row, diff col

A_emp = np.mean(A_vals)
B_emp = np.mean(B_vals)
C_emp = np.mean(C_vals)
F_emp = np.mean(F_vals)
A_th = 24 / poch(D, 4)
B_th = 4 / poch(D, 4)

print(f"  D={D}, N={N:,}")
print(f"  A = E[|U_11|^8]              = {A_emp:.8f}  theory={A_th:.8f}  err={abs(A_emp-A_th):.2e}")
print(f"  B = E[|U_11|^4|U_21|^4]      = {B_emp:.8f}  theory={B_th:.8f}  err={abs(B_emp-B_th):.2e}")
print(f"  C = E[|U_11|^4|U_12|^4]      = {C_emp:.8f}  (should = B)")
print(f"  F = E[|U_11|^4|U_22|^4]      = {F_emp:.8f}  (diff row+col)")
print(f"  B=C? |B-C|={abs(B_emp-C_emp):.2e}")

# ============================================================
# Part 2: D=2 exact variance
# ============================================================
print()
print("=" * 60)
print("Part 2: D=2 exact variance — Var[r_2] = 1/45")
print("=" * 60)

# For D=2: |U_11|^2 = |U_22|^2 always, so F=A=1/5
A2 = 24 / poch(2, 4)   # = 24/120 = 1/5
B2 = 4 / poch(2, 4)    # = 4/120 = 1/30
F2 = A2                  # exact for D=2

Er2sq_exact = A2 + 2 * (2 - 1) * B2 + (2 - 1)**2 * F2  # = 2A+2B = 56/120 = 7/15
Er2_exact = 2 / (2 + 1)  # = 2/3
Var_exact = Er2sq_exact - Er2_exact**2  # = 7/15 - 4/9 = 1/45

print(f"  D=2:  A=1/5={A2:.6f}, B=1/30={B2:.6f}, F=A=1/5={F2:.6f}")
print(f"  E[r_2^2] = 2A + 2B = 7/15 = {Er2sq_exact:.8f}")
print(f"  (E[r_2])^2 = (2/3)^2 = 4/9 = {Er2_exact**2:.8f}")
print(f"  Var[r_2] = 7/15 - 4/9 = 1/45 = {Var_exact:.8f}")

# Verify numerically
D = 2
N = 500000
r2_D2 = [r2(haar_unitary(D, rng), D) for _ in range(N)]
print(f"  Numerical: Var[r_2] = {np.var(r2_D2):.8f}  (theory 1/45={1/45:.8f})")

# ============================================================
# Part 3: Case F for D=2,3,4,8 numerically (single matrix)
# ============================================================
print()
print("=" * 60)
print("Part 3: Case F = E[|U_11|^4|U_22|^4] for various D")
print("=" * 60)

N = 200000
print(f"  {'D':>3}  {'F_emp':>12}  {'A_th':>12}  {'B_th':>12}  {'F/A':>8}  {'F/B':>8}")
F_by_D = {}
for D in [2, 3, 4, 5, 6, 8, 10]:
    F_samp = []
    for _ in range(N):
        U = haar_unitary(D, rng)
        F_samp.append(np.abs(U[0, 0]) ** 4 * np.abs(U[1, 1]) ** 4)
    F_val = np.mean(F_samp)
    F_by_D[D] = F_val
    A_th = 24 / poch(D, 4)
    B_th = 4 / poch(D, 4)
    print(f"  {D:>3}  {F_val:>12.8f}  {A_th:>12.8f}  {B_th:>12.8f}  {F_val/A_th:>8.4f}  {F_val/B_th:>8.4f}")

# ============================================================
# Part 4: Exact Var[r_2] using A, B, F numerically
# ============================================================
print()
print("=" * 60)
print("Part 4: Var[r_2] from decomposition vs direct Monte Carlo")
print("=" * 60)

N_mc = 300000
print(f"  {'D':>3}  {'Var_direct':>12}  {'Var_fromABF':>12}  {'E[r2]':>8}  {'std/mean':>10}  {'std*D':>10}")
for D in [2, 3, 4, 6, 8, 12, 16]:
    # Direct Monte Carlo
    r2_vals = [r2(haar_unitary(D, rng), D) for _ in range(N_mc)]
    var_direct = np.var(r2_vals)
    mean_direct = np.mean(r2_vals)

    # From decomposition using numerical F
    if D in F_by_D:
        F_val = F_by_D[D]
    else:
        F_samp = []
        for _ in range(100000):
            U = haar_unitary(D, rng)
            F_samp.append(np.abs(U[0, 0]) ** 4 * np.abs(U[1, 1]) ** 4)
        F_val = np.mean(F_samp)

    A_th = 24 / poch(D, 4)
    B_th = 4 / poch(D, 4)
    Er2sq_form = A_th + 2 * (D - 1) * B_th + (D - 1)**2 * F_val
    Er2_th = 2.0 / (D + 1)
    var_form = Er2sq_form - Er2_th**2

    print(f"  {D:>3}  {var_direct:>12.8f}  {var_form:>12.8f}  "
          f"{mean_direct:>8.6f}  {np.std(r2_vals)/mean_direct:>10.6f}  "
          f"{np.std(r2_vals)*D:>10.6f}")

# ============================================================
# Part 5: Concentration — std/E[r_2] vs D
# ============================================================
print()
print("=" * 60)
print("Part 5: Concentration — std/E[r_2] vs D")
print("=" * 60)

N = 200000
print(f"  {'D':>4}  {'E[r_2]':>10}  {'std':>10}  {'std/mean':>10}  "
      f"{'std*D^{3/2}':>14}  {'Var*D^3':>12}")
for D in [2, 4, 8, 16, 32, 64]:
    r2_vals = [r2(haar_unitary(D, rng), D) for _ in range(N)]
    m = np.mean(r2_vals)
    s = np.std(r2_vals)
    v = np.var(r2_vals)
    print(f"  {D:>4}  {m:>10.6f}  {s:>10.6f}  {s/m:>10.6f}  "
          f"{s * D**1.5:>14.6f}  {v * D**3:>12.6f}")

# ============================================================
# Part 6: Exact Var[r_2] for D=2 and asymptotic scaling
# ============================================================
print()
print("=" * 60)
print("Part 6: Asymptotic analysis of Var[r_2]")
print("=" * 60)

print("  For large D, using F ~ (E[|U_11|^4])^2 = 4/(D(D+1))^2 ~ 4/D^4:")
print("  E[r_2^2] = A + 2(D-1)*B + (D-1)^2*F")
print("  Leading: A~24/D^4, B~4/D^4, F~4/D^4")
print("  E[r_2^2] ~ (24 + 8D + 4D^2)/D^4 = 4/D^2 + 8/D^3 + 24/D^4")
print("  (E[r_2])^2 = 4/(D+1)^2 = 4/D^2 - 8/D^3 + 12/D^4 - ...")
print("  Var[r_2] ~ 16/D^3 + 12/D^4 + ... => std ~ 4/D^{3/2}")
print("  std/mean ~ (4/D^{3/2}) / (2/D) = 2/D^{1/2} = 2/sqrt(D)")
print()
print("  Checking std/mean * sqrt(D) ~ const:")
N = 200000
print(f"  {'D':>4}  {'std/mean*sqrt(D)':>18}")
for D in [4, 8, 16, 32, 64]:
    r2_vals = [r2(haar_unitary(D, rng), D) for _ in range(N)]
    m = np.mean(r2_vals)
    s = np.std(r2_vals)
    print(f"  {D:>4}  {s/m * D**0.5:>18.6f}")

print()
print("  Exact Var[r_2] for D=2:")
print(f"    Var[r_2] = 1/45 = {1/45:.10f}")
print()
print("  Summary: concentration theorem")
print("    std[r_2] / E[r_2] = O(1/sqrt(D)) -> 0 as D -> infinity")
print("    By Chebyshev: Pr[|r_2 - 2/(D+1)| > eps] <= Var[r_2]/eps^2 ~ C/(eps^2 * D^3)")
