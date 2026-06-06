"""
cnt_exact_variance.py  — Section 47

Exact closed-form formulas for F = E[|U_11|^4 |U_22|^4]
and Var[r_2] for the Haar measure on U(D).

KEY RESULTS (all proved analytically and verified numerically):

  F(D)      = 4 / [(D-1) D^2 (D+3)]
  Var[r_2]  = 4(D-1) / [D^2 (D+1)^2 (D+3)]
  sigma/mu  = sqrt(D-1) / [D sqrt(D+3)]  ~  1/D  as D -> inf

Derivation strategy (Beta-Dirichlet conditional approach):
  E[|U_11|^4 |U_22|^4]
  = E[d^2] * E[x_1^2 (1-x_2)^2]
  where:
    d ~ Beta(1, D-2)  =>  E[d^2] = 2/[(D-1)D]
    (x_1,...,x_D) ~ Dir(1,...,1)  =>  E[x_1^2(1-x_2)^2] = 2/[D(D+3)]
  Product = 4/[(D-1)D^2(D+3)].
"""

import numpy as np
from fractions import Fraction

rng = np.random.default_rng(42)


def haar(D):
    Z = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def poch(D, k):
    return int(np.prod([D + i for i in range(k)]))


# ============================================================
# Part 1: Closed-form F(D) via Beta-Dirichlet approach
# ============================================================
print("=" * 62)
print("Part 1: Closed-form F(D) = 4/[(D-1) D^2 (D+3)]")
print("=" * 62)
print()
print("Proof sketch:")
print("  Given 1st column u_1: |u_2[2]|^2 = (1-|u_1[2]|^2)*d, d ~ Beta(1,D-2)")
print("  => E[d^2] = 2/[(D-1)D]  (2nd moment of Beta(1,D-2))")
print("  (x_1,...,x_D) ~ Dir(1,...,1):")
print("  E[x_1^2(1-x_2)^2] = E[x_1^2] - 2E[x_1^2 x_2] + E[x_1^2 x_2^2]")
print("  = 2/D(D+1) - 4/D(D+1)(D+2) + 4/D(D+1)(D+2)(D+3)")
print("  = 2(D+1-2+2/(D+3)) / [D(D+1)(D+2)(D+3)]  <- telescopes")
print("  = 2/[D(D+3)]  (after simplification)")
print("  F = E[d^2]*E[x_1^2(1-x_2)^2] = 4/[(D-1)D^2(D+3)]  QED")
print()

N_mc = 500_000
print(f"  {'D':>3}  {'F_exact':>14}  {'F_MC':>12}  {'error':>10}")
for D in [2, 3, 4, 5, 6, 8]:
    F_ex = Fraction(4, (D - 1) * D**2 * (D + 3))
    F_mc = np.mean([abs(haar(D)[0, 0])**4 * abs(haar(D)[1, 1])**4
                    for _ in range(N_mc)])  # bug: two matrices
    # correct:
    vals = []
    for _ in range(N_mc):
        U = haar(D)
        vals.append(abs(U[0, 0])**4 * abs(U[1, 1])**4)
    F_mc = np.mean(vals)
    print(f"  {D:>3}  {float(F_ex):>14.8f}  {F_mc:>12.8f}  {abs(float(F_ex)-F_mc):>10.2e}")


# ============================================================
# Part 2: Exact Var[r_2] = 4(D-1)/[D^2(D+1)^2(D+3)]
# ============================================================
print()
print("=" * 62)
print("Part 2: Exact Var[r_2] = 4(D-1)/[D^2*(D+1)^2*(D+3)]")
print("=" * 62)
print()
print("Derivation:")
print("  r_2 = (1/D) sum_{j,k} |U_kj|^4")
print("  E[r_2^2] = A + 2(D-1)B + (D-1)^2 F")
print("  where A=24/poch(D,4), B=4/poch(D,4), F=4/[(D-1)D^2(D+3)]")
print("  Var = E[r_2^2] - (E[r_2])^2 = A+2(D-1)B+(D-1)^2*F - 4/(D+1)^2")
print("  Numerator: 8(D+2)/[D(D+1)(D+2)(D+3)] + 4(D-1)/[D^2(D+3)]")
print("           - 4/(D+1)^2")
print("  = 8/[D(D+1)(D+3)] + 4(D-1)/[D^2(D+3)] - 4/(D+1)^2")
print("  = [8D + 4(D^2-1)] / [D^2(D+1)(D+3)] - 4/(D+1)^2")
print("  = 4(D^2+2D-1)/[D^2(D+1)(D+3)] - 4/(D+1)^2")
print("  = 4/(D+1) * {(D^2+2D-1)(D+1) - D^2(D+3)} / [D^2(D+1)^2(D+3)]")
print("  Bracket: D^3+3D^2+D-1 - D^3-3D^2 = D-1")
print("  => Var[r_2] = 4(D-1)/[D^2(D+1)^2(D+3)]  QED")
print()

N_mc2 = 200_000
print(f"  {'D':>3}  {'Var_exact':>14}  {'Var_MC':>12}  {'error':>10}")
for D in [2, 3, 4, 5, 6, 8]:
    V_ex = Fraction(4 * (D - 1), D**2 * (D + 1)**2 * (D + 3))
    r2_vals = [np.sum(abs(haar(D))**4) / D for _ in range(N_mc2)]
    V_mc = np.var(r2_vals)
    print(f"  {D:>3}  {float(V_ex):>14.8f}  {V_mc:>12.8f}  {abs(float(V_ex)-V_mc):>10.2e}")


# ============================================================
# Part 3: sigma/mean exact formula; large-D asymptotics
# ============================================================
print()
print("=" * 62)
print("Part 3: sigma/mean = sqrt(D-1)/[D*sqrt(D+3)]  ~  1/D")
print("=" * 62)
print()
print("  sigma/mean = sqrt(Var)/E[r_2]")
print("  = sqrt(4(D-1)/[D^2(D+1)^2(D+3)]) / (2/(D+1))")
print("  = sqrt((D-1)/[D^2(D+3)])")
print("  = sqrt(D-1) / [D*sqrt(D+3)]")
print("  Large D: ~ sqrt(D)/(D*sqrt(D)) = 1/D  (NOT 1/sqrt(D))  QED")
print()
print(f"  {'D':>4}  {'sigma/mean exact':>18}  {'1/D':>10}  {'1/sqrt(D)':>12}")
for D in [2, 4, 8, 16, 32, 64]:
    sm = np.sqrt((D - 1) / (D**2 * (D + 3)))
    print(f"  {D:>4}  {sm:>18.8f}  {1/D:>10.6f}  {1/np.sqrt(D):>12.6f}")


# ============================================================
# Part 4: Chebyshev concentration bound
# ============================================================
print()
print("=" * 62)
print("Part 4: Chebyshev bound  Pr[|r_2 - 2/(D+1)| > eps] <= Var/eps^2")
print("=" * 62)
print()
print("  For eps = 0.03 * E[r_2] = 0.03 * 2/(D+1):")
eps_frac = 0.03
print(f"  {'D':>4}  {'E[r_2]':>10}  {'Var_exact':>14}  {'Chebyshev_bound':>18}")
for D in [4, 8, 16, 32, 64]:
    Er2 = 2.0 / (D + 1)
    Var = 4.0 * (D - 1) / (D**2 * (D + 1)**2 * (D + 3))
    eps = eps_frac * Er2
    bound = Var / eps**2
    print(f"  {D:>4}  {Er2:>10.6f}  {Var:>14.10f}  {bound:>18.6f}")


# ============================================================
# Part 5: Section 46 correction — F_code was wrong at D=3
# ============================================================
print()
print("=" * 62)
print("Part 5: Correction — previous code gave wrong F(D=3)")
print("=" * 62)
print()
print("  Old code (Weingarten linear system) gave F(D=3) = 3/32 = 0.09375")
print("  This was wrong: the 5x5 Gram matrix is SINGULAR for D=3,m=4")
print("  (representation (1^4) with ell=4>D=3 causes degeneracy)")
print("  Correct: F(D=3) = 1/27 = 0.037037... (proved via Beta-Dir approach)")
print()
print("  Impact on Section 46 Var[r_2] formula:")
print("  Old (D=3): formula used F=3/32 -> Var=0.236 (WRONG, factor ~25 off)")
print("  New (D=3): F=1/27 -> Var=1/108=0.00926 (correct, matches MC)")
print()

# Final summary table
print("=" * 62)
print("Final summary: F(D) and Var[r_2] exact formulas")
print("=" * 62)
print(f"  {'D':>3}  {'F = 4/(...)':>20}  {'Var = 4(D-1)/...':>20}  {'sigma/mean':>12}")
for D in [2, 3, 4, 5, 6, 8, 16]:
    F_ex = 4.0 / ((D - 1) * D**2 * (D + 3))
    V_ex = 4.0 * (D - 1) / (D**2 * (D + 1)**2 * (D + 3))
    sm = np.sqrt((D - 1) / (D**2 * (D + 3)))
    print(f"  {D:>3}  {F_ex:>20.10f}  {V_ex:>20.10f}  {sm:>12.8f}")
