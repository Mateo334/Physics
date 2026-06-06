"""
cnt_var_alpha.py  — Section 48

Exact closed-form Var[r_alpha] for all alpha > 0 and D >= 2.
Uses the Beta-Dirichlet framework of Section 47.

KEY RESULTS:

  F_alpha(D) = (D-1) * Gamma(alpha+1)^2 * Gamma(D-1)^2
               / [Gamma(D-1+alpha)^2 * (D-1+2*alpha)]

  Var[r_alpha] = [Gamma(2*alpha+1) + 2(D-1)*Gamma(alpha+1)^2] * Gamma(D) / Gamma(D+2*alpha)
               + (D-1)^3 * Gamma(alpha+1)^2 * Gamma(D-1)^2
                 / [Gamma(D-1+alpha)^2 * (D-1+2*alpha)]
               - Gamma(alpha+1)^2 * Gamma(D+1)^2 / Gamma(D+alpha)^2

Special cases:
  alpha=1: Var = 0 (r_1 = 1 always)
  alpha=2: Var = 4(D-1)/[D^2(D+1)^2(D+3)]  (Section 47)
  alpha=3: Var = 48(D-1)^2(D^2+1) / [D^2(D+1)^2(D+2)^2(D+3)(D+5)]  (new)
"""

import numpy as np
from scipy.special import gamma as G
from fractions import Fraction

rng = np.random.default_rng(42)


def haar(D):
    Z = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def A_alpha(alpha, D):
    """E[|U_11|^{4*alpha}] = Gamma(2a+1)*Gamma(D)/Gamma(D+2a)"""
    return G(2 * alpha + 1) * G(D) / G(D + 2 * alpha)


def B_alpha(alpha, D):
    """E[|U_11|^{2a}|U_21|^{2a}] same col = Gamma(a+1)^2*Gamma(D)/Gamma(D+2a)"""
    return G(alpha + 1)**2 * G(D) / G(D + 2 * alpha)


def F_alpha(alpha, D):
    """E[|U_11|^{2a}|U_22|^{2a}] diff row+col.
    = (D-1)*Gamma(a+1)^2*Gamma(D-1)^2 / [Gamma(D-1+a)^2*(D-1+2a)]
    Proof: F = E[d^a]^2 * E[(1-x2)^{2a}]
      E[d^a] = Gamma(a+1)*Gamma(D-1)/Gamma(D-1+a)  (d ~ Beta(1,D-2))
      E[(1-x2)^{2a}] = (D-1)/(D-1+2a)  (1-x2 ~ Beta(D-1,1))
    """
    Ed_a = G(alpha + 1) * G(D - 1) / G(D - 1 + alpha)
    E_w = (D - 1) / (D - 1 + 2 * alpha)
    return Ed_a**2 * E_w


def Er_alpha(alpha, D):
    """E[r_alpha] = Gamma(a+1)*Gamma(D+1)/Gamma(D+a)"""
    return G(alpha + 1) * G(D + 1) / G(D + alpha)


def Var_alpha(alpha, D):
    """Var[r_alpha] closed form."""
    A = A_alpha(alpha, D)
    B = B_alpha(alpha, D)
    F = F_alpha(alpha, D)
    Er = Er_alpha(alpha, D)
    return A + 2 * (D - 1) * B + (D - 1)**2 * F - Er**2


# ============================================================
# Part 1: Verify building blocks A, B, F against MC
# ============================================================
print("=" * 64)
print("Part 1: Verify A_alpha, B_alpha, F_alpha vs MC")
print("=" * 64)
N = 300_000
print(f"  {'alpha':>6}  {'D':>2}  "
      f"{'A_err':>10}  {'B_err':>10}  {'F_err':>10}")
for alpha in [0.5, 1.0, 1.5, 2.0, 3.0]:
    for D in [2, 3, 4]:
        A = A_alpha(alpha, D); B = B_alpha(alpha, D); F = F_alpha(alpha, D)
        As, Bs, Fs = [], [], []
        for _ in range(N):
            U = haar(D)
            As.append(abs(U[0, 0])**(4 * alpha))
            Bs.append(abs(U[0, 0])**(2 * alpha) * abs(U[1, 0])**(2 * alpha))
            Fs.append(abs(U[0, 0])**(2 * alpha) * abs(U[1, 1])**(2 * alpha))
        print(f"  {alpha:>6.1f}  {D}  "
              f"{abs(A-np.mean(As)):>10.2e}  "
              f"{abs(B-np.mean(Bs)):>10.2e}  "
              f"{abs(F-np.mean(Fs)):>10.2e}")
    print()

# ============================================================
# Part 2: Var[r_alpha] formula vs MC
# ============================================================
print("=" * 64)
print("Part 2: Var[r_alpha] formula vs MC")
print("=" * 64)
N2 = 200_000
print(f"  {'alpha':>6}  D  {'Var_form':>14}  {'Var_MC':>12}  {'err':>10}")
for alpha in [0.5, 1.0, 1.5, 2.0, 3.0]:
    for D in [2, 3, 4, 5]:
        V_f = Var_alpha(alpha, D)
        r_vals = [np.sum(abs(haar(D))**(2 * alpha)) / D for _ in range(N2)]
        V_mc = np.var(r_vals)
        print(f"  {alpha:>6.1f}  {D}  {V_f:>14.8f}  {V_mc:>12.8f}  "
              f"{abs(V_f - V_mc):>10.2e}")
    print()

# ============================================================
# Part 3: Exact rational values for alpha=2,3 at small D
# ============================================================
print("=" * 64)
print("Part 3: Exact rational Var[r_alpha] at alpha=2,3")
print("=" * 64)
print()
print("alpha=2 (cross-check with Section 47):")
for D in range(2, 7):
    V = Var_alpha(2.0, D)
    V_sec47 = 4 * (D - 1) / (D**2 * (D + 1)**2 * (D + 3))
    print(f"  D={D}: formula={V:.10f}, 4(D-1)/[D^2(D+1)^2(D+3)]={V_sec47:.10f}, "
          f"match={abs(V-V_sec47)<1e-10}")

print()
print("alpha=3:")
for D in range(2, 7):
    V = Var_alpha(3.0, D)
    # The exact rational for alpha=3: derive by substitution
    # A_3 = Gamma(7)*Gamma(D)/Gamma(D+6) = 720/poch(D,6)
    # B_3 = Gamma(4)^2*Gamma(D)/Gamma(D+6) = 36/poch(D,6)
    # F_3 = Gamma(4)^2*Gamma(D-1)^2*(D-1) / [Gamma(D+2)^2*(D+5)]
    #      = 36*(D-2)!^2*(D-1) / [(D+1)!^2*(D+5)]
    #      = 36*(D-1)/[(D-1)^2*D^2*(D+1)^2*(D+5)]... wait let me use Fraction
    A3 = Fraction(720, 1)
    for i in range(6): A3 /= (D + i)
    B3 = Fraction(36, 1)
    for i in range(6): B3 /= (D + i)
    # F3: Gamma(4)=6; E[d^3] for d~Beta(1,D-2) = 6*(D-2)!/D! = 6/(D*(D-1)*(D-2))...
    # Actually Gamma(4)*Gamma(D-1)/Gamma(D+2) = 6*(D-2)!/(D+1)! = 6/((D-1)*D*(D+1))
    Ed3 = Fraction(6, (D - 1) * D * (D + 1))
    Ew3 = Fraction(D - 1, D + 5)
    F3 = Ed3**2 * Ew3
    Er3 = Fraction(6, 1)
    for i in range(1, 4): Er3 *= i  # 6 = 3!
    Er3_numer = Er3 * Fraction(D * 1, 1)
    # E[r_3] = Gamma(4)*Gamma(D+1)/Gamma(D+3) = 6*D!/(D+2)! = 6/(D+1)(D+2)
    Er3_val = Fraction(6, (D + 1) * (D + 2))
    Var3 = A3 + 2 * (D - 1) * B3 + (D - 1)**2 * F3 - Er3_val**2
    print(f"  D={D}: Var[r_3]={Var3} = {float(Var3):.8f}  (formula: {V:.8f})")

# ============================================================
# Part 4: sigma/mu vs D for multiple alpha
# ============================================================
print()
print("=" * 64)
print("Part 4: sigma[r_alpha]/E[r_alpha] = coefficient of variation")
print("=" * 64)
print(f"  {'D':>4}  {'alpha=0.5':>12}  {'alpha=1.5':>12}  "
      f"{'alpha=2':>12}  {'alpha=3':>12}")
for D in [2, 4, 8, 16, 32, 64]:
    row = [D]
    for alpha in [0.5, 1.5, 2.0, 3.0]:
        V = Var_alpha(alpha, D)
        Er = Er_alpha(alpha, D)
        row.append(np.sqrt(max(V, 0)) / Er)
    print(f"  {D:>4}  {row[1]:>12.6f}  {row[2]:>12.6f}  "
          f"{row[3]:>12.6f}  {row[4]:>12.6f}")

# ============================================================
# Part 5: Large-D asymptotics of sigma/mu
# ============================================================
print()
print("=" * 64)
print("Part 5: Large-D asymptotics of sigma/mu")
print("=" * 64)
print()
print("  For large D: sigma[r_alpha]/E[r_alpha] ~ C_alpha / D^{alpha}")
print()
print(f"  {'alpha':>6}  {'D=16':>10}  {'D=32':>10}  {'D=64':>10}  "
      f"{'sigma/mu * D^alpha (D=64)':>28}")
for alpha in [0.5, 1.0, 1.5, 2.0, 3.0]:
    vals = []
    for D in [16, 32, 64]:
        V = Var_alpha(alpha, D)
        Er = Er_alpha(alpha, D)
        vals.append(np.sqrt(max(V, 0)) / Er)
    # sigma/mu * D^alpha at D=64
    C = vals[-1] * 64**alpha
    print(f"  {alpha:>6.1f}  {vals[0]:>10.6f}  {vals[1]:>10.6f}  {vals[2]:>10.6f}  "
          f"{C:>28.4f}")

print()
print("  Note: sigma[r_2]/mu = sqrt(D-1)/(D*sqrt(D+3)) ~ 1/D (Section 47)")
print("  For general alpha > 1: sigma/mu ~ C_alpha / D^alpha  (D -> inf)")
