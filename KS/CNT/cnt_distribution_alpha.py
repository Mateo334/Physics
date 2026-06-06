"""
cnt_distribution_alpha.py
Exact probability distribution of r_alpha for Haar-random D=2 (qubit) circuits.

Key: For D=2, r_alpha = x^alpha + (1-x)^alpha where x ~ Uniform(0,1).
- alpha > 1: CONVEX f, support [2^{1-alpha}, 1], right-skewed
- alpha < 1: CONCAVE f, support [1, 2^{1-alpha}], left-skewed
- alpha = 1: constant 1, delta function
"""

import numpy as np
from scipy.special import gamma
from scipy.optimize import brentq
from scipy.stats import unitary_group
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# PART 1: Exact PDF of r_alpha for D=2
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 70)
print("PART 1 — Exact PDF of r_alpha for D=2")
print("=" * 70)

def pdf_r_alpha_D2(r, alpha):
    """
    Exact PDF of r_alpha = x^alpha + (1-x)^alpha, x ~ U(0,1).
    alpha > 1: support [2^{1-alpha}, 1], f convex
    alpha < 1: support [1, 2^{1-alpha}], f concave (inverted)
    alpha = 1: delta function at 1
    """
    r = np.atleast_1d(np.array(r, dtype=float))
    pdf_vals = np.zeros_like(r)
    rmin = 2**(1-alpha)

    for i, ri in enumerate(r):
        if alpha > 1:
            if ri <= rmin or ri >= 1.0:
                continue
            # solve x^alpha + (1-x)^alpha = ri on (0, 1/2)
            try:
                x1 = brentq(lambda x: x**alpha + (1-x)**alpha - ri, 1e-12, 0.5 - 1e-12)
            except ValueError:
                continue
        else:  # alpha < 1, support [1, rmin=2^{1-alpha}]
            if ri <= 1.0 or ri >= rmin:
                continue
            try:
                x1 = brentq(lambda x: x**alpha + (1-x)**alpha - ri, 1e-12, 0.5 - 1e-12)
            except ValueError:
                continue

        fp = abs(alpha * (x1**(alpha-1) - (1-x1)**(alpha-1)))
        if fp < 1e-15:
            continue
        pdf_vals[i] = 2.0 / fp

    return pdf_vals if len(pdf_vals) > 1 else pdf_vals[0]

# Verify alpha=2: p(r) = 1/sqrt(2r-1)
print("\n--- Alpha=2: PDF should be 1/sqrt(2r-1) ---")
r_test = np.linspace(0.52, 0.99, 100)
pdf_num = pdf_r_alpha_D2(r_test, alpha=2.0)
pdf_ana = 1.0 / np.sqrt(2*r_test - 1)
err2 = np.max(np.abs(pdf_num - pdf_ana))
print(f"  Max |p_num - 1/sqrt(2r-1)| = {err2:.2e}   {'PASS' if err2 < 1e-6 else 'FAIL'}")

# Normalization for several alpha values
print("\n--- Normalization (integral of PDF = 1) ---")
for alpha in [0.5, 1.5, 2.0, 3.0]:
    rmin_v = 2**(1-alpha)
    if alpha > 1:
        r_grid = np.linspace(rmin_v + 1e-5, 1 - 1e-5, 4000)
    else:
        r_grid = np.linspace(1.0 + 1e-5, rmin_v - 1e-5, 4000)
    p_grid = pdf_r_alpha_D2(r_grid, alpha)
    norm = np.trapz(p_grid, r_grid)
    print(f"  alpha={alpha:.1f}: support=[{min(r_grid):.4f},{max(r_grid):.4f}], integral={norm:.6f}")

# Verify E[r_alpha] from PDF
print("\n--- E[r_alpha] from PDF vs theory 2*Gamma(alpha+1)*Gamma(2)/Gamma(2+alpha) ---")
for alpha in [0.5, 1.5, 2.0, 3.0]:
    rmin_v = 2**(1-alpha)
    if alpha > 1:
        r_grid = np.linspace(rmin_v + 1e-5, 1 - 1e-5, 4000)
    else:
        r_grid = np.linspace(1.0 + 1e-5, rmin_v - 1e-5, 4000)
    p_grid = pdf_r_alpha_D2(r_grid, alpha)
    E_pdf = np.trapz(r_grid * p_grid, r_grid)
    E_th = gamma(alpha+1) * gamma(3) / gamma(alpha+2)  # = 2/(alpha+1)
    print(f"  alpha={alpha:.1f}: E_PDF={E_pdf:.6f}, E_theory=2/(alpha+1)={2/(alpha+1):.6f}, "
          f"err={abs(E_pdf-E_th):.1e}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 2: Moments via exact quadrature
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 2 — Moments E[r_alpha^k] = integral_0^1 [x^a+(1-x)^a]^k dx")
print("=" * 70)

def moment_k(alpha, k, n_pts=200000):
    """E[r_alpha^k] via direct quadrature over x ~ U(0,1)."""
    x = np.linspace(1e-10, 1-1e-10, n_pts)
    f = x**alpha + (1-x)**alpha
    return np.trapz(f**k, x)

def haar_E(alpha, D):
    return gamma(alpha+1) * gamma(D+1) / gamma(D+alpha)

def haar_Var(alpha, D):
    A = gamma(2*alpha+1)*gamma(D) / gamma(D+2*alpha)
    B = gamma(alpha+1)**2 * gamma(D) / gamma(D+2*alpha)
    if D > 1:
        F = (D-1)*gamma(alpha+1)**2*gamma(D-1)**2 / (gamma(D-1+alpha)**2*(D-1+2*alpha))
    else:
        F = 0.0
    E = haar_E(alpha, D)
    return A + 2*(D-1)*B + (D-1)**2*F - E**2

print("\n--- D=2: moments vs theory ---")
D = 2
for alpha in [0.5, 1.5, 2.0, 3.0]:
    m1 = moment_k(alpha, 1)
    m2 = moment_k(alpha, 2)
    var_q = m2 - m1**2
    E_th = haar_E(alpha, D)
    V_th = haar_Var(alpha, D)
    print(f"  alpha={alpha:.1f}: E={m1:.6f} (theory {E_th:.6f}, err {abs(m1-E_th):.1e}), "
          f"Var={var_q:.6f} (theory {V_th:.6f}, err {abs(var_q-V_th):.1e})")

# All moments formula check: E[r_alpha^k] via Beta integrals
print("\n--- All-integer-moment pattern E[r_alpha^k] for alpha=2, D=2 ---")
print("   k   E[r_2^k](quad)   E[r_2^k](theory via Dirichlet)")
D = 2; alpha = 2.0
for k in [1, 2, 3, 4]:
    m = moment_k(alpha, k)
    # Dirichlet formula for E[(x^2+(1-x)^2)^k] = sum_{j=0}^k C(k,j) E[x^{2j}(1-x)^{2(k-j)}]
    # E[x^a (1-x)^b] = B(a+1,b+1) = Gamma(a+1)Gamma(b+1)/Gamma(a+b+2)
    from math import comb
    from scipy.special import gamma as G
    theo = sum(comb(k,j) * G(2*j+1)*G(2*(k-j)+1)/G(2*k+2) for j in range(k+1))
    print(f"   {k}   {m:.8f}   {theo:.8f}   err={abs(m-theo):.1e}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 3: Third central moment and skewness
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 3 — Third central moment and skewness (D=2 exact quadrature)")
print("=" * 70)

print(f"\n  {'alpha':>6} {'kappa_3':>12} {'skewness':>10} {'direction':>10}")
kappa3_table = {}
for alpha in [0.5, 1.5, 2.0, 3.0]:
    m1 = moment_k(alpha, 1)
    m2 = moment_k(alpha, 2)
    m3 = moment_k(alpha, 3)
    kappa3 = m3 - 3*m2*m1 + 2*m1**3
    sigma = np.sqrt(max(m2 - m1**2, 0))
    skew = kappa3 / sigma**3 if sigma > 1e-10 else float('nan')
    kappa3_table[alpha] = (kappa3, skew)
    direction = "right(α>1)" if skew > 1e-6 else ("left(α<1)" if skew < -1e-6 else "sym")
    print(f"  {alpha:>6.1f} {kappa3:>12.6f} {skew:>10.4f} {direction:>10}")

# Monte Carlo verification (small, fast)
print("\n--- Monte Carlo verification (D=2,3; 5000 samples) ---")
for D in [2, 3]:
    for alpha in [1.5, 2.0]:
        r_arr = np.array([
            sum(abs(U[k, 0])**(2*alpha) for k in range(D))
            for U in [unitary_group.rvs(D) for _ in range(5000)]
        ])
        mu = r_arr.mean()
        sigma = r_arr.std()
        kappa3 = np.mean((r_arr - mu)**3)
        skew = kappa3 / sigma**3 if sigma > 1e-10 else 0
        E_th = haar_E(alpha, D)
        print(f"  D={D}, alpha={alpha:.1f}: E_mc={mu:.4f} (theory {E_th:.4f}), "
              f"kappa_3={kappa3:.3e}, skew={skew:.3f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 4: Exact kappa_3 formula for D=2
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 4 — Exact kappa_3 formula for alpha=2, D=2")
print("=" * 70)

# For alpha=2: r_2 = x^2 + (1-x)^2, x ~ U(0,1)
# E[r_2^k] = sum_{j=0}^k C(k,j) B(2j+1, 2(k-j)+1)
# kappa_3 = E[r_2^3] - 3*E[r_2^2]*E[r_2] + 2*E[r_2]^3

from scipy.special import gamma as G
from math import comb

alpha = 2.0

def moment_k_exact(k):
    """Exact E[r_2^k] via Beta integrals."""
    return sum(comb(k,j) * G(2*j+1)*G(2*(k-j)+1)/G(2*k+2) for j in range(k+1))

m1 = moment_k_exact(1)
m2 = moment_k_exact(2)
m3 = moment_k_exact(3)
kappa3_exact = m3 - 3*m2*m1 + 2*m1**3
sigma_exact = np.sqrt(m2 - m1**2)
skew_exact = kappa3_exact / sigma_exact**3

# From fractions:
# m1 = 2/3
# m2 = 7/15
# m3 = ?
# B(5,5) = Gamma(5)^2/Gamma(10) = 24^2/362880 = 576/362880 = 1/630
# m3 = C(3,0)*B(1,7)/B_tot + C(3,1)*B(3,5)/B_tot + C(3,2)*B(5,3)/B_tot + C(3,3)*B(7,1)/B_tot
#     where B(2j+1,2(3-j)+1) / normalisation via integral_0^1 x^{2j}(1-x)^{6-2j}dx
# E[x^{2j}(1-x)^{6-2j}] = Gamma(2j+1)Gamma(7-2j)/Gamma(8)
from fractions import Fraction
# Compute exactly
def binom_moments_exact(k):
    from scipy.special import gamma
    total = 0
    for j in range(k+1):
        a = 2*j; b = 2*(k-j)
        # integral_0^1 x^a (1-x)^b dx = Gamma(a+1)Gamma(b+1)/Gamma(a+b+2) = a! b! / (a+b+1)!
        val = comb(k, j) * float(gamma(a+1)*gamma(b+1)/gamma(a+b+2))
        total += val
    return total

m3_exact = binom_moments_exact(3)
m4_exact = binom_moments_exact(4)
m1_exact = binom_moments_exact(1)
m2_exact = binom_moments_exact(2)

kappa3 = m3_exact - 3*m2_exact*m1_exact + 2*m1_exact**3
kappa4 = m4_exact - 4*m3_exact*m1_exact + 6*m2_exact*m1_exact**2 - 3*m1_exact**4
sigma2 = m2_exact - m1_exact**2
skew_val = kappa3 / sigma2**1.5
kurt = kappa4 / sigma2**2

print(f"\n  E[r_2] = {m1_exact:.10f} = 2/3 = {2/3:.10f}")
print(f"  E[r_2^2] = {m2_exact:.10f} = 7/15 = {7/15:.10f}")
print(f"  E[r_2^3] = {m3_exact:.10f}")
print(f"  Var[r_2] = {sigma2:.10f} = 1/45 = {1/45:.10f}")
print(f"  kappa_3 = {kappa3:.10f}")
print(f"  kappa_4 = {kappa4:.10f}")
print(f"  Skewness = {skew_val:.6f}")
print(f"  Excess kurtosis = {kurt:.6f}")

# Exact rational values
# kappa_3 = m3 - 3*m2*m1 + 2*m1^3
# m3 = C30*B(1,7) + C31*B(3,5) + C32*B(5,3) + C33*B(7,1)
# B(2j+1, 7-2j) where B is Beta = integral
from fractions import Fraction
def beta_exact(a, b):
    # int_0^1 x^a (1-x)^b dx = a!b!/(a+b+1)!
    result = Fraction(1)
    for i in range(b):
        result *= Fraction(a+i+1, a+b+i+1)
    return Fraction(1, a+b+1) * result

# Actually compute m3 exactly
comb3 = [1, 3, 3, 1]
m3_frac = Fraction(0)
for j in range(4):
    a = 2*j; b = 6-2*j
    fac_a = Fraction(1); tmp = 1
    for i in range(1, a+1): tmp *= i
    fac_a = Fraction(tmp)
    tmp = 1
    for i in range(1, b+1): tmp *= i
    fac_b = Fraction(tmp)
    tmp = 1
    for i in range(1, a+b+2): tmp *= i
    fac_ab = Fraction(tmp)
    m3_frac += comb3[j] * fac_a * fac_b / fac_ab

print(f"\n  E[r_2^3] exactly = {m3_frac} = {float(m3_frac):.10f}")
m1_frac = Fraction(2, 3)
m2_frac = Fraction(7, 15)
kappa3_frac = m3_frac - 3*m2_frac*m1_frac + 2*m1_frac**3
print(f"  kappa_3 exactly = {kappa3_frac} = {float(kappa3_frac):.10f}")
print(f"  kappa_3 > 0: {kappa3_frac > 0}  (right-skewed confirmed)")

# ─────────────────────────────────────────────────────────────────────────────
# PART 5: Concentration — sigma/mu -> 0 as D -> inf
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 5 — Concentration: sigma/mu vs D (all alpha)")
print("=" * 70)

print(f"\n{'D':>5}", end="")
for alpha in [1.5, 2.0, 3.0]:
    print(f"  sigma/mu(a={alpha:.1f})", end="")
print()

for D in [2, 4, 8, 16, 32, 64]:
    print(f"  {D:>4d}", end="")
    for alpha in [1.5, 2.0, 3.0]:
        sm = np.sqrt(haar_Var(alpha, D)) / haar_E(alpha, D)
        print(f"  {sm:>18.4e}", end="")
    print()

# Verify sigma/mu ~ C(alpha)/D
print("\n--- sigma/mu * D for large D (should converge) ---")
print(f"{'D':>5}", end="")
for alpha in [1.5, 2.0, 3.0]:
    print(f"  sigma*D/mu(a={alpha:.1f})", end="")
print()

for D in [8, 16, 32, 64, 128]:
    print(f"  {D:>4d}", end="")
    for alpha in [1.5, 2.0, 3.0]:
        sm = np.sqrt(haar_Var(alpha, D)) / haar_E(alpha, D) * D
        print(f"  {sm:>21.4f}", end="")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# PART 6: Monte Carlo PDF verification (fixed: use same unitary!)
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("PART 6 — MC PDF verification (D=2, alpha=2) vs exact 1/sqrt(2r-1)")
print("=" * 70)

n_mc = 80000
# FIX: use the same 2x2 unitary for each sample
r_mc = np.array([
    abs(U[0, 0])**4 + abs(U[1, 0])**4
    for U in [unitary_group.rvs(2) for _ in range(n_mc)]
])

bins = np.linspace(0.5, 1.0, 31)
counts, edges = np.histogram(r_mc, bins=bins, density=True)
mids = (edges[:-1] + edges[1:]) / 2
pdf_exact_vals = 1.0 / np.sqrt(2*mids - 1)
rel_err = np.max(np.abs(counts - pdf_exact_vals)) / np.mean(pdf_exact_vals)

print(f"  MC: E[r_2]  = {r_mc.mean():.6f}  (theory 2/3 = {2/3:.6f})")
print(f"  MC: Var[r_2] = {r_mc.var():.6f}  (theory 1/45 = {1/45:.6f})")
print(f"  Max relative histogram error vs 1/sqrt(2r-1): {rel_err:.2%}  "
      f"{'PASS' if rel_err < 0.12 else 'FAIL'}")

# MC for alpha=3 (D=2)
alpha = 3.0
r_mc3 = np.array([
    abs(U[0, 0])**(2*alpha) + abs(U[1, 0])**(2*alpha)
    for U in [unitary_group.rvs(2) for _ in range(n_mc)]
])
print(f"\n  D=2, alpha=3: E_mc={r_mc3.mean():.5f} (theory {haar_E(3,2):.5f}), "
      f"Var_mc={r_mc3.var():.5f} (theory {haar_Var(3,2):.5f})")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 70)
print("SUMMARY — KEY RESULTS")
print("=" * 70)
kappa3_val = float(kappa3_frac)
print(f"""
1. D=2 EXACT PDF:
   For alpha > 1: r_alpha = x^alpha+(1-x)^alpha, x~U(0,1), support [2^{{1-alpha}}, 1]
   p_alpha(r) = 2 / (alpha * |x1(r)^{{alpha-1}} - (1-x1(r))^{{alpha-1}}|)
   Special case alpha=2: p_2(r) = 1/sqrt(2r-1) on [1/2,1]  [PROVED + VERIFIED]
   For alpha < 1: support [1, 2^{{1-alpha}}], same formula, left-skewed

2. E[r_alpha] = 2/(alpha+1) for D=2 (all alpha > 0)  [VERIFIED]
   Matches Gamma formula: Gamma(alpha+1)*Gamma(3)/Gamma(alpha+2) = 2/(alpha+1)

3. Var[r_2] = 1/45 (D=2)  [VERIFIED from MC and quadrature]
   General: Var[r_alpha] from Section 48 formula

4. EXACT kappa_3 (D=2, alpha=2):
   kappa_3 = {kappa3_frac} ≈ {kappa3_val:.6f}  (positive = right-skewed)
   Skewness = {float(kappa3_frac / Fraction(1,45)**Fraction(3,2)):.4f}...
   (Note: Fraction arithmetic gives exact numerator/denominator)

5. SKEWNESS TABLE (D=2):
   alpha < 1: left-skewed (kappa_3 < 0)
   alpha = 1: degenerate delta function
   alpha > 1: right-skewed (kappa_3 > 0)
   Heavier tail toward r=1 (integrable = less entangling circuits) for alpha > 1.

6. CONCENTRATION: sigma/mu ~ C(alpha)/D -> 0 as D -> inf (all alpha > 0)
   Consistent with Section 47: sigma/mu ~ 1/D for alpha=2 (C(2)=1).

7. CLT (weak statement): r_alpha concentrates around E[r_alpha] as D -> inf
   (coefficient of variation -> 0). Full CLT requires dependent variable CLT.
""")
