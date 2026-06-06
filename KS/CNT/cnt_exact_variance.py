"""
cnt_exact_variance.py

Exact Var[r_2] using the S_4 Weingarten linear system to compute
F = E[|U_11|^4 |U_22|^4] for all D.

Key results:
  D=2: F = 1/5 = A (proved; |U_11|=|U_22| for 2x2 unitaries)
  D>=3: F from Weingarten linear system (numerical, exact for rational D)
  Var[r_2] = (A-mu^2) + 2(D-1)(B-mu^2) + (D-1)^2(F-mu^2)
"""

import numpy as np
from itertools import permutations

rng = np.random.default_rng(42)


def haar_unitary(D, rng):
    Z = (rng.standard_normal((D, D)) + 1j * rng.standard_normal((D, D))) / np.sqrt(2)
    Q, R = np.linalg.qr(Z)
    return Q * (np.diag(R) / np.abs(np.diag(R)))


def poch(a, k):
    return float(np.prod([a + i for i in range(int(k))]))


def n_cycles(perm):
    n = len(perm); vis = [False]*n; c = 0
    for i in range(n):
        if not vis[i]:
            c += 1; j = i
            while not vis[j]: vis[j] = True; j = perm[j]
    return c


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def cycle_type(perm):
    n = len(perm); vis = [False]*n; lens = []
    for i in range(n):
        if not vis[i]:
            l = 0; j = i
            while not vis[j]: vis[j] = True; j = perm[j]; l += 1
            lens.append(l)
    return tuple(sorted(lens, reverse=True))


all_perms = list(permutations(range(4)))
ct_to_perms = {}
for p in all_perms:
    ct = cycle_type(p)
    ct_to_perms.setdefault(ct, []).append(p)
ct_list = sorted(ct_to_perms.keys(), key=lambda c: sum(x-1 for x in c))
ct_idx = {ct: i for i, ct in enumerate(ct_list)}


def compute_F_wg(D):
    """Compute F = E[|U_11|^4|U_22|^4] from the Weingarten linear system for S_4 on U(D).
    M[i,j] = sum_{tau in C_j} D^{cycles(sigma_i * tau)}  (sum over ALL elements of class j)
    """
    reps = [ct_to_perms[ct][0] for ct in ct_list]
    nc = len(ct_list)
    M = np.zeros((nc, nc))
    for i, si in enumerate(reps):
        for j, ctj in enumerate(ct_list):
            # Sum over ALL elements of class j
            s = sum(D**n_cycles(compose(si, tau)) for tau in ct_to_perms[ctj])
            M[i, j] = s
    e0 = np.zeros(nc); e0[0] = 1.0
    try:
        wg = np.linalg.solve(M, e0)
    except np.linalg.LinAlgError:
        return np.nan
    # F = 4*Wg(id) + 8*Wg([21^2]) + 4*Wg([2^2])
    return 4*wg[ct_idx[(1,1,1,1)]] + 8*wg[ct_idx[(2,1,1)]] + 4*wg[ct_idx[(2,2)]]


def var_r2_exact(D, F_val):
    A = 24/poch(D,4); B = 4/poch(D,4)
    mu = 2.0/(D*(D+1))
    return (A-mu**2) + 2*(D-1)*(B-mu**2) + (D-1)**2*(F_val-mu**2)


# ============================================================
# Part 1: F values from Weingarten linear system
# ============================================================
print("="*60)
print("Part 1: F(D) from Weingarten linear system")
print("="*60)
F_by_D = {2: 24/poch(2,4)}  # D=2: F=A=1/5 (proved)
print(f"  D=2: F = A = 1/5 = {F_by_D[2]:.10f} (exact)")
for D in range(3, 12):
    F_by_D[D] = compute_F_wg(D)
    print(f"  D={D}: F = {F_by_D[D]:.10f}")

# ============================================================
# Part 2: Exact Var[r_2] vs Monte Carlo
# ============================================================
print()
print("="*60)
print("Part 2: Exact Var[r_2] from Weingarten F")
print("="*60)
N_mc = 500000
print(f"  {'D':>3}  {'Var_formula':>14}  {'Var_MC':>14}  {'error':>10}  {'std/mean':>10}")
for D in [2, 3, 4, 5, 6, 7, 8]:
    F_val = F_by_D[D]
    var_form = var_r2_exact(D, F_val)
    r2_vals = [np.sum(np.abs(haar_unitary(D,rng))**4)/D for _ in range(N_mc)]
    var_mc = np.var(r2_vals)
    mean_mc = np.mean(r2_vals)
    print(f"  {D:>3}  {var_form:>14.10f}  {var_mc:>14.10f}  {abs(var_form-var_mc):>10.2e}  "
          f"{np.sqrt(var_mc)/mean_mc:>10.6f}")

# ============================================================
# Part 3: F(D) as a fraction — check for rational form
# ============================================================
print()
print("="*60)
print("Part 3: Rational form of F(D)")
print("="*60)
print("  D=2: F = 1/5 (exact)")
print("  D=4: F = 1/84? Check:", 1/84)

# The formula F = (D^2+2) / (D*(D+1)*(D+2)*(D+3)) * something?
# Let's check F * poch(D,4):
print("\n  F * poch(D,4):")
for D in range(2, 12):
    val = F_by_D.get(D, np.nan) * poch(D,4)
    print(f"    D={D}: {val:.8f}")

# Try F * D(D+1)(D+2)(D+3)/2:
# D=4: 10/2 = 5
# D=5: 8.4/2 = 4.2 = 21/5
# Hmm.

# Try: F * poch(D,4) as polynomial in D
# D=2: 24, D=3: 33.75, D=4: 10, D=5: 8.4...
# Doesn't look polynomial. Let's fit as rational P(D)/Q(D).

# Check F * D^2 * (D^2-1):
print("\n  F * D^2*(D^2-1):")
for D in range(3, 12):
    val = F_by_D.get(D, np.nan) * D**2 * (D**2-1)
    print(f"    D={D}: {val:.8f}")

# Check F * D^2 * (D^2-1) * (D^2-4) / const:
# This is the structure from U(D) Weingarten with (D^2-1)(D^2-4) in denominator.
print("\n  F * D^2*(D^2-1)*(D^2-4):")
for D in range(4, 12):
    val = F_by_D.get(D, np.nan) * D**2 * (D**2-1) * (D**2-4)
    print(f"    D={D}: {val:.8f}")

# From the S_4 Weingarten, denominators involve D*(D^2-1)*(D^2-4)*(D+3) etc.
# Let's try D*(D^2-1)*(D^2-4)*(D+3):
print("\n  F * D*(D^2-1)*(D^2-4)*(D+3):")
for D in range(4, 12):
    val = F_by_D.get(D, np.nan) * D*(D**2-1)*(D**2-4)*(D+3)
    print(f"    D={D}: {val:.8f}")

# ============================================================
# Part 4: Exact Var[r_2] for D=2,3 analytically
# ============================================================
print()
print("="*60)
print("Part 4: Exact Var[r_2] for D=2 and D=3")
print("="*60)

# D=2: F = A = 24/120 = 1/5 (proved)
D = 2
A = 24/poch(2,4); B = 4/poch(2,4); mu = 2/(2*3); F = A
var = var_r2_exact(2, F)
print(f"  D=2: Var[r_2] = 7/15 - 4/9 = 1/45 = {1/45:.10f} (exact)")
print(f"       Formula: {var:.10f}")

# D=3: F from Weingarten (approximately 0.037037 = 1/27?)
D = 3
F3 = F_by_D[3]
print(f"\n  D=3: F from Weingarten = {F3:.10f}")
# Check if F3 = 1/27 = 0.037037
print(f"       1/27 = {1/27:.10f}, diff = {abs(F3-1/27):.2e}")
# Check F3 * poch(3,4) = F3 * 3*4*5*6 = F3*360
print(f"       F3 * 360 = {F3*360:.6f}")

# ============================================================
# Part 5: Summary table
# ============================================================
print()
print("="*60)
print("Part 5: Summary — Var[r_2] and F(D)")
print("="*60)
print(f"  {'D':>3}  {'F':>14}  {'Var[r_2]':>14}  {'std/mean':>10}")
for D in [2, 3, 4, 5, 6, 8, 10]:
    F_val = F_by_D.get(D, np.nan)
    if np.isnan(F_val): continue
    var_form = var_r2_exact(D, F_val)
    E_r2 = 2.0/(D+1)
    std = np.sqrt(max(var_form, 0))
    print(f"  {D:>3}  {F_val:>14.8f}  {var_form:>14.10f}  {std/E_r2:>10.6f}")

print()
print("Key exact values:")
print(f"  D=2: Var[r_2] = 1/45 = {1/45:.10f}")
D2_form = var_r2_exact(2, F_by_D[2])
print(f"       Formula gives: {D2_form:.10f}")
print(f"  Concentration: std/mean ~ O(1/sqrt(D)) -> 0")
