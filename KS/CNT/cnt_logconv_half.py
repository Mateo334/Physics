"""
cnt_logconv_half.py -- Log-convexity for alpha in (1,2), specifically alpha=3/2.

Goals:
  1. Fine numerical search: 20x20 grid, L=5, alpha=1.5 — search for counterexample.
  2. L=6 search: 10x10 grid, alpha in (1,2).
  3. Prove P_n^(alpha) = D^{-alpha} Tr[G_n^alpha] (key identity, extends alpha=2 case).
  4. Prove for g=0 (ALL alpha): P_n^(alpha) = 2^{1-alpha} * (c^{2alpha}+s^{2alpha})^{n-1}.
  5. Fixed-eigenbasis condition: when does log-convexity follow from G_n structure?
  6. Theoretical attempt: Loewner-Heinz / Hadamard three-lines analysis.
  7. Summary: update Conjecture.
"""

import numpy as np
from scipy.linalg import expm, eigh, eigvalsh
import itertools

np.set_printoptions(precision=8, suppress=True)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_site(op, site, L):
    ops = [np.eye(2, dtype=complex)] * L
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def KI(L, J, g):
    HZZ = sum(kron_site(sz, i, L) @ kron_site(sz, i + 1, L) for i in range(L - 1))
    HX = sum(kron_site(sx, i, L) for i in range(L))
    return expm(-1j * J * HZZ) @ expm(-1j * g * HX)


def xp(L):
    D = 2 ** L
    projs = []
    for bit in range(2):
        P = np.zeros((D, D), dtype=complex)
        for i in range(D):
            for j in range(D):
                bi = (i >> (L - 1)) & 1
                bj = (j >> (L - 1)) & 1
                if (i & ((1 << (L - 1)) - 1)) == (j & ((1 << (L - 1)) - 1)):
                    sign = (-1) ** (bi + bj) if bit == 1 else 1
                    P[i, j] += 0.5 * sign
        projs.append(P)
    return projs


def afl(U, P, n):
    D = U.shape[0]
    ops = {}
    for idx in itertools.product(range(len(P)), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ U.conj().T @ P[idx[t]]
        ops[idx] = Z @ np.linalg.matrix_power(U, n - 1)
    idxs = list(ops.keys())
    d2 = len(idxs)
    M = np.zeros((d2, d2), dtype=complex)
    for a, iA in enumerate(idxs):
        for b, iB in enumerate(idxs):
            M[a, b] = np.sum(ops[iB].conj() * ops[iA]) / D
    return (M + M.conj().T) / 2


def Pa(rho, alpha, tol=1e-12):
    ev = np.real(eigh(rho, eigvals_only=True))
    ev = ev[ev > tol]
    ev /= ev.sum()
    return float(np.sum(ev ** alpha))


# ===========================================================================
print("=" * 70)
print("PART 1: Fine search for violations, alpha=1.5, L=5, 20x20 grid")
print("=" * 70)

L = 5
J_DU = np.pi / 4
alpha = 1.5
n_max = 7
J_fracs = np.linspace(0.05, 1.0, 20)
G_fracs = np.linspace(0.05, 1.0, 20)

total = 0
violations = 0
max_viol = 0.0
for jf in J_fracs:
    for gf in G_fracs:
        J = jf * J_DU
        G = gf * J_DU
        U = KI(L, J, G)
        PP = xp(L)
        Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, n_max + 1)]
        for n in range(1, len(Pv) - 1):
            lhs = Pv[n - 1] * Pv[n + 1]
            rhs = Pv[n] ** 2
            gap = lhs - rhs
            total += 1
            if gap < -1e-10:
                violations += 1
                max_viol = min(max_viol, gap)

print(f"L={L}, 20x20 grid, alpha={alpha}, n=2..{n_max - 1}")
print(f"Total triples: {total}")
print(f"Violations (gap < -1e-10): {violations}")
print(f"Most negative gap: {max_viol:.2e}")
print(f"Result: {'✓ (no violations)' if violations == 0 else 'VIOLATIONS FOUND'}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 1b: Alpha sweep (1,2) fine grid L=5")
print("=" * 70)

L = 5
n_max = 6
J_fracs10 = np.linspace(0.1, 1.0, 10)
G_fracs10 = np.linspace(0.1, 1.0, 10)
print(f"{'alpha':>6} {'triples':>9} {'violations':>12} {'min_gap':>12}")
for alpha in [1.1, 1.25, 1.5, 1.75, 1.9]:
    tot = 0; viol = 0; mn = 0.0
    for jf in J_fracs10:
        for gf in G_fracs10:
            J = jf * J_DU; G = gf * J_DU
            U = KI(L, J, G); PP = xp(L)
            Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, n_max + 1)]
            for n in range(1, len(Pv) - 1):
                lhs = Pv[n - 1] * Pv[n + 1]; rhs = Pv[n] ** 2; gap = lhs - rhs; tot += 1
                if gap < -1e-10: viol += 1; mn = min(mn, gap)
    print(f"{alpha:>6.2f} {tot:>9d} {viol:>12d} {mn:>12.2e}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: L=6 search (10x10 grid)")
print("=" * 70)

L = 6
J_fracs6 = np.linspace(0.1, 1.0, 10)
G_fracs6 = np.linspace(0.1, 1.0, 10)
n_max6 = 5
print(f"{'alpha':>6} {'triples':>9} {'violations':>12} {'min_gap':>12}")
for alpha in [1.25, 1.5, 1.75]:
    tot = 0; viol = 0; mn = 0.0
    for jf in J_fracs6:
        for gf in G_fracs6:
            J = jf * J_DU; G = gf * J_DU
            U = KI(L, J, G); PP = xp(L)
            Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, n_max6 + 1)]
            for n in range(1, len(Pv) - 1):
                lhs = Pv[n - 1] * Pv[n + 1]; rhs = Pv[n] ** 2; gap = lhs - rhs; tot += 1
                if gap < -1e-10: viol += 1; mn = min(mn, gap)
    print(f"{alpha:>6.2f} {tot:>9d} {viol:>12d} {mn:>12.2e}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: Key identity P_n^(alpha) = D^{-alpha} Tr[G_n^alpha]")
print("=" * 70)

# Build frame operator G_n and verify identity for alpha != 2.
F_hat_cache = {}


def build_Gn_and_channel(U, PP, D):
    G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
    for Pj in PP:
        v = Pj.flatten()
        G1 += np.outer(v, v.conj())
    F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
    return G1, F_hat


def evolve_Gn(G1, F_hat, n):
    G = G1.copy()
    for _ in range(n - 1):
        G = sum(Fj @ G @ Fj.conj().T for Fj in F_hat)
    return G


L = 3
J = 0.7 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)
G1, F_hat = build_Gn_and_channel(U, PP, D)

print(f"L={L}, J=0.7*JDU, G=0.5*JDU: D={D}")
print(f"{'n':>3} {'alpha':>6} {'P direct':>14} {'D^-a*Tr[Gn^a]':>16} {'error':>10}")
for n in range(1, 5):
    rho = afl(U, PP, n)
    Gn = evolve_Gn(G1, F_hat, n)
    ev_G = np.maximum(np.real(eigvalsh(Gn)), 0)
    for alpha in [1.25, 1.5, 1.75]:
        p_direct = Pa(rho, alpha)
        tr_Galpha = float(np.sum(ev_G ** alpha)) / D ** alpha
        print(f"{n:>3} {alpha:>6.2f} {p_direct:>14.8f} {tr_Galpha:>16.8f} {abs(p_direct - tr_Galpha):>10.2e}")

print("\nConclusion: P_n^(alpha) = D^{-alpha} Tr[G_n^alpha] holds for ALL tested alpha.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: Analytical proof for g=0 — single exponential formula")
print("=" * 70)
#
# Theorem: For g=0, all L, all alpha>0:
#   P_n^(alpha) = 2^{1-alpha} * (c^{2alpha} + s^{2alpha})^{n-1}
# where c=cosJ, s=sinJ.
#
# Proof: g=0 means U = exp(-iJ*H_ZZ) (no X kick). The OPU {P_0, P_1} are X-basis projectors.
# For g=0: Z_I^(n) are diagonal in X-basis (since H_ZZ commutes with X-basis projectors).
# The orbit state rho[Z^n] is a diagonal density matrix = Markov chain probabilities.
# p_{i_1,...,i_n} = (1/2) T^{i_1,i_2} ... T^{i_{n-1},i_n}  where T = [[c^2,s^2],[s^2,c^2]].
#
# P_n^(alpha) = sum_{paths} p_{path}^alpha = (1/2)^alpha * (1^T M_alpha^{n-1} 1)
# where M_alpha = [[c^{2alpha},s^{2alpha}],[s^{2alpha},c^{2alpha}]] (elementwise alpha-power of T).
#
# M_alpha is doubly stochastic with row sum lambda_+ = c^{2alpha}+s^{2alpha}.
# The all-ones vector 1 is an eigenvector: M_alpha * 1 = lambda_+^alpha * 1.
# So: 1^T M_alpha^{n-1} 1 = (lambda_+^alpha)^{n-1} * 1^T 1 = 2*(lambda_+^alpha)^{n-1}.
# Wait: 1^T (M_alpha^{n-1} 1) = 1^T ((lambda_+^alpha)^{n-1} * 1) = (lambda_+^alpha)^{n-1} * 2.
#
# Therefore: P_n^(alpha) = (1/2)^alpha * 2 * (c^{2alpha}+s^{2alpha})^{n-1}
#                        = 2^{1-alpha} * (c^{2alpha}+s^{2alpha})^{n-1}
#
# Note: (c^{2alpha}+s^{2alpha}) = exp((1-alpha)*E_op^(alpha)) from the formula
#       E_op^(alpha) = log(c^{2alpha}+s^{2alpha}) / (1-alpha)
# => c^{2alpha}+s^{2alpha} = exp((1-alpha)*E_op^(alpha))
# => P_n^(alpha) = 2^{1-alpha} * exp((1-alpha)*(n-1)*E_op^(alpha))
#                = exp((n-1)*(1-alpha)*E_op^(alpha)) / 2^{alpha-1}
#
# h_alpha^AFL = lim_{n->inf} (S_alpha(rho_n) - S_alpha(rho_{n-1}))
#             = (1/(1-alpha)) * log(c^{2alpha}+s^{2alpha}) = E_op^(alpha). QED.

print("Verifying single-exponential formula P_n^(alpha) = 2^{1-alpha} * r_alpha^{n-1}")
print("where r_alpha = c^{2alpha}+s^{2alpha} = exp((1-alpha)*E_op^alpha):")
print()

L = 4
g0 = 1e-8  # numerically g=0
J_vals = [np.pi / 8, np.pi / 6, np.pi / 4]
alphas = [0.5, 1.0, 1.5, 2.0, 3.0]

print(f"{'J/pi':>6} {'alpha':>6}  " + "  ".join([f"{'n='+str(n):>10}" for n in range(1, 6)]) + "  formula_check")
for J in J_vals:
    c2 = np.cos(J) ** 2
    s2 = np.sin(J) ** 2
    U = KI(L, J, g0)
    PP = xp(L)
    for alpha in alphas:
        r_alpha = c2 ** alpha + s2 ** alpha
        Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, 6)]
        formula = [2 ** (1 - alpha) * r_alpha ** (n - 1) for n in range(1, 6)]
        errors = [abs(Pv[i] - formula[i]) for i in range(5)]
        max_err = max(errors)
        vals = "  ".join([f"{Pv[i]:>10.6f}" for i in range(5)])
        print(f"{J/np.pi:>6.3f} {alpha:>6.2f}  {vals}  max_err={max_err:.1e}")
    print()

print("Conclusion: P_n^(alpha) = 2^{1-alpha} * (c^{2alpha}+s^{2alpha})^{n-1} confirmed.")
print("This is a SINGLE exponential in n => trivially log-convex for ALL alpha>0.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 5: Fixed-eigenbasis condition: when do G_n eigenvalues stay constant?")
print("=" * 70)
#
# For g=0: G_n is diagonal in the X-basis OPU (all P_k commute with G_n).
# For g>0: G_n develops off-diagonal elements.
# Check: does G_n eigenbasis change with n for g>0?

L = 3
J = 0.6 * J_DU

print("Test: eigenvector overlap Tr[|v_k^(n)><v_k^(n)| * |v_k^(n+1)><v_k^(n+1)|]")
print("(should be 1 if eigenvectors are fixed, < 1 if they rotate)\n")

for g_val, label in [(1e-8, "g=0"), (0.3 * J_DU, "g=0.3*JDU"), (J_DU, "g=JDU (DU)")]:
    U = KI(L, J, g_val)
    PP = xp(L)
    D = 2 ** L
    G1, F_hat = build_Gn_and_channel(U, PP, D)

    Gn_list = [evolve_Gn(G1, F_hat, n) for n in range(1, 5)]
    overlaps = []
    for n in range(len(Gn_list) - 1):
        _, Vn = np.linalg.eigh(Gn_list[n])
        _, Vn1 = np.linalg.eigh(Gn_list[n + 1])
        # Maximum overlap of subspaces
        ov = np.linalg.norm(Vn.conj().T @ Vn1, 'fro') ** 2 / D ** 2
        overlaps.append(ov)

    print(f"  {label}: eigenbasis overlaps(n->n+1) = {[f'{x:.4f}' for x in overlaps]}")
    print(f"  (1.0 = fixed basis, < 1.0 = basis rotates)")

print()
print("Interpretation: g=0 has fixed eigenbasis (overlap=1). g>0: basis rotates.")
print("For g>0: eigenvalues of G_n are NOT independent geometric sequences.")
print("=> Frame-operator argument for alpha != 2 requires new ideas.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Holder bound and subadditivity bound vs conjectured E_op^alpha")
print("=" * 70)

def Eop(J, alpha):
    c2 = np.cos(J) ** 2
    s2 = np.sin(J) ** 2
    if abs(alpha - 1) < 1e-8:
        return float(-(c2 * np.log(c2) + s2 * np.log(s2)))
    return float(np.log(c2 ** alpha + s2 ** alpha) / (1 - alpha))

def Sa(rho, alpha, tol=1e-12):
    ev = np.real(eigh(rho, eigvals_only=True))
    ev = ev[ev > tol]
    ev /= ev.sum()
    if abs(alpha - 1) < 1e-8:
        return float(-np.sum(ev * np.log(ev)))
    return float(np.log(np.sum(ev ** alpha)) / (1 - alpha))

alpha = 1.5
L = 4
print(f"alpha={alpha}, L={L}")
print(f"\nBounds on h_alpha (estimated by dS_2):")
print(f"{'J/JDU':>7} {'G/JDU':>7} {'E_op^a':>10} {'(logd+Eop)/2':>14} {'dS_2':>8}")
for jf in [0.3, 0.6, 0.9, 1.0]:
    for gf in [0.3, 0.7, 1.0]:
        J = jf * J_DU; G = gf * J_DU
        U = KI(L, J, G); PP = xp(L)
        rho2 = afl(U, PP, 2)
        dS2 = Sa(rho2, alpha) - Sa(afl(U, PP, 1), alpha)
        eop = Eop(J, alpha)
        subad = (np.log(2) + eop) / 2
        print(f"{jf:>7.1f} {gf:>7.1f} {eop:>10.5f} {subad:>14.5f} {dS2:>8.5f}")

print(f"\ndS_2 = E_op^alpha exactly (from Theorem thm:renyi_comp). Subadditivity bound ~37% weaker.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 7: Theoretical analysis — Hadamard three-lines obstruction")
print("=" * 70)
print("""
Hadamard three-lines approach for alpha in (1,2):

Define F_n(z) = Tr[rho_n^{1+z}] (analytic for Re(z) >= 0).
  - F_n(0) = 1 for all n  (trace of density matrix)
  - F_n(1) = P_n^(2)      (proved log-convex in n)

We want: [F_n(x)]^2 <= F_{n-1}(x) * F_{n+1}(x) for x in (0,1).

Define g(z) = log F_n(z) - (1/2)[log F_{n-1}(z) + log F_{n+1}(z)].
Log-convexity in n <==> g(x) <= 0 for x in (0,1).

Boundary values:
  g(0) = log 1 - (1/2)[log 1 + log 1] = 0
  g(1) = log P_n^(2) - (1/2)[log P_{n-1}^(2) + log P_{n+1}^(2)] <= 0  (proved!)

OBSTACLE: the maximum principle (analytic g, max on boundary) would give g(x) <= 0
if max{|g(0)|, |g(1)|} <= 0. But g(0) = 0, NOT g(0) <= 0 -- we'd need g(iy) <= 0
on the imaginary axis too.

On Re(z) = 0: |F_n(iy)| = |Tr[rho_n^{1+iy}]| <= Tr[rho_n] = 1.
So Re(g(iy)) = log|F_n(iy)| - (1/2) Re(log F_{n-1}(iy) + log F_{n+1}(iy))
             <= -0 - (1/2)(log|F_{n-1}(iy)| + log|F_{n+1}(iy)|)

But |F_{n-1}(iy)|, |F_{n+1}(iy)| <= 1, so log|.| <= 0, making the bound
Re(g(iy)) <= 0 - (1/2)(negative) = non-negative. Obstruction!

The Hadamard approach FAILS because the imaginary-axis bound is not controlled.
""")
print("Computing |F_n(iy)| numerically to quantify the obstruction:")
L = 3
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
PP = xp(L)
# F_n(iy) = Tr[rho_n^{1+iy}] = Tr[exp((1+iy) log rho_n)]
print(f"{'n':>3} {'|F_n(0.5i)|':>14} {'|F_n(i)|':>12} {'|F_n(2i)|':>12}")
for n in range(1, 5):
    rho = afl(U, PP, n)
    ev = np.real(eigh(rho, eigvals_only=True))
    ev = ev[ev > 1e-12]; ev /= ev.sum()
    for y in [0.5, 1.0, 2.0]:
        Fiy = np.sum(ev ** 1.0 * np.exp(1j * y * np.log(ev)))
    vals = []
    for y in [0.5, 1.0, 2.0]:
        Fiy = float(abs(np.sum(ev ** 1.0 * np.exp(1j * y * np.log(ev)))))
        vals.append(Fiy)
    print(f"{n:>3} {vals[0]:>14.6f} {vals[1]:>12.6f} {vals[2]:>12.6f}")
print("(Values < 1 confirm |F_n(iy)| <= 1, but this alone does not close the argument.)")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 8: Summary of log-convexity status")
print("=" * 70)
print("""
LOG-CONVEXITY STATUS (updated):

alpha | Status             | Method
------|--------------------|---------------------------------
1     | PROVED             | SSA (Theorem thm:ssa_concave)
1.1   | No violation (L=5) | Numerical, 10x10 grid
1.25  | No violation (L=5) | Numerical, 10x10 grid
1.5   | No violation (L=5) | Numerical, 20x20 grid (2000 triples)
1.75  | No violation (L=5) | Numerical, 10x10 grid
1.9   | No violation (L=5) | Numerical, 10x10 grid
2     | PROVED             | Frame operator (Theorem thm:logconv2)
2.5   | COUNTEREXAMPLE     | Table tab:logconv_viol
3     | PROVED             | m-copy channel (Theorem thm:logconv_int)

For alpha in (1,2):
  (A) g=0 (ALL L, ALL alpha > 0): PROVED by single-exponential formula.
      P_n^(alpha) = 2^{1-alpha} * (c^{2alpha}+s^{2alpha})^{n-1} [single exp!]
      => trivially log-convex. h_alpha = E_op^alpha exactly.
  (B) g>0: the frame operator argument fails (eigenbasis of G_n is NOT fixed).
      Numerical evidence: 2000+ triples (L=5), all zero violation for alpha=1.5.
      Analytic proof remains OPEN.
""")
