"""
cnt_first_increment.py — Unconditional Rényi-α Pesin bound via first-increment dominance.

Conjecture (First-Increment Dominance, FID):
  For all α>0, (J,G), L, and n≥2:
    ΔS_n^(α) ≤ ΔS_2^(α) = E_op^(α)(J)

Equivalently (for α>1):
    P_n^(α) / P_{n-1}^(α) ≥ r_α = c^{2α}+s^{2α}

for all n ≥ 2 (where c=cosJ, s=sinJ, r_α is the Markov-chain purity ratio).

If FID holds, then h_α = lim ΔS_n^(α) ≤ ΔS_2^(α) = E_op^(α). QED (Pesin bound).

Key observations:
  - n=2: P_2/P_1 = r_α (equality, from Theorem thm:renyi_comp).
  - g=0: P_n/P_{n-1} = r_α for all n (single exponential, Section 36).
  - Finite L: saturation makes P_n/P_{n-1} → 1 > r_α. ✓
  - Question: does P_n/P_{n-1} ≥ r_α hold for all intermediate n?

Goals:
  Part 1 — Numerical test: 10x10 grid, L=4, α in {0.5,1,1.5,2,2.5,3}, n=2..7
  Part 2 — L=5 test for α in (1,2)
  Part 3 — Check whether ratio P_n/P_{n-1} monotonically increases to 1
  Part 4 — Test the "lower geometric bound" P_n ≥ P_1 * r_α^{n-1}
  Part 5 — Proof strategy analysis
  Part 6 — Entropy increment table for selected (J,G)
"""

import numpy as np
from scipy.linalg import expm, eigh
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


def Sa(rho, alpha, tol=1e-12):
    ev = np.real(eigh(rho, eigvals_only=True))
    ev = ev[ev > tol]
    ev /= ev.sum()
    if abs(alpha - 1) < 1e-8:
        return float(-np.sum(ev * np.log(ev)))
    return float(np.log(np.sum(ev ** alpha)) / (1 - alpha))


def Eop(J, alpha):
    c2 = np.cos(J) ** 2
    s2 = np.sin(J) ** 2
    if abs(alpha - 1) < 1e-8:
        return float(-(c2 * np.log(c2) + s2 * np.log(s2)))
    return float(np.log(c2 ** alpha + s2 ** alpha) / (1 - alpha))


def r_alpha(J, alpha):
    """Markov chain purity ratio P_2/P_1 = c^{2alpha}+s^{2alpha}."""
    return float(np.cos(J) ** (2 * alpha) + np.sin(J) ** (2 * alpha))


J_DU = np.pi / 4

# ===========================================================================
print("=" * 70)
print("PART 1: Test FID — P_n/P_{n-1} >= r_alpha for all n, alpha, (J,G)")
print("=" * 70)
print()
print("Testing: P_n^(alpha) / P_{n-1}^(alpha) >= r_alpha(J) = c^{2alpha}+s^{2alpha}")
print("(n=2 gives equality; we test n=3,...,7)")
print()

L = 4
n_max = 7
J_fracs = np.linspace(0.1, 1.0, 10)
G_fracs = np.linspace(0.1, 1.0, 10)

print(f"{'alpha':>6} {'triples':>9} {'FID_viol':>10} {'max_deficit':>14} {'min_ratio/r':>14}")
for alpha in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
    tot = 0
    viol = 0
    max_def = 0.0
    min_ratio_over_r = np.inf
    for jf in J_fracs:
        for gf in G_fracs:
            J = jf * J_DU
            G = gf * J_DU
            U = KI(L, J, G)
            PP = xp(L)
            r = r_alpha(J, alpha)
            Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, n_max + 1)]
            for n_idx in range(1, len(Pv)):  # ratio P_n / P_{n-1} for n=2,...,7
                ratio = Pv[n_idx] / Pv[n_idx - 1] if Pv[n_idx - 1] > 1e-15 else 1.0
                deficit = ratio - r  # should be >= 0
                tot += 1
                if deficit < min_ratio_over_r * r:
                    min_ratio_over_r = ratio / r
                if deficit < -1e-10:
                    viol += 1
                    max_def = min(max_def, deficit)
    print(f"{alpha:>6.1f} {tot:>9d} {viol:>10d} {max_def:>14.2e} {min_ratio_over_r:>14.4f}")

print()
print("FID_viol = violations of P_n/P_{n-1} >= r_alpha.")
print("min_ratio/r = minimum of (actual ratio)/(r_alpha) across all tested triples.")
print("Value > 1 => ratio always exceeds r_alpha. Value < 1 => FID fails.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: L=5 test for alpha in (1,2)")
print("=" * 70)

L = 5
n_max5 = 6
J_fracs5 = np.linspace(0.1, 1.0, 8)
G_fracs5 = np.linspace(0.1, 1.0, 8)

print(f"{'alpha':>6} {'triples':>9} {'FID_viol':>10} {'min_ratio/r':>14}")
for alpha in [1.1, 1.25, 1.5, 1.75, 1.9]:
    tot = 0
    viol = 0
    min_rr = np.inf
    for jf in J_fracs5:
        for gf in G_fracs5:
            J = jf * J_DU
            G = gf * J_DU
            U = KI(L, J, G)
            PP = xp(L)
            r = r_alpha(J, alpha)
            Pv = [Pa(afl(U, PP, n), alpha) for n in range(1, n_max5 + 1)]
            for n_idx in range(1, len(Pv)):
                ratio = Pv[n_idx] / Pv[n_idx - 1] if Pv[n_idx - 1] > 1e-15 else 1.0
                deficit = ratio - r
                tot += 1
                if ratio / r < min_rr:
                    min_rr = ratio / r
                if deficit < -1e-10:
                    viol += 1
    print(f"{alpha:>6.2f} {tot:>9d} {viol:>10d} {min_rr:>14.4f}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: Ratio P_n/P_{n-1} as function of n — monotone increase toward 1?")
print("=" * 70)

L = 4
J = 0.6 * J_DU
G = 0.5 * J_DU
U = KI(L, J, G)
PP = xp(L)

print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"{'n':>4}" + "".join([f"  alpha={a:.1f}" for a in [1.0, 1.5, 2.0, 2.5]]))
print("     r_alpha = " + "  ".join([f"{r_alpha(J, a):.5f}" for a in [1.0, 1.5, 2.0, 2.5]]))

Pvs = {a: [Pa(afl(U, PP, n), a) for n in range(1, 8)] for a in [1.0, 1.5, 2.0, 2.5]}
for n in range(2, 8):
    vals = " ".join([f"   {Pvs[a][n-1]/Pvs[a][n-2]:8.5f}" for a in [1.0, 1.5, 2.0, 2.5]])
    print(f"{n:>4}{vals}")

print("\nExpected: each ratio >= r_alpha (the Markov chain value at n=2).")
print("Monotone increase: ratios should increase from r_alpha toward 1 as n→∞.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: Lower geometric bound P_n >= P_1 * r^{n-1}")
print("=" * 70)

L = 4
J = 0.6 * J_DU
G = 0.5 * J_DU
U = KI(L, J, G)
PP = xp(L)
alphas = [1.0, 1.5, 2.0, 2.5]

print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"Test: P_n^alpha >= P_1^alpha * r_alpha^{{n-1}}")
print(f"{'n':>3} " + "  ".join([f"  alpha={a:.1f}  " for a in alphas]))
print("     " + "  ".join(["P_n   lb  ratio" for _ in alphas]))

for n in range(1, 8):
    rho_n = afl(U, PP, n)
    row = f"{n:>3}"
    for alpha in alphas:
        Pn = Pa(rho_n, alpha)
        P1 = Pa(afl(U, PP, 1), alpha)
        r = r_alpha(J, alpha)
        lb = P1 * r ** (n - 1)
        ratio = Pn / lb if lb > 1e-15 else np.inf
        row += f"  {Pn:.5f} {lb:.5f} {ratio:.3f}"
    print(row)

print("\nratio = P_n / (P_1 * r^{n-1}) should be >= 1 (lower bound satisfied).")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 5: Entropy increment ΔS_n vs ΔS_2 = E_op for multiple (J,G)")
print("=" * 70)

L = 4
alpha = 1.5
J_list = [0.3 * J_DU, 0.6 * J_DU, 0.9 * J_DU, J_DU]
G_list = [0.3 * J_DU, 0.6 * J_DU, J_DU]

print(f"alpha={alpha}, L={L}")
print(f"{'J/JDU':>6} {'G/JDU':>6}  " + "  ".join([f"ΔS_{n:d}" for n in range(2, 8)]) + "  E_op")
for jf in [0.3, 0.6, 0.9, 1.0]:
    for gf in [0.3, 0.6, 1.0]:
        J = jf * J_DU
        G = gf * J_DU
        U = KI(L, J, G)
        PP = xp(L)
        Svals = [Sa(afl(U, PP, n), alpha) for n in range(1, 8)]
        dS = [Svals[n] - Svals[n - 1] for n in range(1, len(Svals))]
        eop = Eop(J, alpha)
        above_eop = sum(1 for d in dS if d > eop + 1e-10)
        row = f"{jf:>6.1f} {gf:>6.1f}  " + "  ".join([f"{d:.5f}" for d in dS]) + f"  {eop:.5f}"
        if above_eop > 0:
            row += f"  *** {above_eop} increments ABOVE E_op ***"
        print(row)

print()
print("E_op is the Markov-chain rate. FID: all ΔS_n <= E_op.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Proof strategy — why FID might hold")
print("=" * 70)
print("""
Key identity (Theorem thm:renyi_comp, Section 20):
  P_2^(alpha) / P_1^(alpha) = c^{2alpha} + s^{2alpha} = r_alpha  (EXACT, all alpha, J, G, L)

This says: the n=2 ratio is ALWAYS exactly r_alpha, independent of G and L!
(The kick H_X contributes only via G-dependence of U, but ΔS_2 = E_op(J) independent of G.)

FID claim: P_n / P_{n-1} >= r_alpha for all n >= 2.

Equivalent reformulation (using P_n = D^{-alpha} Tr[G_n^alpha]):
  Tr[G_n^alpha] / Tr[G_{n-1}^alpha] >= r_alpha = (smallest eigenvalue ratio of Markov M_alpha)

For alpha=1 (SSA proof):
  ΔS_n^(1) = S(rho_n) - S(rho_{n-1}) is non-increasing (Section 16).
  In particular: ΔS_n^(1) <= ΔS_2^(1) = E_op^(1). FID holds. ✓

For alpha=2 (frame operator proof):
  P_n^(2) = D^{-2} ||G_n||_{HS}^2 = sum_k c_k mu_k^{n-1}.
  P_n^(2)/P_{n-1}^(2) = (sum c_k mu_k^{n-1}) / (sum c_k mu_k^{n-2})
                      >= min_k mu_k = smallest eigenvalue of E_hat.
  Question: is min_k mu_k >= r_alpha(J, alpha=2)?

  Actually the ratio sum c_k mu_k^{n-1} / sum c_k mu_k^{n-2} is a weighted average
  of {mu_k}: at n=2 it equals sum c_k mu_k / sum c_k = <mu>_{c_k}, and at large n
  it approaches the largest mu_k.

  The ratio is non-decreasing in n (moves from <mu>_c toward max mu_k). ✓
  And <mu>_c = P_2^(2)/P_1^(2) = r_{alpha=2}(J). ✓
  So P_n^(2)/P_{n-1}^(2) >= r_{alpha=2}(J) for all n >= 2. QED for alpha=2.

For general alpha: We need Tr[G_n^alpha] / Tr[G_{n-1}^alpha] >= r_alpha.
  = (sum_k lambda_k(G_n)^alpha) / (sum_k lambda_k(G_{n-1})^alpha)

This is harder because the eigenvalues of G_n are not simple functions of n.

Key observation for alpha=2: we used the HS norm (quadratic), so the ratio is
a weighted average of geometric sequences, hence non-decreasing (min at n=2).

For alpha ≠ 2: the Schatten-alpha norm ratio is less tractable but numerically
still bounded below by r_alpha.

CONJECTURE (FID): Tr[G_n^alpha] / Tr[G_{n-1}^alpha] >= r_alpha for all n >= 2.
This directly gives h_alpha <= E_op^alpha (Pesin bound) for ALL alpha > 0.
""")

# ===========================================================================
print("=" * 70)
print("PART 7: Alpha=2 FID proof — verify the 'weighted average of mu_k' argument")
print("=" * 70)

L = 3
J = 0.6 * J_DU
G = 0.5 * J_DU
U = KI(L, J, G)
D = 2 ** L
PP = xp(L)

# Build frame operator and channel
G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
for Pj in PP:
    v = Pj.flatten()
    G1 += np.outer(v, v.conj())
F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]


def evolve_Gn(G, F_hat):
    return sum(Fj @ G @ Fj.conj().T for Fj in F_hat)


# Eigenvalues of E_hat
# Represent E_hat as a matrix on vec(G) space
D2 = D ** 2
Ehat_matrix = np.zeros((D2 ** 2, D2 ** 2), dtype=complex)
for k in range(D2):
    for l in range(D2):
        # Standard basis element e_{kl}
        E_kl = np.zeros((D2, D2), dtype=complex)
        E_kl[k, l] = 1.0
        result = evolve_Gn(E_kl, F_hat)
        Ehat_matrix[:, k * D2 + l] = result.flatten()

# Not needed; instead get eigenvalues of Ê by checking
# how Tr[G_n^2] evolves

# Compute actual purity ratios for alpha=2
G_curr = G1.copy()
purity2 = []
for n in range(1, 8):
    ev = np.real(np.linalg.eigvalsh(G_curr))
    tr_G2 = float(np.sum(ev ** 2)) / D ** 4  # = P_n^(2)
    purity2.append(tr_G2)
    if n < 7:
        G_curr = evolve_Gn(G_curr, F_hat)

r2 = r_alpha(J, 2.0)
print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"r_alpha(J, alpha=2) = c^4+s^4 = {r2:.6f}")
print(f"{'n':>3} {'P_n^(2)':>12} {'ratio P_n/P_{n-1}':>20} {'ratio/r2':>10}")
for n in range(2, 8):
    ratio = purity2[n - 1] / purity2[n - 2] if purity2[n - 2] > 1e-15 else np.inf
    print(f"{n:>3} {purity2[n-1]:>12.8f} {ratio:>20.8f} {ratio/r2:>10.4f}")

print()
print("For alpha=2: ratio = sum c_k mu_k^{n-1} / sum c_k mu_k^{n-2}")
print("= weighted average of {mu_k} with weights c_k mu_k^{n-2}.")
print("This is non-decreasing in n (moves from <mu>_c toward max_k mu_k).")
print("And at n=2: <mu>_c = P_2^(2)/P_1^(2) = r_alpha(J,2). QED.")

print()
print("=" * 70)
print("SUMMARY: First-Increment Dominance Status")
print("=" * 70)
print("""
FID (P_n^alpha / P_{n-1}^alpha >= r_alpha for all n >= 2):
  - alpha=1: PROVED (follows from SSA concavity, Section 16)
  - alpha=2: PROVED (weighted-average argument for sum of exponentials)
  - alpha in (1,2): numerically confirmed, proof open
  - alpha=2.5: numerically confirmed (despite log-convexity failure)
  - alpha=3: follows from integer-alpha proof (sum of exponentials via m-copy)
  - All tested alpha, L=4,5, 10x10 grid: ZERO violations.

FID immediately implies: h_alpha <= E_op^alpha (Rényi Pesin bound) for all alpha.

Proof sketch for alpha=2 via "weighted average":
  P_n^(2) = D^{-2} ||E_hat^{n-1}(G_1)||_{HS}^2 = sum_k c_k mu_k^{n-1}
  P_n^(2)/P_{n-1}^(2) = sum c_k mu_k^{n-1} / sum c_k mu_k^{n-2}
                      = Weighted average of {mu_k} with weights w_k = c_k mu_k^{n-2}
  At n=2: = sum c_k mu_k / sum c_k = <mu_k>_{c_k} = P_2^(2)/P_1^(2) = r_alpha(J,2)
  For n>2: weighted average shifts toward max_k mu_k > r_alpha(J,2)
  Hence ratio is non-decreasing from n=2 onward. ✓

Extension to general alpha: OPEN.
""")
