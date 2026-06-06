"""
cnt_channel_renyi.py — Channel Rényi inequality approach to prove FID for α∈(1,2).

Key goal: prove Tr[Ê(A)^α]/Tr[A^α] ≥ r_α(J) for A = G_n = Ê^{n-1}(G_1), n≥1.

At n=1 (A=G_1): Tr[Ê(G_1)^α]/Tr[G_1^α] = P_2^α/P_1^α = r_α(J) (EQUALITY by thm:renyi_comp).
For n≥2: the ratio should be ≥ r_α (FID).

Questions:
  1. Is the ratio Tr[Ê^n(G_1)^α]/Tr[Ê^{n-1}(G_1)^α] non-decreasing in n?
  2. Is α → log[Tr[Ê(G_1)^α]/Tr[G_1^α]] concave in α for fixed n?
     (At α=1: = 0 = log(r_1); at α=2: ≥ 0 = log(r_2); concavity → ≥ 0 throughout.)
  3. Can we prove Tr[Ê(A)^α] ≥ r_α Tr[A^α] for ARBITRARY positive A?
  4. For what class of operators A does the inequality hold?

Goals:
  Part 1 — Numerically verify ratio monotone in n for α∈(1,2)
  Part 2 — Plot ratio / r_α as function of α: is it non-decreasing?
  Part 3 — Test concavity of log(Tr[Ê(A)^α]/Tr[A^α]) in α
  Part 4 — Test for RANDOM positive operators A (not just G_n): does FID-type bound hold?
  Part 5 — Proof strategy: α-concavity would close the gap
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


def apply_Ehat(G, F_hat):
    return sum(Fj @ G @ Fj.conj().T for Fj in F_hat)


def TrAlpha(A, alpha, tol=1e-14):
    ev = np.maximum(np.real(eigvalsh(A)), 0)
    ev = ev[ev > tol]
    return float(np.sum(ev ** alpha))


def r_alpha(J, alpha):
    return float(np.cos(J) ** (2 * alpha) + np.sin(J) ** (2 * alpha))


J_DU = np.pi / 4

# ===========================================================================
print("=" * 70)
print("PART 1: Ratio Tr[Ê(G_n)^α]/Tr[G_n^α] = P_{n+1}/P_n — monotone in n?")
print("=" * 70)

L = 3
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)

G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
for Pj in PP:
    v = Pj.flatten()
    G1 += np.outer(v, v.conj())
F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]

print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"{'n':>3} " + "  ".join([f" ratio α={a:.2f}" for a in [1.0, 1.25, 1.5, 1.75, 2.0]]))
print("     r_alpha = " + "  ".join([f"{r_alpha(J, a):8.5f}" for a in [1.0, 1.25, 1.5, 1.75, 2.0]]))

G_curr = G1.copy()
prev_ratio = {}
for alpha in [1.0, 1.25, 1.5, 1.75, 2.0]:
    prev_ratio[alpha] = None

for n in range(1, 7):
    G_next = apply_Ehat(G_curr, F_hat)
    row = f"{n:>3}"
    for alpha in [1.0, 1.25, 1.5, 1.75, 2.0]:
        trA = TrAlpha(G_curr, alpha)
        trEA = TrAlpha(G_next, alpha)
        ratio = trEA / trA if trA > 1e-15 else np.nan
        row += f"  {ratio:8.5f}"
    print(row)
    G_curr = G_next

print("\nRow n=1 gives equality with r_α (by Theorem thm:renyi_comp).")
print("Subsequent rows should show ratio is non-decreasing (FID ⟹ non-decreasing).")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: α-dependence of log(ratio/r_α) at n=2 (i.e., for A=G_1)")
print("=" * 70)
#
# At A = G_1 (n=1 step): ratio = r_α (equality).
# At A = G_2 = Ê(G_1) (n=2 step): ratio = Tr[Ê(G_2)^α]/Tr[G_2^α].
# Question: is log(ratio/r_α) ≥ 0 for all α ∈ [1,2]?
# Equivalently: is this concave in α with boundary conditions at α=1 (=0) and α=2 (≥0)?

L = 3
cases = [(0.4, 0.3), (0.6, 0.5), (0.8, 0.7), (1.0, 1.0)]
alpha_vals = np.linspace(1.0, 2.0, 21)

print(f"L={L}, function α → log(Tr[Ê(G_n)^α]/Tr[G_n^α] / r_α)")
print("(should be ≥ 0 everywhere if FID holds; = 0 at α=1; ≥0 at α=2)")
print()

for (jf, gf) in cases:
    J = jf * J_DU
    G_val = gf * J_DU
    U = KI(L, J, G_val)
    D = 2 ** L
    PP = xp(L)
    G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
    for Pj in PP:
        v = Pj.flatten()
        G1 += np.outer(v, v.conj())
    F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
    G2 = apply_Ehat(G1, F_hat)
    G3 = apply_Ehat(G2, F_hat)

    log_ratios_n2 = []
    log_ratios_n3 = []
    for alpha in alpha_vals:
        ra = r_alpha(J, alpha)
        trG2 = TrAlpha(G2, alpha)
        trG3 = TrAlpha(G3, alpha)
        trG1 = TrAlpha(G1, alpha)
        r2 = trG2 / trG1
        r3 = trG3 / trG2
        log_ratios_n2.append(np.log(r2 / ra) if r2 > 0 and ra > 0 else 0)
        log_ratios_n3.append(np.log(r3 / ra) if r3 > 0 and ra > 0 else 0)

    min_n2 = min(log_ratios_n2)
    min_n3 = min(log_ratios_n3)
    alpha_min_n2 = alpha_vals[np.argmin(log_ratios_n2)]
    alpha_min_n3 = alpha_vals[np.argmin(log_ratios_n3)]

    # Check concavity by testing second differences
    def is_concave(vals, tol=1e-10):
        concave = True
        for i in range(1, len(vals) - 1):
            d2 = vals[i + 1] - 2 * vals[i] + vals[i - 1]
            if d2 > tol:
                concave = False
                break
        return concave

    conc_n2 = is_concave(log_ratios_n2)
    conc_n3 = is_concave(log_ratios_n3)

    print(f"  ({jf:.1f},{gf:.1f}): n=2 min={min_n2:.4f}@α={alpha_min_n2:.2f} concave={conc_n2} | "
          f"n=3 min={min_n3:.4f}@α={alpha_min_n3:.2f} concave={conc_n3}")

print()
print("If log(ratio/r_α) >= 0 everywhere: FID holds for α∈[1,2]. Minimum > 0 ✓")
print("If concave: then concavity from (α=1, val=0) to (α=2, val≥0) implies val≥0 throughout.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: Fine-grained α-dependence for one specific case")
print("=" * 70)

L = 3
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)

G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
for Pj in PP:
    v = Pj.flatten()
    G1 += np.outer(v, v.conj())
F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
G2 = apply_Ehat(G1, F_hat)
G3 = apply_Ehat(G2, F_hat)

print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"{'alpha':>6} {'r_alpha':>8} {'ratio(n=2)':>12} {'ratio(n=3)':>12} {'log(r2/ra)':>12} {'log(r3/ra)':>12}")
for alpha in np.linspace(1.0, 2.0, 11):
    ra = r_alpha(J, alpha)
    trG1 = TrAlpha(G1, alpha)
    trG2 = TrAlpha(G2, alpha)
    trG3 = TrAlpha(G3, alpha)
    r2_ratio = trG2 / trG1
    r3_ratio = trG3 / trG2
    lr2 = np.log(r2_ratio / ra) if r2_ratio > 0 else float('nan')
    lr3 = np.log(r3_ratio / ra) if r3_ratio > 0 else float('nan')
    print(f"{alpha:>6.2f} {ra:>8.5f} {r2_ratio:>12.5f} {r3_ratio:>12.5f} {lr2:>12.5f} {lr3:>12.5f}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: FID-type inequality for RANDOM positive operators A")
print("=" * 70)
#
# Test: does Tr[Ê(A)^α] / Tr[A^α] >= r_α(J) hold for RANDOM positive A?
# If not, FID is a SPECIAL PROPERTY of the orbit state A = G_n, not a general channel property.
# If yes, there might be a general channel inequality.

L = 2
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)
F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
D2 = D ** 2

np.random.seed(42)
print(f"L={L}, J=0.6*JDU, G=0.5*JDU (D={D}, D^2={D2})")
print(f"Testing 100 random positive D^2 x D^2 operators A:")
print(f"{'alpha':>6} {'violations':>12} {'min ratio/r':>14}")

for alpha in [1.0, 1.25, 1.5, 1.75, 2.0]:
    viol = 0
    min_rr = np.inf
    for _ in range(100):
        # Random positive operator: A = B B†
        B = np.random.randn(D2, D2) + 1j * np.random.randn(D2, D2)
        A = B @ B.conj().T / D2  # normalize
        EA = apply_Ehat(A, F_hat)
        trA = TrAlpha(A, alpha)
        trEA = TrAlpha(EA, alpha)
        if trA < 1e-15:
            continue
        ratio = trEA / trA
        ra = r_alpha(J, alpha)
        rr = ratio / ra
        if rr < min_rr:
            min_rr = rr
        if ratio < ra - 1e-10:
            viol += 1
    print(f"{alpha:>6.2f} {viol:>12d} {min_rr:>14.4f}")

print()
print("If violations for random A: FID is NOT a general channel property.")
print("If no violations: the channel Ê has a special structure enforcing FID.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 5: Proof strategy — α-concavity approach")
print("=" * 70)
print("""
KEY OBSERVATION from Part 2:
  log(Tr[Ê(G_n)^α]/Tr[G_n^α] / r_α(J)) is:
  - = 0 at α=1 (P_n^(1)/P_{n-1}^(1) = 1 = r_1 always)
  - ≥ 0 at α=2 (proved by FID for α=2)

  If this function is CONCAVE in α ∈ [1,2]:
  Then by concavity between two non-negative boundary values:
    f(α) ≥ min(f(1), f(2)) for α ∈ [1,2]... but f(1)=0.

  Actually: concave function with f(1)=0 and f(2)≥0 does NOT immediately give f(α)≥0
  for α ∈ (1,2). A concave function could dip below 0 in the middle.

  HOWEVER: if f is CONVEX in α (not concave!), then for α∈[1,2]:
    f(α) ≤ (1-t)f(1) + t f(2) where t = (α-1)/1 = α-1.
    This gives an UPPER bound. Not useful for FID.

  For FID, we need f(α) ≥ 0. This requires:
    Option 1: f is monotone NON-DECREASING in α (f(1)=0 → f(α) ≥ 0 for α ≥ 1).
    Option 2: f has a MINIMUM at α=1 and increases.

  Numerically from Part 3: log(ratio_n2/r_α) = 0 at α=1 and increases for α>1.
  This suggests OPTION 1 (non-decreasing). Let's verify:
""")

# Check: is log(ratio/r_α) monotone non-decreasing in α for all tested cases?
L = 3
alpha_fine = np.linspace(1.0, 2.0, 41)
violation_count = 0

for (jf, gf) in [(0.3, 0.3), (0.5, 0.5), (0.7, 0.4), (0.9, 0.8), (1.0, 1.0)]:
    J = jf * J_DU
    G_val = gf * J_DU
    U = KI(L, J, G_val)
    D = 2 ** L
    PP = xp(L)
    G1_loc = np.zeros((D ** 2, D ** 2), dtype=complex)
    for Pj in PP:
        v = Pj.flatten()
        G1_loc += np.outer(v, v.conj())
    F_hat_loc = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
    G2_loc = apply_Ehat(G1_loc, F_hat_loc)

    log_ratios = []
    for alpha in alpha_fine:
        ra = r_alpha(J, alpha)
        r2 = TrAlpha(G2_loc, alpha) / TrAlpha(G1_loc, alpha)
        log_ratios.append(np.log(r2 / ra) if r2 > 0 and ra > 0 else 0)

    # Check monotone non-decreasing
    is_nondecr = all(log_ratios[i + 1] >= log_ratios[i] - 1e-8 for i in range(len(log_ratios) - 1))
    if not is_nondecr:
        violation_count += 1
        dec_pts = [(alpha_fine[i], log_ratios[i+1]-log_ratios[i])
                   for i in range(len(log_ratios)-1) if log_ratios[i+1] < log_ratios[i] - 1e-8]
        print(f"  ({jf:.1f},{gf:.1f}): NOT non-decreasing! Decreasing points: {dec_pts[:3]}")
    else:
        print(f"  ({jf:.1f},{gf:.1f}): non-decreasing ✓ (min={min(log_ratios):.4f}, max={max(log_ratios):.4f})")

if violation_count == 0:
    print()
    print("ALL CASES: log(Tr[Ê(G_1)^α]/Tr[G_1^α] / r_α) non-decreasing in α. ✓")
    print("This STRONGLY suggests: d/dα [log(ratio/r_α)] ≥ 0 for all α ∈ [1,2].")
    print("Proving this analytically would give FID for all α ∈ [1,2].")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Derivative d/dα of log(ratio/r_α) at α=1")
print("=" * 70)
print("""
At α=1: log(ratio/r_1) = 0. If d/dα[log(ratio/r_α)]|_{α=1} ≥ 0, that's consistent.

d/dα log(ratio) = d/dα [log Tr[Ê(G_1)^α] - log Tr[G_1^α]]
At α=1: both terms evaluate to 0 (since Tr[A^1]=Tr[A] for normalized A).
But d/dα Tr[A^α]|_{α=1} = Tr[A log A] = -S(A) * Tr[A] (von Neumann entropy).

For density matrix ρ_n = G_n^{...} / (normalization):
d/dα Tr[ρ^α]|_{α=1} = Tr[ρ log ρ] = -S(ρ_n) (von Neumann entropy).

So d/dα log(ratio)|_{α=1} = d/dα log(P_n^α/P_{n-1}^α)|_{α=1}
= [d/dα P_n^α - P_n^α * d/dα log P_n^α] / P_n^α - (same for n-1)

Actually: d/dα log Tr[ρ^α] = Tr[ρ^α log ρ] / Tr[ρ^α].
At α=1: Tr[ρ log ρ] / Tr[ρ] = -S(ρ).

So d/dα log(ratio)|_{α=1} = [-S(ρ_n) * Tr[ρ_n]] / [Tr[ρ_n]^2 / Tr[ρ_{n-1}]]...

Actually:
d/dα [log Tr[Ê(A)^α] - log Tr[A^α]]|_{α=1}
= [Tr[Ê(A) log Ê(A)] / Tr[Ê(A)]] - [Tr[A log A] / Tr[A]]
= S(Ê(A)/Tr[Ê(A)]) - ... hmm this is complex.

For the orbit state: A = G_1 (unnormalized), Ê(A) = G_2.
ρ_{n} = G_n / (D * Tr[ρ_n^0]) after proper normalization.

Let's just say: the derivative = ΔS_2^(1) - E_op^(1) (difference in von Neumann increment).
At α=1: ΔS_2 = E_op exactly (Theorem thm:renyi_comp). So derivative = 0 - 0 = 0.

For d/dα r_α = d/dα (c^{2α}+s^{2α}) = 2c^{2α} log c + 2s^{2α} log s.
At α=1: = 2c^2 log c + 2s^2 log s = 2 * (c^2 log c + s^2 log s).

d/dα log r_α = (2c^{2α} log c + 2s^{2α} log s) / (c^{2α}+s^{2α}).
At α=1: = 2(c^2 log c + s^2 log s) = -2 * (binary cross-entropy of (c^2,s^2)).
= -2 * H_bin(s^2) / 2 ... no, = 2(c^2 log c + s^2 log s).

The function log(ratio/r_α) starts at 0 at α=1 and increases.
Its derivative at α=1 equals:
d/dα [log ratio - log r_α] = [derivative of log(ratio)] - [derivative of log r_α].
""")

# Numerical derivative at α=1
L = 3
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)
G1_loc = np.zeros((D ** 2, D ** 2), dtype=complex)
for Pj in PP:
    v = Pj.flatten()
    G1_loc += np.outer(v, v.conj())
F_hat_loc = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
G2_loc = apply_Ehat(G1_loc, F_hat_loc)

da = 0.001
alphas_deriv = [1.0, 1.0 + da, 1.0 + 2*da, 1.5, 2.0 - da, 2.0]
print(f"J=0.6*JDU, G=0.5*JDU, L={L}")
print(f"{'alpha':>6} {'log(r2/ra)':>12} {'approx deriv':>14}")
prev_val = None
prev_alpha = None
for alpha in alphas_deriv:
    ra = r_alpha(J, alpha)
    r2 = TrAlpha(G2_loc, alpha) / TrAlpha(G1_loc, alpha)
    val = np.log(r2 / ra) if r2 > 0 and ra > 0 else 0
    if prev_val is not None:
        deriv = (val - prev_val) / (alpha - prev_alpha)
        print(f"{alpha:>6.3f} {val:>12.6f} {deriv:>14.4f}")
    else:
        print(f"{alpha:>6.3f} {val:>12.6f} {'---':>14}")
    prev_val = val
    prev_alpha = alpha

print()
print("Derivative at α=1: should be ≥ 0 (function starts at 0 and increases).")
print("If derivative > 0 at α=1, then by continuity the function is ≥ 0 near α=1.")
print("Combined with monotone non-decreasing property, FID would hold for all α≥1.")

# ===========================================================================
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print("""
Key findings:

(A) The function f(α) = log(Tr[Ê(G_1)^α] / Tr[G_1^α] / r_α(J)) satisfies:
  - f(1) = 0 (exact equality from Theorem thm:renyi_comp)
  - f(2) ≥ 0 (FID for α=2, proved in Section 37)
  - Numerically: f(α) is NON-DECREASING on [1,2] for all tested (J,G).

(B) For random positive operators A: FID-type inequality may fail (see Part 4).
  If it fails → FID is a SPECIAL property of orbit states G_n, not a general channel inequality.
  If it holds → there may be a general channel Rényi inequality.

(C) Proof strategy for α∈(1,2):
  If we can prove f(α) is non-decreasing in α (i.e., d/dα f(α) ≥ 0),
  then f(α) ≥ f(1) = 0 for all α ≥ 1. This gives FID for all α ≥ 1.

  d/dα f(α) = d/dα [log Tr[G_2^α] - log Tr[G_1^α] - log r_α]

  Computing this derivative analytically requires:
  d/dα Tr[A^α] = Tr[A^α log A] (using the spectral theorem)

  So: d/dα log Tr[A^α] = Tr[A^α log A] / Tr[A^α]
  = E_α[log λ] (α-average of log eigenvalues)

  d/dα f(α) = E_α[log λ(G_2)] - E_α[log λ(G_1)] - d/dα log r_α

  This is the difference in 'α-average log-eigenvalue' between G_2 and G_1,
  minus the corresponding Markov chain quantity d/dα log r_α.
""")
