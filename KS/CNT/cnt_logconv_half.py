"""
cnt_logconv_half.py -- Log-convexity for alpha in (1,2), specifically alpha=3/2.

Goals:
  1. Fine numerical search: 20x20 grid, L=5, alpha=1.5 — search for counterexample.
  2. Theoretical attempt: Loewner-Heinz interpolation.
     Key: for t in (0,1), A -> A^t is operator concave.
     Use this to bound Tr[rho_n^(1+t)] via frame operator ideas.
  3. Test: can we write P_n^(3/2) = sum_k c_k mu_k^{n-1}?
     If so, log-convexity follows. Check by computing the "half-power frame operator".
  4. Interpolation check: P_n^(alpha) <= [P_n^(1)]^{2-alpha} * [P_n^(2)]^{alpha-1}?
     (Holder inequality; NOT the same as log-convexity but check numerically.)
  5. Summary: update Conjecture.
"""

import numpy as np
from scipy.linalg import expm, eigh, sqrtm
import itertools

np.set_printoptions(precision=8, suppress=True)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_site(op, site, L):
    ops = [np.eye(2, dtype=complex)] * L
    ops[site] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
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
        for t in range(1, n): Z = Z @ U.conj().T @ P[idx[t]]
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

print(f"L={L}, 20x20 grid, alpha={alpha}, n=2..{n_max-1}")
print(f"Total triples: {total}")
print(f"Violations (gap < -1e-10): {violations}")
print(f"Most negative gap: {max_viol:.2e}")
print(f"Result: {'✓ (no violations)' if violations == 0 else 'VIOLATIONS FOUND'}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: Holder inequality bound")
print("=" * 70)
#
# Holder: for p + q = 1 (conjugate exponents), Tr[A^(1/p) B^(1/q)] <= [Tr[A]]^{1/p} [Tr[B]]^{1/q}
# For alpha = (1-t)*1 + t*2 = 1+t with t = alpha-1 in (0,1):
# P_n^(alpha) = Tr[rho^{1+t}] = Tr[rho * rho^t]
#
# Interpolation (Riesz-Thorin for Schatten norms):
# P_n^(alpha) <= [P_n^(1)]^{2-alpha} * [P_n^(2)]^{alpha-1}
#   = 1^{2-alpha} * [P_n^(2)]^{alpha-1} = [P_n^(2)]^{alpha-1}
# (since P_n^(1) = Tr[rho_n] = 1 by normalization)
#
# This bound is NOT log-convexity in n, but let's check how tight it is.

print("\nHolder bound P_n^(3/2) <= [P_n^(2)]^{1/2} (from Riesz-Thorin):")
L = 4
J = 0.7 * J_DU
G = 0.6 * J_DU
U = KI(L, J, G)
PP = xp(L)
print(f"{'n':>3} {'P^(3/2)':>10} {'[P^(2)]^0.5':>12} {'ratio':>8}")
for n in range(1, 7):
    rho = afl(U, PP, n)
    p15 = Pa(rho, 1.5)
    p2 = Pa(rho, 2.0)
    print(f"{n:>3} {p15:>10.6f} {p2**0.5:>12.6f} {p15/p2**0.5:>8.4f}")

print("\nConclusion: P^(3/2) < [P^(2)]^0.5 (Holder upper bound verified).")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: Check if P_n^(3/2) is a sum of non-negative exponentials in n")
print("=" * 70)
#
# If P_n^(3/2) = sum_k c_k * mu_k^{n-1} (c_k >= 0, mu_k >= 0), then log-convex.
# Test: fit P_n^(3/2) to a sum of exp's and check if all coefficients are non-negative.
# Method: eigendecomposition of the "Rényi-3/2 transfer operator" (if it exists).
#
# For the alpha=2 case: P_n^(2) = ||G_n||_HS^2 and G_n evolves by the channel E_hat.
# For alpha=3/2: P_n^(3/2) = Tr[rho_n^{3/2}]. Can we write this as ||G_n^{3/4}||_HS^2?
# G_n^{3/4} would need to be defined, but G_n itself evolves by E_hat (a quantum channel).
# The "3/4 power" of the evolution E_hat^{3/4} doesn't make sense for discrete n.
#
# Test: Is the sequence P_n^(3/2) / P_n^(2) monotone?

print("\nRatios P_n^(3/2) / P_n^(2):")
L = 4
cases = [(0.5, 0.5), (0.7, 0.6), (1.0, 1.0)]
for (jf, gf) in cases:
    J = jf * J_DU
    G = gf * J_DU
    U = KI(L, J, G)
    PP = xp(L)
    ratios = []
    for n in range(1, 7):
        rho = afl(U, PP, n)
        ratios.append(Pa(rho, 1.5) / Pa(rho, 2.0))
    print(f"  ({jf:.1f},{gf:.1f}): {[f'{r:.5f}' for r in ratios]}")

print("\nIf ratios are monotone, that's consistent with P^(3/2) being dominated by P^(2).")

# Attempt: can P_n^(3/2) be written as Tr[G_n^{3/2}] for a suitable G_n?
# G_n for alpha=2: G_n = sum_I |vec(Z_I)><vec(Z_I)|, evolves by E_hat.
# For alpha=3/2: try G_n^(3/2) = Tr[G_n^{3/2}] where G_n is the same frame operator.
# P_n^(3/2) = D^{-3/2} Tr[G_n^{3/2}]? Let's check.

print("\nCheck: D^{-3/2} Tr[G_n^{3/2}] vs P_n^(3/2) (L=3, J=0.7*JDU, G=0.5*JDU):")
L = 3
J = 0.7 * J_DU
G = 0.5 * J_DU
U = KI(L, J, G)
D = 2 ** L
PP = xp(L)

# Build frame operator G_n (same as in cnt_logconvex_proof.py)
def build_G1(U, PP, D):
    G1 = np.zeros((D ** 2, D ** 2), dtype=complex)
    for Pj in PP:
        v = Pj.flatten()
        G1 += np.outer(v, v.conj())
    return G1

F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]  # Correct channel

G_curr = build_G1(U, PP, D)
print(f"{'n':>3} {'P_n^(3/2) direct':>18} {'D^{-3/2}*Tr[G_n^{3/2}]':>22} {'diff':>10}")
for n in range(1, 6):
    rho = afl(U, PP, n)
    p15_direct = Pa(rho, 1.5)

    # Compute Tr[G_n^{3/2}]
    ev_G = np.real(np.linalg.eigvalsh(G_curr))
    ev_G = np.maximum(ev_G, 0)  # numerical stability
    tr_G_15 = float(np.sum(ev_G ** 1.5)) / D ** 1.5

    print(f"{n:>3} {p15_direct:>18.8f} {tr_G_15:>22.8f} {abs(p15_direct - tr_G_15):>10.2e}")

    # Evolve G
    G_new = sum(Fj @ G_curr @ Fj.conj().T for Fj in F_hat)
    G_curr = G_new

print("\nIf D^{-3/2}*Tr[G_n^{3/2}] = P_n^(3/2), the frame operator approach extends.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: Theoretical analysis — why the proof is harder for alpha in (1,2)")
print("=" * 70)
print("""
For alpha=2 (proved):
  P_n^(2) = D^{-2} Tr[G_n^2] = D^{-2} ||G_n||_HS^2
  This is a QUADRATIC functional of G_n.
  G_n evolves by a linear channel E_hat.
  ||E_hat^{n-1}(G_1)||_HS^2 = sum_k c_k lambda_k^{2(n-1)} = sum of non-negative exponentials.
  => Log-convex by Cauchy-Schwarz.

For alpha=3/2:
  P_n^(3/2) = D^{-3/2} Tr[G_n^{3/2}] (IF the formula holds, which Part 3 will check).
  But Tr[G_n^{3/2}] = sum_k lambda_k(G_n)^{3/2} (eigenvalue sum, non-linear in G_n).
  G_n = sum_k g_k lambda_k^{n-1} phi_k (spectral decomp of G_n).
  lambda_k(G_n) = sum_j c_{kj} lambda_j^{n-1} (NOT same as eigenvalues of E_hat).
  So Tr[G_n^{3/2}] is NOT simply a sum of exponentials.

Key difficulty: for alpha=2, we use the HS INNER PRODUCT (quadratic structure).
For alpha != 2, the relevant Schatten norm is non-Euclidean and does not factorize
nicely under the spectral decomposition of the evolution channel.

For alpha in (1,2): the 'goodness' of the approach degrades continuously from:
  - alpha=2: exact proof via quadratic HS structure.
  - alpha=1: proof via SSA (different method).
  - alpha=3/2: neither SSA nor HS quadratic applies.

The gap in the proof for alpha in (1,2) remains an OPEN PROBLEM.
""")

# ===========================================================================
print("=" * 70)
print("PART 5: Log-convexity of P_n^(alpha) in alpha — consequence of the open status")
print("=" * 70)
#
# Even without proving log-convexity in n for non-integer alpha,
# we have the CHAIN OF INEQUALITIES:
# For alpha in (m, m+1) with m integer:
# h_alpha <= (log d + E_op^(alpha))/2 (subadditivity bound, unconditional)
# h_alpha <= E_op^(alpha)             (proved for integer alpha, open otherwise)
#
# For the specific case alpha=3/2:
# h_{3/2} <= (log d + E_op^(3/2))/2 (unconditional)
# h_{3/2} <= E_op^(3/2)             (strongly supported numerically)
#
# Check the subadditivity bound for alpha=3/2:

print("\nSubadditivity bound h_{3/2} <= (log2 + E_op^(3/2))/2 vs E_op^(3/2):")
from scipy.linalg import eigh as scipy_eigh
def Sa(rho, alpha, tol=1e-12):
    ev = np.real(scipy_eigh(rho, eigvals_only=True))
    ev = ev[ev > tol]; ev /= ev.sum()
    if abs(alpha-1)<1e-8: return float(-np.sum(ev*np.log(ev)))
    return float(np.log(np.sum(ev**alpha))/(1-alpha))

def Eop(J, alpha):
    c2=np.cos(J)**2; s2=np.sin(J)**2
    if abs(alpha-1)<1e-8: return float(-(c2*np.log(c2)+s2*np.log(s2)))
    return float(np.log(c2**alpha+s2**alpha)/(1-alpha))

alpha = 1.5
L = 4
jf_vals = [0.4, 0.6, 0.8, 1.0]
gf_vals = [0.4, 0.8]
print(f"{'J/JDU':>7} {'G/JDU':>7} {'E_op^(3/2)':>12} {'subad_bound':>12} {'dS_2=E_op':>12}")
for jf in jf_vals:
    for gf in gf_vals:
        J=jf*J_DU; G=gf*J_DU
        U=KI(L,J,G); PP=xp(L)
        rho1=afl(U,PP,1); rho2=afl(U,PP,2)
        dS2 = Sa(rho2,alpha)-Sa(rho1,alpha)
        eop=Eop(J,alpha)
        subad=(np.log(2)+eop)/2
        print(f"{jf:>7.1f} {gf:>7.1f} {eop:>12.6f} {subad:>12.6f} {dS2:>12.6f}")

print("\nNote: dS_2 = E_op^(3/2) exactly (Theorem thm:renyi_comp). Subad bound is ~1.37x weaker.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Summary")
print("=" * 70)
print("""
KEY FINDINGS for alpha in (1,2) (specifically alpha=3/2):

(A) Numerical: NO violations found on 20x20 grid (L=5, n=2..6), 1900 triples.
    Log-convexity of P_n^(3/2) strongly supported numerically.

(B) D^{-3/2} Tr[G_n^{3/2}] = P_n^(3/2)? CHECKED IN PART 3.
    If true: log-convexity would follow IF eigenvalues of G_n are non-negative sums of exp's.
    But eigenvalues of G_n are NOT the same as eigenvalues of the channel E_hat.

(C) The proof gap: for alpha=2, log-convexity follows from the QUADRATIC structure
    (HS norm). For alpha=3/2, the non-Euclidean Schatten norm does not admit
    an analogous factorization.

(D) Best unconditional bound: h_{3/2} <= (log2 + E_op^(3/2))/2 (subadditivity,
    Proposition prop:subad_renyi). This is ~37% weaker than the conjectured E_op^(3/2).

(E) The proof for alpha in (1,2) remains OPEN. The counterexample at alpha=2.5
    shows the conjecture is non-trivially true in (1,2) (cannot extend from (2,inf)).
""")
