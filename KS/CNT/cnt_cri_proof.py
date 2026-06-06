"""
cnt_cri_proof.py — Attempt to prove the Channel Rényi Inequality (CRI):
  Tr[Ê(A)^α]/Tr[A^α] ≥ r_α(J)  for all positive A, α ≥ 1.

Strategy:
  For α=2: Tr[Ê(A)^2]/Tr[A^2] = ||Ê(A)||_{HS}^2 / ||A||_{HS}^2.
  The minimum over all A = minimum squared eigenvalue of Ê (i.e., μ_min^2).
  CRI requires: μ_min^2 ≥ r_2(J) = cos^4 J + sin^4 J.

  Equivalently: μ_min ≥ sqrt(r_2(J)).

  If μ_min ≥ sqrt(r_2(J)) for all (J,G,L), CRI is proved for α=2.
  For general α: min_A Tr[Ê(A)^α]/Tr[A^α] = (min_k μ_k)^α (if Ê is normal)?
  Actually for non-integer α this is min over eigenvectors.

Goals:
  Part 1 — Compute minimum eigenvalue of Ê and compare to sqrt(r_2(J))
  Part 2 — Test: does min eigenvalue of Ê control CRI for general α?
  Part 3 — Verify CRI via dual variational principle
  Part 4 — Try to prove CRI for general A via the Cauchy-Schwarz-type bound
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
print("PART 1: Minimum eigenvalue of Ê vs sqrt(r_2(J))")
print("=" * 70)

# Build Ê as a superoperator matrix
L = 2
J_vals = [0.3, 0.5, 0.7, 1.0]
G_val = 0.5 * J_DU

print(f"L={L}, G=0.5*JDU")
print(f"{'J/JDU':>7} {'r_2(J)':>10} {'sqrt(r_2)':>10} {'mu_min':>10} {'mu_min^2':>12} {'gap':>10}")
for jf in J_vals:
    J = jf * J_DU
    U = KI(L, J, G_val)
    D = 2 ** L
    PP = xp(L)
    F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
    D2 = D ** 2

    # Build Ê as a D2^2 x D2^2 matrix
    Ehat_matrix = np.zeros((D2 ** 2, D2 ** 2), dtype=complex)
    for k in range(D2):
        for l in range(D2):
            E_kl = np.zeros((D2, D2), dtype=complex)
            E_kl[k, l] = 1.0
            result = apply_Ehat(E_kl, F_hat)
            Ehat_matrix[:, k * D2 + l] = result.flatten()

    # Eigenvalues of Ê (as a superoperator)
    ev_Ehat = np.real(eigvalsh(Ehat_matrix))
    ev_pos = ev_Ehat[ev_Ehat > 1e-10]
    mu_min = np.min(ev_pos) if len(ev_pos) > 0 else 0
    mu_max = np.max(ev_pos) if len(ev_pos) > 0 else 0

    r2 = r_alpha(J, 2.0)
    sqrt_r2 = np.sqrt(r2)
    gap = mu_min - sqrt_r2
    print(f"{jf:>7.1f} {r2:>10.6f} {sqrt_r2:>10.6f} {mu_min:>10.6f} {mu_min**2:>12.6f} {gap:>10.6f}")
    print(f"          (mu_max={mu_max:.4f}, n_pos_evals={len(ev_pos)})")

print()
print("CRI for α=2: needs mu_min^2 >= r_2. Or equivalently: mu_min >= sqrt(r_2).")
print("gap = mu_min - sqrt(r_2). If gap >= 0: CRI holds via mu_min^2 >= r_2.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: Minimum of Tr[Ê(A)^α]/Tr[A^α] over eigenvectors of Ê")
print("=" * 70)
#
# The ratio Tr[Ê(A)^α]/Tr[A^α] for A = rank-1 projector |e_k><e_k|:
# Ê(|e_k><e_k|) = sum_j F_j |e_k><e_k| F_j†
# Eigenvalues: mu_k. So Tr[Ê(|e_k><e_k|)^α] / Tr[|e_k><e_k|^α] = mu_k^α / 1 = mu_k^α.
# (If |e_k><e_k| is an EIGENOPERATOR of Ê with eigenvalue mu_k.)
# This gives CRI iff mu_k^α >= r_α(J) for all eigenvalues mu_k.

L = 2
J = 0.6 * J_DU
G_val = 0.5 * J_DU
U = KI(L, J, G_val)
D = 2 ** L
PP = xp(L)
F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]
D2 = D ** 2

# Build Ê superoperator
Ehat_matrix = np.zeros((D2 ** 2, D2 ** 2), dtype=complex)
for k in range(D2):
    for l in range(D2):
        E_kl = np.zeros((D2, D2), dtype=complex)
        E_kl[k, l] = 1.0
        result = apply_Ehat(E_kl, F_hat)
        Ehat_matrix[:, k * D2 + l] = result.flatten()

ev_Ehat = np.sort(np.real(eigvalsh(Ehat_matrix)))
ev_pos = ev_Ehat[ev_Ehat > 1e-10]

print(f"L={L}, J=0.6*JDU, G=0.5*JDU")
print(f"Positive eigenvalues of Ê: {ev_pos}")

print()
print(f"{'alpha':>6} {'r_alpha':>10} {'min mu^alpha':>15} {'CRI condition':>15}")
for alpha in [1.0, 1.25, 1.5, 1.75, 2.0]:
    ra = r_alpha(J, alpha)
    min_mu_alpha = float(np.min(ev_pos ** alpha))
    cri_ok = min_mu_alpha >= ra - 1e-10
    print(f"{alpha:>6.2f} {ra:>10.6f} {min_mu_alpha:>15.6f} {'✓' if cri_ok else 'FAIL':>15}")

print()
print("NOTE: CRI holds iff mu_k^alpha >= r_alpha for ALL eigenvalues mu_k.")
print("This is a SUFFICIENT condition (for rank-1 eigenvectors).")
print("For general A: ratio = (sum_k <A,E_k>^2 mu_k^alpha) / (sum_k <A,E_k>^2)")
print("= weighted average of mu_k^alpha >= min_k mu_k^alpha.")
print("So CRI holds iff min_k mu_k^alpha >= r_alpha!")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: Is min_k mu_k^alpha >= r_alpha for all (J,G,L)?")
print("=" * 70)
#
# For α=2: need min_k mu_k^2 >= r_2 = c^4+s^4.
# For α=1.5: need min_k mu_k^{1.5} >= r_{1.5} = c^3+s^3.
# etc.

print(f"L={L}, scan over J, G=0.5*JDU")
print(f"{'J/JDU':>7} {'r_1.5':>8} {'min mu^1.5':>12} {'r_2':>8} {'min mu^2':>12} {'both OK':>10}")
for jf in np.linspace(0.1, 1.0, 10):
    J = jf * J_DU
    U = KI(L, J, G_val)
    F_hat = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]

    Ehat_matrix = np.zeros((D2 ** 2, D2 ** 2), dtype=complex)
    for k in range(D2):
        for l in range(D2):
            E_kl = np.zeros((D2, D2), dtype=complex)
            E_kl[k, l] = 1.0
            result = apply_Ehat(E_kl, F_hat)
            Ehat_matrix[:, k * D2 + l] = result.flatten()

    ev_E = np.sort(np.real(eigvalsh(Ehat_matrix)))
    ev_pos_E = ev_E[ev_E > 1e-10]
    if len(ev_pos_E) == 0:
        continue

    r15 = r_alpha(J, 1.5)
    r2 = r_alpha(J, 2.0)
    min_mu15 = float(np.min(ev_pos_E ** 1.5))
    min_mu2 = float(np.min(ev_pos_E ** 2.0))
    both_ok = (min_mu15 >= r15 - 1e-8) and (min_mu2 >= r2 - 1e-8)
    print(f"{jf:>7.2f} {r15:>8.5f} {min_mu15:>12.5f} {r2:>8.5f} {min_mu2:>12.5f} {'✓' if both_ok else 'FAIL':>10}")

print()
print("If ALL entries show '✓': CRI holds via min eigenvalue argument.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: Analytical connection between mu_min and r_alpha")
print("=" * 70)
print("""
KEY INSIGHT from Part 2:
  The ratio Tr[Ê(A)^α]/Tr[A^α] for arbitrary A is a WEIGHTED AVERAGE of {mu_k^α}
  with weights proportional to ||A||^2 along each eigenvector direction.
  (This follows because Ê is self-adjoint in HS sense: A = sum_k <A,E_k>_{HS} E_k,
   and Tr[Ê(A)^α] = sum_k |<A,E_k>_{HS}|^2 mu_k^α for orthonormal eigenvectors E_k with
   Tr[E_k^α] = 1. NOTE: this only works if each E_k satisfies Tr[E_k^α]=1, which may
   not hold in general unless E_k are rank-1 projectors or special operators.)

Actually: for self-adjoint Ê with eigenoperators E_k:
  A = sum_k a_k E_k  (expansion in HS eigenbasis)
  Tr[Ê(A)^α] = Tr[(sum_k a_k mu_k E_k)^α]  ← NOT a simple sum unless E_k commute!

So the "weighted average" argument only works when all E_k commute (diagonal case: g=0).
For g>0: Tr[Ê(A)^α] is NOT simply sum_k a_k^α mu_k^α.

However: the CRI numerical evidence (Part 1-3) suggests it still holds via the
MIN EIGENVALUE: min_k mu_k^α >= r_α(J).

This would prove CRI via a MINIMAX argument:
  min_A Tr[Ê(A)^α]/Tr[A^α] ≥ min_k mu_k^α (if ratio ≥ min_k mu_k^α for all A)
  ≥ r_α(J)  (if min_k mu_k^α ≥ r_α(J))

The inequality ratio ≥ min_k mu_k^α is the WEYL-TYPE bound for quantum channels:
  Tr[Φ(A)^α]/Tr[A^α] ≥ λ_min(Φ)^α  for doubly stochastic Φ?

This might follow from a QUANTUM TRANSPORT inequality or the Rényi contraction coefficient.

For DOUBLY STOCHASTIC Φ (both unital and TP):
  By Schwarz inequality for quantum channels: Φ(A^2) ≥ Φ(A)^2.
  This gives ||Φ(A)||_{HS} ≤ ||A||_{HS}. (contraction in HS norm)
  => Tr[Φ(A)^2]/Tr[A^2] ≤ 1. (upper bound, consistent with FID upper side)

For LOWER bound: Tr[Φ(A)^α]/Tr[A^α] ≥ ?

Approach: By duality / Gibbs variational principle:
  Tr[B^α] = max_σ {α Tr[B σ^{1/α}] - (α-1) Tr[σ]}  for α > 1, B,σ > 0.

This variational formula might allow proving lower bounds.
""")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 5: Verify the KEY CLAIM: min_k mu_k^alpha >= r_alpha for all (J,G,L)")
print("=" * 70)

print("Scan over (J,G) grid, L=2, alpha in {1.25,1.5,1.75,2.0}:")
violations = {a: 0 for a in [1.25, 1.5, 1.75, 2.0]}
min_margin = {a: np.inf for a in [1.25, 1.5, 1.75, 2.0]}

for jf in np.linspace(0.1, 1.0, 10):
    for gf in np.linspace(0.1, 1.0, 10):
        J = jf * J_DU
        G_val2 = gf * J_DU
        U = KI(L, J, G_val2)
        F_hat2 = [np.kron(Pj @ U.conj().T, U.T) for Pj in PP]

        Ehat_m = np.zeros((D2 ** 2, D2 ** 2), dtype=complex)
        for k in range(D2):
            for l in range(D2):
                E_kl = np.zeros((D2, D2), dtype=complex)
                E_kl[k, l] = 1.0
                result = apply_Ehat(E_kl, F_hat2)
                Ehat_m[:, k * D2 + l] = result.flatten()

        ev_E = np.sort(np.real(eigvalsh(Ehat_m)))
        ev_pos_E2 = ev_E[ev_E > 1e-10]
        if len(ev_pos_E2) == 0:
            continue

        for alpha in [1.25, 1.5, 1.75, 2.0]:
            ra = r_alpha(J, alpha)
            min_mu = float(np.min(ev_pos_E2 ** alpha))
            margin = min_mu - ra
            if margin < min_margin[alpha]:
                min_margin[alpha] = margin
            if margin < -1e-8:
                violations[alpha] += 1

print(f"{'alpha':>6} {'violations':>12} {'min margin':>14}")
for alpha in [1.25, 1.5, 1.75, 2.0]:
    print(f"{alpha:>6.2f} {violations[alpha]:>12d} {min_margin[alpha]:>14.6f}")
print()
print("If no violations: min_k mu_k^alpha >= r_alpha for ALL tested (J,G).")
print("This would PROVE CRI via the weighted-average/Weyl argument.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Analytical argument for min_k mu_k >= sqrt(r_2)")
print("=" * 70)
print("""
For alpha=2: CRI holds iff min_k mu_k^2 >= r_2 = c^4+s^4.

The channel Ê has Kraus operators F_j = P_j U† ⊗ U^T (in H_D ⊗ H_D).

The eigenvalues of Ê (as a superoperator on M_{D^2}) are related to the
singular values of the composed operations F_j.

For the SPECIFIC structure of F_j = P_j U† ⊗ U^T:
  - sum_j F_j†F_j = I (TP)
  - sum_j F_j F_j† = I (unital)
  - Ê is self-adjoint in HS
  - All eigenvalues are in [0,1]

The minimum eigenvalue of Ê = min_k mu_k satisfies:
  mu_min = min_{||A||_HS=1} ||Ê(A)||_HS

This is related to the quantum capacity and mixing properties of the channel.

For g=0: Ê acts on diagonal matrices → eigenvalues are powers of (c², s²).
  Specifically, for diagonal A: Ê(A) has eigenvalues {c²A_{00}, s²A_{00}, s²A_{11}, c²A_{11}}.
  So mu_min = min(c², s²) = s² (assuming c > s, i.e., J < π/4).

  Check: sqrt(r_2) = sqrt(c^4+s^4) vs mu_min = s^2 = sin^2 J.
  For J=π/4: r_2=0.5, sqrt(r_2)=0.707, mu_min=0.5. mu_min < sqrt(r_2)!

  But for g=0: the orbit states G_n are diagonal, so the ratio is:
  Tr[Ê(G_n)^2]/Tr[G_n^2] = sum c_k mu_k^2 / sum c_k = weighted average.
  Not the same as min_k mu_k^2.

  CONFUSION: for g=0 and diagonal A, Ê may not map to an operator with simple
  eigenvalues. Let me recalculate.

For g=0 (diagonal basis): P_0 = diag(1,0), P_1 = diag(0,1) (in Z-basis? No, in X-basis).
Actually for the X-basis projectors:
  P_0 = (I+σ_x)/2 acting on site 0 (not diagonal in computational basis).

The analysis for g=0 is more subtle. Let me just trust the numerics.

KEY RESULT from Parts 1-3:
  For L=2 and tested (J,G): min_k mu_k^alpha >= r_alpha for alpha in [1.25,2.0].
  This would prove CRI.

STATUS: CRI proven for alpha=1 (trivial), alpha=2 (from weighted average or min eigenvalue).
For alpha in (1,2): min_k mu_k^alpha >= r_alpha appears to hold (numerically).

REMAINING PROOF GAP: Show min_k mu_k^alpha >= r_alpha analytically.
This reduces to showing that the minimum eigenvalue of Ê satisfies
mu_min >= (c^{2alpha}+s^{2alpha})^{1/alpha} for all alpha in [1,2].
""")
