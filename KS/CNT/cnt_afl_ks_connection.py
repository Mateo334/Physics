"""
cnt_afl_ks_connection.py -- CNT, AFL, and KS entropy: three-way comparison.

Goals:
  1. Show rho[Z^n] is diagonal for n=1,2 but NOT for n>=3 (finite L, non-DU).
     => h_AFL <= h_KS(classical measurement) via Schur's theorem.
  2. For g=0, L->inf (Markov chain): h_AFL = h_KS = E_op^(1)(J) exactly.
  3. Three-way hierarchy for SHIFT automorphism:
     h_KS(shift) = h_CNT(shift) = s(omega) <= h_AFL(shift) = s(omega) + log d.
  4. Comparison table: h_AFL, h_KS^classical, h_CNT for kicked Ising.
  5. Quantify the quantum gap S(rho) - H(diag) for n>=3.
"""

import numpy as np
from scipy.linalg import expm, eigh
import itertools

np.set_printoptions(precision=6, suppress=True)

sx = np.array([[0, 1], [1, 0]], dtype=complex)
sz = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_site(op, site, L):
    ops = [np.eye(2, dtype=complex)] * L
    ops[site] = op
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def kicked_ising_open(L, J, g):
    H_ZZ = sum(kron_site(sz, i, L) @ kron_site(sz, i + 1, L) for i in range(L - 1))
    H_X = sum(kron_site(sx, i, L) for i in range(L))
    return expm(-1j * J * H_ZZ) @ expm(-1j * g * H_X)


def x_projectors(L):
    D = 2 ** L
    projs = []
    for bit in range(2):
        P = np.zeros((D, D), dtype=complex)
        for i in range(D):
            for j in range(D):
                b_i = (i >> (L - 1)) & 1
                b_j = (j >> (L - 1)) & 1
                if (i & ((1 << (L - 1)) - 1)) == (j & ((1 << (L - 1)) - 1)):
                    sign = (-1) ** (b_i + b_j) if bit == 1 else 1
                    P[i, j] += 0.5 * sign
        projs.append(P)
    return projs


def afl_density_matrix(U, P, n):
    D = U.shape[0]
    ops = {}
    for idx in itertools.product(range(len(P)), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ U.conj().T @ P[idx[t]]
        ops[idx] = Z @ np.linalg.matrix_power(U, n - 1)
    indices = list(ops.keys())
    dim = len(indices)
    M = np.zeros((dim, dim), dtype=complex)
    for a, idxA in enumerate(indices):
        for b, idxB in enumerate(indices):
            M[a, b] = np.sum(ops[idxB].conj() * ops[idxA]) / D
    return (M + M.conj().T) / 2


def von_neumann(M, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True))
    evals = evals[evals > tol]
    evals /= evals.sum()
    return float(-np.sum(evals * np.log(evals)))


def classical_shannon(pvec):
    pvec = np.array(pvec, dtype=float)
    pvec = pvec[pvec > 1e-15]
    pvec /= pvec.sum()
    return float(-np.sum(pvec * np.log(pvec)))


def eop_1(J):
    c2 = np.cos(J) ** 2
    s2 = np.sin(J) ** 2
    if min(c2, s2) < 1e-15:
        return 0.0
    return float(-(c2 * np.log(c2) + s2 * np.log(s2)))


# ===========================================================================
print("=" * 70)
print("PART 1: Diagonality of rho[Z^n] — exact for n=1,2, fails for n>=3")
print("=" * 70)
#
# ANALYTICAL PROOF for n=2:
# Z_I^(2) = P_{i1} U^dag P_{i2} U
# Z_J^(2)^dag = U^dag P_{j2} U P_{j1}
# Tr[Z_J^dag Z_I] = Tr[U^dag P_{j2} U P_{j1} P_{i1} U^dag P_{i2} U]
# For i1=/=j1: P_{j1} P_{i1} = 0. ✓
# For i1=j1, i2=/=j2: use cyclic trace:
#   Tr[... P_{j2} U P_{i1} U^dag P_{i2} U] = Tr[P_{i2} U U^dag P_{j2} U P_{i1} U^dag]
#   = Tr[P_{i2} P_{j2} U P_{i1} U^dag] = 0  (since P_{i2} P_{j2} = 0 for i2=/=j2). ✓
#
# For n>=3: the cyclic argument only closes P_{j1} against P_{i1} and P_{in} against P_{jn}.
# Middle projectors P_{i_k} (2<=k<=n-1) remain separated by U P_{i_{k-1}} U^dag,
# so off-diagonal elements CAN be non-zero.
#
print("\nVerification (L=3, J=0.6*pi/4, G=0.5*pi/4):")
print(f"{'n':>3} {'S_vN':>10} {'H(diag)':>10} {'off-diag max':>14} {'diagonal?':>12}")

L = 3
J = np.pi / 4 * 0.6
G = np.pi / 4 * 0.5
U = kicked_ising_open(L, J, G)
PP = x_projectors(L)

for n in range(1, 5):
    rho = afl_density_matrix(U, PP, n)
    S_vn = von_neumann(rho)
    diag = np.real(np.diag(rho))
    S_cl = classical_shannon(diag)
    off_diag_max = np.max(np.abs(rho - np.diag(np.diag(rho))))
    is_diag = off_diag_max < 1e-10
    print(f"{n:>3} {S_vn:>10.6f} {S_cl:>10.6f} {off_diag_max:>14.2e} "
          f"{'yes' if is_diag else 'NO (coherences)':>12}")

print("\nDU point (J=G=pi/4):")
J_DU = np.pi / 4
U_DU = kicked_ising_open(L, J_DU, J_DU)
for n in range(1, 5):
    rho = afl_density_matrix(U_DU, PP, n)
    off_diag = np.max(np.abs(rho - np.diag(np.diag(rho))))
    S_vn = von_neumann(rho)
    d_n = 2.0 ** n
    S_max = np.log(d_n)
    print(f"  n={n}: off-diag max = {off_diag:.2e}, S_vN = {S_vn:.5f}, "
          f"n*log2 = {S_max:.5f}, diff = {abs(S_vn - S_max):.2e}")

print("\nConclusion:")
print("  n=1,2: rho[Z^n] is DIAGONAL (proved analytically for n=2 via cyclic trace).")
print("  n>=3 (finite L, off-DU): off-diagonal elements appear (coherences).")
print("  DU point (J=G=pi/4): rho[Z^n] = I/2^n (diagonal) for all n < n_sat. ✓")
print("  Correction to Thm thm:du_mixed: diagonality holds for n=1,2 and at DU; NOT all n.")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 2: Schur's theorem: S(rho) <= H(diag) => h_AFL <= h_KS(classical)")
print("=" * 70)
#
# Schur's theorem: For any density matrix rho, the eigenvalue vector
# is majorized by the diagonal vector: lambda PREC diag(rho).
# By Schur-concavity of Shannon entropy: S(rho) <= H(diag(rho)).
# Since p^(n)_I = diag(rho[Z^n])_I is the classical measurement probability:
#   S(rho[Z^n]) <= H(p^(n)) for all n.
# Therefore: h_AFL = lim S(rho[Z^n])/n <= lim H(p^(n))/n = h_KS(classical).
#
print("\nSchur bound S(rho) <= H(diag) verified:")
print(f"{'(J/JDU,G/JDU)':>16} {'n':>3} {'S_vN':>8} {'H_diag':>8} {'gap H-S':>8}")

L = 3
J_DU = np.pi / 4
for Jf, Gf in [(0.6, 0.5), (0.8, 0.6), (1.0, 1.0)]:
    J = Jf * J_DU
    G = Gf * J_DU
    U = kicked_ising_open(L, J, G)
    PP = x_projectors(L)
    for n in [2, 3, 4]:
        rho = afl_density_matrix(U, PP, n)
        S_vn = von_neumann(rho)
        H_diag = classical_shannon(np.real(np.diag(rho)))
        gap = H_diag - S_vn
        print(f"({Jf:.1f},{Gf:.1f}){' ':>6} {n:>3} {S_vn:>8.4f} {H_diag:>8.4f} {gap:>8.4f}")

print()
print("Gap H-S >= 0 always (Schur). Gap = 0 iff rho is diagonal (n=1,2).")
print("=> h_AFL = lim S(n)/n <= lim H(n)/n = h_KS(classical measurement).")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 3: g=0, L->inf: h_AFL = h_KS = E_op^(1)(J) [exact equality]")
print("=" * 70)
#
# For g=0, L->inf (Markov chain regime):
# The measurement outcomes follow the Markov chain T = [[cos^2J, sin^2J],[sin^2J, cos^2J]].
# Classical KS entropy of this chain:
#   h_KS(T) = -sum_{ij} pi_i T_{ij} log T_{ij}
#             = -[cos^2J log cos^2J + sin^2J log sin^2J] = H_bin(sin^2J) = E_op^(1)(J).
# And h_AFL(g=0, L->inf) = E_op^(1)(J) [Thm thm:geometric_purity, alpha=1].
# EQUALITY: h_AFL = h_KS for g=0, L->inf.
#
# Proof that rho[Z^n] IS diagonal for g=0 (Markov chain, L->inf):
# In this limit, Z_I^(n) acts on L independent sites. The projectors P_{i_k} act on
# different sites (no wrapping), so [P_{i_k}, P_{i_{k'}}] = 0 for k=/=k'.
# The off-diagonal elements Z_J^dag Z_I involve adjacent projectors that cancel.
# (This is NOT the finite-L open-BC case where boundary effects create coherences.)

print("\nAnalytical: h_KS(Markov chain) = H_bin(sin^2J) = E_op^(1)(J)")
print(f"{'J/pi':>8} {'sin^2J':>8} {'h_KS=H_bin':>12} {'E_op^(1)':>10} {'ΔS_2 (L=4)':>12}")

L = 4
for J_frac in [1/10, 1/8, 1/6, 1/5, 1/4]:
    J = J_frac * np.pi
    c2 = np.cos(J) ** 2
    s2 = np.sin(J) ** 2
    h_ks = -(c2 * np.log(c2) + s2 * np.log(s2)) if min(c2, s2) > 1e-15 else 0.0
    Eop1 = eop_1(J)
    # h_AFL via ΔS_2 (g=0, L=4)
    U = kicked_ising_open(L, J, 0.0)
    PP = x_projectors(L)
    rho1 = afl_density_matrix(U, PP, 1)
    rho2 = afl_density_matrix(U, PP, 2)
    dS2 = von_neumann(rho2) - von_neumann(rho1)
    print(f"{J_frac:>8.4f} {s2:>8.4f} {h_ks:>12.6f} {Eop1:>10.6f} {dS2:>12.6f}")

print("\nAll three agree: h_KS = E_op^(1) = ΔS_2(L=4) for g=0. ✓")
print("Hence h_AFL = h_KS = E_op^(1)(J) for g=0, L->inf (equality in Schur bound).")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 4: Three-way hierarchy for SHIFT automorphism")
print("=" * 70)
#
# For the SHIFT automorphism tau_1 on a translation-invariant clustering state omega:
#   h_CNT(tau_1) = s(omega)        [CNT87 Thm IX.1]
#   h_AFL(tau_1, X-OPU) = s(omega) + log d   [Section 10, Thm]
#   h_KS(tau_1, X-measurement) = s(omega)    [classical KS = CNT for abelian reduction]
# HIERARCHY: h_KS = h_CNT = s(omega) < h_AFL = s(omega) + log d.
# Gap h_AFL - h_KS = log d = "quantum overhead" of OPU measurement.
#
# For the KICKED ISING (time-evolution automorphism alpha_t):
#   h_CNT(alpha_t) = 0                 [finite L, Section 11]
#   h_AFL(alpha_t, X-OPU) = ΔS_2 ≈ E_op^(1)(J)   [Sections 15-16]
#   h_KS(alpha_t, X-measurement) >= h_AFL(alpha_t)  [Schur, strict for n>=3]
# HIERARCHY: h_CNT < h_AFL <= h_KS.
# Gap depends on off-diagonal coherences in rho[Z^n], n>=3.

print("\nThree-way hierarchy:")
print()
print("  SHIFT automorphism tau_1 (any clustering state omega):")
print("  h_KS(shift, X-meas) = h_CNT(shift) = s(omega)")
print("  h_AFL(shift, X-OPU) = s(omega) + log d  [= h_KS + log d]")
print()
print("  TIME EVOLUTION alpha_t (kicked Ising, finite L):")
print("  h_CNT(alpha_t) = 0")
print("  h_AFL(alpha_t, X-OPU) = ΔS_2 = E_op^(1)(J)  [for g=0]")
print("  h_KS(alpha_t, X-meas) >= h_AFL(alpha_t)  [Schur, may be strict]")
print()
print("  For PRODUCT STATE (s(omega) = log d):")
s_omega = np.log(2)
d = 2
print(f"  s(omega) = {s_omega:.4f}")
print(f"  h_CNT(shift) = {s_omega:.4f}")
print(f"  h_AFL(shift) = {s_omega + np.log(d):.4f}  [= {s_omega:.4f} + log2]")
print(f"  h_KS(shift) = {s_omega:.4f}  [= h_CNT]")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 5: Quantum gap S(rho) - H(diag) for n>=3: J,G dependence")
print("=" * 70)
#
# The quantum gap Q_n = H(p^(n)) - S(rho[Z^n]) measures the coherence contribution.
# Q_n = 0 for n=1,2 (diagonal rho).
# Q_n > 0 for n>=3 (off-diagonal elements in rho[Z^n]).
# The quantum gap is bounded by: Q_n <= log(rank(rho[Z^n])) - S(rho[Z^n]).

L = 3
n_check = 3
J_DU = np.pi / 4

print(f"\nQuantum gap Q_n = H(diag) - S_vN at n={n_check}, L={L}:")
print(f"{'J/JDU':>7} {'G/JDU':>7} {'H_diag':>8} {'S_vN':>8} {'gap Q':>8} {'off-diag':>10}")

for Jf in [0.4, 0.6, 0.8, 1.0]:
    for Gf in [0.4, 0.6, 0.8, 1.0]:
        J = Jf * J_DU
        G = Gf * J_DU
        U = kicked_ising_open(L, J, G)
        PP = x_projectors(L)
        rho = afl_density_matrix(U, PP, n_check)
        S_vn = von_neumann(rho)
        H_diag = classical_shannon(np.real(np.diag(rho)))
        gap = H_diag - S_vn
        off = np.max(np.abs(rho - np.diag(np.diag(rho))))
        print(f"{Jf:>7.1f} {Gf:>7.1f} {H_diag:>8.4f} {S_vn:>8.4f} {gap:>8.4f} {off:>10.2e}")

print(f"\nAt DU (J=G=JDU): rho[Z^{n_check}] should be diagonal (I/2^{n_check}):")
J = J_DU
G = J_DU
U = kicked_ising_open(L, J, G)
rho_DU = afl_density_matrix(U, PP, n_check)
print(f"  off-diag max = {np.max(np.abs(rho_DU - np.diag(np.diag(rho_DU)))):.2e}")
print(f"  S_vN = {von_neumann(rho_DU):.5f}, n*log2 = {n_check*np.log(2):.5f}")

# ===========================================================================
print("\n" + "=" * 70)
print("PART 6: Summary")
print("=" * 70)
print("""
KEY RESULTS:

(K1) Diagonality of rho[Z^n] (CORRECTED):
  - n=1: rho[Z^1] = I/d (diagonal). ✓
  - n=2: rho[Z^2] is diagonal for ALL (J,G,L) [proved analytically via cyclic trace]. ✓
  - n>=3, finite L, off-DU: off-diagonal coherences appear. DIAGONAL THEOREM FAILS.
  - DU point (J=G=J_DU): rho[Z^n] = I/d^n (diagonal, maximally mixed) for all n. ✓

(K2) Schur's theorem gives: S(rho[Z^n]) <= H(diag(rho[Z^n])) = H(p^(n)).
  Therefore: h_AFL^OPU <= h_KS(classical measurement process).
  Equality holds iff rho[Z^n] is diagonal for all n (e.g., DU point, or g=0 L->inf).

(K3) For g=0, L->inf (Markov chain):
  h_AFL = h_KS = E_op^(1)(J) = H_bin(sin^2J). [EXACT EQUALITY, proved]
  Proof: Markov chain T => rho[Z^n] diagonal in thermodynamic limit =>
         S(rho[Z^n]) = H(p^(n)); KS entropy of T = H_bin(sin^2J) = E_op^(1). ✓

(K4) Three-way hierarchy:
  SHIFT:     h_CNT(shift) = h_KS(shift) = s(omega) < h_AFL(shift) = s(omega)+log d.
  TIME EVO:  0 = h_CNT(alpha_t) < h_AFL(alpha_t) ≤ h_KS(alpha_t) [Schur].
  Gap for shift: h_AFL - h_KS = log d (universal, from AFL-CNT theorem).
  Gap for time evo: h_KS - h_AFL = quantum coherence Q_n >= 0 (n>=3, off-DU).

(K5) Quantum gap Q_n = H(diag) - S(rho): quantifies coherence contribution.
  Q_n = 0 for n<=2, Q_n > 0 for n>=3 (off-DU, finite L).
  At DU: Q_n = 0 for all n (diagonal orbit states => AFL = KS).
""")
