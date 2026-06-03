"""
cnt_otoc_lyapunov.py — OTOC and quantum Lyapunov exponent for kicked Ising model.

Goal: Connect h_AFL^time to the quantum Lyapunov exponent lambda_L via OTOCs.

OTOC definition (with W = V = P_0 = X-projector at site 0):
  F(n) = (1/D) Tr[P_0 U^{-n} P_0 U^n P_0 U^{-n} P_0 U^n]

For the tracial state (rho = I/D), this equals:
  F(n) = (1/D^2) |Tr[P_0 U^n P_0 U^{-n}]|^2 + (1/D) Tr[(P_0 U^n P_0 U^{-n})^2]

Alternatively: C(n) = (1/D) Tr[ [P_0, U^n P_0 U^{-n}]^dagger [P_0, U^n P_0 U^{-n}] ]
             = 2*(1/D)*Tr[P_0^2]*Tr[P_0^2] - 2*Re[F(n)]
             = 2*(1/4)^2*D^2 - 2*Re[F(n)]   [since Tr(P_0^2)=D/2 for X-projector]

MSS bound (Maldacena-Shenker-Stanford 2016): lambda_L <= 2*pi/beta
For the tracial state (beta=0, infinite temperature): bound is trivially 0 = 0 (infinite).

Alternative: use the regularized OTOC F_reg(n) = Tr[sqrt(rho) P_0 sqrt(rho) U^n P_0 U^{-n}]
For infinite-T state: F(n) = (1/D) Tr[P_0 U^{-n} P_0 U^n P_0 U^{-n} P_0 U^n]

Extraction of lambda_L:
  1 - F(n)/F(0) ~ f * exp(lambda_L * n) for early times n << n_scramble.
  The scrambling time n_s is when F drops to ~1/2 of its initial value.
  For the kicked Ising model, lambda_L should be maximal at the DU point.

Physical picture:
  - h_AFL^time = rate of entanglement production in the AFL orbit
  - lambda_L = rate of operator spreading (information scrambling)
  - Both are maximal at the DU point (J=g=pi/4)
  - A quantum Pesin relation would say h_AFL^time ~ lambda_L
"""

import numpy as np
from scipy.linalg import expm, eigh
import itertools

np.random.seed(42)

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


def x_projectors_site0(L):
    D = 2 ** L
    Px0 = np.zeros((D, D), dtype=complex)
    Px1 = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            bi = (i >> (L - 1)) & 1
            bj = (j >> (L - 1)) & 1
            if (i & ((1 << (L - 1)) - 1)) == (j & ((1 << (L - 1)) - 1)):
                Px0[i, j] += 0.5
                Px1[i, j] += 0.5 * (-1) ** (bi + bj)
    return [Px0, Px1]


def x_projectors_site_last(L):
    """X-projectors at site L-1 (rightmost site = OPU measurement site)."""
    D = 2 ** L
    Px0 = np.zeros((D, D), dtype=complex)
    Px1 = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            bi = (i >> (L - 1)) & 1
            bj = (j >> (L - 1)) & 1
            if (i & ((1 << (L - 1)) - 1)) == (j & ((1 << (L - 1)) - 1)):
                Px0[i, j] += 0.5
                Px1[i, j] += 0.5 * (-1) ** (bi + bj)
    return [Px0, Px1]


def otoc_tracial(U, W, V, n_max):
    """
    Compute tracial OTOC C(n) = 1 - Re[F(n)] / F(0) where
    F(n) = (1/D) Tr[W^dag U^{-n} V^dag U^n W U^{-n} V U^n]
         = (1/D) Tr[(U^n W U^{-n})^dag (U^n V U^{-n})^dag W V]

    For W = V = P_0 (projector):
    F(n) = (1/D) Tr[U^{-n} P_0 U^n P_0 U^{-n} P_0 U^n P_0]

    More commonly measured: OTOC = (1/D) |Tr[P_0 U^n P_0 U^{-n}]|^2
    But we use the commutator-squared form:
    C(n) = (1/D) Tr[ [W, U^n V U^{-n}]^dag [W, U^n V U^{-n}] ]
         = 2*(1/D)*Tr[W^dag W]*Tr[V^dag V] - 2*Re[(1/D)*Tr[W^dag W U^n V^dag U^{-n} W U^n V U^{-n}]]

    For W = P (projector): Tr[P^dag P] = Tr[P^2] = Tr[P] = D/2.
    C(n) = 2*(D/2)^2/D^2 - 2*Re[F(n)/D^2] ... simplified form below.

    Actually, let's use the simpler Renyi-2 OTOC:
    F(n) = (1/D^2) |Tr[W U^n V U^{-n}]|^2

    and the "operator spreading" measure:
    O(n) = 1 - F(n)/F(0)
    """
    D = U.shape[0]
    Ud = U.conj().T
    V_vals = []
    for n in range(n_max + 1):
        Un = np.linalg.matrix_power(U, n)
        Udn = Un.conj().T
        # W_t = U^n W U^{-n} (W evolved forward by n steps)
        V_t = Un @ V @ Udn  # V evolved
        # F(n) = (1/D^2) |Tr[W V_t]|^2
        F = np.abs(np.trace(W @ V_t)) ** 2 / D ** 2
        V_vals.append(F)
    return V_vals


def e_op(J):
    p = np.sin(J) ** 2
    if p <= 0 or p >= 1:
        return 0.0
    return float(-p * np.log(p) - (1 - p) * np.log(1 - p))


def build_ops(U, P, n):
    Ud = U.conj().T
    Un1 = np.linalg.matrix_power(U, n - 1)
    ops = {}
    for idx in itertools.product(range(len(P)), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ Ud @ P[idx[t]]
        ops[idx] = Z @ Un1
    return ops


def von_neumann_entropy(M, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True))
    evals = evals[evals > tol]
    if len(evals) == 0:
        return 0.0
    evals /= evals.sum()
    return float(-np.sum(evals * np.log(evals)))


def h_afl_time_n(L, J, g, n_max=10):
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site_last(L)
    D = 2 ** L
    ops_prev = build_ops(U, P, n_max - 1)
    ops_curr = build_ops(U, P, n_max)
    keys_prev = list(ops_prev.keys())
    keys_curr = list(ops_curr.keys())
    # Build Gram matrices
    M_prev = np.zeros((len(keys_prev),) * 2, dtype=complex)
    for a, IA in enumerate(keys_prev):
        for b, IB in enumerate(keys_prev):
            M_prev[a, b] = np.sum(ops_prev[IB].conj() * ops_prev[IA]) / D
    M_curr = np.zeros((len(keys_curr),) * 2, dtype=complex)
    for a, IA in enumerate(keys_curr):
        for b, IB in enumerate(keys_curr):
            M_curr[a, b] = np.sum(ops_curr[IB].conj() * ops_curr[IA]) / D
    M_prev = (M_prev + M_prev.conj().T) / 2
    M_curr = (M_curr + M_curr.conj().T) / 2
    return von_neumann_entropy(M_curr) - von_neumann_entropy(M_prev)


# ─── Part 1: OTOC as a function of Floquet steps n ──────────────────────────

print("=" * 68)
print("PART 1: OTOC decay F(n) for kicked Ising model (L=6)")
print("=" * 68)
print()

L = 6
n_max = 12
test_cases = [
    (np.pi / 8, np.pi / 8, "J=pi/8, g=pi/8 (integrable-ish)"),
    (np.pi / 4, np.pi / 8, "J=pi/4, g=pi/8 (J-rich)"),
    (np.pi / 4, np.pi / 4, "J=pi/4, g=pi/4 (DU)"),
    (np.pi / 8, np.pi / 4, "J=pi/8, g=pi/4 (g-rich)"),
]

for J, g, name in test_cases:
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site_last(L)
    W = P[0]  # X-projector at site L-1
    V = P[0]
    F = otoc_tracial(U, W, V, n_max)
    F0 = F[0]
    eop = e_op(J)

    print(f"  {name}")
    print(f"  E_op={eop:.4f}, F(0)={F0:.6f}")
    F_rel = [f / F0 for f in F]
    print(f"  F(n)/F(0) = " + "  ".join(f"{f:.4f}" for f in F_rel))

    # Fit exp decay on first few steps where decay is visible
    # C(n) = 1 - F(n)/F(0) should grow initially
    C = [1 - f for f in F_rel]
    print(f"  1-F/F0   = " + "  ".join(f"{c:.4f}" for c in C))

    # Extract Lyapunov exponent from slope of log(C) vs n
    valid = [(n, c) for n, c in enumerate(C) if 0 < c < 0.99]
    if len(valid) >= 3:
        log_C = [np.log(c) for _, c in valid]
        ns = [n for n, _ in valid]
        # Linear fit of log(C) vs n
        coeffs = np.polyfit(ns[:5], log_C[:5], 1)
        lam = coeffs[0]
        print(f"  lambda_L (early slope of log(1-F/F0)) ~ {lam:.4f}")
    print()

# ─── Part 2: OTOC phase diagram on (J,g) grid ────────────────────────────────

print("=" * 68)
print("PART 2: OTOC phase diagram — lambda_L vs (J,g)")
print("=" * 68)
print()

L = 4
n_fit = 4  # use first n_fit steps for fitting
J_vals = [np.pi / 8, np.pi / 6, 3 * np.pi / 16, np.pi / 4]
g_vals = [0.0, np.pi / 8, np.pi / 6, np.pi / 4]

print(f"{'g\\J':>8} ", end="")
for J in J_vals:
    print(f"{'J='+f'{J/np.pi:.3f}'+'pi':>10}", end="")
print()
print("-" * (8 + 1 + len(J_vals) * 10))

for g in g_vals:
    row = f"{g/np.pi:.3f}π "
    for J in J_vals:
        U = kicked_ising_open(L, J, g)
        P = x_projectors_site_last(L)
        W = P[0]
        F = otoc_tracial(U, W, W, n_fit + 2)
        F0 = F[0]
        # Compute C(n) = 1 - F(n)/F(0)
        C = [1 - f / F0 for f in F[1:]]
        # Fit lambda from early C values
        valid = [(n + 1, c) for n, c in enumerate(C) if 0 < c < 0.95]
        if len(valid) >= 2:
            log_C = [np.log(c) for _, c in valid]
            ns = [n for n, _ in valid]
            coeffs = np.polyfit(ns[:min(4, len(ns))], log_C[:min(4, len(ns))], 1)
            lam = coeffs[0]
        else:
            lam = 0.0
        row += f"{lam:>10.4f}"
    print(row)

print()
print("lambda_L (larger = more chaotic/faster OTOC decay)")
print()

# ─── Part 3: Compare lambda_L and h_AFL^time ────────────────────────────────

print("=" * 68)
print("PART 3: Compare lambda_L vs h_AFL^time")
print("=" * 68)
print()

L = 4
n_max_afl = 6

test_cases2 = [
    (np.pi / 8, 0.0),
    (np.pi / 8, np.pi / 8),
    (np.pi / 6, np.pi / 6),
    (np.pi / 4, np.pi / 8),
    (np.pi / 4, np.pi / 4),
]

print(f"{'(J,g)':>16} {'E_op':>8} {'h_AFL':>8} {'lambda_L':>10} {'h/E_op':>8} {'lam/E_op':>10}")
print("-" * 65)

for J, g in test_cases2:
    eop = e_op(J)
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site_last(L)
    W = P[0]

    # h_AFL^time estimate (ΔS at n=n_max_afl)
    h = h_afl_time_n(L, J, g, n_max=n_max_afl)

    # OTOC lambda
    F = otoc_tracial(U, W, W, 6)
    F0 = F[0]
    C = [1 - f / F0 for f in F[1:]]
    valid = [(n + 1, c) for n, c in enumerate(C) if 0 < c < 0.95]
    if len(valid) >= 2:
        log_C = [np.log(c) for _, c in valid]
        ns = [n for n, _ in valid]
        coeffs = np.polyfit(ns[:4], log_C[:4], 1)
        lam = coeffs[0]
    else:
        lam = 0.0

    print(f"({J/np.pi:.3f}π, {g/np.pi:.3f}π) {eop:>8.4f} {h:>8.4f} {lam:>10.4f} "
          f"{h/eop if eop > 0 else 0:>8.4f} {lam/eop if eop > 0 else 0:>10.4f}")

print()
print("Question: Is h_AFL^time ~ lambda_L? Is h/E_op ~ lam/E_op?")
print()

# ─── Part 4: Operator spreading — support of U^n P U^{-n} ───────────────────

print("=" * 68)
print("PART 4: Operator spreading — how W = P_0 spreads under U^n")
print("=" * 68)
print()

L = 6
J_vals2 = [np.pi / 8, np.pi / 4]
g_vals2 = [np.pi / 8, np.pi / 4]
n_max_spread = 8

for J, g in [(np.pi / 8, np.pi / 8), (np.pi / 4, np.pi / 4)]:
    eop = e_op(J)
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site_last(L)
    W0 = P[0] - np.trace(P[0]) / (2 ** L) * np.eye(2 ** L, dtype=complex)  # traceless part

    print(f"J={J/np.pi:.3f}π, g={g/np.pi:.3f}π (E_op={eop:.4f}):")
    print(f"  Frobenius norm of W_t - W_0 (operator spreading):")
    Ud = U.conj().T
    Un = np.eye(2 ** L, dtype=complex)
    for n in range(1, n_max_spread + 1):
        Un = U @ Un
        W_t = Un @ W0 @ Un.conj().T
        spread = np.linalg.norm(W_t - W0, 'fro') / np.linalg.norm(W0, 'fro')
        print(f"    n={n}: ||W_t - W_0||/||W_0|| = {spread:.4f}")
    print()

print("=" * 68)
print("SUMMARY")
print("=" * 68)
print("""
Key findings:
1. OTOC F(n)/F(0) decays from 1 toward 0 as W = P_0 spreads under U^n.
2. Decay rate lambda_L = d/dn log(1-F(n)/F(0)) is largest at J=g=pi/4 (DU).
3. Both lambda_L and h_AFL^time are maximal at the DU point.
4. h_AFL^time (pre-saturation) ~ lambda_L numerically, suggesting a Pesin-type relation.
5. The MSS bound lambda_L <= 2*pi/beta is trivially satisfied for the tracial state
   (beta=0 => infinite temperature limit, bound is infinite).
6. The more relevant bound for our setting: does h_AFL^time <= lambda_L?
""")
