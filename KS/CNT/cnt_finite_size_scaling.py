"""
cnt_finite_size_scaling.py — Finite-size scaling of the Pesin gap.

Goal: compute Delta(J,g,L) = E_op(J) - h_AFL^time(J,g,L) for L = 4,5,6,8
to determine the thermodynamic limit behavior.

Key questions:
  (1) Does Delta(J,0,L) -> 0 for all J as L -> inf? (Markovian limit, Prop 16.4)
  (2) Is the g=0 line the ONLY line where the gap closes?
  (3) What is the finite-size correction exponent: Delta ~ C(J,g)/L^nu?
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


def time_afl_density_matrix(U, P, n):
    D = U.shape[0]
    ops = build_ops(U, P, n)
    indices = list(ops.keys())
    dim = len(indices)
    M = np.zeros((dim, dim), dtype=complex)
    for a, idxA in enumerate(indices):
        for b, idxB in enumerate(indices):
            M[a, b] = np.sum(ops[idxB].conj() * ops[idxA]) / D
    return (M + M.conj().T) / 2


def von_neumann_entropy(M, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True))
    evals = evals[evals > tol]
    if len(evals) == 0:
        return 0.0
    evals /= evals.sum()
    return float(-np.sum(evals * np.log(evals)))


def h_afl_time(L, J, g, n_max=8):
    """Estimate h_AFL^time as lim_{n->inf} (S_n - S_{n-1})."""
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site0(L)
    S = [von_neumann_entropy(time_afl_density_matrix(U, P, n)) for n in range(1, n_max + 1)]
    dS = [S[i] - S[i - 1] for i in range(1, len(S))]
    # Use last delta as best estimate (increments are non-increasing)
    return dS[-1], dS


def e_op(J):
    """E_op(J) = H_bin(sin^2 J)."""
    p = np.sin(J) ** 2
    if p <= 0 or p >= 1:
        return 0.0
    return float(-p * np.log(p) - (1 - p) * np.log(1 - p))


# ─── Part 1: Finite-size scaling along g=0 (Markovian line) ─────────────────

print("=" * 72)
print("PART 1: Finite-size scaling along g=0 (Markovian limit)")
print("=" * 72)
print()
print("Prediction (Prop 16.4): h_AFL^time(J,0,L->inf) = E_op(J)")
print("Expected: Delta(J,0,L) -> 0 as L -> inf.")
print()

Ls = [4, 5, 6]
J_vals = [np.pi / 8, np.pi / 6, 3 * np.pi / 16, np.pi / 4]
g = 0.0

print(f"{'J':>8} {'E_op':>8} ", end="")
for L in Ls:
    print(f"{'h(L='+str(L)+')':>10} {'Delta':>8}", end="")
print()
print("-" * (8 + 8 + 3 + len(Ls) * 18))

results_g0 = {}
for J in J_vals:
    eop = e_op(J)
    row = f"{J/np.pi:.4f}π {eop:.5f} "
    row_data = []
    for L in Ls:
        h, dS = h_afl_time(L, J, g, n_max=8)
        delta = eop - h
        row += f" {h:>9.5f}  {delta:>7.4f}"
        row_data.append((L, h, delta))
    print(row)
    results_g0[J] = (eop, row_data)

print()
print("Observation: Delta(J,0,L) decreasing with L (approaching 0)?")
print()

# ─── Part 2: Finite-size scaling at g=J (diagonal line) ─────────────────────

print("=" * 72)
print("PART 2: Finite-size scaling along g=J (self-dual diagonal)")
print("=" * 72)
print()

print(f"{'J=g':>8} {'E_op':>8} ", end="")
for L in Ls:
    print(f"{'h(L='+str(L)+')':>10} {'Delta':>8}", end="")
print()
print("-" * (8 + 8 + 3 + len(Ls) * 18))

results_diag = {}
for J in J_vals:
    eop = e_op(J)
    row = f"{J/np.pi:.4f}π {eop:.5f} "
    row_data = []
    for L in Ls:
        h, dS = h_afl_time(L, J, J, n_max=8)
        delta = eop - h
        row += f" {h:>9.5f}  {delta:>7.4f}"
        row_data.append((L, h, delta))
    print(row)
    results_diag[J] = (eop, row_data)

print()

# ─── Part 3: Finite-size scaling at fixed J=pi/4 varying g ─────────────────

print("=" * 72)
print("PART 3: Finite-size scaling at J=pi/4 varying g")
print("=" * 72)
print()

J = np.pi / 4
eop = e_op(J)
g_vals = [0.0, np.pi / 8, np.pi / 4 - 0.1, np.pi / 4]

print(f"J=pi/4, E_op={eop:.5f}")
print(f"{'g':>8} ", end="")
for L in Ls:
    print(f"{'h(L='+str(L)+')':>10} {'Delta':>8}", end="")
print()
print("-" * (8 + 3 + len(Ls) * 18))

results_J_pi4 = {}
for g in g_vals:
    row = f"{g/np.pi:.4f}π "
    row_data = []
    for L in Ls:
        h, dS = h_afl_time(L, J, g, n_max=8)
        delta = eop - h
        row += f" {h:>9.5f}  {delta:>7.4f}"
        row_data.append((L, h, delta))
    print(row)
    results_J_pi4[g] = row_data

print()

# ─── Part 4: Power-law fit Delta ~ C/L^nu for g=0 line ─────────────────────

print("=" * 72)
print("PART 4: Power-law fit Delta(J,0,L) ~ C(J)/L^nu(J)")
print("=" * 72)
print()

print("Fitting Delta = C * L^(-nu) using L=4,5,6:")
print()
print(f"{'J':>10} {'nu':>8} {'C':>10} {'Delta(L=4)':>12} {'Delta(L=6)':>12}")
print("-" * 55)

for J in [np.pi / 8, np.pi / 6, 3 * np.pi / 16]:
    eop, row_data = results_g0[J]
    Ls_fit = [r[0] for r in row_data]
    deltas_fit = [r[2] for r in row_data]
    # Filter out near-zero deltas (DU point)
    valid = [(l, d) for l, d in zip(Ls_fit, deltas_fit) if d > 1e-6]
    if len(valid) >= 2:
        log_L = np.log([v[0] for v in valid])
        log_D = np.log([v[1] for v in valid])
        nu, log_C = np.polyfit(log_L, log_D, 1)
        C = np.exp(log_C)
        print(f"{J/np.pi:.4f}π  {-nu:>7.3f}  {C:>9.4f}  "
              f"{deltas_fit[0]:>11.4f}  {deltas_fit[-1]:>11.4f}")
    else:
        print(f"{J/np.pi:.4f}π  (insufficient data)")

print()

# ─── Part 5: Entropy increments at large n (convergence check) ──────────────

print("=" * 72)
print("PART 5: Entropy increment convergence ΔS_n vs n for L=4,6")
print("=" * 72)
print()

test_params = [
    (np.pi / 8, 0.0, "J=pi/8, g=0"),
    (np.pi / 8, np.pi / 8, "J=pi/8, g=pi/8"),
    (np.pi / 4, 0.0, "J=pi/4, g=0"),
]

for J, g, name in test_params:
    eop = e_op(J)
    print(f"  {name} (E_op={eop:.4f})")
    for L in [4, 6]:
        _, dS = h_afl_time(L, J, g, n_max=8)
        dS_str = "  ".join(f"{d:.4f}" for d in dS)
        print(f"    L={L}: ΔS = {dS_str}  (h_est={dS[-1]:.4f}, Δ={eop-dS[-1]:.4f})")
    print()

print("=" * 72)
print("SUMMARY")
print("=" * 72)
print("""
Key findings — CORRECTED INTERPRETATION (see cnt_L_independence.py):

1. L-INDEPENDENCE THEOREM (proved in Section 18 of Output.tex):
   rho[Z^n](L) = rho[Z^n](2) for all L >= 2 and all n <= n_sat(L=2).
   This means ΔS_n(J,g,L) is the same for ALL finite L (for short orbits).
   In particular: h_AFL^time(J,g,L) = 0 for ALL finite L (rank argument, Sec 11).
   Therefore: Delta(J,g,L) = E_op(J) for ALL finite L (except DU point).

2. MARKOVIAN EQUALITY at g=0 (L -> inf only):
   lim_{L->inf} lim_{n->inf} ΔS_n(J,0,L) = E_op(J)  [Markov chain, Prop 16.4]
   lim_{n->inf} ΔS_n(J,0,L_finite) = 0               [rank saturation, Sec 11]
   The limits n->inf and L->inf DO NOT COMMUTE for g=0.

3. DU POINT (J=g=pi/4): h = E_op for ALL L and ALL n. No saturation ever.
   This is the ONLY (J,g) where h > 0 for finite L.

4. For g>0 and n large: h_est(L) appears to grow with L (pre-saturation regime
   is longer for larger L), but this is just the saturation happening at larger n.
   The TRUE h_AFL^time = 0 for all finite L with g != (pi/4, pi/4).

5. Pre-saturation rate: ΔS_n(J,g,L) for n << n_sat(L) is L-independent.
   The rate at n=2 (exactly ΔS_2 = E_op) is the finite-L signature of the
   Markovian rate, valid for ALL L.
""")
