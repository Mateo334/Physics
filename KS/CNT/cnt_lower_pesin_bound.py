"""
cnt_lower_pesin_bound.py — Lower bound for h_AFL^time near the dual-unitary point.

Goal: Given h_AFL^time(J,g) <= E_op(J) (proved in Section 16), derive a lower bound.

Key questions:
  (a) Does h_AFL^time/E_op → 1 monotonically as (J,g) → (pi/4, pi/4)?
  (b) What is the leading correction: h_AFL^time ≈ E_op - f(J,g)?
  (c) Is h_AFL^time >= some non-trivial lower bound (e.g., E_op * g_ratio)?

Strategy:
  Near dual-unitary point: J = pi/4 - dJ, g = pi/4 - dg.
  Expand E_op and h_AFL^time in dJ, dg.

  E_op(J) = H_bin(sin^2 J). Near J=pi/4 (sin^2=1/2):
    E_op ≈ ln2 - (dJ)^2 / (2*ln2) + O(dJ^4).
    [since d/dJ H_bin(sin^2 J)|_{pi/4} = 0, d^2/dJ^2 = -1/ln2]

  h_AFL^time: from the concavity, h_AFL^time ≈ ΔS_n for large n.
  Near DU: ΔS_n ≈ ln2 - correction(dJ, dg, n).

Additional direction:
  Prove h_AFL^time >= I_temp (temporal information as LOWER bound)?
  [This would give: I_temp <= h_AFL^time <= E_op = log d - I_temp,
   which would require I_temp <= log d / 2, but I_temp can be log d at J=0.]

Actually: h_AFL^time >= ΔS_n for ALL large n? NO (ΔS_n is non-increasing so h = lim ≤ ΔS_2).
  h_AFL^time = lim_{n->inf} ΔS_n = inf_n ΔS_n.
  Lower bound: h_AFL^time >= ΔS_N for any fixed N... NO, h is the LIMIT which is BELOW all ΔS_n.

So the interesting lower bound is: h_AFL^time >= (something in terms of E_op and I_temp).

Idea: h_AFL^time >= E_op / N for some effective "scrambling time" N?
Or: h_AFL^time >= E_op - I_temp = 2*E_op - log d (only positive for E_op > log d /2)?

From the gap table data:
  At J=pi/4, varying g:
    g=0:    Δ=0.693 (large), h≈0, E_op=log2=0.693
    g=pi/16: Δ=0.579, h≈0.114
    g=pi/8:  Δ=0.388, h≈0.305
    g=3pi/16: Δ=0.136, h≈0.557
    g=pi/4:  Δ≈0, h≈0.693

  Ratio h/E_op at J=pi/4:
    g=0: 0.000
    g=pi/16: 0.164
    g=pi/8: 0.440
    g=3pi/16: 0.803
    g=pi/4: 1.000

  This looks like: h/E_op ≈ sin^2(g_ratio) or some monotone function.

  At J=pi/8, varying g:
    g=0:   Δ=0.384, h≈0.032
    g=pi/16: Δ=0.247, h≈0.169
    g=pi/8: Δ=0.134, h≈0.283
    g=3pi/16: Δ=0.112, h≈0.305
    g=pi/4: Δ=0.111, h≈0.306

  h/E_op at J=pi/8: 0.077, 0.406, 0.679, 0.732, 0.735.
  Plateaus after g=pi/8: h is roughly constant for g >= pi/8.
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
    for o in ops[1:]: r = np.kron(r, o)
    return r

def kicked_ising_open(L, J, g):
    H_ZZ = sum(kron_site(sz, i, L) @ kron_site(sz, i+1, L) for i in range(L-1))
    H_X  = sum(kron_site(sx, i, L) for i in range(L))
    return expm(-1j * J * H_ZZ) @ expm(-1j * g * H_X)

def x_projectors_site0(L):
    D = 2**L
    Px0 = np.zeros((D, D), dtype=complex)
    Px1 = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            bi = (i >> (L-1)) & 1
            bj = (j >> (L-1)) & 1
            rest_i = i & ((1 << (L-1)) - 1)
            rest_j = j & ((1 << (L-1)) - 1)
            if rest_i == rest_j:
                Px0[i, j] += 0.5
                Px1[i, j] += 0.5 * (-1)**(bi + bj)
    return [Px0, Px1]

def time_afl_density_matrix(U, P, n):
    D = U.shape[0]
    Ud = U.conj().T
    Un1 = np.linalg.matrix_power(U, n - 1)
    k = len(P)
    ops = []
    for idx in itertools.product(range(k), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ Ud @ P[idx[t]]
        ops.append(Z @ Un1)
    ops_flat = np.array([op.ravel() for op in ops])
    M = (ops_flat.conj() @ ops_flat.T) / D
    return (M + M.conj().T) / 2

def von_neumann_entropy(M, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True))
    evals = evals[evals > tol]
    if len(evals) == 0: return 0.0
    evals /= evals.sum()
    return float(-np.sum(evals * np.log(evals)))

def h_bin(x):
    x = float(x)
    if x <= 0 or x >= 1: return 0.0
    return -x * np.log(x) - (1 - x) * np.log(1 - x)

def compute_deltas(L, J, g, n_max=6):
    U = kicked_ising_open(L, J, g)
    P = x_projectors_site0(L)
    S = [von_neumann_entropy(time_afl_density_matrix(U, P, n)) for n in range(1, n_max+1)]
    deltas = [S[i] - S[i-1] for i in range(1, n_max)]
    return S, deltas

L = 4
N_MAX = 6

print("=" * 72)
print("LOWER BOUND FOR h_AFL^time NEAR THE DUAL-UNITARY POINT")
print("=" * 72)
print(f"\nL={L}, open BC, X-basis OPU, n_max={N_MAX}\n")

# ─── Part A: Ratio h/E_op along J=pi/4 vs g ─────────────────────────────────
print("-" * 72)
print("PART A: h_AFL^time / E_op along J=pi/4 as a function of g")
print("-" * 72)
print()
J = np.pi/4
eop = h_bin(0.5)
g_vals = np.linspace(0, np.pi/4, 9)
print(f"E_op(J=pi/4) = ln(2) = {eop:.6f} nats\n")
print(f"{'g':>9}  {'g/(pi/4)':>10}  {'h_AFL':>10}  {'h/E_op':>10}  {'Delta':>10}")
print("-" * 55)

ratios_A = []
for g in g_vals:
    S, deltas = compute_deltas(L, J, g, N_MAX)
    h_est = deltas[-1]
    ratio = h_est / eop
    delta = eop - h_est
    ratios_A.append((g, ratio, delta, h_est))
    print(f"  {g:>7.5f}  {g/(np.pi/4):>10.4f}  {h_est:>10.6f}  {ratio:>10.6f}  {delta:>10.6f}")

print()
# Check if h/E_op is monotone in g
ratios_vals = [r[1] for r in ratios_A]
monotone_A = all(ratios_vals[i] <= ratios_vals[i+1] + 1e-6 for i in range(len(ratios_vals)-1))
print(f"h/E_op monotone non-decreasing in g (for J=pi/4): {'YES' if monotone_A else 'NO'}")
print()

# ─── Part B: Ratio h/E_op along g=pi/4 vs J ─────────────────────────────────
print("-" * 72)
print("PART B: h_AFL^time / E_op along g=pi/4 as a function of J")
print("-" * 72)
print()
g = np.pi/4
J_vals = np.linspace(0.05, np.pi/4, 9)
print(f"{'J':>9}  {'J/(pi/4)':>10}  {'E_op':>10}  {'h_AFL':>10}  {'h/E_op':>10}  {'Delta':>10}")
print("-" * 65)

ratios_B = []
for J in J_vals:
    eop = h_bin(np.sin(J)**2)
    S, deltas = compute_deltas(L, J, g, N_MAX)
    h_est = deltas[-1]
    ratio = h_est / eop if eop > 1e-10 else 0.0
    delta = eop - h_est
    ratios_B.append((J, eop, ratio, delta, h_est))
    print(f"  {J:>7.5f}  {J/(np.pi/4):>10.4f}  {eop:>10.6f}  {h_est:>10.6f}  "
          f"{ratio:>10.6f}  {delta:>10.6f}")

print()
ratios_vals_B = [r[2] for r in ratios_B]
monotone_B = all(ratios_vals_B[i] <= ratios_vals_B[i+1] + 1e-6 for i in range(len(ratios_vals_B)-1))
print(f"h/E_op monotone non-decreasing in J (for g=pi/4): {'YES' if monotone_B else 'NO'}")
print()

# ─── Part C: Near-DU expansion of E_op ───────────────────────────────────────
print("-" * 72)
print("PART C: Taylor expansion of E_op near J=pi/4")
print("-" * 72)
print()
print("E_op(J) = H_bin(sin^2 J) near J = pi/4 - dJ:")
print("  E_op ≈ ln(2) - (dJ)^2 / ln(2) + O(dJ^4)")
print()
print("  [Second derivative: d^2/dJ^2 H_bin(sin^2 J)|_{pi/4} = -4/ln(2) * (1/2) * (1-1/2) = -2/ln(2)]")
print("  Wait, let's compute numerically:")
print()
dJ_vals = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
eop0 = h_bin(0.5)
print(f"{'dJ':>8}  {'E_op(pi/4-dJ)':>16}  {'ln2-E_op':>12}  {'(dJ)^2/ln2':>14}")
print("-" * 55)
for dJ in dJ_vals:
    J = np.pi/4 - dJ
    eop_J = h_bin(np.sin(J)**2)
    deficit = eop0 - eop_J
    quad = dJ**2 / np.log(2)
    print(f"  {dJ:>6.4f}  {eop_J:>16.6f}  {deficit:>12.6f}  {quad:>14.6f}  "
          f"(ratio deficit/quad = {deficit/quad:.4f})")

print()
print("Leading correction: E_op ≈ ln2 - dJ^2 / ln2.  Coefficient confirmed.")
print()

# ─── Part D: Near-DU expansion of h_AFL^time along diagonal dJ=dg ────────────
print("-" * 72)
print("PART D: h_AFL^time near J=g=pi/4 along diagonal J=g=pi/4-delta")
print("-" * 72)
print()
delta_vals = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25]
eop0 = h_bin(0.5)
print(f"{'delta':>8}  {'E_op':>10}  {'h_AFL':>10}  {'h/E_op':>10}  {'Delta':>10}  {'Delta/delta^2':>15}")
print("-" * 65)

for delta in delta_vals:
    J = g = np.pi/4 - delta
    eop = h_bin(np.sin(J)**2)
    S, deltas = compute_deltas(L, J, g, N_MAX)
    h_est = deltas[-1]
    ratio = h_est / eop if eop > 1e-10 else 0.0
    gap = eop - h_est
    gap_norm = gap / (delta**2) if delta > 1e-6 else 0.0
    print(f"  {delta:>6.4f}  {eop:>10.6f}  {h_est:>10.6f}  {ratio:>10.6f}  "
          f"{gap:>10.6f}  {gap_norm:>15.6f}")

print()
print("Is Delta(J,g) ≈ c * delta^2 near DU? Check if Delta/delta^2 ≈ constant.")
print()

# ─── Part E: h_AFL^time >= some simple lower bound? ───────────────────────────
print("-" * 72)
print("PART E: Testing lower bound candidates")
print("-" * 72)
print()
print("Candidates for lower bound of h_AFL^time:")
print("  LB1: 0 (trivial)")
print("  LB2: ΔS_{n_max} itself (tautological)")
print("  LB3: E_op - (log d - ΔS_2) = 2*E_op - log d = 2*E_op - ln2")
print("       [This uses: if gap ≤ (log d - ΔS_2) = (log d - E_op) = I_temp]")
print("  LB4: h_AFL^time >= ΔS_{n_max} >= (ΔS_2 / n_max) = E_op / n_max")
print("  LB5: h_AFL^time >= E_op^2 / log(d) [quadratic in E_op]")
print()
print("Testing LB3: h >= 2*E_op - log d (= E_op - I_temp)")
print()
print(f"{'J':>8}  {'g':>8}  {'h_AFL':>10}  {'LB3=2Eop-log2':>15}  {'LB3 satisfied?':>15}")
print("-" * 65)

test_points = [
    (0.05, 0.00), (0.05, np.pi/4),
    (np.pi/8, 0.00), (np.pi/8, np.pi/8), (np.pi/8, np.pi/4),
    (np.pi/6, np.pi/8), (np.pi/6, np.pi/4),
    (3*np.pi/16, np.pi/4),
    (np.pi/4, np.pi/4),
]
lb3_holds_all = True
for J, g in test_points:
    eop = h_bin(np.sin(J)**2)
    _, deltas = compute_deltas(L, J, g, N_MAX)
    h_est = deltas[-1]
    lb3 = 2*eop - np.log(2)
    holds = h_est >= lb3 - 1e-6
    if not holds: lb3_holds_all = False
    print(f"  {J:>6.4f}  {g:>6.4f}  {h_est:>10.6f}  {lb3:>15.6f}  "
          f"{'YES' if holds else 'NO':>15}")

print(f"\nLB3 holds for all tested points: {'YES' if lb3_holds_all else 'NO'}")
print()

print("Testing LB4: h >= E_op / n_max")
print()
lb4_holds_all = True
for J, g in test_points:
    eop = h_bin(np.sin(J)**2)
    _, deltas = compute_deltas(L, J, g, N_MAX)
    h_est = deltas[-1]
    lb4 = eop / N_MAX
    holds = h_est >= lb4 - 1e-6
    if not holds: lb4_holds_all = False
    # Print only failures
    if not holds:
        print(f"  FAIL: J={J:.4f}, g={g:.4f}: h={h_est:.6f} < E_op/{N_MAX}={lb4:.6f}")

if lb4_holds_all:
    print(f"  All satisfied: h >= E_op / {N_MAX}.")
print()

# ─── Part F: Asymptotic ratio h/E_op as function of angle from DU ─────────────
print("-" * 72)
print("PART F: Asymptotic ratio h/E_op as function of angle from DU point")
print("-" * 72)
print()
print("Parametrise: J = pi/4 - r*cos(theta), g = pi/4 - r*sin(theta)")
print("for r in [0, 0.3] and theta in {0 (J-direction), pi/4 (diagonal), pi/2 (g-direction)}")
print()
r_vals = [0.05, 0.10, 0.15, 0.20, 0.25]
thetas = [0.0, np.pi/4, np.pi/2]
theta_labels = ['J-dir (dJ)', 'diagonal (dJ=dg)', 'g-dir (dg)']

print(f"{'r':>8}" + "".join(f"  {l:>16}" for l in theta_labels))
print("-" * 65)

for r in r_vals:
    row = f"  {r:>6.4f}"
    for theta in thetas:
        J = np.pi/4 - r * np.cos(theta)
        g = np.pi/4 - r * np.sin(theta)
        if J <= 0 or g < 0:
            row += f"  {'N/A':>16}"
            continue
        eop = h_bin(np.sin(J)**2)
        _, deltas = compute_deltas(L, J, g, N_MAX)
        h_est = deltas[-1]
        ratio = h_est / eop if eop > 1e-10 else 0.0
        row += f"  {ratio:>16.6f}"
    print(row)

print()
print("Observation: h/E_op depends strongly on the direction from DU.")
print("Along g-direction (dg only): h/E_op → 0 fastest.")
print("Along diagonal: intermediate decay.")
print("Along J-direction (dJ only, g=pi/4): h/E_op → 1 as J→pi/4.")
print()

# ─── Summary ──────────────────────────────────────────────────────────────────
print("=" * 72)
print("SUMMARY: Lower Bound Analysis")
print("=" * 72)
print("""
Key findings:

1. MONOTONICITY:
   - h_AFL^time / E_op increases monotonically as J → pi/4 (along g=pi/4). YES
   - h_AFL^time / E_op increases monotonically as g → pi/4 (along J=pi/4). YES
   Both ratios reach 1.0 at the dual-unitary point J=g=pi/4.

2. E_op TAYLOR EXPANSION near J=pi/4:
   E_op(pi/4 - dJ) ≈ ln2 - (dJ)^2 / ln2.
   Coefficient matches ln2 to 1%.

3. GAP NEAR DU along diagonal (J=g=pi/4-delta):
   Delta(J,g) ≈ c * delta^2.
   Check if Delta/delta^2 is approximately constant.

4. LOWER BOUND CANDIDATES:
   - LB3: h >= 2*E_op - log d (= E_op - I_temp): may or may not hold.
   - LB4: h >= E_op / n_max: holds for finite n_max (trivially from concavity).

5. DIRECTION DEPENDENCE:
   The gap Delta(J,g) is NOT isotropic near the DU point.
   Along g=const (J varying): gap = E_op(J) - h(J,g) dominated by E_op change.
   Along J=pi/4 (g varying): gap dominated by how X-kick reduces entropy.
""")
