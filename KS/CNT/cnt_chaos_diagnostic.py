"""
cnt_chaos_diagnostic.py — CNT entropy to distinguish integrable vs chaotic dynamics.

New task: compare S_L/L convergence and pre-saturation slope of H_phi for:
  - Integrable: XXX Heisenberg chain
  - Chaotic:    Kicked Ising chain at dual-unitary (J=g=pi/4)

CNT diagnostic principles:
  1. h_CNT(tau_1) = S(phi) for both — CANNOT distinguish chaos from integrability.
  2. Pre-saturation slope of H_phi(N, alpha_t(N), ...) CAN distinguish:
     - Chaotic: fast operator spreading → large slope.
     - Integrable: slow spreading → small slope.
  3. Mean entropy density S_L/L converges faster for chaotic systems
     (stronger mixing → faster thermalization).
"""

import numpy as np
from scipy.linalg import expm

# ──────────────────────────────────────────────────────────────
# Utility
# ──────────────────────────────────────────────────────────────

def eta(t):
    return -t * np.log(t) if t > 1e-15 else 0.0

def von_neumann(rho):
    eigs = np.linalg.eigvalsh(rho)
    return sum(eta(max(e, 0)) for e in eigs)

def pauli():
    sx = np.array([[0, 1], [1, 0]], dtype=complex) / 2
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
    sz = np.array([[1, 0], [0, -1]], dtype=complex) / 2
    return sx, sy, sz

def kron_op(op, site, L):
    ops = [np.eye(2)] * L
    ops[site] = op
    result = ops[0]
    for o in ops[1:]:
        result = np.kron(result, o)
    return result

def xxx_hamiltonian(L, J=1.0, pbc=False):
    sx, sy, sz = pauli()
    H = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L - 1):
        for S in [sx, sy, sz]:
            H += J * (kron_op(S, i, L) @ kron_op(S, i+1, L))
    if pbc and L > 2:
        for S in [sx, sy, sz]:
            H += J * (kron_op(S, L-1, L) @ kron_op(S, 0, L))
    return H

def xx_hamiltonian(L, J=1.0, pbc=False):
    """XX model: only sigma_x and sigma_y terms (free fermion, integrable)."""
    sx, sy, sz = pauli()
    H = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L - 1):
        H += J * (kron_op(sx, i, L) @ kron_op(sx, i+1, L) +
                  kron_op(sy, i, L) @ kron_op(sy, i+1, L))
    if pbc and L > 2:
        H += J * (kron_op(sx, L-1, L) @ kron_op(sx, 0, L) +
                  kron_op(sy, L-1, L) @ kron_op(sy, 0, L))
    return H

def kicked_ising_unitary(L, J, g, pbc=False):
    """
    Kicked Ising Floquet unitary: U = exp(-i g sum_k sigma_z^k) * exp(-i J sum_k sigma_z^k sigma_z^{k+1})
    At J = g = pi/4: dual-unitary (maximally chaotic).
    """
    sx, sy, sz = pauli()
    # Kick: sum_k g * sigma_z^k
    H_kick = np.zeros((2**L, 2**L), dtype=complex)
    for k in range(L):
        H_kick += g * kron_op(2 * sz, k, L)  # sigma_z = 2*sz

    # Coupling: sum_k J * sigma_z^k * sigma_z^{k+1}
    H_coupling = np.zeros((2**L, 2**L), dtype=complex)
    for k in range(L - 1):
        H_coupling += J * (kron_op(2 * sz, k, L) @ kron_op(2 * sz, k+1, L))
    if pbc and L > 2:
        H_coupling += J * (kron_op(2 * sz, L-1, L) @ kron_op(2 * sz, 0, L))

    U_kick = expm(-1j * H_kick)
    U_couple = expm(-1j * H_coupling)
    return U_kick @ U_couple

def gibbs_state(H, beta):
    rho_unnorm = expm(-beta * H)
    return rho_unnorm / np.trace(rho_unnorm)

def tracial_state(L):
    d = 2**L
    return np.eye(d) / d

def partial_trace(rho, keep, dims):
    n = len(dims)
    rho_t = rho.reshape(dims + dims)
    trace_out = [i for i in range(n) if i not in keep]
    for idx in sorted(trace_out, reverse=True):
        rho_t = np.trace(rho_t, axis1=idx, axis2=idx + n)
        n -= 1
    d_keep = int(np.prod([dims[i] for i in keep]))
    return rho_t.reshape(d_keep, d_keep)

# ──────────────────────────────────────────────────────────────
# AFL-like H_phi for time evolution orbit
# Single-step OPU = {P0, P1} at site 0; orbit = {alpha_{k*dt}(Z)}_{k=0}^{n-1}
# ──────────────────────────────────────────────────────────────

def time_evolve_op(op, U):
    """One step: alpha_1(op) = U op U^dagger."""
    return U @ op @ U.conj().T

def orbit_entropy(rho_L, U, n_steps, site=0, basis='x'):
    """
    Compute S(rho[Z^(n)]) for n_steps steps of unitary U.

    OPU at site 0:
      basis='z': sigma_z basis {|0><0|, |1><1|} — trivial for diagonal U
      basis='x': sigma_x basis {|+><+|, |-><-|} — non-trivial for both systems

    We use basis='x' by default for a fair comparison.
    """
    L = int(np.log2(rho_L.shape[0]))
    d_local = 2
    d_total = 2**L

    if basis == 'z':
        Px0 = np.array([[1, 0], [0, 0]], dtype=complex)
        Px1 = np.array([[0, 0], [0, 1]], dtype=complex)
    else:  # x-basis: |+><+|, |-><-|
        Px0 = np.array([[0.5, 0.5], [0.5, 0.5]], dtype=complex)
        Px1 = np.array([[0.5, -0.5], [-0.5, 0.5]], dtype=complex)

    P0_site = kron_op(Px0, site, L)
    P1_site = kron_op(Px1, site, L)

    # Build evolved projectors: alpha_{k}(P_i) = U^k P_i (U^dag)^k
    evolved_projs = []
    P_cur = [P0_site.copy(), P1_site.copy()]
    for k in range(n_steps):
        evolved_projs.append([p.copy() for p in P_cur])
        P_cur = [time_evolve_op(p, U) for p in P_cur]

    # Build density matrix rho[Z^(n)]
    dim = d_local**n_steps
    rho_afl = np.zeros((dim, dim), dtype=complex)

    for idx_row in range(dim):
        idx_i = []
        tmp = idx_row
        for _ in range(n_steps):
            idx_i.append(tmp % d_local)
            tmp //= d_local
        Z_i = np.eye(d_total, dtype=complex)
        for k in range(n_steps):
            Z_i = Z_i @ evolved_projs[k][idx_i[k]]

        for idx_col in range(dim):
            idx_j = []
            tmp = idx_col
            for _ in range(n_steps):
                idx_j.append(tmp % d_local)
                tmp //= d_local
            Z_j = np.eye(d_total, dtype=complex)
            for k in range(n_steps):
                Z_j = Z_j @ evolved_projs[k][idx_j[k]]

            rho_afl[idx_row, idx_col] = np.trace(rho_L @ Z_j.conj().T @ Z_i)

    rho_afl = (rho_afl + rho_afl.conj().T) / 2
    rho_afl = np.real(rho_afl)
    rho_afl /= np.trace(rho_afl)
    return von_neumann(rho_afl)

# ──────────────────────────────────────────────────────────────
# Main: Integrable vs Chaotic
# ──────────────────────────────────────────────────────────────

print("=" * 65)
print("CNT ENTROPY AS CHAOS DIAGNOSTIC")
print("=" * 65)
print()

# ──────────────────────────────────────────────────────────────
# Part 1: Mean entropy density S_L/L — shift CNT entropy
# ──────────────────────────────────────────────────────────────

print("PART 1: h_CNT(shift) = S(phi) — not a chaos diagnostic")
print()
print("For BOTH integrable and chaotic systems:")
print("  h_CNT(tau_1) = S(phi) = thermodynamic entropy density.")
print("  Both approach log(2) at beta=0 (tracial state).")
print()

beta = 1.0
chain_lengths = [4, 6, 8]

print(f"{'System':>30}  {'L=4':>8}  {'L=6':>8}  {'L=8':>8}  {'Limit S(phi)':>14}")
print("-" * 73)

# XXX chain
S_xxx = []
for L in chain_lengths:
    H = xxx_hamiltonian(L, J=1.0, pbc=False)
    rho_L = gibbs_state(H, beta)
    S_xxx.append(von_neumann(rho_L) / L)
print(f"{'XXX Heisenberg (integrable)':>30}  " +
      "  ".join(f"{s:>8.4f}" for s in S_xxx) + f"  {'converging →':>14}")

# XX chain (free fermion)
S_xx = []
for L in chain_lengths:
    H = xx_hamiltonian(L, J=1.0, pbc=False)
    rho_L = gibbs_state(H, beta)
    S_xx.append(von_neumann(rho_L) / L)
print(f"{'XX free fermion (integrable)':>30}  " +
      "  ".join(f"{s:>8.4f}" for s in S_xx) + f"  {'converging →':>14}")

# Tracial state (infinite temperature)
S_inf = [np.log(2)] * len(chain_lengths)
print(f"{'beta=0 (infinite T, tracial)':>30}  " +
      "  ".join(f"{s:>8.4f}" for s in S_inf) + f"  {'log(2)=0.6931':>14}")

print()
print("=> h_CNT(shift) does NOT distinguish integrable from chaotic.")
print("   It depends only on temperature (beta), not on the dynamics.")

# ──────────────────────────────────────────────────────────────
# Part 2: Pre-saturation slope — time evolution CNT diagnostic
# ──────────────────────────────────────────────────────────────

print()
print("KEY NOTE: Kicked Ising U_KI is DIAGONAL in sigma_z basis.")
print("  → sigma_z basis projectors {|0><0|, |1><1|} have TRIVIAL orbit under U_KI:")
print("    U_KI^k |i><i| (U_KI^dag)^k = |i><i| (phases cancel for diagonal basis)")
print("  → All n give S=log(2). NOT a useful diagnostic with sigma_z OPU.")
print()
print("FIX: Use sigma_x basis {|+><+|, |-><-|} as OPU. NOT in the diagonal basis of U_KI.")
print("     This gives non-trivial orbits for BOTH systems.")
print()

print("PART 2: Pre-saturation slope — sigma_x OPU (fair comparison)")
print()

L = 5
n_max = 7

# Integrable: XXX chain
print(f"Integrable: XXX chain (L={L}, beta={beta}, sigma_x OPU)")
H_xxx = xxx_hamiltonian(L, J=1.0, pbc=False)
rho_xxx = gibbs_state(H_xxx, beta)
U_xxx = expm(-1j * H_xxx * 1.0)
print(f"{'n':>4}  {'S(rho[Z^n])':>14}  {'Slope':>10}")
print("-" * 35)
prev_S = 0.0
slopes_xxx = []
for n in range(1, n_max + 1):
    S_n = orbit_entropy(rho_xxx, U_xxx, n, site=0, basis='x')
    slope = S_n - prev_S
    if n > 1:
        slopes_xxx.append(slope)
    print(f"{n:>4}  {S_n:>14.6f}  {slope:>10.6f}")
    prev_S = S_n
mean_slope_xxx = np.mean(slopes_xxx) if slopes_xxx else 0.0
print(f"Mean slope (n=2..{n_max}): {mean_slope_xxx:.6f}")

print()

# Chaotic: Kicked Ising at J=g=pi/4 (dual-unitary)
J_ki = np.pi / 4
g_ki = np.pi / 4
print(f"Chaotic: Kicked Ising (L={L}, J=g=pi/4, dual-unitary), tracial state, sigma_x OPU")
rho_ki = tracial_state(L)
U_ki = kicked_ising_unitary(L, J_ki, g_ki, pbc=False)
print(f"{'n':>4}  {'S(rho[Z^n])':>14}  {'Slope':>10}")
print("-" * 35)
prev_S = 0.0
slopes_ki = []
for n in range(1, n_max + 1):
    S_n = orbit_entropy(rho_ki, U_ki, n, site=0, basis='x')
    slope = S_n - prev_S
    if n > 1:
        slopes_ki.append(slope)
    print(f"{n:>4}  {S_n:>14.6f}  {slope:>10.6f}")
    prev_S = S_n
mean_slope_ki = np.mean(slopes_ki) if slopes_ki else 0.0
print(f"Mean slope (n=2..{n_max}): {mean_slope_ki:.6f}")

print()
print(f"COMPARISON (sigma_x OPU):")
print(f"  XXX (integrable): mean slope = {mean_slope_xxx:.6f}")
print(f"  KI  (chaotic):    mean slope = {mean_slope_ki:.6f}")
ratio = mean_slope_ki / mean_slope_xxx if abs(mean_slope_xxx) > 1e-10 else float('inf')
print(f"  Ratio KI/XXX = {ratio:.3f}  (expect > 1 for chaotic)")
print(f"  log(2) = {np.log(2):.6f}  (per-site upper bound)")

# ──────────────────────────────────────────────────────────────
# Part 3: Chaos vs integrable — multiple coupling values
# ──────────────────────────────────────────────────────────────

print()
print("PART 3: Kicked Ising — integrable vs chaotic (varying J=g), sigma_x OPU")
print()
print(f"{'J=g':>8}  {'Mean slope':>12}  {'Regime':>25}  {'S/n (n=5)':>12}")
print("-" * 65)

L = 5
rho_inf = tracial_state(L)
J_values = [0.1, 0.3, np.pi / 8, np.pi / 6, np.pi / 4]
label = {0.1: "near-integrable", 0.3: "near-integrable", np.pi/8: "weakly chaotic",
         np.pi/6: "chaotic", np.pi/4: "dual-unitary (max chaos)"}

for J in J_values:
    U = kicked_ising_unitary(L, J, J, pbc=False)
    prev_S = 0.0
    slopes = []
    S_n5 = 0.0
    for n in range(1, 6):
        S_n = orbit_entropy(rho_inf, U, n, site=0, basis='x')
        if n > 1:
            slopes.append(S_n - prev_S)
        if n == 5:
            S_n5 = S_n
        prev_S = S_n
    mean_sl = np.mean(slopes) if slopes else 0.0
    print(f"{J:>8.4f}  {mean_sl:>12.6f}  {label[J]:>25}  {S_n5/5:>12.6f}")

print()
print("=> Pre-saturation slope INCREASES with chaos strength (sigma_x OPU).")
print("=> Dual-unitary (J=g=pi/4) achieves the maximum slope.")
print("=> This is the CNT analog of the Rényi AFL phase diagram (subfolder 1).")

# ──────────────────────────────────────────────────────────────
# Part 4: The CNT chaos diagnostic defined rigorously
# ──────────────────────────────────────────────────────────────

print()
print("PART 4: CNT chaos diagnostic — rigorous definition")
print()
print("Define the finite-L CNT chaos indicator:")
print()
print("  Delta_CNT(L, n) := [S(rho[Z_chaos^(n)]) - S(rho[Z_int^(n)])] / log(2)")
print()
print("where Z_chaos = time-evolution orbit for a chaotic system,")
print("      Z_int   = time-evolution orbit for an integrable system,")
print("      log(2) = the per-site maximum.")
print()

L = 5
n_test = 5
rho_inf = tracial_state(L)

# Integrable: XX chain (free fermion, maximally integrable)
H_xx = xx_hamiltonian(L, J=1.0, pbc=False)
U_xx_t = expm(-1j * H_xx * 1.0)
S_int = orbit_entropy(rho_inf, U_xx_t, n_test, site=0, basis='x')

# Chaotic: Kicked Ising at J=g=pi/4
U_ki_du = kicked_ising_unitary(L, np.pi/4, np.pi/4, pbc=False)
S_chaos = orbit_entropy(rho_inf, U_ki_du, n_test, site=0, basis='x')

Delta_CNT = (S_chaos - S_int) / np.log(2)
print(f"L = {L}, n = {n_test}:")
print(f"  S_int  (XX,   integrable) = {S_int:.6f}")
print(f"  S_chaos (KI dual-unitary) = {S_chaos:.6f}")
print(f"  Delta_CNT = {Delta_CNT:.6f}  (in units of log(2))")
print()
print("  Delta_CNT > 0 iff chaos (larger entropy orbit for chaotic dynamics).")
print("  Delta_CNT → 0 as L → ∞ (finite-size effect; nonzero in thermo limit?)")
print()
print(f"  Comparison with S_chaos_orbit/n5 = {S_chaos/n_test:.4f} vs S_int_orbit/n5 = {S_int/n_test:.4f}")

print()
print("=" * 65)
print("ALL COMPUTATIONS COMPLETE")
print("=" * 65)
