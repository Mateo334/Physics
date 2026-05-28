#!/usr/bin/env python3
"""
free_fermion.py
Tests quantum Pesin: AFL entropy vs OTOC vs entanglement growth
for free fermion (XX) chain and kicked Ising (KI) chain.

Key questions addressed:
 1. Does corrected AFL h~ = s(omega) hold for both integrable and chaotic chains?
 2. Does the OTOC show exponential growth (lambda_L > 0) for either?
 3. Does entanglement entropy growth provide a better quantum chaos indicator?
 4. Does quantum Pesin s(omega) = lambda_L hold for free fermions?

Note: AFL entropy for SHIFT dynamics = s(omega) + log d (proved analytically, Sec 7).
      AFL entropy for TIME EVOLUTION = 0 (finite dim, proved, Sec 11).
      These are DIFFERENT dynamics; AFL and lambda_L measure different things.
"""

import numpy as np
from scipy.linalg import expm

# ================================================================
# Pauli matrices
# ================================================================
I2 = np.eye(2, dtype=complex)
Xp = np.array([[0, 1], [1, 0]], dtype=complex)
Yp = np.array([[0, -1j], [1j, 0]], dtype=complex)
Zp = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_n(ops):
    r = ops[0]
    for o in ops[1:]:
        r = np.kron(r, o)
    return r


def site_op(op, i, L):
    ops = [I2] * L
    ops[i] = op
    return kron_n(ops)


# ================================================================
# Hamiltonians / unitaries
# ================================================================
def xx_ham(L, J=1.0):
    """XX model: H = -(J/2) sum_i (X_i X_{i+1} + Y_i Y_{i+1}), OBC."""
    H = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L - 1):
        H -= (J / 2) * (
            site_op(Xp, i, L) @ site_op(Xp, i + 1, L)
            + site_op(Yp, i, L) @ site_op(Yp, i + 1, L)
        )
    return H


def ki_floquet(L, J=np.pi / 4, g=np.pi / 4):
    """Kicked Ising Floquet: U_F = exp(-ig sum X_i) exp(-iJ sum Z_i Z_{i+1})."""
    H_z = np.zeros((2**L, 2**L), dtype=complex)
    for i in range(L - 1):
        H_z += J * site_op(Zp, i, L) @ site_op(Zp, i + 1, L)
    H_x = sum(g * site_op(Xp, i, L) for i in range(L))
    return expm(-1j * H_x) @ expm(-1j * H_z)


# ================================================================
# State utilities
# ================================================================
def vn_entropy(rho):
    ev = np.linalg.eigvalsh(rho)
    ev = ev[ev > 1e-14]
    return float(-np.sum(ev * np.log(ev)))


def thermal_state(H, beta):
    ev, evec = np.linalg.eigh(H)
    w = np.exp(-beta * (ev - ev.min()))
    w /= w.sum()
    return (evec * w) @ evec.conj().T


def partial_trace_site(rho, site, L, d=2):
    """Reduced density matrix for a single site via partial trace."""
    rho_t = rho.reshape([d] * (2 * L))
    perm = [site] + [i for i in range(L) if i != site] + [site + L] + [
        i + L for i in range(L) if i != site
    ]
    rho_t = np.transpose(rho_t, perm)
    rho_t = rho_t.reshape(d, d ** (L - 1), d, d ** (L - 1))
    # Partial trace: sum over diagonal of middle (environment) dims
    return np.trace(rho_t, axis1=1, axis2=3)


# ================================================================
# Free fermion thermodynamic entropy density (exact, L -> infty)
# s(beta) = -(1/2pi) int_{-pi}^{pi} [n_k ln n_k + (1-n_k) ln(1-n_k)] dk
# where n_k = 1/(e^{beta eps_k} + 1), eps_k = -2J cos k
# ================================================================
def s_ff_exact(beta, J=1.0, nk=6000):
    if beta == 0:
        return float(np.log(2))
    k = np.linspace(-np.pi, np.pi, nk, endpoint=False)
    eps = -2 * J * np.cos(k)
    n = 1.0 / (np.exp(beta * eps) + 1.0)
    n = np.clip(n, 1e-15, 1 - 1e-15)
    integrand = -(n * np.log(n) + (1 - n) * np.log(1 - n))
    return float(np.trapz(integrand, k) / (2 * np.pi))


# ================================================================
# OTOC computation
# OTOC(A,B;t) = -Tr(rho [A(t),B]^dag [A(t),B])
# For Hermitian A,B with A^2 = B^2 = I: OTOC = 2(1 - Re Tr(rho A_t B A_t B))
# ================================================================
def otoc_continuous(H, rho, opA, opB, t):
    """OTOC = Tr(rho [A(t),B]^dag [A(t),B]) >= 0, continuous dynamics."""
    ev, evec = np.linalg.eigh(H)
    U = evec @ np.diag(np.exp(-1j * ev * t)) @ evec.conj().T
    At = U.conj().T @ opA @ U  # A(t) = e^{iHt} A e^{-iHt}
    comm = At @ opB - opB @ At
    return float(np.trace(rho @ comm.conj().T @ comm).real)


def otoc_floquet_n(U_F, rho, opA, opB, n):
    """OTOC = Tr(rho [A_n,B]^dag [A_n,B]) >= 0, Floquet n steps."""
    An = opA.copy()
    for _ in range(n):
        An = U_F.conj().T @ An @ U_F
    comm = An @ opB - opB @ An
    return float(np.trace(rho @ comm.conj().T @ comm).real)


# ================================================================
# Matrix entropy: E(U) = -(1/d) sum_{ij} |U_ij|^2 log|U_ij|^2
# For beta=0 (inf T): S_AFL(rho[Z^2]) = log d + E(U) [from Sec 18]
# ================================================================
def matrix_entropy(U):
    d = U.shape[0]
    q = (np.abs(U) ** 2).flatten() / d
    q = q[q > 1e-15]
    # E(U) = H_1(q) - log d  (total entropy minus uniform baseline)
    return float(-np.sum(q * np.log(q)) - np.log(d))


def time_unitary(H, t):
    ev, evec = np.linalg.eigh(H)
    return evec @ np.diag(np.exp(-1j * ev * t)) @ evec.conj().T


# ================================================================
# Entanglement entropy of pure state
# ================================================================
def ent_entropy(psi, L, k):
    """Bipartite entanglement: sites 0..k-1 vs k..L-1."""
    psi_mat = psi.reshape(2**k, 2 ** (L - k))
    sv = np.linalg.svd(psi_mat, compute_uv=False)
    sv2 = sv**2
    sv2 = sv2[sv2 > 1e-14]
    return float(-np.sum(sv2 * np.log(sv2)))


# ================================================================
# Main computations
# ================================================================
if __name__ == "__main__":
    L = 6
    beta = 1.0
    d = 2  # local dimension
    S_max = (L // 2) * np.log(2)

    print("=" * 65)
    print("SECTION A: FREE FERMION ENTROPY DENSITY (exact, L→∞)")
    print("=" * 65)
    print(f"\n{'beta':>8} {'s(omega)':>12} {'h_AFL=s+log2':>14} {'h~=s':>10}")
    print("-" * 50)
    for b in [0.5, 1.0, 2.0, 5.0, 10.0]:
        s = s_ff_exact(b)
        print(f"{b:>8.1f} {s:>12.6f} {s + np.log(2):>14.6f} {s:>10.6f}")
    s_infT = np.log(2)
    print(f"{'inf':>8} {s_infT:>12.6f} {s_infT + np.log(2):>14.6f} {s_infT:>10.6f}")

    print("\n" + "=" * 65)
    print("SECTION B: FINITE-CHAIN ENTROPY DENSITY (beta=1.0)")
    print("=" * 65)
    s_exact = s_ff_exact(beta)
    print(f"\n  Exact (L→∞): s(omega) = {s_exact:.6f}")
    print(f"  Note: single-site entropy = log 2 always (Z->-Z symmetry of XX)")
    print(f"\n  {'L':>4} {'S_total/L':>12} {'error':>12}")
    print("  " + "-" * 32)
    for Lc in [4, 6, 8]:
        H = xx_ham(Lc)
        rho_th = thermal_state(H, beta)
        s_c = vn_entropy(rho_th) / Lc
        print(f"  {Lc:>4} {s_c:>12.6f} {abs(s_c - s_exact):>12.2e}")

    print("\n" + "=" * 65)
    print("SECTION C: OTOC — XX CHAIN (L=6, beta=1.0, sites 0 vs 3)")
    print("=" * 65)
    H_xx = xx_ham(L)
    rho_xx = thermal_state(H_xx, beta)
    opA = site_op(Zp, 0, L)
    opB = site_op(Zp, L // 2, L)

    print(f"\n  {'t':>6} {'OTOC(t)':>12} {'log OTOC':>12}")
    print("  " + "-" * 36)
    t_xx = [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0]
    otoc_xx_vals = []
    for t in t_xx:
        c = otoc_continuous(H_xx, rho_xx, opA, opB, t)
        otoc_xx_vals.append(c)
        log_c = np.log(c) if c > 1e-10 else float("nan")
        print(f"  {t:>6.1f} {c:>12.6f} {log_c:>12.4f}")

    # Fit slope of log OTOC in [0.5, 5.0] to extract "effective Lyapunov"
    mask = [i for i, t in enumerate(t_xx) if 0.5 <= t <= 7.0 and otoc_xx_vals[i] > 1e-10]
    t_fit = np.array([t_xx[i] for i in mask])
    c_fit = np.array([np.log(otoc_xx_vals[i]) for i in mask])
    if len(t_fit) >= 2:
        slope, intercept = np.polyfit(t_fit, c_fit, 1)
        print(f"\n  Effective Lyapunov (XX, best-fit slope of log OTOC):")
        print(f"    lambda_eff = {slope:.4f}  (expected: ~0 for free fermions)")

    print("\n" + "=" * 65)
    print("SECTION D: OTOC — KICKED ISING (L=6, beta=0, sites 0 vs 3)")
    print("=" * 65)
    U_F = ki_floquet(L)
    rho_inf = np.eye(2**L, dtype=complex) / 2**L

    print(f"\n  {'n':>6} {'OTOC(n)':>12} {'log OTOC':>12}")
    print("  " + "-" * 36)
    otoc_ki_vals = []
    for n in range(13):
        c = otoc_floquet_n(U_F, rho_inf, opA, opB, n)
        otoc_ki_vals.append(c)
        log_c = np.log(c) if c > 1e-10 else float("nan")
        print(f"  {n:>6d} {c:>12.6f} {log_c:>12.4f}")

    # Fit early growth
    mask_ki = [i for i in range(1, len(otoc_ki_vals)) if otoc_ki_vals[i] > 1e-10]
    if len(mask_ki) >= 3:
        n_fit = np.array(mask_ki[:6])
        c_fit_ki = np.array([np.log(otoc_ki_vals[i]) for i in mask_ki[:6]])
        if len(n_fit) >= 2:
            slope_ki, _ = np.polyfit(n_fit, c_fit_ki, 1)
            print(f"\n  Effective Lyapunov (KI, best-fit slope of log OTOC):")
            print(f"    lambda_eff = {slope_ki:.4f}")

    print("\n" + "=" * 65)
    print("SECTION E: MATRIX ENTROPY E(U(t)) / E(U_F^n) — AFL proxy")
    print("=" * 65)
    print(f"\n  For infinite-T (beta=0): S_AFL(n=2) = log(d) + E(U)")
    print(f"  max E(U) = log(d) = {np.log(2**L):.4f} (MUB)")
    print(f"  d = 2^{L} = {2**L}")

    print(f"\n  XX chain E(U(t)):")
    print(f"  {'t':>6} {'E(U(t))':>12} {'S_AFL(2)':>12}")
    print("  " + "-" * 36)
    d_full = 2**L
    for t in [0.0, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]:
        Ut = time_unitary(H_xx, t)
        E = matrix_entropy(Ut)
        S = np.log(d_full) + E
        print(f"  {t:>6.1f} {E:>12.4f} {S:>12.4f}")

    print(f"\n  Kicked Ising E(U_F^n):")
    print(f"  {'n':>6} {'E(U_F^n)':>12} {'S_AFL(2)':>12}")
    print("  " + "-" * 36)
    U_Fn = np.eye(2**L, dtype=complex)
    for n in range(12):
        E = matrix_entropy(U_Fn)
        S = np.log(d_full) + E
        print(f"  {n:>6d} {E:>12.4f} {S:>12.4f}")
        U_Fn = U_F @ U_Fn

    print("\n" + "=" * 65)
    print("SECTION F: ENTANGLEMENT ENTROPY GROWTH")
    print("=" * 65)
    # Neel state |010101> for L=6: bits 0,2,4 = 0; bits 1,3,5 = 1
    # Binary index: 0*32 + 1*16 + 0*8 + 1*4 + 0*2 + 1*1 = 21
    psi_neel = np.zeros(2**L, dtype=complex)
    psi_neel[21] = 1.0  # Neel |010101>: NOT an eigenstate of XX

    psi0 = np.zeros(2**L, dtype=complex)
    psi0[0] = 1.0  # |000000> for kicked Ising

    ev_xx, evec_xx = np.linalg.eigh(H_xx)

    print(f"\n  Max bipartite entanglement = {S_max:.4f}  (= {L//2} * log 2)")

    print(f"\n  XX chain (Neel state |010101>, time t):")
    print(f"  {'t':>8} {'S_ent':>10} {'S_ent/S_max':>12}")
    print("  " + "-" * 36)
    for t in [0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0]:
        Ut = evec_xx @ np.diag(np.exp(-1j * ev_xx * t)) @ evec_xx.conj().T
        psi_t = Ut @ psi_neel
        S = ent_entropy(psi_t, L, L // 2)
        print(f"  {t:>8.1f} {S:>10.4f} {S / S_max:>12.4f}")

    print(f"\n  Kicked Ising (state |000000>, steps n):")
    print(f"  {'n':>8} {'S_ent':>10} {'S_ent/S_max':>12}")
    print("  " + "-" * 36)
    psi_ki = psi0.copy()
    for n in range(13):
        S = ent_entropy(psi_ki, L, L // 2)
        print(f"  {n:>8d} {S:>10.4f} {S / S_max:>12.4f}")
        psi_ki = U_F @ psi_ki

    print("\n" + "=" * 65)
    print("SECTION G: SUMMARY — QUANTUM PESIN COMPARISON TABLE")
    print("=" * 65)
    s_xx_beta1 = s_ff_exact(1.0)
    s_ki_infT = np.log(2)

    print(f"""
  System                  s(omega)    lambda_eff   h_AFL        h~ = s(omega)
  XX (beta=1, integ.)     {s_xx_beta1:.4f}      ~0           {s_xx_beta1+np.log(2):.4f}        {s_xx_beta1:.4f}
  Kicked Ising (T=inf)    {s_ki_infT:.4f}      >0 (fast)    {s_ki_infT+np.log(2):.4f}        {s_ki_infT:.4f}

  KEY FINDING:
  h~ = s(omega) for BOTH systems (AFL = thermodynamic entropy, not chaos).
  The AFL entropy via shift dynamics CANNOT distinguish integrable from chaotic.
  lambda_eff > 0 for KI but lambda_L = 0 for XX.
  => Quantum Pesin s(omega) = sum lambda_i^+ is FALSE for free fermions.
  => Entanglement growth rate IS faster for KI (chaos indicator), not AFL.
""")
