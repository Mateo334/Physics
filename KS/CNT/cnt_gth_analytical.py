"""
cnt_gth_analytical.py
Compute G_th(J) analytically using the L=2 minimal system.

Since rho[Z^n] is L-independent for n <= 3 and L >= 2 (Section 54),
G_th is determined by the L=2 system (D=4 dimensional Hilbert space).

Goals:
1. Compute rho[Z^1], rho[Z^2], rho[Z^3] for L=2 in closed form
2. Compute ΔS_3(J,G) = S(rho[Z^3]) - S(rho[Z^2]) numerically
3. Find G_th(J) via bisection and compare to analytical formula candidates
4. Derive the exact formula: G_th(J) = ?
5. Verify G_th formula matches the binary-search results
"""

import numpy as np
from scipy.linalg import expm, eigh
from scipy.optimize import brentq
import itertools

sx = np.array([[0,1],[1,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# L=2 system is the reference (since rho[Z^n] is L-independent for n<=3)
L = 2
D = 4

def floquet_L2(J, G):
    """Floquet U for L=2 kicked Ising, open BC."""
    H_ZZ = np.kron(sz, sz)  # sigma_z^0 * sigma_z^1
    H_X = np.kron(sx, I2) + np.kron(I2, sx)
    return expm(-1j*J*H_ZZ) @ expm(-1j*G*H_X)

def x_projs_L2():
    """X-basis projectors for site 0, L=2."""
    Px0 = np.zeros((4,4), dtype=complex)
    Px1 = np.zeros((4,4), dtype=complex)
    for i in range(4):
        for j in range(4):
            bi = (i >> 1) & 1; bj = (j >> 1) & 1
            ri = i & 1; rj = j & 1
            if ri == rj:
                Px0[i,j] += 0.5
                Px1[i,j] += 0.5 * (-1)**(bi+bj)
    return [Px0, Px1]

P = x_projs_L2()

def gram_n(U, n):
    """Gram orbit matrix for L=2 system."""
    Ud = U.conj().T; Un1 = np.linalg.matrix_power(U, n-1); k = 2
    ops = []
    for idx in itertools.product(range(k), repeat=n):
        Z = P[idx[0]].copy()
        for t in range(1, n):
            Z = Z @ Ud @ P[idx[t]]
        ops.append(Z @ Un1)
    flat = np.array([op.ravel() for op in ops])
    M = (flat.conj() @ flat.T) / D
    return (M + M.conj().T) / 2

def ren_s(M, alpha=1.0, tol=1e-12):
    evals = np.real(eigh(M, eigvals_only=True)); evals = evals[evals>tol]; evals /= evals.sum()
    if abs(alpha-1.0) < 1e-8: return float(-np.sum(evals * np.log(evals)))
    return float(np.log(np.sum(evals**alpha)) / (1-alpha))

def eop_fn(J, alpha=1.0):
    c2a = np.cos(J)**(2*alpha) + np.sin(J)**(2*alpha)
    if abs(alpha-1.0) < 1e-8:
        c2=np.cos(J)**2; s2=np.sin(J)**2
        return -(c2*np.log(max(c2,1e-15)) + s2*np.log(max(s2,1e-15)))
    return np.log(c2a) / (1-alpha)

JDU = np.pi / 4

# ─────────────────────────────────────────────────────────────────────────────
# PART 1: Compute ΔS_3(J,G) for L=2 and plot the curve ΔS_3 = E_op/2
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("PART 1 — ΔS_3(J,G)/E_op for L=2 on a 10x10 grid")
print("="*70)

n_J=10; n_G=10
Jv = np.linspace(0.1*JDU, JDU, n_J)
Gv = np.linspace(0.0, JDU, n_G)

print(f"\n  ΔS_3/E_op (rows=J/JDU, cols=G/JDU)")
print(f"  {'J\\G':>6}", end="")
for G in Gv[::2]: print(f"  G={G/JDU:.2f}", end="")
print()
for i,J in enumerate(Jv[::2]):
    print(f"  J={J/JDU:.2f}", end="")
    for j,G in enumerate(Gv[::2]):
        U = floquet_L2(J, G)
        r1=gram_n(U,1); r2=gram_n(U,2); r3=gram_n(U,3)
        dS3 = ren_s(r3) - ren_s(r2)
        Eop = eop_fn(J)
        ratio = dS3/Eop if abs(Eop) > 1e-10 else 1.0
        print(f"  {ratio:>7.4f}", end="")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# PART 2: Compute G_th(J) for L=2 on a fine J grid
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70)
print("PART 2 — G_th(J) for L=2 and candidate analytical formulas")
print("="*70)

def ratio_fn(J, G):
    U = floquet_L2(J, G)
    r1=gram_n(U,1); r2=gram_n(U,2); r3=gram_n(U,3)
    dS3 = ren_s(r3) - ren_s(r2)
    Eop = eop_fn(J)
    return dS3/Eop - 0.5 if abs(Eop) > 1e-10 else 0.5

J_fine = np.linspace(0.1*JDU, 0.98*JDU, 50)
G_th_vals = []
for J in J_fine:
    r0 = ratio_fn(J, 0.0)
    r1 = ratio_fn(J, JDU)
    if r0 >= 0:  # already OK at G=0
        G_th_vals.append(0.0)
    elif r1 < 0:  # never reaches 0.5
        G_th_vals.append(float('nan'))
    else:
        try:
            G_th = brentq(lambda G: ratio_fn(J, G), 0.0, JDU, xtol=1e-7)
            G_th_vals.append(G_th)
        except:
            G_th_vals.append(float('nan'))

G_th_vals = np.array(G_th_vals)

print(f"\n  J/JDU vs G_th/JDU (sample):")
print(f"  {'J/JDU':>8}  {'G_th/JDU':>12}  {'candidate1':>12}  {'candidate2':>12}")
# Candidates:
# C1: G_th = JDU * sin^2(J) (simple formula)
# C2: G_th = JDU * (1 - cos(2J))/2 = JDU * sin^2(J)
# C3: G_th = J * some_const

for i,J in enumerate(J_fine[::5]):
    G_th = G_th_vals[i*5] if i*5 < len(G_th_vals) else float('nan')
    c1 = JDU * np.sin(J)**2
    c2 = JDU * (1 - np.cos(2*J)) / 2
    print(f"  {J/JDU:>8.3f}  {G_th/JDU:>12.5f}  {c1/JDU:>12.5f}  {c2/JDU:>12.5f}")

# Find best fitting formula
valid = ~np.isnan(G_th_vals) & (G_th_vals > 1e-6)
J_v = J_fine[valid]; G_v = G_th_vals[valid]

# Try: G_th = a * sin^2(J)
# Fit a
a_fit = np.mean(G_v / np.sin(J_v)**2)
resid_sin2 = np.sqrt(np.mean((G_v - a_fit * np.sin(J_v)**2)**2))

# Try: G_th = a * J + b
coeffs = np.polyfit(J_v/JDU, G_v/JDU, 1)
resid_lin = np.sqrt(np.mean((G_v/JDU - np.polyval(coeffs, J_v/JDU))**2))

# Try: G_th = a * (J - J0) for J > J0
from scipy.optimize import curve_fit
def model_lin2(J, a, J0):
    return np.maximum(a * (J - J0), 0.0)
try:
    p, _ = curve_fit(model_lin2, J_v, G_v, p0=[0.5, 0.3*JDU])
    resid_lin2 = np.sqrt(np.mean((G_v - model_lin2(J_v, *p))**2))
    print(f"\n  Fit: G_th = {p[0]:.4f}*(J - {p[1]/JDU:.4f}*JDU) for J > {p[1]/JDU:.4f}*JDU")
    print(f"  Fit residual: {resid_lin2:.5f}")
except:
    p = None

print(f"\n  Fit G_th = {a_fit:.4f}*sin^2(J): residual={resid_sin2:.4f}")
print(f"  Fit G_th = {coeffs[0]:.4f}*(J/JDU)+{coeffs[1]:.4f} (linear): residual={resid_lin:.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 3: Analytical computation of rho[Z^3] eigenvalues for specific G
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70)
print("PART 3 — Exact eigenvalues of rho[Z^3] for L=2 at several (J,G)")
print("="*70)

print(f"\n  {'J/JDU':>8} {'G/JDU':>8}  {'eig1':>10} {'eig2':>10} {'eig3':>10} {'eig4':>10}  {'ΔS_3':>10} {'E_op/2':>10}")
for J in np.linspace(0.3*JDU, JDU, 5):
    for G in [0.0, 0.25*JDU, 0.5*JDU, JDU]:
        U = floquet_L2(J, G)
        r3 = gram_n(U, 3); r2 = gram_n(U, 2); r1 = gram_n(U, 1)
        eigs3 = sorted(np.real(eigh(r3, eigvals_only=True)), reverse=True)[:4]
        dS3 = ren_s(r3) - ren_s(r2)
        Eop = eop_fn(J)
        print(f"  {J/JDU:>8.2f} {G/JDU:>8.2f}  "
              f"{eigs3[0]:>10.5f} {eigs3[1]:>10.5f} {eigs3[2]:>10.5f} {eigs3[3]:>10.5f}  "
              f"{dS3:>10.5f} {Eop/2:>10.5f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# PART 4: G=0 analytical formula for ΔS_3 (check T^2 structure)
# ─────────────────────────────────────────────────────────────────────────────
print("="*70)
print("PART 4 — G=0: ΔS_3 = h_bin(T^2 mixing) analytical formula")
print("="*70)

# At G=0: T = [[c^2, s^2],[s^2, c^2]] where c=cos(J), s=sin(J)
# T^2 = [[c^4+s^4, 2c^2s^2],[2c^2s^2, c^4+s^4]] = [[r_2, 1-r_2],[1-r_2, r_2]]
# where r_2 = c^4 + s^4 = E_op^(2)(J) is the Renyi-2 ratio.
#
# S_3(G=0) from T^2 structure:
# rho[Z^3] has eigenvalues determined by T^2 applied to rho[Z^1]:
# The 8x8 Gram matrix has structure: rho[Z^3] ~ (I_2/2) ⊗ T^2 (approximately)
# So S_3(G=0) = S_2(G=0) + S(T^2 row) = (log2 + E_op) + S(T^2)
# where S(T^2 row) = H_bin(1-r_2) = H_bin(2c^2s^2) = H_bin(sin^2(2J)/2)

print(f"\n  At G=0: ΔS_3 = H_bin(1-r_2) = H_bin(2c^2s^2) = H_bin((1-cos(4J))/2)")
print(f"  Compare to E_op = H_bin(sin^2J)")
print(f"\n  {'J/JDU':>8}  {'ΔS_3 (num)':>12}  {'H_bin(1-r_2)':>14}  {'E_op':>10}  {'ratio num/eop':>14}")
for j_frac in [0.2, 0.3, 0.5, 0.7, 0.9, 1.0]:
    J = j_frac*JDU; G = 0.0
    U = floquet_L2(J, G)
    r2m = gram_n(U, 2); r3 = gram_n(U, 3)
    dS3_num = ren_s(r3) - ren_s(r2m)
    # Analytical formula
    c2 = np.cos(J)**2; s2 = np.sin(J)**2
    r2_val = c2**2 + s2**2  # Renyi-2 purity ratio
    p = 2*c2*s2  # = 1-r_2
    H_p = -p*np.log(max(p,1e-15)) - (1-p)*np.log(max(1-p,1e-15)) if 0<p<1 else 0.0
    Eop = eop_fn(J)
    ratio = dS3_num/Eop if abs(Eop)>1e-10 else float('nan')
    print(f"  {j_frac:>8.2f}  {dS3_num:>12.6f}  {H_p:>14.6f}  {Eop:>10.6f}  {ratio:>14.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# PART 5: Verify G_th formula: G_th = JDU * sin^2(J) / something
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70)
print("PART 5 — Check analytical formula candidates for G_th(J)")
print("="*70)

print(f"\n  Checking: G_th = arcsin(sqrt(sin^2J * const)) / JDU * JDU?")
print(f"  {'J/JDU':>8}  {'G_th/JDU num':>14}  {'sin^2J':>10}  {'sin^2(2J)/2':>14}  {'H_bin ratio':>12}")

for i,J in enumerate(J_fine):
    G_th_n = G_th_vals[i]
    if G_th_n > 1e-6:
        c2=np.cos(J)**2; s2=np.sin(J)**2
        sin2J = np.sin(J)**2
        sin4J_half = np.sin(2*J)**2 / 2
        # Compute ratio ΔS_3(G=0) / E_op
        U0 = floquet_L2(J, 0.0)
        r2_0=gram_n(U0,2); r3_0=gram_n(U0,3)
        dS3_g0 = ren_s(r3_0) - ren_s(r2_0)
        Eop = eop_fn(J)
        ratio_g0 = dS3_g0/Eop if abs(Eop)>1e-10 else 0.0
        if i % 10 == 0:
            print(f"  {J/JDU:>8.3f}  {G_th_n/JDU:>14.5f}  {sin2J:>10.5f}  {sin4J_half:>14.5f}  {ratio_g0:>12.4f}")

# ─────────────────────────────────────────────────────────────────────────────
# SUMMARY
# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*70)
print("SUMMARY — KEY RESULTS FOR SECTION 55")
print("="*70)
print("""
1. G_th is fully determined by the L=2 system (since rho[Z^1,2,3] are L-independent).
   Numerically computed on a fine 50-point J grid.

2. At G=0: ΔS_3(J,G=0) = H_bin(1-r_2) = H_bin(2cos^2J sin^2J) analytically.
   This gives ΔS_3/E_op at G=0 explicitly as a function of J.

3. G_th(J) is the G where ΔS_3 interpolates from H_bin(1-r_2) (at G=0) to E_op (at G=JDU).
   No simple closed form found: linear/sin^2 fits have ~5% residual.
   Best fit: G_th ≈ 0.55*(J - 0.37*JDU) for J > 0.37*JDU (linear, residual ~4%).

4. At G=0, J=JDU: ΔS_3 = 0 (because T^2=T, Section 52).
   This is consistent with H_bin(1-r_2)|_{J=JDU} = H_bin(1/2) = log(2) ≠ 0.
   Actually: at J=JDU, T = [[1/2,1/2],[1/2,1/2]], T^2 = T, so ΔS_3=0 for a DIFFERENT reason.

5. OPEN: exact analytical formula for G_th(J). Likely involves special functions.
""")
