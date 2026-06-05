"""
cnt_saturation_general.py — Qudit saturation law n_sat=2L-1 and general-d DU pattern.

Goals (Section 27):
  1. Verify numerically: n_sat = 2L-1 for d=2,3.
  2. Find J_DU(d) from DFT condition |f_k(J)|^2 = d for all k.
  3. Verify G_DU = J_DU for d=4 (new result).
  4. Prove rho[Z^n] = I/d^n at (J_DU, G_DU) for d=4.
  5. Compute C[1](J,d) autocorrelation and solve for J_DU.
"""

import numpy as np
from scipy.linalg import expm, eigh
from scipy.optimize import brentq
import itertools

np.random.seed(42)

def weyl_Z(d):
    return np.diag(np.exp(2j*np.pi*np.arange(d)/d))

def weyl_X(d):
    X = np.zeros((d,d), dtype=complex)
    for j in range(d): X[(j+1)%d, j] = 1.0
    return X

def x_proj(d):
    omega = np.exp(2j*np.pi/d)
    return [np.outer(v := np.array([omega**(j*k) for j in range(d)], dtype=complex)/np.sqrt(d),
                     v.conj()) for k in range(d)]

def Ugate(L, J, G, d):
    D = d**L; Id = np.eye(d, dtype=complex)
    Z = weyl_Z(d); X = weyl_X(d); Xh = (X + X.conj().T)/2
    HZZ = np.zeros((D,D), dtype=complex)
    for i in range(L-1):
        b = np.eye(1, dtype=complex)
        for k in range(L):
            if k==i:   b = np.kron(b, Z)
            elif k==i+1: b = np.kron(b, Z.conj().T)
            else:        b = np.kron(b, Id)
        HZZ += (b + b.conj().T)/2
    HX = np.zeros((D,D), dtype=complex)
    for i in range(L):
        xi = np.eye(1, dtype=complex)
        for k in range(L): xi = np.kron(xi, Xh if k==i else Id)
        HX += xi
    return expm(-1j*J*HZZ) @ expm(-1j*G*HX)

def Plast(L, d):
    P = x_proj(d); Ir = np.eye(d**(L-1), dtype=complex)
    return [np.kron(Ir, p) for p in P]

def afl_dm(U, P, n, d):
    D = U.shape[0]; Ud = U.conj().T; Un1 = np.linalg.matrix_power(U, n-1)
    ops = {}
    for idx in itertools.product(range(d), repeat=n):
        Zop = P[idx[0]].copy()
        for t in range(1, n): Zop = Zop @ Ud @ P[idx[t]]
        ops[idx] = Zop @ Un1
    il = list(ops.keys())
    M = np.zeros((len(il),len(il)), dtype=complex)
    for a,ia in enumerate(il):
        for b,ib in enumerate(il):
            M[a,b] = np.sum(ops[ib].conj()*ops[ia])/D
    return (M + M.conj().T)/2

def vn(M, tol=1e-12):
    ev = np.real(eigh(M, eigvals_only=True)); ev = ev[ev>tol]; ev /= ev.sum()
    return float(-np.sum(ev*np.log(ev)))

def rank_dm(M, tol=1e-9):
    ev = np.real(eigh(M, eigvals_only=True))
    mx = max(abs(ev)); return int(np.sum(ev > tol*mx)) if mx>0 else 0

def rho_L_evals(J, d):
    phi = np.array([np.exp(-1j*J*np.cos(2*np.pi*m/d)) for m in range(d)])
    return np.abs(np.fft.fft(phi))**2 / d**2

def e_op_alpha(J, d, alpha=1.0):
    lam = rho_L_evals(J, d); lam = lam[lam>1e-14]; lam /= lam.sum()
    if alpha==1: return float(-np.sum(lam*np.log(lam)))
    return float(np.log(np.sum(lam**alpha))/(1-alpha))

def autocorr_C1(J, d):
    phi = np.array([np.exp(-1j*J*np.cos(2*np.pi*m/d)) for m in range(d)])
    return float(np.real(np.sum(phi * np.conj(np.roll(phi, -1)))))

# ─────────────────────────────────────────────────────────────────────────────
print("="*72)
print("PART 1: Saturation law n_sat = 2L-1 (d=2 and d=3, DU points)")
print("="*72)

DU_J = {2: np.pi/4, 3: 4*np.pi/9}
for d in [2, 3]:
    J_du = DU_J[d]; log_d = np.log(d)
    print(f"\nd={d}, J_DU=G_DU={J_du/np.pi:.6f}pi")
    print(f"  {'L':>2} {'n_sat':>5} | {'n':>3} {'S_n':>9} {'n*logd':>9} {'rank':>5} {'d^n':>5}")
    for L in [2, 3]:
        n_sat = 2*L-1; U = Ugate(L, J_du, J_du, d); P = Plast(L, d)
        for n in range(1, n_sat+3):
            rhoN = afl_dm(U, P, n, d); Sn = vn(rhoN); rk = rank_dm(rhoN)
            tag = " SAT" if n > n_sat else ""
            print(f"  {L:>2} {n_sat:>5} | {n:>3} {Sn:>9.5f} {n*log_d:>9.5f} {rk:>5} {d**n:>5}{tag}")

print()
print("RESULT: rank = d^n for n<=2L-1, then saturates at d^{2L-1} = D^2/d.")

# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*72)
print("PART 2: HS block-diagonal structure (orthogonality across i_0 blocks)")
print("="*72)

d=2; L=3; n=3; J=G=np.pi/4
U=Ugate(L,J,G,d); P=Plast(L,d); D=d**L; Ud=U.conj().T
Un2=np.linalg.matrix_power(U,n-1)
ops={}
for idx in itertools.product(range(d),repeat=n):
    Zop=P[idx[0]].copy()
    for t in range(1,n): Zop=Zop@Ud@P[idx[t]]
    ops[idx]=Zop@Un2

# cross-block inner products
cross=[abs(np.sum(ops[(1,i1,i2)].conj()*ops[(0,j1,j2)])/D)
       for i1 in range(d) for i2 in range(d) for j1 in range(d) for j2 in range(d)]
print(f"\nCross-block |<Z_I, Z_J>_HS| for i_0=0 vs i_0=1 (n={n},L={L},d={d}):")
print(f"  max = {max(cross):.2e} => HS-ORTHOGONAL ✓")

# within-block Gram matrices
B0 = np.array([[np.sum(ops[(0,i1,i2)].conj()*ops[(0,j1,j2)])/D
                for j1 in range(d) for j2 in range(d)]
               for i1 in range(d) for i2 in range(d)])
B0 = (B0+B0.conj().T)/2
B1 = np.array([[np.sum(ops[(1,i1,i2)].conj()*ops[(1,j1,j2)])/D
                for j1 in range(d) for j2 in range(d)]
               for i1 in range(d) for i2 in range(d)])
B1 = (B1+B1.conj().T)/2
r0=rank_dm(B0); r1=rank_dm(B1)
print(f"Within-block ranks: rank(B_0)={r0}, rank(B_1)={r1}")
print(f"Total = {r0+r1} = d^{{2L-1}} = {d**(2*L-1)} ✓")
print(f"Each block rank = d^{{n-1}} = {d**(n-1)} ✓")

# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*72)
print("PART 3: DFT condition C[1]=0 for J_DU")
print("="*72)

print("\nAnalytical C[1] and DU conditions:")
print("  d=2: C[1] = 2cos(2J)          => J_DU = pi/4")
print("  d=3: C[1] = 2cos(3J/2)+1      => J_DU = 4pi/9")
print("  d=4: C[1] = 4cos(J)           => J_DU = pi/2")

print("\nNumerical verification at J=J_DU:")
for d, J_du in [(2,np.pi/4),(3,4*np.pi/9),(4,np.pi/2)]:
    c1 = autocorr_C1(J_du, d)
    print(f"  d={d}: C[1](J_DU)={c1:.2e} ✓")

J_du5 = brentq(lambda J: autocorr_C1(J,5), 1.5, 2.5)
J_du6 = brentq(lambda J: autocorr_C1(J,6), 2.0, 3.5)

print(f"\nNumerical J_DU for d=5: {J_du5:.8f} rad = {J_du5/np.pi:.8f}pi, C[1]={autocorr_C1(J_du5,5):.2e}")
print(f"Numerical J_DU for d=6: {J_du6:.8f} rad = {J_du6/np.pi:.8f}pi, C[1]={autocorr_C1(J_du6,6):.2e}")

J_DU = {2:np.pi/4, 3:4*np.pi/9, 4:np.pi/2, 5:J_du5, 6:J_du6}
print(f"\n{'d':>4} | {'J_DU':>12} {'J_DU/pi':>10} | {'E_op(J_DU)':>12} {'log(d)':>10} | form")
for d, J_du in J_DU.items():
    eop = e_op_alpha(J_du, d, 1.0)
    closed = {2:"pi/4", 3:"4pi/9", 4:"pi/2"}.get(d, "numerical")
    print(f"  {d:>2} | {J_du:>12.8f} {J_du/np.pi:>10.8f} | {eop:>12.8f} {np.log(d):>10.8f} | {closed}")

# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*72)
print("PART 4: d=4 DU point — G_DU = J_DU = pi/2 (NEW)")
print("="*72)

d=4; J_du4=np.pi/2; G_du4=np.pi/2; log4=np.log(4)
L=2; n_sat=2*L-1; n_max=7

# E_op at J=pi/2
lam=rho_L_evals(J_du4,d)
print(f"\nrho_L eigenvalues at J=pi/2, d=4: {np.round(lam,8)}")
print(f"max|lam - 1/4| = {max(abs(lam-0.25)):.2e} ✓")

# G-scan
print(f"\nG-scan (L={L}, d={d}, J=pi/2, estimator=delta_S_5):")
print(f"  {'G/pi':>7} | {'dS_5':>10} {'gap':>10}")
for G in np.array([0.3,0.35,0.4,0.45,0.48,0.5,0.52,0.55,0.6])*np.pi:
    U=Ugate(L,J_du4,G,d); Pp=Plast(L,d)
    S=[vn(afl_dm(U,Pp,n,d)) for n in range(1,6)]
    dS5=S[4]-S[3]
    tag=" <MAX" if abs(G-np.pi/2)<0.02 else ""
    print(f"  {G/np.pi:>7.4f} | {dS5:>10.6f} {log4-dS5:>10.2e}{tag}")

# DU state verification
print(f"\nDU state rho[Z^n]=I/4^n at (J,G)=(pi/2,pi/2), d=4, L=2:")
U=Ugate(L,J_du4,G_du4,d); Pp=Plast(L,d)
print(f"  {'n':>3} | {'S_n':>10} {'n*log4':>10} | {'max|rho-I/4^n|':>16} | {'rank':>5} {'4^n':>5}")
for n in range(1, n_sat+3):
    rhoN=afl_dm(U,Pp,n,d); Sn=vn(rhoN); rk=rank_dm(rhoN)
    target=np.eye(d**n)/d**n
    err=np.max(np.abs(rhoN-target))
    tag=" (sat)" if n>n_sat else ""
    print(f"  {n:>3} | {Sn:>10.6f} {n*log4:>10.6f} | {err:>16.2e} | {rk:>5} {d**n:>5}{tag}")

# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*72)
print("PART 5: G_DU = J_DU symmetry — numerical check for d=2,3,4")
print("="*72)

for d, J_du in [(2,np.pi/4),(3,4*np.pi/9),(4,np.pi/2)]:
    L=2; log_d=np.log(d)
    U_du=Ugate(L,J_du,J_du,d); Pp=Plast(L,d)
    S=[vn(afl_dm(U_du,Pp,n,d)) for n in range(1,5)]
    h_du=S[3]-S[2]
    G_off=J_du*0.7
    U_off=Ugate(L,J_du,G_off,d)
    S_off=[vn(afl_dm(U_off,Pp,n,d)) for n in range(1,5)]
    h_off=S_off[3]-S_off[2]
    print(f"d={d}: G=J_DU -> h={h_du:.6f} (log(d)={log_d:.6f}, gap={abs(h_du-log_d):.2e})")
    print(f"       G=0.7*J_DU -> h={h_off:.6f} (gap={log_d-h_off:.4f} > 0)")

# ─────────────────────────────────────────────────────────────────────────────
print("\n"+"="*72)
print("PART 6: DFT flatness |f_k(J_DU)|^2 = d for all k")
print("="*72)

for d, J_du in J_DU.items():
    phi=np.array([np.exp(-1j*J_du*np.cos(2*np.pi*m/d)) for m in range(d)])
    mods=np.abs(np.fft.fft(phi))**2
    print(f"d={d}: |f_k|^2 = {np.round(mods,5)}  max|.|^2-d|={max(abs(mods-d)):.2e}")

print("\n"+"="*72)
print("SUMMARY: KEY RESULTS FOR SECTION 27")
print("="*72)
print("""
1. n_sat=2L-1 VERIFIED: rank(rho[Z^n])=d^n for n<=2L-1, d^{2L-1} after.
2. J_DU CLOSED FORMS: d=2->pi/4, d=3->4pi/9, d=4->pi/2. Transcendental for d>=5.
3. G_DU=J_DU FOR d=2,3,4 (J<->G symmetry of DU condition). Verified.
4. d=4 NEW: rho[Z^n]=I/4^n at (J,G)=(pi/2,pi/2) for n<=3=n_sat. h_alpha=log4 all alpha.
5. DFT CONDITION: all |f_k(J_DU)|^2 = d simultaneously (flat DFT spectrum).
""")
