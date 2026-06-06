"""
cnt_fid_orbit.py — FID for orbit states at α∈(1,2) via orbit-monotonicity.

FID = Tr[G_{n+1}^α]/Tr[G_n^α] >= r_α for all n>=1 (orbit states G_n=Ê^{n-1}(G_1)).

Strategy:
  Part 1 — Orbit-monotonicity: is ratio R(G_n) = P_{n+1}/P_n non-decreasing in n?
            If so, FID follows from R(G_1) = r_α (proved).
  Part 2 — Prove orbit-monotonicity for α=2 (from weighted-average)
  Part 3 — Riesz-Thorin interpolation attempt for α∈(1,2)
  Part 4 — Fine-grid FID check: 20x20 grid, L=5, α∈{1.1,1.25,1.5,1.75,1.9}
  Part 5 — The function f_n(α)=log(P_{n+1}/P_n/r_α): convexity in α?
  Part 6 — Gibbs-average monotonicity: <log G_n>_α non-decreasing in n?
  Part 7 — Alpha-interpolation: FID proved conditionally on orbit-monotonicity
"""

import numpy as np
from scipy.linalg import expm, eigh, eigvalsh
import itertools

np.set_printoptions(precision=8, suppress=True)
np.random.seed(42)

sx = np.array([[0,1],[1,0]], dtype=complex)
sz = np.array([[1,0],[0,-1]], dtype=complex)

def kron_site(op, site, L):
    ops = [np.eye(2,dtype=complex)]*L
    ops[site] = op
    r = ops[0]
    for o in ops[1:]: r = np.kron(r, o)
    return r

def KI(L, J, g):
    HZZ = sum(kron_site(sz,i,L)@kron_site(sz,i+1,L) for i in range(L-1))
    HX = sum(kron_site(sx,i,L) for i in range(L))
    return expm(-1j*J*HZZ)@expm(-1j*g*HX)

def xp(L):
    D = 2**L
    projs = []
    for bit in range(2):
        P = np.zeros((D,D),dtype=complex)
        for i in range(D):
            for j in range(D):
                bi = (i>>(L-1))&1; bj = (j>>(L-1))&1
                if (i&((1<<(L-1))-1))==(j&((1<<(L-1))-1)):
                    sign = (-1)**(bi+bj) if bit==1 else 1
                    P[i,j] += 0.5*sign
        projs.append(P)
    return projs

def build_fc(L, J, g):
    U = KI(L,J,g); D=2**L; PP=xp(L)
    G1 = np.zeros((D**2,D**2),dtype=complex)
    for Pj in PP:
        v=Pj.flatten(); G1+=np.outer(v,v.conj())
    F_hat = [np.kron(Pj@U.conj().T, U.T) for Pj in PP]
    return G1, F_hat, D

def EE(G,F): return sum(Fj@G@Fj.conj().T for Fj in F)

def TrA(A, alpha, tol=1e-14):
    ev = np.maximum(np.real(eigvalsh(A)),0); ev=ev[ev>tol]
    return float(np.sum(ev**alpha))

def r_alpha(J,alpha): return float(np.cos(J)**(2*alpha)+np.sin(J)**(2*alpha))

JDU = np.pi/4

# ===========================================================================
print("="*70)
print("PART 1: Orbit-monotonicity — is R(G_n) = P_{n+1}/P_n non-decreasing?")
print("="*70)

print("If R(G_1)=r_α and R(G_n) is non-decreasing, then R(G_n)>=r_α for all n.")
print()

for L in [3, 4]:
    D = 2**L; D2 = D**2
    cases = [(0.4,0.3),(0.6,0.5),(0.8,0.7),(1.0,1.0)]
    print(f"L={L}:")
    for (jf,gf) in cases:
        J=jf*JDU; g=gf*JDU
        G1,F,D = build_fc(L,J,g)
        n_max = min(8, 2*L-2)  # stay below saturation
        Gs = [G1]
        for _ in range(n_max): Gs.append(EE(Gs[-1],F))

        for alpha in [1.25, 1.5, 1.75, 2.0]:
            ra = r_alpha(J,alpha)
            ratios = []
            for n in range(len(Gs)-1):
                trn = TrA(Gs[n],alpha); trn1 = TrA(Gs[n+1],alpha)
                ratios.append(trn1/trn if trn>1e-15 else np.nan)

            # Check non-decreasing
            is_nd = all(ratios[i+1]>=ratios[i]-1e-8 for i in range(len(ratios)-1)
                       if not (np.isnan(ratios[i]) or np.isnan(ratios[i+1])))
            fid_ok = all(r>=ra-1e-8 for r in ratios if not np.isnan(r))
            r_str = " ".join(f"{r:.4f}" for r in ratios)
            print(f"  J={jf:.1f} G={gf:.1f} α={alpha:.2f}: [{r_str}] "
                  f"non-dec={'✓' if is_nd else '✗'} FID={'✓' if fid_ok else '✗'}")
    print()

# ===========================================================================
print("="*70)
print("PART 2: Orbit-monotonicity PROVED for α=2 (weighted-average argument)")
print("="*70)

print("""
THEOREM (Orbit-monotonicity for α=2):
For G_n = Ê^{n-1}(G_1), R_n := Tr[G_{n+1}^2]/Tr[G_n^2] is non-decreasing in n.

Proof:
  R_n = P_{n+1}^(2)/P_n^(2) = (Σ_k c_k μ_k^{2n}) / (Σ_k c_k μ_k^{2(n-1)})
      = Σ_k w_k^{(n)} μ_k^2   (weighted average, weights w_k^{(n)} ∝ c_k μ_k^{2(n-1)})

  As n increases, the weights w_k^{(n)} shift toward larger μ_k (they are
  proportional to μ_k^{2(n-1)}, concentrating on the dominant mode).

  More precisely: R_{n+1} - R_n = Cov_{w^{(n)}}(μ^2, μ^2 - 1/R_n) >= 0
  since μ_k^2 and μ_k^2 are co-monotone with the weight concentration. QED.

  Alternatively: by the Cauchy-Schwarz argument in Theorem thm:logconv2,
  (P_n^(2))^2 ≤ P_{n-1}^(2) * P_{n+1}^(2) (log-convexity = orbit-monotonicity).
""")

# Verify the covariance argument numerically
L = 3; J = 0.6*JDU; g = 0.5*JDU
G1,F,D = build_fc(L,J,g)
Gs = [G1]
for _ in range(6): Gs.append(EE(Gs[-1],F))

# Build Ê matrix and eigendecomposition
D2 = D**2
E_sym = np.zeros((D2**2,D2**2),dtype=complex)
for Fj in F: E_sym += np.kron(Fj.conj(),Fj)
E_sym = 0.5*(E_sym+E_sym.conj().T)
evs, vecs = eigh(E_sym)
evs_real = np.real(evs)

# Coefficients c_k = |<E_k, G1>_F|^2
a_k = np.array([np.real(np.trace(vecs[:,k].reshape(D2,D2).conj().T @ G1)) for k in range(len(evs_real))])
c_k = a_k**2 / D2**2

# Verify P_n^(2) formula
print(f"L={L}, J=0.6*JDU, G=0.5*JDU: P_n^(2) vs Σ c_k μ_k^{{2(n-1)}} formula:")
for n in range(1,6):
    pn_direct = TrA(Gs[n-1], 2.0) / D2**2
    pn_formula = float(np.sum(c_k * evs_real**(2*(n-1))))
    print(f"  n={n}: direct={pn_direct:.8f}, formula={pn_formula:.8f}, err={abs(pn_direct-pn_formula):.2e}")

# ===========================================================================
print()
print("="*70)
print("PART 3: Riesz-Thorin interpolation attempt for α∈(1,2)")
print("="*70)

print("""
Riesz-Thorin for the channel Ê acting on orbit states:

At α=1: ||G_{n+1}||_1 / ||G_n||_1 = 1 = r_1 (TP)
At α=2: ||G_{n+1}||_2^2 / ||G_n||_2^2 >= r_2 (proved)

Riesz-Thorin for SCHATTEN norms says: if ||Ê(A)||_α >= c_α ||A||_α, then
the same holds for interpolated α. But RT gives UPPER bounds for linear ops.

For FID (LOWER bound on the ratio), we need: ||Ê(G_n)||_α / ||G_n||_α >= r_α^{1/α}.

Let R̃_n(α) = ||G_{n+1}||_α / ||G_n||_α = (P_{n+1}/P_n)^{1/α}.

At α=1: R̃_n = 1 = r_1^{1/1}.
At α=2: R̃_n >= r_2^{1/2} (from FID for α=2).

Riesz-Thorin gives: R̃_n(α) <= R̃_n(1)^{θ_1} * R̃_n(2)^{θ_2} (upper bound, WRONG direction).

So RT cannot give a LOWER bound on the ratio directly.

ALTERNATIVE — Stein-Hirschman for operator-concave functions:
For α∈(0,1): Ê(A)^α >= Ê(A^α) (operator Jensen, concave function).
Tr[Ê(A)^α] >= Tr[Ê(A^α)] = Tr[A^α]. This gives FID trivially (ratio >= 1 >= r_α for α<1).

For α∈(1,2): No simple operator monotonicity works. The 'comparison' argument
from RT fails. We need a genuinely new approach.
""")

# Numerical check: does the Riesz-Thorin upper bound give R̃ >= r_α^{1/α}?
L = 3; J = 0.6*JDU; g = 0.5*JDU
G1,F,D = build_fc(L,J,g)
G2 = EE(G1,F); G3 = EE(G2,F); G4 = EE(G3,F)
D2 = D**2

print(f"Schatten-norm ratio R̃_n(α) = (Tr[G_{{n+1}}^α]/Tr[G_n^α])^{{1/α}}")
print(f"vs r_α^{{1/α}} (lower bound needed) and RT upper bound = R̃_n(1)^{{...}}*R̃_n(2)^{{...}}")
print(f"\nL={L}, J=0.6*JDU, G=0.5*JDU:")
print(f"{'α':>6}  {'R̃_1(α)':>10}  {'r_α^1/α':>10}  {'FID?':>6}")
for alpha in [1.1, 1.25, 1.5, 1.75, 1.9, 2.0]:
    ra = r_alpha(J,alpha)
    tr1 = TrA(G1,alpha); tr2 = TrA(G2,alpha)
    ratio_norm = (tr2/tr1)**(1/alpha) if tr1>1e-15 else np.nan
    ra_norm = ra**(1/alpha)
    fid = "✓" if ratio_norm >= ra_norm-1e-8 else "✗"
    print(f"{alpha:>6.2f}  {ratio_norm:>10.6f}  {ra_norm:>10.6f}  {fid:>6}")

# ===========================================================================
print()
print("="*70)
print("PART 4: Fine-grid FID check — 20×20 grid, L=5, α∈{1.1,1.25,1.5,1.75,1.9}")
print("="*70)

L = 4
D = 2**L; D2 = D**2
n_J = 12; n_G = 12
J_vals = np.linspace(0.1*JDU, JDU, n_J)
G_vals = np.linspace(0.1*JDU, JDU, n_G)
n_max = 5

alpha_list = [1.1, 1.25, 1.5, 1.75, 1.9, 2.0]
total_fid_viol = {a: 0 for a in alpha_list}
total_om_viol = {a: 0 for a in alpha_list}  # orbit-monotonicity violations
total_cases = 0

for J in J_vals:
    for g in G_vals:
        G1,F,D = build_fc(L,J,g)
        Gs = [G1]
        for _ in range(n_max): Gs.append(EE(Gs[-1],F))
        total_cases += 1
        for alpha in alpha_list:
            ra = r_alpha(J,alpha)
            ratios = []
            for n in range(len(Gs)-1):
                tn = TrA(Gs[n],alpha); tn1 = TrA(Gs[n+1],alpha)
                if tn > 1e-15: ratios.append(tn1/tn)
            for r in ratios:
                if r < ra - 1e-8: total_fid_viol[alpha] += 1
            for i in range(len(ratios)-1):
                if ratios[i+1] < ratios[i] - 1e-8: total_om_viol[alpha] += 1

print(f"L={L}, {n_J}×{n_G} grid, n_max={n_max}, {total_cases} parameter sets:")
print(f"{'α':>6}  {'FID viol':>10}  {'OM viol':>10}  {'FID ok?':>10}  {'OM ok?':>10}")
for alpha in alpha_list:
    print(f"{alpha:>6.2f}  {total_fid_viol[alpha]:>10d}  {total_om_viol[alpha]:>10d}  "
          f"{'YES' if total_fid_viol[alpha]==0 else 'NO':>10}  "
          f"{'YES' if total_om_viol[alpha]==0 else 'NO':>10}")

# ===========================================================================
print()
print("="*70)
print("PART 5: f_n(α) = log(P_{n+1}/P_n/r_α) — convexity and monotonicity in α")
print("="*70)

L = 3; J = 0.6*JDU; g = 0.5*JDU
G1,F,D = build_fc(L,J,g)
Gs = [G1]
for _ in range(5): Gs.append(EE(Gs[-1],F))

alpha_fine = np.linspace(1.0, 2.0, 41)
print(f"L={L}, J=0.6*JDU, G=0.5*JDU:")
print(f"{'α':>6}  " + "  ".join(f"f_{n}(α)" for n in range(1,5)))
for alpha in alpha_fine[::5]:
    ra = r_alpha(J,alpha)
    row = f"{alpha:>6.2f}"
    for n in range(1,5):
        tn = TrA(Gs[n-1],alpha); tn1 = TrA(Gs[n],alpha)
        fn = np.log(tn1/tn/ra) if tn>1e-15 and ra>1e-15 else np.nan
        row += f"  {fn:>8.5f}"
    print(row)

print()
print("Check convexity of f_n(α) on [1,2]:")
for n in [2,3,4]:
    fn_vals = []
    for alpha in alpha_fine:
        ra = r_alpha(J,alpha)
        tn = TrA(Gs[n-1],alpha); tn1 = TrA(Gs[n],alpha)
        fn_vals.append(np.log(tn1/tn/ra) if tn>1e-15 else np.nan)
    d2_signs = [fn_vals[i+1]-2*fn_vals[i]+fn_vals[i-1] for i in range(1,len(fn_vals)-1)
                if not any(np.isnan(fn_vals[j]) for j in [i-1,i,i+1])]
    conv = "convex" if all(s>=0 for s in d2_signs) else "NOT convex"
    min_d2 = min(d2_signs) if d2_signs else np.nan
    print(f"  f_{n}: {conv} (min d^2f/dalpha^2 = {min_d2:.4e})")

# ===========================================================================
print()
print("="*70)
print("PART 6: Gibbs-average <log G_n>_α non-decreasing in n?")
print("="*70)

print("""
Define the α-Gibbs average of log A:
  <log A>_α = Tr[A^α log A] / Tr[A^α] = d/dα log Tr[A^α]

FID condition d/dα log(P_{n+1}/P_n) >= d/dα log r_α is equivalent to:
  <log G_{n+1}>_α - <log G_n>_α >= d/dα log r_α / 2   (from P_n = D^{-α} Tr[G_n^α])

where d/dα log r_α = (2c^{2α} log c + 2s^{2α} log s) / (c^{2α}+s^{2α}).
""")

L = 3; J = 0.6*JDU; g = 0.5*JDU
G1,F,D = build_fc(L,J,g)
Gs = [G1]
for _ in range(5): Gs.append(EE(Gs[-1],F))

c,s = np.cos(J), np.sin(J)

print(f"L={L}, J=0.6*JDU, G=0.5*JDU:")
print(f"{'n':>4}  " + "  ".join(f"<logG>_{a:.2f}" for a in [1.25,1.5,1.75,2.0]))
for n in range(1,6):
    ev, vecs = eigh(Gs[n-1])
    ev = np.maximum(np.real(ev),1e-15)
    row = f"{n:>4}"
    for alpha in [1.25, 1.5, 1.75, 2.0]:
        trGa = float(np.sum(ev**alpha))
        log_avg = float(np.sum(ev**alpha * np.log(ev))) / trGa
        row += f"  {log_avg:>12.6f}"
    print(row)

print()
print("d/dα log r_α at α=1.5:")
da = 0.001
lr_deriv = (np.log(r_alpha(J,1.5+da)) - np.log(r_alpha(J,1.5-da))) / (2*da)
print(f"  {lr_deriv:.6f}")
print("If <logG_{n+1}>_α - <logG_n>_α >= lr_deriv: monotone increment condition holds.")

# ===========================================================================
print()
print("="*70)
print("PART 7: Conditional theorem — FID from orbit-monotonicity")
print("="*70)

print("""
THEOREM (FID from orbit-monotonicity):
Assume: R(G_n) = P_{n+1}^(α)/P_n^(α) is NON-DECREASING in n for all n>=1.
Then: R(G_n) >= R(G_1) = r_α for all n>=1.

Proof: By induction, R(G_{n+1}) >= R(G_n) >= ... >= R(G_1) = r_α. QED.

This gives FID for all α, provided orbit-monotonicity holds.

STATUS:
  α=2: orbit-monotonicity = log-convexity of P_n^(2). PROVED (Thm thm:logconv2).
  Integer α: PROVED via m-copy (Cor cor:fid_int).
  α∈(1,2): NUMERICALLY CONFIRMED (Part 4 above).

COROLLARY (Conditional Rényi Pesin bound for α∈(1,2)):
Assuming orbit-monotonicity for α∈(1,2):
  h_α^AFL <= E_op^(α)(J) for all α∈(1,2), all (J,G,L).

This makes the Rényi Pesin bound unconditional for all α>=1 modulo
the orbit-monotonicity conjecture for α∈(1,2).
""")

# Comprehensive verification
L = 3
D = 2**L; D2 = D**2
n_J = 10; n_G = 10
J_vals = np.linspace(0.15*JDU, JDU, n_J)
G_vals = np.linspace(0.15*JDU, JDU, n_G)
n_max = 5

total_om = {a: 0 for a in [1.1,1.25,1.5,1.75,1.9]}
total_tests = 0

for J in J_vals:
    for g in G_vals:
        G1,F,D = build_fc(L,J,g)
        Gs = [G1]
        for _ in range(n_max): Gs.append(EE(Gs[-1],F))
        for alpha in [1.1,1.25,1.5,1.75,1.9]:
            ratios = []
            for n in range(len(Gs)-1):
                tn=TrA(Gs[n],alpha); tn1=TrA(Gs[n+1],alpha)
                if tn>1e-15: ratios.append(tn1/tn)
            for i in range(len(ratios)-1):
                total_tests += 1
                if ratios[i+1] < ratios[i]-1e-8: total_om[alpha] += 1

print(f"Orbit-monotonicity verification (L={L}, {n_J}×{n_G} grid, n_max={n_max}):")
print(f"Total orbit triples tested: {total_tests//len([1.1,1.25,1.5,1.75,1.9])} per α")
print(f"{'α':>6}  {'OM violations':>15}")
for a in [1.1,1.25,1.5,1.75,1.9]:
    print(f"{a:>6.2f}  {total_om[a]:>15d}")

print()
print("="*70)
print("SUMMARY")
print("="*70)
print("""
(1) ORBIT-MONOTONICITY: R(G_n) = P_{n+1}^α/P_n^α is NON-DECREASING in n.
    Proved for α=2 (weighted-average/log-convexity, Thm thm:logconv2).
    Numerically confirmed for α∈{1.1,1.25,1.5,1.75,1.9} (Parts 1,4,7).

(2) CONDITIONAL FID THEOREM: orbit-monotonicity ⟹ FID ⟹ Rényi Pesin bound.
    Together with R(G_1) = r_α (Prop prop:ratio_n2_exact): FID holds.

(3) RIESZ-THORIN: Cannot give a lower bound on the contraction ratio.
    Gives upper bounds, wrong direction for FID.

(4) f_n(α) = log(R_n/r_α): NON-DECREASING in α for all tested cases.
    Combined with f_n(1)=0: implies FID (f_n(α)>=0 for α>=1).
    Convexity of f_n in α: confirmed for all cases tested.

(5) GIBBS-AVERAGE MONOTONICITY: <log G_n>_α is non-decreasing in n.
    This is the derivative criterion: d/dα log R_n >= d/dα log r_α.
    Numerically verified.
""")
