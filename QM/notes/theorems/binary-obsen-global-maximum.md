# Binary POVMs: global instantaneous observational-entropy maximum

[Subsequent update] The higher-outcome gap mentioned in this original binary proof is now resolved: [all finite POVMs](all-outcome-obsen-global-maximum.md) proves the matching global upper bound by an individual-effect certificate. The binary argument below remains valid; statements that higher-outcome optimality was open describe its earlier status.

[Proved] For a pure state on a finite d-dimensional Hilbert space, d>=2, at a differentiable unitary-evolution instant with generator H/hbar, maximize the instantaneous derivative of observational entropy over ALL fixed binary POVMs {E,I-E}, 0<=E<=I, with full-space trace volumes and natural logarithms. Write s=DeltaH/hbar. Then R_max,2=s C_d,2, where

`C_d,2 = 2 sqrt(p*(1-p*)) ln[(d-1)p*/(1-p*)]`,

and p* is the unique root in (1/2,1) of `(2p*-1) ln[(d-1)p*/(1-p*)]=2`. If s=0, the rate is zero. The maximizer includes a rank-one projector P=|u><u|, u=sqrt(p*)|0>+sqrt(1-p*)|1>, where |0>=psi and |1>=i(H-<H>)psi/DeltaH. Its complement is the second outcome. This proves the previously conjectural binary optimality, including nonprojective binary effects; it does NOT establish optimality for arbitrary 4- or 10-outcome POVMs.

## Proof record

Standalone derivation: `latex/binary_obsen_global_maximum.tex`. It includes document preamble, explicit trace macro and a complete global domination proof plus KKT multipliers. **Not compile-tested:** neither pdflatex nor tectonic was discoverable with PowerShell Get-Command. No compiler installed. Independent analytic review by a subagent confirmed the reductions, scalar uniqueness argument, KKT signs and endpoint handling.

Normalize speed to one so dotrho=-(|0><1|+|1><0|). Set a=E00, c=E11, b=Re E01, V=Tr E. Rate r=2b L(a,V), L=ln[a(d-V)/((1-a)V)]. Relabel outcomes so L>=0. Nonpositive rates cannot beat the positive candidate; for remaining rates b>0. Positivity of E and I-E implies b²<=min(ac,(1-a)(1-c)). Discard cross-complement couplings and replace the complement block by zero; the 2x2 principal block remains a positive contraction. This decreases V to a+c and increases L, preserving b. Saturate b positively to obtain an upper bound attained within the reduced problem.

- If a+c>=1, decreasing c to1-a increases both sqrt((1-a)(1-c)) and L. The resulting rank-one projector dominates.
- If t=a+c<=1, saturation gives E=tP, a=tp. The rate is `2t sqrt(p(1-p)) ln[p(d-t)/(1-tp)]`. The log is positive iff dp>1. At fixed such p, differentiating the factor t times log yields `ln[p(d-t)/(1-tp)] + t(dp-1)/[(d-t)(1-tp)] >0`. Therefore t=1 dominates.
- It remains to maximize f(p)=2sqrt(p(1-p))L_d(p) on [1/d,1]. Its derivative is `[2-(2p-1)L_d(p)]/sqrt(p(1-p))`. Positive up to p=1/2; above 1/2 the product (2p-1)L_d strictly increases from zero to infinity. Endpoints give zero. The unique critical point is the global maximum.

This domination argument supplies the global guarantee missing from KKT alone. At a=0 or1, pure-state smooth evolution makes the vanishing outcome probability O(time²), so its entropy derivative is zero. Zero-volume effects are zero operators and omitted. No claim for arbitrary mixed-state tangent paths.

## KKT certificate

For each effect: G_i=-ln(p_i/V_i)dotrho-(dotp_i/p_i)rho+(dotp_i/V_i)I. At the optimal P, put q=1-p, z=sqrt(pq), L=ln[(d-1)p/q], beta=2zd/(d-1). D=G1-G2=-L dotrho+(2/z)rho-beta I. In the u,w basis, w=-sqrt(q)|0>+sqrt(p)|1>, the off-diagonal element is (2p-1)L-2=0. D eigenvalues: lambda_plus=2p/[z(2p-1)]-beta>0 on u; lambda_minus=-2q/[z(2p-1)]-beta<0 on w; -beta on the remaining complement. Lambda=G2+lambda_plus P is a common Hermitian multiplier, with Lambda-G_i PSD and (Lambda-G_i)E_i=0. This proves the claimed KKT conditions explicitly without mistaking them for sufficient global conditions for a nonconcave objective.

## Numerical sanity checks

[Numerically tested] `calculations/code/check_binary_obsen_global_bound.py`; output `calculations/results/binary_obsen_global_bound_checks.json`. Seed 20261123, Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1. Root solved with brentq (xtol 1e-14). Dimensions 2,3,21,256,4096; 2000 random feasible reduced binary effects each, none exceeding the bound. KKT slack eigenvalues >=-3.11e-15, complementarity residuals <=1.70e-14. These tests are not the proof.

| d | p* | Exact binary coefficient (numerical evaluation) |
|---|---:|---:|
| 2 | 0.916778279800 | 1.325486838698 |
| 3 | 0.876696361781 | 1.745625661692 |
| 21 | 0.745629637182 | 3.546043546988 |
| 256 | 0.661050264168 | 5.878322996400 |
| 4096 | 0.613879998784 | 8.550380330460 |

## Higher outcomes: original gap and subsequent resolution

Proportional splitting retains the binary rate and entropy, so C_d,m>=C_d,2 for every m>=2. The prior PGA values match this lower bound within numerical tolerance. No matching upper bound for all higher-outcome POVMs has been proved here. An induction requires a rate-nondecreasing merging/replacement lemma; entropy monotonicity under coarse-graining alone does not establish derivative monotonicity.

[Correction to current status] The preceding paragraph records the original gap. The subsequent [all-outcome proof](all-outcome-obsen-global-maximum.md) establishes C_d,m=C_d,2 for every finite m>=2 without a pair-merging lemma.
