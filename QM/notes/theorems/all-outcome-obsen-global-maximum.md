# All finite POVMs: binary optimality of the instantaneous OE rate

[Proved, analytic derivation in this workspace; no literature-novelty claim] For finite integer d>=2, a pure state at a differentiable unitary-evolution instant, full-space trace volumes, natural logarithms, and fixed POVM effects during differentiation, every finite outcome count m>=2 has the same exact maximum:

`R_max,m=(DeltaH/hbar) C_d`, with `C_d=2sqrt(p(1-p))ln[(d-1)p/(1-p)]`, where the unique p in (1/2,1) solves `(2p-1)ln[(d-1)p/(1-p)]=2`.

This resolves the higher-outcome gap previously recorded in the binary theorem and numerical notes. Their earlier PGA outputs remain numerical results; the new analytic certificate independently identifies the exact optimum they matched. It does not make an unconverged numerical run globally convergent.

## Deliverable and scope

Full standalone proof: `latex/all_outcome_obsen_global_maximum.tex`. **NOT compile-tested:** no pdflatex or tectonic is available. Independent subagent algebra review checked the global certificate, scalar branches, complex effects, zero probabilities and trace budget. No papers were needed or cited.

Binary attainment uses P=|phi><phi|, phi=sqrt(p)|0>+sqrt(1-p)|1>, with |0>=psi and |1>=i(H-<H>)psi/DeltaH. For m>2, split I-P proportionally into m-1 nonzero effects. Zero energy uncertainty gives zero derivative.

## Global certificate, not a merging assumption

Normalize speed to one, so dotrho=-X, X=|0><1|+|1><0|. For every PSD effect set a=E00, b=Re E01, c=E11, V=Tr E. The rate is the sum of f(E)=2b ln(a/V), with zero-probability contributions defined by the actual zero derivative.

Put D=d-1, t=sqrt((1-p)/p), u=t². Stationarity gives ln(D/t²)=2(1+t²)/(1-t²). Define:

- `k=2t/[D(1+t²)]`
- `nu=2t/(1-t²)-2t³d/[D(1+t²)²]`
- `mu=-2ln[D(1+t²)/t²]+2nu/t`
- `lambda=nu/t²-2/t`

The common Hermitian certificate is `Lambda=lambda |0><0|+(mu/2)X+nu |1><1|+kI`. Its trace is `lambda+nu+kd=4t/(1-t²)=C_d`. The crucial global statement is `f(E)<=Tr(Lambda E)` for EVERY PSD effect, not merely first-order KKT stationarity.

For a>0, let x=b/a and y=V/a. PSD gives c/a>=x² and y>=1+x². Since nu>t>0, the effect inequality follows from

`G(x,y)=lambda+mu x+nu x²+ky+2x ln y >=0`, for all real x and y>=1+x².

### Complete scalar proof

1. Identities: `lambda+k=2t³[(3d-2)+(d-2)u]/[D(1-u)(1+u)²]>0`; `mu/2+ln(1+u)=-2u/(1-u)-2ud/[D(1+u)²]<0`; `nu/t-1 >= (1-u+7u²+u³)/[(1-u)(1+u)²]>0`; also 0<k<1.
2. For x>=0 minimize y at1+x². F(x)=G(x,1+x²) satisfies F(t)=F'(t)=0 and `F''=2(nu+k)+4x(x²+3)/(1+x²)²>0`. Thus F>=0.
3. For x=-r<0, the constrained minimizer is y=max(1+r²,2r/k). Let r_-<1<r_+ solve k(1+r²)=2r. Since k(1+t²)=2t/D<=2t, r_-<=t.
4. For r<=r_-, the minimum is `g(r)=lambda+k-r[mu+2ln(1+r²)]+(nu+k)r²>0`, using step1 and r<=t.
5. For every r>0, unconstrained minimization over y>0 gives the lower bound `h(r)=lambda-mu r+nu r²+2r-2r ln(2r/k)`. At r_-, h=g>0. At r*=1/t>r_-, h=h'=0; h''=2nu-2/r changes sign exactly once and h''(r*)>0. Thus r* is the only local minimum; any other stationary point is a maximum. Since h tends to +infinity, its minimum on [r_-,infinity) is at r_- or r*. Both are nonnegative. The actual constrained gap is always >=h, even when the unconstrained minimizing y is infeasible.
6. This handles d2 without any special-case assumption on the sign of lambda. For a=0, PSD gives b=0 and Tr(Lambda E)=nu c+kV>=0. Along a differentiable pure-state lift the outcome probability is O(time²), so its entropy derivative is zero.
7. Summing the effect inequalities and using completeness yields rate<=Tr Lambda=C_d, independently of outcome count. Binary attainment and proportional splitting prove equality for every finite m>=2.

The proof does not require effects to have isotropic complement blocks; arbitrary complex effects and all cross-complement couplings are included. No rate-nondecreasing pair-merging lemma is asserted. This is a direct replacement for the proposed induction.

## Numerical sanity check

[Numerically tested, not a premise] `calculations/code/test_obsen_multioutcome_dual.py`; output `calculations/results/obsen_multioutcome_dual_test.json`. Float64, SciPy brentq xtol1e-14, deterministic dimensions 2,3,4,8,21,256,4096. At each dimension, 30,001 uniform x points on [-30,30] plus both exact contact x values, with y minimized analytically; local bounded minimization near sampled minima. Assert tolerance -1e-9, smallest sampled gap approximately -2.24e-14 (roundoff), with contacts at x=t and x=-1/t. This finite-domain numerical check is not a global proof. Python3.11.9, NumPy2.2.3, SciPy1.15.1.

## Limitations

Not a mixed-state, common-ensemble-POVM, restricted-measurement, finite-time, continuous-outcome, ALF/CNT, or chaos theorem. No derivative at a delta-kick discontinuity is claimed. Effects are selected for a chosen instant and held fixed during differentiation.
