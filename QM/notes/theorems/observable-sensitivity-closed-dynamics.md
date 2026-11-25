# Observable sensitivity with an exact underlying state

2026-11-25. User explicitly asks to explore observables/probabilities while retaining exact closed-system rho dynamics, without forcing autonomous probability dynamics or a pre-Ehrenfest window.

## Setup

Finite dimension, fixed H and finite fixed nonzero-effect POVM E_i; differentiable initial preparation rho_epsilon, full-rank reference for Fisher formulas; V=partial_epsilon rho_epsilon. Evolve rho and V under the same commutator generator L. Readouts p_i=Tr(E_i rho), s_i=Tr(E_i V). This is exact at every time and requires no intermediate measurement or reduced closure. A restricted preparation orbit requires its own tangent space, not every traceless Hermitian V.

## Results

- [Proved] I_E=sum s_i²/p_i has derivative sum(2s_i dot s_i/p_i-s_i² dot p_i/p_i²). Its square-root logarithmic rate is dot I_E/(2I_E) where I_E>0. Readout rates are not identified with chaos exponents.
- [Proved, existing Fisher bound applied here] I_E<=F_Q(rho,V); common unitary conserves F_Q. The fraction eta=sqrt(I_E/F_Q) is in [0,1] for nonzero V. This is an accessibility measure, not a newly growing global information quantity. Self-contained SLD/Cauchy--Schwarz proof in TeX. Earlier local work: `fixed-povm-relative-entropy-derivative.md`.
- [Proved] Write CV=(Tr E_i V), T_t V=U_t V U_t†. For C nonzero on the chosen invariant ambient tangent space, sup_{CV!=0} ||CT_t V||/||CV|| is finite iff ker C is contained in ker(CT_t). Failure produces arbitrarily large ratios from V=Z+epsilon W, with CZ=0, CT_tZ!=0. Normalizing global tangent magnitude cannot fix this homogeneous ratio. For a restricted domain use the kernel inclusion on that actual domain.
- [Proved analytic control] H=hbar omega sigma_y/2, rho0=(I+r sigma_y)/2 with 0<r<1, V0=(a sigma_x+b sigma_z)/2, sigma_x PVM. The reference is stationary; V is a valid unitary-orbit tangent. p±=1/2, s±=±(a cos omega t+b sin omega t)/2, I_E=(a cos+b sin)², F_Q=a²+b². Arbitrary readout amplification from a small a occurs under periodic regular rotation. No numerical test needed.
- [Proved from finite transport duality] For a specified metric d_ij on outcomes, ||s||_tr=max_{|f_i-f_j|<=d_ij} sum f_i s_i is a norm for zero-sum s. W1(p_epsilon,p0)=|epsilon| ||s||_tr+o(|epsilon|). Pullback to V is a seminorm for incomplete POVMs; it is exact for any time, not a heuristic autonomous probability flow. Physical ground metric must be justified; maximizing potentials can switch, so generic smooth Hessian reasoning is unavailable.
- [Proved] With fixed quadratic output weight W, M=C*WC and q=<V,MV>, dot q=<V,(ML-LM)V>, using HS skew-adjointness of L. If M>0 (IC readout and positive weights), length gain <=sqrt(cond M) uniformly in time. For fixed ambient V, this geometry has no reference-state dependence; it can enter through restricted preparations, normalization or state-dependent output geometry.

## Interpretation / next decision

[Unresolved] Which physically fixed POVM, ground metric and initial perturbation family distinguish chaotic dynamics from regular visibility changes? No claim of intrinsic full-state stretching, new autonomous golden braid, or ALF/CNT/KS/Pesin identity. Candidate choices are statistical accessibility (Fisher, already studied locally) and physically structured transport (new direction). Preserve exact rho and V; do not repeat the prior argument that nonclosure makes the programme impossible.

Standalone derivation: `latex/observable_sensitivity_closed_dynamics.tex`. No pdflatex or tectonic found on PATH on 2026-11-25; NOT compile-tested. No numerical simulations performed.
