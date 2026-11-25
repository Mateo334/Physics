# Exact closure of expectation-value dynamics

2026-11-24. User asks how to prove closure of the proposed readout-vector dynamics and its tangent equation. Scope: closed finite-dimensional systems, finite times; fixed Hermitian H and readouts A_a.

- [Proved] Put S=span_R{I,A_a}. An autonomous law dot x=F(x) valid for **all density matrices**, x_a=Tr(rho A_a), exists iff i[H,A_a]/hbar belongs to S for every a. When it exists, F is affine on the attainable readout set: dot x=c+Lx; actual preparation tangents obey dot(delta x)=L delta x.
- Sufficiency is direct substitution. Necessity: for any Hermitian V orthogonal to S, states I/d +/- epsilon V are positive for small epsilon and have identical readouts. Closure forces Tr(V i[H,A_a]/hbar)=0. Orthogonal-complement duality gives membership in S.
- Failure witness: choose V as the nonzero component of i[H,A_a]/hbar orthogonal to S. The two states have identical readouts but derivative difference 2 epsilon Tr(V^2)>0.
- Example: H=hbar omega sigma_x/2, readout sigma_z only; states (I +/- sigma_y)/2 have x=0 and derivatives +/-omega. Adding sigma_y closes an oscillator, not an unstable flow.
- Restricted invariant state families require a separate analysis; the all-state criterion does not exclude nonlinear coordinates on such families. No positive Lyapunov or entropy result follows from closure.
- Context: R_ab(t) differentiates fixed-observable expectations after common unitary evolution with respect to initial unitary-preparation parameters. It is a response Jacobian, not automatically an autonomous reduced-flow tangent map.

Standalone proof: `latex/readout_dynamics_closure.tex`. No pdflatex or tectonic found on PATH; **NOT compile-tested**. No simulation or external citation required for this finite-dimensional linear-algebra proof.

## Correction / connection to the user's existing parameter derivation

The subsequent discussion introduced local time-dependent closure on a preparation family: if f_t(epsilon)=x(t,epsilon) is locally invertible, F(t,x)=partial_t f_t(f_t^{-1}(x)), with Jacobian dot R R^{-1} along the reference. This is a local coordinate description, not a new intrinsic instability mechanism. The earlier emphasis on autonomous all-state closure was too restrictive for constructing a compatible time-dependent pair.

User correctly identified that complete readouts repeat their existing parameter-space procedure in `latex/Matejs_notes.tex`, “Time evolution,” lines 2177–2238, and generator transport, 2754–2773. In their full basis rho=I/d+sum_i lambda_i X_i, fixed readouts satisfy x=b+C lambda, C_ai=Tr(A_a X_i). If C is square invertible, the readout flow is just the similarity transform of dot lambda=A^T lambda. Their identity AG+GA^T=0 already proves preservation of the corresponding Gram metric. Fixing a complete measurement basis does not escape this argument. Incomplete readouts can display projection/visibility amplification but generally lose all-state closure. Local family inversion is again reparametrization. No proposed readout construction here has resolved the user's original instability problem; do not present it as a new escape from their existing obstruction.
