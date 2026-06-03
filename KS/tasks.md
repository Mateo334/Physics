# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the notes.md file in each subfolder you are currently working in, mateo will regularly write what you should focus on. If there is none, create one
- Reproduce all derivations step by step in LaTeX, without skipping any steps.
- Attempt to derive something new or extend the result — a sharper bound, a generalisation, a connection to another entropy measure, or an application to a specific class of quantum channels.
- Try to work with POVM, CP maps, OTOC, anything that is connected even remotely to this concept.
- Use small easy systems for fact checking - include the analytical matrices and the derivation in the latex output.

- [x] CNT open problem: prove rigorously that S(rho[Z^(n)]) is subadditive for the AFL orbit under time evolution (not superadditive), and show when equality S_2=2*S_1 holds; connect to operator entanglement and the dual-unitary property
- [x] CNT extension: Prove analytically that E_op(u) = log(d) iff J=pi/4 (for the kicked Ising gate); found exact formula E_op = H_bin(sin^2(J)), proved independence from g; showed {E_op=log d} = {I_temp=0} = {J=pi/4} x [0,pi/2] (a line, not a single point); verified numerically.
- [x] CNT Pesin gap: Proved the Operator-Entanglement Pesin inequality h_AFL^time(J,g) <= E_op(u_J) = H_bin(sin^2(J)). Key steps: (1) ΔS_2 = E_op from complementarity law; (2) concavity of S_n (ΔS_n non-increasing) proved via SSA + marginal consistency; (3) therefore h_AFL^time ≤ ΔS_2 = E_op. Gap function Delta(J,g) computed on 5x5 grid (all ≥ 0); equality only at J=g=pi/4. Added Section 16 to CNT/Output.tex with full proofs and tables.
- [x] CNT lower Pesin bound: Proved near-DU lower bound h_AFL^time(J,g) >= E_op(J)*(1 - 4*r^2/ln2) where r=dist from DU. Key results: (1) E_op ≈ ln2 - 2*(dJ)^2 Taylor expansion; (2) gap Delta ≈ 4*r^2 is isotropic near DU (verified in all 3 directions, coefficient 4/ln2≈5.77); (3) implies h >= E_op/2 for r ≤ sqrt(ln2/8)≈0.294. No global lower bound found (LB3, LB4 both fail). Section 17 added to CNT/Output.tex.
- [ ] CNT second-marginal consistency: The SSA proof of concavity in Section 16 used that the partial trace of rho[Z^n] over the LAST measurement also gives rho[Z^{n-1}] (by time stationarity). Verify this analytically: prove that Tr_{last}[rho[Z^n]] = rho[Z^{n-1}] for time-stationary (kicked Ising) dynamics. If true, establish the two-sided marginal consistency and give a rigorous SSA proof. If false, find the correction term and its effect on the concavity bound. Add the rigorous proof (or correction) to Section 16 of CNT/Output.tex.
