# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the notes.md file in each subfolder you are currently working in, mateo will regularly write what you should focus on. If there is none, create one
- Reproduce all derivations step by step in LaTeX, without skipping any steps.
- Attempt to derive something new or extend the result — a sharper bound, a generalisation, a connection to another entropy measure, or an application to a specific class of quantum channels.
- Try to work with POVM, CP maps, OTOC, anything that is connected even remotely to this concept.
- Use small easy systems for fact checking - include the analytical matrices and the derivation in the latex output.


- [x] Renyi-alpha Pesin inequality without SSA — prove or find counterexample: The α=1 Pesin bound h_AFL ≤ E_op uses SSA + marginal consistency (Section 16). For α≠1, SSA fails. (1) For α=2 specifically: try to prove h_2^AFL ≤ E_op^(2) = (1/(1-2))*log(cos^4J+sin^4J) rigorously using the Rényi-2 entropy formula S_2^(2) = -log Tr[rho^2]. (2) Use the Rényi-2 conditional entropy: S_2(A|B) = S_2(AB) - S_2(B) (may differ sign from classical case). (3) Attempt to adapt the marginal consistency proof to bound ΔS^(2)_n (Rényi-2 increments) using the matrix inequality Tr[rho^alpha] properties. (4) Alternatively, construct a potential counterexample: find a quantum channel where h_2^AFL > E_op^(2) (if SSA genuinely fails). (5) Numerical search: on a 5x5 (J,g) grid for d=2, L=4, compute h_2^AFL and E_op^(2) and check whether any violations occur. Add rigorous proof or counterexample as Section 32 to CNT/Output.tex. DONE: proved exactly for g=0 (geometric purity, h_alpha=E_op^alpha); verified on 5x5 grid (max violation 5.6e-16); weaker subadditivity bound proved; log-convexity conjecture stated (Section 32).

- [ ] Connecting the CNT entropy with either AFL entropy or KS entropy for some specific cases, if they generalize also. 
