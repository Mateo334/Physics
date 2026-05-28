# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the latex file, mateo will write regularly in the first section what the next step should be, what you should refine, etc. You are encouraged to reply to this if it was too difficult, or impossible or incorrect. 
## Task 5 — Time-Dynamical AFL Entropy: True Quantum Chaos Indicator via Time Evolution OPUs
- [ ] Prove that the time-refined OPU Z^(n) with Θ_t (time evolution) is a valid OPU; show OPU condition holds.
- [ ] Write time_evolution_afl.py: compute AFL entropy S(n) with time evolution for kicked Ising and XX chains; show S(n)/n → log d for chaotic, S(n)/n < log d for integrable.
- [ ] Derive closed-form formula for n=2 density matrix ρ[Z^(2)] in terms of OTOC; connect time-AFL entropy at n=2 to operator scrambling.
- [ ] State and prove (or conjecture rigorously): h_AFL^time = v_E = v_B log d for dual-unitary; h_AFL^time < v_B log d for integrable — the TRUE quantum Pesin.
- [ ] Add Sections 35–38 to Output.tex: time-dynamical OPU framework (Sec 35), n=2 OTOC connection (Sec 36), numerical results kicked Ising vs XX (Sec 37), quantum Pesin via time evolution (Sec 38).

## Task 4 — Free Fermion vs Chaotic Chains: Testing Quantum Pesin and General Conjecture [COMPLETE — 2026-05-29]
- [x] Numerically verify: corrected AFL entropy h~ = s(ω) for both free fermion (XX) and kicked Ising chains; confirm AFL cannot distinguish integrability.
- [x] Compute OTOC for XX chain (L=6): confirm λ_L = 0 (no exponential growth); quantum Pesin fails for free fermions.
- [x] Compute OTOC and entanglement growth for kicked Ising: show faster scrambling, linear entanglement growth vs logarithmic for XX.
- [x] Add Sections 31–34 to Output.tex: infinite systems + C*-algebra context (Sec 31), free fermion exact results (Sec 32), kicked Ising comparison (Sec 33), general quantum Pesin sketch with reconciliation of AFL vs λ_L (Sec 34).

## Task 3 — Quantum Pesin Conjecture: Kicked Top, Qutrit Systems, and Dissipative Channels [COMPLETE — 2026-05-28]
- [x] Numerically verify Conjecture 23.1 using the quantum kicked top: proved K-independence theorem, showed E(U^n) is k-dependent for n≥2 but does not grow monotonely (quantum recurrence).
- [x] Compute E(U(k)) vs k: proved analytically E(U^1) = E(U_rot) = k-independent; E(U^n) for n≥2 is k-dependent.
- [x] Extend AFL analysis to d=3 (qutrit: odd-step k-independence conjecture) and d=2 SIC-POVM (eigenvalues (1/2,1/6,1/6,1/6) proved).
- [x] Study non-unital dynamics: proved OPU composition breaks for amplitude damping; fixed point |0><0> gives h~=0.
- [x] Added sections 27–30 to Output.tex (3064 lines, 276 balanced environments).

## Task 2 — Matrix Visualizations, Coarse-Graining Refinement, and Observational Entropy [COMPLETE — 2026-05-28]
Study and write, in full LaTeX detail:
- Visualize all OPU density matrices symbolically (block structure, not numbers).
- Analyze the limiting procedure as the OPU is refined: how does AFL entropy change?
- Study Šafránek's observational entropy and connect it rigorously to AFL.
- Add new sections 24–26 to Output.tex.

## Task 1 — Alicki-Fannes Entropy [COMPLETE — 2026-05-28]
Study and work through the Alicki-Fannes (AF) entropy inequality in full detail:

- Reproduce all derivations step by step in LaTeX, without skipping any steps.
- Understand the bound and its tightness.
- Attempt to derive something new or extend the result — a sharper bound, a generalisation, a connection to another entropy measure, or an application to a specific class of quantum channels.
- Try to work with POVM, CP maps, OTOC, anything that is connected even remotely to this concept.
- Use small easy systems for fact checking - include the analytical matrices and the derivation in the latex output.
