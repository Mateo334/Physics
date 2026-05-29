# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the latex file, mateo will write regularly in the first section what the next step should be, what you should refine, etc. You are encouraged to reply to this if it was too difficult, or impossible or incorrect. 
## Task 9 — MIPT Order Parameter and Lindbladian Spectral Formula for AFL Entropy [COMPLETE — 2026-05-29]
- [x] Implement hybrid_circuit.py: L-qubit circuit with alternating Haar-random 2-qubit gates + projective measurements at rate p; compute h_AFL^time vs p for L=4,6; identify AFL phase transition at p_c; compare with entanglement entropy transition.
- [x] Derive the Lindbladian spectral formula: for Lindbladian with Lindblad operators {L_k, γ_k}, express h_AFL^open analytically in terms of the Lindblad spectrum; prove h_AFL^open → 0 in the strong-coupling limit and recover the unitary result at γ→0.
- [x] Prove (or rigorously argue) h_AFL^time is an MIPT order parameter: volume-law phase → h_AFL^time ~ log d; area-law phase → h_AFL^time → 0 as L → ∞; show the transition is sharp.
- [x] Add Sections 50–53 to Output.tex: hybrid circuit AFL (Sec 50), MIPT transition (Sec 51), Lindbladian spectral formula (Sec 52), unified quantum chaos phase diagram (Sec 53).

## Task 5 — Time-Dynamical AFL Entropy: True Quantum Chaos Indicator via Time Evolution OPUs [COMPLETE — 2026-05-29]
- [x] Prove that the time-refined OPU Z^(n) with Θ_t (time evolution) is a valid OPU; show OPU condition holds.
- [x] Write time_evolution_afl.py: compute AFL entropy S(n) with time evolution for kicked Ising and XX chains; show S(n)/n → log d for chaotic, S(n)/n < log d for integrable.
- [x] Derive closed-form formula for n=2 density matrix ρ[Z^(2)] in terms of OTOC; connect time-AFL entropy at n=2 to operator scrambling.
- [x] State and prove (or conjecture rigorously): h_AFL^time = v_E = v_B log d for dual-unitary; h_AFL^time < v_B log d for integrable — the TRUE quantum Pesin.
- [x] Add Sections 35–38 to Output.tex: time-dynamical OPU framework (Sec 35), n=2 OTOC connection (Sec 36), numerical results kicked Ising vs XX (Sec 37), quantum Pesin via time evolution (Sec 38).

## Task 6 — Quantum Pesin Inequality: Rigorous Bounds and Phase Diagram [COMPLETE — 2026-05-29]
- [x] Prove h_AFL^time ≤ v_B log d rigorously (or as a rigorous conjecture with full supporting argument); prove the inequality is tight at dual-unitary.
- [x] Map the (coupling strength, h_AFL^time) phase diagram for kicked Ising chain: compute h_AFL^time vs J=g from 0 to π/4; identify the chaotic vs integrable crossover.
- [x] Prove the Qutrit Odd-Step K-Independence Conjecture analytically using Wigner d-matrix symmetry.
- [x] Add Sections 39–42 to Output.tex: rigorous bound proof (Sec 39), phase diagram numerics (Sec 40), qutrit conjecture proof (Sec 41), synthesis and final quantum Pesin picture (Sec 42).

## Task 7 — Open-System Quantum Chaos: Lindbladian Dynamics and Stinespring Dilation AFL [COMPLETE — 2026-05-29]
- [x] Formulate AFL entropy for open quantum systems via Stinespring dilation: define OPU for Lindbladian dynamics using the dilated unitary on system ⊗ environment.
- [x] Implement lindbladian_afl.py: compute time-AFL entropy for amplitude damping, dephasing, and depolarizing channels; show how dissipation reduces h_AFL^time.
- [x] Derive the entropy production formula: h_AFL^time(Lindblad) = h_AFL^time(Stinespring) − (dissipation rate); prove that the integrable-chaotic gap persists at weak noise.
- [x] Add Sections 43–46 to Output.tex: Stinespring OPU framework (Sec 43), Lindbladian entropy formula (Sec 44), numerical results for open-system chaos (Sec 45), connection to quantum error correction (Sec 46).

## Task 8 — Many-Body Noise Threshold and Global Depolarization [COMPLETE — 2026-05-29]
- [x] Derive the critical noise rate γ_c below which the integrable-chaotic AFL gap persists: study all-site depolarizing noise vs site-0 amplitude damping; show γ_c → 0 for global noise.
- [x] Implement noise_threshold.py: sweep noise strength γ for site-0 vs all-site channels; map the gap Δ_n(γ) and identify the threshold; compare KI vs XX chains.
- [x] Prove analytically: for all-site depolarizing at rate p, the gap collapses as Δ_n(p) ≈ Δ_n(0)·(1-p)^{Ln}; show global noise drives XX toward maximal entropy, collapsing the gap (not reducing KI).
- [x] Add Sections 47–49 to Output.tex: many-body threshold theory (Sec 47), numerical phase diagram (Sec 48), connection to topological order and scrambling (Sec 49).

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
