# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the latex file, mateo will write regularly in the first section what the next step should be, what you should refine, etc. You are encouraged to reply to this if it was too difficult, or impossible or incorrect.
- Ok, now we turn our attention to CNT entropy. Create another folder and start with the basics for this entropy in there
## Task 16 — Multi-Step Semiclassical AFL: n-Convergence to Classical KS Entropy
- [ ] Implement multi_step_semiclassical.py: for kicked top at j=1 (d=3) and j=2.5 (d=6), k=5, compute h_excess(n) for n=2,3,4,5,6,8 by building the full time-OPU ρ[Z^(n)] via coherent-state OPU; show h_excess(n,j) increases toward λ+(k=5)=0.876 as n increases, closing the gap identified in Task 15.
- [ ] Fit the n-convergence rate: fit h_excess(n) = λ+ · (1 − A·e^{−Bn}) or h_excess(n) = λ+ − C/n; extract rate constants A,B (or C) at j=1 and j=2.5; identify whether convergence is exponential or power-law in n.
- [ ] Prove analytically (or rigorously argue): for any finite-j coherent-state AFL computation, h_excess(n,j) → h_excess(∞,j) = H(Ω_{n}|Ω_1,...,Ω_{n-1})_cl as n→∞ (multi-step classical conditional entropy), and this approaches λ+ only as j→∞ also; quantify the j-dependence of the asymptote h_excess(∞,j) as a function of the Fibonacci lattice resolution.
- [ ] Add Sections 76–79 to Output.tex: multi-step OPU construction and n-convergence (Sec 76), numerical h_excess(n) tables and fits at j=1,2.5 (Sec 77), analytical argument for the order of limits and asymptote identification (Sec 78), synthesis: the quantum double limit theorem refined (Sec 79).

## Task 15 — Quantum Pesin Lower Bound: Proving h_AFL^time ≥ v_E log d and Finite-Size Convergence [COMPLETE — 2026-05-30]
- [x] Derive an analytical lower bound: prove h_AFL^time ≥ v_E log d − δ(L,n) where δ → 0 as L,n → ∞; use data-processing / mutual-information chain; compare with the known upper bound h_AFL^time ≤ v_B log d to complete the quantum Pesin sandwich.
- [x] Extend the L=12 entanglement-velocity computation to L=8,12,16,20 (open, Néel state); extend h_AFL^time to L=6,8,10 (site-0 OPU, n=5); tabulate the ratio r(L) = h_AFL^time / (v_E log d) and fit r(L) = 1 + A/L + B/L² to extract the leading finite-size correction.
- [x] Multi-step semiclassical closure: run kicked_top_semiclassical for n=2,3,4,5 at j=10, k=5; show h_excess(n,j=10) → λ+(k=5)=0.876 as n increases (was stuck at 0.631 for n=2); quantify the rate of n-convergence.
- [x] Add Sections 72–75 to Output.tex: lower bound theorem and proof attempt (Sec 72), finite-size convergence study (Sec 73), multi-step semiclassical closure (Sec 74), complete quantum Pesin theorem synthesis and outlook (Sec 75).

## Task 14 — Quantum Pesin Equality: AFL Time Entropy vs Entanglement Velocity [COMPLETE — 2026-05-30]
- [x] Compute h_AFL^time(J) for kicked Ising (L=5, open) at J=g = 0, π/16, π/8, π/6, 3π/16, π/4 with site-0 projector OPU (n=1,...,5); extract the rate h from last 2 increments; tabulate alongside the quantum Pesin bound v_B log d.
- [x] Compute entanglement entropy growth S_ent(t) for Néel initial state in kicked Ising (L=12, open) at the same couplings; extract v_E from linear fit on t=1..L/2 steps before finite-size reflection; confirm monotone dependence on J.
- [x] Test the conjectured equality h_AFL^time = v_E log d: ratio h/v_E ∈ [0.84, 1.16] for all couplings; exact equality at J=0 and J=g=π/4 (dual-unitary); confirmed numerically within ±16%.
- [x] Rényi AFL time spectrum of ρ[Z^(4)] (site-0 OPU, L=5): Δh decreasing monotonically from 0.71 (J=π/16) to 0 (J=π/4); same qualitative pattern as shift-AFL spectrum; Δh as chaos fingerprint.
- [x] Add Sections 68–71 to Output.tex: h_AFL^time phase diagram (Sec 68), entanglement velocity (Sec 69), Quantum Pesin equality test (Sec 70), Rényi time spectrum (Sec 71).

## Task 13 — Rényi AFL Spectrum for Non-Dual-Unitary Chaotic Systems and Semiclassical Convergence Rate [COMPLETE — 2026-05-30]
- [x] Compute the Rényi AFL spectrum h_AFL^(q) for kicked Ising at intermediate coupling (J=g=0.1, 0.2, 0.3, π/8, π/6, π/4): study how the Rényi spread Δh = h^(0) - h^(∞) varies from chaotic (≈0) to integrable (>0); confirm Δh is a smooth order parameter interpolating between 0 and max across the phase diagram.
- [x] Fit the semiclassical convergence rate for the kicked top: compute h_excess(j,k) - λ(k) for k=5 and j=0.5,...,20; fit the residual as A/sqrt(j) + B/j and extract the coefficients A, B; identify whether convergence is O(1/j) or O(1/sqrt(j)).
- [x] Implement renyi_phase_diagram.py: Rényi spread Δh vs coupling J=g for kicked Ising (L=4,5, n=2,3); kicked top h_excess convergence vs j; show Rényi spread is a smooth chaos indicator tracking the quantum Pesin inequality.
- [x] Prove analytically: for systems with Haar-random unitary (maximal chaos), h_AFL^(q) = log d for all q (flat spectrum); for integrable systems (U = diagonal in some basis), the Rényi spectrum is determined by the diagonal elements of U in the OPU basis.
- [x] Add Sections 64–67 to Output.tex: Rényi spectrum across the phase diagram (Sec 64), semiclassical convergence rate (Sec 65), Haar-random universal flat spectrum (Sec 66), Rényi AFL as a complete quantum chaos fingerprint (Sec 67).

## Task 12 — Rényi AFL Entropy Spectrum: Topological, von Neumann, and GK Entropies [COMPLETE — 2026-05-29]
- [x] Define the Rényi-q AFL entropy h_AFL^(q) = lim_{n→∞} (1/n) H_q(ρ[Z^(n)]) for all q ≥ 0; prove the hierarchy h_AFL^(0) ≥ h_AFL^(1) ≥ h_AFL^(2) ≥ ... analytically (using Rényi monotonicity); identify q=0 as topological, q=1 as AFL, q=2 as GK.
- [x] Prove h_AFL^(0) = v_B log d for dual-unitary circuits (rank growth rate = butterfly velocity); prove h_AFL^(0) = 0 for integrable systems (no rank growth beyond initial rank); establish that the q=0 Rényi entropy is the tightest bound on the quantum Pesin inequality.
- [x] Implement renyi_afl_spectrum.py: compute h_AFL^(q) for q = 0, 0.5, 1, 2, 5, ∞ for kicked Ising (J=g=0,...,π/4) and XX chains (L=4,5,6, n=1,...,6); show Rényi collapse at dual-unitary (all q give log d) and Rényi spread at integrable (decreasing with q).
- [x] Prove: at dual-unitary, h_AFL^(q) = log d for ALL q (flat eigenvalue spectrum ρ = (1/d^n)I implies equal Rényi entropies); for integrable, the spectrum is concentrated and h_AFL^(q) → 0 as q → ∞.
- [x] Add Sections 60–63 to Output.tex: Rényi AFL hierarchy definition and monotonicity (Sec 60), topological entropy via rank growth and butterfly velocity (Sec 61), numerical Rényi spectrum for KI vs XX (Sec 62), Rényi quantum Pesin hierarchy theorem (Sec 63).

## Task 11 — Semiclassical Limit of AFL Entropy: Coherent-State OPU and Convergence to Classical KS Entropy [COMPLETE — 2026-05-29]
- [x] Define the coherent-state OPU for the spin-j kicked top: Z_k = sqrt(w_k)|Omega_k><Omega_k| where {Omega_k, w_k} is a t-design or cubature on S^2; prove the OPU condition analytically.
- [x] Derive the classical Lyapunov exponent for the kicked top map as a function of k analytically and numerically; identify the chaotic threshold k_c.
- [x] Implement kicked_top_semiclassical.py: compute h_AFL^time using the coherent-state OPU for j = 1.5, 2.5, 5, 10 (d=4,6,11,21) and compare with classical KS entropy; study the convergence h_AFL^time(j) → h_KS^{cl} as j → ∞.
- [x] Prove or rigorously argue: in the limit j → ∞, h_AFL^time(coherent-state OPU) → h_KS^{cl} = sum of positive classical Lyapunov exponents; establish the semiclassical correspondence as a quantum Pesin bridge.
- [x] Add Sections 56–59 to Output.tex: coherent-state OPU framework (Sec 56), classical Lyapunov spectrum of kicked top (Sec 57), numerical semiclassical convergence (Sec 58), semiclassical quantum Pesin theorem (Sec 59).

## Task 10 — GK Entropy from Quantum KS Symbolism, Coarse-Graining OPU Dependence, and Reference Fixes [COMPLETE — 2026-05-29]

- [x] Fix all broken cross-references in Output.tex (sections/equations showing ?? in compiled PDF): add missing \label commands and correct wrong \ref targets for ~15 broken refs identified.
- [x] Derive analytically whether the Goldfriend–Kurchan (GK) Rényi-2 entropy can be obtained as a special case of AFL entropy with a specific OPU choice (noise-OPU or coherent-state OPU); compare GK and AFL structures term-by-term; write new Section 54.
- [x] Numerically study OPU-choice independence: compute h_AFL^time for projector OPU vs SIC-POVM OPU vs matrix-unit OPU on the kicked Ising chain (L=5, n=4); quantify spread Δh across OPU choices; write new Section 55.
- [x] Implement opu_comparison.py: h_AFL^time for three OPU types on KI and XX chains; tabulate the spread; confirm whether h is OPU-independent at dual-unitary point.
- [x] Add Sections 54–55 to Output.tex: GK–AFL unification attempt (Sec 54), OPU independence study (Sec 55); reply to Mateo's notes.

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
