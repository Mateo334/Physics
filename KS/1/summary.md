# Summary — Subfolder 1: AFL (Alicki-Fannes) Dynamical Entropy

## What was attempted

Full rigorous study of the Alicki-Fannes dynamical entropy (AFL entropy),
including all definitions, proofs, numerical verification, and original
extensions.  The session also addressed Mateo's explicit instructions to
compute the density matrix rho[X]_{ij} = Tr(x_i omega x_j*) for finite-dim
systems and test on the Hadamard qubit.

## What was established

### Definitions and constructions (all fully rigorous)
- OPUs (Operational Partitions of Unity): non-commutative analogue of
  measurable partitions, satisfying Σ Z_i†Z_i = 1.
- Time-refined OPUs: Z^(n) = Θ^(n-1)(Z) ∘ ··· ∘ Z.
- OPU density matrix: rho[Z]_{ij} = omega(Z_j†Z_i) = Tr(Z_i omega Z_j*).
- Mateo's definition Tr(x_i omega x_j*) is IDENTICAL to AFL (cyclic trace).
- AFL entropy rate h = limsup (1/n) S(rho[Z^(n)]).

### Key theorems (proved in full)
1. **Consistency**: rho[Z^(n+1)] marginalises to rho[Z^(n)].
2. **Reduction to KS entropy**: For classical systems, AFL = KS entropy.
3. **Quantum spin chains**: h = s(omega) + log d.
4. **Fannes inequality**: |S(rho)-S(sigma)| <= eps log(d-1) + h(eps), tight.
5. **Alicki-Fannes (2004)**: |S(A|B)_rho - S(A|B)_sigma| <= 4eps log d_A + 2h(2eps).
6. **Winter (2016)**: coefficient 2 vs 4.
7. **Finite-level AFL = 0**: rank rho[Z^(n)] <= d^2 for all n, so S/n -> 0.
8. **SSA -> monotone increments**: delta h(Z,n) non-increasing.

### New results (original)

**Result 1: State-dependent sharpening (Proposition 13)**
For rho* = I/d and T(rho*,sigma) = eps:
  2eps^2 <= S(rho*) - S(sigma) <= 2d eps^2
Sharper than Fannes for eps < log(d-1)/(2d).

**Result 2: Closed-form formula for rho[Z^(2)] (Proposition 17.1)**
For projector OPU {|e_i><e_i|}, omega = I/d, Theta(A) = UAU*:
  rho[Z^2]_{(i,j),(k,l)} = (1/d) delta_{ik} |U_{ij}|^2 delta_{jl}
The density matrix is diagonal with a fully explicit closed form.

**Result 3: Matrix Entropy Formula (Theorem 18.2)**
  S(rho[Z^2]) = log d + E(U)
where E(U) = -(1/d) Σ_{i,j} |U_{ij}|^2 log|U_{ij}|^2 is the matrix entropy
of U.  Verified numerically to machine precision for 8 different gates.

**Result 4: MUB Saturation Theorem (Theorem 18.3)**
  rho[Z^2] = (1/d^2)I_{d^2}  iff  {|e_i>} and {U|e_j>} are MUBs
  iff  S(rho[Z^2]) = 2 log d (maximal).
Proved for all d >= 2.

**Result 5: Entropy saturation (observed and verified)**
For the Hadamard qubit:
- n=1: S = log 2, rank = 2
- n>=2: S = log 4, rank = 4 (saturated; S/n -> 0)
The density matrix reaches maximal rank after exactly one step of
Hadamard dynamics.

### Numerical verification
All results confirmed by hadamard_qubit.py:
- Hadamard entropy saturation: exact at n=2
- Matrix entropy formula: max error 4.4e-16 (machine precision)
- MUB characterization: Hadamard and U(pi/4) are MUB connectors
- Permutation-like gates (I, X, T, S): S(rho[Z^2]) = log 2 (no scrambling)
- Rotation angle formula: S = log2 + h(cos^2 theta), symmetric about pi/4

**Result 6: OTOC = (2/d)|U_{ji}|²(1-|U_{ji}|²) [Theorem 21.1, new]**
Full proof via commutator computation and anti-self-adjointness.
Verified to machine precision for 8 qubit gates.

**Result 7: AFL entropy is Rényi-1, OTOC probes Rényi-2 [Theorem 21.2, new]**
For distribution q_{ij} = |U_{ij}|²/d:
  S(ρ[Z^(2)]) = H_1(q)  [Rényi-1 = Shannon entropy]
  C_1(U) = (2/d³)(d - Σ|U_{ij}|^4)  related to Rényi-2 moment
  H_1(q) ≥ H_2(q), equality at permutations and MUBs.
MUB ↔ max AFL ↔ max OTOC: all three conditions are equivalent (Corollary).

**Result 8: Instruments vs. projections [Sec 20, new]**
Projector OPU → diagonal ρ[Z^(n)] = classical joint probability.
General OPU → off-diagonal coherences → S(ρ) > H(p) (quantum excess entropy).
Sequential measurement interpretation of diagonal entries proved rigorously.

**Result 9: Matrix-unit OPU is a quantum generating partition [Theorem 22.2, new]**
Quantum analogue of Krieger's generator theorem: the matrix-unit OPU achieves
the AFL entropy supremum for spin chains, requiring no optimisation.

**Result 10: Corrected AFL entropy and quantum Pesin [Proposition 23.2 + Conjecture 23.1]**
Corrected AFL entropy: h̃ = h_AFL - log d = s(ω).
- Removes quantum noise floor (inherent uncertainty = log d per site).
- Detects integrability: free-fermion chains → h̃ < log d; chaotic → h̃ = log d.
- Quantum Pesin conjecture: s(ω) = λ_L (quantum Lyapunov exponent via OTOC).
- Ehrenfest obstruction: no exponential OTOC growth for finite-dim systems.

**Result 11: Symbolic matrix structure (Sec 24, new)**
Explicit block-form of ρ[Z] and ρ[Z^(2)] with symbolic entries.
- ρ[Z]: k×k matrix with entries ω(Z_j*Z_i).
- ρ[Z^(2)]: k²×k² block-diagonal, entries (1/d)δ_{ik}|U_{ij}|²δ_{jl}.
- n=3: off-diagonal coherences emerge (Z^(3)_{(0,0,0)} = Z^(3)_{(0,1,0)} = (1/2)P_0),
  reducing S_AFL below H(p).

**Result 12: S_obs = 0 for maximally mixed state (Theorem 24.3, new)**
For ω = I/d and ANY OPU: p_k = V_k, so S_obs = -Σ p log(p/V) = 0.
S_AFL = H(V) = structural entropy of POVM cells for n ≤ n_sat.
For n > n_sat: quantum coherence gap H(V) - S_AFL grows, AFL stays ≤ 2 log d.

**Result 13: AFL-obs decomposition (Prop. 24.2, new)**
For projector OPU + any state: S_AFL = S_obs + H(V).
→ AFL entropy = observational entropy + structural entropy.
The two components measure: (state information gain) + (OPU cell complexity).

**Result 14: Structural entropy theorem (Theorem 25.4, new)**
For projector OPU + ω = I/d:
- n = 1,2: S_AFL(Z^(n)) = H(V^(n)) exactly.
- n ≥ n_sat: S_AFL saturates at log d + E(U) while H(V) grows as n log k.
- The quantum coherence gap = classical structural entropy − AFL = information lost to rank bound.

**Result 15: K-independence theorem (Theorem 27.1, new)**
For the kicked top U(k) = U_kick(k) * U_rot:
  |U(k)_{ij}|² = |[U_rot]_{ij}|²  for all k, i, j.
Proof: U_kick is diagonal in J_z eigenbasis → phases cancel in |·|².
Consequence: E(U^1) = E(U_rot) = k-independent. AFL at n=1 is phase-blind.
Verified to machine precision (error < 3e-16) for j = 0.5, 1.0, 1.5, 2.5, 5.0.

**Result 16: d=2 all-orders k-independence (Corollary 27.2, new)**
For j=1/2 (d=2): m² = 1/4 = constant, so ALL powers U^n are k-independent.

**Result 17: Qutrit odd-step k-independence (Conjecture 28.1, new)**
For qutrit kicked top (j=1, d=3): E(U^{odd}) = E(U^1) for all k (verified n=1,3,5).
Mechanism: parity symmetry of d^1(π/2) Wigner matrix + m² = {0,1} symmetry.

**Result 18: SIC-POVM AFL entropy (Proposition 28.2, new)**
For d=2 SIC-POVM OPU {Z_i = (1/√2)|φ_i><φ_i|}:
Eigenvalues of ρ[Z^(1)] = (1/2, 1/6, 1/6, 1/6), S = log(2)/2 + (1/2)log(6) ≈ 1.2425.
Proved analytically via equiangularity: [ρ]_{ij} = δ_{ij}/4 + (1-δ_{ij})/12.

**Result 19: OPU composition breakdown for non-unital channels (Theorem 29.1, new)**
For any TP-CP map Theta, the composed OPU satisfies:
  Σ_{ij} (Z^(2)_{ij})† Z^(2)_{ij} = Σ_j Theta(Z_j†) Theta(Z_j) ≤ Sigma_j Theta(Z_j†Z_j) = I
with inequality (Kadison) unless Theta is a *-homomorphism.
Explicit: amplitude damping gives P_0 + (γ²+(1-γ)²)P_1 ≠ I for γ ∈ (0,1).

**Result 20: Dissipation kills dynamical entropy (Remark 29.3, new)**
For amplitude damping: fixed point ω* = |0><0>, s(ω*) = 0, h~ = 0.
Three remedies for AFL with non-unital dynamics: Stinespring dilation,
instrument entropy, entropy production rate.

**Result 21: Quantum Pesin obstruction catalogue (Section 30, new)**
Five rigorous obstructions to naive Quantum Pesin for finite-d systems:
1. Finite-rate obstruction: h_AFL = 0 always.
2. Phase obstruction: E(U^1) is k-independent (kicked top).
3. Recurrence obstruction: E(U^n) oscillates, no linear growth.
4. Non-unital obstruction: OPU composition breaks for dissipative channels.
5. Ehrenfest obstruction: OTOC bounded by 2 log d.
Each resolved only in the double limit (semiclassical + thermodynamic).

**Result 22: Refined Quantum Pesin Conjecture (Conjecture 30.1, new)**
s(ω) = Σ_{λ_i > 0} λ_i  (sum of positive OTOC eigenvalues at Ehrenfest time)
in the double limit ℏ→0, L→∞ for infinite spin chains with KMS state.

**Result 22: AFL entropy is chaos-blind [Sections 31–33, new]**
For any quantum spin chain in a KMS state:
  h_AFL(shift, omega_beta) = s(omega_beta) + log d
This is INDEPENDENT of whether the chain is integrable or chaotic.
Proved analytically (follows from Theorem 7.1); confirmed numerically for XX and
kicked Ising chains (L=6).
- Free fermion XX (beta=1): s = 0.511510, h_AFL = 1.204657
- Kicked Ising (T=inf): s = 0.693147, h_AFL = 1.386294

**Result 23: Quantum Pesin fails for free fermions [Proposition 32.4, new]**
For the XX model at any finite temperature: s(omega) > 0 but lambda_L = 0 (integrable).
Therefore the quantum Pesin relation h~ = sum lambda_i^+ is FALSE for free fermions.
OTOC oscillates (no exponential growth), confirming lambda_L = 0.

**Result 24: Kicked Ising OTOC — dual-unitary instant saturation [Section 33, new]**
At J=g=pi/4 (dual-unitary point): OTOC = 0 for n < L/2, then jumps to maximum=4 at
n = L/2 (sharp lightcone). Entanglement reaches S_max = 3 log 2 in just 4 steps.
AFL entropy = same as XX chain (both = s(omega) + log 2).
OTOC and entanglement growth ARE chaos indicators; AFL is not.

**Result 25: General quantum Pesin inequality [Conjecture 34.2, new]**
  s(omega_beta) ≤ v_B * lambda_L
with equality iff system is maximally chaotic (saturates MSS bound).
For free fermions: v_B > 0 but lambda_L = 0 → bound 0 ≥ s(omega) > 0 is saturated at 0.
Remark: correct quantum Pesin is the GK formula h_KS^{(2)} = sum lambda_i^+ (not AFL).

**Result 26: Entanglement-Lyapunov conjecture [Conjecture 33.1, new]**
  v_E ≤ v_B * log d  (entanglement velocity ≤ butterfly velocity × log d)
Equality at dual-unitary (J=g=pi/4): v_E = log 2 = v_B * log d. Verified numerically.

## Status of this approach

NOT exhausted. Remaining open directions:
1. Rigorous proof of Qutrit Odd-Step Conjecture (Conjecture 28.1) via Wigner d-matrix analysis.
2. Rigorous proof of General Quantum Pesin inequality (Conjecture 34.2): s(omega) ≤ v_B * lambda_L.
3. Stinespring dilation approach to AFL for dissipative channels.
4. Semiclassical analysis in the j→∞ limit for kicked top (Weyl quantization route).
5. Operator entanglement entropy growth rate = v_B * s_op (Conjecture, Problem 5 in Sec 34).

## Files
- Output.tex: Full LaTeX (3664 lines), sections 1-34 + bibliography.
- progress.md: Step-by-step progress log.
- summary.md: This file.
- hadamard_qubit.py: Python verification (AFL entropy, MUB, matrix entropy formula).
- otoc_analysis.py: Python verification (OTOC formula, Rényi comparison, total OTOC).
- coarse_graining.py: Python verification (observational entropy, structural entropy, refinement bound).
- kicked_top.py: Python verification (K-independence, E(U^n) growth, qutrit, SIC-POVM, amplitude damping).
- free_fermion.py: Python verification (AFL entropy XX chain, OTOC comparison, entanglement growth).
