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

## Status of this approach

NOT exhausted. Remaining open directions:
1. Rigorous proof of Quantum Pesin Conjecture (Conjecture 23.1) for specific models.
2. Extension to d > 2 (qutrit Hadamard, SIC-POVMs).
3. Non-unital CP maps: what happens when the dynamics is dissipative?
4. Direct comparison of s(ω) with known Lyapunov exponents in XXZ chain.

## Files
- Output.tex: Full LaTeX (~2114 lines), sections 1-23 + bibliography.
- progress.md: Step-by-step progress log.
- summary.md: This file.
- hadamard_qubit.py: Python verification (AFL entropy, MUB, matrix entropy formula).
- otoc_analysis.py: Python verification (OTOC formula, Rényi comparison, total OTOC).
