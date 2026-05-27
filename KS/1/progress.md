# Progress Log — Subfolder 1: AFL Entropy

## Session: 2026-05-27

### Step 1 — Literature survey COMPLETE
- Alicki & Fannes (1994): quantum dynamical entropy via OPUs.
- Benatti book Ch. 8.2: AFL entropy construction, OPU density matrices.
- Fannes (1973) continuity inequality.
- Alicki & Fannes (2004) conditional entropy continuity.
- Winter (2016) tight uniform continuity bound.

### Step 2 — Output.tex written INCREMENTALLY (all sections done)

Sections written:
1. Introduction (motivation, classical KS vs quantum generalizations)
2. Mathematical Preliminaries (C*-algebras, states, von Neumann entropy, CP maps)
3. OPUs (definition, composition, refinement, Lemma: composition is OPU)
4. Density matrices from OPUs (positivity proof, consistency/Lemma under marginalisation)
5. AFL dynamical entropy (definition, properties: superadditivity, subadditivity, Fekete)
6. Reduction to classical KS entropy (Theorem 6.1, full proof)
7. AFL entropy for quantum spin chains (Theorem 7.1: h = s(ω) + log d)
8. Fannes inequality (1973) — FULL PROOF via variational worst-case argument
9. Alicki-Fannes (2004) continuity inequality — derivation via flag state + SSA
10. Winter (2016) tight bound
11. Properties: finite-level quantum systems have AFL = 0
12. Numerical verification (Python/NumPy, all tests pass)
13. NEW RESULT: state-dependent quadratic bound near maximally mixed state
14. Entropy increment monotonicity via SSA (Prop. + proof)
15. Connection to quantum Pesin theorem (outlook)

### Step 3 — Numerical verification COMPLETE (all pass)

- Fannes inequality: 0 violations in 5000 trials for d=2,4,8.
- Worst-case example achieves ratio exactly 1.000000 (tight).
- AFL entropy of finite system (d=3 shift): S/n → 0 confirmed.
- SSA monotonicity: non-increasing increments confirmed.
- Quadratic estimate near maximally mixed: matches actual ΔS to <5% for ε≤0.1.

### Step 4 — NEW RESULT proved

**Proposition (state-dependent sharpening)**
For ρ* = (1/d)I and any σ with T(ρ*,σ) = ε:
  2ε² ≤ S(ρ*) - S(σ) ≤ 2dε²

- Lower bound: Pinsker inequality + identity S(σ||ρ*) = log d - S(σ)
- Upper bound: Tr(δ²) ≤ ||δ||₁² for traceless Hermitian δ
- Sharper than Fannes for ε < (log(d-1))/(2d)

### Current state
Previous session complete. New session starting 2026-05-27.

### Step 5 — Literature survey for extensions
- Read Goldfriend-Kurchan (2021) "Quantum KS entropy and Pesin relation": 
  uses entropy production under weak noise coupling; derives quantum Pesin relation
  h_KS^{(2)} = sum of positive eigenvalues of OTOC matrix A-hat.
- Read Slomczynski-Zyczkowski (1994): CS-quantum entropy, correspondence principle.
- Plan: add Sections 16-19 to Output.tex covering GK entropy, OTOC, AFL-GK comparison,
  new rigorous results.

### Step 6 — New sections being written (in progress)
Plan for new sections:
- Sec 16: GK quantum KS entropy (setup, main formula, additivity)
- Sec 17: OTOC matrix, quantum Lyapunov exponents, GK-Pesin relation
- Sec 18: AFL vs GK comparison (table, semiclassical equivalence, Ehrenfest obstruction)
- Sec 19: New result — shift automorphism OTOC = 0 despite h_AFL > 0; implications for quantum Pesin
