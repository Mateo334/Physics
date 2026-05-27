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

### Step 6 — New sections COMPLETE

Sections added to Output.tex (now 1926 lines, 48 theorem environments):

**Sec 16 (GK quantum KS entropy):**
- Setup: H = H₀ + noise coupling
- Main GK formula: Tr{ρ²(t)} = det[1 + (ε²τ/ℏ²)A]^{-1/2} (with derivation outline)
- Additivity proved (block-diagonal A structure)
- Semiclassical limit → classical KS (citing GK Sec IV)

**Sec 17 (OTOC matrix and GK-Pesin relation):**
- OTOC definition: C_ω(A,B;n) = ω([A_n,B]†[A_n,B])
- Lemma: OTOC bounded in finite dimensions, λ_L ≤ 0
- OTOC matrix: Â_{αβ}(t) = ⟨φ₀|[x^α(t),ρ₀][x^β(t),ρ₀]|φ₀⟩
- Quantum Lyapunov exponents: λᵢ = (1/t) ln σᵢ(Â(t))
- Theorem (GK Pesin): h̄_KS^{(2)} = Σ_{λᵢ>0} λᵢ — FULL PROOF via Minkowski inequality
  + Matrix Determinant Lemma + exponential dominance argument
- Remark: Ehrenfest obstruction (valid only for t < t_E)

**Sec 18 (AFL vs GK comparison):**
- Comparison table (6 rows)
- Theorem (semiclassical equivalence): lim h_AFL = lim h_KS^{(2)} = h_KS^cl
- Corollary: AFL-Pesin in semiclassical limit
- Proposition (Ehrenfest obstruction): finite-dim → AFL=0, GK>0 only for t<t_E

**Sec 19 (NEW RESULT — shift OTOC = 0):**
- Theorem: shift automorphism → [Θ^n(A),B]=0 for n>k, so OTOC=0 for all large n
- Corollary: butterfly velocity v_B = -∞ for shift
- Theorem (AFL-Pesin BREAKDOWN): h_AFL = s(ω)+log d > 0, sum of positive λ_L = 0
  → quantum Pesin FAILS for shift automorphism
- Physical diagnosis: information production ≠ sensitivity (Bernoulli shift analogy)
- Proposition: necessary conditions for quantum Pesin to hold
- Conjecture: quantum Pesin for chaotic spin chains (in thermodynamic+semiclassical limit)

All environments balanced (48 theorem-like, 26 equations, etc).
Bibliography extended with 6 new entries (GK2021, MSS2016, SZ1994, Hosur2016, LR1972).

### Current state
Output.tex complete and structurally verified. Approach still not exhausted.
Next directions:
- CNT entropy comparison
- Numerical verification of OTOC growth and AFL for random unitary circuits
- Concrete model where quantum Pesin holds (SYK, random circuits)
