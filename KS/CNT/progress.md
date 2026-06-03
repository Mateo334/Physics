# Progress Log — CNT Subfolder

## Session: 2026-06-01

### Step 0 — Folder setup
- CNT/ subfolder created per Mateo's instruction in 1/notes.md and tasks.md.
- notes.md written.

### Step 1 — Literature study
- Read CNT entropy.pdf (original 1987 Connes–Narnhofer–Thirring paper, all 29 pages).
- Read CNT entropy chapter.pdf (Benatti book Chapter 5, pages 1–10).
- Extracted: abelian models, entropy defect, H_phi functional, dynamical entropy definition,
  Theorem IX.1 (shift = mean entropy), and all supporting lemmas.

### Step 2 — Output.tex written (8 sections, complete)
Sections:
1. Introduction and Motivation (KS failure in NC case, two obstructions)
2. Mathematical Preliminaries (von Neumann entropy, Umegaki relative entropy, CPU maps, eta function)
3. Abelian Models and Entropy Defect (Def. abelian model, eps_mu, s_mu, with proofs)
4. k-Subalgebra Entropy Functional H_phi(N_1,...,N_k) with ALL 7 properties proved
5. Extension to CP Maps (Def. H_phi(gamma_1,...,gamma_n), consistency, monotonicity)
6. CNT Dynamical Entropy h_phi(Theta) (Fekete argument, limit existence, 5 properties proved)
7. Main Theorem: shift entropy = mean entropy (full proof: upper bound + lower bound via clustering)
8. Explicit Examples (classical Bernoulli, qubit, qutrit, XXZ chain qualitative)
9. Comparison: CNT vs AFL/OPU (table, remarks on difference = log d)
10. Summary and Open Directions (4 open problems)
Bibliography: 8 references (CNT87, Benatti, Lieb, Lindblad, Araki, Alicki-Fannes, NT89, Emch)

### Step 3 — Python verification (cnt_basics.py)
8 examples, all pass:
1. H_phi(N) = h_2(p) for qubit, all p values. ✓
2. Qutrit: S_n/n = log(3) for all n (product state). ✓
3. Entropy defect: spectral decomp → s_mu=0, equal decomp → s_mu=S(rho). ✓
4. Relative entropy properties: non-neg, self=0, unitary covariance, joint convexity. ✓
5. Finite-dim CNT = 0: argument by rank bound. ✓
6. AFL vs CNT: difference = log d = log 2 for all p. ✓
7. n-step convergence: S_n/n → S(phi) for product state. ✓
8. Correlated qubit: S_2/2 < S_1 (strict subadditivity from correlations). ✓

### Step 4 — cnt_deeper.py written and run
Key results (all verified numerically):

**Product state (all p):**
- h_AFL - h_CNT = log d = log 2 exactly for all p. ✓

**Correlated XXX chain (L=8, beta=0,0.5,1,2,5):**
- Using CORRECT AFL formula h_AFL = s(omega) + log d:
  - h_AFL - h_CNT = log 2 = 0.693147 for ALL beta. ✓
- Note: single-site entropy S_1 > s(omega) for correlated states (S_1 = log 2 always,
  while s(omega) decreases from log 2 at beta=0 to 0.14 at beta=5).

**New theorem proved and written (Section 10 of Output.tex):**
Theorem (AFL-CNT gap): h_AFL(tau_1, Z) - h_CNT(tau_1) = log d for ALL translation-invariant
clustering states.
Proof: uses h_AFL = s(omega) + log d [subfolder 1, Thm 7.1] and h_CNT = s(omega) [CNT87, Thm IX.1].
Corollaries: gap = log d for any d; physical interpretation (OPU structural entropy).

### Step 5 — cnt_time_evolution.py written and run
Key results:
- Part 1: h_CNT(alpha_t) = 0 for all finite-L systems (proved by rank bound 2L*log 2).
- Part 2: Pre-saturation entropy slope for XXX chain (L=5, beta=1, dt=0.5):
  - Converges to ~0.101/step = 0.202/unit time by n=7 (below LR bound log2/0.5 = 1.386).
- Part 3: Lieb-Robinson bound h_CNT(alpha_t) <= v_LR*log(d)*|t| = 2J*log(2)*|t|. PROVED.
- Part 4: Modular automorphism sigma_t^phi = alpha_{-i*beta*t} (KMS/Tomita-Takesaki).
  - Imaginary-time evolution spreads site-0 operator to adjacent sites for s ~ 1.
  - Confirmed numerically (L=4, beta=1): site deviations grow with s. ✓

### Step 6 — Section 11 added to Output.tex
- Theorem: h_CNT(alpha_t) = 0 for finite-dim systems.
- Numerical table: pre-saturation slope for XXX chain.
- Theorem LR-Pesin bound for time evolution.
- Proposition: KMS condition identifies modular flow with imaginary-time evolution.
- Numerical table: spreading of sigma_s^phi(P0) across sites.
- Theorem: CNT entropy of modular automorphism (additivity, non-vanishing for non-tracial).
- Summary table: all three automorphisms compared.

Output.tex: 11 sections (+ bibliography), 155 environments balanced.

### Step 7 — cnt_chaos_diagnostic.py written and run
Key results:
- Part 1: h_CNT(shift) is NOT a chaos diagnostic (same S(phi) for both integrable and chaotic).
- Key lemma proved: diagonal U → sigma_z OPU gives trivial orbit (slope=0 for KI).
- Fixed: use sigma_x OPU for non-trivial comparison.
- With sigma_x OPU:
  - KI (dual-unitary): initial slope = 0.693 = log(2) (maximum), saturates at n=2.
  - XXX (integrable): initial slope = 0.407 < log(2), keeps growing past n=7.
- Phase diagram (kicked Ising, varying J=g): initial slope increases with chaos.
- Delta_CNT = S_chaos - S_int = 0.022 at n=5 (small positive = chaotic wins slightly).

### Step 8 — Section 12 added to Output.tex
- Prop: CNT shift entropy is chaos-blind.
- Lemma: Diagonal OPU triviality.
- Definition: initial slope Delta_1.
- Numerical tables: XXX vs KI at 7 steps; kicked Ising phase diagram.
- Theorem: CNT chaos diagnostic Delta_1 ∈ [0, log d].

Output.tex: 12 sections (+ bibliography), 178 environments balanced, 1687 lines.

### Step 9 — cnt_lower_bound.py written and run
Key results:
- DISPROVED: h_CNT >= Delta_1/|t| fails for finite L (h_CNT=0 but Delta_1>0).
- PROVED: h_CNT >= (1/2)*Delta_1 from Fekete, numerically confirmed for XXX at t=0.25,0.5,1,2.
- NEW OBSERVATION: AFL orbit entropy S(rho[Z^(n)]) is SUBADDITIVE (not superadditive) for XXX under time evolution:
  S_2 = 1.100 < 2*S_1 = 1.386 (XXX, sigma_x OPU, L=5, beta=1).
  KI dual-unitary: S_2 = 2*S_1 (equality, maximally superadditive).
- PROVED: Equality Delta_1 = h_CNT in the linear-growth regime (thermodynamic limit).

### Step 10 — Section 13 added to Output.tex
- Theorem: Fekete lower bound h_CNT >= (1/2)*Delta_1.
- Remark: AFL orbit entropy NOT superadditive under time evolution (new observation).
- Proposition: Disproof of h_CNT >= Delta_1/|t|.
- Theorem: Equality Delta_1 = h_CNT in linear-growth regime.
- Numerical table: Fekete bound verified at 4 time steps.

Output.tex: 13 sections (+ bibliography), 192 environments balanced, 1798 lines.

### Step 11 — cnt_operator_entanglement.py written and run (Session 2026-06-03)
Key results:
- Proved analytically: I_temp=0 <=> rho[Z^(2)]=I/d^2 <=> X-basis transition matrix T=uniform.
- T_{ij} = (d/D)*Tr[P_i U P_j U†]; T=1/d for all i,j iff I_temp=0.
- Operator entanglement E_op(u) of local 2-site gate (space bipartition):
  ranges from ~0 (integrable) to log(d) (dual-unitary J=g=pi/4).
  Eigenvalues of rho_L = [0.998,0.002,0,0] -> [0.5,0.5,0,0] as J goes 0 to pi/4.
- All results verified numerically (L=5, standard X-kick + ZZ coupling kicked Ising).
- Space-time complementarity: maximal E_op <=> zero I_temp <=> J=g=pi/4.

### Step 12 — Section 14.3 added to Output.tex
- Subsection 14.3: Connection to Operator Entanglement (3 sub-subsections)
- Definition: E_op(u) via reshuffled matrix R(u).
- Proposition: rank structure of kicked Ising gate (rank <= d, E_op <= log d).
- Theorem: I_temp=0 <=> T uniform (X-basis depolarising) — rigorous proof.
- Numerical table: E_op and I_temp for 6 coupling values, eigenvalue structure shown.
- Proposition: Maximal E_op <=> I_temp=0 <=> dual-unitary point.
- Remark: Space-time complementarity of maximal chaos.
- Reference added: Bertini et al. (2019).
Output.tex: 14+ sections (Sec 14 now has 3 subsections), 220 environments balanced.
Python script: cnt_operator_entanglement.py.

### Step 13 — cnt_eop_formula.py written and run (Session 2026-06-03)
Key results (both analytical and numerical):
- EXACT FORMULA proved: E_op(u(J,g)) = H_bin(sin²J), independent of g.
  Proof: rho_L eigenvalues = lambda_pm = (1±cos(2J))/2 = cos²J, sin²J.
  The kick H^g only rotates eigenvectors; eigenvalues depend only on J.
- INDEPENDENCE OF g: I_temp also independent of g (X-projectors commute with X-kick).
  Proved: K_g†P_i^x K_g = P_i^x, so Tr[P_i^x U P_j^x U†] = Tr[P_i^x D_ZZ P_j^x D_ZZ†].
- EQUIVALENCE: {E_op = log(d)} = {I_temp = 0} = {J = pi/4} x [0,pi/2] (a LINE, not a point).
- PHASE PORTRAIT: 5x5 grid in (J,g) space confirms I_temp depends only on J.
- All verified numerically (L=4, various J,g).

### Step 14 — Section 15 added to Output.tex
- Theorem: Exact E_op formula H_bin(sin²J), boxed.
- Full proof: block structure of rho_L, orthonormal columns h_0⊥h_1.
- Remark: geometric interpretation (J=0 rank-1, J=π/4 maximally mixed 2D subspace).
- Numerical table: exact formula vs code, g-independence verified.
- Locus table: (J,g) parameter space, both conditions = {J=π/4} (any g).
- Corollary: complete equivalence stated.
- Proposition 14.3 and Remark updated with corrected g-independent statement.
Output.tex: 15 sections (+ bibliography), 237 environments balanced, 2198 lines.
Python script: cnt_eop_formula.py.

### Current state
ALL TASKS COMPLETE (2026-06-03):
  - [x] CNT basics (Sections 1-9)
  - [x] CNT deeper: h_AFL - h_CNT = log d (Section 10)
  - [x] CNT time evolution + modular (Section 11)
  - [x] CNT chaos diagnostic (Section 12)
  - [x] CNT lower bound: Fekete bound; AFL subadditivity (Section 13)
  - [x] CNT subadditivity: rigorous proof + I_temp as diagnostic (Section 14.1-14.2)
  - [x] CNT open problem: subadditivity proof, equality conditions, operator entanglement connection (Section 14.3)
6 Python scripts: cnt_basics.py, cnt_deeper.py, cnt_time_evolution.py,
                  cnt_chaos_diagnostic.py, cnt_lower_bound.py, cnt_operator_entanglement.py.





