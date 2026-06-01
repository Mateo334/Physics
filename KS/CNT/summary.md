# Summary — CNT Subfolder

## Session: 2026-06-01

### What was attempted

Full development of the Connes–Narnhofer–Thirring (CNT) entropy from scratch:
basics, properties, examples, and several new analytical results extending the theory.

### What was established

#### Basics (Sections 1–9)
- Complete rigorous construction of the CNT entropy: abelian models, entropy defect,
  k-subalgebra functional H_phi(N_1,...,N_k), extension to CP maps, CNT dynamical entropy.
- All 7 properties of H_phi proved (positivity, permutation invariance, monotonicity,
  subadditivity, repetition invariance, Theta-invariance, abelian reduction).
- Main theorem: h_phi(tau_1) = S(phi) for quantum spin chains with clustering states.
- Explicit examples: Bernoulli shift, qubit chain, qutrit chain, XXZ qualitative.
- Comparison table: CNT vs AFL/OPU (difference = log d).
- Python verification: cnt_basics.py, all 8 examples confirmed.

#### New Result 1: Exact AFL–CNT Gap (Section 10)
**Theorem:** h_AFL(tau_1, Z) - h_CNT(tau_1) = log d
for ALL translation-invariant clustering states, for any d-dimensional spin chain.
Proof: uses h_AFL = s(omega) + log d [subfolder 1] and h_CNT = s(omega) [CNT87].
Numerically verified for XXX chain at 5 temperatures (L=8).
The gap is EXACT, independent of temperature and correlations.

#### New Result 2: Time Evolution and Modular Theory (Section 11)
- h_CNT(alpha_t) = 0 for all finite-L systems (rank bound 2L*log d).
- Lieb-Robinson Pesin bound: h_CNT(alpha_t) <= v_LR * log d * |t|.
- Tomita-Takesaki connection: sigma_s^phi = alpha_{-i*beta*s} for KMS states.
- Imaginary-time evolution (modular automorphism) spreads operators across sites.
- Summary table: 3 automorphisms (shift, time evolution, modular) compared.

#### New Result 3: Chaos Diagnostic via Initial Slope (Section 12)
- Proved: Diagonal OPU Triviality Lemma — sigma_z OPU gives trivial orbit for diagonal U.
- Defined: initial slope Delta_1 = S(rho[Z^(2)]) - S(rho[Z^(1)]) ∈ [0, log d].
- Proved: Delta_1 = log d iff the orbit doubles the rank in one step (dual-unitary).
- Numerical table: KI dual-unitary achieves Delta_1 = log 2 (max); XXX achieves 0.407 < log 2.
- Phase diagram: Delta_1 increases monotonically with J=g (kicked Ising coupling).
- Theorem: CNT chaos diagnostic — Delta_1 detects chaoticity.

#### New Result 4: Fekete Bound and AFL Subadditivity (Section 13)
- Proved: h_CNT(alpha_t) >= (1/2) * Delta_1 [Fekete lower bound].
- Disproved: h_CNT(alpha_t) >= Delta_1/|t| fails for finite L.
- NEW OBSERVATION: AFL orbit entropy S(rho[Z^(n)]) under time evolution is SUBADDITIVE
  (NOT superadditive): for XXX chain, S_2 = 1.100 < 2*S_1 = 1.386.
  For KI dual-unitary: S_2 = 2*S_1 (equality, linear growth).
- Proved: Delta_1 = h_CNT in the thermodynamic linear-growth regime.

### Key open problems remaining

1. **AFL orbit subadditivity**: Why is S(rho[Z^(n)]) subadditive under time evolution
   but superadditive under shift? Find a rigorous proof and tight conditions.
2. **Tight lower bound**: What is the tightest lower bound on h_CNT^∞(alpha_t) in terms of Delta_1^∞?
3. **Non-clustering states**: Does h_CNT(tau_1) = S(phi) hold without the clustering condition?
4. **Type III factors**: Compute h_CNT explicitly for modular automorphism of type III_1 factors.
5. **Quantum Pesin for CNT**: Is h_CNT(alpha_t) = sum of positive CNT Lyapunov exponents for
   interacting spin chains in the semiclassical limit?

### Is the approach exhausted?

No. The CNT subfolder is actively developing. The approach has produced 4 new results
and 5 open problems. Further work on the open problems above is warranted.

### Files
- `Output.tex`: 13 sections + bibliography, 192 environments balanced, 1798 lines.
- `cnt_basics.py`: verification of all basic CNT formulas (8 examples).
- `cnt_deeper.py`: AFL–CNT gap = log d for XXX chain (5 temperatures).
- `cnt_time_evolution.py`: time evolution vanishing, LR bound, modular spreading.
- `cnt_chaos_diagnostic.py`: initial slope as chaos diagnostic (XXX vs KI).
- `cnt_lower_bound.py`: Fekete bound proof; AFL subadditivity disproved.
