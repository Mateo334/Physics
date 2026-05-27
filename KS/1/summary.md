# Summary — Subfolder 1: AFL (Alicki-Fannes) Dynamical Entropy

## What was attempted

Full rigorous study of the Alicki-Fannes dynamical entropy (AFL entropy),
including all definitions, proofs, numerical verification, and an original
extension.

## What was established

### Definitions and constructions (all fully rigorous)
- OPUs (Operational Partitions of Unity): the non-commutative analogue of
  measurable partitions, satisfying Σ Z_i†Z_i = 1.
- Time-refined OPUs: Z^(n) = Θ^(n-1)(Z) ∘ ··· ∘ Z.
- OPU density matrix ρ[Z] with entries ρ[Z]_{ij} = ω(Z_j† Z_i).
- AFL entropy rate h_ω^AFL(Θ,Z) = limsup (1/n) S(ρ[Z^(n)]).
- AFL dynamical entropy = supremum over OPUs.

### Key theorems (proved in full)
1. **Consistency**: ρ[Z^(n+1)] marginalises to ρ[Z^(n)] (proved using
   the OPU property and Θ-invariance of ω).
2. **Reduction to KS entropy** (Theorem 6.1): For classical systems,
   AFL = KS entropy.
3. **Quantum spin chains** (Theorem 7.1): h_ω^AFL(Θ_σ) = s(ω) + log d.
4. **Fannes inequality**: |S(ρ) - S(σ)| ≤ ε log(d-1) + h(ε), tight.
5. **Alicki-Fannes (2004)**: |S(A|B)_ρ - S(A|B)_σ| ≤ 4ε log d_A + 2h(2ε).
6. **Winter (2016)**: Coefficient 2 vs 4, using quantum optimal transport.
7. **Finite-level AFL = 0**: S(ρ[Z^(n)]) ≤ 2 log d for all n, so S/n → 0.
8. **SSA → monotone increments**: δh(Z,n) non-increasing (Proposition 14).

### New result (original)
**State-dependent sharpening near the maximally mixed state** (Proposition 13):
For ρ* = I/d and any σ with T(ρ*,σ) = ε:
  2ε² ≤ S(ρ*) - S(σ) ≤ 2dε²

- Quadratic (not linear) behavior near the maximally mixed state.
- Sharper than Fannes for ε < log(d-1)/(2d).
- Lower bound: Pinsker + the identity S(σ||I/d) = log d - S(σ).
- Upper bound: Taylor expansion of relative entropy + Tr(δ²) ≤ ||δ||₁².
- Confirmed numerically to <5% error for ε ≤ 0.1, d=4.

## Status of this approach

NOT exhausted. The subfolder-1 approach (direct study of AFL entropy
construction) has successfully established all foundational results.
Remaining extensions (OTOC, quantum Lyapunov, Pesin) belong in future
subfolders or in this subfolder if Mateo specifies a direction.

## Files
- Output.tex: Full LaTeX with all definitions, proofs, and numerical results.
- progress.md: Step-by-step progress log.
- summary.md: This file.
