# Research context

## Start here
Thesis: **Quantum Kolmogorov–Sinai Entropy**. Main threads: ALF/AFL versus CNT and classical KS entropy; history effects and generating cones; observational entropy; measured relative entropy; kicked-top/SYK chaos diagnostics.

- [Index](INDEX.md): current theorem, literature, calculation and open-problem notes.
- [Calculations guide](../calculations/README.md): commands, code and saved results.
- [Kicked-top summary](calculations/kicked-top-simple-povms-lyapunov.md): latest final-only entropy, random-PVM and classical-exponent computations.
- [Historical handoff](archive/research-context-2026-11-15.md): full earlier context, exact formulas and original pending requests. Its pending final-only entropy plots have since been computed; repeated weak-measurement histories remain uncomputed.

## Standards and conventions
- Distinguish proofs, numerical checks, conjectures and citations. Recheck scripts/source assumptions before citing; finite-model observations are not universal chaos results.
- Search local files recursively before downloading. Store verified findings in notes and index them; preserve papers and thesis files.
- Density operators: ρ, σ. Partition/Kraus operators satisfy Σ Xᵢ†Xᵢ=I; effects Eᵢ=Xᵢ†Xᵢ are not Kraus operators. Check time ordering separately for each history convention.
- Generating-cone claims must specify topology, horizon and quantifiers. Observational entropy uses effect trace volumes and natural logarithms, omitting zero effects.
- Keep final-only measurements, full recorded histories and discarded-history protocols distinct. Full-history refinement and final-only entropy growth are different questions.
- Be concise; keep LaTeX within requested scope. Do not restore removed epsilon-expansion/SLD material to the shortened relative-entropy document without a request. `latex/fixed_povm_relative_entropy_derivative.tex` was not compile-tested.

## Settled results and remaining gaps
- [Proved] Relative-entropy derivatives here vary the initial rho in both preparations, holding evolution, perturbation and measurement fixed; they are not time/epsilon derivatives. See [theorem note](theorems/fixed-povm-relative-entropy-derivative.md).
- [Proved] Optimization over all single-copy POVMs is invariant under common unitary evolution as a function of rho, including its initial-state derivatives. The envelope theorem does not interchange derivative and supremum; the optimized Hessian includes optimizer response. Positivity/full-rank assumptions matter.
- [Proved] Data processing bounds fixed-measurement and common-instrument record divergences by initial quantum relative entropy, not automatically their derivatives. Nested full-record KL is nondecreasing; separate final-measurement values can oscillate.
- [Numerically tested] The j=10 derivative experiment agrees with finite differences through 200 kicks, without convincing exponential growth. [Details](calculations/kicked-top-frechet.md).
- [Not computed] Repeated weak Jz instrument with suggested eta=0.2. Rank-one recorded Jz projections erase preparation dependence of the conditional state after the first outcome; later outcomes then add no history KL distinguishability.
- [Proved] Tracial binary shifts have an explicit ALF partition of rate log 2, including zero-CNT cases. Strictly local partition supremum is log 2; unrestricted partitions have only the established lower bound. [Source and scope](papers/binary-shift-alf.md).
- [Conjecture here] K-property implies positive ALF in a nontrivial faithful tracial system. The conditional-expectation sufficient condition is proved; deriving it or another uniform history-entropy bound from K alone remains a gap. General literature-open status and the other historical source hypotheses are unverified. [Proof programme](open-problems/k-property-positive-alf.md).
