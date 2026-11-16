# Algebraic K-property and positive ALF entropy: proof foundations

## Scope
[Conjecture] For a nontrivial finite von Neumann algebra with faithful normal trace tau and a trace-preserving automorphism alpha satisfying the user's algebraic K-filtration conditions, the supremum of ALF rates over all finite operational partitions in M is positive. This note does not claim novelty or verified open status in the literature. The unrestricted nontracial question is not settled here. A restricted partition algebra must be specified separately.

## Proven groundwork
[Proved] Write N_k=alpha^k(M_0), an increasing filtration. Trace-preserving conditional expectations onto N_{-r} converge in L2 to a -> tau(a)I, by the decreasing-projection/reverse-martingale argument and the trivial tail. This is pointwise and does not establish uniform control of growing history families.

[Proved: sufficient condition] If a unitary u in N_0 satisfies E_{N_{-1}}(u)=0, then X=(I/sqrt(2),u/sqrt(2)) has R_n=2^{-n}I and ALF rate log 2. Proof: in the inner product of distinct binary words, cancel the common latest factors. The latest remaining differing factor is alpha^k(u) or its adjoint, flanked by elements of N_{k-1}. Bimodularity and zero conditional expectation kill its trace. This does not require commutation or self-adjointness of u. The partition class must include X.

[Proved: sufficient estimate] Exponentially small purity Tr(R_n^2)<=C exp(-cn) for a single fixed partition implies S(R_n)>=cn-log C. No such purity estimate has been derived from the K-property here.

## Exact gap and cautions
[Unresolved] Does the K-property alone produce a partition with linear ALF entropy growth? The zero-conditional-expectation unitary is additional structure, not a consequence established here. Strict inclusion alone does not imply its existence: C subset C^2 with trace weights (0.9,0.1) has no centered unitary. This inclusion is NOT a counterexample to the K-conjecture.

- Exclude M=CI, which otherwise trivially satisfies the displayed algebraic conditions and has zero entropy.
- A single positive finite-time entropy is insufficient; one partition must work for arbitrarily long histories.
- Outcome Shannon entropy is an upper bound, not a lower bound, for correlation-matrix entropy.
- Pointwise tail decorrelation does not control exponentially many history overlaps.
- ALF>=CNT, even if available asymptotically, yields no positive lower bound at CNT=0.
- No claim about Lyapunov exponents follows here.

## Source connection and deliverable
Binary-shift source check and exact calculation: [existing note](../papers/binary-shift-alf.md). No new source-specific K theorem was verified. The algebraic K-conditions are those posed by the user; the Narnhofer–Thirring source has not been checked here for additional hypotheses.

Readable standalone draft: `latex/k_property_positive_alf_foundations.tex`. No pdflatex or tectonic available on PATH at the compilation check; NOT compile-tested.
