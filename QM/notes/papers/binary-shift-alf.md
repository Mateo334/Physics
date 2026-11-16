# Binary shifts: exact ALF calculation and source check

## Bibliography / local PDF and URL
- Erling Størmer, *A survey of noncommutative dynamical entropy* (2000), arXiv:math/0007010.
- https://arxiv.org/abs/math/0007010
- Local PDF: `papers/Stormer_2000_Survey_noncommutative_dynamical_entropy.pdf`.
- Checked printed/PDF pp. 31–35 (especially §9, pp. 31–33), using pypdf extraction and the arXiv HTML.
- Neshveyev–Størmer, *Dynamical Entropy in Operator Algebras* (2006), DOI 10.1007/3-540-34673-2: bibliography confirmed on Neshveyev's institutional personal site; book p. 219 NOT accessed. Tuyls's comparison paper NOT verified.

## Source summary relevant to this question
[From source] Section 9 defines bilateral generators u_j, self-adjoint unitaries with bitstream-controlled commutation/anticommutation. The canonical trace vanishes on every nonempty ordered monomial on distinct sites; alpha(u_j)=u_{j+1}. Equations (9.1)–(9.6) describe local algebra dimensions and mean entropies. These mean entropies are not ALF entropy.

[From source] Theorem 9.3, p. 33, gives a sufficient anticommutation condition for zero CNT entropy and asserts existence of bitstreams satisfying it. Lemma 9.1 and Corollary 9.4 give CNT entropy log 2 for the tensor square, despite zero CNT entropy of the factors. These zero-CNT results do not apply to every binary shift.

## Exact ALF result
[Proved] For the canonical trace, X=(I/sqrt(2),u_0/sqrt(2)) is an operational partition. Length-n histories are scaled distinct binary monomials. Trace orthogonality makes the entire ALF correlation matrix 2^{-n}I, not merely its diagonal. Hence S_n=n log 2 and the fixed-partition rate is log 2 for every bitstream in this setup. The projective partition ((I+u_0)/2,(I-u_0)/2) has the same spectrum by Hadamard mixing of outcomes.

[Proved] For the supremum over strictly local partitions, h_ALF^loc=log 2. A partition supported on L sites has all n-step histories in a vector space of dimension at most 2^{n+L-1}; the Gram rank bounds entropy by (n+L-1)log 2. Thus the lower bound is sharp for this partition class.

[Unresolved here] Equality for arbitrary nonlocal partitions was not established. The unrestricted von Neumann-algebra ALF supremum is at least log 2. The general K-implies-positive-ALF implication is neither proved nor disproved here. A one-step comparison cannot by itself control asymptotic rates; even full ALF >= CNT gives no strict positivity when CNT=0.

## Deliverable / limitations
- Derivation: `latex/binary_shift_alf.tex` (standalone).
- No pdflatex or tectonic found on PATH; NOT compile-tested.
- No numerics needed: exact Gram-matrix argument.
- Exact correspondence with the user's book p. 219 remains to be checked against that page.
