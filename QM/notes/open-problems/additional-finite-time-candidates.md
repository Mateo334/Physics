# Additional finite-time state-function candidates

2026-11-25. Scope: closed finite-dimensional quantum systems, finite observation times, state functions and their consistent initial-state variations. Search of `latex/Matejs_notes.tex`, current notes and paper annotations found no developed stellar-root, transport-metric or Krylov-tangent programme. Reduced states have a precursor at notebook lines 1488–1492, but not spectral tangent dynamics. These are candidate directions, not solutions to a quantum golden braid.

## Stellar roots

[From source, institutional summary only] P. Bruno, *Quantum Geometric Phase in Majorana’s Stellar Representation: Mapping onto a Many-Body Aharonov-Bohm Phase*, Phys. Rev. Lett. 108, 240402 (2012). ESRF summary: https://www.esrf.fr/news/spotlight/spotlight163/index_html (read 2026-11-25). It describes pure spin-j states as 2j points and energy, symplectic structure and evolution in these coordinates. Primary paper not checked here.

[Proved, elementary implicit differentiation] For a spin-state polynomial P and a simple finite root z_k, its velocity is -partial_t P(z_k)/P'(z_k); its variation is -delta P(z_k)/P'(z_k). Overall polynomial scale/phase cancels, so the constellation is a function of the pure density matrix once the spin basis is fixed. Exact full-root dynamics is a projective-state reparametrization, not intrinsic stretching. Root matching, infinity and multiple roots need treatment; coherent spin states have coincident roots. Product-sphere distances are not the induced Fubini–Study metric. No mixed-state extension selected.

## Physical transport geometry

[From source] See `../papers/wang-wang-wu-physical-distance.md`. Probabilities on a fixed physically metrized basis or POVM outcome space are state functions. Transport distinguishes moving probability a long physical distance from redistributing it among nearby outcomes. It changes the geometry, not merely the displayed coordinates.

[Unresolved] This gives a finite-time separation measure, not automatically a closed Hamiltonian flow of distributions or a smooth tangent metric. Basis probabilities discard coherences and give a pseudometric on states. A spin-coherent Husimi distribution with a sphere ground metric is a separate candidate; do not attribute its details to the basis construction. Transport growth alone is not chaos.

## Reduced spectra

[Candidate] For a fixed bipartition, use eigenvalues of Tr_B rho and their variations under the same global unitary. Global closed-system assumptions remain satisfied; the marginal need not be autonomous. At a simple reduced eigenvalue p_a with eigenvector a, delta p_a=<a|Tr_B(delta rho)|a>. At degeneracy use spectral subspaces, not a chosen eigenvector derivative. Global correlations/Schmidt frames are needed for an exact lifted dynamics.

This measures entanglement/information redistribution rather than full-state distinguishability growth. Near-zero eigenvalues make -log p_a ill-conditioned. No entropy-production/Pesin identity established.

## Krylov direction: lower priority

[Literature lead, abstract and bibliography only] Parker, Cao, Avdoshkin, Scaffidi and Altman, *A Universal Operator Growth Hypothesis*, Phys. Rev. X 9, 041017 (2019), https://arxiv.org/abs/1812.08657, DOI 10.1103/PhysRevX.9.041017.

Krylov/Lanczos coordinates organize commutator evolution by operator depth; weights measure spreading. Complete amplitudes in a fixed orthonormal operator basis still evolve isometrically in Hilbert–Schmidt norm. A seed-dependent basis must not silently be rebuilt for each perturbation. Finite chains terminate. Abstract's large-system growth claims have not been checked for the user's finite setting. This is not presently a better answer to the requested tangent dynamics.

## Priority

Investigate stellar roots for an exact geometric base-plus-variation description, transport for a physically motivated non-unitarily-invariant distance, and reduced spectra for subsystem sensitivity. None currently supplies an established quantum Hamiltonian/tangent mechanism with intrinsic positive stretching. No simulations or full candidate derivations performed.
