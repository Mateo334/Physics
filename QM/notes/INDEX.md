# Research knowledge index

This is a curated navigation index for reusable mathematical and literature findings. Link new notes here and keep entries short.

## Core research
- [Research context and concise conversation handoff](RESEARCH_CONTEXT.md) — updated 2026-11-15; conventions, completed computations, pending requests, proof gaps, and user scope preferences.

## Definitions and theorem notes
- [Classical chaos diagnostics versus quantum spectral tests](definitions/classical-chaos-diagnostics.md) — Lyapunov exponents, KS/topological entropy, mixing, Pesin assumptions, Bethe ansatz and level statistics; standalone LaTeX overview.
- [Fixed and optimized measured relative entropy differentials](theorems/fixed-povm-relative-entropy-derivative.md) — full-rank derivatives, unitary tangent, SLD Fisher response, and time invariance of the optimized quantity.

## Literature
- [Probing chaos and criticality with OE](papers/probing-chaos-criticality-oe.md) — 2026 numerical phase-space POVM/Lyapunov correspondence, transition diagnostics and finite-shot resolution tradeoff.
- [Observational entropy in the kicked top](papers/observational-entropy-kicked-top.md) — verified Jz coarse-grainings and rate definitions; exponential deficit relaxation is distinct from linear entropy growth.
- [Binary shifts: source check and exact ALF calculation](papers/binary-shift-alf.md) — Størmer §9; zero-CNT examples have an explicit ALF partition of rate log 2; exact local-partition supremum.

## Numerical calculations
- [Kicked-top Husimi phase space and parameter animations](calculations/kicked-top-husimi-k3-k6.md) — same initial coherent state after five kicks, two static densities and k-sweep animations, including a 50-fps 1080p version.
- [Four-model observational entropy, one graph](calculations/four-model-obsen-100-states.md) — SYK, kicked tops at k=2,4,5,6, Bose--Hubbard and GOE; 100 coherent preparations each, five kicks and normalized fixed-PVM observational entropy.
- [Kicked-top Chernoff observational entropy across POVM ranks](calculations/kicked-top-chernoff-observational-povm-ranks.md) — fixed Jx/Jz PVM blocks and an unsharp full-rank Jz POVM with discarded Lüders backaction; k-distinction beyond rank-one Jz.
- [Coupled kicked tops: OE versus positive classical exponents](calculations/coupled-kicked-tops-obsen-lyapunov.md) — product phase-space PGM, 100-state mean, two-top tangent-map QR exponents, and short pre-Ehrenfest slope comparison.
- [Kicked-top outcome relative entropy](calculations/kicked-top-relative-entropy.md) — matched 100-state fixed-POVM comparison with/without unrecorded Lüders backaction, plus exact full-record history comparison.
- [Simple kicked-top code and classical Lyapunov comparison](calculations/kicked-top-simple-povms-lyapunov.md) — paper/random-axis POVMs, 20-state and 20-PVM averages, tangent-map exponents, and window-sensitive short-time entropy slopes.
- [Kicked-top observational entropy](calculations/kicked-top-observational-entropy.md) — j=200, four fixed POVMs; single-state plots and 20-coherent-state averages for k=4,5,6,7 through 20 kicks.
- [SYK fixed-POVM observational entropy](calculations/syk-observational-entropy.md) — four plots: occupation and particle-number PVMs and their unsharp counterparts; N=16, one disorder realization.
- [Kicked-top Frechet check](calculations/kicked-top-frechet.md) — spin 10, fixed Jz measurement; finite-difference verification through 200 kicks, semilog plots and exponential fits show no convincing growth in tested windows.

## Open problems and conjectures
- [Intermediate linear observational-entropy growth](open-problems/observational-entropy-linear-growth.md) — conditional conjecture; commuting-POVM counterexample to a universal chaos implication.
- [K-property and positive ALF: proof foundations](open-problems/k-property-positive-alf.md) — tracial formulation, an exact conditional-expectation sufficient condition, and the remaining gap to the general implication.
