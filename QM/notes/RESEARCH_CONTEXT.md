# Current research context

## Scope and priority
Thesis: **Quantum Kolmogorov–Sinai Entropy**. Current goal: a quantum analogue of classical Hamiltonian-plus-tangent dynamics (the “golden braid”), in a **closed finite-dimensional system over finite times**, connected to classical instability and ultimately entropy. Do not substitute a generic chaos diagnostic or repeatedly object using infinite-time limits.

Original notebook: `latex/Matejs_notes.tex`, containing two concatenated documents. Never overwrite it. [Notebook map](matejs-notes-research-map.md) preserves source lines; [decision ledger](open-problems/quantum-lyapunov-programme.md) preserves the latest discussion and rejected shortcuts.

## Current conclusions (2026-11-25)
- Exact rho and its initial-state tangent obey the same unitary conjugation; invariant state-space lengths do not grow. This is NOT a proof of no quantum chaos.
- Fixed-probability/readout tangents can grow by changing visibility. The proposed normalized response spectrum is **unvalidated**, overlaps existing notebook work, and is not an established golden braid.
- Keeping rho underneath makes readout nonclosure acceptable. Requiring autonomous readout equations is an additional goal, not a prerequisite for observing sensitivity.
- Neither nonlinearity nor noninvertibility supplies instability. CPTP evolution is linear; conditional measurement is generally nonlinear but changes the physical problem. TDVP/mean-field can give classical-style tangent equations, but generally approximates exact dynamics.
- Classical Liouville evolution is also linear. The remaining route is a specified physical preparation/detector geometry with **derivative-controlled finite-time classical correspondence**; none established here.

## Established results / boundaries
- [OE optimum](theorems/all-outcome-obsen-global-maximum.md): binary optimality for every finite outcome count >=2, pure-state instantaneous rate; not a chaos theorem.
- [Relative-entropy derivatives](theorems/fixed-povm-relative-entropy-derivative.md): vary initial rho in both preparations; optimized measurements restore unitary invariance. Value bounds do not automatically bound derivatives.
- [Binary shifts](papers/binary-shift-alf.md): explicit ALF partition with rate log 2, including zero-CNT cases; general [K ⇒ positive ALF](open-problems/k-property-positive-alf.md) remains unresolved here.
- Final-only readouts, discarded backaction and conditioned histories are different experiments. No probability-based Pesin identity established.

## Retrieval and operations
Start with [INDEX](INDEX.md), then only the relevant note. Full inventories: `notes/CATALOG.md`, `calculations/CATALOG.md`. Preserve original paths for papers, code, data and manuscripts. Read archives only for historical specifics.
No pdflatex/tectonic was found on 2026-11-25; generated TeX must not be called compile-tested. Do not restore deliberately removed material to `latex/fixed_povm_relative_entropy_derivative.tex`.
GitHub: `Mateo334/Physics`, branch `main`, workspace maps to `QM/`; no local `.git`. **Before any pull, read [sync checkpoints](workspace-operations.md).**
