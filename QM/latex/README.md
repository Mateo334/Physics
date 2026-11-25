# LaTeX retrieval guide

Do not read every TeX file to recover research context. Start with `../notes/INDEX.md` and the topic note; it identifies the relevant deliverable and assumptions.

## Original manuscript
- `Matejs_notes.tex`: user's primary notebook, unchanged. It contains two concatenated documents; consult `../notes/matejs-notes-research-map.md` for section/line retrieval. Do not repair or overwrite it without a specific request.

## Lyapunov work: current status matters
- `observable_sensitivity_closed_dynamics.tex`: exact lifted readout dynamics, bounds and counterexamples; **not** a solved golden-braid construction.
- `readout_dynamics_closure.tex`: autonomous linear-readout closure theorem.
- `readout_jacobian_method.tex`: earlier method sketch. Read with `../notes/open-problems/quantum-lyapunov-programme.md`; retaining rho does not require autonomous readout closure, and response singular values are not automatically classical Lyapunov exponents.
- `quantum_hamiltonian_tangent_geometry.tex`: exact geometry and limits; user did not accept it as solving the target.

## Other deliverables
- Current OE optimum: `all_outcome_obsen_global_maximum.tex`; earlier binary proof retained for provenance.
- Relative-entropy derivative scope: `fixed_povm_relative_entropy_derivative.tex` deliberately shortened at user request.
- `current_results_summary.tex`: dated snapshot, not a replacement for current research notes.
- Cogno introduction: `soledad_gonzalo_cogno_intro.tex`; sources in `../papers/soledad_gonzalo_cogno/README.md`.

No pdflatex or tectonic found on 2026-11-25. Do not infer compilation success from a `.tex` file's existence. Preserve all original filenames for existing links and citations.
