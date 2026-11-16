# Historical research context — 2026-11-15

Preserved during reorganization on 2026-11-16. Read [current context](../RESEARCH_CONTEXT.md) first. Pending-computation statements below describe this earlier handoff: later final-only observational-entropy runs are indexed in the current notes.

Historical paths below: scripts moved from `calculations/` to `calculations/code/`; generated files moved to `calculations/results/`, with unchanged filenames.

## Thesis
Working title: **Quantum Kolmogorov–Sinai Entropy**. Research on quantum dynamical entropy, particularly the Alicki–Lindblad–Fannes (ALF/AFL) framework and comparison with Connes–Narnhofer–Thirring (CNT), and links to quantum chaos and putative quantum Lyapunov exponents.

## Main topics
- Quantum partitions of unity, POVMs, quantum instruments, positive cones of measurement history effects and their generating properties.
- ALF history correlation matrices and entropy rates; CNT versus classical Kolmogorov–Sinai entropy.
- Observational entropy and its behavior under histories/refinement; von Neumann entropy lower bound.
- Measured quantum relative entropy and Helstrom state discrimination (single-copy versus many-copy distinctions).
- Quantum kicked top and SYK model; unitary closed systems and open/noisy counterparts.
- OTOCs, operator growth, finite-time diagnostics, semiclassical Lyapunov exponents and spectral chaos.

## Notation to verify per source
- `ρ` and `σ`: density operators; `d`: Hilbert-space dimension where finite.
- `X_i`: Kraus/partition operators satisfying Σ_i X_i† X_i = I.
- `E_i = X_i† X_i`: POVM effects.
- `Π_{i^n}`: effects from time-ordered instrument histories. **Check the history convention and order in each derivation.**
- Generating-cone questions concern approximation of density operators by positive combinations of history effects; specify the topology, time horizon and quantifiers.
- `S_obs(ρ) = -Σ_i Tr(ρ E_i) log[Tr(ρ E_i)/Tr(E_i)]` for finite-dimensional POVMs, omitting zero effects.

## Scientific standards
Prior numerical observations and exploratory conjectures are **not proven theorems**. Re-verify specific formulas and parameter values from actual scripts/files before citing them. State when a result applies only to a model, a POVM, a state, or a finite time.

## Conversation handoff — 2026-11-15

### User preferences
- Keep explanations human-readable, concise, and limited to what was asked. Do not add unsolicited derivations or inflate the LaTeX.
- The user explicitly requested shortening the relative-entropy document. Preserve that scope; do not restore the removed epsilon-expansion/SLD material unless asked.
- Current standalone deliverable: `latex/fixed_povm_relative_entropy_derivative.tex`. No LaTeX compiler is available; it is **NOT compile-tested**.

### Relative entropy: settled conventions and results
- The variable is the **initial rho in both arguments**, varied as rho+hA with Hermitian traceless A. Fix V=U_epsilon, U_t, and the measurement. This is NOT differentiation in epsilon or time, nor equal perturbations of the two final states.
- States: U_t rho U_t† and U_t V rho V† U_t†. Final directions: U_t A U_t† and U_t V A V† U_t†. POVM symbols M_i/E_i denote effects, not Kraus operators.
- [Proved] For fixed effects, probabilities p_i,q_i and their directional derivatives a_i,b_i give first differential sum[a_i log(p_i/q_i)-(p_i/q_i)b_i] and Hessian sum p_i(a_i/p_i-b_i/q_i)^2. Dropping +1 requires Tr A=0. Assume finite dimension and positive probabilities; full rank suffices.
- [Proved] Supremum over all single-copy POVMs is invariant under common unitary evolution **as a function of rho**; its first and second initial-state derivatives are also time-independent. The envelope theorem evaluates the first differential at an optimizer; it does NOT interchange derivative and supremum. The optimized Hessian includes optimizer response. Exact implicit formulas are in the current LaTeX and `notes/theorems/fixed-povm-relative-entropy-derivative.md`.
- [Proved] Data processing bounds every fixed final-measurement divergence and every common instrument's full-record divergence by D(rho||V rho V†). Fixed final-measurement values can oscillate: there is generally no classical channel mapping one experiment's final outcome distribution to the next time's. A full nested record's KL divergence is nondecreasing. These are different protocols.
- [Numerically evaluated] For the simulation's initial states, the quantum upper bound is **0.04460889 nats**, independent of time and kappa. This bounds the divergence, **not automatically its Frechet derivatives**. Full-rank assumptions matter; distinct pure states can have infinite optimized relative entropy.

### Kicked-top computations actually completed
- Scripts: `calculations/kicked_top_frechet.py` and `calculations/kicked_top_growth.py`; detailed parameters/results: `notes/calculations/kicked-top-frechet.md`.
- Spin j=10, d=21; U=exp[-i kappa Jz^2/(2j)] exp[-i(pi/2)Jy]; kappa=0.5,6; V=exp[-i0.07Jy]. Initial noisy coherent state: 70% coherent projector (theta=1.1,phi=0.7), 30% I/21. One random traceless Hermitian direction, operator norm 1, seed 20261112.
- **Only unitary evolution followed by a final projective Jz measurement was simulated. No repeated measurements or backaction were included.** Times are separate fresh-state experiments.
- [Numerically tested] Derivatives agree with central finite differences through 200 kicks. At h=3e-5, maximum absolute errors: first 1.89e-8, Hessian 8.10e-7. State/unitarity and epsilon=0 checks pass.
- Semilog and early-time plots: `calculations/kicked_top_growth_semilog.png`, `calculations/kicked_top_growth_early.png`; fit table: `calculations/kicked_top_growth_fits.csv`. No convincing exponential growth in tested windows. For kappa=6, Hessian fit over 0–20 kicks has rate 0.00770 but log-space R²=0.0144; do not report this as a Lyapunov exponent. Original 30-kick outputs were preserved.

### Discussed but NOT computed
- Repeated rank-one Jz projections: after the first recorded outcome the conditional state is the same for both preparations, so subsequent outcomes add no history KL distinguishability. No repeated-measurement script was completed.
- Proposed weak alternative: E_m=eta|m><m|+(1-eta)I/d, Kraus operators sqrt(E_m), suggested eta=0.2. This instrument generally preserves preparation dependence after conditioning. **It has not been simulated.**
- User requested observational-entropy plots for kappa=4–7, saturation, and initial slopes for kappa>6. This remains pending: protocol clarification was dismissed. Full-history observational entropy decreases under refinement; final-only or repeated/discarded-history observational entropy is different and can grow. Do not silently choose among them or claim those plots exist.

### K-property / ALF thread
- [Proved] Tracial binary shifts admit X=(I/sqrt(2),u_0/sqrt(2)) with R_n=2^{-n}I and rate log 2, including zero-CNT bitstreams. Exact supremum over strictly local partitions is log 2; only a lower bound was established for unrestricted von Neumann-algebra partitions.
- Source: downloaded Størmer survey, §9 pp.31–33; see `notes/papers/binary-shift-alf.md`. The user's book p.219, Tuyls comparison theorem, and Narnhofer–Thirring source hypotheses were **not independently checked**.
- [Conjecture here] Algebraic K-property implies positive ALF in a nontrivial faithful tracial system. No general proof or literature-open-status claim established.
- [Proved sufficient condition] For N_k=alpha^k(M_0), a unitary u in N_0 with E_{N_-1}(u)=0 yields orthogonal binary histories and ALF rate log 2. The gap is deriving such a partition, or another uniform entropy lower bound, from K alone. Pointwise tail triviality is not uniform control of exponentially many histories.
- Readable proof programme: `latex/k_property_positive_alf_foundations.tex`; status/gaps: `notes/open-problems/k-property-positive-alf.md`. This and the binary-shift LaTeX were left unchanged when the user shortened the relative-entropy document.
