# Observational entropy, Fisher speed, and a possible Pesin connection

Primary concise deliverable: `latex/observational_entropy_qfi_pesin_concise.tex`. Detailed derivations and caveats: `latex/observational_entropy_qfi_tangent_pesin.tex`. Neither is compile-tested: no pdflatex or tectonic was available. No new simulation.

[Proved] For fixed finite POVM, positive probabilities, fixed trace volumes, and differentiable state curve: dot S = -Cov_p(dot p/p, log(p/v)); hence |dot S| <= sqrt(classical time FI times Var(log(p/v))) <= sqrt(SLD time QFI times that variance). Applies to general observation coordinates, not only density matrices. Zero-probability extensions require regularity; smooth unitary pure trajectories suffice, arbitrary C1 density curves do not.

[Proved] All pure state--tangent pairs with equal Fubini–Study speed are unitarily equivalent. Therefore unrestricted m-outcome instantaneous supremum = C_d,m Delta H/hbar = C_d,m sqrt(time QFI)/2. This proves scaling, not an exact numerical C. Effects are frozen for differentiation; reoptimizing at each point does not give one entropy trajectory attaining all optima. Independently optimized ensemble mean uses mean sqrt(QFI), not sqrt(mean QFI); common-POVM optimum only obeys the resulting upper bound.

[Proved, subsequently extended to every finite m>=2] The [binary global maximum](binary-obsen-global-maximum.md) first established the exact coefficient for two outcomes. The new [all-outcome certificate](all-outcome-obsen-global-maximum.md) proves C_d,m=C_d,2 for every finite m>=2, resolving the previously recorded gap. A rank-one projector/complement pair is globally optimal; proportional splitting gives any larger finite outcome count. The global individual-effect inequality, not KKT stationarity alone, proves the upper bound.

[Proved] Time velocity and initial-condition perturbations are distinct. Classical Hamiltonian tangent expansion can be transverse to the neutral flow direction. Finite-dimensional quantum amplitudes form a Hamiltonian but also orthogonal/unitary flow; initial-label QFI is preserved under a common parameter-independent U. Autonomous time QFI is constant. Positive optimized instantaneous OE rate already occurs in a periodic qubit.

[Proved, illustrative classical model only] Diagonal linear hyperbolic flow, Gaussian preparation, and fixed Gaussian blur give observed covariance sigma² M M^T + epsilon² I. Differential entropy slope tends to sum of positive exponents; classical time FI tends to twice the sum of their squares. Without blur the entropy slope is zero. This unbounded model is not itself a Pesin theorem or a quantum OE theorem.

[From source] See `notes/papers/goldfriend-kurchan-pesin.md` for checked local paper and page/equation details. Its bath/noise-averaged Renyi entropy, response determinant and required limits differ from fixed-POVM OE. [Proved deduction] Its pure endpoint symmetric covariance is hbar²/2 times a multiparameter perturbation-QFI matrix. No scalar-time-QFI identity with a Lyapunov sum follows.

[Unresolved] Transfer a response-volume/coarse-graining mechanism to specified semiclassical POVM OE with controlled errors, resolution, ensemble and pre-saturation time window. Finite-d single-time OE is bounded and has zero asymptotic net rate. No general quantum Pesin theorem established.
