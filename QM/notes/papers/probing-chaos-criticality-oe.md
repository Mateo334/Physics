# Probing Chaos and Criticality with Observational Entropy and Finite-Resolution Measurements

## Bibliography and sources
J. Bharathi Kannan, Sreeram PG, M. S. Santhanam (2026). arXiv:2605.23585v1. Source read: https://arxiv.org/html/2605.23585v1 (Sections II–VI, Eqs. 6–17, Figs. 1–4). PDF downloaded from https://arxiv.org/pdf/2605.23585v1 to `papers/Kannan_Sreeram_Santhanam_2026_Probing_Chaos_Criticality_OE.pdf`. Summary checked against HTML; PDF page numbers not checked.

## Main results [From source: numerical evidence]
- Section III, Figs. 1–2: derivatives of averaged OE locate the standard-map global-transport threshold near K=0.9716 and the Aubry–André localization transition at modulation strength 2J. The former uses entropy curvature (body specifies second derivative with respect to log K), the latter first derivative with respect to modulation strength. Long-time and eigenstate/Floquet ensemble estimates agree approximately; resolution and finite-size dependence are examined. Excessively coarse bins lose accuracy. Body/captions differ on derivative convention and panel labels; consult figures/code before reproducing.
- Section IV, Eqs. 6–12 and Fig. 3: finite-time OE grows approximately linearly after a transient and before Ehrenfest time in standard and regularized singular kicked rotors. Slopes from phase-space coherent-state POVMs normalized by the frame-operator inverse square root (PGM) closely reproduce reference classical Lyapunov exponents. Momentum-only binning underestimates them in these tests. Averaging uses 100 coherent initial conditions, uniformly sampled initial position and initial momentum zero. Phase-space slopes stabilize around 128 cells in reported examples; not a universal threshold.
- Section V, Fig. 4: Naimark dilation supplies an in-principle measurement realization; multinomial shot-noise simulations demonstrate a finite-resolution/statistical-error tradeoff. This is not a reported laboratory experiment.

## Interpretation and limitations
- The entropy/Lyapunov correspondence is a model-dependent numerical result supported by a semiclassical cell-counting heuristic, not a universal theorem for all POVMs or an asymptotic ALF/CNT entropy identity.
- Protocol computes OE of separately evolved states, not an entropy of repeated measurement histories.
- Frame inverse requires spanning coherent states (otherwise completeness holds on the support only). Nonorthogonality alone does not prevent a normalized coherent-state POVM: continuous Husimi POVMs already have a resolution of identity. The normalization issue here concerns discretized frames.
- Finite-size derivative peaks do not by themselves prove singularities or a universal classification by derivative order. The introduction mentions a critical-exponent result, but no explicit supporting fit was identified in Sections III–VI; do not repeat it as verified.
- Relevant to workspace: a concrete positive numerical example of intermediate linear OE growth with a specifically chosen phase-space POVM, unlike an assertion for arbitrary fixed measurements. No numerical reproduction performed.
