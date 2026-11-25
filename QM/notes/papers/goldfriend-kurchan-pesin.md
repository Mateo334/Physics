# Goldfriend–Kurchan: bath entropy, determinants, and the scope of the Pesin claim

Source: Tomer Goldfriend and Jorge Kurchan, *Quantum Kolmogorov-Sinai entropy and Pesin relation*, Phys. Rev. Research **3**, 023234 (2021), DOI 10.1103/PhysRevResearch.3.023234. Local PDF: `../../papers/PhysRevResearch.3.023234-2.pdf`.

Verified 2026-11-22 using local pypdf text extraction: main text pp. 1–9 and Appendix B, p. 12; Appendix A text also inspected. PDF page n equals printed page 023234-n. No web sources used.

## What actually produces entropy

- **Pp. 2–3, Eqs. (7)–(12):** classical white Gaussian noise, correlation scale τ shorter than system timescales, produces a Gaussian cloud around an unperturbed trajectory. Its covariance integrates products of tangent propagators. Noise stops contraction while unstable directions expand; the logarithmic volume growth selects positive classical exponents.
- **P. 4, Sec. III, Eqs. (18)–(19):** quantum evolution is coupled to a weak auxiliary bath. Infinite bath temperature permits neglecting friction and replacing the bath by classical fluctuating Hamiltonian fields coupled to physical observables. The density matrix whose purity is calculated is averaged over noise. Each individual noise realization is still unitary.
- **P. 5, Eq. (20); p. 12, Appendix B1:** expand the logarithm of the two-noise overlap to second order in noise, then perform the Gaussian noise average. The response kernel is the positive Hilbert–Schmidt Gram kernel

```latex
K_{\alpha\beta}(s,u)=
\frac{\operatorname{Tr}
 \{i[x^\alpha(s),\rho_0]i[x^\beta(u),\rho_0]\}}
{\hbar^2\operatorname{Tr}\rho_0^2}.
```

Thus this is a weak-noise approximation, not an exact arbitrary-coupling entropy formula.

- **P. 5, Eqs. (22)–(29):** for a pure initial state, center the Heisenberg observables and form time-integrated outer products of the vectors xᵅ_c(s)|φ₀⟩, including their complex-conjugate sector. The determinant lemma converts the noise-history determinant to a 2L-dimensional determinant (L is Hilbert-space dimension):

```latex
\Delta S_2(t)
=\frac12\log\det
\left(I+\frac{\epsilon^2\tau}{\hbar^2}A(t)\right),
\qquad S_2=-\log\operatorname{Tr}\rho^2.
```

Here A is the positive integrated Gram matrix in Eq. (29). It is **not** just a single OTOC, the classical tangent matrix, or a physical-time derivative of the state. Their finite-time rate is ΔS₂/t.

## What the quantum Pesin step establishes

**Pp. 7–8, Sec. V, Eqs. (42)–(45):** split A into early and late positive contributions. Dropping the early part yields a determinant lower bound.

**P. 8, Eqs. (46)–(48):** approximate the late interval, of duration comparable to an autocorrelation time, by its endpoint. The determinant lemma reduces this bound to a 2N-dimensional matrix (N physical degrees of freedom):

```latex
\begin{aligned}
C_{\alpha\beta}(t)
=\langle\{x_c^\alpha(t),x_c^\beta(t)\}\rangle,\\
t\bar h_{\rm KS}^{(2)}
\simeq\frac12\log\det(I+cC(t)),\\
c=\frac{\bar\tau\epsilon^2\tau}{\hbar^2}.
\end{aligned}
```

This endpoint approximation has assumptions beyond the preceding determinant inequality. The paper connects C to commutators with the initial projector ρ₀, hence an **OTOC-like** object. It explicitly distinguishes its construction from Ref. [12]'s coordinate-commutator Lyapunov-spectrum definition.

For eigenvalues of C scaling as exp(2λᵢt), the one-half log determinant gives the sum of positive λᵢ once the growing contributions dominate the identity. If eigenvalues are instead defined to grow as exp(γᵢt), the contribution is γᵢ/2. Keep this normalization explicit.

The paragraph after Eq. (48) explicitly says the authors have **not shown** that the individual rates deserve to be called Lyapunov exponents. They argue that exponential expansion can saturate the bound and check the semiclassical case, Eq. (49). Do not present a universal equality for arbitrary finite quantum systems as a proved theorem.

### Formula caveat

The PDF text contains internally inconsistent sign/factor expressions: Eq. (21)'s bare commutator trace has the opposite sign to the positive covariance; the second forms of (47)–(48) also should not be copied as a positive covariance identity. For a pure state, the correct elementary identity is

```latex
-\operatorname{Tr}([X,\rho_0][Y,\rho_0])
=\langle\{X_c,Y_c\}\rangle.
```

The text after (48) states eigenvalues exp(λᵢt) and then a rate sum without the one-half present in the determinant. Equations (40)–(41) likewise identify a positive entropy expression with log purity instead of minus log purity. Use Eq. (28)'s purity/entropy convention and the positive covariance form above; the source has not been checked against an erratum.

## Limits and timescales are essential

- **P. 6–7, Sec. III B, Eqs. (32)–(34):** at fixed ħ and fixed observation time, taking noise strength ε to zero gives zero entropy production, not a noise-independent positive KS rate. Amplification must overcome the identity term before the Ehrenfest/scrambling cutoff.
- Their timescale estimate gives the condition below, in their λ₁ convention:

```latex
\frac{\epsilon^2\tau k}{\lambda_1\hbar^2}\gg1.
```

- **P. 7, end of Sec. IV, after Eq. (41):** take ħ → 0 **before** ε → 0 to recover the classical Lyapunov expansion. Reversing these limits fails. The initial wavepacket cell plays the role of classical initial width squared.
- **Pp. 6–7:** an alternative is a large-system hierarchy where k grows with the number of degrees of freedom while ordinary correlation times do not; then ε can vanish while ε²τk remains large. This is a conditional separation of scales, not an unrestricted t → ∞ statement at fixed finite Hilbert-space dimension.
- **P. 6 and p. 9, Sec. VII:** the choice to couple noise to a number of physical observables comparable to the degrees of freedom is physically motivated but not unique. Coupling to exponentially many arbitrary projectors can change the answer.
- **P. 9, Sec. VII:** colored noise can enter the kernel in Eq. (25); finite-temperature friction requires additional Schwinger–Keldysh structure. The present calculation is not a generic finite-temperature dissipative-bath theorem.

## Relevance to OE / QFI exploration — deductions, not claims of this paper

1. **Fixed-POVM observational entropy is different.** It uses final probabilities and effect trace volumes, rather than the purity of a noise-averaged state. A measurement/coarse-graining theorem would be needed to transfer this determinant rate to OE. Neither the determinant lower bound nor its saturation establishes that transfer.
2. **There is a multiparameter QFI connection.** For a pure state with local parameter generators xᵅ(t)/ħ, the symmetric covariance matrix satisfies F = 2C/ħ². This follows directly from the pure-state SLD-QFI formula. It is a displacement/perturbation-parameter metric, not automatically physical-time QFI. The full A additionally integrates independent perturbations throughout the noise history.
3. **Physical-time QFI is a scalar in one direction.** For an autonomous Hamiltonian and a pure initial state it is 4 Var(H)/ħ², constant under that Hamiltonian's evolution. It cannot generally represent the growing volume determinant or the sum over unstable directions. Driven systems require their actual time-local generator.
4. **A defensible exploratory bridge:** compare an explicitly specified OE protocol with the bath-induced entropy and the positive covariance spectrum over a justified pre-Ehrenfest interval; document all noise, resolution, dimension, and initial-state limits. Any identification of the rates remains an additional conjecture.
5. **P. 9, Sec. VI:** the paper proposes coherent-state initial-condition averaging and a spin-j construction using two independent operators per site, with j → ∞ as semiclassical limit. Useful for kicked-top exploration, but not proof that a fixed Jy QFI or Jz-POVM OE has its KS rate.

No numerical reproduction or independent validation of the perturbative remainder or saturation argument was performed.
