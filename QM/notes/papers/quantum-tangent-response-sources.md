# Quantum tangent response: checked sources and exact/semiclassical distinctions

Checked 2026-11-24. Targeted source check, not an exhaustive literature review. PDF page numbers below refer to the saved arXiv PDFs unless specified. Formula derivations in §2–§4 are independent deductions, not quotations.

## 1. Source ledger

### A. Gharibyan–Hanada–Swingle–Tezuka: the closest matrix construction

Hrant Gharibyan, Masanori Hanada, Brian Swingle, Masaki Tezuka, **Quantum Lyapunov Spectrum**, JHEP **04** (2019) 082, DOI [10.1007/JHEP04(2019)082](https://doi.org/10.1007/JHEP04(2019)082). [arXiv:1809.01671v2](https://arxiv.org/abs/1809.01671v2), revised 5 November 2018.

Local PDF: `../../papers/Gharibyan_Hanada_Swingle_Tezuka_2019_Quantum_Lyapunov_spectrum_1809.01671v2.pdf` (40 pages; downloaded this check).

Verified locations:

- **Pp. 4–5, Eqs. (1)–(3):** classical tangent matrix, its singular values, and the positive matrix M†M. Their finite-time exponent is log(sᵢ)/t, not log(sᵢ²)/t.
- **Pp. 7–8, §2.3:** finite-time coordinate/basis dependence; classical Hamiltonian ± pairing need not survive the quantum definition.
- **P. 9, Eqs. (10)–(11), footnotes 2–3:** bosonic operator matrix is i[ẑᵢ(t),Π̂ⱼ(0)], with Π̂ⱼ the canonical conjugate to ẑⱼ, in their units/sign convention. The ordinary numerical matrix is the state expectation of the operator Gram matrix, with a summed output index. Their star means Hilbert-space adjoint. They explicitly reject using only the mean commutator as a reliable growth diagnostic.
- **P. 10, Eq. (12):** fermionic SYK uses an **anticommutator** of Majoranas. Eq. (16) uses raising/lowering spin commutators for XXZ. These are choices adapted to the algebra, not literal canonical-coordinate derivatives.
- **P. 11, Eq. (17), §4.3:** the sum of positive quantum exponents is a **definition** of their “KS entropy”; they acknowledge other definitions and unclear properties away from the classical limit.
- **Pp. 11–12, §4.4:** finite-size perturbations can mix energy scales substantially, especially in their SYK calculations.
- **Pp. 12–14, §4.5, Eqs. (19)–(23):** relation to the largest OTOC growth rate is argued using ETH, a growth window, and dominance of the largest exponent. Eq. (20) expresses the trace as a sum of squared singular values. The Hilbert-space “Lyapunovian” and projected nonlinear-state dynamics are explicitly distinguished from their degrees-of-freedom-indexed matrix.
- **P. 4:** approximately exponential SYK growth and brief exponential then power-law XXZ growth are numerical findings. Entanglement-production comparison is numerical, not a general Pesin theorem.

**Use:** direct precedent for an operator-valued tangent surrogate and a positive, state-dependent response spectrum. Do not replace their Gram expectation by the Gram of expectations, or transfer their ± pairing/entropy language without qualifications.

### B. Rozenbaum–Ganeshan–Galitski: averaging changes the exponent

Efim B. Rozenbaum, Sriram Ganeshan, Victor Galitski, **Lyapunov Exponent and Out-of-Time-Ordered Correlator's Growth Rate in a Chaotic System**, Phys. Rev. Lett. **118**, 086801 (2017), DOI [10.1103/PhysRevLett.118.086801](https://doi.org/10.1103/PhysRevLett.118.086801). [arXiv:1609.01707v3](https://arxiv.org/abs/1609.01707v3), revised 18 January 2017; PDF header dated 20 January 2017.

Local PDF: `../../papers/Rozenbaum_Ganeshan_Galitski_2017_Lyapunov_exponent_OTOC_1609.01707v3.pdf` (7 pages; downloaded this check).

Verified **pp. 1–4**, especially:

- **P. 2, Eq. (2):** momentum–momentum squared commutator in the kicked rotor.
- **P. 2, Eq. (3):** operator derivative followed by a **semiclassical approximation** by the squared classical initial-position derivative, averaged over initial phase space. The latter is not asserted as an exact quantum/classical equality.
- **P. 2, discussion below Eq. (3):** log of averaged squared stretching differs from average of logarithmic stretching. Their correlator growth convention is exp(2 λ̃ t).
- **P. 3, Figs. 2–3 and fitting paragraph:** quantum and classical correlator growth agree in the studied pre-Ehrenfest window; at large K, λ̃ exceeds the usual Lyapunov exponent by approximately log(√2). Fits omit t=0 because this particular commutator vanishes there.
- **Pp. 3–4:** growth slows beyond the Ehrenfest time; small chaotic islands can dominate the correlator even without global chaos.

**Use:** a mandatory warning against identifying a squared-response average with the usual trajectory Lyapunov exponent, even after taking a valid semiclassical limit.

### C. Goldfriend–Kurchan: a different entropy-producing response geometry

Tomer Goldfriend, Jorge Kurchan, **Quantum Kolmogorov-Sinai entropy and Pesin relation**, Phys. Rev. Research **3**, 023234 (2021), DOI [10.1103/PhysRevResearch.3.023234](https://doi.org/10.1103/PhysRevResearch.3.023234).

Existing local PDF: `../../papers/PhysRevResearch.3.023234-2.pdf`. Existing detailed check: [goldfriend-kurchan-pesin.md](goldfriend-kurchan-pesin.md).

Freshly rechecked **printed/PDF pp. 7–9 and reference page 13**:

- **P. 7, Eqs. (35)–(37):** expand around a classical trajectory to a time-dependent quadratic fluctuation Hamiltonian. Its linear fluctuation evolution is the classical tangent propagator R. This is a **semiclassical linearization** of the original nonlinear problem; exact linear evolution applies to the quadratic fluctuation model.
- **P. 7, Eqs. (38)–(41):** localized oscillator-vacuum initial state and doubled one-boson-sector determinant recover a classical integrated tangent Gram matrix. The limit order is explicitly ħ → 0 before noise strength → 0.
- **P. 8, Eqs. (43)–(45):** positive early/late decomposition gives a determinant lower bound. **Eqs. (46)–(48)** additionally approximate the late interval by endpoint correlations.
- Their endpoint matrix is the centered **anticommutator covariance of evolved coordinates**, related to commutators with the initial density projector. It is **not** the coordinate-coordinate commutator Gram in source A.
- **P. 8, paragraph following Eq. (48):** explicitly “We have not shown” that their individual rates deserve the name Lyapunov exponents. Exponential expansion is used to argue possible saturation.
- **P. 8, Eq. (49):** semiclassical endpoint covariance reduces to R Rᵀ times the initial fluctuation scale.
- **P. 9:** coherent initial-state averaging, spin generalization, and arbitrariness in the set of bath-coupled observables.
- **P. 13, Refs. [12], [28]:** A above and D below are explicitly different spectrum constructions.

The pre-existing note checked pp. 4–6 and Appendix B: entropy is Rényi-2 entropy of a **weak-noise-averaged state**, with a time-integrated response determinant. Those pages were not freshly re-extracted here. Printed sign/factor inconsistencies in Eqs. (40)–(41), (47)–(49) are another reason to use the positive covariance and one-half log determinant in the existing note, not copy formulas mechanically. No erratum search completed.

### D. Rozenbaum–Ganeshan–Galitski: Hilbert-space Lyapunovian and Moyal discussion

Same three authors, **Universal level statistics of the out-of-time-ordered operator**, Phys. Rev. B **100**, 035112 (2019), DOI 10.1103/PhysRevB.100.035112; [arXiv:1801.10591v3](https://arxiv.org/html/1801.10591v3), 16 July 2019. Journal metadata corroborated by source C's Ref. [28].

**HTML checked; no local PDF, PDF pagination not checked.**

- **§I, Eq. (1):** logarithm of minus a squared coordinate commutator divided by 2t defines a Hilbert-space operator, the “Lyapunovian.” It is not a 2N-by-2N tangent-direction matrix.
- **§III, Eqs. (4)–(5):** spectral restrictions for logarithms are explicit for odd powers. In general a log of a positive semidefinite operator also requires handling its kernel; positivity alone does not give a finite logarithm.
- **§V, Eqs. (6)–(9):** Moyal evolution and an even-ħ expansion describe quantum trajectories and higher corrections; the phase-space level-correlation interpretation is explicitly **heuristic**.
- **§VI, Eqs. (10)–(12):** Gaussian wavepackets, finite pre-Ehrenfest windows, and phase-space averaged squared derivatives. Available times did not achieve complete self-averaging.

**Ordering caution (independent check):** the exact Weyl symbol of an operator square is a **star square**, not generally the pointwise square of its symbol. Do not read §V Eq. (6)'s square as a universal exact pointwise formula. Our explicit star-product formula below avoids that ambiguity.

### E. Robert: exact linear generators versus general Moyal corrections

Didier Robert, **When Poisson and Moyal Brackets are equal?**, [arXiv:2207.08798v2](https://arxiv.org/html/2207.08798v2), 2 March 2023 (first submitted 18 July 2022).

**HTML checked; no local PDF, PDF pagination not checked; journal publication not verified.**

- **§2.1, Remark 2.1, Eqs. (2.2)–(2.3):** commutator correspondence with a linear phase-space function is exact, via Weyl translations.
- **§2.2, Eqs. (2.4)–(2.9), Lemma 2.2:** exact oscillatory star product, formal ħ expansion, and symbol remainder estimates.
- **Theorem 1.1:** within the stated moderate-growth symbol class (δ < 1/2), equality of Poisson and Moyal brackets with H for **every Schwartz observable**, at fixed ħ, forces H to be at most quadratic. **Remark 1.2:** quadratic H suffices.
- This does **not** prohibit exact correspondence for a particular observable under a particular nonquadratic Hamiltonian.
- Robert uses the opposite Poisson-bracket sign from §2 below and i/ħ times the commutator. Translate conventions rather than copy signs.

**Egorov reference trail:** Robert Ref. [1] identifies A. Bouzouina and D. Robert, *Uniform semi-classical estimates for the propagation of quantum observables*, Duke Math. J. **111**, 223–252 (2002). Its original theorem/proof was not inspected. Search also located Wang Zuoqin's “Lecture 17: Egorov's theorem,” Theorem 2.1 (compactly supported Hamiltonian symbol, fixed-time O(ħ²)); only search-extracted text was inspected, so do not use this check as verification of an arbitrary Hamiltonian/Ehrenfest-time remainder theorem.

## 2. Exact response identity: what C actually differentiates

Use canonical Z=(Q₁,…,Qₙ,P₁,…,Pₙ), with J the standard symplectic matrix, and the convention {qⱼ,pₖ}=δⱼₖ:

```latex
\begin{aligned}
[Z_a,Z_b]&=i\hbar J_{ab}I,\\
J&=\begin{pmatrix}0&I\\-I&0\end{pmatrix},\\
C_{ab}(t)&=\frac{[Z_a(t),Z_b]}{i\hbar}.
\end{aligned}
```

For unbounded operators these identities are on a suitable common invariant domain, with states having the required moments; arbitrary products/derivatives need not exist on all Hilbert space. Canonical CCR cannot hold exactly in finite dimension (take a trace). Spins and finite torus quantizations require their own algebra and local/geometric interpretation.

An impulse generated by Z_b gives an **exact directional derivative**, with no small-ħ approximation:

```latex
\begin{aligned}
V_b(s)&=e^{-isZ_b/\hbar},\\
\left.\frac{d}{ds}
 V_b(s)^\dagger Z_a(t)V_b(s)\right|_{s=0}
 &=C_{ab}(t).
\end{aligned}
```

Here t is held fixed. This is a derivative with respect to the initial perturbation parameter s, not a physical-time derivative. Z_b generates a symplectically rotated displacement, not displacement along coordinate b itself.

For literal initial phase-space translation by ξ choose

```latex
\begin{aligned}
D(\xi)&=\exp[-i\xi^T JZ/\hbar],\\
D(\xi)^\dagger ZD(\xi)&=Z+\xi I,\\
T_{aj}(t)
&=\left.\partial_{\xi_j}
 D(\xi)^\dagger Z_a(t)D(\xi)\right|_{\xi=0}\\
&=\sum_b C_{ab}(t)J_{jb},\\
T(t)&=C(t)J^T=-C(t)J.
\end{aligned}
```

Thus **C(0)=J I, whereas T(0)=I**. This symplectic factor is necessary to call the object a derivative with respect to initial coordinate components. With dimensionally normalized canonical coordinates J is orthogonal, so the associated direction-space Gram spectra coincide up to that basis rotation. Arbitrary coordinate units/rescalings change the finite-time metric.

For displaced initial states the corresponding measurable mean response is exactly

```latex
\begin{aligned}
\rho_\xi&=D(\xi)\rho D(\xi)^\dagger,\\
\left.\partial_{\xi_j}
\operatorname{Tr}(\rho_\xi Z_a(t))\right|_0
&=\operatorname{Tr}(\rho T_{aj}(t)).
\end{aligned}
```

This is not itself the squared-response diagnostic. A positive numerical Gram matrix is

```latex
\begin{aligned}
L_{jk}(t)
 &=\sum_a\operatorname{Tr}
   [\rho T_{aj}(t)^\dagger T_{ak}(t)],\\
v^\dagger L v
 &=\sum_a\left\|
   \sum_jT_{aj}(t)v_j\sqrt{\rho}
   \right\|_{\rm HS}^2\geq0,\\
\lambda_j(t)&=\frac{1}{2t}\log\ell_j(t).
\end{aligned}
```

The last line is a **definition** where eigenvalues ℓⱼ are positive and the normalization is fixed; it does not prove existence of asymptotic Lyapunov exponents. Source A supplies a closely matching construction. L differs from the Gram of the mean response, which discards fluctuations. For real physical displacement vectors the real symmetric part of L gives the same quadratic form, but its spectrum need not equal the complex Hermitian Gram spectrum: specify the choice.

## 3. Weyl/Moyal versus classical tangent flow

Let zᵃ_t be the **exact Weyl symbol** of Z_a(t), not the classical trajectory. Since the other commutator argument z_b is linear, the following is exact in the Weyl calculus:

```latex
\begin{aligned}
\sigma_W(C_{ab}(t))
 &=\{z^a_t,z_b\}
 =\sum_c(\partial_c z^a_t)J_{cb},\\
\sigma_W(T_{aj}(t))&=\partial_j z^a_t.
\end{aligned}
```

This precise identity is stronger than a generic “commutator approximates Poisson bracket” slogan: **the derivative of the quantum symbol is exact; replacing that symbol by the classical flow is not**.

For a Weyl Hamiltonian with a suitable smooth, ħ-independent symbol h, the exact evolution is

```latex
\partial_t z^a_t
=\frac{z^a_t\star h-h\star z^a_t}{i\hbar}.
```

Its leading semiclassical evolution is classical Hamiltonian transport. Under appropriate symbol/flow hypotheses, fixed-time Egorov gives the leading replacement zᵃ_t ≈ Φᵃ_t. For Weyl quantization without subprincipal contributions, the usual formal first correction is O(ħ²); symbol derivatives and long-time remainders need their own control. An operator-norm estimate alone does not justify differentiating its remainder against an unbounded coordinate.

Consequently, only in a controlled semiclassical regime,

```latex
\begin{aligned}
\sigma_W(T_{aj}(t))
 &\simeq \partial_j\Phi_t^a,\\
C(t)&\longleftrightarrow M(t)J,
\qquad M=D\Phi_t.
\end{aligned}
```

For quadratic Hamiltonians this correspondence is exact under the usual canonical/Weyl assumptions. For nonlinear chaotic dynamics, derivative growth degrades the expansion, typically on a logarithmic-in-1/ħ Ehrenfest scale; its coefficient and validity depend on the model and observable. No universal long-time estimate is proved here.

Moreover, taking products introduces star products even when the single commutator identity is exact:

```latex
\sigma_W(C_{ab}^\dagger C_{cd})
=\overline{c_{ab}}\star c_{cd},
\qquad c_{ab}=\sigma_W(C_{ab}).
```

Taking a state expectation then averages with its Wigner distribution. Only further semiclassical localization/averaging assumptions identify this with the Gram matrix of one classical trajectory. Wigner distributions are not generally positive measures. Source B demonstrates that even a legitimate classical positive ensemble average changes the extracted exponent through the order of averaging and taking logarithms.

## 4. Noncommutative Fréchet derivative: retain left and right multiplication

An inner derivation is an exact Leibniz derivation:

```latex
\begin{aligned}
\delta_b A&=[A,Z_b]/(i\hbar),\\
\delta_b(AB)&=(\delta_b A)B+A(\delta_b B).
\end{aligned}
```

For a fixed ordered polynomial F in several operators, its derivative inserts the perturbation into **each position**:

```latex
DF(Z)[V]
=\sum_{r=1}^m
 Z_{i_1}\cdots Z_{i_{r-1}}
 V_{i_r}
 Z_{i_{r+1}}\cdots Z_{i_m},
\quad F(Z)=Z_{i_1}\cdots Z_{i_m}.
```

This is a linear map on an operator tuple, not generally a scalar matrix multiplying a vector. The elementary counterexample to the commuting chain rule is

```latex
D(X^3)[V]=X^2V+XVX+VX^2,
```

which becomes 3X²V only when the needed commutation holds. For an exponential (bounded operators, or suitable differentiable functional calculus),

```latex
D(e^X)[V]=\int_0^1 e^{(1-s)X}Ve^{sX}\,ds.
```

On polynomials the derivation chain rule is exact:

```latex
\delta_bF(Z)=DF(Z)[\delta_b Z].
```

If Heisenberg coordinates obey a closed ordered equation dot Z=F(Z), applying the fixed initial derivation gives the exact operator variational equation

```latex
\dot C_{\cdot b}(t)
=DF(Z(t))[C_{\cdot b}(t)].
```

“Closed” here permits specified ordered polynomials/functional calculus with the necessary domains. It does not mean a finite c-number system. For a quartic potential V(Q)=gQ⁴/4, for example,

```latex
\begin{aligned}
\dot Q&=P/m,\qquad \dot P=-gQ^3,\\
\dot C_{Qb}&=C_{Pb}/m,\\
\dot C_{Pb}
 &=-g(Q^2C_{Qb}+QC_{Qb}Q+C_{Qb}Q^2).
\end{aligned}
```

Replacing the last line by −3gQ²C_Qb is generally wrong. Evaluating a classical Hessian on noncommuting coordinates does not fix this ordering issue.

Finally distinguish this from a derivative of **the state evolution map**:

```latex
D(\rho\mapsto U_t\rho U_t^\dagger)[\Delta]
=U_t\Delta U_t^\dagger.
```

This exact linear map preserves unitarily invariant norms. It is not an exponentially unstable tangent map on density matrices. Response growth concerns selected observables/derivations and their metric, not increasing full-state distinguishability under a common unitary.

## 5. Search/check log and scope

### Independent check of the main-report conventions

The main report's generators and response are consistent with §2 (its D is our T, not our Weyl unitary D(ξ)):

```latex
\begin{aligned}
G_b&=\sum_c(J^{-1})_{cb}R_c
    =\sum_cJ_{bc}R_c,\\
D_{ab}(t)&=[R_a(t),G_b]/(i\hbar),\\
D_{ab}(0)&=\sum_cJ_{ac}(J^{-1})_{cb}I
            =\delta_{ab}I.
\end{aligned}
```

In particular G_q=P and G_p=−Q. Thus initial translations are generated by exp(−i ξᵇG_b/ħ), and their adjoint action shifts R by +ξ. Reversing either the unitary exponent or the commutator order without also reversing the other changes the sign.

For a **fixed scalar-coefficient noncommutative polynomial vector field** F, the report's exact equation dot D=DF(R(t))[D] is correct column by column: differentiate in time, commute F with the fixed initial G_b, and apply the ordered Leibniz rule. Scalar coefficients are material: noncentral operator coefficients would generate additional commutator terms. No replacement by a commuting Jacobian is justified. The quartic force derivative in §4 agrees exactly with the proposed main-report formula.

**Independent inverted-oscillator check (exact quadratic example, not evidence of global chaos):** in dimensionally normalized canonical variables choose H=κ(P²−Q²)/2. Then

```latex
\begin{aligned}
R(t)&=M(t)R(0),\\
M(t)&=\begin{pmatrix}
\cosh(\kappa t)&\sinh(\kappa t)\\
\sinh(\kappa t)&\cosh(\kappa t)
\end{pmatrix},\\
\Delta\mu(t)&=M(t)\xi,\\
\Sigma(t)&=M(t)\Sigma_0M(t)^T.
\end{aligned}
```

For a pure Gaussian with symmetric covariance Σ₀, initially displaced along ξ and then evolved under the same unitary, the Fubini–Study infinitesimal line element (one quarter of SLD QFI) is

```latex
\begin{aligned}
ds_{\rm FS}^2
 &=\frac14\,\xi^T\Sigma_0^{-1}\xi\\
 &=\frac14\,\Delta\mu(t)^T
       \Sigma(t)^{-1}\Delta\mu(t).
\end{aligned}
```

Here ξ is an infinitesimal displacement. Along an unstable direction its raw centroid difference grows as exp(κt), while the covariance grows in the same direction and the FS line element is exactly constant. The constant-distance statement also follows directly for finite separations from unitary invariance of state overlaps. This example supports observable-response growth alongside fixed state distinguishability; the inverted oscillator is an unbounded hyperbolic quadratic system, so it does not prove finite-system chaotic entropy production. Its covariance grows without a finite-dimensional saturation cutoff.

Source connection: A justifies using a Gram of operator responses rather than squared mean response; E supports exact linear-generator and quadratic-Hamiltonian correspondence; C's quadratic fluctuation construction explains why a covariance tangent geometry can emerge semiclassically, but its bath-averaged entropy is absent from this isolated oscillator example. B's average-versus-log warning applies when unstable rates vary across initial conditions; this uniform quadratic example does not test that distinction.

1. Read local `notes/INDEX.md`, `notes/RESEARCH_CONTEXT.md`, `papers/README.md`, and `notes/papers/goldfriend-kurchan-pesin.md`.
2. Recursively listed local `papers/` and searched local notes for `Rozenbaum`, `Egorov`, `Moyal`, `Lyapunov spectrum`, `tangent matrix`. No existing PDFs for A/B were present. Used installed pypdf 5.3.0; installed nothing.
3. Re-extracted source C pp. 7–9 and 13. An initial console-encoding failure was corrected using Python UTF-8 mode.
4. Web queries: “Rozenbaum Ganeshan Galitski Lyapunov exponent out of time ordered correlator 2017 arxiv”; “quantum Lyapunov spectrum matrix commutator Hanada Gharibyan Tezuka”.
5. Verified A/B metadata against arXiv abstract records. Downloaded **exactly two** PDFs from `https://arxiv.org/pdf/1609.01707v3` and `https://arxiv.org/pdf/1809.01671v2`, guarded against overwriting existing paths. Confirmed parseable PDFs and page counts 7 and 40.
6. Read A pp. 3–14 and B pp. 1–4 via local extraction; locations above delimit what supports the claims. Did not reproduce their numerical results.
7. Followed C's reference trail to D and read its HTML §§I, III, V–VI, particularly Eqs. (1), (4), (6)–(12).
8. Queried “Egorov theorem Weyl quantization Moyal bracket Ehrenfest time Bouzouina Robert 2002” and “Didier Robert When Poisson and Moyal brackets are equal arxiv 2022”; read E's HTML theorem and Weyl-calculus sections. The Bouzouina–Robert original remains an unverified primary-source lead.
9. Cross-check: García-Mata, Jalabert, Wisniacki, *Out-of-time-order correlators and quantum chaos*, [arXiv:2209.07965v1](https://arxiv.org/html/2209.07965v1), HTML §I.3 Eqs. (9)–(11), §II.1, and finite-torus Weyl-translation Eqs. (14)–(20). This review corroborates the semiclassical/averaging restrictions; no PDF downloaded or journal metadata checked.

**Bottom line:** C is an exact inner-derivation response matrix; −CJ is the canonical initial-coordinate displacement response. Its connection to a classical tangent matrix is semiclassical except in special exact cases. A state-weighted Gram supplies a defensible quantum response spectrum, but neither classical Pesin equality nor Goldfriend–Kurchan entropy production follows just by naming its eigenvalue rates Lyapunov exponents.

Not checked: general unbounded-operator domain theorems, rigorous long-time Egorov constants, source errata, or numerical reproduction. No universal identification with fixed-POVM observational entropy is established.
