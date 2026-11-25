# Paper annotations

## Soledad Gonzalo Cogno collection
- Dedicated subfolder: [soledad_gonzalo_cogno/](soledad_gonzalo_cogno/README.md).
- Source-grounded summaries cover single-neuron dynamics, balanced visual networks, population codes, neural sequences and model uncertainty; reading/access status is recorded per paper.
- Physicist-oriented introduction: `../latex/soledad_gonzalo_cogno_intro.tex` (standalone; not compile-tested).

## Gharibyan–Hanada–Swingle–Tezuka (2019): Quantum Lyapunov Spectrum
- PDF: `Gharibyan_Hanada_Swingle_Tezuka_2019_Quantum_Lyapunov_spectrum_1809.01671v2.pdf`; https://arxiv.org/pdf/1809.01671v2; JHEP 04 (2019) 082, DOI 10.1007/JHEP04(2019)082.
- Checked pp. 3–14, especially p. 9 Eqs. (10)–(11): canonical commutator matrix and state expectation of its operator Gram matrix. This is not the Gram of mean commutators. Spin/fermion versions use algebra-specific probes.
- Finite-time basis dependence and loss of classical symplectic pairing are explicit. Their sum of positive rates is a proposed “KS entropy” definition, not a theorem identifying ALF/CNT or observational entropy. Numerical SYK/XXZ results do not establish a universal asymptotic spectrum.
- Source details and exact response-derivation conventions: `../notes/papers/quantum-tangent-response-sources.md`.

## Rozenbaum–Ganeshan–Galitski (2017): Lyapunov Exponent and Out-of-Time-Ordered Correlator's Growth Rate in a Chaotic System
- PDF: `Rozenbaum_Ganeshan_Galitski_2017_Lyapunov_exponent_OTOC_1609.01707v3.pdf`; https://arxiv.org/pdf/1609.01707v3; Physical Review Letters 118, 086801, DOI 10.1103/PhysRevLett.118.086801.
- Checked pp. 1–4, especially p. 2 Eqs. (2)–(3): kicked-rotor squared commutator becomes an averaged squared classical derivative semiclassically. Log of an averaged square differs from the average log; the OTOC growth rate need not equal the usual trajectory Lyapunov exponent.
- Pre-Ehrenfest agreement, later slowdown, and domination by small chaotic islands are model-specific numerical findings. Details: `../notes/papers/quantum-tangent-response-sources.md`.

## Ashtekar–Schilling (1997): Geometrical Formulation of Quantum Mechanics
- PDF: `Ashtekar_Schilling_1997_Geometrical_Formulation_Quantum_Mechanics.pdf`; https://arxiv.org/pdf/gr-qc/9706069v1.
- Read pp. 6–8, 11–14 and 33–34 using `pypdf` (fitz unavailable). Eqs. (2.1)–(2.5), (2.14), (2.16)–(2.21), Corollary 1: exact Schrödinger flow is Hamiltonian and Killing on projective Hilbert space. Thus common-unitary FS tangent lengths are preserved, unlike generic classical Hamiltonian stretching.
- Eqs. (4.7)–(4.10): classical motion of oscillator expectation coordinates does not imply invariance of arbitrary coherent-state sections; standard ground-state-generated oscillator coherent states are a special invariant section. Metric normalization and unbounded-generator domains require care.
- Verified locations, conventions and deductions: `../notes/papers/quantum-geometric-mechanics-sources.md`.

## Hackl et al. (2020): Geometry of variational methods: dynamics of closed quantum systems
- PDF: `Hackl_et_al_2020_Geometry_Variational_Methods_Closed_Quantum_Systems.pdf`; published SciPost Physics 9, 048, https://scipost.org/SciPostPhys.9.4.048/pdf; DOI 10.21468/SciPostPhys.9.4.048.
- Read published pp. 26–28 and 35–37 using `pypdf`. Eqs. (81)–(96), Propositions 4–6: projective variational action and tangent projections; equivalence of variational principles for Kähler families, with non-Kähler caveats.
- Eqs. (122)–(128): fixed-point Hessian linearization, symplectic but generally non-isometric projected flow, and essential derivative-of-projector contribution. P. 36 explicitly distinguishes exact unitary Kähler preservation from projected dynamics. Published equation numbers differ from the arXiv HTML rendition.
- Detailed source verification: `../notes/papers/quantum-geometric-mechanics-sources.md`.

## Deffner (2020): Quantum speed limits and the maximal rate of information production
- PDF: `Deffner_2020_Quantum_speed_limits_information_production_1910.06811v3.pdf`; arXiv:1910.06811v3, https://arxiv.org/pdf/1910.06811v3; Physical Review Research 2, 013161 (2020), DOI 10.1103/PhysRevResearch.2.013161.
- Checked pp. 1–7 via text extraction, especially pp. 2–5, Eqs. (6)–(16), (24)–(33). Fannes/trace-speed finite-time entropy-change bounds; separate unitary measurement-Shannon discussion. Its rate notation divides finite changes by a QSL time, not a pointwise time derivative.
- Connection: our finite-d pure-state OE theorem gives the sharp local bound |dot S_obs|<=(C_d/2)||dotrho||_1 and an integrated fixed-POVM speed limit. Rank-one projective OE is the paper's discrete Shannon entropy.
- Caution: Eq. (26) labels an Lp difference norm as Wasserstein distance; Eq. (28) reverses the genuine Wasserstein ordering. Do not import the unitary continuity argument as written. Full summary and limitations: `../notes/papers/deffner-information-production-speed-limit.md`.

## Kannan, Sreeram PG, Santhanam (2026): Probing Chaos and Criticality with Observational Entropy and Finite-Resolution Measurements
- PDF: `Kannan_Sreeram_Santhanam_2026_Probing_Chaos_Criticality_OE.pdf`; arXiv:2605.23585v1, https://arxiv.org/html/2605.23585v1.
- Read HTML Sections II–VI, Figs. 1–4; PDF page numbers not checked. Detailed summary: `../notes/papers/probing-chaos-criticality-oe.md`.
- Numerical findings: OE derivatives locate kicked-rotor and Aubry–André transition points; pre-Ehrenfest linear OE growth with PGM-normalized Husimi phase-space measurements closely matches classical Lyapunov rates in standard/singular kicked rotors. Momentum binning underestimates rates in their examples. Resolution thresholds and finite-shot tradeoffs are demonstrated; not a general entropy/Lyapunov theorem or laboratory implementation.
- PDF verification: printed pp. 5–7, Section IV, Eqs. (6)–(8), (12), and Fig. 3 explicitly study temporal OE growth rates, not only mean entropy. Page 6 specifies averaging over 100 initial coherent states (uniform q0, p0=0). The growth slope is extracted within the Ehrenfest window; ensemble averaging does not remove the time dependence. Page 7 also describes finite-shot slope fits in Fig. 4. These are finite-window numerical comparisons, not an asymptotic entropy-rate theorem.
- Section III separately checked on printed pp. 3–4: mean OE from long-time dynamics and eigenstate/Floquet ensembles is differentiated with respect to control parameters to locate transition points, not to estimate Lyapunov exponents. KR ensembles average quasi-momentum sectors; AA ensembles average phase offsets. The main text uses second derivative with respect to log K for KR and minus the first derivative with respect to potential strength lambda for AA. Lambda here is not a Lyapunov exponent. Captions and main text differ in some panel/derivative labels; do not silently equate their conventions.

## Sreeram PG, Modak, Aravinda: Witnessing quantum chaos using observational entropy
- PDF: `2212.01585v3.pdf`, arXiv:2212.01585v3 (7 June 2023).
- Detailed annotation: `../notes/papers/observational-entropy-kicked-top.md`; checked pp. 1–7.
- Additional checks: p. 8 Conclusions and p. 10 Ref. [28]. The positive-Lyapunov link is indirect through an OTOC comparison and prior classical literature, not a demonstrated equality of OE slope with a classical Lyapunov exponent. Eq. (15) has C(t) approximately exp(2 lambda_q t), requiring care with the factor of two in a log-OTOC slope.
- Jz projective coarse-grainings; dynamics Fig. 3 uses rank-2 cells over half of d=400 and rank-4 cells over the rest, averaging 100 coherent states. Small unequal volumes preserve dynamic range and a state-dependent Boltzmann term. Fig. 3(b) fits log(log(d)-OE), not a sustained affine entropy law. Fig. 4 instead studies the short-time OE slope at the third kick. Our nonprojective POVMs are not this paper's protocol.

## Other PDFs screened for the model/measurement question
- `1803.00665v2.pdf`: Safranek, Deutsch, Aguirre, *Quantum coarse-grained entropy and thermalization in closed systems*. First-page abstract describes observational entropy, multiple possibly noncommuting coarse-grainings, and thermodynamic relaxation. Model details not checked here.
- `PhysRevResearch.3.023234-2.pdf`: Goldfriend, Kurchan, *Quantum Kolmogorov-Sinai entropy and Pesin relation*, Physical Review Research 3, 023234 (2021), DOI 10.1103/PhysRevResearch.3.023234. Verified pp. 1–9 and Appendix B: weak-noise, infinite-temperature auxiliary-bath Rényi-2 production; Eqs. (20), (28)–(29) give a time-integrated commutator/covariance Gram determinant. Sec. V, pp. 7–8, Eqs. (42)–(49), reduces this to a late-time covariance/OTOC-like determinant lower bound; saturation and a positive-expansion-rate sum require additional exponential-growth assumptions. Sec. III B, pp. 6–7, Eqs. (32)–(34), requires a pre-Ehrenfest window and correlation/scrambling separation; take the semiclassical limit before zero noise (p. 7), or an appropriate large-system limit. This is neither physical-time QFI nor fixed-POVM OE; no theorem equating their rates follows. Sign/factor conventions in the printed formulas need care; use the positive covariance form and the one-half log determinant. Detailed verified source note: `../notes/papers/goldfriend-kurchan-pesin.md`.
- `2603.14025v3.pdf`: Nichele, Benatti, *Quantum Dynamical Entropy and Dissipative Information Flows*, PDF dated 4 November 2026. First page: open-system ALF extension and environment-memory witnesses; not fixed-POVM OE.
- `entropy-23-00645.pdf`: Facchi, Gramegna, Konderak, *Entropy of Quantum States*, Entropy 23, 645 (2021), DOI 10.3390/e23060645. First page: algebraic entropy under selection rules; not a model-specific OE-growth prescription.
- `frechet_derivative.pdf`: Murk, Tan, Muller, Safranek, *Connecting Quantum Tomography and Quantum Retrodiction*. First page: Petz map as a log-likelihood gradient and tomography algorithm; not a model-specific OE-growth prescription.
- A pypdf text search of these six PDFs and `2212.01585v3.pdf` for SYK/Sachdev/Kitaev found no matches. This does not replace a broader SYK literature search, which was not done.

## Størmer (2000): A survey of noncommutative dynamical entropy
- PDF: `Stormer_2000_Survey_noncommutative_dynamical_entropy.pdf`
- Source: https://arxiv.org/abs/math/0007010
- Detailed annotation: `../notes/papers/binary-shift-alf.md`.
- Checked §9, printed pp. 31–33, plus pp. 34–35. Defines canonical tracial binary shifts and proves existence of zero-CNT examples (Theorem 9.3). Tensor-square CNT entropy is log 2 (Lemma 9.1 / Corollary 9.4). The survey is not the source of the ALF calculation: the linked note independently proves fixed-partition rate log 2 and exact strictly-local ALF supremum log 2.
# Additional source checked 2026-11-25

- `Wang_Wang_Wu_2020_Physical_Distance_Quantum_Chaos.pdf`: Zhenduo Wang, Yijie Wang and Biao Wu, *Genuine Quantum Chaos and Physical Distance Between Quantum States*, arXiv:2007.06244v2 (2020), journal reference Phys. Rev. E 103, 042209 (2021), DOI 10.1103/PhysRevE.103.042209. Source https://arxiv.org/pdf/2007.06244v2. Checked pp. 1–5: Wasserstein distances between probabilities in a physically metrized basis, explicitly a pseudometric on states (§II); proposed Lyapunov expression (§IV A) explicitly interpreted only before wavepacket spreading, not a strict infinite-time limit. Detailed assessment: `../notes/papers/wang-wang-wu-physical-distance.md`. Numerical sections not independently verified.
