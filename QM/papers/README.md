# Paper annotations

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
- `PhysRevResearch.3.023234-2.pdf`: Goldfriend, Kurchan, *Quantum Kolmogorov-Sinai entropy and Pesin relation*, Physical Review Research 3, 023234 (2021), DOI 10.1103/PhysRevResearch.3.023234. First page: entropy production induced by a weak auxiliary bath, not our fixed-POVM OE protocol; requires separation of Ehrenfest and correlation times.
- `2603.14025v3.pdf`: Nichele, Benatti, *Quantum Dynamical Entropy and Dissipative Information Flows*, PDF dated 4 November 2026. First page: open-system ALF extension and environment-memory witnesses; not fixed-POVM OE.
- `entropy-23-00645.pdf`: Facchi, Gramegna, Konderak, *Entropy of Quantum States*, Entropy 23, 645 (2021), DOI 10.3390/e23060645. First page: algebraic entropy under selection rules; not a model-specific OE-growth prescription.
- `frechet_derivative.pdf`: Murk, Tan, Muller, Safranek, *Connecting Quantum Tomography and Quantum Retrodiction*. First page: Petz map as a log-likelihood gradient and tomography algorithm; not a model-specific OE-growth prescription.
- A pypdf text search of these six PDFs and `2212.01585v3.pdf` for SYK/Sachdev/Kitaev found no matches. This does not replace a broader SYK literature search, which was not done.

## Størmer (2000): A survey of noncommutative dynamical entropy
- PDF: `Stormer_2000_Survey_noncommutative_dynamical_entropy.pdf`
- Source: https://arxiv.org/abs/math/0007010
- Detailed annotation: `../notes/papers/binary-shift-alf.md`.
- Checked §9, printed pp. 31–33, plus pp. 34–35. Defines canonical tracial binary shifts and proves existence of zero-CNT examples (Theorem 9.3). Tensor-square CNT entropy is log 2 (Lemma 9.1 / Corollary 9.4). The survey is not the source of the ALF calculation: the linked note independently proves fixed-partition rate log 2 and exact strictly-local ALF supremum log 2.
