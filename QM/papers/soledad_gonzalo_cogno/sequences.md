# Gonzalo Cogno: slow sequences and their recurrent mechanism

Read/access checked 2026-11-24. This note separates experimental findings, model results, and analytical claims. All paths below are relative to the workspace root.

## Sources and access

1. **Soledad Gonzalo Cogno, Horst A. Obenhaus, Ane Lautrup, R. Irene Jacobsen, Claudia Clopath, Sebastian O. Andersson, Flavio Donato, May-Britt Moser, Edvard I. Moser.** “Minute-scale oscillatory sequences in medial entorhinal cortex.” *Nature* **625**, 338–344 (2024). Published online 20 December 2023; issue 11 January 2024.
   - DOI: https://doi.org/10.1038/s41586-023-06864-1
   - Full text: https://www.nature.com/articles/s41586-023-06864-1
   - Authorized open-access PDF: https://www.nature.com/articles/s41586-023-06864-1.pdf
   - Local: `papers/soledad_gonzalo_cogno/GonzaloCogno_et_al_2024_Nature_Minute_scale_oscillatory_sequences.pdf` (49 PDF pages; includes Methods, Extended Data and reporting material).
   - Data/code links reported by the paper: https://doi.org/10.25493/SKKX-4W3 and https://github.com/soledadgcogno/Ultraslow-oscillatory-sequences.git
2. **Lea Marie Braun, Mathias Karsrud Nordal, Silje-Neneva Hanssen Rambø, Claudia Clopath, Soledad Gonzalo Cogno.** “A common structure in recurrent networks supports neural sequence generation locally and in downstream neurons.” *bioRxiv*, version 1, **posted 8 September 2026** (not 2 September, despite the DOI identifier). Preprint, not certified by peer review.
   - DOI: https://doi.org/10.64898/2026.09.02.748960
   - Versioned record: https://www.biorxiv.org/content/10.64898/2026.09.02.748960v1
   - Authorized PDF: https://www.biorxiv.org/content/10.64898/2026.09.02.748960v1.full.pdf
   - Local: `papers/soledad_gonzalo_cogno/Braun_et_al_2026_Common_structure_recurrent_networks_bioRxiv_v1.pdf` (55 PDF pages).
   - License displayed: CC BY-NC-ND 4.0.
   - Analytical supplement, discovered via the supplementary-material page:
     https://www.biorxiv.org/content/biorxiv/early/2026/09/08/2026.09.02.748960/DC1/embed/media-1.pdf?download=true
   - Local: `papers/soledad_gonzalo_cogno/Braun_et_al_2026_Common_structure_analytical_supplement_v1.pdf` (11 pages).
   - Main p. 11: no new data collected; code promised upon publication. Code availability was not independently investigated.

Downloads returned HTTP 200 and were verified as PDFs with pypdf. Files were created exclusively, without overwriting. Local recursive PDF search preceded downloading. PyMuPDF/fitz was unavailable; pypdf was already installed. No package installation. Same-stem `.txt` files preserve page-labelled extracted text; extraction mangles some mathematical subscripts, so use the PDF/supplement for equations.

### Actual reading ledger

- **Nature:** read the publisher's complete main text, including Discussion and captions to Figs. 1–5. Main text occupies PDF pp. 1–7 / journal pp. 338–344. Direct PDF/text checks: pp. 12–16 (PCA, nonlinear embedding, phase, oscillation score, sequence identification, phase locking and mutual information); p. 19 (ensemble transition and sequence statistics); p. 21 (linear readout model, data/code); pp. 45–46 (Extended Data Fig. 12 and caption). Visually inspected Fig. 2 on PDF p. 3 and Extended Data Fig. 12 on p. 45 using extracted PDF pages passed to the PDF reader. This is not a claim to have inspected every Extended Data panel or all 49 pages.
- **Braun main:** read pp. 2–10 (Introduction, all Results, Discussion), metadata/abstract on publisher record, p. 11 availability, pp. 12–15 main-figure legends; pp. 21–31 (model, targets, preprocessing, FORCE, connectivity profile, sparse-training procedures and Gaussian fits); pp. 32–36 (perturbations and downstream readout). Visually inspected all five main figures on PDF pp. 41–45. Equation checks particularly pp. 21–23, 27, 35–36.
- **Braun analytical supplement:** read all pp. 1–11, including support definition, Fourier truncation, circulant averaging, degeneracy proof, idealized ridge calculation and downstream width/period bounds.
- No simulations or statistical analyses were rerun; no raw recording data or videos were inspected.

## A concise physicist's introduction

These papers concern **collective dynamics and an inverse problem**. Hundreds of neurons can have a complicated instantaneous activity vector, yet its slow evolution is approximately parametrized by one circular phase. The experimental paper finds repeated minute-scale traversal of this population-state loop. The modelling paper asks which features of recurrent coupling are required to generate such a trajectory, how much microscopic connectivity remains undetermined, and which temporal patterns a downstream linear observer can reconstruct.

The useful physics connection is a moving activity bump, Fourier modes of an approximately translation-invariant interaction in **sequence order**, and low-dimensional collective coordinates. It is not evidence of quantum dynamics or deterministic chaos. “Ring” refers to activity space/order, not to a physical ring of neurons.

## 1. What the Nature experiment actually establishes

### Protocol and evidence

Head-fixed mice ran freely on a non-motorized wheel in darkness, without scheduled rewards or changes in navigational position/direction. MEC activity was recorded by two-photon calcium imaging and independently with Neuropixels spiking recordings.

- Single-cell autocorrelations and their spectra show ultraslow oscillations, below 0.1 Hz. In the 15 oscillatory imaging sessions, 5,691/6,231 recorded cell observations (91%) met the single-cell oscillation criterion. Neuropixels: 683/879 units (78%). These are not claims that 91% of cells oscillated in every MEC session.
- **Coordination**, not merely slow periodic firing, is the important result: appropriately sorted cells activate sequentially with a reproducible relative ordering. Fig. 2 shows repeated slanted raster bands, raw fluorescence confirmation, independent spiking confirmation, and ring-like trajectories.
- Fifteen of 27 MEC imaging sessions were classified as oscillatory; 421 individual sequences were identified. Durations range from tens of seconds to minutes, more variable across sessions than within a session. Some sequences are interrupted; 279/406 inter-sequence intervals were zero.
- 5,841/6,231 cell observations (93.7%) in oscillatory sessions were significantly phase locked. Cells can nevertheless skip individual sequences; phase locking and participation are distinct quantities.
- Sequential activation of three or more coactive ensembles was 2.3 times more likely than in shuffled data. Sequence score was significant in all 15 oscillatory sessions and in 5/12 non-oscillatory sessions.
- Sequences persisted through immobility, including continuous epochs up to 258 s, and one sequence could accompany 0–86 wheel laps. They were more prevalent during running: median conditional probabilities 0.93 running versus 0.69 immobility in the synchronized behavioural subset (10 sessions, 3 mice). Thus not a simple wheel-position or distance code, but not independent of behavioural state either.
- Similar oscillatory sequences were absent in sampled parasubiculum (0/25 sessions) and visual cortex (0/19), even though some individual neurons there had ultraslow oscillations.
- Preferred phase and anatomical distance showed no meaningful travelling-wave organization in the recorded fields. This does not rule out every spatial organization at unsampled scales.

### Precisely what “low dimensional” means here

Let X contain neurons as rows and time bins as columns. For sorting, PCA uses the binary deconvolved calcium-event matrix at 129 ms resolution, all movement and immobility epochs, neurons as variables and time bins as observations. For cell i, take its **eigenvector components/loadings**, without multiplying by eigenvalues, on PC1 and PC2. Sort cells by the polar angle of that pair:

```latex
\theta_i=\operatorname{atan2}(\ell_{i2},\ell_{i1}).
```

This is a fixed ordering of cells, not a time coordinate. The instantaneous population phase instead uses the **PC scores of the population activity vector**:

```latex
\phi(t)=\operatorname{atan2}(s_2(t),s_1(t)).
```

The paper writes arctg of a ratio but specifies the full interval [−π, π); atan2 here makes the required quadrant-aware interpretation explicit, not a claim about inspected source code. PCA sign ambiguity allows reversal of sequence direction in plots.

The proposed collective trajectory has approximately one intrinsic phase coordinate, embedded in two displayed dimensions. **Two PCs for sorting/visualization do not mean that all recorded variance is rank two.** The Methods explicitly retain full activity matrices/all components for analyses except visualization; sequence-to-sequence variability can blur the PCA ring. No Lyapunov exponents or rigorous attractor-topology proof are reported.

Visualization has additional preprocessing: Gaussian smoothing with width four times the session's “oscillation bin size” for the PCA manifold; for the Laplacian-eigenmaps view, resampling to 0.52 s, PCA denoising to five PCs, then LEM with k=15 and σ=2. LEM and the alternative sorting methods have hyperparameters. The unsupervised PCA sorting itself avoids those nonlinear-embedding choices. First-half/second-half sorting comparisons and shuffled controls test whether ordering is merely an artefact of fitting the displayed data.

### Quantification is more than seeing a ring

- The **oscillation score is not literally fraction of recording time oscillating**, despite shorthand main-text language. First require a prominent nonzero spectral peak in sin φ(t). Then compute the joint distribution of cross-correlation peak lag and cell-pair angular separation. For each of 11 angular-distance bins, test whether the conditional lag distribution has periodic spectral structure. The score is the fraction passing the specified peak tests. Threshold 0.72 selects the high group in a bimodal distribution; thresholds between 0.27 and 0.72 give the same classification (Methods p. 14).
- The session oscillation bin is one tenth of the period inferred from the peak of the spectrum of sin φ, ranging 3–17 s across the oscillatory sessions.
- Individual sequences are identified by a smoothed phase traversing a full 2π with limited backwards fluctuations; beginnings/endings are arbitrary choices of phase origin (p. 15).
- Cell locking is the length of the average unit phasor over that cell's event times, assessed against 1,000 temporal shuffles and their 99th percentile. A bias-corrected mutual information between event count and phase is a separate consistency check (p. 16).
- Ten ensembles are formed from nearby cells in PCA order. The most active ensemble at each coarse time bin gives an ensemble trajectory. Reported transition weights are counts divided by **all transitions**, not row-normalized conditional Markov probabilities. Repeated self-transitions are discarded. Sequence score counts increasing runs/subsequences of at least three ensembles and is compared with shuffled controls (p. 19).

### Mechanism and function remain open

The observations are compatible with ring-attractor-like dynamics, but the paper explicitly discusses recurrent networks, feedforward propagation and slow plasticity as alternatives. They do not experimentally identify the generating connectivity, establish a specific cellular source of the minute timescale, or demonstrate a causal memory function.

Its illustrative readout is 500 prescribed sequential Gaussian input units feeding **one linear output**, R=WX. Trained weights reproduce a ramp or an Ornstein–Uhlenbeck-generated target over a roughly 100 s interval if the input sequence repeats slowly enough. A fast repeated input cannot distinguish different cycles for a memoryless readout. This is a model-based plausibility argument for a temporal template, not a recording of a downstream memory circuit.

Source inconsistency worth preserving: the Extended Data Fig. 12 caption says the longest input sequence is 400 s, whereas Methods specify a 240 s template and the error plot extends to approximately 240 s. Avoid quoting a unique longest duration without checking code.

## 2. Recurrent structure in the 2026 preprint

### Model, fitted quantities, and numerical evidence

The rate network has sigmoid outputs and noisy leaky activation dynamics (main Eq. 1, p. 21):

```latex
\begin{aligned}
r_i&=\phi(h_i),&
\phi(h)&=(1+e^{-h})^{-1},\\
\tau\dot h_i
&=-h_i+\sum_j W_{ij}r_j+\zeta_i.
\end{aligned}
```

Wᵢⱼ is from presynaptic j to postsynaptic i; τ=0.3224 s. Synthetic networks have 501 neurons, target period 8.5 s, and phase-shifted wrapped Gaussian firing profiles. Width means the duration of each neuron's activation bump, **not the full population sequence period**. Amplitude, baseline/offset and width are varied. Fifty random initializations/snippets are trained with a FORCE variant, 300 training iterations followed by 10 fixed-weight test iterations.

The connectivity profile averages all weights at the same circular separation d in sequence ordering:

```latex
c(d)=\frac1N\sum_{i=0}^{N-1}
W_{i,(i+d)\bmod N}.
```

This is a reduction from N² microscopic weights to N phase-distance averages. It is not averaging by anatomical distance.

- Fig. 1: two example learned matrices correlate only 0.07, whereas their profiles correlate 0.98. The across-realization mean matrix itself generates a sequence.
- Fig. 2: increasing baseline shifts profile asymmetry; increasing bump width reduces profile range; larger amplitude reveals stronger/higher-frequency structure in the profile. Asymmetry supports directed bump progression.
- Fig. 3 compares networks trained on recordings with networks trained on Gaussian approximations to those recordings: one MEC session with 484 cells; a bat CA1 dataset with 168 time cells pooled across sessions and averaged across hanging intervals. Profile correlations are 0.95 (mouse) and 0.83 (bat), while matrix correlations are −0.09 and 0.13.
- These are **model-inferred couplings**, not measured synapses. Both compared connectivities are obtained by training the same model class. No new experimental perturbation or connectivity dataset was collected.
- Dale's law was not imposed. Band-like weights in experimental-data fits are suggestive but do not establish biological excitatory/inhibitory connectivity.

### Exact content of the analytical result

Supplement pp. 1–7 gives the clearest and more qualified statement. For prescribed rates f(t), set x*(t)=logit f(t). In the noiseless model, exact support requires

```latex
\begin{aligned}
y^*(t)&=x^*(t)+\tau\dot x^*(t),\\
Wf(t)&=y^*(t).
\end{aligned}
```

Assumptions include uniformly phase-shifted identical periodic bump shapes, rates bounded away from sigmoid saturation, and a retained Fourier cutoff D chosen to approximate **both** f and y* at tolerance ε. Require D<N/2 to avoid aliasing of retained modes. A smooth Gaussian is not exactly finite rank: all its Fourier coefficients are nonzero; “low rank” is accuracy dependent.

Cyclically average W using the shift operator S:

```latex
W_{\rm cp}=\frac1N\sum_{r=0}^{N-1}S^{-r}WS^r.
```

This circulant matrix has the profile as its cyclic diagonals and acts like W on the supported truncated sequence. On retained Fourier vectors vₖ, the required eigenvalues are fixed:

```latex
\lambda_k=
\frac{(1+ik\omega\tau)\widehat{x^*}_k}
{\widehat f_k},
\qquad |k|\le D,\quad \omega=2\pi/T.
```

Writing P_D for the active-subspace projector and W_D for the canonical low-pass circulant extension, all matrices supporting the **D-mode constraint** are

```latex
W=W_D+\widetilde W,
\qquad \widetilde W P_D=0.
```

Rank of W_D is at most 2D+1, equal when every retained eigenvalue is nonzero. The free part includes high-frequency circulant modes as well as noncirculant structure. Thus **only the active/low-pass profile is fixed**, not the full connectivity profile in arbitrary supporting networks. Similarity of the full learned profiles is a numerical observation, interpreted through regularization.

Supplement pp. 7–9 studies an idealized ridge regression: each required eigenvalue is multiplied by the shrinkage factor below:

```latex
\lambda_k^{(\alpha)}
=\frac{|\widehat f_k|^2}
{|\widehat f_k|^2+\alpha}\lambda_k.
```

This explains suppression of weak high-frequency modes. It is an idealized learning calculation, not a proof that finite FORCE training reaches this global optimum.

**Scope caution from the mathematics:** satisfying the on-trajectory linear constraint does not prove the trajectory is a stable attractor, nor that an arbitrarily large residual preserves accurate full nonlinear dynamics away from the retained subspace. Perturbation stability is tested numerically in selected trained networks, not supplied by this degeneracy theorem.

### Perturbation mechanism

Fig. 4/main pp. 6–7: silencing currently active neighbouring neurons in sequence order can reset activity to another ensemble. Removing their negative outgoing contributions disinhibits recipients; connectivity-profile minima help predict the new bump location. Multiple minima can support transient parallel bumps. Silencing inactive ensembles has little effect in the synthetic examples.

Sufficiently strong injected current instead starts a new sequence near the stimulated ensemble. Distributed random noise often leaves subsequent sequence progression intact. These claims concern simulated networks and tested parameter ranges.

The mouse-recording-trained network is an important exception to a simple robustness narrative: after perturbation it can reset and then settle into a non-progressing attractor state (Fig. 4H, legend p. 14). The paper's tests on “experimental data” perturb the fitted RNN, not the animal.

### Downstream computation and its real limitations

Fig. 5/main pp. 34–36: 501 recurrent upstream units feed ten downstream units, with **no downstream recurrence**, via a learned linear matrix A:

```latex
z(t)=A\,u(t).
```

Training uses 1,000 gradient-like iterations of averaged error times input (learning rate 0.015). The same upstream sequence supports many downstream amplitudes and offsets. Wider downstream bumps are fitted accurately; narrower ones are less accurate (displayed examples: upstream width 1 s, downstream 2 s gives correlation 0.99; downstream 0.5 s gives 0.87). This is degradation, not complete absence of a sequence.

Supplement pp. 9–11 sharpens two distinct restrictions:

1. A time-independent linear output of a Tᵤ-periodic input is itself Tᵤ-periodic. For nonconstant outputs with fundamental period T_d, Tᵤ=qT_d for an integer q≥1. It cannot create a new, longer fundamental period without extra state or input.
2. At finite approximation accuracy, a narrower Gaussian output is **not mathematically excluded**. It requires amplifying the input's weak high-frequency Fourier modes. In the exact continuum/increasing-population limit, with fixed population-size ratio and uniformly bounded readout operator norm, Gaussian widths must satisfy σ_d≥σᵤ. The supplement explicitly states the finite-accuracy qualification.

This width limit is separate from the period limit. At equal periods it is the familiar ill-conditioning of deconvolution: sharpening a Gaussian requires exponentially increasing high-frequency gain. No universal biological impossibility theorem follows for nonlinear, recurrent, time-dependent or externally driven downstream circuits.

### Small source discrepancies

- Main Methods p. 21 says 170 bat neurons; detailed preprocessing pp. 23–24, Results p. 5, Fig. 3 and Fig. 4 consistently say **168**.
- Fig. 5D legend on p. 15 accidentally repeats width 2 s; the actual Fig. 5D on p. 45 and Methods p. 35 specify **0.5 s**.
- Main training error Eq. 7 is printed as S−H_target, whereas exact dynamical support in the analytical supplement includes the derivative term x*+τ dx*/dt. The supplement's idealized inverse problem should not silently be substituted for the implementation described in main Methods.
- Some perturbation parameter listings/figure references in p. 33 are inconsistent with main figure legends. Mechanistic summaries above do not rely on a single supposedly exact noise/current setting.

## Suggested short synthesis

Gonzalo Cogno's sequence work asks how large neural populations organize activity over behavioural timescales. In MEC she and collaborators observed repeated slow sequences whose collective activity traces a ring-like trajectory, approximately tracked by one phase. Later modelling connects these trajectories to a shared low-frequency, circulant component of recurrent coupling: very different microscopic networks can implement the same sequence. A downstream linear readout can use that sequence as a temporal basis, but its period and usable bandwidth limit the patterns it can generate. The experiments establish the phenomenon; the proposed recurrent mechanisms and memory-scaffold function remain model-based explanations and hypotheses.
