# Single-neuron dynamics and coding: Gonzalo Cogno

Checked 2026-11-24. Full main texts read for the 2014 and 2016 papers. The 2013 paper remains **metadata-only, not full-text read**. No paywall bypass or package installation was used. PDF extraction tried `fitz` first (unavailable), then used installed `pypdf`; extracted working text is in the approved temporary directory.

## 1. Dynamics and reliability of bistable neurons driven with time-dependent stimuli

**Citation:** Soledad Gonzalo Cogno, Susanne Schreiber, Inés Samengo, *Neural Computation* **26**(12), 2798–2826 (2014).

- DOI: https://doi.org/10.1162/NECO_a_00671
- PubMed metadata: https://pubmed.ncbi.nlm.nih.gov/25248087/
- Open institutional manuscript: https://ri.conicet.gov.ar/handle/11336/21388
- Downloaded PDF: https://ri.conicet.gov.ar/server/api/core/bitstreams/bb0e8e7c-fc0d-4d58-aed2-17bd4b3228ad/content
- Local file: `papers/soledad_gonzalo_cogno/cogno_schreiber_samengo_2014_bistable_reliability.pdf`
- Version: 30-page author-formatted institutional manuscript, **not journal pagination**. OpenAlex identifies this repository deposit as submittedVersion, CC BY-NC-SA; exact equivalence to the published version was not checked.
- Publisher PDF request returned HTTP 403; the public CONICET DSpace record and its public bitstream supplied the full text.
- **Read:** manuscript pp. 1–28, §§1–8, equations (1)–(9), captions of Figures 1–12; bibliography pp. 29–30 inspected. Particularly pp. 4–5 (model/reliability), 6–7 (bifurcations), 8–15 (noise/transitions), 15–21 (periodic forcing), 21–28 (combined forcing and conclusions).

### Question and exact model specification recoverable here

How does the frequency of a repeated stimulus control spike reproducibility near a bistable firing threshold, where subthreshold and spiking frequencies differ?

The model is a **two-dimensional type-II Morris–Lecar conductance-based neuron**, with membrane voltage V and potassium activation variable w:

```latex
\begin{aligned}
C\dot V&=-I_{\rm ion}(V,w)+I(t),\\
\dot w&=\frac{w_\infty(V)-w}{\tau_w(V)},\\
I(t)&=I_0+I_1\sin(2\pi ft)+I_{\rm noise}(t),\\
\langle I_{\rm noise}(t)I_{\rm noise}(t')\rangle
&=\sigma^2\delta(t-t').
\end{aligned}
```

The manuscript writes the ionic current as leak, calcium and potassium terms, with the calcium term printed/extracted as `gCa(V−VCa)` and the potassium term as `gK w(V−VK)`. **Do not silently turn this into a fully specified standard Morris–Lecar implementation:** its voltage-dependent calcium activation and the parameter/gating functions are not explicitly supplied here; §2 delegates “other parameters” to Schreiber, Samengo & Herz (2009), *J. Neurophysiol.* 101, 2239–2251. That source was not independently read in this task. Thus model identity and equations of motion are verified, but a complete replication-ready parameter table is not.

Numerics: stochastic second-order Runge–Kutta (Honeycutt), step 0.1 ms. Typically 500 or 1000 independent-noise trials with identical initial state and deterministic stimulus. Representative operating current is 25.4 µA/cm².

For the reported parameter set (§3), a stable/unstable limit-cycle pair appears near 24.82 µA/cm²; the unstable cycle shrinks into a stable spiral fixed point near 26.23 µA/cm² in a subcritical Hopf bifurcation. Between these currents, the stable resting point and stable spiking cycle coexist, separated by an unstable cycle. Subthreshold oscillations are about 68 Hz; spiking is about 50 Hz. These are model-specific numbers, not universal physiological constants.

### Reliability is not mutual information

Each trial is represented by a binary spike sequence rᵢ(t), binned at 1 ms. The PSTH is its trial average. Equation (4):

```latex
\chi^2=
\frac{\operatorname{Var}_t[\langle r_i(t)\rangle_i]}
{\langle\operatorname{Var}_t[r_i(t)]\rangle_i}.
```

This is a variance-ratio synchronization/reliability statistic, with an independent-trial floor of order inverse trial count. It measures reproducibility under independent noise with common forcing; it is **neither Shannon mutual information nor a Lyapunov exponent**. All-silent output is displayed as zero reliability in their maps; the raw ratio would need a convention when its denominator vanishes.

### Mechanism and verified results

1. **Two sources of unreliability:** timing jitter within a spiking episode, and switches into the resting region that remove entire stretches of spikes (§7).
2. Input near the spiking frequency (and harmonics) locks the circulating trajectory, minimizes jitter and, for sufficient forcing, favors entering/staying in the spiking state. Input near the subthreshold frequency more easily destabilizes rest and generates mixed trajectories with complicated locking patterns (§§5–7).
3. Consequently there is **no smooth interpolation rule** between the optimum frequency below threshold and above threshold. Initial condition, stimulus phase, amplitude, noise and observation duration all matter.
4. If weak noise and weak forcing start a trajectory near rest, the nominally low-jitter spiking-frequency input may fail to elicit spikes during the recording. The subthreshold-resonant input can then be the more useful/reliable drive (§7, manuscript p. 27).
5. Noise is not simply detrimental: it can promote transitions, but it can also suppress transitions already generated by periodic forcing. This direction-dependent effect appears clearly for stronger sinusoidal forcing (§6, Figure 11).
6. The autonomous two-dimensional model does not intrinsically burst in this regime. Its apparent bursts arise from noise-induced alternation between rest and spiking (§4.2). This differs from the intrinsic slow–fast bursting model in the 2016 study.
7. At constant mean input, weak-noise firing-rate curves can exhibit apparent hysteresis and inverse stochastic resonance over finite windows. The manuscript argues that initial-state dependence and the reported inverse-resonance structure disappear in the infinite-time limit for its parameter regime, but stationarity may require unrealistically long recordings (§4.1). This is not established as a theorem for all bistable neurons.

### Useful mathematics: ISIs, residence weighting, transitions

The interspike-interval distribution P(τ) has a narrow first peak near the firing period and subsequent peaks separated roughly by the subthreshold period, followed by a long exponential tail. A threshold τlim at the minimum between the first two peaks separates short “suprathreshold” and long “subthreshold” intervals.

The fraction of intervals q and the fraction of elapsed time p are different. Equations (5)–(8):

```latex
\begin{aligned}
p_{\rm supra}
&=\frac{\int_0^{\tau_{\rm lim}}\tau P(\tau)\,d\tau}
{\langle\tau\rangle},\\
p_{\rm sub}
&=\frac{\int_{\tau_{\rm lim}}^\infty\tau P(\tau)\,d\tau}
{\langle\tau\rangle},\\
\nu_{\rm supra\to sub}
&=\frac{q_{\rm sub}}{p_{\rm supra}\langle\tau\rangle},\\
\nu_{\rm sub\to supra}
&=\frac{q_{\rm sub}}{p_{\rm sub}\langle\tau\rangle},\\
p_{\rm supra}\nu_{\rm supra\to sub}
&=p_{\rm sub}\nu_{\rm sub\to supra}.
\end{aligned}
```

These are operational ISI-based coarse-grained estimates: classifying a whole long interval as subthreshold is not exact continuous-state basin tracking. The same qsub appears in both rates because each long silent episode has both an entry and an exit. Rare long intervals may dominate elapsed time even when most observed intervals are short. Stationarity requires recordings much longer than **both inverse transition rates**, not merely many spikes.

### Source cautions

- Some manuscript details are internally inconsistent: Figure 9 defines n:m as n stimulus cycles to m spikes but its caption prints a reversed frequency ratio; use the verbal event-count definition, giving firing rate m f/n. The periods listed after equation (9) appear reversed relative to the approximately 50 Hz firing and 68 Hz subthreshold frequencies. Neither ambiguity is needed for the qualitative synthesis.
- Results are mainly numerical dynamical-system analysis, with local linearization and counting identities; they are not experimental demonstrations in this paper.
- “Crossing trajectories” in the V,w projection of the forced system does not establish deterministic chaos; forcing phase/time is missing from that projection. The paper supplies no chaos diagnostic.

## 2. Bursting neurons in the hippocampal formation encode features of LFP rhythms

**Citation:** Maria Constantinou, Soledad Gonzalo Cogno, Daniel H. Elijah, Emilio Kropff, John Gigg, Inés Samengo, Marcelo A. Montemurro, *Frontiers in Computational Neuroscience* **10**, 133 (2016). Constantinou and Gonzalo Cogno contributed equally.

- DOI: https://doi.org/10.3389/fncom.2016.00133
- Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC5183636/
- Publisher: https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2016.00133/full
- Successful PDF source: https://www.frontiersin.org/journals/computational-neuroscience/articles/10.3389/fncom.2016.00133/pdf
- Local file: `papers/soledad_gonzalo_cogno/constantinou_cogno_2016_bursting_LFP_rhythms.pdf` (18 pages, published version, CC BY).
- **Read:** complete main HTML §§1–4, including §§2.1–2.10 and §§3.1–3.6, equations (1)–(4), all main figure captions. PDF downloaded and text extraction verified.
- Supplement link: https://pmc.ncbi.nlm.nih.gov/articles/instance/5183636/bin/DataSheet1.pdf . This returned HTML rather than PDF; Europe PMC supplementary-files service timed out. **Supplement not obtained or read**, including its conductance equations and parameter tables. The ORO mirror PDF request returned HTTP 403.

### Model and experimental comparison

The model is a **two-compartment conductance-based pyramidal bursting neuron**, with a dendritic input compartment and somatic spike output (§2.3). The paper describes its lineage from the Traub CA3 model through Pinsky–Rinzel reduction and Kepecs–Wang minimal bursting conductances. It is not a Morris–Lecar model.

The actual current/gating equations are Supplementary Equations 1–2 and parameters Supplementary Tables 1–2. **Exact ionic equations and numerical parameters remain unverified here**; they must not be reconstructed by assuming that every model in this lineage uses identical currents.

The input is a phase-randomized surrogate of a measured LFP: Fourier amplitudes are preserved, Fourier phases randomized, and the inverse transform gives a signal with the same power spectrum but altered temporal realization. Signals are scaled to mean zero and standard deviation 0.7 nA (anesthetized/subicular spectral profile) or 0.4 nA (awake/MEC profile). Current enters the dendrite. This models a spectrally structured effective input; it does not identify measured extracellular voltage with an actual injected current through a calibrated physical conversion.

Experimental comparison:

- 28 bursting units in subiculum, four urethane-anesthetized rats; dominant approximately 1 Hz delta/slow rhythm.
- 42 bursting units in MEC layers III/V, **one awake foraging rat**; dominant approximately 8 Hz theta.
- Main feature filters: 0.5–3 Hz (delta), 6–12 Hz (theta). Voltage and derivative give instantaneous value and slope; analytic-signal/Hilbert phase and modulus give phase and amplitude. Phase zero is the oscillation maximum.
- Burst segmentation joins consecutive spikes with ISIs below 8 ms (subiculum), 5 ms (MEC), or 16 ms (slower simulated spike patterns). Single spikes are called one-spike bursts.

### Precisely what “information” means

Time bins have width δt = 5 ms. Each bin is labeled at burst onset with N = 0, 1, 2, or 3, where 3 pools all bursts of at least three spikes. An LFP feature X, possibly at a time lag, is discretized into four equiprobable bins.

Three codes (§2.8):

- **Full burst code:** silence plus distinct burst sizes.
- **Burst rate code:** silence versus any burst onset, ignoring size. “Rate” here is the binary occurrence representation, not a long-window spike-count rate.
- **Burst distinction code:** size only, conditioned on a burst having occurred; silent bins are discarded.

Information is Shannon uncertainty reduction about X by observing the neural symbol. In the short-bin rate formulation used by the paper:

```latex
\begin{aligned}
I_{\rm bin}
&\simeq \delta t\sum_n\sum_x p(x)r_n(x)
\log_2\frac{r_n(x)}{\bar r_n},\\
\bar r_n&=\sum_xp(x)r_n(x).
\end{aligned}
```

For an unambiguous first-order interpretation, the sum is over nonzero event types; silence contributes only higher-order terms in the usual small-bin expansion. Exact finite-bin mutual information instead uses all response probabilities, including silence. The paper’s surrounding text lists the silence symbol as part of its code alphabet without separately clarifying this rate-notation distinction.

Dividing information per bin by r δt (expected burst count per bin) gives bits/burst. Let B indicate burst occurrence, and K its size when B = 1. The chain-rule decomposition is:

```latex
\begin{aligned}
I(X;N)
&=I(X;B)+p(B=1)I(X;K\mid B=1),\\
I_{\rm full}^{\rm per\ burst}
&=I_{\rm rate}^{\rm per\ burst}
+I_{\rm distinction}^{\rm per\ burst}.
\end{aligned}
```

The conditional feature distribution in the distinction term is the feature distribution **at burst times**, not the unconditional four-bin distribution. This matters when comparing quantities with different conditioning/units. The identity explains the extra information retained by distinguishing burst sizes; it does not prove that the brain decodes it.

Estimation uses 100 shuffles to estimate sampling bias, subtracted separately for each cell. A cell is counted as significantly informative when real-data information exceeds the maximum shuffled value over the tested time window (§2.8). These are empirical bias-corrected estimates, not exact information capacities.

### Verified findings

1. Both simulated and recorded bursts carry information about the dominant rhythm’s voltage, slope and phase; amplitude is generally less strongly encoded (§§3.1–3.3).
2. **Most information is in burst onset timing.** Size adds a smaller, cell-dependent contribution; the discussion summarizes it as roughly 10–15% additional information. Detailed population mean ratios of distinction/full information vary by feature: subiculum 13.8%, 12.6%, 13.1%, 21.2%; MEC 16%, 16%, 23%, 13%, respectively for voltage, slope, phase, amplitude. Do not replace all those values with a universal 10–15% law.
3. Using the stronger threshold of significant information at least 0.1 bits/burst, 57.1% of subicular cells and about 40% of MEC cells encoded at least one feature. These are not the fractions passing significance alone.
4. Preferred phase shifts with burst size in some cells; others phase-lock without distinguishing sizes. Phase-locking alone therefore does not establish a phase-to-burst-size code (§3.5).
5. Most informative subicular cells and all informative MEC cells encoded maximally in the dominant-power frequency band; five subicular cells additionally encoded amplitude at higher frequencies (§3.6).
6. Information about voltage/slope versus lag can oscillate at twice the rhythm frequency. The proposed explanation is that burst-size-conditioned average waveforms cross twice per cycle, reducing distinguishability at those lags (§3.4).
7. Information about future LFP values is compatible with causality: temporally correlated input allows present responses to predict future features. The simulations explicitly demonstrate this without neural output causing the future drive (§4).

### Limits important for an introductory account

- Recorded LFP–burst associations are statistical, not a causal intervention on extracellular fields.
- Anatomical region and behavioral/anesthetic state differ together; this design does not isolate their separate contributions.
- Awake evidence comes from one rat; 42 cells are not 42 independent animals.
- Four feature values, pooled ≥3-spike bursts, discarded within-burst timing, and fixed lag/filters define a deliberately coarse code. It is not information in the complete spike train.
- Features are correlated. Information about several features must not be added as independent channels; the authors discuss a phase-centered explanation for MEC.
- The paper explicitly does not demonstrate whether/how downstream neurons use this information.
- Supplement-dependent biophysical details were inaccessible in this reading.

## 3. Effect of input noise on neuronal firing rate

**Citation:** S. Gonzalo-Cogno and I. Samengo, *AIP Conference Proceedings* **1510**, 247–250 (2013), in *Physics, Computation, and the Mind—Advances and Challenges at Interfaces*.

- DOI verified through Crossref title matching: https://doi.org/10.1063/1.4776516
- Author bibliography: https://gonzalocognolab.com/publications/
- Proceedings contents independently list the paper starting at p. 247: https://www.proceedings.com/content/020/020931webtoc.pdf
- OpenAlex lookup: https://api.openalex.org/works/https://doi.org/10.1063/1.4776516
- Repository metadata location reported by OpenAlex: https://hdl.handle.net/11250/4703110
- **Access status:** no lawful downloadable full text located. OpenAlex reported no open PDF and only publisher/NVA metadata locations. DOI fetch returned a citation, not article content. NVA API retrieval did not produce usable metadata JSON/full text. Exact-title searches and author/coauthor bibliographies did not yield an author manuscript.
- **Local PDF: none. Pages/sections read: none. Even an abstract was not treated as a full-text source.**

No model equations, definition of firing rate, noise normalization, reliability/information statistic, or scientific conclusion is verified for this specific paper. Do not attribute the 2014 inverse-stochastic-resonance analysis to this 2013 proceedings paper without obtaining it.

## Synthesis for a physicist’s introduction

The verified bridge is from **nonlinear stochastic dynamics** to **operational coding statistics**. In the 2014 work, a neuron is a forced bistable dynamical system: the stable focus and spiking limit cycle introduce different characteristic frequencies, while noise controls timing precision and switching. Finite recording windows can retain memory of preparation and produce qualitatively different measured rates and reliability. In the 2016 work, intrinsic bursting is examined as a coarse-grained output alphabet: burst timing indicates when a rhythmic feature occurs and burst size can add finer distinctions. Information is explicitly Shannon mutual information about measured/simulated signal features, not thermodynamic entropy production, KS entropy, or evidence of chaos.

Safe introductory formulation: “Gonzalo Cogno’s early single-neuron work connects the geometry of noisy, driven neuronal dynamics with the reproducibility and information content of spike patterns. Bistability creates competing resonances and noise-sensitive transitions, while burst timing and size can represent features of dominant hippocampal rhythms.”

Scope of that formulation: supported by full reading of 2014 and 2016; the 2013 citation is bibliographic context only. No simulations or reported statistical results were independently reproduced.
