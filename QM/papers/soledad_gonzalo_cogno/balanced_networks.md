# Gonzalo Cogno–Mato: balanced networks, orientation selectivity and plasticity

Verified 2026-11-24 against the complete texts, not just their abstracts. Equations below are transcribed with explicit sign conventions and minor typographical cleanup. Numerical results are the authors' simulations, not independently reproduced results.

## Sources, files and reading record

1. **Soledad Gonzalo Cogno and Germán Mato (2015), “The effect of synaptic plasticity on orientation selectivity in a balanced model of primary visual cortex.” Frontiers in Neural Circuits 9:42**, published 20 August 2015; received 7 May, accepted 4 August. DOI: **10.3389/fncir.2015.00042**. PMCID: PMC4542321; PMID: 26347615. Original research, CC BY.
   - DOI: https://doi.org/10.3389/fncir.2015.00042
   - Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4542321/
   - Publisher: https://www.frontiersin.org/journals/neural-circuits/articles/10.3389/fncir.2015.00042/full
   - Successful PDF download: https://www.frontiersin.org/journals/neural-circuits/articles/10.3389/fncir.2015.00042/pdf
   - Local PDF: `papers/soledad_gonzalo_cogno/cogno_mato_2015_front_neural_circuits_9_42.pdf` (7,087,476 bytes, 13 pages).
   - **Actually read:** PDF-extracted text of all pp. 1–13, including Introduction (pp. 1–2), Materials and Methods and all model/plasticity/reconnection definitions (pp. 2–5), reduced-model Results (pp. 6–7), network Results and captions of Figs. 3–9 (pp. 7–10), Discussion (pp. 10–12), supplementary-material notice and references (pp. 12–13). Cross-checked the main equations and context against PMC HTML. Figure captions and associated result text were read; plotted data were not digitized.
2. **Soledad Gonzalo Cogno and Germán Mato (2015), “Orientation selectivity in a model of primary visual cortex with and without orientation map.” BMC Neuroscience 16(Suppl 1):P160**, published 18 December 2015. DOI: **10.1186/1471-2202-16-S1-P160**. PMCID: PMC4697541. Open-access, CC BY 4.0, **one-page poster presentation/conference abstract**, CNS*2015, Prague, 18–23 July 2015.
   - DOI/publisher: https://doi.org/10.1186/1471-2202-16-S1-P160
   - Full text: https://pmc.ncbi.nlm.nih.gov/articles/PMC4697541/
   - Successful PDF download: https://pmc.ncbi.nlm.nih.gov/articles/PMC4697541/pdf/1471-2202-16-S1-P160.pdf?download=1
   - Local PDF: `papers/soledad_gonzalo_cogno/cogno_mato_2015_bmc_neuroscience_16_s1_p160.pdf` (92,098 bytes, 1 page).
   - **Actually read:** complete PDF p. 1, all four body paragraphs, publication/conference metadata and all three references; also complete PMC HTML. This is its entire published text, not an unseen longer research article. It gives no numbered equations and no separate Methods section.

Download/extraction: local recursive searches found no existing matching sources. PDFs were saved with exclusive creation, preserving existing files. `fitz` was unavailable; installed `pypdf` read both PDFs, with no installations. Frontiers extraction reported duplicate `/MediaBox` keys but extracted all 13 pages. Publisher/PMC requests for the BMC PDF initially returned HTML rather than PDF; the successful PMC response was checked for the `%PDF` signature and parsed. Temporary Frontiers extraction is outside the workspace.

## Main physical picture

The model studies recurrent excitation and inhibition in visual-cortical layer 2/3, driven by weakly orientation-tuned layer-4 input. Strong excitatory and inhibitory currents dynamically cancel at leading order. Thus recurrent connectivity must reproduce and cancel the spatial/orientation profile of the external drive: it acts as a convolution constrained by balance. In the tractable stable regime, increasing the relative orientation modulation of the recurrent kernel reduces, rather than amplifies, the response modulation. This is distinct from a recurrent-amplification account of orientation tuning.

Orientation-map and salt-and-pepper refer to the relation between **physical position and preferred orientation**, not to the presence or absence of orientation selectivity. Local physical connectivity implies functional connectivity in an orientation map, but need not do so in a salt-and-pepper layout.

## Network and scaling [Frontiers pp. 2–5, Eqs. 1–6, 13–19]

- Two-dimensional periodic square patch, approximately 1 mm per side; 8,100 excitatory and 2,025 inhibitory leaky integrate-and-fire neurons by default. Each neuron receives on average K connections **from each presynaptic population**, with K = 500 by default; this is not 500 total recurrent inputs.
- Membrane time constant 20 ms, resistance 38.3 MΩ, threshold 30 mV, reset 0 mV. Recurrent synaptic currents have exponentially decaying kernels, normalized by the synaptic time constant; excitatory/inhibitory time constants 25/4 ms.
- Physical connection probability is a product of periodic Gaussians of the two coordinate differences, normalized to the desired indegree; default width 0.2 times the patch side.
- Layer-4 input uses 250 synapses, efficacy 1.65 nA ms, and average rate 15 Hz. Its orientation profile is:

```latex
I_{L4,iA}(\vartheta)
=g_{AL}K_{AL}f_L
\left[1+2\rho_A\cos 2(\vartheta-\vartheta_{iA})\right].
```

Here default ρE = ρI = 0.06. Preferred feed-forward angles are either independent uniform variables on [0,π), or an explicitly prescribed periodic pinwheel map (Eq. 3). The map is imposed, not learned.

The recurrent-coupling and background-current scaling is:

```latex
\begin{aligned}
g_{AB}&=\frac{G_{AB}}{\sqrt{K}},\\
I_{\mathrm{back},iA}&=\sqrt{K}\,I_{\mathrm{back},A}.
\end{aligned}
```

In the **reduced model**, the external current explicitly scales as:

```latex
i_{\mathrm{ext},A}(\theta)
=\sqrt{K}\,I_{\mathrm{ext},A}(\theta).
```

Consequently mean excitatory and inhibitory components are of order √K, whereas their balanced residual is order one. The reduced-model limit assumes 1 ≪ K ≪ N. The full model specifies its feed-forward parameters separately; do not silently claim that Eq. 2 explicitly imposes the reduced model's external-current scaling for every K variation.

**Sign convention:** the full simulation uses signed GEE = 32, GEI = −96, GIE = 96, GII = −128. In the reduced equations inhibition is explicitly subtracted, so JAI and the corresponding GAI must be interpreted as **positive inhibitory magnitudes**. Mixing these conventions gives wrong stability conditions.

## Fourier-space balance and stability [pp. 5–7]

The reduced coordinate θ lies on a unit-period circle of preferred orientations (physical angle divided by π), not generally on the physical cortical surface. Define the Fourier transform and circular convolution in this coordinate. Eqs. 17–19 give:

```latex
\begin{aligned}
\widetilde f(n)&=\int_0^1 f(\theta)e^{-2\pi in\theta}\,d\theta,\\
\mu_A(\theta)
&=\sqrt K\left[
J_{AE}*\nu_E-J_{AI}*\nu_I+I_{\mathrm{ext},A}
\right](\theta),\\
J_{AE}*\nu_E-J_{AI}*\nu_I+I_{\mathrm{ext},A}
&=O(K^{-1/2}).
\end{aligned}
```

The leading-order solution, with all Fourier coefficients evaluated at n, is [Eqs. 20–21]:

```latex
\begin{aligned}
D_n&=\widetilde J_{IE}\widetilde J_{EI}
-\widetilde J_{EE}\widetilde J_{II},\\
\widetilde\nu_E
&=\frac{\widetilde J_{II}\widetilde I_E
-\widetilde J_{EI}\widetilde I_I}{D_n},\\
\widetilde\nu_I
&=\frac{\widetilde J_{IE}\widetilde I_E
-\widetilde J_{EE}\widetilde I_I}{D_n}.
\end{aligned}
```

Here IE and II denote external inputs, not synaptic indices. The authors' mean-field rate stability test [Eqs. 22–24] requires negative-real-part eigenvalues of:

```latex
\begin{aligned}
A_n&=
\begin{pmatrix}
\widetilde J_{EE}(n)&-\widetilde J_{EI}(n)\\
\widetilde J_{IE}(n)&-\widetilde J_{II}(n)
\end{pmatrix},\\
\operatorname{tr}A_n
&=\widetilde J_{EE}(n)-\widetilde J_{II}(n)<0,\\
\det A_n&=D_n>0.
\end{aligned}
```

These are conditions within the reduced mean-field description, not a full stability theorem for the finite spiking network. For the zero mode, the positive-rate branch considered also satisfies Eq. 25:

```latex
\frac{\widetilde I_E(0)}{\widetilde I_I(0)}
>
\frac{\widetilde J_{EI}(0)}{\widetilde J_{II}(0)}
>
\frac{\widetilde J_{EE}(0)}{\widetilde J_{IE}(0)}.
```

## Orientation selectivity: the central analytic result

OSI is the modulus of the first angular Fourier coefficient divided by the mean [p. 4, Eq. 7]:

```latex
\mathrm{OSI}_i
=\frac{\left|\sum_{\vartheta}
f_i(\vartheta)e^{2i\vartheta}\right|}
{\sum_{\vartheta}f_i(\vartheta)}.
```

Simulations sample nine angles, 0° through 160° at 20° spacing. The range is [0,1] for nonnegative firing rates with positive total response.

For kernels whose shape depends only on the **presynaptic** population, and E/I external inputs having the same orientation profile, Eqs. 26–30 yield:

```latex
\begin{aligned}
\widetilde J_{AB}(n)&=G_{AB}\widetilde J_B(n),\\
I_{\mathrm{ext},A}(\theta)
&=I_A[1+2\rho\cos 2\pi(\theta-\theta_0)],\\
\mathrm{OSI}_A
&=\rho\left|
\frac{\widetilde J_A(0)}{\widetilde J_A(1)}
\right|.
\end{aligned}
```

The paper writes the aligned, positive-coefficient ratio without absolute values. This is the **inverse** of relative kernel modulation. Broadening an orientation-space kernel reduces its relative first harmonic and can sharpen the response needed to cancel the input.

**Singular-limit warning [p. 6]:** a perfectly flat salt-and-pepper mean kernel has zero first harmonic. The formula cannot be extrapolated to infinite OSI. Finite-connectivity fluctuations, which the authors describe as order K⁻¹ᐟ², then matter; their asymptotic balanced salt-and-pepper discussion also requires input modulation ρ to scale as K⁻¹ᐟ². This is not an invertible deterministic first-harmonic balance problem at fixed nonzero ρ. A nonnegative pure first-harmonic tuning curve additionally has OSI at most 1/2 (direct algebraic consequence, not a separate result claimed by the paper); numerical OSI above 1/2 is outside that simple harmonic ansatz.

## Plasticity [pp. 4, 6–9]

Only synapses **from excitatory neurons**, onto either E or I cells, undergo STDP. Initial normalized weights are 1. Inhibitory-source weights stay 1. With Δt = tpost − tpre:

```latex
\Delta w=
\begin{cases}
a_+e^{-\Delta t/\tau_+}(2-w),&\Delta t>0,\\
a_-e^{\Delta t/\tau_-}w,&\Delta t<0.
\end{cases}
```

Here a+ = 0.0128, a− = −0.0045, τ+ = 30 ms, τ− = 40 ms. Multiplicative factors soften updates near bounds 0 and 2. This changes existing weights, not the adjacency matrix. Learning runs for the first 13.33 s of a 20-s run, then weights are frozen for 6.66 s of measurements. Each run has one fixed stimulus orientation; the orientation ensemble is assembled across runs. Do not describe this as an explicitly rotating-stimulus training protocol.

The early-learning theory neglects spike-train cross-correlations, weight-dependent nonlinear corrections about w = 1, and feedback of weight changes onto rates. In the paper's Eq. 31, “Δw” means change per unit time; written explicitly as a rate:

```latex
\begin{aligned}
\alpha&=a_+\tau_++a_-\tau_-,\\
\dot w_{Ak,Ej}&\simeq\alpha\nu_{Ak}\nu_{Ej},\\
\overline{\dot w}_{Ak,Ej}
&\simeq\alpha\nu^\circ_{Ak}\nu^\circ_{Ej}
\left[1+2\xi_{Ak}\xi_{Ej}
\cos 2\pi(\theta_{Ak}-\theta_{Ej})\right].
\end{aligned}
```

The last expression is Eq. 34, averaged over stimulus orientation; ξ denotes initial response selectivity. The same pair-difference average follows for a fixed stimulus under the stated homogeneity assumptions. The default α is positive (0.204 ms with these units). The nonmodulated increment is order one in selectivity, the modulated increment order ξ². Accordingly STDP can strengthen connections between similarly tuned cells while sharpening responses because mean coupling grows sufficiently relative to its modulation. The analytic argument applies to early learning; a universal claim that STDP must improve tuning is not justified.

## Reconnection and instability [pp. 4–5, 7, 9–10]

Reconnection is a separate imposed mechanism, based on preferred orientations measured from responses:

```latex
\widehat p_{iA,jB}
=p_{iA,jB}
\left[1+\epsilon_{cB}
\cos 2(\mathrm{PO}_{jB}-\mathrm{PO}_{iA})\right].
```

It regenerates connectivity rather than merely changing weights. In the ideal homogeneous orientation-space model it preserves mean connectivity while modifying the first harmonic. With preferred response angles aligned to preferred inputs:

```latex
\begin{aligned}
J_B(\theta)&=1+\epsilon_{cB}\cos(2\pi\theta),\\
\widetilde J_B(0)&=1,\qquad
\widetilde J_B(1)=\epsilon_{cB}/2,\\
\mathrm{OSI}_A&=\frac{2\rho}{\epsilon_{cA}}.
\end{aligned}
```

The last line is an algebraic specialization of Eqs. 29–30 and 37–40, valid only on an admissible stable branch. Increasing own-population functional modulation therefore decreases own-population OSI there. Finite-network results also show cross-population changes not predicted by the separable reduced model; authors explicitly leave those unexplained.

The first-harmonic trace loses stability when [Eq. 41, positive inhibitory magnitudes]:

```latex
\frac{\epsilon_{cE}}{\epsilon_{cI}}
>\frac{G_{II}}{G_{EE}}=4.
```

At equality the trace is zero, so strict negative-trace stability is already lost. At εcE = 0.44 the boundary is εcI = 0.11. Once unstable, the inverse-modulation prediction cannot be used; nonlinear spiking dynamics can instead increase selectivity. This is not a proof that every unstable network becomes more selective.

## Reported numerical findings

| Finding | Reported result | Source |
|---|---|---|
| Balanced-current check | K = 250 vs 1000: larger cancelling E/I currents, similar net current; excitatory rates 3.02 vs 2.61 Hz, inhibitory rates 7.81 vs 5.71 Hz | pp. 7–8, Fig. 3 |
| Architecture, otherwise fixed parameters | Mean excitatory OSI 0.27 (map) vs 0.57 (salt-and-pepper); mean excitatory rates 2.51 vs 2.62 Hz | p. 8, Fig. 4 |
| STDP | Map OSI 0.27 → 0.36; salt-and-pepper 0.57 → 0.63 | p. 8, Fig. 5 |
| Stable reconnection | Increasing excitatory modulation reduces excitatory OSI; increasing inhibitory modulation reduces inhibitory OSI; mean rates change less than 5% in these comparisons | pp. 9–10, Fig. 8 |
| First-mode unstable reconnection | εcE = 0.44, εcI = 0.03 gives excitatory OSI 0.67 vs baseline 0.57; mean rates need not change much | p. 10, Fig. 9 |

Numerical method: Euler step 0.05 ms, checked by authors at 0.025 ms; typical runs 20 s, initial 0.1-s transient; some doubled-size comparisons (16,129 E and 4,096 I); STDP results compared between 20-s and 40-s runs. These are reported robustness checks, not new computations here.

Use the explicit OSI values rather than the Discussion's approximate “50%” map improvement: 0.27 → 0.36 is about 33%. Figure 6's caption calls its quantity final normalized efficacy, while the adjacent prose interprets fitted coefficients as efficacy **changes**; this note does not use that ambiguous label to infer absolute final weights.

## What the BMC item adds, and citation guidance

The complete one-page BMC text confirms the two architectures, balanced layer-2/3 and reduced models, excitatory-only STDP, OSI comparison 0.57 vs 0.27, broader-kernel sharpening, and the relative-modulation explanation of learning. It is a concise related conference account, not an independent replication or a second detailed derivation. Its phrase “ratio between the first Fourier component ... and its mean value” does not state the direction of dependence: use Frontiers Eqs. 29–30 to specify that response selectivity is **inverse** in that relative modulation. Reconnection and its instability threshold should be sourced to Frontiers.

### Compact introductory synthesis

Gonzalo Cogno and Mato model orientation-selective visual cortex as a strongly coupled excitatory–inhibitory network in which synaptic strengths scale as K⁻¹ᐟ² and order-√K inputs cancel to leave an order-one residual. A reduced balance equation diagonalizes in orientation-space Fourier modes. For stable, separable connectivity kernels and a common input profile, response selectivity is proportional to input selectivity and inversely proportional to the relative first harmonic of recurrent connectivity. Hence stronger functional connectivity need not sharpen tuning. Excitatory STDP can nevertheless increase selectivity when it reinforces mean coupling more than its orientation modulation; reconnection at fixed mean coupling typically decreases selectivity on the stable branch, but can enhance it after a first-harmonic instability.

These papers support an introduction about balance, structured disorder, sensory computation and plasticity. They do **not** calculate Lyapunov exponents, KS entropy, or an information-production rate; citing their balanced irregular activity as a demonstrated chaos/entropy result would overstate the evidence.

**Not checked:** separate Frontiers supplementary file, linked C implementation, independent simulation reproduction, or cited experimental papers. The finite-connectivity correction is reported from the main text, not independently verified from its supplementary derivation. Normalization of intermediate discrete kernels in Eqs. 13–16 is terse and inconsistent-looking across the source; use the explicit continuum balance and scaling equations above rather than copying those discrete formulas uncritically.
