# Deffner: quantum speed limits and information production

## Bibliography
- Sebastian Deffner, *Quantum speed limits and the maximal rate of information production*, Physical Review Research 2, 013161 (2020).
- Verified metadata: https://arxiv.org/abs/1910.06811 ; arXiv:1910.06811v3; journal DOI 10.1103/PhysRevResearch.2.013161.
- Local PDF: `papers/Deffner_2020_Quantum_speed_limits_information_production_1910.06811v3.pdf`; downloaded from https://arxiv.org/pdf/1910.06811v3.
- Checked PDF/printed pp. 1–7 by pypdf extraction, concentrating on pp. 2–5; dedicated PDF reader returned text but the viewer reported an unsupported-format error. No reliance on figure interpretation. PyMuPDF unavailable; no package installed.

## Setup and main results
[From source] Eq. (6), p. 2 defines trace distance ell=||rho-sigma||_1/2. Eqs. (8)–(10) bound its speed by ||dotrho||_1/2 and define Lambda_tau=(1/tau) integral ||dotrho||_1 dt and tau_QSL=2ell/Lambda_tau.

[From source] For finite-dimensional nonunitary dynamics, p. 3, Eqs. (13)–(16), the paper combines a Fannes bound |Delta H_vN|<=2ell ln d+1/e with the QSL. The quantity denoted dot I is |Delta H_vN|/(tau_QSL ln2), not a pointwise derivative or generally |Delta H_vN|/(actual elapsed time ln2). Eq. (15) drops the additive term and uses an approximate inequality. The exact retained inequality has an extra 1/(e tau_QSL ln2). Large d alone does not justify neglect at arbitrarily small endpoint distance.

[From source] P. 4, Eqs. (22)–(23) discuss energy-constrained infinite-dimensional extensions; no extension of our finite-d pure-state theorem was inferred from these.

[From source] Pp. 4–5, Eqs. (24)–(33) discuss Shannon entropy of measurement distributions in an observable eigenbasis under unitary evolution. Eq. (33) uses a distribution-dependent prefactor and a time-averaged marginal speed; its dot I_X is again |Delta S_X|/(tau_QSL^X ln2).

## Direct connection to our result
[Proved from our all-outcome theorem] Finite-d pure-state unitary dynamics, fixed finite POVM, full trace volumes, natural logs:

`|dot S_obs| <= C_d DeltaH/hbar = (C_d/2)||dotrho||_1 = (C_d/2)sqrt(F_Q,time)`.

Here C_d is the sharp binary coefficient in `../theorems/all-outcome-obsen-global-maximum.md`. The norm identity follows because dotrho has two nonzero eigenvalues +/-DeltaH/hbar. Applying the theorem to the reversed tangent bounds the negative derivative too. Sharpness is instantaneous, not generally sustained on an interval.

For a single fixed POVM along a smooth pure-state unitary trajectory:

`|S_obs(tau)-S_obs(0)| <= (C_d/2) integral_0^tau ||dotrho||_1 dt = (C_d/2)tau Lambda_tau`.

Hence, if Lambda_tau>0, `tau >= 2|Delta S_obs|/(C_d Lambda_tau)`; for constant H, `tau >= hbar|Delta S_obs|/(C_d DeltaH)`. If the speed vanishes throughout, Delta S_obs=0 instead of dividing by zero. Our bound is an entropy-change speed limit in the paper's terminology, but it is not obtained by identifying actual time with tau_QSL.

Rank-one PVMs have V_i=1, so OE equals the Shannon entropy in Eq. (24). Thus our theorem directly bounds its finite-d pure-state case; general POVMs additionally contain sum p_i ln V_i. A binary projector/complement maximizes OE, not generally binary Shannon entropy when d>2. No unrestricted POVM outcome Shannon-entropy bound follows just by dropping the volume correction.

## Mathematical caution on the unitary section
[Checked mathematical issue] P. 5, Eq. (26) calls an Lp norm of density differences a Wasserstein-p distance. The genuine Wasserstein distance minimizes transport cost over couplings and is not that Lp norm. Its monotonicity is W1<=W2 (when finite), opposite to Eq. (28). Therefore the chain (25)–(29), and its use in (33), cannot be adopted as written. Ref. [65] itself has NOT been separately audited; this observation concerns this paper's stated formulas, not a claim about the cited theorem. General finite-simplex Shannon entropy is not uniformly Lipschitz in L1 near zero probabilities: for p=(1,0), q=(1-epsilon,epsilon), H(q)/||p-q||_1 diverges. Our pure-state amplitude-speed bound avoids that invalid inference.

## Interpretation and limitations
Describe our result as a sharp finite-dimensional pure-state observational-entropy speed limit, in the spirit of this paper's information-speed-limit programme. Do not claim a general strengthening of all its results: entropy, state class, instantaneous versus finite-time rate, and normalization differ. Neither result alone proves channel capacity, thermodynamic entropy production, or a chaos/Lyapunov relation. Literature priority for our exact coefficient has not been checked.

Standalone connection: `latex/obsen_deffner_speed_limit_connection.tex`, not compile-tested.
