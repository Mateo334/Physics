# Wang–Wang–Wu: physical transport distance

## Bibliography and files

Zhenduo Wang, Yijie Wang, Biao Wu. Preprint title: *Genuine Quantum Chaos and Physical Distance Between Quantum States*, arXiv:2007.06244v2 (2020), https://arxiv.org/abs/2007.06244. Journal reference verified on arXiv: Phys. Rev. E 103, 042209 (2021); DOI https://doi.org/10.1103/PhysRevE.103.042209.

Downloaded PDF: `papers/Wang_Wang_Wu_2020_Physical_Distance_Quantum_Chaos.pdf`, from https://arxiv.org/pdf/2007.06244v2. Checked printed/PDF pp. 1–5 with pypdf on 2026-11-25; later numerical sections not independently checked.

## Setup and claims

[From source] P. 2, §II, Eqs. (1)–(3): choose an orthonormal basis with a ground metric and compare basis-outcome probabilities using classical Wasserstein distance. Mixed states use p_i=Tr(|xi_i><xi_i|rho). Basis and metric are physical choices, not unitarily invariant constructions.

P. 2 explicitly acknowledges this is not a metric on Hilbert-space states: distinct phase-superposition states can have identical probabilities and zero distance. The source also distinguishes prior Husimi/Monge constructions from its basis-probability construction.

P. 4, §IV A, Eq. (11): writes a quantum Lyapunov expression with infinitesimal and long-time limits. Crucial caveat directly after Eq. (12): t→infinity is NOT intended strictly, but as a sufficiently long time before wavepacket spreading; discussion invokes Ehrenfest time and semiclassical correspondence.

Pp. 1–2 and 5 introduce kicked rotor, three-site Bose–Hubbard and XXZ examples. No numerical result independently reproduced.

## Relation and limitations

[Assessment] A relevant lead for finite-time sensitivity in physically structured probability geometry, not a proof of intrinsic unitary-state stretching or an autonomous Hamiltonian flow of probabilities. Hidden coherences can produce nonclosure; denominator zero and nonsmooth optimal transport require care in a differential construction. The paper's use of “genuine chaos” is its interpretation, not a theorem verified here. No ALF/CNT/KS/Pesin bridge established by this reading.
