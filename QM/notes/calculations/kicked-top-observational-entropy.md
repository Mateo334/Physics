# Fixed-POVM observational entropy: larger kicked top

## Twenty-state average

[Numerically tested] `python calculations/kicked_top_entropy_state_average.py` evolves 20 independent spin-coherent initial states sampled uniformly on the sphere (cos(theta) uniform on [-1,1], phi uniform on [0,2pi)), seed 20261114. These are not Haar-random vectors in dimension 401. The same 20 states are used for k=4,5,6,7, j=200, 0–20 kicks; all other dynamics and POVMs below are unchanged. Plotted quantity is the arithmetic mean of the individual entropies, not the entropy of the mixture. Outputs have prefix `calculations/kicked_top_observational_entropy_mean20`: PNG/PDF, mean CSV, `_sem.csv` (sample standard deviation divided by sqrt(20)), `_individual.csv`, and `_summary.json` including all sampled angles and package versions. Shared plotting code is `calculations/plot_kicked_top_four_k.py`; means only are shown to avoid visual clutter, with independently scaled panel y axes.

Passed POVM and entropy checks, coherent-state normalization and Jz expectation, comparison of the binomial coherent-state formula with matrix exponentiation for one sampled state (1e-11 absolute tolerance), common initial entropy across k, and averaging check. Unitarity residual 5.11e-15; maximum state norm error 2.62e-13 (threshold 1e-10). Initial mean entropies: 3.24294, 3.99620, 4.39097, 4.90309 nats. Fine-PVM means at kick 20 for k=4,5,6,7 are 5.45568, 5.58057, 5.50966, 5.55096, with SEM 0.09085, 0.00590, 0.05610, 0.01239. Sampling error is not negligible; no convergence claim from only 20 states. The PNG was visually inspected.

Numerical correction discovered during this run: evaluating p*ln(p/v) through the ratio can underflow p/v for tiny nonzero probabilities in coarse bins, returning infinity. The shared entropy helper in `syk_observational_entropy.py` now evaluates -xlogy(p,p)+p*ln(v), an algebraically identical, more stable form. Added a smallest-positive-float regression check; it and the existing self-check pass. Previous output files were not overwritten. Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0; float64/complex128.

Presentation update: `calculations/plot_kicked_top_four_k.py` replots the existing CSV for k=4,5,6,7 and kicks 0–20, with distinct colors, dash patterns and hollow markers. Outputs `calculations/kicked_top_observational_entropy_k4_7.png` and `.pdf`. Each panel has an independently scaled vertical axis for readability. Data completeness/finiteness assertions passed; PNG visually checked. No dynamics rerun or data changed.

[Numerically tested] Run `python calculations/kicked_top_observational_entropy.py`. Reuses the entropy function from `calculations/syk_observational_entropy.py`. Outputs: `calculations/kicked_top_observational_entropy.png`, `.pdf`, `.csv`, `_early.png` (0–20 kicks), and `_summary.json`.

## Parameters and conventions

j=200, d=401; k=linspace(3,10,8); 100 kicks. U=exp(-i k Jz^2/(2j)) exp(-i pi Jy/2), so rotation acts first, matching `kicked_top_frechet.py`. Pure initial coherent state exp(-i phi Jz) exp(-i theta Jy)|j,j>, theta=1.1, phi=0.7. Unlike the prior derivative experiment, no initial-state mixing. Deterministic computation, no random seed. State vectors for all eight k values evolve together; no repeated-measurement backaction.

Four fixed POVMs: rank-one Jz projectors; ten contiguous bins in ascending m (first rank 41, others rank 40); and each unsharpened as E_a=0.75 P_a+0.25 Tr(P_a)I/401. Bin boundaries are recorded in JSON. Entropy is -sum p_a ln(p_a/Tr(E_a)), in nats; plotted as absolute entropy, not initial-baseline-subtracted growth. All volumes refer to the full spin space.

## Checks and observations

Float64/complex128: rotation unitarity residual 5.11e-15 (threshold 1e-11), maximum state normalization error 1.35e-13 (threshold 1e-10). Passed initial coherent-state Jz expectation check, effect positivity/completeness, trace volumes, sharpness checks, maximally mixed entropy, identity-evolution toy check, binary entropy toy check, entropy bounds, common initial entropy across k, and entropy ordering under coarse-graining and added noise. Density-matrix positivity is automatic from pure state vectors. The full-time PNG was visually inspected.

Initial entropies in panel order: 3.60613, 4.07068, 4.64104, 4.94444 nats. Curves can initially decrease and fluctuate; no monotonic-growth assertion. For this state the k=3 curve approaches its fluctuating plateau more slowly than higher k values. Coarse POVMs approach ln(401), while fine POVMs remain lower. No Lyapunov fit, initial-state averaging, or j-scaling check; these curves alone do not establish chaos.

Environment: Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0. No new packages or sources.
