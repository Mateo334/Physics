# Four-outcome instantaneous OE constant versus dimension

[Numerically tested] User requested four outcomes and a dimension sweep feasible on this laptop. Script: `calculations/code/plot_obsen_constant_m4_dimension.py`; run with `py -3.11` and OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1. Outputs: `calculations/results/obsen_constant_m4_dimension.{png,pdf,csv,json,npz}`. Runtime 45.47 seconds. Plot visually inspected.

Pure states, unrestricted four-effect POVMs, natural-log observational entropy, trace volumes on the full d-dimensional space. Plotted constant is the attained instantaneous rate divided by DeltaH/hbar, equivalently 2R/sqrt(physical-time QFI). Exact two-dimensional tangent-space reduction avoids full matrices and Hamiltonians; the complement is scalar and absent for d=2. Thus d=4096 is a conservative feasible stopping point, not a measured laptop memory limit. This is one model-independent coefficient sweep, not three Hamiltonian simulations.

Dimensions: 2,3,4,8,16,21,32,64,128,256,512,1024,2048,4096. Reused existing d64 and d512 four-outcome results; otherwise three starts (split binary seed plus two random starts), 200-iteration cap, seed 20261123+d, using existing reduced Dykstra/simplex PGA. Saved JSON retains all optimizer diagnostics. PSD roundoff repair and common normalization slightly reduce some raw best rates (below 3e-10); CSV and plot use repaired feasible candidates. Selected coefficients:

| d | Attained coefficient |
|---|---:|
| 2 | 1.32548683868 |
| 16 | 3.29771142247 |
| 64 | 4.57380836022 |
| 256 | 5.87832299640 |
| 1024 | 7.20628849939 |
| 4096 | 8.55038033046 |

Markers agree with the attainable split-binary candidate curve within numerical tolerance. Four outcomes permit proportional splitting of one binary effect; four independent measurement directions are not required. The plot uses a base-two logarithmic dimension axis. Observed growth is approximately logarithmic over this range; no asymptotic theorem about the true optimum was established.

Checks at every dimension: reduced effect PSD to 1e-12, completeness to 1e-12, scalar simplex normalization/nonnegativity, probability and volume normalization, directional analytic-gradient versus central-difference relative error below 1e-6, entropy-derivative absolute error below 1e-7. No full d-dimensional effects constructed. NumPy 2.2.3, SciPy 1.15.1; Python/Matplotlib versions recorded in JSON.

Limitations: attained lower bounds, not certified global optima. Both random starts at d2048 and one at d4096 hit the iteration cap substantially below the seeded candidate; the seed supplies the best feasible value there. All other fresh runs stopped with no improving step. This sweep does not show that four outcomes attain the unrestricted global maximum or that extra outcomes cannot improve it.
