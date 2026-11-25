# Calculations

- `code/`: Python scripts; shared imports remain alongside their callers.
- `code/four_model_setup.py` holds the four-model observational-entropy preparations; `code/chernoff_utils.py` evaluates classical outcome Chernoff entropy.
- `results/`: generated PNG/PDF figures, NPZ data, JSON summaries and CSV tables. Existing filenames are retained.
- [Calculation notes](../notes/INDEX.md#numerical-calculations): current summaries, interpretation and historical details.

## Entry points
Run from the workspace root with NumPy, SciPy and Matplotlib installed:

- `py -3.11 calculations/code/quantum_tangent_geometry.py` — exact inverted-oscillator Gaussian covariance/centroid/QFI comparison and seeded ordered polynomial derivative checks; JSON/CSV/PNG; refuses overwrite.

- `py -3.11 calculations/code/check_binary_obsen_global_bound.py` — sanity checks for the analytic binary-POVM global bound: scalar optimum, KKT semidefinite slacks/complementarity and 2000 random effects at each of five dimensions; proof in `latex/binary_obsen_global_maximum.tex`.
- `py -3.11 calculations/code/test_obsen_multioutcome_dual.py` — deterministic finite-grid sanity check of the global individual-effect certificate in seven dimensions; refuses overwrite. The analytic all-outcome proof is `latex/all_outcome_obsen_global_maximum.tex`, independent of this numerical check.
- `py -3.11 calculations/code/rotor_energy_100_coherent.py` — mean post-first-kick DeltaH times ln2 over 100 random-position, zero-momentum Gaussian packets; K4–12, d512, d1024 cutoff verification at every K, plot and standard-error band.
- `py -3.11 calculations/code/kicked_rotor_d2_energy_uncertainty.py` — d2 rotor toy at fixed hbar=2pi/64, K8–13, post-first-kick physical energy uncertainty and generator spread; analytic self-check.
- `py -3.11 calculations/code/plot_obsen_attained_rate_models_m4.py` — fixed four-outcome OE curves for saved d256 SYK/top/rotor states, with best attained initial-rate tangents; uses the saved dimension-sweep POVM and performs no new optimization.
- `py -3.11 calculations/code/plot_obsen_constant_m4_dimension.py` — four-outcome pure-state instantaneous OE constant for d=2 through 4096, exact small-block PGA, reuses saved d64/d512 results; PNG/PDF plot, CSV, JSON and compact POVMs, no full Hamiltonians.
- `py -3.11 calculations/code/kicked_models_d256_m10_compare.py` — d=256, ten-outcome instantaneous PGA comparison of kicked top K=3,6 and rotor K=6,10 against saved SYK; physical-time QFI, energy spread, full-effect checks and d=512 rotor cutoff check at fixed hbar.
- `py -3.11 calculations/code/syk_instantaneous_pga_pair.py` — two N=16 quartic-SYK disorder realizations, ten-outcome instantaneous OE PGA using the exact pure-state reduction, energy uncertainty and physical-time QFI; full-space effect and finite-difference checks, refuses overwrite.

- `py -3.11 calculations/code/syk_majorana_probe_qfi.py` — N=16 quartic SYK vacuum, final single-Majorana probe QFI through Jt=40; exact parity-based constant checked numerically.

- `py -3.11 calculations/code/kicked_top_preparation_qfi_pair.py` — one coherent preparation family, initial Jy parameter 0 and 0.07, two normalized-QFI panels under common evolution; exact invariance tested using propagated tangents.
- `py -3.11 calculations/code/kicked_top_time_translation_qfi_100_states.py` — same fixed-Jy QFI plot averaged over 100 sphere-uniform coherent states; preserves single-state outputs and saves individual QFIs and standard errors.
- `py -3.11 calculations/code/kicked_top_physical_time_qfi_100_states.py` — physical-time QFI for H=pi Jy then H=K Jz²/j, each lasting half a period; reuses the preceding saved 100-state ensemble.
- `py -3.11 calculations/code/plot_physical_time_qfi_phases.py` — readable two-panel view of the saved physical-time QFI, separating rotation and twist phase samples without overlapping vertical jumps.
- `py -3.13 calculations/code/kicked_top_time_translation_qfi_j10.py` — eight K values and 50 complete kicks for the j=10 pure-state fixed-Jy time-translation QFI; single plot, full CSV, NPZ and summary JSON. Not kicking-strength QFI.
- `py -3.13 calculations/code/kicked_top_instantaneous_povm_pga.py --outcomes 10` — j=10 instantaneous ten-outcome POVM optimization for K=3,6; omit flag for the earlier four-outcome run. Optional CVXPY/SCS; Dykstra projection fallback if CVXPY is unavailable. Saves distinct `results/kicked_top_instantaneous_povm_j10_m10_*` files.
- `py -3.13 calculations/code/kicked_top_ensemble_mixed_povm.py` — common four-outcome POVM on ten pure and ten mixed preparations; then `py -3.13 calculations/code/refine_kicked_top_ensemble_mixed.py` and `py -3.13 calculations/code/finalize_kicked_top_ensemble_mixed.py` for refined, strictly PSD results. Finalization refuses to overwrite existing final results.
- `py -3.13 calculations/code/kicked_top_30coherent_m10_highdim.py --spin 50 --restarts 3 --budget 180` — d=101, 30 coherent states, 10-outcome shared POVM; first saved pilot reached its time limit **without convergence**, see [calculation note](../notes/calculations/kicked-top-30-coherent-ten-outcome-j50.md). Script refuses to overwrite results for the same spin.
- `py -3.13 calculations/code/canonical_instantaneous_obsen_coefficient.py --dimension 1024 --outcomes 10` — exact 2D reduction for the unrestricted pure-state instantaneous coefficient; saves compact POVM and 16-start diagnostics in `results/canonical_obsen_coefficient_d1024_m10.{json,npz}`; a 1024-dimensional rotor Hamiltonian is not simulated.
- `py -3.13 calculations/code/canonical_instantaneous_obsen_coefficient.py --dimension 64 --outcomes 4 --restarts 12 --iterations 800` then `py -3.13 calculations/code/kicked_rotor_gaussian_free_obsen_d64.py` — d=64 torus Gaussian packet, measured free kinetic-energy uncertainty, and explicit four-outcome POVM slope. Existing results are not overwritten.
- `py -3.13 calculations/code/plot_kicked_rotor_pre_first_kick_obsen.py` — plot fixed-POVM OE for the d=64 Gaussian during free evolution from preparation at t=0 up to (not including) the first kick at t=1; saves a PNG and sampled NPZ.
- `py -3.13 calculations/code/canonical_instantaneous_obsen_coefficient.py --dimension 512 --outcomes 4 --restarts 12 --iterations 800` then `py -3.13 calculations/code/kicked_rotor_K10_postkick_instantaneous.py` — K=10 post-kick instantaneous four-POVM rates at t=1+,2+, comparing badly momentum-aliased d=64 with extended-cutoff d=512 at fixed hbar_eff=2pi/64. Saves PNG/JSON/NPZ and refuses overwrite.
- `py -3.13 calculations/code/kicked_rotor_fixed_povm_interval_slopes.py` — compute endpoint secant and 401-point OLS slopes with each saved POVM held fixed through both K=10 free intervals and the second kick; saves `results/kicked_rotor_K10_postkick_d64_d512_m4_fixed_interval_slopes.json`.
- `py -3.13 calculations/code/kicked_rotor_fixed_povm_0_to_3_optimize.py` — optimize a **single fixed four-effect POVM** for the discrete-time OLS OE slope on kicks 0–3 for K10; exact four-dimensional span reduction, d64 versus d768 at fixed hbar, plot/JSON/compact NPZ.
- `py -3.13 calculations/code/kicked_rotor_K10_linearity_penalty_d64.py --mu 1 --restarts 4 --steps 400` — fixed four-effect POVM, maximize first-to-second-kick secant minus 50-point derivative variance; analytical gradient, Dykstra, 501-point check and entropy/derivative plot. Refuses overwriting saved results.
- `py -3.13 calculations/code/kicked_rotor_K10_linearity_best_povm_full_curve.py` — plot S_obs from preparation through kick 3 using the already saved highest-objective POVM for the mu=1 linearity run, without reoptimizing.
- `py -3.13 calculations/code/kicked_rotor_d21_low_initial_fixed_povm.py --restarts 4 --steps 900 --beta 20 --margin .35 --warm-start calculations/results/kicked_rotor_d21_K10_K0_lowS0_T3_margin035_warm_r4s400.npz --tag margin035_refine_r4s900` — d21 periodic kicked rotor K10/K0 fixed four-effect POVM, three complete kicks, low-S0 and non-saturation penalties, analytical gradient, full-effect Dykstra; results/plot/POVM. For fresh, independent starts, omit `--warm-start` and use a new tag. Refuses overwrite.

- `py -3.13 calculations/code/compare_obsen_four_models_100_states.py` — 100-state observational entropy for SYK, kicked tops, Bose--Hubbard and GOE.
- `py -3.13 calculations/code/kicked_top_chernoff_observational_povm_ranks.py` — j=30, 100-state Chernoff observational entropy for several fixed PVMs and a full-rank unsharp Jz POVM. Saves `results/kicked_top_chernoff_observational_povm_ranks_j30_100_states.png` and `.npz`; the plot does not display a difference statistic.
- `py -3.13 calculations/code/kicked_top_fixed_povm_backaction_comparison.py` — fixed-POVM outcome KL with/without discarded-outcome backaction.
- `py -3.13 calculations/code/kicked_top_husimi_parameter_animation.py` — five-kick Husimi MP4 animation; requires ffmpeg. Options: `--resolution 540p --k-start 3 --k-stop 6 --fps 30`.
- `python calculations/code/kicked_top_relative_entropy.py` — fixed-Jz final-only outcome KL for 200 kicks.

Scripts contain numerical checks; saved plots and finite-time fits do not establish universal chaos or quantum Pesin relations. Read the linked notes before interpreting slopes.
