# Calculations

- `code/`: Python scripts; shared imports remain alongside their callers.
- `code/four_model_setup.py` holds the four-model observational-entropy preparations; `code/chernoff_utils.py` evaluates classical outcome Chernoff entropy.
- `results/`: generated PNG/PDF figures, NPZ data, JSON summaries and CSV tables. Existing filenames are retained.
- [Calculation notes](../notes/INDEX.md#numerical-calculations): current summaries, interpretation and historical details.

## Entry points
Run from the workspace root with NumPy, SciPy and Matplotlib installed:

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
