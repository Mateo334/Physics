# Direct-evolution relative entropy

[Numerically tested] `calculations/code/kicked_top_relative_entropy.py` runs successfully and saves `calculations/results/kicked_top_relative_entropy.png`.

Both density matrices evolve by the same kicked-top Floquet unitary. The perturbed preparation is V rho V† with V=exp(-i epsilon Jy). The plotted quantity is D(p_reference||p_perturbed) for a fixed final-only rank-one Jz measurement, using natural logarithms. It is not quantum relative entropy or a measurement-history divergence. Quantum relative entropy is invariant under this common unitary evolution.

Parameters match the earlier derivative setup: j=10, dimension 21, epsilon=0.07, maximally mixed weight 0.3, coherent-state angles theta=1.1 and phi=0.7, rotation pi/2, kappa=0.5 and 6, times 0 through 200 kicks. Deterministic; no random seed required. Complex128 arithmetic. Assertions check unitary error at absolute tolerance 1e-12, state Hermiticity and trace at 1e-11, positive eigenvalues, finite divergence and nonnegativity to 1e-12, and exact zero for identical probability arrays. No fits or chaos inference.

Test environment: Python 3.11.9, NumPy 2.2.3, SciPy 1.15.1, Matplotlib 3.10.0. The script saves a PNG without requiring an interactive display.
