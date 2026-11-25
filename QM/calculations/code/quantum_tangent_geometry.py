"""Exact Gaussian instability and noncommutative tangent sanity checks.

Run: py -3.11 calculations/code/quantum_tangent_geometry.py
No Hilbert-space truncation is used for the Gaussian calculation.
"""
from pathlib import Path
import json
import platform

import numpy as np
import scipy
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    root = Path(__file__).resolve().parents[1] / "results"
    stem = root / "quantum_tangent_geometry"
    paths = [stem.with_suffix(s) for s in (".json", ".csv", ".png")]
    if any(p.exists() for p in paths):
        raise FileExistsError("Refusing to overwrite quantum_tangent_geometry results")
    root.mkdir(parents=True, exist_ok=True)
    # Dimensionless units: hbar=m=lambda=1; H=(p²-q²)/2.
    J = np.array([[0., 1.], [-1., 0.]])
    A = np.array([[0., 1.], [1., 0.]])
    sigma0 = np.eye(2) / 2
    v = np.ones(2) / np.sqrt(2)
    rows = []
    symplectic_error = purity_error = 0.
    for t in np.linspace(0., 3., 121):
        M = expm(t * A)
        sigma = M @ sigma0 @ M.T
        delta = M @ v
        qfi = delta @ np.linalg.solve(sigma, delta)
        assert np.linalg.eigvalsh(sigma).min() > 0
        symplectic_error = max(symplectic_error, np.max(np.abs(M.T @ J @ M - J)))
        purity_error = max(purity_error, abs(np.linalg.det(sigma) - .25))
        rows.append((t, np.linalg.norm(delta), qfi, np.exp(-.01 * qfi / 4)))
    rows = np.array(rows)
    assert np.max(np.abs(rows[:, 1] / np.exp(rows[:, 0]) - 1)) < 1e-12
    assert np.max(np.abs(rows[:, 2] - 2)) < 1e-9
    assert symplectic_error < 1e-10 and purity_error < 1e-10

    rng = np.random.default_rng(20261124)
    def hermitian():
        x = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
        return (x + x.conj().T) / 2
    Q, G = hermitian(), hermitian()
    comm = lambda x, y: x @ y - y @ x
    dQ = comm(Q, G) / 1j
    ordered = Q @ Q @ dQ + Q @ dQ @ Q + dQ @ Q @ Q
    direct = comm(Q @ Q @ Q, G) / 1j
    naive = 3 * Q @ Q @ dQ
    scale = np.linalg.norm(direct)
    leibniz_error = np.linalg.norm(ordered - direct) / scale
    naive_error = np.linalg.norm(naive - direct) / scale
    eps = 1e-5
    plus, minus = Q + eps * dQ, Q - eps * dQ
    finite_difference = (plus @ plus @ plus - minus @ minus @ minus) / (2 * eps)
    fd_error = np.linalg.norm(finite_difference - ordered) / scale
    assert leibniz_error < 1e-12 and fd_error < 1e-7
    assert naive_error > .01
    summary = {
        "python": platform.python_version(), "numpy": np.__version__,
        "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
        "gaussian": {"hbar": 1, "mass": 1, "lambda": 1, "time_max": 3,
                     "samples": 121, "covariance_initial": "I/2",
                     "initial_direction": "(1,1)/sqrt(2)",
                     "overlap_displacement": .1,
                     "max_qfi_error": float(np.max(abs(rows[:, 2] - 2))),
                     "max_symplectic_error": float(symplectic_error),
                     "max_covariance_determinant_error": float(purity_error)},
        "polynomial": {"dimension": 5, "seed": 20261124, "fd_step": eps,
                       "relative_leibniz_error": float(leibniz_error),
                       "relative_fd_error": float(fd_error),
                       "relative_naive_error": float(naive_error)},
        "scope": "Gaussian result exact analytically; matrix test checks Leibniz rule, not finite-dimensional CCR. No chaos or asymptotic inference."
    }
    np.savetxt(paths[1], rows, delimiter=",",
               header="t,centroid_tangent_norm,displacement_qfi,squared_overlap", comments="")
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    axes[0].semilogy(rows[:, 0], rows[:, 1], label="Centroid tangent")
    axes[0].semilogy(rows[:, 0], np.exp(rows[:, 0]), "--", label="exp(t)")
    axes[0].set(xlabel="t", ylabel="Euclidean tangent norm")
    axes[0].legend()
    axes[1].plot(rows[:, 0], rows[:, 2], label="Displacement QFI")
    axes[1].plot(rows[:, 0], rows[:, 3], label="Squared overlap, displacement 0.1")
    axes[1].set(xlabel="t", ylim=(0, 2.3), ylabel="Unitary-invariant quantities")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(paths[2], dpi=160)
    plt.close(fig)
    paths[0].write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
