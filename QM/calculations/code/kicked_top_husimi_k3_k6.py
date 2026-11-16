"""Husimi phase-space densities for the kicked top at k=3 and k=6."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.special import binom

# Keep the same spin size and five-kick horizon as the recent comparison.
j, kicks = 20, 5
theta0, phi0 = 1.1, 0.7
d = 2 * j + 1
m = np.arange(-j, j + 1)
raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy = (raising - raising.T) / (2j)
Jz = np.diag(m)
rotation = expm(-1j * np.pi * Jy / 2)

# Coordinates on the spin sphere: q is azimuth, p = cos(theta).
# The coherent-state resolution of identity makes Q(q,p) a probability
# density with respect to dq dp, not a density on a flat Cartesian plane.
q_grid = np.linspace(0, 2 * np.pi, 256, endpoint=False)
p_grid = np.linspace(-1, 1, 181)
cos_half = np.sqrt((1 + p_grid) / 2)
sin_half = np.sqrt((1 - p_grid) / 2)
coherent_amplitudes = np.sqrt(binom(2 * j, j + m))[:, None] * (
    cos_half[None, :] ** (j + m[:, None])) * (
    sin_half[None, :] ** (j - m[:, None]))
azimuth_phases = np.exp(1j * np.outer(q_grid, m))


def coherent(theta, phi):
    coefficients = np.sqrt(binom(2 * j, j + m))
    return coefficients * np.cos(theta / 2) ** (j + m) * np.sin(theta / 2) ** (j - m) * np.exp(-1j * m * phi)


psi0 = coherent(theta0, phi0)
reference = expm(-1j * phi0 * Jz) @ expm(-1j * theta0 * Jy) @ np.eye(d)[:, -1]
assert np.allclose(psi0, reference, atol=1e-11)
assert np.isclose(np.linalg.norm(psi0), 1, atol=1e-12)
assert np.allclose(rotation.conj().T @ rotation, np.eye(d), atol=1e-11)


def husimi(psi):
    overlaps = azimuth_phases @ (coherent_amplitudes * psi[:, None])
    density = d / (4 * np.pi) * abs(overlaps.T) ** 2
    normalization = 2 * np.pi * np.trapezoid(density.mean(axis=1), p_grid)
    assert abs(normalization - 1) < 5e-3
    return density, normalization


results = []
for kappa in (3, 6):
    twist = np.exp(-1j * kappa * m**2 / (2 * j))
    psi = psi0.copy()
    for _ in range(kicks):
        psi = twist * (rotation @ psi)
    assert np.isclose(np.linalg.norm(psi), 1, atol=1e-11)
    density, normalization = husimi(psi)
    results.append((kappa, density, normalization))

maximum = max(density.max() for _, density, _ in results)
output_directory = Path(__file__).resolve().parent.parent / "results"
output_directory.mkdir(parents=True, exist_ok=True)
for kappa, density, normalization in results:
    fig, ax = plt.subplots(figsize=(8, 4.5))
    image = ax.imshow(density, extent=[0, 2 * np.pi, -1, 1], origin="lower",
                      aspect="auto", cmap="viridis", vmin=0, vmax=maximum)
    ax.set(xlabel=r"$q=\phi$ (azimuth)", ylabel=r"$p=\cos\theta$",
           title=f"Kicked-top Husimi density: j={j}, k={kappa}, {kicks} kicks")
    ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
                  ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
    fig.colorbar(image, ax=ax, label=r"$Q(q,p)=(2j+1)|\langle q,p|\psi\rangle|^2/(4\pi)$")
    fig.tight_layout()
    output = output_directory / f"kicked_top_husimi_k{kappa}_t{kicks}.png"
    fig.savefig(output, dpi=160)
    plt.close(fig)
    print(f"{output} (integral approximately {normalization:.6f})")
