"""Animate the five-kick Husimi density while sweeping the kick parameter."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter, writers
import numpy as np
from scipy.linalg import expm
from scipy.special import binom

# Defaults: k=2..10, 50 fps, four seconds; options also reproduce earlier movies.
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--k-start", type=float, default=2)
parser.add_argument("--k-stop", type=float, default=10)
parser.add_argument("--fps", type=int, default=50)
parser.add_argument("--seconds", type=float, default=4)
parser.add_argument("--resolution", choices=("540p", "1080p"), default="1080p")
args = parser.parse_args()
assert args.k_start < args.k_stop and args.fps > 0 and args.seconds > 0
fps, frames, kicks, j = args.fps, round(args.fps * args.seconds), 5, 20
assert frames >= 2
k_values = np.linspace(args.k_start, args.k_stop, frames)
d = 2 * j + 1
m = np.arange(-j, j + 1)
raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy = (raising - raising.T) / (2j)
rotation = expm(-1j * np.pi * Jy / 2)

theta0, phi0 = 1.1, 0.7
psi0 = (np.sqrt(binom(2 * j, j + m)) * np.cos(theta0 / 2) ** (j + m)
        * np.sin(theta0 / 2) ** (j - m) * np.exp(-1j * m * phi0))
assert np.isclose(np.linalg.norm(psi0), 1, atol=1e-12)
assert np.allclose(rotation.conj().T @ rotation, np.eye(d), atol=1e-11)

# q=azimuth, p=cos(polar angle); Q is normalized with respect to dq dp.
high_resolution = args.resolution == "1080p"
q_grid = np.linspace(0, 2 * np.pi, 512 if high_resolution else 256, endpoint=False)
p_grid = np.linspace(-1, 1, 361 if high_resolution else 181)
polar_coefficients = (np.sqrt(binom(2 * j, j + m))[:, None]
                      * np.sqrt((1 + p_grid) / 2)[None, :] ** (j + m[:, None])
                      * np.sqrt((1 - p_grid) / 2)[None, :] ** (j - m[:, None]))
azimuth_phases = np.exp(1j * np.outer(q_grid, m))

# Propagate every k value in one batch.
states = np.tile(psi0, (frames, 1))
twists = np.exp(-1j * k_values[:, None] * m[None, :] ** 2 / (2 * j))
for _ in range(kicks):
    states = twists * (rotation @ states.T).T
assert np.allclose(np.sum(abs(states) ** 2, axis=1), 1, atol=1e-10)

# Verify the first and last frames against an independent, single-state loop.
for index in (0, frames - 1):
    independent = psi0.copy()
    for _ in range(kicks):
        independent = twists[index] * (rotation @ independent)
    assert np.allclose(states[index], independent, atol=1e-11)

def husimi_batch(block):
    """Compute densities in batches rather than storing all 1080p frames."""
    coefficients = (polar_coefficients[:, :, None] * block.T[:, None, :]).reshape(d, -1)
    overlaps = azimuth_phases @ coefficients
    overlaps = overlaps.reshape(len(q_grid), len(p_grid), len(block)).transpose(2, 1, 0)
    return d / (4 * np.pi) * abs(overlaps) ** 2


# First pass finds the global color range and checks normalization.
maximum, max_error = 0.0, 0.0
for first in range(0, frames, 12):
    block = husimi_batch(states[first:first + 12])
    if first == 0:
        first_density = block[0].copy()
    maximum = max(maximum, float(block.max()))
    integrals = 2 * np.pi * np.trapezoid(block.mean(axis=2), p_grid, axis=1)
    max_error = max(max_error, float(np.max(abs(integrals - 1))))
assert max_error < 5e-3
assert writers.is_available("ffmpeg"), "Saving MP4 requires ffmpeg on PATH."

figsize, dpi = ((12.8, 7.2), 150) if high_resolution else ((8, 4.5), 120)
fig, ax = plt.subplots(figsize=figsize)
image = ax.imshow(first_density, extent=[0, 2 * np.pi, -1, 1], origin="lower",
                  aspect="auto", cmap="viridis", vmin=0, vmax=maximum)
ax.set(xlabel=r"$q=\phi$ (azimuth)", ylabel=r"$p=\cos\theta$")
ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi],
              ["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$"])
title = ax.set_title("")
fig.colorbar(image, ax=ax, label=r"$Q(q,p)=(2j+1)|\langle q,p|\psi\rangle|^2/(4\pi)$")
fig.tight_layout()

output = Path(__file__).resolve().parent.parent / "results" / (
    f"kicked_top_husimi_k{args.k_start:g}_to_k{args.k_stop:g}_t{kicks}_{fps}fps"
    f"{'_1080p' if high_resolution else ''}.mp4")
output.parent.mkdir(parents=True, exist_ok=True)
writer = FFMpegWriter(fps=fps, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "18"])
with writer.saving(fig, str(output), dpi=dpi):
    for first in range(0, frames, 12):
        block = husimi_batch(states[first:first + 12])
        for kappa, density in zip(k_values[first:first + len(block)], block):
            image.set_data(density)
            title.set_text(f"Kicked-top Husimi density: j={j}, k={kappa:.2f}, {kicks} kicks")
            writer.grab_frame()
plt.close(fig)
print(f"{output} ({frames} frames at {fps} fps; max grid-integral error {max_error:.6g})")
