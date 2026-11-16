"""Measure in a new random basis at each time, keeping the backaction."""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from scipy.special import rel_entr

# Use the same two preparations as in kicked_top_relative_entropy.py.
j, epsilon, noise, kicks = 10, 0.07, 0.3, 5
d = 2 * j + 1
identity = np.eye(d)
m = np.arange(-j, j + 1, dtype=float)
raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy = (raising - raising.T) / (2j)
Jz = np.diag(m)
psi = expm(-1j * 0.7 * Jz) @ expm(-1j * 1.1 * Jy) @ identity[:, -1]
rho0 = (1 - noise) * np.outer(psi, psi.conj()) + noise * identity / d
V = expm(-1j * epsilon * Jy)
sigma0 = V @ rho0 @ V.conj().T

# Choose the measurements once, without looking at either state or kappa.
# At each time both preparations use exactly the same random basis.
rng = np.random.default_rng(20261116)
bases = []
for _ in range(kicks + 1):
    random_matrix = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    basis, _ = np.linalg.qr(random_matrix)
    bases.append(basis)

fig, ax = plt.subplots()
for kappa in (1.3, 4, 6.0, 10):
    U = expm(-1j * kappa * Jz @ Jz / (2 * j)) @ expm(-1j * np.pi / 2 * Jy)
    rho, sigma = rho0.copy(), sigma0.copy()
    divergences = []

    for t, basis in enumerate(bases):
        # Measure both states in this time's common random basis.
        p = np.diag(basis.conj().T @ rho @ basis).real
        q = np.diag(basis.conj().T @ sigma @ basis).real
        divergences.append(np.sum(rel_entr(p, q)))

        # The outcome is not retained: a projective measurement removes
        # off-diagonal terms in this basis before the next kick.
        rho = basis @ np.diag(p) @ basis.conj().T
        sigma = basis @ np.diag(q) @ basis.conj().T
        if t < kicks:
            rho = U @ rho @ U.conj().T
            sigma = U @ sigma @ U.conj().T

    ax.plot(range(kicks + 1), divergences, label=f"kappa={kappa}")

ax.set(xlabel="Number of kicks", ylabel="Outcome relative entropy (nats)",
       title="Repeated, unrecorded random-basis measurements")
ax.legend()
ax.grid(alpha=0.3)
fig.tight_layout()
output = Path(__file__).resolve().parent.parent / "results" / "kicked_top_random_repeated_measurement.png"
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=160)
print(output)
plt.show()
