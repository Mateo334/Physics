"""Plot fixed-Jz outcome relative entropy after evolving both states."""

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm, logm
from scipy.special import rel_entr

# Same preparation and perturbation as kicked_top_frechet.py.
j, epsilon, noise, kicks = 10, 0.07, 0.3, 200
theta, phi = 1.1, 0.7
d = 2 * j + 1
eye = np.eye(d)
m = np.arange(-j, j + 1, dtype=float)
plus = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy, Jz = (plus - plus.T) / (2j), np.diag(m)
psi = expm(-1j * phi * Jz) @ expm(-1j * theta * Jy) @ eye[:, -1]
rho0 = (1 - noise) * np.outer(psi, psi.conj()) + noise * eye / d
V = expm(-1j * epsilon * Jy)
sigma0 = V @ rho0 @ V.conj().T

fig, ax = plt.subplots()
for kappa in (4, 6.0):
    # Rightmost rotation acts first.
    U = expm(-1j * kappa * (Jz @ Jz) / (2 * j)) @ expm(-1j * np.pi / 2 * Jy)
    assert np.allclose(U.conj().T @ U, eye, atol=1e-12, rtol=0)
    rho, sigma = rho0.copy(), sigma0.copy()
    values = []
    for _ in range(kicks + 1):
        for state in (rho, sigma):
            assert np.allclose(state, state.conj().T, atol=1e-11, rtol=0)
            assert abs(np.trace(state) - 1) < 1e-11
            assert np.linalg.eigvalsh(state).min() > 0
        # Final-only Jz measurement: no collapse is fed into the evolution.
        p, q = np.diag(rho).real, np.diag(sigma).real
        value = np.sum(rel_entr(p, q))  # D(reference || perturbed), in nats
        assert np.isfinite(value) and value >= -1e-12
        assert np.sum(rel_entr(p, p)) == 0
        values.append(value)
        rho = U @ rho @ U.conj().T
        sigma = U @ sigma @ U.conj().T
    ax.plot(range(kicks + 1), values, label=f"kappa={kappa}")

ax.set(xlabel="Number of kicks", ylabel="Relative entropy (nats)",
       title="Fixed Jz measurement: reference vs perturbed state")
ax.legend()
ax.grid(alpha=0.3)
bound = np.trace(rho0 @ (logm(rho0) - logm(sigma0))).real
ax.axhline(bound, linestyle="--", label="Data-processing bound")
fig.tight_layout()
plt.show()
