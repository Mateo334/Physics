"""Compare final-only and unrecorded-measurement protocols for one fixed POVM."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.special import rel_entr

j, epsilon, noise, kicks = 10, 0.07, 0.3, 20
d, rank, outcomes, number_of_states = 2 * j + 1, 7, 4, 100
identity = np.eye(d)
m = np.arange(-j, j + 1, dtype=float)
raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy = (raising - raising.T) / (2j)
Jz = np.diag(m)
V = expm(-1j * epsilon * Jy)

# Use the same 100 reference/perturbed state pairs in both protocols.
state_rng = np.random.default_rng(20261117)
jy_values, jy_vectors = np.linalg.eigh(Jy)
north = jy_vectors.conj().T @ identity[:, -1]
reference, perturbed = [], []
for _ in range(number_of_states):
    theta = np.arccos(state_rng.uniform(-1, 1))
    phi = state_rng.uniform(0, 2 * np.pi)
    psi = np.exp(-1j * phi * m) * (jy_vectors @ (np.exp(-1j * theta * jy_values) * north))
    rho = (1 - noise) * np.outer(psi, psi.conj()) + noise * identity / d
    reference.append(rho)
    perturbed.append(V @ rho @ V.conj().T)
initial_states = np.array([reference, perturbed])

# Choose one random rank-7 POVM, then hold it fixed for all states and kicks.
rng = np.random.default_rng(20261116)
projectors = []
for _ in range(outcomes):
    random_matrix = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    basis, _ = np.linalg.qr(random_matrix)
    subspace = basis[:, :rank]
    projectors.append(subspace @ subspace.conj().T)

total = sum(projectors)
eigenvalues, eigenvectors = np.linalg.eigh(total)
assert eigenvalues.min() > 1e-12
inverse_root = (eigenvectors / np.sqrt(eigenvalues)) @ eigenvectors.conj().T
effects, kraus = [], []
for projector in projectors:
    effect = inverse_root @ projector @ inverse_root
    values, vectors = np.linalg.eigh(effect)
    effects.append(effect)
    kraus.append((vectors * np.sqrt(np.clip(values, 0, None))) @ vectors.conj().T)
assert np.allclose(sum(effects), identity, atol=1e-10)


def outcome_divergence(rho, sigma):
    """Mean of the 100 individual outcome KL divergences, not KL of the mean."""
    p = np.einsum("aij,sji->sa", effects, rho).real
    q = np.einsum("aij,sji->sa", effects, sigma).real
    assert np.all(p > 0) and np.all(q > 0)
    for probabilities in (p, q):
        assert np.allclose(probabilities.sum(axis=1), 1, atol=1e-10)
    return np.sum(rel_entr(p, q), axis=1).mean()


def unrecorded_measurement(states):
    """Apply the same Lüders instrument and discard its outcome."""
    return sum(K @ states @ K.conj().T for K in kraus)


fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for kappa in (1.3, 4, 6.0, 10):
    U = expm(-1j * kappa * Jz @ Jz / (2 * j)) @ expm(-1j * np.pi / 2 * Jy)
    free, measured = initial_states.copy(), initial_states.copy()
    without_backaction, with_backaction = [], []

    for _ in range(kicks):
        # Evolve both state pairs, then evaluate the same POVM before backaction.
        free = U @ free @ U.conj().T
        measured = U @ measured @ U.conj().T
        without_backaction.append(outcome_divergence(*free))
        with_backaction.append(outcome_divergence(*measured))
        # Only the measured pair feeds discarded-outcome states into the next kick.
        measured = unrecorded_measurement(measured)

    assert np.isclose(without_backaction[0], with_backaction[0], atol=1e-12)
    axes[0].plot(range(1, kicks + 1), without_backaction, label=f"kappa={kappa}")
    axes[1].plot(range(1, kicks + 1), with_backaction, label=f"kappa={kappa}")

axes[0].set_title("Fixed POVM, no backaction")
axes[1].set_title("Fixed POVM, unrecorded backaction")
for ax in axes:
    ax.set_xlabel("Number of kicks")
    ax.grid(alpha=0.3)
    ax.legend()
axes[0].set_ylabel("Mean single-time outcome KL (nats)")
fig.suptitle(f"Same {number_of_states} state pairs and same rank-{rank} POVM")
fig.tight_layout()
output = Path(__file__).resolve().parent.parent / "results" / "kicked_top_fixed_povm_backaction_comparison_100_states.png"
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=160)
plt.close(fig)
print(output)
