"""One plot of mean observational entropy for four finite quantum models."""

from itertools import combinations
from math import factorial
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.special import xlogy

COUNT, TIMES = 100, np.arange(6)  # Initial value and five kicks.
SEED = 20261116
SQRT2 = np.sqrt(2.0)


def entropy_curve(states, evolve, dimension):
    """Fixed rank-one basis PVM: S_obs is the Shannon entropy of |psi_i|^2."""
    result = []
    for t in TIMES:
        probabilities = abs(states) ** 2
        assert np.allclose(probabilities.sum(axis=0), 1, atol=1e-10)
        values = -np.sum(xlogy(probabilities, probabilities), axis=0)
        assert np.all((values >= -1e-10) & (values <= np.log(dimension) + 1e-10))
        result.append(values.mean() / np.log(dimension))
        if t < TIMES[-1]:
            states = evolve(states, t)
    return np.array(result)


def spectral_evolution(hamiltonian, states):
    """Normalize rms eigenenergy to one; return an exact unitary time step."""
    assert np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-11)
    energies, vectors = np.linalg.eigh(hamiltonian)
    width = np.std(energies)
    assert width > 1e-12
    # Removing the mean changes only an irrelevant global phase.
    phases = np.exp(-1j * (energies - energies.mean()) / width)
    def step(states, _):
        return vectors @ (phases[:, None] * (vectors.conj().T @ states))
    return step


# The same spin-coherent preparations are used for both kicked tops and GOE.
j = 20
d_spin = 2 * j + 1
m = np.arange(-j, j + 1, dtype=float)
raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
Jy = (raising - raising.T) / (2j)
jy_values, jy_vectors = np.linalg.eigh(Jy)
north = jy_vectors.conj().T[:, -1]
spin_rng = np.random.default_rng(SEED + 1)
spin_states = []
for _ in range(COUNT):
    theta = np.arccos(spin_rng.uniform(-1, 1))
    phi = spin_rng.uniform(0, 2 * np.pi)
    psi = np.exp(-1j * phi * m) * (jy_vectors @ (np.exp(-1j * theta * jy_values) * north))
    spin_states.append(psi)
spin_states = np.column_stack(spin_states)
rotation = expm(-1j * np.pi * Jy / 2)
assert np.allclose(rotation.conj().T @ rotation, np.eye(d_spin), atol=1e-11)

curves = {}
for kappa in (2.0, 4.0, 5.0, 6.0):
    twist = np.exp(-1j * kappa * m**2 / (2 * j))
    def kick(states, _, twist=twist):
        return twist[:, None] * (rotation @ states)
    curves[f"Kicked top (k={kappa:g})"] = entropy_curve(spin_states.copy(), kick, d_spin)


# q=4 Majorana SYK with 16 Majoranas, restricted to its even-parity sector.
# A Slater determinant of four complex modes is a physical fermionic
# generalized coherent state, unlike a Grassmann-valued coherent state.
orbitals = 8
syk_basis = [b for b in range(1 << orbitals) if b.bit_count() % 2 == 0]
syk_index = {bits: index for index, bits in enumerate(syk_basis)}
d_syk = len(syk_basis)


def majorana(index, bits):
    site = index // 2
    occupied = (bits >> site) & 1
    parity = -1 if (bits & ((1 << site) - 1)).bit_count() % 2 else 1
    factor = parity / SQRT2
    if index % 2:
        factor *= -1j if occupied else 1j
    return bits ^ (1 << site), factor


# Check {chi_i,chi_j}=delta_ij on every occupation-basis vector.
for first in range(2 * orbitals):
    for second in range(first, 2 * orbitals):
        for bits in range(1 << orbitals):
            middle, a = majorana(second, bits)
            _, b = majorana(first, middle)
            middle, c = majorana(first, bits)
            _, e = majorana(second, middle)
            assert abs(a * b + c * e - (1 if first == second else 0)) < 1e-12


syk_rng = np.random.default_rng(SEED + 2)
syk_h = np.zeros((d_syk, d_syk), dtype=complex)
for sites in combinations(range(2 * orbitals), 4):
    coupling = syk_rng.normal(scale=np.sqrt(6 / (2 * orbitals)**3))
    for column, bits in enumerate(syk_basis):
        amplitude, final_bits = coupling, bits
        for site in reversed(sites):
            final_bits, factor = majorana(site, final_bits)
            amplitude *= factor
        syk_h[syk_index[final_bits], column] += amplitude

slater_rng = np.random.default_rng(SEED + 3)
slater_basis = list(combinations(range(orbitals), 4))
slater_states = np.zeros((d_syk, COUNT), dtype=complex)
for n in range(COUNT):
    random_matrix = slater_rng.normal(size=(orbitals, 4)) + 1j * slater_rng.normal(size=(orbitals, 4))
    occupied_orbitals, _ = np.linalg.qr(random_matrix)
    for subset in slater_basis:
        bits = sum(1 << site for site in subset)
        slater_states[syk_index[bits], n] = np.linalg.det(occupied_orbitals[list(subset), :])
assert np.allclose(np.sum(abs(slater_states) ** 2, axis=0), 1, atol=1e-10)
curves["Majorana SYK (even parity)"] = entropy_curve(
    slater_states, spectral_evolution(syk_h, slater_states), d_syk)


# Three-site Bose--Hubbard ring at fixed boson number N=8. An SU(3)
# coherent state puts all eight bosons into the same random single-particle mode.
bosons, sites, hopping, interaction = 8, 3, 1.0, 0.8
boson_basis = [(a, b, bosons - a - b) for a in range(bosons + 1)
               for b in range(bosons + 1 - a)]
boson_index = {occupation: index for index, occupation in enumerate(boson_basis)}
d_bh = len(boson_basis)
bh_h = np.zeros((d_bh, d_bh), dtype=complex)
for column, occupation in enumerate(boson_basis):
    bh_h[column, column] += interaction * sum(n * (n - 1) for n in occupation) / 2
    for left, right in [(0, 1), (1, 2), (2, 0)]:
        for target, source in [(left, right), (right, left)]:
            if occupation[source] == 0:
                continue
            new = list(occupation)
            new[target] += 1
            new[source] -= 1
            bh_h[boson_index[tuple(new)], column] -= hopping * np.sqrt(
                (occupation[target] + 1) * occupation[source])

bh_rng = np.random.default_rng(SEED + 4)
bh_states = np.zeros((d_bh, COUNT), dtype=complex)
for n in range(COUNT):
    z = bh_rng.normal(size=sites) + 1j * bh_rng.normal(size=sites)
    z /= np.linalg.norm(z)
    for index, occupation in enumerate(boson_basis):
        coefficient = np.sqrt(factorial(bosons) / np.prod([factorial(x) for x in occupation]))
        bh_states[index, n] = coefficient * np.prod(z ** np.array(occupation))
assert np.allclose(np.sum(abs(bh_states) ** 2, axis=0), 1, atol=1e-10)
curves["Bose--Hubbard (3 sites)"] = entropy_curve(
    bh_states, spectral_evolution(bh_h, bh_states), d_bh)


# A single GOE Hamiltonian in the spin-j basis, using the SAME 100 spin
# coherent states as the two kicked tops (only the Hamiltonian changes).
goe_rng = np.random.default_rng(SEED + 5)
random_matrix = goe_rng.normal(size=(d_spin, d_spin))
goe_h = (random_matrix + random_matrix.T) / np.sqrt(2 * d_spin)
curves["Random-matrix Hamiltonian (GOE)"] = entropy_curve(
    spin_states.copy(), spectral_evolution(goe_h, spin_states), d_spin)
assert np.isclose(curves["Random-matrix Hamiltonian (GOE)"][0],
                  curves["Kicked top (k=2)"][0], atol=1e-12)

fig, ax = plt.subplots(figsize=(9, 5))
for label, curve in curves.items():
    ax.plot(TIMES, curve, label=label)
ax.set(xlabel="Kicks / normalized Hamiltonian time (different physical clocks)",
       ylabel="Mean observational entropy / log(sector dimension)",
       title="Four models, kicked-top k=2,4,5,6: 100 coherent states each")
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
output = Path(__file__).resolve().parent.parent / "results" / "obsen_four_models_k2_4_5_6_5kicks_100_states.png"
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=160)
plt.close(fig)
print(output)
print(f"Dimensions: kicked top and GOE={d_spin}, SYK parity={d_syk}, Bose--Hubbard={d_bh}")
