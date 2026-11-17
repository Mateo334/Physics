"""Preparations and one-step unitaries shared by the four-model comparisons."""

from itertools import combinations
from math import factorial

import numpy as np
from scipy.linalg import expm

COUNT = 100
SEED = 20261116
TIMES = np.arange(6)  # Initial value and five steps.


def spectral_unitary(hamiltonian):
    """Set the rms eigenenergy to one; discard the irrelevant mean phase."""
    assert np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-11)
    energies, vectors = np.linalg.eigh(hamiltonian)
    width = np.std(energies)
    assert width > 1e-12
    phases = np.exp(-1j * (energies - energies.mean()) / width)
    return (vectors * phases) @ vectors.conj().T


def spin_models():
    """Spin-j coherent states and four kicked-top steps in the fixed Jz basis."""
    j = 20
    m = np.arange(-j, j + 1, dtype=float)
    raising = np.diag(np.sqrt(j * (j + 1) - m[:-1] * (m[:-1] + 1)), -1)
    jy = (raising - raising.T) / (2j)
    values, vectors = np.linalg.eigh(jy)
    north = vectors.conj().T[:, -1]
    rng = np.random.default_rng(SEED + 1)
    states = []
    for _ in range(COUNT):
        theta = np.arccos(rng.uniform(-1, 1))
        phi = rng.uniform(0, 2 * np.pi)
        states.append(np.exp(-1j * phi * m) *
                      (vectors @ (np.exp(-1j * theta * values) * north)))
    states = np.column_stack(states)
    check_states(states)
    rotation = expm(-1j * np.pi * jy / 2)
    assert np.allclose(rotation.conj().T @ rotation, np.eye(len(m)), atol=1e-11)
    steps = {}
    for kappa in (2.0, 4.0, 5.0, 6.0):
        twist = np.exp(-1j * kappa * m**2 / (2 * j))
        steps[f"Kicked top (k={kappa:g})"] = twist[:, None] * rotation
    return states, steps


def majorana(index, bits):
    site = index // 2
    occupied = (bits >> site) & 1
    parity = -1 if (bits & ((1 << site) - 1)).bit_count() % 2 else 1
    factor = parity / np.sqrt(2.0)
    if index % 2:
        factor *= -1j if occupied else 1j
    return bits ^ (1 << site), factor


def syk_model():
    """N=16, q=4 Majorana SYK, even parity; four-mode Slater states."""
    orbitals = 8
    basis = [bits for bits in range(1 << orbitals) if bits.bit_count() % 2 == 0]
    index = {bits: row for row, bits in enumerate(basis)}
    for first in range(2 * orbitals):
        for second in range(first, 2 * orbitals):
            for bits in range(1 << orbitals):
                middle, a = majorana(second, bits)
                _, b = majorana(first, middle)
                middle, c = majorana(first, bits)
                _, e = majorana(second, middle)
                assert abs(a * b + c * e - (1 if first == second else 0)) < 1e-12

    rng = np.random.default_rng(SEED + 2)
    hamiltonian = np.zeros((len(basis), len(basis)), dtype=complex)
    for sites in combinations(range(2 * orbitals), 4):
        coupling = rng.normal(scale=np.sqrt(6 / (2 * orbitals)**3))
        for column, bits in enumerate(basis):
            amplitude, final_bits = coupling, bits
            for site in reversed(sites):
                final_bits, factor = majorana(site, final_bits)
                amplitude *= factor
            hamiltonian[index[final_bits], column] += amplitude

    rng = np.random.default_rng(SEED + 3)
    states = np.zeros((len(basis), COUNT), dtype=complex)
    for n in range(COUNT):
        random_matrix = rng.normal(size=(orbitals, 4)) + 1j * rng.normal(size=(orbitals, 4))
        occupied_orbitals, _ = np.linalg.qr(random_matrix)
        for subset in combinations(range(orbitals), 4):
            bits = sum(1 << site for site in subset)
            states[index[bits], n] = np.linalg.det(occupied_orbitals[list(subset), :])
    check_states(states)
    return states, spectral_unitary(hamiltonian)


def bose_hubbard_model():
    """Three-site, eight-boson ring with SU(3) coherent preparations."""
    bosons, sites, hopping, interaction = 8, 3, 1.0, 0.8
    basis = [(a, b, bosons - a - b) for a in range(bosons + 1)
             for b in range(bosons + 1 - a)]
    index = {occupation: row for row, occupation in enumerate(basis)}
    hamiltonian = np.zeros((len(basis), len(basis)), dtype=complex)
    for column, occupation in enumerate(basis):
        hamiltonian[column, column] += interaction * sum(n * (n - 1) for n in occupation) / 2
        for left, right in ((0, 1), (1, 2), (2, 0)):
            for target, source in ((left, right), (right, left)):
                if occupation[source] == 0:
                    continue
                new = list(occupation)
                new[target] += 1
                new[source] -= 1
                hamiltonian[index[tuple(new)], column] -= hopping * np.sqrt(
                    (occupation[target] + 1) * occupation[source])

    rng = np.random.default_rng(SEED + 4)
    states = np.zeros((len(basis), COUNT), dtype=complex)
    for n in range(COUNT):
        z = rng.normal(size=sites) + 1j * rng.normal(size=sites)
        z /= np.linalg.norm(z)
        for row, occupation in enumerate(basis):
            coefficient = np.sqrt(factorial(bosons) / np.prod([factorial(x) for x in occupation]))
            states[row, n] = coefficient * np.prod(z ** np.array(occupation))
    check_states(states)
    return states, spectral_unitary(hamiltonian)


def check_states(states):
    assert np.allclose(np.sum(abs(states)**2, axis=0), 1, atol=1e-10)


def models():
    """Ordered labels, initial-state columns and unitaries for all seven curves."""
    spin_states, spin_steps = spin_models()
    result = {label: (spin_states, step) for label, step in spin_steps.items()}
    result["Majorana SYK (even parity)"] = syk_model()
    result["Bose--Hubbard (3 sites)"] = bose_hubbard_model()
    rng = np.random.default_rng(SEED + 5)
    random_matrix = rng.normal(size=(len(spin_states), len(spin_states)))
    goe_h = (random_matrix + random_matrix.T) / np.sqrt(2 * len(spin_states))
    result["Random-matrix Hamiltonian (GOE)"] = spin_states, spectral_unitary(goe_h)
    return result
