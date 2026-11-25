"""Two quartic-SYK realizations: instantaneous OE PGA, energy spread, time QFI."""
from itertools import combinations
from pathlib import Path
import json
import platform

import numpy as np
import scipy
from four_model_setup import majorana
from canonical_instantaneous_obsen_coefficient import run, objective


def main():
    stem = Path(__file__).resolve().parents[1] / "results/syk_N16_instantaneous_pga_pair"
    for suffix in (".json", ".npz", ".csv"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    n, d, m = 16, 256, 10
    optimization, A, alpha = run(d, m, 4, 300, 20261123)
    # Remove roundoff PSD violations by a common normalization.
    values, vectors = np.linalg.eigh(A)
    A = (vectors * np.maximum(values, 0)[:, None, :]) @ vectors.conj().transpose(0, 2, 1)
    values, vectors = np.linalg.eigh(A.sum(axis=0))
    inverse = (vectors * values**(-.5)) @ vectors.conj().T
    A = inverse[None] @ A @ inverse[None]
    coefficient = objective(A, alpha, d)[0]
    rows, saved = [], dict(A=A, alpha=alpha)
    for seed in (20261113, 20261123):
        rng = np.random.default_rng(seed)
        H = np.zeros((d, d), complex)
        for sites in combinations(range(n), 4):
            coupling = rng.normal(scale=np.sqrt(6 / n**3))
            for col in range(d):
                row, amplitude = col, coupling
                for site in reversed(sites):
                    row, factor = majorana(site, row)
                    amplitude *= factor
                H[row, col] += amplitude
        assert np.max(abs(H-H.conj().T)) < 1e-12
        psi = np.eye(d, dtype=complex)[:, 0]
        mean = np.vdot(psi, H @ psi).real
        tangent = -1j * (H @ psi - mean * psi)
        delta = np.linalg.norm(tangent)
        # Canonical dotrho has negative off-diagonals.
        B = np.column_stack((psi, -tangent/delta))
        complement = np.eye(d)-B @ B.conj().T
        effects = B[None] @ A @ B.conj().T[None] + alpha[:, None, None]*complement
        p = np.einsum("a,iab,b->i", psi.conj(), effects, psi).real
        dp = 2*np.einsum("a,iab,b->i", psi.conj(), effects, tangent).real
        volumes = np.trace(effects, axis1=1, axis2=2).real
        rate = float(-dp @ np.log(p/volumes))
        completeness = float(np.linalg.norm(effects.sum(axis=0)-np.eye(d)))
        minimum = float(np.linalg.eigvalsh(effects).min())
        assert completeness < 1e-10 and minimum > -1e-12
        assert p.min() > 0 and abs(p.sum()-1) < 1e-10
        assert abs(rate-coefficient*delta) < 1e-10
        energies, vectors = np.linalg.eigh(H)
        times = np.linspace(0, 40, 81)
        states = vectors @ (vectors[0].conj()[:, None]*np.exp(-1j*energies[:, None]*times))
        acted = H @ states
        means = np.sum(states.conj()*acted, axis=0)
        qfi = 4*np.sum(abs(acted-states*means)**2, axis=0)
        assert np.max(abs(np.sum(abs(states)**2, axis=0)-1)) < 1e-11
        qfi_error = float(np.max(abs(qfi-4*delta**2)))
        assert qfi_error < 1e-10
        def entropy(t):
            state = vectors @ (vectors[0].conj()*np.exp(-1j*energies*t))
            probs = np.einsum("a,iab,b->i", state.conj(), effects, state).real
            assert probs.min() > 0
            return float(-probs @ np.log(probs/volumes))
        finite_difference = (entropy(1e-5)-entropy(-1e-5))/2e-5
        assert abs(finite_difference-rate) < 1e-7
        parity = np.array([(-1)**b.bit_count() for b in range(d)])
        assert np.max(abs((parity[:, None]-parity[None, :])*H)) < 1e-12
        assert np.max(abs((parity[None, :, None]-parity[None, None, :])*effects)) < 1e-11
        even = np.flatnonzero(parity == 1)
        ev = np.linalg.eigvalsh(H[np.ix_(even, even)])[16:-16]
        gaps = np.diff(ev)
        gap_ratio = float(np.mean(np.minimum(gaps[1:], gaps[:-1])/np.maximum(gaps[1:], gaps[:-1])))
        rows.append(dict(seed=seed, energy_mean=float(mean), energy_uncertainty=float(delta),
                         physical_time_qfi=float(4*delta**2), attained_rate=rate,
                         rate_over_deltaH=rate/delta, finite_difference_rate=finite_difference,
                         completeness_error=completeness, minimum_effect_eigenvalue=minimum,
                         qfi_constancy_error=qfi_error, even_sector_central_gap_ratio=gap_ratio))
        saved[f"H_{seed}"] = H
        saved[f"basis_{seed}"] = B
        saved[f"qfi_{seed}"] = qfi
    summary = dict(N=n, dimension=d, outcomes=m, J=1, hbar=1,
                   convention="{chi_i,chi_j}=delta_ij; H=sum_{i<j<k<l} J_ijkl chi_i chi_j chi_k chi_l",
                   coupling_variance=6/n**3, initial_state="occupation vacuum", evaluation_time=0,
                   volume_convention="full Fock-space trace, d=256; effects commute with parity",
                   optimization=optimization, feasible_coefficient=coefficient, results=rows,
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   caveat="Attained lower bounds, not certified global maxima. Two finite disorder samples do not establish chaos.")
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    np.savez_compressed(stem.with_suffix(".npz"), times=times, **saved)
    np.savetxt(stem.with_suffix(".csv"),
               [[r[k] for k in ("seed", "energy_uncertainty", "physical_time_qfi", "attained_rate")] for r in rows],
               delimiter=",", header="seed,DeltaH,FQ_time,attained_OE_rate", comments="")
    print(json.dumps({"coefficient": coefficient, "results": rows}, indent=2))


if __name__ == "__main__":
    main()
