"""Final single-Majorana probe QFI, full Fock space, no applied probe."""
from itertools import combinations
from pathlib import Path
import json
import platform

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from four_model_setup import majorana


def main():
    n, seed = 16, 20261113
    d = 2 ** (n // 2)
    stem = Path(__file__).resolve().parents[1] / "results/syk_N16_majorana_probe_qfi"
    for suffix in (".png", ".csv", ".json"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    h = np.zeros((d, d), complex)
    rng = np.random.default_rng(seed)
    for sites in combinations(range(n), 4):
        coupling = rng.normal(scale=np.sqrt(6 / n**3))
        for col in range(d):
            row, amplitude = col, coupling
            for site in reversed(sites):
                row, factor = majorana(site, row)
                amplitude *= factor
            h[row, col] += amplitude
    chi = np.zeros_like(h)
    for col in range(d):
        row, factor = majorana(0, col)
        chi[row, col] = factor
    parity = np.array([(-1) ** b.bit_count() for b in range(d)])
    assert np.max(abs(h - h.conj().T)) < 1e-12
    assert np.max(abs(chi - chi.conj().T)) < 1e-12
    assert np.max(abs(chi @ chi - np.eye(d) / 2)) < 1e-12
    assert np.max(abs((parity[:, None] - parity[None, :]) * h)) < 1e-12
    assert np.max(abs((parity[:, None] + parity[None, :]) * chi)) < 1e-12
    energies, vectors = np.linalg.eigh(h)
    times = np.linspace(0, 40, 801)
    states = vectors @ (vectors[0].conj()[:, None] *
                        np.exp(-1j * energies[:, None] * times))
    norm_error = float(np.max(abs(np.sum(abs(states)**2, axis=0) - 1)))
    acted = chi @ states
    means = np.sum(states.conj() * acted, axis=0)
    tangent = acted - states * means
    qfi = 4 * np.sum(abs(tangent)**2, axis=0)
    error = float(np.max(abs(qfi - 2)))
    assert norm_error < 1e-11
    assert np.max(abs(means)) < 1e-11
    assert np.max(abs(states[parity == -1])) < 1e-11
    assert error < 1e-10
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(times, qfi, lw=2, label=r"$F_Q(t)=4\,\mathrm{Var}_{\psi(t)}(\chi_1)$")
    ax.set(xlabel=r"$Jt$", ylabel="Final-probe QFI",
           title="Quartic SYK: single-Majorana probe\nN=16, occupation vacuum, one disorder realization",
           ylim=(0, 2.5), xlim=(0, 40))
    ax.legend()
    ax.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    plt.close(fig)
    np.savetxt(stem.with_suffix(".csv"), np.column_stack((times, qfi)),
               delimiter=",", header="Jt,FQ_chi1", comments="")
    summary = dict(N=n, dimension=d, q=4, J=1, hbar=1, seed=seed,
                   coupling_variance=6/n**3, initial_state="occupation vacuum",
                   convention="{chi_i,chi_j}=delta_ij",
                   protocol="exp(-i theta chi_1) after exp(-i H t); probe not applied",
                   time_points=len(times), time_interval=[0, 40],
                   norm_error=norm_error, max_error_from_exact_2=error,
                   python=platform.python_version(), numpy=np.__version__,
                   matplotlib=matplotlib.__version__, dtype="complex128",
                   caveat="Odd probe is formal unless parity exchange/reference is allowed.")
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
