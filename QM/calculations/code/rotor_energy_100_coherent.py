"""Mean post-first-kick energy uncertainty times ln(2), resolved rotor ensemble."""
from pathlib import Path
import json
import platform
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def uncertainties(d, centers, strengths):
    hbar = 2*np.pi/64
    sigma = np.sqrt(hbar/2)
    q = 2*np.pi*np.arange(d)/d
    p = hbar*np.fft.fftfreq(d, d=1/d)
    energy = p*p/2
    wave = sum(np.exp(-(q[:, None]-centers[None, :]+2*np.pi*s)**2/(4*sigma**2))
               for s in range(-3, 4)).astype(complex)
    wave /= np.linalg.norm(wave, axis=0)
    momentum = np.fft.fft(wave, axis=0, norm="ortho")
    initial_prob = abs(momentum)**2
    initial_mean = energy @ initial_prob
    initial_delta = np.sqrt(np.sum(initial_prob*(energy[:, None]-initial_mean)**2, axis=0))
    assert np.max(abs(initial_delta-hbar/(2*np.sqrt(2)))) < 1e-11
    pre_kick = np.fft.ifft(np.exp(-1j*energy[:, None]/hbar)*momentum, axis=0, norm="ortho")
    result, edges = [], []
    for k in strengths:
        state = np.fft.fft(np.exp(-1j*k*np.cos(q[:, None])/hbar)*pre_kick,
                           axis=0, norm="ortho")
        prob = abs(state)**2
        assert np.max(abs(prob.sum(axis=0)-1)) < 1e-11
        mean = energy @ prob
        delta = np.sqrt(np.sum(prob*(energy[:, None]-mean)**2, axis=0))
        result.append(delta)
        edges.append(float(prob[abs(p)>.8*max(abs(p))].sum(axis=0).max()))
    return np.array(result), max(edges)


def main():
    stem = Path(__file__).resolve().parents[1]/"results/rotor_d512_100coherent_energy_ln2_K4to12"
    for suffix in (".png", ".pdf", ".npz", ".csv", ".json"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    seed = 20261123
    centers = np.random.default_rng(seed).uniform(0, 2*np.pi, 100)
    strengths = np.linspace(4, 12, 401)
    delta, edge = uncertainties(512, centers, strengths)
    reference, _ = uncertainties(1024, centers, strengths)
    relative = float(np.max(abs(delta-reference)/reference))
    assert relative < 1e-6
    values = delta*np.log(2)
    mean = values.mean(axis=1)
    sem = values.std(axis=1, ddof=1)/10
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(strengths, mean, label="Mean over 100 coherent packets")
    ax.fill_between(strengths, mean-sem, mean+sem, alpha=.2, label="±1 standard error")
    ax.set(xlabel="Kick strength K", ylabel="Mean energy uncertainty × ln(2)",
           title="Kicked rotor: immediately after first kick\n"
                 "d=512, ℏ_eff=2π/64; random position centers, initial p₀=0",
           xlim=(4, 12))
    ax.grid(alpha=.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    np.savez_compressed(stem.with_suffix(".npz"), K=strengths, centers=centers,
                        deltaH=delta, deltaH_d1024=reference, mean=mean, sem=sem)
    np.savetxt(stem.with_suffix(".csv"), np.column_stack((strengths, mean, sem)),
               delimiter=",", header="K,mean_DeltaH_ln2,standard_error", comments="")
    summary = dict(d=512, cutoff_check_d=1024, state_count=100, seed=seed,
                   K_interval=[4, 12], K_step=.02, hbar=2*np.pi/64,
                   initial_sigma_q=float(np.sqrt(np.pi/64)),
                   ensemble="Uniform independent q0 in [0,2pi), fixed p0=0; same packets for every K",
                   protocol="Prepare Gaussian; free one time unit; kick; DeltaH for H_phys=p^2/2 at t=1+",
                   maximum_relative_cutoff_error=relative, maximum_outer20percent_probability=edge,
                   mean_increases_at_every_sample=bool(np.all(np.diff(mean)>0)),
                   selected_values=[dict(K=float(strengths[i]), mean=float(mean[i]), sem=float(sem[i]))
                                    for i in (0, 100, 200, 300, 400)],
                   python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__,
                   caveat="Mean of individual pure-state energy uncertainties, not uncertainty of the ensemble mixture; no monotonicity theorem or chaos inference.")
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
