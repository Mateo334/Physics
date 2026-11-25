"""Fixed four-outcome OE curves realizing the saved best attained initial slope."""
from pathlib import Path
import json
import platform

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from canonical_instantaneous_obsen_coefficient import objective


def main():
    out = Path(__file__).resolve().parents[1]/"results"
    stem = out/"obsen_attained_rate_models_d256_m4"
    for suffix in (".png", ".pdf", ".npz", ".json"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    with np.load(out/"obsen_constant_m4_dimension.npz") as data:
        A, alpha = data["A_d256"], data["alpha_d256"]
    coefficient = objective(A, alpha, 256)[0]
    assert np.linalg.eigvalsh(A).min() > -1e-12
    assert np.linalg.norm(A.sum(axis=0)-np.eye(2)) < 1e-12
    assert alpha.min() >= 0 and abs(alpha.sum()-1) < 1e-12
    with np.load(out/"syk_N16_instantaneous_pga_pair.npz") as data:
        H = data["H_20261113"]
        B = data["basis_20261113"]
    models = [("SYK, seed 20261113", H, B[:, 0], B, "0")]
    with np.load(out/"kicked_models_d256_m10_compare.npz") as data:
        j = 127.5
        z = np.arange(-j, j+1)
        raising = np.diag(np.sqrt(j*(j+1)-z[:-1]*(z[:-1]+1)), -1)
        jy = (raising-raising.T)/(2j)
        models.append(("Kicked top, K=6", np.pi*jy, data["top_K6_state"],
                       data["top_K6_frame"], "1.1"))
        p = data["rotor_momentum"]
        models.append(("Kicked rotor, K=10", np.diag(p*p/(2*(2*np.pi/64))),
                       data["rotor_K10_state"], data["rotor_K10_frame"], "1+"))
    volumes = np.trace(A, axis1=1, axis2=2).real+254*alpha
    saved, rows = dict(A=A, alpha=alpha), []
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
    for index, (ax, (name, generator, psi, B, t0)) in enumerate(zip(axes, models)):
        assert np.linalg.norm(B.conj().T @ B-np.eye(2)) < 1e-10
        mean = np.vdot(psi, generator @ psi).real
        tangent = -1j*(generator @ psi-mean*psi)
        spread = np.linalg.norm(tangent)
        assert np.linalg.norm(tangent+spread*B[:, 1]) < 1e-10
        values, vectors = np.linalg.eigh(generator)
        coords = vectors.conj().T @ psi
        def entropy(times):
            states = vectors @ (coords[:, None]*np.exp(-1j*values[:, None]*times))
            assert np.max(abs(np.sum(abs(states)**2, axis=0)-1)) < 1e-11
            projected = B.conj().T @ states
            remaining = 1-np.sum(abs(projected)**2, axis=0)
            assert remaining.min() > -1e-11
            probs = np.einsum("at,iab,bt->it", projected.conj(), A, projected).real
            probs += alpha[:, None]*remaining
            assert probs.min() > 0
            assert np.max(abs(probs.sum(axis=0)-1)) < 1e-11
            return -np.sum(probs*np.log(probs/volumes[:, None]), axis=0)
        times = np.linspace(0, 2/spread, 801)
        # Windows remain inside the same physical generator phase.
        if index == 1:
            assert times[-1] < .4
        if index == 2:
            assert times[-1] < 1
        curve = entropy(times)
        rate = coefficient*spread
        step = 1e-5/spread
        check = entropy(np.array([-step, step]))
        fd = (check[1]-check[0])/(2*step)
        assert abs(fd/rate-1) < 1e-7
        assert curve.min() >= -1e-10 and curve.max() <= np.log(256)+1e-10
        line_times = np.linspace(0, .4/spread, 51)
        ax.plot(times, curve, lw=2, label="Fixed optimized POVM")
        ax.plot(line_times, curve[0]+rate*line_times, "--", color="C1",
                label="Initial-rate tangent")
        ax.axhline(np.log(256), color="gray", ls=":", label="ln(256)")
        ax.set(title=f"{name}\nInitial rate = {rate:.6f} nats/time",
               xlabel=f"Elapsed time from t₀={t0}", xlim=(0, times[-1]), ylim=(0, 5.8))
        ax.grid(alpha=.2)
        rows.append(dict(model=name, t0=t0, generator_uncertainty=float(spread),
                         rate=float(rate), finite_difference_rate=float(fd),
                         initial_entropy=float(curve[0]), final_elapsed_time=float(times[-1]),
                         entropy_range=[float(curve.min()), float(curve.max())]))
        saved[f"time_{index}"], saved[f"entropy_{index}"] = times, curve
        saved[f"frame_{index}"] = B
    axes[0].set_ylabel("Observational entropy (nats)")
    axes[0].legend(fontsize=8, loc="lower right")
    fig.suptitle("d=256, four outcomes: best attained initial rate, then POVM held fixed", fontsize=12)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    np.savez_compressed(stem.with_suffix(".npz"), **saved)
    summary = dict(d=256, outcomes=4, coefficient=coefficient, results=rows,
                   python=platform.python_version(), numpy=np.__version__,
                   matplotlib=matplotlib.__version__,
                   protocol="One different fixed optimized POVM per model; no measurement backaction",
                   window="0 <= elapsed time * Delta(generator) <= 2; all within fixed-generator phases",
                   caveat="Best attained initial slopes, not globally certified maxima or sustained maximal rates. Rotor finite-d cutoff.")
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
