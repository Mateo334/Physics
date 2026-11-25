"""Laptop-friendly four-outcome pure-state OE coefficient versus dimension."""
from pathlib import Path
import json
import platform
import time

import numpy as np
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from canonical_instantaneous_obsen_coefficient import (
    run, objective, random_povm, herm, binary_seed,
)


def main():
    out = Path(__file__).resolve().parents[1]/"results"
    stem = out/"obsen_constant_m4_dimension"
    for suffix in (".json", ".npz", ".csv", ".png", ".pdf"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    start_time = time.perf_counter()
    dimensions = [2, 3, 4, 8, 16, 21, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    rows, saved = [], {}
    for d in dimensions:
        old = out/f"canonical_obsen_coefficient_d{d}_m4"
        if old.with_suffix(".json").exists() and old.with_suffix(".npz").exists():
            summary = json.loads(old.with_suffix(".json").read_text())
            with np.load(old.with_suffix(".npz")) as data:
                A, alpha = data["two_by_two_effects"], data["complement_weights"]
            source = str(old.name)
        else:
            summary, A, alpha = run(d, 4, 3, 200, 20261123+d)
            source = "fresh reduced PGA: binary seed plus two random starts"
        values, vectors = np.linalg.eigh(A)
        A = (vectors*np.maximum(values, 0)[:, None, :]) @ vectors.conj().transpose(0, 2, 1)
        values, vectors = np.linalg.eigh(A.sum(axis=0))
        inverse = (vectors*values**(-.5)) @ vectors.conj().T
        A = inverse[None] @ A @ inverse[None]
        coefficient, _, _, p, dp, V = objective(A, alpha, d)
        assert np.linalg.eigvalsh(A).min() > -1e-12
        assert np.linalg.norm(A.sum(axis=0)-np.eye(2)) < 1e-12
        assert alpha.min() >= 0 and abs(alpha.sum()-1) < 1e-12
        assert abs(p.sum()-1) < 1e-12 and abs(V.sum()-d) < 1e-9
        assert abs(dp.sum()) < 1e-12
        # Directional gradient check in the full HS metric.
        rng = np.random.default_rng(20261124+d)
        X, a = random_povm(4, rng)
        DX = herm(rng.normal(size=X.shape)+1j*rng.normal(size=X.shape))
        DX /= np.linalg.norm(DX)
        da = rng.normal(size=4)
        da /= np.linalg.norm(da)
        _, G, g, *_ = objective(X, a, d)
        analytic = np.einsum("iab,iba->", G, DX).real+(d-2)*g @ da
        h = 1e-7
        fd = (objective(X+h*DX, a+h*da, d)[0]-objective(X-h*DX, a-h*da, d)[0])/(2*h)
        gradient_error = float(abs(fd-analytic)/max(1, abs(analytic)))
        assert gradient_error < 1e-6
        def entropy(t):
            psi = np.array([np.cos(t), -np.sin(t)])
            probs = np.einsum("a,iab,b->i", psi, A, psi).real
            return -probs @ np.log(probs/V)
        fd_rate = (entropy(1e-6)-entropy(-1e-6))/2e-6
        assert abs(fd_rate-coefficient) < 1e-7
        rows.append(dict(d=d, coefficient=coefficient, source=source,
                         gradient_relative_error=gradient_error,
                         entropy_difference_error=float(abs(fd_rate-coefficient)),
                         optimizer=summary))
        saved[f"A_d{d}"], saved[f"alpha_d{d}"] = A, alpha
        print("DIMENSION", d, "COEFFICIENT", coefficient, flush=True)
    coefficients = np.array([r["coefficient"] for r in rows])
    smooth_d = np.unique(np.geomspace(2, dimensions[-1], 240).astype(int))
    binary = np.array([binary_seed(int(d), 4)[2] for d in smooth_d])
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(smooth_d, binary, color="C0", alpha=.65, label="Split binary candidate (attainable)")
    ax.plot(dimensions, coefficients, "o", color="C1", ms=5, label="Four-outcome PGA: best attained")
    ax.set_xscale("log", base=2)
    ax.set_xticks([2, 8, 32, 128, 512, 2048, 4096])
    ax.set_xticklabels(["2", "8", "32", "128", "512", "2048", "4096"])
    ax.set(xlabel="Hilbert-space dimension d",
           ylabel=r"$C_{\mathrm{att}}=R/(\Delta H/\hbar)=2R/\sqrt{F_Q}$",
           title="Pure-state instantaneous observational entropy\nFour POVM outcomes; numerical lower bounds")
    ax.grid(alpha=.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=180)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)
    metadata = dict(outcomes=4, dimensions=dimensions, results=rows,
                    elapsed_seconds=time.perf_counter()-start_time,
                    python=platform.python_version(), numpy=np.__version__,
                    scipy=scipy.__version__, matplotlib=matplotlib.__version__,
                    method="Exact 2D pure-state reduction; scalar complement absent for d=2",
                    caveat="Attained lower bounds; no certified global optimum. Proportional outcomes allowed.")
    stem.with_suffix(".json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    np.savez_compressed(stem.with_suffix(".npz"), dimensions=dimensions,
                        coefficients=coefficients, **saved)
    np.savetxt(stem.with_suffix(".csv"), np.column_stack((dimensions, coefficients)),
               delimiter=",", header="dimension,attained_four_outcome_constant", comments="")
    print("ELAPSED", metadata["elapsed_seconds"], "PLOT", stem.with_suffix(".png"))


if __name__ == "__main__":
    main()
