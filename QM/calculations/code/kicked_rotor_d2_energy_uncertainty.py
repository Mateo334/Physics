"""Two-state rotor cutoff at the previously used fixed effective hbar."""
from pathlib import Path
import json
import platform

import numpy as np
from kicked_models_d256_m10_compare import rotor_state


def main():
    stem = Path(__file__).resolve().parents[1]/"results/kicked_rotor_d2_K8to13_energy"
    for suffix in (".json", ".csv"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    hbar = 2*np.pi/64
    rows = []
    for k in range(8, 14):
        p, generator, state = rotor_state(2, k)
        energy = p*p/2
        probabilities = abs(state)**2
        assert abs(probabilities.sum()-1) < 1e-12
        mean = float(probabilities @ energy)
        delta = float(np.linalg.norm((energy-mean)*state))
        # The sampled symmetric Gaussian is constant on q=(0,pi).
        analytic = hbar**2/4*abs(np.sin(2*k/hbar))
        assert abs(delta-analytic) < 1e-12
        rows.append(dict(K=k, energy_mean=mean, physical_energy_uncertainty=delta,
                         generator_uncertainty=delta/hbar,
                         momentum_minus_hbar_probability=float(probabilities[1]),
                         analytic_check_absolute_error=float(abs(delta-analytic))))
    summary = dict(d=2, hbar=hbar, evaluation_time="1+",
                   preparation="Previous periodized Gaussian q0=pi/2, p0=0, sigma_q=sqrt(hbar/2), sampled at q=0,pi",
                   protocol="Free evolution for one unit, then exp(-i K cos(q)/hbar); FFT momentum indices 0,-1",
                   energy="H_phys=p^2/2; evaluate after kick, not during delta impulse",
                   caveat="Severely truncated two-state grid model, not a converged rotor or chaos diagnostic",
                   python=platform.python_version(), numpy=np.__version__, results=rows)
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    np.savetxt(stem.with_suffix(".csv"),
               [[r["K"], r["physical_energy_uncertainty"], r["generator_uncertainty"]] for r in rows],
               delimiter=",", header="K,DeltaH_phys,DeltaH_over_hbar", comments="")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
