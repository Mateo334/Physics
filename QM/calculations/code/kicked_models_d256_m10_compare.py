"""Matched-dimension pure-state instantaneous OE rates: top, rotor and saved SYK."""
from pathlib import Path
import json
import platform

import numpy as np
import scipy
from scipy.linalg import expm
from canonical_instantaneous_obsen_coefficient import run, objective
from kicked_rotor_K10_postkick_instantaneous import setup, free


def check_model(name, psi, generator, hbar, A, alpha):
    d = len(psi)
    assert np.linalg.norm(generator-generator.conj().T) < 1e-11
    assert abs(np.linalg.norm(psi)-1) < 1e-11
    mean = np.vdot(psi, generator @ psi).real
    tangent = -1j*(generator @ psi-mean*psi)
    spread = np.linalg.norm(tangent)
    B = np.column_stack((psi, -tangent/spread))
    assert np.linalg.norm(B.conj().T @ B-np.eye(2)) < 1e-10
    effects = B[None] @ A @ B.conj().T[None]+alpha[:, None, None]*(np.eye(d)-B @ B.conj().T)
    completeness = float(np.linalg.norm(effects.sum(axis=0)-np.eye(d)))
    minimum = float(np.linalg.eigvalsh(effects).min())
    assert completeness < 1e-10 and minimum > -1e-12
    p = np.einsum("a,iab,b->i", psi.conj(), effects, psi).real
    dp = 2*np.einsum("a,iab,b->i", psi.conj(), effects, tangent).real
    volumes = np.trace(effects, axis1=1, axis2=2).real
    assert p.min() > 0 and abs(p.sum()-1) < 1e-10
    rate = float(-dp @ np.log(p/volumes))
    qfi = float(4*np.vdot(tangent, tangent).real)
    values, vectors = np.linalg.eigh(generator)
    coordinates = vectors.conj().T @ psi
    def entropy(t):
        state = vectors @ (np.exp(-1j*values*t)*coordinates)
        probabilities = np.einsum("a,iab,b->i", state.conj(), effects, state).real
        assert probabilities.min() > 0
        return float(-probabilities @ np.log(probabilities/volumes))
    dt = 1e-5/max(1, spread)
    fd = (entropy(dt)-entropy(-dt))/(2*dt)
    assert abs(fd-rate)/max(1, abs(rate)) < 1e-7
    assert abs(rate/spread-objective(A, alpha, d)[0]) < 1e-10
    return dict(model=name, physical_energy_uncertainty=float(hbar*spread),
                hbar=hbar, generator_uncertainty=float(spread), physical_time_qfi=qfi,
                attained_rate=rate, coefficient=rate/spread,
                coefficient_from_qfi=2*rate/np.sqrt(qfi),
                completeness_error=completeness, minimum_effect_eigenvalue=minimum,
                finite_difference_rate=fd, finite_difference_step=dt), B


def rotor_state(d, k):
    q, p, g, psi = setup(d)
    evolved = free(psi, g, 1)
    psi = np.fft.fft(np.exp(-1j*k*np.cos(q)/(2*np.pi/64)) *
                     np.fft.ifft(evolved, norm="ortho"), norm="ortho")
    return p, g, psi


def main():
    out = Path(__file__).resolve().parents[1]/"results"
    stem = out/"kicked_models_d256_m10_compare"
    for suffix in (".json", ".npz", ".csv"):
        if stem.with_suffix(suffix).exists():
            raise FileExistsError(stem.with_suffix(suffix))
    d, m, j = 256, 10, 127.5
    optimization, A, alpha = run(d, m, 4, 300, 20261124)
    values, vectors = np.linalg.eigh(A)
    A = (vectors*np.maximum(values, 0)[:, None, :]) @ vectors.conj().transpose(0, 2, 1)
    values, vectors = np.linalg.eigh(A.sum(axis=0))
    inverse = (vectors*values**(-.5)) @ vectors.conj().T
    A = inverse[None] @ A @ inverse[None]
    saved, rows = dict(A=A, alpha=alpha), []
    z = np.arange(-j, j+1)
    raising = np.diag(np.sqrt(j*(j+1)-z[:-1]*(z[:-1]+1)), -1)
    jy = (raising-raising.T)/(2j)
    north = np.eye(d, dtype=complex)[:, -1]
    initial = np.exp(-1j*.7*z)*(expm(-1j*1.1*jy) @ north)
    rotation = expm(-1j*np.pi*jy/2)
    for k in (3., 6.):
        state = expm(-1j*.1*np.pi*jy) @ (np.exp(-1j*k*z*z/(2*j))*(rotation @ initial))
        row, B = check_model(f"top_K{k:g}", state, np.pi*jy, 1., A, alpha)
        row.update(k=k, spin=j, evaluation_time=1.1,
                   phase="rotation: H=pi Jy; twist half-period H=K Jz^2/j")
        rows.append(row)
        saved[f"top_K{k:g}_state"], saved[f"top_K{k:g}_frame"] = state, B
    for k in (6., 10.):
        p, g, state = rotor_state(d, k)
        row, B = check_model(f"rotor_K{k:g}", state, np.diag(g), 2*np.pi/64, A, alpha)
        # Fixed-hbar larger cutoff checks dynamics, not the d-dependent OE coefficient.
        p2, g2, state2 = rotor_state(512, k)
        spread2 = np.linalg.norm((g2-np.vdot(state2, g2*state2).real)*state2)
        indices = (np.fft.fftfreq(d, d=1/d).astype(int)) % 512
        overlap = np.vdot(state, state2[indices])
        outside = float(1-np.sum(abs(state2[indices])**2))
        cutoff_error = float(abs(spread2/row["generator_uncertainty"]-1))
        row.update(k=k, evaluation_time="1+",
                   phase="free immediately after first delta kick; H_phys=p^2/2",
                   cutoff_relative_spread_error=cutoff_error,
                   cutoff_overlap_magnitude=float(abs(overlap)),
                   probability_outside_d256_at_d512=outside,
                   cutoff_converged_at_1e_minus_5=bool(cutoff_error < 1e-5))
        rows.append(row)
        saved[f"rotor_K{k:g}_state"], saved[f"rotor_K{k:g}_frame"] = state, B
        saved["rotor_momentum"] = p
    syk = json.loads((out/"syk_N16_instantaneous_pga_pair.json").read_text())
    summary = dict(d=d, outcomes=m, optimizer=optimization, results=rows,
                   syk_reference=syk["results"], syk_coefficient=syk["feasible_coefficient"],
                   normalized_constant="rate/(DeltaH/hbar)=2*rate/sqrt(FQ_time)",
                   top_initial="theta=1.1, phi=0.7 spin coherent",
                   rotor_initial="q0=pi/2, p0=0, Gaussian sigma_q=sqrt(hbar/2)",
                   rotor_hbar=2*np.pi/64, rotor_period=1,
                   rotor_cutoff="d=256 momentum truncation, NOT hbar=2pi/256 torus",
                   python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                   caveat="Attained lower bounds, not global certificates or chaos indicators. Different physical energy units.")
    stem.with_suffix(".json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    np.savez_compressed(stem.with_suffix(".npz"), **saved)
    with stem.with_suffix(".csv").open("w", encoding="utf-8") as f:
        f.write("model,DeltaH,hbar,Delta_generator,FQ_time,attained_rate,coefficient\n")
        for r in rows:
            f.write(",".join(str(r[key]) for key in ("model", "physical_energy_uncertainty", "hbar",
                         "generator_uncertainty", "physical_time_qfi", "attained_rate", "coefficient"))+"\n")
    print(json.dumps({"results": rows, "syk_coefficient": syk["feasible_coefficient"]}, indent=2))


if __name__ == "__main__":
    main()
