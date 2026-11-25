"""Numerical checks accompanying the analytic binary-POVM global proof."""
from pathlib import Path
import json
import platform
import numpy as np
import scipy
from scipy.optimize import brentq


def main():
    target = Path(__file__).resolve().parents[1]/"results/binary_obsen_global_bound_checks.json"
    if target.exists():
        raise FileExistsError(target)
    rng = np.random.default_rng(20261123)
    rows = []
    for d in (2, 3, 21, 256, 4096):
        p = brentq(lambda p: (2*p-1)*np.log((d-1)*p/(1-p))-2,
                   .5, 1-1e-12, xtol=1e-14)
        q, s = 1-p, np.sqrt(p*(1-p))
        coefficient = 2*s*np.log((d-1)*p/q)
        size = min(d, 3)
        rho = np.zeros((size, size))
        rho[0, 0] = 1
        drho = np.zeros_like(rho)
        drho[0, 1] = drho[1, 0] = -1
        u = np.zeros(size)
        u[:2] = [np.sqrt(p), np.sqrt(q)]
        E = np.outer(u, u)
        I = np.eye(size)
        G1 = -np.log(p)*drho+2*s/p*rho-2*s*I
        G2 = -np.log(q/(d-1))*drho-2*s/q*rho+2*s/(d-1)*I
        D = G1-G2
        lam = float(u @ D @ u)
        multiplier = G2+lam*E
        slack1, slack2 = multiplier-G1, multiplier-G2
        slack_min = float(min(np.linalg.eigvalsh(slack1).min(),
                              np.linalg.eigvalsh(slack2).min()))
        residual = float(max(np.linalg.norm(slack1 @ E),
                             np.linalg.norm(slack2 @ (I-E))))
        assert slack_min > -1e-10 and residual < 1e-10
        assert abs((2*p-1)*np.log((d-1)*p/q)-2) < 1e-10
        highest = -np.inf
        for _ in range(2000):
            Z = rng.normal(size=(2, 2))+1j*rng.normal(size=(2, 2))
            U, _ = np.linalg.qr(Z)
            A = (U*rng.uniform(0, 1, 2)) @ U.conj().T
            a = A[0, 0].real
            V = np.trace(A).real+(d-2)*rng.uniform()
            rate = 2*A[0, 1].real*np.log(a*(d-V)/((1-a)*V))
            highest = max(highest, rate)
            assert rate <= coefficient+1e-10
        rows.append(dict(d=d, optimal_p=p, binary_global_coefficient=coefficient,
                         KKT_minimum_slack_eigenvalue=slack_min,
                         KKT_complementarity_residual=residual,
                         highest_random_candidate=float(highest), random_candidates=2000))
    summary = dict(seed=20261123, python=platform.python_version(),
                   numpy=np.__version__, scipy=scipy.__version__, results=rows,
                   scope="Numerical sanity checks only; the global guarantee comes from the analytic proof.")
    target.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
