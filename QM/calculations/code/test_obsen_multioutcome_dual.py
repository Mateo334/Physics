"""Test a proposed per-effect dual inequality, not a proof."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import brentq, minimize_scalar


def main():
    target = Path(__file__).resolve().parents[1]/"results/obsen_multioutcome_dual_test.json"
    if target.exists():
        raise FileExistsError(target)
    rows = []
    for d in (2, 3, 4, 8, 21, 256, 4096):
        p = brentq(lambda p: (2*p-1)*np.log((d-1)*p/(1-p))-2, .5, 1-1e-12, xtol=1e-14)
        q, z = 1-p, np.sqrt(p*(1-p))
        lam = 2*p/(z*(2*p-1))-2*z*d/(d-1)
        k = 2*z/(d-1)
        a = k-2*z/q+lam*p
        b = np.log(q/(d-1))+lam*z
        c = lam*q
        def gap(x):
            y = max(1+x*x, -2*x/k)
            return a+2*b*x+c*x*x+k*(y-1)+2*x*np.log(y)
        grid = np.unique(np.r_[np.linspace(-30, 30, 30001), -np.sqrt(p/q), np.sqrt(q/p)])
        gaps = np.array([gap(x) for x in grid])
        minima = []
        for i in range(1, len(grid)-1):
            if gaps[i] < gaps[i-1] and gaps[i] < gaps[i+1]:
                result = minimize_scalar(gap, bounds=(grid[i-1], grid[i+1]), method="bounded")
                minima.append(dict(x=float(result.x), gap=float(result.fun)))
        assert gaps.min() > -1e-9
        rows.append(dict(d=d, p=p, dual_a=a, dual_b=b, dual_c=c, dual_k=k,
                         negative_branch_value_at_zero=a-k,
                         minimum_sampled_gap=float(gaps.min()), local_minima=minima))
    target.write_text(json.dumps(rows, indent=2))
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
