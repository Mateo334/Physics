"""Classical Chernoff observational entropy of rows of POVM probabilities."""

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import logsumexp


def mean_chernoff_entropy(probabilities, reference, dimension):
    """Average log(d) + log(min_s sum_i p_i^s q_i^(1-s)).

    At s=0 zero-probability outcomes contribute zero by continuous extension.
    Remove any zero-trace effects before calling (all reference entries > 0).
    """
    p, q = np.asarray(probabilities), np.asarray(reference)
    assert p.ndim == 2 and p.shape[1] == len(q)
    assert np.all(q > 0) and np.isclose(q.sum(), 1)
    assert np.all(p >= -1e-12) and np.allclose(p.sum(axis=1), 1, atol=1e-10)
    values = []
    for row in p:
        positive = row > 0
        log_p, log_q = np.log(row[positive]), np.log(q[positive])

        def log_overlap(s):
            return logsumexp(s*log_p + (1-s)*log_q)

        optimum = minimize_scalar(log_overlap, bounds=(0, 1), method="bounded")
        values.append(np.log(dimension) + min(log_overlap(0), log_overlap(1), optimum.fun))
    mean = np.mean(values)
    assert -1e-9 <= mean <= np.log(dimension) + 1e-9
    return mean
