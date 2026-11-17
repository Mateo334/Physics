"""Fixed-POVM Chernoff observational entropy, with nonselective backaction.

Uses the same preparation recipe as the four-model experiment, but at j=30.
Different PVM axes/block ranks and one full-rank unsharp POVM.
At t=0 only the outcome distribution is evaluated; each later step is U,
unrecorded Lüders measurement, then evaluation of its outcome distribution.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from chernoff_utils import mean_chernoff_entropy as entropy

J, COUNT, TIMES, SEED = 30, 100, np.arange(6), 20261116
KAPPAS = (2., 4., 5., 6.)


def initial_states():
    m = np.arange(-J, J+1, dtype=float)
    raising = np.diag(np.sqrt(J*(J+1)-m[:-1]*(m[:-1]+1)), -1)
    jy = (raising-raising.T)/(2j)
    jx = (raising+raising.T)/2
    eig, vectors = np.linalg.eigh(jy)
    north = vectors.conj().T[:, -1]  # Exactly the existing comparison's convention.
    rng = np.random.default_rng(SEED+1)
    states = np.column_stack([
        np.exp(-1j*phi*m) * (vectors @ (np.exp(-1j*theta*eig)*north))
        for theta, phi in [(np.arccos(rng.uniform(-1, 1)), rng.uniform(0, 2*np.pi))
                           for _ in range(COUNT)]
    ])
    assert np.allclose(np.sum(abs(states)**2, axis=0), 1, atol=1e-10)
    rotation = expm(-1j*np.pi*jy/2)
    return m, states, rotation, np.linalg.eigh(jx)[1]


def pvm_curve(states, u, basis, rank, first_probs=None):
    d = len(states)
    u = basis.conj().T @ u @ basis
    psi = basis.conj().T @ states
    group = np.arange(d)//rank  # The last cell has rank d % rank if nonzero.
    sizes = np.bincount(group)
    q = sizes/d
    mask = group[:, None] == group[None, :]

    def outcomes(rho):
        diag = rho.diagonal(axis1=1, axis2=2).real
        return np.stack([diag[:, group == i].sum(axis=1) for i in range(len(q))], axis=1)

    rho = np.einsum("an,bn->nab", psi, psi.conj())
    curve = [entropy(outcomes(rho), q, d)]
    for t in TIMES[1:]:
        rho = (u @ rho @ u.conj().T) * mask
        assert np.allclose(np.trace(rho, axis1=1, axis2=2), 1, atol=1e-10)
        assert np.min(np.linalg.eigvalsh(rho[0])) >= -1e-10
        curve.append(entropy(first_probs if t == 1 and first_probs is not None else outcomes(rho), q, d))
    return np.array(curve)


def unsharp_z_curve(states, u, eta, first_diag=None):
    """E_m=eta P_m+(1-eta) I/d, Kraus K_m=sqrt(E_m)."""
    d = len(states)
    a = np.sqrt((1-eta)/d)
    b = np.sqrt(eta+(1-eta)/d)
    lam = d*a*a+2*a*(b-a)
    assert np.isclose(1-lam, (b-a)**2)
    assert 0 < lam < 1
    rho = np.einsum("an,bn->nab", states, states.conj())
    q = np.ones(d)/d

    def outcomes(matrix):
        return eta * matrix.diagonal(axis1=1, axis2=2).real + (1-eta)/d

    curve = [entropy(outcomes(rho), q, d)]
    for t in TIMES[1:]:
        rho = u @ rho @ u.conj().T
        diagonal = np.zeros_like(rho)
        ix = np.arange(d)
        diagonal[:, ix, ix] = rho[:, ix, ix]
        rho = lam*rho+(1-lam)*diagonal
        assert np.allclose(np.trace(rho, axis1=1, axis2=2), 1, atol=1e-10)
        assert np.min(np.linalg.eigvalsh(rho[0])) >= -1e-10
        p = eta*first_diag+(1-eta)/d if t == 1 and first_diag is not None else outcomes(rho)
        curve.append(entropy(p, q, d))
    return np.array(curve)


def main():
    # Analytic two-outcome checks including the continuous-extension endpoint.
    assert np.isclose(entropy(np.array([[1., 0.]]), np.array([.5, .5]), 2), 0)
    assert np.isclose(entropy(np.array([[.5, .5]]), np.array([.5, .5]), 2), np.log(2))
    m, states, rotation, x_basis = initial_states()
    d = len(m)
    first_diag = abs(rotation @ states).T**2
    choices = {"Jz rank 1 (control)": (np.eye(d), 1),
               "Jz rank <=4": (np.eye(d), 4),
               "Jx rank 1": (x_basis, 1),
               "Jx rank <=4": (x_basis, 4),
               "Jx rank <=8": (x_basis, 8)}
    curves = {name: {} for name in (*choices, "unsharp Jz (eta=0.35, full-rank effects)")}
    for k in KAPPAS:
        u = np.exp(-1j*k*m**2/(2*J))[:, None]*rotation
        assert np.allclose(u.conj().T @ u, np.eye(d), atol=1e-11)
        for name, (basis, rank) in choices.items():
            # Analytically D_z(D_k X D_k^dagger)=D_z(X). Evaluating the
            # equivalent rotation-only channel also avoids spurious support
            # changes of tiny probabilities from finite-precision cancellation.
            effective_u = rotation if name == "Jz rank 1 (control)" else u
            first_p = (np.stack([first_diag[:, b*rank:(b+1)*rank].sum(axis=1)
                                 for b in range((d+rank-1)//rank)], axis=1)
                       if name.startswith("Jz") else None)
            curves[name][k] = pvm_curve(states, effective_u, basis, rank, first_p)
        curves["unsharp Jz (eta=0.35, full-rank effects)"][k] = unsharp_z_curve(states, u, .35, first_diag)

    control = np.array(list(curves["Jz rank 1 (control)"].values()))
    assert np.max(np.ptp(control, axis=0)) < 1e-12
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), sharex=True, sharey=True)
    data = {"times": TIMES}
    for ax, (name, group) in zip(axes.flat, curves.items()):
        for k, values in group.items():
            ax.plot(TIMES, values, marker="o", label=f"k={k:g}")
            data[f"{name} k={k:g}"] = values
        ax.set_title(name)
        ax.grid(alpha=.3)
    axes[0, 0].legend(fontsize=8)
    for ax in axes[1, :]:
        ax.set_xlabel("Kicks (t=0 before measurement)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Mean $S_C^{\mathcal{C}}$ (nats)")
    fig.suptitle(f"Chernoff observational entropy; fixed POVM, discarded Lüders outcomes; j={J}, {COUNT} states")
    fig.tight_layout()
    output = Path(__file__).resolve().parent.parent/"results"
    output.mkdir(parents=True, exist_ok=True)
    name = f"kicked_top_chernoff_observational_povm_ranks_j{J}_{COUNT}_states"
    fig.savefig(output/(name+".png"), dpi=160)
    plt.close(fig)
    np.savez(output/(name+".npz"), **data)
    print(output/(name+".png"))


if __name__ == "__main__":
    main()
