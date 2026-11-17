"""One plot of mean observational entropy for four finite quantum models."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import xlogy

from four_model_setup import TIMES, models


def entropy_curve(states, unitary, dimension):
    """Fixed rank-one basis PVM: S_obs is the Shannon entropy of |psi_i|^2."""
    result = []
    for t in TIMES:
        probabilities = abs(states) ** 2
        assert np.allclose(probabilities.sum(axis=0), 1, atol=1e-10)
        values = -np.sum(xlogy(probabilities, probabilities), axis=0)
        assert np.all((values >= -1e-10) & (values <= np.log(dimension) + 1e-10))
        result.append(values.mean() / np.log(dimension))
        if t < TIMES[-1]:
            states = unitary @ states
    return np.array(result)


setups = models()
curves = {label: entropy_curve(states, unitary, len(states))
          for label, (states, unitary) in setups.items()}
assert np.isclose(curves["Random-matrix Hamiltonian (GOE)"][0],
                   curves["Kicked top (k=2)"][0], atol=1e-12)

fig, ax = plt.subplots(figsize=(9, 5))
for label, curve in curves.items():
    ax.plot(TIMES, curve, label=label)
ax.set(xlabel="Kicks / normalized Hamiltonian time (different physical clocks)",
       ylabel="Mean observational entropy / log(sector dimension)",
       title="Four models, kicked-top k=2,4,5,6: 100 coherent states each")
ax.set_ylim(0, 1.05)
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.tight_layout()
output = Path(__file__).resolve().parent.parent / "results" / "obsen_four_models_k2_4_5_6_5kicks_100_states.png"
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=160)
plt.close(fig)
print(output)
print(f"Dimensions: kicked top and GOE={len(setups['Kicked top (k=2)'][0])}, "
      f"SYK parity={len(setups['Majorana SYK (even parity)'][0])}, "
      f"Bose--Hubbard={len(setups['Bose--Hubbard (3 sites)'][0])}")
