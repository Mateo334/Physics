"""Continuous-K two-state rotor energy uncertainty times natural log(2)."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

stem = Path(__file__).resolve().parents[1]/"results/kicked_rotor_d2_DeltaH_ln2_K4to12"
if any(stem.with_suffix(s).exists() for s in (".png", ".csv")):
    raise FileExistsError("Refusing to overwrite existing curve")
hbar = 2*np.pi/64
k = np.linspace(4, 12, 4001)
values = hbar**2/4*np.abs(np.sin(2*k/hbar))*np.log(2)
probability = np.sin(k/hbar)**2
assert np.allclose(values, hbar**2/2*np.sqrt(probability*(1-probability))*np.log(2),
                   atol=1e-12, rtol=1e-10)
np.savetxt(stem.with_suffix(".csv"), np.column_stack((k, values)),
           delimiter=",", header="K,DeltaH_times_ln2", comments="")
fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(k, values, lw=1.4)
ax.set(xlabel="K", ylabel="Energy uncertainty × ln(2)",
       title="Two-state kicked rotor: immediately after first kick", xlim=(4, 12))
ax.grid(alpha=.25)
fig.tight_layout()
fig.savefig(stem.with_suffix(".png"), dpi=180)
