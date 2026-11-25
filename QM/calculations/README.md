# Calculations: choose a topic first

- `code/`: runnable scripts and shared helpers; existing imports and filenames are unchanged.
- `results/`: preserved data, figures, tables and logs. Read the associated note/JSON summary before opening arrays or full logs.
- [Command catalog](CATALOG.md): all existing commands, dependencies, fixtures, parameters and overwrite warnings. Search for the model/task; do not load it routinely.
- [Calculation-note catalog](../notes/CATALOG.md#numerical-calculations): interpretation, convergence and numerical limitations.

## Common entry points
| Task | Read first |
| --- | --- |
| Classical Lyapunov / simple kicked top | [matched tests](../notes/calculations/kicked-top-simple-povms-lyapunov.md) |
| Closed-system tangent / isometry | [analytic and numerical controls](../notes/calculations/quantum-tangent-geometry.md) |
| Probability/relative-entropy derivatives | [Fréchet check](../notes/calculations/kicked-top-frechet.md) |
| Pure-state OE optimization | [current all-outcome proof](../notes/theorems/all-outcome-obsen-global-maximum.md) |
| Husimi phase-space plots | [static and animation notes](../notes/calculations/kicked-top-husimi-k3-k6.md) |
| Other top / rotor / SYK / coupled models | [full catalog](CATALOG.md) |

Run commands from the workspace root. Reuse existing scripts. Do not infer convergence from a saved plot or overwrite old results. `kicked_top_random_pvms.py` needs saved simple data **and** the existing k=4..7 regression fixture; see [exact prerequisites](../notes/calculations/kicked-top-simple-povms-lyapunov.md#run-and-find-results).
