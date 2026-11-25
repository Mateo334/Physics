# Quantum chaos research workspace

Research notes, source papers, LaTeX fragments and reproducible calculations for **Quantum Kolmogorov–Sinai Entropy**.

## Start here
- [Research index](notes/INDEX.md) and [current context](notes/RESEARCH_CONTEXT.md) — short startup reading; follow only the task-relevant branch.
- [Current Lyapunov decision ledger](notes/open-problems/quantum-lyapunov-programme.md) — proposed constructions, notebook overlap, corrections and unresolved classical correspondence.
- [Calculations](calculations/README.md) — scripts in `calculations/code/`, generated data and figures in `calculations/results/`.
- [Papers](papers/README.md) — source annotations; [LaTeX](latex/README.md) — originals and generated deliverables, with status warnings.
- [Full research catalog](notes/CATALOG.md) and [command catalog](calculations/CATALOG.md) — load selectively, not as startup context.

## Main calculations
From the workspace root, with Python, NumPy, SciPy and Matplotlib installed:

```powershell
python calculations/code/kicked_top_simple.py
python calculations/code/kicked_top_random_pvms.py
```

The random-PVM run needs saved simple data/summary **and the existing k=4..7 random-PVM regression fixture**; it is not a standalone fresh-fixture generator. See the [authoritative calculation note](notes/calculations/kicked-top-simple-povms-lyapunov.md#run-and-find-results) for filenames.

## OpenCode
Run `opencode` here and select `quantum-research`. Project instructions are in `AGENTS.md`; specialist configuration is in `.opencode/agents/quantum-research.md`. A configured model is required; downloading papers requires network access. Review research output before incorporating it into the thesis.
