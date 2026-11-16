# Quantum chaos research workspace

Research notes, source papers, LaTeX fragments and reproducible calculations for **Quantum Kolmogorov–Sinai Entropy**.

## Start here
- [Research index](notes/INDEX.md) — concise current findings and links to detailed archives.
- [Research context](notes/RESEARCH_CONTEXT.md) — conventions, established results and open gaps.
- [Calculations](calculations/README.md) — scripts in `calculations/code/`, generated data and figures in `calculations/results/`.
- `papers/` — source PDFs; `latex/` — thesis-ready fragments and standalone presentations.

## Main calculations
From the workspace root, with Python, NumPy, SciPy and Matplotlib installed:

```powershell
python calculations/code/kicked_top_simple.py
python calculations/code/kicked_top_random_pvms.py
```

The random-PVM run needs saved simple data/summary **and the existing k=4..7 random-PVM regression fixture**; it is not a standalone fresh-fixture generator. See the calculations guide for filenames.

## OpenCode
Run `opencode` here and select `quantum-research`. Project instructions are in `AGENTS.md`; specialist configuration is in `.opencode/agents/quantum-research.md`. A configured model is required; downloading papers requires network access. Review research output before incorporating it into the thesis.
