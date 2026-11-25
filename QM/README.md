# Quantum chaos research workspace

Research notes, source papers, LaTeX fragments and reproducible calculations for **Quantum Kolmogorov–Sinai Entropy**.

## Start here
- [Research index](notes/INDEX.md) and [current context](notes/RESEARCH_CONTEXT.md) — short startup reading; follow only the task-relevant branch.
- [Current Lyapunov decision ledger](notes/open-problems/quantum-lyapunov-programme.md) — proposed constructions, notebook overlap, corrections and unresolved classical correspondence.
- [Calculations](calculations/README.md) — scripts in `calculations/code/`, generated data and figures in `calculations/results/`.
- [Papers](papers/README.md) — source annotations; [LaTeX](latex/README.md) — originals and generated deliverables, with status warnings.
- [Full research catalog](notes/CATALOG.md) and [command catalog](calculations/CATALOG.md) — load selectively, not as startup context.

## Main calculations
Choose a topic in the [calculation guide](calculations/README.md), then inspect its prerequisites and overwrite warnings before running from the workspace root.

The previously advertised `kicked_top_simple.py` and `kicked_top_random_pvms.py` are absent from this checkout (2026-11-25). Their [historical results and prerequisites](notes/calculations/kicked-top-simple-povms-lyapunov.md#run-and-find-results) remain documented; recover verified originals before rerunning.

## OpenCode
Run `opencode` here and select `quantum-research`. Shared rules: `AGENTS.md`; primary and three specialist prompts: `.opencode/agents/`; on-demand GitHub sync: `.opencode/skills/github-sync/SKILL.md`. A configured model is required; downloading papers requires network access. Review research output before incorporating it into the thesis.
