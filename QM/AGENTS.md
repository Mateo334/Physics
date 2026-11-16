# Quantum chaos research workspace

This directory is a physics research workspace. The primary specialist agent is `quantum-research` in `.opencode/agents/quantum-research.md`. Its research standards apply whenever it is selected.

- Search locally and recursively before searching/downloading sources.
- Keep useful, verified findings in Markdown under `notes/`; index them in `notes/INDEX.md`.
- Start with `notes/INDEX.md` and `notes/RESEARCH_CONTEXT.md`; read only the relevant summary next. Open `notes/archive/` only for historical details, not as routine context.
- Store downloaded open-access PDFs under `papers/` without overwriting existing papers.
- Store Overleaf-ready LaTeX section fragments under `latex/` (no document preamble by default).
- Store Python source in `calculations/code/` and generated data/figures in `calculations/results/`. Read `calculations/README.md` for entry points; do not load result arrays or full logs into conversation unless needed.
- Treat proofs, numerical checks, conjectures and citations as distinct.
- Never silently alter existing thesis files or delete user data.
