# Quantum chaos research workspace

This directory is a physics research workspace. The primary specialist agent is `quantum-research` in `.opencode/agents/quantum-research.md`. Its research standards apply whenever it is selected.

## GitHub repository and pull convention
- Default repository: https://github.com/Mateo334/Physics.git; branch: `main`; folder: `QM` (https://github.com/Mateo334/Physics/tree/main/QM).
- This workspace root maps directly to `QM/`. When asked to pull from GitHub, use this source without asking for the URL again.
- Read `notes/RESEARCH_CONTEXT.md` for the last sync checkpoint. Without a local `.git`, fetch through a temporary clone and compare against that checkpoint before copying updates. Preserve local edits and unrelated files; report conflicts and never delete files automatically.

- Search locally and recursively before searching/downloading sources.
- Keep useful, verified findings in Markdown under `notes/`; index them in `notes/INDEX.md`.
- Start with `notes/INDEX.md` and `notes/RESEARCH_CONTEXT.md`, then one relevant summary. Full inventories live in `notes/CATALOG.md` and `calculations/CATALOG.md`; search them only when needed. Open `notes/archive/` only for historical details, not routine context.
- Store downloaded open-access PDFs under `papers/` without overwriting existing papers.
- Store Overleaf-ready LaTeX section fragments under `latex/` (no document preamble by default).
- Store Python source in `calculations/code/` and generated data/figures in `calculations/results/`. Read `calculations/README.md` for entry points; do not load result arrays or full logs into conversation unless needed.
- Treat proofs, numerical checks, conjectures and citations as distinct.
- Never silently alter existing thesis files or delete user data.
