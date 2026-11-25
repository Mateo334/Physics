# Memory-preserving reorganization — 2026-11-25

User authorized workspace cleanup while retaining all memory. Strategy: shorten automatic/startup retrieval, preserve scientific filenames and all originals, and record latest conceptual corrections.

## Preserved originals
Subpaths mirror their original workspace locations except the root README, stored as `workspace-README.md`: `AGENTS.md`, `notes/INDEX.md`, `notes/RESEARCH_CONTEXT.md`, `calculations/README.md`, `papers/README.md`. The root README was reconstructed from its complete pre-edit read after a backup-name collision; integrity is checked against the pre-edit SHA256. Other files were copied before cleanup. Their internal relative links refer to the original layout; treat them as historical source text, not active navigation.

`baseline.csv` records SHA256, byte size and workspace-relative path for every existing file outside this backup directory before edits. This permits checking that manuscripts, source PDFs, code, results and existing detailed notes remain untouched. It is an integrity baseline, not a backup of the whole workspace.

## New retrieval layers
- Short `notes/INDEX.md` and `notes/RESEARCH_CONTEXT.md`.
- Full old research index retained as `notes/CATALOG.md`; full command guide retained as `calculations/CATALOG.md`, with active relative links.
- Short calculations guide and a LaTeX guide; papers annotations stay where existing instructions expect them.
- `notes/workspace-operations.md` routes historical technical details out of research startup context.
- `notes/open-problems/quantum-lyapunov-programme.md` saves the latest conversation's proposed spectrum, classical bridge, corrections, failed shortcuts and exact notebook overlap.

## Deletion policy
No scientific originals, PDFs, manuscripts, code, result arrays or unique research notes are deleted. Eight regenerable Python bytecode caches with verified adjacent source files were deleted; five caches without source at the expected path were conservatively retained. Exact paths are in `removed-caches.json`. Superseded navigation prose is removed from active startup files but retained here. Obsolete claims are labelled/corrected rather than erased, because failed approaches are research memory. Physical mass-renaming would break references without reducing reading tokens; existing scientific paths are intentionally preserved.

Validation and measured startup-size changes will be recorded in `verification.md`.
