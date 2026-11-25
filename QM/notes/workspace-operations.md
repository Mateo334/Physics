# Workspace operations — read only when needed

## GitHub synchronization
- Default repository https://github.com/Mateo334/Physics.git, branch `main`, subfolder `QM`. Workspace root maps directly to `QM/`.
- No local `.git` as of the last check. Pull through a temporary clone, compare against the recorded checkpoint, and preserve local edits and unrelated files. Never blindly mirror or delete.
- 2026-11-18 checkpoint: `8ac98d6d6bf5fbf258c563cc3700111310f9d8db`.
- 2026-11-22 fetched `5431c5d379013a086d897a72b40f587d25ddd73d`: 95 added files, updated calculations README and notes INDEX. Research context conflicted locally and was preserved; **its comparison base remains the November 18 commit**. Other synchronized files use the November 22 checkpoint. Subsequent local reorganizations are local edits, not upstream changes.
- The 2026-11-25 cleanup retains scientific paths and makes old index/command catalogs available as `notes/CATALOG.md` and `calculations/CATALOG.md`. Do not overwrite these changes during a future sync without comparison.

## Technical history
- [Working LaTeX terminal workaround](opencode-latex-fix.md): user-confirmed Unicode-cell renderer fix; restoring the old cache backup or refreshing/reinstalling the cache can remove it. Graphical rendering remained affected by the legacy API.
- Full plugin versions, diagnostic paths, temporary instrumentation, failure causes, test scope and superseded statuses are preserved verbatim in [pre-cleanup context](archive/reorganization-2026-11-25/notes/RESEARCH_CONTEXT.md#opencode-plugins-verified-2026-11-22). These are historical checkpoints, not current verification or authorization to modify external files.
- `opencode.cmd` avoids the blocked PowerShell shim. No service restart was needed for configuration changes; the terminal workaround required reopening the client.
- Compiler checks on 2026-11-25 found neither pdflatex nor tectonic. Prior standalone TeX is not compile-tested.

## Preservation
See [cleanup record](archive/reorganization-2026-11-25/README.md). Original archives are historical snapshots with original relative paths, not active navigation trees. Do not routinely read them or repeat already resolved diagnostic steps.
