# Workspace operations — read only when needed

## GitHub synchronization
- 2026-11-25 published the workspace to `main` in `7f953567c574d544789283e111e2ff428af4b177` (154 changed files); remote hash verified. This supersedes the older comparison bases below, including research context. Python bytecode was excluded, remote-only files preserved, and no repository files deleted. The following bookkeeping commit records this checkpoint; scientific tests and TeX compilation were not rerun.
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

## Workflow audit — 2026-11-25
- Preserved primary `quantum-research` plus programmer/skeptic/writer, scientific validation requirements, result statuses and standalone TeX default. Shortened primary/writer prompts; bounded delegation now uses supplied context first and returns evidence/paths rather than repeated background or full artifacts.
- Fixed programmer's writer-role typo and canonical-note ownership conflict; reconciled AGENTS.md's fragment-only rule with writer's standalone default. Writer now defines the nonstandard trace macro and separates TeX-file from chat formatting.
- Renamed `SKILL,md` to `SKILL.md`; OpenCode advertised `github-sync` after the change. Corrected its checkpoint pointer to this file. No synchronization performed.
- Navigation audit checked 135 local Markdown targets across 66 nonarchive Markdown files, including six fragments, with no broken targets reported. Independently confirmed the missing `kicked_top_simple.py` / `kicked_top_random_pvms.py`; marked historical commands unavailable, without reconstructing scripts or changing recorded scientific findings. Corrected stale catalog descriptions and added QFI/prototype retrieval pointers.
- [Numerically tested] Text-size proxy, not tokenizer/billing measurement: primary prompt 4282 → 3236 characters (24.4% reduction); shared rules plus primary 5422 → 4607 (15.0%); shared rules plus all four agent files 9911 → 8703 (12.2%). The last total is not a per-request prompt size. Retrieval/handoff savings remain unmeasured.
- Validation: four agent frontmatter headers and corrected skill filename checked; SHA-256 baseline of 307 files confirmed no changes to papers, manuscripts, Python scripts or result artifacts. Only workflow/navigation documents and the calculation note's reproducibility warning changed. No scientific runs, TeX compilation, global/ancestor configuration audit or end-to-end research-output equivalence test.
- Original edited prompts/startup files and full hash baseline: `C:\Users\matej\AppData\Local\Temp\opencode\workflow-audit-20261125-202349`. This is a local temporary backup, not portable research memory. Existing sync checkpoint above is unchanged.
