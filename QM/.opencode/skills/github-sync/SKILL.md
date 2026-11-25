---
name: github-sync
description: Push or pull the whole research workspace on GitHub when requested.
---

# GitHub synchronization

Repository: https://github.com/Mateo334/Physics.git
Branch: main
Remote directory: QM/
The local workspace root maps directly to QM/.

## Bulk sync

- “Push” authorizes staging and committing the whole workspace to QM/, then
  pushing main. “Pull” authorizes fetching QM/ into the whole workspace.
  No per-topic selection, content review, scientific tests or file-by-file reports.
- Use Git directly when configured for this mapping; otherwise reuse a clean
  temporary clone or create one. Keep repository content outside QM/ untouched.
- Transfer in bulk, including .opencode/. Exclude .git/, secrets (.env files,
  credentials/private keys), caches, temporary files and ignored build output.
  Do not force-add ignored files.
- Retain only mechanical safety checks: correct remote/branch/mapping, command
  exit codes, conflict detection and push/pull success. Use Git or an automated
  baseline comparison, not model inspection of individual files.
- Stop on conflicts or ambiguous overwrites; never force-push, reset, blindly
  mirror deletions or discard uncommitted work. If a copy-based pull has no valid
  baseline, ask before replacing differing local files.
- No research-note/checkpoint reading or bookkeeping commits on the normal path.
  Keep any needed sync baseline machine-readable outside the research notes.
- Reply with success and commit hash, or the blocking error; no inventory.
