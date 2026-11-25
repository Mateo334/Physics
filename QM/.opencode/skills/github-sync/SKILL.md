---
name: github-sync
description: Synchronize this research workspace with its GitHub repository when the user requests a pull or sync.
---

# GitHub synchronization

Repository: https://github.com/Mateo334/Physics.git
Branch: main
Remote directory: QM/
The local workspace root maps directly to QM/.

## Procedure

All workspace paths below are relative to the workspace root, not this skill directory.

1. Read the sync checkpoint in notes/workspace-operations.md.
2. Inspect local modifications and repository status.
3. Fetch remote changes without overwriting local files.
   If no local .git exists, use a temporary clone.
4. Compare local and remote changes against the checkpoint.
5. Apply nonconflicting remote updates.
6. Preserve local edits and report conflicts.
7. Update the checkpoint after successful synchronization.
8. Report added, updated, skipped and conflicting files.

## Restrictions

- Synchronize only when explicitly requested.
- Never automatically delete or overwrite conflicting files.
- Never push, commit or reset without authorization.
- If the checkpoint is missing or invalid, do not assume
  that differing local files are safe to overwrite.
- Clean up only temporary files created during synchronization.
