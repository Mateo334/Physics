# Cleanup verification

2026-11-25; checks run with local Python standard library. Machine-readable results: `verification.json`.

- Baseline: **295 existing files**; **282 remain byte-identical** by SHA256.
- Exactly five existing navigation/instruction files changed: root README, AGENTS, calculations README, notes INDEX, notes RESEARCH_CONTEXT.
- Exactly eight generated `.pyc` files removed, each with verified source at its expected path. Five caches without adjacent source retained. No scientific source, manuscript, PDF, data, figure, detailed note or executable script deleted or modified.
- All six preserved navigation originals match their baseline SHA256. Root README was recovered byte-for-byte after a backup-name collision; its hash matches `7E90316FC7EDE7C950D133FBA0B8D3089260899D18034694B16AA34BC55719AF`.
- Every original research-index bullet remains in `notes/CATALOG.md`; all **32 original Python command bullets** remain in `calculations/CATALOG.md`.
- Checked **135 active Markdown relative file-link targets**: zero missing. Historical archives were excluded because their copied relative paths intentionally reflect the old location. The parser checks explicit Markdown file links, not every inline-code filename, URL, image or heading anchor.
- Required `notes/INDEX.md#numerical-calculations` heading retained. Command catalog now has selective-read headings.
- Default research startup (`INDEX.md` + `RESEARCH_CONTEXT.md`): **23,622 → 5,112 characters, a 78.4% reduction**. This measures characters, not exact model tokens; AGENTS/model instructions and selected task notes are additional context.
- Calculations landing page: **9,283 → 1,611 characters**, while original detailed commands remain in the catalog.

All comparisons/assertions passed after the backup filename was corrected. No research computation, package installation, external configuration edit or TeX compilation performed. Latest mathematical discussion is summarized in `notes/open-problems/quantum-lyapunov-programme.md`; this organization check does not certify every historical scientific claim.
