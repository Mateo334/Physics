# OpenCode LaTeX rendering fix

Confirmed working by user on Windows, OpenCode 2.3.2, `opencode-latex@0.2.0` (2026-11-22).

**Symptom:** fenced `latex` blocks show raw TeX, even with `/latex-render-mode` → **Unicode cells**.

**Cause:** `opentui-math` uses an obsolete OpenTUI layout API. Captured error: `TypeError: undefined is not an object (evaluating 'this.yogaNode.setMeasureFunc')` in `LatexRenderable.setupMeasureFunction`. Not a PowerShell or missing TeX-compiler problem.

## Apply on another PC

1. Locate the active cached package, typically `~/.cache/opencode/npm/opencode-latex@0.2.0/<cache-id>/node_modules/opencode-latex/src/tui.ts`. Cache IDs differ. Back up this file before editing.
2. Add imports:

```ts
import { TextRenderable } from "@opentui/core"
import { renderLatexToString } from "opentui-math"
```

3. Inside the existing `registerCodeBlockRenderer` callback's `try`, **before** constructing `GraphicalLatexRenderable`, add:

```ts
if (effectiveMode() === "cells") {
  return new TextRenderable(context.renderer, {
    content: renderLatexToString(normalizeLatexSource(token.text), {
      displayMode: true,
    }),
    fg: color,
    wrapMode: "none",
  })
}
```

These names already exist in version 0.2.0. Preserve the surrounding error fallback.

4. Choose **Unicode cells**, then close and reopen the terminal client to load the cache edit. No background-service restart needed; `/reload` alone did not activate our earlier cache patch.
5. Verify a fenced `latex` expression `\frac{a+b}{2}` displays a stacked fraction, not source.

**Limits:** fixes cell rendering only; graphical paths retain the incompatible API. Plugin reinstall/cache refresh may erase the patch. Prefer an upstream compatible release when available; do not apply blindly to newer code. Restore the backup to undo.

Local regression test: `node calculations/code/test_latex_cell_fix.cjs` (adjust its cache path on another PC). Tests the actual formatter with a mocked host renderer; user separately confirmed real terminal rendering.
