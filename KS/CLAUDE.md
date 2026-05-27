# KS — Working Rules

Read `tasks.md` in this folder to find what to work on.

- All equations are written in LaTeX. Never skip steps — every derivation must be written out in full detail.
- For computation, use SageMath, NumPy, SciPy, or similar. Numerical verification of analytical results is encouraged.
- Only write results that are rigorous and foolproof. No hand-waving.
- Experimentation is allowed and encouraged — explore, conjecture, test. But only commit results that are solid.
- You may create up to 10 subfolders in this directory, each corresponding to a distinct approach or method of calculation.
- Subfolders follow a strict hierarchy: always begin work in the lowest-numbered existing subfolder. Only open a new subfolder when the current approach has reached a genuine computational dead end — meaning no further progress is possible with the tools available. Do not open a new subfolder just because progress is slow or difficult.
- Write to `Output.tex` incrementally as you work — do not accumulate results in memory and write only at the end. After completing each definition, lemma, or derivation block, append it to `Output.tex` immediately. The session may be cut off at any time.
- Keep a running `progress.md` in the active subfolder. Update it briefly after each significant step (a definition pinned down, a proof step completed, a dead end hit). This is your checkpoint file — it must reflect the current state of understanding at all times.
- Before finishing each session, write a `summary.md` in the subfolder you worked in. It should describe what was attempted, what was established, and — if the approach is exhausted — why it is a dead end.
- All LaTeX output files must be named or marked `Output` (e.g. `Output.tex`) inside the relevant subfolder.
- Each `Output.tex` begins with a `\section*{Notes}` written by Mateo. Read it carefully before working — it contains directions, focus areas, and strategy hints that override general priorities. Never delete or modify this section. Append all your work below it.
- You may use any PDF files found in this folder or any of its subfolders as reference material.
- You may also use online resources — especially arXiv, Scholarpedia, and Wikipedia.
