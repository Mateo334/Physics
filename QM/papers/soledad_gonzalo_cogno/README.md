# Soledad Gonzalo Cogno: reading collection

The author is **Soledad Gonzalo Cogno**, computational neuroscientist at NTNU. This folder supports the user's requested concise introduction for a theoretical physicist.

## Start here
- Standalone introduction: [`../../latex/soledad_gonzalo_cogno_intro.tex`](../../latex/soledad_gonzalo_cogno_intro.tex).
- Publication-list source, checked 2026-11-24: https://gonzalocognolab.com/publications/
- Completed 2026-11-25: eight listed works covered from accessible full texts (including the short 2015 conference abstract and the 2026 preprint); the ninth, the 2013 AIP contribution, is metadata-only and explicitly not summarized scientifically. Preprint/final versions are not counted twice.
- Summaries distinguish full-text reading, abstract-only checks, and inaccessible PDFs. Author-listed does not imply peer-reviewed: the 2026 recurrent-network work is a preprint.

## Source-grounded annotations
- [Single neurons, noise and bursting](single_neurons.md).
- [Balanced visual networks and plasticity](balanced_networks.md).
- [Slow sequences and recurrent-network mechanisms](sequences.md).
- [CA1 population coding and uncertainty methodology](population_coding_and_uncertainty.md).

Each annotation records verified bibliographic information, download URLs, files, pages/sections checked, findings and limitations. Extracted text is a reading aid, not a replacement for the original equations and figures. No downloaded code was executed; no experiments were reproduced.

## Access and compilation
- CA1: full 2022 preprint HTML read and saved as `Nagelhus_et_al_2022_Object_centered_preprint_web.txt`; published 2023 metadata/abstract checked. PDF requests failed with HTTP 403. No substitute HTML is labelled as a PDF.
- Uncertainty audit: published PDF saved with an adjacent extracted text file; main text pp. 1–14 read.
- Balanced-network work: full 13-page Frontiers article and complete one-page BMC conference contribution downloaded and read; see `balanced_networks.md` for exact filenames and source locations.
- Single neurons: 2014 institutional manuscript and 2016 published article obtained and read. The 2016 biophysical supplement was inaccessible. The 2013 AIP full text was not obtained; its citation is retained only as bibliographic context.
- Sequences: published Nature PDF and 2026 bioRxiv v1 plus analytical supplement obtained; detailed reading ledger in `sequences.md`. Experimental findings and model-based mechanisms are distinguished.
- No `pdflatex`, `tectonic`, `lualatex` or `xelatex` command was available when checked. The standalone introduction is **not compile-tested**; no compiler was installed.
- Final checks (2026-11-25): all eight saved PDFs have PDF signatures and parse with pypdf (one publisher PDF emits duplicate MediaBox warnings). The introduction has nine references; environment nesting, brace/display delimiters, and citation keys pass basic structural checks. These are not a substitute for compilation or visual layout inspection.
