# Theory Documentation Setup

Use this one-time setup rule when constructing `docs/theory` for a new project
from original paper TeX sources.

- Start from the original paper sources.
  - Put the full main-paper TeX source in `docs/theory/main.tex`.
  - Put the full supplement/end-matter TeX source in `docs/theory/supp.tex`
    or another clearly named source file.
  - Keep the original PDF only as a fallback reference when TeX sources are
    insufficient.
- Build targeted theory-specific files from the original TeX sources.
  - Extract each major appendix or self-contained theory topic into its own
    standalone TeX file in `docs/theory`.
  - Preserve the original scientific content verbatim whenever extraction is the
    goal.
  - Only add the minimal preamble, packages, and `\begin{document}` /
    `\end{document}` wrapper needed to make the extracted file renderable.
  - Do not rewrite, summarize, or reinterpret the paper text during extraction
    unless the user explicitly asks for a summary.
- Create or update `docs/theory/theory_overview.md`.
  - Use it as the theory entry point for Codex.
  - Include a short description of each theory-specific file.
  - State when each file should be read, using concrete trigger topics such as
    "bad-cavity elimination", "weak-drive limit", or "single-particle
    decoherence".
  - Prefer a small label map when useful, for example
    `sec:weak -> docs/theory/appendix_weak_drive_limit.tex`.
  - Make clear that targeted theory files should be read before the full
    supplement or PDF.
- Update the theory hierarchy in `AGENTS.md`.
  - `AGENTS.md` should point to `docs/theory/theory_overview.md` as the theory
    entry point.
  - The hierarchy should list only theory-specific files that are actually
    referenced by `docs/theory/theory_overview.md`.
  - The hierarchy should include file names explicitly so Codex can route to the
    correct source without searching the full paper first.
- Validate renderability when TeX files are created.
  - Compile standalone TeX files with `docs/build_tex.sh` when possible.
  - Fix render-blocking errors.
  - Unresolved citations or cross-references are acceptable for extracted
    standalone sections if the PDF still builds.
