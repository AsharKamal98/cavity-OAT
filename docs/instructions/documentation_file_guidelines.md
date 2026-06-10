# Documentation File Guidelines

Read this before creating a new file in `docs/theory` or `docs/instructions`.

## Core Rules

- If the user asks for a new theory or instruction file in TeX format, make it renderable by itself.
- For standalone TeX files, use only the minimal LaTeX preamble needed to compile cleanly, for example `article`, `geometry`, `amsmath`, `amssymb`, and `hyperref` when useful.
- After creating or editing a standalone TeX file, compile it with `docs/build_tex.sh` when possible and fix render errors before finishing.
- Files may reference other documentation files. Prefer references over re-explaining material that already has a dedicated source-of-truth file.
- If a topic already has a theory file and an instruction file, keep the split clear: theory files explain equations and physics assumptions, while instruction files explain how the repository should implement or use them.
- Always update AGENTS.md file accordingly.

## Reference Style

When pointing to another document, use the repository path, for example:

```text
docs/instructions/bloch_vector_averaging.tex
```

For TeX files, use escaped underscores inside `\texttt{...}`:

```tex
\texttt{docs/instructions/bloch\_vector\_averaging.tex}
```

## Avoid Duplication

Before adding a long explanation, check whether an existing file already covers
the idea. If it does, briefly state the local relevance and refer to that file.
This keeps the documentation maintainable and avoids conflicting versions of
the same rule or formula.
