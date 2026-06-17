# Documentation File Guidelines

Read this before creating or modifying files in `docs/theory` or
`docs/instructions`.

## How to write instruction files

Instruction files should explain how the repository should implement or use a
piece of theory. They should not be full derivations, and they should not try to
document every line of code.

Use this structure when it fits the task:

1. Method in theory.
   - State the physical or mathematical method in the same order the code should
     follow.
   - Include only the key equations needed to define the implementation
     contract.
   - If a derivation already exists in a theory file, refer to that file instead
     of re-deriving it.
   - If the code needs a derived quantity, state the minimal operation clearly.
     For example, if the theory file defines `LHS = RHS`, the instruction file
     can say that the diagnostic should compute `residual = LHS - RHS` and
     minimize or plot that residual.
   - Briefly state any approximation, convention, normalization, gauge choice etc.
     that changes how the theory should be used in code.

2. Method in pseudo-code.
   - Describe the implementation flow in pseudo-code or function-level steps.
   - Prefer data-flow descriptions over long prose.
   - Before proposing new helpers, scan the existing codebase for functions that
     already implement the same or closely related logic.
   - If an existing helper fits the required behavior exactly, mention that it
     should be reused.
   - If an existing helper could be reused only after modifying its behavior,
     pause the task and ask the user before continuing. Do not silently broaden
     or repurpose existing helper semantics.
   - For each important function or helper, state:
     - what data should go in;
     - which part of the theory method it implements;
     - what data should come out.
   - If multiple functions are needed, make the output of one function clearly
     feed into the input of the next.
   - Do not prescribe every internal line of code unless that detail is an
     intentional invariant that future refactors must preserve.

3. Invariants and edge cases.
   - Explicitly list conventions that must not change silently, such as
     averaging order, normalization, sector labels, seed behavior, or whether a
     quantity should be computed before or after trajectory averaging.
   - State important invalid inputs or parameter regimes and how code should
     fail or warn.
   - Mention when an existing helper or instruction file is authoritative for a
     convention.

## Core Rules

- If the user asks for a new theory or instruction file in TeX format, make it renderable by itself.
- For standalone TeX files, use only the minimal LaTeX preamble needed to compile cleanly, for example `article`, `geometry`, `amsmath`, `amssymb`, and `hyperref` when useful.
- After creating or editing a standalone TeX file, compile it with `docs/build_tex.sh` when possible and fix render errors before finishing.
- Files may reference other documentation files. Prefer references over re-explaining material that already has a dedicated source-of-truth file.
- If a topic already has a theory file and an instruction file, keep the split clear: theory files explain equations and physics assumptions, while instruction files explain how the repository should implement or use them.
- Always update AGENTS.md file accordingly.
- Before adding a long explanation, check whether an existing file already covers
  the idea. If it does, briefly state the local relevance and refer to that file.
  This keeps the documentation maintainable and avoids conflicting versions of
  the same rule or formula.

## Core Rules: Instruction files

- Instruction files should be concise and straight to the point. 
- Use numbered sections and subsections deliberately so the file is easy to
  scan and easy to reference. Each section should cover one clear part of the
  implementation logic, and subsections should be used when a section contains
  multiple distinct cases or steps.
- Use notation consistent with the relevant theory files whenever possible.
  If implementation notation must differ from theory notation, state the mapping
  explicitly.
- Use instruction-style wording for implementation rules: prefer "should be" or
  "should" when stating intended behavior. Use "is" for descriptive statements,
  such as where a rule or reference file is found.

## Core Rules: Theory files

TDOO: add instructions on how to write theory files

## Reference Style

When pointing to another document, use the repository path, for example:

```text
docs/instructions/bloch_vector_averaging.tex
```

For TeX files, use escaped underscores inside `\texttt{...}`:

```tex
\texttt{docs/instructions/bloch\_vector\_averaging.tex}
```
