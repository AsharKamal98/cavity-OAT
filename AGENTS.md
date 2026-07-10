# AGENTS.md — Project Instructions for Coding Agents

This repository studies and simulates the protocol from the paper
*Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong Symmetries*.

The core topic is engineered one-axis twisting (OAT) from dissipative Berry
phases in a driven-dissipative three-level atom-cavity system with a strong
symmetry.

## Documentation Structure

The `docs` directory is split by purpose. Use this split before deciding what
to read.

### `docs/theory`

This directory contains the project theory. Use these files to understand
the theory when constructing instruction files, or to get the physics context
while implementing code changes. Paper-derived theory files are the ultimate
source of truth. Project theory notes are useful working derivations, but should
defer to the actual paper files if the two disagree. Instruction files should be
checked against theory files in case of inconsistencies or doubts about
implementation.

Use `docs/theory/theory_overview.md` as the entry point. That file gives the
high-level theory background and points to other theory-specific files for details.

The theory-document hierarchy should be:

```text
AGENTS.md
  -> docs/theory/theory_overview.md
      -> docs/theory/main.tex
      -> docs/theory/end_matter.tex
      -> docs/theory/appendix_cavity_model.tex
      -> docs/theory/appendix_mean_field_theory.tex
      -> docs/theory/appendix_coherence_preservation.tex
      -> docs/theory/appendix_weak_drive_limit.tex
      -> docs/theory/appendix_holstein_primakoff.tex
      -> docs/theory/appendix_single_particle_decoherence.tex
      -> docs/theory/notes_inhomogeneous_coupling.tex
      -> docs/theory/supp.tex
      -> docs/theory/Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong.pdf
```

For one-time setup or later editing of standalone theory TeX files, use the
`theory-writer` skill.

### `docs/instructions`

This directory contains high-level implementation rules for Codex and future coding tasks.
Use these files when changing code behavior, APIs, diagnostics, averaging logic,
post-processing logic, or requiring implementation-specific knowledge for a task.

For implementation-specific tasks, use
`docs/instructions/implementation_overview.md` as the entry point. That file gives the
high-level intended code structure and points to task-specific instruction files
for details. Use this for:

- the intended top-level workflow of custom MCWF runs;
- how parameter setup, phase construction, initial states, validation,
  simulation, observables, diagnostics, and plotting should fit together;
- deciding which task-specific instruction file to read next;
- keeping notebook and helper code aligned with the repository's intended
  structure.

The instruction-document hierarchy should be:

```text
AGENTS.md
      -> docs/instructions/implementation_overview.md
      -> docs/instructions/common/utils.typ
      -> docs/instructions/common/plotting.typ
      -> docs/instructions/model_parameters.typ
      -> docs/instructions/parser.typ
      -> docs/instructions/solvers/mcwf/method_parameters.typ
      -> docs/instructions/solvers/mcwf/initial_sector_state.typ
      -> docs/instructions/parameter_validation.md (planned/missing)
      -> docs/instructions/solvers/mcwf/simulation_precompute.typ
      -> docs/instructions/solvers/mcwf/sector_operators.typ
      -> docs/instructions/paper_inhomogeneous_couplings.typ
      -> docs/instructions/solvers/mfe/mfe-solver.typ
      -> docs/instructions/solvers/mcwf/ensemble_simulation.typ
      -> docs/instructions/solvers/mcwf/single_trajectory_simulation.typ
      -> docs/instructions/observable_moment_pipeline.typ
      -> docs/instructions/j_moments.typ
      -> docs/instructions/post_analysis/mfe_residuals.typ
      -> docs/instructions/bloch_vector_averaging.typ
      -> docs/instructions/plot_spin_components.typ
      -> docs/instructions/post_analysis/squeezing.typ
      -> docs/instructions/dephasing_diagnostics.typ
      -> docs/instructions/post_analysis/post_analysis.md
      -> docs/instructions/plotting_workflows.md
```

Keep task-specific routing in `docs/instructions/implementation_overview.md`,
not directly in `AGENTS.md`.

All future changes to intended code behavior should be reflected in the relevant instruction file.
The goal is that, given the implementation overview and the
task-specific instruction files it references, Codex can reconstruct the
high-level code logic without relying on hidden assumptions.





### `docs/output_texts`

This directory contains generated writeups, method summaries, captions, and
result descriptions requested by the user. These are outputs, not source-of-truth
instructions.

Use these files only when the task asks for prose consistency, figure captions,
result summaries, or wording that should match previous text. Do not use them as
the primary source for formulas or implementation behavior.

Read relevant theory and instruction files before writing text. 

## General Rules
- For implementation changes, or any task requiring knowledge about implementation, read
  `docs/instructions/implementation_overview.md` first.
- For editing instructions, or any task requiring knowledge about theory, read
  `docs/theory/theory_overview.md` first.
- Always read the relevant overview file before reading a task-specific
  instruction or theory-specific file. 
- Before creating a new theory or instruction file, use the relevant writing
  skill when available, such as `theory-writer` for standalone theory TeX
  files or `agent-instruction-writer` for instruction files.
- If both a theory file and an instruction file apply, use the theory file for
  equations and the instruction file for repository-specific implementation
  behavior.
- Keep `AGENTS.md` as the clean hierarchy/front-door map. It should point to
  overview files, and the hierarchy should list task-specific or theory-specific
  files by name under the overview file that references them.
  Example:
  ```text
  AGENTS.md
    -> docs/instructions/implementation_overview.md
        -> docs/instructions/task_specific_file.md
        -> docs/instructions/other_task_specific_file.md
          -> ...
        -> ...
    -> docs/theory/theory_overview.md
        -> docs/theory/task_specific_theory_file.tex
  ```
