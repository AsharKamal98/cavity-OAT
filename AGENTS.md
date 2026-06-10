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

This directory contains physics/theory context. Use these files to understand
what equations, approximations, and physical assumptions should be implemented.

Theory files answer questions like:

- What is the Hamiltonian or jump operator?
- What approximation is being used?
- What are the relevant paper formulas?
- What physics should a diagnostic or plot represent?

Theory files should not be treated as implementation instructions unless they
explicitly say so. After reading theory, check `docs/instructions` for how the
code should realize that theory.

### `docs/instructions`

This directory contains implementation rules for Codex and future coding tasks.
Use these files when changing code behavior, APIs, diagnostics, averaging logic,
or post-processing logic.

Instruction files answer questions like:

- How should a helper function behave?
- What ordering or averaging convention should the code use?
- What logic must be preserved across future refactors?
- What implementation details are intentional and should not be changed
  accidentally?

If theory and implementation instructions both exist for a topic, read both:
the theory explains the equations, while the instruction file explains how the
repository should implement them.

### `docs/output_texts`

This directory contains generated writeups, method summaries, captions, and
result descriptions requested by the user. These are outputs, not source-of-truth
instructions.

Use these files only when the task asks for prose consistency, figure captions,
result summaries, or wording that should match previous text. Do not use them as
the primary source for formulas or implementation behavior.

## Theory Files

### `docs/theory/paper_summary.md`

Read this for high-level paper context:

- the main idea of the protocol;
- the three phases of the protocol;
- why strong symmetry matters;
- how Berry-phase accumulation leads to OAT;
- what the simulations are broadly trying to reproduce.

Use this for orientation, README-style summaries, plotting interpretation, and
short explanations that do not require detailed equations.

### `docs/theory/paper_key_equations.md`

Read this whenever the task involves paper formulas, observables, Hamiltonians,
jump operators, squeezing, or benchmark parameters.

Use this for:

- master equation and effective spin model;
- collective operators;
- critical drive and mean-field steady state;
- shifted jump operator;
- effective OAT and collective dephasing expressions;
- generalized squeezing parameter;
- Fig. 4-style benchmark parameters.

For any task involving squeezing, jump operators, Berry phase, \(N_J\) sectors,
or Fig. 4 reproduction, read this file first unless a more specific theory file
below is clearly the better starting point.

### `docs/theory/paper_appendix_cavity_model.tex`

Read this whenever the task mentions cavity parameters, the cavity model,
bad-cavity elimination, or converting experimental/cavity parameters into the
effective spin parameters.

Use this for:

- full driven cavity-spin model;
- cavity mean-field equation;
- bad-cavity condition;
- identification of \(\Omega\) and \(\Gamma\);
- relation between cavity coherence and the shifted jump operator.

### `docs/theory/paper_inhomogeneous_couplings.tex`

Read this whenever the task involves inhomogeneous couplings, weighted collective
jumps, split active-manifold groups, or comparing homogeneous vs inhomogeneous
MCWF simulations.

Use this for:

- two-group active-manifold sectors \((N_{J,1},N_{J,2})\);
- product Dicke basis states \(|n_{e,1},n_{e,2}\rangle\);
- group-resolved operators \(J_{1,\pm}\), \(J_{2,\pm}\), \(N_{e,1}\),
  \(N_{e,2}\);
- weighted drive Hamiltonian;
- weighted collective jump operator;
- shifted weighted jump operator;
- inhomogeneous mean-field residual equations.

Before modifying Hamiltonian construction, jump operators, sector keys, basis
construction, precomputation, observables, or trajectory propagation for
inhomogeneous coupling, read this file and then check the relevant files in
`docs/instructions`.

### `docs/theory/main.tex`

This is the full main-paper TeX source. Read it only when the existing theory
summary files do not contain enough context for the task, or when the user
explicitly asks to inspect the paper text.

Prefer the targeted theory summaries above first. They are faster to use and are
intended to avoid unnecessary full-paper reading.

### `docs/theory/supp.tex`

This is the full supplemental TeX source. Read it only when the existing theory
summary files do not contain enough context for a supplemental/appendix-level
task, or when the user explicitly asks to inspect the supplement.

Prefer targeted theory summaries, such as
`docs/theory/paper_appendix_cavity_model.tex`, before reading the full
supplement.

### `docs/theory/Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong.pdf`

This is the original paper PDF. Use the TeX and Markdown files first. Open the
PDF only if the TeX/Markdown sources are insufficient or if the user explicitly
requests PDF-level verification.

## Instruction Files

### `docs/instructions/documentation_file_guidelines.md`

Read this before creating a new file in `docs/theory` or `docs/instructions`.

Use this for:

- deciding whether a new documentation file should be Markdown or standalone
  TeX;
- keeping standalone TeX files minimal and renderable;
- avoiding duplicated explanations by referencing existing documentation files;
- preserving the split between theory files and implementation-instruction
  files.

### `docs/instructions/bloch_vector_averaging.tex`

Read this whenever a task involves averaging or plotting active-manifold Bloch
vectors, Bloch angles, dressed-state directions, or group/sector/trajectory
averages of \(\theta_J,\phi_J\).

Use this for:

- why raw collective vectors from different \(N_J\) sectors should not be
  compared directly as directions;
- how `active_manifold_angles(...)` normalizes by \(N_{\rm active}\);
- how to average sector, group, and trajectory moments before constructing
  Bloch directions;
- when this logic should not be used, e.g. genuinely extensive observables such
  as total \(\langle J_x\rangle\), jump rates, or jump counts.

### `docs/instructions/ensemble_simulation_implementation.md`

Read this whenever a task involves `run_trajectory_ensemble(...)`,
`simulate_single_trajectory(...)`, `build_precomputed_trajectory_data(...)`,
the shared `t_eval` grid, MCWF step splitting, multiprocessing, precomputed
propagators, or full-step versus partial-step propagation.

Use this for:

- how ensemble-level precomputation is structured;
- what data are shared across trajectories and worker processes;
- when precomputed full-step propagators can be reused;
- why phase boundaries, `t_eval` boundaries, and jump bisection require
  variable-step propagation;
- how homogeneous and inhomogeneous sector structures differ at the
  implementation level.

### `docs/instructions/squeezing.tex`

Read this whenever a task involves implementing or modifying the generalized
three-level squeezing parameter.

Use this for:

- how the dressed states \(|1\rangle\), \(|c\rangle\), \(|j\rangle\), and
  \(|s\rangle\) should be constructed;
- what moments are required;
- how ensemble squeezing should be computed from averaged moments rather than
  averaging per-trajectory squeezing values;
- what data should be saved or post-processed.

When this file refers to active-manifold angle averaging, use
`docs/instructions/bloch_vector_averaging.tex` as the authoritative convention.

### `docs/instructions/dephasing_diagnostics.tex`

Read this whenever a task involves visualizing dephasing, plotting effective
\(S\)-Bloch vector lengths, or comparing total and group-resolved coherence
loss.

Use this for:

- why Bloch-vector components must be ensemble-averaged before taking the
  vector length;
- how to plot total effective-\(S\) Bloch-vector length for homogeneous results;
- how to plot total, group-1, and group-2 lengths for inhomogeneous results;
- how to interpret total shrinkage versus group-resolved shrinkage.

## General Rules

- Do not guess paper-specific formulas from memory. Check the relevant theory
  file first.
- For pure coding tasks unrelated to physics formulas, do not read all theory
  notes unless needed.
- Before creating a new theory or instruction file, read
  `docs/instructions/documentation_file_guidelines.md`.
- For implementation changes, check whether a relevant instruction file exists
  before editing code.
- If both a theory file and an instruction file apply, use the theory file for
  equations and the instruction file for repository-specific implementation
  behavior.
- Keep implementations consistent with the reduced strong-symmetry sector basis
  used in the code.
- Prefer sparse/reduced-basis constructions over full \(3^N\) tensor-product
  operators.
- Do not use files in `docs/output_texts` as implementation requirements unless
  the user explicitly asks to follow one of those writeups.
