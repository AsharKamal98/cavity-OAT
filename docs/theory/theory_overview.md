# theory_overview.md — Project theory for Coding Agents

This is the main entry point for theory-specific guidance. Use these files to understand
what equations, approximations, and physical assumptions when constructing
instruction files, or to get the physics context while implementing code changes.

TODO: add a high-level explanation of what the paper is doing. 

## Source-of-Truth Paper Files

These files are extracted from, or are direct copies of, the paper TeX sources.
They are the ultimate source of truth for paper equations, approximations, and
physics assumptions.

### `docs/theory/main.tex`

This is the main-paper TeX source after splitting the end matter into
`docs/theory/end_matter.tex`.

TODO: list topics it covers and when to read. 

### `docs/theory/end_matter.tex`

This is the end-matter section extracted from `docs/theory/main.tex`.

Read this whenever the task involves main-paper end-matter derivations or
benchmarks.

Use this for:

- Berry-phase Hamiltonian derivation from \(N_J\) fluctuations;
- collective dephasing from the shifted jump operator;
- quantum-trajectory benchmarking of squeezing dynamics.


### `docs/theory/appendix_cavity_model.tex`

This is the cavity-model appendix extracted from `docs/theory/supp.tex`.

Read this whenever the task mentions cavity parameters, the cavity model,
bad-cavity elimination, or converting experimental/cavity parameters into the
effective spin parameters.

Use this for:

- full driven cavity-spin model;
- cavity mean-field equation;
- bad-cavity condition;
- identification of \(\Omega\) and \(\Gamma\);
- relation between cavity coherence and the shifted jump operator.

### `docs/theory/appendix_mean_field_theory.tex`

This is the mean-field-theory appendix extracted from `docs/theory/supp.tex`.

Read this whenever the task involves the two-level mean-field steady state,
the Berry-frame frequency, or the steady-state angles \(\theta_J,\phi_J\).

Use this for:

- Schwinger-boson mean-field equations;
- identifying \(\omega_B\);
- steady-state relations for \(\theta_J,\phi_J\);
- the \(\delta=0\) critical-drive relation.

### `docs/theory/appendix_coherence_preservation.tex`

This is the coherence-preservation appendix extracted from
`docs/theory/supp.tex`.

Read this whenever the task involves weak versus strong symmetries, preserved
coherences, dissipative freezing, or coherence transfer after the drive is
turned off.

Use this for:

- weak and strong symmetry distinctions;
- Liouvillian block structure;
- mean-field preservation of \(u^\dagger d\) and \(u^\dagger e\);
- analytic coherence transfer when \(\Omega=0\).

### `docs/theory/appendix_weak_drive_limit.tex`

This is the weak-drive-limit appendix extracted from `docs/theory/supp.tex`.

Read this whenever the task involves adiabatic elimination of \(|e\rangle\),
ground-manifold effective Hamiltonians, weak-drive jump operators, or
collective interactions from cavity detuning.

Use this for:

- non-Hermitian perturbative elimination;
- effective \(|\uparrow\rangle,|\downarrow\rangle\)-manifold Hamiltonian;
- effective weak-drive jump operator;
- weak-drive OAT and collective dephasing rates.

### `docs/theory/appendix_holstein_primakoff.tex`

This is the Holstein-Primakoff appendix extracted from `docs/theory/supp.tex`.

Read this whenever the task involves the generalized HP expansion, dressed
modes \(c,s,\jmath\), effective \(S\)-dynamics, or adiabatic elimination of the
\(J\)-fluctuation mode.

Use this for:

- the \(c,s,\jmath\) mode transformation;
- HP quadrature-to-spin mappings;
- quadratic HP Hamiltonian and linear jump operator;
- adiabatic elimination leading to effective OAT and dephasing.

### `docs/theory/appendix_single_particle_decoherence.tex`

This is the single-particle-decoherence appendix extracted from
`docs/theory/supp.tex`.

Read this whenever the task involves squeezing limits from spontaneous
emission, single-particle dephasing, competing decoherence sources, or
comparison with other cavity squeezing protocols.

Use this for:

- optimal OAT squeezing scalings;
- collective-decoherence-modified squeezing scalings;
- spontaneous-emission and dephasing limits;
- comparison table for common cavity squeezing protocols.

### `docs/theory/supp.tex`

This is the full supplemental TeX source. Read it only when the existing theory
summary files do not contain enough context for a supplemental/appendix-level
task, or when the user explicitly asks to inspect the supplement.

Prefer targeted theory summaries, such as
`docs/theory/appendix_cavity_model.tex`, before reading the full
supplement.

### `docs/theory/Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong.pdf`

This is the original paper PDF. Use the TeX and Markdown files first. Open the
PDF only if the TeX/Markdown sources are insufficient or if the user explicitly
requests PDF-level verification.

## Project Theory Notes

These files contain theory notes derived or written within this project. They
are useful for project-specific extensions and working derivations, but they are
not the ultimate source of truth when they overlap with the actual paper files
above.

If a project note conflicts with a source-of-truth paper file, trust the paper
file first and flag the inconsistency before implementing changes.

### `docs/theory/notes_inhomogeneous_coupling.tex`

This is the project theory note for inhomogeneous couplings.

Read this whenever the task involves the physics of inhomogeneous couplings,
weighted drives, split active-manifold groups, or group-resolved mean-field
equations.

Use this for:

- inhomogeneous Hamiltonian structure;
- group-resolved Dicke/occupation basis;
- split Schwinger-boson operators;
- weighted collective jump operator;
- group-resolved mean-field equations and residuals.


## General rules
- Theory files should not be treated as implementation instructions unless they
  explicitly say so.
- Source-of-truth paper files override project theory notes if the two disagree.
