# AGENTS.md — Project Instructions for Coding Agents

This repository studies and simulates the protocol from:

`docs/papers/Engineering One Axis Twisting via a Dissipative Berry Phase Using Strong.pdf`

The core topic is engineered one-axis twisting (OAT) from dissipative Berry phases in a driven-dissipative three-level atom-cavity system with a strong symmetry.

## How to use the paper notes

Use the Markdown notes below instead of trying to read the full PDF for every task.

### `docs/paper_summary.md`

Read this for high-level context:
- What the protocol does.
- The three phases of the protocol.
- Why strong symmetry matters.
- What Figs. 3 and 4 are about.
- What the main numerical implementation is trying to reproduce.

Use this for general orientation, summaries, plotting interpretation, and README-style text.

### `docs/paper_key_equations.md`

Read this whenever the task involves formulas, implementation, observables, or physics logic:
- Master equation.
- Collective operators.
- Mean-field steady state.
- Critical drive.
- Effective `S`-spin basis.
- Berry-phase/OAT Hamiltonian.
- Collective dephasing jump operator.
- Generalized squeezing parameter.
- Fig. 4 benchmark parameters.

Use this for coding tasks involving Hamiltonians, jump operators, phase protocols, squeezing, or comparison to analytics.

### `docs/paper_theory_notes.md`

Read this when writing longer theory explanations or when the task requires conceptual accuracy:
- Strong vs weak symmetry.
- Why coherences between `N_J` sectors survive.
- Why detuning causes collective dephasing.
- How the effective OAT picture emerges.
- How single-particle decoherence enters.

Use this for drafting text, documentation, theory summaries, or explaining results.

### `docs/paper_appendix_notes.md`

This file is intentionally referenced but not yet created.

When it exists, use it only for detailed Appendix-level derivations:
- Berry phase derivation.
- Collective dephasing derivation.
- Benchmarking squeezing dynamics.
- Supplemental derivations.
- Holstein-Primakoff approximation.
- Single-particle decoherence scaling.

Do not assume this file exists yet.

## General rules

- Do not guess paper-specific formulas from memory. Check the relevant Markdown notes first.
- For pure coding tasks unrelated to physics formulas, do not read all paper notes unless needed.
- For any task involving squeezing, jump operators, Berry phase, `N_J` sectors, or Fig. 4 reproduction, read `docs/paper_key_equations.md` first.
- Keep implementations consistent with the reduced strong-symmetry sector basis used in the code.
- Prefer sparse/reduced-basis constructions over full `3^N` tensor-product operators.
