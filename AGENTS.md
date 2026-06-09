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

### `docs/paper_inhomogeneous_couplings.md`

Read this whenever the task involves inhomogeneous couplings, weighted collective jumps, split active-manifold groups, or comparing homogeneous vs inhomogeneous MCWF simulations:

* Two-group active-manifold sectors ((N_{J,1},N_{J,2})).
* Product Dicke basis states (|n_{e,1},n_{e,2}\rangle).
* Group-resolved operators (J_{1,\pm}), (J_{2,\pm}), (N_{e,1}), (N_{e,2}).
* Drive replacement
  [
  \Omega J_x \rightarrow \Omega(\omega_1J_{1,x}+\omega_2J_{2,x}).
  ]
* Input convention: pass (N_1), (N_2), and (\omega_1), then set (\omega_2=(N-N_1\omega_1)/N_2) so \(N_1\omega_1+N_2\omega_2=N\).
* Weighted collective jump
  [
  A=\omega_1J_{1,-}+\omega_2J_{2,-}.
  ]
* Shifted weighted jump
  [
  l=A+i\frac{\Omega}{\Gamma}.
  ]
* Performance requirements: precompute two-group operators, jump operators, generators, and propagators; do not rebuild sparse Kronecker products in the MCWF time loop.
* Validation requirements: when one group is empty, recover the current homogeneous one-group implementation within numerical tolerance. The homogenous code path can be changed to run the inhomogenous code path, only if it can be guranteed that the runtime won't increase for e.g. NJ1 = 0.

Use this before modifying Hamiltonian construction, jump operators, sector keys, basis construction, precomputation, observables, or trajectory propagation for inhomogeneous coupling. Keep the homogeneous code path working unchanged.

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
