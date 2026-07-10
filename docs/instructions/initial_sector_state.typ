#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Initial Sector State: Implementation Instructions]]

= Purpose

This file describes the initial strong-symmetry sector expansion used by the
MCWF ensemble run. Use it when editing `solvers/mcwf/sim.py` or helper
functions that construct the initial sector coefficients and sector blocks.

= Method

The initial wavefunction should first be defined over strong-symmetry sectors.
The MCWF simulation should then evolve one internal array per sector in the
active-manifold excitation basis.

== Sector Keys

Single-group runs use scalar sector keys:

```python
{NJ: coeff}
```

with selected sectors

$
N_(J) in {floor(N / 2) - Delta N, ..., floor(N / 2) + Delta N}.
$

Group-resolved runs use tuple sector keys:

```python
{(NJ1, NJ2): coeff}
```

where

$
N_(J,1) + N_(J,2) = N_(J),
0 <= N_(J,1) <= N_1,
0 <= N_(J,2) <= N_2.
$

For future multi-group runs, the tuple key should generalize to
$(N_(J,1), ..., N_(J,G))$.

== Sector Coefficients

High-level initialization helpers should be preferred over manually writing
low-level sector dictionaries. Supported sector coefficient choices should be:

```python
sector_distribution="square"    # equal total weight over selected NJ sectors
sector_distribution="binomial"  # product-state binomial weights, truncated to the selected window
```

== Sector Blocks

After the sector coefficients are chosen, the propagated state should be stored
as sector blocks:

```python
{NJ: psi_NJ}
{(NJ1, NJ2): psi_NJ1_NJ2}
```

For a scalar sector key, `psi_NJ` should have shape `(NJ + 1,)` over
$n_(e) = 0, ..., N_(J)$. For a two-group tuple key, `psi_NJ1_NJ2` should have
shape `((NJ1 + 1) * (NJ2 + 1),)` over $(n_(e,1), n_(e,2))$.

By default, each internal sector state should start in the all-active-atoms-down
state, i.e. $n_(e)=0$ or $(n_(e,1), n_(e,2))=(0,0)$.

= Implementation Flow

Initial sector construction should be part of the MCWF ensemble run, not a
separate notebook-level step:

```python
parameters = MCWFSolverParameters(...)

run_trajectory_ensemble(parameters, *, t_eval, ntraj, ...)
    -> centered_sector_initial_coeffs(...)
    -> check_initial_sector_omega_ratio(...)
    -> build_precomputed_trajectory_data(...)
    -> _simulate_single_trajectory(...)
    -> TrajectoryEnsemble(...)
```

The ensemble entry point should choose the sector coefficients, validate the
initial sector regime, build reusable precomputed data, and pass the resulting
sector metadata into each trajectory.

= Invariants

- Initial sector construction should be deterministic for fixed `Ni`, `dN`,
  `omega_i`, and `sector_distribution`.
- The represented sector keys should be the same for every trajectory in an
  ensemble.
- The internal sector-block basis should be fixed by the sector key and should
  not change during a trajectory.
