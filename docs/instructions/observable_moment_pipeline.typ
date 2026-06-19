#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Observable Moment Pipeline]
]

= Purpose

Use this note when changing trajectory observable extraction, ensemble
post-processing, diagnostic data builders, or plotting functions that consume
trajectory moments.

The goal is to keep plotting thin and to avoid recomputing trajectory moments
inside every plot. Plotting functions should preferably consume already-built
observable or diagnostic data.

= Core Order

For ensemble diagnostics whose target is the unconditioned MCWF state, the
implementation order should be:

```python
TrajectoryEnsemble
-> compute required moments per trajectory per timestep
-> average those moments across trajectories on the shared t_eval grid
-> construct nonlinear derived quantities from the averaged moments
-> plot already-computed data
```

Do not average trajectory wavefunctions or sector-block states. Do not average
nonlinear per-trajectory outputs such as angles, normalized directions,
covariance eigenvalues, squeezing parameters, or vector lengths unless the
diagnostic is explicitly intended to describe typical conditioned trajectories.

= First-Order Observable Pipeline

Existing first-order trajectory observables include

$
J_x, J_y, J_z, N_e, N_J, "jump rate",
$

with optional two-group moments

$
J_(x,g), J_(y,g), J_(z,g), N_(e,g), quad g=1,2.
$

The intended function split is:

```python
collect_trajectory_observables(traj) -> ObservableSeries
```

Compute first-order moments for one trajectory on its saved `t_eval` grid.
Current implementation reference:

```python
quantum_trajectories.aggregator.trajectory_observables(...)
```

```python
collect_ensemble_observable_samples(ensemble) -> list[ObservableSeries]
```

Compute the same first-order moment series for every trajectory. All entries
must share the same `t_eval` grid.

```python
average_observable_samples(samples) -> ObservableSeries
```

Average first-order moments across trajectories. The averaged
`ObservableSeries` is the right input for active-manifold angles, spin
directions, group-angle diagnostics, and most first-order plots.

Current implementation reference:

```python
quantum_trajectories.aggregator.ensemble_observables(...)
```

The current implementation combines sample collection, averaging, and derived
Bloch quantities in one helper. Future refactors may split these steps, but the
averaging order should not change.

= Bloch Directions and Angles

After first-order moments have been averaged, construct active-manifold
quantities:

```python
build_active_bloch_quantities(averaged_observables)
    -> N_active, sx, sy, sz, theta_J, phi_J
```

The normalization and angle convention should follow:

```text
docs/instructions/bloch_vector_averaging.typ
```

Do not compute `sx`, `sy`, `sz`, `theta`, or `phi` separately for each
trajectory and then average them for an ensemble plot.

= Squeezing Moment Pipeline

Generalized squeezing requires staged moment extraction because the operators
used later depend on angles built from earlier averaged moments.

```python
average first-order J moments
-> theta_J(t), phi_J(t)
```

Use averaged active-manifold moments to construct the dressed $J$ direction.

```python
collect_s_moments(ensemble, theta_J, phi_J)
    -> per-trajectory Sx, Sy, Sz, N_c
average_s_moments(...)
    -> averaged S moments, theta_S, phi_S
```

The $S$ moments are moments of an effective dressed pseudospin. The $S$
operators depend on the $J$ angles because the dressed basis depends on
$theta_J(t), phi_J(t)$.

```python
collect_covariance_moments(ensemble, theta_J, phi_J, theta_S, phi_S)
    -> per-trajectory <O_a>, <O_a O_b>, N_c
average_covariance_moments(...)
    -> averaged <O_a>, <O_a O_b>, N_c
build_squeezing_data(...)
    -> covariance matrix, eigenvalues, xi^2
```

Second-order moments such as $chevron.l O_a O_b chevron.r$ must be computed
from trajectory states or sector blocks. They cannot be reconstructed from only
the first-order moments $chevron.l O_a chevron.r$.

Current implementation references:

```python
quantum_trajectories.squeezing.generalized_squeezing_for_ensemble(...)
quantum_trajectories.squeezing.generalized_squeezing_for_inhomogeneous(...)
```

= Plotting Inputs

Plotting functions should prefer one of these inputs:

- an averaged `ObservableSeries` for first-order plots;
- a diagnostic data dictionary for derived diagnostics;
- an explicit moment bundle for diagnostics that need reusable first- or
  second-order moments.

Plotting functions may accept a full `TrajectoryResult` or `TrajectoryEnsemble`
as a convenience wrapper, but they should not hide expensive recomputation if a
matching precomputed observable or diagnostic object is supplied.

= Invariants

- Ensemble moment averaging should happen on the shared saved `t_eval` grid.
- Moment extraction should be per trajectory; wavefunctions should not be
  averaged into one wavefunction.
- Nonlinear quantities should be built after the moments they depend on have
  been averaged.
- Moment extraction should be demand-driven. Do not compute second-order
  moments for ordinary first-order plots.
- Group-resolved inhomogeneous data should use the same order, with tuple
  sector keys and group moments handled before nonlinear group quantities are
  built.
