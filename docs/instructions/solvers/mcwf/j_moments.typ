#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[MCWF J-Moment Extraction]
]

= Purpose
This file specifies how to compute first-order $J$-moments in
`solvers/mcwf/j_moments.py`. Use this file for tasks related to
MCWF J-moment extraction or when writing diagnostics that consume its outputs.

Shared spin-vector direction and angle conventions are defined in
`docs/instructions/common/spin_moment_conventions.typ`.

The main function structure should be:

```python
_compute_snapshot_j_moments(snapshot, ...)
    -> JMomentSnapshot:
    # Compute J moments for one saved snapshot of one trajectory.
    return JMomentSnapshot(
        t, integration_phase_index,
        x, y, z, N_e, N_j,
        aggregated=AggregatedJMomentSnapshot(
            sum(x), sum(y), sum(z), sum(N_e), sum(N_j)
        ),
        jump_rate,
    )

compute_trajectory_j_moments(
    trajectory: TrajectoryResult,
    metadata: TrajectoryEnsembleMetadata,
    *,
    tol=1e-12,
)
    -> JMomentSeries:
    
    j_moment_snapshots = [
        _compute_snapshot_j_moments(
            snapshot,
            metadata.phase_protocol.integration_phases,
            metadata.Gamma,
            ...
            )
        for snapshot in trajectory.snapshots
    ]
    return stack into JMomentSeries

compute_mcwf_j_moments(ensemble: TrajectoryEnsemble, ...)
    -> JMomentSeries

    samples = map_with_optional_pool(
        compute_trajectory_j_moments(traj, ensemble.metadata)
        for traj in ensemble.trajectories
    )
    averaged = compute_average_j_moments(samples)
    JMomentSeries.attach_norm_spin_components_from_spin_components(averaged, tol=tol)
    JMomentSeries.attach_angles_from_norm_spin_components(averaged, tol=tol)
    JMomentSeries.attach_aggregated_from_group_fields(averaged, tol=tol)
    return averaged
  
compute_average_j_moments(samples: list[JMomentSeries])
    -> JMomentSeries
    "average raw J-moment samples across trajectories"
```

`_compute_snapshot_j_moments(...)` computes the J moments for one snapshot of
one trajectory and returns group tuples for both scalar and tuple sector keys.


`compute_trajectory_j_moments(...)` calls `_compute_snapshot_j_moments(...)` once per saved snapshot.

`compute_mcwf_j_moments(...)` calls `compute_trajectory_j_moments(...)` once
per trajectory, then uses `compute_average_j_moments(...)` to average the raw
moment fields across the shared saved-time grid.

`compute_mcwf_j_moments(...)` uses the `JMomentSeries` class methods to attach
group spin directions and angles from the averaged spin components, then
constructs all aggregated fields from those averaged raw group fields. The
`JMomentSeries` class methods are covered in `docs/instructions/parser.typ`.

= J-Moment Definitions

This section specifies the definitions of the $J$-sphere moments which should be computed for each trajectory snapshot.

At a saved time $t_(k)$, the trajectory stores sector blocks

$
psi(t_(k)) = sum_(alpha) psi_(alpha)(t_(k)),
$

where $alpha$ labels a sector (e.g. $N_(J)$ for a homogeneous sector or
$(N_(J,1), ..., N_(J,G))$ for a group-resolved sector). The implementation should use the block norm

$
Z(t_(k)) = sum_(alpha) chevron.l psi_(alpha)(t_(k)) | psi_(alpha)(t_(k)) chevron.r
$

to normalize all expectation values. For each requested moment operator
$O_(alpha)$,

$
chevron.l O(t_(k)) chevron.r =
frac(
  sum_(alpha) chevron.l psi_(alpha)(t_(k)) | O_(alpha) | psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
).
$

Apply this rule to the group-resolved fields. For active-sector atom number
$N_(J,g)$, do not construct an operator. Use the scalar group sector number and
replace each numerator term by
$N_(J,g)(alpha) chevron.l psi_(alpha) | psi_(alpha) chevron.r$. For
strong-symmetry-preserving dynamics, this weighted sector number is expected
to remain constant in time.

== Jump Rate

The jump-rate field should use the same phase-local jump operator used by the
MCWF solver. Let $l_(alpha)(t_(k))$ be this unscaled sector jump operator, built
from the precomputed sector operators in
`docs/instructions/solvers/mcwf/simulation_precompute.typ` by
`build_phase_jump_operator_for_sector(...)`. The physical rate is

$
r(t_(k)) =
frac(
  Gamma sum_(alpha)
    chevron.l l_(alpha)(t_(k)) psi_(alpha)(t_(k)) |
    l_(alpha)(t_(k)) psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
).
$

Tiny negative rates caused by floating-point noise should be clipped to zero.
The jump rate should not be normalized by $N$, $N_(J)$, or the number of
trajectories.

== Group and Aggregated Fields

Unsuffixed fields are the primary group-resolved representation and should be
returned for every group $g$, including for $G=1$. For a group-local operator
$O_(g,alpha)$,

$
chevron.l O_g(t_(k)) chevron.r =
frac(
  sum_(alpha) chevron.l psi_(alpha)(t_(k)) |
    O_(g,alpha) |
    psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
).
$
Use this rule for $J_(i,g)$ with $i in {x,y,z}$ and for $N_(e,g)$ with the
group-local operators from `build_sector_ops_for_key(...)`.

The group-resolved active-sector atom number `N_j[g]` should use the same
sector-weighted average, with $N_(J,g)(alpha)$ for each group.

For homogeneous or single-group results, every group field should be a
one-element tuple. The current implementation supports homogeneous scalar
sector keys and group-resolved tuple sector keys; mixed or malformed sector
structures should not be introduced silently.

Full-system fields live in the nested `aggregated` container. Additive fields
are constructed from the group-resolved fields:

$
q_"aggregated"(t_(k)) = sum_g q_g(t_(k)),
quad q in {J_x, J_y, J_z, N_e, N_J}.
$

The snapshot subcontainer stores these sums directly. For an ensemble series,
average the raw group fields first and then call
`JMomentSeries.attach_aggregated_from_group_fields(...)`. Do not independently
evaluate or average full-system operators. Nonlinear fields such as
`aggregated.nx` and `aggregated.theta` are derived afterward from the summed
aggregated spin vector.

== J-Vector Direction and Angles

Apply `docs/instructions/common/spin_moment_conventions.typ` to the
trajectory-averaged J-sphere components, using $S_a=J_a$ and the J-sphere pole
choice $p=-1$.

For each group $g$, compute `length[g]`, `nx[g]`, `ny[g]`, `nz[g]`,
`theta[g]`, and `phi[g]`, corresponding to
$J_("len",g)$, $n_(x,g)$, $n_(y,g)$, $n_(z,g)$,
$theta_(J,g)$, and $phi_(J,g)$. Compute the same six fields in the nested
`aggregated` container from the already-summed aggregated J vector.

Their mathematical definitions, averaging order, pole convention, and
zero-length behavior are defined in
`docs/instructions/common/spin_moment_conventions.typ`.

Use the corresponding `JMomentSeries` class methods documented in
`docs/instructions/parser.typ`.

= Output

The snapshot helper should return a `JMomentSnapshot` with tuple-valued group
fields and scalar aggregated fields:

```python
JMomentSnapshot(
    t=t,
    integration_phase_index=integration_phase_index,
    x=tuple[float, ...],
    y=tuple[float, ...],
    z=tuple[float, ...],
    N_e=tuple[float, ...],
    N_j=tuple[float, ...],
    aggregated=AggregatedJMomentSnapshot(
        x=float,
        y=float,
        z=float,
        N_e=float,
        N_j=float,
    ),
    jump_rate=float,
)
```

The trajectory helper should stack one `JMomentSnapshot` per saved snapshot and
return `JMomentSeries`:

```python
JMomentSeries(
    t,
    integration_phase_index=None or array,
    x=None or tuple[array, ...],
    y=None or tuple[array, ...],
    z=None or tuple[array, ...],
    length=None or tuple[array, ...],
    nx=None or tuple[array, ...],
    ny=None or tuple[array, ...],
    nz=None or tuple[array, ...],
    N_e=None or tuple[array, ...],
    N_j=None or tuple[array, ...],
    theta=None or tuple[array, ...],
    phi=None or tuple[array, ...],
    aggregated=AggregatedJMomentSeries(
        x=None or array,
        y=None or array,
        z=None or array,
        length=None or array,
        nx=None or array,
        ny=None or array,
        nz=None or array,
        N_e=None or array,
        N_j=None or array,
        theta=None or array,
        phi=None or array,
    ),
    jump_rate=None or array,
)
```

All trajectory arrays should be defined on the trajectory's saved `t_eval`
grid. `JMomentSeries` returned by `compute_mcwf_j_moments(...)` has the
same fields, but each numeric moment field is averaged across trajectories.
For averaged outputs, the J-vector length and normalized direction fields
should be attached inside `compute_mcwf_j_moments(...)` after
`compute_average_j_moments(...)` returns the raw averaged series. Angle fields
should then be attached there from those normalized directions. Finally,
construct the complete aggregated subcontainer from the averaged group fields.
These helpers should not compute squeezing or covariance matrices.

Legacy note: older field names included `Jx`, `Jy`, `Jz`, `Jx_groups`,
`Jy_groups`, `Jz_groups`, `J_len`, `sx`, `sy`, and `sz`. The immediately
preceding API used a `_groups` suffix for group tuples and unsuffixed names for
aggregated fields.


= Invariants

- `compute_trajectory_j_moments(...)` should return per-trajectory moments only.
- `compute_average_j_moments(...)` should average moments, not nonlinear derived
  quantities.
- Aggregate raw group moments only after trajectory averaging; do not average
  per-trajectory normalized directions or angles.
- `compute_mcwf_j_moments(...)` should return trajectory-averaged moments and
  should be the preferred input for plots that do not need per-trajectory
  samples.
- `compute_mcwf_j_moments(...)` should attach the group `length`, `nx`, `ny`,
  `nz`, `theta`, and `phi` fields, plus their counterparts in `aggregated`
  from the aggregated J components.
- `compute_mcwf_j_moments(...)` should require all internally computed
  samples to share the same `t` and `integration_phase_index` grids.
