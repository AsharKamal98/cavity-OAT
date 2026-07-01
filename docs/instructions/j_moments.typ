#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[J-Moment Extraction]
]

= Purpose
This file specifies how to compute first-order $J$-moments in
`quantum_trajectories/j_moments.py`. Use this file for tasks related to
J-moment extraction or when writing diagnostics that consume its outputs.

The main function structure should be:

```python
_compute_snapshot_j_moments(snapshot, ...)
    -> JMomentSnapshot:
    # Compute J moments for one saved snapshot of one trajectory.
    return JMomentSnapshot(
        t, phase_index, x, y, z, N_e, N_j, jump_rate,
        J_drive, x_groups, y_groups, z_groups, N_e_groups,
        N_j_groups,
    )

compute_trajectory_j_moments(trajectory: TrajectoryResult, *, tol=1e-12)
    -> JMomentSeries:
    
    j_moment_snapshots = [
        _compute_snapshot_j_moments(
            snapshot,
            trajectory.phases,
            trajectory.Gamma,
            ...
            )
        for snapshot in trajectory.snapshots
    ]
    return stack into JMomentSeries

compute_ensemble_j_moments(ensemble: TrajectoryEnsemble, ...)
    -> JMomentSeries

    samples = map_with_optional_pool(
        compute_trajectory_j_moments(traj)
        for traj in ensemble.trajectories
    )
    averaged = compute_average_j_moments(samples)
    _attach_spin_direction_fields(averaged, tol=tol)
    _attach_spin_angles(averaged, tol=tol)
    return averaged
  
compute_average_j_moments(samples: list[JMomentSeries])
    -> JMomentSeries
    "average raw J-moment samples across trajectories"
```

`_compute_snapshot_j_moments(...)` computes the J-moments for one snapshot of
one trajectory. If the trajectory uses tuple sector keys, it should also return
group-resolved moment fields. It is called once per saved snapshot in
`compute_trajectory_j_moments(...)`. `compute_average_j_moments(...)` averages
the raw per-trajectory moment fields across the shared saved time grid.


Ensemble and plotting code should call `compute_ensemble_j_moments(...)`,
which collects trajectory samples, averages them, and attaches derived
direction and angle fields to the averaged result.
`_attach_spin_direction_fields(...)` attaches `length`, `nx`, `ny`, and `nz`,
plus group-resolved versions when group fields exist.
`_attach_spin_angles(...)` attaches `theta` and `phi`, plus group-resolved
versions when group direction fields exist.
MFE residuals are computed separately; use
`docs/instructions/mfe_residuals.typ` for that diagnostic.

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

Use this rule for $J_x$, $J_y$, $J_z$ and $N_(e)$ with the sector-local
operators from the precomputed `ops_list` described in
`docs/instructions/simulation_precompute.typ`. For the full wavefunction
average active-sector atom number field $N_(J)$, do not construct an operator.
Use the scalar active-sector atom number for each sector and replace each
numerator term by
$N_(J)(alpha) chevron.l psi_(alpha) | psi_(alpha) chevron.r$. For strong-symmetry-preserving dynamics, this weighted sector number is expected
to remain constant in time.

== Jump Rate

The jump-rate field should use the same phase-local jump operator used by the
MCWF solver. Let $l_(alpha)(t_(k))$ be this unscaled sector jump operator, built
from the precomputed sector operators in
`docs/instructions/simulation_precompute.typ` by
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

== Group Fields

Group-resolved fields should be returned for each group $g$ (if two or more groups used). For a group-local operator $O_(g,alpha)$,

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

The group-resolved active-sector atom number `N_j_groups[g]` should be the
same sector-weighted average as `N_j`, but using $N_(J,g)(alpha)$ for each
group.

For homogeneous / single-group results, group fields should be `None`. The current implementation supports homogeneous scalar sector keys and group-resolved tuple sector keys; mixed or malformed sector structures should not be introduced silently.

The full-system fields (`Ji`,`N_e`, `N_j` etc.) should still be
computed by the full-system rules above. Group-resolved fields are additional diagnostic outputs and should not replace the full fields.

== Drive Term

The drive field `J_drive` should use the drive-coupled operator defined by
`build_sector_ops_for_key(...)`; its construction is covered in
`docs/instructions/simulation_precompute.typ`. This field was previously called
`J_x_drive` or `Jx_drive`; future code should use `J_drive` instead. In this
file it is just another normalized sector expectation value following the rule
above.

== Normalized Spin Components

After the trajectory-averaged fields $J_i(t_(k))$ have been computed, define
the J-vector length

$
J_("len")(t_(k)) =
sqrt(J_x(t_(k))^2 + J_y(t_(k))^2 + J_z(t_(k))^2).
$

The normalized spin component for $i in {x,y,z}$ is

$
n_i(t_(k)) = frac(J_i(t_(k)), J_("len")(t_(k))).
$

If $J_("len")(t_(k))$ is below the numerical tolerance, set the normalized
components to zero. For group-resolved fields, apply the same rule separately
to each group using $J_(i,g)$ and $J_("len",g)$.
This is implemented by `_attach_spin_direction_fields(...)`.

== Angles

Angles should be computed from the normalized spin components:

$
theta(t_(k)) = arccos(-n_z(t_(k))),
$

The minus sign keeps the same active-manifold convention as the old angle
pipeline: $theta=0$ points along the south/down pole.

$
phi(t_(k)) = "atan2"(n_y(t_(k)), n_x(t_(k))).
$

For group-resolved fields, apply the same rule separately to each group's
normalized components $n_(i,g)$.
This is implemented by `_attach_spin_angles(...)`.

= Output

The snapshot helper should return a `JMomentSnapshot` with scalar fields:

```python
JMomentSnapshot(
    t, phase_index,
    x, y, z,
    N_e,
    N_j,
    jump_rate,
    J_drive,
    x_groups=None or tuple[float, ...],
    y_groups=None or tuple[float, ...],
    z_groups=None or tuple[float, ...],
    N_e_groups=None or tuple[float, ...],
    N_j_groups=None or tuple[float, ...],
)
```

The trajectory helper should stack one `JMomentSnapshot` per saved snapshot and
return `JMomentSeries`:

```python
JMomentSeries(
    t, phase_index,
    x, y, z,
    x_groups=None or tuple[array, ...],
    y_groups=None or tuple[array, ...],
    z_groups=None or tuple[array, ...],
    length=None or array,
    length_groups=None or tuple[array, ...],
    nx=None or array,
    ny=None or array,
    nz=None or array,
    nx_groups=None or tuple[array, ...],
    ny_groups=None or tuple[array, ...],
    nz_groups=None or tuple[array, ...],
    N_e,
    N_j,
    N_e_groups=None or tuple[array, ...],
    N_j_groups=None or tuple[array, ...],
    theta=None or array,
    phi=None or array,
    theta_groups=None or tuple[array, ...],
    phi_groups=None or tuple[array, ...],
    jump_rate,
    J_drive,
)
```

All trajectory arrays should be defined on the trajectory's saved `t_eval`
grid. `JMomentSeries` returned by `compute_ensemble_j_moments(...)` has the
same fields, but each numeric moment field is averaged across trajectories.
For averaged outputs, the J-vector length and normalized direction fields
should be attached inside `compute_ensemble_j_moments(...)` after
`compute_average_j_moments(...)` returns the raw averaged series. Angle fields
should then be attached there from those normalized directions. These helpers
should not compute active-manifold angles, squeezing, or covariance matrices.

Legacy note: the previous field names were `Jx`, `Jy`, `Jz`,
`Jx_groups`, `Jy_groups`, `Jz_groups`, `J_len`, and `sx`, `sy`, `sz`.


= Invariants

- `compute_trajectory_j_moments(...)` should return per-trajectory moments only.
- `compute_average_j_moments(...)` should average moments, not nonlinear derived
  quantities.
- `compute_ensemble_j_moments(...)` should return trajectory-averaged moments and
  should be the preferred input for plots that do not need per-trajectory
  samples.
- `compute_ensemble_j_moments(...)` should attach `length`, `nx`, `ny`, `nz`,
  `theta`, and `phi` from the averaged J components, plus group-resolved
  versions when group fields exist.
- `compute_ensemble_j_moments(...)` should require all internally computed
  samples to share the same `t` and `phase_index` grids.
