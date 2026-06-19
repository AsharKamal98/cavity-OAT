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
        t, phase_index, Jx, Jy, Jz, N_e, jump_rate, N_j,
        Jx_drive, Jx_groups, Jy_groups, Jz_groups, N_e_groups,
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

    samples = [
        compute_trajectory_j_moments(traj)
        for traj in ensemble.trajectories
    ]
      
    return compute_average_j_moments(samples)
  
compute_average_j_moments(samples: list[JMomentSeries])
    -> JMomentSeries
    "average J-moment samples across trajectories"
```

`_compute_snapshot_j_moments(...)` computes the J-moments for one snapshot of
one trajectory. If the trajectory uses tuple sector keys, it should also return
group-resolved moment fields. It is called once per saved snapshot in
`compute_trajectory_j_moments(...)`. `compute_average_j_moments(...)` averages
the per-trajectory moment fields across the shared saved time grid. Ensemble and
plotting code should call `compute_ensemble_j_moments(...)`, which performs
both steps and returns trajectory-averaged moments.

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

Use this rule for $J_x$, $J_y$, and $N_(e)$ with the sector-local operators from
the precomputed `ops_list` described in `docs/instructions/simulation_precompute.typ`.

Implementation caveat: $J_z$ should follow the same expectation-value rule, but
the current code evaluates it directly from the reduced-basis diagonal
$n_(e) - N_(J) / 2$ rather than constructing a separate sparse operator.

The full wavefunction average active-sector atom number field `N_j` is the sector-weighted active
manifold size:

$
N_(J)(t_(k)) =
frac(
  sum_(alpha) N_(J)(alpha) chevron.l psi_(alpha)(t_(k)) | psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
).
$

For homogeneous sectors, $N_(J)(alpha)=N_(J)$. For group-resolved sectors,
$N_(J)(alpha)=sum_g N_(J,g)$.
For strong-symmetry-preserving dynamics, this weighted sector number is expected
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

== Drive and Group Fields

The drive field `Jx_drive` should use the drive-coupled operator defined by
`build_sector_ops_for_key(...)`; its construction is covered in
`docs/instructions/simulation_precompute.typ`. In this file it is just another
normalized sector expectation value following the rule above.

For group-resolved inhomogeneous results, group-resolved fields should be
returned for each group $g$:

$
chevron.l J_(i,g)(t_(k)) chevron.r =
frac(
  sum_(alpha) chevron.l psi_(alpha)(t_(k)) |
    J_(i,g,alpha) |
    psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
),
quad i in {x,y,z},
$

and

$
chevron.l N_(e,g)(t_(k)) chevron.r =
frac(
  sum_(alpha) chevron.l psi_(alpha)(t_(k)) |
    N_(e,g,alpha) |
    psi_(alpha)(t_(k)) chevron.r,
  Z(t_(k))
).
$

For homogeneous results, group fields should be `None`. The current
implementation supports homogeneous scalar sector keys and group-resolved tuple
sector keys; mixed or malformed sector structures should not be introduced
silently.

= Output

The snapshot helper should return a `JMomentSnapshot` with scalar fields:

```python
JMomentSnapshot(
    t, phase_index,
    Jx, Jy, Jz,
    N_e,
    jump_rate,
    N_j,
    Jx_drive,
    Jx_groups=None or tuple[float, ...],
    Jy_groups=None or tuple[float, ...],
    Jz_groups=None or tuple[float, ...],
    N_e_groups=None or tuple[float, ...],
)
```

The trajectory helper should stack one `JMomentSnapshot` per saved snapshot and
return `JMomentSeries`:

```python
JMomentSeries(
    t, phase_index,
    Jx, Jy, Jz,
    N_e,
    jump_rate,
    N_j,
    Jx_drive,
    Jx_groups=None or tuple[array, ...],
    Jy_groups=None or tuple[array, ...],
    Jz_groups=None or tuple[array, ...],
    N_e_groups=None or tuple[array, ...],
)
```

All trajectory arrays should be defined on the trajectory's saved `t_eval`
grid. `JMomentSeries` returned by `compute_ensemble_j_moments(...)` has the
same fields, but each numeric field is averaged across trajectories. These
helpers should not compute active-manifold angles, normalized spin directions,
squeezing, or covariance matrices.


= Invariants

- `compute_trajectory_j_moments(...)` should return per-trajectory moments only.
- `compute_average_j_moments(...)` should average moments, not nonlinear derived
  quantities.
- `compute_ensemble_j_moments(...)` should return trajectory-averaged moments and
  should be the preferred input for plots that do not need per-trajectory
  samples.
- `compute_ensemble_j_moments(...)` should require all internally computed
  samples to share the same `t` and `phase_index` grids.
