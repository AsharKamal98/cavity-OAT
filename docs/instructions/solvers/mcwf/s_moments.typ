#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#let ket(x) = math.equation(block: false, $|#x chevron.r$)
#let bra(x) = math.equation(block: false, $chevron.l #x|$)

#align(center)[
  #text(size: 1.6em, weight: "bold")[
    MCWF S-Moment Extraction
  ]
]

= Purpose

This file specifies how `solvers/mcwf/s_moments.py` should convert saved MCWF
sector-block wavefunctions and an ensemble-averaged `JMomentSeries` into a
group-resolved `SMomentSeries`. Store the result as `moments.S`; it supplies
the S-sphere directions required by later squeezing-moment extraction.

Define `SMomentSeries` in `parser/s_moments.py` with spin-vector fields parallel
to `JMomentSeries`. Shared spin-vector direction and angle conventions are
defined in
`docs/instructions/common/spin_moment_conventions.typ`.

= S-Moment Definitions

== Effective S Qubit

At each saved time, group $g$ uses its ensemble-averaged J-sphere angles
$theta_(J,g), phi_(J,g)$ to define

$
|1_g chevron.r
=
cos(theta_(J,g)/2)|d chevron.r
+
e^(-i phi_(J,g)) sin(theta_(J,g)/2)|e chevron.r.
$

The effective S qubit is the group-local subspace spanned by
$|u chevron.r$ and $|1_g chevron.r$. Population in the orthogonal J-fluctuation
state lies outside this S sphere.

== Single-Particle and Collective Operators

For atom $i$ in group $g$, retain the S-sphere convention

$
s_(x,g)^((i))
=
frac(ket(1_g) bra(u) + ket(u) bra(1_g), 2),
quad
s_(y,g)^((i))
=
frac(ket(1_g) bra(u) - ket(u) bra(1_g), 2i),
$

$
s_(z,g)^((i))
=
frac(ket(u) bra(u) - ket(1_g) bra(1_g), 2).
$

Lift these to collective group operators:

$
S_(a,g) = sum_(i in g) s_(a,g)^((i)),
quad a in {x,y,z}.
$

== Trajectory and Ensemble Moments

MCWF snapshots may be unnormalized. For trajectory $r$ at saved time $t_k$,
use the complete sector-block norm

$
Z_r(t_k)
=
sum_(alpha)
chevron.l psi_(r,alpha)(t_k) | psi_(r,alpha)(t_k) chevron.r
$

and calculate

$
S_(a,g)^((r))(t_k)
=
frac(
  chevron.l psi_r(t_k) | S_(a,g) | psi_r(t_k) chevron.r,
  Z_r(t_k)
).
$

Average these raw expectation values equally over the $R$ trajectories:

$
bar(S)_(a,g)(t_k)
=
frac(1, R) sum_(r=1)^R S_(a,g)^((r))(t_k).
$

Do not average per-trajectory directions or angles.

== Normalized Direction and Angles

Apply `docs/instructions/common/spin_moment_conventions.typ` to the
trajectory-averaged S-sphere components, using the S-sphere pole choice
$p=+1$.

For each group $g$, compute `length[g]`, `nx[g]`, `ny[g]`, `nz[g]`,
`theta[g]`, and `phi[g]`, corresponding to
$S_("len",g)$, $n_(x,g)$, $n_(y,g)$, $n_(z,g)$,
$theta_(S,g)$, and $phi_(S,g)$.

Their mathematical definitions, averaging order, pole convention, and
zero-length behavior are defined in
`docs/instructions/common/spin_moment_conventions.typ`.

Use parallel `SMomentSeries` class methods to attach normalized components and
angles from the averaged raw fields. These methods should reuse the helpers in
`common/utils/moments.py` and encapsulate the S-sphere $p=+1$ sign mapping.
Do not construct an aggregated S vector because different groups may use
different local $|1_g chevron.r$ frames.

= Function Flow

Keep the extraction lean while returning an output parallel to
`JMomentSeries`:

```python
compute_trajectory_s_moments(
    trajectory: TrajectoryResult,
    metadata: TrajectoryEnsembleMetadata,
    j_moments: JMomentSeries,
    *,
    tol=1e-12,
) -> SMomentSeries
    # Return one raw group-resolved S-moment series for one trajectory.

compute_average_s_moments(
    samples: list[SMomentSeries],
) -> SMomentSeries
    # Average only raw x, y, and z fields across trajectories.

compute_mcwf_s_moments(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    *,
    tol=1e-12,
    n_processes=None,
) -> SMomentSeries

    trajectory_moments = map_with_optional_pool(
        compute_trajectory_s_moments(...)
        for trajectory in ensemble.trajectories
    )
    averaged = compute_average_s_moments(trajectory_moments)
    SMomentSeries.attach_norm_spin_components_from_spin_components(
        averaged,
        tol=tol,
    )
    SMomentSeries.attach_angles_from_norm_spin_components(
        averaged,
        tol=tol,
    )
    return averaged
```

`compute_trajectory_s_moments(...)` constructs $|1_g chevron.r$ from the
supplied ensemble-averaged J angles and returns one raw `SMomentSeries` for
each trajectory. `compute_average_s_moments(...)` averages those raw
per-trajectory series. `compute_mcwf_s_moments(...)` coordinates optional
multiprocessing and then attaches nonlinear direction and angle fields only to
the single averaged `SMomentSeries`.

Reuse the reduced collective one-body operator application logic currently
used by `solvers/mcwf/squeezing_moments.py`; do not duplicate its transition
rules or construct full tensor-product operators.

= Output

`SMomentSeries` should mirror the shared spin-vector fields of
`JMomentSeries`. Physical fields are tuples ordered by group, including for a
single-group simulation:

```python
SMomentSeries(
    t,
    integration_phase_index=None or array,
    x=tuple[array, ...],
    y=tuple[array, ...],
    z=tuple[array, ...],
    length=None or tuple[array, ...],
    nx=None or tuple[array, ...],
    ny=None or tuple[array, ...],
    nz=None or tuple[array, ...],
    theta=None or tuple[array, ...],
    phi=None or tuple[array, ...],
)
```

All series arrays have shape `(len(t),)`. Store the ensemble-averaged object as
`moments.S`.

= Data Requirements

- A non-empty `TrajectoryEnsemble` with metadata and saved `sector_blocks`.
- An ensemble-averaged `JMomentSeries` on the same `t_eval` grid, with one
  `theta` and `phi` series per physical group.
- Scalar single-group sectors or the current two-group tuple sectors.

= Invariants and Edge Cases

- Normalize every trajectory expectation by that snapshot's complete block
  norm before ensemble averaging.
- Follow `docs/instructions/common/spin_moment_conventions.typ`: average raw S
  moments before calculating nonlinear fields and normalize directions by the
  S-vector length.
- Preserve $S_z=(N_u-N_1)/2$ and the corresponding S-sphere pole choice
  $p=+1$.
- Return group tuples even when there is only one physical group.
- Do not add an `N_u1` field at this stage.
- Do not compute full-system, cross-group, covariance, or squeezing quantities.
- Moment extraction must not affect MCWF propagation.
- Do not construct full $3^N$ tensor-product operators.
