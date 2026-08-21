#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Generic J-Moment Post-Analysis]
]

= Purpose

This file specifies `compute_j_moments(...)` in
`post_analysis/j_moments.py`. Use it when converting already-averaged,
group-resolved SU(3) generator moments into the shared `JMomentSeries` output.

The generator normalization and storage conventions originate in
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`. Shared
spin-vector direction and angle conventions are defined in
`docs/instructions/common/spin_moment_conventions.typ`, and the output classes
are documented in `docs/instructions/parser.typ`.

= J-Moment Definitions

== Averaged Generator Input

For group $g$ and one-based Gell--Mann label $a$, let
$bar(mu)_(g,a)(t)$ denote the normalized first generator moment after averaging
the selected complete trajectories. The input
`TrajectoryResolvedMomentSeries` retains this mean on a singleton trajectory
axis, so its `first_moments` shape is `(1, len(t), G, 8)`.

The implementation uses $g_a=lambda_a/2$ in the
$(|u chevron.r, |d chevron.r, |e chevron.r)$ basis. Equations in this file use
one-based labels $a in {1, dots.c, 8}$; Python indexes the final array axis from
zero.

== Generator-to-J Mapping

For every group and saved time, calculate

$
J_(x,g) = bar(mu)_(g,6),
quad
J_(y,g) = -bar(mu)_(g,7),
$

$
J_(z,g)
=
frac(1, 2) bar(mu)_(g,3)
-
frac(sqrt(3), 2) bar(mu)_(g,8).
$

The minus sign in $J_y=-G_7$ is required by the repository's convention
$J_y=(J_+-J_-)/(2i)$.

For physical group size $N_g$, calculate the atom-number fields

$
N_(e,g)
=
frac(N_g, 3)
-
frac(2, sqrt(3)) bar(mu)_(g,8),
$

$
N_(J,g)
=
frac(2 N_g, 3)
- bar(mu)_(g,3)
- frac(1, sqrt(3)) bar(mu)_(g,8).
$

Only first generator moments enter these expressions. This function does not
transform the stored symmetric second moments.

== Group and Aggregated Fields

The primary physical fields are tuples ordered by group, including for $G=1$.
Construct the additive full-system fields from the group fields:

$
q_"aggregated"(t) = sum_g q_g(t),
quad
q in {J_x, J_y, J_z, N_e, N_J}.
$

Call `JMomentSeries.attach_aggregated_from_group_fields(...)` after the raw
group fields have been constructed. It sums the additive fields first and then
derives the nonlinear aggregated direction and angles from the summed J
vector.

== Directions, Angles, and Phase Indices

Apply `docs/instructions/common/spin_moment_conventions.typ` with the J-sphere
pole choice $p=-1$. First call
`JMomentSeries.attach_norm_spin_components_from_spin_components(...)`, then
`JMomentSeries.attach_angles_from_norm_spin_components(...)`.

Compute `integration_phase_index` from the output time grid and
`metadata.phase_protocol` through
`integration_phase_indices_at_times(...)`. Do not take the phase-index series
from the generator-moment input, which stores only its time grid and moments.

= Method in Pseudo-Code

```python
def _group_series(values: array) -> tuple[array, ...]:
    # values has shape (len(t), G)
    return one copied time series per physical group

def compute_j_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *,
    metadata: SimulationMetadata,
    tol: float = 1e-12,
) -> JMomentSeries:
    require tol > 0
    require mean_generator_moments.num_trajectories == 1

    first = mean_generator_moments.first_moments[0]
    require first.shape[1] == len(metadata.Ni)

    calculate group x, y, z, N_e, and N_j from first and metadata.Ni
    result = JMomentSeries(
        t=copy of mean_generator_moments.t,
        integration_phase_index=integration_phase_indices_at_times(
            mean_generator_moments.t,
            metadata.phase_protocol,
        ),
        x=_group_series(x),
        y=_group_series(y),
        z=_group_series(z),
        N_e=_group_series(N_e),
        N_j=_group_series(N_j),
    )

    JMomentSeries.attach_norm_spin_components_from_spin_components(
        result, tol=tol
    )
    JMomentSeries.attach_angles_from_norm_spin_components(
        result, tol=tol
    )
    JMomentSeries.attach_aggregated_from_group_fields(
        result, tol=tol
    )
    return result
```

`_group_series(...)` converts each dense `(time, group)` result into the tuple
representation used by the parser class. `compute_j_moments(...)` validates the
central generator input, applies the fixed generator mapping, constructs the
raw output, and then delegates all nonlinear and aggregated fields to the
`JMomentSeries` class methods.

= Output

Return one `JMomentSeries`:

```python
JMomentSeries(
    t=array,
    integration_phase_index=array,
    x=tuple[array, ...],
    y=tuple[array, ...],
    z=tuple[array, ...],
    length=tuple[array, ...],
    nx=tuple[array, ...],
    ny=tuple[array, ...],
    nz=tuple[array, ...],
    N_e=tuple[array, ...],
    N_j=tuple[array, ...],
    theta=tuple[array, ...],
    phi=tuple[array, ...],
    aggregated=AggregatedJMomentSeries(
        x=array,
        y=array,
        z=array,
        length=array,
        nx=array,
        ny=array,
        nz=array,
        N_e=array,
        N_j=array,
        theta=array,
        phi=array,
    ),
)
```

Every numeric series uses the input saved-time grid. Store this output as
`moments.J`.

= Data Requirements

- `mean_generator_moments`: a `TrajectoryResolvedMomentSeries` with a required
  singleton trajectory axis and first moments shaped `(1, len(t), G, 8)`.
- `metadata.Ni`: one physical group size per generator group.
- `metadata.phase_protocol`: the protocol used to derive phase indices.
- A strictly positive numerical tolerance.

Although `SimulationMetadata` also contains coupling and decay information,
the current function consumes only `Ni` and `phase_protocol`.

= Invariants and Edge Cases

- Average or resample complete trajectories before calling this function. It
  does not average internally.
- Preserve $g_a=lambda_a/2$ and the exact $J_y=-G_7$ sign convention.
- Keep the singleton trajectory axis on the mean generator input.
- Return group tuples for every group count, including $G=1$.
- Derive group directions and angles only after the raw group fields exist.
- Sum raw additive group fields before deriving aggregated nonlinear fields.
- Use the J-sphere $p=-1$ pole convention and the shared zero-length fallback.
- Do not add jump rate to the J-moment output.
- Do not transform symmetric second moments or compute covariance or squeezing.

