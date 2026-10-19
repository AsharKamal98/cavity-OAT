#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Generic J-Moment Post-Analysis]
]

= Purpose

This file specifies `compute_j_components(...)` and `compute_j_moments(...)`
in `post_analysis/j_moments.py`. Use it when converting already-averaged,
group-resolved SU(3) generator moments into `JMomentSeries` with optional SEMs.

The generator normalization and storage conventions originate in
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`. Shared
spin-vector direction and angle conventions are defined in
`docs/instructions/common/spin_moment_conventions.typ`, and the output classes
are documented in `docs/instructions/parser.typ`.

= J-Moment Definitions

== Averaged Generator Input

For group $g$ and one-based Gell--Mann label $a$, let
$overline(mu)_(g,a)(t)$ denote the normalized first generator moment after averaging
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
J_(x,g) = overline(mu)_(g,6),
quad
J_(y,g) = -overline(mu)_(g,7),
$

$
J_(z,g)
=
frac(1, 2) overline(mu)_(g,3)
-
frac(sqrt(3), 2) overline(mu)_(g,8).
$

The minus sign in $J_y=-G_7$ is required by the repository's convention
$J_y=(J_+-J_-)/(2i)$.

Define the shared, read-only `J_FROM_GENERATORS` matrix in
`post_analysis/j_moments.py`, and export it from `post_analysis` as well.
Its shape is `(3, 8)`, with rows ordered `(J_x, J_y, J_z)` and columns
ordered `(G_1, ..., G_8)`. It applies identically to every time and group.
Define separate named coefficient vectors for excited and active populations.
Combine them with `J_FROM_GENERATORS` to form the read-only `(5, 8)` covariance
map, ordered `(x, y, z, N_e, N_j)`. Compute the central spin components and
populations as explicitly named arrays using these same coefficients.
Population offsets depend only on the fixed physical group sizes.

Keep extraction code readable: name the physical population being computed
instead of modifying an unnamed packed column. Use distinct names for group
and aggregate covariance, standard deviations, and SEMs, and attach output
fields explicitly. Packed numeric axes belong to covariance propagation.

For physical group size $N_g$, calculate the atom-number fields

$
N_(e,g)
=
frac(N_g, 3)
-
frac(2, sqrt(3)) overline(mu)_(g,8),
$

$
N_(J,g)
=
frac(2 N_g, 3)
- overline(mu)_(g,3)
- frac(1, sqrt(3)) overline(mu)_(g,8).
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
group fields and their available sampling covariance have been attached. It
owns aggregate values, covariance, and SEMs. It sums the
additive fields first and derives nonlinear aggregated direction and angles
from the summed J vector.

== Directions, Angles, and Phase Indices

Apply `docs/instructions/common/spin_moment_conventions.typ` with the J-sphere
pole choice $p=-1$. First call
`JMomentSeries.attach_norm_spin_components_from_spin_components(...)`, then
`JMomentSeries.attach_angles_from_norm_spin_components(...)`.

Compute `integration_phase_index` from the output time grid and
`metadata.phase_protocol` through
`integration_phase_indices_at_times(...)`. Do not take the phase-index series
from the generator-moment input, which does not store phase indices.

== SEM Construction

`compute_j_components(...)` constructs `x_sem`, `y_sem`, `z_sem`, `N_e_sem`,
and `N_j_sem` alongside the corresponding central values. It follows the
generic covariance propagation, standard-deviation extraction, and SEM
normalization in `docs/instructions/post_analysis/uncertainty.typ`, using
`propagate_covariance(...)` and `std_from_covariance(...)`.

Use the full generator covariance from the mean input's optional `uncertainty`.
Store the propagated `(time, 5G, 5G)` sample covariance and original sample
count in `JMomentSeries.uncertainty`, with group-major field order
`(x, y, z, N_e, N_j)`. Convert its diagonal standard deviations to SEMs once.
Store SEMs as group tuples; do not store additional standard-deviation arrays.

The aggregation class method constructs aggregate SEMs from the sum of
matching fields across groups, including every group-pair covariance block.
It calls the shared uncertainty helpers in `common/utils/uncertainty.py`,
saves aggregate covariance as `(time, 5, 5)` with the same field order, and
saves one array per aggregate SEM. Summing group SEMs or variances would omit
correlations. Repeated calls recompute aggregate statistics from the current
group covariance; absent group uncertainty clears aggregate statistics.

If generator uncertainty is absent, leave SEMs and J uncertainty unset.

== Length, Direction, and Angle SEMs

`attach_norm_spin_components_from_spin_components(...)` owns `length_sem`,
`nx_sem`, `ny_sem`, and `nz_sem`. It calls
`spin_length_and_direction_sems(...)` from `common/utils/moments.py`, using
each group's Cartesian $3 times 3$ sample-covariance block. For averaged
Cartesian vector $bold(J)$, let $L=norm(bold(J))$ and $bold(n)=bold(J)/L$.
The Jacobian rows for length and normalized direction are

$
D_L = bold(n)^T,
quad D_n = (I - bold(n) bold(n)^T) / L.
$

`attach_angles_from_norm_spin_components(...)` owns `theta_sem` and
`phi_sem`. Its shared `spin_angle_sems(...)` helper composes angle and
normalization derivatives directly with respect to Cartesian components,
avoiding another stored covariance. With $rho=sqrt(n_x^2+n_y^2)$, the rows are

$
D_theta = vec(-n_x n_z, -n_y n_z, rho^2)^T / (L rho),
quad D_phi = vec(-n_y, n_x, 0)^T / (L rho^2).
$

The generic propagation and sample-count normalization remain those in
`docs/instructions/post_analysis/uncertainty.typ`. These nonlinear SEMs are
first-order approximations for functions of the averaged components. Keep
central values unchanged; do not average per-trajectory nonlinear quantities.
Aggregation calls the same helpers with its own summed vector and covariance,
including cross-group correlations before nonlinear propagation.

At $L <= "tol"$, all six derived SEMs are `NaN`. On the polar axis
($rho <= "tol"$), angle SEMs are `NaN`; length/direction SEMs remain available
when $L > "tol"$. Central fallback values remain unchanged. Angle SEMs are
local widths in radians; a plotting consumer crossing the azimuthal branch
cut must handle wrapping. Missing input uncertainty leaves all derived SEMs
`None`, and repeated class-method calls recompute or clear their owned SEMs.

= Method in Pseudo-Code

```python
def _group_series(values: array) -> tuple[array, ...]:
    # values has shape (len(t), G)
    return one copied time series per physical group

def compute_j_components(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *, Ni: tuple[int, ...], tol: float = 1e-12,
) -> JMomentSeries:
    require tol > 0 and a singleton input trajectory axis
    require one finite nonnegative Ni entry per physical group
    calculate named Cartesian x, y, z using J_FROM_GENERATORS
    calculate excited_population and active_population with named coefficients
    convert each (time, group) field with _group_series(...)
    result = JMomentSeries(x=x, y=y, z=z, N_e=excited_population,
                          N_j=active_population, ...)
    if source sampling uncertainty is available:
        construct the group-block-diagonal (5G, 8G) transformation
        group_covariance = propagate_covariance(source covariance, transformation)
        group_std = std_from_covariance(group_covariance)
        group_sem = group_std / sqrt(source sample count)
        attach each named group SEM and SamplingUncertainty explicitly

    JMomentSeries.attach_aggregated_from_group_fields(result, tol=tol)
        # This class method owns aggregate values, covariance, and SEMs.
    return result

def compute_j_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *,
    metadata: SimulationMetadata,
    tol: float = 1e-12,
) -> JMomentSeries:
    result = compute_j_components(
        mean_generator_moments, Ni=metadata.Ni, tol=tol,
    )
    result.integration_phase_index = integration_phase_indices_at_times(
        mean_generator_moments.t, metadata.phase_protocol,
    )

    JMomentSeries.attach_norm_spin_components_from_spin_components(
        result, tol=tol
    )
        # Attach length/direction values and their available SEMs.
    JMomentSeries.attach_angles_from_norm_spin_components(
        result, tol=tol
    )
        # Attach angle values and their available SEMs.
    return result
```

`_group_series(...)` converts each dense `(time, group)` result into the tuple
representation used by the parser class. `compute_j_components(...)` owns the
group affine values and optional SEMs, then delegates full-system values and
uncertainties to the aggregation class method, as defined in
`docs/instructions/parser.typ`. `compute_j_moments(...)` calls it
and attaches phase indices and group directions and angles through the shared
class methods. Both functions are exported from `post_analysis`.

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
    x_sem=None or tuple[array, ...],
    y_sem=None or tuple[array, ...],
    z_sem=None or tuple[array, ...],
    N_e_sem=None or tuple[array, ...],
    N_j_sem=None or tuple[array, ...],
    length_sem=None or tuple[array, ...],
    nx_sem=None or tuple[array, ...],
    ny_sem=None or tuple[array, ...],
    nz_sem=None or tuple[array, ...],
    theta_sem=None or tuple[array, ...],
    phi_sem=None or tuple[array, ...],
    uncertainty=None or SamplingUncertainty(...),
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
        x_sem=None or array,
        y_sem=None or array,
        z_sem=None or array,
        N_e_sem=None or array,
        N_j_sem=None or array,
        length_sem=None or array,
        nx_sem=None or array,
        ny_sem=None or array,
        nz_sem=None or array,
        theta_sem=None or array,
        phi_sem=None or array,
        uncertainty=None or SamplingUncertainty(...),
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
- Optional first-moment sampling covariance and sample count attached to the
  mean input; these are required only for SEMs.

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
- Preserve cross-group correlations and the original sample count in SEMs.
- A raw singleton or deterministic output has no SEM estimate; use `None`.
- Keep group and aggregate SEM shapes aligned with their central fields.
- Derived SEMs use first-order Jacobians of the averaged spin; preserve `NaN`
  at zero-length and angular singularities instead of reporting zero.
- Do not add jump rate to the J-moment output.
- Do not transform symmetric second moments or construct quantum covariance
  or squeezing during J extraction.

