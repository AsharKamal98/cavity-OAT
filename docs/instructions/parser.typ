#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Parser Containers]
]

= Purpose

This file specifies parser/output-container conventions used by the current
post-processing pipeline. Use it when editing Pydantic classes in `parser/`,
especially `parser/j_moments.py`, `parser/s_moments.py`, and
`parser/moments.py`.

General parser-class style rules live in the `parser` skill. This file records
the repository-specific J- and S-moment container contracts.

= J-Moment Containers

`JMomentSnapshot` stores first-order J-sphere moments for one saved MCWF
trajectory snapshot:

```python
AggregatedJMomentSnapshot(
    x, y, z,
    N_e, N_j,
)

JMomentSnapshot(
    t,
    integration_phase_index,
    x, y, z,                 # tuples ordered by group
    N_e, N_j,                # tuples ordered by group
    aggregated=AggregatedJMomentSnapshot(...),
    jump_rate,
)
```

`JMomentSeries` stores one time series, either for one trajectory, an ensemble
average, or a non-MCWF method converted into the shared J-sphere format:

```python
AggregatedJMomentSeries(
    x=None, y=None, z=None,
    length=None,
    nx=None, ny=None, nz=None,
    N_e=None, N_j=None,
    theta=None, phi=None,
)

JMomentSeries(
    t,
    integration_phase_index=None,
    x=None, y=None, z=None,                   # tuples ordered by group
    length=None,
    nx=None, ny=None, nz=None,
    N_e=None, N_j=None,
    theta=None, phi=None,
    aggregated=AggregatedJMomentSeries(...),
    jump_rate=None,
)
```

The unsuffixed physical fields are the primary representation and should be
tuples ordered by group, including for $G=1$. The nested `aggregated` container
is the secondary full-system representation. Additive fields are summed over
groups; nonlinear fields are derived from the aggregated spin vector. Keep
`t`, `integration_phase_index`, and `jump_rate` on the outer container.

`SMomentSnapshot` and `SMomentSeries` use the parallel group-resolved
S-sphere structure. `SBasisSeries` stores the effective zero and one states:

```python
SMomentSnapshot(
    t,
    integration_phase_index,
    x, y, z,                 # tuples ordered by group
)

SMomentSeries(
    t,
    integration_phase_index=None,
    x=None, y=None, z=None,
    length=None,
    nx=None, ny=None, nz=None,
    theta=None, phi=None,
    basis=None or SBasisSeries(
        zero_state, one_state,  # tuples ordered by group
    ),
)
```

Each `zero_state[g]` and `one_state[g]` array has shape `(len(t), 3)` in the
single-particle $(|u chevron.r, |d chevron.r, |e chevron.r)$ basis. The state
series must be normalized, mutually orthogonal, and match the `SMomentSeries`
time grid and group count.

Store the ensemble-averaged series as
`moments.S: SMomentSeries | None`. Do not add J-specific atom-number,
jump-rate, or aggregated fields to the S containers.

`BlochVectorSeries` is the reusable Cartesian-vector container defined in
`parser/bloch_vector.py`:

```python
BlochVectorSeries(
    x, y, z,
    length,     # derived automatically
    xy_length,  # derived automatically
)
```

Every construction should derive

$
"length" = sqrt(x^2+y^2+z^2),
quad
"xy_length" = sqrt(x^2+y^2).
$

Callers should supply matching one-dimensional x/y/z arrays and should not
compute either derived field themselves.

`JModeSeries` stores post-analysis combinations separately from raw J moments:

```python
JModeSeries(
    t,
    common=BlochVectorSeries(x, y, z, length, xy_length),
    contrast=BlochVectorSeries(x, y, z, length, xy_length),
    bright=BlochVectorSeries(x, y, z, length, xy_length),
    dark=BlochVectorSeries(x, y, z, length, xy_length),
)
```

Store it as nullable `moments.J_modes`; do not add derived mode fields to
`JMomentSeries`.

`HarmonicAnalysisSeries` stores solver-independent harmonic analyses for one or
more curves. Each `HarmonicSignalAnalysis` contains its fitted offset and
fundamental frequency, observed cycle count, harmonic indices, frequencies,
amplitudes, RMS oscillation amplitude, phase offsets, and total harmonic
distortion. The output does not retain the input time grid. Keep this diagnostic
separate from `JMomentSeries` and `JModeSeries` so it can consume any compatible
real series.

`SqueezingMomentSeries` stores the ensemble-averaged raw moments required to
compute generalized squeezing for one physical group:

```python
SqueezingMomentSeries(
    t,
    N_group,
    mean_operators,          # shape (n_t, 4)
    second_moments,          # shape (n_t, 4, 4)
    N_c,
    excited_fraction_active,
    theta_J, phi_J,
    theta_S, phi_S,
)
```

Store group-local squeezing moments as
`moments.squeezing_moments: tuple[SqueezingMomentSeries, ...] | None`, ordered
by physical group. These are moments of dressed fluctuation operators and
should not be added to `JMomentSeries`.

`SqueezingSeries` stores the derived generalized squeezing result for one
physical group:

```python
SqueezingSeries(
    t,
    N_group,
    xi2,
    xi2_db,
    lambda_min,
    covariance_eigvals,     # shape (n_t, 4)
    N_c,
    excited_fraction_active,
    theta_J, phi_J,
    theta_S, phi_S,
)
```

Store derived group results as
`moments.squeezing: tuple[SqueezingSeries, ...] | None`, with the same group
ordering as `moments.squeezing_moments`.

= Simulation Metadata

`SimulationMetadata` stores shared physical model inputs and the supplied phase
protocol for one run:

```python
SimulationMetadata(
    Ni, omega_i,
    Gamma,
    phase_protocol,
)
```

`omega_i` contains the first $G-1$ independent couplings and `Ni` contains
$G$ group sizes. Its validator constructs the final weighted-average coupling
and saves the full-length `omega_groups`. Construct `phase_protocol` separately,
as defined in `docs/instructions/phases.typ`, and provide it explicitly. Store
this context as `moments.metadata`.

When `MomentSeries` receives `metadata` and `num_snapshots`, its validator
constructs `t` from `metadata.phase_protocol.total_duration`; do not pass
`total_time` separately.

Use this object to construct solver-input containers, but do not store it as a
solver-parameter field. Each solver container should explicitly receive only
the physical fields it consumes, such as `Ni`, completed `omega_i`, `Gamma`,
and `phase_protocol`, alongside its method-specific controls.

= Derived-Field Class Methods

The `JMomentSeries` class methods should attach missing equivalent
representations in place:

```python
JMomentSeries.attach_norm_spin_components_from_spin_components(j_moments, tol=...)
    -> attach length, nx, ny, nz from x, y, z

JMomentSeries.attach_angles_from_norm_spin_components(j_moments, tol=...)
    -> attach theta, phi from length, nx, ny, nz

JMomentSeries.attach_norm_spin_components_from_angles(j_moments)
    -> attach nx, ny, nz from theta, phi

JMomentSeries.attach_spin_components_from_norm_spin_components(j_moments)
    -> attach x, y, z from length, nx, ny, nz

JMomentSeries.attach_aggregated_from_group_fields(j_moments, tol=...)
    -> attach aggregated.x, aggregated.y, aggregated.z,
       aggregated.N_e, aggregated.N_j by summing their group fields
    -> attach aggregated.length, aggregated.nx, aggregated.ny,
       aggregated.nz, aggregated.theta, aggregated.phi from the
       summed spin components

SMomentSeries.attach_norm_spin_components_from_spin_components(s_moments, tol=...)
    -> attach group length, nx, ny, nz from x, y, z

SMomentSeries.attach_angles_from_norm_spin_components(s_moments, tol=...)
    -> attach group theta, phi using the S-sphere positive-z pole convention
```

The spin-component and angle conversion methods operate on the primary
group-resolved fields. `attach_aggregated_from_group_fields(...)` is the sole
series-level aggregation method and constructs a new complete
`AggregatedJMomentSeries`. Each method should raise a clear `ValueError` when
its required group inputs are unavailable.

The conversion helpers used by these class methods live in
`common/utils/moments.py`.
Shared mathematical definitions and the J/S pole mapping live in
`docs/instructions/common/spin_moment_conventions.typ`.

= Conversion Logic

Given spin components, compute vector length and normalized direction:

$
L(t) = sqrt(x(t)^2 + y(t)^2 + z(t)^2),
quad
n_i(t) = frac(i(t), L(t)).
$

Given normalized direction, compute angles:

$
theta(t) = arccos(-n_z(t)),
quad
phi(t) = "atan2"(n_y(t), n_x(t)).
$

The minus sign in $theta$ preserves the active-manifold convention used by the
J-moment pipeline. If the vector length is below tolerance, normalized
components and angles should use the safe fallback behavior implemented by
`common/utils/moments.py`.

Angles alone determine only the normalized direction. To recover spin
components, the series must also contain a vector length.

For additive observables $q in {x, y, z, N_e, N_J}$, construct the aggregated
field from the primary group fields as

$
q(t) = sum_g q_g(t).
$

Do not apply this sum to vector lengths, normalized directions, or angles.

= Invariants

- Parser classes should stay as Pydantic containers with explicit typed fields.
- `BlochVectorSeries` must derive `length` and `xy_length` from x/y/z on every
  construction so stored values cannot disagree with the Cartesian components.
- `SqueezingMomentSeries` must validate `(n_t, 4)` operator means,
  `(n_t, 4, 4)` second moments, and one-dimensional fields matching `t`.
- `SqueezingSeries` must validate `(n_t, 4)` covariance eigenvalues and
  one-dimensional result fields matching `t`.
- Derived-field class methods should mutate the supplied `JMomentSeries`
  instance in place.
- `SMomentSeries` should attach derived fields only after raw per-trajectory
  S moments have been averaged.
- `SMomentSeries` should use the S-sphere $p=+1$ pole convention and should not
  construct aggregated fields.
- `SMomentSeries.basis` should hold the shared effective S-qubit zero and one
  states used to compute its S moments.
- Aggregated and group conversions should call the same pure conversion
  helpers.
- Only additive fields `x`, `y`, `z`, `N_e`, and `N_j` should be summed across
  groups; save the sums in `aggregated` and derive aggregated lengths,
  directions, and angles from `aggregated.x/y/z`.
- Construct series-level aggregated fields only through
  `attach_aggregated_from_group_fields(...)`, after raw group fields have been
  computed or averaged.
- Unsuffixed fields are always group tuples in public backend outputs, including
  for $G=1$.
- J-moment extraction code should call these class methods instead of
  duplicating spin-component, direction, or angle conversion logic.
- Legacy field names such as `Jx`, `Jy`, `Jz`, `J_len`, and `sx`, `sy`, `sz`
  should not be reintroduced.
