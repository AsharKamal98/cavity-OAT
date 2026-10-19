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
especially `parser/generator_moments.py`, `parser/mcwf.py`,
`parser/j_moments.py`, `parser/s_moments.py`, `parser/squeezing.py`, and
`parser/moments.py`, together with `parser/uncertainty.py`.

General parser-class style rules live in the `parser` skill. This file records
the repository-specific J- and S-moment container contracts.

= MCWF Initial-Sector Containers

`MCWFSolverParameters` should expose independent optional initial-sector filters:

```python
MCWFSolverParameters(
    ...,
    dN=2,
    max_omega_over_omega_c=None,
    linear_constraints=None,
    joint_probability=None,
)
```

Any combination, including all `None`, is valid. Apply non-`None` inputs in
the order `dN`, `max_omega_over_omega_c`, `linear_constraints`, then
`joint_probability`; the latter is always the final tuple-wise filter rather
than an alternative selection branch. Validate
`0 < max_omega_over_omega_c <= 1`. A supplied
`linear_constraints=(A, b)` should be converted into
`LinearSectorConstraints(matrix=A, offset=b)`. Its integer matrix must have one
column per physical group and one offset per row.

`InitialSectorFilterSummary` stores `filter_name`, `input_sector_count`,
`remaining_sector_count`, `input_joint_probability`, and
`remaining_joint_probability` for one applied filter. Both probabilities are
relative to the complete initial wavefunction. `InitialSectorSelectionSummary`
stores these ordered rows in `filter_steps` together with the retained,
pre-probability-filter, and full tuple-sector counts and the corresponding
full-wavefunction probabilities. It is a transient setup diagnostic printed
once by the ensemble runner; it is not copied into every trajectory result.

= Raw Trajectory-Resolved Generator-Moment Container

`TrajectoryResolvedMomentSeries` in `parser/generator_moments.py` stores the
complete SU(3) moment output for one or more individually retained trajectories
in one caller-selected generator frame:

```python
TrajectoryResolvedMomentSeries(
    t,                       # shape (n_t,)
    first_moments,           # shape (R, n_t, G, 8)
    second_moments,          # shape (R, n_t, K)
    uncertainty=None,        # optional statistics attached to an ensemble mean
)
```

Here $R$ is the explicit trajectory axis and must always be present. A single
MCWF trajectory and a deterministic solver output use $R=1$; an MCWF ensemble
uses $R=n_"traj"$. The first moments use one explicit physical-group axis
followed by the eight Gell--Mann generators. Let $D=8G$. The second moments
store the packed upper triangle of the real symmetric $D times D$
anticommutator-moment matrix, so $K=D(D+1)/2$. Its packed axis uses the
row-major ordering returned by `numpy.triu_indices(D)`.

The container does not encode frame provenance. MCWF simulation output and its
trajectory mean use the lab frame; `rotate_generator_moments(...)` may return
the same container type in a dressed frame. Keep those series distinct, and do
not pass a dressed series into analyses whose generator mapping assumes the lab
frame. The values are generator moments rather than per-trajectory covariance
matrices. A selected ensemble mean reuses the class with $R=1$.
The saved-time grid must be nonempty, finite, and strictly increasing. Both
moment arrays must be finite and real; construction may discard only an
imaginary residual that is zero within the parser tolerance.

`uncertainty` is an optional `SamplingUncertainty` from
`parser/uncertainty.py`. When present, require a singleton trajectory axis
and `sample_covariance` shaped `(n_t, 8G, 8G)` in group-major generator order.
Its `n_samples` records the selected source-trajectory count, independently of
the singleton stored axis. The generic container validates finite real
symmetric covariance with nonnegative variances and an integer count of at
least two. Keep this container limited to `sample_covariance` and `n_samples`;
compute mean covariance and SEM in downstream analysis when needed.

`compute_mean_generator_moments(...)` automatically attaches first-moment
sampling statistics, including cross-group correlations. Raw and
deterministic singleton data leave the field unset. These statistics do not
describe the quantum covariance built from `second_moments`. The complete
contract is in `docs/instructions/post_analysis/uncertainty.typ`.

= Analysis Moment Container

`MomentSeries` should keep the raw generator moments together with derived J,
S, squeezing, and diagnostic results:

```python
MomentSeries(
    t,                                    # shape (n_t,)
    metadata=None or SimulationMetadata(...),
    generator_moments=None or TrajectoryResolvedMomentSeries(...),
    mean_generator_moments=None or TrajectoryResolvedMomentSeries(...),
    J=None,
    S=None,
    squeezing=None,
    ...,
)
```

When `generator_moments` is present, `MomentSeries.t` and
`generator_moments.t` must have equal shape and agree to the saved-time grid
tolerance. Its physical-group count must equal `len(metadata.Ni)` when metadata
is present. These invariants must also be checked when `generator_moments`,
`mean_generator_moments`, `t`, or `metadata` is assigned after construction.

If no explicit `t` is supplied, `MomentSeries` should first use the time grid
from `generator_moments`, then from `mean_generator_moments`. Otherwise, the
existing `metadata` plus `num_snapshots` construction should build the grid.
Keep `first_moments` and `second_moments` nested inside their generator-moment
containers; do not copy them into loose top-level `MomentSeries` array fields.

`mean_generator_moments` stores the trajectory mean in the same container
type, with a required singleton trajectory axis $R=1$. It must match
`generator_moments` in time, group, generator, and packed-second-moment axes.
It may also be stored without the raw trajectories when loading an artifact
that contains only the central ensemble result.
Keep its sampling statistics nested in `mean_generator_moments.uncertainty`
so linear-observable SEMs remain recoverable even without the raw trajectories.


= J-Moment Containers

`JMomentSeries` stores one selected trajectory average or a deterministic
method output converted into the shared J-sphere format:

```python
AggregatedJMomentSeries(
    x=None, y=None, z=None,
    length=None,
    nx=None, ny=None, nz=None,
    N_e=None, N_j=None,
    theta=None, phi=None,
    x_sem=None, y_sem=None, z_sem=None,       # one array per aggregate field
    N_e_sem=None, N_j_sem=None,
    length_sem=None, nx_sem=None, ny_sem=None, nz_sem=None,
    theta_sem=None, phi_sem=None,
    uncertainty=None,
)

JMomentSeries(
    t,
    integration_phase_index=None,
    x=None, y=None, z=None,                   # tuples ordered by group
    length=None,
    nx=None, ny=None, nz=None,
    N_e=None, N_j=None,
    theta=None, phi=None,
    x_sem=None, y_sem=None, z_sem=None,       # tuples ordered by group
    N_e_sem=None, N_j_sem=None,
    length_sem=None, nx_sem=None, ny_sem=None, nz_sem=None,
    theta_sem=None, phi_sem=None,
    uncertainty=None,
    aggregated=AggregatedJMomentSeries(...),
)
```

The unsuffixed physical fields are the primary representation and should be
tuples ordered by group, including for $G=1$. The nested `aggregated` container
is the secondary full-system representation. Additive fields are summed over
groups; nonlinear fields are derived from the aggregated spin vector. Keep
`t` and `integration_phase_index` on the outer container. Jump rate is not a
J-moment field.

The eleven nullable `*_sem` fields store sampling standard errors of the
corresponding observables, never trajectory standard deviations. Arrays must
be one-dimensional, nonnegative, and match their central fields and saved-time
grid; group tuples must match group order and count. Affine-field SEMs must
be finite. Derived `length/nx/ny/nz/theta/phi` SEMs may contain `NaN` where
the local derivative is undefined, but must not contain infinities. Missing
source statistics leave every SEM field and `uncertainty` unset, including
for deterministic MFE output.

`JMomentSeries.uncertainty` uses the minimal `SamplingUncertainty` container.
Its sample covariance has shape `(len(t), 5G, 5G)`, with group-major order
`(x, y, z, N_e, N_j)` per group. The nested aggregate uncertainty stores
`(len(t), 5, 5)` covariance in the same field order. Both retain the original
sample count; neither stores an extra standard-deviation array. Preserve the
full covariance as the source for first-order nonlinear error propagation;
do not store additional derived covariance or standard-deviation arrays. See
`docs/instructions/post_analysis/j_moments.typ` for construction and
`docs/instructions/post_analysis/uncertainty.typ` for the generic method.

The generic post-analysis constructs `SMomentSeries`; `SBasisSeries` stores its
effective zero and one states:

```python
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
`moments.S: SMomentSeries | None`. Do not add J-specific atom-number or
aggregated fields to the S containers.

`BlochVectorSeries` is the reusable Cartesian-vector container defined in
`parser/bloch_vector.py`:

```python
BlochVectorSeries(
    x, y, z,
    length,     # derived automatically
    xy_length,  # derived automatically
    uncertainty=None,
    x_sem=None, y_sem=None, z_sem=None,
    length_sem=None, xy_length_sem=None,
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

Optional `uncertainty` is a `SamplingUncertainty` with `(len(x), 3, 3)`
Cartesian sample covariance and the original sample count. The container
automatically derives component SEMs and first-order length/xy-length SEMs.
It calls the generic helpers in `common/utils/uncertainty.py` and
`vector_length_sem(...)` from `common/utils/moments.py`; follow
`docs/instructions/post_analysis/uncertainty.typ` for sample normalization.
Zero-length derivatives give `NaN` norm SEMs. Without covariance, all SEMs
are reset to `None`. These fields are derived outputs, not independent inputs.
Retain the Cartesian covariance for subsequent observables.

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
`JMomentSeries`. Each mode carries its own Cartesian covariance and SEMs when
constructed with sampling statistics. Construction and population-dependent
Jacobians are defined in `docs/instructions/post_analysis/j_modes.typ`.

`HarmonicAnalysisSeries` stores solver-independent harmonic analyses for one or
more curves. Each `HarmonicSignalAnalysis` contains its fitted offset and
fundamental frequency, observed cycle count, harmonic indices, frequencies,
amplitudes, RMS oscillation amplitude, phase offsets, and total harmonic
distortion. The output does not retain the input time grid. Keep this diagnostic
separate from `JMomentSeries` and `JModeSeries` so it can consume any compatible
real series.

`WindowedHarmonicAnalysisSeries` in the same parser module stores
`window_centers` with shape `(windows,)` and `offsets`,
`fundamental_amplitudes`, `rms_amplitudes`, `fundamental_frequencies`, and
`total_harmonic_distortions` with
shape `(windows, curves)`, including a singleton column for one curve.
Unresolved curve/window fits use `NaN`; constant curves have zero amplitudes
and frequency with their fitted offset retained and `NaN` THD. Window construction and
fitting conventions live in
`docs/instructions/post_analysis/harmonic_analysis.typ`.

`SqueezingSeries` stores the complete generalized squeezing result for either
one physical group or the full system:

```python
SqueezingSeries(
    t,
    N_atoms,
    covariance,              # shape (n_t, 4, 4)
    covariance_eigvals,      # shape (n_t, 4), ascending
    lambda_min,
    xi2,
    xi2_db,
    N_c,
    excited_fraction_active,
)
```

`SqueezingResult` keeps the mandatory group-resolved series and optional
full-system series distinct:

```python
SqueezingResult(
    group_resolved: tuple[SqueezingSeries, ...],
    full_system: SqueezingSeries | None = None,
)
```

Store this object as `moments.squeezing: SqueezingResult | None`. The group
tuple follows `metadata.Ni` order and its atom counts must match those group
sizes. The full-system atom count, when present, is their sum. All series must
share `MomentSeries.t`. Dressed moments may remain a transient analysis value
between explicit rotation and squeezing and do not need a separate parser
container. J and S angles remain in
their source series rather than being duplicated here.

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
    -> attach length_sem, nx_sem, ny_sem, nz_sem when covariance is available

JMomentSeries.attach_angles_from_norm_spin_components(j_moments, tol=...)
    -> attach theta, phi from length, nx, ny, nz
    -> attach theta_sem, phi_sem when covariance is available

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
    -> attach aggregate covariance and x/y/z/N_e/N_j SEMs when
       j_moments.uncertainty is available
    -> attach first-order aggregate length, direction, and angle SEMs

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

Aggregation owns aggregate values and their sampling uncertainties. Attach
group covariance before calling it. Propagate the full group covariance,
including cross-group blocks, using the shared helpers in
`common/utils/uncertainty.py`; follow
`docs/instructions/post_analysis/uncertainty.typ` for the generic SEM flow.
Recompute these statistics on every aggregation call. Without group
uncertainty, the new aggregate leaves covariance and SEMs unset, clearing any
stale aggregate statistics.

Length/direction and angle attachment methods own their corresponding SEMs.
They delegate to `spin_length_and_direction_sems(...)` and
`spin_angle_sems(...)` in `common/utils/moments.py`, using each group's
Cartesian covariance block. Aggregation uses the same helpers after summing
the full covariance. The angle helper composes derivatives directly with
respect to Cartesian components. Without input statistics, reset the owned
SEMs to `None`. At length below tolerance, length and direction SEMs are
`NaN`; angle SEMs are also `NaN` on the polar axis. Keep existing central-value
fallbacks. Angle SEMs are local widths in radians, so plotting across the
azimuthal branch cut requires wrap-aware handling.

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
- `TrajectoryResolvedMomentSeries` should require a nonempty finite increasing
  time grid and finite real first- and second-moment arrays.
- `BlochVectorSeries` must derive `length` and `xy_length` from x/y/z on every
  construction so stored values cannot disagree with the Cartesian components.
- `SqueezingSeries` must validate an `(n_t, 4, 4)` symmetric covariance,
  `(n_t, 4)` covariance eigenvalues, and one-dimensional result fields matching
  `t`.
- `SqueezingResult` must contain at least one group result on a shared grid;
  any full-system result must use that grid and the sum of the group atom
  counts.
- Derived-field class methods should mutate the supplied `JMomentSeries`
  instance in place.
- `SMomentSeries` should attach derived fields only after the selected mean
  first generator moments have been converted into raw S components.
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
  computed or averaged and their available covariance has been attached.
- Unsuffixed fields are always group tuples in public backend outputs, including
  for $G=1$.
- J-moment extraction code should call these class methods instead of
  duplicating spin-component, direction, or angle conversion logic.
