#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#let ket(x) = math.equation(block: false, $|#x chevron.r$)
#let bra(x) = math.equation(block: false, $chevron.l #x|$)

#align(center)[
  #text(size: 1.6em, weight: "bold")[Generic S-Moment Post-Analysis]
]

= Purpose

This file specifies `compute_s_moments(...)` in
`post_analysis/s_moments.py`. Use it when converting already-averaged SU(3)
generator moments into the group-local `SMomentSeries` defined by the
corresponding group-resolved J angles.

The effective S-qubit, generator mapping, and pole conventions are specified
below. Shared spin-vector conversions are defined in
`docs/instructions/common/spin_moment_conventions.typ`, and the output classes
are documented in `docs/instructions/parser.typ`.

= S-Moment Definitions

== Effective S Qubit

At each saved time, group $g$ uses its J-sphere angles
$theta_(J,g), phi_(J,g)$ to define

$
|0_g chevron.r = |u chevron.r,
quad
|1_g chevron.r
=
cos(theta_(J,g)/2)|d chevron.r
+
e^(-i phi_(J,g)) sin(theta_(J,g)/2)|e chevron.r.
$

Store this shared basis as `SMomentSeries.basis`. Each group entry in
`zero_state` and `one_state` has shape `(len(t), 3)` in the
$(|u chevron.r, |d chevron.r, |e chevron.r)$ basis. `SBasisSeries` validates
the shapes, normalization, and orthogonality of these states.

The corresponding single-particle S operators retain the definitions

$
s_(x,g)^((i))
=
frac(ket(1_g) bra(0_g) + ket(0_g) bra(1_g), 2),
$

$
s_(y,g)^((i))
=
frac(ket(1_g) bra(0_g) - ket(0_g) bra(1_g), 2i),
quad
s_(z,g)^((i))
=
frac(ket(0_g) bra(0_g) - ket(1_g) bra(1_g), 2),
$

with collective group operators $S_(a,g)=sum_(i in g)s_(a,g)^((i))$.

== Generator-to-S Mapping

Let

$
c_g = cos(theta_(J,g)/2),
quad
s_g = sin(theta_(J,g)/2),
$

and let $bar(mu)_(g,a)$ be the already-averaged first generator moments using
$g_a=lambda_a/2$. The generic implementation evaluates the S components as
the following linear combinations:

$
S_(x,g)
=
c_g bar(mu)_(g,1)
+ s_g cos(phi_(J,g)) bar(mu)_(g,4)
- s_g sin(phi_(J,g)) bar(mu)_(g,5),
$

$
S_(y,g)
=
-c_g bar(mu)_(g,2)
- s_g sin(phi_(J,g)) bar(mu)_(g,4)
- s_g cos(phi_(J,g)) bar(mu)_(g,5),
$

$
S_(z,g)
=
frac(1+c_g^2, 2) bar(mu)_(g,3)
- frac(sin(theta_(J,g)) cos(phi_(J,g)), 2) bar(mu)_(g,6)
+ frac(sin(theta_(J,g)) sin(phi_(J,g)), 2) bar(mu)_(g,7)
+ frac(sqrt(3) s_g^2, 2) bar(mu)_(g,8).
$

Equivalently, at each time and group,

$
bold(S)_g
=
C_(S,g)(theta_(J,g), phi_(J,g)) bar(bold(mu))_g,
$

where the implementation stores $C_S$ with shape `(len(t), G, 3, 8)` and
contracts it with first moments of shape `(len(t), G, 8)`. Equations use
one-based generator labels; Python uses zero-based indexes on the final axis.

The J angles define one shared frame for the selected trajectory average, so
the linear transformation acts directly on the averaged first moments. The
function does not calculate per-trajectory S vectors and does not transform
symmetric second moments.

== Normalized Direction and Angles

Apply `docs/instructions/common/spin_moment_conventions.typ` with the S-sphere
pole choice $p=+1$. After constructing the raw group components, call
`SMomentSeries.attach_norm_spin_components_from_spin_components(...)` and then
`SMomentSeries.attach_angles_from_norm_spin_components(...)`.

The current angle class method implements $p=+1$ by passing `-nz` into the
shared helper that natively uses the J-sphere $p=-1$ convention. Do not
construct an aggregated S vector because different groups generally use
different local $|1_g chevron.r$ frames.

= Method in Pseudo-Code

```python
def _validate_reference_data(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    j_moments: JMomentSeries,
) -> tuple[array, array, array, array | None]:
    require mean_generator_moments.num_trajectories == 1
    require mean_generator_moments.t and j_moments.t to match within 1e-12
    require group-resolved j_moments.theta and j_moments.phi
    require both J-angle arrays to have shape (G, len(t))
    if integration_phase_index exists, require shape (len(t),)
    return mean first moments, J angles, and optional phase indices

def _build_s_basis(theta_J_groups, phi_J_groups) -> SBasisSeries:
    build |0_g> = |u> and the J-angle-dependent |1_g> for every group and time
    return validated SBasisSeries

def _build_s_generator_coefficients(
    theta_J_groups,
    phi_J_groups,
) -> array:
    return C_S with shape (len(t), G, 3, 8)

def _group_series(values: array) -> tuple[array, ...]:
    # values has shape (len(t), G)
    return one copied time series per physical group

def compute_s_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    j_moments: JMomentSeries,
    *,
    tol: float = 1e-12,
) -> SMomentSeries:
    require tol > 0
    first, theta_J, phi_J, phase_index = _validate_reference_data(
        mean_generator_moments,
        j_moments,
    )
    coefficients = _build_s_generator_coefficients(theta_J, phi_J)
    components = einsum(
        "tgsa,tga->tgs", coefficients, first, optimize=True
    )

    result = SMomentSeries(
        t=copy of mean_generator_moments.t,
        integration_phase_index=copy of phase_index or None,
        x=_group_series(components[:, :, 0]),
        y=_group_series(components[:, :, 1]),
        z=_group_series(components[:, :, 2]),
        basis=_build_s_basis(theta_J, phi_J),
    )
    SMomentSeries.attach_norm_spin_components_from_spin_components(
        result, tol=tol
    )
    SMomentSeries.attach_angles_from_norm_spin_components(
        result, tol=tol
    )
    return result
```

`_validate_reference_data(...)` enforces the current time-, group-, and angle
shape requirements. `_build_s_generator_coefficients(...)` implements the
analytic generator mapping, while `_build_s_basis(...)` stores the equivalent
effective-qubit basis. `compute_s_moments(...)` contracts the averaged first
moments, constructs the raw output, and delegates nonlinear fields to the
`SMomentSeries` class methods.

= Output

Return one group-resolved `SMomentSeries`:

```python
SMomentSeries(
    t=array,
    integration_phase_index=None or array,
    x=tuple[array, ...],
    y=tuple[array, ...],
    z=tuple[array, ...],
    length=tuple[array, ...],
    nx=tuple[array, ...],
    ny=tuple[array, ...],
    nz=tuple[array, ...],
    theta=tuple[array, ...],
    phi=tuple[array, ...],
    basis=SBasisSeries(
        zero_state=tuple[array, ...],  # each shape (len(t), 3)
        one_state=tuple[array, ...],   # each shape (len(t), 3)
    ),
)
```

All group fields are tuples even for $G=1$. Store the result as `moments.S`.

= Data Requirements

- `mean_generator_moments`: a `TrajectoryResolvedMomentSeries` with singleton
  trajectory axis and first moments shaped `(1, len(t), G, 8)`.
- `j_moments`: a `JMomentSeries` on the same time grid with `theta` and `phi`
  tuples that form arrays of shape `(G, len(t))`.
- If supplied, `j_moments.integration_phase_index` must have shape `(len(t),)`.
- A strictly positive numerical tolerance.

The intended pipeline computes `j_moments` from the same selected mean
generator moments before calling this function. The current validation checks
shape and time-grid consistency but cannot verify this provenance.

= Invariants and Edge Cases

- Average or resample complete trajectories before computing J and S moments.
- Recompute J angles and then S moments for every resampled trajectory set.
- Preserve the effective basis phase $e^(-i phi_J)$ and all signs in the
  generator-to-S mapping.
- Preserve $S_z=(N_u-N_1)/2$ and the S-sphere $p=+1$ pole convention.
- Return group tuples for every group count; the generic array implementation
  is not restricted to one or two groups.
- Copy the optional integration-phase index from `j_moments`; leave it `None`
  when the J series does not contain it.
- Do not aggregate S components across groups with different local frames.
- Do not transform symmetric second moments or compute covariance or squeezing.
- Require matching time grids within absolute tolerance $10^(-12)$ and reject
  nonpositive `tol`.
