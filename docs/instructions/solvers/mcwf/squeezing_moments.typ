#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#let ket(x) = math.equation(block: false, $|#x chevron.r$)
#let bra(x) = math.equation(block: false, $chevron.l #x|$)

#align(center)[
  #text(size: 1.6em, weight: "bold")[
    MCWF Squeezing-Moment Extraction
  ]
]

= Purpose

The `SqueezingMomentSeries` objects produced here are passed as input to the
solver-independent squeezing code described in
`docs/instructions/post_analysis/squeezing.typ`.

This file specifies how `solvers/mcwf/squeezing_moments.py` combines saved
MCWF sector-block wavefunctions with precomputed `JMomentSeries` and
`SMomentSeries` inputs to produce one ensemble-averaged
`SqueezingMomentSeries` per physical group.

= Squeezing Moment Definitions

== Group-Local Dressed Basis

At each saved time, group $g$ reuses the effective states
$|0_g chevron.r$ and $|1_g chevron.r$ stored in
`SMomentSeries.basis.zero_state` and `SMomentSeries.basis.one_state`. Their
construction from the ensemble-averaged J angles is defined in
`docs/instructions/solvers/mcwf/s_moments.typ`; do not reconstruct
$|1_g chevron.r$ here.

Reuse the group-resolved ensemble-averaged $S$ angles
$theta_(S,g), phi_(S,g)$ from `SMomentSeries`, as computed according to
`docs/instructions/solvers/mcwf/s_moments.typ`. Then define

$
|c_g chevron.r
=
cos(theta_(S,g)/2)|0_g chevron.r
+
e^(-i phi_(S,g)) sin(theta_(S,g)/2)|1_g chevron.r,
$

$
|j_g chevron.r
=
-sin(theta_(J,g)/2)|d chevron.r
+
e^(-i phi_(J,g)) cos(theta_(J,g)/2)|e chevron.r,
$

$
|s_g chevron.r
=
-sin(theta_(S,g)/2)|0_g chevron.r
+
e^(-i phi_(S,g)) cos(theta_(S,g)/2)|1_g chevron.r.
$

Construct these once per group and saved time.

== Fluctuation Moments

Construct the four group-local single-particle operators for atom $i$

$
o_(1,g)^((i)) =
frac(ket(c_g) bra(j_g) + ket(j_g) bra(c_g), 2),
quad
o_(2,g)^((i)) =
frac(ket(c_g) bra(j_g) - ket(j_g) bra(c_g), 2i),
$

$
o_(3,g)^((i)) =
frac(ket(c_g) bra(s_g) + ket(s_g) bra(c_g), 2),
quad
o_(4,g)^((i)) =
frac(ket(c_g) bra(s_g) - ket(s_g) bra(c_g), 2i),
$

Construct the fluctuation operators once per group and saved time. Lift them to collective group operators

$
O_(a,g) = sum_(i in g) o_(a,g)^((i)).
$

Then, for every trajectory $r$ and saved time $t_k$, extract normalized moments

$
mu_(a,g)^((r))(t_k)
=
frac(
  chevron.l psi_r(t_k)|O_(a,g)|psi_r(t_k) chevron.r,
  Z_r(t_k)
),
quad
M_(a b,g)^((r))(t_k)
=
frac(
  "Re" chevron.l O_(a,g) psi_r(t_k)
  | O_(b,g) psi_r(t_k) chevron.r,
  Z_r(t_k)
),
$

and

$
N_(c,g)^((r))(t_k)
=
frac(
  chevron.l psi_r(t_k) |
  sum_(i in g) |c_g chevron.r_i chevron.l c_g|_i
  | psi_r(t_k) chevron.r,
  Z_r(t_k)
).
$

These are raw per-trajectory moments. After normalizing each trajectory by its
own $Z_r(t_k)$, average equally over the $R$ trajectories:

$
mu_(a,g)(t_k)
=
frac(1, R) sum_(r=1)^R mu_(a,g)^((r))(t_k),
quad
M_(a b,g)(t_k)
=
frac(1, R) sum_(r=1)^R M_(a b,g)^((r))(t_k),
$

$
N_(c,g)(t_k)
=
frac(1, R) sum_(r=1)^R N_(c,g)^((r))(t_k).
$

Store these ensemble averages as `mean_operators`, `second_moments`, and
`N_c`, respectively, in the group `SqueezingMomentSeries`.

Use the four vectors $O_(a,g)|psi chevron.r$ to obtain $M_(a b,g)$; do not
construct products $O_(a,g)O_(b,g)$ explicitly.

The effective basis and all J/S angles are ensemble-averaged inputs. Therefore
$|c_g chevron.r$, $|j_g chevron.r$, $|s_g chevron.r$, the four fluctuation
operators, and the $|c_g chevron.r chevron.l c_g|$ projector are identical for
all trajectories at a given group and saved time. Construct these shared
series before starting trajectory workers.

= Pseudo-code

Keep the extraction compact and separate shared construction from
trajectory-dependent expectation values:

```python
_build_dressed_state_series(
    zero_state_groups,
    one_state_groups,
    theta_J_groups,
    phi_J_groups,
    theta_S_groups,
    phi_S_groups,
) -> tuple[Array, Array, Array]
    # Return shared c, j, and s arrays, each shaped (G, n_t, 3).

_build_fluctuation_operator_series(
    dressed_state_groups: tuple[Array, Array, Array],
) -> tuple[Array, Array]
    # Return shared fluctuation operators shaped (G, n_t, 4, 3, 3)
    # and N_c projectors shaped (G, n_t, 3, 3).

_trajectory_fluctuation_moments_worker(args) -> tuple[Array, Array, Array]
    # Apply the shared operators to one trajectory and return raw
    # mu, M, and N_c arrays.

compute_mcwf_squeezing_moments(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol=1e-12,
    n_processes=None,
) -> tuple[SqueezingMomentSeries, ...]
    _validate_inputs(ensemble, j_moments, s_moments, tol=tol)

    dressed_state_groups = _build_dressed_state_series(
        s_moments.basis.zero_state,
        s_moments.basis.one_state,
        j_moments.theta,
        j_moments.phi,
        s_moments.theta,
        s_moments.phi,
    )
    fluctuation_operators, N_c_operators = (
        _build_fluctuation_operator_series(dressed_state_groups)
    )

    trajectory_moments = map_with_optional_pool(
        _trajectory_fluctuation_moments_worker(
            trajectory,
            ensemble.metadata.t_eval,
            fluctuation_operators,
            N_c_operators,
            ensemble.metadata.Ni,
            tol,
        )
        for trajectory in ensemble.trajectories
    )
    averaged = average raw mu, M, and N_c over trajectory_moments
    return one SqueezingMomentSeries per group from averaged and input moments
```

`compute_mcwf_squeezing_moments(...)` builds the shared dressed states and
single-particle operators before optional multiprocessing. Each worker applies
those operators to one trajectory. The raw trajectory moments are then
averaged before constructing one `SqueezingMomentSeries` per group.

The reduced-basis transition rules remain in
`solvers/mcwf/squeezing_moments.py`. They support scalar single-group sectors
and the current two-group tuple sectors.

= Output

The returned tuple is ordered by physical group. Each item is

```python
SqueezingMomentSeries(
    t,
    N_group,
    mean_operators,          # shape (n_t, 4)
    second_moments,          # shape (n_t, 4, 4)
    N_c,
    excited_fraction_active,
    theta_J,
    phi_J,
    theta_S,
    phi_S,
)
```

Store the tuple as `moments.squeezing_moments`.

= Data Requirements

- A non-empty `TrajectoryEnsemble` with metadata.
- Saved `sector_blocks` on `ensemble.metadata.t_eval`.
- An ensemble-averaged `JMomentSeries` with group-local `theta`, `phi`, `N_e`,
  and `N_j`.
- An ensemble-averaged `SMomentSeries` on the same time grid with group-local
  `theta`, `phi`, and the stored effective `basis.zero_state` and
  `basis.one_state`.
- Scalar one-group sectors or two-group tuple sectors.

= Invariants

- Moment extraction must not affect MCWF propagation.
- Return one independent moment series per physical group.
- Average raw moments over trajectories before returning them.
- Reuse `SMomentSeries.theta` and `SMomentSeries.phi`; do not recompute S
  components or S angles inside squeezing-moment extraction.
- Reuse `SMomentSeries.basis`; do not reconstruct the effective S basis inside
  squeezing-moment extraction.
- Construct dressed states and their single-particle operators once from
  ensemble-level moment data, before mapping over trajectories.
- Trajectory workers must apply the shared precomputed operators; they must not
  reconstruct dressed states or fluctuation operators.
- Do not construct covariance eigenvalues or squeezing parameters here.
- Do not compute full-system or cross-group squeezing.
- Do not add these moments to `JMomentSeries`.
- Do not construct full $3^N$ tensor-product operators.
