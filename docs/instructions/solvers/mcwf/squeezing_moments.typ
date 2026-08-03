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

At each saved time, group $g$ uses its group-resolved ensemble-averaged
$J$ angles $theta_(J,g), phi_(J,g)$ from `JMomentSeries`. Define

$
|1_g chevron.r
=
cos(theta_(J,g)/2)|d chevron.r
+
e^(-i phi_(J,g)) sin(theta_(J,g)/2)|e chevron.r.
$

Reuse the group-resolved ensemble-averaged $S$ angles
$theta_(S,g), phi_(S,g)$ from `SMomentSeries`, as computed according to
`docs/instructions/solvers/mcwf/s_moments.typ`. Then define

$
|c_g chevron.r
=
cos(theta_(S,g)/2)|u chevron.r
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
-sin(theta_(S,g)/2)|u chevron.r
+
e^(-i phi_(S,g)) cos(theta_(S,g)/2)|1_g chevron.r.
$

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

and lift them to collective group operators

$
O_(a,g) = sum_(i in g) o_(a,g)^((i)).
$

For every trajectory and saved time, extract

$
mu_(a,g) = chevron.l O_(a,g) chevron.r,
quad
M_(a b,g)
=
"Re" chevron.l O_(a,g) psi | O_(b,g) psi chevron.r,
$

and

$
N_(c,g)
=
chevron.l sum_(i in g) |c_g chevron.r_i
chevron.l c_g|_i chevron.r.
$

Use the four vectors $O_(a,g)|psi chevron.r$ to obtain $M_(a b,g)$; do not
construct products $O_(a,g)O_(b,g)$ explicitly.

= Function Flow

The public entry point is

```python
compute_mcwf_squeezing_moments(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol=1e-12,
    n_processes=None,
) -> tuple[SqueezingMomentSeries, ...]
```

Its workflow is:

```python
# 1. Validate the shared t_eval grid and one- or two-group sector structure.
# 2. Reuse theta_J[g], phi_J[g], N_e[g], and N_j[g] from j_moments.
# 3. Reuse theta_S[g] and phi_S[g] from s_moments.
# 4. Extract mu_a,g, M_ab,g, and N_c,g for every trajectory.
# 5. Average these raw moments over trajectories.
# 6. Return one SqueezingMomentSeries per group.
```

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
  `theta` and `phi`.
- Scalar one-group sectors or two-group tuple sectors.

= Invariants

- Moment extraction must not affect MCWF propagation.
- Return one independent moment series per physical group.
- Average raw moments over trajectories before returning them.
- Reuse `SMomentSeries.theta` and `SMomentSeries.phi`; do not recompute S
  components or S angles inside squeezing-moment extraction.
- Do not construct covariance eigenvalues or squeezing parameters here.
- Do not compute full-system or cross-group squeezing.
- Do not add these moments to `JMomentSeries`.
- Do not construct full $3^N$ tensor-product operators.
