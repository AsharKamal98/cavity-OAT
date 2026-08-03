#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[
    Solver-Independent Generalized Squeezing
  ]
]

= Purpose

This file specifies how `post_analysis/squeezing.py` converts precomputed
group-local `SqueezingMomentSeries` objects into one `SqueezingSeries` per
physical group.

Solver-specific extraction of MCWF squeezing moments is defined separately in
`docs/instructions/solvers/mcwf/squeezing_moments.typ`.

= Calculation

For one group, let $mu_a$ be the four mean fluctuation operators and let
$M_(a b)$ be their symmetric second-moment matrix (see one of solver-specific `squeezing_moments.typ`instruction files for definitionss). Construct

$
C_(a b) = M_(a b) - mu_a mu_b.
$

Use `numpy.linalg.eigvalsh` to obtain the ordered covariance eigenvalues and
take the smallest,

$
lambda_"min" = lambda_1(C).
$

The generalized squeezing parameter is

$
xi^2
=
frac(
  N_g lambda_"min",
  (N_c / 2)^2
).
$

Where the denominator is below tolerance, return `NaN` for $xi^2$. Convert
positive values to decibels using

$
xi_"dB"^2 = 10 log_10(xi^2),
$

and leave non-positive values as `NaN` in the dB series.

= Function Flow

The public entry point is

```python
compute_squeezing(
    squeezing_moments: tuple[SqueezingMomentSeries, ...],
    *,
    tol=1e-12,
) -> tuple[SqueezingSeries, ...]
```

It should:

```python
# 1. Require at least one input and a positive tolerance.
# 2. Validate that all groups share the same time grid.
# 3. For each group:
#       construct C from second_moments and mean_operators;
#       compute covariance_eigvals and lambda_min;
#       compute xi2 and xi2_db;
#       copy the group metadata and reference-angle fields.
# 4. Return the group tuple in the original order.
```

= Output

Each returned group item is

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
    theta_J,
    phi_J,
    theta_S,
    phi_S,
)
```

Store the tuple as `moments.squeezing`.

= Data Requirements

`compute_squeezing(...)` requires one or more ensemble-averaged
`SqueezingMomentSeries` objects on a common time grid.

= Invariants

- `post_analysis/squeezing.py` must remain solver-independent.
- Do not import MCWF result containers or inspect saved wavefunctions.
- Do not average covariance matrices or squeezing values across trajectories;
  the input raw moments must already be ensemble averaged.
- Preserve the physical-group ordering of the input tuple.
- Do not compute full-system or cross-group squeezing.
