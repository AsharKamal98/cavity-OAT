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

`post_analysis/squeezing.py` converts already-averaged dressed generator
moments into generalized squeezing. It always computes one result per physical
group and optionally computes the squeezing of the unweighted sum of all group
fluctuation operators. It consumes dressed generator moments, J moments, and
shared simulation metadata. The caller supplies the explicitly rotated dressed
moments.

= Theory

Let $g_a=lambda_a/2$ be the normalized Gell--Mann generators. The J and S
angles define the time- and group-dependent dressed basis
$(|c chevron.r, |j chevron.r, |s chevron.r)$. Before this function is called,
construct $U$ and $R$ and rotate both moment orders as specified in
`docs/instructions/post_analysis/dressed_basis.typ` and
`docs/instructions/post_analysis/generator_rotation.typ`. Denote the resulting
dressed moments by $tilde(mu)_(g alpha)$ and
$tilde(M)_(g alpha, h beta)$. The four squeezing fluctuation generators in
each group are

$
(T_1, T_2, T_4, T_5).
$

For $a_i in (1, 2, 4, 5)$, retain every group block of the covariance matrix,

$
C^(g h)_(i j)
=
tilde(M)_(g a_i, h a_j)
-
tilde(mu)_(g a_i) tilde(mu)_(h a_j).
$

The group-resolved covariance is the diagonal block $C^(g g)$. The dressed
$c$-state population is recovered without another operator evaluation:

$
N_(c,g)
=
N_g/3 + tilde(mu)_(g 3) + tilde(mu)_(g 8)/sqrt(3).
$

Let $lambda_"min"$ be the smallest eigenvalue of the relevant four-by-four
covariance. For a subsystem containing $N$ atoms, define

$
xi^2
=
frac(N lambda_"min", (N_c/2)^2),
quad
xi_"dB"^2 = 10 log_10(xi^2).
$

Return `NaN` for $xi^2$ when its denominator is below tolerance and for the dB
value when $xi^2$ is non-positive.

When full-system squeezing is requested, define each full-system fluctuation
operator as the unweighted group sum. Its covariance and dressed population
are therefore

$
C^"tot" = sum_(g,h) C^(g h),
quad
N_c^"tot" = sum_g N_(c,g),
quad
N^"tot" = sum_g N_g.
$

Cross-group covariance blocks are essential in this sum. The corresponding
active-manifold excited fraction is
$sum_g N_(e,g) / sum_g N_(j,g)$.

= Function Flow

The public entry point is

```python
compute_squeezing(
    mean_dressed_generator_moments,
    j_moments,
    *,
    metadata,
    include_full_system=False,
    tol=1e-12,
) -> SqueezingResult
```

It should:

```python
# 1. Require averaged dressed moments with a singleton trajectory axis.
# 2. Validate shared time and physical-group axes.
# 3. Form the full 4G x 4G covariance for (T1, T2, T4, T5).
# 4. Diagonalize each 4 x 4 diagonal block and compute group squeezing.
# 5. If include_full_system is True only:
#       sum all covariance blocks and group N_c values;
#       diagonalize the resulting 4 x 4 covariance;
#       compute the full-system squeezing series.
# 6. Return the group tuple and optional full-system series.
```

= Output

`SqueezingSeries` represents either one group or the full system:

```python
SqueezingSeries(
    t,
    N_atoms,
    covariance,                 # shape (n_t, 4, 4)
    covariance_eigvals,         # shape (n_t, 4), ascending
    lambda_min,
    xi2,
    xi2_db,
    N_c,
    excited_fraction_active,
)
```

The public result is

```python
SqueezingResult(
    group_resolved=(...),
    full_system=None or SqueezingSeries(...),
)
```

Store it as `moments.squeezing`. J and S angles remain in `moments.J` and
`moments.S`; do not duplicate them in the squeezing result.

= Invariants

- The function remains solver-independent and never imports MCWF containers.
- The caller must supply moments already rotated into the dressed
  $(|c chevron.r, |j chevron.r, |s chevron.r)$ frame. The function must not
  infer or apply that rotation internally.
- Average complete trajectories before calling it; do not average covariance
  matrices, eigenvalues, or squeezing values across trajectories.
- Preserve all cross-group second moments until the covariance blocks have
  been formed.
- Always compute group-resolved results in metadata group order.
- Do not construct or diagonalize the full-system covariance when
  `include_full_system` is false.
- The full-system definition is the unweighted sum of corresponding local
  fluctuation operators. Weighted collective modes require a separate API.
