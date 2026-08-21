#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Gell--Mann Moment Rotation]
]

= Purpose

This file specifies `compute_gell_mann_rotation(...)` and
`rotate_generator_moments(...)` in
`post_analysis/generator_rotation.py`. They convert stored lab-basis SU(3)
moments into the dressed basis by transforming the moment arrays directly.
The generator normalization and packed storage convention are defined in
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`.

= Theory

Use $g_a=lambda_a/2$ with

$
"Tr"(g_a g_b) = frac(delta_(a b), 2).
$

For the transformation from
`docs/instructions/post_analysis/dressed_basis.typ`, define

$
T_(alpha,g) = U_g g_alpha U_g^dagger
= sum_(a=1)^8 R_(alpha a,g) g_a,
quad
R_(alpha a,g) = 2 "Tr"(U_g g_alpha U_g^dagger g_a).
$

The adjoint matrix $R_g$ is real and orthogonal. For $G$ groups, let
$cal(R) = "diag"(R_1, dots.c, R_G)$ be block diagonal, using the flattened
generator index $(g,a)$. For every retained trajectory $r$,

$
tilde(mu)^((r)) = cal(R) mu^((r)),
quad
tilde(M)^((r)) = cal(R) M^((r)) cal(R)^T.
$

Here $M$ is the full symmetric anticommutator-moment matrix. These are moment
rotations, not covariance rotations; no product of first moments is
subtracted.

= Method in Pseudo-Code

```python
def compute_gell_mann_rotation(
    transformation: Array,
    *,
    tol: float = 1e-12,
) -> Array:
    require U with shape (len(t), G, 3, 3) and U_dagger @ U == identity
    construct g_a = lambda_a / 2 in the lab basis
    T_alpha = U @ g_alpha @ U_dagger
    R[alpha, a] = 2 * trace(T_alpha @ g_a)
    require negligible imaginary residual
    require R @ swapaxes(R, -1, -2) == identity
    return real R with shape (len(t), G, 8, 8)

def rotate_generator_moments(
    moments: TrajectoryResolvedMomentSeries,
    rotation: Array,
    *,
    tol: float = 1e-12,
) -> TrajectoryResolvedMomentSeries:
    require matching time/group axes and orthogonal local rotations
    rotate first_moments on each group-local generator axis
    unpack second_moments into the full symmetric matrix
    rotate with the block-diagonal global matrix and repack its upper triangle
    return a new TrajectoryResolvedMomentSeries on a copied time grid
```

`compute_gell_mann_rotation(...)` converts the three-level basis change into
the eight-generator adjoint representation. `rotate_generator_moments(...)`
then applies it to every retained trajectory and both stored moment orders.

= Invariants

- Preserve the trajectory axis; this function does not average trajectories.
- Preserve all cross-group second-moment blocks.
- Preserve the group-major flattened order `8 * group + generator` and the
  `np.triu_indices(8 * G)` packing convention.
- Use one shared time- and group-resolved rotation for first and second moments.
- Do not store `U` or `R` in the returned moment series.
- Compute J and S moments from the lab-frame mean before rotating. Do not pass
  the dressed result back to lab-mapping functions such as
  `compute_j_moments(...)` or `compute_s_moments(...)`.
