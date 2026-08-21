#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[MCWF Squeezing Post-Analysis]
]

= Purpose

MCWF squeezing is computed from selected mean generator moments through the
same solver-independent post-analysis pipeline used by any simulator that
provides the shared moment format.

The dressed-state construction, generator rotation, and squeezing calculation
are specified in:

- `docs/instructions/post_analysis/dressed_basis.typ`,
- `docs/instructions/post_analysis/generator_rotation.typ`, and
- `docs/instructions/post_analysis/squeezing.typ`.

= Theory

J and S angles define the group- and time-dependent dressed basis
$(|c chevron.r, |j chevron.r, |s chevron.r)$. If $U_g$ contains those dressed
states and $g_a=lambda_a/2$, the dressed generators and their adjoint rotation
are

$
T_(g alpha) = U_g g_alpha U_g^dagger,
quad
T_(g alpha) = sum_a R_(g alpha,a) g_a,
$

$
R_(g alpha,a) = 2 "Tr"(U_g g_alpha U_g^dagger g_a).
$

Apply the corresponding block-diagonal rotation to both selected mean moment
orders. The squeezing covariance is then formed from the dressed first and
symmetric second moments. Group squeezing uses each diagonal four-by-four
fluctuation block; optional full-system squeezing sums all group-pair blocks
before diagonalization.

= Method in Pseudo-Code

```python
U = build_dressed_basis_transformation(moments.J, moments.S)
R = compute_gell_mann_rotation(U)

mean_dressed_generator_moments = rotate_generator_moments(
    moments.mean_generator_moments,
    R,
)

moments.squeezing = compute_squeezing(
    mean_dressed_generator_moments,
    moments.J,
    metadata=moments.metadata,
    include_full_system=include_full_system,
)
```

= Output

Store the `SqueezingResult` as `MomentSeries.squeezing`. Its
`group_resolved` tuple always contains one `SqueezingSeries` per physical
group. Its `full_system` field is populated only when
`include_full_system=True`.

= Invariants

- Average or resample complete trajectories before constructing J and S
  angles, the dressed basis, covariance matrices, or squeezing values.
- Rotate both the first and symmetric second selected mean moments with the
  same group- and time-dependent adjoint rotation.
- Preserve cross-group symmetric second moments through the rotation and
  covariance construction.
- Supply dressed moments explicitly to `compute_squeezing(...)`.
- Compute the full-system covariance only when requested.
