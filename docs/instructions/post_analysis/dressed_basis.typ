#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Dressed-Basis Transformation]
]

= Purpose

This file specifies `build_dressed_basis_transformation(...)` in
`post_analysis/dressed_basis.py`. It constructs the group- and time-resolved
dressed basis from the J and S angles used by post-analysis.
`docs/theory/appendix_holstein_primakoff.tex` provides physical context for the
dressed modes; the exact numerical gauge is defined below.

= Theory

In the lab order $(|u chevron.r, |d chevron.r, |e chevron.r)$, define

$
|1_g chevron.r
= cos(theta_(J,g)/2)|d chevron.r
+ e^(-i phi_(J,g)) sin(theta_(J,g)/2)|e chevron.r,
$

$
|j_g chevron.r
= -sin(theta_(J,g)/2)|d chevron.r
+ e^(-i phi_(J,g)) cos(theta_(J,g)/2)|e chevron.r,
$

$
|c_g chevron.r
= cos(theta_(S,g)/2)|u chevron.r
+ e^(-i phi_(S,g)) sin(theta_(S,g)/2)|1_g chevron.r,
$

$
|s_g chevron.r
= -sin(theta_(S,g)/2)|u chevron.r
+ e^(-i phi_(S,g)) cos(theta_(S,g)/2)|1_g chevron.r.
$

The transformation is

$
U_g = (|c_g chevron.r quad |j_g chevron.r quad |s_g chevron.r),
quad U_g^dagger U_g = I_3.
$

Its columns are dressed states in lab coordinates. Therefore $U_g$ maps
dressed coordinates to lab coordinates, while $U_g^dagger$ maps lab
coordinates to dressed coordinates. Preserve the stated phases and signs;
they fix the dressed-state gauge used by later observables.

= Method in Pseudo-Code

```python
def build_dressed_basis_transformation(
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol: float = 1e-12,
) -> Array:
    require matching nonempty time grids and group-resolved J/S angles
    construct |1>, |j>, |c>, and |s> from the equations above
    U = stack (|c>, |j>, |s>) as columns
    require U_dagger @ U == identity within tol
    return U with shape (len(t), G, 3, 3)
```

The function validates the shared frame data, constructs the exact dressed
states, and returns one unitary transformation for every time and group.

= Invariants

- Use the lab-state order $(|u chevron.r, |d chevron.r, |e chevron.r)$.
- Use the dressed column order $(|c chevron.r, |j chevron.r, |s chevron.r)$.
- Do not rephase or orthogonalize the analytically constructed columns.
- Keep `U` as a transient post-analysis object rather than simulation output.
