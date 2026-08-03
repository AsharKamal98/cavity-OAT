#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[
    Spin-Moment Conventions
  ]
]

= Purpose

This file defines shared conversions between Cartesian spin components,
normalized spin directions, and Bloch-sphere angles. Use these conventions for
any spin sphere after specifying which z pole corresponds to $theta=0$.

= Raw Spin Components

For a spin vector with Cartesian components $S_x$, $S_y$, and $S_z$, calculate
all nonlinear fields only after any required trajectory or ensemble averaging.
Do not average per-sample vector lengths, normalized directions, or angles.

When group spin operators are additive and share the same coordinate frame,
construct an aggregated spin vector by summing the raw components:

$
S_(a,"aggregated") = sum_g S_(a,g),
quad a in {x,y,z}.
$

Derive the aggregated direction and angles afterward. Do not sum vector
lengths, normalized directions, or angles. A sphere-specific instruction must
state whether aggregation is physically meaningful.

= Normalized Spin Direction

Define

$
S_"len" = sqrt(S_x^2 + S_y^2 + S_z^2),
quad
n_a = frac(S_a, S_"len"),
quad a in {x,y,z}.
$

Use `norm_spin_components_from_spin_components(...)` from
`common/utils/moments.py`. If $S_"len"$ is below the numerical tolerance, set
the normalized components to zero.

= Bloch-Sphere Angles

Let $p$ identify the z pole at $theta=0$:

$
p =
cases(
  +1 & "for the +z pole",
  -1 & "for the -z pole",
).
$

Then

$
theta = arccos(p n_z),
quad
phi = "atan2"(n_y, n_x).
$

The inverse relations are

$
n_x = sin(theta) cos(phi),
quad
n_y = sin(theta) sin(phi),
quad
n_z = p cos(theta).
$

The current `angles_from_norm_spin_components(...)` helper implements
$p=-1$ directly. For a sphere with $p=+1$, pass `-nz` to that helper while
retaining the physical `nz` in the output container.

The repository conventions are:

- J sphere: $p=-1$, so $theta_J=0$ points along the down or south pole.
- S sphere defined by $S_z=(N_u-N_1)/2$: $p=+1$, so $theta_S=0$ points along
  $|u chevron.r$.

If the spin length is below tolerance, set both angles to zero.

= Invariants

- Average raw Cartesian components before calculating nonlinear fields.
- Normalize directions by the actual spin-vector length, not by an atom
  population.
- State the sphere's pole convention wherever these shared conversions are
  applied.
- Aggregate only additive vectors expressed in a shared coordinate frame.
- Use the helpers in `common/utils/moments.py` rather than duplicating their
  numerical tolerance logic.
