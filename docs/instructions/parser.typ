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
especially `parser/j_moments.py` and `parser/moments.py`.

General parser-class style rules live in the `parser` skill. This file records
the repository-specific `JMomentSeries` contract.

= J-Moment Containers

`JMomentSnapshot` stores first-order J-sphere moments for one saved MCWF
trajectory snapshot:

```python
JMomentSnapshot(
    t,
    phase_index,
    x, y, z,
    N_e,
    N_j,
    jump_rate,
    x_groups=None,
    y_groups=None,
    z_groups=None,
    N_e_groups=None,
    N_j_groups=None,
)
```

`JMomentSeries` stores one time series, either for one trajectory, an ensemble
average, or a non-MCWF method converted into the shared J-sphere format:

```python
JMomentSeries(
    t,
    phase_index=None,
    x=None, y=None, z=None,
    x_groups=None, y_groups=None, z_groups=None,
    length=None, length_groups=None,
    nx=None, ny=None, nz=None,
    nx_groups=None, ny_groups=None, nz_groups=None,
    N_e=None, N_j=None,
    N_e_groups=None, N_j_groups=None,
    theta=None, phi=None,
    theta_groups=None, phi_groups=None,
    jump_rate=None,
)
```

Group-resolved fields should be tuples ordered by group. If no group-resolved
data exist, the corresponding fields should remain `None`.

= Derived-Field Class Methods

The `JMomentSeries` class methods should attach missing equivalent
representations in place:

```python
JMomentSeries.attatch_norm_spin_components_from_spin_components(j_moments, tol=...)
    -> attach length, nx, ny, nz from x, y, z

JMomentSeries.attatch_angles_from_norm_spin_components(j_moments, tol=...)
    -> attach theta, phi from length, nx, ny, nz

JMomentSeries.attatch_norm_spin_components_from_angles(j_moments)
    -> attach nx, ny, nz from theta, phi

JMomentSeries.attatch_spin_components_from_norm_spin_components(j_moments)
    -> attach x, y, z from length, nx, ny, nz
```

The same methods should also attach group-resolved fields when the required
group inputs are present. If neither full-system nor group-resolved inputs are
available for the requested conversion, the method should raise a clear
`ValueError`.

The conversion helpers used by these class methods live in
`common/utils/moments.py`.

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

= Invariants

- Parser classes should stay as Pydantic containers with explicit typed fields.
- Derived-field class methods should mutate the supplied `JMomentSeries`
  instance in place.
- Full-system and group-resolved conversions should use the same formulas.
- J-moment extraction code should call these class methods instead of
  duplicating spin-component, direction, or angle conversion logic.
- Legacy field names such as `Jx`, `Jy`, `Jz`, `J_len`, and `sx`, `sy`, `sz`
  should not be reintroduced.
