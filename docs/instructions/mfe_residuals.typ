#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[MFE Residual Diagnostics]
]

= Purpose
This file specifies how to compute two-group mean-field-equation (MFE)
residuals in `quantum_trajectories/mfe_residuals.py`. Use this file for tasks
related to residual diagnostics or plots that consume `moments.MFE_residuals`.

The MFE residual diagnostic consumes averaged J moments from
`quantum_trajectories/j_moments.py`. It should not recompute J moments or old
observable-series fields.

The main function structure should be:

```python
compute_mfe_residuals(
    j_moments: JMomentSeries,
    *,
    parameters: MomentParameters,
    tol=1e-12,
) -> MFEResidualSeries | None
    "compute two-group MFE residuals from stored J-vector group angles"

attach_mfe_residuals(moments: MomentSeries, *, tol=1e-12)
    -> MFEResidualSeries | None
    moments.MFE_residuals = compute_mfe_residuals(
        moments.J,
        parameters=moments.parameters,
        tol=tol,
    )
    return moments.MFE_residuals
```

= MFE Residual Definitions

For two-group inhomogeneous results, define the weighted collective transverse
sum using the J-vector group angles stored on `moments.J`:

$
C_(J)(t_(k)) =
omega_(1) N_(J,1)(t_(k)) e^(i phi_(J,1)(t_(k))) sin(theta_(J,1)(t_(k)))
+ omega_(2) N_(J,2)(t_(k)) e^(i phi_(J,2)(t_(k))) sin(theta_(J,2)(t_(k))).
$

The residual for group $a in {1,2}$ is the left side minus the right side of
the two-group steady-state equation:

$
R_(a)(t_(k)) =
frac(Omega(t_(k)) omega_(a), 2)
e^(-i phi_(J,a)(t_(k))) sin(theta_(J,a)(t_(k)))
- frac(delta(t_(k)), 2) sin(theta_(J,a)(t_(k))) tan(theta_(J,a)(t_(k)))
+ frac(i Gamma omega_(a), 4)
e^(-i phi_(J,a)(t_(k))) sin(theta_(J,a)(t_(k))) C_(J)(t_(k)).
$

Here $Omega(t_(k))$ and $delta(t_(k))$ are selected from the saved
`phase_index`, while $omega_(a)$ is read from `MomentParameters.omega_groups`.
The group atom-number weights are `moments.J.N_j_groups[a]`, not the fixed
total group sizes.

The L2 norm used by residual plots is

$
norm(R)_2(t_(k)) = sqrt(abs(R_(1)(t_(k)))^2 + abs(R_(2)(t_(k)))^2).
$

= Output

`compute_mfe_residuals(...)` should return `None` when the input does not
contain exactly two group-resolved J-angle and atom-number fields. Otherwise it
should return:

```python
MFEResidualSeries(
    t,
    phase_index,
    residuals_groups=(R_1, R_2),
)
```

The top-level notebook container should store this result as:

```python
moments.MFE_residuals = compute_mfe_residuals(
    moments.J,
    parameters=moments.parameters,
)
```

= Invariants

- MFE residuals should be computed after J moments have been ensemble-averaged
  and after J-vector group angles have been attached.
- The residual calculation should use `theta_groups`, `phi_groups`, and
  `N_j_groups` from `moments.J`.
- The residual calculation should use shared protocol metadata from
  `moments.parameters`.
- The residual calculation should not live in `quantum_trajectories/j_moments.py`
  and should not attach residual fields to `JMomentSeries`.
