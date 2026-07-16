#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Numerical MFE Solver]
]

= Purpose
This file specifies the preferred structure for numerically solving the
mean-field equations (MFEs). Use it when implementing the standalone MFE solver
outside the MCWF backend, for example in `solvers/mfe/sim.py`, or when writing
diagnostics that consume its outputs.

This file is a structural implementation guide. For the theory derivation, use
`docs/theory/notes_inhomogeneous_coupling.tex`. For residual diagnostics that
consume solved or simulated angles, use `docs/instructions/post_analysis/mfe_residuals.typ`.

= MFE Definitions

The numerical solver should integrate the complex mean-field amplitudes
$D_(a)(t)$ and $E_(a)(t)$ for groups $a=1,dots,G$. The group count should be an
input; homogeneous runs are the special case $G=1$ with $omega_(1)=1$.

For phase-local drive $Omega(t)$ and detuning $delta(t)$, solve

$
partial_t D_(a) =
- frac(i Omega(t) omega_(a), 2) E_(a)
+ frac(Gamma omega_(a), 2)
  (sum_(b=1)^G omega_(b) E_(b)^* D_(b)) E_(a),
$

$
partial_t E_(a) =
- frac(i Omega(t) omega_(a), 2) D_(a)
+ i delta(t) E_(a)
- frac(Gamma omega_(a), 2) D_(a)
  (sum_(b=1)^G omega_(b) D_(b)^* E_(b)).
$

At saved time $t_(k)$, the amplitudes are related to the group-resolved J-sphere
angles by

$
D_(a)(t_(k)) = sqrt(N_(J,a)) cos(frac(theta_(J,a)(t_(k)), 2)),
$

$
E_(a)(t_(k)) =
sqrt(N_(J,a)) e^(-i phi_(J,a)(t_(k))) sin(frac(theta_(J,a)(t_(k)), 2)).
$

The inverse relations used for output should be

$
N_(J,a)(t_(k)) = abs(D_(a)(t_(k)))^2 + abs(E_(a)(t_(k)))^2,
$

$
theta_(J,a)(t_(k)) =
arccos(frac(abs(D_(a)(t_(k)))^2 - abs(E_(a)(t_(k)))^2, N_(J,a)(t_(k)))),
$

$
phi_(J,a)(t_(k)) = arg(D_(a)(t_(k))) - arg(E_(a)(t_(k))).
$

After taking this phase difference, wrap it back to the principal interval
$[-pi, pi]$ before storing or plotting it.

The angle convention matches the J-vector convention used by
`docs/instructions/j_moments.typ`: all atoms in $|d chevron.r$ correspond to
$theta_(J,a)=0$.

The corresponding group-resolved spin components should be

$
J_(x,a)(t_(k)) =
frac(N_(J,a)(t_(k)), 2) sin(theta_(J,a)(t_(k))) cos(phi_(J,a)(t_(k))),
$

$
J_(y,a)(t_(k)) =
+ frac(N_(J,a)(t_(k)), 2) sin(theta_(J,a)(t_(k))) sin(phi_(J,a)(t_(k))),
$

$
J_(z,a)(t_(k)) =
- frac(N_(J,a)(t_(k)), 2) cos(theta_(J,a)(t_(k))),
$

with

$
|J_(a)|(t_(k)) =
sqrt(J_(x,a)(t_(k))^2 + J_(y,a)(t_(k))^2 + J_(z,a)(t_(k))^2).
$

The minus sign in $J_(z,a)$ follows the same repository convention as the
J-moment pipeline: $theta_(J,a)=0$ corresponds to the all-$|d chevron.r$
state.

The group excited-state population is $N_(e,a)(t_(k)) =
abs(E_(a)(t_(k)))^2$. Full additive fields are sums over groups, following the
`JMomentSeries` convention in `docs/instructions/parser.typ`.

= Method in Pseudo-code

The solver should be split into small functions with pure data flow. The core
solver should not import `solvers.mcwf`. The shared `PhaseProtocol` convention
is defined in `docs/instructions/phases.typ`.

```python
def mfe_rhs(t, y, parameters, integration_phase) -> Array:
    G = parameters.group_count
    D, E = y[:G], y[G:]
    omega = np.asarray(parameters.omega_i)
    Omega_t = integration_phase.omega
    delta_t = integration_phase.delta

    ED = sum(omega_b * np.conj(E_b) * D_b for omega_b, E_b, D_b in zip(omega, E, D))
    DE = sum(omega_b * np.conj(D_b) * E_b for omega_b, D_b, E_b in zip(omega, D, E))

    dD = -0.5j * Omega_t * omega * E + 0.5 * parameters.Gamma * omega * ED * E
    dE = -0.5j * Omega_t * omega * D + 1j * delta_t * E - 0.5 * parameters.Gamma * omega * D * DE
    return np.concatenate([dD, dE])
```

```python
def solve_mfe(parameters, *, t_eval, rtol=1e-9, atol=1e-11) -> MFEResult:
    zero_angles = (0.0,) * parameters.group_count
    y0 = amplitudes_from_initial_state(
        zero_angles, zero_angles, parameters
    )
    integration_phases = parameters.phase_protocol.integration_phases
    integration_boundaries = phase_boundary_times(integration_phases)

    def rhs(t, y):
        phase_index = searchsorted(integration_boundaries, t, side="left")
        return mfe_rhs(t, y, parameters, integration_phases[phase_index])

    solution = solve_ivp(
        rhs,
        (t_eval[0], t_eval[-1]),
        y0,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
    )
    G = parameters.group_count
    D_groups, E_groups = solution.y[:G], solution.y[G:]
    return MFEResult(t=t_eval, D_groups=D_groups, E_groups=E_groups, ...)
```

```python
def compute_mfe_j_moments(result: MFEResult, *, tol=1e-12) -> JMomentSeries:
    N_j_groups, theta_groups, phi_groups = angles_from_amplitudes(
        result.D_groups,
        result.E_groups,
    )
    j_moments = JMomentSeries(
        result.t,
        integration_phase_index=integration_phase_indices_at_times(...),
        N_e_groups=tuple(abs(E_g)**2 for E_g in result.E_groups),
        N_j_groups=N_j_groups,
        theta_groups=theta_groups,
        phi_groups=phi_groups,
        ...,
    )
    JMomentSeries.attatch_norm_spin_components_from_angles(j_moments)
    JMomentSeries.attatch_spin_components_from_norm_spin_components(j_moments)
    JMomentSeries.attatch_additive_full_fields_from_group_fields(j_moments)
    JMomentSeries.attatch_norm_spin_components_from_spin_components(j_moments, tol=tol)
    JMomentSeries.attatch_angles_from_norm_spin_components(j_moments, tol=tol)
    return j_moments
```

Run one adaptive solve over the complete protocol. Flatten the integration
`Phase` sequence and compute its boundaries once before calling `solve_ivp`;
the RHS selects the active integration `Phase` from the current time.

Undefined helper notes:

- `amplitudes_from_initial_state(...)` is a new local helper that constructs
  `(D_groups, E_groups)` from `(theta_groups, phi_groups, N_(J))`, with
  `N_(J,a) = N_(a) / 2`.
- `angles_from_amplitudes(...)` is a new local helper that applies the inverse
  relations in the MFE definitions section.

Function flow: `solve_mfe(...)` is the main entry point. It calls
`amplitudes_from_initial_state(...)` to build the initial solver vector,
`solve_ivp(...)` to integrate `mfe_rhs(...)`.
`compute_mfe_j_moments(...)` is the post-processing step that converts an
`MFEResult` into a `JMomentSeries` containing group-resolved and full-system
fields.
`mfe_rhs(...)` receives the current integration `Phase` explicitly.

If a top-level `MomentSeries` container is already in use,
store the solved observable series explicitly as `moments.J`, for example via
`moments.J = compute_mfe_j_moments(result)`.

= Data Requirements

The solver needs:

- `metadata`: `SimulationMetadata` with `Ni`, the first `G-1` `omega_i`
  inputs, `Gamma`, and a supplied `phase_protocol`. Its validator supplies the
  full `omega_groups` vector used by the MFE equations;
- `t_eval`: saved output times.

The metadata validator completes the final coupling before `solve_mfe(...)`
runs, using the same weighted-coupling convention as the MCWF code.

= Output

Important solver-facing data should be stored in Pydantic classes. In the
current layout, solver-input and raw-solution classes live in `parser/mfe.py`,
while the observable output container is `parser/j_moments.py`. A suggested
minimal output structure is:

```python
MFESolverParameters(
    Gamma,
    phase_protocol,
    omega_i,
    Ni,
)

MFEResult(
    t,
    D_groups,
    E_groups,
    success,
    message,
    parameters,
)

JMomentSeries(
    t,
    integration_phase_index,
    N_e,
    N_j,
    theta,
    phi,
    x,
    y,
    z,
    length,
    nx,
    ny,
    nz,
    N_e_groups,
    N_j_groups,
    theta_groups,
    phi_groups,
    x_groups,
    y_groups,
    z_groups,
    length_groups,
    nx_groups,
    ny_groups,
    nz_groups,
)
```

The observable series should be the object consumed by plotting or residual
diagnostics. Do not require plotting functions to unpack raw `solve_ivp`
objects.

= Invariants

- The standalone MFE solver should not import from `solvers.mcwf`.
- Use the shared `PhaseProtocol` and integration `Phase` classes; do not define
  solver-specific phase classes.
- The solver should support arbitrary group count $G$ when the equations are
  written as sums over groups. Two-group-only logic belongs in residual
  diagnostics unless the theory explicitly requires two groups.
- Do not recompute coupling weights from instantaneous $N_(J,a)(t_(k))$. Coupling
  weights are fixed simulation inputs.
- The per-group norm $abs(D_(a))^2 + abs(E_(a))^2$ should remain constant up to
  solver tolerance. Large drift should be treated as a numerical warning.
- Phase and angle conventions should match `docs/instructions/j_moments.typ`
  and `docs/instructions/post_analysis/mfe_residuals.typ`.
