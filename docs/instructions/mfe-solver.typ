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
outside `quantum_trajectories`, for example in `mfe/solver.py`, or when writing
diagnostics that consume its outputs.

This file is a structural implementation guide. For the theory derivation, use
`docs/theory/notes_inhomogeneous_coupling.tex`. For residual diagnostics that
consume solved or simulated angles, use `docs/instructions/mfe_residuals.typ`.

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

The angle convention matches the J-vector convention used by
`docs/instructions/j_moments.typ`: all atoms in $|d chevron.r$ correspond to
$theta_(J,a)=0$.

= Method in Pseudo-code

The solver should be split into small functions with pure data flow. The core
solver should not import `quantum_trajectories`. `Phase` and
`phase_values_at_time(...)` are defined in `common` and follow the convention
described in `docs/instructions/simulation_parameters.typ`.

```python
def mfe_rhs(t: float, y: Array, parameters: MFESolverParameters) -> Array:
    G = parameters.group_count
    D, E = y[:G], y[G:]
    omega = np.asarray(parameters.omega_groups)
    Omega_t, delta_t = phase_values_at_time(t, parameters.phases)

    ED = sum(omega_b * np.conj(E_b) * D_b for omega_b, E_b, D_b in zip(omega, E, D))
    DE = sum(omega_b * np.conj(D_b) * E_b for omega_b, D_b, E_b in zip(omega, D, E))

    dD = -0.5j * Omega_t * omega * E + 0.5 * parameters.Gamma * omega * ED * E
    dE = -0.5j * Omega_t * omega * D + 1j * delta_t * E - 0.5 * parameters.Gamma * omega * D * DE
    return np.concatenate([dD, dE])
```

```python
def solve_mfe(parameters, initial_state, *, t_eval, rtol=1e-9, atol=1e-11) -> MFEResult:
    y0 = amplitudes_from_initial_state(initial_state, parameters)
    solution = solve_ivp(lambda t, y: mfe_rhs(t, y, parameters), ...)
    G = parameters.group_count
    D_groups, E_groups = solution.y[:G], solution.y[G:]
    result = MFEResult(t=t_eval, D_groups=D_groups, E_groups=E_groups, ...)
    result.observables = compute_mfe_observables(result)
    return result
```

```python
def compute_mfe_observables(result: MFEResult, *, tol=1e-12) -> MFEObservableSeries:
    N_j_groups, theta_groups, phi_groups = angles_from_amplitudes(
        result.D_groups,
        result.E_groups,
    )
    return MFEObservableSeries(result.t, result.D_groups, result.E_groups, ...)
```

If the phase protocol has discontinuities, solving phase-by-phase and using the
end state of one phase as the initial state of the next is preferred over
forcing one adaptive solve across discontinuous coefficients.

Undefined helper notes:

- `amplitudes_from_initial_state(...)` is a new local helper that constructs
  `(D_groups, E_groups)` from `(theta_groups, phi_groups, N_j_groups)`.
- `angles_from_amplitudes(...)` is a new local helper that applies the inverse
  relations in the MFE definitions section.

Function flow: `solve_mfe(...)` is the main entry point. It calls
`amplitudes_from_initial_state(...)` to build the initial solver vector,
`solve_ivp(...)` to integrate `mfe_rhs(...)`, and
`compute_mfe_observables(...)` to attach angle and atom-number outputs.
`mfe_rhs(...)` calls `phase_values_at_time(...)` to evaluate the phase-local
equations of motion.

If a top-level `MomentSeries` container is already in use,
`attach_mfe_observables(moments, result)` should store the solved observable
series as `moments.mfe`.

= Data Requirements

The solver needs:

- `phases`: piecewise-constant `Phase` objects with `duration`, `omega`, and
  `delta`;
- `Gamma`: collective decay scale;
- `omega_groups`: one coupling weight per group;
- `N_j_groups`: active atom number per group for the mean-field state being
  solved;
- `initial_state`: initial J angles `(theta_groups, phi_groups)`;
- `t_eval`: saved output times.

Parameter validation may live outside the solver. In particular,
`omega_groups` should already follow the same weighted-coupling convention as
the MCWF code.

= Output

Important solver-facing data should be stored in Pydantic classes in
`parser/mfe.py`. A suggested minimal output structure is:

```python
MFESolverParameters(
    Gamma,
    phases,
    omega_groups,
    N_j_groups,
)

MFEInitialState(
    theta_groups,
    phi_groups,
)

MFEResult(
    t,
    D_groups,
    E_groups,
    success,
    message,
    parameters,
)

MFEObservableSeries(
    t,
    D_groups,
    E_groups,
    N_j_groups,
    theta_groups,
    phi_groups,
)
```

The observable series should be the object consumed by plotting or residual
diagnostics. Do not require plotting functions to unpack raw `solve_ivp`
objects.

= Invariants

- The standalone MFE solver should not import from `quantum_trajectories`.
- Use `parser.common.Phase` for phase metadata unless a stronger reason exists
  to define a solver-specific phase class.
- The solver should support arbitrary group count $G$ when the equations are
  written as sums over groups. Two-group-only logic belongs in residual
  diagnostics unless the theory explicitly requires two groups.
- Do not recompute coupling weights from instantaneous $N_(J,a)(t_(k))$. Coupling
  weights are fixed simulation inputs.
- The per-group norm $abs(D_(a))^2 + abs(E_(a))^2$ should remain constant up to
  solver tolerance. Large drift should be treated as a numerical warning.
- Phase and angle conventions should match `docs/instructions/j_moments.typ`
  and `docs/instructions/mfe_residuals.typ`.
