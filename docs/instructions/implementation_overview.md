# Implementation Overview

This is the main entry point for implementation-specific guidance. Read it
before changing code.

Instruction files describe how the repository implements the physics from
`docs/theory`. They should focus on intended logic, data flow,
performance-relevant choices, and invariants that must survive refactors. They
should not duplicate detailed theory or document every line of code.

Shared repository utility helpers are summarized in

- `docs/instructions/common/utils.typ` for shared utility helpers.
- `docs/instructions/common/plotting.typ` for shared plotting helpers and
  common plotting-function contracts.


## High-Level Workflow

Most runs follow this order:

1. Define model and method parameters.
2. Define the phase protocol.
3. Run the chosen simulation backend.
4. Convert simulation output into moment series.
5. Build derived diagnostics.
6. Plot or export results.

Notebook functions should follow this same order when possible, so the main
analysis cells stay predictable.

## 1. Model and Method Parameters

Runs should define the effective physical/model parameters such as `Omega`,
`Gamma`, `delta`, coupling choices, and the phase protocol.

These may be defined explicitly or through shared helpers in
`common/utils/parameters.py`.

Detailed implementation conventions live in:
  `docs/instructions/model_parameters.typ`

Numerical method parameters such as timesteps, saved-time grids, and solver tolerances live in:
  `docs/instructions/solvers/mcwf/method_parameters.typ`

## 2. Phase Protocol

For the standard protocol, construct the phase protocol first and provide it to
`SimulationMetadata`:

```python
phase_protocol = default_three_phase_protocol(
    durations=(T1, T2, T3),
    delta0=delta0,
    Omega0=Omega0,
)
metadata = SimulationMetadata(
    Ni=Ni,
    omega_i=omega_i,
    Gamma=Gamma,
    phase_protocol=phase_protocol,
)
```

Each `FamilyPhase` carries a duration and target `Omega` and `delta`. Its
integration `Phase` objects are piecewise constant, which enables operator and
propagator precomputation. Optional ramp durations and ramp-segment counts
control a piecewise-constant ramp followed by a hold at the target values.

The resulting standard phase protocol is:

```text
Phase 1: T = T1, Omega = Omega0, delta = 0
Phase 2: T = T2, Omega = Omega0, delta = delta0
Phase 3: T = T3, Omega = 0,      delta = 0
```

Detailed phase-protocol conventions live in
`docs/instructions/phases.typ`.

## 3. Simulation Backends

The repository currently supports multiple simulation backends that should feed
into the same post-processing pipeline.

High-level method flow should look like:

```python
# define shared model/protocol data and result container
metadata = SimulationMetadata(...)
moments = MomentSeries(num_snapshots=100, metadata=metadata, ...)

# define method-specific parameters from shared metadata
parameters = <Method>SolverParameters(
    Ni=metadata.Ni,
    omega_i=metadata.omega_groups,
    Gamma=metadata.Gamma,
    phase_protocol=metadata.phase_protocol,
    ...,
)

# run method
mcwf_ensemble = run_trajectory_ensemble(parameters, ...)
mfe_result = solve_mfe(parameters, ...)
qutip_result = simulate_fixed_nj_mc_trajectory(parameters, ...) or simulate_fixed_nj_me_trajectory(parameters, ...)
```

The backend choice changes how the raw physical evolution is computed, but the
analysis flow after that should be shared as much as possible.

### 3.1 MCWF

The main custom MCWF entry point is:

```python
run_trajectory_ensemble(...)
```

The ensemble function should:

- create trajectory seeds consistently;
- build reusable precomputed data once;
- run single trajectories serially or with multiprocessing;
- return a `TrajectoryEnsemble` containing the individual trajectory results
  and shared simulation metadata.

High-level data flow should look like:

```python
parameters = MCWFSolverParameters(
    Ni=metadata.Ni,
    omega_i=metadata.omega_groups,
    Gamma=metadata.Gamma,
    phase_protocol=metadata.phase_protocol,
    dN=dN,
    dt=dt,
    ...,
)

def run_trajectory_ensemble(
    parameters, *, t_eval, ntraj, seed=None, n_processes=None, ...
) -> TrajectoryEnsemble:
    sector_coeffs = centered_sector_initial_coeffs(...)
    precomputed = build_precomputed_trajectory_data(...)
    trajectories = [_simulate_single_trajectory(..., precomputed) for _ in range(ntraj)]
    return TrajectoryEnsemble(trajectories=trajectories, metadata=metadata)
```

Shared simulation-level fields such as `Ni`, `omega_i`, `Gamma`,
`phase_protocol`, `shifted_jump_operator`, `t_eval`, and sector metadata should
live in `TrajectoryEnsemble.metadata`.

Detailed initialization conventions live in:

- `docs/instructions/solvers/mcwf/initial_sector_state.typ`

For precomputation rules, use:

- `docs/instructions/solvers/mcwf/simulation_precompute.typ`
- `docs/instructions/solvers/mcwf/sector_operators.typ`

For ensemble-level simulation flow, use:

- `docs/instructions/solvers/mcwf/ensemble_simulation.typ` 

For single-trajectory simulation flow, use:

- `docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`


### 3.2 MFE

The deterministic MFE backend should solve the group-resolved mean-field
equations and return an `MFEResult` on the requested saved-time grid:

```python
mfe_result = solve_mfe(parameters, t_eval=t_eval)
```

The MFE solver structure and the conversion from solved amplitudes to
`JMomentSeries` are defined in:

- `docs/instructions/solvers/mfe/mfe-solver.typ`

### 3.3 QuTiP

The QuTiP backends should run either `mcsolve` or `mesolve` and return the raw
QuTiP result together with the metadata needed for later J-moment extraction.

The preferred pattern is still that QuTiP outputs are converted into the same
shared `JMomentSeries` representation as the other methods.

## 4. Moment Series

The primary post-processing pipeline should be moment-first. Outputs from MCWF,
MFE, and QuTiP should all be converted into shared moment containers before
diagnostics or plotting.

For parser output-container conventions, use:

- `docs/instructions/parser.typ`
- the `parser` skill

### 4.1 Moment Container

The notebook-level container is `MomentSeries`, defined in
`parser/moments.py`. It should be initialized on the shared
`t_eval` grid, from the metadata phase durations and `num_snapshots`, and then filled in as
post-processing steps are run:

```python
raw_result = run_backend(...)
moments = MomentSeries(
    num_snapshots=num_snapshots,
    metadata=SimulationMetadata(
        Ni=raw_result.metadata.Ni,
        omega_i=raw_result.metadata.omega_i[:-1],
        Gamma=raw_result.metadata.Gamma,
        phase_protocol=raw_result.metadata.phase_protocol,
    ),
)
moments.J = compute_method_j_moments(raw_result)
```

Current top-level fields are:

- `moments.t`: the shared saved-time grid.
- `moments.metadata`: validated physical model and standard-protocol context,
  including `Ni`, independent `omega_i`, completed `omega_groups`, `Gamma`,
  and the supplied `phase_protocol`.
- `moments.J`: a `JMomentSeries` containing first-order J-sphere moments plus
  derived J-vector direction fields and angles.
- `moments.mfe_residuals`: an `MFEResidualSeries` containing two-group MFE
  residual diagnostics when computed from `moments.J`.
- `moments.S`: placeholder for future S-moment or spin-direction data.

Parser container conventions live in `docs/instructions/parser.typ` and the
`parser` skill.

### 4.2 J-Sphere Moments

J-sphere moments are the common post-processing representation used across the
simulation backends. The solvers should stay responsible for producing raw
dynamics, while method-specific converters turn those outputs into a shared
`JMomentSeries`.

Use the method-specific converters as the main entry points:

```python
moments.J = compute_mcwf_j_moments(ensemble, n_processes=n_processes)
moments.J = compute_mfe_j_moments(mfe_result)
moments.J = compute_qutip_j_moments(qutip_result)
```

The shared series may contain spin components, normalized directions, angles,
atom-number fields, jump rates, and group-resolved counterparts when those
quantities are available. Derived fields such as spin directions and angles
should be attached inside the J-moment layer, so plotting and diagnostics can
consume one common structure independent of whether the data came from MCWF,
MFE, or QuTiP. Additive full-system fields may be constructed by summing their
group-resolved counterparts before nonlinear derived fields are attached.

Detailed definitions and averaging rules live in
`docs/instructions/j_moments.typ`.

Numerical MFE solving and MFE-to-J-moment conversion are defined in:

- `docs/instructions/solvers/mfe/mfe-solver.typ`

QuTiP-to-J-moment conversion is defined in  

- `docs/instructions/solvers/qutip/....typ` (planned)


## 5. Diagnostics

Diagnostics should consume already-computed moment series unless they genuinely
need lower-level saved simulation data. Reusable post-simulation physics
diagnostics should live in the root-level `post_analysis/` package.

Current task-specific diagnostic instructions include:

- `docs/instructions/post_analysis/mfe_residuals.typ`
- `docs/instructions/post_analysis/squeezing.typ`

summarized in
  `docs/instructions/post_analysis/post_analysis.md`

## 6. Plotting
The shared spin-component plot lives in `common/plotting/j_spin.py`:
- `plot_spin_components(series, ...)`: plots stored `x`, `y`, `z`, `length`,
  and the matching group-resolved fields when present.

The shared angle plot lives in `common/plotting/j_spin.py`:
- `plot_bloch_angles(series, ...)`: plots whatever stored `theta`, `phi`,
  `theta_groups`, and `phi_groups` fields are available on the input series,
  using the selected `colour_family_index` / `shade_index` palette and
  `linestyle`.

The shared MFE residual diagnostic plotting function lives in
`common/plotting/mfe_residuals.py`:
- `plot_mfe_residuals(moments.mfe_residuals, ...)`: plots stored two-group
  residuals in a single residual panel with the L2 norm.

Repo-specific plotting-function contracts live in
`docs/instructions/plotting_workflows.md`.

Future diagnostics and plots should consume `MomentSeries` or `JMomentSeries`
when the required data are already present, instead of rerunning older observable extraction steps.

## TODOs

Open implementation cleanup notes are tracked in `docs/instructions/TODOs.md`.
