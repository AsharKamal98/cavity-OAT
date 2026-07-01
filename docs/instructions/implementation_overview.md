# Implementation Overview

This is the main entry point for implementation-specific guidance. Read it
before changing code.

Instruction files describe how the repository implements the physics from
`docs/theory`. They should focus on intended logic, data flow,
performance-relevant choices, and invariants that must survive refactors. They
should not duplicate detailed theory or document every line of code.

Generic reusable instruction files live in `docs/instructions/generic/` and can
be copied between projects.

Shared repository utility helpers are summarized in
`docs/instructions/common_utils.typ`.

When implementing a feature:

- Read the relevant theory file for equations and physical meaning.
- Read this overview for the expected code structure.
- Read the task-specific instruction file for conventions and edge cases.
- If no task-specific instruction file exists, implement the minimal consistent
  behavior and add or update an instruction file when the convention should be
  preserved.

## High-Level Workflow

Most custom-code runs should follow this order:

1. Define physical parameters.
2. Define the phase protocol.
3. Construct the initial sector expansion.
4. Validate parameter regimes.
5. Run the custom MCWF ensemble.
6. Compute trajectory or ensemble observables.
7. Build derived diagnostics.
8. Plot or export results.

Notebook functions should follow this same order when possible, so the main
analysis cells stay predictable.

## 1. Physical Parameters

Runs should define the effective parameters `Omega`, `Gamma`, and `delta`
explicitly, or derive them from cavity parameters when using the cavity-model
helpers.

Detailed implementation conventions should live in:

- `docs/instructions/common_utils.typ`
- `docs/instructions/simulation_parameters.typ` (needs cleaning!)

## 2. Phase Protocol

Protocol phases should be built through the shared phase-construction helpers,
such as `default_three_phase_protocol(...)`, rather than duplicated manually in
many places:

```python
phases = default_three_phase_protocol(
    T1=T1,
    T2=T2,
    T3=T3,
    delta0=delta0,
    Omega0=Omega0,
)
```

Each phase should carry duration, `Omega`, and `delta`. Phases are piecewise
constant, which enables phase-level operator and propagator precomputation.
Shared phase helpers are summarized in `docs/instructions/common_utils.typ`.

## 3. Initial Sector State

The initial wavefunction should first be defined over strong-symmetry sectors.
The simulation should then evolve one internal array per sector in the
active-manifold excitation basis.

Homogeneous runs use scalar sector keys:

```python
{Nj: coeff}
```

with selected sectors

```python
Nj in range(N//2 - dN, N//2 + dN + 1)
```

Inhomogeneous runs use group-resolved sector keys:

```python
{(Nj1, Nj2): coeff}
```

where

```python
Nj1 + Nj2 = Nj
0 <= Nj1 <= N1
0 <= Nj2 <= N2
```

High-level initialization helpers should be preferred over manually writing
low-level sector dictionaries. Supported sector coefficient choices should be:

```python
sector_distribution="square"    # equal total weight over selected Nj sectors
sector_distribution="binomial"  # product-state binomial weights, truncated to the selected window
```

After the sector coefficients are chosen, the propagated state should be stored
as sector blocks:

```python
{Nj: psi_Nj}              # psi_Nj has shape (Nj + 1,) over n_e = 0,...,Nj
{(Nj1, Nj2): psi_Nj1_Nj2} # shape ((Nj1 + 1) * (Nj2 + 1),) over (n_e1, n_e2)
```

By default, each internal sector state should start in the all-active-atoms-down
state, i.e. `n_e=0` or `(n_e1, n_e2)=(0, 0)`.

Detailed initialization conventions should live in:

- `docs/instructions/initial_sector_state.md` (planned/missing)

For inhomogeneous coupling conventions, use:

- `docs/instructions/paper_inhomogeneous_couplings.typ` (needs cleaning!)

## 4. Parameter Validation

Validation should happen before long simulations begin. Existing validation
checks include:

```python
Omega_Gamma_from_cavity_parameters(...)   # validates bad-cavity limit
validated_mcwf_dt(...)                    # enforces the MCWF dt rule
check_initial_sector_omega_ratio(...)     # checks Omega / Omega_c for initial sectors
validate_sector_distribution(...)         # checks "square" / "binomial"
omega2_from_weighted_average(...)         # validates inhomogeneous group sizes
```

Initial-state helpers should validate sector ranges, even `N` assumptions,
group sizes, and internal sector-state shapes. Simulation entry points should
validate `dt`, `num_snapshots`, `ntraj`, `n_processes`, and homogeneous versus
inhomogeneous sector metadata.

Validation helpers should fail early with clear messages. Notebook workflows may
stop immediately after a failed validation rather than entering a long loop that
later crashes.

Detailed validation conventions should live in:

- `docs/instructions/parameter_validation.md` (planned/missing)

## 5. Ensemble Simulation

The main custom MCWF entry point is:

```python
run_trajectory_ensemble(...)
```

The ensemble function should:

- create trajectory seeds consistently;
- build reusable precomputed data once;
- run single trajectories serially or with multiprocessing;
- return a `TrajectoryEnsemble` containing the individual trajectory results
  and shared simulation parameters.

High-level data flow should look like:

```python
def run_trajectory_ensemble(
    N, Gamma, phases, sector_coeffs, dt, num_snapshots, ntraj, ...
) -> TrajectoryEnsemble:
    precomputed = build_precomputed_trajectory_data(...)
    trajectories = [simulate_single_trajectory(..., precomputed) for _ in range(ntraj)]
    return TrajectoryEnsemble(trajectories=trajectories, parameters=parameters)
```

`parameters` should contain shared metadata from the ensemble inputs, including
computed `omega_groups` and `N_groups`; `TrajectoryEnsemble` should validate it
into `MomentParameters`.

For precomputation rules, use:

- `docs/instructions/simulation_precompute.typ`
- `docs/instructions/sector_operators.typ`

These files should cover reusable sector operators, phase generators,
full-`dt` propagators, and when precomputed data can or cannot be used.

They should be read before changing build precomputed trajectory data(...),
phase-dependent jump operators, non-Hermitian generators, full-dt propagators, 
or the logic that chooses between precomputed and variable-step propagation.


For ensemble-level simulation flow, use:

- `docs/instructions/ensemble_simulation.typ` (needs cleaning!)

This file should cover seed construction, multiprocessing, worker state, shared
precomputed data, and collecting `TrajectoryResult`s into a
`TrajectoryEnsemble`.

### Single Trajectories

Each trajectory should evolve the same physical model as the ensemble run and
save snapshots on the common `t_eval` grid implied by `num_snapshots`.

Trajectory evolution may use adaptive internal logic such as jump refinement or
partial steps, but saved snapshots should be aligned across trajectories so
ensemble observables can be averaged without interpolation.

High-level data flow should look like:

```python
simulate_single_trajectory(
    N, Gamma, phases, sector_coeffs, dt, num_snapshots, seed, precomputed, ...
) -> TrajectoryResult(
    snapshots=[TrajectorySnapshot, ...],
    final_sector_blocks={sector_key: psi_sector},
    jump_times=[...],
)
```

Each `TrajectorySnapshot` should contain the saved time, phase index, and sector
blocks on the internal `n_e` or `(n_e1, n_e2)` basis.

For single-trajectory simulation flow, use:

- `docs/instructions/single_trajectory_simulation.typ` (needs cleaning!)

This file should cover `simulate_single_trajectory(...)`, `t_eval` saving,
full-`dt` steps versus partial steps, jump detection, jump bisection, and
snapshot/result construction.

## 6. Observables

Observable code should convert saved trajectory snapshots into expectation-value
time series. For ensembles, average the underlying observables or moments before
constructing nonlinear derived quantities when the target object is the
unconditioned state.

For parser output-container conventions, use:

- `docs/instructions/generic/parser.md`

High-level data flow should look like:

```python
trajectory_observables(TrajectoryResult) -> ObservableSeries
ensemble_observables(TrajectoryEnsemble) -> ObservableSeries
```

For reusable first-order J-moment extraction, use:

- `docs/instructions/j_moments.typ`

This file should cover `compute_trajectory_j_moments(...)`,
`compute_average_j_moments(...)`, `compute_ensemble_j_moments(...)`,
`JMomentSnapshot`, `JMomentSeries`, and the meaning of each saved J-moment
field.

Use task-specific rules for nonlinear diagnostics:

- Numerical MFE solving:
  `docs/instructions/mfe-solver.typ`
- MFE residuals:
  `docs/instructions/mfe_residuals.typ`
- Bloch angles and active-manifold directions:
  `docs/instructions/bloch_vector_averaging.typ`
- Normalized active-manifold spin-component plots:
  `docs/instructions/plot_spin_components.typ`
- Generalized squeezing:
  `docs/instructions/squeezing.typ`
- Dephasing Bloch-vector lengths:
  `docs/instructions/dephasing_diagnostics.typ` (needs cleaning!)

Extensive observables such as total atom numbers, jump rates, and jump counts
should not be normalized using Bloch-direction conventions unless a diagnostic
explicitly asks for that.

## 7. Diagnostics

Diagnostics should usually be standalone post-processing functions. They should
read `TrajectoryResult`, `TrajectoryEnsemble`, or already-computed observables
without changing the MCWF propagation path unless new saved data are genuinely
required.

High-level data flow should look like:

```python
diagnostic(result_or_observables, ...) -> diagnostic_data
plot_diagnostic(diagnostic_data or result, ...) -> matplotlib figure/axes
```

Current task-specific diagnostic instructions include:

- `docs/instructions/squeezing.typ`
- `docs/instructions/dephasing_diagnostics.typ`
- `docs/instructions/bloch_vector_averaging.typ`
- `docs/instructions/plot_spin_components.typ`

Future diagnostics should get their own instruction files when they introduce
new averaging rules, new physical conventions, or nontrivial plotting logic.

## 8. Plotting and Notebook Workflows

Plotting code should be thin: it should visualize already-computed observables
or diagnostics rather than hiding heavy physics calculations inside plotting
calls.

Notebook functions should:

- define their own inputs instead of relying on hidden globals;
- construct phases locally;
- validate parameters before simulations;
- expose runtime-heavy options explicitly;
- reuse shared helpers instead of duplicating logic.

Plot labels should reflect the physical quantity being plotted, especially when
comparing simulation, theory, and approximations.

Detailed plotting conventions live in:

- `docs/instructions/plotting_workflows.md`

## 9. New pipeline

Sections 6-8 describe the established observable, diagnostic, and plotting
pipeline. This section documents the newer moment-first pipeline that is being
built to replace parts of the old observable flow. Use this section when the
task explicitly mentions the new pipeline, `MomentSeries`, J moments, or
`quantum_trajectories/plotting_j_moments.py`.

### 9.1 Moment Container

The notebook-level container is `MomentSeries`, defined in
`parser/moments.py`. It should be initialized on the shared
`t_eval` grid, from `phases` and `num_snapshots`, and then filled in as
post-processing steps are run:

```python
ensemble = run_trajectory_ensemble(...)
moments = MomentSeries(
    phases=phases,
    num_snapshots=num_snapshots,
    parameters=ensemble.parameters,
)
moments.J = compute_ensemble_j_moments(ensemble)
```

Current top-level fields are:

- `moments.t`: the shared saved-time grid.
- `moments.parameters`: simulation metadata needed by moment-level diagnostics,
  such as `Gamma`, `phases`, `omega_groups`, and `N_groups`; this can mirror
  `ensemble.parameters`.
- `moments.J`: a `JMomentSeries` containing first-order J-sphere moments plus
  derived J-vector direction fields and angles when produced by
  `compute_ensemble_j_moments(...)`.
- `moments.mfe`: an `MFEObservableSeries` containing the solved MFE observable
  series when attached from the standalone MFE pipeline.
- `moments.mfe_residuals`: an `MFEResidualSeries` containing two-group MFE
  residual diagnostics when computed from `moments.J`.
- `moments.S`: placeholder for future S-moment or spin-direction data.

Parser container conventions live in `docs/instructions/generic/parser.md`.

### 9.2 J-Sphere Moments

Use `compute_ensemble_j_moments(...)` from `quantum_trajectories/j_moments.py`
as the main new-pipeline entry point for trajectory-averaged J moments:

```python
moments.J = compute_ensemble_j_moments(ensemble, n_processes=n_processes)
```

The returned `JMomentSeries` contains arrays on the saved `t_eval` grid,
including `x`, `y`, `z`, `N_e`, `N_j`, `jump_rate`, `J_drive`, and optional
group-resolved fields such as `x_groups`, `y_groups`, `z_groups`,
`N_e_groups`, and `N_j_groups`. It also includes derived fields attached after
ensemble averaging: `length`, `nx`, `ny`, `nz`, `theta`, and `phi`, plus
group-resolved versions when group fields exist.

Detailed definitions and averaging rules live in
`docs/instructions/j_moments.typ`.

### 9.3 J-Vector Direction Fields

After the raw first-order moments are averaged across trajectories, derived
direction fields include `length`, `nx`, `ny`, and `nz`, plus group-resolved
counterparts when present. The current J-moment pipeline uses the Euclidean
direction of the averaged J vector:

```python
moments.J.x
moments.J.y
moments.J.z
```

rather than the older active-manifold normalization by `N_active`. The derived
direction fields are attached inside `compute_ensemble_j_moments(...)` after
`compute_average_j_moments(...)` returns the raw ensemble average. Angle fields
such as `theta`, `phi`, `theta_groups`, and `phi_groups` are then attached from
those normalized directions before returning the final `JMomentSeries`.

Legacy note: these fields were previously named `Jx`, `Jy`, `Jz`, `Jx_groups`,
`Jy_groups`, `Jz_groups`, `J_len`, and `sx`, `sy`, `sz`.

### 9.4 MFE Residuals

Use `compute_mfe_residuals(...)` from `quantum_trajectories/mfe_residuals.py`
after `moments.J` has been computed:

```python
moments.mfe_residuals = compute_mfe_residuals(
    moments.J,
    parameters=moments.parameters,
)
```

Detailed residual definitions live in `docs/instructions/mfe_residuals.typ`.

### 9.5 New-Pipeline Plotting

New-pipeline plotting functions should take moment series objects directly,
usually `moments.J`, and should visualize already-computed fields rather than
recomputing moments.

Current moment plotting functions live in `quantum_trajectories/plotting_j_moments.py`:

- `plot_j_spin_components(moments.J, ...)`: plots `x`, `y`, and `z`, plus
  group-resolved curves when present.
- `plot_j_angles(moments.J, ...)`: plots stored `theta` and `phi`, plus stored
  group-resolved angle curves when present.

Current diagnostic plotting functions live in
`quantum_trajectories/plotting_mfe_residuals.py`:

- `plot_mfe_residuals(moments.mfe_residuals, ...)`: plots stored two-group
  residuals in a single residual panel with the L2 norm.

General diagnostic plotting functions live in
`quantum_trajectories/plotting_diagnostics.py`:

- `plot_sector_probabilities(result, ...)`: plots normalized represented-sector
  probabilities `p_alpha(t)` computed directly from saved snapshot sector
  blocks.

Detailed plotting conventions live in `docs/instructions/plotting_workflows.md`.
Future new-pipeline diagnostics should consume `MomentSeries` or
`JMomentSeries` when the required data are already present, instead of rerunning
old observable extraction.

## General Rules

- This implementation overview should mainly point to task-specific
  instruction files, rather than directly to theory files.
- Theory files may be referenced directly when this overview gives an
  implementation instruction that depends on a specific theory source.
- During development, this overview may refer to task-specific instruction
  files that do not yet exist, as long as the reference clearly says that the
  file is planned or missing.
- Instruction files should not include every coding detail. They should
  document the logic, conventions, and performance-relevant choices that matter
  for correctly implementing the physics.
- When a code change establishes or changes intended behavior, update the relevant
  instruction file. If no relevant file exists, either add one or add a concise
  placeholder reference here.
- Instruction files should be concise and to the point. Avoid text that is not
  necessary for implementation. When theory context is needed, point to the
  relevant theory file instead of repeating the theory.
- If an instruction can be explained clearly with pseudo-code or a short code
  sketch, prefer that over a lengthy prose explanation.
- Instruction files also serve as documentation. When it clarifies a function
  or feature, include:
  1. the physics expected from that functionality;
  2. pseudo-code-style data input shape;
  3. pseudo-code-style data output shape;
  4. high-level pseudo-code for long or complex functions.
- Do not automatically modify instruction files or parts of instruction files
  unrelated to the given task, unless it is a smaller typo. If a larger
  inconsistency is found in unrelated instructions to the task, let the user know.

## TODOs

Open implementation cleanup notes are tracked in `docs/instructions/TODOs.md`.
