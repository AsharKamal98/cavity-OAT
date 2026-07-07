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
5. Run the chosen simulation backend.
6. Convert simulation output into moment series.
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
validate `dt`, `t_eval`, `ntraj`, `n_processes`, and homogeneous versus
inhomogeneous sector metadata.

Validation helpers should fail early with clear messages. Notebook workflows may
stop immediately after a failed validation rather than entering a long loop that
later crashes.

Detailed validation conventions should live in:

- `docs/instructions/parameter_validation.md` (planned/missing)

## 5. Simulation Backends

The repository currently supports multiple simulation backends that should feed
into the same post-processing pipeline.

High-level method flow should look like:

```python
mcwf_ensemble = run_trajectory_ensemble(...)
mfe_result = solve_mfe(...)
qutip_result = simulate_fixed_nj_mc_trajectory(...) or simulate_fixed_nj_me_trajectory(...)
```

The backend choice changes how the raw physical evolution is computed, but the
analysis flow after that should be shared as much as possible.

### 5.1 MCWF

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
parameters = MCWFSolverParameters(...)

def run_trajectory_ensemble(
    parameters, *, t_eval, ntraj, seed=None, n_processes=None, ...
) -> TrajectoryEnsemble:
    precomputed = build_precomputed_trajectory_data(...)
    trajectories = [_simulate_single_trajectory(..., precomputed) for _ in range(ntraj)]
    return TrajectoryEnsemble(trajectories=trajectories, metadata=metadata)
```

Shared simulation-level fields such as `Ni`, `omega_i`, `Gamma`, `phases`,
`shifted_jump_operator`, `t_eval`, and sector metadata should live in
`TrajectoryEnsemble.metadata` rather than on each `TrajectoryResult`.

For precomputation rules, use:

- `docs/instructions/simulation_precompute.typ`
- `docs/instructions/sector_operators.typ`

For ensemble-level simulation flow, use:

- `docs/instructions/ensemble_simulation.typ` (needs cleaning!)

For single-trajectory simulation flow, use:

- `docs/instructions/single_trajectory_simulation.typ` (needs cleaning!)

Each `TrajectorySnapshot` should contain the saved time, phase index, and sector
blocks on the internal `n_e` or `(n_e1, n_e2)` basis. Saved snapshots should be
aligned across trajectories so trajectory-level moment samples can be averaged
without interpolation.

### 5.2 MFE

The deterministic MFE backend should solve the group-resolved mean-field
equations and return an `MFEResult` on the requested saved-time grid:

```python
mfe_result = solve_mfe(parameters, t_eval=t_eval)
```

The MFE solver structure and the conversion from solved amplitudes to
`JMomentSeries` are defined in:

- `docs/instructions/mfe-solver.typ`

### 5.3 QuTiP

The QuTiP backends should run either `mcsolve` or `mesolve` and return the raw
QuTiP result together with the metadata needed for later J-moment extraction.

The preferred pattern is still that QuTiP outputs are converted into the same
shared `JMomentSeries` representation as the other methods.

## 6. Moment Series

The primary post-processing pipeline should be moment-first. Outputs from MCWF,
MFE, and QuTiP should all be converted into shared moment containers before
diagnostics or plotting.

For parser output-container conventions, use:

- `docs/instructions/generic/parser.md`

### 6.1 Moment Container

The notebook-level container is `MomentSeries`, defined in
`parser/moments.py`. It should be initialized on the shared
`t_eval` grid, from `total_time` and `num_snapshots`, and then filled in as
post-processing steps are run:

```python
raw_result = run_backend(...)
moments = MomentSeries(
    total_time=raw_result.metadata.t_eval[-1],
    num_snapshots=num_snapshots,
    parameters=MomentParameters(
        Gamma=raw_result.metadata.Gamma,
        phases=raw_result.metadata.phases,
        omega_groups=raw_result.metadata.omega_i,
        N_groups=raw_result.metadata.Ni,
    ),
)
moments.J = compute_method_j_moments(raw_result)
```

Current top-level fields are:

- `moments.t`: the shared saved-time grid.
- `moments.parameters`: simulation metadata needed by moment-level diagnostics,
  such as `Gamma`, `phases`, `omega_groups`, and `N_groups`; this can mirror
  selected fields from `ensemble.metadata`.
- `moments.J`: a `JMomentSeries` containing first-order J-sphere moments plus
  derived J-vector direction fields and angles when produced by
  `compute_mcwf_j_moments(...)`, or group-resolved MFE moment fields when
  produced by `compute_mfe_j_moments(...)`.
- `moments.mfe_residuals`: an `MFEResidualSeries` containing two-group MFE
  residual diagnostics when computed from `moments.J`.
- `moments.S`: placeholder for future S-moment or spin-direction data.

Parser container conventions live in `docs/instructions/generic/parser.md`.

### 6.2 J-Sphere Moments

Use the method-specific J-moment converters as the main post-processing entry
points:

```python
moments.J = compute_mcwf_j_moments(ensemble, n_processes=n_processes)
moments.J = compute_mfe_j_moments(mfe_result)
moments.J = compute_qutip_j_moments(qutip_result)
```

The returned `JMomentSeries` should be the shared spin-series representation
used by analysis and plotting. Depending on the method, it may contain full and
group-resolved `x`, `y`, `z`, `N_e`, `N_j`, `jump_rate`, `length`, `nx`, `ny`,
`nz`, `theta`, and `phi` fields when those quantities are available or
attached.

Detailed definitions and averaging rules live in
`docs/instructions/j_moments.typ`.

Numerical MFE solving and MFE-to-J-moment conversion are defined in:

- `docs/instructions/mfe-solver.typ`

QuTiP-to-J-moment conversion should follow the same shared output contract,
even if the internal extraction logic differs from the MCWF path.

When saving Slurm or batch-run J moments, store a self-contained artifact with
both the averaged `JMomentSeries` and the corresponding `phases`. Plotting
notebooks should read the saved `phases` from that artifact instead of
reconstructing the protocol manually.

### 6.3 J-Vector Direction Fields

After the raw first-order moments are available, derived direction fields
include `length`, `nx`, `ny`, and `nz`, plus group-resolved counterparts when
present. The current J-moment pipeline uses the Euclidean direction of the
stored J vector:

```python
moments.J.x
moments.J.y
moments.J.z
```

rather than the older active-manifold normalization by `N_active`. For MCWF,
the derived direction fields are attached after `compute_average_j_moments(...)`
returns the raw ensemble average. For other methods, the same fields should be
attached from the method's stored or reconstructed J components before
returning the final `JMomentSeries`.

Legacy note: these fields were previously named `Jx`, `Jy`, `Jz`, `Jx_groups`,
`Jy_groups`, `Jz_groups`, `J_len`, and `sx`, `sy`, `sz`.

## 7. Diagnostics

Diagnostics should consume already-computed moment series unless they genuinely
need lower-level saved simulation data. Reusable post-simulation physics
diagnostics should preferably live in the root-level `post_analysis/` package
rather than in `common/`.

Current task-specific diagnostic instructions include:

- `docs/instructions/mfe_residuals.typ`
- `docs/instructions/squeezing.typ`
- `docs/instructions/dephasing_diagnostics.typ`
- `docs/instructions/post_analysis.md`

Future diagnostics should get their own instruction files when they introduce
new averaging rules, new physical conventions, or nontrivial plotting logic.

### 7.1 MFE Residuals

Use `compute_mfe_residuals(...)` from `post_analysis/mfe_residuals.py`
after `moments.J` has been computed:

```python
moments.mfe_residuals = compute_mfe_residuals(
    moments.J,
    parameters=moments.parameters,
)
```

Detailed residual definitions live in `docs/instructions/mfe_residuals.typ`.

## 8. Plotting and Notebook Workflows

Plotting code should be thin: it should visualize already-computed moment
series or diagnostics rather than hiding moment extraction or other heavy
calculations inside plotting calls.

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

### 8.1 Shared Plotting

The shared spin-component plot now lives in `common/plotting/j_spin.py`:

- `plot_spin_components(series, ...)`: plots stored `x`, `y`, `z`, `length`,
  and the matching group-resolved fields when present.

The shared angle plot now lives in `common/plotting/j_spin.py`:

- `plot_bloch_angles(series, ...)`: plots whatever stored `theta`, `phi`,
  `theta_groups`, and `phi_groups` fields are available on the input series,
  using the selected `colour_index` palette and `linestyle`.

Current diagnostic plotting functions live in
`common/plotting/mfe_residuals.py`:

- `plot_mfe_residuals(moments.mfe_residuals, ...)`: plots stored two-group
  residuals in a single residual panel with the L2 norm.

General diagnostic plotting functions currently live in
`legacy/plotting_diagnostics.py`:

- `plot_sector_probabilities(result, ...)`: plots normalized represented-sector
  probabilities `p_alpha(t)` computed directly from saved snapshot sector
  blocks.

Detailed plotting conventions live in `docs/instructions/plotting_workflows.md`.
Future diagnostics and plots should consume `MomentSeries` or `JMomentSeries`
when the required data are already present, instead of rerunning older
observable extraction steps.

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
