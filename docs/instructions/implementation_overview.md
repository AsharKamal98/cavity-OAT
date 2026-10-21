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
4. Attach the simulator's stored moments to the shared moment series.
5. Build derived diagnostics.
6. Plot or export results.

Notebook functions should follow this same order when possible, so the main
analysis cells stay predictable.

In comparison notebooks, write each method's setup, analysis, and plotting
calls explicitly, with MCWF followed by MFE. Keep analysis and plotting cells
free of loops and comprehensions. Pass multiple curves to existing plotters
and move repeated diagnostic work into shared helpers.

## 1. Model and Method Parameters

Runs should define the effective physical/model parameters such as `Omega`,
`Gamma`, `delta`, coupling choices, and the phase protocol.

These may be defined explicitly or through shared helpers in
`common/utils/parameters.py`.

Detailed implementation conventions live in:
  `docs/instructions/model_parameters.typ`

Numerical method parameters such as timesteps, saved-time grids, and solver tolerances live in:
  `docs/instructions/solvers/mcwf/method_parameters.typ`

When asked to estimate two-group MCWF runtime as a function of `Ni` and `dN`,
use the SURFcua baseline and factors in:
  `docs/instructions/scaling.typ`
Treat `Ni`-dependent runtime estimates as provisional until the additional
SURFcua calibration described there is available.

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

# retain shared raw moment output before constructing derived observables
moments.generator_moments = mcwf_ensemble.moments

# form the central ensemble moments before J, S, and squeezing analysis
moments.mean_generator_moments = compute_mean_generator_moments(
    moments.generator_moments
)

# derive solver-independent J moments from the central generator moments
moments.J = compute_j_moments(
    moments.mean_generator_moments,
    metadata=moments.metadata,
)

# derive group-local S moments in the J-angle-defined frame
moments.S = compute_s_moments(
    moments.mean_generator_moments,
    moments.J,
)

# rotate both averaged generator-moment orders into the dressed frame
U = build_dressed_basis_transformation(moments.J, moments.S)
R = compute_gell_mann_rotation(U)
mean_dressed_generator_moments = rotate_generator_moments(
    moments.mean_generator_moments,
    R,
)

# compute group-resolved squeezing and, when requested, the full system
moments.squeezing = compute_squeezing(
    mean_dressed_generator_moments,
    moments.J,
    metadata=moments.metadata,
    include_full_system=False,
)
```

The backend choice changes how the raw physical evolution is computed, but the
analysis flow after that should be shared as much as possible. Resampling for
error bars should select complete rows from `generator_moments`, recompute the
mean generator moments, and then repeat the full derived analysis.

### 3.1 MCWF

The main custom MCWF entry point is:

```python
run_trajectory_ensemble(...)
```

Slurm jobs should use the unified command-line entry point
`slurm/run_solvers/run_mcwf.py`. Pass the physical groups directly through
`--Ni`: a one-group run omits `--omega-i`, while a two-group run supplies one
independent coupling through `--omega-i`. The runner must not silently collapse
multiple supplied groups into one. An explicit `--dt` is optional; when omitted,
the runner should use the standard `mcwf_dt_from_scales(...)` result.

Submit each analysis SLURM script from its own directory, for example with
`cd slurm/analyses/squeezing_omega1-scan && sbatch submit.slurm`. Runtime paths
should resolve from `SLURM_SUBMIT_DIR`; analysis artifacts should use that
analysis's local `output` directory, while scheduler logs remain in
`slurm/logs`.

For post-analysis parameter sweeps, use one SLURM array task per parameter
value instead of processing the full sweep sequentially within one task. Each
task should run the complete requested trajectory ensemble using its internal
multiprocessing and write one trajectory-resolved moment artifact for that
parameter value. This gives every parameter value its own wall-time allocation,
precomputation, and output file without a separate combination step.

An analysis-local notebook may load those parameter-resolved artifacts from
the analysis's `output` directory and run the shared mean-generator, J/S,
dressed-rotation, squeezing, and plotting pipeline locally. Keep each raw
trajectory axis intact until the notebook explicitly computes its ensemble
mean.

The bright/dark workflow lives in
`slurm/analyses/bright_dark_modes_mcwf`. Run `bash submit.sh` there to submit
four trajectory chunks and an `afterok`-dependent combination job. `NTRAJ`
in `submit.slurm` is the count per task, so the combined total is `4 * NTRAJ`.
Each task uses the unified runner with its own array index and seed, saving under
`output/<array-job-id>/`. The combination job validates and concatenates raw
moments using `slurm/trajectory_ensemble_io.py`, then saves `ensemble.pkl`. The local
`post_analysis.ipynb` computes J moments and modes, solves MFE with the saved
physical context, and plots J angles/moments and bright/dark diagnostics with
MCWF SEMs. Keep its simulation settings and saved-time resolution aligned with
`notebooks/bright_dark_modes_mcwf.ipynb`. Retain raw trajectory moments for the
matching windowed harmonics and early-time bootstrap analysis, including the
dashed MFE overlays. Follow `docs/instructions/post_analysis/j_modes.typ` for
this workflow.

When an ensemble is split across multiple SLURM tasks, each task must write a
separate chunk. A dependent combination step should concatenate first- and
second-moment arrays along the trajectory axis before
generic analysis. Multiple tasks must never append to or rewrite one shared
pickle artifact concurrently. Do not average per-task derived moments.

The ensemble function should:

- create trajectory seeds consistently;
- build reusable precomputed data once;
- run single trajectories serially or with multiprocessing;
- stack the raw trajectory moments without averaging and discard the
  single-trajectory diagnostics after any runtime summary;
- return a `TrajectoryEnsemble` containing one trajectory-resolved moment
  series and shared simulation metadata.

For a homogeneous grouping-equivalence check, use the same total atom number,
phase protocol, total-$N_J$ sector window, timestep, trajectory seeds, and unit
coupling in both the one-group and explicitly split two-group runs. Compare the
one-group observables with additive two-group J and S quantities and with
two-group full-system squeezing; a local subgroup is not the counterpart of
the complete one-group system.

High-level data flow should look like:

```python
parameters = MCWFSolverParameters(
    Ni=metadata.Ni,
    omega_i=metadata.omega_groups,
    Gamma=metadata.Gamma,
    phase_protocol=metadata.phase_protocol,
    dN=dN,
    max_omega_over_omega_c=None,
    linear_constraints=None,
    joint_probability=None,
    dt=dt,
    ...,
)

def run_trajectory_ensemble(
    parameters, *, t_eval, ntraj, seed=None, n_processes=None, ...
) -> TrajectoryEnsemble:
    sector_coeffs, sector_summary = group_generic_initial_sector_coeffs(...)
    print_initial_sector_summary(sector_summary)
    precomputed = build_precomputed_trajectory_data(...)
    trajectory_results = (
        _simulate_single_trajectory(..., precomputed) for _ in range(ntraj)
    )
    moments, total_steps, non_precomputed_steps = _collect_trajectory_moments(
        trajectory_results, ntraj=ntraj, t_eval=t_eval, num_groups=len(parameters.Ni)
    )
    return TrajectoryEnsemble(moments=moments, seeds=seed_keys, metadata=metadata)
```

Shared simulation-level fields such as `Ni`, `omega_i`, `Gamma`,
`phase_protocol`, `shifted_jump_operator`, `t_eval`, and sector metadata should
live in `TrajectoryEnsemble.metadata`.

The custom MCWF solver should use the shifted jump-operator picture by default.
Pass `shifted_jump_operator=False` explicitly when the unshifted picture is
required.

Detailed initial-sector selection and state-construction conventions live in:

- `docs/instructions/solvers/mcwf/initial_sector_selection.typ`

Apply initial-sector filters independently in a clear sequence: `dN` first,
`max_omega_over_omega_c` second, exact $N_J$-ratio equations in
`linear_constraints` third, and `joint_probability` last. Execute each filter
only when its input is non-`None`; all may be omitted to retain the complete
support. The drive-ratio filter removes complete total-$N_J$ shells below its
derived strict minimum. `joint_probability` caps
exact retained binomial probability relative to the complete unconstrained
wavefunction and always ranks individual group-resolved tuples, not complete
total-$N_J$ sectors.
The ensemble setup should print the initial tuple-sector count followed by one
aligned row per applied filter. Each row reports sectors remaining, retention
relative to that filter's input, and retention relative to the initial tuple
space, for both sector count and exact binomial probability.

For precomputation rules, use:

- `docs/instructions/solvers/mcwf/simulation_precompute.typ`
- `docs/instructions/solvers/mcwf/sector_operators.typ`

The precompute builder may construct aligned `SectorOperators` objects
transiently, but its returned worker state should retain only the matrices and
metadata consumed by trajectory propagation.

For ensemble-level simulation flow, use:

- `docs/instructions/solvers/mcwf/ensemble_simulation.typ` 

For single-trajectory simulation flow, use:

- `docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`

For how MCWF moment output enters dressed-frame squeezing post-analysis, use:

- `docs/instructions/solvers/mcwf/squeezing_moments.typ`


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

The primary post-processing pipeline should be moment-first. MCWF output should
be averaged as trajectory-resolved generator moments and converted into shared
observable containers before diagnostics or plotting. MFE and QuTiP retain
their method-specific converters into the same downstream observable
containers.

For parser output-container conventions, use:

- `docs/instructions/parser.typ`
- the `parser` skill

### 4.1 Moment Container

The notebook-level container is `MomentSeries`, defined in
`parser/moments.py`. It should be initialized on the shared
`t_eval` grid, from the metadata phase durations and `num_snapshots`, and then filled in as
post-processing steps are run:

```python
raw_result = run_trajectory_ensemble(...)
moments = MomentSeries(
    num_snapshots=num_snapshots,
    metadata=SimulationMetadata(
        Ni=raw_result.metadata.Ni,
        omega_i=raw_result.metadata.omega_i[:-1],
        Gamma=raw_result.metadata.Gamma,
        phase_protocol=raw_result.metadata.phase_protocol,
    ),
)
moments.generator_moments = raw_result.moments
moments.mean_generator_moments = compute_mean_generator_moments(
    moments.generator_moments
)
moments.J = compute_j_moments(
    moments.mean_generator_moments,
    metadata=moments.metadata,
)
moments.S = compute_s_moments(
    moments.mean_generator_moments,
    moments.J,
)
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
    include_full_system=False,
)
```

Current top-level fields are:

- `moments.t`: the shared saved-time grid.
- `moments.metadata`: validated physical model and standard-protocol context,
  including `Ni`, independent `omega_i`, completed `omega_groups`, `Gamma`,
  and the supplied `phase_protocol`.
- `moments.generator_moments`: trajectory-resolved first and symmetric second
  SU(3) generator moments, with complete trajectories retained along axis zero.
  Its shared saved-time grid and stored real moment arrays must be finite.
- `moments.mean_generator_moments`: the selected trajectory average of the
  same generator moments, retaining a singleton trajectory axis. Averaging
  automatically attaches first-moment sampling covariance and the selected
  sample count in its nullable `uncertainty` field, including cross-group
  correlations. A raw singleton leaves uncertainty unset.
- `moments.J`: a `JMomentSeries` containing first-order J-sphere moments plus
  derived J-vector direction fields and angles. Unsuffixed physical fields are
  tuples ordered by group, including for $G=1$; full-system fields live in
  `moments.J.aggregated`. Cartesian, population, length, normalized-direction,
  and angle fields have nullable `*_sem` counterparts, computed automatically
  when generator sampling statistics are present; nonlinear SEMs use
  first-order propagation and mark singular points with `NaN`. Both levels
  retain propagated affine-field sample covariance
  and sample count in `uncertainty`.
- `moments.J_modes`: an optional `JModeSeries` containing post-analysis common,
  population-weighted contrast, drive-bright, and drive-dark one/two-group
  vectors. Each vector can retain Cartesian sampling covariance and sample
  count, plus component and first-order length/xy-length SEMs.
- `moments.mfe_residuals`: an `MFEResidualSeries` containing two-group MFE
  residual diagnostics when computed from `moments.J`.
- `moments.squeezing`: an optional `SqueezingResult`. Its `group_resolved`
  tuple contains one covariance and generalized-squeezing series per physical
  group. Its `full_system` field is present only when the caller requests
  squeezing for the unweighted sum of all group fluctuation operators.
- `moments.S`: an optional `SMomentSeries` containing group-local S-sphere
  components, normalized directions, angles, and the effective zero/one basis
  used to define the S sphere.

Parser container conventions live in `docs/instructions/parser.typ` and the
`parser` skill.

Reusable Cartesian vector outputs should use `BlochVectorSeries` from
`parser/bloch_vector.py`, which derives its full and xy-plane lengths from the
stored x/y/z components. With optional Cartesian sampling uncertainty it also
derives the matching component and norm SEMs.

### 4.2 J-Sphere Moments

J-sphere moments are the common post-processing representation used across the
simulation backends. The solvers should stay responsible for producing raw
dynamics, while post-analysis converters turn those outputs into a shared
`JMomentSeries`. MCWF uses the solver-independent generator-moment converter;
MFE and QuTiP keep their method-specific converters because they do not need
the generator-moment path for S moments or squeezing.

Use these converter entry points:

```python
moments.J = compute_j_moments(
    moments.mean_generator_moments,
    metadata=moments.metadata,
)
moments.J = compute_mfe_j_moments(mfe_result)
moments.J = compute_qutip_j_moments(qutip_result)
```

The shared series may contain group spin components, normalized directions,
angles, atom-number fields, and counterparts in the nested `aggregated`
  container when those quantities are available. Jump rate is not part of the
J-moment output contract. Derived fields
such as spin directions and angles should be attached inside the J-moment
layer, so plotting and diagnostics can consume one common structure independent
of whether the data came from MCWF, MFE, or QuTiP. Additive aggregated fields
should be constructed through
`JMomentSeries.attach_aggregated_from_group_fields(...)` by summing their group
counterparts before nonlinear aggregated fields are attached. This class
method also owns aggregate covariance and SEM attachment, including
cross-group correlations. Attach group covariance before calling it; absent
input uncertainty leaves aggregate covariance and SEMs unset.

Shared spin-vector direction and angle conventions live in
`docs/instructions/common/spin_moment_conventions.typ`. The generic conversion
from averaged generator moments is defined in
`docs/instructions/post_analysis/j_moments.typ` and summarized in
`docs/instructions/post_analysis/post_analysis.md`.
The Cartesian conversion uses the shared, read-only `J_FROM_GENERATORS`
matrix from `post_analysis/j_moments.py`, also exported by `post_analysis`,
so later uncertainty calculations can reuse the same generator coefficients.
`compute_j_components(...)` extracts Cartesian components and populations,
attaches their group sampling covariance and SEMs, then delegates aggregate
values, covariance, and SEMs to the aggregation class method. Keep physical
values as named arrays and share their named coefficients
with covariance propagation, as described in the J-moment instructions.
`compute_j_moments(...)` calls it before attaching phase indices and
group directions and angles. The direction and angle class methods also
attach their SEMs with first-order Jacobian propagation from the stored
Cartesian covariance; aggregation uses the same shared helpers. Undefined
derivatives produce `NaN` SEMs. Deterministic converters leave SEMs unset.

Numerical MFE solving and MFE-to-J-moment conversion are defined in:

- `docs/instructions/solvers/mfe/mfe-solver.typ`

A task-specific QuTiP-to-J-moment instruction has not yet been added.

### 4.3 S-Sphere Moments

For generator-moment workflows, construct group-local S moments directly from
the averaged first generator moments and the corresponding group-resolved J
angles:

```python
moments.S = compute_s_moments(
    moments.mean_generator_moments,
    moments.J,
)
```

The converter uses one shared J-angle-defined S frame per group and saved time,
then delegates lengths, normalized directions, and S angles to the
`SMomentSeries` class methods. Symmetric second generator moments remain the
input for later covariance and squeezing analysis. The generic converter is
defined in `docs/instructions/post_analysis/s_moments.typ` and summarized in
`docs/instructions/post_analysis/post_analysis.md`.

### 4.4 Dressed Generator Moments

Construct the dressed basis from the group-resolved J and S angles, convert it
to the real adjoint rotation of the normalized Gell--Mann generators, and
rotate both first and symmetric second moments:

```python
U = build_dressed_basis_transformation(moments.J, moments.S)
R = compute_gell_mann_rotation(U)
mean_dressed_generator_moments = rotate_generator_moments(
    moments.mean_generator_moments,
    R,
)
```

The rotation preserves the input trajectory axis and all cross-group
second-moment blocks. `U` and `R` are transient derived arrays rather than
fields of `MomentSeries`; the rotated result reuses
`TrajectoryResolvedMomentSeries`. Basis construction is defined in
`docs/instructions/post_analysis/dressed_basis.typ`, and generator/moment
rotation is defined in
`docs/instructions/post_analysis/generator_rotation.typ`. Both are summarized
in `docs/instructions/post_analysis/post_analysis.md`. Keep the lab and dressed
series distinct: the fixed generator mappings in `compute_j_moments(...)` and
`compute_s_moments(...)` require the lab-frame series.

### 4.5 Generalized Squeezing

Compute squeezing from the averaged dressed generator moments and the shared J
series:

```python
moments.squeezing = compute_squeezing(
    mean_dressed_generator_moments,
    moments.J,
    metadata=moments.metadata,
    include_full_system=False,
)
```

The caller must construct `mean_dressed_generator_moments` explicitly as in
Section 4.4; `compute_squeezing(...)` does not infer or apply a basis rotation.
It forms the full `4G x 4G` covariance of `(T1, T2, T4, T5)` and always
computes the diagonal `4 x 4` block for every physical group. With
`include_full_system=True`, it additionally sums every group-pair covariance
block and computes squeezing for the unweighted sum of corresponding group
fluctuation operators. This option must not change the group-resolved results.
The detailed contract is in
`docs/instructions/post_analysis/squeezing.typ`.


## 5. Diagnostics

Diagnostics should consume already-computed moment series. If a new diagnostic
requires additional simulation data, extend the stored moment contract
explicitly. Reusable post-simulation physics diagnostics should live in the
root-level `post_analysis/` package.

Current task-specific diagnostic instructions include:

- `docs/instructions/post_analysis/j_modes.typ`
- `docs/instructions/post_analysis/harmonic_analysis.typ`
- `docs/instructions/post_analysis/mfe_residuals.typ`
- `docs/instructions/post_analysis/squeezing.typ`

The generic `compute_sample_covariance(...)`, `propagate_covariance(...)`,
and `std_from_covariance(...)` helpers live in `common/utils/uncertainty.py`,
with compatibility exports through `post_analysis`. Their contracts and the
minimal `SamplingUncertainty` container are specified in
`docs/instructions/post_analysis/uncertainty.typ`, which includes the matrix
form, Jacobian evaluation at the mean, and a worked two-variable example.
`compute_mean_generator_moments(...)` calls the helper automatically for
trajectory first moments and stores only sample covariance and sample count,
separately from quantum second moments. Mean covariance and SEM are calculated
downstream when needed. Observable extraction should construct associated
standard deviations using the shared helpers and store SEMs where requested.
The J affine-field and first-order nonlinear integration is implemented in
`docs/instructions/post_analysis/j_moments.typ`. For J modes, prefer
`compute_j_modes_from_j_moments(moments.J, omega_groups=...)`, which owns active
population normalization and propagates the full J covariance, including
population-dependent weights and cross-group correlations. Its contract is
in `docs/instructions/post_analysis/j_modes.typ`. The array-based
`compute_j_modes(...)` retains its central-value interface. S/squeezing
uncertainties remain future work.
Apply sampling calculations to eligible stochastic ensembles;
deterministic MFE outputs remain unchanged.

The public post-analysis entry points are summarized in
`docs/instructions/post_analysis/post_analysis.md`.

For offsets and oscillation amplitudes over time, use
`compute_windowed_harmonic_analysis(...)` from `post_analysis/harmonic_analysis.py`.
It reuses the existing fitter in complete caller-sized windows and returns
window-centre times with per-curve metric arrays, including THD for reuse of
the existing harmonic plotter. The bright/dark MCWF notebook uses this
windowed analysis for phase-2 traces and plots the results against time.
Use adequate saved-time resolution and flag poor periodic fits before
interpreting changes in amplitude or frequency. The notebook prints failure
reasons and R-squared values through `compute_harmonic_quality(...)`, which
returns the boolean display masks. It shows flagged candidate fits as crosses, and
leaves gaps only where no candidate exists. Fit quality alone is not
statistical evidence for an oscillation. Early-time damping is orchestrated
separately by `compute_short_time_dark_mode_analysis(...)` in
`post_analysis/short_time_oscillations.py`. It selects the early phase-2 data,
resamples whole cached trajectories through the established averaging-to-J-to-mode
pipeline, and calls the existing bootstrap fitter. The result contains the
curves, optional SEMs, and fits for direct plotting; notebook cells only
configure, call, and plot this helper. Inspect fit quality and
the model-conditional intervals before claiming damping. Both analyses reuse
cached simulation results; do not rerun parameter initialization or simulation
to refresh their imports. Validation notes are in
`docs/investigations/windowed_dark_mode_harmonics.md`.
Window and output conventions are defined in
`docs/instructions/post_analysis/harmonic_analysis.typ`.

## 6. Plotting
The shared spin-component plot lives in `common/plotting/j_spin.py`:
- `plot_spin_components(t, x_components, y_components, z_components, lengths,
  labels=..., x_sem=None, y_sem=None, z_sem=None, length_sem=None,
  spin_symbol="J", ...)`: plots caller-selected vectors with optional matching
  mean-plus/minus-SEM bands with darker colors and stronger opacity, and uses
  the supplied spin symbol in the component labels. MCWF notebook calls
  plotting raw J components should pass their corresponding SEM fields;
  derived observables require their own errors.

The transverse spin-component plot lives in `common/plotting/specific.py`:
- `plot_transverse_spin_components(t, x_components, y_components,
  z_components, lengths, labels=..., ...)`: plots a selected transverse
  combination, z, and full length.

The shared angle plot lives in `common/plotting/j_spin.py`:
- `plot_bloch_angles(t, theta_curves, phi_curves, labels=..., title=..., ...)`:
  plots caller-selected angle curves using the selected palette and
  `linestyle`.

The shared MFE residual diagnostic plotting function lives in
`common/plotting/mfe_residuals.py`:
- `plot_mfe_residuals(moments.mfe_residuals, ...)`: plots stored two-group
  residuals in a single residual panel with the L2 norm.

The shared harmonic-sweep plot lives in
`common/plotting/harmonic_analysis.py`:
- `plot_harmonic_sweep(...)`: plots frequency, RMS amplitude, offset, and THD
  against a parameter or window-centre time; optional `quality_masks` mark
  flagged candidates as crosses.
- `plot_damped_oscillations(...)`: overlays early MCWF curves and SEM shading
  with already-computed damped fits, backgrounds, and envelopes. Optional
  reference inputs add dashed cached MFE curves on the same interval.

The shared generalized-squeezing plot lives in
`common/plotting/squeezing.py`:
- `plot_squeezing(t, xi2, lambda_min, N_c, excited_fraction_active,
  labels=..., covariance_eigvals=..., plot_all_eigenvalues=False, ...)`: plots
  one or more squeezing results on the shared `2x2` diagnostic grid, with an
  opt-in view of all four stored covariance eigenvalues.

Repo-specific plotting-function contracts live in
`docs/instructions/plotting_workflows.md`.

Future diagnostics and plots should consume `MomentSeries` or `JMomentSeries`
when the required data are already present, instead of rerunning older observable extraction steps.

## TODOs

Open implementation cleanup notes are tracked in `docs/TODOs.md`.
