#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Ensemble Simulation: Implementation Instructions]
]

= Purpose

This file describes the ensemble-level custom MCWF runner implemented in
`solvers/mcwf/ensamble_sim.py`. Use it when editing
`run_trajectory_ensemble(...)`, ensemble seeding, multiprocessing dispatch, or
the `TrajectoryEnsemble` output contract.

The ensemble layer coordinates repeated single-trajectory runs. It should not
redefine the trajectory physics, sector operators, or precomputed propagators;
those rules live in the MCWF instruction files referenced below.

= Method

For fixed `MCWFSolverParameters(...)`, a saved-time grid `t_eval`, and `ntraj`
trajectories, the ensemble is a collection of independent MCWF unravelings:

$
cal(E) = {T_1, dots.c, T_(n_"traj")}.
$

All trajectories should share the same completed model inputs, precomputed
sector data, and saved-time grid so that post-analysis can combine or average
their stored moment rows by index without interpolation.

The ensemble runner is responsible for:

- constructing and validating the initial sector coefficients;
- building reusable precomputed data once;
- creating reproducible child seeds;
- dispatching `_simulate_single_trajectory(...)` serially or with
  multiprocessing;
- copying each singleton trajectory moment result into its corresponding row
  of the ensemble moment arrays without averaging;
- consuming single-trajectory diagnostics for optional runtime summaries and
  then discarding them;
- returning a `TrajectoryEnsemble` with shared metadata.

= Method in Pseudo-Code

```python
def run_trajectory_ensemble(parameters, *, t_eval, ntraj, seed=None,
                            n_processes=None, chunksize=1, verbose=False):
    validate runtime inputs
    integration_phases = parameters.phase_protocol.integration_phases
    sector_coeffs = centered_sector_initial_coeffs(...)
    check_initial_sector_omega_ratio(...)

    precomputed = build_precomputed_trajectory_data(...)
    seed_sequences = np.random.SeedSequence(seed).spawn(ntraj)

    allocate first_moments with shape (ntraj, n_t, G, 8)
    allocate second_moments with shape (ntraj, n_t, K)

    if serial:
        results = (
            _simulate_single_trajectory(..., seed_sequence=child,
                                        precomputed=precomputed)
            for child in seed_sequences
        )
    else:
        initialize worker-local state once
        results = pool.imap(_simulate_single_trajectory_worker,
                            seed_sequences, chunksize=chunksize)

    for r, result in enumerate(results):
        verify result.moments uses the shared t_eval and has R=1
        first_moments[r] = result.moments.first_moments[0]
        second_moments[r] = result.moments.second_moments[0]
        accumulate optional step-count summary from result.diagnostics

    return TrajectoryEnsemble(
        moments=TrajectoryResolvedMomentSeries(
            t=t_eval,
            first_moments=first_moments,
            second_moments=second_moments,
        ),
        seeds=[tuple(child.spawn_key) for child in seed_sequences],
        metadata=metadata,
    )
```

`run_trajectory_ensemble(...)` is the public MCWF entry point. It creates the
shared setup, then delegates trajectory physics to
`_simulate_single_trajectory(...)`. The multiprocessing helpers should only
reuse worker-local copies of the same setup; they should not change the
simulation semantics.

= Data Requirements

The ensemble entry point should receive:

```python
MCWFSolverParameters(
    Ni=metadata.Ni,
    omega_i=metadata.omega_groups,
    Gamma=metadata.Gamma,
    phase_protocol=metadata.phase_protocol,
    dN,
    sector_distribution="binomial",
    dt=1e-3,
    shifted_jump_operator=True,
)

run_trajectory_ensemble(
    parameters,
    *,
    t_eval,
    ntraj,
    seed=None,
    n_processes=None,
    chunksize=1,
    verbose=False,
) -> TrajectoryEnsemble
```

The public MCWF parameters default to the shifted jump-operator picture. Set
`shifted_jump_operator=False` explicitly for an unshifted ensemble.

`SimulationMetadata.omega_i` should contain the first `G-1` group couplings.
Its validator constructs the final weighted-average coupling as
`metadata.omega_groups`; pass that completed list as the solver's `omega_i`.
For single-group runs, use `Ni=[N]` and `omega_i=[]`.

Initial-sector conventions live in
`docs/instructions/solvers/mcwf/initial_sector_state.typ`.

Precompute contents and sector-operator conventions live in:

- `docs/instructions/solvers/mcwf/simulation_precompute.typ`
- `docs/instructions/solvers/mcwf/sector_operators.typ`

Single-trajectory propagation and moment saving live in
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`.

= Output

The returned object should be:

```python
TrajectoryEnsemble(
    moments=TrajectoryResolvedMomentSeries(
        t=...,                       # shape (n_t,)
        first_moments=...,           # shape (ntraj, n_t, G, 8)
        second_moments=...,          # shape (ntraj, n_t, K)
    ),
    seeds=seed_keys,
    metadata=metadata,
)
```

The first axis of `TrajectoryEnsemble.moments` should preserve submission
order. `TrajectoryEnsemble.seeds[r]` must identify the same trajectory as
moment row `r`.
`TrajectoryEnsemble.metadata` should store shared simulation metadata such as
`Ni`, completed `omega_i`, `Gamma`, `phase_protocol`, `shifted_jump_operator`,
`t_eval`, sector keys, sector multiplicities, and sector dimensions.

Each private single-trajectory call should return a complete
`TrajectoryResult`; see
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ` for its
fields. The ensemble should retain only its moment arrays. Jump times, jump
counts, and step counts are MCWF runtime diagnostics and should not be returned
by the parent ensemble.

For array simulations that require ensemble-wide derived bases, save each
`TrajectoryEnsemble` moment chunk temporarily in shared storage. After every
array job has finished, load the chunks in array-index order, verify their
shared metadata and saved-time grid, and concatenate their moment arrays along
the trajectory axis into one `TrajectoryEnsemble`. Compute J, S, and squeezing
from that combined ensemble rather than averaging separately derived array
outputs.

= Invariants and Edge Cases

- Build `precomputed` exactly once per `run_trajectory_ensemble(...)` call.
- Serial and multiprocessing paths should produce the same seed semantics,
  trajectory ordering, and output structure.
- Each trajectory should receive exactly one child `SeedSequence`.
- The trajectory axis must never be averaged inside the ensemble runner.
- Single-trajectory results must be copied into their ensemble rows as they are
  received so their diagnostics need not remain resident in the parent.
- `TrajectoryEnsemble.seeds` and the first moment axis must have equal length.
- Worker-local precomputed state should be treated as read-only.
- All trajectories in one ensemble should use the same explicit `t_eval` grid.
- `t_eval` should contain only finite values, and `chunksize` should be positive
  even when the serial path is selected.
- The physical-group axis of returned moments should match both `metadata.Ni`
  and the completed coupling list in `metadata.omega_i`.
- Concatenated array ensembles must share physical, phase, sector, and time-grid
  metadata.
- `n_processes=None` and `n_processes=1` should use the serial path;
  `n_processes=-1` should mean all available CPUs.
- Runtime diagnostics such as mean total steps and mean non-precomputed steps
  may be printed, but must not affect trajectory physics.
