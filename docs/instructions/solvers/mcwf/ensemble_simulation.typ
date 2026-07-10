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

For fixed `MCWFSolverParameters`, a saved-time grid `t_eval`, and `ntraj`
trajectories, the ensemble is a collection of independent MCWF unravelings:

$
cal(E) = {T_1, dots.c, T_(n_"traj")}.
$

All trajectories should share the same completed model inputs, precomputed
sector data, and saved-time grid so that later moment extraction can average
snapshots by index without interpolation.

The ensemble runner is responsible for:

- completing group couplings once;
- constructing and validating the initial sector coefficients;
- building reusable precomputed data once;
- creating reproducible child seeds;
- dispatching `_simulate_single_trajectory(...)` serially or with
  multiprocessing;
- returning a `TrajectoryEnsemble` with shared metadata.

= Method in Pseudo-Code

```python
def run_trajectory_ensemble(parameters, *, t_eval, ntraj, seed=None,
                            n_processes=None, chunksize=1, verbose=False):
    validate runtime inputs
    omega_i = complete final group coupling
    sector_coeffs = centered_sector_initial_coeffs(...)
    check_initial_sector_omega_ratio(...)

    precomputed = build_precomputed_trajectory_data(...)
    seed_sequences = np.random.SeedSequence(seed).spawn(ntraj)

    if serial:
        trajectories = [
            _simulate_single_trajectory(..., seed_sequence=child,
                                        precomputed=precomputed)
            for child in seed_sequences
        ]
    else:
        initialize worker-local state once
        trajectories = pool.imap(_simulate_single_trajectory_worker,
                                 seed_sequences, chunksize=chunksize)

    return TrajectoryEnsemble(
        trajectories=trajectories,
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
    Ni,
    dN,
    omega_i,
    Gamma,
    phases,
    sector_distribution="binomial",
    dt=1e-3,
    shifted_jump_operator=False,
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

`omega_i` should contain the first `G-1` group couplings. The ensemble runner
should append the final coupling with `omega_G_from_weighted_average(...)`
before constructing precomputed data. For single-group runs, use `Ni=[N]` and
`omega_i=[]`.

Initial-sector conventions live in
`docs/instructions/solvers/mcwf/initial_sector_state.typ`.

Precompute contents and sector-operator conventions live in:

- `docs/instructions/solvers/mcwf/simulation_precompute.typ`
- `docs/instructions/solvers/mcwf/sector_operators.typ`

Single-trajectory propagation and snapshot saving live in
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ`.

= Output

The returned object should be:

```python
TrajectoryEnsemble(
    trajectories=trajectories,
    seeds=seed_keys,
    metadata=metadata,
)
```

`TrajectoryEnsemble.trajectories` should preserve submission order.
`TrajectoryEnsemble.metadata` should store shared simulation metadata such as
`Ni`, completed `omega_i`, `Gamma`, `phases`, `shifted_jump_operator`, `t_eval`,
sector keys, sector multiplicities, and sector dimensions.

Each trajectory should already be a complete `TrajectoryResult`; see
`docs/instructions/solvers/mcwf/single_trajectory_simulation.typ` for its
fields.

= Invariants and Edge Cases

- Build `precomputed` exactly once per `run_trajectory_ensemble(...)` call.
- Serial and multiprocessing paths should produce the same seed semantics,
  trajectory ordering, and output structure.
- Each trajectory should receive exactly one child `SeedSequence`.
- Worker-local precomputed state should be treated as read-only.
- All trajectories in one ensemble should use the same explicit `t_eval` grid.
- `n_processes=None` and `n_processes=1` should use the serial path;
  `n_processes=-1` should mean all available CPUs.
- Runtime diagnostics such as mean total steps and mean non-precomputed steps
  may be printed, but must not affect trajectory physics.
