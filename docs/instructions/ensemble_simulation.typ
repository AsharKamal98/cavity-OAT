#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Ensemble Simulation: Implementation Instructions]
]

= Scope

This file describes the intended implementation of the ensemble-level custom
MCWF runner. It should be read together with
`docs/instructions/simulation_precompute.typ` and
`docs/instructions/single_trajectory_simulation.typ`.

The ensemble layer should not change the physical model. It should:

- construct reproducible trajectory seeds;
- build reusable phase/sector data once;
- run many independent `_simulate_single_trajectory(...)` calls;
- return a `TrajectoryEnsemble` on one common saved-time grid.

= Method

For fixed
$
N_i, omega_i, Gamma, {"phases"}, d t,
$
the ensemble should consist of `ntraj` independent unravelings of the same
piecewise-constant non-Hermitian evolution and jump process:

$
T_k = "_simulate_single_trajectory"(..., {"seed_sequence"}_k, {"precomputed"}),
quad k = 1, dots.c, n_"traj".
$

The ensemble output should therefore be
$
cal(E) = {T_1, dots.c, T_(n_"traj")},
$
with all trajectories sharing the same `t_eval` convention so that ensemble
observables can be averaged snapshot-by-snapshot without interpolation.

The strong-symmetry decomposition is already handled inside the single-trajectory
solver and precompute layer. The ensemble code should only coordinate repeated,
independent calls to that solver.

= Data In

The ensemble entry point should be:

```python
MCWFSolverParameters(
    Ni,
    dN=0,
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
    seed=None,
    ntraj,
    n_processes=None,
    chunksize=1,
    verbose=False,
) -> TrajectoryEnsemble
```

The same homogeneous and inhomogeneous sector-key conventions from
`docs/instructions/simulation_precompute.typ` should be supported here.

`MCWFSolverParameters` should validate at least:

- `Ni` is nonempty and all group sizes are non-negative;
- `len(omega_i) = len(Ni) - 1`;
- `shifted_jump_operator=True` implies `Gamma > 0`;
- `dt > 0`.

The ensemble entry point should validate at least:

- `ntraj > 0`;
- `t_eval` is one-dimensional, strictly increasing, starts at `0`, and ends at
  the total protocol time;
- `n_processes` is one of `None`, `1`, `-1`, or a positive integer.

The ensemble layer should take `Ni` and the first `G-1` couplings in
`omega_i`, complete the final coupling once through
`omega_G_from_weighted_average(...)`, and then pass the completed `omega_i`
list together with `Ni` to precompute and single-trajectory helpers. For
single-group runs, this means `Ni=[N]` and `omega_i=[]`.

If `sector_coeffs is None`, the ensemble runner should construct centered
initial coefficients itself through `centered_sector_initial_coeffs(...)` using
`dN` and `sector_distribution`. It should then run the same
`check_initial_sector_omega_ratio(...)` validation that previously lived in the
deleted `run_h_sim(...)` / `run_inh_sim(...)` wrappers.

= Seed Construction and Reproducibility

The ensemble seed convention should be based on one parent `SeedSequence`:

```python
parent_seed_sequence = np.random.SeedSequence(seed)
seed_sequences = parent_seed_sequence.spawn(ntraj)
```

Each child `SeedSequence` should be used for exactly one trajectory.

This convention should preserve the current reproducibility rule:

- trajectory `0` in `run_trajectory_ensemble(..., seed=seed, ntraj>=1)` should
  use the same child seed as a direct
  `_simulate_single_trajectory(..., seed_sequence=child)` call;
- serial and multiprocessing ensemble runs should use the same child seed list;
- the returned `TrajectoryEnsemble.seeds` should store
  `tuple(child.spawn_key)` for each trajectory in output order.

The ensemble layer should not draw random numbers itself beyond seed
construction. All physical randomness should remain inside
`_simulate_single_trajectory(...)`.

= Precompute Reuse

The ensemble code should build reusable phase/sector data once:

```python
precomputed = build_precomputed_trajectory_data(
    Ni=Ni,
    omega_i=omega_i,
    Gamma=Gamma,
    phases=phases,
    sector_coeffs=sector_coeffs,
    dt=dt,
    shifted_jump_operator=shifted_jump_operator,
)
```

This should happen before entering either the serial or multiprocessing path.
The same `precomputed` dictionary should then be reused by every trajectory in
the ensemble run.

For the contents and indexing conventions of `precomputed`, use
`docs/instructions/simulation_precompute.typ`.

= Execution Paths

== Serial Path

If `n_processes is None` or `n_processes == 1`, the ensemble should run
serially:

```python
trajectories = [
    _simulate_single_trajectory(..., seed_sequence=child, precomputed=precomputed)
    for child in seed_sequences
]
```

This path should be kept as the simplest reference behavior for debugging and
reproducibility checks.

The current code shows a progress bar with `tqdm` in serial mode; the
instruction should preserve that behavior unless the user explicitly asks for a
quiet path.

== Multiprocessing Path

If `n_processes == -1`, the ensemble should use `mp.cpu_count()`.

If `n_processes > 1`, the code should create a multiprocessing pool and
initialize worker-local read-only state once through:

```python
_init_trajectory_worker(
    Ni,
    Gamma,
    phases,
    sector_coeffs,
    dt,
    t_eval,
    shifted_jump_operator,
    precomputed,
    omega_i,
)
```

The worker initializer should store these objects in a process-local global
dictionary such as `_WORKER_STATE`.

Per-task payloads should stay minimal:

```python
pool.imap(_simulate_single_trajectory_worker, seed_sequences, chunksize=chunksize)
```

Each worker task should receive only one child `SeedSequence`. Large objects
such as `precomputed` should be reused from `_WORKER_STATE` rather than sent
with every task.

The multiprocessing path should preserve the same trajectory order as the input
`seed_sequences` list so the returned `TrajectoryEnsemble.seeds` stays aligned
with `TrajectoryEnsemble.trajectories`.

== Chunksize

`chunksize` should be passed directly to `pool.imap(...)`.

Small `chunksize`, especially `1`, should remain the default because the number
of jumps can vary strongly across trajectories and better load balancing is more
important than minimizing task-dispatch overhead for the current workflow.

= Method in Pseudo-Code

```python
def run_trajectory_ensemble(...):
    validate ensemble runtime inputs
    complete omega_i by appending omega_G_from_weighted_average(...)
    sector_coeffs = centered_sector_initial_coeffs(...)
    run check_initial_sector_omega_ratio(...)

    parent_seed_sequence = np.random.SeedSequence(seed)
    seed_sequences = parent_seed_sequence.spawn(ntraj)
    seed_keys = [tuple(child.spawn_key) for child in seed_sequences]

    precomputed = build_precomputed_trajectory_data(...)

    if n_processes is None or n_processes == 1:
        trajectories = [
            _simulate_single_trajectory(..., seed_sequence=child, precomputed=precomputed)
            for child in seed_sequences
        ]
    else:
        if n_processes == -1:
            n_processes = mp.cpu_count()
        initialize worker-local state once
        trajectories = list(pool.imap(worker, seed_sequences, chunksize=chunksize))

    total_steps = sum(traj.total_step_count for traj in trajectories)
    non_precomputed_steps = sum(traj.non_precomputed_step_count for traj in trajectories)
    print average step diagnostics

    return TrajectoryEnsemble(
        trajectories=trajectories,
        seeds=seed_keys,
        metadata=metadata,
    )
```

The ensemble layer should reuse the existing helpers:

- `build_precomputed_trajectory_data(...)` for shared phase/sector operators;
- `_simulate_single_trajectory(...)` for trajectory physics and snapshot saving;
- `_init_trajectory_worker(...)` and `_simulate_single_trajectory_worker(...)`
  for multiprocessing state reuse.

Do not duplicate the single-trajectory propagation logic inside the ensemble
runner.

= Data Out

The returned object should be:

```python
TrajectoryEnsemble(
    trajectories=trajectories,
    seeds=seed_keys,
    metadata=metadata,
)
```

`TrajectoryEnsemble.trajectories` should preserve submission order, not sorted
post hoc by jump count, runtime, or seed.
`TrajectoryEnsemble.metadata` should store shared simulation metadata such as
`Ni`, `omega_i`, `Gamma`, `phases`, `shifted_jump_operator`, `t_eval`,
`sectors`, `sector_multiplicities`, and `sector_dimensions`.

Each element should already be a complete `TrajectoryResult`, including:

- saved `snapshots`;
- `jump_times` and `jump_count`;
- `final_sector_blocks`;
- runtime step counters.

= Invariants and Edge Cases

- The ensemble layer should build `precomputed` exactly once per call to
  `run_trajectory_ensemble(...)`, not once per trajectory.
- Serial and multiprocessing paths should differ only in execution strategy, not
  in seed semantics or returned data structure.
- All trajectories in one ensemble should use the same explicit `t_eval` grid,
  because downstream ensemble averaging assumes aligned saved times.
- Worker-local global state should be treated as read-only during simulation.
  Do not mutate the shared `precomputed` dictionary inside worker tasks.
- The current code prints a step-summary diagnostic at the end of both serial
  and multiprocessing runs:

  ```python
  total steps = mean_k traj_k.total_step_count
  non-precomputed steps = mean_k traj_k.non_precomputed_step_count
  ```

  These counters are runtime diagnostics only and should not affect trajectory
  physics.
- `verbose=True` should control setup timing prints such as precompute and pool
  startup. It should not disable the trajectory progress bar in the current
  implementation.
