#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Single-Trajectory Simulation: Implementation Instructions]
]

= Scope

This file describes the intended implementation of
`simulate_single_trajectory(...)` in the custom strong-symmetry MCWF solver.
It should be read together with
`docs/instructions/simulation_precompute.typ`.

The single-trajectory layer should:

- build one trajectory-local RNG from one child seed;
- evolve the direct-sum sector wavefunction through the piecewise-constant
  protocol;
- save snapshots exactly on the common `t_eval` grid;
- record jump times and final sector blocks in a `TrajectoryResult`.

= Method

For one trajectory, the reduced sector-block wavefunction should evolve between
jumps under the non-Hermitian effective generator
$
H_{"eff"},
$
sector by sector:

$
tilde(psi)_s(t + Delta t)
=
exp[-i H_{"eff"},s Delta t] psi_s(t).
$

The MCWF survival probability over the attempted step should be computed from
the norm of the unnormalized propagated state:

$
P_"surv"(Delta t) = sum_s ||tilde(psi)_s(t + Delta t)||^2.
$

With one random threshold
$
r in [0,1),
$
the step should be accepted as a no-jump step when
$
P_"surv"(Delta t) > r.
$
Otherwise, a jump should be placed inside the step, the state should be
propagated to the refined jump time, renormalized, acted on by the phase-local
jump operator, and renormalized again.

The sector-resolved effective generators, jump operators, and full-`dt`
propagators used here are defined in
`docs/instructions/simulation_precompute.typ`.

= Data In

The single-trajectory entry point should be:

```python
simulate_single_trajectory(
    N,
    Gamma,
    phases,
    sector_coeffs,
    *,
    internal_sector_states=None,
    dt=1e-3,
    num_snapshots=101,
    seed=1234,
    seed_sequence=None,
    shifted_jump_operator=False,
    precomputed=None,
    omega_1=None,
    N1=None,
    N2=None,
) -> TrajectoryResult
```

The function should validate at least:

- `N > 0`;
- `Gamma >= 0`;
- `dt > 0`;
- `num_snapshots >= 2`;
- `shifted_jump_operator=True` implies `Gamma > 0`.

If any sector key is a tuple, the solver should require:

- `omega_1 is not None`;
- `N1 is not None` and `N2 is not None`;
- `N1 + N2 = N`.

= Data Out

The returned `TrajectoryResult` should contain:

```python
TrajectoryResult(
    N=N,
    Gamma=Gamma,
    phases=list(phases),
    shifted_jump_operator=shifted_jump_operator,
    t_eval=t_eval,
    sectors=sector_list,
    sector_multiplicities=multiplicities,
    final_sector_blocks=...,
    snapshots=...,
    jump_times=...,
    jump_count=len(jump_times),
    sector_dimensions=dims,
    omega_1=omega_1,
    omega_2=omega_2,
    N1=N1,
    N2=N2,
    total_step_count=...,
    non_precomputed_step_count=...,
)
```

For inhomogeneous runs, `omega_2` should be derived once at the end through
`omega2_from_weighted_average(...)`. For homogeneous runs it should remain
`None`.

= Initialization

== RNG

If `seed_sequence is None`, the current code should continue to build the
trajectory-local RNG through:

```python
seed_seq = np.random.SeedSequence(seed).spawn(1)[0]
rng = np.random.default_rng(seed_seq)
```

If `seed_sequence` is provided, it should override `seed`.

This should stay consistent with the ensemble seeding convention documented in
`docs/instructions/ensemble_simulation.typ`.

== Precompute and Sector Ordering

If `precomputed is None`, the solver should call
`build_precomputed_trajectory_data(...)` itself. Otherwise it should reuse the
supplied dictionary without rebuilding it.

The hot loop should read the following aligned objects from `precomputed`:

- `sector_list`;
- `multiplicities`;
- `dims`;
- `phase_jump_operators`;
- `phase_generators`;
- `phase_propagators`.

The sector-block state inside the trajectory loop should therefore be kept as an
ordered list:

```python
psi_blocks = [psi_sector for sector in sector_list]
```

while snapshots and result outputs should convert back to dictionaries through
`blocks_list_to_dict(...)`.

== Initial Wavefunction

The initial state should be built through the existing helper:

```python
initial_blocks = build_initial_sector_state(N, sector_coeffs, internal_sector_states)
psi_blocks = renormalize_psi_blocks([initial_blocks[key] for key in sector_list])
```

`build_initial_sector_state(...)` should remain authoritative for:

- normalized sector coefficients;
- default `|n_e = 0>` or `|(n_(e,1), n_(e,2)) = (0,0)>` internal states;
- validation of supplied internal block shapes.

The common saved-time grid should be constructed once through:

```python
t_eval = build_t_eval_from_phases(phases, num_snapshots)
```

which currently uses `np.linspace(0.0, total_time, num_snapshots)`.

= Snapshot Construction

The first snapshot should always be saved at time `0.0`:

```python
TrajectorySnapshot(
    time=0.0,
    sector_blocks=blocks_list_to_dict(sector_list, psi_blocks),
    norm=1.0,
    phase_index=0,
)
```

Later snapshots should be saved only when the internal evolution lands exactly
on the next requested `t_eval` value up to the solver tolerance.

The helper logic should continue to enforce:

- save exactly one snapshot per requested `t_eval` point;
- raise an error if the step-splitting logic misses a requested `t_eval` time;
- store the current `phase_index` with each saved snapshot.

An important current convention is that the saved `sector_blocks` need not be
renormalized between jumps. On no-jump evolution segments, the code currently
saves the unnormalized propagated state and stores
$
"norm" = sqrt(sum_s ||psi_s||^2).
$
Downstream observable code then divides by the saved block norm when computing
expectation values. This convention should not change silently.

= Time-Stepping Logic

For each nonzero-duration phase, the solver should load:

```python
jump_operators_list = phase_jump_operators[phase_index]
generators_list = phase_generators[phase_index]
full_step_propagators = phase_propagators[phase_index]
```

The attempted step length should be:

```python
step = min(dt, phase_end - current_time, next_eval_time - current_time)
```

This ensures the solver lands exactly on:

- the base internal grid when possible;
- phase boundaries;
- the next requested saved time.

== Full-Step Path

If
```python
abs(step - dt) <= 1e-15
```
the solver should use the precomputed full-step propagators:

```python
trial = propagate_blocks_with_propagators(psi_blocks, full_step_propagators)
```

This is the fast path for ordinary interior full-`dt` evolution.

== Variable-Step Path

If `step != dt` within tolerance, the solver should use:

```python
trial = propagate_blocks(psi_blocks, generators_list, step)
```

This path should be used for:

- phase-boundary steps;
- `t_eval`-boundary steps;
- jump bisection midpoint propagations;
- propagation to the refined jump time;
- the optional post-jump remainder step.

= Jump Detection and Refinement

At the beginning of the trajectory, the solver should draw one threshold:

```python
threshold = rng.random()
```

For each attempted step, it should compute:

```python
trial_norm2 = total_norm2_list(trial)
```

If `trial_norm2 > threshold`, the step is accepted as a no-jump step.
The solver should then:

- set `psi_blocks = trial`;
- advance `current_time += step`;
- call the snapshot-saving helper.

If `trial_norm2 <= threshold`, the current code should resolve one jump inside
that attempted step by bisection.

== Jump Bisection

The current implementation should keep the existing fixed-depth bisection:

```python
lo, hi = 0.0, step
for _ in range(5):
    mid = 0.5 * (lo + hi)
    mid_state = propagate_blocks(pre_blocks, generators_list, mid)
    if total_norm2_list(mid_state) > threshold:
        lo = mid
    else:
        hi = mid
tau = hi
```

This refines the jump location to one of the size-`step / 2^5` subintervals.
If a more accurate root finder is ever desired, that should be introduced as an
explicit behavior change rather than silently replacing the current five-step
bisection rule.

== Applying the Jump

After finding `tau`, the solver should:

```python
psi_blocks = propagate_blocks(pre_blocks, generators_list, tau)
current_time += tau
psi_blocks = renormalize_psi_blocks(psi_blocks)
psi_blocks = apply_jump(psi_blocks, jump_operators_list)
jump_times.append(current_time)
threshold = rng.random()
maybe_save_snapshot()
```

`apply_jump(...)` should continue to apply the same phase-local jump operators
used for the jump-rate and precompute logic. If the jumped state has zero norm,
the current helper keeps the pre-jump state instead of crashing. This fallback
should remain an explicit runtime-safety convention.

== Post-Jump Remainder

After a jump, the current code attempts to use the remainder
$
"remainder" = "step" - tau
$
only if it does not immediately contain another jump:

```python
trial = propagate_blocks(psi_blocks, generators_list, remainder)
if total_norm2_list(trial) <= threshold:
    pass
else:
    psi_blocks = trial
    current_time += remainder
    maybe_save_snapshot()
```

This means the current implementation resolves at most one explicitly refined
jump per attempted outer-loop step. If the post-jump remainder would already
contain another jump, the solver leaves the state at the refined jump time and
lets the next outer-loop iteration handle subsequent evolution. The instruction
file should preserve this current behavior unless the user explicitly asks to
change it.

= Phase and Final-Time Handling

- Negative phase durations should raise an error.
- Zero-duration phases should be skipped.
- After all phases, the solver should renormalize the final state before
  constructing `final_sector_blocks`.
- If the last `t_eval` point still needs to be saved, the solver should require
  that `current_time` matches it within tolerance and then save it.
- If not all requested `t_eval` points were saved, the solver should raise an
  error.
- The current code may append one extra final snapshot only if the last saved
  snapshot time does not already match `current_time`.

= Step Counters

The solver should continue to track:

```python
total_step_count
non_precomputed_step_count
```

`total_step_count` should count every actual propagation call, including:

- ordinary full-`dt` propagation;
- partial-step propagation;
- each bisection midpoint propagation;
- refined jump-time propagation;
- post-jump remainder propagation.

`non_precomputed_step_count` should count the subset of these calls that use
`propagate_blocks(...)` rather than `propagate_blocks_with_propagators(...)`.

These counters are runtime diagnostics only. They should not change the
trajectory state-update rules.

= Method in Pseudo-Code

```python
def simulate_single_trajectory(...):
    validate inputs
    build one trajectory-local RNG
    load or build precomputed data
    build t_eval and initial psi_blocks
    save snapshot at t = 0

    threshold = rng.random()
    for each phase:
        skip zero-duration phases
        load phase-local jump operators, generators, propagators
        while current_time has not reached this phase end:
            step = min(dt, phase boundary, next t_eval boundary)
            trial = precomputed full-step path or variable-step path
            if total_norm2_list(trial) > threshold:
                accept no-jump step
                maybe save snapshot
                continue

            refine one jump time by five-step bisection
            propagate to refined jump time
            renormalize
            apply jump
            record jump time
            draw new threshold
            maybe save snapshot

            if remainder exists:
                attempt one remainder propagation
                if it stays above threshold:
                    accept it

    renormalize final state
    ensure all t_eval snapshots were saved
    derive omega_2 when needed
    return TrajectoryResult(...)
```

= Invariants and Edge Cases

- The sector ordering in `precomputed["sector_list"]` should define the ordering
  used by `psi_blocks`, all phase-resolved operator lists, saved snapshots, and
  final output dictionaries.
- Snapshot times should stay identical across trajectories because ensemble
  post-processing assumes aligned `t_eval`.
- The solver should split steps to hit `t_eval` and phase boundaries exactly,
  rather than saving by interpolation afterward.
- The single-trajectory code should not rebuild sparse operators inside the hot
  loop when `precomputed` data are available.
- The shifted-jump picture should still validate `Gamma > 0` because the jump
  operator contains `omega / Gamma`.
- The saved `phase_index` should correspond to the phase active when the
  snapshot was written. The initial snapshot should use `0`, and the final
  fallback snapshot should use `max(len(phases) - 1, 0)` as in the current
  code.
