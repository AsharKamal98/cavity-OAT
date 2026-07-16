#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Single-Trajectory Simulation: Implementation Instructions]
]

= Purpose

This file describes `_simulate_single_trajectory(...)` in the custom MCWF
solver, currently implemented in `solvers/mcwf/sim.py`. Use it when editing
trajectory-local evolution, jump handling, snapshot saving, or the
`TrajectoryResult` output contract.

The single-trajectory layer should consume already-validated inputs and
precomputed sector data from `run_trajectory_ensemble(...)`. Ensemble seeding
and multiprocessing are defined in
`docs/instructions/solvers/mcwf/ensemble_simulation.typ`.

= Method

For one trajectory, the sector-block wavefunction evolves between jumps under
the precomputed non-Hermitian generators:

$
tilde(psi)_s(t + Delta t)
=
exp[-i H_("eff",s) Delta t] psi_s(t).
$

The survival probability for an attempted step is the norm of the propagated
unnormalized state:

$
P_("surv")(Delta t) = sum_s ||tilde(psi)_s(t + Delta t)||^2.
$

If $P_("surv")$ stays above the trajectory threshold, the step is accepted as a
no-jump step. If it crosses the threshold, the code refines one jump time
inside that attempted step, applies the phase-local jump operator, records the
jump, and draws a new threshold.

The current jump-time refinement uses a fixed ten-step bisection inside the
attempted step.

The effective generators, jump operators, and full-`dt` propagators are defined
in `docs/instructions/solvers/mcwf/simulation_precompute.typ`. Initial
sector-block conventions are defined in
`docs/instructions/solvers/mcwf/initial_sector_state.typ`.

= Method in Pseudo-Code

```python
def _simulate_single_trajectory(..., seed_sequence, precomputed):
    rng = np.random.default_rng(seed_sequence)
    psi_blocks = build initial sector blocks in precomputed sector order
    save snapshot at t = 0

    threshold = rng.random()
    for each nonzero-duration integration_phase:
        load integration-Phase-local operators, generators, and propagators
        while current_time is before the integration Phase end:
            step = min(dt, phase boundary, next t_eval boundary)
            trial = full-step propagator path if step == dt else variable-step path

            if total_norm2_list(trial) > threshold:
                accept no-jump step
                maybe save snapshot
                continue

            refine jump time by ten-step bisection
            propagate to the refined jump time
            renormalize, apply jump, renormalize
            record jump time and draw a new threshold
            maybe save snapshot

            attempt the post-jump remainder only if it does not already cross
            the new threshold

    renormalize final state
    ensure every requested t_eval snapshot was saved
    return TrajectoryResult(...)
```

`_simulate_single_trajectory(...)` should keep the trajectory physics local to
one child seed. It should reuse the supplied precomputed objects and should not
rebuild sector operators, phase generators, or full-step propagators.

= Data Requirements

The private entry point should receive the completed model inputs and the shared
precompute object from the ensemble layer:

```python
_simulate_single_trajectory(
    Ni,
    omega_i,
    Gamma,
    integration_phases,
    sector_coeffs,
    *,
    dt,
    t_eval,
    seed_sequence,
    shifted_jump_operator=False,
    precomputed,
) -> TrajectoryResult
```

The ensemble layer should validate public inputs before calling this helper.
This function may use simple internal assertions for assumptions such as
`len(omega_i) = len(Ni)` and `precomputed is not None`.

= Output

The returned object should be:

```python
TrajectoryResult(
    final_sector_blocks=...,
    snapshots=...,
    jump_times=...,
    jump_count=len(jump_times),
    total_step_count=...,
    non_precomputed_step_count=...,
)
```

Shared run metadata such as `Ni`, `omega_i`, `Gamma`, `phase_protocol`,
`t_eval`, and sector dimensions should live on `TrajectoryEnsemble.metadata`,
not on each `TrajectoryResult`.

= Invariants and Edge Cases

- The sector order in `precomputed["sector_list"]` should define the order of
  `psi_blocks`, phase-local operator lists, snapshots, and final output blocks.
- Save exactly one snapshot per requested `t_eval` point, with the first
  snapshot at `t=0`.
- Each snapshot stores the corresponding `integration_phase_index`.
- Split steps so the solver lands on phase boundaries and saved-time points
  exactly, rather than interpolating snapshots afterward.
- Saved no-jump `sector_blocks` may be unnormalized; snapshot `norm` stores
  the corresponding block norm for later expectation-value normalization.
- Refine jump times with the current fixed ten-step bisection rule unless the
  user explicitly requests a behavior change.
- Resolve at most one explicitly refined jump per attempted outer-loop step;
  if the post-jump remainder would cross the new threshold, leave further
  evolution to the next loop iteration.
- Track `total_step_count` and `non_precomputed_step_count` as runtime
  diagnostics only; they must not affect trajectory physics.
- If all requested `t_eval` snapshots are not saved by the end, raise an error
  instead of returning a partially aligned trajectory.
