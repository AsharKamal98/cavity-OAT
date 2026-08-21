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

= Theory

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

Use the normalized single-particle Gell--Mann generators

$
g_a = lambda_a / 2,
quad a in {1, dots.c, 8}.
$

For physical group $q$, define the corresponding collective generator by

$
hat(G)_(q,a) = sum_(i in q) g_a^((i)).
$

For a one-group run, omit the group label and write $hat(G)_a$. At each saved
time, the normalized first-order moment for trajectory $r$ is

$
mu_(q,a)^((r))(t)
=
frac(
  ⟨psi_r(t) | hat(G)_(q,a) | psi_r(t)⟩,
  ⟨psi_r(t) | psi_r(t)⟩
).
$

The normalized symmetric second-order moment is

$
M_(q a, h b)^((r))(t)
=
frac(
  ⟨psi_r(t) |
  {hat(G)_(q,a), hat(G)_(h,b)}
  | psi_r(t)⟩,
  2 ⟨psi_r(t) | psi_r(t)⟩
).
$

Here $q,h$ run over all physical groups and $a,b in {1, dots.c, 8}$. Together,
the complete first- and second-order moment sets contain enough information to
compute the observables in the post-analysis pipeline without retaining the
wavefunction.

= Method in Pseudo-Code

```python
def _simulate_single_trajectory(..., seed_sequence, precomputed):
    rng = np.random.default_rng(seed_sequence)
    psi_blocks = build initial sector blocks in precomputed sector order
    compute and save one raw generator-moment row at t = 0

    threshold = rng.random()
    for each nonzero-duration integration_phase:
        load integration-Phase-local operators, generators, and propagators
        while current_time is before the integration Phase end:
            step = min(dt, phase boundary, next t_eval boundary)
            trial = full-step propagator path if step == dt else variable-step path

            if total_norm2_list(trial) > threshold:
                accept no-jump step
                maybe compute and save a generator-moment row
                continue

            refine jump time with four bisection iterations
            propagate to the refined jump time
            renormalize, apply jump, renormalize
            record jump time and draw a new threshold
            maybe compute and save a generator-moment row

            attempt the post-jump remainder only if it does not already cross
            the new threshold

    ensure every requested t_eval moment row was saved
    return TrajectoryResult(moments=..., diagnostics=...)
```

The current jump-time refinement uses four bisection iterations inside the
attempted step. This reduces an attempted interval `step` to `step / 16`, which
is at most `dt / 16`.

The effective generators, jump operators, and full-`dt` propagators are defined
in `docs/instructions/solvers/mcwf/simulation_precompute.typ`. Initial
sector-block conventions are defined in
`docs/instructions/solvers/mcwf/initial_sector_state.typ`.

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
    shifted_jump_operator=True,
    precomputed,
) -> TrajectoryResult
```

The shifted jump-operator picture is the default. Pass
`shifted_jump_operator=False` explicitly to use the unshifted picture.

The ensemble layer should validate public inputs before calling this helper.
This function may use simple internal assertions for assumptions such as
`len(omega_i) = len(Ni)` and `precomputed is not None`.

= Output

The returned object should be:

```python
TrajectoryResult(
    moments=TrajectoryResolvedMomentSeries(
        t=...,                       # shape (n_t,)
        first_moments=...,           # shape (1, n_t, G, 8)
        second_moments=...,          # shape (1, n_t, K), K=(8G)(8G+1)/2
    ),
    diagnostics=MCWFTrajectoryDiagnostics(
        jump_times=...,
        jump_count=len(jump_times),
        total_step_count=...,
        non_precomputed_step_count=...,
    ),
)
```

Shared run metadata such as `Ni`, `omega_i`, `Gamma`, `phase_protocol`,
`t_eval`, and sector dimensions should live on `TrajectoryEnsemble.metadata`,
not on each `TrajectoryResult`.

= Invariants and Edge Cases

- The sector order in `precomputed["sector_list"]` should define the order of
  `psi_blocks` and phase-local operator lists throughout propagation.
- Save exactly one raw moment row per requested `t_eval` point, with the first
  row at `t=0`.
- Retain the explicit trajectory axis in `TrajectoryResolvedMomentSeries`; it
  has length one for this single-trajectory function.
- Split steps so the solver lands on phase boundaries and saved-time points
  exactly, rather than interpolating moments afterward.
- Normalize every saved first and second moment by the complete sector-block
  norm at that saved time.
- Do not retain saved sector-block wavefunctions in `TrajectoryResult`.
- Refine jump times with the current fixed four-iteration bisection rule unless
  the user explicitly requests a behavior change.
- Resolve at most one explicitly refined jump per attempted outer-loop step;
  if the post-jump remainder would cross the new threshold, leave further
  evolution to the next loop iteration.
- Track `total_step_count` and `non_precomputed_step_count` as runtime
  diagnostics only; they must not affect trajectory physics.
- Require `jump_count` to equal the number of stored `jump_times`.
- If all requested `t_eval` moment rows are not saved by the end, raise an
  error instead of returning a partially aligned trajectory.
