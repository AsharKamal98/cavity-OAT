# Simulation Precompute

Use this file when changing `build_precomputed_trajectory_data(...)`,
phase-dependent operators, full-step propagators, or the logic deciding whether
a propagation step can use precomputed data.

## Expected Physics

The MCWF simulator evolves sector blocks under a piecewise-constant protocol.
Within a fixed phase, `Omega`, `delta`, `Gamma`, the sector basis, and the base
timestep `dt` are fixed. Therefore the reduced jump operator, non-Hermitian
generator, and full-`dt` propagator should be built once and reused.

The precompute layer should not change the physics. It should only cache objects
that are identical across trajectories and repeated full-`dt` steps.

## Data In

High-level input should look like:

```python
build_precomputed_trajectory_data(
    N=N,
    Gamma=Gamma,
    phases=phases,
    sector_coeffs={sector_key: coeff},
    dt=dt,
    shifted_jump_operator=shifted_jump_operator,
    omega_1=omega_1,  # only for inhomogeneous tuple sector keys
    N1=N1,
    N2=N2,
)
```

Homogeneous sector keys should be integers:

```python
sector_coeffs = {Nj: coeff}
```

Inhomogeneous sector keys should be two-group tuples:

```python
sector_coeffs = {(Nj1, Nj2): coeff}
```

Tuple sector keys require `omega_1`, `N1`, and `N2`. The group-2 coupling should
be derived by the operator helpers from the fixed physical-group normalization.

## Data Out

The returned object should be a dictionary with aligned sector ordering:

```python
precomputed = {
    "sector_list": [sector_key_0, sector_key_1, ...],
    "ops_list": [SectorOperators_0, SectorOperators_1, ...],
    "multiplicities": {sector_key: multiplicity},
    "dims": {sector_key: reduced_dimension},
    "phase_jump_operators": [[l_phase_sector, ...], ...],
    "phase_generators": [[H_eff_phase_sector, ...], ...],
    "phase_propagators": [[U_eff_dt_phase_sector, ...], ...],
}
```

The list-valued entries should use the same sector order as `sector_list`.
Simulation hot loops should use these aligned lists rather than repeatedly
indexing dictionaries.

## Construction Logic

The high-level precompute flow should be:

```python
sector_list = sorted(sector_coeffs.keys(), key=split_sector_key)
ops_list = [build_sector_ops_for_key(sector_key, omega_1, N1, N2) for sector_key in sector_list]
multiplicities = {sector_key: sector_multiplicity_or_two_group_multiplicity(...)}
dims = {sector_key: ops.J_minus.shape[0]}

phase_jump_operators = [
    [
        build_phase_jump_operator_for_sector(
            ops, phase.omega, Gamma, shifted_jump_operator=shifted_jump_operator
        )
        for ops in ops_list
    ]
    for phase in phases
]

phase_generators = [
    [
        heff_for_sector(
            ops,
            phase.omega,
            phase.delta,
            Gamma,
            shifted_jump_operator=shifted_jump_operator,
            jump_operator=jump_operator,
        )
        for ops, jump_operator in zip(ops_list, jump_operators_for_phase)
    ]
    for phase, jump_operators_for_phase in zip(phases, phase_jump_operators)
]

phase_propagators = [
    [expm((-1j * H_eff) * dt).tocsc() for H_eff in phase_generators_for_phase]
    for phase_generators_for_phase in phase_generators
]
```

## Jump Operators

For each phase and sector, the jump operator should be:

```python
unshifted_l = ops.A_weighted if inhomogeneous else ops.J_minus
```

If `shifted_jump_operator=False`, use:

```python
l = unshifted_l
```

If `shifted_jump_operator=True`, use:

```python
l = unshifted_l + 1j * phase.omega / Gamma * identity
```

The shifted jump operator requires `Gamma > 0`.

## Effective Generators

For the regular jump operator, the effective generator should be:

```python
drive_op = ops.J_x_drive if inhomogeneous else ops.J_x
decay_term = ops.AdagA_weighted if inhomogeneous else ops.JpJm
H_eff = phase.omega * drive_op - phase.delta * ops.N_e - 0.5j * Gamma * decay_term
```

For the shifted jump operator, the current implementation should use:

```python
H_eff = -phase.delta * ops.N_e - 0.5j * Gamma * (l.conj().T @ l)
```

This keeps the shifted-jump picture consistent with the current custom MCWF
implementation.

## When Precompute Should Be Used

A step should use `phase_propagators[phase_index]` only when the actual step
length is exactly the base `dt`:

```python
if abs(step - dt) <= tolerance:
    psi_blocks = propagate_blocks_with_propagators(psi_blocks, full_step_propagators)
else:
    psi_blocks = propagate_blocks(psi_blocks, generators_list, step)
```

This is the fast path:

```python
psi_sector(t + dt) = exp(-1j * H_eff * dt) @ psi_sector(t)
```

## When Precompute Should Not Be Used

Precomputed full-step propagators should not be used for any step whose length
is not exactly `dt`. The variable-step path should use `expm_multiply(...)`.

This includes:

1. phase-boundary steps;
2. `t_eval` boundary steps;
3. jump-time bisection midpoint steps;
4. propagation to the refined jump time;
5. propagation of the post-jump remainder.

These cases need:

```python
psi_sector(t + step) = expm_multiply((-1j * H_eff) * step, psi_sector(t))
```

## Ensemble Reuse

`run_trajectory_ensemble(...)` should build this precomputed dictionary once per
ensemble run and pass it to every `simulate_single_trajectory(...)` call.

In multiprocessing mode, worker initialization should store the precomputed data
once per worker process. Individual trajectory tasks should pass only the
trajectory seed sequence and reuse the worker-local precomputed objects.

## Diagnostics

The simulator tracks:

```python
total_step_count
non_precomputed_step_count
```

`total_step_count` should include full-step propagation calls, partial-step
propagation calls, and jump-bisection propagation calls.

`non_precomputed_step_count` should count propagation calls that use
`propagate_blocks(...)` instead of `propagate_blocks_with_propagators(...)`.

These counters are runtime diagnostics only; they should not affect trajectory
physics.
