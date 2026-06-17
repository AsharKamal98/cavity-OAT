#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Simulation Precompute: Implementation Instructions]
]

= Scope

This file describes the intended implementation of the precompute layer used by
the custom MCWF simulator. The precompute layer should not change the physical
model. It should only cache objects that are identical across trajectories and
repeated full-`dt` steps.

= Method

The simulator evolves sector blocks under a piecewise-constant protocol. Within
one phase, the following quantities are fixed:

- the drive $Omega$;
- the detuning $delta$;
- the decay rate $Gamma$;
- the reduced sector basis;
- the base internal timestep `dt`.

Therefore, for each phase and each sector, the code should build once:

$
l_(alpha,s), quad
H_("eff",alpha,s), quad
U_(alpha,s)(d t) =
exp[-i H_("eff",alpha,s) d t],
$

where $alpha$ labels the phase and $s$ labels the sector key.

For a given sector $s$, let $A_s$ denote the sector lowering operator used by
the collapse channel. In homogeneous sectors 
$
A_s = J_s^-.
$
In inhomogeneous
sectors $A_s$ is the weighted group lowering operator stored as
`ops.A_weighted`
$
A_s = omega_1 J_(1,s)^- + omega_2 J_(2,s)^-.
$
See `docs/inststructions/paper_inhomogeneous_couplings.typ` for further details about the inhomogeneous jump operator. 

The phase-resolved jump operator should therefore look like
the following.

Regular jump picture:
$
l_(alpha,s) = A_s.
$

Shifted jump picture:

$
l_(alpha,s) = A_s + i Omega_alpha bb(1)_s / Gamma.
$

Here $bb(1)_s$ is the identity matrix on sector $s$. The corresponding
non-Hermitian generator should look like the following.

Regular jump picture:

$
H_("eff",alpha,s)
=
Omega_alpha J_(x,"drive",s)
- delta_alpha N_(e,s)
- (i Gamma)/2 A_s^dagger A_s.
$

Shifted jump picture:

$
H_("eff",alpha,s)
=
-delta_alpha N_(e,s)
- (i Gamma)/2 l_(alpha,s)^dagger l_(alpha,s).
$

These matrices act only inside the reduced basis of sector $s$: dimension
$N_J + 1$ for homogeneous sectors, or $(N_(J,1)+1)(N_(J,2)+1)$ for two-group
inhomogeneous sectors.

These operator forms should be read as the sector-resolved implementation of
the master-equation Hamiltonian and jump operators in `docs/theory/main.tex`.

= Data In

The precompute entry point is:

```python
build_precomputed_trajectory_data(
    N,
    Gamma,
    phases,
    sector_coeffs,
    dt,
    shifted_jump_operator=False,
    omega_1=None,
    N1=None,
    N2=None,
) -> dict
```

Homogeneous sector keys should be integers:

```python
sector_coeffs = {Nj: coeff}
```

Inhomogeneous two-group sector keys should be tuples:

```python
sector_coeffs = {(Nj1, Nj2): coeff}
```
Tuple sector keys require `N1` and `N2`. The function should
validate that $N_1 + N_2 = N$. Tuple sector keys  also require `omega_1`, from which the group-2 coupling should be fixed by `omega2_from_weighted_average(...)` in
`quantum_trajectories/operator_helpers.py`; it should not be recomputed per sector.

= Data Out

The returned dictionary should have one shared sector ordering:

```python
precomputed = {
    "sector_list": [...],
    "ops_list": [...],
    "multiplicities": {...},
    "dims": {...},
    "phase_jump_operators": [[...], ...],
    "phase_generators": [[...], ...],
    "phase_propagators": [[...], ...],
}
```

The list-valued objects should all use the same sector order as `sector_list`.
The hot simulation loop should use these aligned lists instead of repeatedly
indexing dictionaries.

= Construction Flow

The implementation should follow this high-level flow:

```python
sector_list = sorted(sector_coeffs.keys(), key=split_sector_key)

ops_list = [
    build_sector_ops_for_key(key, omega_1=omega_1, N1=N1, N2=N2)
    for key in sector_list
]

multiplicities = {
    key: sector_multiplicity(...) or two_group_sector_multiplicity(...)
    for key in sector_list
}

dims = {key: ops.Jm.shape[0] for key, ops in zip(sector_list, ops_list)}

phase_jump_operators = [
    [build_phase_jump_operator_for_sector(ops, phase.omega, Gamma, ...)]
    for phase in phases
]

phase_generators = [
    [heff_for_sector(ops, phase.omega, phase.delta, Gamma, jump_operator=l, ...)]
    for phase, l in zip(phases, phase_jump_operators)
]

phase_propagators = [
    [expm((-1j * H_eff) * dt).tocsc() for H_eff in generators_for_phase]
    for generators_for_phase in phase_generators
]
```

= Sector Operators

Sector operators should be constructed through:

```python
build_sector_ops_for_key(sector_key, omega_1=omega_1, N1=N1, N2=N2)
```

For homogeneous integer keys, this returns the usual single-sector operators on
the `|n_e>` basis. For inhomogeneous tuple keys, this returns product
Dicke-basis operators on `|n_e1, n_e2>`.

The simulation code should prefer the generic fields on `SectorOperators`:

```python
ops.A_weighted
ops.AdagA_weighted
ops.J_x_drive
```

For homogeneous sectors these fields are aliases of the ordinary unweighted
operators. For inhomogeneous sectors they contain the weighted group operators.
This convention keeps the propagation code mostly independent of whether the
sector key is homogeneous or group-resolved.

= Jump Operators

For each sector and phase, the unshifted jump operator should be:

$
l = A,
$

where $A$ is represented by `ops.A_weighted`. In homogeneous runs this is
$J_-$. In inhomogeneous runs this is the weighted collective lowering operator.

If `shifted_jump_operator=False`, the returned jump operator should be the
unshifted operator.

If `shifted_jump_operator=True`, the returned jump operator should be:

$
l = A + i frac(Omega, Gamma) bb(1).
$

The shifted jump operator requires $Gamma > 0$. The implementation should fail
with a clear `ValueError` if this condition is not satisfied.

= Effective Generators

For the regular jump picture, the effective generator should be:

$
H_("eff")
=
Omega J_(x,"drive")
- delta N_e
- frac(i Gamma, 2) A^dagger A.
$

In code this corresponds to:

```python
H = omega * ops.J_x_drive - delta * ops.N_e
H_eff = H - 0.5j * Gamma * ops.AdagA_weighted
```

For the shifted jump picture, the current custom MCWF implementation should use:

$
H_("eff")
=
-delta N_e
- frac(i Gamma, 2) l^dagger l,
$

where $l = A + i Omega bb(1) / Gamma$. The explicit drive term should not be
added separately in this shifted picture, because it is already represented
through the shifted collapse operator convention used by the simulator.

= Fast Full-Step Path

A propagation step should use precomputed propagators only when the actual step
length is exactly the base `dt`, up to the small tolerance used by the
simulation loop:

```python
if abs(step - dt) <= 1e-15:
    psi_blocks = propagate_blocks_with_propagators(
        psi_blocks,
        full_step_propagators,
    )
```

This path applies:

$
psi_s(t + d t) = U_s(d t) psi_s(t)
$

for each sector block $s$. It avoids recomputing a matrix exponential during
ordinary full-`dt` evolution.

= Variable-Step Path

Precomputed full-`dt` propagators should not be used when the required step
length differs from `dt`. In that case the simulator should call:

```python
propagate_blocks(psi_blocks, generators_list, step)
```

which uses `expm_multiply` sector by sector:

$
psi_s(t + Delta t)
=
exp[-i H_("eff",s) Delta t] psi_s(t),
quad
Delta t != d t.
$

The variable-step path should be used for:

- phase-boundary steps;
- `t_eval` boundary steps;
- jump-time bisection midpoint propagations;
- propagation to the refined jump time;
- propagation of the post-jump remainder.

= Ensemble Reuse

`run_trajectory_ensemble(...)` should call
`build_precomputed_trajectory_data(...)` once per ensemble run and pass the
resulting dictionary to each `simulate_single_trajectory(...)` call.

In multiprocessing mode, `_init_trajectory_worker(...)` should store the
precomputed dictionary in worker-local global state. Individual worker tasks
should receive only the trajectory seed sequence and should reuse the
worker-local precomputed objects.

= Step Counters

The simulator should track:

```python
total_step_count
non_precomputed_step_count
```

`total_step_count` should count actual propagation calls, including ordinary
full-`dt` calls, partial steps, bisection midpoint propagations, refined
jump-time propagation, and post-jump remainder propagation.

`non_precomputed_step_count` should count propagation calls that use
`propagate_blocks(...)` instead of `propagate_blocks_with_propagators(...)`.

These counters should be runtime diagnostics only. They should not change
trajectory physics.

= Invariants

- Precomputation should happen outside the per-trajectory hot loop whenever the
  data are shared by all trajectories.
- All phase- and sector-resolved lists should remain aligned with `sector_list`.
- Full-`dt` propagators should only be used for exact full `dt` steps.
- Variable $Delta t$ propagation should use `expm_multiply`.
- Shifted and unshifted jump pictures should use the generator conventions
  listed above.
- Inhomogeneous tuple sectors should require `omega_1`, `N1`, and `N2`.
