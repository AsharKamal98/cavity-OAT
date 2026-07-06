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
sectors $A_s$ is the weighted group lowering operator is
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
    Ni,
    omega_i,
    Gamma,
    phases,
    sector_coeffs,
    dt,
    shifted_jump_operator=False,
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
`Ni` should contain the group sizes and `omega_i` should contain the completed
group couplings. For homogeneous runs, use `Ni=[N]` and `omega_i=[1.0]`. For
two-group tuple sectors, use `Ni=[N1, N2]` together with the completed
two-entry `omega_i` list. The precompute layer should receive this completed
list and should not reconstruct couplings per sector.

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

The keys have the following meanings (details further down):

- `sector_list` is the sorted list of populated sector keys, defines the
  canonical sector ordering used everywhere else in the precomputed dictionary. 
- `ops_list` contains the `SectorOperators` object for each sector in
  `sector_list`; each object packages the reduced-basis operators needed for
  propagation and observables in that sector.    
- `multiplicities` maps each sector key to its combinatorial multiplicity in the
  full atom ensemble.
- `dims` maps each sector key (`NJ` in homogeneous case, `(NJ1,NJ2)` in two-group inhomogeneous case) to the reduced Hilbert-space dimension of that
  sector.
- `phase_jump_operators` stores the jump operator for every phase and sector,
  indexed as `[phase_index][sector_index]`.
- `phase_generators` stores the non-Hermitian effective generator
  `H_eff` for every phase and sector, using the same two-level indexing.
- `phase_propagators` stores the precomputed full-`dt` propagator
  `exp(-i H_eff dt)` for every phase and sector, used only on exact full-`dt`
  steps. 


== Sector List (`sector_list`)
Sector list should be constructed through:
```Python sector_list = sorted(sector_coeffs.keys(), key=lambda key: split_sector_key(key))
```
For homogeneous cases, sorts first by $N_J$ and then by $n_e$. For two-group inhomogeneous sectors, the keys are tuples and the sorting is first by the sum $N_(J,1) + N_(J,2)$ and then by $N_(J,1)$.

E.g. in the $N=4$ homogeneous case, `sector_list` is 
$
[0, 1, 2, 3, 4];
$
in the two-group $N_1=N_2=1$ inhomogeneous case it is
$
[(0,0), (0,1), (1,0), (1,1)].
$



== Sector Operators (`ops_list`)

Detailed construction rules for `build_sector_ops_for_key(...)`,
`build_sector_ops(...)`, and `build_two_group_sector_ops(...)` are in
`docs/instructions/sector_operators.typ`. Below is a summary focusing mostly on the operator forms.

Sector operators should be constructed through:

```python ops_list = [
    build_sector_ops_for_key(key, Ni=Ni, omega_i=omega_i)
    for key in sector_list
]
```

For homogeneous integer keys, `Ni` and `omega_i` should be one-entry lists and
the function should return the usual single-sector operators on the `|n_e>`
basis. For inhomogeneous tuple keys, the current implementation expects
two-entry `Ni` and `omega_i` lists and returns product Dicke-basis operators on
`|n_e1, n_e2>`.

In both cases, the function returns a `SectorOperators` object for a given sector key. Each object is one element of `ops_list`. In pseudo-code, `ops` means the `SectorOperators` object for the current sector. This object should contain the following fields:
- `Jp`, `Jm`, `JpJm`, `J_x`, `J_y`, `N_e` --- unweighted sector operators, used for physical observables.

  $
  "inhomogeneous:" quad J^+ = J_1^+ + J_2^+, quad
  N_e = N_(e,1) + N_(e,2).
  $
- `N_e_groups`, `J_x_groups`, `J_y_groups` --- unweighted group-resolved operators, used for
    group-resolved observables.

  $
  "inhomogeneous:" quad J_x^("groups") = (J_(x,1), J_(x,2)), quad
  N_e^("groups") = (N_(e,1), N_(e,2)).
  $

- `J_drive`, `A_weighted`, `AdagA_weighted` --- drive/jump operators, used for propagation.

  $
  "homogeneous:" quad J_(x,"drive") = J_x, quad A = J^-,
  $

  $
  "inhomogeneous:" quad
  J_(x,"drive") = omega_1 J_(x,1) + omega_2 J_(x,2), quad
  A = omega_1 J_1^- + omega_2 J_2^-.
  $

The simulation code should use the generic fields on `SectorOperators`:

```python
ops.A_weighted
ops.AdagA_weighted
ops.J_drive
```

For homogeneous sectors these fields are aliases of the ordinary unweighted
operators. For inhomogeneous sectors they contain the weighted group operators.
This convention keeps the propagation code mostly independent of whether the
sector key is homogeneous or group-resolved.

== Sector Multiplicities (`multiplicities`)
The multiplicity of each sector should count how many atom-label configurations
are represented by the reduced sector.

E.g. in the homogeneous $N=4$ and $N_J=2$ case, the reduced
sector represents the atom-label configurations
$(1,1,0,0)$, $(1,0,1,0)$,
..., $(0,0,1,1)$. Thus the multiplicity is $binom(N, N_J)$ in homogeneous
cases and $binom(N_1, N_(J,1)) binom(N_2, N_(J,2))$ in two-group inhomogeneous
cases.


== Sector Dimensions (`dims`)
The reduced Hilbert-space dimension of each sector should be computed through:

```python
dims = {
    sector_key: ops.Jm.shape[0]
    for sector_key, ops in zip(sector_list, ops_list)
}
```

In homogeneous cases it is `Nj + 1`, since for a
given sector `Nj`
$
n_e = 0, 1, ... N_J;
$ 
in two-group inhomogeneous cases it is `(Nj1 + 1) * (Nj2 + 1)`, since for a
given sector `(Nj1, Nj2)`
$
n_(e,1) = 0, 1, ... N_(J,1),quad
n_(e,2) = 0, 1, ... N_(J,2).
$

== Jump Operators (`phase_jump_operators`)
jump operators should be constructed through:
```python phase_jump_operators = [
    [build_phase_jump_operator_for_sector(
        ops,
        phase.omega, 
        Gamma,
        build_phase_jump_operator_for_sector,
    )
    ]
    for phase in phases
]
```

If `shifted_jump_operator=False`, the `omega` and `Gamma` arguments are nullable, and for each sector and phase, the jump operator should be:
$
l = A,
$
where $A$ is the sector specific `ops.A_weighted` inside `ops_list`.

If `shifted_jump_operator=True`, the returned jump operator should be:
$
l = A + i frac(Omega, Gamma) bb(1),
$
where the phase-specific $Omega$ is given by `phase.omega`.



== Effective Generators (`phase_generators`)

Effective generators should be constructed through:

```python
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
```

This produces one generator for every phase and every sector. The indexing
should match the other phase-resolved lists: `phase_generators[phase_index][sector_index]`
is the $H_("eff")$ matrix for that phase and sector.

The effective generator should always be written in the MCWF form

$
H_("eff",s)
=
H_s
- frac(i Gamma, 2) l_s^dagger l_s.
$

The sector lowering operator $A_s$ should be defined by the sector type:

For homogeneous sectors,

$
A_s = J_s^-.
$

For two-group inhomogeneous sectors,

$
A_s = omega_1 J_(1,s)^- + omega_2 J_(2,s)^-.
$

For the regular jump picture,

$
H_s =
Omega J_(x,"drive",s)
- delta N_(e,s),
$

$
l_s = A_s.
$

In code, the regular picture corresponds to:

```python
H = omega * ops.J_drive - delta * ops.N_e
H_eff = H - 0.5j * Gamma * ops.AdagA_weighted
```

For the shifted jump picture, the current custom MCWF implementation should use:

$
H_s =
-delta N_(e,s),
$

$
l_s = A_s + i frac(Omega, Gamma) bb(1)_s.
$

The explicit drive term should not be added separately in this shifted picture,
because it is already represented through the shifted collapse operator
convention used by the simulator.

== Phase Propagators (`phase_propagators`)

Full-step propagators should be constructed from the precomputed effective
generators:

```python
phase_propagators = [
    [
        expm((-1j * H_eff) * dt).tocsc()
        for H_eff in generators_for_phase
    ]
    for generators_for_phase in phase_generators
]
```

This produces one full-`dt` propagator for every phase and every sector, with
the same indexing convention as the other phase-resolved lists:
`phase_propagators[phase_index][sector_index]`.

Mathematically, each stored propagator is

$
U_(alpha,s)(d t)
=
exp[-i H_("eff",alpha,s) d t].
$

= Fast Full-Step Path

A propagation step should use precomputed propagators only when the actual step
length is exactly the base `dt`, up to the small tolerance used by the
simulation loop:

For a given phase, `full_step_propagators` should be the corresponding entry
from `phase_propagators`:

```python
full_step_propagators = phase_propagators[phase_index]
```

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


/*
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
*/
