#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Sector Operators: Implementation Instructions]
]

= Purpose

This file specifies how reduced-basis sector operators should be constructed in
`quantum_trajectories/operator_helpers.py`. Use it for changes involving
`build_sector_ops_for_key(...)`, `build_sector_ops(...)`, or
`build_two_group_sector_ops(...)`.

The precompute layer should collect these objects into `ops_list`; see
`docs/instructions/simulation_precompute.typ` for how `ops_list` is used by
phase generators, jump operators, and full-`dt` propagators.

= Function Structure

The public dispatch helper should be:

```python
build_sector_ops_for_key(key, omega_1=None, omega_2=None, N1=None, N2=None)
    -> SectorOperators
```

For homogeneous integer keys, `omega_1`, `omega_2`, `N1`, and `N2` are nullable
and the function should call:

```python
build_sector_ops(Nj) -> SectorOperators
```

For inhomogeneous tuple keys, `omega_1`, `omega_2`, `N1`, and `N2` are required
and the function should call:

```python
build_two_group_sector_ops(Nj1, Nj2, omega_1, omega_2, N1, N2)
    -> SectorOperators
```

= Homogeneous Sector Operators

For a homogeneous integer sector key $N_(J)$, `build_sector_ops(Nj)` should
construct collective two-level operators on the symmetric `|n_e>` basis, with
$n_e = 0, ..., N_(J)$.

The returned `SectorOperators` object should include:

- `Jp`, `Jm`, `JpJm`, `J_x`, `J_y`, `J_z`, and `N_e`;
- `J_drive = J_x`;
- `A_weighted = Jm`;
- `AdagA_weighted = JpJm`;
- single-group fields such as `J_x_groups=(J_x,)` and `N_e_groups=(N_e,)`.

= Two-Group Sector Operators

For an inhomogeneous tuple sector key $(N_(J,1), N_(J,2))$,
`build_two_group_sector_ops(...)` should construct product Dicke-basis
operators on `|n_e1, n_e2>`.

The unweighted full-system observables should be:

```python
Jp = J1p + J2p
Jm = J1m + J2m
N_e = N_e1 + N_e2
```

The unweighted group-resolved fields should include:

```python
J_x_groups = (J1_x, J2_x)
J_y_groups = (J1_y, J2_y)
J_z_groups = (J1_z, J2_z)
N_e_groups = (N_e1, N_e2)
```

The drive and jump operators used for propagation should be:

```python
J_drive = omega_1 * J_x1 + omega_2 * J_x2
A_weighted = omega_1 * J1m + omega_2 * J2m
```

The corresponding decay operator should be stored as:

```python
AdagA_weighted = A_weighted.conjugate().transpose() @ A_weighted
```

= SectorOperators Contract

Both homogeneous and two-group builders should return a `SectorOperators`
object. In pseudo-code, `ops` means this object for the current sector.

The generic fields should be used by downstream code:

```python
ops.A_weighted
ops.AdagA_weighted
ops.J_drive
```

Code that computes physical observables should use the unweighted operators
`J_x`, `J_y`, `J_z`, and `N_e`. Code that computes group-resolved diagnostics
should use the corresponding `*_groups` fields.

= Caching

`build_sector_ops(...)` and `build_two_group_sector_ops(...)` should be cached
with `functools.lru_cache(maxsize=None)`.

The cache key is the exact function argument tuple. For example, the first call
to

```python
build_two_group_sector_ops(Nj1, Nj2, omega_1, omega_2, N1, N2)
```

constructs the sparse Kronecker-product operators, while later calls with the
same arguments return the already-built `SectorOperators` object.

This cache avoids repeated sparse-matrix construction during post-processing and
diagnostics. The cache is process-local, so multiprocessing workers build and
hold their own cached operator objects.
