#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#let ket(x) = math.equation(block: false, $|#x chevron.r$)
#let bra(x) = math.equation(block: false, $chevron.l #x|$)

#align(center)[
  #text(size: 1.6em, weight: "bold")[Generalized Three-Level Squeezing Parameter]
  ]

= Purpose
This instruction specifies how to construct the generalized three-level squeezing parameter. The goal is to compute $xi^2_("gen")(t)$ at each saved timestep after the simulation has run.

The squeezing calculation should be implemented as a standalone post-processing function that takes the simulation result or observable output as input.

= Theory

== Per-timestep calculation

At each timestep, do the following.

=== Define the effective dressed state $|1⟩$
Step 1: Identify the $theta_J, phi_J$ effective qubit angles.


Work in the single-particle basis
$
(|u⟩, |d⟩, |e⟩).
$

Use the ansatz

$
|1⟩ = 0 |u⟩ + cos(theta_J / 2) |d⟩ + e^(-i phi_J) sin(theta_J / 2) |e⟩.
$

Here $theta_J,phi_J$ should be found by comparing this ansatz to the simulation state at the current timestep. Follow the active-manifold Bloch-vector averaging convention in:

```
docs/instructions/bloch_vector_averaging.typ
```

=== Define the instantaneous mean direction $|c⟩$
Step 2: Using the $J$ angles from step 1, identify the $theta_S, phi_S$ effective $S$-Bloch angles.

Use the ansatz
$
|c⟩ = cos(theta_S / 2) |u⟩ + e^(-i phi_S) sin(theta_S / 2) |1⟩.
$

Here $theta_S,phi_S$ should be found by comparing this ansatz to the instantaneous effective $S$-Bloch vector of the simulation state. The angles $theta_J,phi_J$ appearing inside $|1⟩$ are the angles found in step 1.

In the explicit $(u,d,e)$ basis,

$
|c⟩ = mat(cos(theta_S/2); e^(-i phi_S) sin(theta_S/2) cos(theta_J/2); e^(-i(phi_S+phi_J)) sin(theta_S/2) sin(theta_J/2)).
$

=== Define the $J$-fluctuation direction $|j⟩$
Step 3: Using the $J$ angles from step 1and $S$ angles from step 2, construct the state orthogonal to $|c⟩$ on the $J$-sphere ($(d,e)$ manifold):

$
|j⟩ = 0 |u⟩ - sin(theta_J / 2) |d⟩ + e^(-i phi_J) cos(theta_J / 2) |e⟩.
$

In vector form:

$
|j⟩ = mat(0; -sin(theta_J/2); e^(-i phi_J) cos(theta_J/2)).
$

=== Define the $S$-fluctuation direction $|s⟩$

Step 4: Using the $J$-angles from step 1 and the $S$-angles from step 2, construct the state orthogonal to $|c⟩$ on the $S$-sphere ($(u,1)$ manifold):

$
|s⟩ = -sin(theta_S / 2) |u⟩ + e^(-i phi_S) cos(theta_S / 2) |1⟩.
$

In the $(u,d,e)$ basis,

$
|s⟩ = mat(-sin(theta_S/2); e^(-i phi_S) cos(theta_S/2) cos(theta_J/2); e^(-i(phi_S+phi_J)) cos(theta_S/2) sin(theta_J/2)).
$

The three states $|c⟩,|j⟩,|s⟩$ should form an orthonormal single-particle basis up to numerical precision.

=== Construct the four local fluctuation operators

Construct the four $3 times 3$ single-particle operators:

$
o_1 = frac(ket(c) bra(j) + ket(j) bra(c), 2),
$

$
o_2 = frac(ket(c) bra(j) - ket(j) bra(c), 2 i),
$

$
o_3 = frac(ket(c) bra(s) + ket(s) bra(c), 2),
$

$
o_4 = frac(ket(c) bra(s) - ket(s) bra(c), 2 i).
$

Then construct the corresponding collective operators

$
O_a = sum_i o_a^((i)), quad a=1,2,3,4.
$

=== Construct the covariance matrix

For the current state $|psi(t)⟩$, compute

$
mu_a = ⟨ O_a⟩.
$

Then construct the $4 times 4$ covariance matrix

$
C_(a,b) = frac(1, 2) ⟨ O_a O_b+O_b O_a⟩ - mu_a mu_b.
$

=== Minimum fluctuation direction

Diagonalize the $4 times 4$ covariance matrix and take

$
lambda_("min")(C).
$

This is the minimum generalized transverse variance.

=== Generalized squeezing parameter

Compute
$
xi^2_("gen")(t) = frac(N lambda_("min")(C), ⟨ N_c/2⟩^2),
$

where
$
N_c = sum_i ket(c_i) bra(c_i).
$

Compute $⟨ N_c⟩$ using the same collective-operator construction:
$
N_c = sum_(mu,nu) (ket(c) bra(c))_(mu nu) A_(mu nu).
$
If the state is well polarized along $|c⟩$, then $⟨ N_c/2⟩ approx N/2$, but the code should compute it explicitly rather than assuming this.

== Inhomogeneous Coupling Extension

This section should be used for two-group inhomogeneous simulations whose
sector keys are $(N_(J,1),N_(J,2))$. The homogeneous squeezing instruction above
remains the base algorithm. The inhomogeneous extension should only change how
the local dressed bases and fluctuation operators are grouped and combined.

=== Per-group squeezing

For each group $g in {1,2}$, compute a group-local squeezing parameter by
applying the homogeneous squeezing construction to that subgroup only.

For group $g$, first compute group-resolved active-manifold moments and obtain
one group-local effective qubit direction,

$
theta_(J,g)(t), quad phi_(J,g)(t).
$

Then compute the group-local effective $S$-Bloch direction,

$
theta_(S,g)(t), quad phi_(S,g)(t),
$

and construct the corresponding group-local states

$
|1_g⟩, quad |c_g⟩, quad |j_g⟩, quad |s_g⟩.
$

The group-local fluctuation operators should be

$
O_(a,g) = sum_(i in g) o_(a,g)^((i)), quad a=1,2,3,4,
$

where $o_(a,g)$ is constructed from $|c_g⟩,|j_g⟩,|s_g⟩$ in the same way as in
the homogeneous instruction.

For each group, construct

$
C_(a,b,g)
= frac(1,2) ⟨ O_(a,g) O_(b,g) + O_(b,g) O_(a,g) ⟩
- ⟨O_(a,g)⟩ ⟨O_(b,g)⟩,
$

and compute

$
xi_g^2(t) =
frac(N_g lambda_("min")(C_g), ⟨N_(c,g) / 2⟩^2),
quad g=1,2.
$

=== Full-system inhomogeneous squeezing

The full inhomogeneous squeezing should keep the group-local dressed bases, but
combine the group fluctuation operators before constructing the covariance
matrix. First compute the separate group operators,

$
O_(a,1), quad O_(a,2),
$

then define the full-system fluctuation operator

$
O_a = O_(a,1) + O_(a,2), quad a=1,2,3,4.
$

The full covariance matrix should then be constructed exactly as in the
homogeneous case:

$
C_(a,b)
= frac(1,2) ⟨ O_a O_b + O_b O_a ⟩
- ⟨O_a⟩ ⟨O_b⟩.
$

This construction keeps cross-group covariance terms automatically, because
$O_a O_b$ contains products of both group-1 and group-2 operators.

The full-system active mean population should be
$
N_c = N_(c,1) + N_(c,2),
$

and the full squeezing is as usual
$
xi^2(t) =
frac(N lambda_("min")(C), ⟨N_c / 2⟩^2).
$

= Data Requirements

These data requirements support both the regular homogeneous squeezing plot and
the inhomogeneous extension plot. The current squeezing and squeezing-plot
helpers should be treated as post-processing on saved snapshots.

Required per-timestep data points are:

- time $t_k$ of snapshot;
- the full sector-block wavefunction at that time `snapshot.sector_blocks`, in
  the reduced homogeneous basis $(N_J,n_e)$ or the reduced two-group basis
  $(N_(J,1),N_(J,2),n_(e,1),n_(e,2))$.

The required trajectory- or ensemble-level metadata is:

- `N`, `Gamma`, and the full phase list;
- the sector-key list;
- inhomogeneous metadata such as `N1`, `N2`, `omega_1`, and `omega_2` when
  tuple sectors are used.

From this saved data, the current post-processing code recomputes the quantities
used by plotting:

- homogeneous plotting:
  - $J$-Bloch and $S$-Bloch angles $theta_J(t)$, $phi_J(t)$, $theta_S(t)$, $phi_S(t)$;
  - ordered covariance eigenvalues $lambda(t)$;
  - mean dressed active population $⟨N_c⟩(t)$;
  - active-manifold excited fraction $frac(⟨N_e⟩(t), ⟨N_J⟩(t))$;
  - Wineland squeezing parameter $xi^2_("gen")(t)$.
- inhomogeneous plotting, with full and group-resolved quantities:
  - $J$-Bloch and $S$-Bloch angles $theta_J(t)$, $phi_J(t)$, $theta_S(t)$, $phi_S(t)$ and $theta_(J,g)(t)$, $phi_(J,g)(t)$, $theta_(S,g)(t)$, $phi_(S,g)(t)$;
  - smallest covariance eigenvalues $lambda_("min")(t)$ and $lambda_g("min")(t)$;
  - mean dressed active populations $⟨N_c⟩(t)$ and $⟨N_(c,g)⟩(t)$;
  - active-manifold excited fractions $frac(⟨N_e⟩(t), ⟨N_J⟩(t))$ and $frac(⟨N_(e,g)⟩(t), ⟨N_(J,g)⟩(t))$;
  - Wineland squeezing parameters $xi^2(t)$ and $xi_g^2(t)$.

Therefore, the implementation should prefer saving the snapshot wavefunction and
metadata once, then recomputing these derived squeezing series during
post-processing rather than storing redundant per-timestep angle or covariance
arrays.

= Pseudo-code

The new-pipeline squeezing code should keep the old pipeline physics order, but
split it into small post-processing functions.

```python
compute_squeezing_reference_angles(result, moments=None, ...)
    -> SqueezingReferenceAngles:
    # 1. Compute averaged J moments on the saved t grid.
    #    Reuse moments.J when available.
    # 2. Convert the averaged J moments into theta_J, phi_J.
    # 3. Expand the S operators into the (u,d,e) basis using theta_J, phi_J.
    # 4. Extract S-Bloch moments per trajectory, average those moments, then
    #    construct theta_S, phi_S.
    # 5. For tuple-sector runs, repeat steps 1-4 for each group and return
    #    group-resolved angle arrays as well as full-system angle arrays.
```

For the S-Bloch extraction, rewrite the $S$ operators into the $(u,d,e)$ basis.
For example, with fixed $theta_J,phi_J$,

$
S_x =
frac(1,2) (ket(1) bra(u) + ket(u) bra(1)),
$

where

$
|1⟩ =
cos(theta_J / 2) |d⟩
+ e^(-i phi_J) sin(theta_J / 2) |e⟩.
$

Therefore

$
ket(1) bra(u) =
cos(theta_J / 2) ket(d) bra(u)
+ e^(-i phi_J) sin(theta_J / 2) ket(e) bra(u).
$

The same construction gives $S_y$, $S_z$, and $N_("u1")$. The resulting
trajectory moments should be averaged before constructing $theta_S,phi_S$.

```python
compute_trajectory_squeezing_moments(trajectory, reference_angles, ...)
    -> SqueezingMomentSample:
    for each saved snapshot:
        build |c>, |j>, |s> from the reference angles
        build local o_a for a = 1,2,3,4
        build O_a on the reduced simulation basis
        v_a = O_a @ psi
        mu_a = <psi | O_a | psi>
        second_ab = Re <v_a | v_b>
        N_c = <psi | N_c | psi>
        excited_fraction_active = <N_e> / <N_J>
    return raw moment arrays
```

Do not construct full $3^N$-dimensional tensor-product operators. Use the
existing reduced sector basis. When possible, construct each collective
fluctuation operator from precomputed collective one-body transition operators:

$
A_(mu nu)=sum_i ket(mu_i) bra(nu_i), quad mu,nu in (u,d,e),
$

$
O_a = sum_(mu,nu) (o_a)_(mu nu) A_(mu nu).
$

Do not explicitly construct $O_a O_b$. Use the vectors $v_a=O_a |psi⟩$ and set

$
"second"_(a,b) = "Re" [⟨ v_a | v_b ⟩].
$

For tuple-sector runs, the group-local angle extraction should supply the local
operators for the relevant group, e.g.
$J_(x,g),J_(y,g),J_(z,g),N_(e,g)$, to the same Bloch-vector averaging logic as
the single-group case. The corresponding group-local fluctuation operators
should act only on the chosen group factor of the product basis.

```python
compute_average_squeezing_moments(samples: list[SqueezingMomentSample])
    -> SqueezingMomentSample:
    average mu_a, second_ab, N_c, and excited_fraction_active over trajectories

finalize_squeezing_series(averaged_moments, reference_angles, ...)
    -> SqueezingSeries:
    for each saved timestep:
        C = averaged_second - outer(averaged_mu, averaged_mu)
        covariance_eigvals = eigvalsh(C)
        lambda_min = min(covariance_eigvals)
        xi2 = N * lambda_min / (N_c / 2)^2
        xi2_db = 10 * log10(xi2) where xi2 > 0
    attach full-system fields
    attach group-resolved fields when present
    return SqueezingSeries

compute_ensemble_squeezing(ensemble, moments=None, ...)
    -> SqueezingSeries:
    reference_angles = compute_squeezing_reference_angles(ensemble, moments)
    samples = map_with_optional_pool(
        compute_trajectory_squeezing_moments(traj, reference_angles)
        for traj in ensemble.trajectories
    )
    averaged = compute_average_squeezing_moments(samples)
    return finalize_squeezing_series(averaged, reference_angles)
```

The plotting entry points should be thin wrappers:

```python
plot_generalized_xi(squeezing: SqueezingSeries, ...):
    make_2x2_grid_plot(
        squeezing_db=squeezing.xi2_db,
        covariance_eigenvalues=squeezing.covariance_eigvals,
        dressed_population=squeezing.N_c,
        excited_fraction=squeezing.excited_fraction_active,
    )

plot_inhomogeneous_generalized_xi(squeezing: SqueezingSeries, ...):
    make_2x2_grid_plot(
        squeezing_db=squeezing.xi2_db and squeezing.xi2_db_groups,
        smallest_covariance_eigenvalue=squeezing.lambda_min and squeezing.lambda_min_groups,
        dressed_population=squeezing.N_c and squeezing.N_c_groups,
        excited_fraction=squeezing.excited_fraction_active and squeezing.excited_fraction_active_groups,
    )
```

= Output

The planned squeezing output should be a Pydantic class:

```python
class SqueezingSeries(BaseModel):
    t: Array

    # full-system squeezing fields
    xi2: Array
    xi2_db: Array
    lambda_min: Array
    covariance_eigvals: Array  # shape (n_t, 4)
    N_c: Array
    excited_fraction_active: Array

    # full-system reference angles
    theta_J: Array
    phi_J: Array
    theta_S: Array
    phi_S: Array

    # optional group-resolved squeezing fields
    xi2_groups: tuple[Array, ...] | None = None
    xi2_db_groups: tuple[Array, ...] | None = None
    lambda_min_groups: tuple[Array, ...] | None = None
    covariance_eigvals_groups: tuple[Array, ...] | None = None
    N_c_groups: tuple[Array, ...] | None = None
    excited_fraction_active_groups: tuple[Array, ...] | None = None

    # optional group-resolved reference angles
    theta_J_groups: tuple[Array, ...] | None = None
    phi_J_groups: tuple[Array, ...] | None = None
    theta_S_groups: tuple[Array, ...] | None = None
    phi_S_groups: tuple[Array, ...] | None = None
```

Internal helper outputs such as `SqueezingReferenceAngles` and
`SqueezingMomentSample` may also be Pydantic classes if they are passed between
public helpers. At minimum, they should make clear whether a field is a raw
trajectory moment, an ensemble-averaged moment, or a nonlinear derived field.

= Invariants

- Squeezing should be computed as post-processing on saved snapshots; it should
  not change the MCWF propagation path.
- For ensemble results, average the required raw operator moments over
  trajectories before constructing the covariance matrix and $xi^2$.
- Do not compute one covariance matrix or one $xi^2$ per trajectory and then
  average those final nonlinear quantities.
- Homogeneous squeezing should define one effective qubit for the full saved
  wavefunction or full ensemble-averaged state at each timestep. Do not select
  only the central $N_J=N/2$ sector and do not define one separate effective
  qubit per $N_J$ sector.
- Inhomogeneous full-system squeezing should compute full and group-resolved
  quantities separately. The full-system covariance should use
  $O_a = O_(a,1) + O_(a,2)$ so cross-group covariance terms are retained.
