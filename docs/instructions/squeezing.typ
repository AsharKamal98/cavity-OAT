#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")
#let ket(x) = math.equation(block: false, $|#x chevron.r$)
#let bra(x) = math.equation(block: false, $chevron.l #x|$)

#align(center)[
  #text(size: 1.6em, weight: "bold")[Generalized Three-Level Squeezing Parameter]
  ]

This instruction specifies how to construct the generalized three-level squeezing parameter. The goal is to compute $xi^2_("gen")(t)$ at each saved timestep after the simulation has run.

The squeezing calculation should be implemented as a standalone post-processing function that takes the simulation result or observable output as input.
= Per-timestep calculation

At each timestep, do the following.

== Define the effective dressed state $|1⟩$
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

Reuse `active_manifold_angles` where appropriate.

*Implementation note:* the current homogeneous squeezing code should define one
effective qubit for the full saved wavefunction, or for the full
ensemble-averaged state, at each timestep. Do not select the central
$N_J=N/2$ sector and do not define one separate effective qubit per
$N_J$ sector. Instead, it should compute the collective moments

$
⟨J_x⟩, quad ⟨J_y⟩, quad ⟨J_z⟩, quad ⟨N_e⟩
$

using the full sector-block wavefunction and the Bloch-vector averaging rule,
then convert those moments into one $theta_J(t),phi_J(t)$.

== Define the instantaneous mean direction $|c⟩$
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

*Implementation note:* Rewrite the $S$-Bloch sphere operators $S_i$ into the $(u,d,e)$ basis.

So when we want

$
S_x =
frac(1,2) (ket(1) bra(u) + ket(u) bra(1)),
$

the code expands $|1⟩$ in the original basis:

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
with $theta_J,phi_J$ fixed from the effective qubit construction.

Using the above basis for the spin operators allows $theta_S,phi_S$ to be
computed in the regular way in
```
docs/instructions/bloch_vector_averaging.typ
```



== Define the $J$-fluctuation direction $|j⟩$
Step 3: Using the $J$ angles from step 1and $S$ angles from step 2, construct the state orthogonal to $|c⟩$ on the $J$-sphere ($(d,e)$ manifold):

$
|j⟩ = 0 |u⟩ - sin(theta_J / 2) |d⟩ + e^(-i phi_J) cos(theta_J / 2) |e⟩.
$

In vector form:

$
|j⟩ = mat(0; -sin(theta_J/2); e^(-i phi_J) cos(theta_J/2)).
$

== Define the $S$-fluctuation direction $|s⟩$

Step 4: Using the $J$-angles from step 1 and the $S$-angles from step 2, construct the state orthogonal to $|c⟩$ on the $S$-sphere ($(u,1)$ manifold):

$
|s⟩ = -sin(theta_S / 2) |u⟩ + e^(-i phi_S) cos(theta_S / 2) |1⟩.
$

In the $(u,d,e)$ basis,

$
|s⟩ = mat(-sin(theta_S/2); e^(-i phi_S) cos(theta_S/2) cos(theta_J/2); e^(-i(phi_S+phi_J)) cos(theta_S/2) sin(theta_J/2)).
$

The three states $|c⟩,|j⟩,|s⟩$ should form an orthonormal single-particle basis up to numerical precision.

== Construct the four local fluctuation operators

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

Do not construct full $3^N$-dimensional tensor-product operators unless absolutely necessary. Use the existing reduced $(N_J,n_e)$ basis and exploit symmetry. Each collective operator should be represented as a sparse operator acting on the reduced simulation basis.

If possible, construct $O_a$ from precomputed collective one-body transition operators

$
A_(mu nu)=sum_i ket(mu_i) bra(nu_i), quad mu,nu in (u,d,e).
$

Then

$
O_a = sum_(mu,nu) (o_a)_(mu nu) A_(mu nu).
$

This avoids building large tensor-product matrices.

== Construct the covariance matrix

For the current state $|psi(t)⟩$, compute

$
mu_a = ⟨ O_a⟩.
$

Then construct the $4 times 4$ covariance matrix

$
C_(a,b) = frac(1, 2) ⟨ O_a O_b+O_b O_a⟩ - mu_a mu_b.
$

For efficiency, avoid explicitly constructing $O_a O_b$. Instead use

$
C_(a,b) = "Re" [ ⟨ O_a psi | O_b psi⟩ ] - mu_a mu_b.
$

== Minimum fluctuation direction

Diagonalize the $4 times 4$ covariance matrix and take

$
lambda_("min")(C).
$

This is the minimum generalized transverse variance.

== Generalized squeezing parameter

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

For ensemble results, the required operator moments should be averaged over
trajectories before constructing the covariance matrix and squeezing
parameter. Do not compute one covariance matrix or one $xi^2$ per trajectory
and then average those final nonlinear quantities.

= Inhomogeneous Coupling Extension

This section should be used for two-group inhomogeneous simulations whose
sector keys are $(N_(J,1),N_(J,2))$. The homogeneous squeezing instruction above
remains the base algorithm. The inhomogeneous extension should only change how
the local dressed bases and fluctuation operators are grouped and combined.

Angle and Bloch-vector averages should follow:

```
  docs/instructions/bloch_vector_averaging.typ
```

== 1. Per-group squeezing

For each group $g in {1,2}$, compute a group-local squeezing parameter by
applying the homogeneous squeezing construction to that subgroup only.

When the saved wavefunction is organized by tuple sectors
$(N_(J,1),N_(J,2))$, the group-local angles should be extracted by supplying the
local operators for the relevant group, e.g. $J_(x,g),J_(y,g),J_(z,g),N_(e,g)$ to the bloch-vector averaging routine.
This should follow the same active-manifold angle logic as
in the single-group case, except that the operators act only on the chosen
group factor of the product basis.

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

== 2. Full-system inhomogeneous squeezing

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

The plotting entry points should use this data flow:

```python
ensemble = run_trajectory_ensemble(...)

plot_generalized_xi(result: TrajectoryEnsemble, ...):
    xi_data = generalized_squeezing_for_trajectory_or_ensemble(result)
    make_2x2_grid_plot(
        squeezing_db=10 * log10(xi_data["xi2_gen"]),
        covariance_eigenvalues=xi_data["covariance_eigvals"],
        dressed_population=xi_data["N_c"],
        excited_fraction=xi_data["excited_fraction_active"],
    )

plot_inhomogeneous_generalized_xi(result: TrajectoryEnsemble, ...):

    xi_data = generalized_squeezing_for_inhomogeneous(result)
    make_2x2_grid_plot(
        squeezing_db=10 * log10(xi_data["xi2"]),
        smallest_covariance_eigenvalue=xi_data["lambda_min"],
        dressed_population=xi_data["N_c"],
        excited_fraction=xi_data["excited_fraction_active"],
    )
```
