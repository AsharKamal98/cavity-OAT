#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Instruction: Implement Generalized Three-Level Squeezing Parameter]]

Construct the generalized three-level squeezing parameter for the custom MCWF code. The goal is to compute $xi^2_("gen")(t)$ at each saved timestep after the simulation has run.

The squeezing calculation should preferably be implemented as a standalone post-processing function that takes the simulation result or observable output as input. If it is simpler, it may be integrated into the observable function, but it must be optional via a boolean flag, since this calculation may be expensive and is not always needed.

If additional per-timestep data must be saved for this calculation, add a boolean flag such as `save_squeezing_data` or `compute_generalized_squeezing`.

= Per-timestep calculation

At each timestep, do the following.

== 1. Define the effective dressed state $|1⟩$

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
docs/instructions/bloch_vector_averaging.tex
```

Reuse `active_manifold_angles` where appropriate.

== 2. Define the instantaneous mean direction $|c⟩$

Use the ansatz

$
|c⟩ = cos(theta_S / 2) |u⟩ + e^(-i phi_S) sin(theta_S / 2) |1⟩.
$

Here $theta_S,phi_S$ should be found by comparing this ansatz to the instantaneous effective $S$-Bloch vector of the simulation state. The angles $theta_J,phi_J$ appearing inside $|1⟩$ are the angles found in step 1.

In the explicit $(u,d,e)$ basis,

$
|c⟩ = mat(cos(theta_S/2); e^(-i phi_S) sin(theta_S/2) cos(theta_J/2); e^(-i(phi_S+phi_J)) sin(theta_S/2) sin(theta_J/2)).
$

== 3. Define the $J$-fluctuation direction $|j⟩$

Using the $J$-angles found in step 1, construct the state orthogonal to $|1⟩$ inside the $(d,e)$ manifold:

$
|j⟩ = 0 |u⟩ - sin(theta_J / 2) |d⟩ + e^(-i phi_J) cos(theta_J / 2) |e⟩.
$

In vector form:

$
|j⟩ = mat(0; -sin(theta_J/2); e^(-i phi_J) cos(theta_J/2)).
$

== 4. Define the $S$-fluctuation direction $|s⟩$

Using the $S$-angles from step 2 and the $J$-angles from step 1, construct

$
|s⟩ = -sin(theta_S / 2) |u⟩ + e^(-i phi_S) cos(theta_S / 2) |1⟩.
$

In the $(u,d,e)$ basis,

$
|s⟩ = mat(-sin(theta_S/2); e^(-i phi_S) cos(theta_S/2) cos(theta_J/2); e^(-i(phi_S+phi_J)) cos(theta_S/2) sin(theta_J/2)).
$

The three states $|c⟩,|j⟩,|s⟩$ should form an orthonormal single-particle basis up to numerical precision.

== 5. Construct the four local fluctuation operators

Construct the four $3 times 3$ single-particle operators:

$
o_1 = frac("ket"(c) "bra"(j) + "ket"(j) "bra"(c), 2),
$

$
o_2 = frac("ket"(c) "bra"(j) - "ket"(j) "bra"(c), 2 i),
$

$
o_3 = frac("ket"(c) "bra"(s) + "ket"(s) "bra"(c), 2),
$

$
o_4 = frac("ket"(c) "bra"(s) - "ket"(s) "bra"(c), 2 i).
$

Then construct the corresponding collective operators

$
O_a = sum_i o_a^((i)), quad a=1,2,3,4.
$

Do not construct full $3^N$-dimensional tensor-product operators unless absolutely necessary. Use the existing reduced $(N_J,n_e)$ basis and exploit symmetry. Each collective operator should be represented as a sparse operator acting on the reduced simulation basis.

If possible, construct $O_a$ from precomputed collective one-body transition operators

$
A_(mu nu)=sum_i "ket"(mu_i) "bra"(nu_i), quad mu,nu in (u,d,e).
$

Then

$
O_a = sum_(mu,nu) (o_a)_(mu nu) A_(mu nu).
$

This avoids building large tensor-product matrices.

== 6. Construct the covariance matrix

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

== 7. Minimum fluctuation direction

Diagonalize the $4 times 4$ covariance matrix and take

$
lambda_("min")(C).
$

This is the minimum generalized transverse variance.

== 8. Generalized squeezing parameter

Compute

$
xi^2_("gen")(t) = frac(N lambda_("min")(C), ⟨ N_c/2⟩^2),
$

where

$
N_c = sum_i "ket"(c_i) "bra"(c_i).
$

Compute $⟨ N_c⟩$ using the same collective-operator construction:

$
N_c = sum_(mu,nu) ("ket"(c) "bra"(c))_(mu nu) A_(mu nu).
$

If the state is well polarized along $|c⟩$, then $⟨ N_c/2⟩ approx N/2$, but the code should compute it explicitly rather than assuming this.

= Data that may need to be saved per timestep

Save enough data to reconstruct the squeezing calculation after the simulation:

- the full state vector $|psi(t)⟩$ in the reduced $(N_J,n_e)$
    basis, or equivalent data sufficient to reconstruct it;
- the timestep values $t$;
- $N$, $Gamma$, $Omega$, $delta$, and any sector truncation
    metadata;
- the basis/index mapping for $(N_J,n_e)$;
- if not recomputed during post-processing, save
    $theta_J(t),phi_J(t),theta_S(t),phi_S(t)$.

Prefer recomputing angles from the state during post-processing if this avoids storing redundant data.
