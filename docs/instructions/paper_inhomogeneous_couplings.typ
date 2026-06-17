#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Inhomogeneous Couplings: Implementation Instructions]]

This file gives repository-specific implementation rules for inhomogeneous couplings in the custom MCWF code. For the physics derivation and mean-field equations, use `docs/theory/notes_inhomogeneous_coupling.tex`. If this instruction file and the theory note disagree on equations, check the theory note first and flag the inconsistency.

= Scope

Inhomogeneous coupling should be implemented as a two-group extension of the existing strong-symmetry sector code. The code should keep the same top-level MCWF workflow as homogeneous simulations:

```
sector_coeffs -> build_precomputed_trajectory_data(...)
              -> simulate_single_trajectory(...)
              -> TrajectoryResult / TrajectoryEnsemble
```

Only the sector keys, sector operators, Hamiltonian pieces, and jump operator should change.

= User-Facing Inputs

High-level inhomogeneous setup should use:

```
N, dN, N1, N2, omega1, sector_distribution
```

with

$
N=N_1+N_2.
$

The second coupling weight should be computed once from the physical group sizes:

$
omega_2=frac(N-N_1 omega_1, N_2), quad N_1 omega_1+N_2 omega_2=N.
$

This convention should be used so homogeneous and inhomogeneous runs have the same atom-number weighted average coupling. Do not recompute $omega_2$ from instantaneous, trajectory-averaged, or sector-averaged active populations.

If $N_2=0$, the code should avoid dividing by $N_2$; the second group is empty and $omega_2$ is irrelevant.

= Sector Keys and State Blocks

Homogeneous runs use scalar sector keys:

```
{Nj: coeff}
```

Inhomogeneous runs should use group-resolved sector keys:

```
{(Nj1, Nj2): coeff}
```

The physical active-manifold count should be

$
N_J=N_(J,1)+N_(J,2).
$

The coupling weights $omega_1,omega_2$ should affect operators, not the definition of the sector label. In particular, do not define sectors by $omega_1 N_(J,1)+omega_2 N_(J,2)$.

For a fixed group-resolved sector, the internal state block should be an array over the product active-manifold excitation basis:

```
psi[(ne1, ne2)]  with  ne1 = 0,...,Nj1 and ne2 = 0,...,Nj2
```

Flattened arrays should therefore have shape:

$
(N_(J,1)+1)(N_(J,2)+1).
$

The default internal state should be the all-active-atoms-down state:

$
(n_(e,1),n_(e,2))=(0,0).
$

= High-Level Sector Construction

User-facing code should not require manually writing every $(N_(J,1),N_(J,2))$ coefficient. A helper should construct the low-level dictionary from the high-level inputs.

The selected total active sectors should be:

$
N_J in {N/2-d N, dots.h, N/2+d N},
$

after removing invalid values. For each selected $N_J$, include every pair satisfying:

$
N_(J,1)+N_(J,2)=N_J, quad 0<= N_(J,1)<= N_1, quad 0<= N_(J,2)<= N_2.
$

The helper output should be:

```
{(Nj1, Nj2): coeff}
```

== Sector Distribution

The same high-level distribution names should be supported as in homogeneous runs.

`sector_distribution="binomial"` should use product-state binomial weights over group-resolved sectors, restricted to the selected window and renormalized.

`sector_distribution="square"` should preserve the old homogeneous meaning: equal total probability over selected total $N_J$ sectors. For each fixed $N_J$, that total probability should be split across valid $(N_(J,1),N_(J,2))$ pairs using the conditional binomial distribution.

Do not make all group-resolved pairs equally probable by default. If that behavior is ever needed, it should be added as an explicit new option, for example:

```
sector_distribution="square_group_resolved"
```

= Operator Construction

For a sector $(N_(J,1),N_(J,2))$, group operators should act only on their own group and as identity on the other group:

$
J_(1,-)=J_-(N_(J,1))⊗ I_2, quad J_(2,-)=I_1⊗ J_-(N_(J,2)),
$

$
N_(e,1)=N_e(N_(J,1))⊗ I_2, quad N_(e,2)=I_1⊗ N_e(N_(J,2)).
$

The drive Hamiltonian should use the weighted group operators:

$
H_Omega = Omega(omega_1J_(1,x)+omega_2J_(2,x)).
$

The detuning should remain homogeneous across groups:

$
H_delta=-delta(N_(e,1)+N_(e,2)).
$

All group-sector operators should be built during precomputation when possible. Do not rebuild sparse Kronecker products, $l^dagger l$, effective generators, or full-`dt` propagators inside the trajectory time loop.

For general precompute rules, use `docs/instructions/simulation_precompute.typ`.

= Collective Jump Convention

The inhomogeneous collective jump should be one shared weighted collective channel:

$
A=omega_1J_(1,-)+omega_2J_(2,-).
$

This is not two independent jump channels. Do not replace it by separate $l_1=omega_1J_(1,-)$ and $l_2=omega_2J_(2,-)$ unless independent group-resolved decay is explicitly requested.

In the shifted-jump picture, the jump operator should be:

$
l=A+i frac(Omega, Gamma).
$

This shifted form should absorb the weighted drive Hamiltonian consistently with the theory note.

= Simulation Behavior

The same ensemble and single-trajectory machinery should be used for homogeneous and inhomogeneous runs. The code should dispatch based on sector-key type:

```
int key        -> homogeneous sector path
(int, int) key -> inhomogeneous two-group sector path
```

The simulation result should continue to store snapshots as sector blocks:

```
snapshot.sector_blocks[(Nj1, Nj2)] = psi_Nj1_Nj2
```

The homogeneous limit should be recoverable numerically when the weights are equal or when one group is empty, although using two explicit groups can still be slower than the scalar homogeneous path.

= Diagnostics and Observables

Group-resolved diagnostics should compute observables from saved sector blocks where possible. Add only the minimum extra support needed to compute:

$
J_(x,1),J_(y,1),J_(z,1),N_(e,1), quad J_(x,2),J_(y,2),J_(z,2),N_(e,2).
$

Angle and Bloch-vector averages across sectors or trajectories should follow `docs/instructions/bloch_vector_averaging.typ`. In particular, do not average angles directly. Average the appropriate normalized components or moments first, then construct $theta$ and $phi$.

= Mean-Field Residual Diagnostic

The residual check should be standalone post-processing. It should not change core trajectory propagation unless missing data must be exposed.

The residual equation and definitions are given in `docs/theory/notes_inhomogeneous_coupling.tex`. The implementation should:

- use the same fixed $omega_1,omega_2$ as the simulation;
- use the phase-local $Omega(t)$ and $delta(t)$;
- use the simulation $Gamma$;
- compute group-resolved $theta_1,phi_1,theta_2,phi_2$;
- evaluate $R_1(t)$ and $R_2(t)$;
- print or report $|R_1|+|R_2|$ at phase ends when requested;
- plot $|"Re" R_1|$, $|"Im" R_1|$,
    $|"Re" R_2|$, and $|"Im" R_2|$ in a $2 times 2$
    grid.

For multiple group-resolved sectors, residuals should be evaluated per $(N_(J,1),N_(J,2))$ sector before any probability-weighted aggregation. Do not first average nonlinear inputs globally and then plug those averages into the residual equation.

= Plotting Helpers

Inhomogeneous plotting helpers should be standalone diagnostics. They should not duplicate core observable code unless needed.

The group-angle plot should show:

- $theta_1(t)$, $theta_2(t)$, and the correctly averaged
    $theta_("avg")(t)$;
- $phi_1(t)$, $phi_2(t)$, and the correctly averaged
    $phi_("avg")(t)$.

The average angle convention should follow `docs/instructions/bloch_vector_averaging.typ`.
