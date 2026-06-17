#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Dephasing Diagnostics from Bloch-Vector Length]]

= Purpose

Use this note whenever implementing or modifying diagnostics intended to visualize dephasing in custom MCWF results.

For any averaging over different $N_J$ sectors, group-resolved sectors, or trajectory ensembles, follow `docs/instructions/bloch_vector_averaging.tex`. This file only adds the dephasing-specific rule that the averaged vector length is computed after the relevant vector components have been averaged.

= Goal

Visualize dephasing by plotting the length of the ensemble-averaged effective $S$-Bloch vector:

$
|⟨ bold(S)(t)⟩_("ens") | = sqrt( ⟨ S_x(t)⟩_("ens")^2+ ⟨ S_y(t)⟩_("ens")^2+ ⟨ S_z(t)⟩_("ens")^2 ).
$

Dephasing is loss of ensemble phase coherence. In a collective-spin picture, it appears as shrinkage of this averaged vector.

= Core Averaging Rule

Always average the vector components first, then compute the length:

$
|⟨ bold(S)(t)⟩_("ens") | = sqrt( ⟨ S_x(t)⟩_("ens")^2+ ⟨ S_y(t)⟩_("ens")^2+ ⟨ S_z(t)⟩_("ens")^2 ).
$

Do not average per-trajectory lengths:

$
frac(1, N_("traj")) sum_r |bold(S)_r(t) |.
$

The latter measures the typical conditioned-trajectory spin length, not ensemble dephasing.

= Recommended API

Implement this as standalone post-processing, for example:

```
plot_dephasing_bloch_lengths(
    result,
    *,
    normalize=True,
    axes=None,
    output_path=None,
)
```

The function should accept both `TrajectoryResult` and `TrajectoryEnsemble`.

For a single trajectory, compute the vector length from that trajectory's saved observables. For an ensemble, average $S_x,S_y,S_z$ over trajectories on the shared `t_eval` grid before computing the length.

Do not modify MCWF propagation unless the required observables are genuinely unavailable.

= Homogeneous Results

For homogeneous simulations, plot one curve:

$
|⟨ bold(S)_("total")(t)⟩|.
$

Use the existing effective-$S$ observable or moment-extraction path if one exists. Do not invent a second convention for $S_x,S_y,S_z$.

If effective-$S$ observables are not available, add only the minimal support needed to compute and post-process them.

= Inhomogeneous Results

For two-group inhomogeneous simulations, plot one figure with three curves: total, group 1, and group 2.

For groups $a=1,2$,

$
|⟨ bold(S)_a(t)⟩| = sqrt( ⟨ S_(x,a)(t)⟩^2+ ⟨ S_(y,a)(t)⟩^2+ ⟨ S_(z,a)(t)⟩^2 ).
$

For the total curve, sum components first:

$
S_(x,"total")=S_(x,1)+S_(x,2), quad S_(y,"total")=S_(y,1)+S_(y,2), quad S_(z,"total")=S_(z,1)+S_(z,2).
$

Then compute

$
|⟨ bold(S)_("total")(t)⟩| = sqrt( ⟨ S_(x,"total")(t)⟩^2+ ⟨ S_(y,"total")(t)⟩^2+ ⟨ S_(z,"total")(t)⟩^2 ).
$

If group-resolved effective-$S$ observables are unavailable, add the minimal support needed to compute them. Keep the diagnostic standalone.

= Normalization

If `normalize=True`, divide by

$
S_(max)=frac(N, 2)
$

for total homogeneous or effective-spin length, and

$
S_((max),a)=frac(N_a, 2)
$

for group-resolved inhomogeneous curves.

= Optional Relative-Phase Diagnostic

If group angles are already available, it is useful to also plot

$
Delta phi(t)=phi_1(t)-phi_2(t).
$

This can help distinguish relative dephasing between groups from dephasing within each group.

= Interpretation

- If the total length shrinks, ensemble coherence is being lost.
- If the total length shrinks while group-resolved lengths remain large, the dominant effect is likely relative dephasing between groups.
- If group-resolved lengths also shrink, there is dephasing within the groups as well.
