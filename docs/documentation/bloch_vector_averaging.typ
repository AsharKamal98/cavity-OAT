#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 0.7em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Bloch-Vector Averaging Across $N_J$ Sectors]]

= Purpose

Use this note whenever changing code that computes, averages, or plots active-manifold Bloch vectors, Bloch angles, dressed-state directions, or quantities derived from $theta_J, phi_J$.

= Core Rule

Do not compare or average raw collective Bloch vectors from different $N_J$ sectors as directions.

Raw collective components such as

$
J_x, quad J_y, quad J_z
$

scale with the number of active atoms in the sector. A sector with larger $N_J$ can have a larger raw vector simply because it contains more atoms, not because its single-particle Bloch direction is different. Therefore, direction-like quantities must be computed from normalized single-particle components.

For the active manifold ${|↓ chevron.r, |e chevron.r}$, the code uses

$
N_("active") = chevron.l N_↓ + N_e chevron.r = 2 (chevron.l N_e chevron.r - chevron.l J_z chevron.r),
$

because

$
J_z = frac(N_e - N_↓, 2).
$

The normalized active-manifold Bloch direction is then

$
s_x = frac(2 chevron.l J_x chevron.r, N_("active")), quad
s_y = frac(2 chevron.l J_y chevron.r, N_("active")), quad
s_z = -frac(2 chevron.l J_z chevron.r, N_("active")).
$

Angles are computed only after this normalization:

$
theta = arccos(s_z), quad
phi = "atan2"(s_y, s_x).
$

This is implemented in:

```
common.utils.active_manifold_angles(...)
```

= How Sector Averaging Should Work

For a direct-sum wavefunction over sectors,

$
|psi chevron.r = ⊕_alpha |psi_alpha chevron.r,
$

where $alpha$ may be a homogeneous sector $N_J$ or an inhomogeneous sector $(N_(J,1), N_(J,2))$, first compute expectation values with the correct sector probabilities:

$
chevron.l J_a chevron.r = frac(sum_alpha chevron.l psi_alpha | J_(a,alpha) | psi_alpha chevron.r, sum_alpha chevron.l psi_alpha | psi_alpha chevron.r),
$

and similarly for $chevron.l N_e chevron.r$. Then compute $N_("active")$, $s_x, s_y, s_z$, $theta$, and $phi$ from those averaged expectation values.

Equivalently, this is an active-population-weighted average of the normalized sector directions. It is not an unweighted average over sectors and it is not an average of raw collective vector lengths.

Do not do:

$
theta_("avg") = frac(1, N_("sectors")) sum_alpha theta_alpha, quad
phi_("avg") = frac(1, N_("sectors")) sum_alpha phi_alpha.
$

Do not do:

$
bold(s)_("avg") = frac(1, N_("sectors")) sum_alpha
( frac(2 J_(x,alpha), N_("active",alpha)),
  frac(2 J_(y,alpha), N_("active",alpha)),
 -frac(2 J_(z,alpha), N_("active",alpha)) )
$

unless an explicitly unweighted diagnostic is desired and clearly labeled as such.

= Inhomogeneous Group Averages

For two-group inhomogeneous coupling, group-resolved angles $(theta_1, phi_1)$ and $(theta_2, phi_2)$ are computed by applying the same normalization separately to each group:

$
N_("active",a) = 2 (chevron.l N_(e,a) chevron.r - chevron.l J_(z,a) chevron.r),
$

$
s_(x,a) = frac(2 chevron.l J_(x,a) chevron.r, N_("active",a)), quad
s_(y,a) = frac(2 chevron.l J_(y,a) chevron.r, N_("active",a)), quad
s_(z,a) = -frac(2 chevron.l J_(z,a) chevron.r, N_("active",a)).
$

The combined or average angle should be computed by summing the group expectation values first,

$
chevron.l J_x chevron.r = chevron.l J_(x,1) chevron.r + chevron.l J_(x,2) chevron.r, quad
chevron.l J_y chevron.r = chevron.l J_(y,1) chevron.r + chevron.l J_(y,2) chevron.r,
$

$
chevron.l J_z chevron.r = chevron.l J_(z,1) chevron.r + chevron.l J_(z,2) chevron.r, quad
chevron.l N_e chevron.r = chevron.l N_(e,1) chevron.r + chevron.l N_(e,2) chevron.r,
$

and then calling `active_manifold_angles(...)` on the summed values. This avoids arithmetic angle averaging such as $(theta_1 + theta_2) / 2$, which is generally wrong.

Current implementation reference:

```
quantum_trajectories.inhomogeneous_diagnostics.inhomogeneous_group_angles(...)
```

= Trajectory-Ensemble Averages

For ensemble observables, average the required expectation values over trajectories on the common `t_eval` grid, then compute the normalized Bloch direction from the averaged moments.

Current implementation reference:

```
quantum_trajectories.aggregator.ensemble_observables(...)
```

This is the same physical logic used for squeezing: direction and covariance information should be constructed from ensemble-averaged moments when the target object is the unconditioned density matrix.

= When To Use This Logic

Use this logic for:

- active-manifold angles $theta_J, phi_J$;
- dressed-state directions such as $|1 chevron.r$, $|c chevron.r$, $|j chevron.r$, and $|s chevron.r$;
- plots comparing Bloch angles across sectors, groups, or trajectories;
- any diagnostic where the intended object is a single-particle direction rather than a collective vector length.

Do not use this normalization for genuinely extensive observables where the total collective quantity is the object of interest, such as total $chevron.l J_x chevron.r$, total $chevron.l N_e chevron.r$, total $chevron.l N_J chevron.r$, jump rates, or jump counts.
