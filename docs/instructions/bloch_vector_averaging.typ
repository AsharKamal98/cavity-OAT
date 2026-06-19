#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[#text(size: 1.6em, weight: "bold")[Bloch-Vector Averaging Across $N_J$ Sectors]]

= Purpose

Use this note whenever changing code that computes, averages, or plots active-manifold Bloch vectors, Bloch angles, dressed-state directions, or quantities derived from $theta_J, phi_J$.

These Bloch-vector averaging rules apply after trajectory averaging has already
been resolved. At each time $t$, the input should be one effective state or
one already-averaged moment record.

For a sector-block state, write

$
|psi(t) chevron.r = sum_alpha |psi_alpha (t) chevron.r,
$

with
$
alpha = N_J
quad "for homogeneous sectors",
\
alpha = (N_(J,1), N_(J,2))
quad "for two-group inhomogeneous sectors".
$


= Within-Sector Averaging

Within a sector, expand the sector block in its internal excitation basis:

$
|psi_alpha chevron.r =
sum_(nu in cal(B)_alpha) c_(alpha,nu) |nu; alpha chevron.r,
$

where $nu$ is the internal basis label. For homogeneous sectors $nu=n_e$,
while for two-group inhomogeneous sectors $nu=(n_(e,1),n_(e,2))$.

For a collective component $J_(alpha,i)$ in sector $alpha$, compute the sector
expectation value as

$
chevron.l J_(alpha,i) chevron.r =
frac(chevron.l psi_alpha | J_(alpha,i) | psi_alpha chevron.r, chevron.l psi_alpha | psi_alpha chevron.r) =

(
  sum_mu sum_nu c_(alpha,mu)^* c_(alpha,nu)
  chevron.l mu; alpha | J_(alpha,i) | nu; alpha chevron.r
)
/
(
  sum_nu |c_(alpha,nu)|^2
).
$

Here both sums run over the internal excitation basis labels for sector $alpha$.

The normalized active-manifold Bloch direction is then
$
s_(alpha,i) = frac(2 chevron.l J_(alpha,i) chevron.r, N_("active",alpha)).
$

For the active manifold ${|↓ chevron.r, |e chevron.r}$,

$
N_("active",alpha) = chevron.l N_(↓,alpha) + N_(e,alpha) chevron.r
= 2 (chevron.l N_(e,alpha) chevron.r - chevron.l J_(alpha,z) chevron.r).
$


= Averaging Across Sectors

Never directly compare or average raw collective Bloch vector components $J_(alpha,i)$ from different sectors $alpha$. Such components scale with the number of active atoms in the sector. Compare only direction-like quantities $s_(alpha,i)$, or use an active-population-weighted average.

Averaged spin components can be found in one of two ways:
1. Weighted average where only the sector probabilities are used as weights:
$
s_i^("sector") =
frac(sum_alpha chevron.l psi_alpha|psi_alpha chevron.r s_(alpha,i), sum_alpha chevron.l psi_alpha|psi_alpha chevron.r)
$
2. Weighted average where the sector probabilities and active population of each sector are used as weights:
$
s_i =
frac(sum_alpha chevron.l psi_alpha|psi_alpha chevron.r N_("active",alpha) s_(alpha,i), sum_alpha chevron.l psi_alpha|psi_alpha chevron.r N_("active",alpha)).
$
Use the second method as the default for the full active-manifold Bloch direction. It is equivalent to summing the collective moments first and normalizing,

$
s_i =
frac(2 sum_alpha chevron.l psi_alpha | J_(alpha,i) | psi_alpha chevron.r, sum_alpha chevron.l psi_alpha|psi_alpha chevron.r N_("active",alpha)).
$

Use the first method only for explicitly sector-weighted diagnostics.

Once the averaged spin direction $s_i$ is computed, the angles can be found by
$
theta = arccos(-s_z), quad
phi = "atan2"(s_y, s_x).
$

The negative sign in $theta$ is a convention to make the $|d chevron.r$ state point towards the north pole in the $J$-bloch sphere.

This convention should be implemented in:
```
  common.utils.active_manifold_angles(...)
```


= Implementation note for Group-resolved Cases

Averaging across group-resolved sectors is already covered above by taking $alpha = (N_(J,1), N_(J,2))$.

For each group $a=1,2$, compute group-resolved angles from local group moments using the same active-manifold normalization:

$
N_("active",a) = 2 (chevron.l N_(e,a) chevron.r - chevron.l J_(z,a) chevron.r),
$

$
s_(i,a) = frac(2 chevron.l J_(i,a) chevron.r, N_("active",a)), quad
i in {x,y,z}.
$

Then use the same angle convention as above: $theta_a = arccos(-s_(z,a))$ and $phi_a = "atan2"(s_(y,a), s_(x,a))$.

The full inhomogeneous angle is not an arithmetic average of group angles. Compute it by summing group moments first,

$
chevron.l J_x chevron.r = chevron.l J_(x,1) chevron.r + chevron.l J_(x,2) chevron.r, quad
chevron.l J_y chevron.r = chevron.l J_(y,1) chevron.r + chevron.l J_(y,2) chevron.r,
$

$
chevron.l J_z chevron.r = chevron.l J_(z,1) chevron.r + chevron.l J_(z,2) chevron.r, quad
chevron.l N_e chevron.r = chevron.l N_(e,1) chevron.r + chevron.l N_(e,2) chevron.r,
$

and then calling `active_manifold_angles(...)` on the summed values.

Current implementation reference:

```
quantum_trajectories.inhomogeneous_diagnostics.inhomogeneous_group_angles(...)
```
