#set page(margin: 1in)
#set text(size: 11pt)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 17pt, weight: "bold")[Spin Direction Plot: Implementation Instructions]
]

= Scope

This file describes the intended post-processing plot for the normalized
active-manifold spin components (spin directions)

$
s_x(t), quad s_y(t), quad s_z(t).
$

The plot should follow the averaging convention in
`
    docs/instructions/bloch_vector_averaging.typ
`. 

= Recommended API

The plotting function should be used for both homogeneous and
inhomogeneous simulations. The homogeneous case should plot only the full-system
normalized spin components. The inhomogeneous case should plot the full-system components and
the group-resolved components on the same figure.

```python
plot_spin_direction(
    result: TrajectoryResult | TrajectoryEnsemble,
    *,
    group_resolved="auto",
    tol=1e-12,
    n_processes=None,
    averaged_observables=None,
    ax=None,
    output_path=None,
) -> tuple[dict, Figure, Axes]
```

Use `trajectory_observables(...)` for a `TrajectoryResult` and
`ensemble_observables(...)` for a `TrajectoryEnsemble`, unless a matching
`averaged_observables` object is supplied by the caller.

The function should have a clear two-stage structure:

1. Build a data dictionary from averaged observables.
2. Plot the data dictionary without recomputing observables.

The shared data-building flow should be:

```python
if averaged_observables is not None:
    obs = averaged_observables
elif isinstance(result, TrajectoryEnsemble):
    obs = ensemble_observables(result, tol=tol, n_processes=n_processes)
else:
    obs = trajectory_observables(result, tol=tol)

full = spin_component_bundle(
    label="full",
    sx=obs.sx,
    sy=obs.sy,
    sz=obs.sz,
)

data = {
    "t": obs.t,
    "full": full,
    "groups": [],
}

if group_resolved is True or (
    group_resolved == "auto" and obs.Jx_groups is not None
):
    for group_index in range(len(obs.Jx_groups)):
        _, _, _, sx_g, sy_g, sz_g = active_manifold_angles(
            obs.Jx_groups[group_index],
            obs.Jy_groups[group_index],
            obs.Jz_groups[group_index],
            obs.N_e_groups[group_index],
            tol=tol,
        )
        data["groups"].append(
            spin_component_bundle(
                label=f"group {group_index + 1}",
                sx=sx_g,
                sy=sy_g,
                sz=sz_g,
            )
        )
```

where the small shared helper is:

```python
def spin_component_bundle(label, sx, sy, sz):
    spin_length = sqrt(sx**2 + sy**2 + sz**2)
    return {
        "label": label,
        "sx": sx,
        "sy": sy,
        "sz": sz,
        "spin_length": spin_length,
    }
```

This single-function design is appropriate because the existing
`ObservableSeries.sx`, `ObservableSeries.sy`, and `ObservableSeries.sz` fields
are already the global active-manifold components for both scalar $N_J$ sectors
and tuple sectors $(N_(J,1), N_(J,2))$. The same normalization helper can be
reused for group-resolved raw moments when tuple-sector group observables are
present.

= Plotting

The plot should be a thin visualization of the observable data:

```python
make 2x2 axes with panels for sx, sy, sz, and spin_length
plot data["full"] on each panel
for each entry in data["groups"], plot that group on each panel
add phase-boundary vertical lines from the reference result phases when available
label each panel clearly
return data, fig, axes
```

Use clear panel labels, for example `$s_x$`, `$s_y$`, `$s_z$`, and `$|s|$`.
Use curve labels such as `full`, `group 1`, and `group 2`. The components are
already normalized by their corresponding $N_("active")$, so the plotting
function should not divide them by $N$, $N_J$, group size, or the number of
trajectories.

= Homogeneous and Inhomogeneous Results

Homogeneous and inhomogeneous results should share the same data-building and
plotting code wherever possible. The only branch should be whether
group-resolved raw moments are available.

For homogeneous results, `data["groups"]` should be empty and only the full
curve should be plotted.

For inhomogeneous results, the figure should include:

- full-system curves computed from `obs.sx`, `obs.sy`, and `obs.sz`;
- group 1 curves computed from
  `obs.Jx_groups[0]`, `obs.Jy_groups[0]`, `obs.Jz_groups[0]`,
  and `obs.N_e_groups[0]`;
- group 2 curves computed from
  `obs.Jx_groups[1]`, `obs.Jy_groups[1]`, `obs.Jz_groups[1]`,
  and `obs.N_e_groups[1]`.

The full-system inhomogeneous components are formed by summing group-resolved
moments first and then applying the same normalization rule, matching
`docs/instructions/bloch_vector_averaging.typ`. The group-resolved components
use the same rule separately for each group:

$
N_("active",g)(t) =
2 (chevron.l N_(e,g)(t) chevron.r - chevron.l J_(z,g)(t) chevron.r),
$

$
s_(x,g)(t) = frac(2 chevron.l J_(x,g)(t) chevron.r, N_("active",g)(t)), quad
s_(y,g)(t) = frac(2 chevron.l J_(y,g)(t) chevron.r, N_("active",g)(t)), quad
s_(z,g)(t) = frac(2 chevron.l J_(z,g)(t) chevron.r, N_("active",g)(t)).
$

If `group_resolved=True` is requested for homogeneous data, fail clearly because
there are no group-resolved moments to plot. If `group_resolved="auto"`, omit
group curves when group moments are absent.

= Data Requirements and Edge Cases

The plot should require no new simulation data. It should be reconstructible
from a saved `TrajectoryResult`, a `TrajectoryEnsemble`, or a matching
`ObservableSeries` produced from those objects.

All trajectories in an ensemble should share the same saved `t_eval` grid, as
already required by `ensemble_observables(...)`.

When $N_("active")(t)$ is zero or numerically below tolerance, use the behavior
of `active_manifold_angles(...)`; the returned components should remain finite
and the spin length should be computed from those returned components.
