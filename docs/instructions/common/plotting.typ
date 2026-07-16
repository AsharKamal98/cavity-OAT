#set page(paper: "a4", margin: 1in)
#set text(size: 11pt)
#set par(first-line-indent: 0pt, spacing: 1em)
#set heading(numbering: "1.")

#align(center)[
  #text(size: 1.6em, weight: "bold")[Common Plotting Instructions]
]

= Purpose

This file lists shared plotting conventions for repository-level plotting code,
especially helpers and functions under `common/plotting/`. Use it when adding
or revising plotting functions that should be reusable across simulation
backends.

For shared overlay, axes, palette, and line-style behavior, use the
plotting-workflows skill.

= Plot Functions

== `fig, axes = plot_spin_components(...)`

This function lives in `common/plotting/j_spin.py`.

- The function should take a series object directly, for example `moments.J`
  or another series with matching `x/y/z/length` field names.
- The default call should plot the stored `x`, `y`, `z`, and `length`
  fields.
- If `normalized=True`, the function should instead plot `nx`, `ny`, `nz`,
  and `length`.
- The function should read `x`, `y`, `z`, `length`, and optional
  `x_groups`, `y_groups`, `z_groups`, `length_groups`. For normalized plots it
  should analogously read `nx`, `ny`, `nz`, and optional group-resolved
  normalized fields.
- The output should be a `2x2` panel showing x, y, z, and length.
- The current common implementation plots only group-resolved curves when the
  corresponding `x_groups`, `y_groups`, `z_groups`, or `length_groups` fields
  exist. It currently does not draw the full-system curves even if `x`, `y`,
  `z`, and `length` are also present.
- The function should support `axes`, `output_path`, `label`, `phases`,
  `colour_family_index`, `shade_index`, and `linestyle`.
- For shared overlay and palette behavior, use the plotting-workflows skill.

== `fig, axes = plot_bloch_angles(...)`

This function lives in `common/plotting/j_spin.py`.

- The function should take a series object as input.
- The input must contain `t`, and may contain full-system `theta`/`phi`,
  group-resolved `theta_groups`/`phi_groups`, or both.
- The output should be a `2x1` panel showing `theta` and `phi`.
- The plotting function should read whatever stored angle fields are available
  on the input series. Angle construction should happen upstream; the plotting
  function should not recompute angles from spin components.
- If group-resolved fields exist, plot those first using the selected palette
  and `linestyle`.
- The current common implementation plots only group-resolved angles when
  `theta_groups` and `phi_groups` exist. It currently does not draw the
  full-system `theta` and `phi` curves even if those fields are also present.
- The function should support `axes`, `output_path`, `label`, `phases`,
  `colour_family_index`, `shade_index`, and `linestyle`.
- For shared overlay and palette behavior, use the plotting-workflows skill.

== `fig, axes = plot_mfe_residuals(...)`

This function lives in `common/plotting/mfe_residuals.py`.

- The function should take `moments.mfe_residuals` as input.
- The function should only plot already-computed residuals from
  `moments.mfe_residuals.residuals_groups`. It should not recompute MFE
  residuals inside the plotting layer.
- The current common implementation should assume exactly two residual groups.
- The output should be a single panel showing `Re R_1`, `Im R_1`, `Re R_2`,
  `Im R_2`, and the L2 norm `sqrt(|R_1|^2 + |R_2|^2)`.
- With `show_components=False`, plot only the nonnegative L2 norm.
- With `symlog=True` (the default), use a symmetric logarithmic y-axis with
  `linthresh=1e-5` for signed residual components. For an L2-only plot, use a
  regular logarithmic y-axis when all values are positive; fall back to symlog
  if any L2 value is zero.
- The function should support `axes`, `output_path`, `label`, `phases`,
  `colour_family_index`, `shade_index`, `linestyle`, `symlog`, and
  `show_components`.
- If `phases` are provided, the function may also print a phase-end residual
  summary.
- Resolve residual colors through the shared `colour_family_index` and
  `shade_index` palette controls. The signed residuals and L2 norm cycle
  through that palette in plotted-curve order.
- For shared overlay and line-style behavior, use the plotting-workflows
  skill.

= Shared Plotting Rules

== Plotting Scope

- Plotting functions should stay thin: they should visualize already-computed
  observables, moments, or diagnostics rather than recomputing expensive
  physics.
- Shared plotting helpers should live in `common/plotting/utils.py`.
- Common plotting functions should support both single-group and
  group-resolved data when the input container provides the relevant group
  fields.

== Time and Phase Handling

- Time axes should use scientific notation when useful and should disable
  additive offset notation.
- Plotting functions for time-series data should accept optional `phases`.
- If `phases` are provided, show protocol phases with subtle background bands
  and phase boundaries.
- Phase-boundary logic should be reused from shared plotting helpers rather
  than recomputed inside each plotting function.

== Styling

- Use a white figure background, white axes, light gray grids, hidden top/right
  spines, and small x-margins.
- New figures should use constrained layout when possible.
- Figure titles should use slightly larger fonts, for example with
  `fig.suptitle(..., y=1.02)`.
- Labels should name the physical quantity being plotted, using math labels
  when appropriate.

== Skill-Routed Behavior

For the following shared plotting behaviors, use the plotting-workflows skill
instead of restating local rules in each file:

- optional `fig, axes` overlay behavior;
- shared color-palette selection, including `colour_palette(...)`,
  `colour_family_index`, and `shade_index`;
- `linestyle` conventions.


= Invariants

- Shared plotting behavior should be implemented once in `common/plotting/`
  helpers where possible, not redefined separately in each plotting function.
- Backend-specific plotting assumptions should not be hidden inside common
  plotting helpers.
- If a plotting rule is reused across multiple plot types, it should live in
  this file or in the plotting-workflows skill, not only in one plot-specific
  instruction.
