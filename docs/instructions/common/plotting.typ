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

- The function should take `t` and x, y, z, and length inputs. Each input may
  be one array or a sequence of arrays.
- Each position across the four sequences describes one Bloch vector, and the
  matching `labels` entry identifies that vector in every panel.
- The caller chooses whether the supplied curves are full-system,
  group-resolved, normalized, or derived combinations. The plotting function
  should not inspect a moment class or transform the supplied data.
- The output should be a `2x2` panel showing x, y, z, and length.
- The function should support `axes`, `output_path`, `phase_protocol`, `title`,
  `colour_family_index`, `shade_index`, and `linestyle`.
- For shared overlay and palette behavior, use the plotting-workflows skill.

== `fig, axes = plot_bloch_angles(...)`

This function lives in `common/plotting/j_spin.py`.

- The function should take `t` and matching theta and phi inputs. Each input
  may be one array or a sequence of arrays.
- Each position across the two sequences describes one Bloch vector, and the
  matching `labels` entry identifies that vector in both panels.
- The output should be a `2x1` panel showing `theta` and `phi`.
- The caller chooses whether the supplied curves are full-system,
  group-resolved, or derived combinations. The plotting function should not
  inspect a moment class or recompute angles.
- The function should support `axes`, `output_path`, `phase_protocol`,
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
- The function should support `axes`, `output_path`, `label`, `phase_protocol`,
  `colour_family_index`, `shade_index`, `linestyle`, `symlog`, and
  `show_components`.
- If `phase_protocol` is provided, the function may also print a family-phase-end residual
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
- Common plotting functions should accept caller-selected curve collections so
  single-group, group-resolved, full-system, and derived data can use the same
  plotting code.

== Time and Phase Handling

- Time axes should use scientific notation when useful and should disable
  additive offset notation.
- Plotting functions for time-series data should accept an optional
  `phase_protocol`.
- If `phase_protocol` is provided, show its `family_phases` with subtle
  background bands and family-phase boundaries. Do not shade individual
  integration `Phase` segments.
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
- Multi-panel plots should show shared curve labels once in a figure legend
  centered below the panel grid.

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
