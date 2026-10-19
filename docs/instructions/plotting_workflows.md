# Plotting Workflows

This file defines repo-specific contracts for the current plotting functions.
For shared plotting conventions under `common/plotting/`, use
`docs/instructions/common/plotting.typ`. For shared overlay, `fig, axes`,
palette, and line-style behavior, use the plotting-workflows skill.

## `fig, axes = plot_spin_components(...)`

1. `plot_spin_components` lives in `common/plotting/j_spin.py`.

2. The function should take explicit `t`, `x_components`, `y_components`,
   `z_components`, and `lengths` inputs. Each component input may be one array
   or a sequence of arrays.

3. The caller chooses the curves and passes one matching label per curve. For
   example, group-resolved moments can be plotted with:

   ```python
   fig, axes = plot_spin_components(
       moments.J.t,
       moments.J.x,
       moments.J.y,
       moments.J.z,
       moments.J.length,
       x_sem=moments.J.x_sem,
       y_sem=moments.J.y_sem,
       z_sem=moments.J.z_sem,
       labels=group_labels,
   )
   ```

4. The plotting function should not inspect `JMomentSeries` or automatically
   add group-resolved or aggregated fields. The caller may instead pass
   `moments.J.aggregated.*` explicitly when plotting the aggregated curve.
   Match aggregate curves with `moments.J.aggregated.x_sem/y_sem/z_sem`.

5. The output should be a `2x2` panel showing the x, y, z, and length
   components.

6. The function should support `axes`, `output_path`, `labels`,
   `phase_protocol`, `title`,
   `colour_family_index`, `shade_index`, `linestyle`, and `marker` so multiple
   results can be plotted into the same figure.

7. Optional `x_sem`, `y_sem`, `z_sem`, and `length_sem` inputs draw shaded
   `value ± SEM` bands beneath the central lines. Darken the matching line's
   RGB channels by a factor of `0.75` and use alpha `0.45`. They accept one
   array or one entry per curve; `None` omits all errors for that panel, and
   a sequence may contain `None` for individual curves. Preserve central
   line styles and a legend containing only the central curves. Validate
   real nonnegative errors with the same shapes as their values and `t`.
   Allow `NaN` to omit shading at undefined derived points; reject infinities.

8. MCWF notebooks plotting raw J components should pass their matching SEMs,
   preserving group selection, curve order, and any time slicing. Normalize
   errors only when the same already-defined linear scaling is applied to the
   central values. J lengths and normalized directions have their own matching
   `length_sem/nx_sem/ny_sem/nz_sem` fields from nonlinear propagation; use
   those rather than Cartesian SEMs. J-mode vectors expose their own matching
   `x_sem/y_sem/z_sem/length_sem` fields; pass these for MCWF bright/dark plots,
   preserving the same time selection in phase-specific plots. MFE mode SEMs
   remain `None`. S components still require their own uncertainty calculation.
   Legacy loaded results without SEM fields may use `getattr(..., None)`.

## `fig, axes = plot_bloch_angles(...)`

1. `plot_bloch_angles` lives in `common/plotting/j_spin.py`.

2. The function should take explicit `t`, `theta_curves`, and `phi_curves`
   inputs. Each angle input may be one array or a sequence of arrays, with one
   matching label per curve.

3. The output should be a `2x1` panel showing `theta` and `phi`.

4. The caller chooses whether to pass group-resolved or aggregated angles.
   The plotting function should not inspect `JMomentSeries` or recompute angles
   from spin components.

5. The function should support `axes`, `output_path`, `labels`,
   `phase_protocol`,
   `colour_family_index`, `shade_index`, `linestyle`, and `marker`.

## `fig, axes = plot_transverse_spin_components(...)`

1. `plot_transverse_spin_components` lives in `common/plotting/specific.py`.

2. It accepts x, y, z, and full-length curves on a shared time grid. The
   selected transverse quantity is either `sqrt(x^2+y^2)` or `x+y`.

3. The output is a `3x1` panel showing `xy_length`, z, and full length.

4. It follows the shared axes, phase-region, palette, line-style, marker,
   legend, and output-path conventions.

## `fig, axes = plot_mfe_residuals(...)`

1. `plot_mfe_residuals` lives in
   `common/plotting/mfe_residuals.py`.

2. The function should take `moments.mfe_residuals` as input.

3. The function should only plot already-computed residuals from
   `moments.mfe_residuals.residuals_groups`. It should not recompute MFE
   residuals.

4. The output should be a single panel showing `Re R_1`, `Im R_1`,
   `Re R_2`, `Im R_2`, and the L2 norm `sqrt(|R_1|^2 + |R_2|^2)`.

5. The signed residual components and L2 norm should use the selected
   `colour_family_index` / `shade_index` palette, cycling by plotted-curve
   order. The L2 norm should use a solid line.

6. The function should support `axes`, `output_path`, `label`, `phases`,
   `colour_family_index`, `shade_index`, `linestyle`, and `marker`.
   If `phases` are provided, optionally print the phase-end residual summary.

## `fig, axes = plot_harmonic_sweep(...)`

1. `plot_harmonic_sweep` lives in
   `common/plotting/harmonic_analysis.py`.

2. It accepts already-extracted parameter values and one or more matching
   curves of fundamental frequencies, total harmonic distortions, RMS
   oscillation amplitudes, and fitted offsets. It must not inspect simulation
   outputs or recompute harmonic analysis.

3. The output is a `2x2` panel with fundamental frequency and RMS oscillation
   amplitude in the top row, followed by fitted offset and THD in the bottom
   row. All panels share the parameter x-axis.
   The figure title must have a clear vertical gap above the panel titles.
   If `family_phase_index` is supplied, shade all panel backgrounds using that
   shared phase color without changing the outer figure background.

4. The function supports one `labels` entry per curve pair, plus
   `parameter_label`, `family_phase_index`, `title`, `axes`,
   `output_path`, `colour_family_index`, `shade_index`, `linestyle`, and
   `marker`.
   Repeated calls with supplied `axes` overlay additional curves and preserve
   the existing layout. Its figure legend uses at most three columns and
   preserves row-major call order, so each three-curve call occupies one row.
   Each call restarts color assignment from the beginning of its selected
   palette.

5. Each curve must be non-empty, one-dimensional, and match the varied
   parameter array. Parameters, frequency curves, RMS amplitudes, and offsets
   must be finite. THD curves may contain `NaN` for non-oscillating series but
   must not contain infinities. All four metric groups and labels must have
   equal counts. Their supplied order is preserved.

6. For windowed analysis, pass window-centre times as the parameter values and
   use a time-axis label. Select finite frequency, RMS amplitude, and offset
   entries from each curve before plotting; select the same entries for THD.
   Constant signals may retain `NaN` THD. In the bright/dark MCWF notebook,
   use one color per dark component and solid MCWF/dashed MFE overlays.
