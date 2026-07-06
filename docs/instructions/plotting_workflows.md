# Plotting Workflows

This file defines repository-level plotting conventions. Plotting functions
should stay thin: they should visualize already-computed observables, moments,
or diagnostics rather than recomputing expensive physics.

## General Plotting Rules

1. Time axes should use scientific notation when useful and should disable
   additive offset notation:

   ```python
   ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)
   ```

2. Plotting functions should support both single-group and group-resolved data
   when the input container provides the relevant group fields.

3. When both group-resolved and full-system quantities are available, plot the
   group-resolved curves first and the full-system curve last.

4. Group-resolved curves should use dashed lines. Each group should keep its
   own color consistently across panels.

5. Full-system curves should use a solid line and accept an optional
   `colour_index`. The default `colour_index=0` should be gray, and later
   indices should cycle through other dark colors for overlaid runs.

6. Plotting functions should accept an optional `axes` argument. If `axes is
   None`, create a new figure and axes. If `axes` is provided, plot into those
   axes and return the same figure and axes. This allows overlaying multiple
   results:

   ```python
   fig, axes = plot_some_quantity(result_a, label="run A")
   fig, axes = plot_some_quantity(result_b, axes=axes, label="run B")
   ```

7. Plotting functions that save figures should accept `output_path`, create the
   parent directory if needed, and save with `dpi=200` and `bbox_inches="tight"`.

8. Labels should name the physical quantity being plotted, using math labels
   when appropriate.

9. Plotting functions should accept optional `phases` when the x-axis is time.

## `fig, axes = plot_spin_components(...)`

1. `plot_spin_components` lives in `common/plotting/j_spin.py`.

2. The function should take a series object directly, for example `moments.J`
   or another series with matching `x/y/z/length` field names.

3. The default call plots stored `x`, `y`, `z`, and `length` fields:

   ```python
   fig, axes = plot_spin_components(moments.J)
   ```

4. The function should read `x`, `y`, `z`, `length`, and optional
   `x_groups`, `y_groups`, `z_groups`, `length_groups`.

5. The output should be a `2x2` panel showing the x, y, z, and length
   components.

6. Each panel should include the full-system curve. If group fields exist, each
   panel should also include all group-resolved curves.

7. The function should support `axes`, `output_path`, `label`, `phases`,
   `colour_index`, and `linestyle` so multiple results can be plotted into the
   same figure.

## `fig, axes = plot_bloch_angles(...)`

1. `plot_bloch_angles` lives in `common/plotting/j_spin.py`.

2. The function should take a series object as input. The input must contain
   `t`, and may contain full-system `theta`/`phi`, group-resolved
   `theta_groups`/`phi_groups`, or both.

3. The output should be a `2x1` panel showing `theta` and `phi`.

4. `plot_bloch_angles` should read whatever stored angle fields are available on
   the input series. Angle construction should happen upstream; the plotting
   function should not recompute angles from spin components.

5. If group-resolved `theta_groups` and `phi_groups` exist, plot group angles
   first using the selected `colour_index` palette and selected `linestyle`.
   Plot the full-system angles last using the next color from the same palette
   when those fields exist.

6. The function should support `axes`, `output_path`, `label`, `phases`,
   `colour_index`, and `linestyle`.

## `fig, axes = plot_mfe_residuals(...)`

1. `plot_mfe_residuals` lives in
   `common/plotting/mfe_residuals.py`.

2. The function should take `moments.mfe_residuals` as input.

3. The function should only plot already-computed residuals from
   `moments.mfe_residuals.residuals_groups`. It should not recompute MFE residuals or
   call old observable extraction.

4. The output should be a single panel showing `Re R_1`, `Im R_1`,
   `Re R_2`, `Im R_2`, and the L2 norm `sqrt(|R_1|^2 + |R_2|^2)`.

5. The signed residual components should use the selected `colour_index`
   palette and `linestyle`. The L2 norm should use a solid gray curve.

6. The function should support `axes`, `output_path`, `label`, `phases`,
   `colour_index`, and `linestyle`.
   If `phases` are provided, optionally print the same phase-end residual
   summary as the old pipeline.

## `fig, axes = plot_sector_probabilities(...)`

1. This plot should show the normalized represented-sector probabilities
   `p_alpha(t_k)`, not raw sector amplitudes. For each represented sector
   `alpha`, define

   ```python
   w_alpha(t_k) = ||psi_alpha(t_k)||^2
   p_alpha(t_k) = w_alpha(t_k) / sum_beta w_beta(t_k)
   ```

2. `plot_sector_probabilities` may compute the normalized sector probabilities
   directly from `TrajectoryResult` or `TrajectoryEnsemble` snapshot data,
   using `TrajectoryEnsemble.trajectories[i].snapshots[k].sector_blocks`,
   which has the form
   ```python
   sector_blocks = {alpha: psi_alpha(t_k)}
   ```

3. The default output should be a single panel with one curve per represented
   sector, plotting `p_alpha(t_k)`. For homogeneous runs, the labels should use `N_J`. For
   two-group inhomogeneous runs, the labels should use the tuple key
   `(N_(J,1), N_(J,2))`.

4. This plot should be used to monitor how the relative weights/probabilities of the
   represented sectors change in time. It is a sector-ratio diagnostic, not a
   direct test for leakage into sectors that were not included in the
   simulation basis.

5. `plot_sector_probabilities` should live in
   `legacy/plotting_diagnostics.py`.

6. The function should support `linestyle`, for example `"-"` for solid,
   `"--"` for dashed, or `"-."` for dash-dot overlay comparisons.


Legacy note: the previous J-moment field names were `Jx`, `Jy`, `Jz`,
`Jx_groups`, `Jy_groups`, `Jz_groups`, `J_len`, and `sx`, `sy`, `sz`.

## Global Styling Rules

Shared styling helpers should live in `common/plotting/utils.py`.
Indexed color and line-pattern helpers may live in
`common/plotting/utils.py`.

1. Use a colorblind-friendly manual palette, such as Okabe-Ito, rather than
   changing Matplotlib's global style.

2. If `phases` are provided, show protocol phases with subtle background bands:

   ```python
   ax.axvspan(phase_start, phase_end, alpha=0.35)
   ```

3. Use a white figure background, white axes, light gray grids, hidden
   top/right spines, and small x-margins.

4. New figures should use constrained layout when possible, and figure titles
   should use slightly larger fonts with `fig.suptitle(..., y=1.02)`.
