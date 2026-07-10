# Plotting Workflows

This file defines repo-specific contracts for the current plotting functions.
For shared plotting conventions under `common/plotting/`, use
`docs/instructions/common/plotting.typ`. For shared overlay, `fig, axes`,
palette, and line-style behavior, use the plotting-workflows skill.

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
   `colour_family_index`, `shade_index`, and `linestyle` so multiple results
   can be plotted into the same figure.

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
   first using the selected `colour_family_index` / `shade_index` palette and
   selected `linestyle`. Plot the full-system angles last using the next color
   from the same palette when those fields exist.

6. The function should support `axes`, `output_path`, `label`, `phases`,
   `colour_family_index`, `shade_index`, and `linestyle`.

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
