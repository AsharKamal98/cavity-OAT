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

## `fig, axes = plot_j_spin_components(...)`

1. `plot_j_spin_components` lives in `quantum_trajectories/plotting_j_moments.py`.

2. The function should take the moment series directly, e.g. `moments.J`, not
   the full `MomentSeries` container.

3. The default call plots J moments:

   ```python
   fig, axes = plot_j_spin_components(moments.J)
   ```

4. The function should support `spin_component="j"` and `spin_component="s"`.
   For `spin_component="j"`, it reads `x`, `y`, `z`, `length`, and optional
   `x_groups`, `y_groups`, `z_groups`, `length_groups`. For
   `spin_component="s"`, it should read the analogous normalized-direction
   fields `nx`, `ny`, `nz`, plus `length`.

5. The output should be a `2x2` panel showing the x, y, z, and length
   components.

6. Each panel should include the full-system curve. If group fields exist, each
   panel should also include all group-resolved curves.

7. The function should support `axes`, `output_path`, `label`, and `phases` so
   multiple results can be plotted into the same figure.

## `fig, axes = plot_j_angles(...)`

1. `plot_j_angles` lives in `quantum_trajectories/plotting_j_moments.py`.

2. The function should take `moments.J` as input.

3. The output should be a `2x1` panel showing `theta` and `phi`.

4. `plot_j_angles` should read stored `theta` and `phi` fields from
   `moments.J`. Angle construction should happen upstream in the moment
   pipeline; the plotting function should not recompute angles from
   `x/y/z`, `length`, or `nx/ny/nz`.

5. If group-resolved `theta_groups` and `phi_groups` exist, plot group angles
   first using dashed group-colored curves. Plot the full-system angles last
   using a solid gray curve.

6. The function should support `axes`, `output_path`, `label`, and `phases`.

## `fig, axes = plot_mfe_residuals(...)`

1. `plot_mfe_residuals` lives in
   `quantum_trajectories/plotting_diagnostics.py`.

2. The function should take `moments.J` as input.

3. The function should only plot already-computed residuals from
   `moments.J.mfe_residuals_groups`. It should not recompute MFE residuals or
   call old observable extraction.

4. The output should be a single panel showing `Re R_1`, `Im R_1`,
   `Re R_2`, `Im R_2`, and the L2 norm `sqrt(|R_1|^2 + |R_2|^2)`.

5. The signed residual components should be dashed. The two `R_1` components
   should use different blue shades, and the two `R_2` components should use
   different orange shades. The L2 norm should use a solid gray curve.

6. The function should support `axes`, `output_path`, `label`, and `phases`.
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
   `quantum_trajectories/plotting_diagnostics.py`.


Legacy note: the previous J-moment field names were `Jx`, `Jy`, `Jz`,
`Jx_groups`, `Jy_groups`, `Jz_groups`, `J_len`, and `sx`, `sy`, `sz`.

## Global Styling Rules

Shared styling helpers should live in `quantum_trajectories/plotting_utils.py`.

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
