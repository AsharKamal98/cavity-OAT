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

5. Full-system curves should use a solid gray line.

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

## `fig, axes = plot_j_spin_components(...)`

1. `plot_j_spin_components` lives in `quantum_trajectories/plotting_j_moments.py`.

2. The function should take the moment series directly, e.g. `moments.J`, not
   the full `MomentSeries` container.

3. The default call plots J moments:

   ```python
   fig, axes = plot_j_spin_components(moments.J)
   ```

4. The function should support `spin_component="j"` and `spin_component="s"`.
   For `spin_component="j"`, it reads `Jx`, `Jy`, `Jz` and optional
   `Jx_groups`, `Jy_groups`, `Jz_groups`. For `spin_component="s"`, it should
   read the analogous `sx`, `sy`, `sz` fields once they exist.

5. The output should be a `3x1` panel showing the x, y, and z components.

6. Each panel should include the full-system curve. If group fields exist, each
   panel should also include all group-resolved curves.

7. The function should support `axes`, `output_path`, and `label` so multiple
   results can be plotted into the same figure.

## `fig, axes = plot_j_angles(...)`

1. `plot_j_angles` lives in `quantum_trajectories/plotting_j_moments.py`.

2. The function should take `moments.J` as input.

3. The output should be a `2x1` panel showing `theta` and `phi`.

4. Angles should be computed from the direction of the mean J vector. For each
   full-system or group-resolved curve, first compute

   ```python
   J_len = sqrt(Jx**2 + Jy**2 + Jz**2)
   sx = Jx / J_len
   sy = Jy / J_len
   sz = Jz / J_len
   ```

   for points with `J_len > tol`.

5. The angle convention should match the old active-manifold angle convention:

   ```python
   theta = arccos(-sz)
   phi = atan2(sy, sx)
   ```

   If the transverse length is below tolerance, set `phi = 0`.

6. If group-resolved J fields exist, plot group angles first using dashed
   group-colored curves. Plot the full-system angles last using a solid gray
   curve.

7. The function should support `axes`, `output_path`, `label`, and `tol`.
