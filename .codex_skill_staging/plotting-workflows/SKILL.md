---
name: plotting-workflows
description: Create or modify plotting functions, shared plotting helpers, or plotting instructions. Use when the task involves plot styling, overlay conventions, panel layouts, legends, phase boundaries, scientific axes, color handling, or requests like "use colour index for this plot".
---

# Plotting Workflows

Use this skill when implementing or revising plotting code.

## Core Rules

- Keep plotting functions thin: plot already-computed data; do not recompute
  expensive physics inside plotting unless the user explicitly wants that.
- Prefer shared plotting helpers for backend-neutral styling, colors, axis
  formatting, phase shading, and save logic.
- Plot group-resolved curves first and full-system curves last.
- Group-resolved curves should usually be dashed; full-system curves should
  usually be solid.
- Use scientific notation on time axes when labels would crowd, and disable
  additive offset notation.
- If a plot supports `axes`, preserve overlay behavior:
  - create new axes only when `axes is None`;
  - otherwise plot into the provided axes and return the same figure/axes.
- Keep panel titles, labels, and legends short and physically meaningful.

## Shared-vs-Local Split

- Put backend-neutral plotting style helpers in a shared module such as
  `common/plotting_utils.py`.
- Put generic plotting functions that can work across multiple data series in a
  shared plotting module such as `common/plotting.py`.
- Keep backend-specific plotting composition local to the relevant package.

## When Asked To Use `colour_index`

If the user asks to use `colour_index` for a plot:

- apply `colour_index` only to the full-system curve color, unless the user
  explicitly asks for group curves to vary too;
- keep group-resolved curves on their fixed group palette by default;
- implement the full-curve color lookup through one shared helper, not by
  hardcoding colors in each plotting function.

If the plotted series has no full-system curve, `colour_index` will have no
effect unless the user explicitly asks for alternate behavior.

## Phase Styling

When `phases` are available:

- use shared helpers to add phase shading and vertical phase boundaries;
- avoid duplicating manual boundary logic across plotting functions.

## Documentation

If the plotting behavior changes intentionally, update the relevant local
plotting instruction file as well.
