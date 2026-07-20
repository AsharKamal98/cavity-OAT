from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Literal, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.utils import (
    colour_palette,
    finish_time_plot,
    format_time_axis,
    get_axes,
    palette_curve_color,
    set_bottom_figure_legend,
    style_axis,
    validated_linestyle,
    validated_marker,
)
from common.utils.moments import as_series_tuple
from parser.common import Array, PhaseProtocol


def plot_transverse_spin_components(
    t,
    x_components: Sequence[Array] | Array,
    y_components: Sequence[Array] | Array,
    z_components: Sequence[Array] | Array,
    lengths: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    transverse_definition: Literal["magnitude", "sum"] = "magnitude",
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str | None = "-",
    marker: str | None = None,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
    title: str = "Spin components",
):
    """Plot a selected transverse component, longitudinal component, and spin length."""

    fig, axes = get_axes(
        axes,
        n_axes=3,
        create_figure=lambda: plt.subplots(
            3,
            1,
            figsize=(8, 8),
            sharex=True,
            constrained_layout=True,
        ),
        error_message="axes must contain exactly three axes for the 3x1 spin-component grid.",
    )
    line_style = validated_linestyle(linestyle)
    marker_style = validated_marker(marker)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

    t = np.asarray(t, dtype=float)
    x_components, y_components, z_components, lengths = (
        as_series_tuple(components)
        for components in (x_components, y_components, z_components, lengths)
    )
    labels = tuple(labels)
    component_sets = (x_components, y_components, z_components, lengths)
    if not labels:
        raise ValueError(
            "plot_transverse_spin_components requires at least one spin vector."
        )
    if any(len(components) != len(labels) for components in component_sets):
        raise ValueError("Each component input requires one curve per label.")

    if transverse_definition == "magnitude":
        transverse_components = tuple(
            np.hypot(x, y) for x, y in zip(x_components, y_components)
        )
        transverse_label = r"$J_\perp$"
    elif transverse_definition == "sum":
        transverse_components = tuple(
            x + y for x, y in zip(x_components, y_components)
        )
        transverse_label = r"$J_x+J_y$"
    else:
        raise ValueError("transverse_definition must be 'magnitude' or 'sum'.")

    for curve_index, (transverse, z, length, label) in enumerate(
        zip(transverse_components, z_components, lengths, labels)
    ):
        color = palette_curve_color(palette, curve_index)
        curve_style = {
            "linewidth": 1.8,
            "color": color,
            "linestyle": line_style,
            "marker": marker_style,
            "label": label,
        }
        axes[0].plot(t, transverse, **curve_style)
        axes[1].plot(t, z, **curve_style)
        axes[2].plot(t, length, **curve_style)

    panel_specs = (
        (axes[0], transverse_label, transverse_label),
        (axes[1], r"$J_z$", r"$J_z$"),
        (axes[2], r"$|J|$", r"$|J|$"),
    )
    for ax, ylabel, panel_title in panel_specs:
        ax.set_ylabel(ylabel)
        ax.set_title(panel_title, fontsize=11)
        style_axis(ax)
        format_time_axis(ax)

    set_bottom_figure_legend(fig, axes[0])
    finish_time_plot(
        fig,
        axes,
        phase_protocol=phase_protocol,
        title=title,
        output_path=output_path,
    )

    return fig, axes


__all__ = ["plot_transverse_spin_components"]
