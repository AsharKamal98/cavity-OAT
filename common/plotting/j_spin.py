from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.utils import (
    colour_palette,
    finish_time_plot,
    format_time_axis,
    get_axes,
    palette_curve_color,
    style_axis,
    validated_linestyle,
)
from common.utils.moments import as_series_tuple
from parser.common import Array, PhaseProtocol


def plot_spin_components(
    t,
    x_components: Sequence[Array] | Array,
    y_components: Sequence[Array] | Array,
    z_components: Sequence[Array] | Array,
    lengths: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str = "-",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
    title: str = "Spin components",
):
    """
    Plot supplied x/y/z/length time series on a 2x2 grid.
    """
    fig, axes = get_axes(
        axes,
        n_axes=4,
        create_figure=lambda: plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly four axes for the 2x2 spin-component grid.",
    )
    line_style = validated_linestyle(linestyle)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

    t = np.asarray(t, dtype=float)
    component_sets = tuple(
        as_series_tuple(components)
        for components in (x_components, y_components, z_components, lengths)
    )
    labels = tuple(labels)
    if not labels:
        raise ValueError("plot_spin_components requires at least one spin vector.")
    if any(len(components) != len(labels) for components in component_sets):
        raise ValueError("Each component input requires one curve per label.")

    axis_names = ("x", "y", "z", "len")
    panel_titles = (r"$J_x$", r"$J_y$", r"$J_z$", r"$|J|$")

    for ax, components, axis_name, panel_title in zip(
        axes,
        component_sets,
        axis_names,
        panel_titles,
    ):
        component_label = r"$|J|$" if axis_name == "len" else rf"$J_{{{axis_name}}}$"
        for curve_index, (values, label) in enumerate(zip(components, labels)):
            ax.plot(
                t,
                values,
                linewidth=1.8,
                color=palette_curve_color(palette, curve_index),
                linestyle=line_style,
                label=label,
            )

        ax.set_ylabel(component_label)
        ax.set_title(panel_title, fontsize=11)
        style_axis(ax)
        ax.legend()
        format_time_axis(ax)
    finish_time_plot(
        fig,
        axes,
        phase_protocol=phase_protocol,
        title=title,
        output_path=output_path,
    )

    return fig, axes


def plot_bloch_angles(
    t,
    theta_curves: Sequence[Array] | Array,
    phi_curves: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str = "-",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
):
    """
    Plot supplied polar and azimuthal angle time series on a 2x1 grid.
    """
    fig, axes = get_axes(
        axes,
        n_axes=2,
        create_figure=lambda: plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly two axes for the 2x1 angle grid.",
    )
    line_style = validated_linestyle(linestyle)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

    t = np.asarray(t, dtype=float)
    theta_curves = as_series_tuple(theta_curves)
    phi_curves = as_series_tuple(phi_curves)
    labels = tuple(labels)
    if not labels:
        raise ValueError("plot_bloch_angles requires at least one Bloch vector.")
    if len(theta_curves) != len(labels) or len(phi_curves) != len(labels):
        raise ValueError("Theta and phi inputs require one curve per label.")

    for curve_index, (theta, phi, label) in enumerate(
        zip(theta_curves, phi_curves, labels)
    ):
        color = palette_curve_color(palette, curve_index)
        axes[0].plot(
            t,
            theta,
            linewidth=1.8,
            color=color,
            linestyle=line_style,
            label=label,
        )
        axes[1].plot(
            t,
            phi,
            linewidth=1.8,
            color=color,
            linestyle=line_style,
            label=label,
        )

    angle_specs = (
        (axes[0], r"$\theta$", r"Polar angle $\theta(t)$"),
        (axes[1], r"$\phi$", r"Azimuthal angle $\phi(t)$"),
    )
    for ax, ylabel, title in angle_specs:
        ax.set_ylabel(ylabel)
        ax.set_title(title, fontsize=11)
        style_axis(ax)
        ax.legend()
        format_time_axis(ax)

    finish_time_plot(
        fig,
        axes,
        phase_protocol=phase_protocol,
        title="Bloch angles",
        output_path=output_path,
    )

    return fig, axes
