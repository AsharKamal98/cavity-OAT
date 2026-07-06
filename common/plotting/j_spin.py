from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.utils import (
    curve_label,
    finish_time_plot,
    format_time_axis,
    get_axes,
    indexed_curve_color,
    style_axis,
    validated_linestyle,
)


def plot_spin_components(
    spin_series: Any,
    normalized: bool = False,
    *,
    colour_index: int = 0,
    linestyle: str = "-",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot any available full and group-resolved x/y/z/length series on a 2x2 grid.
    """
    fig, axes = get_axes(
        axes,
        n_axes=4,
        create_figure=lambda: plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly four axes for the 2x2 spin-component grid.",
    )
    line_style = validated_linestyle(linestyle)

    t = np.asarray(spin_series.t, dtype=float)
    if not normalized:
        component_fields = ("x", "y", "z", "length")
        group_fields = ("x_groups", "y_groups", "z_groups", "length_groups")
    else:
        component_fields = ("nx", "ny", "nz", "length")
        group_fields = ("nx_groups", "ny_groups", "nz_groups", "length_groups")
    axis_names = ("x", "y", "z", "len")
    panel_titles = (r"$J_x$", r"$J_y$", r"$J_z$", r"$|J|$")

    for ax, field_name, group_field_name, axis_name, panel_title in zip(
        axes,
        component_fields,
        group_fields,
        axis_names,
        panel_titles,
    ):
        values = getattr(spin_series, field_name, None)
        group_values = getattr(spin_series, group_field_name, None)

        if values is None and group_values is None:
            raise ValueError(
                f"plot_spin_components requires {field_name} or {group_field_name} on the input series."
            )

        component_label = r"$|J|$" if axis_name == "len" else rf"$J_{{{axis_name}}}$"

        if group_values is not None:
            for group_index, group_data in enumerate(group_values, start=1):
                group_label = (
                    rf"$|J_{group_index}|$"
                    if axis_name == "len"
                    else rf"$J_{{{axis_name},{group_index}}}$"
                )
                group_color = indexed_curve_color(colour_index, group_index - 1)
                ax.plot(
                    t,
                    np.asarray(group_data, dtype=float),
                    linewidth=1.8,
                    color=group_color,
                    linestyle=line_style,
                    label=curve_label(group_label, label=label),
                )

        # if values is not None:
        if False:
            full_color = indexed_curve_color(
                colour_index,
                len(group_values) if group_values is not None else 0,
            )
            ax.plot(
                t,
                np.asarray(values, dtype=float),
                linewidth=1.8,
                color=full_color,
                linestyle=line_style,
                label=curve_label(component_label, label=label),
            )

        ax.set_ylabel(component_label)
        ax.set_title(panel_title, fontsize=11)
        style_axis(ax)
        ax.legend()
        format_time_axis(ax)
    title = "Normalized spin components" if normalized else "Spin components"
    finish_time_plot(fig, axes, phases=phases, title=title, output_path=output_path)

    return fig, axes


def plot_bloch_angles(
    angle_series: Any,
    *,
    colour_index: int = 0,
    linestyle: str = "-",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot any available full and group-resolved Bloch-angle series stored on an input object.
    """
    fig, axes = get_axes(
        axes,
        n_axes=2,
        create_figure=lambda: plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly two axes for the 2x1 angle grid.",
    )
    line_style = validated_linestyle(linestyle)

    t = np.asarray(angle_series.t, dtype=float)
    theta = getattr(angle_series, "theta", None)
    phi = getattr(angle_series, "phi", None)
    theta_groups = getattr(angle_series, "theta_groups", None)
    phi_groups = getattr(angle_series, "phi_groups", None)

    has_full = theta is not None and phi is not None
    has_groups = theta_groups is not None and phi_groups is not None
    if not has_full and not has_groups:
        raise ValueError(
            "plot_bloch_angles requires full theta/phi fields or group-resolved theta_groups/phi_groups."
        )

    if has_groups:
        for group_index, (theta_g, phi_g) in enumerate(zip(theta_groups, phi_groups), start=1):
            group_color = indexed_curve_color(colour_index, group_index - 1)
            axes[0].plot(
                t,
                np.asarray(theta_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle=line_style,
                label=curve_label(rf"$\theta_{group_index}$", label=label),
            )
            axes[1].plot(
                t,
                np.asarray(phi_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle=line_style,
                label=curve_label(rf"$\phi_{group_index}$", label=label),
            )

    if False:
        full_color = indexed_curve_color(
            colour_index,
            len(theta_groups) if has_groups else 0,
        )
        axes[0].plot(
            t,
            np.asarray(theta, dtype=float),
            linewidth=1.8,
            color=full_color,
            linestyle=line_style,
            label=curve_label(r"$\theta$", label=label),
        )
        axes[1].plot(
            t,
            np.asarray(phi, dtype=float),
            linewidth=1.8,
            color=full_color,
            linestyle=line_style,
            label=curve_label(r"$\phi$", label=label),
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

    finish_time_plot(fig, axes, phases=phases, title="Bloch angles", output_path=output_path)

    return fig, axes
