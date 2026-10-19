from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgb

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


def _validated_sem_curves(
    sem: Sequence[Array | None] | Array | None,
    *,
    curves: tuple[Array, ...],
    t: Array,
    name: str,
) -> tuple[Array | None, ...]:
    """Match nonnegative SEM curves; NaN omits undefined uncertainty points."""
    if sem is None:
        return (None,) * len(curves)
    try:
        entries = tuple(sem)
    except TypeError as error:
        raise ValueError(f"{name} must be a time series or sequence of SEM curves.") from error
    if any(np.iscomplexobj(entry) for entry in entries if entry is not None):
        raise ValueError(f"{name} must be real.")
    errors = entries if any(entry is None for entry in entries) else as_series_tuple(sem)
    if len(errors) != len(curves):
        raise ValueError(f"{name} requires one SEM curve per plotted curve.")
    result = []
    for error, values in zip(errors, curves):
        if error is None:
            result.append(None)
            continue
        error = np.asarray(error)
        if np.iscomplexobj(error):
            raise ValueError(f"{name} must be real.")
        error = np.asarray(error, dtype=float)
        if error.ndim != 1 or error.shape != values.shape or error.shape != t.shape:
            raise ValueError(f"{name} must match the central curves and time grid.")
        if np.any(np.isinf(error)) or np.any(error < 0.0):
            raise ValueError(f"{name} must be nonnegative, with NaN only for undefined points.")
        result.append(error)
    return tuple(result)


def plot_spin_components(
    t,
    x_components: Sequence[Array] | Array,
    y_components: Sequence[Array] | Array,
    z_components: Sequence[Array] | Array,
    lengths: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    x_sem: Sequence[Array | None] | Array | None = None,
    y_sem: Sequence[Array | None] | Array | None = None,
    z_sem: Sequence[Array | None] | Array | None = None,
    length_sem: Sequence[Array | None] | Array | None = None,
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str | None = "-",
    marker: str | None = None,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
    title: str = "Spin components",
    spin_symbol: str = "J",
):
    """
    Plot supplied x/y/z/length time series with optional mean +/- SEM bands.

    SEM inputs follow the corresponding curve order and may be omitted for
    the whole panel or individual curves. No uncertainties are computed here.
    """
    fig, axes = get_axes(
        axes,
        n_axes=4,
        create_figure=lambda: plt.subplots(2, 2, figsize=(10, 6.5), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly four axes for the 2x2 spin-component grid.",
    )
    line_style = validated_linestyle(linestyle)
    marker_style = validated_marker(marker)
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
    sem_sets = tuple(
        _validated_sem_curves(errors, curves=curves, t=t, name=name)
        for errors, curves, name in zip(
            (x_sem, y_sem, z_sem, length_sem),
            component_sets,
            ("x_sem", "y_sem", "z_sem", "length_sem"),
        )
    )

    axis_names = ("x", "y", "z", "len")
    panel_titles = (
        rf"${spin_symbol}_x$",
        rf"${spin_symbol}_y$",
        rf"${spin_symbol}_z$",
        rf"$|{spin_symbol}|$",
    )

    for ax, components, errors, axis_name, panel_title in zip(
        axes,
        component_sets,
        sem_sets,
        axis_names,
        panel_titles,
    ):
        component_label = (
            rf"$|{spin_symbol}|$"
            if axis_name == "len"
            else rf"${spin_symbol}_{{{axis_name}}}$"
        )
        for curve_index, (values, error, label) in enumerate(zip(components, errors, labels)):
            line, = ax.plot(
                t,
                values,
                linewidth=1.8,
                color=palette_curve_color(palette, curve_index),
                linestyle=line_style,
                marker=marker_style,
                label=label,
            )
            if error is not None:
                band_color = 0.75 * np.asarray(to_rgb(line.get_color()))
                ax.fill_between(
                    t,
                    values - error,
                    values + error,
                    color=band_color,
                    alpha=0.45,
                    linewidth=0.0,
                    zorder=line.get_zorder() - 0.5,
                    label="_nolegend_",
                )

        ax.set_ylabel(component_label)
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


def plot_bloch_angles(
    t,
    theta_curves: Sequence[Array] | Array,
    phi_curves: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str | None = "-",
    marker: str | None = None,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
    title: str = "Bloch angles",
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
    marker_style = validated_marker(marker)
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
            marker=marker_style,
            label=label,
        )
        axes[1].plot(
            t,
            phi,
            linewidth=1.8,
            color=color,
            linestyle=line_style,
            marker=marker_style,
            label=label,
        )

    angle_specs = (
        (axes[0], r"$\theta$", r"Polar angle $\theta(t)$"),
        (axes[1], r"$\phi$", r"Azimuthal angle $\phi(t)$"),
    )
    for ax, ylabel, panel_title in angle_specs:
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
