from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.utils import phase1_ss_angles_for_nj

from quantum_trajectories.plotting_utils import (
    GROUP_CURVE_COLORS,
    curve_label,
    finish_time_plot,
    format_time_axis,
    full_curve_color,
    get_axes,
    style_axis,
)


def _spin_component_fields(spin_component: str) -> tuple[tuple[str, str, str], tuple[str, str, str], str]:
    spin_component = spin_component.lower()
    if spin_component == "j":
        return (
            ("x", "y", "z", "length"),
            ("x_groups", "y_groups", "z_groups", "length_groups"),
            "J spin components",
        )
    if spin_component == "s":
        return (
            ("nx", "ny", "nz", "length"),
            ("nx_groups", "ny_groups", "nz_groups", "length_groups"),
            "Normalized spin directions",
        )
    raise ValueError("spin_component must be 'j' or 's'.")


def _spin_symbol(spin_component: str) -> str:
    return "J" if spin_component.lower() == "j" else "n"


def plot_j_spin_components(
    spin_moments: Any,
    *,
    spin_component: str = "j",
    colour_index: int = 0,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot full and optional group-resolved spin components on a 2x2 grid.

    Pass ``moments.J`` for the current J-moment pipeline. The ``spin_component``
    option is kept generic so the same function can later plot S-direction
    fields once those exist on the input container.
    """
    component_fields, group_fields, title = _spin_component_fields(spin_component)
    spin_symbol = _spin_symbol(spin_component)
    full_color = full_curve_color(colour_index)

    fig, axes = get_axes(
        axes,
        n_axes=4,
        create_figure=lambda: plt.subplots(2, 2, figsize=(10, 7), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly four axes for the 2x2 spin-component grid.",
    )

    t = np.asarray(spin_moments.t, dtype=float)
    axis_names = ("x", "y", "z", "len")
    panel_titles = (
        rf"${spin_symbol}_x$",
        rf"${spin_symbol}_y$",
        rf"${spin_symbol}_z$",
        r"$|J|$",
    )

    for ax, field_name, group_field_name, axis_name, panel_title in zip(
        axes,
        component_fields,
        group_fields,
        axis_names,
        panel_titles,
    ):
        component_label = rf"$|{spin_symbol}|$" if axis_name == "len" else rf"${spin_symbol}_{{{axis_name}}}$"
        group_values = getattr(spin_moments, group_field_name)
        if group_values is not None:
            for group_index, values in enumerate(group_values, start=1):
                group_label = (
                    rf"$|{spin_symbol}_{group_index}|$"
                    if axis_name == "len"
                    else rf"${spin_symbol}_{{{axis_name},{group_index}}}$"
                )
                group_color = GROUP_CURVE_COLORS[(group_index - 1) % len(GROUP_CURVE_COLORS)]
                ax.plot(
                    t,
                    np.asarray(values, dtype=float),
                    linewidth=1.8,
                    color=group_color,
                    linestyle="--",
                    label=curve_label(group_label, label=label),
                )

        ax.plot(
            t,
            np.asarray(getattr(spin_moments, field_name), dtype=float),
            linewidth=1.8,
            color=full_color,
            linestyle="-",
            label=curve_label(component_label, label=label),
        )

        ax.set_ylabel(component_label)
        ax.set_title(panel_title, fontsize=11)
        style_axis(ax)
        ax.legend()
        format_time_axis(ax)

    finish_time_plot(fig, axes, phases=phases, title=title, output_path=output_path)

    return fig, axes


def plot_j_angles(
    j_moments: Any,
    *,
    colour_index: int = 0,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
    show_phase1_ss: bool = False,
    Gamma: Optional[float] = None,
):
    """
    Plot stored polar and azimuthal angles.

    The angle arrays must already be stored on the input moment series.
    """
    fig, axes = get_axes(
        axes,
        n_axes=2,
        create_figure=lambda: plt.subplots(2, 1, figsize=(8, 6), sharex=True, constrained_layout=True),
        error_message="axes must contain exactly two axes for the 2x1 angle grid.",
    )
    full_color = full_curve_color(colour_index)

    t = np.asarray(j_moments.t, dtype=float)
    theta = j_moments.theta
    phi = j_moments.phi
    if theta is None or phi is None:
        raise ValueError("plot_j_angles requires j_moments.theta and j_moments.phi.")

    theta_groups = j_moments.theta_groups
    phi_groups = j_moments.phi_groups
    if theta_groups is not None and phi_groups is not None:
        for group_index, (theta_g, phi_g) in enumerate(
            zip(theta_groups, phi_groups),
            start=1,
        ):
            group_color = GROUP_CURVE_COLORS[(group_index - 1) % len(GROUP_CURVE_COLORS)]
            axes[0].plot(
                t,
                np.asarray(theta_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=curve_label(rf"$\theta_{group_index}$", label=label),
            )
            axes[1].plot(
                t,
                np.asarray(phi_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=curve_label(rf"$\phi_{group_index}$", label=label),
            )

    axes[0].plot(
        t,
        np.asarray(theta, dtype=float),
        linewidth=1.8,
        color=full_color,
        linestyle="-",
        label=curve_label(r"$\theta$", label=label),
    )
    axes[1].plot(
        t,
        np.asarray(phi, dtype=float),
        linewidth=1.8,
        color=full_color,
        linestyle="-",
        label=curve_label(r"$\phi$", label=label),
    )

    # FIXME: should be moved to moment claculations in the future
    if show_phase1_ss:
        if phases is None:
            raise ValueError("plot_j_angles requires phases when show_phase1_ss=True.")
        if Gamma is None:
            raise ValueError("plot_j_angles requires Gamma when show_phase1_ss=True.")
        n_j = np.asarray(j_moments.N_j, dtype=float)
        finite_nj = n_j[np.isfinite(n_j)]
        if finite_nj.size == 0:
            raise ValueError("plot_j_angles requires finite j_moments.N_j values.")
        theta_ss, _ = phase1_ss_angles_for_nj(float(finite_nj[0]), phases[0].omega, Gamma)
        if np.isfinite(theta_ss):
            axes[0].hlines(
                y=theta_ss,
                xmin=0.0,
                xmax=phases[0].duration,
                linewidth=1.8,
                color="#009E73",
                linestyle=":",
                label=curve_label(r"phase 1 ss", label=label),
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

    finish_time_plot(fig, axes, phases=phases, title="J-vector angles", output_path=output_path)

    return fig, axes
