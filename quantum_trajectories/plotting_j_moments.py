from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

_FULL_CURVE_COLOR = "gray"
_GROUP_CURVE_COLORS = ("tab:blue", "tab:orange", "tab:green", "tab:red")


def _spin_component_fields(spin_component: str) -> tuple[tuple[str, str, str], tuple[str, str, str], str]:
    spin_component = spin_component.lower()
    if spin_component == "j":
        return (
            ("x", "y", "z"),
            ("x_groups", "y_groups", "z_groups"),
            "J spin components",
        )
    if spin_component == "s":
        return (
            ("nx", "ny", "nz"),
            ("nx_groups", "ny_groups", "nz_groups"),
            "Normalized spin directions",
        )
    raise ValueError("spin_component must be 'j' or 's'.")


def _curve_label(base_label: str, *, label: Optional[str]) -> str:
    if label is None:
        return base_label
    return f"{label} {base_label}"


def _add_phase_boundaries(axes, phases) -> None:
    if phases is None:
        return

    boundaries = np.cumsum([phase.duration for phase in phases], dtype=float)[:-1]
    for ax in np.asarray(axes).ravel():
        existing = [
            line.get_xdata()[0]
            for line in ax.lines
            if line.get_gid() == "phase_boundary"
        ]
        for boundary in boundaries:
            if any(np.isclose(boundary, x) for x in existing):
                continue
            line = ax.axvline(boundary, linestyle="--", color="black", alpha=0.6)
            line.set_gid("phase_boundary")


def plot_j_spin_components(
    spin_moments: Any,
    *,
    spin_component: str = "j",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot full and optional group-resolved spin components on a 1x3 grid.

    Pass ``moments.J`` for the current J-moment pipeline. The ``spin_component``
    option is kept generic so the same function can later plot S-direction
    fields once those exist on the input container.
    """
    component_fields, group_fields, title = _spin_component_fields(spin_component)
    spin_symbol = spin_component.upper() if spin_component.lower() == "j" else spin_component.lower()

    if axes is None:
        fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes).ravel()
    if axes.size != 3:
        raise ValueError("axes must contain exactly three axes for the 3x1 spin-component grid.")

    t = np.asarray(spin_moments.t, dtype=float)
    axis_names = ("x", "y", "z")

    for ax, field_name, group_field_name, axis_name in zip(
        axes,
        component_fields,
        group_fields,
        axis_names,
    ):
        component_label = rf"${spin_symbol}_{axis_name}$"
        group_values = getattr(spin_moments, group_field_name)
        if group_values is not None:
            for group_index, values in enumerate(group_values, start=1):
                group_label = rf"${spin_symbol}_{{{axis_name},{group_index}}}$"
                group_color = _GROUP_CURVE_COLORS[(group_index - 1) % len(_GROUP_CURVE_COLORS)]
                ax.plot(
                    t,
                    np.asarray(values, dtype=float),
                    linewidth=1.8,
                    color=group_color,
                    linestyle="--",
                    label=_curve_label(group_label, label=label),
                )

        ax.plot(
            t,
            np.asarray(getattr(spin_moments, field_name), dtype=float),
            linewidth=1.8,
            color=_FULL_CURVE_COLOR,
            linestyle="-",
            label=_curve_label(component_label, label=label),
        )

        ax.set_xlabel(r"$\Gamma t$")
        ax.set_ylabel(component_label)
        ax.set_title(component_label)
        ax.grid(alpha=0.3)
        ax.legend()
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)

    _add_phase_boundaries(axes, phases)

    fig.suptitle(title)
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes


def plot_j_angles(
    j_moments: Any,
    *,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot stored polar and azimuthal angles.

    The angle arrays must already be stored on the input moment series.
    """
    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes).ravel()
    if axes.size != 2:
        raise ValueError("axes must contain exactly two axes for the 1x2 angle grid.")

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
            group_color = _GROUP_CURVE_COLORS[(group_index - 1) % len(_GROUP_CURVE_COLORS)]
            axes[0].plot(
                t,
                np.asarray(theta_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=_curve_label(rf"$\theta_{group_index}$", label=label),
            )
            axes[1].plot(
                t,
                np.asarray(phi_g, dtype=float),
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=_curve_label(rf"$\phi_{group_index}$", label=label),
            )

    axes[0].plot(
        t,
        np.asarray(theta, dtype=float),
        linewidth=1.8,
        color=_FULL_CURVE_COLOR,
        linestyle="-",
        label=_curve_label(r"$\theta$", label=label),
    )
    axes[1].plot(
        t,
        np.asarray(phi, dtype=float),
        linewidth=1.8,
        color=_FULL_CURVE_COLOR,
        linestyle="-",
        label=_curve_label(r"$\phi$", label=label),
    )

    angle_specs = (
        (axes[0], r"$\theta$", r"Polar $\theta(t)$"),
        (axes[1], r"$\phi$", r"Azimuthal $\phi(t)$"),
    )
    for ax, ylabel, title in angle_specs:
        ax.set_xlabel(r"$\Gamma t$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(alpha=0.3)
        ax.legend()
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)

    _add_phase_boundaries(axes, phases)

    fig.suptitle("J-vector angles")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes
