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
            ("Jx", "Jy", "Jz"),
            ("Jx_groups", "Jy_groups", "Jz_groups"),
            "J spin components",
        )
    if spin_component == "s":
        return (
            ("sx", "sy", "sz"),
            ("sx_groups", "sy_groups", "sz_groups"),
            "S spin directions",
        )
    raise ValueError("spin_component must be 'j' or 's'.")


def _curve_label(base_label: str, *, label: Optional[str]) -> str:
    if label is None:
        return base_label
    return f"{label} {base_label}"


def _angles_from_j_components(
    jx: Any,
    jy: Any,
    jz: Any,
    *,
    tol: float,
) -> tuple[np.ndarray, np.ndarray]:
    jx = np.asarray(jx, dtype=float)
    jy = np.asarray(jy, dtype=float)
    jz = np.asarray(jz, dtype=float)

    j_len = np.sqrt(jx**2 + jy**2 + jz**2)
    valid = j_len > tol

    sx = np.zeros_like(jx, dtype=float)
    sy = np.zeros_like(jy, dtype=float)
    sz = np.zeros_like(jz, dtype=float)
    sx[valid] = jx[valid] / j_len[valid]
    sy[valid] = jy[valid] / j_len[valid]
    sz[valid] = jz[valid] / j_len[valid]

    theta = np.zeros_like(sz, dtype=float)
    theta[valid] = np.arccos(np.clip(-sz[valid], -1.0, 1.0))

    phi = np.arctan2(sy, sx)
    r_perp = np.sqrt(sx**2 + sy**2)
    phi[r_perp < tol] = 0.0

    return theta, phi


def plot_j_spin_components(
    spin_moments: Any,
    *,
    spin_component: str = "j",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
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
        group_values = getattr(spin_moments, group_field_name, None)
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
    tol: float = 1e-12,
):
    """
    Plot polar and azimuthal angles from the direction of the mean J vector.

    The normalized direction is computed internally as J_i / |J|, so the angles
    are tied directly to the Jx, Jy, Jz moments supplied by the input container.
    """
    if axes is None:
        fig, axes = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes).ravel()
    if axes.size != 2:
        raise ValueError("axes must contain exactly two axes for the 2x1 angle grid.")

    t = np.asarray(j_moments.t, dtype=float)
    group_fields = (
        getattr(j_moments, "Jx_groups", None),
        getattr(j_moments, "Jy_groups", None),
        getattr(j_moments, "Jz_groups", None),
    )

    if all(field is not None for field in group_fields):
        jx_groups, jy_groups, jz_groups = group_fields
        for group_index, (jx_g, jy_g, jz_g) in enumerate(
            zip(jx_groups, jy_groups, jz_groups),
            start=1,
        ):
            theta_g, phi_g = _angles_from_j_components(jx_g, jy_g, jz_g, tol=tol)
            group_color = _GROUP_CURVE_COLORS[(group_index - 1) % len(_GROUP_CURVE_COLORS)]
            axes[0].plot(
                t,
                theta_g,
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=_curve_label(rf"$\theta_{group_index}$", label=label),
            )
            axes[1].plot(
                t,
                phi_g,
                linewidth=1.8,
                color=group_color,
                linestyle="--",
                label=_curve_label(rf"$\phi_{group_index}$", label=label),
            )

    theta, phi = _angles_from_j_components(j_moments.Jx, j_moments.Jy, j_moments.Jz, tol=tol)
    axes[0].plot(
        t,
        theta,
        linewidth=1.8,
        color=_FULL_CURVE_COLOR,
        linestyle="-",
        label=_curve_label(r"$\theta$", label=label),
    )
    axes[1].plot(
        t,
        phi,
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

    fig.suptitle("J-vector angles")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes
