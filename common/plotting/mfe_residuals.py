from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.utils import (
    colour_palette,
    curve_label,
    finish_time_plot,
    format_time_axis,
    get_axes,
    palette_curve_color,
    style_axis,
    validated_linestyle,
)
from common.utils.phases import phase_boundary_times
from parser.common import PhaseProtocol

RESIDUAL_LINTHRESH = 1e-5


def _set_symmetric_ylim_from_lines(ax) -> None:
    finite_values = []
    for line in ax.lines:
        if line.get_gid() == "phase_boundary":
            continue
        ydata = np.asarray(line.get_ydata(), dtype=float)
        ydata = ydata[np.isfinite(ydata)]
        if ydata.size:
            finite_values.append(ydata)

    if not finite_values:
        ax.set_ylim(-1.0, 1.0)
        return

    max_abs = float(np.max(np.abs(np.concatenate(finite_values))))
    if max_abs <= 0.0:
        max_abs = 1.0
    ax.set_ylim(-max_abs, max_abs)


def _set_residual_scale(ax, *, show_components: bool, symlog: bool, residual_l2) -> None:
    if not symlog:
        ax.set_yscale("linear")
    elif show_components or np.any(np.asarray(residual_l2) <= 0.0):
        ax.set_yscale("symlog", linthresh=RESIDUAL_LINTHRESH)
    else:
        ax.set_yscale("log")


def _print_phase_end_residuals(
    mfe_residuals: Any,
    phase_protocol: PhaseProtocol | None,
) -> None:
    if phase_protocol is None:
        return
    residuals = mfe_residuals.residuals_groups
    if len(residuals) != 2:
        return

    t = np.asarray(mfe_residuals.t, dtype=float)
    phase_end_times = phase_boundary_times(phase_protocol.family_phases)
    if t.size == 0:
        return

    r1 = np.asarray(residuals[0], dtype=complex)
    r2 = np.asarray(residuals[1], dtype=complex)
    residual_l1 = (
        np.abs(np.real(r1))
        + np.abs(np.imag(r1))
        + np.abs(np.real(r2))
        + np.abs(np.imag(r2))
    )
    print("Residual L1 norm")
    phase_start = 0.0
    for phase_idx, phase_end in enumerate(phase_end_times, start=1):
        snapshot_idx = int(np.argmin(np.abs(t - phase_end)))
        phase_mask = (t >= phase_start) & (t <= phase_end)
        phase_min = float(np.nanmin(residual_l1[phase_mask]))
        print(f"Phase {phase_idx} - end: {residual_l1[snapshot_idx]:.3e}, min: {phase_min:.3e}")
        phase_start = phase_end


def plot_mfe_residuals(
    mfe_residuals: Any,
    *,
    show_components: bool = True,
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str = "--",
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phase_protocol: PhaseProtocol | None = None,
    print_phase_end_summary: bool = True,
    symlog: bool = True,
):
    """
    Plot stored two-group MFE residuals and their L2 norm.
    """
    if mfe_residuals is None:
        raise ValueError("plot_mfe_residuals requires moments.mfe_residuals.")
    residuals = mfe_residuals.residuals_groups
    if len(residuals) != 2:
        raise ValueError("plot_mfe_residuals currently requires exactly two residual groups.")
    line_style = validated_linestyle(linestyle)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

    fig, axes = get_axes(
        axes,
        n_axes=1,
        create_figure=lambda: plt.subplots(1, 1, figsize=(9, 4), constrained_layout=True),
        error_message="axes must contain exactly one axis for the residual plot.",
    )

    t = np.asarray(mfe_residuals.t, dtype=float)
    r1 = np.asarray(residuals[0], dtype=complex)
    r2 = np.asarray(residuals[1], dtype=complex)
    residual_l2 = np.sqrt(np.abs(r1) ** 2 + np.abs(r2) ** 2)
    residual_specs = (
        (np.real(r1), r"$\mathrm{Re}\,R_1$"),
        (np.imag(r1), r"$\mathrm{Im}\,R_1$"),
        (np.real(r2), r"$\mathrm{Re}\,R_2$"),
        (np.imag(r2), r"$\mathrm{Im}\,R_2$"),
    )

    if show_components:
        for curve_index, (values, residual_label) in enumerate(residual_specs):
            axes[0].plot(
                t,
                values,
                linewidth=1.8,
                color=palette_curve_color(palette, curve_index),
                linestyle=line_style,
                label=curve_label(residual_label, label=label),
            )
    axes[0].plot(
        t,
        residual_l2,
        linewidth=1.8,
        color=palette_curve_color(palette, len(residual_specs) if show_components else 0),
        linestyle=line_style,
        label=curve_label("L2 norm", label=label) if show_components else label,
    )
    if show_components:
        axes[0].axhline(0.0, linestyle=":", color="black", alpha=0.7)
        axes[0].set_ylabel("Residual")
        _set_symmetric_ylim_from_lines(axes[0])
    else:
        axes[0].set_ylabel("L2 norm")
    _set_residual_scale(
        axes[0],
        show_components=show_components,
        symlog=symlog,
        residual_l2=residual_l2,
    )
    style_axis(axes[0])
    if show_components or label is not None:
        axes[0].legend()
    format_time_axis(axes[0])

    finish_time_plot(
        fig,
        axes,
        phase_protocol=phase_protocol,
        title="MFE steady-state residuals" if show_components else "MFE residual L2 norm",
        output_path=output_path,
    )

    if print_phase_end_summary:
        _print_phase_end_residuals(mfe_residuals, phase_protocol)

    return fig, axes


__all__ = [
    "plot_mfe_residuals",
]
