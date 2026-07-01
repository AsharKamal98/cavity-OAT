from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting_utils import (
    curve_label,
    finish_time_plot,
    format_time_axis,
    full_curve_color,
    get_axes,
    style_axis,
)


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


def _print_phase_end_residuals(mfe_residuals: Any, phases) -> None:
    if phases is None:
        return
    residuals = mfe_residuals.residuals_groups
    if len(residuals) != 2:
        return

    t = np.asarray(mfe_residuals.t, dtype=float)
    phase_end_times = np.cumsum([phase.duration for phase in phases], dtype=float)
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
    colour_index: int = 0,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
    print_phase_end_summary: bool = True,
):
    """
    Plot stored two-group MFE residuals and their L2 norm.
    """
    if mfe_residuals is None:
        raise ValueError("plot_mfe_residuals requires moments.mfe_residuals.")
    residuals = mfe_residuals.residuals_groups
    if len(residuals) != 2:
        raise ValueError("plot_mfe_residuals currently requires exactly two residual groups.")
    full_color = full_curve_color(colour_index)

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

    residual_specs = [
        (np.real(r1), r"$\mathrm{Re}\,R_1$", "#0072B2"),
        (np.imag(r1), r"$\mathrm{Im}\,R_1$", "#56B4E9"),
        (np.real(r2), r"$\mathrm{Re}\,R_2$", "#D55E00"),
        (np.imag(r2), r"$\mathrm{Im}\,R_2$", "#E69F00"),
    ]
    for values, residual_label, color in residual_specs:
        axes[0].plot(
            t,
            values,
            linewidth=1.8,
            color=color,
            linestyle="--",
            label=curve_label(residual_label, label=label),
        )
    axes[0].plot(
        t,
        residual_l2,
        linewidth=1.8,
        color=full_color,
        linestyle="-",
        label=curve_label("L2 norm", label=label),
    )
    axes[0].axhline(0.0, linestyle=":", color="black", alpha=0.7)
    axes[0].set_ylabel("Residual")
    style_axis(axes[0])
    axes[0].legend()
    format_time_axis(axes[0])
    _set_symmetric_ylim_from_lines(axes[0])

    finish_time_plot(
        fig,
        axes,
        phases=phases,
        title="MFE steady-state residuals",
        output_path=output_path,
        title_y=1.08,
    )

    if print_phase_end_summary:
        _print_phase_end_residuals(mfe_residuals, phases)

    return fig, axes


__all__ = [
    "plot_mfe_residuals",
]
