from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from quantum_trajectories.plotting_utils import (
    curve_label,
    finish_time_plot,
    format_time_axis,
    full_curve_color,
    get_axes,
    sector_curve_color,
    style_axis,
)
from quantum_trajectories.parser import TrajectoryEnsemble, TrajectoryResult


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


def _print_phase_end_residuals(j_moments: Any, phases) -> None:
    if phases is None:
        return
    residuals = j_moments.mfe_residuals_groups
    if residuals is None or len(residuals) != 2:
        return

    t = np.asarray(j_moments.t, dtype=float)
    phase_end_times = np.cumsum([phase.duration for phase in phases], dtype=float)
    if t.size == 0:
        return

    r1 = np.asarray(residuals[0], dtype=complex)
    r2 = np.asarray(residuals[1], dtype=complex)
    print("Residual L1 norm")
    for phase_idx, phase_end in enumerate(phase_end_times, start=1):
        snapshot_idx = int(np.argmin(np.abs(t - phase_end)))
        residual_sum = (
            abs(np.real(r1[snapshot_idx]))
            + abs(np.imag(r1[snapshot_idx]))
            + abs(np.real(r2[snapshot_idx]))
            + abs(np.imag(r2[snapshot_idx]))
        )
        print("Phase " f"{phase_idx}: {residual_sum:.3e}")


def plot_mfe_residuals(
    j_moments: Any,
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
    residuals = j_moments.mfe_residuals_groups
    if residuals is None:
        raise ValueError("plot_mfe_residuals requires j_moments.mfe_residuals_groups.")
    if len(residuals) != 2:
        raise ValueError("plot_mfe_residuals currently requires exactly two residual groups.")
    full_color = full_curve_color(colour_index)

    fig, axes = get_axes(
        axes,
        n_axes=1,
        create_figure=lambda: plt.subplots(1, 1, figsize=(9, 4), constrained_layout=True),
        error_message="axes must contain exactly one axis for the residual plot.",
    )

    t = np.asarray(j_moments.t, dtype=float)
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
        _print_phase_end_residuals(j_moments, phases)

    return fig, axes


def _sector_label(sector: Any) -> str:
    if isinstance(sector, tuple):
        return rf"$({sector[0]}, {sector[1]})$"
    return rf"$N_J={sector}$"


def _sector_probability_data(result: TrajectoryResult | TrajectoryEnsemble):
    if isinstance(result, TrajectoryEnsemble):
        trajectories = result.trajectories
        reference = trajectories[0]
    else:
        trajectories = [result]
        reference = result

    sectors = tuple(reference.sectors)
    t = np.asarray([snap.time for snap in reference.snapshots], dtype=float)
    probabilities = []
    for traj in trajectories:
        weights = np.asarray(
            [
                [
                    float(np.vdot(snap.sector_blocks[sector], snap.sector_blocks[sector]).real)
                    if sector in snap.sector_blocks
                    else 0.0
                    for snap in traj.snapshots
                ]
                for sector in sectors
            ],
            dtype=float,
        )
        totals = np.sum(weights, axis=0, keepdims=True)
        probabilities.append(
            np.divide(weights, totals, out=np.zeros_like(weights), where=totals > 0.0)
        )

    return t, sectors, np.mean(probabilities, axis=0), reference.phases


def plot_sector_probabilities(
    result: TrajectoryResult | TrajectoryEnsemble,
    *,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    label: Optional[str] = None,
    phases=None,
):
    """
    Plot normalized represented-sector probabilities from saved sector blocks.
    """
    fig, axes = get_axes(
        axes,
        n_axes=1,
        create_figure=lambda: plt.subplots(1, 1, figsize=(9, 4), constrained_layout=True),
        error_message="axes must contain exactly one axis for the sector-probability plot.",
    )

    t, sectors, probabilities, default_phases = _sector_probability_data(result)
    if phases is None:
        phases = default_phases

    for sector_index, (sector, values) in enumerate(zip(sectors, probabilities)):
        axes[0].plot(
            t,
            np.asarray(values, dtype=float),
            linewidth=1.8,
            color=sector_curve_color(sector_index),
            linestyle="-",
            label=curve_label(_sector_label(sector), label=label),
        )

    axes[0].set_ylabel(r"$p_\alpha$")
    axes[0].set_title("Sector probabilities", fontsize=11)
    style_axis(axes[0])
    axes[0].legend()
    format_time_axis(axes[0])

    finish_time_plot(
        fig,
        axes,
        phases=phases,
        title="Represented-sector probabilities",
        output_path=output_path,
    )

    return fig, axes


__all__ = [
    "plot_mfe_residuals",
    "plot_sector_probabilities",
]
