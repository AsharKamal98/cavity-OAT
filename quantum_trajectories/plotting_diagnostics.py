from __future__ import annotations

from pathlib import Path
from typing import Any, Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from parser.quantum_trajectories import TrajectoryEnsemble, TrajectoryResult
from common.plotting_utils import (
    curve_label,
    finish_time_plot,
    format_time_axis,
    get_axes,
    sector_curve_color,
    style_axis,
)


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
    style_index: int = 0,
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
    line_styles = ("-", "--")
    line_style = line_styles[style_index % len(line_styles)]

    for sector_index, (sector, values) in enumerate(zip(sectors, probabilities)):
        axes[0].plot(
            t,
            np.asarray(values, dtype=float),
            linewidth=1.8,
            color=sector_curve_color(sector_index),
            linestyle=line_style,
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
    "plot_sector_probabilities",
]
