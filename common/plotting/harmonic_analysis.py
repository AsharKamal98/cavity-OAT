from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np

from common.plotting.utils import (
    colour_palette,
    get_axes,
    palette_curve_color,
    save_figure,
    set_bottom_figure_legend,
    shade_axes_for_family_phase,
    style_axis,
    validated_linestyle,
    validated_marker,
)
from common.utils.moments import as_series_tuple
from parser.common import Array


def _validated_sweep_arrays(
    parameter_values: Sequence[float] | Array,
    fundamental_frequencies: Sequence[Array] | Array,
    total_harmonic_distortions: Sequence[Array] | Array,
    rms_amplitudes: Sequence[Array] | Array,
    offsets: Sequence[Array] | Array,
    labels: Sequence[str],
) -> tuple[
    Array,
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[str, ...],
]:
    """Return matching finite one-dimensional sweep curves."""
    parameter_values = np.asarray(parameter_values, dtype=float)
    frequency_curves = as_series_tuple(fundamental_frequencies)
    distortion_curves = as_series_tuple(total_harmonic_distortions)
    rms_amplitude_curves = as_series_tuple(rms_amplitudes)
    offset_curves = as_series_tuple(offsets)
    labels = tuple(labels)

    if parameter_values.ndim != 1:
        raise ValueError("Parameter values must be one-dimensional.")
    if parameter_values.size == 0:
        raise ValueError("Harmonic sweep inputs must not be empty.")
    curve_counts = tuple(
        len(curves)
        for curves in (
            frequency_curves,
            distortion_curves,
            rms_amplitude_curves,
            offset_curves,
        )
    )
    if any(count != len(labels) for count in curve_counts):
        raise ValueError("Each harmonic sweep curve requires one label.")

    curves = (
        frequency_curves
        + distortion_curves
        + rms_amplitude_curves
        + offset_curves
    )
    if any(values.ndim != 1 for values in curves):
        raise ValueError("Harmonic sweep curves must be one-dimensional.")
    if any(values.size != parameter_values.size for values in curves):
        raise ValueError("Harmonic sweep inputs must have matching lengths.")
    finite_curves = frequency_curves + rms_amplitude_curves + offset_curves
    if not np.all(np.isfinite(parameter_values)) or any(
        not np.all(np.isfinite(values)) for values in finite_curves
    ):
        raise ValueError(
            "Parameters, frequencies, RMS amplitudes, and offsets must be finite."
        )
    if any(np.any(np.isinf(values)) for values in distortion_curves):
        raise ValueError("THD curves must not contain infinite values.")
    return (
        parameter_values,
        frequency_curves,
        distortion_curves,
        rms_amplitude_curves,
        offset_curves,
        labels,
    )


def plot_harmonic_sweep(
    parameter_values: Sequence[float] | Array,
    fundamental_frequencies: Sequence[Array] | Array,
    total_harmonic_distortions: Sequence[Array] | Array,
    rms_amplitudes: Sequence[Array] | Array,
    offsets: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    parameter_label: str,
    family_phase_index: int | None = None,
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str | None = "-",
    marker: str | None = None,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    title: str = "Harmonic analysis",
):
    """Plot four harmonic-analysis metrics against a varied parameter."""
    (
        parameter_values,
        frequency_curves,
        distortion_curves,
        rms_amplitude_curves,
        offset_curves,
        labels,
    ) = _validated_sweep_arrays(
        parameter_values,
        fundamental_frequencies,
        total_harmonic_distortions,
        rms_amplitudes,
        offsets,
        labels,
    )
    line_style = validated_linestyle(linestyle)
    marker_style = validated_marker(marker)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

    creating_figure = axes is None
    fig, axes = get_axes(
        axes,
        n_axes=4,
        create_figure=lambda: plt.subplots(
            2,
            2,
            figsize=(10, 7),
            sharex=True,
            constrained_layout=True,
        ),
        error_message="axes must contain exactly four axes for the harmonic sweep.",
    )

    curve_sets = zip(
        frequency_curves,
        distortion_curves,
        rms_amplitude_curves,
        offset_curves,
        labels,
    )
    for curve_index, curve_set in enumerate(curve_sets):
        frequencies, distortions, rms_values, offset_values, label = curve_set
        color = palette_curve_color(palette, curve_index)
        curve_style = {
            "linewidth": 1.8,
            "color": color,
            "linestyle": line_style,
            "marker": marker_style,
            "label": label,
        }
        axes[0].plot(parameter_values, frequencies, **curve_style)
        axes[1].plot(parameter_values, rms_values, **curve_style)
        axes[2].plot(parameter_values, offset_values, **curve_style)
        axes[3].plot(parameter_values, distortions, **curve_style)

    if creating_figure:
        panel_specs = (
            (axes[0], r"$f_0$", "Fundamental frequency"),
            (axes[1], r"$A_{\mathrm{RMS}}$", "RMS oscillation amplitude"),
            (axes[2], r"$c$", "Offset"),
            (axes[3], "THD", "Total harmonic distortion"),
        )
        for ax, ylabel, panel_title in panel_specs:
            ax.set_ylabel(ylabel)
            ax.set_title(panel_title, fontsize=11)
            style_axis(ax)
            ax.ticklabel_format(
                axis="x",
                style="sci",
                scilimits=(0, 0),
                useOffset=False,
            )
        for ax in axes[2:]:
            ax.set_xlabel(parameter_label)
        if family_phase_index is not None:
            shade_axes_for_family_phase(axes, family_phase_index)
        fig.suptitle(title, y=1.06, fontsize=14)

    set_bottom_figure_legend(
        fig,
        axes[0],
        max_columns=3,
        row_major=True,
    )
    save_figure(fig, output_path)
    return fig, axes


__all__ = ["plot_harmonic_sweep"]
