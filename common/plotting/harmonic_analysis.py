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
from parser.harmonic_analysis import DampedOscillationSeries


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
    """Validate matching curves, allowing NaNs to mark unresolved metrics."""
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
    if not np.all(np.isfinite(parameter_values)):
        raise ValueError("Parameter values must be finite.")
    if any(np.any(np.isinf(values)) for values in curves):
        raise ValueError("Harmonic sweep curves must not contain infinite values.")
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
    quality_masks: Sequence[Array] | Array | None = None,
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
    masks = None
    if quality_masks is not None:
        mask_values = np.asarray(quality_masks)
        if mask_values.ndim not in (1, 2):
            raise ValueError("Quality masks must be one-dimensional curves.")
        masks = (mask_values,) if mask_values.ndim == 1 else tuple(mask_values)
        if len(masks) != len(labels) or any(
            mask.shape != parameter_values.shape for mask in masks
        ):
            raise ValueError("Each quality mask must match its harmonic curve.")
        if any(mask.dtype != np.dtype(bool) for mask in masks):
            raise ValueError("Quality masks must contain booleans.")
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
        for axis, values in (
            (axes[0], frequencies), (axes[1], rms_values), (axes[3], distortions)
        ):
            if masks is None:
                axis.plot(parameter_values, values, **curve_style)
            else:
                accepted = masks[curve_index]
                axis.plot(
                    parameter_values, np.where(accepted, values, np.nan), **curve_style
                )
                flagged = ~accepted & np.isfinite(values)
                axis.scatter(
                    parameter_values[flagged], values[flagged],
                    color=color, marker="x", s=28, alpha=0.65,
                )
        axes[2].plot(parameter_values, offset_values, **curve_style)

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


def plot_damped_oscillations(
    t: Array,
    series: Sequence[Array] | Array,
    analysis: DampedOscillationSeries,
    *,
    labels: Sequence[str],
    sem: Sequence[Array] | Array | None = None,
    reference_t: Array | None = None,
    reference_series: Sequence[Array] | Array | None = None,
    reference_label: str = "Reference",
    time_scale: float = 1.0,
    time_label: str = r"Time since fit start",
    output_path: Optional[Union[str, Path]] = None,
):
    """Plot early traces, damped fits, and an optional dashed reference per panel."""
    t = np.asarray(t, dtype=float)
    curves = as_series_tuple(series)
    labels = tuple(labels)
    if len(curves) != len(labels) or len(curves) != len(analysis.signals):
        raise ValueError("Each damped fit and input curve requires one label.")
    if t.ndim != 1 or t.size == 0 or any(curve.shape != t.shape for curve in curves):
        raise ValueError("Damped-fit curves must match the supplied time grid.")
    if not np.isfinite(time_scale) or time_scale <= 0.0:
        raise ValueError("time_scale must be finite and positive.")
    sem_curves = (None,) * len(curves) if sem is None else as_series_tuple(sem)
    if len(sem_curves) != len(curves) or any(
        values is not None and (
            values.shape != t.shape
            or np.any(values < 0.0)
            or np.any(np.isinf(values))
        )
        for values in sem_curves
    ):
        raise ValueError(
            "Each nonnegative SEM curve must match t and contain no infinities."
        )
    reference_curves = (None,) * len(curves)
    if reference_t is not None or reference_series is not None:
        if reference_t is None or reference_series is None:
            raise ValueError("Provide both reference_t and reference_series.")
        reference_t = np.asarray(reference_t, dtype=float)
        reference_curves = as_series_tuple(reference_series)
        if (
            reference_t.ndim != 1
            or not np.all(np.isfinite(reference_t))
            or np.any(np.diff(reference_t) <= 0.0)
            or len(reference_curves) != len(curves)
            or any(values.shape != reference_t.shape for values in reference_curves)
        ):
            raise ValueError("Each reference curve must match an increasing reference grid.")
        reference_mask = (reference_t >= t[0]) & (reference_t <= t[-1])
        if not np.any(reference_mask):
            raise ValueError("The reference must overlap the fit interval.")
        reference_plot_t = time_scale * (reference_t[reference_mask] - t[0])
    fig, axes = plt.subplots(
        len(curves), 1, figsize=(10, 2.6 * len(curves)),
        sharex=True, squeeze=False, constrained_layout=True,
    )
    axes = axes[:, 0]
    plot_t = time_scale * (t - t[0])
    palette = colour_palette(shade_index=3)
    mean_color, fit_color, background_color = palette[0], palette[1], palette[4]
    for ax, curve, fit, label, errors, reference in zip(
        axes, curves, analysis.signals, labels, sem_curves, reference_curves
    ):
        ax.plot(plot_t, curve, color=mean_color, linewidth=1.2, label="MCWF mean")
        if errors is not None:
            ax.fill_between(
                plot_t, curve - errors, curve + errors,
                color=mean_color, alpha=0.22, label="±1 SEM",
            )
        if fit is not None:
            if not np.array_equal(fit.t, t):
                raise ValueError("Each damped fit must use the supplied time grid.")
            ax.plot(
                plot_t, fit.fitted_values, color=fit_color,
                label="Damped-sinusoid fit",
            )
            ax.plot(
                plot_t, fit.background, color=background_color,
                linestyle=":", label="Linear background",
            )
            ax.plot(
                plot_t, fit.background + fit.envelope, color=fit_color,
                linestyle="--", alpha=0.7, label="Fitted envelope",
            )
            ax.plot(
                plot_t, fit.background - fit.envelope, color=fit_color,
                linestyle="--", alpha=0.7,
            )
            ax.set_title(f"{label}; fit R²={fit.fit_r_squared:.3f}")
        else:
            ax.set_title(f"{label}; no oscillatory fit")
        if reference is not None:
            ax.plot(
                reference_plot_t, reference[reference_mask],
                color=palette[2], linestyle="--", label=reference_label,
            )
        style_axis(ax)
        ax.set_ylabel(label)
    axes[-1].set_xlabel(time_label)
    fig.suptitle("Early-time decay fit", fontsize=14)
    axes[0].legend(loc="best", fontsize=8)
    save_figure(fig, output_path)
    return fig, axes


__all__ = ["plot_damped_oscillations", "plot_harmonic_sweep"]
