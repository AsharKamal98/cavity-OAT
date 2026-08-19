from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np

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


def _validated_squeezing_curves(
    t: Array,
    xi2: Sequence[Array] | Array,
    lambda_min: Sequence[Array] | Array,
    N_c: Sequence[Array] | Array,
    excited_fraction_active: Sequence[Array] | Array,
    labels: Sequence[str],
) -> tuple[
    Array,
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[Array, ...],
    tuple[str, ...],
]:
    """Return matching one-dimensional squeezing curves."""
    t = np.asarray(t, dtype=float)
    curve_sets = tuple(
        as_series_tuple(values)
        for values in (
            xi2,
            lambda_min,
            N_c,
            excited_fraction_active,
        )
    )
    labels = tuple(labels)

    if t.ndim != 1 or t.size == 0:
        raise ValueError("t must be a non-empty one-dimensional array.")
    if not np.all(np.isfinite(t)):
        raise ValueError("t must contain only finite values.")
    if not labels:
        raise ValueError("plot_squeezing requires at least one curve label.")
    if any(len(curves) != len(labels) for curves in curve_sets):
        raise ValueError("Each squeezing metric requires one curve per label.")

    for curve in sum(curve_sets, ()):
        if curve.shape != t.shape:
            raise ValueError("All squeezing curves must match the t grid.")
        if np.any(np.isinf(curve)):
            raise ValueError("Squeezing curves must not contain infinities.")

    return (t, *curve_sets, labels)


def _validated_covariance_eigenvalues(
    covariance_eigvals: Sequence[Array] | Array | None,
    *,
    t: Array,
    num_curves: int,
    plot_all_eigenvalues: bool,
) -> tuple[Array, ...] | None:
    """Return matching four-eigenvalue series when requested."""
    if covariance_eigvals is None:
        if plot_all_eigenvalues:
            raise ValueError(
                "covariance_eigvals is required when plot_all_eigenvalues=True."
            )
        return None

    values = np.asarray(covariance_eigvals, dtype=float)
    if values.ndim == 2:
        curves = (values,)
    elif values.ndim == 3:
        curves = tuple(values)
    else:
        raise ValueError(
            "covariance_eigvals must have shape (len(t), 4) or "
            "(num_curves, len(t), 4)."
        )

    if len(curves) != num_curves:
        raise ValueError("covariance_eigvals requires one array per label.")
    for curve in curves:
        if curve.shape != (t.size, 4):
            raise ValueError(
                "Each covariance_eigvals array must have shape (len(t), 4)."
            )
        if np.any(np.isinf(curve)):
            raise ValueError("covariance_eigvals must not contain infinities.")
    return curves


def _xi2_to_db(xi2: Array) -> Array:
    """Convert positive squeezing values to dB and mask invalid values."""
    xi2 = np.asarray(xi2, dtype=float)
    xi2_db = np.full_like(xi2, np.nan, dtype=float)
    valid = np.isfinite(xi2) & (xi2 > 0.0)
    xi2_db[valid] = 10.0 * np.log10(xi2[valid])
    return xi2_db


def plot_squeezing(
    t: Array,
    xi2: Sequence[Array] | Array,
    lambda_min: Sequence[Array] | Array,
    N_c: Sequence[Array] | Array,
    excited_fraction_active: Sequence[Array] | Array,
    *,
    labels: Sequence[str],
    covariance_eigvals: Sequence[Array] | Array | None = None,
    plot_all_eigenvalues: bool = False,
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
    linestyle: str | None = "-",
    marker: str | None = None,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    phase_protocol: PhaseProtocol | None = None,
    title: str = "Generalized squeezing",
):
    """Plot one or more generalized-squeezing results on a 2x2 grid."""
    (
        t,
        xi2_curves,
        lambda_min_curves,
        N_c_curves,
        excited_fraction_curves,
        labels,
    ) = _validated_squeezing_curves(
        t,
        xi2,
        lambda_min,
        N_c,
        excited_fraction_active,
        labels,
    )
    eigenvalue_curves = _validated_covariance_eigenvalues(
        covariance_eigvals,
        t=t,
        num_curves=len(labels),
        plot_all_eigenvalues=plot_all_eigenvalues,
    )
    line_style = validated_linestyle(linestyle)
    marker_style = validated_marker(marker)
    palette = colour_palette(
        colour_family_index=colour_family_index,
        shade_index=shade_index,
    )

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
        error_message="axes must contain exactly four axes for the squeezing plot.",
    )

    curve_sets = zip(
        xi2_curves,
        lambda_min_curves,
        N_c_curves,
        excited_fraction_curves,
        labels,
    )
    for curve_index, curve_set in enumerate(curve_sets):
        xi2_values, lambda_values, N_c_values, excited_values, label = curve_set
        curve_style = {
            "linewidth": 1.8,
            "color": palette_curve_color(palette, curve_index),
            "linestyle": line_style,
            "marker": marker_style,
            "label": label,
        }
        axes[0].plot(t, _xi2_to_db(xi2_values), **curve_style)
        if plot_all_eigenvalues:
            assert eigenvalue_curves is not None
            eigenvalues = eigenvalue_curves[curve_index]
            for eigenvalue_index in range(eigenvalues.shape[1]):
                axes[1].plot(
                    t,
                    np.where(
                        eigenvalues[:, eigenvalue_index] > 0.0,
                        eigenvalues[:, eigenvalue_index],
                        np.nan,
                    ),
                    **{
                        **curve_style,
                        "alpha": 1.0 - 0.15 * eigenvalue_index,
                        "label": "_nolegend_",
                    },
                )
        else:
            axes[1].plot(
                t,
                np.where(lambda_values > 0.0, lambda_values, np.nan),
                **curve_style,
            )
        axes[2].plot(t, N_c_values, **curve_style)
        axes[3].plot(t, excited_values, **curve_style)

    panel_specs = (
        (axes[0], r"$\xi^2$ [dB]", "Generalized squeezing"),
        (
            axes[1],
            (
                r"$\lambda_i(C)$"
                if plot_all_eigenvalues
                else r"$\lambda_{\min}(C)$"
            ),
            (
                "Covariance eigenvalues"
                if plot_all_eigenvalues
                else "Smallest covariance eigenvalue"
            ),
        ),
        (axes[2], r"$\langle N_c\rangle$", "Mean population"),
        (
            axes[3],
            r"$\langle N_e\rangle/\langle N_J\rangle$",
            "Active excited fraction",
        ),
    )
    for ax, ylabel, panel_title in panel_specs:
        ax.set_ylabel(ylabel)
        ax.set_title(panel_title, fontsize=11)
        style_axis(ax)
        format_time_axis(ax)
    axes[1].set_yscale("log")

    set_bottom_figure_legend(fig, axes[0])
    finish_time_plot(
        fig,
        axes,
        phase_protocol=phase_protocol,
        title=title,
        output_path=output_path,
    )
    return fig, axes


__all__ = ["plot_squeezing"]
