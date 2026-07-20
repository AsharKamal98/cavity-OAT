from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from common.utils.moments import as_series_tuple
from parser.bloch_vector import BlochVectorSeries
from parser.common import Array
from parser.j_modes import JModeSeries


def _bloch_vector_series(vector: Array) -> BlochVectorSeries:
    return BlochVectorSeries(
        x=vector[0],
        y=vector[1],
        z=vector[2],
    )


def compute_j_modes(
    t: Array,
    x_components: Sequence[Array] | Array,
    y_components: Sequence[Array] | Array,
    z_components: Sequence[Array] | Array,
    *,
    populations: Sequence[float | Array],
    omega_groups: Sequence[float],
) -> JModeSeries:
    """Compute common/contrast and drive-bright/dark modes for one or two groups."""
    x_groups = as_series_tuple(x_components)
    y_groups = as_series_tuple(y_components)
    z_groups = as_series_tuple(z_components)
    group_count = len(x_groups)
    if (
        group_count not in (1, 2)
        or len(y_groups) != group_count
        or len(z_groups) != group_count
        or len(populations) != group_count
        or len(omega_groups) != group_count
    ):
        raise ValueError("J modes require matching one- or two-group inputs.")

    t = np.asarray(t, dtype=float)
    vector_1 = np.stack((x_groups[0], y_groups[0], z_groups[0]))
    if vector_1.shape != (3, t.size):
        raise ValueError("Each input vector must contain x/y/z arrays matching t.")

    population_1 = np.asarray(populations[0], dtype=float)
    if np.any(population_1 <= 0.0):
        raise ValueError("Group populations must be non-negative with a positive sum.")

    if group_count == 1:
        drive_weight = population_1 * float(omega_groups[0])
        if np.any(drive_weight == 0.0):
            raise ValueError("Bright/dark modes require at least one nonzero drive weight.")
        zero = np.zeros_like(vector_1)
        return JModeSeries(
            t=t,
            common=_bloch_vector_series(vector_1),
            contrast=_bloch_vector_series(zero),
            bright=_bloch_vector_series(np.sign(drive_weight) * vector_1),
            dark=_bloch_vector_series(zero),
        )

    vector_2 = np.stack((x_groups[1], y_groups[1], z_groups[1]))
    if vector_2.shape != (3, t.size):
        raise ValueError("Each input vector must contain x/y/z arrays matching t.")

    population_2 = np.asarray(populations[1], dtype=float)
    total_population = population_1 + population_2
    if np.any(population_2 < 0.0) or np.any(total_population <= 0.0):
        raise ValueError("Group populations must be non-negative with a positive sum.")

    common = (population_1 * vector_1 + population_2 * vector_2) / total_population
    contrast = 0.5 * (vector_1 - vector_2)

    drive_weight_1 = population_1 * float(omega_groups[0])
    drive_weight_2 = population_2 * float(omega_groups[1])
    drive_scale = np.abs(drive_weight_1) + np.abs(drive_weight_2)
    if np.any(drive_scale == 0.0):
        raise ValueError("Bright/dark modes require at least one nonzero drive weight.")
    bright = (drive_weight_1 * vector_1 + drive_weight_2 * vector_2) / drive_scale
    dark = (drive_weight_2 * vector_1 - drive_weight_1 * vector_2) / drive_scale

    return JModeSeries(
        t=t,
        common=_bloch_vector_series(common),
        contrast=_bloch_vector_series(contrast),
        bright=_bloch_vector_series(bright),
        dark=_bloch_vector_series(dark),
    )


__all__ = ["compute_j_modes"]
