from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries


def _validate_trajectory_indices(
    trajectory_indices: Array,
    *,
    num_trajectories: int,
) -> Array:
    """Return a non-empty, one-dimensional integer trajectory selection."""
    indices = np.asarray(trajectory_indices)
    if indices.ndim != 1:
        raise ValueError("trajectory_indices must be one-dimensional.")
    if indices.size == 0:
        raise ValueError("trajectory_indices must select at least one trajectory.")
    if not np.issubdtype(indices.dtype, np.integer):
        raise TypeError("trajectory_indices must contain integers.")

    indices = indices.astype(np.intp, copy=False)
    if np.any(indices < 0) or np.any(indices >= num_trajectories):
        raise IndexError("trajectory_indices contains an out-of-range index.")
    return indices


def compute_mean_generator_moments(
    moments: TrajectoryResolvedMomentSeries,
    *,
    trajectory_indices: Array | None = None,
) -> TrajectoryResolvedMomentSeries:
    """Average complete trajectories while retaining a singleton trajectory axis.

    Repeated ``trajectory_indices`` are allowed so callers can use the same
    operation for bootstrap resampling. First and symmetric second moments are
    averaged independently; covariance construction belongs to later analysis.
    """
    if trajectory_indices is None:
        mean_first_moments = np.mean(
            moments.first_moments,
            axis=0,
            keepdims=True,
            dtype=float,
        )
        mean_second_moments = np.mean(
            moments.second_moments,
            axis=0,
            keepdims=True,
            dtype=float,
        )
    else:
        indices = _validate_trajectory_indices(
            trajectory_indices,
            num_trajectories=moments.num_trajectories,
        )
        weights = np.bincount(
            indices,
            minlength=moments.num_trajectories,
        ).astype(float)
        weights /= indices.size
        mean_first_moments = np.tensordot(
            weights,
            moments.first_moments,
            axes=(0, 0),
        )[np.newaxis, ...]
        mean_second_moments = np.tensordot(
            weights,
            moments.second_moments,
            axes=(0, 0),
        )[np.newaxis, ...]

    return TrajectoryResolvedMomentSeries(
        t=moments.t.copy(),
        first_moments=mean_first_moments,
        second_moments=mean_second_moments,
    )


__all__ = ["compute_mean_generator_moments"]
