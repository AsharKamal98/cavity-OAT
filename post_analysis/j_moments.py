from __future__ import annotations

import numpy as np

from common.utils.phases import integration_phase_indices_at_times
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.j_moments import JMomentSeries
from parser.moments import SimulationMetadata


def _group_series(values: np.ndarray) -> tuple[np.ndarray, ...]:
    """Split a ``(time, group)`` array into independent group series."""
    return tuple(values[:, group_index].copy() for group_index in range(values.shape[1]))


def compute_j_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *,
    metadata: SimulationMetadata,
    tol: float = 1e-12,
) -> JMomentSeries:
    """Construct averaged, group-resolved J moments from mean SU(3) moments.

    The input must already have been averaged over the desired complete
    trajectories and retain the resulting singleton trajectory axis. Generator
    labels follow the one-based Gell-Mann convention with
    ``g_a = lambda_a / 2``; the stored array uses the corresponding zero-based
    final axis.

    Only first generator moments are required for this conversion. Symmetric
    second moments remain available to later S-moment, covariance, and
    squeezing analyses.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")
    if mean_generator_moments.num_trajectories != 1:
        raise ValueError(
            "compute_j_moments requires already-averaged generator moments "
            "with a singleton trajectory axis."
        )

    first_moments = mean_generator_moments.first_moments[0]
    group_count = first_moments.shape[1]
    if group_count != len(metadata.Ni):
        raise ValueError(
            "Generator moments must contain one group-axis entry for each "
            "group in metadata.Ni."
        )

    group_sizes = np.asarray(metadata.Ni, dtype=float)[np.newaxis, :]
    sqrt_three = np.sqrt(3.0)

    # The MCWF sector operators define J_y=(J_+-J_-)/(2i), which is -G_7
    # in the stored (|u>, |d>, |e>) Gell-Mann basis.
    x = first_moments[:, :, 5]
    y = -first_moments[:, :, 6]
    z = 0.5 * first_moments[:, :, 2] - 0.5 * sqrt_three * first_moments[:, :, 7]
    N_e = group_sizes / 3.0 - (2.0 / sqrt_three) * first_moments[:, :, 7]
    N_j = (
        2.0 * group_sizes / 3.0
        - first_moments[:, :, 2]
        - first_moments[:, :, 7] / sqrt_three
    )

    j_moments = JMomentSeries(
        t=mean_generator_moments.t.copy(),
        integration_phase_index=integration_phase_indices_at_times(
            mean_generator_moments.t,
            metadata.phase_protocol,
        ),
        x=_group_series(x),
        y=_group_series(y),
        z=_group_series(z),
        N_e=_group_series(N_e),
        N_j=_group_series(N_j),
    )
    JMomentSeries.attach_norm_spin_components_from_spin_components(
        j_moments,
        tol=tol,
    )
    JMomentSeries.attach_angles_from_norm_spin_components(
        j_moments,
        tol=tol,
    )
    JMomentSeries.attach_aggregated_from_group_fields(
        j_moments,
        tol=tol,
    )
    return j_moments


__all__ = ["compute_j_moments"]
