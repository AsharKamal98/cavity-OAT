from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.j_moments import JMomentSeries
from parser.moments import SimulationMetadata
from parser.squeezing import SqueezingResult, SqueezingSeries


# In the dressed (|c>, |j>, |s>) basis these normalized Gell--Mann
# generators are the four fluctuation operators used for squeezing.
_FLUCTUATION_GENERATOR_INDICES = np.array((0, 1, 3, 4), dtype=np.intp)


def _unpack_second_moments(
    packed_moments: Array,
    *,
    generator_count: int,
) -> Array:
    """Unpack symmetric moments into ``(time, generator, generator)``."""
    packed_moments = np.asarray(packed_moments, dtype=float)
    upper = np.triu_indices(generator_count)
    expected_shape = (packed_moments.shape[0], upper[0].size)
    if packed_moments.shape != expected_shape:
        raise ValueError(
            f"packed_moments must have shape {expected_shape}."
        )

    full_moments = np.zeros(
        (packed_moments.shape[0], generator_count, generator_count),
        dtype=float,
    )
    full_moments[:, upper[0], upper[1]] = packed_moments
    full_moments[:, upper[1], upper[0]] = packed_moments
    return full_moments


def _validated_populations(
    j_moments: JMomentSeries,
    *,
    num_times: int,
    num_groups: int,
    tol: float,
) -> tuple[Array, Array, Array]:
    """Return group populations and their active excited fractions."""
    if j_moments.N_e is None or j_moments.N_j is None:
        raise ValueError("Squeezing requires group-resolved N_e and N_j.")

    N_e = np.asarray(j_moments.N_e, dtype=float)
    N_j = np.asarray(j_moments.N_j, dtype=float)
    expected_shape = (num_groups, num_times)
    if N_e.shape != expected_shape or N_j.shape != expected_shape:
        raise ValueError(
            "J-moment N_e and N_j must have shape (num_groups, len(t))."
        )
    if not np.all(np.isfinite(N_e)) or not np.all(np.isfinite(N_j)):
        raise ValueError("J-moment populations must be finite.")

    N_e = N_e.T.copy()
    N_j = N_j.T.copy()
    excited_fraction = np.full_like(N_e, np.nan, dtype=float)
    valid = N_j > tol
    excited_fraction[valid] = N_e[valid] / N_j[valid]
    return N_e, N_j, excited_fraction


def _compute_squeezing_series(
    *,
    t: Array,
    N_atoms: int,
    covariance: Array,
    N_c: Array,
    excited_fraction_active: Array,
    tol: float,
) -> SqueezingSeries:
    """Compute covariance eigenvalues and normalized squeezing."""
    covariance = 0.5 * (
        np.asarray(covariance, dtype=float)
        + np.swapaxes(covariance, -1, -2)
    )
    covariance_eigvals = np.linalg.eigvalsh(covariance)
    lambda_min = covariance_eigvals[:, 0]

    N_c = np.asarray(N_c, dtype=float).copy()
    small_negative = (N_c < 0.0) & (np.abs(N_c) <= 100.0 * tol)
    N_c[small_negative] = 0.0

    denominator = (N_c / 2.0) ** 2
    xi2 = np.full(t.shape, np.nan, dtype=float)
    valid_denominator = denominator > tol
    xi2[valid_denominator] = (
        N_atoms
        * lambda_min[valid_denominator]
        / denominator[valid_denominator]
    )

    xi2_db = np.full_like(xi2, np.nan, dtype=float)
    valid_xi2 = np.isfinite(xi2) & (xi2 > 0.0)
    xi2_db[valid_xi2] = 10.0 * np.log10(xi2[valid_xi2])

    return SqueezingSeries(
        t=np.asarray(t, dtype=float).copy(),
        N_atoms=int(N_atoms),
        covariance=covariance,
        covariance_eigvals=covariance_eigvals,
        lambda_min=lambda_min,
        xi2=xi2,
        xi2_db=xi2_db,
        N_c=N_c,
        excited_fraction_active=excited_fraction_active,
    )


def compute_squeezing(
    mean_dressed_generator_moments: TrajectoryResolvedMomentSeries,
    j_moments: JMomentSeries,
    *,
    metadata: SimulationMetadata,
    include_full_system: bool = False,
    tol: float = 1e-12,
) -> SqueezingResult:
    """Compute group-resolved and optional full-system generalized squeezing.

    The input generator moments must already be averaged over the selected
    complete trajectories, rotated into the ``(|c>, |j>, |s>)`` frame, and
    retain a singleton trajectory axis. The four dressed generators
    ``(T_1, T_2, T_4, T_5)`` per group define the squeezing covariance.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")
    if not isinstance(include_full_system, (bool, np.bool_)):
        raise TypeError("include_full_system must be a boolean.")
    if mean_dressed_generator_moments.num_trajectories != 1:
        raise ValueError(
            "compute_squeezing requires averaged dressed generator moments "
            "with a singleton trajectory axis."
        )

    num_times = mean_dressed_generator_moments.t.size
    num_groups = mean_dressed_generator_moments.first_moments.shape[2]
    if num_groups != len(metadata.Ni):
        raise ValueError(
            "Generator moments must contain one group for each metadata.Ni entry."
        )
    if np.asarray(j_moments.t).shape != (num_times,) or not np.allclose(
        j_moments.t,
        mean_dressed_generator_moments.t,
        atol=tol,
        rtol=0.0,
    ):
        raise ValueError(
            "J moments and dressed generator moments must share a time grid."
        )

    N_e, N_j, excited_fraction = _validated_populations(
        j_moments,
        num_times=num_times,
        num_groups=num_groups,
        tol=tol,
    )
    generator_count = 8 * num_groups
    dressed_first = mean_dressed_generator_moments.first_moments[0]
    dressed_second = _unpack_second_moments(
        mean_dressed_generator_moments.second_moments[0],
        generator_count=generator_count,
    )
    selected_indices = np.concatenate(
        [
            8 * group_index + _FLUCTUATION_GENERATOR_INDICES
            for group_index in range(num_groups)
        ]
    )
    selected_first = dressed_first.reshape(num_times, generator_count)[
        :, selected_indices
    ]
    selected_second = dressed_second[
        :, selected_indices[:, None], selected_indices[None, :]
    ]
    covariance = selected_second - (
        selected_first[:, :, None] * selected_first[:, None, :]
    )
    covariance = 0.5 * (covariance + np.swapaxes(covariance, -1, -2))

    N_groups = np.asarray(metadata.Ni, dtype=int)
    N_c_groups = (
        N_groups[np.newaxis, :] / 3.0
        + dressed_first[:, :, 2]
        + dressed_first[:, :, 7] / np.sqrt(3.0)
    )
    group_resolved = tuple(
        _compute_squeezing_series(
            t=mean_dressed_generator_moments.t,
            N_atoms=int(N_groups[group_index]),
            covariance=covariance[
                :,
                4 * group_index : 4 * (group_index + 1),
                4 * group_index : 4 * (group_index + 1),
            ],
            N_c=N_c_groups[:, group_index],
            excited_fraction_active=excited_fraction[:, group_index],
            tol=tol,
        )
        for group_index in range(num_groups)
    )

    full_system = None
    if include_full_system:
        covariance_blocks = covariance.reshape(
            num_times,
            num_groups,
            4,
            num_groups,
            4,
        )
        full_covariance = np.sum(covariance_blocks, axis=(1, 3))
        total_N_j = np.sum(N_j, axis=1)
        full_excited_fraction = np.full(num_times, np.nan, dtype=float)
        valid_active_population = total_N_j > tol
        full_excited_fraction[valid_active_population] = (
            np.sum(N_e, axis=1)[valid_active_population]
            / total_N_j[valid_active_population]
        )
        full_system = _compute_squeezing_series(
            t=mean_dressed_generator_moments.t,
            N_atoms=int(np.sum(N_groups)),
            covariance=full_covariance,
            N_c=np.sum(N_c_groups, axis=1),
            excited_fraction_active=full_excited_fraction,
            tol=tol,
        )

    return SqueezingResult(
        group_resolved=group_resolved,
        full_system=full_system,
    )


__all__ = ["compute_squeezing"]
