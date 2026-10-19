from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from common.utils.moments import as_series_tuple
from common.utils.uncertainty import propagate_covariance
from parser.bloch_vector import BlochVectorSeries
from parser.common import Array
from parser.j_modes import JModeSeries
from parser.j_moments import JMomentSeries
from parser.uncertainty import SamplingUncertainty


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
    contrast = (population_2 * vector_1 - population_1 * vector_2) / total_population

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


def _j_mode_jacobians(
    normalized_vectors: Array,
    populations: Array,
    omega_groups: Array,
) -> dict[str, Array]:
    """Differentiate each mode with respect to raw group J fields and N_j.

    Source columns follow (x, y, z, N_e, N_j) within each group. The
    derivatives include both 2 J / N_j and the population-dependent weights.
    """
    num_times, group_count, _ = normalized_vectors.shape
    identity = np.eye(3)
    if group_count == 1:
        common = np.zeros((num_times, 3, 5), dtype=float)
        common[:, :, :3] = (2.0 / populations[:, 0])[:, None, None] * identity
        common[:, :, 4] = -normalized_vectors[:, 0] / populations[:, 0, None]
        zero = np.zeros_like(common)
        return {
            "common": common,
            "contrast": zero,
            "bright": np.sign(omega_groups[0]) * common,
            "dark": zero,
        }

    # Numerators of the two vector coefficients are linear in populations.
    # Each pair contains that linear map and the denominator's gradient.
    omega_1, omega_2 = omega_groups
    coefficient_maps = {
        "common": (np.eye(2), np.ones(2)),
        "contrast": (np.array([[0.0, 1.0], [-1.0, 0.0]]), np.ones(2)),
        "bright": (np.diag(omega_groups), np.abs(omega_groups)),
        "dark": (np.array([[0.0, omega_2], [-omega_1, 0.0]]), np.abs(omega_groups)),
    }
    jacobians = {}
    for name, (numerator_map, denominator_gradient) in coefficient_maps.items():
        denominator = populations @ denominator_gradient
        coefficients = (populations @ numerator_map.T) / denominator[:, None]
        # Quotient rule: rows identify vector coefficients, columns populations.
        coefficient_derivatives = (
            numerator_map[None, :, :]
            - coefficients[:, :, None] * denominator_gradient[None, None, :]
        ) / denominator[:, None, None]

        jacobian = np.zeros((num_times, 3, 10), dtype=float)
        for group in range(2):
            coefficient = coefficients[:, group]
            vector = normalized_vectors[:, group]
            component_columns = slice(5 * group, 5 * group + 3)
            own_population_column = 5 * group + 4
            jacobian[:, :, component_columns] = (
                2.0 * coefficient / populations[:, group]
            )[:, None, None] * identity

            for population_group in range(2):
                population_column = 5 * population_group + 4
                jacobian[:, :, population_column] += (
                    coefficient_derivatives[:, group, population_group, None] * vector
                )
            # The normalized input vector itself changes with its population.
            jacobian[:, :, own_population_column] -= (
                coefficient / populations[:, group]
            )[:, None] * vector
        jacobians[name] = jacobian
    return jacobians


def compute_j_modes_from_j_moments(
    j_moments: JMomentSeries,
    *,
    omega_groups: Sequence[float],
    tol: float = 1e-12,
) -> JModeSeries:
    """Build active-population-normalized modes and their available SEMs.

    Propagate the full J sample covariance, including population and cross-group
    correlations. Population-dependent ratios use first-order Jacobians.
    Without source sampling statistics, return the same central modes with
    uncertainty and SEMs unset. The array-based compute_j_modes remains valid
    for callers that already supply normalized vectors.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")
    group_fields = (j_moments.x, j_moments.y, j_moments.z, j_moments.N_j)
    if any(field is None for field in group_fields):
        raise ValueError("J modes require group x/y/z components and N_j populations.")
    group_count = len(j_moments.x)
    if group_count not in (1, 2) or any(
        len(field) != group_count for field in group_fields
    ):
        raise ValueError("J modes require matching one- or two-group inputs.")

    t = np.asarray(j_moments.t, dtype=float)
    if t.ndim != 1 or t.size == 0 or not np.all(np.isfinite(t)):
        raise ValueError("J-mode times must be a nonempty finite one-dimensional array.")
    group_vectors = [
        np.stack((x, y, z), axis=-1)
        for x, y, z in zip(j_moments.x, j_moments.y, j_moments.z)
    ]
    vectors = np.stack(group_vectors, axis=1)
    populations = np.stack(j_moments.N_j, axis=-1)
    if (
        vectors.shape != (t.size, group_count, 3)
        or populations.shape != (t.size, group_count)
    ):
        raise ValueError("Group components and populations must match the saved-time grid.")
    if not np.all(np.isfinite(vectors)) or not np.all(np.isfinite(populations)):
        raise ValueError("J-mode components and populations must be finite.")
    if np.any(populations <= tol):
        raise ValueError("Active-population-normalized J modes require every N_j > tol.")
    couplings = np.asarray(omega_groups, dtype=float)
    if couplings.shape != (group_count,) or not np.all(np.isfinite(couplings)):
        raise ValueError("omega_groups must contain one finite coupling per group.")

    normalized_vectors = vectors * (2.0 / populations)[..., None]
    modes = compute_j_modes(
        t,
        tuple(normalized_vectors[:, group, 0] for group in range(group_count)),
        tuple(normalized_vectors[:, group, 1] for group in range(group_count)),
        tuple(normalized_vectors[:, group, 2] for group in range(group_count)),
        populations=tuple(populations[:, group] for group in range(group_count)),
        omega_groups=couplings,
    )
    source_uncertainty = getattr(j_moments, "uncertainty", None)
    if source_uncertainty is None:
        return modes
    covariance_size = 5 * group_count
    expected_covariance_shape = (t.size, covariance_size, covariance_size)
    if source_uncertainty.sample_covariance.shape != expected_covariance_shape:
        raise ValueError("J sample covariance must have shape (len(t), 5 * G, 5 * G).")

    jacobians = _j_mode_jacobians(normalized_vectors, populations, couplings)
    for name, jacobian in jacobians.items():
        covariance = propagate_covariance(
            source_uncertainty.sample_covariance, jacobian, tol=tol
        )
        vector = getattr(modes, name)
        vector_with_uncertainty = BlochVectorSeries(
            x=vector.x,
            y=vector.y,
            z=vector.z,
            uncertainty=SamplingUncertainty(
                sample_covariance=covariance,
                n_samples=source_uncertainty.n_samples,
            ),
        )
        setattr(modes, name, vector_with_uncertainty)
    return modes


__all__ = ["compute_j_modes", "compute_j_modes_from_j_moments"]
