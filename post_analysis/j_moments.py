from __future__ import annotations

import numpy as np

from common.utils.phases import integration_phase_indices_at_times
from common.utils.uncertainty import propagate_covariance, std_from_covariance
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.j_moments import JMomentSeries
from parser.moments import SimulationMetadata
from parser.uncertainty import SamplingUncertainty


# Fixed lab-frame map: (J_x, J_y, J_z)^T = A (G_1, ..., G_8)^T.
# Generators use g_a = lambda_a / 2 in the (|u>, |d>, |e>) basis.
# J_y=(J_+-J_-)/(2i) gives the minus sign on G_7.
J_FROM_GENERATORS = np.array(
    [
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0],
        [0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, -np.sqrt(3.0) / 2.0],
    ],
    dtype=float,
)
J_FROM_GENERATORS.setflags(write=False)

# N_e = N_group / 3 - 2 G_8 / sqrt(3): number of excited atoms.
_EXCITED_POPULATION_FROM_GENERATORS = np.array(
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -2.0 / np.sqrt(3.0)]
)
_EXCITED_POPULATION_FROM_GENERATORS.setflags(write=False)

# N_j = 2 N_group / 3 - G_3 - G_8 / sqrt(3): active-manifold population.
_ACTIVE_POPULATION_FROM_GENERATORS = np.array(
    [0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0 / np.sqrt(3.0)]
)
_ACTIVE_POPULATION_FROM_GENERATORS.setflags(write=False)

# Covariance rows follow this same physical-field order within each group.
_J_FIELD_NAMES = ("x", "y", "z", "N_e", "N_j")
_J_FIELDS_FROM_GENERATORS = np.vstack(
    (
        J_FROM_GENERATORS,
        _EXCITED_POPULATION_FROM_GENERATORS,
        _ACTIVE_POPULATION_FROM_GENERATORS,
    )
)
_J_FIELDS_FROM_GENERATORS.setflags(write=False)


def _group_series(values: np.ndarray) -> tuple[np.ndarray, ...]:
    """Split a ``(time, group)`` array into independent group series."""
    return tuple(
        values[:, group_index].copy()
        for group_index in range(values.shape[1])
    )


def compute_j_components(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *,
    Ni: tuple[int, ...],
    tol: float = 1e-12,
) -> JMomentSeries:
    """Construct affine J fields and their group/full-system sampling SEMs.

    Input first moments are already averaged and use ``g_a = lambda_a / 2``
    in the lab basis. Optional sampling covariance is propagated with the same
    coefficients used for the values. The result retains full sample covariance
    and sample count, plus SEMs (not standard deviations) for x/y/z/N_e/N_j.
    Group direction/angle fields and phase indices are added by
    ``compute_j_moments``. Quantum second moments are not consumed.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")
    if mean_generator_moments.num_trajectories != 1:
        raise ValueError(
            "compute_j_components requires already-averaged generator moments "
            "with a singleton trajectory axis."
        )

    first_moments = mean_generator_moments.first_moments[0]
    group_count = first_moments.shape[1]
    group_sizes = np.asarray(Ni, dtype=float)
    if group_sizes.shape != (group_count,):
        raise ValueError(
            "Generator moments must contain one group-axis entry for each "
            "group in Ni."
        )
    if not np.all(np.isfinite(group_sizes)) or np.any(group_sizes < 0.0):
        raise ValueError("Ni must contain finite, nonnegative group sizes.")

    # Compute the Cartesian spin vector and the two physical populations.
    spin_components = first_moments @ J_FROM_GENERATORS.T
    x = spin_components[:, :, 0]
    y = spin_components[:, :, 1]
    z = spin_components[:, :, 2]
    excited_population = (
        group_sizes / 3.0
        + first_moments @ _EXCITED_POPULATION_FROM_GENERATORS
    )
    active_population = (
        2.0 * group_sizes / 3.0
        + first_moments @ _ACTIVE_POPULATION_FROM_GENERATORS
    )

    j_moments = JMomentSeries(
        t=mean_generator_moments.t.copy(),
        x=_group_series(x),
        y=_group_series(y),
        z=_group_series(z),
        N_e=_group_series(excited_population),
        N_j=_group_series(active_population),
    )
    source_uncertainty = getattr(mean_generator_moments, "uncertainty", None)
    if source_uncertainty is not None:
        n_samples = source_uncertainty.n_samples
        sqrt_sample_count = np.sqrt(n_samples)
        num_times = len(mean_generator_moments.t)
        num_fields = len(_J_FIELD_NAMES)

        # Transform every group's generators into (x, y, z, N_e, N_j), keeping
        # correlations between different groups in the full covariance matrix.
        group_transformation = np.kron(np.eye(group_count), _J_FIELDS_FROM_GENERATORS)
        group_sample_covariance = propagate_covariance(
            source_uncertainty.sample_covariance, group_transformation, tol=tol
        )
        group_std = std_from_covariance(group_sample_covariance, tol=tol)
        group_sem = (group_std / sqrt_sample_count).reshape(
            num_times, group_count, num_fields
        )

        j_moments.x_sem = _group_series(group_sem[:, :, 0])
        j_moments.y_sem = _group_series(group_sem[:, :, 1])
        j_moments.z_sem = _group_series(group_sem[:, :, 2])
        j_moments.N_e_sem = _group_series(group_sem[:, :, 3])
        j_moments.N_j_sem = _group_series(group_sem[:, :, 4])
        j_moments.uncertainty = SamplingUncertainty(
            sample_covariance=group_sample_covariance,
            n_samples=n_samples,
        )

    JMomentSeries.attach_aggregated_from_group_fields(j_moments, tol=tol)
    return j_moments


def compute_j_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    *,
    metadata: SimulationMetadata,
    tol: float = 1e-12,
) -> JMomentSeries:
    """Construct J fields, optional affine-field SEMs, directions and angles."""
    j_moments = compute_j_components(
        mean_generator_moments, Ni=metadata.Ni, tol=tol
    )
    j_moments.integration_phase_index = integration_phase_indices_at_times(
        mean_generator_moments.t, metadata.phase_protocol
    )
    JMomentSeries.attach_norm_spin_components_from_spin_components(
        j_moments,
        tol=tol,
    )
    JMomentSeries.attach_angles_from_norm_spin_components(
        j_moments,
        tol=tol,
    )
    return j_moments


__all__ = ["J_FROM_GENERATORS", "compute_j_components", "compute_j_moments"]
