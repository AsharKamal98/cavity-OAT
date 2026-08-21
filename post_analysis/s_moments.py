from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.j_moments import JMomentSeries
from parser.s_moments import SBasisSeries, SMomentSeries


def _validate_reference_data(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    j_moments: JMomentSeries,
    *,
    tol: float,
) -> tuple[Array, Array, Array, Array | None]:
    """Validate and return mean generators, J angles, and phase indices."""
    if mean_generator_moments.num_trajectories != 1:
        raise ValueError(
            "compute_s_moments requires already-averaged generator moments "
            "with a singleton trajectory axis."
        )

    t_ref = mean_generator_moments.t
    if np.asarray(j_moments.t).shape != t_ref.shape or not np.allclose(
        j_moments.t,
        t_ref,
        atol=tol,
        rtol=0.0,
    ):
        raise ValueError(
            "mean_generator_moments and j_moments must share a time grid."
        )
    if j_moments.theta is None or j_moments.phi is None:
        raise ValueError("S moments require group-resolved J angles.")

    mean_first_moments = mean_generator_moments.first_moments[0]
    num_times, num_groups, _ = mean_first_moments.shape
    theta_J_groups = np.asarray(j_moments.theta, dtype=float)
    phi_J_groups = np.asarray(j_moments.phi, dtype=float)
    expected_angle_shape = (num_groups, num_times)
    if theta_J_groups.shape != expected_angle_shape or (
        phi_J_groups.shape != expected_angle_shape
    ):
        raise ValueError(
            "J-angle fields must match the generator-moment time and group axes."
        )

    if j_moments.integration_phase_index is None:
        integration_phase_index = None
    else:
        integration_phase_index = np.asarray(
            j_moments.integration_phase_index,
            dtype=int,
        )
        if integration_phase_index.shape != (num_times,):
            raise ValueError(
                "j_moments.integration_phase_index must match the time grid."
            )

    return (
        mean_first_moments,
        theta_J_groups,
        phi_J_groups,
        integration_phase_index,
    )


def _build_s_basis(
    theta_J_groups: Array,
    phi_J_groups: Array,
) -> SBasisSeries:
    """Construct the group- and time-resolved effective S-qubit basis."""
    num_times = theta_J_groups.shape[1]
    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    zero_states = tuple(
        np.broadcast_to(ket_u, (num_times, 3)).copy()
        for _ in range(theta_J_groups.shape[0])
    )
    one_states = tuple(
        np.column_stack(
            (
                np.zeros(num_times, dtype=np.complex128),
                np.cos(theta_J / 2.0),
                np.exp(-1j * phi_J) * np.sin(theta_J / 2.0),
            )
        )
        for theta_J, phi_J in zip(theta_J_groups, phi_J_groups)
    )
    return SBasisSeries(
        zero_state=zero_states,
        one_state=one_states,
    )


def _build_s_generator_coefficients(
    theta_J_groups: Array,
    phi_J_groups: Array,
) -> Array:
    """Return ``(time, group, S component, generator)`` coefficients."""
    theta_J = np.asarray(theta_J_groups, dtype=float).T
    phi_J = np.asarray(phi_J_groups, dtype=float).T
    cos_half_theta = np.cos(theta_J / 2.0)
    sin_half_theta = np.sin(theta_J / 2.0)
    cos_phi = np.cos(phi_J)
    sin_phi = np.sin(phi_J)

    coefficients = np.zeros(
        theta_J.shape + (3, 8),
        dtype=float,
    )

    coefficients[:, :, 0, 0] = cos_half_theta
    coefficients[:, :, 0, 3] = sin_half_theta * cos_phi
    coefficients[:, :, 0, 4] = -sin_half_theta * sin_phi

    coefficients[:, :, 1, 1] = -cos_half_theta
    coefficients[:, :, 1, 3] = -sin_half_theta * sin_phi
    coefficients[:, :, 1, 4] = -sin_half_theta * cos_phi

    coefficients[:, :, 2, 2] = 0.5 * (1.0 + cos_half_theta**2)
    coefficients[:, :, 2, 5] = -0.5 * np.sin(theta_J) * cos_phi
    coefficients[:, :, 2, 6] = 0.5 * np.sin(theta_J) * sin_phi
    coefficients[:, :, 2, 7] = (
        0.5 * np.sqrt(3.0) * sin_half_theta**2
    )
    return coefficients


def _group_series(values: Array) -> tuple[Array, ...]:
    """Split a ``(time, group)`` array into independent group series."""
    return tuple(
        values[:, group_index].copy()
        for group_index in range(values.shape[1])
    )


def compute_s_moments(
    mean_generator_moments: TrajectoryResolvedMomentSeries,
    j_moments: JMomentSeries,
    *,
    tol: float = 1e-12,
) -> SMomentSeries:
    """Construct group-resolved S moments from mean SU(3) generator moments.

    The J angles define one shared S frame for each group and time. Because the
    resulting generator transformation is linear and shared across the selected
    trajectories, it acts directly on their already-averaged first moments.
    Symmetric second moments are intentionally left for later covariance and
    squeezing analysis.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")

    (
        mean_first_moments,
        theta_J_groups,
        phi_J_groups,
        integration_phase_index,
    ) = _validate_reference_data(
        mean_generator_moments,
        j_moments,
        tol=tol,
    )
    coefficients = _build_s_generator_coefficients(
        theta_J_groups,
        phi_J_groups,
    )
    components = np.einsum(
        "tgsa,tga->tgs",
        coefficients,
        mean_first_moments,
        optimize=True,
    )

    s_moments = SMomentSeries(
        t=mean_generator_moments.t.copy(),
        integration_phase_index=(
            None
            if integration_phase_index is None
            else integration_phase_index.copy()
        ),
        x=_group_series(components[:, :, 0]),
        y=_group_series(components[:, :, 1]),
        z=_group_series(components[:, :, 2]),
        basis=_build_s_basis(theta_J_groups, phi_J_groups),
    )
    SMomentSeries.attach_norm_spin_components_from_spin_components(
        s_moments,
        tol=tol,
    )
    SMomentSeries.attach_angles_from_norm_spin_components(
        s_moments,
        tol=tol,
    )
    return s_moments


__all__ = ["compute_s_moments"]
