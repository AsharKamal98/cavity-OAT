from __future__ import annotations

import numpy as np

from parser.j_moments import JMomentSeries
from parser.mfe_residuals import MFEResidualSeries
from parser.moments import SimulationMetadata


def compute_mfe_residuals(
    j_moments: JMomentSeries,
    *,
    metadata: SimulationMetadata | None,
    tol: float = 1e-12,
) -> MFEResidualSeries | None:
    """
    Compute two-group MFE residuals from averaged J-vector angles.
    """
    if (
        j_moments.theta_groups is None
        or j_moments.phi_groups is None
        or j_moments.N_j_groups is None
    ):
        return None

    group_count = len(j_moments.theta_groups)
    if group_count != 2:
        return None
    if (
        len(j_moments.phi_groups) != group_count
        or len(j_moments.N_j_groups) != group_count
    ):
        raise ValueError("MFE residuals require matching two-group moment fields.")
    if metadata is None:
        raise ValueError("MFE residuals require moments.metadata.")
    omega_groups = metadata.omega_groups
    if len(omega_groups) != group_count:
        raise ValueError("MFE residuals require two inhomogeneous coupling weights.")

    phase_indices = np.asarray(j_moments.phase_index, dtype=int)
    omega_t = np.asarray([metadata.phases[idx].omega for idx in phase_indices], dtype=float)
    delta_t = np.asarray([metadata.phases[idx].delta for idx in phase_indices], dtype=float)

    theta_groups = tuple(np.asarray(theta, dtype=float) for theta in j_moments.theta_groups)
    phi_groups = tuple(np.asarray(phi, dtype=float) for phi in j_moments.phi_groups)
    nj_groups = tuple(np.asarray(nj, dtype=float) for nj in j_moments.N_j_groups)
    omega_groups = tuple(float(omega_g) for omega_g in omega_groups)

    weighted_collective_transverse_sum = sum(
        omega_g * nj_g * np.exp(1j * phi_g) * np.sin(theta_g)
        for theta_g, phi_g, nj_g, omega_g in zip(
            theta_groups,
            phi_groups,
            nj_groups,
            omega_groups,
        )
    )

    residuals = []
    for theta_g, phi_g, omega_g in zip(theta_groups, phi_groups, omega_groups):
        sin_theta = np.sin(theta_g)
        cos_theta = np.cos(theta_g)

        with np.errstate(divide="ignore", invalid="ignore"):
            detuning_factor = np.where(
                np.abs(cos_theta) > tol,
                sin_theta * np.tan(theta_g),
                np.nan,
            )

        drive_term = 0.5 * omega_t * omega_g * np.exp(-1j * phi_g) * sin_theta
        detuning_term = -0.5 * delta_t * detuning_factor
        decay_term = (
            0.25j
            * metadata.Gamma
            * omega_g
            * np.exp(-1j * phi_g)
            * sin_theta
            * weighted_collective_transverse_sum
        )
        residuals.append(drive_term + detuning_term + decay_term)

    return MFEResidualSeries(
        t=j_moments.t,
        phase_index=j_moments.phase_index,
        residuals_groups=tuple(residuals),
    )


__all__ = [
    "compute_mfe_residuals",
]
