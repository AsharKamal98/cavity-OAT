from __future__ import annotations

import numpy as np

from common.utils.phases import integration_phase_indices_at_times
from parser.j_moments import JMomentSeries
from parser.mfe import MFEResult

from solvers.mfe.utils import angles_from_amplitudes


def compute_mfe_j_moments(
    result: MFEResult,
    *,
    tol: float = 1e-12,
) -> JMomentSeries:
    """
    Build observable time series from an MFE result.
    """
    N_j, theta, phi = angles_from_amplitudes(
        result.D_groups,
        result.E_groups,
        tol=tol,
    )
    j_moments = JMomentSeries(
        t=result.t,
        integration_phase_index=integration_phase_indices_at_times(
            result.t,
            result.parameters.phase_protocol,
        ),
        N_e=tuple(np.abs(E_g) ** 2 for E_g in result.E_groups),
        N_j=N_j,
        length=tuple(0.5 * np.asarray(N_j_g, dtype=float) for N_j_g in N_j),
        theta=theta,
        phi=phi,
    )
    JMomentSeries.attach_norm_spin_components_from_angles(j_moments)
    JMomentSeries.attach_spin_components_from_norm_spin_components(j_moments)
    JMomentSeries.attach_aggregated_from_group_fields(j_moments, tol=tol)
    return j_moments
