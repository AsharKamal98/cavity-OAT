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
    N_j_groups, theta_groups, phi_groups = angles_from_amplitudes(
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
        N_e_groups=tuple(np.abs(E_g) ** 2 for E_g in result.E_groups),
        N_j_groups=N_j_groups,
        length_groups=tuple(0.5 * np.asarray(N_j, dtype=float) for N_j in N_j_groups),
        theta_groups=theta_groups,
        phi_groups=phi_groups,
    )
    JMomentSeries.attatch_norm_spin_components_from_angles(j_moments)
    JMomentSeries.attatch_spin_components_from_norm_spin_components(j_moments)
    JMomentSeries.attatch_additive_full_fields_from_group_fields(j_moments)
    JMomentSeries.attatch_norm_spin_components_from_spin_components(j_moments, tol=tol)
    JMomentSeries.attatch_angles_from_norm_spin_components(j_moments, tol=tol)
    return j_moments
