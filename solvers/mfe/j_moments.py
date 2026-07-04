from __future__ import annotations

import numpy as np

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
        N_j_groups=N_j_groups,
        length_groups=tuple(0.5 * np.asarray(N_j, dtype=float) for N_j in N_j_groups),
        theta_groups=theta_groups,
        phi_groups=phi_groups,
    )
    JMomentSeries.attatch_norm_spin_components_from_angles(j_moments)
    JMomentSeries.attatch_spin_components_from_norm_spin_components(j_moments)
    return j_moments
