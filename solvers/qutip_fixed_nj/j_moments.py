from __future__ import annotations

import numpy as np

from common.utils.phases import integration_phase_indices_at_times
from parser.j_moments import JMomentSeries


def _mean_series(series) -> np.ndarray:
    """
    Convert a QuTiP expectation-value series into one real-valued time series.

    mesolve typically returns a 1D array per observable. mcsolve may return a
    2D array with one row per trajectory. In the latter case, average over
    trajectories.
    """
    arr = np.real(np.asarray(series, dtype=float))
    if arr.ndim == 2:
        return np.mean(arr, axis=0)
    return arr


def compute_qutip_j_moments(
    sim_data,
    *,
    tol: float = 1e-12,
) -> JMomentSeries:
    """
    Convert QuTiP mesolve or mcsolve output into a JMomentSeries.

    Only quantities already present in the QuTiP expectation-value output or in
    the fixed-sector model metadata are attached here. No derived fields such
    as angles, normalized directions, or jump-rate series are computed.
    """
    _ = tol  # Kept for API symmetry with the other moment extractors.

    result = sim_data["result"]
    model = sim_data["model"]
    t = np.asarray(sim_data["tlist"], dtype=float)
    NJi = tuple(int(NJ) for NJ in model.NJi)
    group_count = len(NJi)

    x = _mean_series(result.expect[0])
    y = _mean_series(result.expect[1])
    z = _mean_series(result.expect[2])
    N_e = _mean_series(result.expect[3])
    N_j = np.full_like(t, float(sum(NJi)), dtype=float)

    if group_count == 1:
        x_groups = (x,)
        y_groups = (y,)
        z_groups = (z,)
        N_e_groups = (N_e,)
    else:
        expected_group_fields = 4 * group_count
        if len(result.expect) < 4 + expected_group_fields:
            raise ValueError("QuTiP result does not contain the expected group-resolved observables.")

        offset = 4
        x_groups = tuple(_mean_series(result.expect[offset + g]) for g in range(group_count))
        offset += group_count
        y_groups = tuple(_mean_series(result.expect[offset + g]) for g in range(group_count))
        offset += group_count
        z_groups = tuple(_mean_series(result.expect[offset + g]) for g in range(group_count))
        offset += group_count
        N_e_groups = tuple(_mean_series(result.expect[offset + g]) for g in range(group_count))

    N_j_groups = tuple(
        np.full_like(t, float(NJ), dtype=float)
        for NJ in NJi
    )

    j_moments = JMomentSeries(
        t=t,
        integration_phase_index=integration_phase_indices_at_times(
            t,
            model.phase_protocol,
        ),
        x=x,
        y=y,
        z=z,
        x_groups=x_groups,
        y_groups=y_groups,
        z_groups=z_groups,
        N_e=N_e,
        N_j=N_j,
        N_e_groups=N_e_groups,
        N_j_groups=N_j_groups,
    )
    JMomentSeries.attatch_norm_spin_components_from_spin_components(j_moments, tol=tol)
    JMomentSeries.attatch_angles_from_norm_spin_components(j_moments, tol=tol)
    return j_moments


__all__ = [
    "compute_qutip_j_moments",
]
