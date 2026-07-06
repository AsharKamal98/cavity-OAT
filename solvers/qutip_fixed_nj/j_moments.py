from __future__ import annotations

import numpy as np

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

    x = _mean_series(result.expect[0])
    y = _mean_series(result.expect[1])
    z = _mean_series(result.expect[2])
    N_e = _mean_series(result.expect[3])
    N_j = np.full_like(t, float(getattr(model, "NJ")), dtype=float)

    x_groups = y_groups = z_groups = N_e_groups = N_j_groups = None
    if hasattr(model, "Jx_groups") and len(result.expect) >= 12:
        x_groups = (
            _mean_series(result.expect[4]),
            _mean_series(result.expect[5]),
        )
        y_groups = (
            _mean_series(result.expect[6]),
            _mean_series(result.expect[7]),
        )
        z_groups = (
            _mean_series(result.expect[8]),
            _mean_series(result.expect[9]),
        )
        N_e_groups = (
            _mean_series(result.expect[10]),
            _mean_series(result.expect[11]),
        )
        N_j_groups = (
            np.full_like(t, float(getattr(model, "NJ1")), dtype=float),
            np.full_like(t, float(getattr(model, "NJ2")), dtype=float),
        )

    j_moments = JMomentSeries(
        t=t,
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
