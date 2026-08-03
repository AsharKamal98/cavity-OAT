from __future__ import annotations

import numpy as np

from parser.squeezing import SqueezingMomentSeries, SqueezingSeries


def _compute_group_squeezing(
    moments: SqueezingMomentSeries,
    *,
    tol: float,
) -> SqueezingSeries:
    """Compute generalized squeezing from one group's averaged raw moments."""
    covariance = moments.second_moments - (
        moments.mean_operators[:, :, None]
        * moments.mean_operators[:, None, :]
    )
    covariance_eigvals = np.linalg.eigvalsh(covariance)
    lambda_min = covariance_eigvals[:, 0]

    denominator = (moments.N_c / 2.0) ** 2
    xi2 = np.full_like(moments.t, np.nan, dtype=float)
    valid_denominator = denominator > tol
    xi2[valid_denominator] = (
        moments.N_group
        * lambda_min[valid_denominator]
        / denominator[valid_denominator]
    )

    xi2_db = np.full_like(xi2, np.nan, dtype=float)
    valid_xi2 = xi2 > 0.0
    xi2_db[valid_xi2] = 10.0 * np.log10(xi2[valid_xi2])

    return SqueezingSeries(
        t=moments.t,
        N_group=moments.N_group,
        xi2=xi2,
        xi2_db=xi2_db,
        lambda_min=lambda_min,
        covariance_eigvals=covariance_eigvals,
        N_c=moments.N_c,
        excited_fraction_active=moments.excited_fraction_active,
        theta_J=moments.theta_J,
        phi_J=moments.phi_J,
        theta_S=moments.theta_S,
        phi_S=moments.phi_S,
    )


def compute_squeezing(
    squeezing_moments: tuple[SqueezingMomentSeries, ...],
    *,
    tol: float = 1e-12,
) -> tuple[SqueezingSeries, ...]:
    """Compute one generalized squeezing series per physical group."""
    if not squeezing_moments:
        raise ValueError("At least one SqueezingMomentSeries is required.")
    if tol <= 0.0:
        raise ValueError("tol must be positive.")

    t_ref = np.asarray(squeezing_moments[0].t, dtype=float)
    for moments in squeezing_moments[1:]:
        if moments.t.shape != t_ref.shape or not np.allclose(
            moments.t,
            t_ref,
            atol=1e-12,
            rtol=0.0,
        ):
            raise ValueError(
                "All squeezing moment series must share the same time grid."
            )

    return tuple(
        _compute_group_squeezing(moments, tol=tol)
        for moments in squeezing_moments
    )


__all__ = ["compute_squeezing"]
