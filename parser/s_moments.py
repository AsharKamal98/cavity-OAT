from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict

from parser.common import Array


class SMomentSnapshot(BaseModel):
    """First-order S-sphere moments for one saved MCWF trajectory snapshot."""

    t: float
    integration_phase_index: int
    x: tuple[float, ...]
    y: tuple[float, ...]
    z: tuple[float, ...]


class SMomentSeries(BaseModel):
    """Group-resolved first-order S-sphere moment time series."""

    t: Array
    integration_phase_index: Array | None = None
    x: tuple[Array, ...] | None = None
    y: tuple[Array, ...] | None = None
    z: tuple[Array, ...] | None = None
    length: tuple[Array, ...] | None = None
    nx: tuple[Array, ...] | None = None
    ny: tuple[Array, ...] | None = None
    nz: tuple[Array, ...] | None = None
    theta: tuple[Array, ...] | None = None
    phi: tuple[Array, ...] | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @classmethod
    def attach_norm_spin_components_from_spin_components(
        cls,
        s_moments: "SMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """Attach S-vector lengths and normalized group directions."""
        from common.utils.moments import norm_spin_components_from_spin_components

        if (
            s_moments.x is None
            or s_moments.y is None
            or s_moments.z is None
        ):
            raise ValueError(
                "Spin-direction attachment requires group spin components."
            )

        group_results = [
            norm_spin_components_from_spin_components(x_g, y_g, z_g, tol=tol)
            for x_g, y_g, z_g in zip(
                s_moments.x,
                s_moments.y,
                s_moments.z,
            )
        ]
        s_moments.length = tuple(result[0] for result in group_results)
        s_moments.nx = tuple(result[1] for result in group_results)
        s_moments.ny = tuple(result[2] for result in group_results)
        s_moments.nz = tuple(result[3] for result in group_results)

    @classmethod
    def attach_angles_from_norm_spin_components(
        cls,
        s_moments: "SMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """Attach S-sphere angles using the positive-z pole convention."""
        from common.utils.moments import angles_from_norm_spin_components

        if (
            s_moments.length is None
            or s_moments.nx is None
            or s_moments.ny is None
            or s_moments.nz is None
        ):
            raise ValueError(
                "Spin angle attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            angles_from_norm_spin_components(
                nx_g,
                ny_g,
                -nz_g,
                valid=length_g > tol,
                tol=tol,
            )
            for length_g, nx_g, ny_g, nz_g in zip(
                s_moments.length,
                s_moments.nx,
                s_moments.ny,
                s_moments.nz,
            )
        ]
        s_moments.theta = tuple(result[0] for result in group_results)
        s_moments.phi = tuple(result[1] for result in group_results)


__all__ = ["SMomentSeries", "SMomentSnapshot"]
