from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from parser.common import Array


class AggregatedJMomentSnapshot(BaseModel):
    """Full-system J moments for one saved trajectory snapshot."""

    x: float
    y: float
    z: float
    N_e: float
    N_j: float


class JMomentSnapshot(BaseModel):
    """First-order J-sphere moments for one saved trajectory snapshot."""

    # Legacy note: these fields were previously named Jx/Jy/Jz and
    # Jx_groups/Jy_groups/Jz_groups.
    t: float
    integration_phase_index: int
    x: tuple[float, ...]
    y: tuple[float, ...]
    z: tuple[float, ...]
    N_e: tuple[float, ...]
    N_j: tuple[float, ...]
    aggregated: AggregatedJMomentSnapshot
    jump_rate: float


class AggregatedJMomentSeries(BaseModel):
    """Full-system J-moment time series."""

    x: Array | None = None
    y: Array | None = None
    z: Array | None = None
    length: Array | None = None
    nx: Array | None = None
    ny: Array | None = None
    nz: Array | None = None
    N_e: Array | None = None
    N_j: Array | None = None
    theta: Array | None = None
    phi: Array | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)


class JMomentSeries(BaseModel):
    """Per-timestep first-order J-sphere moments for one (single or averaged) trajectory."""
    # Legacy note: these fields were previously named Jx/Jy/Jz, Jx_groups/
    # Jy_groups/Jz_groups, J_len, and sx/sy/sz.

    t: Array
    integration_phase_index: Array | None = None
    # spin components
    x: tuple[Array, ...] | None = None
    y: tuple[Array, ...] | None = None
    z: tuple[Array, ...] | None = None
    length: tuple[Array, ...] | None = None
    # normalized spin components / directions
    nx: tuple[Array, ...] | None = None
    ny: tuple[Array, ...] | None = None
    nz: tuple[Array, ...] | None = None
    # atom numbers
    N_e: tuple[Array, ...] | None = None
    N_j: tuple[Array, ...] | None = None
    # angles
    theta: tuple[Array, ...] | None = None
    phi: tuple[Array, ...] | None = None
    # full-system fields
    aggregated: AggregatedJMomentSeries = Field(
        default_factory=AggregatedJMomentSeries
    )
    # other
    jump_rate: Array | None = None

    @classmethod
    def attach_aggregated_from_group_fields(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach all aggregated fields from the primary group-resolved fields.
        """
        from common.utils.moments import (
            angles_from_norm_spin_components,
            norm_spin_components_from_spin_components,
        )

        field_names = ("x", "y", "z", "N_e", "N_j")
        group_fields = {
            field_name: getattr(j_moments, field_name)
            for field_name in field_names
        }
        missing_fields = [
            field_name
            for field_name, groups in group_fields.items()
            if groups is None
        ]
        if missing_fields:
            raise ValueError(
                "Additive aggregated-field attachment requires "
                + ", ".join(missing_fields)
                + "."
            )

        group_count = len(group_fields["x"])
        if group_count == 0 or any(
            len(groups) != group_count for groups in group_fields.values()
        ):
            raise ValueError(
                "All additive group fields must have the same nonzero group count."
            )

        aggregated_fields = {}
        for field_name, groups in group_fields.items():
            group_arrays = [np.asarray(group, dtype=float) for group in groups]
            aggregated_fields[field_name] = np.sum(
                np.stack(group_arrays),
                axis=0,
            )

        length, nx, ny, nz = norm_spin_components_from_spin_components(
            aggregated_fields["x"],
            aggregated_fields["y"],
            aggregated_fields["z"],
            tol=tol,
        )
        theta, phi = angles_from_norm_spin_components(
            nx,
            ny,
            nz,
            valid=length > tol,
            tol=tol,
        )
        j_moments.aggregated = AggregatedJMomentSeries(
            **aggregated_fields,
            length=length,
            nx=nx,
            ny=ny,
            nz=nz,
            theta=theta,
            phi=phi,
        )

    @classmethod
    def attach_norm_spin_components_from_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach Euclidean vector length and normalized spin-direction fields.
        """
        from common.utils.moments import norm_spin_components_from_spin_components

        if (
            j_moments.x is None
            or j_moments.y is None
            or j_moments.z is None
        ):
            raise ValueError(
                "Spin-direction attachment requires group spin components."
            )

        group_results = [
            norm_spin_components_from_spin_components(x_g, y_g, z_g, tol=tol)
            for x_g, y_g, z_g in zip(
                j_moments.x,
                j_moments.y,
                j_moments.z,
            )
        ]
        j_moments.length = tuple(result[0] for result in group_results)
        j_moments.nx = tuple(result[1] for result in group_results)
        j_moments.ny = tuple(result[2] for result in group_results)
        j_moments.nz = tuple(result[3] for result in group_results)

    @classmethod
    def attach_angles_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach theta and phi from already-normalized spin-direction fields.
        """
        from common.utils.moments import angles_from_norm_spin_components

        if (
            j_moments.length is None
            or j_moments.nx is None
            or j_moments.ny is None
            or j_moments.nz is None
        ):
            raise ValueError(
                "Spin angle attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            angles_from_norm_spin_components(
                nx_g,
                ny_g,
                nz_g,
                valid=length_g > tol,
                tol=tol,
            )
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            )
        ]
        j_moments.theta = tuple(result[0] for result in group_results)
        j_moments.phi = tuple(result[1] for result in group_results)

    @classmethod
    def attach_norm_spin_components_from_angles(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach normalized spin-direction fields from theta and phi.
        """
        from common.utils.moments import norm_spin_components_from_angles

        if j_moments.theta is None or j_moments.phi is None:
            raise ValueError(
                "Normalized spin-component attachment requires group angles."
            )

        group_results = [
            norm_spin_components_from_angles(theta_g, phi_g)
            for theta_g, phi_g in zip(
                j_moments.theta,
                j_moments.phi,
            )
        ]
        j_moments.nx = tuple(result[0] for result in group_results)
        j_moments.ny = tuple(result[1] for result in group_results)
        j_moments.nz = tuple(result[2] for result in group_results)

    @classmethod
    def attach_spin_components_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach spin-component fields from vector length and normalized direction.
        """
        from common.utils.moments import spin_components_from_norm_spin_components

        if (
            j_moments.length is None
            or j_moments.nx is None
            or j_moments.ny is None
            or j_moments.nz is None
        ):
            raise ValueError(
                "Spin-component attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            spin_components_from_norm_spin_components(length_g, nx_g, ny_g, nz_g)
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            )
        ]
        j_moments.x = tuple(result[0] for result in group_results)
        j_moments.y = tuple(result[1] for result in group_results)
        j_moments.z = tuple(result[2] for result in group_results)

    class Config:
        arbitrary_types_allowed = True
