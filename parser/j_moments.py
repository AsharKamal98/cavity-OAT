from __future__ import annotations

from pydantic import BaseModel

from parser.common import Array


class JMomentSnapshot(BaseModel):
    """First-order J-sphere moments for one saved trajectory snapshot."""

    # Legacy note: these fields were previously named Jx/Jy/Jz and
    # Jx_groups/Jy_groups/Jz_groups.
    t: float
    phase_index: int
    x: float
    y: float
    z: float
    N_e: float
    N_j: float
    jump_rate: float
    x_groups: tuple[float, ...] | None = None
    y_groups: tuple[float, ...] | None = None
    z_groups: tuple[float, ...] | None = None
    N_e_groups: tuple[float, ...] | None = None
    N_j_groups: tuple[float, ...] | None = None


class JMomentSeries(BaseModel):
    """Per-timestep first-order J-sphere moments for one (single or averaged) trajectory."""
    # Legacy note: these fields were previously named Jx/Jy/Jz, Jx_groups/
    # Jy_groups/Jz_groups, J_len, and sx/sy/sz.

    t: Array
    phase_index: Array | None = None
    # spin components
    x: Array | None = None
    y: Array | None = None
    z: Array | None = None
    x_groups: tuple[Array, ...] | None = None
    y_groups: tuple[Array, ...] | None = None
    z_groups: tuple[Array, ...] | None = None
    length: Array | None = None
    length_groups: tuple[Array, ...] | None = None
    # normalized spin components / directions
    nx: Array | None = None
    ny: Array | None = None
    nz: Array | None = None
    nx_groups: tuple[Array, ...] | None = None
    ny_groups: tuple[Array, ...] | None = None
    nz_groups: tuple[Array, ...] | None = None
    # atom numbers
    N_e: Array | None = None
    N_j: Array | None = None
    N_e_groups: tuple[Array, ...] | None = None
    N_j_groups: tuple[Array, ...] | None = None
    # angles
    theta: Array | None = None
    phi: Array | None = None
    theta_groups: tuple[Array, ...] | None = None
    phi_groups: tuple[Array, ...] | None = None
    # other
    jump_rate: Array | None = None

    @classmethod
    def attatch_norm_spin_components_from_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach Euclidean vector length and normalized spin-direction fields.
        """
        from common.utils_moments import norm_spin_components_from_spin_components

        has_full = (
            j_moments.x is not None
            and j_moments.y is not None
            and j_moments.z is not None
        )
        has_groups = (
            j_moments.x_groups is not None
            and j_moments.y_groups is not None
            and j_moments.z_groups is not None
        )

        if not has_full and not has_groups:
            raise ValueError("Spin-direction attachment requires full spin components or group-resolved spin components.")

        if has_full:
            (
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            ) = norm_spin_components_from_spin_components(
                j_moments.x,
                j_moments.y,
                j_moments.z,
                tol=tol,
            )

        if not has_groups:
            return

        group_results = [
            norm_spin_components_from_spin_components(x_g, y_g, z_g, tol=tol)
            for x_g, y_g, z_g in zip(
                j_moments.x_groups,
                j_moments.y_groups,
                j_moments.z_groups,
            )
        ]
        j_moments.length_groups = tuple(result[0] for result in group_results)
        j_moments.nx_groups = tuple(result[1] for result in group_results)
        j_moments.ny_groups = tuple(result[2] for result in group_results)
        j_moments.nz_groups = tuple(result[3] for result in group_results)

    @classmethod
    def attatch_angles_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach theta and phi from already-normalized spin-direction fields.
        """
        from common.utils_moments import angles_from_norm_spin_components

        if (
            j_moments.length is None
            or j_moments.nx is None
            or j_moments.ny is None
            or j_moments.nz is None
        ):
            raise ValueError("Spin angle attachment requires normalized full spin-direction fields.")

        valid = j_moments.length > tol
        j_moments.theta, j_moments.phi = angles_from_norm_spin_components(
            j_moments.nx,
            j_moments.ny,
            j_moments.nz,
            valid=valid,
            tol=tol,
        )

        if (
            j_moments.length_groups is None
            or j_moments.nx_groups is None
            or j_moments.ny_groups is None
            or j_moments.nz_groups is None
        ):
            return

        group_results = [
            angles_from_norm_spin_components(
                nx_g,
                ny_g,
                nz_g,
                valid=length_g > tol,
                tol=tol,
            )
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length_groups,
                j_moments.nx_groups,
                j_moments.ny_groups,
                j_moments.nz_groups,
            )
        ]
        j_moments.theta_groups = tuple(result[0] for result in group_results)
        j_moments.phi_groups = tuple(result[1] for result in group_results)

    @classmethod
    def attatch_norm_spin_components_from_angles(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach normalized spin-direction fields from theta and phi.
        """
        from common.utils_moments import norm_spin_components_from_angles

        has_full = (
            j_moments.theta is not None
            and j_moments.phi is not None
        )
        has_groups = (
            j_moments.theta_groups is not None
            and j_moments.phi_groups is not None
        )

        if not has_full and not has_groups:
            raise ValueError("Normalized spin-component attachment requires full angles or group-resolved angles.")

        if has_full:
            (
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            ) = norm_spin_components_from_angles(
                j_moments.theta,
                j_moments.phi,
            )

        if not has_groups:
            return

        group_results = [
            norm_spin_components_from_angles(theta_g, phi_g)
            for theta_g, phi_g in zip(
                j_moments.theta_groups,
                j_moments.phi_groups,
            )
        ]
        j_moments.nx_groups = tuple(result[0] for result in group_results)
        j_moments.ny_groups = tuple(result[1] for result in group_results)
        j_moments.nz_groups = tuple(result[2] for result in group_results)

    @classmethod
    def attatch_spin_components_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach spin-component fields from vector length and normalized direction.
        """
        from common.utils_moments import spin_components_from_norm_spin_components

        has_full = (
            j_moments.length is not None
            and j_moments.nx is not None
            and j_moments.ny is not None
            and j_moments.nz is not None
        )
        has_groups = (
            j_moments.length_groups is not None
            and j_moments.nx_groups is not None
            and j_moments.ny_groups is not None
            and j_moments.nz_groups is not None
        )

        if not has_full and not has_groups:
            raise ValueError("Spin-component attachment requires full normalized spin-direction fields or group-resolved normalized spin-direction fields.")

        if has_full:
            (
                j_moments.x,
                j_moments.y,
                j_moments.z,
            ) = spin_components_from_norm_spin_components(
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            )

        if not has_groups:
            return

        group_results = [
            spin_components_from_norm_spin_components(length_g, nx_g, ny_g, nz_g)
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length_groups,
                j_moments.nx_groups,
                j_moments.ny_groups,
                j_moments.nz_groups,
            )
        ]
        j_moments.x_groups = tuple(result[0] for result in group_results)
        j_moments.y_groups = tuple(result[1] for result in group_results)
        j_moments.z_groups = tuple(result[2] for result in group_results)

    class Config:
        arbitrary_types_allowed = True
