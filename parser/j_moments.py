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
    def attach_spin_direction_fields(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach Euclidean vector length and normalized spin-direction fields.
        """
        from common.moment_utils import norm_spin_components

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
            ) = norm_spin_components(j_moments.x, j_moments.y, j_moments.z, tol=tol)

        if not has_groups:
            return

        group_results = [
            norm_spin_components(x_g, y_g, z_g, tol=tol)
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

    class Config:
        arbitrary_types_allowed = True
