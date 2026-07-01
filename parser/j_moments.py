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
    J_drive: float
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
    phase_index: Array
    # spin components
    x: Array
    y: Array
    z: Array
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
    N_e: Array
    N_j: Array
    N_e_groups: tuple[Array, ...] | None = None
    N_j_groups: tuple[Array, ...] | None = None
    # angles
    theta: Array | None = None
    phi: Array | None = None
    theta_groups: tuple[Array, ...] | None = None
    phi_groups: tuple[Array, ...] | None = None
    # other
    jump_rate: Array
    J_drive: Array

    class Config:
        arbitrary_types_allowed = True
