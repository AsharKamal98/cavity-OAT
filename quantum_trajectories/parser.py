from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from common.parser import Array, AveragedResult, ObservableSeries, Phase
from pydantic import BaseModel, root_validator
from scipy.sparse import csc_matrix

SectorKey = int | tuple[int, int]


# ----------------------------------------------------
# Classes
# ----------------------------------------------------


@dataclass
class SectorOperators:
    Nj: int
    Jp: csc_matrix
    Jm: csc_matrix
    J_x: csc_matrix
    J_y: csc_matrix
    J_z: csc_matrix
    N_e: csc_matrix
    JpJm: csc_matrix
    sector_key: SectorKey | None = None
    Nj_groups: tuple[int, ...] | None = None
    omega_groups: tuple[float, ...] | None = None
    J_drive: csc_matrix | None = None
    A_weighted: csc_matrix | None = None
    AdagA_weighted: csc_matrix | None = None
    N_e_groups: tuple[csc_matrix, ...] | None = None
    J_x_groups: tuple[csc_matrix, ...] | None = None
    J_y_groups: tuple[csc_matrix, ...] | None = None
    J_z_groups: tuple[csc_matrix, ...] | None = None


@dataclass
class SectorWavefunction:
    """Wavefunction in one fixed-Nj block on the symmetric |n_e> basis."""

    Nj: int
    amplitudes: Array  # shape (Nj+1,), basis |n_e>, n_e = 0..Nj


@dataclass
class TrajectorySnapshot:
    time: float
    sector_blocks: dict[SectorKey, Array]
    norm: float
    phase_index: int


@dataclass
class TrajectoryResult:
    N: int
    Gamma: float
    phases: list[Phase]
    shifted_jump_operator: bool
    t_eval: Array
    sectors: list[SectorKey]
    sector_multiplicities: dict[SectorKey, int]
    final_sector_blocks: dict[SectorKey, Array]
    snapshots: list[TrajectorySnapshot]
    jump_times: list[float]
    jump_count: int
    sector_dimensions: dict[SectorKey, int]
    omega_1: float | None = None
    omega_2: float | None = None
    N1: int | None = None
    N2: int | None = None
    total_step_count: int = 0
    non_precomputed_step_count: int = 0



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
    # residual diagnostics
    mfe_residuals_groups: tuple[Array, ...] | None = None
    # other
    jump_rate: Array
    J_drive: Array

    class Config:
        arbitrary_types_allowed = True


class MomentParameters(BaseModel):
    """Shared simulation parameters needed by moment-level diagnostics."""

    Gamma: float
    phases: list[Phase]
    omega_groups: tuple[float, ...] | None = None
    N_groups: tuple[int, ...] | None = None

    class Config:
        arbitrary_types_allowed = True


class MomentSeries(BaseModel):
    """Container for moment series computed on a shared time grid."""

    t: Array
    parameters: MomentParameters | None = None
    J: JMomentSeries | None = None
    S: Any | None = None

    @root_validator(pre=True)
    def build_t_eval_from_num_snapshots(cls, values: dict[str, Any]) -> dict[str, Any]:
        if values.get("t") is not None:
            return values

        phases = values.get("phases")
        num_snapshots = values.get("num_snapshots")
        if phases is None or num_snapshots is None:
            return values

        if num_snapshots < 2:
            raise ValueError("num_snapshots must be at least 2.")
        total_time = float(sum(phase.duration for phase in phases))
        values["t"] = np.linspace(0.0, total_time, num_snapshots, dtype=float)
        return values

    class Config:
        arbitrary_types_allowed = True

# -----------------------------------------------------------------------------
# Ensambles
# -----------------------------------------------------------------------------

@dataclass
class TrajectoryEnsemble:
    trajectories: list[TrajectoryResult]
    seeds: list[tuple[int, ...]]
    parameters: MomentParameters | dict[str, Any] | None = None

    def __post_init__(self) -> None:
        self.parameters = self._validate_parameters(self.parameters)

    @staticmethod
    def _validate_parameters(
        parameters: MomentParameters | dict[str, Any] | None,
    ) -> MomentParameters | None:
        if parameters is None or isinstance(parameters, MomentParameters):
            return parameters

        values = dict(parameters)
        if "phases" in values:
            values["phases"] = list(values["phases"])
        if values.get("omega_groups") is not None:
            values["omega_groups"] = tuple(float(value) for value in values["omega_groups"])
        if values.get("N_groups") is not None:
            values["N_groups"] = tuple(int(value) for value in values["N_groups"])

        return MomentParameters(**values)
