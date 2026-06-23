from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
from common.parser import Array, AveragedResult, ObservableSeries, Phase
from pydantic import BaseModel, root_validator
from scipy.sparse import csc_matrix

SectorKey = Union[int, Tuple[int, int]]


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
    sector_key: Optional[SectorKey] = None
    Nj_groups: Optional[Tuple[int, ...]] = None
    omega_groups: Optional[Tuple[float, ...]] = None
    J_drive: Optional[csc_matrix] = None
    A_weighted: Optional[csc_matrix] = None
    AdagA_weighted: Optional[csc_matrix] = None
    N_e_groups: Optional[Tuple[csc_matrix, ...]] = None
    J_x_groups: Optional[Tuple[csc_matrix, ...]] = None
    J_y_groups: Optional[Tuple[csc_matrix, ...]] = None
    J_z_groups: Optional[Tuple[csc_matrix, ...]] = None


@dataclass
class SectorWavefunction:
    """Wavefunction in one fixed-Nj block on the symmetric |n_e> basis."""

    Nj: int
    amplitudes: Array  # shape (Nj+1,), basis |n_e>, n_e = 0..Nj


@dataclass
class TrajectorySnapshot:
    time: float
    sector_blocks: Dict[SectorKey, Array]
    norm: float
    phase_index: int


@dataclass
class TrajectoryResult:
    N: int
    Gamma: float
    phases: List[Phase]
    shifted_jump_operator: bool
    t_eval: Array
    sectors: List[SectorKey]
    sector_multiplicities: Dict[SectorKey, int]
    final_sector_blocks: Dict[SectorKey, Array]
    snapshots: List[TrajectorySnapshot]
    jump_times: List[float]
    jump_count: int
    sector_dimensions: Dict[SectorKey, int]
    omega_1: Optional[float] = None
    omega_2: Optional[float] = None
    N1: Optional[int] = None
    N2: Optional[int] = None
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
    x_groups: Tuple[float, ...] | None = None
    y_groups: Tuple[float, ...] | None = None
    z_groups: Tuple[float, ...] | None = None
    N_e_groups: Tuple[float, ...] | None = None
    N_j_groups: Tuple[float, ...] | None = None


class JMomentSeries(BaseModel):
    """Per-timestep first-order J-sphere moments for one (single or averaged) trajectory."""

    # Legacy note: these fields were previously named Jx/Jy/Jz, Jx_groups/
    # Jy_groups/Jz_groups, J_len, and sx/sy/sz.
    t: Array
    phase_index: Array
    x: Array
    y: Array
    z: Array
    N_e: Array
    N_j: Array
    jump_rate: Array
    J_drive: Array
    x_groups: Tuple[Array, ...] | None = None
    y_groups: Tuple[Array, ...] | None = None
    z_groups: Tuple[Array, ...] | None = None
    N_e_groups: Tuple[Array, ...] | None = None
    N_j_groups: Tuple[Array, ...] | None = None
    length: Array | None = None
    nx: Array | None = None
    ny: Array | None = None
    nz: Array | None = None
    length_groups: Tuple[Array, ...] | None = None
    nx_groups: Tuple[Array, ...] | None = None
    ny_groups: Tuple[Array, ...] | None = None
    nz_groups: Tuple[Array, ...] | None = None
    theta: Array | None = None
    phi: Array | None = None
    theta_groups: Tuple[Array, ...] | None = None
    phi_groups: Tuple[Array, ...] | None = None

    class Config:
        arbitrary_types_allowed = True


class MomentSeries(BaseModel):
    """Container for moment series computed on a shared time grid."""

    t: Array
    J: JMomentSeries | None = None
    S: Any | None = None

    @root_validator(pre=True)
    def build_t_eval_from_num_snapshots(cls, values: Dict[str, Any]) -> Dict[str, Any]:
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
    trajectories: List[TrajectoryResult]
    seeds: List[Tuple[int, ...]]
