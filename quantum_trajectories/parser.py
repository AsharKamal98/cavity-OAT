from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

from common.parser import Array, AveragedResult, ObservableSeries, Phase
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
    N_e: csc_matrix
    JpJm: csc_matrix
    sector_key: Optional[SectorKey] = None
    Nj_groups: Optional[Tuple[int, ...]] = None
    omega_groups: Optional[Tuple[float, ...]] = None
    J_x_drive: Optional[csc_matrix] = None
    A_weighted: Optional[csc_matrix] = None
    AdagA_weighted: Optional[csc_matrix] = None
    N_e_groups: Optional[Tuple[csc_matrix, ...]] = None
    J_x_groups: Optional[Tuple[csc_matrix, ...]] = None
    J_y_groups: Optional[Tuple[csc_matrix, ...]] = None


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

# -----------------------------------------------------------------------------
# Ensambles
# -----------------------------------------------------------------------------

@dataclass
class TrajectoryEnsemble:
    trajectories: List[TrajectoryResult]
    seeds: List[Tuple[int, ...]]
