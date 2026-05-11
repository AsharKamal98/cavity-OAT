from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

from common.parser import Array, AveragedResult, ObservableSeries, Phase
from scipy.sparse import csc_matrix

# ----------------------------------------------------
# Classes
# ----------------------------------------------------


@dataclass
class SectorOperators:
    Nj: int
    J_plus: csc_matrix
    J_minus: csc_matrix
    J_x: csc_matrix
    J_y: csc_matrix
    N_e: csc_matrix
    JpJm: csc_matrix


@dataclass
class SectorWavefunction:
    """Wavefunction in one fixed-Nj block on the symmetric |n_e> basis."""

    Nj: int
    amplitudes: Array  # shape (Nj+1,), basis |n_e>, n_e = 0..Nj


@dataclass
class TrajectorySnapshot:
    time: float
    sector_blocks: Dict[int, Array]
    norm: float
    phase_index: int


@dataclass
class TrajectoryResult:
    N: int
    gamma: float
    phases: List[Phase]
    shifted_jump_operator: bool
    sectors: List[int]
    sector_multiplicities: Dict[int, int]
    final_sector_blocks: Dict[int, Array]
    snapshots: List[TrajectorySnapshot]
    jump_times: List[float]
    jump_count: int
    sector_dimensions: Dict[int, int]

# -----------------------------------------------------------------------------
# Ensambles
# -----------------------------------------------------------------------------

@dataclass
class TrajectoryEnsemble:
    trajectories: List[TrajectoryResult]
    seeds: List[int]
