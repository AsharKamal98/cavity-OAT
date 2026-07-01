from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from scipy.sparse import csc_matrix

from parser.common import Array, Phase
from parser.moments import MomentParameters

SectorKey = int | tuple[int, int]


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
