from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from pydantic import BaseModel, ConfigDict, model_validator
from scipy.sparse import csc_matrix

from parser.common import Array, PhaseProtocol

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
    integration_phase_index: int


class MCWFSolverParameters(BaseModel):
    """Validated inputs for the MCWF ensemble solver."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phase_protocol: PhaseProtocol
    dN: int = 0
    sector_distribution: str = "binomial"
    dt: float = 1e-3
    shifted_jump_operator: bool = True

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "MCWFSolverParameters":
        if len(self.Ni) != len(self.omega_i):
            raise ValueError("Ni and omega_i must contain the same number of groups.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if self.dN < 0:
            raise ValueError("dN must be non-negative.")
        if self.dt <= 0.0:
            raise ValueError("dt must be positive.")
        if self.sector_distribution not in {"square", "binomial"}:
            raise ValueError("sector_distribution must be 'square' or 'binomial'.")
        return self


@dataclass
class TrajectoryResult:
    final_sector_blocks: dict[SectorKey, Array]
    snapshots: list[TrajectorySnapshot]
    jump_times: list[float]
    jump_count: int
    total_step_count: int = 0
    non_precomputed_step_count: int = 0


@dataclass
class TrajectoryEnsembleMetadata:
    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phase_protocol: PhaseProtocol
    shifted_jump_operator: bool
    t_eval: Array
    sectors: list[SectorKey]
    sector_multiplicities: dict[SectorKey, int]
    sector_dimensions: dict[SectorKey, int]


@dataclass
class TrajectoryEnsemble:
    trajectories: list[TrajectoryResult]
    seeds: list[tuple[int, ...]]
    metadata: TrajectoryEnsembleMetadata | dict[str, Any] | None = None

    def __post_init__(self) -> None:
        self.metadata = self._validate_metadata(self.metadata)

    @staticmethod
    def _validate_metadata(
        metadata: TrajectoryEnsembleMetadata | dict[str, Any] | None,
    ) -> TrajectoryEnsembleMetadata | None:
        if metadata is None or isinstance(metadata, TrajectoryEnsembleMetadata):
            return metadata

        values = dict(metadata)
        values["Ni"] = tuple(int(value) for value in values["Ni"])
        values["omega_i"] = tuple(float(value) for value in values["omega_i"])
        if isinstance(values["phase_protocol"], dict):
            values["phase_protocol"] = PhaseProtocol.model_validate(
                values["phase_protocol"]
            )
        values["sectors"] = list(values["sectors"])

        return TrajectoryEnsembleMetadata(**values)
