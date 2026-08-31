from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    StrictInt,
    field_validator,
    model_validator,
)
from scipy.sparse import csc_matrix

from parser.common import Array, PhaseProtocol
from parser.generator_moments import TrajectoryResolvedMomentSeries

SectorKey = int | tuple[int, ...]


class LinearSectorConstraints(BaseModel):
    """Exact integer constraints ``matrix @ sector = offset``."""

    matrix: tuple[tuple[StrictInt, ...], ...]
    offset: tuple[StrictInt, ...]

    model_config = ConfigDict(frozen=True)

    @model_validator(mode="after")
    def validate_constraints(self) -> "LinearSectorConstraints":
        if not self.matrix:
            raise ValueError("linear constraint matrix must contain at least one row.")
        row_size = len(self.matrix[0])
        if row_size == 0 or any(len(row) != row_size for row in self.matrix):
            raise ValueError(
                "linear constraint matrix rows must have one shared non-zero length."
            )
        if len(self.offset) != len(self.matrix):
            raise ValueError(
                "linear constraint offset must contain one value per matrix row."
            )
        return self


class InitialSectorSelectionSummary(BaseModel):
    """Counts and full-wavefunction weights for one sector selection."""

    retained_sector_count: int
    eligible_sector_count: int
    total_sector_count: int
    retained_percentage_of_eligible: float
    retained_percentage_of_full_space: float
    eligible_joint_probability: float
    retained_joint_probability: float

    model_config = ConfigDict(frozen=True)

    @model_validator(mode="after")
    def validate_summary(self) -> "InitialSectorSelectionSummary":
        if not (
            0 < self.retained_sector_count
            <= self.eligible_sector_count
            <= self.total_sector_count
        ):
            raise ValueError(
                "Initial-sector counts must satisfy 0 < retained <= eligible <= total."
            )
        if not 0.0 < self.eligible_joint_probability <= 1.0:
            raise ValueError("eligible_joint_probability must lie in (0, 1].")
        if not (
            0.0
            < self.retained_joint_probability
            <= self.eligible_joint_probability
        ):
            raise ValueError(
                "Initial-sector probabilities must satisfy "
                "0 < retained <= eligible <= 1."
            )
        return self


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
    dN: int | None = 0
    joint_probability: float | None = None
    linear_constraints: LinearSectorConstraints | None = None
    sector_distribution: str = "binomial"
    dt: float = 1e-3
    shifted_jump_operator: bool = True

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("linear_constraints", mode="before")
    @classmethod
    def convert_linear_constraints(cls, value: object) -> object:
        if value is None or isinstance(value, LinearSectorConstraints):
            return value
        if isinstance(value, (tuple, list)) and len(value) == 2:
            return {"matrix": value[0], "offset": value[1]}
        return value

    @model_validator(mode="after")
    def validate_inputs(self) -> "MCWFSolverParameters":
        if not self.Ni or sum(self.Ni) <= 0:
            raise ValueError("Ni must contain at least one atom group.")
        if any(group_size < 0 for group_size in self.Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        if len(self.Ni) != len(self.omega_i):
            raise ValueError("Ni and omega_i must contain the same number of groups.")
        if not all(np.isfinite(coupling) for coupling in self.omega_i):
            raise ValueError("omega_i must contain finite group couplings.")
        if not np.isfinite(self.Gamma) or self.Gamma <= 0.0:
            raise ValueError("Gamma must be finite and positive.")
        if (self.dN is None) == (self.joint_probability is None):
            raise ValueError("Exactly one of dN and joint_probability must be provided.")
        if self.dN is not None:
            total_atoms = sum(self.Ni)
            if self.dN < 0:
                raise ValueError("dN must be non-negative.")
            if total_atoms % 2 != 0:
                raise ValueError("dN selection requires an even total atom number.")
            if self.dN > total_atoms // 2:
                raise ValueError("dN cannot exceed half the total atom number.")
        if self.joint_probability is not None:
            if (
                not np.isfinite(self.joint_probability)
                or not 0.0 < self.joint_probability <= 1.0
            ):
                raise ValueError("joint_probability must lie in (0, 1].")
        if not np.isfinite(self.dt) or self.dt <= 0.0:
            raise ValueError("dt must be finite and positive.")
        if self.sector_distribution not in {"square", "binomial"}:
            raise ValueError("sector_distribution must be 'square' or 'binomial'.")
        if self.joint_probability is not None and self.sector_distribution != "binomial":
            raise ValueError(
                "joint_probability selection requires sector_distribution='binomial'."
            )
        if self.linear_constraints is not None:
            expected_columns = len(self.Ni)
            if any(
                len(row) != expected_columns
                for row in self.linear_constraints.matrix
            ):
                raise ValueError(
                    "Each linear constraint row must contain one coefficient per group."
                )
        return self


class MCWFTrajectoryDiagnostics(BaseModel):
    """MCWF-only diagnostics retained by a single-trajectory call."""

    jump_times: Array
    jump_count: int
    total_step_count: int
    non_precomputed_step_count: int

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("jump_times", mode="before")
    @classmethod
    def convert_jump_times(cls, value: object) -> Array:
        return np.asarray(value, dtype=float)

    @model_validator(mode="after")
    def validate_diagnostics(self) -> "MCWFTrajectoryDiagnostics":
        if self.jump_times.ndim != 1:
            raise ValueError("jump_times must be one-dimensional.")
        if self.jump_count != self.jump_times.size:
            raise ValueError("jump_count must equal len(jump_times).")
        if self.total_step_count < 0 or self.non_precomputed_step_count < 0:
            raise ValueError("MCWF step counts must be non-negative.")
        if self.non_precomputed_step_count > self.total_step_count:
            raise ValueError(
                "non_precomputed_step_count cannot exceed total_step_count."
            )
        return self


class TrajectoryResult(BaseModel):
    """Moment output and transient diagnostics for one MCWF trajectory."""

    moments: TrajectoryResolvedMomentSeries
    diagnostics: MCWFTrajectoryDiagnostics

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_single_trajectory(self) -> "TrajectoryResult":
        if self.moments.num_trajectories != 1:
            raise ValueError("TrajectoryResult.moments must contain one trajectory.")
        return self


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
    moments: TrajectoryResolvedMomentSeries
    seeds: list[tuple[int, ...]]
    metadata: TrajectoryEnsembleMetadata | dict[str, Any] | None = None

    def __post_init__(self) -> None:
        self.metadata = self._validate_metadata(self.metadata)
        if len(self.seeds) != self.moments.num_trajectories:
            raise ValueError(
                "TrajectoryEnsemble.seeds must contain one seed per trajectory."
            )
        if self.metadata is not None and not np.array_equal(
            np.asarray(self.metadata.t_eval, dtype=float), self.moments.t
        ):
            raise ValueError(
                "TrajectoryEnsemble metadata and moments must use the same t_eval grid."
            )
        if self.metadata is not None:
            group_count = self.moments.first_moments.shape[2]
            if len(self.metadata.Ni) != group_count:
                raise ValueError(
                    "TrajectoryEnsemble metadata.Ni must match the moment group axis."
                )
            if len(self.metadata.omega_i) != group_count:
                raise ValueError(
                    "TrajectoryEnsemble metadata.omega_i must contain one "
                    "completed coupling per physical group."
                )

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
