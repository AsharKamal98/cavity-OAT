from __future__ import annotations

import numpy as np
from pydantic import (
    BaseModel,
    ConfigDict,
    ValidationInfo,
    field_validator,
    model_validator,
)

from parser.common import Array, PhaseProtocol
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.j_modes import JModeSeries
from parser.j_moments import JMomentSeries
from parser.mfe_residuals import MFEResidualSeries
from parser.s_moments import SMomentSeries
from parser.squeezing import SqueezingResult
from common.utils.parameters import omega_G_from_weighted_average


class SimulationMetadata(BaseModel):
    """Shared physical model and supplied protocol data for one simulation."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phase_protocol: PhaseProtocol
    omega_groups: tuple[float, ...] = ()

    @model_validator(mode="after")
    def complete_couplings(self) -> "SimulationMetadata":
        if any(N_g < 0 for N_g in self.Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        if len(self.Ni) != len(self.omega_i) + 1:
            raise ValueError("Ni must contain exactly one more element than omega_i.")
        if not all(np.isfinite(coupling) for coupling in self.omega_i):
            raise ValueError("omega_i must contain finite couplings.")
        if not np.isfinite(self.Gamma) or self.Gamma <= 0.0:
            raise ValueError("Gamma must be finite and positive.")

        self.omega_groups = self.omega_i + (
            omega_G_from_weighted_average(self.omega_i, self.Ni),
        )
        return self

    model_config = ConfigDict(arbitrary_types_allowed=True)


class MomentSeries(BaseModel):
    """Container for moment series computed on a shared time grid."""

    t: Array
    metadata: SimulationMetadata | None = None
    generator_moments: TrajectoryResolvedMomentSeries | None = None
    mean_generator_moments: TrajectoryResolvedMomentSeries | None = None
    J: JMomentSeries | None = None
    J_modes: JModeSeries | None = None
    mfe_residuals: MFEResidualSeries | None = None
    squeezing: SqueezingResult | None = None
    S: SMomentSeries | None = None

    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=True,
    )

    @model_validator(mode="before")
    @classmethod
    def build_time_grid(cls, values: object) -> object:
        if not isinstance(values, dict):
            return values
        if values.get("t") is not None:
            return values

        generator_moments = values.get("generator_moments")
        if generator_moments is not None:
            if not isinstance(generator_moments, TrajectoryResolvedMomentSeries):
                generator_moments = TrajectoryResolvedMomentSeries.model_validate(
                    generator_moments
                )
                values["generator_moments"] = generator_moments
            values["t"] = generator_moments.t.copy()
            return values

        mean_generator_moments = values.get("mean_generator_moments")
        if mean_generator_moments is not None:
            if not isinstance(
                mean_generator_moments,
                TrajectoryResolvedMomentSeries,
            ):
                mean_generator_moments = (
                    TrajectoryResolvedMomentSeries.model_validate(
                        mean_generator_moments
                    )
                )
                values["mean_generator_moments"] = mean_generator_moments
            values["t"] = mean_generator_moments.t.copy()
            return values

        num_snapshots = values.get("num_snapshots")
        metadata = values.get("metadata")
        if metadata is None or num_snapshots is None:
            return values

        if num_snapshots < 2:
            raise ValueError("num_snapshots must be at least 2.")
        if isinstance(metadata, dict):
            metadata = SimulationMetadata.model_validate(metadata)
            values["metadata"] = metadata
        total_time = metadata.phase_protocol.total_duration
        values["t"] = np.linspace(0.0, total_time, num_snapshots, dtype=float)
        return values

    @staticmethod
    def _validate_generator_context(
        *,
        t: Array,
        metadata: SimulationMetadata | None,
        generator_moments: TrajectoryResolvedMomentSeries | None,
        field_name: str = "generator_moments",
    ) -> None:
        if generator_moments is None:
            return

        generator_time = generator_moments.t
        if t.shape != generator_time.shape or not np.allclose(
            t,
            generator_time,
            atol=1e-12,
            rtol=0.0,
        ):
            raise ValueError(
                f"MomentSeries.t and {field_name}.t must use the same grid."
            )

        if metadata is not None:
            metadata_group_count = len(metadata.Ni)
            moment_group_count = generator_moments.first_moments.shape[2]
            if moment_group_count != metadata_group_count:
                raise ValueError(
                    f"{field_name} must contain one group axis entry for "
                    "each group in metadata.Ni."
                )

    @classmethod
    def _validate_mean_generator_context(
        cls,
        *,
        t: Array,
        metadata: SimulationMetadata | None,
        generator_moments: TrajectoryResolvedMomentSeries | None,
        mean_generator_moments: TrajectoryResolvedMomentSeries | None,
    ) -> None:
        if mean_generator_moments is None:
            return

        cls._validate_generator_context(
            t=t,
            metadata=metadata,
            generator_moments=mean_generator_moments,
            field_name="mean_generator_moments",
        )
        if mean_generator_moments.num_trajectories != 1:
            raise ValueError(
                "mean_generator_moments must retain a singleton trajectory axis."
            )

        if generator_moments is not None and (
            mean_generator_moments.first_moments.shape[1:]
            != generator_moments.first_moments.shape[1:]
            or mean_generator_moments.second_moments.shape[1:]
            != generator_moments.second_moments.shape[1:]
        ):
            raise ValueError(
                "generator_moments and mean_generator_moments must have "
                "matching time, group, and packed-moment dimensions."
            )

    @staticmethod
    def _validate_squeezing_context(
        *,
        t: Array,
        metadata: SimulationMetadata | None,
        squeezing: SqueezingResult | None,
    ) -> None:
        if squeezing is None:
            return

        for group in squeezing.group_resolved:
            if group.t.shape != t.shape or not np.allclose(
                group.t,
                t,
                atol=1e-12,
                rtol=0.0,
            ):
                raise ValueError(
                    "MomentSeries.t and squeezing results must use the same grid."
                )

        if metadata is not None:
            if len(squeezing.group_resolved) != len(metadata.Ni):
                raise ValueError(
                    "squeezing must contain one group-resolved result for "
                    "each group in metadata.Ni."
                )
            if tuple(
                group.N_atoms for group in squeezing.group_resolved
            ) != tuple(metadata.Ni):
                raise ValueError(
                    "Group-resolved squeezing atom counts must match metadata.Ni."
                )

    @field_validator("t", mode="before")
    @classmethod
    def validate_time_grid(cls, value: object, info: ValidationInfo) -> Array:
        t = np.asarray(value, dtype=float)
        if t.ndim != 1:
            raise ValueError("MomentSeries.t must be one-dimensional.")
        if t.size == 0:
            raise ValueError("MomentSeries.t must contain at least one time point.")
        cls._validate_generator_context(
            t=t,
            metadata=info.data.get("metadata"),
            generator_moments=info.data.get("generator_moments"),
        )
        cls._validate_mean_generator_context(
            t=t,
            metadata=info.data.get("metadata"),
            generator_moments=info.data.get("generator_moments"),
            mean_generator_moments=info.data.get("mean_generator_moments"),
        )
        return t

    @field_validator("metadata")
    @classmethod
    def validate_metadata_context(
        cls,
        metadata: SimulationMetadata | None,
        info: ValidationInfo,
    ) -> SimulationMetadata | None:
        t = info.data.get("t")
        if t is not None:
            cls._validate_generator_context(
                t=t,
                metadata=metadata,
                generator_moments=info.data.get("generator_moments"),
            )
            cls._validate_mean_generator_context(
                t=t,
                metadata=metadata,
                generator_moments=info.data.get("generator_moments"),
                mean_generator_moments=info.data.get(
                    "mean_generator_moments"
                ),
            )
        return metadata

    @field_validator("generator_moments")
    @classmethod
    def validate_generator_moment_context(
        cls,
        generator_moments: TrajectoryResolvedMomentSeries | None,
        info: ValidationInfo,
    ) -> TrajectoryResolvedMomentSeries | None:
        t = info.data.get("t")
        if t is not None:
            cls._validate_generator_context(
                t=t,
                metadata=info.data.get("metadata"),
                generator_moments=generator_moments,
            )
            cls._validate_mean_generator_context(
                t=t,
                metadata=info.data.get("metadata"),
                generator_moments=generator_moments,
                mean_generator_moments=info.data.get(
                    "mean_generator_moments"
                ),
            )
        return generator_moments

    @field_validator("mean_generator_moments")
    @classmethod
    def validate_mean_generator_moment_context(
        cls,
        mean_generator_moments: TrajectoryResolvedMomentSeries | None,
        info: ValidationInfo,
    ) -> TrajectoryResolvedMomentSeries | None:
        t = info.data.get("t")
        if t is not None:
            cls._validate_mean_generator_context(
                t=t,
                metadata=info.data.get("metadata"),
                generator_moments=info.data.get("generator_moments"),
                mean_generator_moments=mean_generator_moments,
            )
        return mean_generator_moments

    @model_validator(mode="after")
    def validate_analysis_context(self) -> "MomentSeries":
        self._validate_generator_context(
            t=self.t,
            metadata=self.metadata,
            generator_moments=self.generator_moments,
        )
        self._validate_mean_generator_context(
            t=self.t,
            metadata=self.metadata,
            generator_moments=self.generator_moments,
            mean_generator_moments=self.mean_generator_moments,
        )
        self._validate_squeezing_context(
            t=self.t,
            metadata=self.metadata,
            squeezing=self.squeezing,
        )
        return self
