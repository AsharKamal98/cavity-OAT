from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator


Array = np.ndarray


class Phase(BaseModel):
    """One piecewise-constant integration segment."""

    duration: float
    omega: float
    delta: float
    label: str = ""

    model_config = ConfigDict(frozen=True)

    @model_validator(mode="after")
    def validate_duration(self) -> "Phase":
        if self.duration < 0.0:
            raise ValueError("Phase duration must be non-negative.")
        return self


class FamilyPhase(BaseModel):
    """One physical protocol phase containing constant integration segments."""

    duration: float
    omega: float
    delta: float
    ramp_duration: float = 0.0
    num_ramp_segments: int = 0
    label: str = ""
    segments: tuple[Phase, ...] = ()

    model_config = ConfigDict(frozen=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "FamilyPhase":
        if self.duration < 0.0:
            raise ValueError("FamilyPhase duration must be non-negative.")
        if not 0.0 <= self.ramp_duration <= self.duration:
            raise ValueError("FamilyPhase ramp_duration must lie within its duration.")
        if self.num_ramp_segments < 0:
            raise ValueError("FamilyPhase num_ramp_segments must be non-negative.")
        if (self.ramp_duration == 0.0) != (self.num_ramp_segments == 0):
            raise ValueError(
                "FamilyPhase ramp_duration and num_ramp_segments must either both "
                "be zero or both be positive."
            )
        return self


class PhaseProtocol(BaseModel):
    """Ordered physical phases and their flattened integration segments."""

    family_phases: tuple[FamilyPhase, ...]

    @staticmethod
    def _build_segments(
        start_values: tuple[float, float],
        target_values: tuple[float, float],
        duration: float,
        ramp_duration: float,
        num_ramp_segments: int,
        label: str,
    ) -> tuple[Phase, ...]:
        start_omega, start_delta = start_values
        target_omega, target_delta = target_values
        if num_ramp_segments == 0:
            return (
                Phase(
                    duration=duration,
                    omega=target_omega,
                    delta=target_delta,
                    label=label,
                ),
            )

        segment_duration = ramp_duration / num_ramp_segments
        ramp_segments = tuple(
            Phase(
                duration=segment_duration,
                omega=start_omega + fraction * (target_omega - start_omega),
                delta=start_delta + fraction * (target_delta - start_delta),
                label=label,
            )
            for fraction in (
                index / num_ramp_segments
                for index in range(1, num_ramp_segments + 1)
            )
        )
        hold_duration = duration - ramp_duration
        if hold_duration == 0.0:
            return ramp_segments
        return ramp_segments + (
            Phase(
                duration=hold_duration,
                omega=target_omega,
                delta=target_delta,
                label=label,
            ),
        )

    @model_validator(mode="after")
    def build_family_segments(self) -> "PhaseProtocol":
        if not self.family_phases:
            raise ValueError("PhaseProtocol requires at least one FamilyPhase.")

        start_values = (0.0, 0.0)
        completed_family_phases = []
        for family_phase in self.family_phases:
            target_values = (family_phase.omega, family_phase.delta)
            segments = self._build_segments(
                start_values,
                target_values,
                family_phase.duration,
                family_phase.ramp_duration,
                family_phase.num_ramp_segments,
                family_phase.label,
            )
            completed_family_phases.append(
                family_phase.model_copy(update={"segments": segments})
            )
            start_values = target_values

        self.family_phases = tuple(completed_family_phases)
        return self

    @property
    def integration_phases(self) -> tuple[Phase, ...]:
        return tuple(
            integration_phase
            for family_phase in self.family_phases
            for integration_phase in family_phase.segments
        )

    @property
    def integration_to_family_index(self) -> tuple[int, ...]:
        return tuple(
            family_index
            for family_index, family_phase in enumerate(self.family_phases)
            for _ in family_phase.segments
        )

    @property
    def total_duration(self) -> float:
        return float(
            sum(family_phase.duration for family_phase in self.family_phases)
        )
