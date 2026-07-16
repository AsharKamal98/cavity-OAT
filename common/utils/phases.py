from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from parser.common import FamilyPhase, Phase, PhaseProtocol


def default_three_phase_protocol(
    durations: tuple[float, float, float],
    delta0: float,
    Omega0: float,
    ramp_durations: tuple[float, float, float] = (0.0, 0.0, 0.0),
    ramp_segment_counts: tuple[int, int, int] = (0, 0, 0),
) -> PhaseProtocol:
    """Return the standard three-family-phase protocol."""
    if (
        len(durations) != 3
        or len(ramp_durations) != 3
        or len(ramp_segment_counts) != 3
    ):
        raise ValueError(
            "durations, ramp_durations, and ramp_segment_counts must contain "
            "one value per family phase."
        )
    T1, T2, T3 = (float(duration) for duration in durations)
    return PhaseProtocol(
        family_phases=(
            FamilyPhase(
                duration=T1,
                omega=Omega0,
                delta=0.0,
                ramp_duration=ramp_durations[0],
                num_ramp_segments=ramp_segment_counts[0],
                label="phase1",
            ),
            FamilyPhase(
                duration=T2,
                omega=Omega0,
                delta=delta0,
                ramp_duration=ramp_durations[1],
                num_ramp_segments=ramp_segment_counts[1],
                label="phase2",
            ),
            FamilyPhase(
                duration=T3,
                omega=0.0,
                delta=0.0,
                ramp_duration=ramp_durations[2],
                num_ramp_segments=ramp_segment_counts[2],
                label="phase3",
            ),
        )
    )


def phase_boundary_times(phases: Sequence[Phase | FamilyPhase]) -> np.ndarray:
    """
    Return the cumulative phase-end times for a piecewise-constant protocol.
    """
    if not phases:
        raise ValueError("Need at least one phase.")
    durations = np.asarray([phase.duration for phase in phases], dtype=float)
    return np.cumsum(durations)

def phase_values_at_time(t: float, phases: Sequence[Phase]) -> tuple[float, float]:
    """
    Return phase-local (Omega, delta) values for a piecewise-constant protocol.
    """
    boundaries = phase_boundary_times(phases)
    t_value = float(t)
    total_time = float(boundaries[-1])
    t_value = min(max(t_value, 0.0), total_time)

    for index, boundary in enumerate(boundaries):
        if t_value <= boundary or index == len(phases) - 1:
            phase = phases[index]
            return float(phase.omega), float(phase.delta)

    phase = phases[-1]
    return float(phase.omega), float(phase.delta)


def integration_phase_indices_at_times(
    times: Sequence[float] | np.ndarray,
    phase_protocol: PhaseProtocol,
) -> np.ndarray:
    """Return the integration-Phase index at each supplied time."""
    boundaries = phase_boundary_times(phase_protocol.integration_phases)
    clipped_times = np.clip(np.asarray(times, dtype=float), 0.0, boundaries[-1])
    return np.searchsorted(boundaries, clipped_times, side="left").astype(int)
