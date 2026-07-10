from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from parser.common import Phase


def default_three_phase_protocol(
    T1: float,
    T2: float,
    T3: float,
    delta0: float,
    Omega0: float,
) -> List[Phase]:
    """Return the standard three-phase piecewise-constant protocol."""
    return [
        Phase(duration=T1, omega=Omega0, delta=0.0, label="phase1"),
        Phase(duration=T2, omega=Omega0, delta=delta0, label="phase2"),
        Phase(duration=T3, omega=0.0, delta=0.0, label="phase3"),
    ]


def phase_boundary_times(phases: Sequence[Phase]) -> np.ndarray:
    """
    Return the cumulative phase-end times for a piecewise-constant protocol.
    """
    if not phases:
        raise ValueError("Need at least one phase.")
    durations = np.asarray([phase.duration for phase in phases], dtype=float)
    return np.cumsum(durations)

def phase_values_at_time(t: float, phases: Sequence[Phase]) -> Tuple[float, float]:
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
