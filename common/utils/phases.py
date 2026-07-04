from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from parser.common import Phase


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
    if t_value < 0.0 or t_value > total_time:
        raise ValueError(f"t must lie in [0, {total_time}], got {t_value}.")

    for index, boundary in enumerate(boundaries):
        if t_value <= boundary or index == len(phases) - 1:
            phase = phases[index]
            return float(phase.omega), float(phase.delta)

    phase = phases[-1]
    return float(phase.omega), float(phase.delta)
