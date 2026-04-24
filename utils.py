from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class Phase:
    """One piecewise-constant stage of the protocol."""

    duration: float
    omega: float
    delta: float
    label: str = ""

# ----------------------------------------------------
# Helper functions
# ----------------------------------------------------

def omega_c(N_J: int, Gamma: float) -> float:
    """Critical drive for the polarized-to-mixed transition at delta = 0."""

    return 0.5 * N_J * Gamma

def default_three_phase_protocol(
    T1: float,
    T2: float,
    T3: float,
    delta0: float,
    Omega0: float,
) -> List[Phase]:
    """Three-phase protocol"""
    return [
        Phase(duration=T1, omega=Omega0, delta=0.0, label="phase1"),
        Phase(duration=T2, omega=Omega0, delta=delta0, label="phase2"),
        Phase(duration=T3, omega=0.0, delta=0.0, label="phase3"),
    ]