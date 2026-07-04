from __future__ import annotations

from dataclasses import dataclass

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class Phase:
    """One piecewise-constant stage of the protocol."""

    duration: float
    omega: float
    delta: float
    label: str = ""
