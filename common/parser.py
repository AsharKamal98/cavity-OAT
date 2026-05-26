from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


Array = np.ndarray


@dataclass(frozen=True)
class Phase:
    """One piecewise-constant stage of the protocol."""

    duration: float
    omega: float
    delta: float
    label: str = ""


@dataclass
class ObservableSeries:
    """Time series of observables on a common saved-time grid."""

    t: Array

    Jx: Array
    Jy: Array
    Jz: Array
    N_e: Array
    # Physical jump rate r(t) evaluated with the same collapse operator used by
    # the underlying solver or MCWF trajectory.
    jump_rate: Array
    N_j: Array
    N_active: Array

    theta: Array
    phi: Array

    sx: Array
    sy: Array
    sz: Array

    Jx_std: Optional[Array] = None
    Jy_std: Optional[Array] = None
    Jz_std: Optional[Array] = None
    N_e_std: Optional[Array] = None
    # Optional spread across trajectories for ensemble or mcsolve data.
    jump_rate_std: Optional[Array] = None
    N_j_std: Optional[Array] = None
    N_active_std: Optional[Array] = None


@dataclass
class AveragedResult:
    N: int
    Gamma: float
    ntraj: int
    observables: ObservableSeries
