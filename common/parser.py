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
    t: Array

    Jx: Array
    Jy: Array
    Jz: Array
    N_e: Array
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
    N_j_std: Optional[Array] = None
    N_active_std: Optional[Array] = None


@dataclass
class AveragedResult:
    N: int
    gamma: float
    ntraj: int
    observables: ObservableSeries
