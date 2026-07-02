from __future__ import annotations

from typing import Any

import numpy as np
from pydantic import BaseModel, root_validator

from parser.common import Array, Phase
from parser.j_moments import JMomentSeries
from parser.mfe_residuals import MFEResidualSeries


class MomentParameters(BaseModel):
    """Shared simulation parameters needed by moment-level diagnostics."""

    Gamma: float
    phases: list[Phase]
    omega_groups: tuple[float, ...] | None = None
    N_groups: tuple[int, ...] | None = None

    class Config:
        arbitrary_types_allowed = True


class MomentSeries(BaseModel):
    """Container for moment series computed on a shared time grid."""

    t: Array
    parameters: MomentParameters | None = None
    J: JMomentSeries | None = None
    mfe_residuals: MFEResidualSeries | None = None
    S: Any | None = None

    @root_validator(pre=True)
    def build_t_eval_from_num_snapshots(cls, values: dict[str, Any]) -> dict[str, Any]:
        if values.get("t") is not None:
            return values

        phases = values.get("phases")
        num_snapshots = values.get("num_snapshots")
        if phases is None or num_snapshots is None:
            return values

        if num_snapshots < 2:
            raise ValueError("num_snapshots must be at least 2.")
        total_time = float(sum(phase.duration for phase in phases))
        values["t"] = np.linspace(0.0, total_time, num_snapshots, dtype=float)
        return values

    class Config:
        arbitrary_types_allowed = True
