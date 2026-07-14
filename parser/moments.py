from __future__ import annotations

from typing import Any

import numpy as np
from pydantic import BaseModel, Field, root_validator

from parser.common import Array, Phase
from parser.j_moments import JMomentSeries
from parser.mfe_residuals import MFEResidualSeries
from common.utils.parameters import omega_G_from_weighted_average
from common.utils.phases import default_three_phase_protocol


class SimulationMetadata(BaseModel):
    """Shared physical model and standard protocol data for one simulation."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    Omega0: float
    delta0: float
    T1: float
    T2: float
    T3: float
    omega_groups: tuple[float, ...] = ()
    phases: list[Phase] = Field(default_factory=list)

    @root_validator(skip_on_failure=True)
    def complete_couplings_and_phases(cls, values: dict[str, Any]) -> dict[str, Any]:
        Ni = values["Ni"]
        omega_i = values["omega_i"]
        if any(N_g < 0 for N_g in Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        if len(Ni) != len(omega_i) + 1:
            raise ValueError("Ni must contain exactly one more element than omega_i.")
        if values["Gamma"] <= 0.0:
            raise ValueError("Gamma must be positive.")

        values["omega_groups"] = omega_i + (
            omega_G_from_weighted_average(omega_i, Ni),
        )
        values["phases"] = default_three_phase_protocol(
            T1=values["T1"],
            T2=values["T2"],
            T3=values["T3"],
            delta0=values["delta0"],
            Omega0=values["Omega0"],
        )
        return values

    class Config:
        arbitrary_types_allowed = True


class MomentSeries(BaseModel):
    """Container for moment series computed on a shared time grid."""

    t: Array
    metadata: SimulationMetadata | None = None
    J: JMomentSeries | None = None
    mfe_residuals: MFEResidualSeries | None = None
    S: Any | None = None

    @root_validator(pre=True)
    def build_t_eval_from_num_snapshots(cls, values: dict[str, Any]) -> dict[str, Any]:
        if values.get("t") is not None:
            return values

        num_snapshots = values.get("num_snapshots")
        metadata = values.get("metadata")
        if metadata is None or num_snapshots is None:
            return values

        if num_snapshots < 2:
            raise ValueError("num_snapshots must be at least 2.")
        if isinstance(metadata, dict):
            total_time = float(metadata["T1"] + metadata["T2"] + metadata["T3"])
        else:
            total_time = float(metadata.T1 + metadata.T2 + metadata.T3)
        values["t"] = np.linspace(0.0, total_time, num_snapshots, dtype=float)
        return values

    class Config:
        arbitrary_types_allowed = True
