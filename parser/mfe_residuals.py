from __future__ import annotations

from pydantic import BaseModel

from parser.common import Array


class MFEResidualSeries(BaseModel):
    """Per-timestep mean-field-equation residual diagnostics."""

    t: Array
    integration_phase_index: Array
    residuals_groups: tuple[Array, ...]

    class Config:
        arbitrary_types_allowed = True
