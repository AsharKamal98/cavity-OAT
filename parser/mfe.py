from __future__ import annotations

from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array, Phase


class MFESolverParameters(BaseModel):
    """Parameters needed by the standalone MFE solver."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phases: list[Phase]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @property
    def group_count(self) -> int:
        return len(self.Ni)

    @model_validator(mode="after")
    def validate_inputs(self) -> "MFESolverParameters":
        if len(self.Ni) != len(self.omega_i):
            raise ValueError("Ni and omega_i must contain the same number of groups.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        return self


class MFEResult(BaseModel):
    """Raw MFE solution and derived observables."""

    t: Array
    D_groups: tuple[Array, ...]
    E_groups: tuple[Array, ...]
    success: bool
    message: str
    parameters: MFESolverParameters

    model_config = ConfigDict(arbitrary_types_allowed=True)
