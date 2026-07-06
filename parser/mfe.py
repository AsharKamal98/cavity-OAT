from __future__ import annotations

from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array, Phase


class MFESolverParameters(BaseModel):
    """Parameters needed by the standalone MFE solver."""

    Gamma: float
    phases: list[Phase]
    omega_i: tuple[float, ...]
    Ni: tuple[float, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @property
    def group_count(self) -> int:
        return len(self.Ni)

    @model_validator(mode="after")
    def validate_groups(self) -> "MFESolverParameters":
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if len(self.Ni) == 0:
            raise ValueError("Ni must contain at least one group.")
        if len(self.omega_i) not in (len(self.Ni) - 1, len(self.Ni)):
            raise ValueError(
                "omega_i must contain either one coupling per group or all but the final group coupling."
            )
        if any(N_g < 0.0 for N_g in self.Ni):
            raise ValueError("Ni must be non-negative.")
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
