from __future__ import annotations

from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array, Phase


class MFESolverParameters(BaseModel):
    """Parameters needed by the standalone MFE solver."""

    Gamma: float
    phases: list[Phase]
    omega_groups: tuple[float, ...]
    N_j_groups: tuple[float, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @property
    def group_count(self) -> int:
        return len(self.N_j_groups)

    @model_validator(mode="after")
    def validate_groups(self) -> "MFESolverParameters":
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.omega_groups:
            raise ValueError("omega_groups must contain at least one group.")
        if len(self.omega_groups) not in (len(self.N_j_groups) - 1, len(self.N_j_groups)):
            raise ValueError(
                "omega_groups must contain either one coupling per group or all but the final group coupling."
            )
        if any(N_j < 0.0 for N_j in self.N_j_groups):
            raise ValueError("N_j_groups must be non-negative.")
        return self


class MFEInitialState(BaseModel):
    """Initial group-resolved J-sphere angles."""

    theta_groups: tuple[float, ...]
    phi_groups: tuple[float, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_groups(self) -> "MFEInitialState":
        if len(self.theta_groups) != len(self.phi_groups):
            raise ValueError("theta_groups and phi_groups must have matching lengths.")
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
