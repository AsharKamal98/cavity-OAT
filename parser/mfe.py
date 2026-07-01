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
        return len(self.omega_groups)

    @model_validator(mode="after")
    def validate_groups(self) -> "MFESolverParameters":
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.omega_groups:
            raise ValueError("omega_groups must contain at least one group.")
        if len(self.N_j_groups) != len(self.omega_groups):
            raise ValueError("N_j_groups and omega_groups must have matching lengths.")
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


class MFEObservableSeries(BaseModel):
    """Observable time series derived from the solved MFE amplitudes."""

    t: Array
    D_groups: tuple[Array, ...]
    E_groups: tuple[Array, ...]
    N_j_groups: tuple[Array, ...]
    theta_groups: tuple[Array, ...]
    phi_groups: tuple[Array, ...]
    x_groups: tuple[Array, ...]
    y_groups: tuple[Array, ...]
    z_groups: tuple[Array, ...]
    length_groups: tuple[Array, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)


class MFEResult(BaseModel):
    """Raw MFE solution and derived observables."""

    t: Array
    D_groups: tuple[Array, ...]
    E_groups: tuple[Array, ...]
    success: bool
    message: str
    parameters: MFESolverParameters

    model_config = ConfigDict(arbitrary_types_allowed=True)
