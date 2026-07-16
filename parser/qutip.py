from __future__ import annotations

from dataclasses import dataclass

import qutip as qt
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import PhaseProtocol


class QutipMESolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mesolve benchmark."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phase_protocol: PhaseProtocol
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMESolverParameters":
        if len(self.Ni) != len(self.omega_i):
            raise ValueError("Ni and omega_i must contain the same number of groups.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        return self


class QutipMCSolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mcsolve benchmark."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phase_protocol: PhaseProtocol
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMCSolverParameters":
        if len(self.Ni) != len(self.omega_i):
            raise ValueError("Ni and omega_i must contain the same number of groups.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        return self

@dataclass(frozen=True)
class QutipGroupedFixedNjModel:
    NJi: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    shifted_jump_operator: bool
    phase_protocol: PhaseProtocol
    Jp: qt.Qobj
    Jm: qt.Qobj
    Jx: qt.Qobj
    Jy: qt.Qobj
    Jz: qt.Qobj
    N_e: qt.Qobj
    Jp_groups: tuple[qt.Qobj, ...]
    Jm_groups: tuple[qt.Qobj, ...]
    Jx_groups: tuple[qt.Qobj, ...]
    Jy_groups: tuple[qt.Qobj, ...]
    Jz_groups: tuple[qt.Qobj, ...]
    N_e_groups: tuple[qt.Qobj, ...]
    J_drive: qt.Qobj
    A_weighted: qt.Qobj
    H: list
    c_ops: list
    psi0: qt.Qobj
