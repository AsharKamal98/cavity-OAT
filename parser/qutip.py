from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import qutip as qt
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Phase


class QutipMESolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mesolve benchmark."""

    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    phases: list[Phase]
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMESolverParameters":
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if not self.Ni:
            raise ValueError("Ni must contain at least one group size.")
        if len(self.omega_i) != len(self.Ni) - 1:
            raise ValueError("omega_i must contain the first G-1 group couplings.")
        if any(group_size < 0 for group_size in self.Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        return self


class QutipMCSolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mcsolve benchmark."""

    Gamma: float
    phases: list[Phase]
    Ni: tuple[int, ...]
    omega_i: tuple[float, ...]
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMCSolverParameters":
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if not self.Ni:
            raise ValueError("Ni must contain at least one group size.")
        if len(self.omega_i) != len(self.Ni) - 1:
            raise ValueError("omega_i must contain the first G-1 group couplings.")
        if any(group_size < 0 for group_size in self.Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        return self

@dataclass(frozen=True)
class QutipGroupedFixedNjModel:
    NJi: tuple[int, ...]
    omega_i: tuple[float, ...]
    Gamma: float
    shifted_jump_operator: bool
    omega0: float
    delta0: float
    phases: Sequence[Phase]
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
    t_step1_end: float
    t_step2_end: float
    t_final: float
