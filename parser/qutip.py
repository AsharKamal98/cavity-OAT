from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import qutip as qt
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Phase


class QutipMESolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mesolve benchmark."""

    N: int
    Gamma: float
    phases: list[Phase]
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMESolverParameters":
        if self.N <= 0:
            raise ValueError("N must be positive.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        return self


class QutipMCSolverParameters(BaseModel):
    """Validated inputs for the fixed-NJ QuTiP mcsolve benchmark."""

    N: int
    Gamma: float
    phases: list[Phase]
    Ni: tuple[int, int]
    omega_i: tuple[float, ...]
    shifted_jump_operator: bool = False

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_inputs(self) -> "QutipMCSolverParameters":
        if self.N <= 0:
            raise ValueError("N must be positive.")
        if self.Gamma <= 0.0:
            raise ValueError("Gamma must be positive.")
        if not self.phases:
            raise ValueError("phases must contain at least one phase.")
        if len(self.Ni) != 2:
            raise ValueError("Ni must contain exactly two group sizes.")
        if len(self.omega_i) != 1:
            raise ValueError("omega_i must contain the first group coupling only.")
        if any(group_size < 0 for group_size in self.Ni):
            raise ValueError("Ni must contain non-negative group sizes.")
        return self


@dataclass(frozen=True)
class QutipFixedNjModel:
    N: int
    NJ: int
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
    H: list
    c_ops: list
    psi0: qt.Qobj
    t_step1_end: float
    t_step2_end: float
    t_final: float


@dataclass(frozen=True)
class QutipTwoGroupFixedNjModel:
    N: int
    N1: int
    N2: int
    NJ1: int
    NJ2: int
    NJ: int
    Gamma: float
    shifted_jump_operator: bool
    omega0: float
    delta0: float
    phases: Sequence[Phase]
    omega_1: float
    omega_2: float
    Jp: qt.Qobj
    Jm: qt.Qobj
    Jx: qt.Qobj
    Jy: qt.Qobj
    Jz: qt.Qobj
    N_e: qt.Qobj
    Jp_groups: tuple[qt.Qobj, qt.Qobj]
    Jm_groups: tuple[qt.Qobj, qt.Qobj]
    Jx_groups: tuple[qt.Qobj, qt.Qobj]
    Jy_groups: tuple[qt.Qobj, qt.Qobj]
    Jz_groups: tuple[qt.Qobj, qt.Qobj]
    N_e_groups: tuple[qt.Qobj, qt.Qobj]
    J_drive: qt.Qobj
    A_weighted: qt.Qobj
    H: list
    c_ops: list
    psi0: qt.Qobj
    t_step1_end: float
    t_step2_end: float
    t_final: float
