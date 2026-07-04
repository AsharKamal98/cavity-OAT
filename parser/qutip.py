from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import qutip as qt

from parser.common import Phase


@dataclass(frozen=True)
class QutipFixedNjModel:
    N: int
    N_J: int
    j: float
    Gamma: float
    shifted_jump_operator: bool
    unraveling_picture: str
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
    N_J1: int
    N_J2: int
    N_J: int
    j1: float
    j2: float
    Gamma: float
    shifted_jump_operator: bool
    unraveling_picture: str
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
