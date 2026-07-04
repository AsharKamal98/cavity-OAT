from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

import numpy as np
import qutip as qt

from parser.qutip import QutipFixedNjModel


def _omega_coeff(t, args):
    if t < args["t_step1_end"]:
        return args["Omega0"]
    if t < args["t_step2_end"]:
        return args["Omega0"]
    return 0.0


def _delta_coeff(t, args):
    if t < args["t_step1_end"]:
        return 0.0
    if t < args["t_step2_end"]:
        return args["delta0"]
    return 0.0


@dataclass(frozen=True)
class OmegaCoeffFromPhases:
    """
    Pickle-safe omega(t) coefficient for time-dependent QuTiP operators.
    """

    omega0: float
    t_step1_end: float
    t_step2_end: float

    def __call__(self, t, args=None):
        if args is not None and "t_step1_end" in args:
            return _omega_coeff(t, args)
        if t < self.t_step1_end:
            return self.omega0
        if t < self.t_step2_end:
            return self.omega0
        return 0.0


def build_tlist_from_phases(phases: Sequence, num_points: int) -> np.ndarray:
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")
    t_final = float(sum(p.duration for p in phases))
    return np.linspace(0.0, t_final, num_points)


def _solver_args(model: QutipFixedNjModel) -> Dict[str, float]:
    return {
        "Omega0": model.omega0,
        "delta0": model.delta0,
        "t_step1_end": model.t_step1_end,
        "t_step2_end": model.t_step2_end,
    }


def _observable_e_ops(model) -> List[qt.Qobj]:
    """
    Return the standard observable list for QuTiP benchmarks.
    """
    e_ops = [model.Jx, model.Jy, model.Jz, model.N_e]
    if hasattr(model, "Jx_groups") and hasattr(model, "N_e_groups"):
        e_ops.extend(
            [
                model.Jx_groups[0],
                model.Jx_groups[1],
                model.Jy_groups[0],
                model.Jy_groups[1],
                model.Jz_groups[0],
                model.Jz_groups[1],
                model.N_e_groups[0],
                model.N_e_groups[1],
            ]
        )
    return e_ops
