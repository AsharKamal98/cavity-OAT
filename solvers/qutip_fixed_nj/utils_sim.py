from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import qutip as qt

from common.utils.phases import phase_boundary_times, phase_values_at_time
from parser.qutip import QutipFixedNjModel


def _omega_coeff(t, args):
    omega_t, _ = phase_values_at_time(t, args["phases"])
    return omega_t


def _delta_coeff(t, args):
    _, delta_t = phase_values_at_time(t, args["phases"])
    return delta_t


@dataclass(frozen=True)
class OmegaCoeffFromPhases:
    """
    Pickle-safe omega(t) coefficient for time-dependent QuTiP operators.
    """

    phases: list

    def __call__(self, t, args=None):
        if args is not None and "phases" in args:
            return _omega_coeff(t, args)
        omega_t, _ = phase_values_at_time(t, self.phases)
        return omega_t


def build_tlist_from_phases(phases: list, num_points: int) -> np.ndarray:
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")
    t_final = float(phase_boundary_times(phases)[-1])
    return np.linspace(0.0, t_final, num_points)


def _solver_args(model: QutipFixedNjModel) -> Dict[str, float]:
    return {
        "phases": model.phases,
    }


def _observable_e_ops(model) -> List[qt.Qobj]:
    """
    Return the standard observable list for QuTiP benchmarks.
    """
    e_ops = [model.Jx, model.Jy, model.Jz, model.N_e]
    if hasattr(model, "Jx_groups") and hasattr(model, "N_e_groups"):
        e_ops.extend(model.Jx_groups)
        e_ops.extend(model.Jy_groups)
        e_ops.extend(model.Jz_groups)
        e_ops.extend(model.N_e_groups)
    return e_ops
