from dataclasses import dataclass
from typing import List

import numpy as np
import qutip as qt

from common.utils.phases import phase_boundary_times, phase_values_at_time
from parser.common import Phase, PhaseProtocol


def _omega_coeff(t, args):
    omega_t, _ = phase_values_at_time(t, args["integration_phases"])
    return omega_t


def _delta_coeff(t, args):
    _, delta_t = phase_values_at_time(t, args["integration_phases"])
    return delta_t


@dataclass(frozen=True)
class OmegaCoeffFromIntegrationPhases:
    """
    Pickle-safe omega(t) coefficient for time-dependent QuTiP operators.
    """

    integration_phases: tuple[Phase, ...]

    def __call__(self, t, args=None):
        if args is not None and "integration_phases" in args:
            return _omega_coeff(t, args)
        omega_t, _ = phase_values_at_time(t, self.integration_phases)
        return omega_t


def build_tlist_from_protocol(
    phase_protocol: PhaseProtocol,
    num_points: int,
) -> np.ndarray:
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")
    integration_boundaries = phase_boundary_times(
        phase_protocol.integration_phases
    )
    requested_times = np.linspace(0.0, phase_protocol.total_duration, num_points)
    return np.unique(np.concatenate((requested_times, integration_boundaries)))


def _solver_args(model) -> dict[str, object]:
    return {
        "integration_phases": model.phase_protocol.integration_phases,
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
