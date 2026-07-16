from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from common.utils.phases import phase_boundary_times
from parser.common import Array, Phase
from parser.mfe import (
    MFEResult,
    MFESolverParameters,
)
from solvers.mfe.utils import amplitudes_from_initial_state


def mfe_rhs(
    t: float,
    y: Array,
    parameters: MFESolverParameters,
    integration_phase: Phase,
) -> Array:
    """
    Right-hand side of the group-resolved mean-field equations.
    """
    _ = t
    G = parameters.group_count
    D, E = y[:G], y[G:]
    omega = np.asarray(parameters.omega_i, dtype=float)
    Omega_t = integration_phase.omega
    delta_t = integration_phase.delta

    ED = sum(omega_b * np.conj(E_b) * D_b for omega_b, E_b, D_b in zip(omega, E, D))
    DE = sum(omega_b * np.conj(D_b) * E_b for omega_b, D_b, E_b in zip(omega, D, E))

    dD = -0.5j * Omega_t * omega * E + 0.5 * parameters.Gamma * omega * ED * E
    dE = -0.5j * Omega_t * omega * D + 1j * delta_t * E - 0.5 * parameters.Gamma * omega * D * DE
    return np.concatenate([dD, dE])


def solve_mfe(
    parameters: MFESolverParameters,
    *,
    t_eval: Array,
    rtol: float = 1e-11,
    atol: float = 1e-13,
    method: str = "RK45",
) -> MFEResult:
    """
    Solve the group-resolved MFEs on the requested saved-time grid.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional array with at least two points.")
    if np.any(np.diff(t_eval) <= 0.0):
        raise ValueError("t_eval must be strictly increasing.")

    phase_protocol = parameters.phase_protocol
    integration_phases = phase_protocol.integration_phases
    if abs(float(t_eval[0])) > 1e-12:
        raise ValueError("The first t_eval point must be 0.0.")
    if abs(float(t_eval[-1]) - phase_protocol.total_duration) > 1e-9:
        raise ValueError("The last t_eval point must match the total protocol time.")

    zero_angles = (0.0,) * parameters.group_count
    y0 = amplitudes_from_initial_state(
        zero_angles,
        zero_angles,
        parameters,
    )
    integration_boundaries = phase_boundary_times(integration_phases)

    def rhs(t: float, y: Array) -> Array:
        phase_index = min(
            int(np.searchsorted(integration_boundaries, t, side="left")),
            len(integration_phases) - 1,
        )
        return mfe_rhs(t, y, parameters, integration_phases[phase_index])

    solution = solve_ivp(
        rhs,
        (float(t_eval[0]), float(t_eval[-1])),
        y0,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        method=method,
    )

    G = parameters.group_count
    D_groups = tuple(solution.y[:G])
    E_groups = tuple(solution.y[G:])
    result = MFEResult(
        t=t_eval,
        D_groups=D_groups,
        E_groups=E_groups,
        success=bool(solution.success),
        message=str(solution.message),
        parameters=parameters,
    )
    return result
