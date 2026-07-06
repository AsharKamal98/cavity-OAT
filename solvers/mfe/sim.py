from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from common.utils.phases import phase_values_at_time
from common.utils.parameters import omega_G_from_weighted_average
from parser.common import Array
from parser.mfe import (
    MFEInitialState,
    MFEResult,
    MFESolverParameters,
)
from solvers.mfe.utils import amplitudes_from_initial_state


def mfe_rhs(t: float, y: Array, parameters: MFESolverParameters) -> Array:
    """
    Right-hand side of the group-resolved mean-field equations.
    """
    G = parameters.group_count
    D, E = y[:G], y[G:]
    omega = np.asarray(parameters.omega_groups, dtype=float)
    Omega_t, delta_t = phase_values_at_time(t, parameters.phases)

    ED = sum(omega_b * np.conj(E_b) * D_b for omega_b, E_b, D_b in zip(omega, E, D))
    DE = sum(omega_b * np.conj(D_b) * E_b for omega_b, D_b, E_b in zip(omega, D, E))

    dD = -0.5j * Omega_t * omega * E + 0.5 * parameters.Gamma * omega * ED * E
    dE = -0.5j * Omega_t * omega * D + 1j * delta_t * E - 0.5 * parameters.Gamma * omega * D * DE
    return np.concatenate([dD, dE])


def solve_mfe(
    parameters: MFESolverParameters,
    initial_state: MFEInitialState,
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

    completed_parameters = parameters
    if len(parameters.omega_groups) == parameters.group_count - 1:
        omega_groups = tuple(parameters.omega_groups) + (
            omega_G_from_weighted_average(parameters.omega_groups, parameters.N_j_groups),
        )
        completed_parameters = parameters.model_copy(update={"omega_groups": omega_groups})

    y0 = amplitudes_from_initial_state(initial_state, completed_parameters)
    solution = solve_ivp(
        lambda t, y: mfe_rhs(t, y, completed_parameters),
        (float(t_eval[0]), float(t_eval[-1])),
        y0,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        method=method,
    )

    G = completed_parameters.group_count
    D_groups = tuple(solution.y[:G])
    E_groups = tuple(solution.y[G:])
    result = MFEResult(
        t=t_eval,
        D_groups=D_groups,
        E_groups=E_groups,
        success=bool(solution.success),
        message=str(solution.message),
        parameters=completed_parameters,
    )
    return result
