from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

from common.parser import Array
from common.utils import phase_values_at_time
from mfe.parser import (
    MFEInitialState,
    MFEObservableSeries,
    MFEResult,
    MFESolverParameters,
)


def amplitudes_from_initial_state(
    initial_state: MFEInitialState,
    parameters: MFESolverParameters,
) -> Array:
    """
    Build the flat complex solver vector from initial J-sphere angles.
    """
    if len(initial_state.theta_groups) != parameters.group_count:
        raise ValueError("Initial-state group count must match parameters.")

    theta = np.asarray(initial_state.theta_groups, dtype=float)
    phi = np.asarray(initial_state.phi_groups, dtype=float)
    N_j = np.asarray(parameters.N_j_groups, dtype=float)

    D = np.sqrt(N_j) * np.cos(0.5 * theta)
    E = np.sqrt(N_j) * np.exp(-1j * phi) * np.sin(0.5 * theta)
    return np.concatenate([D, E]).astype(np.complex128)


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


def angles_from_amplitudes(
    D_groups: tuple[Array, ...],
    E_groups: tuple[Array, ...],
    *,
    tol: float = 1e-12,
) -> tuple[tuple[Array, ...], tuple[Array, ...], tuple[Array, ...]]:
    """
    Convert solved amplitudes into N_j, theta_J, and phi_J time series.
    """
    N_j_groups = []
    theta_groups = []
    phi_groups = []

    for D, E in zip(D_groups, E_groups):
        D = np.asarray(D, dtype=np.complex128)
        E = np.asarray(E, dtype=np.complex128)
        N_j = np.abs(D) ** 2 + np.abs(E) ** 2

        ratio = np.zeros_like(N_j, dtype=float)
        valid = N_j > tol
        ratio[valid] = (np.abs(D[valid]) ** 2 - np.abs(E[valid]) ** 2) / N_j[valid]
        theta = np.zeros_like(N_j, dtype=float)
        theta[valid] = np.arccos(np.clip(ratio[valid], -1.0, 1.0))

        phi = np.angle(D) - np.angle(E)
        phi[~valid] = 0.0

        N_j_groups.append(N_j)
        theta_groups.append(theta)
        phi_groups.append(phi)

    return tuple(N_j_groups), tuple(theta_groups), tuple(phi_groups)


def compute_mfe_observables(
    result: MFEResult,
    *,
    tol: float = 1e-12,
) -> MFEObservableSeries:
    """
    Build observable time series from an MFE result.
    """
    N_j_groups, theta_groups, phi_groups = angles_from_amplitudes(
        result.D_groups,
        result.E_groups,
        tol=tol,
    )
    return MFEObservableSeries(
        t=result.t,
        D_groups=result.D_groups,
        E_groups=result.E_groups,
        N_j_groups=N_j_groups,
        theta_groups=theta_groups,
        phi_groups=phi_groups,
    )


def solve_mfe(
    parameters: MFESolverParameters,
    initial_state: MFEInitialState,
    *,
    t_eval: Array,
    rtol: float = 1e-9,
    atol: float = 1e-11,
    method: str = "RK45",
) -> MFEResult:
    """
    Solve the group-resolved MFEs on the requested saved-time grid.
    """
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional array with at least two points.")

    y0 = amplitudes_from_initial_state(initial_state, parameters)
    solution = solve_ivp(
        lambda t, y: mfe_rhs(t, y, parameters),
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
    result.observables = compute_mfe_observables(result)
    return result
