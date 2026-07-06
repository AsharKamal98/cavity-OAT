from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.mfe import MFESolverParameters


def amplitudes_from_initial_state(
    theta_groups: tuple[float, ...],
    phi_groups: tuple[float, ...],
    parameters: MFESolverParameters,
) -> Array:
    """
    Build the flat complex solver vector from initial J-sphere angles.
    """
    if len(theta_groups) != parameters.group_count or len(phi_groups) != parameters.group_count:
        raise ValueError("Initial-state group count must match parameters.")

    theta = np.asarray(theta_groups, dtype=float)
    phi = np.asarray(phi_groups, dtype=float)
    N_j = 0.5 * np.asarray(parameters.Ni, dtype=float)

    D = np.sqrt(N_j) * np.cos(0.5 * theta)
    E = np.sqrt(N_j) * np.exp(-1j * phi) * np.sin(0.5 * theta)
    return np.concatenate([D, E]).astype(np.complex128)


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
        phi = (phi + np.pi) % (2.0 * np.pi) - np.pi
        phi[~valid] = 0.0

        N_j_groups.append(N_j)
        theta_groups.append(theta)
        phi_groups.append(phi)

    return tuple(N_j_groups), tuple(theta_groups), tuple(phi_groups)
