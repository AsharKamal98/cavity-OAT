from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

from parser.common import Array, Phase
from common.utils_moments import angles_from_norm_spin_components


def phase_change_times(phases: Sequence[Phase]) -> Tuple[float, float]:
    if len(phases) < 2:
        raise ValueError("Need at least two phases to define change times.")
    t1 = phases[0].duration
    t2 = phases[0].duration + phases[1].duration
    return t1, t2


def phase_values_at_time(t: float, phases: Sequence[Phase]) -> Tuple[float, float]:
    """
    Return phase-local (Omega, delta) values for a piecewise-constant protocol.
    """
    if not phases:
        raise ValueError("Need at least one phase.")

    t_value = float(t)
    total_time = float(sum(phase.duration for phase in phases))
    if t_value < 0.0 or t_value > total_time:
        raise ValueError(f"t must lie in [0, {total_time}], got {t_value}.")

    phase_end = 0.0
    for index, phase in enumerate(phases):
        phase_end += phase.duration
        if t_value <= phase_end or index == len(phases) - 1:
            return float(phase.omega), float(phase.delta)

    # The loop always returns for non-empty phases, but keep type checkers happy.
    phase = phases[-1]
    return float(phase.omega), float(phase.delta)

def phase1_ss_angles_for_nj(Nj: int, Omega: float, Gamma: float):
    Omega_c = 0.5 * Nj * Gamma
    if Omega_c <= 0:
        raise ValueError("Omega_c must be positive.")
    ratio = Omega / Omega_c
    if abs(ratio) > 1.0:
        return np.nan, np.nan
    cos_theta = np.sqrt(1.0 - ratio**2)
    theta_ss = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    phi_ss = 0.5 * np.pi
    return theta_ss, phi_ss


def active_manifold_angles(
    Jx: Array,
    Jy: Array,
    Jz: Array,
    N_e: Array,
    *,
    tol: float = 1e-12,
) -> Tuple[Array, Array, Array, Array, Array, Array]:
    """
    
    Compute Bloch-sphere angles inside the active {|down>, |e>} manifold.

    The active-manifold population is
        N_active = <N_down + N_e> = 2 (N_e - J_z),
    since J_z = (N_e - N_down) / 2 in each fixed-Nj sector.

    We then normalize the collective spin components by N_active so that the
    returned (sx, sy, sz) correspond to the averaged single-particle Bloch
    vector within the active manifold, with sz using the same sign convention
    as J_z. The polar angle uses theta = arccos(-sz), so the |down> state
    points to the north pole in the J-Bloch convention.
    """
    Jx = np.asarray(Jx, dtype=float)
    Jy = np.asarray(Jy, dtype=float)
    Jz = np.asarray(Jz, dtype=float)
    N_e = np.asarray(N_e, dtype=float)

    N_active = 2.0 * (N_e - Jz)

    sx = np.zeros_like(Jx, dtype=float)
    sy = np.zeros_like(Jy, dtype=float)
    sz = np.zeros_like(Jz, dtype=float)

    valid = N_active > tol
    sx[valid] = 2.0 * Jx[valid] / N_active[valid]
    sy[valid] = 2.0 * Jy[valid] / N_active[valid]
    sz[valid] = 2.0 * Jz[valid] / N_active[valid]

    if np.any(sz < -1.0 - tol) or np.any(sz > 1.0 + tol):
        raise ValueError("sz values must lie in [-1, 1] to compute angles.")
    sz = np.clip(sz, -1.0, 1.0)

    theta, phi = angles_from_norm_spin_components(sx, sy, sz, valid=valid, tol=tol)
    return theta, phi, N_active, sx, sy, sz
