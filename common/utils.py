from __future__ import annotations

from typing import List, Sequence, Tuple

import numpy as np

from common.parser import Array, Phase


def phase_change_times(phases: Sequence[Phase]) -> Tuple[float, float]:
    if len(phases) < 2:
        raise ValueError("Need at least two phases to define change times.")
    t1 = phases[0].duration
    t2 = phases[0].duration + phases[1].duration
    return t1, t2


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
    vector within the active manifold.
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
    sz[valid] = -2.0 * Jz[valid] / N_active[valid]

    if np.any(sz < -1.0 - tol) or np.any(sz > 1.0 + tol):
        raise ValueError("sz values must lie in [-1, 1] to compute angles.")
    sz = np.clip(sz, -1.0, 1.0)

    theta = np.zeros_like(sz, dtype=float)
    theta[valid] = np.arccos(sz[valid])

    phi = np.arctan2(sy, sx)
    r_perp = np.sqrt(sx**2 + sy**2)
    phi[r_perp < tol] = 0.0

    return theta, phi, N_active, sx, sy, sz


def omega_c(N_J: int, Gamma: float) -> float:
    """Critical drive for the polarized-to-mixed transition at delta = 0."""

    return 0.5 * N_J * Gamma


def default_three_phase_protocol(
    T1: float,
    T2: float,
    T3: float,
    delta0: float,
    Omega0: float,
) -> List[Phase]:
    """Three-phase protocol."""
    return [
        Phase(duration=T1, omega=Omega0, delta=0.0, label="phase1"),
        Phase(duration=T2, omega=Omega0, delta=delta0, label="phase2"),
        Phase(duration=T3, omega=0.0, delta=0.0, label="phase3"),
    ]
