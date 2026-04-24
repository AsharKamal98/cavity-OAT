from __future__ import annotations

from quantum_trajectories.parser import (
    Phase,
)
from quantum_trajectories.state_helpers import (
    total_norm2,
)
from quantum_trajectories.operator_helpers import (
    build_sector_ops,
)

from utils import phase_change_times

from typing import Dict, Mapping, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt

Array = np.ndarray

# -----------------------------------------------------------------------------
# Pretty-print / inspection helpers
# -----------------------------------------------------------------------------

def expected_collective_components(blocks: Mapping[int, Array]) -> Tuple[float, float, float, float]:
    """
    Return normalized expectations (Jx, Jy, Jz, Ne) summed over all Nj blocks.
    """
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0, 0.0, 0.0, 0.0

    jx_total = 0.0
    jy_total = 0.0
    jz_total = 0.0
    ne_total = 0.0

    # Loop over Nj secotrs and sum expectations weighted by the probability of being in each sector.
    for Nj, psi in blocks.items():
        if psi.size == 0:
            continue
        ops = build_sector_ops(Nj)
        # For certain Nj sector, we have n=0,1,...,Nj states
        ne = np.arange(Nj + 1, dtype=float)

        jx_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)
        jy_total += float(np.vdot(psi, ((ops.J_plus - ops.J_minus) / (2.0j)).dot(psi)).real)
        # J_z value for state ∣n_e⟩ = (N_e - N_down) / 2 = (n_e - (Nj - n_e)) / 2 = n_e - Nj/2
        # J_z for N_j = Sum_n_e (probability of basis state ∣n_e⟩ in N_j sector times J_z value for state ∣n_e⟩)
        jz_total += float(np.dot(np.abs(psi) ** 2, ne - 0.5 * Nj))
        # ne for N_j = Sum_n_e (probability of basis state ∣n_e⟩ in N_j sector times n_e value for state ∣n_e⟩)
        ne_total += float(np.dot(np.abs(psi) ** 2, ne))

    # FIXME: correct?
    return (
        jx_total / norm2,
        jy_total / norm2,
        jz_total / norm2,
        ne_total / norm2,
    )


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

    # e,d manifold population that can actually contribute to the active-manifold angles.
    # N_active = <N_down + N_e> = 2 (N_e - J_z), since J_z = (N_e - N_down) / 2 in each fixed-Nj sector.
    N_active = 2.0 * (N_e - Jz)

    sx = np.zeros_like(Jx, dtype=float)
    sy = np.zeros_like(Jy, dtype=float)
    sz = np.zeros_like(Jz, dtype=float)

    # Compute angles
    valid = N_active > tol
    sx[valid] = 2.0 * Jx[valid] / N_active[valid]
    sy[valid] = 2.0 * Jy[valid] / N_active[valid]
    # Fit paper convention of polar angle theta = 0 at |down> and theta = pi at |e>.
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


def trajectory_observables(result: TrajectoryResult, *, tol: float = 1e-12) -> Dict[str, Array]:
    """
    Convert saved snapshots into time series for Jx, Jy, Jz, Ne, Nj, and the
    active-manifold angles.

    The returned theta and phi are computed in the {|down>, |e>} manifold using
    the expected active-manifold population N_active = <N_down + N_e>.
    """
    t = np.array([snap.time for snap in result.snapshots], dtype=float)
    jx = np.zeros_like(t)
    jy = np.zeros_like(t)
    jz = np.zeros_like(t)
    ne = np.zeros_like(t)
    nj = np.zeros_like(t)

    # Loop over snapshots and compute expectations for each one.
    for k, snap in enumerate(result.snapshots):
        # Average Jx, Jy, Jz, Ne over all sectors for this snapshot.
        jx_k, jy_k, jz_k, ne_k = expected_collective_components(snap.sector_blocks)
        jx[k] = jx_k
        jy[k] = jy_k
        jz[k] = jz_k
        ne[k] = ne_k
        nj[k] = sum(Nj * float(np.vdot(psi, psi).real) for Nj, psi in snap.sector_blocks.items())

    theta, phi, n_active, sx, sy, sz = active_manifold_angles(jx, jy, jz, ne, tol=tol)

    return {
        "t": t,
        "Jx": jx,
        "Jy": jy,
        "Jz": jz,
        "theta": theta,
        "phi": phi,
        "N_e": ne,
        "N_j": nj,
        "N_active": n_active,
        "sx": sx,
        "sy": sy,
        "sz": sz,
    }


def phase1_ss_angles_for_nj(Nj: int, Omega: float, Gamma: float):
    Omega_c = 0.5 * Nj * Gamma
    if Omega_c <= 0:
        raise ValueError("Omega_c must be positive.")
    ratio = Omega / Omega_c
    if abs(ratio) > 1.0:
        return np.nan, np.nan  # outside polarized phase
    cos_theta = np.sqrt(1.0 - ratio**2)
    theta_ss = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    phi_ss = 0.5 * np.pi
    return theta_ss, phi_ss


def _extract_observables_like_single(result, *, tol: float = 1e-12) -> Dict[str, Array]:
    """
    Accept either a real single-trajectory result or an averaged result already
    carrying observable arrays.
    """
    if hasattr(result, "observables"):
        return result.observables
    if isinstance(result, Mapping) and "t" in result:
        return result
    return trajectory_observables(result, tol=tol)


def plot_trajectory_angles_and_excitation(
    result,
    phases,
    *,
    output_path=None,
    show_phase1_ss=True,
    show_spread: bool = False,
):

    obs = _extract_observables_like_single(result)
    tlist = obs["t"]
    theta_mc = obs["theta"]
    phi_mc = obs["phi"]
    ne_mc = obs["N_e"]
    t_step1_end, t_step2_end = phase_change_times(phases)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=False)

    label = "qt"
    if hasattr(result, "ntraj") and getattr(result, "ntraj") is not None:
        label = f"MC avg ({result.ntraj} traj)"

    axes[0].plot(tlist, theta_mc, label=label, linewidth=1.8)
    axes[1].plot(tlist, phi_mc, label=label, linewidth=1.8)
    axes[2].plot(tlist, ne_mc, label=label, linewidth=1.8)

    if show_spread and hasattr(result, "std") and result.std is not None:
        for ax, key, mean in zip(
            axes,
            ["theta", "phi", "N_e"],
            [theta_mc, phi_mc, ne_mc],
        ):
            std = result.std.get(key)
            if std is not None:
                ax.fill_between(tlist, mean - std, mean + std, alpha=0.2)

    for ax in axes:
        ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.6)
        ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.6)
        ax.grid(alpha=0.3)

    if show_phase1_ss:
        Nj_ref = result.N // 2
        Omega1 = phases[0].omega
        theta_ss, phi_ss = phase1_ss_angles_for_nj(Nj_ref, Omega1, result.gamma)

        if np.isfinite(theta_ss):
            axes[0].hlines(
                y=theta_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($N_J=N/2$)",
            )
            axes[1].hlines(
                y=phi_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($\phi=\pi/2$)",
            )

    axes[0].set_xlabel(r"$\Gamma t$")
    axes[0].set_ylabel(r"Polar $\theta(t)$")
    axes[0].legend()

    axes[1].set_xlabel(r"$\Gamma t$")
    axes[1].set_ylabel(r"Azimuthal $\phi(t)$")
    axes[1].legend()

    axes[2].set_xlabel(r"$\Gamma t$")
    axes[2].set_ylabel(r"$\langle N_e(t)\rangle$")
    axes[2].legend()

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
    return fig, axes
