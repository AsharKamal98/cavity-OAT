from __future__ import annotations

from typing import Dict, Mapping

import numpy as np

from common.utils import active_manifold_angles
from quantum_trajectories.parser import (
    Array,
    AveragedResult,
    ObservableSeries,
    TrajectoryResult,
    TrajectoryEnsemble,
)
from quantum_trajectories.sim import build_phase_jump_operator_for_sector
from quantum_trajectories.state_helpers import (
    total_norm2,
)
from quantum_trajectories.operator_helpers import (
    build_sector_ops,
)

# -----------------------------------------------------------------------------
# Single-trajectory observables
# -----------------------------------------------------------------------------

_sector_ne_arrays: Dict[int, Array] = {}
_sector_jz_diagonals: Dict[int, Array] = {}


def _cached_ne_array(Nj: int) -> Array:
    ne = _sector_ne_arrays.get(Nj)
    if ne is None:
        ne = np.arange(Nj + 1, dtype=float)
        _sector_ne_arrays[Nj] = ne
    return ne


def _cached_jz_diag(Nj: int) -> Array:
    jz_diag = _sector_jz_diagonals.get(Nj)
    if jz_diag is None:
        jz_diag = _cached_ne_array(Nj) - 0.5 * Nj
        _sector_jz_diagonals[Nj] = jz_diag
    return jz_diag


def expected_collective_components(blocks: Mapping[int, Array]) -> tuple[float, float, float, float]:
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
        ne = _cached_ne_array(Nj)
        jz_diag = _cached_jz_diag(Nj)
        psi_prob = np.abs(psi) ** 2

        jx_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)
        jy_total += float(np.vdot(psi, ops.J_y.dot(psi)).real)
        # J_z value for state ∣n_e⟩ = (N_e - N_down) / 2 = (n_e - (Nj - n_e)) / 2 = n_e - Nj/2
        # J_z for N_j = Sum_n_e (probability of basis state ∣n_e⟩ in N_j sector times J_z value for state ∣n_e⟩)
        jz_total += float(np.dot(psi_prob, jz_diag))
        # ne for N_j = Sum_n_e (probability of basis state ∣n_e⟩ in N_j sector times n_e value for state ∣n_e⟩)
        ne_total += float(np.dot(psi_prob, ne))

    # FIXME: correct?
    return (
        jx_total / norm2,
        jy_total / norm2,
        jz_total / norm2,
        ne_total / norm2,
    )


def jump_rate_for_blocks(
    blocks: Mapping[int, Array],
    *,
    gamma: float,
    omega: float,
    shifted_jump_operator: bool,
) -> float:
    """
    Evaluate the physical jump rate for one saved MCWF state.

    The rate is computed from the full block wavefunction and the same
    phase-dependent jump operator used by the solver:
        r(t) = gamma * sum_Nj <psi_Nj| l_Nj^dagger(t) l_Nj(t) |psi_Nj>.
    """
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0

    jump_rate = 0.0
    for Nj, psi in blocks.items():
        if psi.size == 0:
            continue
        ops = build_sector_ops(Nj)
        jump_operator = build_phase_jump_operator_for_sector(
            ops,
            omega,
            gamma,
            shifted_jump_operator=shifted_jump_operator,
        )
        # The MCWF code stores the unscaled jump operator l, so the physical
        # rate is gamma * <l^dagger l>.
        jumped = jump_operator.dot(psi)
        jump_rate += float(gamma * np.vdot(jumped, jumped).real)

    # Guard against tiny negative values from floating-point noise.
    if jump_rate < 0.0 and abs(jump_rate) < 1e-10:
        return 0.0
    return jump_rate / norm2


def trajectory_observables(result: TrajectoryResult, *, tol: float = 1e-12) -> ObservableSeries:
    """
    Convert saved snapshots into time series for Jx, Jy, Jz, Ne, jump rate,
    Nj, and the active-manifold angles.

    The returned theta and phi are computed in the {|down>, |e>} manifold using
    the expected active-manifold population N_active = <N_down + N_e>.
    """
    t = np.array([snap.time for snap in result.snapshots], dtype=float)
    jx = np.zeros_like(t)
    jy = np.zeros_like(t)
    jz = np.zeros_like(t)
    ne = np.zeros_like(t)
    jump_rate = np.zeros_like(t)
    nj = np.zeros_like(t)

    # Loop over snapshots and compute expectations for each one.
    for k, snap in enumerate(result.snapshots):
        # Average Jx, Jy, Jz, Ne over all sectors for this snapshot.
        jx_k, jy_k, jz_k, ne_k = expected_collective_components(snap.sector_blocks)
        jx[k] = jx_k
        jy[k] = jy_k
        jz[k] = jz_k
        ne[k] = ne_k
        # Recover the phase omega from the saved phase index so the shifted
        # picture uses the same time-local jump operator as the trajectory.
        jump_rate[k] = jump_rate_for_blocks(
            snap.sector_blocks,
            gamma=result.gamma,
            omega=result.phases[snap.phase_index].omega,
            shifted_jump_operator=result.shifted_jump_operator,
        )
        norm2 = total_norm2(snap.sector_blocks)
        if norm2 <= 1e-15:
            nj[k] = 0.0
        else:
            nj[k] = sum(
                Nj * float(np.vdot(psi, psi).real) for Nj, psi in snap.sector_blocks.items()
            ) / norm2

    theta, phi, n_active, sx, sy, sz = active_manifold_angles(jx, jy, jz, ne, tol=tol)

    return ObservableSeries(
        t=t,
        Jx=jx,
        Jy=jy,
        Jz=jz,
        N_e=ne,
        jump_rate=jump_rate,
        N_j=nj,
        N_active=n_active,
        theta=theta,
        phi=phi,
        sx=sx,
        sy=sy,
        sz=sz,
    )

def single_trajectory_to_averaged_result(
    result: TrajectoryResult,
    observables: ObservableSeries,
) -> AveragedResult:

    return AveragedResult(
        N=result.N,
        gamma=result.gamma,
        ntraj=None,   # change to 1 if you want the plot label to show "MC avg (1 traj)"
        observables=observables,
    )

# -----------------------------------------------------------------------------
# Ensemble trajectory observables and averaging
# -----------------------------------------------------------------------------

def _interp_series(t_src: Array, y_src: Array, t_ref: Array) -> Array:
    """
    Linear interpolation onto a common time grid.
    """
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    t_ref = np.asarray(t_ref, dtype=float)

    if t_src.ndim != 1 or y_src.ndim != 1:
        raise ValueError("t_src and y_src must be 1D arrays.")
    if len(t_src) != len(y_src):
        raise ValueError("t_src and y_src must have the same length.")

    # Remove duplicate times if present
    t_unique, idx = np.unique(t_src, return_index=True)
    y_unique = y_src[idx]

    return np.interp(t_ref, t_unique, y_unique)


def ensemble_observables(
    ensemble: TrajectoryEnsemble,
    *,
    tol: float = 1e-12,
    reference: str = "first",
) -> ObservableSeries:
    """
    Compute per-trajectory observables, interpolate them to a common time grid,
    average Jx/Jy/Jz/N_e/jump_rate/N_j across trajectories, then compute
    theta/phi and sx/sy/sz from those averaged observables.

    Returns
    -------
    ObservableSeries
        Mean observables in the standard single-trajectory field names, with
        standard deviations filled for:
        Jx, Jy, Jz, N_e, N_j, N_active

        For now, theta_std and phi_std are not included.
    """
    if len(ensemble.trajectories) == 0:
        raise ValueError("Ensemble is empty.")

    if reference != "first":
        raise ValueError("Currently only reference='first' is supported.")

    t_ref = np.asarray(
        [snap.time for snap in ensemble.trajectories[0].snapshots],
        dtype=float,
    )

    Jx_list = []
    Jy_list = []
    Jz_list = []
    Ne_list = []
    jump_rate_list = []
    Nj_list = []

    for traj in ensemble.trajectories:
        obs = trajectory_observables(traj, tol=tol)

        Jx_list.append(_interp_series(obs.t, obs.Jx, t_ref))
        Jy_list.append(_interp_series(obs.t, obs.Jy, t_ref))
        Jz_list.append(_interp_series(obs.t, obs.Jz, t_ref))
        Ne_list.append(_interp_series(obs.t, obs.N_e, t_ref))
        # Jump rate is averaged on the same reference grid as the other
        # observables so ensemble means/stds line up in time.
        jump_rate_list.append(_interp_series(obs.t, obs.jump_rate, t_ref))
        Nj_list.append(_interp_series(obs.t, obs.N_j, t_ref))

    Jx_arr = np.asarray(Jx_list, dtype=float)
    Jy_arr = np.asarray(Jy_list, dtype=float)
    Jz_arr = np.asarray(Jz_list, dtype=float)
    Ne_arr = np.asarray(Ne_list, dtype=float)
    jump_rate_arr = np.asarray(jump_rate_list, dtype=float)
    Nj_arr = np.asarray(Nj_list, dtype=float)

    Jx_mean = np.mean(Jx_arr, axis=0)
    Jy_mean = np.mean(Jy_arr, axis=0)
    Jz_mean = np.mean(Jz_arr, axis=0)
    N_e_mean = np.mean(Ne_arr, axis=0)
    jump_rate_mean = np.mean(jump_rate_arr, axis=0)
    N_j_mean = np.mean(Nj_arr, axis=0)

    Jx_std = np.std(Jx_arr, axis=0, ddof=0)
    Jy_std = np.std(Jy_arr, axis=0, ddof=0)
    Jz_std = np.std(Jz_arr, axis=0, ddof=0)
    N_e_std = np.std(Ne_arr, axis=0, ddof=0)
    jump_rate_std = np.std(jump_rate_arr, axis=0, ddof=0)
    N_j_std = np.std(Nj_arr, axis=0, ddof=0)

    theta, phi, N_active, sx, sy, sz = active_manifold_angles(
        Jx_mean,
        Jy_mean,
        Jz_mean,
        N_e_mean,
        tol=tol,
    )

    # # Compute N_active per trajectory too, so its std matches the same ensemble logic
    # Nactive_list = []
    # for i in range(len(Jx_arr)):
    #     _, _, N_active_i, _, _, _ = active_manifold_angles(
    #         Jx_arr[i],
    #         Jy_arr[i],
    #         Jz_arr[i],
    #         Ne_arr[i],
    #         tol=tol,
    #     )
    #     Nactive_list.append(N_active_i)

    # Nactive_arr = np.asarray(Nactive_list, dtype=float)
    # N_active_std = np.std(Nactive_arr, axis=0, ddof=0)

    return ObservableSeries(
        t=t_ref,
        Jx=Jx_mean,
        Jy=Jy_mean,
        Jz=Jz_mean,
        N_e=N_e_mean,
        jump_rate=jump_rate_mean,
        N_j=N_j_mean,
        N_active=N_active,
        theta=theta,
        phi=phi,
        sx=sx,
        sy=sy,
        sz=sz,
        Jx_std=Jx_std,
        Jy_std=Jy_std,
        Jz_std=Jz_std,
        N_e_std=N_e_std,
        jump_rate_std=jump_rate_std,
        N_j_std=N_j_std,
        N_active_std=0,
    )


def make_averaged_result(
    ensemble: TrajectoryEnsemble,
    observables: ObservableSeries,
) -> AveragedResult:
    first = ensemble.trajectories[0]
    return AveragedResult(
        N=first.N,
        gamma=first.gamma,
        ntraj=len(ensemble.trajectories),
        observables=observables,
    )
