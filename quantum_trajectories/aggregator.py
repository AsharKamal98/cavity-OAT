from __future__ import annotations

from typing import Dict, Mapping, Optional

import numpy as np

from common.utils import active_manifold_angles
from quantum_trajectories.utils import map_with_optional_pool
from quantum_trajectories.parser import (
    Array,
    AveragedResult,
    ObservableSeries,
    SectorKey,
    TrajectoryResult,
    TrajectoryEnsemble,
)
from quantum_trajectories.sim import build_phase_jump_operator_for_sector
from quantum_trajectories.state_helpers import (
    total_norm2,
)
from quantum_trajectories.operator_helpers import (
    build_sector_ops_for_key,
    total_active_atoms_in_sector,
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


def _expected_collective_components_detailed(
    blocks: Mapping[SectorKey, Array],
    *,
    omega_1: Optional[float] = None,
    N1: Optional[int] = None,
    N2: Optional[int] = None,
) -> dict:
    """
    Return normalized total and group-resolved collective expectations.
    """
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return {
            "Jx": 0.0,
            "Jy": 0.0,
            "Jz": 0.0,
            "N_e": 0.0,
            "Jx_drive": 0.0,
            "Jx_groups": (0.0, 0.0) if any(isinstance(key, tuple) for key in blocks) else None,
            "Jy_groups": (0.0, 0.0) if any(isinstance(key, tuple) for key in blocks) else None,
            "Jz_groups": (0.0, 0.0) if any(isinstance(key, tuple) for key in blocks) else None,
            "N_e_groups": (0.0, 0.0) if any(isinstance(key, tuple) for key in blocks) else None,
        }

    jx_total = 0.0
    jy_total = 0.0
    jz_total = 0.0
    ne_total = 0.0
    jx_drive_total = 0.0
    any_inhomogeneous = any(isinstance(key, tuple) for key in blocks)
    if any_inhomogeneous:
        jx_groups = np.zeros(2, dtype=float)
        jy_groups = np.zeros(2, dtype=float)
        jz_groups = np.zeros(2, dtype=float)
        ne_groups = np.zeros(2, dtype=float)
    else:
        jx_groups = jy_groups = jz_groups = ne_groups = None

    for sector_key, psi in blocks.items():
        if psi.size == 0:
            continue
        ops = build_sector_ops_for_key(sector_key, omega_1=omega_1, N1=N1, N2=N2)
        psi_prob = np.abs(psi) ** 2

        jx_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)
        jy_total += float(np.vdot(psi, ops.J_y.dot(psi)).real)
        ne_total += float(np.vdot(psi, ops.N_e.dot(psi)).real)
        if ops.J_x_drive is not None:
            jx_drive_total += float(np.vdot(psi, ops.J_x_drive.dot(psi)).real)
        else:
            jx_drive_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)

        if isinstance(sector_key, tuple):
            Nj1, Nj2 = sector_key
            ne1_diag = np.repeat(np.arange(Nj1 + 1, dtype=float), Nj2 + 1)
            ne2_diag = np.tile(np.arange(Nj2 + 1, dtype=float), Nj1 + 1)
            jz1_diag = ne1_diag - 0.5 * Nj1
            jz2_diag = ne2_diag - 0.5 * Nj2
            jz_total += float(np.dot(psi_prob, jz1_diag + jz2_diag))
            if ops.J_x_groups is not None and ops.J_y_groups is not None:
                jx_groups[0] += float(np.vdot(psi, ops.J_x_groups[0].dot(psi)).real)
                jx_groups[1] += float(np.vdot(psi, ops.J_x_groups[1].dot(psi)).real)
                jy_groups[0] += float(np.vdot(psi, ops.J_y_groups[0].dot(psi)).real)
                jy_groups[1] += float(np.vdot(psi, ops.J_y_groups[1].dot(psi)).real)
            jz_groups[0] += float(np.dot(psi_prob, jz1_diag))
            jz_groups[1] += float(np.dot(psi_prob, jz2_diag))
            ne_groups[0] += float(np.dot(psi_prob, ne1_diag))
            ne_groups[1] += float(np.dot(psi_prob, ne2_diag))
        else:
            Nj = int(sector_key)
            ne = _cached_ne_array(Nj)
            jz_diag = _cached_jz_diag(Nj)
            jz_total += float(np.dot(psi_prob, jz_diag))
            if any_inhomogeneous:
                jx_groups[0] += float(np.vdot(psi, ops.J_x.dot(psi)).real)
                jy_groups[0] += float(np.vdot(psi, ops.J_y.dot(psi)).real)
                jz_groups[0] += float(np.dot(psi_prob, jz_diag))
                ne_groups[0] += float(np.dot(psi_prob, ne))

    out = {
        "Jx": jx_total / norm2,
        "Jy": jy_total / norm2,
        "Jz": jz_total / norm2,
        "N_e": ne_total / norm2,
        "Jx_drive": jx_drive_total / norm2,
        "Jx_groups": None if jx_groups is None else tuple(jx_groups / norm2),
        "Jy_groups": None if jy_groups is None else tuple(jy_groups / norm2),
        "Jz_groups": None if jz_groups is None else tuple(jz_groups / norm2),
        "N_e_groups": None if ne_groups is None else tuple(ne_groups / norm2),
    }
    return out


def expected_collective_components(
    blocks: Mapping[SectorKey, Array],
    *,
    omega_1: Optional[float] = None,
    N1: Optional[int] = None,
    N2: Optional[int] = None,
) -> tuple[float, float, float, float]:
    detailed = _expected_collective_components_detailed(
        blocks,
        omega_1=omega_1,
        N1=N1,
        N2=N2,
    )
    return detailed["Jx"], detailed["Jy"], detailed["Jz"], detailed["N_e"]


def jump_rate_for_blocks(
    blocks: Mapping[SectorKey, Array],
    *,
    Gamma: float,
    omega: float,
    shifted_jump_operator: bool,
    omega_1: Optional[float] = None,
    N1: Optional[int] = None,
    N2: Optional[int] = None,
) -> float:
    """
    Evaluate the physical jump rate for one saved MCWF state.

    The rate is computed from the full block wavefunction and the same
    phase-dependent jump operator used by the solver:
        r(t) = Gamma * sum_Nj <psi_Nj| l_Nj^dagger(t) l_Nj(t) |psi_Nj>.
    """
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0

    jump_rate = 0.0
    for sector_key, psi in blocks.items():
        if psi.size == 0:
            continue
        ops = build_sector_ops_for_key(sector_key, omega_1=omega_1, N1=N1, N2=N2)
        jump_operator = build_phase_jump_operator_for_sector(
            ops,
            omega,
            Gamma,
            shifted_jump_operator=shifted_jump_operator,
        )
        # The MCWF code stores the unscaled jump operator l, so the physical
        # rate is Gamma * <l^dagger l>.
        jumped = jump_operator.dot(psi)
        jump_rate += float(Gamma * np.vdot(jumped, jumped).real)

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
    jx_drive = np.zeros_like(t)
    jump_rate = np.zeros_like(t)
    nj = np.zeros_like(t)
    group_count = 2 if any(isinstance(key, tuple) for key in result.sectors) else 0
    jx_groups = np.zeros((group_count, len(t)), dtype=float) if group_count else None
    jy_groups = np.zeros((group_count, len(t)), dtype=float) if group_count else None
    jz_groups = np.zeros((group_count, len(t)), dtype=float) if group_count else None
    ne_groups = np.zeros((group_count, len(t)), dtype=float) if group_count else None

    # Loop over snapshots and compute expectations for each one.
    for k, snap in enumerate(result.snapshots):
        # Average Jx, Jy, Jz, Ne over all sectors for this snapshot.
        comp = _expected_collective_components_detailed(
            snap.sector_blocks,
            omega_1=result.omega_1,
            N1=result.N1,
            N2=result.N2,
        )
        jx[k] = comp["Jx"]
        jy[k] = comp["Jy"]
        jz[k] = comp["Jz"]
        ne[k] = comp["N_e"]
        jx_drive[k] = comp["Jx_drive"]
        if group_count:
            for g in range(group_count):
                jx_groups[g, k] = comp["Jx_groups"][g]
                jy_groups[g, k] = comp["Jy_groups"][g]
                jz_groups[g, k] = comp["Jz_groups"][g]
                ne_groups[g, k] = comp["N_e_groups"][g]
        # Recover the phase omega from the saved phase index so the shifted
        # picture uses the same time-local jump operator as the trajectory.
        jump_rate[k] = jump_rate_for_blocks(
            snap.sector_blocks,
            Gamma=result.Gamma,
            omega=result.phases[snap.phase_index].omega,
            shifted_jump_operator=result.shifted_jump_operator,
            omega_1=result.omega_1,
            N1=result.N1,
            N2=result.N2,
        )
        norm2 = total_norm2(snap.sector_blocks)
        if norm2 <= 1e-15:
            nj[k] = 0.0
        else:
            nj[k] = sum(
                total_active_atoms_in_sector(sector_key) * float(np.vdot(psi, psi).real)
                for sector_key, psi in snap.sector_blocks.items()
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
        Jx_drive=jx_drive,
        Jx_groups=None if jx_groups is None else tuple(jx_groups[g] for g in range(group_count)),
        Jy_groups=None if jy_groups is None else tuple(jy_groups[g] for g in range(group_count)),
        Jz_groups=None if jz_groups is None else tuple(jz_groups[g] for g in range(group_count)),
        N_e_groups=None if ne_groups is None else tuple(ne_groups[g] for g in range(group_count)),
    )

def single_trajectory_to_averaged_result(
    result: TrajectoryResult,
    observables: ObservableSeries,
) -> AveragedResult:

    return AveragedResult(
        N=result.N,
        Gamma=result.Gamma,
        ntraj=None,   # change to 1 if you want the plot label to show "MC avg (1 traj)"
        observables=observables,
    )

# -----------------------------------------------------------------------------
# Ensemble trajectory observables and averaging
# -----------------------------------------------------------------------------

def _trajectory_observables_worker(args: tuple[TrajectoryResult, float]) -> ObservableSeries:
    """
    Top-level worker used to parallelize per-trajectory observable extraction.
    """
    traj, tol = args
    return trajectory_observables(traj, tol=tol)


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
    n_processes: Optional[int] = None,
) -> ObservableSeries:
    """
    Compute per-trajectory observables on the shared MCWF t_eval grid,
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
    Jx_drive_list = []
    Ne_list = []
    jump_rate_list = []
    Nj_list = []
    collect_groups = False
    Jx_group_lists = None
    Jy_group_lists = None
    Jz_group_lists = None
    Ne_group_lists = None

    per_traj_obs = map_with_optional_pool(
        _trajectory_observables_worker,
        [(traj, tol) for traj in ensemble.trajectories],
        n_processes=n_processes,
        progress_desc="ensemble_observables",
    )

    for obs in per_traj_obs:
        if len(obs.t) != len(t_ref) or not np.allclose(obs.t, t_ref, atol=1e-12, rtol=0.0):
            raise ValueError(
                "All trajectories must share the same t_eval snapshot grid. "
                "Run the ensemble through the common num_snapshots API."
            )

        Jx_list.append(obs.Jx)
        Jy_list.append(obs.Jy)
        Jz_list.append(obs.Jz)
        if obs.Jx_drive is not None:
            Jx_drive_list.append(obs.Jx_drive)
        Ne_list.append(obs.N_e)
        jump_rate_list.append(obs.jump_rate)
        Nj_list.append(obs.N_j)
        if obs.Jx_groups is not None:
            if not collect_groups:
                group_count = len(obs.Jx_groups)
                Jx_group_lists = [[] for _ in range(group_count)]
                Jy_group_lists = [[] for _ in range(group_count)]
                Jz_group_lists = [[] for _ in range(group_count)]
                Ne_group_lists = [[] for _ in range(group_count)]
                collect_groups = True
            for g in range(len(obs.Jx_groups)):
                Jx_group_lists[g].append(obs.Jx_groups[g])
                Jy_group_lists[g].append(obs.Jy_groups[g])
                Jz_group_lists[g].append(obs.Jz_groups[g])
                Ne_group_lists[g].append(obs.N_e_groups[g])

    Jx_arr = np.asarray(Jx_list, dtype=float)
    Jy_arr = np.asarray(Jy_list, dtype=float)
    Jz_arr = np.asarray(Jz_list, dtype=float)
    Ne_arr = np.asarray(Ne_list, dtype=float)
    Jx_drive_arr = np.asarray(Jx_drive_list, dtype=float) if Jx_drive_list else None
    jump_rate_arr = np.asarray(jump_rate_list, dtype=float)
    Nj_arr = np.asarray(Nj_list, dtype=float)

    Jx_mean = np.mean(Jx_arr, axis=0)
    Jy_mean = np.mean(Jy_arr, axis=0)
    Jz_mean = np.mean(Jz_arr, axis=0)
    N_e_mean = np.mean(Ne_arr, axis=0)
    Jx_drive_mean = np.mean(Jx_drive_arr, axis=0) if Jx_drive_arr is not None else None
    jump_rate_mean = np.mean(jump_rate_arr, axis=0)
    N_j_mean = np.mean(Nj_arr, axis=0)

    Jx_std = np.std(Jx_arr, axis=0, ddof=0)
    Jy_std = np.std(Jy_arr, axis=0, ddof=0)
    Jz_std = np.std(Jz_arr, axis=0, ddof=0)
    N_e_std = np.std(Ne_arr, axis=0, ddof=0)
    Jx_drive_std = np.std(Jx_drive_arr, axis=0, ddof=0) if Jx_drive_arr is not None else None
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

    Jx_groups_mean = None
    Jy_groups_mean = None
    Jz_groups_mean = None
    N_e_groups_mean = None
    if collect_groups:
        Jx_groups_mean = tuple(np.mean(np.asarray(group_vals, dtype=float), axis=0) for group_vals in Jx_group_lists)
        Jy_groups_mean = tuple(np.mean(np.asarray(group_vals, dtype=float), axis=0) for group_vals in Jy_group_lists)
        Jz_groups_mean = tuple(np.mean(np.asarray(group_vals, dtype=float), axis=0) for group_vals in Jz_group_lists)
        N_e_groups_mean = tuple(np.mean(np.asarray(group_vals, dtype=float), axis=0) for group_vals in Ne_group_lists)

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
        Jx_drive=Jx_drive_mean,
        Jx_groups=Jx_groups_mean,
        Jy_groups=Jy_groups_mean,
        Jz_groups=Jz_groups_mean,
        N_e_groups=N_e_groups_mean,
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
        Gamma=first.Gamma,
        ntraj=len(ensemble.trajectories),
        observables=observables,
    )
