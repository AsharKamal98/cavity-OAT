from __future__ import annotations

from typing import Dict, Optional

import numpy as np

from quantum_trajectories.utils import map_with_optional_pool
from common.parser import Phase
from quantum_trajectories.parser import (
    Array,
    JMomentSnapshot,
    JMomentSeries,
    TrajectoryEnsemble,
    TrajectoryResult,
    TrajectorySnapshot,
)
from quantum_trajectories.sim import build_phase_jump_operator_for_sector
from quantum_trajectories.state_helpers import total_norm2
from quantum_trajectories.operator_helpers import (
    build_sector_ops_for_key,
    total_active_atoms_in_sector,
)


# -----------------------------------------------------------------------------
# Single-trajectory J moments
# -----------------------------------------------------------------------------

_sector_ne_arrays: Dict[int, Array] = {}


def _cached_ne_array(Nj: int) -> Array:
    ne = _sector_ne_arrays.get(Nj)
    if ne is None:
        ne = np.arange(Nj + 1, dtype=float)
        _sector_ne_arrays[Nj] = ne
    return ne


def _compute_snapshot_j_moments(
    snapshot: TrajectorySnapshot,
    *,
    has_groups: bool,
    phases: list[Phase],
    Gamma: float,
    shifted_jump_operator: bool,
    omega_1: Optional[float],
    N1: Optional[int],
    N2: Optional[int],
) -> JMomentSnapshot:
    """
    Compute first-order J moments and jump rate for one saved snapshot.
    """
    blocks = snapshot.sector_blocks
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        zero_groups = (0.0, 0.0) if has_groups else None
        return JMomentSnapshot(
            t=snapshot.time,
            phase_index=snapshot.phase_index,
            Jx=0.0,
            Jy=0.0,
            Jz=0.0,
            N_e=0.0,
            N_j=0.0,
            jump_rate=0.0,
            J_drive=0.0,
            Jx_groups=zero_groups,
            Jy_groups=zero_groups,
            Jz_groups=zero_groups,
            N_e_groups=zero_groups,
            N_j_groups=zero_groups,
        )

    jx_total = 0.0
    jy_total = 0.0
    jz_total = 0.0
    ne_total = 0.0
    j_drive_total = 0.0
    jump_rate = 0.0
    nj_total = 0.0
    omega = phases[snapshot.phase_index].omega
    if has_groups:
        jx_groups = np.zeros(2, dtype=float)
        jy_groups = np.zeros(2, dtype=float)
        jz_groups = np.zeros(2, dtype=float)
        ne_groups = np.zeros(2, dtype=float)
        nj_groups = np.zeros(2, dtype=float)
    else:
        jx_groups = jy_groups = jz_groups = ne_groups = nj_groups = None

    for sector_key, psi in blocks.items():
        if psi.size == 0:
            continue

        # FIXME: use cache operators. Either compute new local cache operators
        # and store them in j_moments, or save cache operators from simulation,
        # e.g. by attatching to TrajectoryEnsamble.
        ops = build_sector_ops_for_key(
            sector_key,
            omega_1=omega_1,
            N1=N1,
            N2=N2,
        )
        psi_prob = np.abs(psi) ** 2
        psi_norm2 = float(np.vdot(psi, psi).real)

        jx_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)
        jy_total += float(np.vdot(psi, ops.J_y.dot(psi)).real)
        jz_total += float(np.vdot(psi, ops.J_z.dot(psi)).real)
        ne_total += float(np.vdot(psi, ops.N_e.dot(psi)).real)
        nj_total += total_active_atoms_in_sector(sector_key) * psi_norm2
        if ops.J_drive is not None:
            j_drive_total += float(np.vdot(psi, ops.J_drive.dot(psi)).real)
        else:
            j_drive_total += float(np.vdot(psi, ops.J_x.dot(psi)).real)

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

        if isinstance(sector_key, tuple):
            Nj1, Nj2 = sector_key
            ne1_diag = np.repeat(np.arange(Nj1 + 1, dtype=float), Nj2 + 1)
            ne2_diag = np.tile(np.arange(Nj2 + 1, dtype=float), Nj1 + 1)
            if (
                ops.J_x_groups is not None
                and ops.J_y_groups is not None
                and ops.J_z_groups is not None
            ):
                jx_groups[0] += float(np.vdot(psi, ops.J_x_groups[0].dot(psi)).real)
                jx_groups[1] += float(np.vdot(psi, ops.J_x_groups[1].dot(psi)).real)
                jy_groups[0] += float(np.vdot(psi, ops.J_y_groups[0].dot(psi)).real)
                jy_groups[1] += float(np.vdot(psi, ops.J_y_groups[1].dot(psi)).real)
                jz_groups[0] += float(np.vdot(psi, ops.J_z_groups[0].dot(psi)).real)
                jz_groups[1] += float(np.vdot(psi, ops.J_z_groups[1].dot(psi)).real)
            ne_groups[0] += float(np.dot(psi_prob, ne1_diag))
            ne_groups[1] += float(np.dot(psi_prob, ne2_diag))
            nj_groups[0] += Nj1 * psi_norm2
            nj_groups[1] += Nj2 * psi_norm2
        else:
            Nj = int(sector_key)
            ne = _cached_ne_array(Nj)
            if has_groups:
                jx_groups[0] += float(np.vdot(psi, ops.J_x.dot(psi)).real)
                jy_groups[0] += float(np.vdot(psi, ops.J_y.dot(psi)).real)
                jz_groups[0] += float(np.vdot(psi, ops.J_z.dot(psi)).real)
                ne_groups[0] += float(np.dot(psi_prob, ne))
                nj_groups[0] += Nj * psi_norm2

    # Guard against tiny negative values from floating-point noise.
    if jump_rate < 0.0 and abs(jump_rate) < 1e-10:
        jump_rate = 0.0

    return JMomentSnapshot(
        t=snapshot.time,
        phase_index=snapshot.phase_index,
        Jx=jx_total / norm2,
        Jy=jy_total / norm2,
        Jz=jz_total / norm2,
        N_e=ne_total / norm2,
        N_j=nj_total / norm2,
        jump_rate=jump_rate / norm2,
        J_drive=j_drive_total / norm2,
        Jx_groups=None if jx_groups is None else tuple(jx_groups / norm2),
        Jy_groups=None if jy_groups is None else tuple(jy_groups / norm2),
        Jz_groups=None if jz_groups is None else tuple(jz_groups / norm2),
        N_e_groups=None if ne_groups is None else tuple(ne_groups / norm2),
        N_j_groups=None if nj_groups is None else tuple(nj_groups / norm2),
    )


def compute_trajectory_j_moments(
    trajectory: TrajectoryResult,
    *,
    tol: float = 1e-12,
) -> JMomentSeries:
    """
    Convert saved snapshots into time series for Jx, Jy, Jz, Ne, jump rate,
    Nj, and the saved phase index.

    This uses the same first-order moment extraction as
    `trajectory_observables(...)`, but deliberately stops before constructing
    active-manifold angles or normalized Bloch components.
    """
    _ = tol  # Kept for API symmetry with later averaging/diagnostic functions.

    group_count = 2 if any(isinstance(key, tuple) for key in trajectory.sectors) else 0
    has_groups = group_count > 0
    j_moment_snapshots = [
        _compute_snapshot_j_moments(
            snap,
            has_groups=has_groups,
            phases=trajectory.phases,
            Gamma=trajectory.Gamma,
            shifted_jump_operator=trajectory.shifted_jump_operator,
            omega_1=trajectory.omega_1,
            N1=trajectory.N1,
            N2=trajectory.N2,
        )
        for snap in trajectory.snapshots
    ]

    def series(field_name: str, *, dtype=float) -> Array:
        return np.asarray(
            [getattr(snapshot, field_name) for snapshot in j_moment_snapshots],
            dtype=dtype,
        )

    def group_series(field_name: str) -> Optional[tuple[Array, ...]]:
        if not group_count:
            return None
        return tuple(
            np.asarray(
                [getattr(snapshot, field_name)[g] for snapshot in j_moment_snapshots],
                dtype=float,
            )
            for g in range(group_count)
        )

    return JMomentSeries(
        t=series("t"),
        phase_index=series("phase_index", dtype=int),
        Jx=series("Jx"),
        Jy=series("Jy"),
        Jz=series("Jz"),
        N_e=series("N_e"),
        N_j=series("N_j"),
        jump_rate=series("jump_rate"),
        J_drive=series("J_drive"),
        Jx_groups=group_series("Jx_groups"),
        Jy_groups=group_series("Jy_groups"),
        Jz_groups=group_series("Jz_groups"),
        N_e_groups=group_series("N_e_groups"),
        N_j_groups=group_series("N_j_groups"),
    )


# -----------------------------------------------------------------------------
# Ensemble trajectory J moments
# -----------------------------------------------------------------------------

def _compute_trajectory_j_moments_worker(args: tuple[TrajectoryResult, float]) -> JMomentSeries:
    """
    Top-level worker used to parallelize per-trajectory J-moment extraction.
    """
    trajectory, tol = args
    return compute_trajectory_j_moments(trajectory, tol=tol)


def compute_average_j_moments(samples: list[JMomentSeries]) -> JMomentSeries:
    """
    Average already-computed per-trajectory J moments on a shared time grid.
    """
    if len(samples) == 0:
        raise ValueError("No J-moment samples to average.")

    t_ref = np.asarray(samples[0].t, dtype=float)
    phase_ref = np.asarray(samples[0].phase_index, dtype=int)
    for sample in samples:
        if len(sample.t) != len(t_ref) or not np.allclose(sample.t, t_ref, atol=1e-12, rtol=0.0):
            raise ValueError("All J-moment samples must share the same t grid.")
        if len(sample.phase_index) != len(phase_ref) or not np.array_equal(sample.phase_index, phase_ref):
            raise ValueError("All J-moment samples must share the same phase_index grid.")

    def mean_series(field_name: str) -> Array:
        return np.mean(
            np.asarray([getattr(sample, field_name) for sample in samples], dtype=float),
            axis=0,
        )

    def mean_group_series(field_name: str) -> Optional[tuple[Array, ...]]:
        first_groups = getattr(samples[0], field_name)
        if first_groups is None:
            return None
        group_count = len(first_groups)
        for sample in samples:
            groups = getattr(sample, field_name)
            if groups is None or len(groups) != group_count:
                raise ValueError("All J-moment samples must have matching group fields.")
        return tuple(
            np.mean(
                np.asarray([getattr(sample, field_name)[g] for sample in samples], dtype=float),
                axis=0,
            )
            for g in range(group_count)
        )

    return JMomentSeries(
        t=t_ref,
        phase_index=phase_ref,
        Jx=mean_series("Jx"),
        Jy=mean_series("Jy"),
        Jz=mean_series("Jz"),
        N_e=mean_series("N_e"),
        N_j=mean_series("N_j"),
        jump_rate=mean_series("jump_rate"),
        J_drive=mean_series("J_drive"),
        Jx_groups=mean_group_series("Jx_groups"),
        Jy_groups=mean_group_series("Jy_groups"),
        Jz_groups=mean_group_series("Jz_groups"),
        N_e_groups=mean_group_series("N_e_groups"),
        N_j_groups=mean_group_series("N_j_groups"),
    )


def compute_ensemble_j_moments(
    ensemble: TrajectoryEnsemble,
    *,
    tol: float = 1e-12,
    reference: str = "first",
    n_processes: Optional[int] = None,
) -> JMomentSeries:
    """
    Compute trajectory-averaged first-order J moments for an ensemble.

    This is the moment-only counterpart of `ensemble_observables(...)`: it
    averages per-trajectory moment series without constructing theta/phi or
    sx/sy/sz.

    Returns
    -------
    JMomentSeries
        Trajectory-averaged first-order J moments on the common saved `t_eval`
        grid.
    """
    if len(ensemble.trajectories) == 0:
        raise ValueError("Ensemble is empty.")
    if reference != "first":
        raise ValueError("Currently only reference='first' is supported.")

    t_ref = np.asarray(
        [snap.time for snap in ensemble.trajectories[0].snapshots],
        dtype=float,
    )

    samples = map_with_optional_pool(
        _compute_trajectory_j_moments_worker,
        [(traj, tol) for traj in ensemble.trajectories],
        n_processes=n_processes,
        progress_desc="compute_ensemble_j_moments",
    )

    for sample in samples:
        if len(sample.t) != len(t_ref) or not np.allclose(sample.t, t_ref, atol=1e-12, rtol=0.0):
            raise ValueError(
                "All trajectories must share the same t_eval snapshot grid. "
                "Run the ensemble through the common num_snapshots API."
            )

    return compute_average_j_moments(samples)


__all__ = [
    "JMomentSnapshot",
    "JMomentSeries",
    "compute_average_j_moments",
    "compute_ensemble_j_moments",
    "compute_trajectory_j_moments",
]
