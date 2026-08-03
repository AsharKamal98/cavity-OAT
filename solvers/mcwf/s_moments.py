from __future__ import annotations

from typing import Optional

import numpy as np

from parser.common import Array
from parser.j_moments import JMomentSeries
from parser.mcwf import (
    TrajectoryEnsemble,
    TrajectoryEnsembleMetadata,
    TrajectoryResult,
    TrajectorySnapshot,
)
from parser.s_moments import SMomentSeries, SMomentSnapshot
from solvers.mcwf.squeezing_moments import (
    _apply_group_one_body_operator,
    _block_inner_product,
    _s_bloch_operators,
)
from solvers.mcwf.state_helpers import total_norm2
from solvers.mcwf.utils import map_with_optional_pool


def _validated_reference_data(
    metadata: TrajectoryEnsembleMetadata,
    j_moments: JMomentSeries,
) -> tuple[Array, tuple[int, ...], Array, Array]:
    """Validate and return the shared time grid, group sizes, and J angles."""
    N_groups = tuple(int(value) for value in metadata.Ni)
    group_count = len(N_groups)
    if group_count not in (1, 2):
        raise ValueError("MCWF S moments currently support one or two groups.")
    if j_moments.theta is None or j_moments.phi is None:
        raise ValueError("S moments require group-resolved J angles.")
    if len(j_moments.theta) != group_count or len(j_moments.phi) != group_count:
        raise ValueError("S moments require J-angle fields matching metadata.Ni.")

    t_ref = np.asarray(metadata.t_eval, dtype=float)
    if np.asarray(j_moments.t).shape != t_ref.shape or not np.allclose(
        j_moments.t,
        t_ref,
        atol=1e-12,
        rtol=0.0,
    ):
        raise ValueError(
            "j_moments and ensemble.metadata.t_eval must share a time grid."
        )

    theta_J_groups = np.asarray(j_moments.theta, dtype=float)
    phi_J_groups = np.asarray(j_moments.phi, dtype=float)
    expected_shape = (group_count, t_ref.size)
    if theta_J_groups.shape != expected_shape or (
        phi_J_groups.shape != expected_shape
    ):
        raise ValueError("Each group J-angle series must match t_eval.")

    return t_ref, N_groups, theta_J_groups, phi_J_groups


def _compute_snapshot_s_moments(
    snapshot: TrajectorySnapshot,
    *,
    theta_J: tuple[float, ...],
    phi_J: tuple[float, ...],
    N_groups: tuple[int, ...],
    tol: float,
) -> SMomentSnapshot:
    """Compute group-local S moments for one saved trajectory snapshot."""
    group_count = len(N_groups)
    if len(theta_J) != group_count or len(phi_J) != group_count:
        raise ValueError("Snapshot J angles must match the physical groups.")

    norm2 = total_norm2(snapshot.sector_blocks)
    if norm2 <= tol:
        zero_groups = tuple(0.0 for _ in N_groups)
        return SMomentSnapshot(
            t=snapshot.time,
            integration_phase_index=snapshot.integration_phase_index,
            x=zero_groups,
            y=zero_groups,
            z=zero_groups,
        )

    components = np.zeros((3, group_count), dtype=float)
    for group_index in range(group_count):
        S_x, S_y, S_z = _s_bloch_operators(
            theta_J[group_index],
            phi_J[group_index],
        )
        for component_index, operator in enumerate((S_x, S_y, S_z)):
            applied = _apply_group_one_body_operator(
                snapshot.sector_blocks,
                N_groups,
                group_index,
                operator,
            )
            components[component_index, group_index] = float(
                np.real(
                    _block_inner_product(
                        snapshot.sector_blocks,
                        applied,
                    )
                    / norm2
                )
            )

    return SMomentSnapshot(
        t=snapshot.time,
        integration_phase_index=snapshot.integration_phase_index,
        x=tuple(components[0]),
        y=tuple(components[1]),
        z=tuple(components[2]),
    )


def compute_trajectory_s_moments(
    trajectory: TrajectoryResult,
    metadata: TrajectoryEnsembleMetadata,
    j_moments: JMomentSeries,
    *,
    tol: float = 1e-12,
) -> SMomentSeries:
    """Convert one trajectory's saved snapshots into raw S-moment series."""
    if tol <= 0.0:
        raise ValueError("tol must be positive.")

    (
        t_ref,
        N_groups,
        theta_J_groups,
        phi_J_groups,
    ) = _validated_reference_data(metadata, j_moments)

    snapshot_times = np.asarray(
        [snapshot.time for snapshot in trajectory.snapshots],
        dtype=float,
    )
    if snapshot_times.shape != t_ref.shape or not np.allclose(
        snapshot_times,
        t_ref,
        atol=1e-12,
        rtol=0.0,
    ):
        raise ValueError(
            "Trajectory snapshots and metadata.t_eval must share a time grid."
        )

    s_moment_snapshots = [
        _compute_snapshot_s_moments(
            snapshot,
            theta_J=tuple(theta_J_groups[:, time_index]),
            phi_J=tuple(phi_J_groups[:, time_index]),
            N_groups=N_groups,
            tol=tol,
        )
        for time_index, snapshot in enumerate(trajectory.snapshots)
    ]

    def series(field_name: str, *, dtype=float) -> Array:
        return np.asarray(
            [getattr(snapshot, field_name) for snapshot in s_moment_snapshots],
            dtype=dtype,
        )

    def group_series(field_name: str) -> tuple[Array, ...]:
        return tuple(
            np.asarray(
                [
                    getattr(snapshot, field_name)[group_index]
                    for snapshot in s_moment_snapshots
                ],
                dtype=float,
            )
            for group_index in range(len(N_groups))
        )

    return SMomentSeries(
        t=series("t"),
        integration_phase_index=series("integration_phase_index", dtype=int),
        x=group_series("x"),
        y=group_series("y"),
        z=group_series("z"),
    )


def _compute_trajectory_s_moments_worker(
    args: tuple[
        TrajectoryResult,
        TrajectoryEnsembleMetadata,
        JMomentSeries,
        float,
    ],
) -> SMomentSeries:
    """Top-level worker for per-trajectory S-moment extraction."""
    trajectory, metadata, j_moments, tol = args
    return compute_trajectory_s_moments(
        trajectory,
        metadata,
        j_moments,
        tol=tol,
    )


def compute_average_s_moments(
    samples: list[SMomentSeries],
) -> SMomentSeries:
    """Average raw per-trajectory S moments on a shared saved-time grid."""
    if not samples:
        raise ValueError("No S-moment samples to average.")
    if samples[0].integration_phase_index is None:
        raise ValueError(
            "S-moment samples require an integration_phase_index grid."
        )

    t_ref = np.asarray(samples[0].t, dtype=float)
    integration_phase_ref = np.asarray(
        samples[0].integration_phase_index,
        dtype=int,
    )
    for sample in samples:
        if np.asarray(sample.t).shape != t_ref.shape or not np.allclose(
            sample.t,
            t_ref,
            atol=1e-12,
            rtol=0.0,
        ):
            raise ValueError("All S-moment samples must share the same t grid.")
        if sample.integration_phase_index is None or not np.array_equal(
            sample.integration_phase_index,
            integration_phase_ref,
        ):
            raise ValueError(
                "All S-moment samples must share the same "
                "integration_phase_index grid."
            )

    def mean_group_series(field_name: str) -> Optional[tuple[Array, ...]]:
        first_groups = getattr(samples[0], field_name)
        if first_groups is None:
            return None
        group_count = len(first_groups)
        for sample in samples:
            groups = getattr(sample, field_name)
            if groups is None or len(groups) != group_count:
                raise ValueError("All S-moment samples must have matching groups.")
        return tuple(
            np.mean(
                np.asarray(
                    [
                        getattr(sample, field_name)[group_index]
                        for sample in samples
                    ],
                    dtype=float,
                ),
                axis=0,
            )
            for group_index in range(group_count)
        )

    return SMomentSeries(
        t=t_ref,
        integration_phase_index=integration_phase_ref,
        x=mean_group_series("x"),
        y=mean_group_series("y"),
        z=mean_group_series("z"),
    )


def compute_mcwf_s_moments(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    *,
    tol: float = 1e-12,
    n_processes: int | None = None,
) -> SMomentSeries:
    """Compute ensemble-averaged, group-local S-sphere moments."""
    if tol <= 0.0:
        raise ValueError("tol must be positive.")
    if not ensemble.trajectories:
        raise ValueError("TrajectoryEnsemble is empty.")
    if ensemble.metadata is None:
        raise ValueError(
            "TrajectoryEnsemble.metadata is required to compute S moments."
        )

    _validated_reference_data(ensemble.metadata, j_moments)
    samples = map_with_optional_pool(
        _compute_trajectory_s_moments_worker,
        [
            (
                trajectory,
                ensemble.metadata,
                j_moments,
                tol,
            )
            for trajectory in ensemble.trajectories
        ],
        n_processes=n_processes,
        progress_desc="compute_mcwf_s_moments",
    )

    averaged = compute_average_s_moments(samples)
    SMomentSeries.attach_norm_spin_components_from_spin_components(
        averaged,
        tol=tol,
    )
    SMomentSeries.attach_angles_from_norm_spin_components(
        averaged,
        tol=tol,
    )
    return averaged


__all__ = [
    "SMomentSeries",
    "SMomentSnapshot",
    "compute_average_s_moments",
    "compute_mcwf_s_moments",
    "compute_trajectory_s_moments",
]
