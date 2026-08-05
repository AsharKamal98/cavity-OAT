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
from parser.s_moments import SBasisSeries, SMomentSeries, SMomentSnapshot
from solvers.mcwf.squeezing_moments import (
    _apply_group_one_body_operator,
    _block_inner_product,
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


def _build_s_basis(
    theta_J_groups: Array,
    phi_J_groups: Array,
) -> SBasisSeries:
    """Construct the shared effective S-qubit basis for all groups and times."""
    if theta_J_groups.shape != phi_J_groups.shape:
        raise ValueError("J theta and phi arrays must have matching shapes.")
    if theta_J_groups.ndim != 2:
        raise ValueError("J-angle arrays must have shape (num_groups, num_times).")

    num_times = theta_J_groups.shape[1]
    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    zero_states = tuple(
        np.broadcast_to(ket_u, (num_times, 3)).copy()
        for _ in range(theta_J_groups.shape[0])
    )
    one_states = tuple(
        np.column_stack(
            (
                np.zeros(num_times, dtype=np.complex128),
                np.cos(theta_J / 2.0),
                np.exp(-1j * phi_J) * np.sin(theta_J / 2.0),
            )
        )
        for theta_J, phi_J in zip(theta_J_groups, phi_J_groups)
    )
    return SBasisSeries(
        zero_state=zero_states,
        one_state=one_states,
    )


def _s_bloch_operators(
    zero_state: Array,
    one_state: Array,
) -> tuple[Array, Array, Array]:
    """Construct the three S-Bloch operators from one effective basis."""
    projector_zero = np.outer(zero_state, zero_state.conjugate())
    projector_one = np.outer(one_state, one_state.conjugate())
    zero_one = np.outer(zero_state, one_state.conjugate())
    one_zero = np.outer(one_state, zero_state.conjugate())

    S_x = 0.5 * (one_zero + zero_one)
    S_y = (one_zero - zero_one) / (2.0j)
    S_z = 0.5 * (projector_zero - projector_one)
    return S_x, S_y, S_z


def _build_s_bloch_operator_series(
    basis: SBasisSeries,
) -> tuple[Array, ...]:
    """Construct shared S-Bloch operators once per group and saved time."""
    return tuple(
        np.asarray(
            [
                _s_bloch_operators(zero_state, one_state)
                for zero_state, one_state in zip(
                    zero_state_series,
                    one_state_series,
                )
            ],
            dtype=np.complex128,
        )
        for zero_state_series, one_state_series in zip(
            basis.zero_state,
            basis.one_state,
        )
    )


def _compute_snapshot_s_moments(
    snapshot: TrajectorySnapshot,
    *,
    s_bloch_operators: tuple[Array, ...],
    N_groups: tuple[int, ...],
    tol: float,
) -> SMomentSnapshot:
    """Compute group-local S moments for one saved trajectory snapshot."""
    group_count = len(N_groups)
    if len(s_bloch_operators) != group_count:
        raise ValueError("S-Bloch operators must match the physical groups.")

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
        S_x, S_y, S_z = s_bloch_operators[group_index]
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
    s_bloch_operators: tuple[Array, ...],
    *,
    tol: float = 1e-12,
) -> SMomentSeries:
    """Convert one trajectory's saved snapshots into raw S-moment series."""
    if tol <= 0.0:
        raise ValueError("tol must be positive.")

    t_ref = np.asarray(metadata.t_eval, dtype=float)
    N_groups = tuple(int(value) for value in metadata.Ni)
    expected_operator_shape = (t_ref.size, 3, 3, 3)
    if len(s_bloch_operators) != len(N_groups) or any(
        np.asarray(operators).shape != expected_operator_shape
        for operators in s_bloch_operators
    ):
        raise ValueError(
            "Each group S-Bloch operator series must have shape "
            "(len(t_eval), 3, 3, 3)."
        )

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
            s_bloch_operators=tuple(
                operators[time_index]
                for operators in s_bloch_operators
            ),
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
        tuple[Array, ...],
        float,
    ],
) -> SMomentSeries:
    """Top-level worker for per-trajectory S-moment extraction."""
    trajectory, metadata, s_bloch_operators, tol = args
    return compute_trajectory_s_moments(
        trajectory,
        metadata,
        s_bloch_operators,
        tol=tol,
    )


def compute_average_s_moments(
    samples: list[SMomentSeries],
    *,
    basis: SBasisSeries | None = None,
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
        basis=basis,
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

    (
        _,
        _,
        theta_J_groups,
        phi_J_groups,
    ) = _validated_reference_data(ensemble.metadata, j_moments)
    basis = _build_s_basis(theta_J_groups, phi_J_groups)
    s_bloch_operators = _build_s_bloch_operator_series(basis)
    samples = map_with_optional_pool(
        _compute_trajectory_s_moments_worker,
        [
            (
                trajectory,
                ensemble.metadata,
                s_bloch_operators,
                tol,
            )
            for trajectory in ensemble.trajectories
        ],
        n_processes=n_processes,
        progress_desc="compute_mcwf_s_moments",
    )

    averaged = compute_average_s_moments(samples, basis=basis)
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
