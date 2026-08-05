from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from parser.common import Array
from parser.j_moments import JMomentSeries
from parser.mcwf import SectorKey, TrajectoryEnsemble, TrajectoryResult
from parser.s_moments import SMomentSeries
from parser.squeezing import SqueezingMomentSeries
from solvers.mcwf.operator_helpers import build_sector_ops
from solvers.mcwf.state_helpers import total_norm2
from solvers.mcwf.utils import map_with_optional_pool


# -----------------------------------------------------------------------------
# Reduced-basis collective one-body operators
# -----------------------------------------------------------------------------


def _block_inner_product(
    left: Mapping[SectorKey, Array],
    right: Mapping[SectorKey, Array],
) -> complex:
    """Return the inner product of two reduced sector-block states."""
    return sum(
        np.vdot(left[sector_key], right[sector_key])
        for sector_key in set(left).intersection(right)
    )


def _accumulate_block(
    out: dict[SectorKey, Array],
    sector_key: SectorKey,
    contribution: Array,
) -> None:
    """Add one contribution to a reduced sector-block state."""
    if contribution.size == 0:
        return
    if sector_key in out:
        out[sector_key] = out[sector_key] + contribution
    else:
        out[sector_key] = np.asarray(
            contribution,
            dtype=np.complex128,
        ).copy()


def _apply_single_group_transition(
    blocks: Mapping[SectorKey, Array],
    N_group: int,
    mu: str,
    nu: str,
) -> dict[SectorKey, Array]:
    """Apply one collective transition to a single-group MCWF state."""
    out: dict[SectorKey, Array] = {}

    for sector_key, psi in blocks.items():
        if isinstance(sector_key, tuple):
            raise ValueError("Single-group squeezing requires scalar sector keys.")

        Nj = int(sector_key)
        psi = np.asarray(psi, dtype=np.complex128)
        ne = np.arange(Nj + 1, dtype=float)
        n_u = N_group - Nj

        if mu == "u" and nu == "u":
            _accumulate_block(out, Nj, n_u * psi)
        elif mu == "d" and nu == "d":
            _accumulate_block(out, Nj, (Nj - ne) * psi)
        elif mu == "e" and nu == "e":
            _accumulate_block(out, Nj, ne * psi)
        elif mu == "d" and nu == "e":
            _accumulate_block(out, Nj, build_sector_ops(Nj).Jm.dot(psi))
        elif mu == "e" and nu == "d":
            _accumulate_block(out, Nj, build_sector_ops(Nj).Jp.dot(psi))
        elif mu == "u" and nu == "d":
            if Nj > 0:
                target = (
                    np.sqrt((N_group - Nj + 1.0) * (Nj - ne[:-1]))
                    * psi[:-1]
                )
                _accumulate_block(out, Nj - 1, target)
        elif mu == "d" and nu == "u":
            if n_u > 0:
                target = np.zeros(Nj + 2, dtype=np.complex128)
                target[:-1] = np.sqrt(n_u * (Nj - ne + 1.0)) * psi
                _accumulate_block(out, Nj + 1, target)
        elif mu == "u" and nu == "e":
            if Nj > 0:
                target = (
                    np.sqrt((N_group - Nj + 1.0) * ne[1:])
                    * psi[1:]
                )
                _accumulate_block(out, Nj - 1, target)
        elif mu == "e" and nu == "u":
            if n_u > 0:
                target = np.zeros(Nj + 2, dtype=np.complex128)
                target[1:] = np.sqrt(n_u * (ne + 1.0)) * psi
                _accumulate_block(out, Nj + 1, target)
        else:
            raise ValueError(f"Unsupported transition A_{{{mu}{nu}}}.")

    return out


def _apply_two_group_transition(
    blocks: Mapping[SectorKey, Array],
    group_index: int,
    N_group: int,
    mu: str,
    nu: str,
) -> dict[SectorKey, Array]:
    """Apply one collective transition to one factor of a two-group state."""
    if group_index not in (0, 1):
        raise ValueError("Two-group squeezing requires group_index 0 or 1.")

    out: dict[SectorKey, Array] = {}
    for sector_key, psi in blocks.items():
        if not isinstance(sector_key, tuple) or len(sector_key) != 2:
            raise ValueError(
                "Two-group squeezing requires tuple sector keys (Nj1, Nj2)."
            )

        Nj1, Nj2 = (int(value) for value in sector_key)
        psi_matrix = np.asarray(psi, dtype=np.complex128).reshape(
            Nj1 + 1,
            Nj2 + 1,
        )
        Nj_group = (Nj1, Nj2)[group_index]
        n_u = N_group - Nj_group

        def flatten(values: Array) -> Array:
            return np.asarray(values, dtype=np.complex128).reshape(-1)

        def target_key(Nj_new: int) -> tuple[int, int]:
            if group_index == 0:
                return (int(Nj_new), Nj2)
            return (Nj1, int(Nj_new))

        if group_index == 0:
            ne = np.arange(Nj1 + 1, dtype=float)[:, None]
        else:
            ne = np.arange(Nj2 + 1, dtype=float)[None, :]

        if mu == "u" and nu == "u":
            _accumulate_block(out, sector_key, n_u * psi)
        elif mu == "d" and nu == "d":
            _accumulate_block(
                out,
                sector_key,
                flatten((Nj_group - ne) * psi_matrix),
            )
        elif mu == "e" and nu == "e":
            _accumulate_block(out, sector_key, flatten(ne * psi_matrix))
        elif mu == "d" and nu == "e":
            if group_index == 0:
                target = build_sector_ops(Nj1).Jm.dot(psi_matrix)
            else:
                target = (build_sector_ops(Nj2).Jm.dot(psi_matrix.T)).T
            _accumulate_block(out, sector_key, flatten(target))
        elif mu == "e" and nu == "d":
            if group_index == 0:
                target = build_sector_ops(Nj1).Jp.dot(psi_matrix)
            else:
                target = (build_sector_ops(Nj2).Jp.dot(psi_matrix.T)).T
            _accumulate_block(out, sector_key, flatten(target))
        elif mu == "u" and nu == "d":
            if Nj_group > 0:
                if group_index == 0:
                    factors = np.sqrt(
                        (N_group - Nj1 + 1.0)
                        * (Nj1 - np.arange(Nj1, dtype=float))
                    )[:, None]
                    target = factors * psi_matrix[:-1, :]
                else:
                    factors = np.sqrt(
                        (N_group - Nj2 + 1.0)
                        * (Nj2 - np.arange(Nj2, dtype=float))
                    )[None, :]
                    target = psi_matrix[:, :-1] * factors
                _accumulate_block(
                    out,
                    target_key(Nj_group - 1),
                    flatten(target),
                )
        elif mu == "d" and nu == "u":
            if n_u > 0:
                if group_index == 0:
                    target = np.zeros(
                        (Nj1 + 2, Nj2 + 1),
                        dtype=np.complex128,
                    )
                    factors = np.sqrt(
                        n_u
                        * (
                            Nj1
                            - np.arange(Nj1 + 1, dtype=float)
                            + 1.0
                        )
                    )[:, None]
                    target[:-1, :] = factors * psi_matrix
                else:
                    target = np.zeros(
                        (Nj1 + 1, Nj2 + 2),
                        dtype=np.complex128,
                    )
                    factors = np.sqrt(
                        n_u
                        * (
                            Nj2
                            - np.arange(Nj2 + 1, dtype=float)
                            + 1.0
                        )
                    )[None, :]
                    target[:, :-1] = psi_matrix * factors
                _accumulate_block(
                    out,
                    target_key(Nj_group + 1),
                    flatten(target),
                )
        elif mu == "u" and nu == "e":
            if Nj_group > 0:
                if group_index == 0:
                    factors = np.sqrt(
                        (N_group - Nj1 + 1.0)
                        * np.arange(1, Nj1 + 1, dtype=float)
                    )[:, None]
                    target = factors * psi_matrix[1:, :]
                else:
                    factors = np.sqrt(
                        (N_group - Nj2 + 1.0)
                        * np.arange(1, Nj2 + 1, dtype=float)
                    )[None, :]
                    target = psi_matrix[:, 1:] * factors
                _accumulate_block(
                    out,
                    target_key(Nj_group - 1),
                    flatten(target),
                )
        elif mu == "e" and nu == "u":
            if n_u > 0:
                if group_index == 0:
                    target = np.zeros(
                        (Nj1 + 2, Nj2 + 1),
                        dtype=np.complex128,
                    )
                    factors = np.sqrt(
                        n_u * (np.arange(Nj1 + 1, dtype=float) + 1.0)
                    )[:, None]
                    target[1:, :] = factors * psi_matrix
                else:
                    target = np.zeros(
                        (Nj1 + 1, Nj2 + 2),
                        dtype=np.complex128,
                    )
                    factors = np.sqrt(
                        n_u * (np.arange(Nj2 + 1, dtype=float) + 1.0)
                    )[None, :]
                    target[:, 1:] = psi_matrix * factors
                _accumulate_block(
                    out,
                    target_key(Nj_group + 1),
                    flatten(target),
                )
        else:
            raise ValueError(f"Unsupported transition A_{{{mu}{nu}}}.")

    return out


def _apply_group_one_body_operator(
    blocks: Mapping[SectorKey, Array],
    N_groups: tuple[int, ...],
    group_index: int,
    operator: Array,
) -> dict[SectorKey, Array]:
    """Apply a collective one-body operator to one physical group."""
    operator = np.asarray(operator, dtype=np.complex128)
    if operator.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 operator, got {operator.shape}.")
    if not 0 <= group_index < len(N_groups):
        raise ValueError("group_index is outside N_groups.")

    labels = ("u", "d", "e")
    out: dict[SectorKey, Array] = {}
    for row, mu in enumerate(labels):
        for column, nu in enumerate(labels):
            coefficient = operator[row, column]
            if abs(coefficient) <= 1e-15:
                continue

            if len(N_groups) == 1:
                transitioned = _apply_single_group_transition(
                    blocks,
                    N_groups[0],
                    mu,
                    nu,
                )
            elif len(N_groups) == 2:
                transitioned = _apply_two_group_transition(
                    blocks,
                    group_index,
                    N_groups[group_index],
                    mu,
                    nu,
                )
            else:
                raise ValueError(
                    "MCWF squeezing moments currently support one or two groups."
                )

            for sector_key, contribution in transitioned.items():
                _accumulate_block(
                    out,
                    sector_key,
                    coefficient * contribution,
                )

    return out


# -----------------------------------------------------------------------------
# Dressed bases and local fluctuation operators
# -----------------------------------------------------------------------------


def _build_dressed_state_series(
    zero_state_groups: Array,
    one_state_groups: Array,
    theta_J_groups: Array,
    phi_J_groups: Array,
    theta_S_groups: Array,
    phi_S_groups: Array,
) -> tuple[Array, Array, Array]:
    """Construct shared |c>, |j>, and |s> series for all groups and times."""
    state_shape = zero_state_groups.shape
    if len(state_shape) != 3 or state_shape[-1] != 3:
        raise ValueError(
            "S-basis state arrays must have shape (num_groups, num_times, 3)."
        )
    if one_state_groups.shape != state_shape:
        raise ValueError(
            "S-basis zero and one states must have matching shapes."
        )
    angle_shape = state_shape[:-1]
    if any(
        angles.shape != angle_shape
        for angles in (
            theta_J_groups,
            phi_J_groups,
            theta_S_groups,
            phi_S_groups,
        )
    ):
        raise ValueError(
            "J and S angle arrays must match the group and time axes of the "
            "S basis."
        )

    phase_J = np.exp(-1j * phi_J_groups)
    dressed_j = np.stack(
        (
            np.zeros(angle_shape, dtype=np.complex128),
            -np.sin(theta_J_groups / 2.0),
            phase_J * np.cos(theta_J_groups / 2.0),
        ),
        axis=-1,
    )
    phase_S = np.exp(-1j * phi_S_groups)[..., None]

    mean_state = (
        np.cos(theta_S_groups / 2.0)[..., None] * zero_state_groups
        + phase_S
        * np.sin(theta_S_groups / 2.0)[..., None]
        * one_state_groups
    )
    s_state = (
        -np.sin(theta_S_groups / 2.0)[..., None] * zero_state_groups
        + phase_S
        * np.cos(theta_S_groups / 2.0)[..., None]
        * one_state_groups
    )

    mean_state /= np.linalg.norm(mean_state, axis=-1, keepdims=True)
    dressed_j /= np.linalg.norm(dressed_j, axis=-1, keepdims=True)
    s_state /= np.linalg.norm(s_state, axis=-1, keepdims=True)
    return mean_state, dressed_j, s_state


def _build_fluctuation_operator_series(
    dressed_state_groups: tuple[Array, Array, Array],
) -> tuple[Array, Array]:
    """Construct shared fluctuation operators and |c><c| projectors."""
    if len(dressed_state_groups) != 3:
        raise ValueError("Dressed states must contain |c>, |j>, and |s>.")
    mean_state, dressed_j, s_state = dressed_state_groups
    if mean_state.shape != dressed_j.shape or mean_state.shape != s_state.shape:
        raise ValueError("Dressed-state series must have matching shapes.")
    if len(mean_state.shape) != 3 or mean_state.shape[-1] != 3:
        raise ValueError(
            "Dressed-state arrays must have shape (num_groups, num_times, 3)."
        )

    c_j = mean_state[..., :, None] * dressed_j.conjugate()[..., None, :]
    j_c = dressed_j[..., :, None] * mean_state.conjugate()[..., None, :]
    c_s = mean_state[..., :, None] * s_state.conjugate()[..., None, :]
    s_c = s_state[..., :, None] * mean_state.conjugate()[..., None, :]

    O_1 = 0.5 * (c_j + j_c)
    O_2 = (c_j - j_c) / (2.0j)
    O_3 = 0.5 * (c_s + s_c)
    O_4 = (c_s - s_c) / (2.0j)
    fluctuation_operators = np.stack((O_1, O_2, O_3, O_4), axis=2)
    N_c_operators = (
        mean_state[..., :, None] * mean_state.conjugate()[..., None, :]
    )
    return fluctuation_operators, N_c_operators


# -----------------------------------------------------------------------------
# Per-trajectory raw moment extraction
# -----------------------------------------------------------------------------


def _validate_trajectory_grid(
    trajectory: TrajectoryResult,
    t_ref: Array,
) -> None:
    """Require every trajectory to use the shared saved-time grid."""
    t = np.asarray(
        [snapshot.time for snapshot in trajectory.snapshots],
        dtype=float,
    )
    if t.shape != t_ref.shape or not np.allclose(
        t,
        t_ref,
        atol=1e-12,
        rtol=0.0,
    ):
        raise ValueError(
            "All trajectories must share ensemble.metadata.t_eval."
        )


def _trajectory_fluctuation_moments_worker(
    args: tuple[
        TrajectoryResult,
        Array,
        Array,
        Array,
        tuple[int, ...],
        float,
    ],
) -> tuple[Array, Array, Array]:
    """Extract group-local fluctuation moments for one trajectory."""
    (
        trajectory,
        t_ref,
        fluctuation_operator_groups,
        N_c_operator_groups,
        N_groups,
        tol,
    ) = args
    _validate_trajectory_grid(trajectory, t_ref)

    group_count = len(N_groups)
    num_times = t_ref.size
    mean_operators = np.zeros(
        (group_count, num_times, 4),
        dtype=float,
    )
    second_moments = np.zeros(
        (group_count, num_times, 4, 4),
        dtype=float,
    )
    N_c = np.zeros((group_count, num_times), dtype=float)

    for time_index, snapshot in enumerate(trajectory.snapshots):
        norm2 = total_norm2(snapshot.sector_blocks)
        if norm2 <= tol:
            continue

        for group_index in range(group_count):
            operators = fluctuation_operator_groups[group_index, time_index]
            applied = [
                _apply_group_one_body_operator(
                    snapshot.sector_blocks,
                    N_groups,
                    group_index,
                    operator,
                )
                for operator in operators
            ]

            mean_operators[group_index, time_index] = [
                np.real(
                    _block_inner_product(
                        snapshot.sector_blocks,
                        applied_operator,
                    )
                    / norm2
                )
                for applied_operator in applied
            ]
            for a in range(4):
                for b in range(4):
                    second_moments[
                        group_index,
                        time_index,
                        a,
                        b,
                    ] = np.real(
                        _block_inner_product(applied[a], applied[b])
                        / norm2
                    )

            N_c_applied = _apply_group_one_body_operator(
                snapshot.sector_blocks,
                N_groups,
                group_index,
                N_c_operator_groups[group_index, time_index],
            )
            N_c[group_index, time_index] = float(
                np.real(
                    _block_inner_product(
                        snapshot.sector_blocks,
                        N_c_applied,
                    )
                    / norm2
                )
            )

    return mean_operators, second_moments, N_c


# -----------------------------------------------------------------------------
# Ensemble-level group moment construction
# -----------------------------------------------------------------------------


def _validate_inputs(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol: float,
) -> tuple[
    Array,
    tuple[int, ...],
    Array,
    Array,
    Array,
    Array,
    Array,
    Array,
    Array,
]:
    """Validate and return the shared group-resolved reference data."""
    if not ensemble.trajectories:
        raise ValueError("TrajectoryEnsemble is empty.")
    if ensemble.metadata is None:
        raise ValueError(
            "TrajectoryEnsemble.metadata is required for squeezing moments."
        )

    N_groups = tuple(int(value) for value in ensemble.metadata.Ni)
    group_count = len(N_groups)
    if group_count not in (1, 2):
        raise ValueError(
            "MCWF squeezing moments currently support one or two groups."
        )

    if (
        j_moments.theta is None
        or j_moments.phi is None
        or j_moments.N_e is None
        or j_moments.N_j is None
    ):
        raise ValueError(
            "Squeezing moments require group-resolved J angles, N_e, and N_j."
        )
    if any(
        len(groups) != group_count
        for groups in (
            j_moments.theta,
            j_moments.phi,
            j_moments.N_e,
            j_moments.N_j,
        )
    ):
        raise ValueError(
            "Squeezing moments require J fields matching metadata.Ni."
        )

    t_ref = np.asarray(ensemble.metadata.t_eval, dtype=float)
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
    if theta_J_groups.shape != (group_count, t_ref.size) or (
        phi_J_groups.shape != (group_count, t_ref.size)
    ):
        raise ValueError("Each group J-angle series must match t_eval.")

    if s_moments.theta is None or s_moments.phi is None:
        raise ValueError(
            "Squeezing moments require group-resolved S angles."
        )
    if len(s_moments.theta) != group_count or (
        len(s_moments.phi) != group_count
    ):
        raise ValueError(
            "Squeezing moments require S-angle fields matching metadata.Ni."
        )
    if np.asarray(s_moments.t).shape != t_ref.shape or not np.allclose(
        s_moments.t,
        t_ref,
        atol=1e-12,
        rtol=0.0,
    ):
        raise ValueError(
            "s_moments and ensemble.metadata.t_eval must share a time grid."
        )

    theta_S_groups = np.asarray(s_moments.theta, dtype=float)
    phi_S_groups = np.asarray(s_moments.phi, dtype=float)
    if theta_S_groups.shape != (group_count, t_ref.size) or (
        phi_S_groups.shape != (group_count, t_ref.size)
    ):
        raise ValueError("Each group S-angle series must match t_eval.")

    if s_moments.basis is None:
        raise ValueError("Squeezing moments require the effective S basis.")
    zero_state_groups = np.asarray(
        s_moments.basis.zero_state,
        dtype=np.complex128,
    )
    one_state_groups = np.asarray(
        s_moments.basis.one_state,
        dtype=np.complex128,
    )
    expected_basis_shape = (group_count, t_ref.size, 3)
    if zero_state_groups.shape != expected_basis_shape or (
        one_state_groups.shape != expected_basis_shape
    ):
        raise ValueError(
            "S-basis state series must have shape "
            "(num_groups, len(t_eval), 3)."
        )

    excited_fraction_active = np.full(
        (group_count, t_ref.size),
        np.nan,
        dtype=float,
    )
    for group_index, (N_e, N_j) in enumerate(
        zip(j_moments.N_e, j_moments.N_j)
    ):
        N_e = np.asarray(N_e, dtype=float)
        N_j = np.asarray(N_j, dtype=float)
        if N_e.shape != t_ref.shape or N_j.shape != t_ref.shape:
            raise ValueError("Each group N_e and N_j series must match t_eval.")
        valid = N_j > tol
        excited_fraction_active[group_index, valid] = (
            N_e[valid] / N_j[valid]
        )

    return (
        t_ref,
        N_groups,
        zero_state_groups,
        one_state_groups,
        theta_J_groups,
        phi_J_groups,
        theta_S_groups,
        phi_S_groups,
        excited_fraction_active,
    )


def compute_mcwf_squeezing_moments(
    ensemble: TrajectoryEnsemble,
    j_moments: JMomentSeries,
    s_moments: SMomentSeries,
    *,
    tol: float = 1e-12,
    n_processes: int | None = None,
) -> tuple[SqueezingMomentSeries, ...]:
    """
    Compute ensemble-averaged, group-local moments needed for squeezing.

    The supplied J and S moments define the ensemble-averaged dressed bases.
    Raw fluctuation-operator moments are then extracted from each trajectory
    and averaged. This function does not construct covariance matrices,
    covariance eigenvalues, or squeezing parameters.
    """
    if tol <= 0.0:
        raise ValueError("tol must be positive.")

    (
        t_ref,
        N_groups,
        zero_state_groups,
        one_state_groups,
        theta_J_groups,
        phi_J_groups,
        theta_S_groups,
        phi_S_groups,
        excited_fraction_active,
    ) = _validate_inputs(
        ensemble,
        j_moments,
        s_moments,
        tol=tol,
    )

    dressed_state_groups = _build_dressed_state_series(
        zero_state_groups,
        one_state_groups,
        theta_J_groups,
        phi_J_groups,
        theta_S_groups,
        phi_S_groups,
    )
    (
        fluctuation_operator_groups,
        N_c_operator_groups,
    ) = _build_fluctuation_operator_series(dressed_state_groups)

    fluctuation_samples = map_with_optional_pool(
        _trajectory_fluctuation_moments_worker,
        [
            (
                trajectory,
                t_ref,
                fluctuation_operator_groups,
                N_c_operator_groups,
                N_groups,
                tol,
            )
            for trajectory in ensemble.trajectories
        ],
        n_processes=n_processes,
        progress_desc="compute_mcwf_squeezing_fluctuation_moments",
    )
    mean_operators = np.mean(
        np.asarray(
            [sample[0] for sample in fluctuation_samples],
            dtype=float,
        ),
        axis=0,
    )
    second_moments = np.mean(
        np.asarray(
            [sample[1] for sample in fluctuation_samples],
            dtype=float,
        ),
        axis=0,
    )
    N_c = np.mean(
        np.asarray(
            [sample[2] for sample in fluctuation_samples],
            dtype=float,
        ),
        axis=0,
    )
    N_c[(N_c < 0.0) & (np.abs(N_c) < 1e-10)] = 0.0

    return tuple(
        SqueezingMomentSeries(
            t=t_ref,
            N_group=N_group,
            mean_operators=mean_operators[group_index],
            second_moments=second_moments[group_index],
            N_c=N_c[group_index],
            excited_fraction_active=excited_fraction_active[group_index],
            theta_J=theta_J_groups[group_index],
            phi_J=phi_J_groups[group_index],
            theta_S=theta_S_groups[group_index],
            phi_S=phi_S_groups[group_index],
        )
        for group_index, N_group in enumerate(N_groups)
    )


__all__ = ["compute_mcwf_squeezing_moments"]
