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


def _dressed_basis_states(
    theta_J: float,
    phi_J: float,
) -> tuple[Array, Array]:
    """Construct the local dressed states |1> and |j>."""
    cosine = np.cos(theta_J / 2.0)
    sine = np.sin(theta_J / 2.0)
    phase = np.exp(-1j * phi_J)

    dressed_one = np.array(
        [0.0, cosine, phase * sine],
        dtype=np.complex128,
    )
    dressed_j = np.array(
        [0.0, -sine, phase * cosine],
        dtype=np.complex128,
    )
    return dressed_one, dressed_j


def _s_bloch_operators(theta_J: float, phi_J: float) -> tuple[Array, ...]:
    """Construct the three dressed S-Bloch operators."""
    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    dressed_one, _ = _dressed_basis_states(theta_J, phi_J)
    projector_u = np.outer(ket_u, ket_u.conjugate())
    projector_one = np.outer(dressed_one, dressed_one.conjugate())

    S_x = 0.5 * (
        np.outer(dressed_one, ket_u.conjugate())
        + np.outer(ket_u, dressed_one.conjugate())
    )
    S_y = (
        np.outer(dressed_one, ket_u.conjugate())
        - np.outer(ket_u, dressed_one.conjugate())
    ) / (2.0j)
    S_z = 0.5 * (projector_u - projector_one)
    return S_x, S_y, S_z


def _mean_and_fluctuation_basis(
    theta_J: float,
    phi_J: float,
    theta_S: float,
    phi_S: float,
) -> tuple[Array, Array, Array]:
    """Construct the local dressed states |c>, |j>, and |s>."""
    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    dressed_one, dressed_j = _dressed_basis_states(theta_J, phi_J)
    phase_S = np.exp(-1j * phi_S)

    mean_state = (
        np.cos(theta_S / 2.0) * ket_u
        + phase_S * np.sin(theta_S / 2.0) * dressed_one
    )
    s_state = (
        -np.sin(theta_S / 2.0) * ket_u
        + phase_S * np.cos(theta_S / 2.0) * dressed_one
    )

    mean_state /= np.linalg.norm(mean_state)
    dressed_j /= np.linalg.norm(dressed_j)
    s_state /= np.linalg.norm(s_state)
    return mean_state, dressed_j, s_state


def _fluctuation_operators(
    mean_state: Array,
    dressed_j: Array,
    s_state: Array,
) -> tuple[Array, ...]:
    """Construct the four local generalized-squeezing quadratures."""
    O_1 = 0.5 * (
        np.outer(mean_state, dressed_j.conjugate())
        + np.outer(dressed_j, mean_state.conjugate())
    )
    O_2 = (
        np.outer(mean_state, dressed_j.conjugate())
        - np.outer(dressed_j, mean_state.conjugate())
    ) / (2.0j)
    O_3 = 0.5 * (
        np.outer(mean_state, s_state.conjugate())
        + np.outer(s_state, mean_state.conjugate())
    )
    O_4 = (
        np.outer(mean_state, s_state.conjugate())
        - np.outer(s_state, mean_state.conjugate())
    ) / (2.0j)
    return O_1, O_2, O_3, O_4


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
        theta_J_groups,
        phi_J_groups,
        theta_S_groups,
        phi_S_groups,
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
            basis = _mean_and_fluctuation_basis(
                theta_J_groups[group_index, time_index],
                phi_J_groups[group_index, time_index],
                theta_S_groups[group_index, time_index],
                phi_S_groups[group_index, time_index],
            )
            operators = _fluctuation_operators(*basis)
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

            mean_state = basis[0]
            N_c_operator = np.outer(
                mean_state,
                mean_state.conjugate(),
            )
            N_c_applied = _apply_group_one_body_operator(
                snapshot.sector_blocks,
                N_groups,
                group_index,
                N_c_operator,
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
) -> tuple[Array, tuple[int, ...], Array, Array, Array, Array, Array]:
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

    fluctuation_samples = map_with_optional_pool(
        _trajectory_fluctuation_moments_worker,
        [
            (
                trajectory,
                t_ref,
                theta_J_groups,
                phi_J_groups,
                theta_S_groups,
                phi_S_groups,
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
