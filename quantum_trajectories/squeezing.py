from __future__ import annotations

import multiprocessing as mp
import time
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import numpy as np

from common.utils import active_manifold_angles
from common.utils import phase_change_times, phase1_ss_angles_for_nj
from common.parser import ObservableSeries
from quantum_trajectories.aggregator import expected_collective_components, trajectory_observables
from quantum_trajectories.operator_helpers import build_sector_ops
from quantum_trajectories.parser import Array, TrajectoryEnsemble, TrajectoryResult
from quantum_trajectories.state_helpers import total_norm2


def _block_inner_product(left: Mapping[int, Array], right: Mapping[int, Array]) -> complex:
    """
    Inner product between two reduced-basis states written as {N_J: psi_NJ}.

    Missing sectors are interpreted as zero, so only sectors present in both
    operands contribute.
    """
    shared_sectors = set(left).intersection(right)
    return sum(np.vdot(left[Nj], right[Nj]) for Nj in shared_sectors)


def _interp_series(t_src: Array, y_src: Array, t_ref: Array) -> Array:
    """
    Linear interpolation onto a common reference grid.

    This mirrors the ensemble observable strategy used elsewhere in the code:
    each trajectory may have slightly different saved times, so we interpolate
    scalar moment series before averaging across trajectories.
    """
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    t_ref = np.asarray(t_ref, dtype=float)

    if t_src.ndim != 1 or y_src.ndim != 1:
        raise ValueError("t_src and y_src must be 1D arrays.")
    if len(t_src) != len(y_src):
        raise ValueError("t_src and y_src must have the same length.")

    t_unique, idx = np.unique(t_src, return_index=True)
    y_unique = y_src[idx]
    return np.interp(t_ref, t_unique, y_unique)


def _map_with_optional_pool(
    worker,
    items: Iterable,
    *,
    n_processes: Optional[int],
    progress_desc: str,
):
    """
    Run a trajectory-wise worker either serially or through multiprocessing.

    A tqdm progress bar is shown in both cases so long post-processing stages
    stay visible when profiling ensemble squeezing.
    """
    from tqdm.auto import tqdm

    items = list(items)
    if n_processes is None or n_processes == 1:
        return [worker(item) for item in tqdm(items, desc=progress_desc)]

    if n_processes == -1:
        n_processes = mp.cpu_count()
    if n_processes <= 0:
        raise ValueError("n_processes must be None, 1, -1, or a positive integer.")

    ctx = mp.get_context()
    with ctx.Pool(processes=n_processes) as pool:
        return list(tqdm(pool.imap(worker, items), total=len(items), desc=progress_desc))


def _trajectory_observables_worker(args: tuple[TrajectoryResult, float]) -> ObservableSeries:
    """
    Multiprocessing-safe wrapper for trajectory_observables(...).

    This is used by the ensemble squeezing pipeline when it needs per-trajectory
    Jx, Jy, Jz, and N_e series and no already-averaged observable container was
    supplied by the caller.
    """
    traj, tol = args
    return trajectory_observables(traj, tol=tol)


def _accumulate_block(out: Dict[int, Array], Nj: int, contribution: Array) -> None:
    """
    Add one sector contribution into an output block dictionary in-place.
    """
    if contribution.size == 0:
        return
    if Nj in out:
        out[Nj] = out[Nj] + contribution
    else:
        out[Nj] = np.asarray(contribution, dtype=np.complex128).copy()


def _apply_collective_transition(
    blocks: Mapping[int, Array],
    N: int,
    mu: str,
    nu: str,
) -> Dict[int, Array]:
    """
    Apply A_{mu,nu} = sum_i |mu_i><nu_i| to the reduced MCWF state.

    The reduced basis stores each sector as a symmetric Dicke vector labeled by
    (N_J, n_e), with:
        N_u = N - N_J,
        N_d = N_J - n_e,
        N_e = n_e.

    In this representation the collective one-body transitions are equivalent
    to bosonic occupation transfers between u, d, and e. The formulas below are
    the exact symmetric-basis action of those transfers.
    """
    out: Dict[int, Array] = {}

    for Nj, psi in blocks.items():
        psi = np.asarray(psi, dtype=np.complex128)
        ne = np.arange(Nj + 1, dtype=float)
        n_u = N - Nj

        if mu == "u" and nu == "u":
            _accumulate_block(out, Nj, n_u * psi)
        elif mu == "d" and nu == "d":
            _accumulate_block(out, Nj, (Nj - ne) * psi)
        elif mu == "e" and nu == "e":
            _accumulate_block(out, Nj, ne * psi)
        elif mu == "d" and nu == "e":
            _accumulate_block(out, Nj, build_sector_ops(Nj).J_minus.dot(psi))
        elif mu == "e" and nu == "d":
            _accumulate_block(out, Nj, build_sector_ops(Nj).J_plus.dot(psi))
        elif mu == "u" and nu == "d":
            if Nj > 0:
                target = np.sqrt((N - Nj + 1.0) * (Nj - ne[:-1])) * psi[:-1]
                _accumulate_block(out, Nj - 1, target)
        elif mu == "d" and nu == "u":
            if n_u > 0:
                target = np.zeros(Nj + 2, dtype=np.complex128)
                target[:-1] = np.sqrt(n_u * (Nj - ne + 1.0)) * psi
                _accumulate_block(out, Nj + 1, target)
        elif mu == "u" and nu == "e":
            if Nj > 0:
                target = np.sqrt((N - Nj + 1.0) * ne[1:]) * psi[1:]
                _accumulate_block(out, Nj - 1, target)
        elif mu == "e" and nu == "u":
            if n_u > 0:
                target = np.zeros(Nj + 2, dtype=np.complex128)
                target[1:] = np.sqrt(n_u * (ne + 1.0)) * psi
                _accumulate_block(out, Nj + 1, target)
        else:
            raise ValueError(f"Unsupported transition A_{{{mu}{nu}}}.")

    return out


def _apply_collective_one_body_operator(
    blocks: Mapping[int, Array],
    N: int,
    operator: Array,
) -> Dict[int, Array]:
    """
    Apply a collective one-body operator O = sum_{mu,nu} o_{mu,nu} A_{mu,nu}.

    The 3x3 single-particle matrix uses the basis ordering (u, d, e).
    """
    operator = np.asarray(operator, dtype=np.complex128)
    if operator.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 operator, got shape {operator.shape}.")

    labels = ("u", "d", "e")
    out: Dict[int, Array] = {}
    for row, mu in enumerate(labels):
        for col, nu in enumerate(labels):
            coeff = operator[row, col]
            if abs(coeff) <= 1e-15:
                continue
            transitioned = _apply_collective_transition(blocks, N, mu, nu)
            for Nj, psi in transitioned.items():
                _accumulate_block(out, Nj, coeff * psi)
    return out


def _expect_collective_one_body_operator(
    blocks: Mapping[int, Array],
    N: int,
    operator: Array,
    *,
    norm2: Optional[float] = None,
) -> complex:
    """
    Return the normalized expectation value <O> for a one-body collective O.
    """
    if norm2 is None:
        norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0j
    applied = _apply_collective_one_body_operator(blocks, N, operator)
    return _block_inner_product(blocks, applied) / norm2


def _single_snapshot_j_angles(
    blocks: Mapping[int, Array],
    *,
    tol: float,
) -> tuple[float, float]:
    """
    Reuse the existing active-manifold normalization to extract theta_J, phi_J.
    """
    jx, jy, jz, ne = expected_collective_components(blocks)
    theta, phi, _, _, _, _ = active_manifold_angles(
        np.asarray([jx]),
        np.asarray([jy]),
        np.asarray([jz]),
        np.asarray([ne]),
        tol=tol,
    )
    return float(theta[0]), float(phi[0])


def _j_angles_from_collective_components(
    jx: float,
    jy: float,
    jz: float,
    ne: float,
    *,
    tol: float,
) -> tuple[float, float]:
    """
    Compute theta_J, phi_J from already-averaged collective moments.

    This is the ensemble counterpart of _single_snapshot_j_angles(...): the
    same active-manifold normalization is used, but the inputs are the
    ensemble-averaged moments rather than one pure trajectory snapshot.
    """
    theta, phi, _, _, _, _ = active_manifold_angles(
        np.asarray([jx]),
        np.asarray([jy]),
        np.asarray([jz]),
        np.asarray([ne]),
        tol=tol,
    )
    return float(theta[0]), float(phi[0])


def _dressed_basis_states(theta_j: float, phi_j: float) -> tuple[Array, Array]:
    """
    Construct |1> and |j> from the instruction-sheet formulas.
    """
    c = np.cos(theta_j / 2.0)
    s = np.sin(theta_j / 2.0)
    phase = np.exp(-1j * phi_j)

    dressed_one = np.array([0.0, c, phase * s], dtype=np.complex128)
    dressed_j = np.array([0.0, -s, phase * c], dtype=np.complex128)

    return dressed_one, dressed_j


def _s_bloch_angles(
    blocks: Mapping[int, Array],
    N: int,
    dressed_one: Array,
    *,
    tol: float,
    norm2: Optional[float] = None,
) -> tuple[float, float, float]:
    """
    Compute theta_S, phi_S and <N_{u1}> using the dressed |u>,|1> pseudospin.
    """
    if norm2 is None:
        norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0, 0.0, 0.0

    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)

    projector_u = np.outer(ket_u, ket_u.conjugate())
    projector_one = np.outer(dressed_one, dressed_one.conjugate())

    sx_op = 0.5 * (
        np.outer(dressed_one, ket_u.conjugate())
        + np.outer(ket_u, dressed_one.conjugate())
    )
    sy_op = (
        np.outer(dressed_one, ket_u.conjugate())
        - np.outer(ket_u, dressed_one.conjugate())
    ) / (2.0j)
    sz_op = 0.5 * (projector_u - projector_one)
    n_u1_op = projector_u + projector_one

    n_u1 = float(np.real(_expect_collective_one_body_operator(blocks, N, n_u1_op, norm2=norm2)))
    if n_u1 <= tol:
        return 0.0, 0.0, 0.0

    sx = float(np.real(_expect_collective_one_body_operator(blocks, N, sx_op, norm2=norm2)))
    sy = float(np.real(_expect_collective_one_body_operator(blocks, N, sy_op, norm2=norm2)))
    sz = float(np.real(_expect_collective_one_body_operator(blocks, N, sz_op, norm2=norm2)))

    sx = 2.0 * sx / n_u1
    sy = 2.0 * sy / n_u1
    sz = 2.0 * sz / n_u1
    sz = float(np.clip(sz, -1.0, 1.0))

    theta_s = float(np.arccos(sz))
    phi_s = float(np.arctan2(sy, sx))
    if np.hypot(sx, sy) < tol:
        phi_s = 0.0

    return theta_s, phi_s, n_u1


def _mean_and_fluctuation_basis(
    theta_j: float,
    phi_j: float,
    theta_s: float,
    phi_s: float,
) -> tuple[Array, Array, Array]:
    """
    Construct |c>, |j>, |s> as defined in the squeezing instructions.
    """
    dressed_one, dressed_j = _dressed_basis_states(theta_j, phi_j)
    phase_s = np.exp(-1j * phi_s)

    mean_state = (
        np.cos(theta_s / 2.0) * np.array([1.0, 0.0, 0.0], dtype=np.complex128)
        + phase_s * np.sin(theta_s / 2.0) * dressed_one
    )
    s_state = (
        -np.sin(theta_s / 2.0) * np.array([1.0, 0.0, 0.0], dtype=np.complex128)
        + phase_s * np.cos(theta_s / 2.0) * dressed_one
    )

    # Roundoff can accumulate over long runs, so explicitly normalize the
    # basis states before using them to build fluctuation operators.
    mean_state /= np.linalg.norm(mean_state)
    dressed_j /= np.linalg.norm(dressed_j)
    s_state /= np.linalg.norm(s_state)

    return mean_state, dressed_j, s_state


def _trajectory_s_bloch_series_worker(args: tuple) -> tuple[Array, Array, Array, Array]:
    """
    Worker for trajectory-wise extraction of S-Bloch moment series.

    For a fixed ensemble-wide dressed |1>(t), this returns the interpolated
    per-trajectory series needed to build the ensemble-averaged S direction:
        <N_{u1}>(t), <S_x>(t), <S_y>(t), <S_z>(t).
    """
    traj, t_ref, theta_j_series, phi_j_series, tol = args
    t_src = np.asarray([snap.time for snap in traj.snapshots], dtype=float)
    if len(t_src) != len(t_ref) or not np.allclose(t_src, t_ref, atol=1e-12, rtol=0.0):
        raise ValueError(
            "All trajectories must share the same t_eval snapshot grid. "
            "Run the ensemble through the common num_snapshots API."
        )
    n_u1_src = np.zeros(len(t_src), dtype=float)
    sx_src = np.zeros(len(t_src), dtype=float)
    sy_src = np.zeros(len(t_src), dtype=float)
    sz_src = np.zeros(len(t_src), dtype=float)

    ket_u = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
    for j, snap in enumerate(traj.snapshots):
        norm2 = total_norm2(snap.sector_blocks)
        if norm2 <= tol:
            continue

        dressed_one, _ = _dressed_basis_states(theta_j_series[j], phi_j_series[j])
        projector_u = np.outer(ket_u, ket_u.conjugate())
        projector_one = np.outer(dressed_one, dressed_one.conjugate())
        sx_op = 0.5 * (
            np.outer(dressed_one, ket_u.conjugate())
            + np.outer(ket_u, dressed_one.conjugate())
        )
        sy_op = (
            np.outer(dressed_one, ket_u.conjugate())
            - np.outer(ket_u, dressed_one.conjugate())
        ) / (2.0j)
        sz_op = 0.5 * (projector_u - projector_one)
        n_u1_op = projector_u + projector_one

        n_u1_src[j] = float(
            np.real(
                _expect_collective_one_body_operator(
                    snap.sector_blocks,
                    traj.N,
                    n_u1_op,
                    norm2=norm2,
                )
            )
        )
        sx_src[j] = float(
            np.real(
                _expect_collective_one_body_operator(
                    snap.sector_blocks,
                    traj.N,
                    sx_op,
                    norm2=norm2,
                )
            )
        )
        sy_src[j] = float(
            np.real(
                _expect_collective_one_body_operator(
                    snap.sector_blocks,
                    traj.N,
                    sy_op,
                    norm2=norm2,
                )
            )
        )
        sz_src[j] = float(
            np.real(
                _expect_collective_one_body_operator(
                    snap.sector_blocks,
                    traj.N,
                    sz_op,
                    norm2=norm2,
                )
            )
        )

    return (n_u1_src, sx_src, sy_src, sz_src)


def _trajectory_covariance_series_worker(args: tuple) -> tuple[Array, Array, Array]:
    """
    Worker for trajectory-wise extraction of covariance-building moment series.

    For fixed ensemble-wide dressed directions |c>(t), |j>(t), |s>(t), this
    returns the interpolated per-trajectory series for
        <O_a>(t), <O_a O_b>(t), <N_c>(t).
    """
    traj, t_ref, theta_j_series, phi_j_series, theta_s_series, phi_s_series, tol = args
    t_src = np.asarray([snap.time for snap in traj.snapshots], dtype=float)
    if len(t_src) != len(t_ref) or not np.allclose(t_src, t_ref, atol=1e-12, rtol=0.0):
        raise ValueError(
            "All trajectories must share the same t_eval snapshot grid. "
            "Run the ensemble through the common num_snapshots API."
        )
    mu_src = np.zeros((len(t_src), 4), dtype=float)
    second_src = np.zeros((len(t_src), 4, 4), dtype=float)
    n_c_src = np.zeros(len(t_src), dtype=float)

    for j, snap in enumerate(traj.snapshots):
        norm2 = total_norm2(snap.sector_blocks)
        if norm2 <= tol:
            continue

        mean_state, dressed_j, s_state = _mean_and_fluctuation_basis(
            theta_j_series[j],
            phi_j_series[j],
            theta_s_series[j],
            phi_s_series[j],
        )
        o1 = 0.5 * (
            np.outer(mean_state, dressed_j.conjugate())
            + np.outer(dressed_j, mean_state.conjugate())
        )
        o2 = (
            np.outer(mean_state, dressed_j.conjugate())
            - np.outer(dressed_j, mean_state.conjugate())
        ) / (2.0j)
        o3 = 0.5 * (
            np.outer(mean_state, s_state.conjugate())
            + np.outer(s_state, mean_state.conjugate())
        )
        o4 = (
            np.outer(mean_state, s_state.conjugate())
            - np.outer(s_state, mean_state.conjugate())
        ) / (2.0j)
        operators = (o1, o2, o3, o4)

        applied = [
            _apply_collective_one_body_operator(snap.sector_blocks, traj.N, operator)
            for operator in operators
        ]
        mu_src[j] = [
            np.real(_block_inner_product(snap.sector_blocks, out) / norm2)
            for out in applied
        ]
        for a in range(4):
            for b in range(4):
                second_src[j, a, b] = np.real(
                    _block_inner_product(applied[a], applied[b]) / norm2
                )

        n_c_op = np.outer(mean_state, mean_state.conjugate())
        n_c_src[j] = float(
            np.real(
                _expect_collective_one_body_operator(
                    snap.sector_blocks,
                    traj.N,
                    n_c_op,
                    norm2=norm2,
                )
            )
        )

    return mu_src, second_src, n_c_src


def generalized_squeezing_for_trajectory(
    result: TrajectoryResult,
    *,
    tol: float = 1e-12,
) -> dict:
    """
    Compute the generalized three-level squeezing parameter xi_gen^2(t).

    This is a standalone post-processing routine: it consumes the saved MCWF
    snapshots and reconstructs the dressed directions |1>, |c>, |j>, |s> at
    each saved time. No additional simulation-time data are required.

    Returns
    -------
    dict
        Time series containing:
            t,
            xi2_gen,
            lambda_min,
            N_c,
            excited_fraction_active,
            theta_J,
            phi_J,
            theta_S,
            phi_S.
    """
    t = np.asarray([snap.time for snap in result.snapshots], dtype=float)
    xi2_gen = np.full_like(t, np.nan, dtype=float)
    lambda_min = np.full_like(t, np.nan, dtype=float)
    covariance_eigvals = np.full((len(t), 4), np.nan, dtype=float)
    n_c_expect = np.full_like(t, np.nan, dtype=float)
    excited_fraction_active = np.full_like(t, np.nan, dtype=float)
    theta_j = np.zeros_like(t, dtype=float)
    phi_j = np.zeros_like(t, dtype=float)
    theta_s = np.zeros_like(t, dtype=float)
    phi_s = np.zeros_like(t, dtype=float)

    for k, snap in enumerate(result.snapshots):
        blocks = snap.sector_blocks
        norm2 = total_norm2(blocks)
        if norm2 <= tol:
            continue

        _, _, _, ne_expect = expected_collective_components(blocks)
        n_j_expect = sum(
            Nj * float(np.vdot(psi, psi).real) for Nj, psi in blocks.items()
        ) / norm2
        if n_j_expect > tol:
            excited_fraction_active[k] = ne_expect / n_j_expect

        theta_j[k], phi_j[k] = _single_snapshot_j_angles(blocks, tol=tol)
        dressed_one, _ = _dressed_basis_states(theta_j[k], phi_j[k])
        theta_s[k], phi_s[k], _ = _s_bloch_angles(
            blocks,
            result.N,
            dressed_one,
            tol=tol,
            norm2=norm2,
        )

        mean_state, dressed_j, s_state = _mean_and_fluctuation_basis(
            theta_j[k],
            phi_j[k],
            theta_s[k],
            phi_s[k],
        )

        # Build the four single-particle fluctuation operators from the
        # instruction sheet, then lift them into collective operators through
        # O = sum_{mu,nu} o_{mu,nu} A_{mu,nu}.
        o1 = 0.5 * (
            np.outer(mean_state, dressed_j.conjugate())
            + np.outer(dressed_j, mean_state.conjugate())
        )
        o2 = (
            np.outer(mean_state, dressed_j.conjugate())
            - np.outer(dressed_j, mean_state.conjugate())
        ) / (2.0j)
        o3 = 0.5 * (
            np.outer(mean_state, s_state.conjugate())
            + np.outer(s_state, mean_state.conjugate())
        )
        o4 = (
            np.outer(mean_state, s_state.conjugate())
            - np.outer(s_state, mean_state.conjugate())
        ) / (2.0j)
        operators = (o1, o2, o3, o4)

        applied = [
            _apply_collective_one_body_operator(blocks, result.N, operator)
            for operator in operators
        ]
        mu = np.array(
            [np.real(_block_inner_product(blocks, out) / norm2) for out in applied],
            dtype=float,
        )

        covariance = np.zeros((4, 4), dtype=float)
        for a in range(4):
            for b in range(4):
                covariance[a, b] = (
                    np.real(_block_inner_product(applied[a], applied[b]) / norm2)
                    - mu[a] * mu[b]
                )

        eigvals = np.linalg.eigvalsh(covariance)
        covariance_eigvals[k, : len(eigvals)] = np.sort(np.real(eigvals))
        lambda_min[k] = float(np.min(eigvals))

        n_c_op = np.outer(mean_state, mean_state.conjugate())
        n_c = float(
            np.real(
                _expect_collective_one_body_operator(
                    blocks,
                    result.N,
                    n_c_op,
                    norm2=norm2,
                )
            )
        )
        if n_c < 0.0 and abs(n_c) < 1e-10:
            n_c = 0.0
        n_c_expect[k] = n_c

        denominator = (n_c / 2.0) ** 2
        if denominator > tol:
            xi2_gen[k] = result.N * lambda_min[k] / denominator

    return {
        "t": t,
        "xi2_gen": xi2_gen,
        "lambda_min": lambda_min,
        "covariance_eigvals": covariance_eigvals,
        "N_c": n_c_expect,
        "excited_fraction_active": excited_fraction_active,
        "theta_J": theta_j,
        "phi_J": phi_j,
        "theta_S": theta_s,
        "phi_S": phi_s,
    }



def generalized_squeezing_for_ensemble(
    ensemble: TrajectoryEnsemble,
    *,
    tol: float = 1e-12,
    n_processes: Optional[int] = None,
    verbose: bool = True,
    averaged_observables: Optional[ObservableSeries] = None,
) -> dict:
    """
    Compute xi_gen^2(t) from ensemble-averaged MCWF moments.

    Important
    ---------
    This routine does *not* compute xi_gen^2 separately for each trajectory and
    then average the resulting values. Instead, for each saved time it first
    averages the required operator moments over trajectories:

        <O_a>, <O_a O_b>, <N_c>.

    The covariance matrix is then constructed from those ensemble-averaged
    moments, which is equivalent to evaluating the squeezing on the
    ensemble-averaged density matrix

        rho(t) = (1 / N_traj) sum_r |psi_r(t)><psi_r(t)|.

    Parameters
    ----------
    n_processes
        Parallelism for the trajectory-wise moment extraction.
        - None or 1: serial execution.
        - -1: use all available CPU cores.
        - >1: use that many worker processes.
    verbose
        If True, print timers for the major stages of the squeezing
        calculation so bottlenecks are easy to identify.
    averaged_observables
        Optional precomputed ensemble observable averages on the same t_eval
        grid. If supplied and its time grid matches the ensemble snapshots
        exactly, the squeezing code reuses its Jx, Jy, Jz, N_e, and N_j series
        instead of recomputing trajectory_observables(...) for every
        trajectory.
    """
    if len(ensemble.trajectories) == 0:
        raise ValueError("Ensemble is empty.")

    reference = ensemble.trajectories[0]
    t = np.asarray([snap.time for snap in reference.snapshots], dtype=float)
    nsnaps = len(t)

    xi2_gen = np.full_like(t, np.nan, dtype=float)
    lambda_min = np.full_like(t, np.nan, dtype=float)
    covariance_eigvals = np.full((len(t), 4), np.nan, dtype=float)
    n_c_expect = np.full_like(t, np.nan, dtype=float)
    excited_fraction_active = np.full_like(t, np.nan, dtype=float)
    theta_j = np.zeros_like(t, dtype=float)
    phi_j = np.zeros_like(t, dtype=float)
    theta_s = np.zeros_like(t, dtype=float)
    phi_s = np.zeros_like(t, dtype=float)

    t0_total = time.perf_counter()

    t0 = time.perf_counter()
    if averaged_observables is not None:
        if len(averaged_observables.t) != len(t) or not np.allclose(
            averaged_observables.t,
            t,
            atol=1e-12,
            rtol=0.0,
        ):
            raise ValueError(
                "averaged_observables must be defined on the same t_eval grid "
                "as the ensemble trajectories."
            )
        jx_mean = np.asarray(averaged_observables.Jx, dtype=float)
        jy_mean = np.asarray(averaged_observables.Jy, dtype=float)
        jz_mean = np.asarray(averaged_observables.Jz, dtype=float)
        ne_mean = np.asarray(averaged_observables.N_e, dtype=float)
        nj_mean = np.asarray(averaged_observables.N_j, dtype=float)
    else:
        per_traj_obs = _map_with_optional_pool(
            _trajectory_observables_worker,
            [(traj, tol) for traj in ensemble.trajectories],
            n_processes=n_processes,
            progress_desc="xi phase 0/2: regular observables",
        )
        for obs in per_traj_obs:
            if len(obs.t) != len(t) or not np.allclose(obs.t, t, atol=1e-12, rtol=0.0):
                raise ValueError(
                    "All trajectories must share the same t_eval snapshot grid. "
                    "Run the ensemble through the common num_snapshots API."
                )
        jx_mean = np.mean(np.asarray([obs.Jx for obs in per_traj_obs], dtype=float), axis=0)
        jy_mean = np.mean(np.asarray([obs.Jy for obs in per_traj_obs], dtype=float), axis=0)
        jz_mean = np.mean(np.asarray([obs.Jz for obs in per_traj_obs], dtype=float), axis=0)
        ne_mean = np.mean(np.asarray([obs.N_e for obs in per_traj_obs], dtype=float), axis=0)
        nj_mean = np.mean(np.asarray([obs.N_j for obs in per_traj_obs], dtype=float), axis=0)
    if verbose:
        print(f"xi timing: basic observable alignment finished in {time.perf_counter() - t0:.3f} s")

    valid_nj = nj_mean > tol
    excited_fraction_active[valid_nj] = ne_mean[valid_nj] / nj_mean[valid_nj]

    t0 = time.perf_counter()
    for k in range(nsnaps):
        theta_j[k], phi_j[k] = _j_angles_from_collective_components(
            float(jx_mean[k]),
            float(jy_mean[k]),
            float(jz_mean[k]),
            float(ne_mean[k]),
            tol=tol,
        )
    if verbose:
        print(f"xi timing: theta_J/phi_J construction finished in {time.perf_counter() - t0:.3f} s")

    t0 = time.perf_counter()
    s_bloch_results = _map_with_optional_pool(
        _trajectory_s_bloch_series_worker,
        [(traj, t, theta_j, phi_j, tol) for traj in ensemble.trajectories],
        n_processes=n_processes,
        progress_desc="xi phase 1/2: S-Bloch moments",
    )
    n_u1_arr = np.asarray([item[0] for item in s_bloch_results], dtype=float)
    sx_arr = np.asarray([item[1] for item in s_bloch_results], dtype=float)
    sy_arr = np.asarray([item[2] for item in s_bloch_results], dtype=float)
    sz_arr = np.asarray([item[3] for item in s_bloch_results], dtype=float)
    if verbose:
        print(f"xi timing: trajectory S-Bloch extraction finished in {time.perf_counter() - t0:.3f} s")

    t0 = time.perf_counter()
    n_u1_mean = np.mean(n_u1_arr, axis=0)
    valid_s = n_u1_mean > tol
    sx_mean = np.zeros_like(t, dtype=float)
    sy_mean = np.zeros_like(t, dtype=float)
    sz_mean = np.zeros_like(t, dtype=float)
    sx_mean[valid_s] = 2.0 * np.mean(sx_arr, axis=0)[valid_s] / n_u1_mean[valid_s]
    sy_mean[valid_s] = 2.0 * np.mean(sy_arr, axis=0)[valid_s] / n_u1_mean[valid_s]
    sz_mean[valid_s] = 2.0 * np.mean(sz_arr, axis=0)[valid_s] / n_u1_mean[valid_s]
    sz_mean = np.clip(sz_mean, -1.0, 1.0)
    theta_s[valid_s] = np.arccos(sz_mean[valid_s])
    phi_s = np.arctan2(sy_mean, sx_mean)
    phi_s[np.hypot(sx_mean, sy_mean) < tol] = 0.0
    if verbose:
        print(f"xi timing: theta_S/phi_S construction finished in {time.perf_counter() - t0:.3f} s")

    t0 = time.perf_counter()
    covariance_results = _map_with_optional_pool(
        _trajectory_covariance_series_worker,
        [(traj, t, theta_j, phi_j, theta_s, phi_s, tol) for traj in ensemble.trajectories],
        n_processes=n_processes,
        progress_desc="xi phase 2/2: covariance moments",
    )
    mu_arr = np.asarray([item[0] for item in covariance_results], dtype=float)
    second_arr = np.asarray([item[1] for item in covariance_results], dtype=float)
    n_c_arr = np.asarray([item[2] for item in covariance_results], dtype=float)
    if verbose:
        print(f"xi timing: trajectory covariance extraction finished in {time.perf_counter() - t0:.3f} s")

    t0 = time.perf_counter()
    for k in range(nsnaps):
        if not valid_s[k]:
            continue

        mu_mean = np.mean(mu_arr[:, k, :], axis=0)
        second_moment_mean = np.mean(second_arr[:, k, :, :], axis=0)
        covariance = second_moment_mean - np.outer(mu_mean, mu_mean)

        eigvals = np.linalg.eigvalsh(covariance)
        covariance_eigvals[k, : len(eigvals)] = np.sort(np.real(eigvals))
        lambda_min[k] = float(np.min(eigvals))

        n_c = float(np.mean(n_c_arr[:, k]))
        if n_c < 0.0 and abs(n_c) < 1e-10:
            n_c = 0.0
        n_c_expect[k] = n_c

        denominator = (n_c / 2.0) ** 2
        if denominator > tol:
            xi2_gen[k] = reference.N * lambda_min[k] / denominator
    if verbose:
        print(f"xi timing: final covariance assembly finished in {time.perf_counter() - t0:.3f} s")
        print(f"xi timing: total ensemble squeezing finished in {time.perf_counter() - t0_total:.3f} s")

    return {
        "t": t,
        "xi2_gen": xi2_gen,
        "lambda_min": lambda_min,
        "covariance_eigvals": covariance_eigvals,
        "N_c": n_c_expect,
        "excited_fraction_active": excited_fraction_active,
        "theta_J": theta_j,
        "phi_J": phi_j,
        "theta_S": theta_s,
        "phi_S": phi_s,
    }


def plot_generalized_xi(
    result: TrajectoryResult | TrajectoryEnsemble,
    *,
    tol: float = 1e-12,
    n_processes: Optional[int] = None,
    verbose: bool = True,
    averaged_observables: Optional[ObservableSeries] = None,
    ax: Optional[Any] = None,
    output_path: Optional[str | Path] = None,
) -> tuple[dict, Any, Any]:
    """
    Convenience plot for the generalized squeezing parameter in dB.

    The plotted quantity is

        xi_gen,dB^2(t) = 10 log10(xi_gen^2(t)).

    Non-positive xi_gen^2 values are masked before taking the logarithm so the
    plot does not show infinities from invalid numerical points.
    """
    import matplotlib.pyplot as plt

    if isinstance(result, TrajectoryEnsemble):
        xi_data = generalized_squeezing_for_ensemble(
            result,
            tol=tol,
            n_processes=n_processes,
            verbose=verbose,
            averaged_observables=averaged_observables,
        )
        title = "Generalized three-level squeezing (ensemble)"
        phases = result.trajectories[0].phases
    else:
        xi_data = generalized_squeezing_for_trajectory(result, tol=tol)
        title = "Generalized three-level squeezing"
        phases = result.phases

    xi2 = np.asarray(xi_data["xi2_gen"], dtype=float)
    xi2_db = np.full_like(xi2, np.nan, dtype=float)
    valid = xi2 > 0.0
    xi2_db[valid] = 10.0 * np.log10(xi2[valid])
    xi_data["xi2_gen_db"] = xi2_db
    eigvals = np.asarray(xi_data["covariance_eigvals"], dtype=float)
    eigvals_plot = np.where(eigvals > 0.0, eigvals, np.nan)
    t_step1_end, t_step2_end = phase_change_times(phases)
    n_total = result.trajectories[0].N if isinstance(result, TrajectoryEnsemble) else result.N
    sector_list = result.trajectories[0].sectors if isinstance(result, TrajectoryEnsemble) else result.sectors
    dN = max(abs(int(Nj) - (n_total // 2)) for Nj in sector_list) if sector_list else 0
    Gamma_ref = result.trajectories[0].Gamma if isinstance(result, TrajectoryEnsemble) else result.Gamma
    theta_ss, _ = phase1_ss_angles_for_nj(n_total // 2, phases[0].omega, Gamma_ref)

    if ax is None:
        fig, axes = plt.subplots(2, 2, figsize=(14, 9.0), sharex=False)
        axes = axes.ravel()
    else:
        axes = np.asarray(ax)
        if axes.size != 4:
            raise ValueError(
                "plot_generalized_xi(..., ax=...) expects four axes for the "
                "squeezing, eigenvalue, N_c, and excited-fraction panels."
            )
        fig = axes.flat[0].figure
        axes = axes.ravel()

    axes[0].plot(
        xi_data["t"],
        xi_data["xi2_gen_db"],
        linewidth=1.8,
        label=(
            r"$10\log_{10}(\xi^2)$, "
            r"$\xi^2 = \frac{N\lambda_{\min}(C)}{\langle N_c/2\rangle^2}$"
        ),
    )
    if np.isfinite(theta_ss):
        xi2_ss = np.cos(theta_ss)
        if xi2_ss > 0.0:
            xi2_ss_db = 10.0 * np.log10(xi2_ss)
            axes[0].hlines(
                y=xi2_ss_db,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"$\delta = 0,\ \xi^2 = \cos\tilde{\Theta}_J$",
            )
    axes[0].set_xlabel("time")
    axes[0].set_ylabel(r"$\xi^2$ [dB]")
    axes[0].set_title(title)
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    eigen_labels = [
        r"$\lambda_1$ (smallest)",
        r"$\lambda_2$",
        r"$\lambda_3$",
        r"$\lambda_4$ (largest)",
    ]
    for idx, curve_label in enumerate(eigen_labels):
        axes[1].plot(xi_data["t"], eigvals_plot[:, idx], linewidth=1.6, label=curve_label)
    axes[1].axhline(
        y=n_total / 4.0,
        linestyle="--",
        color="black",
        alpha=0.8,
        label=r"squeezed state ($\lambda_a = N/4$)",
    )
    axes[1].set_xlabel("time")
    axes[1].set_ylabel("Covariance eigenvalues")
    axes[1].set_title("Ordered covariance eigenvalues")
    axes[1].set_yscale("log")
    axes[1].grid(alpha=0.3)
    axes[1].legend()

    axes[2].plot(xi_data["t"], np.asarray(xi_data["N_c"], dtype=float), linewidth=1.8, label=r"$\langle N_c\rangle$")
    axes[2].axhline(
        y=n_total,
        linestyle="--",
        color="black",
        alpha=0.8,
        label=r"squeezed state ($N_c = N$)",
    )
    axes[2].set_xlabel("time")
    axes[2].set_ylabel(r"$\langle N_c\rangle$")
    axes[2].set_title(r"Mean active population $\langle N_c\rangle$")
    axes[2].grid(alpha=0.3)
    axes[2].ticklabel_format(axis="y", style="plain", useOffset=False)
    axes[2].yaxis.get_major_locator().set_params(integer=True)
    axes[2].legend()

    axes[3].plot(
        xi_data["t"],
        np.asarray(xi_data["excited_fraction_active"], dtype=float),
        linewidth=1.8,
        label=r"$\langle e\rangle=\langle N_e\rangle/\langle N_J\rangle$",
    )
    if np.isfinite(theta_ss):
        excited_fraction_ss = 0.5 * (1.0 - np.cos(theta_ss))
        axes[3].hlines(
            y=excited_fraction_ss,
            xmin=0.0,
            xmax=t_step1_end,
            linestyle=":",
            alpha=0.9,
            label=r"$\delta=0,\ P_e=\sin^2(\tilde{\Theta}_J/2)$",
        )
    axes[3].set_xlabel("time")
    axes[3].set_ylabel(r"$\langle e\rangle$")
    axes[3].set_title(r"Active-manifold excited fraction $\langle e\rangle$")
    axes[3].grid(alpha=0.3)
    axes[3].legend()

    for axis in axes:
        axis.axvline(t_step1_end, linestyle="--", color="black", alpha=0.6)
        axis.axvline(t_step2_end, linestyle="--", color="black", alpha=0.6)

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.suptitle(f"Squeezing Dynamics (N={n_total}, dN={dN})", y=0.965)
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        return xi_data, fig, axes

    fig.suptitle(f"Squeezing Dynamics (N={n_total}, dN={dN})", y=0.965)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    return xi_data, fig, axes
