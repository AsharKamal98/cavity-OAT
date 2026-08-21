from __future__ import annotations

import multiprocessing as mp
import sys
import time
from typing import Iterable, Optional

import numpy as np
from tqdm.auto import tqdm

from common.utils.parameters import check_initial_sector_omega_ratio
from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.mcwf import (
    MCWFSolverParameters,
    TrajectoryEnsemble,
    TrajectoryEnsembleMetadata,
    TrajectoryResult,
)
from solvers.mcwf.sim import (
    _simulate_single_trajectory,
    build_precomputed_trajectory_data,
)
from solvers.mcwf.state_helpers import centered_sector_initial_coeffs


# -----------------------------------------------------------------------------
# Multiprocessing worker state
# -----------------------------------------------------------------------------
# Each worker process gets its own copy of this global dictionary.
# We use this to avoid passing large objects, especially `precomputed`, to every
# single trajectory task.
_WORKER_STATE = None


def _init_trajectory_worker(
    Ni,
    Gamma,
    integration_phases,
    sector_coeffs,
    dt,
    t_eval,
    shifted_jump_operator,
    precomputed,
    omega_i,
):
    """
    Initialize one multiprocessing worker.

    This runs once per worker process. The goal is to store large/read-only data
    such as precomputed operators and propagators in a process-local global
    variable, instead of sending them again for every trajectory.
    """
    global _WORKER_STATE

    _WORKER_STATE = {
        "Ni": Ni,
        "Gamma": Gamma,
        "integration_phases": integration_phases,
        "sector_coeffs": sector_coeffs,
        "dt": dt,
        "t_eval": t_eval,
        "shifted_jump_operator": shifted_jump_operator,
        "precomputed": precomputed,
        "omega_i": omega_i,
    }


def _simulate_single_trajectory_worker(seed_sequence: np.random.SeedSequence) -> TrajectoryResult:
    """
    Worker-side wrapper around _simulate_single_trajectory(...).

    The only per-task input is the child SeedSequence. Everything else is read from
    _WORKER_STATE, which was initialized once when the worker process started.
    """
    if _WORKER_STATE is None:
        raise RuntimeError("Worker state was not initialized.")

    return _simulate_single_trajectory(
        Ni=_WORKER_STATE["Ni"],
        Gamma=_WORKER_STATE["Gamma"],
        integration_phases=_WORKER_STATE["integration_phases"],
        sector_coeffs=_WORKER_STATE["sector_coeffs"],
        dt=_WORKER_STATE["dt"],
        t_eval=_WORKER_STATE["t_eval"],
        seed_sequence=seed_sequence,
        shifted_jump_operator=_WORKER_STATE["shifted_jump_operator"],
        precomputed=_WORKER_STATE["precomputed"],
        omega_i=_WORKER_STATE["omega_i"],
    )


def _collect_trajectory_moments(
    results: Iterable[TrajectoryResult],
    *,
    ntraj: int,
    t_eval: Array,
    num_groups: int,
) -> tuple[TrajectoryResolvedMomentSeries, int, int]:
    """Stack single-trajectory moments and consume transient diagnostics."""
    generator_count = 8 * num_groups
    packed_count = generator_count * (generator_count + 1) // 2
    first_moments = np.empty((ntraj, t_eval.size, num_groups, 8), dtype=float)
    second_moments = np.empty((ntraj, t_eval.size, packed_count), dtype=float)
    total_step_count = 0
    non_precomputed_step_count = 0
    result_count = 0

    for trajectory_index, result in enumerate(results):
        if trajectory_index >= ntraj:
            raise RuntimeError("Received more trajectory results than requested.")

        moments = result.moments
        if not np.array_equal(moments.t, t_eval):
            raise RuntimeError("A trajectory returned a different saved-time grid.")
        if moments.first_moments.shape[2] != num_groups:
            raise RuntimeError("A trajectory returned a different physical-group count.")

        first_moments[trajectory_index] = moments.first_moments[0]
        second_moments[trajectory_index] = moments.second_moments[0]
        total_step_count += result.diagnostics.total_step_count
        non_precomputed_step_count += (
            result.diagnostics.non_precomputed_step_count
        )
        result_count += 1

    if result_count != ntraj:
        raise RuntimeError(
            f"Received {result_count} trajectory results for ntraj={ntraj}."
        )

    return (
        TrajectoryResolvedMomentSeries(
            t=np.asarray(t_eval, dtype=float).copy(),
            first_moments=first_moments,
            second_moments=second_moments,
        ),
        total_step_count,
        non_precomputed_step_count,
    )


def run_trajectory_ensemble(
    parameters: MCWFSolverParameters,
    *,
    t_eval: Array,
    seed: Optional[int] = None,
    ntraj: int,
    n_processes: Optional[int] = None,
    chunksize: int = 1,
    verbose: bool = False,
) -> TrajectoryEnsemble:
    """
    Run many independent MC trajectories.

    If n_processes is None or 1, this runs serially.
    If n_processes > 1, trajectories are distributed over multiple processes.

    Parameters
    ----------
    n_processes
        Number of worker processes.
        - None or 1: serial execution.
        - -1: use all available CPU cores.
        - >1: use that many processes.

    chunksize
        Number of seeds sent to a worker at once. Larger chunksize can reduce
        multiprocessing overhead, but chunksize=1 gives better load balancing
        if trajectories have very different numbers of jumps.

    t_eval
        Explicit saved-time grid shared by all trajectories in the ensemble, in
        the same style as the MFE solver.

    parameters
        Validated physical inputs and MCWF-specific controls.

    verbose
        If True, print setup timing such as precompute and multiprocessing pool
        startup. Defaults to False so batch analyses stay quiet.
    """

    # Basic input validation
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")
    if n_processes is not None and n_processes != -1 and n_processes <= 0:
        raise ValueError(
            "n_processes must be None, 1, -1, or a positive integer."
        )
    if chunksize <= 0:
        raise ValueError("chunksize must be positive.")
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional array with at least two points.")
    if not np.all(np.isfinite(t_eval)):
        raise ValueError("t_eval must contain only finite values.")
    if np.any(np.diff(t_eval) <= 0.0):
        raise ValueError("t_eval must be strictly increasing.")
    phase_protocol = parameters.phase_protocol
    integration_phases = phase_protocol.integration_phases
    total_time = phase_protocol.total_duration
    if abs(float(t_eval[0])) > 1e-12:
        raise ValueError("The first t_eval point must be 0.0.")
    if abs(float(t_eval[-1]) - total_time) > 1e-9:
        raise ValueError("The last t_eval point must match the total protocol time.")
    if parameters.shifted_jump_operator and parameters.Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    # Spawn one child SeedSequence per trajectory from a single parent. This
    # makes trajectory 0 use the same child seed as the private
    # _simulate_single_trajectory(...) helper called with the same child seed.
    Ni = [int(group_size) for group_size in parameters.Ni]
    omega_i = [float(coupling) for coupling in parameters.omega_i]

    # initial sector coefficients 
    sector_coeffs = centered_sector_initial_coeffs(
        Ni,
        dN=parameters.dN,
        sector_distribution=parameters.sector_distribution,
    )
    ratio_check = check_initial_sector_omega_ratio(
        sector_coeffs,
        Omega=max(abs(phase.omega) for phase in integration_phases),
        Gamma=parameters.Gamma,
    )
    if not ratio_check["is_valid"]:
        if len(Ni) == 1:
            sys.exit(
                "Omega/Omega_c check not valid for homogeneous run: "
                f"Omega={ratio_check['omega']}, Omega_c={ratio_check['omega_c']}, "
                f"smallest Nj={ratio_check['min_nj']}, ratio={ratio_check['ratio']}"
            )
        sys.exit(
            "Omega/Omega_c check not valid for inhomogeneous run: "
            f"Omega={ratio_check['omega']}, Omega_c={ratio_check['omega_c']}, "
            f"smallest Nj={ratio_check['min_nj']}, ratio={ratio_check['ratio']}"
        )

    parent_seed_sequence = np.random.SeedSequence(seed)
    seed_sequences = parent_seed_sequence.spawn(ntraj)
    seed_keys = [tuple(child.spawn_key) for child in seed_sequences]
    # Build all sector operators, effective generators, jump operators and
    # propagators once. These are reused by every trajectory.
    t0 = time.perf_counter()
    precomputed = build_precomputed_trajectory_data(
        Ni=Ni,
        Gamma=parameters.Gamma,
        integration_phases=integration_phases,
        sector_coeffs=sector_coeffs,
        dt=parameters.dt,
        shifted_jump_operator=parameters.shifted_jump_operator,
        omega_i=omega_i,
    )
    t1 = time.perf_counter()
    if verbose:
        print(f"Precompute: {t1 - t0:.2f} seconds.")

    metadata = TrajectoryEnsembleMetadata(
        Ni=tuple(Ni),
        omega_i=tuple(omega_i),
        Gamma=parameters.Gamma,
        phase_protocol=phase_protocol,
        shifted_jump_operator=parameters.shifted_jump_operator,
        t_eval=t_eval.copy(),
        sectors=list(precomputed["sector_list"]),
        sector_multiplicities=dict(precomputed["multiplicities"]),
        sector_dimensions=dict(precomputed["dims"]),
    )

    # -------------------------------------------------------------------------
    # Serial fallback
    # -------------------------------------------------------------------------
    # Useful for debugging, profiling, and checking exact reproducibility.
    if n_processes in (None, 1) or ntraj == 1:
        trajectory_results = (
            _simulate_single_trajectory(
                Ni=Ni,
                Gamma=parameters.Gamma,
                integration_phases=integration_phases,
                sector_coeffs=sector_coeffs,
                dt=parameters.dt,
                t_eval=t_eval,
                seed_sequence=child_seed_sequence,
                shifted_jump_operator=parameters.shifted_jump_operator,
                precomputed=precomputed,
                omega_i=omega_i,
            )
            for child_seed_sequence in tqdm(
                seed_sequences, desc="simulate trajectories"
            )
        )
        moments, total_steps, non_precomputed_steps = _collect_trajectory_moments(
            trajectory_results,
            ntraj=ntraj,
            t_eval=t_eval,
            num_groups=len(Ni),
        )
    else:
        # ---------------------------------------------------------------------
        # Parallel execution
        # ---------------------------------------------------------------------
        if n_processes == -1:
            n_processes = mp.cpu_count()

        t0 = time.perf_counter()
        # Use Python's default multiprocessing context.
        ctx = mp.get_context()
        with ctx.Pool(
            processes=n_processes,
            # Initializer runs once per worker process.
            initializer=_init_trajectory_worker,
            initargs=(
                Ni,
                parameters.Gamma,
                integration_phases,
                sector_coeffs,
                parameters.dt,
                t_eval,
                parameters.shifted_jump_operator,
                precomputed,
                omega_i,
            ),
        ) as worker_pool:
            t1 = time.perf_counter()
            if verbose:
                print(
                    f"Pool startup: {n_processes} processes in {t1 - t0:.2f} seconds."
                )

            trajectory_results = tqdm(
                worker_pool.imap(
                    _simulate_single_trajectory_worker,
                    seed_sequences,
                    chunksize=chunksize,
                ),
                total=len(seed_sequences),
                desc="simulate trajectories",
            )
            moments, total_steps, non_precomputed_steps = (
                _collect_trajectory_moments(
                    trajectory_results,
                    ntraj=ntraj,
                    t_eval=t_eval,
                    num_groups=len(Ni),
                )
            )

    avg_total_steps = total_steps / ntraj
    avg_non_precomputed_steps = non_precomputed_steps / ntraj
    print(
        "Simulation step summary (avg per trajectory): "
        f"total steps={avg_total_steps:.2f}, "
        f"steps without precompute={avg_non_precomputed_steps:.2f}"
    )
    return TrajectoryEnsemble(
        moments=moments,
        seeds=seed_keys,
        metadata=metadata,
    )
