from __future__ import annotations
from multiprocessing import pool
from time import time

from quantum_trajectories.sim import (
    build_precomputed_trajectory_data,
    simulate_single_trajectory,
)
from quantum_trajectories.operator_helpers import omega2_from_weighted_average
from parser.common import Array, Phase
from parser.quantum_trajectories import SectorKey, TrajectoryEnsemble, TrajectoryResult

from dataclasses import dataclass
from math import comb
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import expm_multiply
import multiprocessing as mp
import time
from tqdm.auto import tqdm


# -----------------------------------------------------------------------------
# Multiprocessing worker state
# -----------------------------------------------------------------------------
# Each worker process gets its own copy of this global dictionary.
# We use this to avoid passing large objects, especially `precomputed`, to every
# single trajectory task.
_WORKER_STATE = None


def _init_trajectory_worker(
    N,
    Gamma,
    phases,
    sector_coeffs,
    internal_sector_states,
    dt,
    num_snapshots,
    shifted_jump_operator,
    precomputed,
    omega_1,
    omega_2,
    N1,
    N2,
):
    """
    Initialize one multiprocessing worker.

    This runs once per worker process. The goal is to store large/read-only data
    such as precomputed operators and propagators in a process-local global
    variable, instead of sending them again for every trajectory.
    """
    global _WORKER_STATE

    _WORKER_STATE = {
        "N": N,
        "Gamma": Gamma,
        "phases": phases,
        "sector_coeffs": sector_coeffs,
        "internal_sector_states": internal_sector_states,
        "dt": dt,
        "num_snapshots": num_snapshots,
        "shifted_jump_operator": shifted_jump_operator,
        "precomputed": precomputed,
        "omega_1": omega_1,
        "omega_2": omega_2,
        "N1": N1,
        "N2": N2,
    }


def _simulate_single_trajectory_worker(seed_sequence: np.random.SeedSequence) -> TrajectoryResult:
    """
    Worker-side wrapper around simulate_single_trajectory(...).

    The only per-task input is the child SeedSequence. Everything else is read from
    _WORKER_STATE, which was initialized once when the worker process started.
    """
    if _WORKER_STATE is None:
        raise RuntimeError("Worker state was not initialized.")

    return simulate_single_trajectory(
        N=_WORKER_STATE["N"],
        Gamma=_WORKER_STATE["Gamma"],
        phases=_WORKER_STATE["phases"],
        sector_coeffs=_WORKER_STATE["sector_coeffs"],
        internal_sector_states=_WORKER_STATE["internal_sector_states"],
        dt=_WORKER_STATE["dt"],
        num_snapshots=_WORKER_STATE["num_snapshots"],
        seed_sequence=seed_sequence,
        shifted_jump_operator=_WORKER_STATE["shifted_jump_operator"],
        precomputed=_WORKER_STATE["precomputed"],
        omega_1=_WORKER_STATE["omega_1"],
        omega_2=_WORKER_STATE["omega_2"],
        N1=_WORKER_STATE["N1"],
        N2=_WORKER_STATE["N2"],
    )


def run_trajectory_ensemble(
    N: int,
    Gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[SectorKey, complex],
    *,
    internal_sector_states: Optional[Mapping[SectorKey, Array]] = None,
    dt: float = 1e-3,
    num_snapshots: int = 101,
    seed: Optional[int] = None,
    shifted_jump_operator: bool = False,
    omega_1: Optional[float] = None,
    N1: Optional[int] = None,
    N2: Optional[int] = None,
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

    num_snapshots
        Number of saved snapshots per trajectory. All trajectories use the same
        internally constructed t_eval grid, so ensemble observables and
        squeezing moments are aligned at identical times without interpolation.

    omega_1, N1, N2
        Inhomogeneous-coupling metadata for tuple sector keys.  The group-2
        coupling is derived once from N1 * omega_1 + N2 * omega_2 = N1 + N2.

    verbose
        If True, print setup timing such as precompute and multiprocessing pool
        startup. Defaults to False so batch analyses stay quiet.
    """

    # Basic input validation
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")
    if num_snapshots < 2:
        raise ValueError("num_snapshots must be at least 2.")
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    # Spawn one child SeedSequence per trajectory from a single parent. This
    # makes trajectory 0 use the same child seed as a direct
    # simulate_single_trajectory(..., seed=seed) call.
    parent_seed_sequence = np.random.SeedSequence(seed)
    seed_sequences = parent_seed_sequence.spawn(ntraj)
    seed_keys = [tuple(child.spawn_key) for child in seed_sequences]
    omega_2 = (
        omega2_from_weighted_average(float(omega_1), int(N1), int(N2))
        if omega_1 is not None and N1 is not None and N2 is not None
        else None
    )
    omega_groups = None if omega_2 is None else (float(omega_1), omega_2)
    N_groups = None if N1 is None or N2 is None else (int(N1), int(N2))
    parameters = {
        "Gamma": Gamma,
        "phases": phases,
        "omega_groups": omega_groups,
        "N_groups": N_groups,
    }

    # Build all sector operators, effective generators, jump operators and
    # propagators once. These are reused by every trajectory.
    t0 = time.perf_counter()
    precomputed = build_precomputed_trajectory_data(
        N=N,
        Gamma=Gamma,
        phases=phases,
        sector_coeffs=sector_coeffs,
        dt=dt,
        shifted_jump_operator=shifted_jump_operator,
        omega_1=omega_1,
        omega_2=omega_2,
        N1=N1,
        N2=N2,
    )
    t1 = time.perf_counter()
    if verbose:
        print(f"Precompute: {t1 - t0:.2f} seconds.")

    # -------------------------------------------------------------------------
    # Serial fallback
    # -------------------------------------------------------------------------
    # Useful for debugging, profiling, and checking exact reproducibility.
    if n_processes is None or n_processes == 1:
        trajectories: List[TrajectoryResult] = []

        for child_seed_sequence in tqdm(seed_sequences, desc="simulate trajectories"):
            result = simulate_single_trajectory(
                N=N,
                Gamma=Gamma,
                phases=phases,
                sector_coeffs=sector_coeffs,
                internal_sector_states=internal_sector_states,
                dt=dt,
                num_snapshots=num_snapshots,
                seed_sequence=child_seed_sequence,
                shifted_jump_operator=shifted_jump_operator,
                precomputed=precomputed,
                omega_1=omega_1,
                omega_2=omega_2,
                N1=N1,
                N2=N2,
            )
            trajectories.append(result)

        total_steps = sum(traj.total_step_count for traj in trajectories)
        non_precomputed_steps = sum(traj.non_precomputed_step_count for traj in trajectories)
        avg_total_steps = total_steps / len(trajectories)
        avg_non_precomputed_steps = non_precomputed_steps / len(trajectories)
        print(
            "Simulation step summary (avg per trajectory): "
            f"total steps={avg_total_steps:.2f}, "
            f"steps without precompute={avg_non_precomputed_steps:.2f}"
        )
        return TrajectoryEnsemble(
            trajectories=trajectories,
            seeds=seed_keys,
            parameters=parameters,
        )

    # -------------------------------------------------------------------------
    # Parallel execution
    # -------------------------------------------------------------------------
    if n_processes == -1:
        n_processes = mp.cpu_count()

    if n_processes <= 0:
        raise ValueError("n_processes must be None, 1, -1, or a positive integer.")

    t0 = time.perf_counter()
    # Use Python's default multiprocessing context
    ctx = mp.get_context()
    with ctx.Pool(
        processes=n_processes,
        # initializer runs once per worker process
        initializer=_init_trajectory_worker,
        initargs=(
            N,
            Gamma,
            phases,
            sector_coeffs,
            internal_sector_states,
            dt,
            num_snapshots,
            shifted_jump_operator,
            precomputed,
            omega_1,
            omega_2,
            N1,
            N2,
        ),
    ) as pool:
        t1 = time.perf_counter()
        if verbose:
            print(f"Pool startup: {n_processes} processes in {t1 - t0:.2f} seconds.")

        trajectories = list(
            tqdm(
                pool.imap(
                    _simulate_single_trajectory_worker,
                    seed_sequences,
                    chunksize=chunksize,
                ),
                total=len(seed_sequences),
                desc="simulate trajectories",
            )
        )

    total_steps = sum(traj.total_step_count for traj in trajectories)
    non_precomputed_steps = sum(traj.non_precomputed_step_count for traj in trajectories)
    avg_total_steps = total_steps / len(trajectories)
    avg_non_precomputed_steps = non_precomputed_steps / len(trajectories)
    print(
        "Simulation step summary (avg per trajectory): "
        f"total steps={avg_total_steps:.2f}, "
        f"steps without precompute={avg_non_precomputed_steps:.2f}"
    )
    return TrajectoryEnsemble(
        trajectories=trajectories,
        seeds=seed_keys,
        parameters=parameters,
    )
