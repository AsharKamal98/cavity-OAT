from __future__ import annotations
from multiprocessing import pool
from time import time

from quantum_trajectories.sim import (
    build_precomputed_trajectory_data,
    simulate_single_trajectory,
)
from quantum_trajectories.parser import (
    Array,
    Phase,
    TrajectoryEnsemble,
    TrajectoryResult
)

from dataclasses import dataclass
from math import comb
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import numpy as np
from scipy.sparse import csc_matrix, diags
from scipy.sparse.linalg import expm_multiply
import multiprocessing as mp
import time


# -----------------------------------------------------------------------------
# Multiprocessing worker state
# -----------------------------------------------------------------------------
# Each worker process gets its own copy of this global dictionary.
# We use this to avoid passing large objects, especially `precomputed`, to every
# single trajectory task.
_WORKER_STATE = None


def _init_trajectory_worker(
    N,
    gamma,
    phases,
    sector_coeffs,
    internal_sector_states,
    dt,
    save_every,
    shifted_jump_operator,
    precomputed,
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
        "gamma": gamma,
        "phases": phases,
        "sector_coeffs": sector_coeffs,
        "internal_sector_states": internal_sector_states,
        "dt": dt,
        "save_every": save_every,
        "shifted_jump_operator": shifted_jump_operator,
        "precomputed": precomputed,
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
        gamma=_WORKER_STATE["gamma"],
        phases=_WORKER_STATE["phases"],
        sector_coeffs=_WORKER_STATE["sector_coeffs"],
        internal_sector_states=_WORKER_STATE["internal_sector_states"],
        dt=_WORKER_STATE["dt"],
        save_every=_WORKER_STATE["save_every"],
        seed_sequence=seed_sequence,
        shifted_jump_operator=_WORKER_STATE["shifted_jump_operator"],
        precomputed=_WORKER_STATE["precomputed"],
    )

def run_trajectory_ensemble(
    N: int,
    gamma: float,
    phases: Sequence[Phase],
    sector_coeffs: Mapping[int, complex],
    *,
    internal_sector_states: Optional[Mapping[int, Array]] = None,
    dt: float = 1e-3,
    save_every: int = 1,
    seed: Optional[int] = None,
    shifted_jump_operator: bool = False,
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

    verbose
        If True, print setup timing such as precompute and multiprocessing pool
        startup. Defaults to False so batch analyses stay quiet.
    """

    # Basic input validation
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")
    if shifted_jump_operator and gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires gamma > 0 because the shifted jump "
            "operator contains omega / gamma."
        )

    # Spawn one child SeedSequence per trajectory from a single parent. This
    # makes trajectory 0 use the same child seed as a direct
    # simulate_single_trajectory(..., seed=seed) call.
    parent_seed_sequence = np.random.SeedSequence(seed)
    seed_sequences = parent_seed_sequence.spawn(ntraj)
    seed_keys = [tuple(child.spawn_key) for child in seed_sequences]

    # Build all sector operators, effective generators, jump operators and
    # propagators once. These are reused by every trajectory.
    t0 = time.perf_counter()
    precomputed = build_precomputed_trajectory_data(
        N=N,
        gamma=gamma,
        phases=phases,
        sector_coeffs=sector_coeffs,
        dt=dt,
        shifted_jump_operator=shifted_jump_operator,
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

        for child_seed_sequence in seed_sequences:
            result = simulate_single_trajectory(
                N=N,
                gamma=gamma,
                phases=phases,
                sector_coeffs=sector_coeffs,
                internal_sector_states=internal_sector_states,
                dt=dt,
                save_every=save_every,
                seed_sequence=child_seed_sequence,
                shifted_jump_operator=shifted_jump_operator,
                precomputed=precomputed,
            )
            trajectories.append(result)

        return TrajectoryEnsemble(trajectories=trajectories, seeds=seed_keys)

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
            gamma,
            phases,
            sector_coeffs,
            internal_sector_states,
            dt,
            save_every,
            shifted_jump_operator,
            precomputed,
        ),
    ) as pool:
        t1 = time.perf_counter()
        if verbose:
            print(f"Pool startup: {n_processes} processes in {t1 - t0:.2f} seconds.")

        trajectories = pool.map(
            _simulate_single_trajectory_worker,
            seed_sequences,
            chunksize=chunksize,
        )

    return TrajectoryEnsemble(trajectories=trajectories, seeds=seed_keys)
