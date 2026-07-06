from __future__ import annotations
import sys
from multiprocessing import pool
from time import time

from solvers.mcwf.sim import (
    _simulate_single_trajectory,
    build_precomputed_trajectory_data,
)
from common.utils.parameters import omega_G_from_weighted_average
from common.utils.parameters import check_initial_sector_omega_ratio
from parser.common import Array, Phase
from parser.mcwf import (
    MCWFSolverParameters,
    SectorKey,
    TrajectoryEnsemble,
    TrajectoryEnsembleMetadata,
    TrajectoryResult,
)
from solvers.mcwf.state_helpers import centered_sector_initial_coeffs

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
    Ni,
    Gamma,
    phases,
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
        "phases": phases,
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
        phases=_WORKER_STATE["phases"],
        sector_coeffs=_WORKER_STATE["sector_coeffs"],
        dt=_WORKER_STATE["dt"],
        t_eval=_WORKER_STATE["t_eval"],
        seed_sequence=seed_sequence,
        shifted_jump_operator=_WORKER_STATE["shifted_jump_operator"],
        precomputed=_WORKER_STATE["precomputed"],
        omega_i=_WORKER_STATE["omega_i"],
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
        Validated MCWF solver inputs. This includes `Ni`, `dN`, the first
        `G-1` group couplings in `omega_i`, `Gamma`, `phases`, the sector
        distribution, the timestep, and whether to use the shifted jump
        operator. The final coupling is completed internally from the
        weighted-average condition.

    verbose
        If True, print setup timing such as precompute and multiprocessing pool
        startup. Defaults to False so batch analyses stay quiet.
    """

    # Basic input validation
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")
    t_eval = np.asarray(t_eval, dtype=float)
    if t_eval.ndim != 1 or t_eval.size < 2:
        raise ValueError("t_eval must be a one-dimensional array with at least two points.")
    if np.any(np.diff(t_eval) <= 0.0):
        raise ValueError("t_eval must be strictly increasing.")
    total_time = float(sum(phase.duration for phase in parameters.phases))
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
    omega_i = omega_i + [omega_G_from_weighted_average(omega_i, Ni)]

    # initial sector coefficients 
    sector_coeffs = centered_sector_initial_coeffs(
        Ni,
        dN=parameters.dN,
        sector_distribution=parameters.sector_distribution,
    )
    ratio_check = check_initial_sector_omega_ratio(
        sector_coeffs,
        Omega=max(abs(phase.omega) for phase in parameters.phases),
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
        phases=parameters.phases,
        sector_coeffs=sector_coeffs,
        dt=parameters.dt,
        shifted_jump_operator=parameters.shifted_jump_operator,
        omega_i=omega_i,
    )
    t1 = time.perf_counter()
    if verbose:
        print(f"Precompute: {t1 - t0:.2f} seconds.")

    # -------------------------------------------------------------------------
    # Serial fallback
    # -------------------------------------------------------------------------
    # Useful for debugging, profiling, and checking exact reproducibility.
    if n_processes in (None, 1) or ntraj == 1:
        trajectories: List[TrajectoryResult] = []

        for child_seed_sequence in tqdm(seed_sequences, desc="simulate trajectories"):
            result = _simulate_single_trajectory(
                Ni=Ni,
                Gamma=parameters.Gamma,
                phases=parameters.phases,
                sector_coeffs=sector_coeffs,
                dt=parameters.dt,
                t_eval=t_eval,
                seed_sequence=child_seed_sequence,
                shifted_jump_operator=parameters.shifted_jump_operator,
                precomputed=precomputed,
                omega_i=omega_i,
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
        metadata = TrajectoryEnsembleMetadata(
            Ni=tuple(Ni),
            omega_i=tuple(omega_i),
            Gamma=parameters.Gamma,
            phases=list(parameters.phases),
            shifted_jump_operator=parameters.shifted_jump_operator,
            t_eval=t_eval.copy(),
            sectors=list(precomputed["sector_list"]),
            sector_multiplicities=dict(precomputed["multiplicities"]),
            sector_dimensions=dict(precomputed["dims"]),
        )
        return TrajectoryEnsemble(
            trajectories=trajectories,
            seeds=seed_keys,
            metadata=metadata,
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
            Ni,
            parameters.Gamma,
            parameters.phases,
            sector_coeffs,
            parameters.dt,
            t_eval,
            parameters.shifted_jump_operator,
            precomputed,
            omega_i,
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
    metadata = TrajectoryEnsembleMetadata(
        Ni=tuple(Ni),
        omega_i=tuple(omega_i),
        Gamma=parameters.Gamma,
        phases=list(parameters.phases),
        shifted_jump_operator=parameters.shifted_jump_operator,
        t_eval=t_eval.copy(),
        sectors=list(precomputed["sector_list"]),
        sector_multiplicities=dict(precomputed["multiplicities"]),
        sector_dimensions=dict(precomputed["dims"]),
    )
    return TrajectoryEnsemble(
        trajectories=trajectories,
        seeds=seed_keys,
        metadata=metadata,
    )
