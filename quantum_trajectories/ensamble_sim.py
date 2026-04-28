from __future__ import annotations

from quantum_trajectories.sim import simulate_single_trajectory
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
    ntraj: int,
) -> TrajectoryEnsemble:
    """
    Run many independent MC trajectories and return them as a list.

    This is intentionally simple so it can be parallelized later by replacing
    the for-loop with multiprocessing / joblib / concurrent.futures.
    """
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")

    master_rng = np.random.default_rng(seed)
    seeds = master_rng.integers(0, 2**32 - 1, size=ntraj, dtype=np.uint64).tolist()

    trajectories: List[TrajectoryResult] = []
    for s in seeds:
        rng = np.random.default_rng(int(s))
        result = simulate_single_trajectory(
            N=N,
            gamma=gamma,
            phases=phases,
            sector_coeffs=sector_coeffs,
            internal_sector_states=internal_sector_states,
            dt=dt,
            save_every=save_every,
            seed=int(s),
        )
        trajectories.append(result)

    return TrajectoryEnsemble(trajectories=trajectories, seeds=[int(s) for s in seeds])
