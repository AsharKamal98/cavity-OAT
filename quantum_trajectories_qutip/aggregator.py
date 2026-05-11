from __future__ import annotations

from typing import Optional, Sequence

import numpy as np

from common.parser import AveragedResult, ObservableSeries
from common.utils import active_manifold_angles
from quantum_trajectories_qutip.sim import jump_rate_from_state


def qutip_fixed_nj_observables(
    sim_data,
    *,
    tol: float = 1e-12,
) -> AveragedResult:
    """
    Convert fixed-N_J QuTiP solver output into AveragedResult.

    The returned observable series now also includes the physical jump rate.
    For mcsolve, the rate is averaged over runs when per-run states are present.
    The same function also accepts the corresponding mesolve output shape.
    """
    result = sim_data["result"]
    N = sim_data["N"]
    gamma = sim_data["gamma"]
    ntraj = sim_data["ntraj"]
    t = np.asarray(sim_data["tlist"], dtype=float)
    model = sim_data["model"]

    Jx = np.real(np.asarray(result.expect[0], dtype=float))
    Jy = np.real(np.asarray(result.expect[1], dtype=float))
    Jz = np.real(np.asarray(result.expect[2], dtype=float))
    N_e = np.real(np.asarray(result.expect[3], dtype=float))

    if Jx.ndim == 2:
        Jx_mean = Jx.mean(axis=0)
        Jy_mean = Jy.mean(axis=0)
        Jz_mean = Jz.mean(axis=0)
        N_e_mean = N_e.mean(axis=0)

        Jx_std = Jx.std(axis=0, ddof=0)
        Jy_std = Jy.std(axis=0, ddof=0)
        Jz_std = Jz.std(axis=0, ddof=0)
        N_e_std = N_e.std(axis=0, ddof=0)
    else:
        Jx_mean = Jx
        Jy_mean = Jy
        Jz_mean = Jz
        N_e_mean = N_e

        Jx_std = None
        Jy_std = None
        Jz_std = None
        N_e_std = None

    jump_rate_mean = None
    jump_rate_std = None

    states = sim_data.get("states")
    runs_states = sim_data.get("runs_states")

    if runs_states is not None:
        # mcsolve can expose one saved state list per trajectory; use that to
        # compute a mean and standard deviation for the jump rate.
        jump_rate_runs = []
        for run_states in runs_states:
            jump_rate_runs.append(
                [jump_rate_from_state(model, state, t_k) for state, t_k in zip(run_states, t)]
            )
        jump_rate_arr = np.asarray(jump_rate_runs, dtype=float)
        jump_rate_mean = jump_rate_arr.mean(axis=0)
        jump_rate_std = jump_rate_arr.std(axis=0, ddof=0)
    elif states is not None:
        # mesolve returns one state per saved time, so only a mean curve exists.
        jump_rate_mean = np.asarray(
            [jump_rate_from_state(model, state, t_k) for state, t_k in zip(states, t)],
            dtype=float,
        )
    else:
        raise ValueError(
            "QuTiP jump-rate extraction requires stored states or runs_states. "
            "Rerun the solver with state storage enabled."
        )

    N_j = np.full_like(t, float(N // 2), dtype=float)
    N_j_std = np.zeros_like(t, dtype=float) if Jx_std is not None else None

    # theta, phi, N_active, sx, sy, sz = active_manifold_angles(
    #     Jx_mean,
    #     Jy_mean,
    #     Jz_mean,
    #     N_e_mean,
    #     tol=tol,
    # )

    N_j = np.full_like(t, float(N // 2), dtype=float)

    sx = 2.0 * Jx_mean / N_j
    sy = 2.0 * Jy_mean / N_j

    # QuTiP convention: initial all-down has Jz = -Nj/2,
    # and we want theta = 0 there.
    sz = -2.0 * Jz_mean / N_j
    sz = np.clip(sz, -1.0, 1.0)

    theta = np.arccos(sz)

    phi = np.arctan2(sy, sx)
    r_perp = np.sqrt(sx**2 + sy**2)
    phi[r_perp < tol] = 0.0

    N_active = N_j

    obs = ObservableSeries(
        t=t,
        Jx=Jx_mean,
        Jy=Jy_mean,
        Jz=Jz_mean,
        N_e=N_e_mean,
        jump_rate=jump_rate_mean,
        N_j=N_j,
        N_active=N_active,
        theta=theta,
        phi=phi,
        sx=sx,
        sy=sy,
        sz=sz,
        Jx_std=Jx_std,
        Jy_std=Jy_std,
        Jz_std=Jz_std,
        N_e_std=N_e_std,
        jump_rate_std=jump_rate_std,
        N_j_std=N_j_std,
        N_active_std=None,
    )

    return AveragedResult(
        N=N,
        gamma=gamma,
        ntraj=ntraj,
        observables=obs,
    )


def qutip_fixed_nj_mcsolve_observables(
    sim_data,
    *,
    tol: float = 1e-12,
) -> AveragedResult:
    """
    Backward-compatible alias for qutip_fixed_nj_observables(...).

    The old name is kept so existing notebook cells keep working, but the new
    neutral name should be preferred because the function also handles mesolve
    output.
    """
    return qutip_fixed_nj_observables(sim_data, tol=tol)
