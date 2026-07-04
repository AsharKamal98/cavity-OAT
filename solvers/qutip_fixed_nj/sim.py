from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import qutip as qt

from solvers.qutip_fixed_nj.models import (
    build_qutip_fixed_nj_model_from_phases,
    build_qutip_two_group_fixed_nj_model_from_phases,
)
from solvers.qutip_fixed_nj.utils_sim import (
    _observable_e_ops,
    _solver_args,
    build_tlist_from_phases,
)

def simulate_fixed_nj_me_trajectory(
    N: int,
    Gamma: float,
    phases: Sequence,
    *,
    num_points: int = 600,
    store_states: bool = True,
    shifted_jump_operator: bool = False,
    sector_distribution: str = "square",
):
    """
    Run fixed-N_J mesolve benchmark and return raw solver output together with
    the metadata needed for later observable extraction.

    The output is intentionally shaped like simulate_fixed_nj_mc_trajectory(...),
    so the same QuTiP post-processing helpers, such as
    `compute_qutip_j_moments(...)`, can consume either solver path.
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    model = build_qutip_fixed_nj_model_from_phases(
        N=N,
        Gamma=Gamma,
        phases=phases,
        shifted_jump_operator=shifted_jump_operator,
        sector_distribution=sector_distribution,
    )
    tlist = build_tlist_from_phases(phases, num_points=num_points)
    tlist = np.asarray(tlist, dtype=float)

    options = {
        "store_states": store_states,
    }

    result = qt.mesolve(
        model.H,
        model.psi0,
        tlist,
        c_ops=model.c_ops,
        e_ops=_observable_e_ops(model),
        args=_solver_args(model),
        options=options,
    )

    return {
        "result": result,
        "model": model,
        "N": N,
        "Gamma": Gamma,
        "ntraj": None,
        "tlist": tlist,
        "num_points": num_points,
        "states": result.states if store_states else None,
    }


def simulate_fixed_nj_mc_trajectory(
    N: int,
    Gamma: float,
    phases: Sequence,
    *,
    num_points: int = 600,
    ntraj: int = 200,
    seed: Optional[int] = None,
    keep_runs_results: bool = True,
    shifted_jump_operator: bool = False,
    sector_distribution: str = "square",
    N1: Optional[int] = None,
    N2: Optional[int] = None,
    omega_1: float = 1.0,
    N_J1: Optional[int] = None,
    N_J2: Optional[int] = None,
    n_processes: Optional[int] = None,
    norm_steps: Optional[int] = None,
    norm_tol: Optional[float] = None,
    norm_t_tol: Optional[float] = None,
    norm_min_step: Optional[float] = None,
):
    """
    Run fixed-N_J mcsolve benchmark and return raw solver output together with
    the metadata needed for later post-processing.

    The time grid is built internally as num_points equally spaced samples over
    the full protocol duration.

    Parameters
    ----------
    sector_distribution
        Initial N_J-sector distribution choice shared with the custom MCWF
        helpers. In the fixed-N_J QuTiP benchmark this only validates and
        documents the initialization choice, because the simulation keeps the
        single sector N_J = N/2.
    n_processes
        QuTiP trajectory parallelism setting. ``None`` preserves the default
        serial behavior. ``-1`` uses all available CPU cores. Positive values
        request that many worker processes.
    norm_steps, norm_tol, norm_t_tol, norm_min_step
        Optional QuTiP collapse-time search controls forwarded to ``mcsolve``.
        Leave them as ``None`` to use QuTiP's defaults.
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    if N1 is None:
        N1 = N // 2
    if N2 is None:
        N2 = N - N1

    model = build_qutip_two_group_fixed_nj_model_from_phases(
        N=N,
        Gamma=Gamma,
        phases=phases,
        N1=N1,
        N2=N2,
        omega_1=omega_1,
        N_J1=N_J1,
        N_J2=N_J2,
        shifted_jump_operator=shifted_jump_operator,
    )
    print(
        "Using QuTiP fixed two-group sector "
        f"(N_J1, N_J2)=({model.N_J1}, {model.N_J2}) "
        f"with N1={model.N1}, N2={model.N2}, omega_1={model.omega_1}."
    )
    tlist = build_tlist_from_phases(phases, num_points=num_points)
    tlist = np.asarray(tlist, dtype=float)

    options = {
        "progress_bar": "",
        "store_states": True,
        "keep_runs_results": keep_runs_results,
    }
    if n_processes is None or n_processes == 1:
        options["map"] = "serial"
    else:
        options["map"] = "parallel"
        if n_processes > 1:
            options["num_cpus"] = n_processes
    if norm_steps is not None:
        options["norm_steps"] = norm_steps
    if norm_tol is not None:
        options["norm_tol"] = norm_tol
    if norm_t_tol is not None:
        options["norm_t_tol"] = norm_t_tol
    if norm_min_step is not None:
        options["norm_min_step"] = norm_min_step

    result = qt.mcsolve(
        model.H,
        model.psi0,
        tlist,
        c_ops=model.c_ops,
        e_ops=_observable_e_ops(model),
        ntraj=ntraj,
        seeds=seed,
        args=_solver_args(model),
        options=options,
    )

    return {
        "result": result,
        "model": model,
        "N": N,
        "Gamma": Gamma,
        "ntraj": ntraj,
        "tlist": tlist,
        "num_points": num_points,
        "states": getattr(result, "states", None),
        "runs_states": getattr(result, "runs_states", None),
    }
