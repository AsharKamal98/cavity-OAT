from __future__ import annotations

import numpy as np
import qutip as qt

from parser.qutip import QutipMCSolverParameters, QutipMESolverParameters
from solvers.qutip_fixed_nj.models import (
    build_qutip_grouped_fixed_nj_model_from_phases,
)
from solvers.qutip_fixed_nj.utils_sim import (
    _observable_e_ops,
    _solver_args,
    build_tlist_from_phases,
)

def simulate_fixed_nj_me_trajectory(
    parameters: QutipMESolverParameters,
    *,
    num_points: int = 600,
    store_states: bool = True,
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

    Ni = [int(group_size) for group_size in parameters.Ni]
    omega_i = [float(coupling) for coupling in parameters.omega_i]
    NJi = [group_size // 2 for group_size in Ni]

    model = build_qutip_grouped_fixed_nj_model_from_phases(
        Gamma=parameters.Gamma,
        phases=parameters.phases,
        omega_i=omega_i,
        NJi=NJi,
        shifted_jump_operator=parameters.shifted_jump_operator,
    )
    print(
        f"Using QuTiP master equation solver fixed {len(Ni)}-group sector "
        f"NJi={NJi} with Ni={Ni} and omega_i={model.omega_i}."
    )
    tlist = build_tlist_from_phases(parameters.phases, num_points=num_points)
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
        "N": sum(parameters.Ni),
        "Gamma": parameters.Gamma,
        "ntraj": None,
        "tlist": tlist,
        "num_points": num_points,
        "states": result.states if store_states else None,
    }


def simulate_fixed_nj_mc_trajectory(
    parameters: QutipMCSolverParameters,
    *,
    num_points: int = 600,
    ntraj: int = 200,
    seed: int | None = None,
    keep_runs_results: bool = True,
    n_processes: int | None = None,
    norm_steps: int | None = None,
    norm_tol: float | None = None,
    norm_t_tol: float | None = None,
    norm_min_step: float | None = None,
):
    """
    Run fixed-N_J mcsolve benchmark and return raw solver output together with
    the metadata needed for later post-processing.

    The time grid is built internally as num_points equally spaced samples over
    the full protocol duration.

    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")
    if ntraj <= 0:
        raise ValueError("ntraj must be positive.")
    if n_processes is not None and n_processes != -1 and n_processes <= 0:
        raise ValueError("n_processes must be None, -1, or a positive integer.")

    Ni = [int(group_size) for group_size in parameters.Ni]
    omega_i = [float(coupling) for coupling in parameters.omega_i]
    NJi = [group_size // 2 for group_size in Ni]

    model = build_qutip_grouped_fixed_nj_model_from_phases(
        Gamma=parameters.Gamma,
        phases=parameters.phases,
        omega_i=omega_i,
        NJi=NJi,
        shifted_jump_operator=parameters.shifted_jump_operator,
    )
    print(
        f"Using QuTiP quantum trajectories fixed {len(Ni)}-group sector "
        f"NJi={NJi} with Ni={Ni} and omega_i={model.omega_i}."
    )
    tlist = build_tlist_from_phases(parameters.phases, num_points=num_points)
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
        "N": sum(parameters.Ni),
        "Gamma": parameters.Gamma,
        "ntraj": ntraj,
        "tlist": tlist,
        "num_points": num_points,
        "states": getattr(result, "states", None),
        "runs_states": getattr(result, "runs_states", None),
    }
