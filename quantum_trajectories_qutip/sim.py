from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import qutip as qt
#from Old.nj_sector_mc import phase_change_times, phase1_ss_angles_for_nj
from .utils import phase_change_times



@dataclass(frozen=True)
class QutipFixedNjModel:
    N: int
    N_J: int
    j: float
    gamma: float
    omega0: float
    delta0: float
    Jp: qt.Qobj
    Jm: qt.Qobj
    Jx: qt.Qobj
    Jy: qt.Qobj
    Jz: qt.Qobj
    N_e: qt.Qobj
    H: list
    c_ops: list
    psi0: qt.Qobj
    t_step1_end: float
    t_step2_end: float
    t_final: float


def _omega_coeff(t, args):
    if t < args["t_step1_end"]:
        return args["Omega0"]
    elif t < args["t_step2_end"]:
        return args["Omega0"]
    return 0.0


def _delta_coeff(t, args):
    if t < args["t_step1_end"]:
        return 0.0
    elif t < args["t_step2_end"]:
        return args["delta0"]
    return 0.0


def build_qutip_fixed_nj_model_from_phases(
    N: int,
    gamma: float,
    phases: Sequence,
) -> QutipFixedNjModel:
    """
    Build a fixed-N_J two-level collective model in QuTiP using the Dicke basis.

    Parameters
    ----------
    N
        Total atom number in the full three-level problem. The fixed two-level
        benchmark uses N_J = N//2.
    gamma
        Collective decay rate Gamma.
    phases
        Sequence with attributes .duration, .omega, .delta.
        Expected convention:
            phase 1: delta = 0, omega = Omega0
            phase 2: delta = delta0, omega = Omega0
            phase 3: delta = 0, omega = 0

    Returns
    -------
    QutipFixedNjModel
        QuTiP model object ready for mesolve or mcsolve.
    """
    if len(phases) < 3:
        raise ValueError("Need at least 3 phases.")
    if N <= 0:
        raise ValueError("N must be positive.")
    if N % 2 != 0:
        raise ValueError("Need even N so that N_J = N/2 is integer.")

    N_J = N // 2
    j = N_J / 2.0

    Omega0 = float(phases[0].omega)
    delta0 = float(phases[1].delta)
    t_step1_end = float(phases[0].duration)
    t_step2_end = float(phases[0].duration + phases[1].duration)
    t_final = float(sum(p.duration for p in phases))

    Jp = qt.jmat(j, "+")
    Jm = qt.jmat(j, "-")
    Jx = qt.jmat(j, "x")
    Jy = qt.jmat(j, "y")
    Jz = qt.jmat(j, "z")
    N_e = Jz + j * qt.qeye(int(2 * j + 1))

    H = [
        [Jx, _omega_coeff],
        [-N_e, _delta_coeff],
    ]
    c_ops = [np.sqrt(gamma) * Jm]

    # psi0 = qt.basis(int(2 * j + 1), 0)  # |m=-j> i.e. all active atoms in |down>
    dim = int(2 * j + 1)
    psi0 = qt.basis(dim, dim - 1)  # correct: all in |down>

    return QutipFixedNjModel(
        N=N,
        N_J=N_J,
        j=j,
        gamma=gamma,
        omega0=Omega0,
        delta0=delta0,
        Jp=Jp,
        Jm=Jm,
        Jx=Jx,
        Jy=Jy,
        Jz=Jz,
        N_e=N_e,
        H=H,
        c_ops=c_ops,
        psi0=psi0,
        t_step1_end=t_step1_end,
        t_step2_end=t_step2_end,
        t_final=t_final,
    )


def build_tlist_from_phases(phases: Sequence, num_points: int) -> np.ndarray:
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")
    t_final = float(sum(p.duration for p in phases))
    return np.linspace(0.0, t_final, num_points)


def _solver_args(model: QutipFixedNjModel) -> Dict[str, float]:
    return {
        "Omega0": model.omega0,
        "delta0": model.delta0,
        "t_step1_end": model.t_step1_end,
        "t_step2_end": model.t_step2_end,
    }


def simulate_fixed_nj_me_trajectory(
    N: int,
    gamma: float,
    phases: Sequence,
    *,
    tlist: Optional[np.ndarray] = None,
    num_points: int = 600,
    store_states: bool = True,
):
    """
    Run fixed-N_J mesolve benchmark using the same high-level inputs as your
    current main file.

    Returns a dictionary with keys:
        model, t, result, states, Jx, Jy, Jz, N_e
    """
    model = build_qutip_fixed_nj_model_from_phases(N=N, gamma=gamma, phases=phases)
    if tlist is None:
        tlist = build_tlist_from_phases(phases, num_points=num_points)

    options = {"store_states": store_states}
    result = qt.mesolve(
        model.H,
        model.psi0,
        tlist,
        c_ops=model.c_ops,
        e_ops=[model.Jx, model.Jy, model.Jz, model.N_e],
        args=_solver_args(model),
        options=options,
    )

    return {
        "model": model,
        "t": np.asarray(tlist, dtype=float),
        "result": result,
        "states": result.states if store_states else None,
        "Jx": np.real(np.asarray(result.expect[0], dtype=float)),
        "Jy": np.real(np.asarray(result.expect[1], dtype=float)),
        "Jz": np.real(np.asarray(result.expect[2], dtype=float)),
        "N_e": np.real(np.asarray(result.expect[3], dtype=float)),
    }


def simulate_fixed_nj_mc_trajectory(
    N: int,
    gamma: float,
    phases: Sequence,
    *,
    num_points: int = 600,
    ntraj: int = 200,
    seed: Optional[int] = None,
    keep_runs_results: bool = True,
):
    """
    Run fixed-N_J mcsolve benchmark and return raw solver output together with
    the metadata needed for later observable extraction.

    The time grid is built internally as num_points equally spaced samples over
    the full protocol duration.
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    model = build_qutip_fixed_nj_model_from_phases(N=N, gamma=gamma, phases=phases)
    tlist = build_tlist_from_phases(phases, num_points=num_points)
    tlist = np.asarray(tlist, dtype=float)

    options = {
        "progress_bar": "",
        "map": "serial",
        "store_states": True,
        "keep_runs_results": keep_runs_results,
    }

    result = qt.mcsolve(
        model.H,
        model.psi0,
        tlist,
        c_ops=model.c_ops,
        e_ops=[model.Jx, model.Jy, model.Jz, model.N_e],
        ntraj=ntraj,
        seeds=seed,
        args=_solver_args(model),
        options=options,
    )

    return {
        "result": result,
        "model": model,
        "N": N,
        "gamma": gamma,
        "ntraj": ntraj,
        "tlist": tlist,
        "num_points": num_points,
    }
