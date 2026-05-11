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
    shifted_jump_operator: bool
    unraveling_picture: str
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


def _make_omega_coeff_from_phases(Omega0: float, t_step1_end: float, t_step2_end: float):
    """
    Build an omega(t) coefficient that works both inside and outside a solver.

    QuTiP may probe coefficients before solver args are attached, so this local
    closure can fall back to the phase boundaries captured at model-build time.
    """
    def _omega_coeff_local(t, args=None):
        if args is not None and "t_step1_end" in args:
            return _omega_coeff(t, args)
        if t < t_step1_end:
            return Omega0
        elif t < t_step2_end:
            return Omega0
        return 0.0

    return _omega_coeff_local


def build_qutip_fixed_nj_model_from_phases(
    N: int,
    gamma: float,
    phases: Sequence,
    *,
    shifted_jump_operator: bool = False,
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
    if shifted_jump_operator and gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires gamma > 0 because the shifted jump "
            "operator contains omega / gamma."
        )

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
    identity = qt.qeye(int(2 * j + 1))
    N_e = Jz + j * identity
    omega_coeff_local = _make_omega_coeff_from_phases(Omega0, t_step1_end, t_step2_end)

    if shifted_jump_operator:
        H = [
            [-N_e, _delta_coeff],
        ]
        shifted_c_op = qt.QobjEvo([
            np.sqrt(gamma) * Jm,
            [1j / np.sqrt(gamma) * identity, omega_coeff_local],
        ])
        c_ops = [shifted_c_op]
        unraveling_picture = "shifted"
    else:
        H = [
            [Jx, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(gamma) * Jm]
        unraveling_picture = "regular"

    # psi0 = qt.basis(int(2 * j + 1), 0)  # |m=-j> i.e. all active atoms in |down>
    dim = int(2 * j + 1)
    psi0 = qt.basis(dim, dim - 1)  # correct: all in |down>

    return QutipFixedNjModel(
        N=N,
        N_J=N_J,
        j=j,
        gamma=gamma,
        shifted_jump_operator=shifted_jump_operator,
        unraveling_picture=unraveling_picture,
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


def collapse_ops_at_time(model: QutipFixedNjModel, t: float) -> List[qt.Qobj]:
    """
    Evaluate the model's actual collapse operators at time t.

    Regular-picture collapse operators are plain Qobj instances.
    Shifted-picture collapse operators are stored as QobjEvo and can be
    evaluated directly at time t in the installed QuTiP API.
    """
    args = _solver_args(model)
    evaluated = []
    for c_op in model.c_ops:
        # if time-dependent, c_op should be a QobjEvo; evaluate at time t
        if isinstance(c_op, qt.QobjEvo):
            evaluated.append(c_op(t, args=args))
        # else t-independent, c_op should be a Qobj; use as-is
        else:
            evaluated.append(c_op)
    return evaluated


def jump_rate_operator_at_time(model: QutipFixedNjModel, t: float) -> qt.Qobj:
    """
    Return sum_c c^\dagger(t) c(t) from the model's actual collapse operators.
    """
    collapse_ops = collapse_ops_at_time(model, t)
    rate_operator = 0 * model.Jm
    # loop through decay channels
    for c_op in collapse_ops:
        rate_operator = rate_operator + c_op.dag() * c_op
    return rate_operator


def jump_rate_from_state(model: QutipFixedNjModel, state: qt.Qobj, t: float) -> float:
    """Evaluate the physical jump rate for a ket or density matrix at time t."""
    rate_operator = jump_rate_operator_at_time(model, t)
    if state.isket:
        value = qt.expect(rate_operator, state)
    else:
        value = (state * rate_operator).tr()
    jump_rate = float(np.real(value))
    # Numerical solvers can produce tiny negative roundoff errors near zero.
    if jump_rate < 0.0 and abs(jump_rate) < 1e-10:
        return 0.0
    return jump_rate


def simulate_fixed_nj_me_trajectory(
    N: int,
    gamma: float,
    phases: Sequence,
    *,
    num_points: int = 600,
    store_states: bool = True,
    shifted_jump_operator: bool = False,
):
    """
    Run fixed-N_J mesolve benchmark and return raw solver output together with
    the metadata needed for later observable extraction.

    The output is intentionally shaped like simulate_fixed_nj_mc_trajectory(...),
    so qutip_fixed_nj_observables(...) can also consume it.
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    model = build_qutip_fixed_nj_model_from_phases(
        N=N,
        gamma=gamma,
        phases=phases,
        shifted_jump_operator=shifted_jump_operator,
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
        e_ops=[model.Jx, model.Jy, model.Jz, model.N_e],
        args=_solver_args(model),
        options=options,
    )

    return {
        "result": result,
        "model": model,
        "N": N,
        "gamma": gamma,
        "ntraj": None,
        "tlist": tlist,
        "num_points": num_points,
        "states": result.states if store_states else None,
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
    shifted_jump_operator: bool = False,
):
    """
    Run fixed-N_J mcsolve benchmark and return raw solver output together with
    the metadata needed for later observable extraction.

    The time grid is built internally as num_points equally spaced samples over
    the full protocol duration.
    """
    if num_points < 2:
        raise ValueError("num_points must be at least 2.")

    model = build_qutip_fixed_nj_model_from_phases(
        N=N,
        gamma=gamma,
        phases=phases,
        shifted_jump_operator=shifted_jump_operator,
    )
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
        "states": getattr(result, "states", None),
        "runs_states": getattr(result, "runs_states", None),
    }
