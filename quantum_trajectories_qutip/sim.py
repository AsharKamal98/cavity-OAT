from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
import qutip as qt
#from Old.nj_sector_mc import phase_change_times, phase1_ss_angles_for_nj
from quantum_trajectories.state_helpers import centered_sector_initial_coeffs
from common.utils import omega2_from_weighted_average
from .utils import phase_change_times



@dataclass(frozen=True)
class QutipFixedNjModel:
    N: int
    N_J: int
    j: float
    Gamma: float
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


@dataclass(frozen=True)
class QutipTwoGroupFixedNjModel:
    N: int
    N1: int
    N2: int
    N_J1: int
    N_J2: int
    N_J: int
    j1: float
    j2: float
    Gamma: float
    shifted_jump_operator: bool
    unraveling_picture: str
    omega0: float
    delta0: float
    omega_1: float
    omega_2: float
    Jp: qt.Qobj
    Jm: qt.Qobj
    Jx: qt.Qobj
    Jy: qt.Qobj
    Jz: qt.Qobj
    N_e: qt.Qobj
    Jp_groups: tuple[qt.Qobj, qt.Qobj]
    Jm_groups: tuple[qt.Qobj, qt.Qobj]
    Jx_groups: tuple[qt.Qobj, qt.Qobj]
    Jy_groups: tuple[qt.Qobj, qt.Qobj]
    Jz_groups: tuple[qt.Qobj, qt.Qobj]
    N_e_groups: tuple[qt.Qobj, qt.Qobj]
    J_drive: qt.Qobj
    A_weighted: qt.Qobj
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


@dataclass(frozen=True)
class OmegaCoeffFromPhases:
    """
    Pickle-safe omega(t) coefficient for time-dependent QuTiP operators.

    QuTiP may probe coefficients before solver args are attached, so this
    callable falls back to the phase boundaries captured at model-build time.
    Keeping it as a top-level dataclass makes it safe to serialize for
    multiprocessing-based mcsolve runs.
    """

    omega0: float
    t_step1_end: float
    t_step2_end: float

    def __call__(self, t, args=None):
        if args is not None and "t_step1_end" in args:
            return _omega_coeff(t, args)
        if t < self.t_step1_end:
            return self.omega0
        if t < self.t_step2_end:
            return self.omega0
        return 0.0


def build_qutip_fixed_nj_model_from_phases(
    N: int,
    Gamma: float,
    phases: Sequence,
    *,
    shifted_jump_operator: bool = False,
    sector_distribution: str = "square",
) -> QutipFixedNjModel:
    """
    Build a fixed-N_J two-level collective model in QuTiP using the Dicke basis.

    Parameters
    ----------
    N
        Total atom number in the full three-level problem. The fixed two-level
        benchmark uses N_J = N//2.
    Gamma
        Collective decay rate Gamma.
    phases
        Sequence with attributes .duration, .omega, .delta.
        Expected convention:
            phase 1: delta = 0, omega = Omega0
            phase 2: delta = delta0, omega = Omega0
            phase 3: delta = 0, omega = 0
    sector_distribution
        Initial N_J-sector distribution choice shared with the custom MCWF
        initialization helpers. The fixed-N_J QuTiP benchmark only keeps the
        single central sector N_J = N/2, so both "square" and "binomial"
        reduce to the same one-sector initial state. The option is accepted for
        API consistency and validation.

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
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    initial_sector_coeffs = centered_sector_initial_coeffs(
        N,
        half_width=0,
        sector_distribution=sector_distribution,
    )
    N_J = next(iter(initial_sector_coeffs))
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
    omega_coeff_local = OmegaCoeffFromPhases(Omega0, t_step1_end, t_step2_end)

    if shifted_jump_operator:
        H = [
            [-N_e, _delta_coeff],
        ]
        shifted_c_op = qt.QobjEvo([
            np.sqrt(Gamma) * Jm,
            [1j / np.sqrt(Gamma) * identity, omega_coeff_local],
        ])
        c_ops = [shifted_c_op]
        unraveling_picture = "shifted"
    else:
        H = [
            [Jx, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(Gamma) * Jm]
        unraveling_picture = "regular"

    # psi0 = qt.basis(int(2 * j + 1), 0)  # |m=-j> i.e. all active atoms in |down>
    dim = int(2 * j + 1)
    psi0 = qt.basis(dim, dim - 1)  # correct: all in |down>

    return QutipFixedNjModel(
        N=N,
        N_J=N_J,
        j=j,
        Gamma=Gamma,
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


def _default_two_group_fixed_sector(N: int, N1: int, N2: int) -> tuple[int, int]:
    """Choose the central fixed active-sector split closest to group sizes."""
    if N % 2 != 0:
        raise ValueError("Need even N so that central N_J = N/2 is integer.")
    N_J = N // 2
    candidates = [
        (Nj1, N_J - Nj1)
        for Nj1 in range(max(0, N_J - N2), min(N1, N_J) + 1)
    ]
    if not candidates:
        raise ValueError("No valid two-group fixed sector exists for the given group sizes.")
    target_Nj1 = N_J * N1 / N
    return min(candidates, key=lambda pair: (abs(pair[0] - target_Nj1), pair[0]))


def build_qutip_two_group_fixed_nj_model_from_phases(
    N: int,
    Gamma: float,
    phases: Sequence,
    *,
    N1: int,
    N2: int,
    omega_1: float = 1.0,
    N_J1: Optional[int] = None,
    N_J2: Optional[int] = None,
    shifted_jump_operator: bool = False,
) -> QutipTwoGroupFixedNjModel:
    """
    Build a fixed two-group inhomogeneous collective model in QuTiP.

    This is the inhomogeneous counterpart of
    `build_qutip_fixed_nj_model_from_phases(...)`, but it keeps one fixed
    group-resolved sector `(N_J1, N_J2)` rather than a superposition of sectors.
    """
    if len(phases) < 3:
        raise ValueError("Need at least 3 phases.")
    if N <= 0:
        raise ValueError("N must be positive.")
    if N1 < 0 or N2 < 0:
        raise ValueError("N1 and N2 must be non-negative.")
    if N1 + N2 != N:
        raise ValueError(f"Expected N1 + N2 = N, got N1={N1}, N2={N2}, N={N}.")
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    if N_J1 is None or N_J2 is None:
        if N_J1 is not None or N_J2 is not None:
            raise ValueError("Provide both N_J1 and N_J2, or neither.")
        N_J1, N_J2 = _default_two_group_fixed_sector(N, N1, N2)
    if N_J1 < 0 or N_J2 < 0 or N_J1 > N1 or N_J2 > N2:
        raise ValueError(
            f"Invalid fixed sector ({N_J1}, {N_J2}) for group sizes N1={N1}, N2={N2}."
        )

    N_J = int(N_J1 + N_J2)
    j1 = N_J1 / 2.0
    j2 = N_J2 / 2.0
    omega_2 = omega2_from_weighted_average(float(omega_1), int(N1), int(N2))

    Omega0 = float(phases[0].omega)
    delta0 = float(phases[1].delta)
    t_step1_end = float(phases[0].duration)
    t_step2_end = float(phases[0].duration + phases[1].duration)
    t_final = float(sum(p.duration for p in phases))

    dim1 = int(2 * j1 + 1)
    dim2 = int(2 * j2 + 1)
    I1 = qt.qeye(dim1)
    I2 = qt.qeye(dim2)
    identity = qt.tensor(I1, I2)

    Jp1 = qt.tensor(qt.jmat(j1, "+"), I2)
    Jm1 = qt.tensor(qt.jmat(j1, "-"), I2)
    Jx1 = qt.tensor(qt.jmat(j1, "x"), I2)
    Jy1 = qt.tensor(qt.jmat(j1, "y"), I2)
    Jz1 = qt.tensor(qt.jmat(j1, "z"), I2)

    Jp2 = qt.tensor(I1, qt.jmat(j2, "+"))
    Jm2 = qt.tensor(I1, qt.jmat(j2, "-"))
    Jx2 = qt.tensor(I1, qt.jmat(j2, "x"))
    Jy2 = qt.tensor(I1, qt.jmat(j2, "y"))
    Jz2 = qt.tensor(I1, qt.jmat(j2, "z"))

    Jp = Jp1 + Jp2
    Jm = Jm1 + Jm2
    Jx = Jx1 + Jx2
    Jy = Jy1 + Jy2
    Jz = Jz1 + Jz2

    N_e1 = Jz1 + j1 * identity
    N_e2 = Jz2 + j2 * identity
    N_e = N_e1 + N_e2
    J_drive = omega_1 * Jx1 + omega_2 * Jx2
    A_weighted = omega_1 * Jm1 + omega_2 * Jm2

    omega_coeff_local = OmegaCoeffFromPhases(Omega0, t_step1_end, t_step2_end)
    if shifted_jump_operator:
        H = [
            [-N_e, _delta_coeff],
        ]
        shifted_c_op = qt.QobjEvo([
            np.sqrt(Gamma) * A_weighted,
            [1j / np.sqrt(Gamma) * identity, omega_coeff_local],
        ])
        c_ops = [shifted_c_op]
        unraveling_picture = "shifted"
    else:
        H = [
            [J_drive, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(Gamma) * A_weighted]
        unraveling_picture = "regular"

    psi0 = qt.tensor(
        qt.basis(dim1, dim1 - 1),
        qt.basis(dim2, dim2 - 1),
    )

    return QutipTwoGroupFixedNjModel(
        N=N,
        N1=N1,
        N2=N2,
        N_J1=int(N_J1),
        N_J2=int(N_J2),
        N_J=N_J,
        j1=j1,
        j2=j2,
        Gamma=Gamma,
        shifted_jump_operator=shifted_jump_operator,
        unraveling_picture=unraveling_picture,
        omega0=Omega0,
        delta0=delta0,
        omega_1=float(omega_1),
        omega_2=omega_2,
        Jp=Jp,
        Jm=Jm,
        Jx=Jx,
        Jy=Jy,
        Jz=Jz,
        N_e=N_e,
        Jp_groups=(Jp1, Jp2),
        Jm_groups=(Jm1, Jm2),
        Jx_groups=(Jx1, Jx2),
        Jy_groups=(Jy1, Jy2),
        Jz_groups=(Jz1, Jz2),
        N_e_groups=(N_e1, N_e2),
        J_drive=J_drive,
        A_weighted=A_weighted,
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


def _observable_e_ops(model) -> List[qt.Qobj]:
    """
    Return the standard observable list for QuTiP benchmarks.

    The first four entries are always the full collective observables
    (Jx, Jy, Jz, N_e). If the model also carries group-resolved operators,
    append those in the same order for each group.
    """
    e_ops = [model.Jx, model.Jy, model.Jz, model.N_e]
    if hasattr(model, "Jx_groups") and hasattr(model, "N_e_groups"):
        e_ops.extend(
            [
                model.Jx_groups[0],
                model.Jx_groups[1],
                model.Jy_groups[0],
                model.Jy_groups[1],
                model.Jz_groups[0],
                model.Jz_groups[1],
                model.N_e_groups[0],
                model.N_e_groups[1],
            ]
        )
    return e_ops


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
    so qutip_fixed_nj_observables(...) can also consume it.
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
    the metadata needed for later observable extraction.

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
