from __future__ import annotations

from typing import Optional

import numpy as np
import qutip as qt

from parser.qutip import QutipFixedNjModel, QutipGroupedFixedNjModel
from solvers.qutip_fixed_nj.utils_sim import OmegaCoeffFromPhases, _delta_coeff, _omega_coeff


def build_qutip_fixed_nj_model_from_phases(
    N: int,
    Gamma: float,
    phases: list,
    *,
    shifted_jump_operator: bool = False,
) -> QutipFixedNjModel:
    """
    Build a fixed-NJ two-level collective model in QuTiP using the Dicke basis.
    """
    if len(phases) < 3:
        raise ValueError("Need at least 3 phases.")
    if N <= 0:
        raise ValueError("N must be positive.")
    if N % 2 != 0:
        raise ValueError("Need even N so that NJ = N/2 is integer.")
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    NJ = N // 2
    j = NJ / 2.0

    # Initial state
    dim = int(2 * j + 1)
    psi0 = qt.basis(dim, dim - 1)

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
    omega_coeff_local = OmegaCoeffFromPhases(phases)

    if shifted_jump_operator:
        H = [
            [-N_e, _delta_coeff],
        ]
        shifted_c_op = qt.QobjEvo([
            np.sqrt(Gamma) * Jm,
            [1j / np.sqrt(Gamma) * identity, omega_coeff_local],
        ])
        c_ops = [shifted_c_op]
    else:
        H = [
            [Jx, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(Gamma) * Jm]

    return QutipFixedNjModel(
        NJi=(NJ,),
        Gamma=Gamma,
        shifted_jump_operator=shifted_jump_operator,
        omega0=Omega0,
        delta0=delta0,
        phases=phases,
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
        raise ValueError("Need even N so that central NJ = N/2 is integer.")
    NJ = N // 2
    candidates = [
        (Nj1, NJ - Nj1)
        for Nj1 in range(max(0, NJ - N2), min(N1, NJ) + 1)
    ]
    if not candidates:
        raise ValueError("No valid two-group fixed sector exists for the given group sizes.")
    target_Nj1 = NJ * N1 / N
    return min(candidates, key=lambda pair: (abs(pair[0] - target_Nj1), pair[0]))


def build_qutip_grouped_fixed_nj_model_from_phases(
    Gamma: float,
    phases: list,
    *,
    omega_i: list[float],
    NJi: list[int],
    shifted_jump_operator: bool = False,
) -> QutipGroupedFixedNjModel:
    """
    Build a fixed grouped inhomogeneous collective model in QuTiP.
    """
    if len(omega_i) != len(NJi):
        raise ValueError("omega_i must contain exactly one coupling per group.")
    if not NJi:
        raise ValueError("NJi must contain at least one group active-atom number.")
    if len(phases) < 3:
        raise ValueError("Need at least 3 phases.")
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    dims = tuple(int(NJ_g + 1) for NJ_g in NJi)

    # Initial state
    psi0 = qt.tensor(*(qt.basis(dim, dim - 1) for dim in dims))

    Omega0 = float(phases[0].omega)
    delta0 = float(phases[1].delta)
    t_step1_end = float(phases[0].duration)
    t_step2_end = float(phases[0].duration + phases[1].duration)
    t_final = float(sum(p.duration for p in phases))

    I = tuple(qt.qeye(dim) for dim in dims)
    identity = qt.tensor(*I)
    ji = tuple(NJ_g / 2.0 for NJ_g in NJi)

    def _group_operator(j: float, op: str, group_index: int) -> qt.Qobj:
        factors = list(I)
        factors[group_index] = qt.jmat(j, op)
        return qt.tensor(*factors)

    Jp_groups = tuple(_group_operator(j, "+", group_index) for group_index, j in enumerate(ji))
    Jm_groups = tuple(_group_operator(j, "-", group_index) for group_index, j in enumerate(ji))
    Jx_groups = tuple(_group_operator(j, "x", group_index) for group_index, j in enumerate(ji))
    Jy_groups = tuple(_group_operator(j, "y", group_index) for group_index, j in enumerate(ji))
    Jz_groups = tuple(_group_operator(j, "z", group_index) for group_index, j in enumerate(ji))

    Jp = sum(Jp_groups[1:], Jp_groups[0])
    Jm = sum(Jm_groups[1:], Jm_groups[0])
    Jx = sum(Jx_groups[1:], Jx_groups[0])
    Jy = sum(Jy_groups[1:], Jy_groups[0])
    Jz = sum(Jz_groups[1:], Jz_groups[0])

    N_e_groups = tuple(Jz_g + j * identity for Jz_g, j in zip(Jz_groups, ji))
    N_e = sum(N_e_groups[1:], N_e_groups[0])

    J_drive = sum(
        (float(omega) * Jx_g for omega, Jx_g in zip(omega_i, Jx_groups)),
        0 * Jx_groups[0],
    )
    A_weighted = sum(
        (float(omega) * Jm_g for omega, Jm_g in zip(omega_i, Jm_groups)),
        0 * Jm_groups[0],
    )

    omega_coeff_local = OmegaCoeffFromPhases(phases)
    if shifted_jump_operator:
        H = [
            [-N_e, _delta_coeff],
        ]
        shifted_c_op = qt.QobjEvo([
            np.sqrt(Gamma) * A_weighted,
            [1j / np.sqrt(Gamma) * identity, omega_coeff_local],
        ])
        c_ops = [shifted_c_op]
    else:
        H = [
            [J_drive, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(Gamma) * A_weighted]

    return QutipGroupedFixedNjModel(
        NJi=tuple(int(NJ) for NJ in NJi),
        omega_i=tuple(float(omega) for omega in omega_i),
        Gamma=Gamma,
        shifted_jump_operator=shifted_jump_operator,
        omega0=Omega0,
        delta0=delta0,
        phases=phases,
        Jp=Jp,
        Jm=Jm,
        Jx=Jx,
        Jy=Jy,
        Jz=Jz,
        N_e=N_e,
        Jp_groups=Jp_groups,
        Jm_groups=Jm_groups,
        Jx_groups=Jx_groups,
        Jy_groups=Jy_groups,
        Jz_groups=Jz_groups,
        N_e_groups=N_e_groups,
        J_drive=J_drive,
        A_weighted=A_weighted,
        H=H,
        c_ops=c_ops,
        psi0=psi0,
        t_step1_end=t_step1_end,
        t_step2_end=t_step2_end,
        t_final=t_final,
    )
