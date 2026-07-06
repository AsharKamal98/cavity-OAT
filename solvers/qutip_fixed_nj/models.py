from __future__ import annotations

from typing import Optional

import numpy as np
import qutip as qt

from parser.qutip import QutipFixedNjModel, QutipTwoGroupFixedNjModel
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
        N=N,
        NJ=NJ,
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


def build_qutip_two_group_fixed_nj_model_from_phases(
    N: int,
    Gamma: float,
    phases: list,
    *,
    omega_i: list[float],
    NJi: list[int],
    shifted_jump_operator: bool = False,
) -> QutipTwoGroupFixedNjModel:
    """
    Build a fixed two-group inhomogeneous collective model in QuTiP.
    """
    if len(omega_i) != 2:
        raise ValueError("omega_i must contain exactly two group couplings.")
    if len(NJi) != 2:
        raise ValueError("NJi must contain exactly two group active-atom numbers.")
    if len(phases) < 3:
        raise ValueError("Need at least 3 phases.")
    if N <= 0:
        raise ValueError("N must be positive.")
    if shifted_jump_operator and Gamma <= 0.0:
        raise ValueError(
            "shifted_jump_operator=True requires Gamma > 0 because the shifted jump "
            "operator contains omega / Gamma."
        )

    N1 = N // 2
    N2 = N - N1
    NJ1, NJ2 = int(NJi[0]), int(NJi[1])
    omega_1, omega_2 = float(omega_i[0]), float(omega_i[1])

    if NJ1 < 0 or NJ2 < 0 or NJ1 > N1 or NJ2 > N2:
        raise ValueError(
            f"Invalid fixed sector ({NJ1}, {NJ2}) for group sizes N1={N1}, N2={N2}."
        )

    NJ = int(NJ1 + NJ2)

    dim1 = int(NJ1 + 1)
    dim2 = int(NJ2 + 1)

    # Initial state
    psi0 = qt.tensor(
        qt.basis(dim1, dim1 - 1),
        qt.basis(dim2, dim2 - 1),
    )

    Omega0 = float(phases[0].omega)
    delta0 = float(phases[1].delta)
    t_step1_end = float(phases[0].duration)
    t_step2_end = float(phases[0].duration + phases[1].duration)
    t_final = float(sum(p.duration for p in phases))


    I1 = qt.qeye(dim1)
    I2 = qt.qeye(dim2)
    identity = qt.tensor(I1, I2)

    j1 = NJ1 / 2.0
    j2 = NJ2 / 2.0
    
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

    return QutipTwoGroupFixedNjModel(
        N=N,
        N1=N1,
        N2=N2,
        NJ1=int(NJ1),
        NJ2=int(NJ2),
        NJ=NJ,
        Gamma=Gamma,
        shifted_jump_operator=shifted_jump_operator,
        omega0=Omega0,
        delta0=delta0,
        phases=phases,
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
