from __future__ import annotations

from typing import Optional, Sequence

import numpy as np
import qutip as qt

from common.utils.parameters import omega2_from_weighted_average
from parser.qutip import QutipFixedNjModel, QutipTwoGroupFixedNjModel
from solvers.mcwf.state_helpers import centered_sector_initial_coeffs
from solvers.qutip_fixed_nj.utils_sim import OmegaCoeffFromPhases, _delta_coeff, _omega_coeff


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
        dN=0,
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
        unraveling_picture = "shifted"
    else:
        H = [
            [Jx, _omega_coeff],
            [-N_e, _delta_coeff],
        ]
        c_ops = [np.sqrt(Gamma) * Jm]
        unraveling_picture = "regular"

    dim = int(2 * j + 1)
    psi0 = qt.basis(dim, dim - 1)

    return QutipFixedNjModel(
        N=N,
        N_J=N_J,
        j=j,
        Gamma=Gamma,
        shifted_jump_operator=shifted_jump_operator,
        unraveling_picture=unraveling_picture,
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
