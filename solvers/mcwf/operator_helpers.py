from __future__ import annotations

from functools import lru_cache
from math import comb
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy.sparse import csc_matrix, diags, eye, kron

from parser.mcwf import SectorKey, SectorOperators


def is_inhomogeneous_sector_key(sector_key: SectorKey) -> bool:
    return isinstance(sector_key, tuple)


def split_sector_key(sector_key: SectorKey) -> Tuple[int, ...]:
    if isinstance(sector_key, tuple):
        return tuple(int(v) for v in sector_key)
    return (int(sector_key),)


def total_active_atoms_in_sector(sector_key: SectorKey) -> int:
    return int(sum(split_sector_key(sector_key)))


def sector_multiplicity(N: int, Nj: int) -> int:
    """
    Degeneracy of the Nj sector coming from choosing which Nj atoms live in
    {|down>, |e>} while the remaining N-Nj are in |up>.
    """
    if not (0 <= Nj <= N):
        raise ValueError(f"Nj must lie in [0, N], got Nj={Nj}, N={N}.")
    return comb(N, Nj)


def two_group_sector_multiplicity(N1: int, N2: int, Nj1: int, Nj2: int) -> int:
    """
    Degeneracy of the two-group sector (Nj1, Nj2).

    The physical groups are fixed, so this counts how many atoms are active
    inside each group separately.
    """
    if Nj1 < 0 or Nj2 < 0 or Nj1 > N1 or Nj2 > N2:
        raise ValueError(
            f"Invalid two-group sector ({Nj1}, {Nj2}) for group sizes N1={N1}, N2={N2}."
        )
    return comb(N1, Nj1) * comb(N2, Nj2)


@lru_cache(maxsize=None)
def build_sector_ops(Nj: int) -> SectorOperators:
    """
    Collective two-level operators on the permutationally symmetric Dicke basis
    |n_e>, n_e = 0, ..., Nj, for the {|down>, |e>} subsystem.
    """
    dim = Nj + 1
    ne = np.arange(dim, dtype=float)

    jplus_vals = np.sqrt((Nj - ne[:-1]) * (ne[:-1] + 1.0))
    jminus_vals = np.sqrt(ne[1:] * (Nj - ne[1:] + 1.0))

    Jp = diags(jplus_vals, offsets=-1, shape=(dim, dim), dtype=np.complex128).tocsc()
    Jm = diags(jminus_vals, offsets=+1, shape=(dim, dim), dtype=np.complex128).tocsc()

    J_x = (Jp + Jm) * 0.5
    J_y = ((Jp - Jm) / (2.0j)).tocsc()
    N_e = diags(ne, 0, shape=(dim, dim), dtype=np.complex128).tocsc()
    J_z = (N_e - 0.5 * Nj * eye(dim, format="csc", dtype=np.complex128)).tocsc()
    JpJm = (Jp @ Jm).tocsc()

    return SectorOperators(
        Nj=Nj,
        Jp=Jp,
        Jm=Jm,
        J_x=J_x,
        J_y=J_y,
        J_z=J_z,
        N_e=N_e,
        JpJm=JpJm,
        sector_key=Nj,
        Nj_groups=(Nj,),
        omega_groups=(1.0,),
        J_drive=J_x,
        A_weighted=Jm,
        AdagA_weighted=JpJm,
        N_e_groups=(N_e,),
        J_x_groups=(J_x,),
        J_y_groups=(J_y,),
        J_z_groups=(J_z,),
    )


@lru_cache(maxsize=None)
def build_two_group_sector_ops(
    Nj1: int,
    Nj2: int,
    omega1: float,
    omega2: float,
    N1: int,
    N2: int,
) -> SectorOperators:
    """
    Build cached reduced operators for one two-group active-manifold sector.

    The basis is the product Dicke basis |n_{e,1}, n_{e,2}> with dimension
    (Nj1 + 1)(Nj2 + 1). Group couplings are fixed by the ensemble layer and
    passed in here.
    """
    if Nj1 < 0 or Nj2 < 0:
        raise ValueError("Two-group sector sizes must be non-negative.")
    if Nj1 > N1 or Nj2 > N2:
        raise ValueError(
            f"Sector ({Nj1}, {Nj2}) exceeds physical group sizes N1={N1}, N2={N2}."
        )

    ops1 = build_sector_ops(Nj1)
    ops2 = build_sector_ops(Nj2)

    I1 = eye(Nj1 + 1, format="csc", dtype=np.complex128)
    I2 = eye(Nj2 + 1, format="csc", dtype=np.complex128)

    J1p = kron(ops1.Jp, I2, format="csc")
    J1m = kron(ops1.Jm, I2, format="csc")
    J1_x = kron(ops1.J_x, I2, format="csc")
    J1_y = kron(ops1.J_y, I2, format="csc")
    J1_z = kron(ops1.J_z, I2, format="csc")
    N_e1 = kron(ops1.N_e, I2, format="csc")

    J2p = kron(I1, ops2.Jp, format="csc")
    J2m = kron(I1, ops2.Jm, format="csc")
    J2_x = kron(I1, ops2.J_x, format="csc")
    J2_y = kron(I1, ops2.J_y, format="csc")
    J2_z = kron(I1, ops2.J_z, format="csc")
    N_e2 = kron(I1, ops2.N_e, format="csc")

    Jp = (J1p + J2p).tocsc()
    Jm = (J1m + J2m).tocsc()
    J_x = (J1_x + J2_x).tocsc()
    J_y = (J1_y + J2_y).tocsc()
    J_z = (J1_z + J2_z).tocsc()
    N_e = (N_e1 + N_e2).tocsc()
    A_weighted = (omega1 * J1m + omega2 * J2m).tocsc()
    J_drive = (omega1 * J1_x + omega2 * J2_x).tocsc()
    AdagA_weighted = (A_weighted.conjugate().transpose() @ A_weighted).tocsc()

    return SectorOperators(
        Nj=Nj1 + Nj2,
        Jp=Jp,
        Jm=Jm,
        J_x=J_x,
        J_y=J_y,
        J_z=J_z,
        N_e=N_e,
        JpJm=(Jp @ Jm).tocsc(),
        sector_key=(Nj1, Nj2),
        Nj_groups=(Nj1, Nj2),
        omega_groups=(omega1, omega2),
        J_drive=J_drive,
        A_weighted=A_weighted,
        AdagA_weighted=AdagA_weighted,
        N_e_groups=(N_e1, N_e2),
        J_x_groups=(J1_x, J2_x),
        J_y_groups=(J1_y, J2_y),
        J_z_groups=(J1_z, J2_z),
    )


def build_sector_ops_for_key(
    sector_key: SectorKey,
    *,
    Ni: Optional[Sequence[int]] = None,
    omega_i: Optional[Sequence[float]] = None,
) -> SectorOperators:
    if isinstance(sector_key, tuple):
        if len(sector_key) != 2:
            raise ValueError("Only two-group inhomogeneous sectors are currently supported.")
        if Ni is None or omega_i is None:
            raise ValueError("Ni and omega_i must be provided for inhomogeneous sector keys.")
        if len(Ni) != 2 or len(omega_i) != 2:
            raise ValueError("Only two-group operator construction is currently supported.")
        N1, N2 = (int(group_size) for group_size in Ni)
        omega_1, omega_2 = (float(coupling) for coupling in omega_i)
        return build_two_group_sector_ops(
            int(sector_key[0]),
            int(sector_key[1]),
            omega_1,
            omega_2,
            N1,
            N2,
        )
    return build_sector_ops(int(sector_key))
