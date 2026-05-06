from quantum_trajectories.parser import (
    SectorOperators,
)

import numpy as np
from math import comb
from scipy.sparse import diags


# -----------------------------------------------------------------------------
# Basic sector construction
# -----------------------------------------------------------------------------

def sector_multiplicity(N: int, Nj: int) -> int:
    """
    Degeneracy of the Nj sector coming from choosing which Nj atoms live in
    {|down>, |e>} while the remaining N-Nj are in |up>.
    """
    if not (0 <= Nj <= N):
        raise ValueError(f"Nj must lie in [0, N], got Nj={Nj}, N={N}.")
    return comb(N, Nj)


def build_sector_ops(Nj: int) -> SectorOperators:
    """
    Collective two-level operators on the permutationally symmetric Dicke basis
    |n_e>, n_e = 0, ..., Nj, for the {|down>, |e>} subsystem.
    """
    dim = Nj + 1
    ne = np.arange(dim, dtype=float)

    # J_+ |n_e> = sqrt[(Nj - n_e)(n_e + 1)] |n_e + 1>
    jplus_vals = np.sqrt((Nj - ne[:-1]) * (ne[:-1] + 1.0))

    # J_- |n_e> = sqrt[n_e (Nj - n_e + 1)] |n_e - 1>
    jminus_vals = np.sqrt(ne[1:] * (Nj - ne[1:] + 1.0))

    # IMPORTANT: offsets are opposite of what you had
    J_plus = diags(jplus_vals, offsets=-1, shape=(dim, dim), dtype=np.complex128).tocsc()
    J_minus = diags(jminus_vals, offsets=+1, shape=(dim, dim), dtype=np.complex128).tocsc()

    J_x = (J_plus + J_minus) * 0.5
    J_y = ((J_plus - J_minus) / (2.0j)).tocsc()
    N_e = diags(ne, 0, shape=(dim, dim), dtype=np.complex128).tocsc()
    JpJm = (J_plus @ J_minus).tocsc()

    return SectorOperators(
        Nj=Nj,
        J_plus=J_plus,
        J_minus=J_minus,
        J_x=J_x,
        J_y=J_y,
        N_e=N_e,
        JpJm=JpJm,
    )
