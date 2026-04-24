from __future__ import annotations

from quantum_trajectories.parser import (
    Array,
)

from typing import Dict, Mapping, Optional

import numpy as np


# -----------------------------------------------------------------------------
# Initial state helpers
# -----------------------------------------------------------------------------

def down_state_in_sector(Nj: int) -> Array:
    """
    All Nj active atoms start in |down>, i.e. |n_e = 0>.
    Returns arrar (1,0,0,...,0) of shape (Nj + 1,).
    """
    psi = np.zeros(Nj + 1, dtype=np.complex128)
    psi[0] = 1.0
    return psi


def normalize_sector_coefficients(coeffs: Mapping[int, complex]) -> Dict[int, complex]:
    """
    Normalize a dictionary of sector coefficients to have unit norm. 
    Assume each sector already has its own internal state normalized to 1, so the total norm is just the sum of |coeff|^2.
    """
    norm2 = float(sum(abs(c) ** 2 for c in coeffs.values()))
    if norm2 <= 0.0:
        raise ValueError("At least one sector coefficient must be non-zero.")
    norm = np.sqrt(norm2)
    return {Nj: c / norm for Nj, c in coeffs.items()}


def centered_sector_initial_coeffs(
    N: int,
    half_width: int,
    *,
    phase_fn=None,
) -> Dict[int, complex]:
    """
    Build a normalized superposition of Nj sectors centered around N/2.

    Parameters
    ----------
    N
        Total atom number. Assumes even N if you want the center exactly at N/2.
    half_width
        How many sectors on each side of N/2 to include.

        half_width = 0  -> {N/2}
        half_width = 1  -> {N/2 - 1, N/2, N/2 + 1}
        half_width = 2  -> {N/2 - 2, ..., N/2 + 2}
    phase_fn
        Optional function phase_fn(Nj) returning the phase to apply to sector Nj.
        If omitted, all included sectors have equal real amplitude.
    
    Center at N/2. For a binomial distribution, the width is about σ∼N/2. Good choice is NJ ​∈ [N/2−3σ,N/2+3σ]

    Returns
    -------
    Dict[int, complex]
        Normalized sector coefficients.
    """
    if N < 0:
        raise ValueError("N must be non-negative.")
    if half_width < 0:
        raise ValueError("half_width must be >= 0.")
    if N % 2 != 0:
        raise ValueError("This helper assumes even N so the center is exactly N/2.")

    center = N // 2
    sector_list = list(range(center - half_width, center + half_width + 1))

    # Keep only physical sectors
    sector_list = [Nj for Nj in sector_list if 0 <= Nj <= N]
    if not sector_list:
        raise ValueError("No valid sectors selected.")

    if phase_fn is None:
        coeffs = {Nj: 1.0 for Nj in sector_list}
    else:
        coeffs = {Nj: np.exp(1j * phase_fn(Nj)) for Nj in sector_list}

    return normalize_sector_coefficients(coeffs)

# -----------------------------------------------------------------------------

def total_norm2(blocks: Mapping[int, Array]) -> float:
    # FIXME: how is this caluclated
    return float(sum(np.vdot(psi, psi).real for psi in blocks.values()))



def down_state_in_sector(Nj: int) -> Array:
    """
    All Nj active atoms start in |down>, i.e. |n_e = 0>.
    Returns arrar (1,0,0,...,0) of shape (Nj + 1,).
    """
    psi = np.zeros(Nj + 1, dtype=np.complex128)
    psi[0] = 1.0
    return psi



def normalize_sector_coefficients(coeffs: Mapping[int, complex]) -> Dict[int, complex]:
    """
    Normalize a dictionary of sector coefficients to have unit norm. 
    Assume each sector already has its own internal state normalized to 1, so the total norm is just the sum of |coeff|^2.
    """
    norm2 = float(sum(abs(c) ** 2 for c in coeffs.values()))
    if norm2 <= 0.0:
        raise ValueError("At least one sector coefficient must be non-zero.")
    norm = np.sqrt(norm2)
    return {Nj: c / norm for Nj, c in coeffs.items()}



def build_initial_sector_state(
    N: int,
    sector_coeffs: Mapping[int, complex],
    internal_sector_states: Optional[Mapping[int, Array]] = None,
) -> Dict[int, Array]:
    """
    Build the initial wavefunction written as a dictionary {Nj: psi_Nj}, where
    psi_Nj is the symmetric Dicke vector on the |n_e> basis for that sector.

    The total norm is
        sum_Nj ||psi_Nj||^2 = 1.

    Returns a dictionary of key: sector Nj, value: normalized state vector in that sector's symmetric |n_e> basis (1,0,0,...,0), multiplied by the normalized sector coefficient.
    """
    # Normalize the sector coefficients so that the total norm is 1
    coeffs = normalize_sector_coefficients(sector_coeffs)
    blocks: Dict[int, Array] = {}

    for Nj, coeff in coeffs.items():
        if Nj < 0 or Nj > N:
            raise ValueError(f"Invalid sector Nj={Nj} for N={N}.")

        if internal_sector_states is None or Nj not in internal_sector_states:
            # Sector starts in |n_e=0> (all active atoms in |down>)
            # Returns array (1,0,0,...,0) of shape (Nj + 1,).
            local = down_state_in_sector(Nj)
        else:
            local = np.asarray(internal_sector_states[Nj], dtype=np.complex128).copy()
            if local.shape != (Nj + 1,):
                raise ValueError(
                    f"Internal state for Nj={Nj} must have shape ({Nj+1},), got {local.shape}."
                )
            local_norm = np.linalg.norm(local)
            if local_norm == 0.0:
                raise ValueError(f"Internal state for Nj={Nj} has zero norm.")
            local /= local_norm

        blocks[Nj] = coeff * local

    # blocks = {Nj: (1,0,0,...,0) * coeff for Nj, coeff in coeffs.items()}
    # (1,0,0,...,0) is the down state in each sector
    return blocks

