from __future__ import annotations

from math import comb, sqrt
from typing import Dict, Mapping, Optional

import numpy as np

from quantum_trajectories.parser import Array


SUPPORTED_SECTOR_DISTRIBUTIONS = {"square", "binomial"}


# -----------------------------------------------------------------------------
# Initial state helpers
# -----------------------------------------------------------------------------

def validate_sector_distribution(sector_distribution: str) -> str:
    """
    Validate the requested initial N_J-sector distribution.

    Supported options
    -----------------
    square
        Equal amplitudes for every included N_J sector. This preserves the
        project's original behavior exactly.
    binomial
        Weights matched to the product state
            ((|u> + |d>) / sqrt(2))^N,
        for which the probability of finding exactly N_J atoms in the active
        manifold is
            p(N_J) = binom(N, N_J) / 2^N.
        Since probabilities are amplitudes squared, the sector amplitudes are
        chosen proportional to sqrt(binom(N, N_J)) before renormalizing over
        the included sectors.
    """
    if sector_distribution not in SUPPORTED_SECTOR_DISTRIBUTIONS:
        allowed = ", ".join(sorted(SUPPORTED_SECTOR_DISTRIBUTIONS))
        raise ValueError(
            f"Unsupported sector_distribution={sector_distribution!r}. "
            f"Expected one of: {allowed}."
        )
    return sector_distribution


def down_state_in_sector(Nj: int) -> Array:
    """
    All Nj active atoms start in |down>, i.e. |n_e = 0>.

    Returns the vector (1, 0, 0, ..., 0) of shape (Nj + 1,).
    """
    psi = np.zeros(Nj + 1, dtype=np.complex128)
    psi[0] = 1.0
    return psi


def normalize_sector_coefficients(coeffs: Mapping[int, complex]) -> Dict[int, complex]:
    """
    Normalize a dictionary of sector coefficients to have unit norm.

    Each sector is assumed to carry an internally normalized state, so the full
    norm is simply sum_NJ |c_NJ|^2.
    """
    norm2 = float(sum(abs(c) ** 2 for c in coeffs.values()))
    if norm2 <= 0.0:
        raise ValueError("At least one sector coefficient must be non-zero.")
    norm = np.sqrt(norm2)
    return {Nj: c / norm for Nj, c in coeffs.items()}


def _sector_distribution_weight(N: int, Nj: int, sector_distribution: str) -> float:
    """
    Return the unnormalized amplitude weight for one N_J sector.

    The returned value is a real non-negative amplitude before any optional
    phase factor is applied.
    """
    sector_distribution = validate_sector_distribution(sector_distribution)
    if sector_distribution == "square":
        return 1.0
    # `math.comb` returns a Python integer. For larger N, using NumPy's sqrt
    # directly on that object can fail because the ufunc does not always cast
    # arbitrary-precision ints the way we want. `math.sqrt` handles the
    # combinatorial integer cleanly before we cast to float.
    return float(sqrt(comb(N, Nj)))


def centered_sector_initial_coeffs(
    N: int,
    dN: int,
    *,
    phase_fn=None,
    sector_distribution: str = "square",
) -> Dict[int, complex]:
    """
    Build a normalized superposition of N_J sectors centered around N/2.

    Parameters
    ----------
    N
        Total atom number. Assumes even N if you want the center exactly at
        N/2.
    dN
        How many sectors on each side of N/2 to include.

        dN = 0  -> {N/2}
        dN = 1  -> {N/2 - 1, N/2, N/2 + 1}
        dN = 2  -> {N/2 - 2, ..., N/2 + 2}
    phase_fn
        Optional function phase_fn(Nj) returning the phase to apply to sector
        Nj. If omitted, all included sectors have real non-negative amplitudes.
    sector_distribution
        Choice of how amplitudes are assigned across the included N_J sectors.

        "square"
            Equal amplitudes for all included sectors, normalized afterward.
            This is the historical default used by the codebase.

        "binomial"
            Amplitudes proportional to sqrt(binom(N, N_J)), corresponding to
            the product state ((|u> + |d>) / sqrt(2))^N. If the range of
            included sectors is truncated, the amplitudes are renormalized over
            only the retained sectors.

    Returns
    -------
    Dict[int, complex]
        Normalized sector coefficients.
    """
    if N < 0:
        raise ValueError("N must be non-negative.")
    if dN < 0:
        raise ValueError("dN must be >= 0.")
    if N % 2 != 0:
        raise ValueError("This helper assumes even N so the center is exactly N/2.")

    sector_distribution = validate_sector_distribution(sector_distribution)

    center = N // 2
    sector_list = list(range(center - dN, center + dN + 1))
    sector_list = [Nj for Nj in sector_list if 0 <= Nj <= N]
    if not sector_list:
        raise ValueError("No valid sectors selected.")

    coeffs: Dict[int, complex] = {}
    for Nj in sector_list:
        phase = 1.0 if phase_fn is None else np.exp(1j * phase_fn(Nj))
        coeffs[Nj] = _sector_distribution_weight(N, Nj, sector_distribution) * phase

    return normalize_sector_coefficients(coeffs)


# -----------------------------------------------------------------------------

def total_norm2(blocks: Mapping[int, Array]) -> float:
    return float(sum(np.vdot(psi, psi).real for psi in blocks.values()))


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

    Returns a dictionary whose values are the normalized local |n_e>-basis
    states multiplied by the normalized sector coefficient for that N_J.
    """
    coeffs = normalize_sector_coefficients(sector_coeffs)
    blocks: Dict[int, Array] = {}

    for Nj, coeff in coeffs.items():
        if Nj < 0 or Nj > N:
            raise ValueError(f"Invalid sector Nj={Nj} for N={N}.")

        if internal_sector_states is None or Nj not in internal_sector_states:
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

    return blocks
