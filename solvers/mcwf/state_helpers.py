from __future__ import annotations

from math import comb, sqrt
from typing import Dict, Mapping, Optional, Tuple

import numpy as np

from solvers.mcwf.operator_helpers import split_sector_key, total_active_atoms_in_sector
from parser.common import Array
from parser.mcwf import SectorKey


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


def normalize_sector_coefficients(coeffs: Mapping[SectorKey, complex]) -> Dict[SectorKey, complex]:
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


def single_group_centered_sector_initial_coeffs(
    N: int,
    dN: int,
    *,
    phase_fn=None,
    sector_distribution: str = "binomial",
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

def _valid_group_resolved_pairs_for_total_sector(
    total_nj: int,
    N1: int,
    N2: int,
) -> list[tuple[int, int]]:
    """
    Return all valid (Nj1, Nj2) pairs with Nj1 + Nj2 = total_nj.
    """
    pairs: list[tuple[int, int]] = []
    nj1_min = max(0, total_nj - N2)
    nj1_max = min(N1, total_nj)
    for Nj1 in range(nj1_min, nj1_max + 1):
        Nj2 = total_nj - Nj1
        if 0 <= Nj2 <= N2:
            pairs.append((Nj1, Nj2))
    return pairs


def _group_resolved_binomial_probability_weight(N1: int, N2: int, Nj1: int, Nj2: int) -> float:
    """
    Unnormalized probability weight for one group-resolved sector.

    For the product state ((|u> + |d>) / sqrt(2))^N this is proportional to
    binom(N1, Nj1) * binom(N2, Nj2).
    """
    return float(comb(N1, Nj1) * comb(N2, Nj2))


def two_group_centered_sector_initial_coeffs(
    N: int,
    dN: int,
    N1: int,
    N2: int,
    *,
    sector_distribution: str = "binomial",
) -> Dict[tuple[int, int], complex]:
    """
    Build a normalized superposition of two-group sectors centered around N/2.

    The returned dictionary uses the low-level inhomogeneous representation

        {(Nj1, Nj2): coeff}

    where Nj1 and Nj2 are the active-manifold atom numbers in groups 1 and 2.

    The selected total active-manifold sectors are

        N_J in {N/2 - dN, ..., N/2 + dN},

    and for each selected total N_J the helper includes every valid pair
    (Nj1, Nj2) satisfying

        Nj1 + Nj2 = N_J,
        0 <= Nj1 <= N1,
        0 <= Nj2 <= N2.

    Distribution options
    --------------------
    binomial
        Use amplitudes proportional to
            sqrt(binom(N1, Nj1) * binom(N2, Nj2)),
        then normalize over only the selected window.

    square
        Keep equal total probability for each selected total N_J sector, but
        split that total-sector probability across the valid (Nj1, Nj2) pairs
        using the conditional binomial distribution within the fixed total N_J.
    """
    if N < 0:
        raise ValueError("N must be non-negative.")
    if N1 < 0 or N2 < 0:
        raise ValueError("N1 and N2 must be non-negative.")
    if N1 + N2 != N:
        raise ValueError(f"Expected N1 + N2 = N, got N1={N1}, N2={N2}, N={N}.")
    if dN < 0:
        raise ValueError("dN must be >= 0.")
    if N % 2 != 0:
        raise ValueError("This helper assumes even N so the center is exactly N/2.")

    sector_distribution = validate_sector_distribution(sector_distribution)

    center = N // 2
    total_sector_list = list(range(center - dN, center + dN + 1))
    total_sector_list = [total_nj for total_nj in total_sector_list if 0 <= total_nj <= N]
    if not total_sector_list:
        raise ValueError("No valid total N_J sectors selected.")

    valid_pairs_by_total: dict[int, list[tuple[int, int]]] = {
        total_nj: _valid_group_resolved_pairs_for_total_sector(total_nj, N1, N2)
        for total_nj in total_sector_list
    }
    valid_pairs_by_total = {
        total_nj: pairs for total_nj, pairs in valid_pairs_by_total.items() if pairs
    }
    if not valid_pairs_by_total:
        raise ValueError("No valid group-resolved sectors found in the selected total-N_J window.")

    # Keep only one (Nj1, Nj2) tuple for each selected total N_J sector while
    # debugging single group-resolved sectors.
    valid_pairs_by_total = {
        total_nj: [
            min(
                pairs,
                key=lambda pair: abs(pair[0] - (total_nj * N1 / N)),
            )
        ]
        for total_nj, pairs in valid_pairs_by_total.items()
    }

    coeffs: Dict[tuple[int, int], complex] = {}

    if sector_distribution == "binomial":
        for pairs in valid_pairs_by_total.values():
            for Nj1, Nj2 in pairs:
                coeffs[(Nj1, Nj2)] = float(
                    sqrt(_group_resolved_binomial_probability_weight(N1, N2, Nj1, Nj2))
                )
        return normalize_sector_coefficients(coeffs)

    # "square": equal total probability over the selected total-N_J sectors,
    # then conditional binomial splitting within each total sector.
    total_sector_probability = 1.0 / len(valid_pairs_by_total)
    for total_nj, pairs in valid_pairs_by_total.items():
        conditional_weights = {
            pair: _group_resolved_binomial_probability_weight(N1, N2, pair[0], pair[1])
            for pair in pairs
        }
        conditional_norm = float(sum(conditional_weights.values()))
        if conditional_norm <= 0.0:
            raise ValueError(f"Conditional weight for total N_J={total_nj} is zero.")
        for pair, weight in conditional_weights.items():
            probability = total_sector_probability * (weight / conditional_norm)
            coeffs[pair] = float(np.sqrt(probability))

    return normalize_sector_coefficients(coeffs)


def centered_sector_initial_coeffs(
    Ni: list[int],
    dN: int,
    *,
    sector_distribution: str = "binomial",
) -> Dict[SectorKey, complex]:
    """
    This is a temprary wrapper to dispatch to the correct initial-sector coefficient helper based on the number of groups.
    This will later be replaced by a proper probability-based initial-sector coefficient generator that can handle arbitrary numbers of groups.
    """
    if not Ni:
        raise ValueError("Ni must contain at least one group size.")
    if any(Ng < 0 for Ng in Ni):
        raise ValueError("All group sizes in Ni must be non-negative.")

    if len(Ni) == 1:
        return single_group_centered_sector_initial_coeffs(
            Ni[0],
            dN=dN,
            sector_distribution=sector_distribution,
        )
    if len(Ni) == 2:
        N1, N2 = Ni
        return two_group_centered_sector_initial_coeffs(
            N1 + N2,
            dN=dN,
            N1=N1,
            N2=N2,
            sector_distribution=sector_distribution,
        )
    raise ValueError(
        f"centered_sector_initial_coeffs currently supports 1 or 2 groups, got {len(Ni)}."
    )


# -----------------------------------------------------------------------------

def total_norm2(blocks: Mapping[SectorKey, Array]) -> float:
    return float(sum(np.vdot(psi, psi).real for psi in blocks.values()))


def _two_group_down_state_in_sector(Nj1: int, Nj2: int) -> Array:
    """
    Product state |n_{e,1}=0, n_{e,2}=0> in the two-group basis.
    """
    psi = np.zeros((Nj1 + 1) * (Nj2 + 1), dtype=np.complex128)
    psi[0] = 1.0
    return psi


def _expected_internal_shape(sector_key: SectorKey) -> Tuple[int, ...]:
    groups = split_sector_key(sector_key)
    if len(groups) == 1:
        return (groups[0] + 1,)
    if len(groups) == 2:
        return ((groups[0] + 1) * (groups[1] + 1),)
    raise ValueError("Only one-group and two-group sectors are currently supported.")


def build_initial_sector_state(
    N: int,
    sector_coeffs: Mapping[SectorKey, complex],
    internal_sector_states: Optional[Mapping[SectorKey, Array]] = None,
) -> Dict[SectorKey, Array]:
    """
    Build the initial wavefunction written as a dictionary {Nj: psi_Nj}, where
    psi_Nj is the symmetric Dicke vector on the |n_e> basis for that sector.

    The total norm is
        sum_Nj ||psi_Nj||^2 = 1.

    Returns a dictionary whose values are the normalized local |n_e>-basis
    states multiplied by the normalized sector coefficient for that N_J.
    """
    coeffs = normalize_sector_coefficients(sector_coeffs)
    blocks: Dict[SectorKey, Array] = {}

    for sector_key, coeff in coeffs.items():
        groups = split_sector_key(sector_key)
        Nj_total = total_active_atoms_in_sector(sector_key)
        if Nj_total < 0 or Nj_total > N:
            raise ValueError(f"Invalid sector {sector_key} for N={N}.")

        if internal_sector_states is None or sector_key not in internal_sector_states:
            if len(groups) == 1:
                local = down_state_in_sector(groups[0])
            elif len(groups) == 2:
                local = _two_group_down_state_in_sector(groups[0], groups[1])
            else:
                raise ValueError("Only one-group and two-group sectors are currently supported.")
        else:
            local = np.asarray(internal_sector_states[sector_key], dtype=np.complex128).copy()
            expected_shape = _expected_internal_shape(sector_key)
            if local.shape != expected_shape:
                raise ValueError(
                    f"Internal state for sector {sector_key} must have shape {expected_shape}, "
                    f"got {local.shape}."
                )
            local_norm = np.linalg.norm(local)
            if local_norm == 0.0:
                raise ValueError(f"Internal state for sector {sector_key} has zero norm.")
            local /= local_norm

        blocks[sector_key] = coeff * local

    return blocks
