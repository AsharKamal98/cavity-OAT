from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from itertools import product as cartesian_product
from math import comb, prod, sqrt

import numpy as np

from parser.common import Array
from parser.mcwf import (
    InitialSectorSelectionSummary,
    LinearSectorConstraints,
    SectorKey,
)


SUPPORTED_SECTOR_DISTRIBUTIONS = {"square", "binomial"}
CanonicalSectorKey = tuple[int, ...]
LinearConstraintInput = (
    LinearSectorConstraints
    | tuple[Sequence[Sequence[int]], Sequence[int]]
    | None
)


def validate_sector_distribution(sector_distribution: str) -> str:
    """Validate the requested initial-sector amplitude distribution."""
    if sector_distribution not in SUPPORTED_SECTOR_DISTRIBUTIONS:
        allowed = ", ".join(sorted(SUPPORTED_SECTOR_DISTRIBUTIONS))
        raise ValueError(
            f"Unsupported sector_distribution={sector_distribution!r}. "
            f"Expected one of: {allowed}."
        )
    return sector_distribution


def _validated_group_sizes(Ni: Sequence[int]) -> tuple[int, ...]:
    groups = tuple(int(group_size) for group_size in Ni)
    if not groups or sum(groups) <= 0:
        raise ValueError("Ni must contain at least one atom group.")
    if any(group_size < 0 for group_size in groups):
        raise ValueError("Ni must contain non-negative group sizes.")
    return groups


def _canonical_sector_key(sector_key: SectorKey) -> CanonicalSectorKey:
    if isinstance(sector_key, tuple):
        sector = tuple(int(value) for value in sector_key)
    else:
        sector = (int(sector_key),)
    if not sector or any(value < 0 for value in sector):
        raise ValueError(f"Invalid sector key: {sector_key!r}.")
    return sector


def _coerce_linear_constraints(
    linear_constraints: LinearConstraintInput,
    *,
    num_groups: int,
) -> LinearSectorConstraints | None:
    if linear_constraints is None:
        return None
    if isinstance(linear_constraints, LinearSectorConstraints):
        constraints = linear_constraints
    else:
        if len(linear_constraints) != 2:
            raise ValueError("linear_constraints must contain (matrix, offset).")
        constraints = LinearSectorConstraints(
            matrix=linear_constraints[0],
            offset=linear_constraints[1],
        )
    if any(len(row) != num_groups for row in constraints.matrix):
        raise ValueError(
            "Each linear constraint row must contain one coefficient per group."
        )
    return constraints


def centered_total_nj_values(total_N: int, dN: int) -> tuple[int, ...]:
    """Return the complete centered total-NJ window selected by ``dN``."""
    if total_N < 0 or total_N % 2 != 0:
        raise ValueError("dN selection requires a non-negative even total atom number.")
    if dN < 0 or dN > total_N // 2:
        raise ValueError("dN must lie between zero and half the total atom number.")
    center = total_N // 2
    return tuple(range(center - dN, center + dN + 1))


def exact_group_sector_weights(
    Ni: Sequence[int],
    sectors: Iterable[CanonicalSectorKey],
) -> dict[CanonicalSectorKey, int]:
    """Compute exact product-binomial integer weights for sector tuples."""
    groups = _validated_group_sizes(Ni)
    binomial_tables = tuple(
        tuple(comb(group_size, active_size) for active_size in range(group_size + 1))
        for group_size in groups
    )
    weights: dict[CanonicalSectorKey, int] = {}
    for sector in sectors:
        canonical = _canonical_sector_key(sector)
        if len(canonical) != len(groups):
            raise ValueError("Each sector must contain one entry per physical group.")
        if any(
            active_size > group_size
            for active_size, group_size in zip(canonical, groups)
        ):
            raise ValueError(f"Sector {canonical} exceeds group sizes {groups}.")
        weights[canonical] = prod(
            table[active_size]
            for table, active_size in zip(binomial_tables, canonical)
        )
    return weights


def constraint_eligible_weights(
    all_weights: Mapping[CanonicalSectorKey, int],
    linear_constraints: LinearConstraintInput,
) -> dict[CanonicalSectorKey, int]:
    """Restrict exact sector weights to integer solutions of ``A @ sector = b``."""
    if not all_weights:
        raise ValueError("all_weights must contain at least one sector.")
    num_groups = len(next(iter(all_weights)))
    constraints = _coerce_linear_constraints(
        linear_constraints,
        num_groups=num_groups,
    )
    if constraints is None:
        return dict(all_weights)

    eligible = {
        sector: weight
        for sector, weight in all_weights.items()
        if all(
            sum(coefficient * value for coefficient, value in zip(row, sector))
            == offset
            for row, offset in zip(constraints.matrix, constraints.offset)
        )
    }
    if not eligible:
        raise ValueError("linear_constraints leave no eligible sectors.")
    return eligible


def _joint_probability_weights(
    eligible_weights: Mapping[CanonicalSectorKey, int],
    joint_probability: float,
    *,
    full_weight_sum: int,
) -> dict[CanonicalSectorKey, int]:
    """Retain the largest complete weight groups without exceeding the cap."""
    if not np.isfinite(joint_probability) or not 0.0 < joint_probability <= 1.0:
        raise ValueError("joint_probability must lie in (0, 1].")

    eligible_weight_sum = sum(eligible_weights.values())
    if eligible_weight_sum / full_weight_sum <= joint_probability:
        return dict(eligible_weights)

    ordered = sorted(
        eligible_weights.items(),
        key=lambda item: (-item[1], item[0]),
    )
    retained: dict[CanonicalSectorKey, int] = {}
    retained_weight_sum = 0
    index = 0
    while index < len(ordered):
        cutoff_weight = ordered[index][1]
        group_end = index
        while group_end < len(ordered) and ordered[group_end][1] == cutoff_weight:
            group_end += 1
        weight_group = ordered[index:group_end]
        group_weight_sum = sum(weight for _, weight in weight_group)
        candidate_probability = (
            retained_weight_sum + group_weight_sum
        ) / full_weight_sum
        if candidate_probability > joint_probability:
            break
        for sector, weight in weight_group:
            retained[sector] = weight
        retained_weight_sum += group_weight_sum
        index = group_end

    if not retained:
        highest_tie_probability = sum(
            weight for _, weight in ordered if weight == ordered[0][1]
        ) / full_weight_sum
        raise ValueError(
            "joint_probability is too small to retain the complete "
            "highest-weight tie group: "
            f"cap={joint_probability:.6g}, "
            f"tie-group probability={highest_tie_probability:.6g}."
        )
    return retained


def select_initial_sector_weights(
    Ni: Sequence[int],
    *,
    dN: int | None = None,
    joint_probability: float | None = None,
    linear_constraints: LinearConstraintInput = None,
) -> tuple[dict[CanonicalSectorKey, int], dict[CanonicalSectorKey, int]]:
    """Select exact group-sector weights using one supported cutoff method."""
    groups = _validated_group_sizes(Ni)
    if (dN is None) == (joint_probability is None):
        raise ValueError("Exactly one of dN and joint_probability must be provided.")

    all_sectors = cartesian_product(
        *(range(group_size + 1) for group_size in groups)
    )
    all_weights = exact_group_sector_weights(groups, all_sectors)
    eligible_weights = constraint_eligible_weights(
        all_weights,
        linear_constraints,
    )

    if dN is not None:
        totals = set(centered_total_nj_values(sum(groups), dN))
        retained_weights = {
            sector: weight
            for sector, weight in eligible_weights.items()
            if sum(sector) in totals
        }
    else:
        retained_weights = _joint_probability_weights(
            eligible_weights,
            float(joint_probability),
            full_weight_sum=sum(all_weights.values()),
        )

    if not retained_weights:
        raise ValueError("The requested initial-sector selection retained no sectors.")
    return retained_weights, eligible_weights


def compute_initial_sector_summary(
    Ni: Sequence[int],
    selected_weights: Mapping[CanonicalSectorKey, int],
    eligible_weights: Mapping[CanonicalSectorKey, int],
) -> InitialSectorSelectionSummary:
    """Compute counts and full-wavefunction weights before renormalization."""
    groups = _validated_group_sizes(Ni)
    total_sector_count = prod(group_size + 1 for group_size in groups)
    eligible_sector_count = len(eligible_weights)
    retained_sector_count = len(selected_weights)
    full_weight_sum = 2 ** sum(groups)
    return InitialSectorSelectionSummary(
        retained_sector_count=retained_sector_count,
        eligible_sector_count=eligible_sector_count,
        total_sector_count=total_sector_count,
        retained_percentage_of_eligible=(
            100.0 * retained_sector_count / eligible_sector_count
        ),
        retained_percentage_of_full_space=(
            100.0 * retained_sector_count / total_sector_count
        ),
        eligible_joint_probability=(
            sum(eligible_weights.values()) / full_weight_sum
        ),
        retained_joint_probability=(
            sum(selected_weights.values()) / full_weight_sum
        ),
    )


def print_initial_sector_summary(summary: InitialSectorSelectionSummary) -> None:
    """Print the required once-per-ensemble sector-selection summary."""
    if summary.eligible_sector_count == summary.total_sector_count:
        print(
            "Initial sector selection: "
            f"{summary.retained_sector_count}/{summary.total_sector_count} "
            "sectors retained "
            f"({summary.retained_percentage_of_full_space:.2f}%); "
            "full-wavefunction probability retained: "
            f"{100.0 * summary.retained_joint_probability:.2f}%."
        )
        return

    retained_fraction_of_eligible = (
        summary.retained_joint_probability / summary.eligible_joint_probability
    )
    print(
        "Initial sector selection: "
        f"{summary.retained_sector_count}/{summary.eligible_sector_count} "
        "eligible sectors retained "
        f"({summary.retained_percentage_of_eligible:.2f}%; "
        f"{summary.retained_sector_count}/{summary.total_sector_count} "
        "of full space, "
        f"{summary.retained_percentage_of_full_space:.2f}%); "
        "full-wavefunction probability: "
        f"{100.0 * summary.retained_joint_probability:.2f}% retained from "
        f"{100.0 * summary.eligible_joint_probability:.2f}% eligible "
        f"({100.0 * retained_fraction_of_eligible:.2f}% of eligible mass)."
    )


def normalize_sector_coefficients(
    coeffs: Mapping[SectorKey, complex],
) -> dict[CanonicalSectorKey, complex]:
    """Canonicalize tuple keys and normalize coefficients to unit norm."""
    canonical: dict[CanonicalSectorKey, complex] = {}
    for sector_key, coefficient in coeffs.items():
        sector = _canonical_sector_key(sector_key)
        if sector in canonical:
            raise ValueError(f"Duplicate sector after key normalization: {sector}.")
        canonical[sector] = complex(coefficient)
    norm2 = float(sum(abs(coefficient) ** 2 for coefficient in canonical.values()))
    if not np.isfinite(norm2) or norm2 <= 0.0:
        raise ValueError("At least one finite sector coefficient must be non-zero.")
    norm = sqrt(norm2)
    return {
        sector: coefficient / norm
        for sector, coefficient in canonical.items()
    }


def _coefficients_from_selected_weights(
    selected_weights: Mapping[CanonicalSectorKey, int],
    *,
    sector_distribution: str,
) -> dict[CanonicalSectorKey, complex]:
    """Convert exact selected weights into normalized initial amplitudes."""
    sector_distribution = validate_sector_distribution(sector_distribution)
    if sector_distribution == "binomial":
        total_weight = sum(selected_weights.values())
        return {
            sector: complex(sqrt(weight / total_weight))
            for sector, weight in selected_weights.items()
        }

    weights_by_total: dict[int, dict[CanonicalSectorKey, int]] = {}
    for sector, weight in selected_weights.items():
        weights_by_total.setdefault(sum(sector), {})[sector] = weight
    total_sector_probability = 1.0 / len(weights_by_total)
    coefficients: dict[CanonicalSectorKey, complex] = {}
    for sector_weights in weights_by_total.values():
        conditional_weight = sum(sector_weights.values())
        for sector, weight in sector_weights.items():
            probability = total_sector_probability * weight / conditional_weight
            coefficients[sector] = complex(sqrt(probability))
    return coefficients


def group_generic_initial_sector_coeffs(
    Ni: Sequence[int],
    *,
    dN: int | None = None,
    joint_probability: float | None = None,
    linear_constraints: LinearConstraintInput = None,
    sector_distribution: str = "binomial",
) -> tuple[dict[CanonicalSectorKey, complex], InitialSectorSelectionSummary]:
    """Build group-generic initial amplitudes and their truncation summary."""
    if joint_probability is not None and sector_distribution != "binomial":
        raise ValueError(
            "joint_probability selection requires sector_distribution='binomial'."
        )
    selected_weights, eligible_weights = select_initial_sector_weights(
        Ni,
        dN=dN,
        joint_probability=joint_probability,
        linear_constraints=linear_constraints,
    )
    summary = compute_initial_sector_summary(
        Ni,
        selected_weights,
        eligible_weights,
    )
    coefficients = _coefficients_from_selected_weights(
        selected_weights,
        sector_distribution=sector_distribution,
    )
    return normalize_sector_coefficients(coefficients), summary


def _expected_internal_shape(sector_key: SectorKey) -> tuple[int, ...]:
    sector = _canonical_sector_key(sector_key)
    return (prod(active_size + 1 for active_size in sector),)


def build_initial_sector_state(
    sector_coeffs: Mapping[SectorKey, complex],
    internal_sector_states: Mapping[SectorKey, Array] | None = None,
) -> dict[CanonicalSectorKey, Array]:
    """Build normalized tuple-keyed direct-sum blocks for the initial state."""
    coefficients = normalize_sector_coefficients(sector_coeffs)
    supplied_states: dict[CanonicalSectorKey, Array] = {}
    if internal_sector_states is not None:
        for sector_key, state in internal_sector_states.items():
            sector = _canonical_sector_key(sector_key)
            if sector in supplied_states:
                raise ValueError(
                    f"Duplicate internal state after key normalization: {sector}."
                )
            supplied_states[sector] = np.asarray(
                state,
                dtype=np.complex128,
            ).copy()

    blocks: dict[CanonicalSectorKey, Array] = {}
    for sector, coefficient in coefficients.items():
        expected_shape = _expected_internal_shape(sector)
        if sector not in supplied_states:
            local_state = np.zeros(expected_shape, dtype=np.complex128)
            local_state[0] = 1.0
        else:
            local_state = supplied_states[sector]
            if local_state.shape != expected_shape:
                raise ValueError(
                    f"Internal state for sector {sector} must have shape "
                    f"{expected_shape}, got {local_state.shape}."
                )
            local_norm = float(np.linalg.norm(local_state))
            if not np.isfinite(local_norm) or local_norm <= 0.0:
                raise ValueError(
                    f"Internal state for sector {sector} must have finite non-zero norm."
                )
            local_state /= local_norm
        blocks[sector] = coefficient * local_state
    return blocks


def total_norm2(blocks: Mapping[SectorKey, Array]) -> float:
    """Return the norm squared of a direct-sum sector state."""
    return float(sum(np.vdot(block, block).real for block in blocks.values()))


__all__ = [
    "build_initial_sector_state",
    "centered_total_nj_values",
    "compute_initial_sector_summary",
    "constraint_eligible_weights",
    "exact_group_sector_weights",
    "group_generic_initial_sector_coeffs",
    "normalize_sector_coefficients",
    "print_initial_sector_summary",
    "select_initial_sector_weights",
    "total_norm2",
    "validate_sector_distribution",
]
