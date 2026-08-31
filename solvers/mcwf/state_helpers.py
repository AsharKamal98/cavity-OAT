from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from itertools import product as cartesian_product
from math import comb, floor, prod, sqrt

import numpy as np

from parser.common import Array
from parser.mcwf import (
    InitialSectorFilterSummary,
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


def filter_sector_weights_by_total_nj(
    candidate_weights: Mapping[CanonicalSectorKey, int],
    *,
    total_N: int,
    dN: int,
) -> dict[CanonicalSectorKey, int]:
    """Retain complete total-NJ sectors in the centered ``dN`` window."""
    totals = set(centered_total_nj_values(total_N, dN))
    return {
        sector: weight
        for sector, weight in candidate_weights.items()
        if sum(sector) in totals
    }


def minimum_total_nj_for_omega_ratio(
    *,
    max_abs_omega: float,
    Gamma: float,
    max_omega_over_omega_c: float,
) -> int:
    """Return the smallest total NJ satisfying the strict drive-ratio limit."""
    if not np.isfinite(max_abs_omega) or max_abs_omega < 0.0:
        raise ValueError("max_abs_omega must be finite and non-negative.")
    if not np.isfinite(Gamma) or Gamma <= 0.0:
        raise ValueError("Gamma must be finite and positive.")
    if (
        not np.isfinite(max_omega_over_omega_c)
        or not 0.0 < max_omega_over_omega_c <= 1.0
    ):
        raise ValueError("max_omega_over_omega_c must lie in (0, 1].")
    threshold = 2.0 * max_abs_omega / (Gamma * max_omega_over_omega_c)
    return floor(threshold) + 1


def filter_sector_weights_by_omega_ratio(
    candidate_weights: Mapping[CanonicalSectorKey, int],
    *,
    max_abs_omega: float,
    Gamma: float,
    max_omega_over_omega_c: float,
) -> dict[CanonicalSectorKey, int]:
    """Remove total-NJ sectors that violate the strict drive-ratio limit."""
    minimum_total_nj = minimum_total_nj_for_omega_ratio(
        max_abs_omega=max_abs_omega,
        Gamma=Gamma,
        max_omega_over_omega_c=max_omega_over_omega_c,
    )
    return {
        sector: weight
        for sector, weight in candidate_weights.items()
        if sum(sector) >= minimum_total_nj
    }


def filter_sector_weights_by_joint_probability(
    candidate_weights: Mapping[CanonicalSectorKey, int],
    joint_probability: float,
    *,
    full_weight_sum: int,
) -> dict[CanonicalSectorKey, int]:
    """Apply the final full-wavefunction probability cap to individual tuples."""
    if not np.isfinite(joint_probability) or not 0.0 < joint_probability <= 1.0:
        raise ValueError("joint_probability must lie in (0, 1].")
    if not candidate_weights:
        raise ValueError(
            "joint_probability cannot filter an empty candidate-sector set."
        )

    candidate_weight_sum = sum(candidate_weights.values())
    if candidate_weight_sum / full_weight_sum <= joint_probability:
        return dict(candidate_weights)

    ordered = sorted(
        candidate_weights.items(),
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


def _filter_summary(
    filter_name: str,
    input_weights: Mapping[CanonicalSectorKey, int],
    remaining_weights: Mapping[CanonicalSectorKey, int],
) -> InitialSectorFilterSummary:
    """Describe the sector-count change from one applied filter."""
    return InitialSectorFilterSummary(
        filter_name=filter_name,
        input_sector_count=len(input_weights),
        remaining_sector_count=len(remaining_weights),
    )


def _select_initial_sector_weights_with_steps(
    Ni: Sequence[int],
    *,
    dN: int | None = None,
    max_omega_over_omega_c: float | None = None,
    max_abs_omega: float | None = None,
    Gamma: float | None = None,
    linear_constraints: LinearConstraintInput = None,
    joint_probability: float | None = None,
) -> tuple[
    dict[CanonicalSectorKey, int],
    dict[CanonicalSectorKey, int],
    tuple[InitialSectorFilterSummary, ...],
]:
    """Return filtered weights and count summaries from the ordered pipeline."""
    groups = _validated_group_sizes(Ni)

    all_sectors = cartesian_product(
        *(range(group_size + 1) for group_size in groups)
    )
    all_weights = exact_group_sector_weights(groups, all_sectors)

    # Apply each optional support filter once, in its defined order.
    retained_weights = all_weights
    filter_steps: list[InitialSectorFilterSummary] = []

    if dN is not None:
        input_weights = retained_weights
        retained_weights = filter_sector_weights_by_total_nj(
            retained_weights,
            total_N=sum(groups),
            dN=dN,
        )
        if not retained_weights:
            raise ValueError("The dN filter retained no initial sectors.")
        filter_steps.append(
            _filter_summary("dN filter", input_weights, retained_weights)
        )

    if max_omega_over_omega_c is not None:
        if max_abs_omega is None or Gamma is None:
            raise ValueError(
                "max_abs_omega and Gamma are required when "
                "max_omega_over_omega_c is provided."
            )
        input_weights = retained_weights
        retained_weights = filter_sector_weights_by_omega_ratio(
            retained_weights,
            max_abs_omega=max_abs_omega,
            Gamma=Gamma,
            max_omega_over_omega_c=max_omega_over_omega_c,
        )
        if not retained_weights:
            raise ValueError(
                "The max_omega_over_omega_c filter retained no initial sectors."
            )
        filter_steps.append(
            _filter_summary(
                "Omega/Omega_c filter",
                input_weights,
                retained_weights,
            )
        )

    if linear_constraints is not None:
        input_weights = retained_weights
        retained_weights = constraint_eligible_weights(
            retained_weights,
            linear_constraints,
        )
        filter_steps.append(
            _filter_summary("NJ-ratio filter", input_weights, retained_weights)
        )

    pre_probability_weights = retained_weights

    if joint_probability is not None:
        input_weights = retained_weights
        retained_weights = filter_sector_weights_by_joint_probability(
            retained_weights,
            joint_probability,
            full_weight_sum=sum(all_weights.values()),
        )
        filter_steps.append(
            _filter_summary(
                "Joint-probability filter",
                input_weights,
                retained_weights,
            )
        )

    if not retained_weights:
        raise ValueError("The requested initial-sector selection retained no sectors.")
    return retained_weights, pre_probability_weights, tuple(filter_steps)


def select_initial_sector_weights(
    Ni: Sequence[int],
    *,
    dN: int | None = None,
    max_omega_over_omega_c: float | None = None,
    max_abs_omega: float | None = None,
    Gamma: float | None = None,
    linear_constraints: LinearConstraintInput = None,
    joint_probability: float | None = None,
) -> tuple[dict[CanonicalSectorKey, int], dict[CanonicalSectorKey, int]]:
    """Return final and pre-probability weights from the ordered filters."""
    selected_weights, pre_probability_weights, _ = (
        _select_initial_sector_weights_with_steps(
            Ni,
            dN=dN,
            max_omega_over_omega_c=max_omega_over_omega_c,
            max_abs_omega=max_abs_omega,
            Gamma=Gamma,
            linear_constraints=linear_constraints,
            joint_probability=joint_probability,
        )
    )
    return selected_weights, pre_probability_weights


def compute_initial_sector_summary(
    Ni: Sequence[int],
    selected_weights: Mapping[CanonicalSectorKey, int],
    eligible_weights: Mapping[CanonicalSectorKey, int],
    filter_steps: Sequence[InitialSectorFilterSummary] = (),
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
        filter_steps=tuple(filter_steps),
    )


def print_initial_sector_summary(summary: InitialSectorSelectionSummary) -> None:
    """Print aligned sector-count retention for every applied filter."""
    print(
        "Initial sector construction: "
        f"{summary.total_sector_count:,} sectors."
    )
    if not summary.filter_steps:
        return

    row_format = "{:<26}   {:>17}   {:>19}   {:>21}"
    header = row_format.format(
        "Filter",
        "Sectors remaining",
        "Retained from input",
        "Retained from initial",
    )
    print()
    print(header)
    print("-" * len(header))
    for step in summary.filter_steps:
        retained_from_input = (
            100.0 * step.remaining_sector_count / step.input_sector_count
        )
        retained_from_initial = (
            100.0
            * step.remaining_sector_count
            / summary.total_sector_count
        )
        print(
            row_format.format(
                step.filter_name,
                f"{step.remaining_sector_count:,}",
                f"{retained_from_input:.2f}%",
                f"{retained_from_initial:.2f}%",
            )
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
    max_omega_over_omega_c: float | None = None,
    max_abs_omega: float | None = None,
    Gamma: float | None = None,
    linear_constraints: LinearConstraintInput = None,
    joint_probability: float | None = None,
    sector_distribution: str = "binomial",
) -> tuple[dict[CanonicalSectorKey, complex], InitialSectorSelectionSummary]:
    """Build group-generic initial amplitudes and their truncation summary."""
    if joint_probability is not None and sector_distribution != "binomial":
        raise ValueError(
            "joint_probability filtering requires sector_distribution='binomial'."
        )
    selected_weights, pre_probability_weights, filter_steps = (
        _select_initial_sector_weights_with_steps(
            Ni,
            dN=dN,
            max_omega_over_omega_c=max_omega_over_omega_c,
            max_abs_omega=max_abs_omega,
            Gamma=Gamma,
            linear_constraints=linear_constraints,
            joint_probability=joint_probability,
        )
    )
    summary = compute_initial_sector_summary(
        Ni,
        selected_weights,
        pre_probability_weights,
        filter_steps,
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
    "filter_sector_weights_by_joint_probability",
    "filter_sector_weights_by_omega_ratio",
    "filter_sector_weights_by_total_nj",
    "group_generic_initial_sector_coeffs",
    "minimum_total_nj_for_omega_ratio",
    "normalize_sector_coefficients",
    "print_initial_sector_summary",
    "select_initial_sector_weights",
    "total_norm2",
    "validate_sector_distribution",
]
