from __future__ import annotations

import heapq
from itertools import product
from math import lgamma

import numpy as np


def _normalize_sector_coefficients(
    coeffs: dict[tuple[int, ...], complex],
) -> dict[tuple[int, ...], complex]:
    """Normalize sector coefficients so sum_alpha |c_alpha|^2 = 1."""
    norm2 = float(sum(abs(coeff) ** 2 for coeff in coeffs.values()))
    if norm2 <= 0.0:
        raise ValueError("At least one sector coefficient must be non-zero.")
    norm = np.sqrt(norm2)
    return {sector: coeff / norm for sector, coeff in coeffs.items()}


def _log_binomial_table(N_groups: list[int]) -> list[list[float]]:
    """Precompute log binomial weights log(binomial(Ng, njg)) for each group."""
    tables: list[list[float]] = []
    for Ng in N_groups:
        tables.append(
            [
                lgamma(Ng + 1) - lgamma(njg + 1) - lgamma(Ng - njg + 1)
                for njg in range(Ng + 1)
            ]
        )
    return tables


def _log_tuple_probability(
    sector: tuple[int, ...],
    log_tables: list[list[float]],
) -> float:
    """Return log probability up to an additive normalization constant."""
    return float(sum(table[njg] for table, njg in zip(log_tables, sector)))


def _center_tuples(N_groups: list[int]) -> list[tuple[int, ...]]:
    """Return all mode tuples at the center of the product binomial distribution."""
    center_options = [
        [Ng // 2] if Ng % 2 == 0 else [Ng // 2, Ng // 2 + 1]
        for Ng in N_groups
    ]
    return [tuple(sector) for sector in product(*center_options)]


def _neighbor_tuples(
    sector: tuple[int, ...],
    N_groups: list[int],
) -> list[tuple[int, ...]]:
    """Return all valid tuples obtained by changing one group by +/- 1."""
    neighbors: list[tuple[int, ...]] = []
    for group_idx, (njg, Ng) in enumerate(zip(sector, N_groups)):
        if njg > 0:
            lowered = list(sector)
            lowered[group_idx] -= 1
            neighbors.append(tuple(lowered))
        if njg < Ng:
            raised = list(sector)
            raised[group_idx] += 1
            neighbors.append(tuple(raised))
    return neighbors


def _top_probability_group_sectors(
    N_groups: list[int],
    num_sectors: int,
) -> list[tuple[tuple[int, ...], float]]:
    """
    Return the exact highest-probability group-resolved sectors.

    The probability is proportional to product_g binomial(Ng, NJg). The search
    starts from all central mode tuples and expands outward by best-first
    traversal in log-probability order.
    """
    log_tables = _log_binomial_table(N_groups)

    heap: list[tuple[float, tuple[int, ...]]] = []
    queued: set[tuple[int, ...]] = set()
    accepted: set[tuple[int, ...]] = set()
    top_sectors: list[tuple[tuple[int, ...], float]] = []

    for sector in _center_tuples(N_groups):
        log_prob = _log_tuple_probability(sector, log_tables)
        heapq.heappush(heap, (-log_prob, sector))
        queued.add(sector)

    while heap and len(top_sectors) < num_sectors:
        neg_log_prob, sector = heapq.heappop(heap)
        if sector in accepted:
            continue

        log_prob = -neg_log_prob
        accepted.add(sector)
        top_sectors.append((sector, log_prob))

        for neighbor in _neighbor_tuples(sector, N_groups):
            if neighbor in queued:
                continue
            neighbor_log_prob = _log_tuple_probability(neighbor, log_tables)
            heapq.heappush(heap, (-neighbor_log_prob, neighbor))
            queued.add(neighbor)

    if len(top_sectors) != num_sectors:
        raise ValueError(
            f"Requested {num_sectors} sectors but found only {len(top_sectors)}."
        )
    return top_sectors


def centered_group_sector_coeffs(
    N_groups: list[int],
    num_sectors: int,
) -> dict[tuple[int, ...], complex]:
    """
    Build normalized coefficients for the highest-probability group sectors.

    For a product-state binomial distribution, the tuple-sector probability is

        P(NJ1, ..., NJG) proportional to product_g binomial(Ng, NJg).

    The returned amplitudes are proportional to sqrt(P), normalized over only
    the retained ``num_sectors`` highest-probability tuples.
    """
    if not N_groups:
        raise ValueError("N_groups must contain at least one group.")
    if any(Ng < 0 for Ng in N_groups):
        raise ValueError("All group sizes must be non-negative.")
    if num_sectors <= 0:
        raise ValueError("num_sectors must be >= 1.")

    total_num_sectors = 1
    for Ng in N_groups:
        total_num_sectors *= Ng + 1
    if num_sectors > total_num_sectors:
        raise ValueError(
            f"Requested num_sectors={num_sectors}, but only {total_num_sectors} "
            "group-resolved sectors exist."
        )

    top_sectors = _top_probability_group_sectors(N_groups, num_sectors)
    max_log_prob = max(log_prob for _, log_prob in top_sectors)
    coeffs = {
        sector: complex(np.exp(0.5 * (log_prob - max_log_prob)))
        for sector, log_prob in top_sectors
    }
    return _normalize_sector_coefficients(coeffs)


__all__ = [
    "centered_group_sector_coeffs",
]
