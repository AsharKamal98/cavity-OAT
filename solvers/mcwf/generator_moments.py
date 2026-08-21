from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

from parser.common import Array


SectorLabel = int | tuple[int, ...]
CanonicalSectorLabel = tuple[int, ...]
BlockState = Mapping[CanonicalSectorLabel, Array]


def _build_normalized_gell_mann_generators() -> Array:
    """Return g_a = lambda_a / 2 in the (|u>, |d>, |e>) basis."""
    generators = np.zeros((8, 3, 3), dtype=np.complex128)

    generators[0, 0, 1] = generators[0, 1, 0] = 0.5
    generators[1, 0, 1] = -0.5j
    generators[1, 1, 0] = 0.5j
    generators[2, 0, 0] = 0.5
    generators[2, 1, 1] = -0.5

    generators[3, 0, 2] = generators[3, 2, 0] = 0.5
    generators[4, 0, 2] = -0.5j
    generators[4, 2, 0] = 0.5j
    generators[5, 1, 2] = generators[5, 2, 1] = 0.5
    generators[6, 1, 2] = -0.5j
    generators[6, 2, 1] = 0.5j

    generators[7, 0, 0] = 1.0 / (2.0 * np.sqrt(3.0))
    generators[7, 1, 1] = 1.0 / (2.0 * np.sqrt(3.0))
    generators[7, 2, 2] = -1.0 / np.sqrt(3.0)
    return generators


_NORMALIZED_GELL_MANN_GENERATORS = _build_normalized_gell_mann_generators()
_LEVEL_LABELS = ("u", "d", "e")


def _canonicalize_blocks(
    sector_blocks: Mapping[SectorLabel, Array],
    N_groups: tuple[int, ...],
) -> dict[CanonicalSectorLabel, Array]:
    """Validate sector labels and return tuple-keyed, flat block arrays."""
    group_count = len(N_groups)
    canonical: dict[CanonicalSectorLabel, Array] = {}

    for sector_label, amplitudes in sector_blocks.items():
        if isinstance(sector_label, tuple):
            sector = tuple(int(value) for value in sector_label)
        else:
            sector = (int(sector_label),)

        if len(sector) != group_count:
            raise ValueError(
                "Each sector label must contain one active-manifold size per "
                "physical group."
            )
        if any(Nj < 0 or Nj > N_group for Nj, N_group in zip(sector, N_groups)):
            raise ValueError(
                f"Sector {sector} is incompatible with group sizes {N_groups}."
            )
        if sector in canonical:
            raise ValueError(f"Duplicate sector label after normalization: {sector}.")

        block = np.asarray(amplitudes, dtype=np.complex128).reshape(-1)
        expected_size = int(np.prod([Nj + 1 for Nj in sector], dtype=int))
        if block.size != expected_size:
            raise ValueError(
                f"Sector {sector} has block size {block.size}; expected "
                f"{expected_size}."
            )
        canonical[sector] = block

    if not canonical:
        raise ValueError("sector_blocks must contain at least one sector block.")
    return canonical


def _accumulate_block(
    output: dict[CanonicalSectorLabel, Array],
    sector: CanonicalSectorLabel,
    contribution: Array,
) -> None:
    """Accumulate one flattened contribution into a sector-block state."""
    flattened = np.asarray(contribution, dtype=np.complex128).reshape(-1)
    if sector in output:
        output[sector] = output[sector] + flattened
    else:
        output[sector] = flattened.copy()


def _replace_group_sector(
    sector: CanonicalSectorLabel,
    group_index: int,
    Nj: int,
) -> CanonicalSectorLabel:
    values = list(sector)
    values[group_index] = Nj
    return tuple(values)


def _apply_collective_transition(
    sector_blocks: BlockState,
    N_groups: tuple[int, ...],
    group_index: int,
    target_level: str,
    source_level: str,
) -> dict[CanonicalSectorLabel, Array]:
    """Apply sum_i |target><source| for one physical group."""
    output: dict[CanonicalSectorLabel, Array] = {}

    for sector, amplitudes in sector_blocks.items():
        Nj = sector[group_index]
        n_u = N_groups[group_index] - Nj
        dimensions = tuple(value + 1 for value in sector)
        tensor = amplitudes.reshape(dimensions)
        front = np.moveaxis(tensor, group_index, 0)
        broadcast_tail = (1,) * (front.ndim - 1)
        n_e = np.arange(Nj + 1, dtype=float).reshape(
            (Nj + 1,) + broadcast_tail
        )
        target_sector = sector

        if target_level == source_level == "u":
            target = n_u * front
        elif target_level == source_level == "d":
            target = (Nj - n_e) * front
        elif target_level == source_level == "e":
            target = n_e * front
        elif target_level == "e" and source_level == "d":
            target = np.zeros_like(front)
            factors = np.sqrt((n_e[:-1] + 1.0) * (Nj - n_e[:-1]))
            target[1:] = factors * front[:-1]
        elif target_level == "d" and source_level == "e":
            target = np.zeros_like(front)
            factors = np.sqrt(n_e[1:] * (Nj - n_e[1:] + 1.0))
            target[:-1] = factors * front[1:]
        elif target_level == "u" and source_level == "d":
            if Nj == 0:
                continue
            factors = np.sqrt((n_u + 1.0) * (Nj - n_e[:-1]))
            target = factors * front[:-1]
            target_sector = _replace_group_sector(sector, group_index, Nj - 1)
        elif target_level == "u" and source_level == "e":
            if Nj == 0:
                continue
            factors = np.sqrt((n_u + 1.0) * n_e[1:])
            target = factors * front[1:]
            target_sector = _replace_group_sector(sector, group_index, Nj - 1)
        elif target_level == "d" and source_level == "u":
            if n_u == 0:
                continue
            target = np.zeros(
                (Nj + 2,) + front.shape[1:],
                dtype=np.complex128,
            )
            factors = np.sqrt(n_u * (Nj - n_e + 1.0))
            target[:-1] = factors * front
            target_sector = _replace_group_sector(sector, group_index, Nj + 1)
        elif target_level == "e" and source_level == "u":
            if n_u == 0:
                continue
            target = np.zeros(
                (Nj + 2,) + front.shape[1:],
                dtype=np.complex128,
            )
            factors = np.sqrt(n_u * (n_e + 1.0))
            target[1:] = factors * front
            target_sector = _replace_group_sector(sector, group_index, Nj + 1)
        else:
            raise ValueError(
                f"Unsupported collective transition {target_level}<-{source_level}."
            )

        contribution = np.moveaxis(target, 0, group_index)
        _accumulate_block(output, target_sector, contribution)

    return output


def _apply_collective_one_body_operator(
    sector_blocks: BlockState,
    N_groups: tuple[int, ...],
    group_index: int,
    operator: Array,
) -> dict[CanonicalSectorLabel, Array]:
    """Apply the collective lift of one 3x3 single-particle operator."""
    operator = np.asarray(operator, dtype=np.complex128)
    if operator.shape != (3, 3):
        raise ValueError(f"Expected a 3x3 operator, got {operator.shape}.")

    output: dict[CanonicalSectorLabel, Array] = {}
    for row, target_level in enumerate(_LEVEL_LABELS):
        for column, source_level in enumerate(_LEVEL_LABELS):
            coefficient = operator[row, column]
            if coefficient == 0.0:
                continue
            transitioned = _apply_collective_transition(
                sector_blocks,
                N_groups,
                group_index,
                target_level,
                source_level,
            )
            for sector, contribution in transitioned.items():
                _accumulate_block(output, sector, coefficient * contribution)
    return output


def _block_inner_product(left: BlockState, right: BlockState) -> complex:
    """Return the inner product of two direct-sum sector-block states."""
    return sum(
        (
            np.vdot(left[sector], right[sector])
            for sector in left
            if sector in right
        ),
        0.0j,
    )


def compute_snapshot_generator_moments(
    sector_blocks: Mapping[SectorLabel, Array],
    Ni: Sequence[int],
) -> tuple[Array, Array]:
    """Compute normalized fixed-basis SU(3) moments for one MCWF snapshot.

    The first returned array has shape ``(num_groups, 8)``. The second stores
    the row-major ``np.triu_indices(8 * num_groups)`` packing of the full
    symmetric second-moment matrix. Both arrays use the complete direct-sum
    state norm, and the second moments include cross-group entries.
    """
    N_groups = tuple(int(group_size) for group_size in Ni)
    if not N_groups:
        raise ValueError("Ni must contain at least one physical group.")
    if any(group_size < 0 for group_size in N_groups):
        raise ValueError("Physical group sizes in Ni must be non-negative.")

    blocks = _canonicalize_blocks(sector_blocks, N_groups)
    norm2 = float(np.real(_block_inner_product(blocks, blocks)))
    if not np.isfinite(norm2) or norm2 <= 0.0:
        raise ValueError("The complete snapshot norm must be finite and positive.")

    group_count = len(N_groups)
    first_moments = np.empty((group_count, 8), dtype=float)
    applied_generators: list[dict[CanonicalSectorLabel, Array]] = []

    for group_index in range(group_count):
        for generator_index, generator in enumerate(
            _NORMALIZED_GELL_MANN_GENERATORS
        ):
            applied = _apply_collective_one_body_operator(
                blocks,
                N_groups,
                group_index,
                generator,
            )
            applied_generators.append(applied)
            first_moments[group_index, generator_index] = float(
                np.real(_block_inner_product(blocks, applied)) / norm2
            )

    generator_count = 8 * group_count
    upper_rows, upper_columns = np.triu_indices(generator_count)
    second_moments = np.empty(upper_rows.size, dtype=float)
    for packed_index, (row, column) in enumerate(
        zip(upper_rows, upper_columns)
    ):
        second_moments[packed_index] = float(
            np.real(
                _block_inner_product(
                    applied_generators[int(row)],
                    applied_generators[int(column)],
                )
            )
            / norm2
        )

    return first_moments, second_moments


__all__ = ["compute_snapshot_generator_moments"]
