from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries


def _build_normalized_gell_mann_generators() -> Array:
    """Return ``g_a = lambda_a / 2`` in the ``(|u>, |d>, |e>)`` basis."""
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


def _validate_orthogonal_rotation(rotation: Array, *, tol: float) -> None:
    """Require an orthogonal 8-generator rotation at every time and group."""
    identity = np.eye(8, dtype=float)
    gram = np.einsum(
        "...ac,...bc->...ab",
        rotation,
        rotation,
        optimize=True,
    )
    if not np.allclose(gram, identity, atol=tol, rtol=0.0):
        raise ValueError("Gell-Mann rotations must be orthogonal.")


def compute_gell_mann_rotation(
    transformation: Array,
    *,
    tol: float = 1e-12,
) -> Array:
    """Return ``R`` such that ``U g_alpha U^dagger = R[alpha,a] g_a``.

    ``transformation`` must have shape ``(num_times, num_groups, 3, 3)``.
    The result has shape ``(num_times, num_groups, 8, 8)`` and is real.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")

    transformation = np.asarray(transformation, dtype=np.complex128)
    if (
        transformation.ndim != 4
        or transformation.shape[2:] != (3, 3)
        or transformation.shape[0] == 0
        or transformation.shape[1] == 0
    ):
        raise ValueError(
            "transformation must have shape (num_times, num_groups, 3, 3) "
            "with nonempty time and group axes."
        )
    if not np.all(np.isfinite(transformation)):
        raise ValueError("transformation must contain only finite values.")

    identity = np.eye(3, dtype=np.complex128)
    gram = np.einsum(
        "...ji,...jk->...ik",
        transformation.conjugate(),
        transformation,
        optimize=True,
    )
    if not np.allclose(gram, identity, atol=tol, rtol=0.0):
        raise ValueError("transformation must be unitary at every time and group.")

    generators = _NORMALIZED_GELL_MANN_GENERATORS
    dressed_generators = np.einsum(
        "tgij,ajk,tglk->tgail",
        transformation,
        generators,
        transformation.conjugate(),
        optimize=True,
    )
    rotation_complex = 2.0 * np.einsum(
        "tgxij,yji->tgxy",
        dressed_generators,
        generators,
        optimize=True,
    )
    if not np.allclose(rotation_complex.imag, 0.0, atol=tol, rtol=0.0):
        raise ValueError("The Gell-Mann rotation has a non-negligible imaginary part.")

    rotation = rotation_complex.real.copy()
    _validate_orthogonal_rotation(rotation, tol=tol)
    return rotation


def rotate_generator_moments(
    moments: TrajectoryResolvedMomentSeries,
    rotation: Array,
    *,
    tol: float = 1e-12,
) -> TrajectoryResolvedMomentSeries:
    """Rotate first and symmetric second moments without averaging trajectories.

    The local rotations act independently on each physical group. Cross-group
    blocks of the full symmetric second-moment matrix are retained.
    """
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")

    rotation_complex = np.asarray(rotation, dtype=np.complex128)
    num_times = moments.t.size
    num_groups = moments.first_moments.shape[2]
    expected_shape = (num_times, num_groups, 8, 8)
    if rotation_complex.shape != expected_shape:
        raise ValueError(f"rotation must have shape {expected_shape}.")
    if not np.all(np.isfinite(rotation_complex)):
        raise ValueError("rotation must contain only finite values.")
    if not np.allclose(rotation_complex.imag, 0.0, atol=tol, rtol=0.0):
        raise ValueError("rotation must be real within tol.")
    rotation = rotation_complex.real.copy()
    _validate_orthogonal_rotation(rotation, tol=tol)

    rotated_first = np.einsum(
        "tgab,rtgb->rtga",
        rotation,
        moments.first_moments,
        optimize=True,
    )

    generator_count = 8 * num_groups
    upper = np.triu_indices(generator_count)
    full_second = np.zeros(
        (
            moments.num_trajectories,
            num_times,
            generator_count,
            generator_count,
        ),
        dtype=float,
    )
    full_second[..., upper[0], upper[1]] = moments.second_moments
    full_second[..., upper[1], upper[0]] = moments.second_moments

    global_rotation = np.zeros(
        (num_times, generator_count, generator_count),
        dtype=float,
    )
    for group_index in range(num_groups):
        group_slice = slice(8 * group_index, 8 * (group_index + 1))
        global_rotation[:, group_slice, group_slice] = rotation[:, group_index]

    rotated_second_full = np.einsum(
        "tac,rtcd,tbd->rtab",
        global_rotation,
        full_second,
        global_rotation,
        optimize=True,
    )
    rotated_second_full = 0.5 * (
        rotated_second_full
        + np.swapaxes(rotated_second_full, -1, -2)
    )
    rotated_second = rotated_second_full[..., upper[0], upper[1]].copy()

    return TrajectoryResolvedMomentSeries(
        t=moments.t.copy(),
        first_moments=rotated_first,
        second_moments=rotated_second,
    )


__all__ = ["compute_gell_mann_rotation", "rotate_generator_moments"]
