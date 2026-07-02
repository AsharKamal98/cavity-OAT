from __future__ import annotations

from typing import Tuple

import numpy as np

from parser.common import Array


def norm_spin_components(
    x: Array,
    y: Array,
    z: Array,
    *,
    tol: float = 1e-12,
) -> Tuple[Array, Array, Array, Array]:
    """
    Compute Euclidean vector length and normalized direction of a spin vector.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    length = np.sqrt(x**2 + y**2 + z**2)
    valid = length > tol

    sx = np.zeros_like(x, dtype=float)
    sy = np.zeros_like(y, dtype=float)
    sz = np.zeros_like(z, dtype=float)
    sx[valid] = x[valid] / length[valid]
    sy[valid] = y[valid] / length[valid]
    sz[valid] = z[valid] / length[valid]

    return length, sx, sy, sz
