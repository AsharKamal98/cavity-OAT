from __future__ import annotations

from collections.abc import Sequence
from typing import Tuple

import numpy as np

from parser.common import Array


def as_series_tuple(series: Sequence[Array] | Array) -> tuple[Array, ...]:
    """Return one time series or a sequence of time series as a tuple."""
    values = np.asarray(series, dtype=float)
    if values.ndim == 1:
        return (values,)
    if values.ndim == 2:
        return tuple(values)
    raise ValueError("Expected one time series or a sequence of time series.")


def norm_spin_components_from_spin_components(
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

    nx = np.zeros_like(x, dtype=float)
    ny = np.zeros_like(y, dtype=float)
    nz = np.zeros_like(z, dtype=float)
    nx[valid] = x[valid] / length[valid]
    ny[valid] = y[valid] / length[valid]
    nz[valid] = z[valid] / length[valid]

    return length, nx, ny, nz


def angles_from_norm_spin_components(
    nx: Array,
    ny: Array,
    nz: Array,
    valid: Array,
    tol: float = 1e-12,
) -> Tuple[Array, Array]:
    """
    Compute J-sphere angles from normalized spin components.
    """
    nx = np.asarray(nx, dtype=float)
    ny = np.asarray(ny, dtype=float)
    nz = np.asarray(nz, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    theta = np.zeros_like(nz, dtype=float)
    theta[valid] = np.arccos(np.clip(-nz[valid], -1.0, 1.0))

    phi = np.arctan2(ny, nx)
    r_perp = np.sqrt(nx**2 + ny**2)
    phi[r_perp < tol] = 0.0

    return theta, phi


def norm_spin_components_from_angles(
    theta: Array,
    phi: Array,
) -> Tuple[Array, Array, Array]:
    """
    Compute normalized spin components from J-sphere angles.
    """
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)

    sin_theta = np.sin(theta)
    nx = sin_theta * np.cos(phi)
    ny = sin_theta * np.sin(phi)
    nz = -np.cos(theta)

    return nx, ny, nz


def spin_components_from_norm_spin_components(
    length: Array,
    nx: Array,
    ny: Array,
    nz: Array,
) -> Tuple[Array, Array, Array]:
    """
    Compute spin components from vector length and normalized direction.
    """
    length = np.asarray(length, dtype=float)
    nx = np.asarray(nx, dtype=float)
    ny = np.asarray(ny, dtype=float)
    nz = np.asarray(nz, dtype=float)

    x = length * nx
    y = length * ny
    z = length * nz

    return x, y, z
