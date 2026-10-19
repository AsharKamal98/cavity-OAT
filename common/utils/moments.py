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


def _vector_sems_from_jacobian(
    jacobian: Array,
    covariance: Array,
    valid: Array,
    *,
    n_samples: int,
    tol: float,
) -> Array:
    """Propagate vector sample covariance and mask undefined derivatives."""
    from common.utils.uncertainty import propagate_covariance, std_from_covariance

    if (
        isinstance(n_samples, (bool, np.bool_))
        or not isinstance(n_samples, (int, np.integer))
        or n_samples < 2
    ):
        raise ValueError("n_samples must be an integer of at least two.")
    covariance = np.asarray(covariance)
    dimension = jacobian.shape[-1]
    if covariance.shape != valid.shape + (dimension, dimension):
        raise ValueError(
            "Vector covariance must match the batch axes and input dimension."
        )
    propagated = propagate_covariance(covariance, jacobian, tol=tol)
    sem = std_from_covariance(propagated, tol=tol) / np.sqrt(n_samples)
    sem[~valid] = np.nan
    return sem


def vector_length_sem(
    components: Array,
    covariance: Array,
    *,
    n_samples: int,
    tol: float = 1e-12,
) -> Array:
    """Return first-order norm SEMs for (..., D) vectors and sample covariance.

    Use two components for an xy-plane norm or three for the full length.
    The norm derivative is undefined at zero length, where the SEM is NaN.
    """
    components = np.asarray(components, dtype=float)
    if components.ndim < 1 or components.shape[-1] == 0:
        raise ValueError("components must have nonempty shape (..., D).")
    length = np.linalg.norm(components, axis=-1)
    valid = length > tol
    jacobian = np.zeros(length.shape + (1, components.shape[-1]), dtype=float)
    jacobian[..., 0, :][valid] = components[valid] / length[valid, np.newaxis]
    return _vector_sems_from_jacobian(
        jacobian, covariance, valid, n_samples=n_samples, tol=tol
    )[..., 0]


def spin_length_and_direction_sems(
    length: Array,
    nx: Array,
    ny: Array,
    nz: Array,
    covariance: Array,
    *,
    n_samples: int,
    tol: float = 1e-12,
) -> tuple[Array, Array, Array, Array]:
    """Return first-order SEMs of length, nx, ny, nz from Cartesian covariance.

    Inputs are the length and direction of the averaged Cartesian spin.
    Covariance is the trajectory sample covariance, not covariance of the mean.
    Undefined derivatives at length <= tol give NaN rather than zero SEMs.
    """
    length = np.asarray(length, dtype=float)
    direction = np.stack((nx, ny, nz), axis=-1)
    if direction.shape != length.shape + (3,):
        raise ValueError("Spin length and normalized components must have matching shapes.")
    valid = length > tol
    jacobian = np.zeros(length.shape + (4, 3), dtype=float)
    jacobian[..., 0, :][valid] = direction[valid]

    # d|J|/dJ = n; d(J/|J|)/dJ = (I - n n.T) / |J|.
    normalized_jacobian = jacobian[..., 1:, :]
    valid_direction = direction[valid]
    direction_outer = (
        valid_direction[..., :, np.newaxis] * valid_direction[..., np.newaxis, :]
    )
    normalized_jacobian[valid] = (
        np.eye(3) - direction_outer
    ) / length[valid, np.newaxis, np.newaxis]

    sem = _vector_sems_from_jacobian(
        jacobian, covariance, valid, n_samples=n_samples, tol=tol
    )
    return sem[..., 0], sem[..., 1], sem[..., 2], sem[..., 3]


def spin_angle_sems(
    length: Array,
    nx: Array,
    ny: Array,
    nz: Array,
    covariance: Array,
    *,
    n_samples: int,
    tol: float = 1e-12,
) -> tuple[Array, Array]:
    """Return local first-order theta/phi SEMs from Cartesian sample covariance.

    Compose angle and normalization derivatives directly, avoiding another
    stored covariance. Angles use theta=acos(-nz) and phi=atan2(ny,nx).
    At zero spin length or the polar axis, both SEMs are NaN. Phi SEM is a
    local angular width in radians, including near its +/-pi branch cut.
    """
    length, nx, ny, nz = (
        np.asarray(values, dtype=float) for values in (length, nx, ny, nz)
    )
    if any(values.shape != length.shape for values in (nx, ny, nz)):
        raise ValueError("Spin length and normalized components must have matching shapes.")
    transverse = np.hypot(nx, ny)
    valid = (length > tol) & (transverse > tol)
    jacobian = np.zeros(length.shape + (2, 3), dtype=float)
    theta_gradient = jacobian[..., 0, :]
    phi_gradient = jacobian[..., 1, :]

    # These are the combined derivatives with respect to the raw J vector.
    theta_denominator = length[valid] * transverse[valid]
    phi_denominator = length[valid] * transverse[valid] ** 2
    theta_gradient[valid, 0] = -nx[valid] * nz[valid] / theta_denominator
    theta_gradient[valid, 1] = -ny[valid] * nz[valid] / theta_denominator
    theta_gradient[valid, 2] = transverse[valid] / length[valid]
    phi_gradient[valid, 0] = -ny[valid] / phi_denominator
    phi_gradient[valid, 1] = nx[valid] / phi_denominator

    sem = _vector_sems_from_jacobian(
        jacobian, covariance, valid, n_samples=n_samples, tol=tol
    )
    return sem[..., 0], sem[..., 1]
