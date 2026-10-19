from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.uncertainty import SamplingUncertainty


def _validated_covariance(covariance: Array, *, tol: float) -> Array:
    """Validate real covariance matrices on the final two array axes."""
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive.")
    values = np.asarray(covariance)
    if np.iscomplexobj(values):
        raise ValueError("covariance must be real.")
    values = np.asarray(values, dtype=float)
    if (
        values.ndim < 2
        or values.shape[-2] != values.shape[-1]
        or any(size == 0 for size in values.shape)
    ):
        raise ValueError("covariance must have nonempty shape (..., D, D).")
    if not np.all(np.isfinite(values)):
        raise ValueError("covariance must contain only finite values.")
    if not np.allclose(
        values, np.swapaxes(values, -1, -2), atol=tol, rtol=0.0
    ):
        raise ValueError("covariance must be symmetric within tol.")
    if np.any(np.diagonal(values, axis1=-2, axis2=-1) < -tol):
        raise ValueError("covariance must have nonnegative variances within tol.")
    return values


def compute_sample_covariance(
    samples: Array,
    *,
    mean: Array | None = None,
) -> SamplingUncertainty | None:
    """Estimate covariance across independent, real vector samples.

    Samples have shape ``(R, ..., D)``: sample axis first, observable axis last.
    A one-dimensional sequence is treated as ``R`` scalar samples (``D=1``).
    An optional precomputed sample mean must match ``samples.shape[1:]``.
    Return ``None`` for one sample, since its sampling variance is unknown.
    The ``R - 1`` denominator gives the unbiased sample covariance. Conversion
    to covariance of the mean or SEM belongs to downstream analysis.
    """
    values = np.asarray(samples)
    if np.iscomplexobj(values):
        raise ValueError("samples must be real.")
    values = np.asarray(values, dtype=float)
    if values.ndim == 0 or any(size == 0 for size in values.shape):
        raise ValueError("samples must have nonempty shape (R, ..., D).")
    scalar_samples = values.ndim == 1
    if scalar_samples:
        values = values[:, np.newaxis]
    if not np.all(np.isfinite(values)):
        raise ValueError("samples must contain only finite values.")

    if mean is None:
        sample_mean = np.mean(values, axis=0)
    else:
        sample_mean = np.asarray(mean)
        if np.iscomplexobj(sample_mean):
            raise ValueError("mean must be real.")
        sample_mean = np.asarray(sample_mean, dtype=float)
        if scalar_samples and sample_mean.ndim == 0:
            sample_mean = sample_mean[np.newaxis]
        if sample_mean.shape != values.shape[1:]:
            raise ValueError("mean must match the sample batch and observable axes.")
        if not np.all(np.isfinite(sample_mean)):
            raise ValueError("mean must contain only finite values.")

    n_samples = values.shape[0]
    if n_samples < 2:
        return None
    centered = values - sample_mean
    covariance = np.einsum(
        "r...i,r...j->...ij", centered, centered, optimize=True
    ) / (n_samples - 1)
    return SamplingUncertainty(
        sample_covariance=covariance,
        n_samples=n_samples,
    )


def propagate_covariance(
    covariance: Array,
    transformation: Array,
    *,
    tol: float = 1e-12,
) -> Array:
    """Return ``A @ C @ A.T`` for a real linear map or local Jacobian.

    ``covariance`` has shape ``(..., D, D)`` and ``transformation`` has shape
    ``(..., M, D)``. Leading batch axes follow NumPy broadcasting rules; the
    output has shape ``(..., M, M)``. This propagation is exact for a fixed
    affine map; a nonlinear observable requires its Jacobian and gives a
    first-order approximation. No sample-count normalization is applied.
    """
    values = _validated_covariance(covariance, tol=tol)
    matrix = np.asarray(transformation)
    if np.iscomplexobj(matrix):
        raise ValueError("transformation must be real.")
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim < 2 or any(size == 0 for size in matrix.shape):
        raise ValueError("transformation must have nonempty shape (..., M, D).")
    if matrix.shape[-1] != values.shape[-1]:
        raise ValueError("transformation must have one column per input observable.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("transformation must contain only finite values.")
    try:
        np.broadcast_shapes(values.shape[:-2], matrix.shape[:-2])
    except ValueError as error:
        raise ValueError(
            "Covariance and transformation batch axes must be broadcast-compatible."
        ) from error

    propagated = matrix @ values @ np.swapaxes(matrix, -1, -2)
    propagated = 0.5 * (propagated + np.swapaxes(propagated, -1, -2))
    return _validated_covariance(propagated, tol=tol)


def std_from_covariance(covariance: Array, *, tol: float = 1e-12) -> Array:
    """Return square roots of diagonal variances, with shape ``(..., D)``.

    Supply sample covariance for trajectory standard deviations, or covariance
    of the mean for SEMs. This helper never divides by a sample count. Negative
    diagonal roundoff no larger than ``tol`` is clipped to zero; larger
    negative variances raise ``ValueError``.
    """
    values = _validated_covariance(covariance, tol=tol)
    variances = np.diagonal(values, axis1=-2, axis2=-1)
    return np.sqrt(np.maximum(variances, 0.0))


__all__ = [
    "compute_sample_covariance",
    "propagate_covariance",
    "std_from_covariance",
]
