from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from parser.common import Array


class SamplingUncertainty(BaseModel):
    """Sample covariance and sample count for an estimated vector mean.

    The final two covariance axes index observables; preceding axes may index
    saved times or other batches. This is sampling scatter across independent
    samples, rather than a quantum operator covariance.
    """

    sample_covariance: Array
    n_samples: int = Field(ge=2, strict=True)

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("sample_covariance", mode="before")
    @classmethod
    def convert_covariance(cls, value: object) -> Array:
        values = np.asarray(value)
        if np.iscomplexobj(values):
            raise ValueError("sample_covariance must be real.")
        return np.asarray(values, dtype=float)

    @model_validator(mode="after")
    def validate_covariance(self) -> "SamplingUncertainty":
        covariance = self.sample_covariance
        if (
            covariance.ndim < 2
            or covariance.shape[-2] != covariance.shape[-1]
            or any(size == 0 for size in covariance.shape)
        ):
            raise ValueError(
                "sample_covariance must have nonempty shape (..., D, D)."
            )
        if not np.all(np.isfinite(covariance)):
            raise ValueError("sample_covariance must contain only finite values.")
        if not np.allclose(
            covariance, np.swapaxes(covariance, -1, -2), atol=1e-12, rtol=0.0
        ):
            raise ValueError("sample_covariance must be symmetric.")
        if np.any(np.diagonal(covariance, axis1=-2, axis2=-1) < -1e-12):
            raise ValueError("sample_covariance must have nonnegative variances.")
        return self


__all__ = ["SamplingUncertainty"]
