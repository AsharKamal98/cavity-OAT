from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array
from parser.uncertainty import SamplingUncertainty


class BlochVectorSeries(BaseModel):
    """Cartesian components and derived lengths of one vector series."""

    x: Array
    y: Array
    z: Array
    length: Array
    xy_length: Array
    uncertainty: SamplingUncertainty | None = None
    x_sem: Array | None = None
    y_sem: Array | None = None
    z_sem: Array | None = None
    length_sem: Array | None = None
    xy_length_sem: Array | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="before")
    @classmethod
    def derive_lengths(cls, values: Any) -> Any:
        """Validate Cartesian components and derive their lengths."""
        if not isinstance(values, Mapping) or not all(
            component in values for component in ("x", "y", "z")
        ):
            return values

        data = dict(values)
        x = np.asarray(data["x"], dtype=float)
        y = np.asarray(data["y"], dtype=float)
        z = np.asarray(data["z"], dtype=float)
        if x.ndim != 1 or x.shape != y.shape or x.shape != z.shape:
            raise ValueError(
                "BlochVectorSeries requires matching one-dimensional x/y/z arrays."
            )

        data.update(
            x=x,
            y=y,
            z=z,
            length=np.sqrt(x**2 + y**2 + z**2),
            xy_length=np.sqrt(x**2 + y**2),
        )
        return data

    @model_validator(mode="after")
    def derive_sems(self) -> "BlochVectorSeries":
        """Derive component and norm SEMs from optional sample covariance."""
        from common.utils.moments import vector_length_sem
        from common.utils.uncertainty import std_from_covariance

        self.uncertainty = getattr(self, "uncertainty", None)
        self.x_sem = self.y_sem = self.z_sem = None
        self.length_sem = self.xy_length_sem = None
        if self.uncertainty is None:
            return self

        covariance = self.uncertainty.sample_covariance
        if covariance.shape != (len(self.x), 3, 3):
            raise ValueError("Vector sample covariance must have shape (len(x), 3, 3).")
        sem = std_from_covariance(covariance) / np.sqrt(self.uncertainty.n_samples)
        self.x_sem = sem[:, 0].copy()
        self.y_sem = sem[:, 1].copy()
        self.z_sem = sem[:, 2].copy()
        components = np.stack((self.x, self.y, self.z), axis=-1)
        self.length_sem = vector_length_sem(
            components, covariance, n_samples=self.uncertainty.n_samples
        )
        self.xy_length_sem = vector_length_sem(
            components[:, :2], covariance[:, :2, :2],
            n_samples=self.uncertainty.n_samples,
        )
        return self


__all__ = ["BlochVectorSeries"]
