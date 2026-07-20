from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array


class BlochVectorSeries(BaseModel):
    """Cartesian components and derived lengths of one vector series."""

    x: Array
    y: Array
    z: Array
    length: Array
    xy_length: Array

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


__all__ = ["BlochVectorSeries"]
