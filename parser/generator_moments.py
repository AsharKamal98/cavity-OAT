from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, field_validator, model_validator

from parser.common import Array


class TrajectoryResolvedMomentSeries(BaseModel):
    """Trajectory-resolved SU(3) moments in one caller-selected frame.

    ``second_moments`` uses the row-major ordering returned by
    ``np.triu_indices(8 * num_groups)``. The container does not encode whether
    its generator frame is the lab frame or a derived dressed frame, so callers
    must preserve that provenance in the surrounding analysis flow.
    """

    t: Array
    first_moments: Array
    second_moments: Array

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @field_validator("t", "first_moments", "second_moments", mode="before")
    @classmethod
    def convert_arrays(cls, value: object) -> Array:
        values = np.asarray(value)
        if np.iscomplexobj(values):
            if not np.allclose(values.imag, 0.0, atol=1e-12, rtol=0.0):
                raise ValueError("Generator-moment arrays must be real.")
            values = values.real
        return np.asarray(values, dtype=float)

    @model_validator(mode="after")
    def validate_shapes(self) -> "TrajectoryResolvedMomentSeries":
        """Validate the trajectory, time, group, and packed-moment axes."""
        if self.t.ndim != 1:
            raise ValueError(
                "TrajectoryResolvedMomentSeries.t must be one-dimensional."
            )
        if self.t.size == 0:
            raise ValueError(
                "TrajectoryResolvedMomentSeries.t must contain at least one point."
            )
        if not np.all(np.isfinite(self.t)):
            raise ValueError("TrajectoryResolvedMomentSeries.t must be finite.")
        if np.any(np.diff(self.t) <= 0.0):
            raise ValueError(
                "TrajectoryResolvedMomentSeries.t must be strictly increasing."
            )
        if self.first_moments.ndim != 4 or self.first_moments.shape[3] != 8:
            raise ValueError(
                "first_moments must have shape "
                "(num_trajectories, len(t), num_groups, 8)."
            )
        if self.first_moments.shape[0] == 0:
            raise ValueError("first_moments must contain at least one trajectory.")
        if self.first_moments.shape[1] != self.t.size:
            raise ValueError("first_moments must use the same time axis as t.")
        if self.first_moments.shape[2] == 0:
            raise ValueError("first_moments must contain at least one group.")
        if not np.all(np.isfinite(self.first_moments)):
            raise ValueError("first_moments must contain only finite values.")

        trajectory_count = self.first_moments.shape[0]
        generator_count = 8 * self.first_moments.shape[2]
        packed_count = generator_count * (generator_count + 1) // 2
        if self.second_moments.shape != (
            trajectory_count,
            self.t.size,
            packed_count,
        ):
            raise ValueError(
                "second_moments must have shape "
                f"(num_trajectories, len(t), {packed_count}) for the packed "
                f"upper triangle of the {generator_count}x{generator_count} "
                "symmetric moment matrix."
            )
        if not np.all(np.isfinite(self.second_moments)):
            raise ValueError("second_moments must contain only finite values.")
        return self

    @property
    def num_trajectories(self) -> int:
        """Number of individually retained trajectories."""
        return int(self.first_moments.shape[0])


__all__ = ["TrajectoryResolvedMomentSeries"]
