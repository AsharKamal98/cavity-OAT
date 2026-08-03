from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array


class SqueezingMomentSeries(BaseModel):
    """Group-local moments required to compute generalized squeezing."""

    t: Array
    N_group: int
    mean_operators: Array
    second_moments: Array
    N_c: Array
    excited_fraction_active: Array
    theta_J: Array
    phi_J: Array
    theta_S: Array
    phi_S: Array

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_shapes(self) -> "SqueezingMomentSeries":
        """Validate the shared time axis and fluctuation-moment dimensions."""
        self.t = np.asarray(self.t, dtype=float)
        if self.t.ndim != 1:
            raise ValueError("SqueezingMomentSeries.t must be one-dimensional.")
        if self.N_group < 0:
            raise ValueError("SqueezingMomentSeries.N_group must be non-negative.")

        num_times = self.t.size
        self.mean_operators = np.asarray(self.mean_operators, dtype=float)
        self.second_moments = np.asarray(self.second_moments, dtype=float)
        if self.mean_operators.shape != (num_times, 4):
            raise ValueError(
                "mean_operators must have shape (len(t), 4)."
            )
        if self.second_moments.shape != (num_times, 4, 4):
            raise ValueError(
                "second_moments must have shape (len(t), 4, 4)."
            )

        for field_name in (
            "N_c",
            "excited_fraction_active",
            "theta_J",
            "phi_J",
            "theta_S",
            "phi_S",
        ):
            values = np.asarray(getattr(self, field_name), dtype=float)
            if values.shape != (num_times,):
                raise ValueError(f"{field_name} must have shape (len(t),).")
            setattr(self, field_name, values)

        return self


class SqueezingSeries(BaseModel):
    """Generalized squeezing calculated for one physical group."""

    t: Array
    N_group: int
    xi2: Array
    xi2_db: Array
    lambda_min: Array
    covariance_eigvals: Array
    N_c: Array
    excited_fraction_active: Array
    theta_J: Array
    phi_J: Array
    theta_S: Array
    phi_S: Array

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_shapes(self) -> "SqueezingSeries":
        """Validate the shared time axis and squeezing-result dimensions."""
        self.t = np.asarray(self.t, dtype=float)
        if self.t.ndim != 1:
            raise ValueError("SqueezingSeries.t must be one-dimensional.")
        if self.N_group < 0:
            raise ValueError("SqueezingSeries.N_group must be non-negative.")

        num_times = self.t.size
        self.covariance_eigvals = np.asarray(
            self.covariance_eigvals,
            dtype=float,
        )
        if self.covariance_eigvals.shape != (num_times, 4):
            raise ValueError(
                "covariance_eigvals must have shape (len(t), 4)."
            )

        for field_name in (
            "xi2",
            "xi2_db",
            "lambda_min",
            "N_c",
            "excited_fraction_active",
            "theta_J",
            "phi_J",
            "theta_S",
            "phi_S",
        ):
            values = np.asarray(getattr(self, field_name), dtype=float)
            if values.shape != (num_times,):
                raise ValueError(f"{field_name} must have shape (len(t),).")
            setattr(self, field_name, values)

        return self


__all__ = ["SqueezingMomentSeries", "SqueezingSeries"]
