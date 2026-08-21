from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array


class SqueezingSeries(BaseModel):
    """Generalized squeezing for one group or the full physical system."""

    t: Array
    N_atoms: int
    covariance: Array
    covariance_eigvals: Array
    lambda_min: Array
    xi2: Array
    xi2_db: Array
    N_c: Array
    excited_fraction_active: Array

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_shapes(self) -> "SqueezingSeries":
        """Validate the common time axis and four-observable covariance."""
        self.t = np.asarray(self.t, dtype=float)
        if self.t.ndim != 1 or self.t.size == 0:
            raise ValueError(
                "SqueezingSeries.t must be a nonempty one-dimensional array."
            )
        if self.N_atoms < 0:
            raise ValueError("SqueezingSeries.N_atoms must be non-negative.")

        num_times = self.t.size
        self.covariance = np.asarray(self.covariance, dtype=float)
        if self.covariance.shape != (num_times, 4, 4):
            raise ValueError(
                "covariance must have shape (len(t), 4, 4)."
            )
        if not np.allclose(
            self.covariance,
            np.swapaxes(self.covariance, -1, -2),
            atol=1e-12,
            rtol=0.0,
            equal_nan=True,
        ):
            raise ValueError("covariance must be symmetric.")

        self.covariance_eigvals = np.asarray(
            self.covariance_eigvals,
            dtype=float,
        )
        if self.covariance_eigvals.shape != (num_times, 4):
            raise ValueError(
                "covariance_eigvals must have shape (len(t), 4)."
            )

        for field_name in (
            "lambda_min",
            "xi2",
            "xi2_db",
            "N_c",
            "excited_fraction_active",
        ):
            values = np.asarray(getattr(self, field_name), dtype=float)
            if values.shape != (num_times,):
                raise ValueError(f"{field_name} must have shape (len(t),).")
            setattr(self, field_name, values)
        return self


class SqueezingResult(BaseModel):
    """Group-resolved squeezing and an optional full-system result."""

    group_resolved: tuple[SqueezingSeries, ...]
    full_system: SqueezingSeries | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_results(self) -> "SqueezingResult":
        """Require nonempty group results on one shared time grid."""
        if not self.group_resolved:
            raise ValueError(
                "SqueezingResult.group_resolved must contain at least one group."
            )

        reference_t = self.group_resolved[0].t
        for group in self.group_resolved[1:]:
            if group.t.shape != reference_t.shape or not np.allclose(
                group.t,
                reference_t,
                atol=1e-12,
                rtol=0.0,
            ):
                raise ValueError(
                    "All group-resolved squeezing results must share a time grid."
                )

        if self.full_system is not None:
            if (
                self.full_system.t.shape != reference_t.shape
                or not np.allclose(
                    self.full_system.t,
                    reference_t,
                    atol=1e-12,
                    rtol=0.0,
                )
            ):
                raise ValueError(
                    "Full-system and group squeezing must share a time grid."
                )
            expected_atoms = sum(
                group.N_atoms for group in self.group_resolved
            )
            if self.full_system.N_atoms != expected_atoms:
                raise ValueError(
                    "Full-system N_atoms must equal the sum over groups."
                )
        return self


__all__ = ["SqueezingResult", "SqueezingSeries"]
