from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, model_validator

from parser.common import Array


class SMomentSnapshot(BaseModel):
    """First-order S-sphere moments for one saved MCWF trajectory snapshot."""

    t: float
    integration_phase_index: int
    x: tuple[float, ...]
    y: tuple[float, ...]
    z: tuple[float, ...]


class SBasisSeries(BaseModel):
    """Group-resolved effective S-qubit basis on a shared time grid."""

    zero_state: tuple[Array, ...]
    one_state: tuple[Array, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_states(self) -> "SBasisSeries":
        """Validate matching normalized three-level basis-state series."""
        if not self.zero_state or len(self.zero_state) != len(self.one_state):
            raise ValueError(
                "S-basis zero_state and one_state must contain the same "
                "nonzero number of groups."
            )

        zero_states = tuple(
            np.asarray(states, dtype=np.complex128)
            for states in self.zero_state
        )
        one_states = tuple(
            np.asarray(states, dtype=np.complex128)
            for states in self.one_state
        )
        reference_shape = zero_states[0].shape
        if len(reference_shape) != 2 or reference_shape[1] != 3:
            raise ValueError(
                "Each S-basis state series must have shape (num_times, 3)."
            )
        if any(
            zero.shape != reference_shape or one.shape != reference_shape
            for zero, one in zip(zero_states, one_states)
        ):
            raise ValueError(
                "All S-basis state series must share shape (num_times, 3)."
            )

        zero_norms = np.linalg.norm(np.asarray(zero_states), axis=-1)
        one_norms = np.linalg.norm(np.asarray(one_states), axis=-1)
        overlaps = np.sum(
            np.asarray(zero_states).conjugate() * np.asarray(one_states),
            axis=-1,
        )
        if not np.allclose(zero_norms, 1.0, atol=1e-12, rtol=0.0):
            raise ValueError("Every S-basis zero_state must be normalized.")
        if not np.allclose(one_norms, 1.0, atol=1e-12, rtol=0.0):
            raise ValueError("Every S-basis one_state must be normalized.")
        if not np.allclose(overlaps, 0.0, atol=1e-12, rtol=0.0):
            raise ValueError(
                "S-basis zero_state and one_state must be orthogonal."
            )

        self.zero_state = zero_states
        self.one_state = one_states
        return self


class SMomentSeries(BaseModel):
    """Group-resolved first-order S-sphere moment time series."""

    t: Array
    integration_phase_index: Array | None = None
    x: tuple[Array, ...] | None = None
    y: tuple[Array, ...] | None = None
    z: tuple[Array, ...] | None = None
    length: tuple[Array, ...] | None = None
    nx: tuple[Array, ...] | None = None
    ny: tuple[Array, ...] | None = None
    nz: tuple[Array, ...] | None = None
    theta: tuple[Array, ...] | None = None
    phi: tuple[Array, ...] | None = None
    basis: SBasisSeries | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_basis(self) -> "SMomentSeries":
        """Require an attached S basis to match the series time and groups."""
        if self.basis is None:
            return self

        num_times = np.asarray(self.t).size
        if self.basis.zero_state[0].shape[0] != num_times:
            raise ValueError("S-moment basis state series must match t.")

        group_fields = (
            self.x,
            self.y,
            self.z,
            self.length,
            self.nx,
            self.ny,
            self.nz,
            self.theta,
            self.phi,
        )
        populated_group_fields = [
            field for field in group_fields if field is not None
        ]
        if populated_group_fields and any(
            len(field) != len(self.basis.one_state)
            for field in populated_group_fields
        ):
            raise ValueError(
                "S-moment basis and group-resolved fields must have matching "
                "groups."
            )
        return self

    @classmethod
    def attach_norm_spin_components_from_spin_components(
        cls,
        s_moments: "SMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """Attach S-vector lengths and normalized group directions."""
        from common.utils.moments import norm_spin_components_from_spin_components

        if (
            s_moments.x is None
            or s_moments.y is None
            or s_moments.z is None
        ):
            raise ValueError(
                "Spin-direction attachment requires group spin components."
            )

        group_results = [
            norm_spin_components_from_spin_components(x_g, y_g, z_g, tol=tol)
            for x_g, y_g, z_g in zip(
                s_moments.x,
                s_moments.y,
                s_moments.z,
            )
        ]
        s_moments.length = tuple(result[0] for result in group_results)
        s_moments.nx = tuple(result[1] for result in group_results)
        s_moments.ny = tuple(result[2] for result in group_results)
        s_moments.nz = tuple(result[3] for result in group_results)

    @classmethod
    def attach_angles_from_norm_spin_components(
        cls,
        s_moments: "SMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """Attach S-sphere angles using the positive-z pole convention."""
        from common.utils.moments import angles_from_norm_spin_components

        if (
            s_moments.length is None
            or s_moments.nx is None
            or s_moments.ny is None
            or s_moments.nz is None
        ):
            raise ValueError(
                "Spin angle attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            angles_from_norm_spin_components(
                nx_g,
                ny_g,
                -nz_g,
                valid=length_g > tol,
                tol=tol,
            )
            for length_g, nx_g, ny_g, nz_g in zip(
                s_moments.length,
                s_moments.nx,
                s_moments.ny,
                s_moments.nz,
            )
        ]
        s_moments.theta = tuple(result[0] for result in group_results)
        s_moments.phi = tuple(result[1] for result in group_results)


__all__ = ["SBasisSeries", "SMomentSeries", "SMomentSnapshot"]
