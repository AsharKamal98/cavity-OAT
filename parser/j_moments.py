from __future__ import annotations

import numpy as np
from pydantic import BaseModel, ConfigDict, Field, model_validator

from parser.common import Array
from parser.uncertainty import SamplingUncertainty


_AFFINE_FIELD_NAMES = ("x", "y", "z", "N_e", "N_j")
_DERIVED_FIELD_NAMES = ("length", "nx", "ny", "nz", "theta", "phi")


def _validated_sem(
    values: Array,
    reference: Array | None,
    *,
    field_name: str,
    allow_undefined: bool = False,
) -> Array:
    """Validate SEMs, allowing NaN for undefined nonlinear uncertainties."""
    array = np.asarray(values)
    if np.iscomplexobj(array):
        raise ValueError(f"{field_name} must be real.")
    array = np.asarray(array, dtype=float)
    if reference is None or array.ndim != 1 or array.shape != np.shape(reference):
        raise ValueError(f"{field_name} must match its central moment time series.")
    invalid = np.isinf(array) if allow_undefined else ~np.isfinite(array)
    if np.any(invalid) or np.any(array < 0.0):
        suffix = "; NaN may mark undefined points" if allow_undefined else ""
        raise ValueError(f"{field_name} must be finite and nonnegative{suffix}.")
    return array


def _group_cartesian_covariances(
    uncertainty: SamplingUncertainty, *, num_times: int, num_groups: int
) -> tuple[Array, ...]:
    """Extract x/y/z sample-covariance blocks in physical group order."""
    num_fields = len(_AFFINE_FIELD_NAMES)
    size = num_fields * num_groups
    covariance = uncertainty.sample_covariance
    if covariance.shape != (num_times, size, size):
        raise ValueError("J sample covariance must have shape (len(t), 5 * G, 5 * G).")
    return tuple(
        covariance[
            :, num_fields * group : num_fields * group + 3,
            num_fields * group : num_fields * group + 3,
        ]
        for group in range(num_groups)
    )


class AggregatedJMomentSnapshot(BaseModel):
    """Full-system J moments for one saved trajectory snapshot."""

    x: float
    y: float
    z: float
    N_e: float
    N_j: float


class JMomentSnapshot(BaseModel):
    """First-order J-sphere moments for one saved trajectory snapshot."""

    # Legacy note: these fields were previously named Jx/Jy/Jz and
    # Jx_groups/Jy_groups/Jz_groups.
    t: float
    integration_phase_index: int
    x: tuple[float, ...]
    y: tuple[float, ...]
    z: tuple[float, ...]
    N_e: tuple[float, ...]
    N_j: tuple[float, ...]
    aggregated: AggregatedJMomentSnapshot


class AggregatedJMomentSeries(BaseModel):
    """Full-system J-moment time series."""

    x: Array | None = None
    y: Array | None = None
    z: Array | None = None
    length: Array | None = None
    nx: Array | None = None
    ny: Array | None = None
    nz: Array | None = None
    N_e: Array | None = None
    N_j: Array | None = None
    theta: Array | None = None
    phi: Array | None = None
    x_sem: Array | None = None
    y_sem: Array | None = None
    z_sem: Array | None = None
    N_e_sem: Array | None = None
    N_j_sem: Array | None = None
    length_sem: Array | None = None
    nx_sem: Array | None = None
    ny_sem: Array | None = None
    nz_sem: Array | None = None
    theta_sem: Array | None = None
    phi_sem: Array | None = None
    uncertainty: SamplingUncertainty | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @model_validator(mode="after")
    def validate_uncertainty(self) -> "AggregatedJMomentSeries":
        for name in (*_AFFINE_FIELD_NAMES, *_DERIVED_FIELD_NAMES):
            field_name = f"{name}_sem"
            values = getattr(self, field_name, None)
            setattr(
                self, field_name,
                None if values is None else _validated_sem(
                    values, getattr(self, name), field_name=field_name,
                    allow_undefined=name in _DERIVED_FIELD_NAMES,
                ),
            )
        self.uncertainty = getattr(self, "uncertainty", None)
        if self.uncertainty is not None:
            if self.x is None or self.uncertainty.sample_covariance.shape != (
                len(self.x), 5, 5
            ):
                raise ValueError(
                    "Aggregated J sample covariance must have shape (len(x), 5, 5)."
                )
        return self


class JMomentSeries(BaseModel):
    """Per-timestep first-order J-sphere moments for one (single or averaged) trajectory."""
    # Legacy note: these fields were previously named Jx/Jy/Jz, Jx_groups/
    # Jy_groups/Jz_groups, J_len, and sx/sy/sz.

    t: Array
    integration_phase_index: Array | None = None
    # spin components
    x: tuple[Array, ...] | None = None
    y: tuple[Array, ...] | None = None
    z: tuple[Array, ...] | None = None
    length: tuple[Array, ...] | None = None
    # normalized spin components / directions
    nx: tuple[Array, ...] | None = None
    ny: tuple[Array, ...] | None = None
    nz: tuple[Array, ...] | None = None
    # atom numbers
    N_e: tuple[Array, ...] | None = None
    N_j: tuple[Array, ...] | None = None
    # angles
    theta: tuple[Array, ...] | None = None
    phi: tuple[Array, ...] | None = None
    # Sampling SEMs; absent for deterministic outputs. Nonlinear SEMs may
    # contain NaN at singular times where first-order propagation is undefined.
    x_sem: tuple[Array, ...] | None = None
    y_sem: tuple[Array, ...] | None = None
    z_sem: tuple[Array, ...] | None = None
    N_e_sem: tuple[Array, ...] | None = None
    N_j_sem: tuple[Array, ...] | None = None
    length_sem: tuple[Array, ...] | None = None
    nx_sem: tuple[Array, ...] | None = None
    ny_sem: tuple[Array, ...] | None = None
    nz_sem: tuple[Array, ...] | None = None
    theta_sem: tuple[Array, ...] | None = None
    phi_sem: tuple[Array, ...] | None = None
    # Full group-major covariance, ordered (x, y, z, N_e, N_j) per group.
    uncertainty: SamplingUncertainty | None = None
    # full-system fields
    aggregated: AggregatedJMomentSeries = Field(
        default_factory=AggregatedJMomentSeries
    )

    @model_validator(mode="after")
    def validate_uncertainty(self) -> "JMomentSeries":
        self.aggregated = AggregatedJMomentSeries.model_validate(self.aggregated)
        for name in (*_AFFINE_FIELD_NAMES, *_DERIVED_FIELD_NAMES):
            field_name = f"{name}_sem"
            groups = getattr(self, field_name, None)
            if groups is None:
                setattr(self, field_name, None)
                continue
            reference = getattr(self, name)
            if reference is None or not groups or len(groups) != len(reference):
                raise ValueError(f"{field_name} must match its central moment groups.")
            validated = tuple(
                _validated_sem(
                    values, central, field_name=field_name,
                    allow_undefined=name in _DERIVED_FIELD_NAMES,
                )
                for values, central in zip(groups, reference)
            )
            if any(values.shape != np.shape(self.t) for values in validated):
                raise ValueError(f"{field_name} must use the same time grid as t.")
            setattr(self, field_name, validated)
        self.uncertainty = getattr(self, "uncertainty", None)
        if self.uncertainty is not None:
            count = 0 if self.x is None else 5 * len(self.x)
            if count == 0 or self.uncertainty.sample_covariance.shape != (
                len(self.t), count, count
            ):
                raise ValueError(
                    "J sample covariance must have shape (len(t), 5 * G, 5 * G)."
                )
            aggregate_uncertainty = getattr(self.aggregated, "uncertainty", None)
            if (
                aggregate_uncertainty is not None
                and aggregate_uncertainty.n_samples != self.uncertainty.n_samples
            ):
                raise ValueError("Group and aggregated J sample counts must agree.")
        return self

    @classmethod
    def attach_aggregated_from_group_fields(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach aggregated values, covariance, and available sampling SEMs.
        """
        from common.utils.moments import (
            angles_from_norm_spin_components,
            norm_spin_components_from_spin_components,
            spin_angle_sems,
            spin_length_and_direction_sems,
        )
        from common.utils.uncertainty import (
            propagate_covariance,
            std_from_covariance,
        )

        field_names = _AFFINE_FIELD_NAMES
        group_fields = {
            field_name: getattr(j_moments, field_name)
            for field_name in field_names
        }
        missing_fields = [
            field_name
            for field_name, groups in group_fields.items()
            if groups is None
        ]
        if missing_fields:
            raise ValueError(
                "Additive aggregated-field attachment requires "
                + ", ".join(missing_fields)
                + "."
            )

        group_count = len(group_fields["x"])
        if group_count == 0 or any(
            len(groups) != group_count for groups in group_fields.values()
        ):
            raise ValueError(
                "All additive group fields must have the same nonzero group count."
            )

        aggregated_fields = {}
        for field_name, groups in group_fields.items():
            group_arrays = [np.asarray(group, dtype=float) for group in groups]
            aggregated_fields[field_name] = np.sum(
                np.stack(group_arrays),
                axis=0,
            )

        length, nx, ny, nz = norm_spin_components_from_spin_components(
            aggregated_fields["x"],
            aggregated_fields["y"],
            aggregated_fields["z"],
            tol=tol,
        )
        theta, phi = angles_from_norm_spin_components(
            nx,
            ny,
            nz,
            valid=length > tol,
            tol=tol,
        )
        aggregated = AggregatedJMomentSeries(
            **aggregated_fields,
            length=length,
            nx=nx,
            ny=ny,
            nz=nz,
            theta=theta,
            phi=phi,
        )

        source_uncertainty = getattr(j_moments, "uncertainty", None)
        if source_uncertainty is not None:
            num_fields = len(field_names)
            covariance_size = num_fields * group_count
            if source_uncertainty.sample_covariance.shape != (
                len(j_moments.t), covariance_size, covariance_size
            ):
                raise ValueError(
                    "J sample covariance must have shape (len(t), 5 * G, 5 * G)."
                )

            # Sum matching fields, including all cross-group correlations.
            sum_over_groups = np.tile(np.eye(num_fields), (1, group_count))
            aggregate_covariance = propagate_covariance(
                source_uncertainty.sample_covariance, sum_over_groups, tol=tol
            )
            aggregate_std = std_from_covariance(aggregate_covariance, tol=tol)
            aggregate_sem = aggregate_std / np.sqrt(source_uncertainty.n_samples)

            aggregated.x_sem = aggregate_sem[:, 0].copy()
            aggregated.y_sem = aggregate_sem[:, 1].copy()
            aggregated.z_sem = aggregate_sem[:, 2].copy()
            aggregated.N_e_sem = aggregate_sem[:, 3].copy()
            aggregated.N_j_sem = aggregate_sem[:, 4].copy()
            aggregated.uncertainty = SamplingUncertainty(
                sample_covariance=aggregate_covariance,
                n_samples=source_uncertainty.n_samples,
            )
            cartesian_covariance = aggregate_covariance[:, :3, :3]
            (
                aggregated.length_sem,
                aggregated.nx_sem,
                aggregated.ny_sem,
                aggregated.nz_sem,
            ) = spin_length_and_direction_sems(
                length, nx, ny, nz, cartesian_covariance,
                n_samples=source_uncertainty.n_samples, tol=tol,
            )
            aggregated.theta_sem, aggregated.phi_sem = spin_angle_sems(
                length, nx, ny, nz, cartesian_covariance,
                n_samples=source_uncertainty.n_samples, tol=tol,
            )

        j_moments.aggregated = aggregated

    @classmethod
    def attach_norm_spin_components_from_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach Euclidean length, normalized directions, and their SEMs.
        """
        from common.utils.moments import (
            norm_spin_components_from_spin_components,
            spin_length_and_direction_sems,
        )

        if (
            j_moments.x is None
            or j_moments.y is None
            or j_moments.z is None
        ):
            raise ValueError(
                "Spin-direction attachment requires group spin components."
            )

        group_results = [
            norm_spin_components_from_spin_components(x_g, y_g, z_g, tol=tol)
            for x_g, y_g, z_g in zip(
                j_moments.x,
                j_moments.y,
                j_moments.z,
            )
        ]
        j_moments.length = tuple(result[0] for result in group_results)
        j_moments.nx = tuple(result[1] for result in group_results)
        j_moments.ny = tuple(result[2] for result in group_results)
        j_moments.nz = tuple(result[3] for result in group_results)

        j_moments.length_sem = None
        j_moments.nx_sem = None
        j_moments.ny_sem = None
        j_moments.nz_sem = None
        uncertainty = getattr(j_moments, "uncertainty", None)
        if uncertainty is not None:
            covariances = _group_cartesian_covariances(
                uncertainty, num_times=len(j_moments.t), num_groups=len(group_results)
            )
            group_sems = [
                spin_length_and_direction_sems(
                    *values, covariance=covariance,
                    n_samples=uncertainty.n_samples, tol=tol,
                )
                for values, covariance in zip(group_results, covariances)
            ]
            j_moments.length_sem = tuple(result[0] for result in group_sems)
            j_moments.nx_sem = tuple(result[1] for result in group_sems)
            j_moments.ny_sem = tuple(result[2] for result in group_sems)
            j_moments.nz_sem = tuple(result[3] for result in group_sems)

    @classmethod
    def attach_angles_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
        *,
        tol: float = 1e-12,
    ) -> None:
        """
        Attach theta, phi, and their SEMs from normalized spin directions.
        """
        from common.utils.moments import angles_from_norm_spin_components, spin_angle_sems

        if (
            j_moments.length is None
            or j_moments.nx is None
            or j_moments.ny is None
            or j_moments.nz is None
        ):
            raise ValueError(
                "Spin angle attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            angles_from_norm_spin_components(
                nx_g,
                ny_g,
                nz_g,
                valid=length_g > tol,
                tol=tol,
            )
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            )
        ]
        j_moments.theta = tuple(result[0] for result in group_results)
        j_moments.phi = tuple(result[1] for result in group_results)

        j_moments.theta_sem = None
        j_moments.phi_sem = None
        uncertainty = getattr(j_moments, "uncertainty", None)
        if uncertainty is not None:
            covariances = _group_cartesian_covariances(
                uncertainty, num_times=len(j_moments.t), num_groups=len(group_results)
            )
            group_sems = [
                spin_angle_sems(
                    length_g, nx_g, ny_g, nz_g, covariance,
                    n_samples=uncertainty.n_samples, tol=tol,
                )
                for length_g, nx_g, ny_g, nz_g, covariance in zip(
                    j_moments.length, j_moments.nx, j_moments.ny,
                    j_moments.nz, covariances,
                )
            ]
            j_moments.theta_sem = tuple(result[0] for result in group_sems)
            j_moments.phi_sem = tuple(result[1] for result in group_sems)

    @classmethod
    def attach_norm_spin_components_from_angles(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach normalized spin-direction fields from theta and phi.
        """
        from common.utils.moments import norm_spin_components_from_angles

        if j_moments.theta is None or j_moments.phi is None:
            raise ValueError(
                "Normalized spin-component attachment requires group angles."
            )

        group_results = [
            norm_spin_components_from_angles(theta_g, phi_g)
            for theta_g, phi_g in zip(
                j_moments.theta,
                j_moments.phi,
            )
        ]
        j_moments.nx = tuple(result[0] for result in group_results)
        j_moments.ny = tuple(result[1] for result in group_results)
        j_moments.nz = tuple(result[2] for result in group_results)

    @classmethod
    def attach_spin_components_from_norm_spin_components(
        cls,
        j_moments: "JMomentSeries",
    ) -> None:
        """
        Attach spin-component fields from vector length and normalized direction.
        """
        from common.utils.moments import spin_components_from_norm_spin_components

        if (
            j_moments.length is None
            or j_moments.nx is None
            or j_moments.ny is None
            or j_moments.nz is None
        ):
            raise ValueError(
                "Spin-component attachment requires normalized group "
                "spin-direction fields."
            )

        group_results = [
            spin_components_from_norm_spin_components(length_g, nx_g, ny_g, nz_g)
            for length_g, nx_g, ny_g, nz_g in zip(
                j_moments.length,
                j_moments.nx,
                j_moments.ny,
                j_moments.nz,
            )
        ]
        j_moments.x = tuple(result[0] for result in group_results)
        j_moments.y = tuple(result[1] for result in group_results)
        j_moments.z = tuple(result[2] for result in group_results)

    model_config = ConfigDict(arbitrary_types_allowed=True)
