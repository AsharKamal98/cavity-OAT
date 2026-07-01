from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from parser.common import ObservableSeries
from common.utils import active_manifold_angles, phase_change_times
from quantum_trajectories.aggregator import ensemble_observables, trajectory_observables
from common.utils import omega2_from_weighted_average
from parser.quantum_trajectories import TrajectoryEnsemble, TrajectoryResult


def _observable_series_for_result(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
) -> tuple[ObservableSeries, list, float, np.ndarray]:
    """
    Return averaged observables, phases, Gamma, and phase indices for snapshots.

    This is intentionally a post-processing helper: it reads the saved snapshots
    and does not modify the trajectory simulation path.
    """
    if isinstance(result, TrajectoryEnsemble):
        if not result.trajectories:
            raise ValueError("TrajectoryEnsemble is empty.")
        reference = result.trajectories[0]
        obs = ensemble_observables(result)
        phase_indices = np.asarray([snap.phase_index for snap in reference.snapshots], dtype=int)
        return obs, reference.phases, reference.Gamma, phase_indices

    obs = trajectory_observables(result)
    phase_indices = np.asarray([snap.phase_index for snap in result.snapshots], dtype=int)
    return obs, result.phases, result.Gamma, phase_indices


def _omega_groups_for_result(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
) -> Optional[tuple[float, float]]:
    """
    Read the two inhomogeneous coupling weights stored on the trajectory result.
    """
    reference = result.trajectories[0] if isinstance(result, TrajectoryEnsemble) else result
    if reference.omega_1 is None or reference.omega_2 is None:
        return None
    return float(reference.omega_1), float(reference.omega_2)


def _require_two_group_observables(obs: ObservableSeries) -> tuple[Tuple[np.ndarray, ...], ...]:
    """
    Extract group-resolved observables and fail clearly for homogeneous data.
    """
    if (
        obs.Jx_groups is None
        or obs.Jy_groups is None
        or obs.Jz_groups is None
        or obs.N_e_groups is None
        or len(obs.Jx_groups) != 2
    ):
        raise ValueError(
            "Inhomogeneous MFE residuals require two-group observables. "
            "Run this on a trajectory or ensemble with tuple sectors (Nj1, Nj2)."
        )

    return (
        tuple(np.asarray(arr, dtype=float) for arr in obs.Jx_groups),
        tuple(np.asarray(arr, dtype=float) for arr in obs.Jy_groups),
        tuple(np.asarray(arr, dtype=float) for arr in obs.Jz_groups),
        tuple(np.asarray(arr, dtype=float) for arr in obs.N_e_groups),
    )


def _print_phase_end_residuals(data: dict) -> None:
    """
    Print component-wise residual sums at the saved snapshot closest to each phase end.

    The residuals are computed only on the saved t_eval grid, so the printed
    value is evaluated at the nearest saved snapshot to each phase-end time.
    The summary matches the plotted component curves by summing
    |Re R1| + |Im R1| + |Re R2| + |Im R2| at that snapshot.
    """
    t = np.asarray(data["t"], dtype=float)
    phase_end_times = data.get("phase_end_times")
    if t.size == 0 or phase_end_times is None:
        return

    for phase_idx, phase_end in enumerate(phase_end_times, start=1):
        snapshot_idx = int(np.argmin(np.abs(t - phase_end)))
        residual_sum = (
            abs(np.real(data["R1"][snapshot_idx]))
            + abs(np.imag(data["R1"][snapshot_idx]))
            + abs(np.real(data["R2"][snapshot_idx]))
            + abs(np.imag(data["R2"][snapshot_idx]))
        )
        print(
            "Phase "
            f"{phase_idx}: |Re R1| + |Im R1| + |Re R2| + |Im R2| = {residual_sum:.6e}"
        )


def inhomogeneous_group_angles(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    tol: float = 1e-12,
) -> dict:
    """
    Compute group-resolved and component-averaged active-manifold angles.

    The average angle is computed by first summing the group observables and
    then converting the summed Bloch vector to angles.  This avoids the
    incorrect arithmetic average (theta_1 + theta_2) / 2.
    """
    obs, phases, _, _ = _observable_series_for_result(result)
    jx_groups, jy_groups, jz_groups, ne_groups = _require_two_group_observables(obs)

    theta = []
    phi = []
    n_active = []
    for group_idx in range(2):
        theta_g, phi_g, n_active_g, _, _, _ = active_manifold_angles(
            jx_groups[group_idx],
            jy_groups[group_idx],
            jz_groups[group_idx],
            ne_groups[group_idx],
            tol=tol,
        )
        theta.append(theta_g)
        phi.append(phi_g)
        n_active.append(n_active_g)

    theta_avg, phi_avg, n_active_avg, _, _, _ = active_manifold_angles(
        jx_groups[0] + jx_groups[1],
        jy_groups[0] + jy_groups[1],
        jz_groups[0] + jz_groups[1],
        ne_groups[0] + ne_groups[1],
        tol=tol,
    )

    # TEMP DEBUG: print scaled group angles at the saved t_eval points.
    # Uncomment this block if the omega_1 != omega_2 angle drift diagnostic is needed again.
    # omega_groups = _omega_groups_for_result(result)
    # if omega_groups is not None:
    #     omega1, omega2 = omega_groups
    #     scaled_theta1 = theta[0] / np.sqrt(omega1) if omega1 > 0 else np.full_like(theta[0], np.nan)
    #     scaled_theta2 = theta[1] / np.sqrt(omega2) if omega2 > 0 else np.full_like(theta[1], np.nan)
    #     print("TEMP DEBUG: theta_i / sqrt(omega_i) at saved t_eval points")
    #     print("idx          t    theta_1/sqrt(omega_1)    theta_2/sqrt(omega_2)")
    #     for idx, (time, value1, value2) in enumerate(
    #         zip(np.asarray(obs.t, dtype=float), scaled_theta1, scaled_theta2)
    #     ):
    #         print(f"{idx:3d} {time:10.6f} {value1:24.12e} {value2:24.12e}")

    return {
        "t": np.asarray(obs.t, dtype=float),
        "theta1": theta[0],
        "theta2": theta[1],
        "theta_avg": theta_avg,
        "phi1": phi[0],
        "phi2": phi[1],
        "phi_avg": phi_avg,
        "N_active1": n_active[0],
        "N_active2": n_active[1],
        "N_active_avg": n_active_avg,
        "phase_boundaries": phase_change_times(phases) if len(phases) >= 2 else None,
    }


def inhomogeneous_mfe_residuals(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    omega1: float,
    N1: int,
    N2: int,
    tol: float = 1e-12,
) -> dict:
    """
    Compute the two-group mean-field steady-state residuals from saved data.

    The residual is the algebraic steady-state phase-equation residual described
    in docs/paper_inhomogeneous_couplings.md. It is evaluated from saved
    trajectory snapshots through group-resolved observables; it does not alter
    the MCWF solver or any stored trajectory data.
    """
    if N1 < 0 or N2 < 0:
        raise ValueError("N1 and N2 must be non-negative.")

    obs, phases, Gamma, phase_indices = _observable_series_for_result(result)
    jx_groups, jy_groups, jz_groups, ne_groups = _require_two_group_observables(obs)

    theta = []
    phi = []
    n_active = []
    for group_idx in range(2):
        theta_g, phi_g, n_active_g, _, _, _ = active_manifold_angles(
            jx_groups[group_idx],
            jy_groups[group_idx],
            jz_groups[group_idx],
            ne_groups[group_idx],
            tol=tol,
        )
        theta.append(theta_g)
        phi.append(phi_g)
        n_active.append(n_active_g)

    theta1, theta2 = theta
    phi1, phi2 = phi
    nj1, nj2 = n_active
    omega2 = omega2_from_weighted_average(omega1, N1, N2)

    omega_t = np.asarray([phases[idx].omega for idx in phase_indices], dtype=float)
    delta_t = np.asarray([phases[idx].delta for idx in phase_indices], dtype=float)
    phase_end_times = np.cumsum([phase.duration for phase in phases], dtype=float)
    phase_boundaries = phase_change_times(phases) if len(phases) >= 2 else None

    residuals = []
    weighted_collective_transverse_sum = (
        omega1 * nj1 * np.exp(1j * phi1) * np.sin(theta1)
        + omega2 * nj2 * np.exp(1j * phi2) * np.sin(theta2)
    )

    for theta_g, phi_g, omega_g in ((theta1, phi1, omega1), (theta2, phi2, omega2)):
        sin_theta = np.sin(theta_g)
        cos_theta = np.cos(theta_g)

        with np.errstate(divide="ignore", invalid="ignore"):
            detuning_factor = np.where(
                np.abs(cos_theta) > tol,
                sin_theta * np.tan(theta_g),
                np.nan,
            )

        # Residual from docs/paper_inhomogeneous_couplings.md:
        # R_a = LHS_a - RHS_a.  The theta_a angle stays group-resolved.
        drive_term = 0.5 * omega_t * omega_g * np.exp(-1j * phi_g) * sin_theta
        detuning_term = -0.5 * delta_t * detuning_factor
        decay_term = (
            0.25j
            * Gamma
            * omega_g
            * np.exp(-1j * phi_g)
            * sin_theta
            * weighted_collective_transverse_sum
        )
        residuals.append(drive_term + detuning_term + decay_term)

    return {
        "t": np.asarray(obs.t, dtype=float),
        "R1": residuals[0],
        "R2": residuals[1],
        "theta1": theta1,
        "theta2": theta2,
        "phi1": phi1,
        "phi2": phi2,
        "N_J1": nj1,
        "N_J2": nj2,
        "omega": omega_t,
        "delta": delta_t,
        "phase_end_times": phase_end_times,
        "phase_boundaries": phase_boundaries,
        "omega1": float(omega1),
        "omega2": float(omega2),
        "N1": int(N1),
        "N2": int(N2),
        "Gamma": float(Gamma),
    }


def plot_inhomogeneous_mfe_residuals(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    omega1: float,
    N1: int,
    N2: int,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    tol: float = 1e-12,
):
    """
    Plot signed real and imaginary parts of the two inhomogeneous MFE residuals.

    The grid is
        Re R1    Im R1
        Re R2    Im R2
    and the curves should stay near zero when the saved states are close to the
    inhomogeneous mean-field steady state.
    """
    data = inhomogeneous_mfe_residuals(
        result,
        omega1=omega1,
        N1=N1,
        N2=N2,
        tol=tol,
    )

    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes)
    if axes.shape != (2, 2):
        raise ValueError("axes must have shape (2, 2) for the residual grid.")

    t = data["t"]
    specs = [
        (axes[0, 0], np.real(data["R1"]), r"$\mathrm{Re}\,R_1$"),
        (axes[0, 1], np.imag(data["R1"]), r"$\mathrm{Im}\,R_1$"),
        (axes[1, 0], np.real(data["R2"]), r"$\mathrm{Re}\,R_2$"),
        (axes[1, 1], np.imag(data["R2"]), r"$\mathrm{Im}\,R_2$"),
    ]

    for ax, values, label in specs:
        ax.plot(t, values, linewidth=1.8, label=label)
        ax.axhline(0.0, linestyle=":", color="black", alpha=0.7)
        ax.set_ylabel(label)
        ax.grid(alpha=0.3)
        ax.legend()
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)
        finite_values = np.asarray(values, dtype=float)
        finite_values = finite_values[np.isfinite(finite_values)]
        if finite_values.size == 0:
            ax.set_ylim(-1.0, 1.0)
        else:
            max_abs = float(np.max(np.abs(finite_values)))
            if max_abs <= 0.0:
                max_abs = 1.0
            ax.set_ylim(-max_abs, max_abs)

    if data["phase_boundaries"] is not None:
        t_step1_end, t_step2_end = data["phase_boundaries"]
        for ax in axes.ravel():
            ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.5)
            ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.5)

    axes[1, 0].set_xlabel(r"$\Gamma t$")
    axes[1, 1].set_xlabel(r"$\Gamma t$")
    fig.suptitle("Inhomogeneous MFE steady-state residuals")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    _print_phase_end_residuals(data)

    return data, fig, axes


def plot_inhomogeneous_group_angles(
    result: Union[TrajectoryResult, TrajectoryEnsemble],
    *,
    axes=None,
    output_path: Optional[Union[str, Path]] = None,
    tol: float = 1e-12,
):
    """
    Plot group-resolved and component-averaged active-manifold angles.

    The grid is
        theta_1, theta_2, theta_avg    phi_1, phi_2, phi_avg
    where theta_avg and phi_avg are computed from summed group observables,
    not from arithmetic averaging of the angles.
    """
    data = inhomogeneous_group_angles(result, tol=tol)

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes).ravel()
    if axes.size != 2:
        raise ValueError("axes must contain exactly two axes for the 1x2 angle grid.")

    t = data["t"]
    angle_specs = [
        (
            axes[0],
            (
                (data["theta1"], r"$\theta_1$"),
                (data["theta2"], r"$\theta_2$"),
                (data["theta_avg"], r"$\theta_{\rm avg}$"),
            ),
            r"$\theta$",
        ),
        (
            axes[1],
            (
                (data["phi1"], r"$\phi_1$"),
                (data["phi2"], r"$\phi_2$"),
                (data["phi_avg"], r"$\phi_{\rm avg}$"),
            ),
            r"$\phi$",
        ),
    ]

    for ax, curves, ylabel in angle_specs:
        for values, label in curves:
            ax.plot(t, values, linewidth=1.8, label=label)
        if data["phase_boundaries"] is not None:
            t_step1_end, t_step2_end = data["phase_boundaries"]
            ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.5)
            ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.5)
        ax.set_xlabel(r"$\Gamma t$")
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.3)
        ax.legend()
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)

    fig.suptitle("Inhomogeneous group active-manifold angles")
    fig.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return data, fig, axes
