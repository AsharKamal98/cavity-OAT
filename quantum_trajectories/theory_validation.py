from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from common.utils import omega_c
from quantum_trajectories.parser import Array, TrajectoryEnsemble, TrajectoryResult
from quantum_trajectories.state_helpers import total_norm2


def _build_theoretical_jump_rate_coefficients(
    *,
    delta: float,
    Gamma: float,
    omega_over_omega_c: float,
    N: int,
) -> dict:
    """
    Build the shared small-detuning coefficients for both theory paths.

        l_eff ≈ C + A S_z

    where

        C = (delta / Gamma) tan(theta_tilde)
        A = 2 delta sin(theta_tilde) / (N Gamma cos(theta_tilde)^3)

    with

        sin(theta_tilde) = Omega / Omega_c.

    Both the theory-based path and the theory-approximation path use the same
    coefficient construction. They differ only in how the S_z moments are
    supplied to the final jump-rate formula.
    """
    if Gamma <= 0:
        raise ValueError("Gamma must be positive.")
    if N <= 0:
        raise ValueError("N must be positive.")

    sin_theta = float(omega_over_omega_c)
    if abs(sin_theta) >= 1.0:
        raise ValueError(
            f"Invalid approximation: |Omega/Omega_c| = {abs(sin_theta):.6g} >= 1. "
            "The small-detuning spin-polarized approximation assumes Omega < Omega_c."
        )

    cos_theta = np.sqrt(1.0 - sin_theta**2)
    tan_theta = sin_theta / cos_theta

    C = (delta / Gamma) * tan_theta
    A = 2.0 * delta * sin_theta / (N * Gamma * cos_theta**3)

    return {
        "sin_theta": sin_theta,
        "cos_theta": cos_theta,
        "tan_theta": tan_theta,
        "C": C,
        "A": A,
        "R0": Gamma * abs(C) ** 2,
        "B": Gamma * abs(A) ** 2,
    }


def theoretical_approx_jump_rate_coefficients(
    *,
    delta: float,
    Gamma: float,
    omega_over_omega_c: float,
    N: int,
) -> dict:
    """
    Public wrapper for the shared coefficient builder.

    Kept as the user-facing entry point for the theory-approximation path.
    """
    return _build_theoretical_jump_rate_coefficients(
        delta=delta,
        Gamma=Gamma,
        omega_over_omega_c=omega_over_omega_c,
        N=N,
    )


def _normalized_sector_probabilities(blocks: dict[int, Array]) -> tuple[float, dict[int, float]]:
    """Return the total norm and normalized Nj-sector probabilities for one snapshot."""
    norm2 = total_norm2(blocks)
    if norm2 <= 1e-15:
        return 0.0, {}
    probs = {
        Nj: float(np.vdot(psi, psi).real) / norm2
        for Nj, psi in blocks.items()
        if psi.size != 0
    }
    return norm2, probs


def _sz_moments_from_sector_probabilities(
    sector_probs: dict[int, float],
    *,
    N: int,
) -> dict:
    """
    Compute <S_z>, <S_z^2>, and Var(S_z) from normalized sector probabilities.

    This is the exact same moment definition used by both theory paths; only
    the source of the probabilities changes.
    """
    if not sector_probs:
        return {
            "S_z_mean": 0.0,
            "S_z2_mean": 0.0,
            "S_z_var": 0.0,
        }

    s_mean = 0.0
    s2_mean = 0.0
    for Nj, prob in sector_probs.items():
        s_nj = 0.5 * N - Nj
        s_mean += prob * s_nj
        s2_mean += prob * (s_nj**2)

    s_var = s2_mean - s_mean**2
    return {
        "S_z_mean": float(np.real(s_mean)),
        "S_z2_mean": float(np.real(s2_mean)),
        "S_z_var": float(np.real(s_var)),
    }


def _theoretical_jump_rate_from_moments(
    *,
    coeffs: dict,
    Gamma: float,
    sz_mean: float,
    sz2_mean: float,
) -> float:
    """
    Evaluate the shared theory jump-rate formula from precomputed S_z moments.

    This is the one place where the actual rate expression lives:

        R = Gamma * ( |C|^2 + 2 Re(C* A) <S_z> + |A|^2 <S_z^2> )
    """
    C = coeffs["C"]
    A = coeffs["A"]
    rate = Gamma * (
        abs(C) ** 2
        + 2.0 * np.real(np.conjugate(C) * A) * sz_mean
        + abs(A) ** 2 * sz2_mean
    )
    return float(np.real(rate))


def _sz_moments_from_uniform_half_width(half_width: int) -> dict:
    """
    Compute S_z moments for the theory-approximation path.

    For the equal-sector approximation, S_z = N/2 - N_J is just the negative
    of the symmetric sector offset, so the moments can be computed from the
    same shared machinery after constructing uniform probabilities.
    """
    if half_width < 0:
        raise ValueError("half_width must be >= 0.")

    offsets = list(range(-half_width, half_width + 1))
    if not offsets:
        return {
            "S_z_mean": 0.0,
            "S_z2_mean": 0.0,
            "S_z_var": 0.0,
        }

    prob = 1.0 / len(offsets)
    sz_mean = 0.0
    sz2_mean = 0.0
    for offset in offsets:
        sz_value = -float(offset)
        sz_mean += prob * sz_value
        sz2_mean += prob * (sz_value**2)

    return {
        "S_z_mean": float(sz_mean),
        "S_z2_mean": float(sz2_mean),
        "S_z_var": float(sz2_mean - sz_mean**2),
    }


def theoretical_jump_rates_for_blocks(
    blocks: dict[int, Array],
    *,
    N: int,
    Gamma: float,
    omega: float,
    delta: float,
) -> dict:
    """
    Evaluate the theory-based jump-rate formulas using the simulated sector distribution.

    This path extracts S_z moments from the actual saved sector weights of one
    MCWF snapshot, then feeds those moments into the shared jump-rate formula.
    """
    _, sector_probs = _normalized_sector_probabilities(blocks)
    if not sector_probs:
        return {
            "R_paper_raw": 0.0,
            "S_z_mean": 0.0,
            "S_z2_mean": 0.0,
            "S_z_var": 0.0,
        }

    central_nj = N // 2
    omega_crit = omega_c(central_nj, Gamma)
    if omega_crit <= 0.0:
        raise ValueError("Omega_c must be positive for the theoretical jump-rate calculation.")

    sin_theta_tilde = omega / omega_crit
    if abs(sin_theta_tilde) >= 1.0:
        return {
            "R_paper_raw": np.nan,
            "S_z_mean": np.nan,
            "S_z2_mean": np.nan,
            "S_z_var": np.nan,
        }

    cos_theta_tilde = np.sqrt(1.0 - sin_theta_tilde**2)
    if abs(cos_theta_tilde) <= 1e-15:
        return {
            "R_paper_raw": np.nan,
            "S_z_mean": np.nan,
            "S_z2_mean": np.nan,
            "S_z_var": np.nan,
        }

    coeffs = _build_theoretical_jump_rate_coefficients(
        delta=delta,
        Gamma=Gamma,
        omega_over_omega_c=sin_theta_tilde,
        N=N,
    )
    moments = _sz_moments_from_sector_probabilities(sector_probs, N=N)
    paper_raw = _theoretical_jump_rate_from_moments(
        coeffs=coeffs,
        Gamma=Gamma,
        sz_mean=moments["S_z_mean"],
        sz2_mean=moments["S_z2_mean"],
    )

    return {
        "R_paper_raw": paper_raw,
        "S_z_mean": moments["S_z_mean"],
        "S_z2_mean": moments["S_z2_mean"],
        "S_z_var": moments["S_z_var"],
    }


def theoretical_jump_rate_for_trajectory(result: TrajectoryResult) -> dict:
    """
    Return time series for the theory-based jump-rate formulas.

    The flow is:
    1. read one saved snapshot
    2. extract sector probabilities
    3. convert them into S_z moments
    4. evaluate the shared jump-rate formula
    """
    t = np.array([snap.time for snap in result.snapshots], dtype=float)
    paper_raw = np.zeros_like(t)
    sz_mean = np.zeros_like(t)
    sz2_mean = np.zeros_like(t)
    sz_var = np.zeros_like(t)

    for k, snap in enumerate(result.snapshots):
        phase = result.phases[snap.phase_index]
        approx = theoretical_jump_rates_for_blocks(
            snap.sector_blocks,
            N=result.N,
            Gamma=result.Gamma,
            omega=phase.omega,
            delta=phase.delta,
        )
        paper_raw[k] = approx["R_paper_raw"]
        sz_mean[k] = approx["S_z_mean"]
        sz2_mean[k] = approx["S_z2_mean"]
        sz_var[k] = approx["S_z_var"]

    return {
        "t": t,
        "R_paper_raw": paper_raw,
        "S_z_mean": sz_mean,
        "S_z2_mean": sz2_mean,
        "S_z_var": sz_var,
    }


def _interp_series(t_src: Array, y_src: Array, t_ref: Array) -> Array:
    """Linear interpolation onto a common time grid."""
    t_src = np.asarray(t_src, dtype=float)
    y_src = np.asarray(y_src, dtype=float)
    t_ref = np.asarray(t_ref, dtype=float)

    if t_src.ndim != 1 or y_src.ndim != 1:
        raise ValueError("t_src and y_src must be 1D arrays.")
    if len(t_src) != len(y_src):
        raise ValueError("t_src and y_src must have the same length.")

    t_unique, idx = np.unique(t_src, return_index=True)
    y_unique = y_src[idx]
    return np.interp(t_ref, t_unique, y_unique)


def theoretical_jump_rate_for_ensemble(
    ensemble: TrajectoryEnsemble,
    *,
    reference: str = "first",
) -> dict:
    """
    Average the theory-based jump-rate calculation over an ensemble on the
    shared MCWF t_eval grid.
    """
    if len(ensemble.trajectories) == 0:
        raise ValueError("Ensemble is empty.")
    if reference != "first":
        raise ValueError("Currently only reference='first' is supported.")

    t_ref = np.asarray(
        [snap.time for snap in ensemble.trajectories[0].snapshots],
        dtype=float,
    )

    raw_list = []
    sz_mean_list = []
    sz2_mean_list = []
    sz_var_list = []

    for traj in ensemble.trajectories:
        approx = theoretical_jump_rate_for_trajectory(traj)
        if len(approx["t"]) != len(t_ref) or not np.allclose(approx["t"], t_ref, atol=1e-12, rtol=0.0):
            raise ValueError(
                "All trajectories must share the same t_eval snapshot grid. "
                "Run the ensemble through the common num_snapshots API."
            )
        raw_list.append(np.asarray(approx["R_paper_raw"], dtype=float))
        sz_mean_list.append(np.asarray(approx["S_z_mean"], dtype=float))
        sz2_mean_list.append(np.asarray(approx["S_z2_mean"], dtype=float))
        sz_var_list.append(np.asarray(approx["S_z_var"], dtype=float))

    raw_arr = np.asarray(raw_list, dtype=float)
    sz_mean_arr = np.asarray(sz_mean_list, dtype=float)
    sz2_mean_arr = np.asarray(sz2_mean_list, dtype=float)
    sz_var_arr = np.asarray(sz_var_list, dtype=float)

    return {
        "t": t_ref,
        "R_paper_raw": np.mean(raw_arr, axis=0),
        "S_z_mean": np.mean(sz_mean_arr, axis=0),
        "S_z2_mean": np.mean(sz2_mean_arr, axis=0),
        "S_z_var": np.mean(sz_var_arr, axis=0),
        "R_paper_raw_std": np.std(raw_arr, axis=0, ddof=0),
    }


def theoretical_jump_rate_ensemble_summary(ensemble: TrajectoryEnsemble) -> dict:
    """Compare integrated exact and theory-based jump rates across trajectories."""
    from quantum_trajectories.aggregator import trajectory_observables

    exact_integrals = []
    theoretical_integrals = []
    jump_counts = []

    for traj in ensemble.trajectories:
        obs = trajectory_observables(traj)
        approx = theoretical_jump_rate_for_trajectory(traj)
        exact_integrals.append(float(np.trapezoid(obs.jump_rate, obs.t)))
        theoretical_integrals.append(float(np.trapezoid(approx["R_paper_raw"], approx["t"])))
        jump_counts.append(float(traj.jump_count))

    def summarize(values: list[float]) -> dict:
        arr = np.asarray(values, dtype=float)
        return {
            "mean": float(np.nanmean(arr)),
            "std": float(np.nanstd(arr, ddof=0)),
        }

    return {
        "integrated_exact_rate": summarize(exact_integrals),
        "integrated_theoretical_jump_rate": summarize(theoretical_integrals),
        "jump_count": summarize(jump_counts),
    }


def uniform_sector_variance(half_width: int) -> float:
    """
    Variance of N_J for equal weights over

        N_J = N/2 - dN, ..., N/2, ..., N/2 + dN

    where ``half_width = dN``.
    """
    if half_width < 0:
        raise ValueError("half_width must be >= 0.")

    return _sz_moments_from_uniform_half_width(half_width)["S_z_var"]


def theoretical_approx_jump_count_vs_half_width_data(
    *,
    delta: float,
    Gamma: float,
    N: int,
    omega_over_omega_c: float,
    max_dN: int,
    simulation_time: float = 10.0,
) -> dict:
    """
    Compute the theoretical-approximation jump count vs number of included N_J sectors.

    This path uses the same shared coefficient builder and the same shared
    moment-to-rate formula as the theory-based path. The only difference is
    that the S_z moments come from the closed-form equal-sector approximation
    instead of from simulated trajectory sector weights.
    """
    if max_dN < 0:
        raise ValueError("max_dN must be >= 0.")
    if simulation_time < 0.0:
        raise ValueError("simulation_time must be non-negative.")

    coeffs = theoretical_approx_jump_rate_coefficients(
        delta=delta,
        Gamma=Gamma,
        omega_over_omega_c=omega_over_omega_c,
        N=N,
    )

    dN_values = np.arange(max_dN + 1, dtype=int)
    n_sectors = 2 * dN_values + 1
    sz_mean = np.zeros_like(dN_values, dtype=float)
    sz2_mean = np.zeros_like(dN_values, dtype=float)
    sz_var = np.zeros_like(dN_values, dtype=float)
    jump_rate = np.zeros_like(dN_values, dtype=float)

    for i, dN in enumerate(dN_values):
        moments = _sz_moments_from_uniform_half_width(int(dN))
        sz_mean[i] = moments["S_z_mean"]
        sz2_mean[i] = moments["S_z2_mean"]
        sz_var[i] = moments["S_z_var"]
        jump_rate[i] = _theoretical_jump_rate_from_moments(
            coeffs=coeffs,
            Gamma=Gamma,
            sz_mean=moments["S_z_mean"],
            sz2_mean=moments["S_z2_mean"],
        )

    jump_count = simulation_time * jump_rate

    return {
        "dN_values": dN_values,
        "half_widths": dN_values,
        "n_sectors": n_sectors,
        "S_z_mean": sz_mean,
        "S_z2_mean": sz2_mean,
        "S_z_var": sz_var,
        "jump_rate": jump_rate,
        "var_Nj": sz_var,
        "jump_count": jump_count,
        "coeffs": coeffs,
    }


def plot_theoretical_approx_jump_count_vs_half_width(
    *,
    delta: float,
    Gamma: float,
    N: int,
    omega_over_omega_c: float,
    max_dN: int,
    simulation_time: float = 10.0,
    ax: Optional[plt.Axes] = None,
    label: str | None = None,
):
    """
    Plot the theoretical-approximation jump count vs number of included N_J sectors.

    Uses

        Count_approx(dN) = T * (R0 + B Var(N_J))

    with equal sector weights, so

        Var(N_J) = dN(dN+1)/3.
    """
    if max_dN < 0:
        raise ValueError("max_dN must be >= 0.")
    if simulation_time < 0.0:
        raise ValueError("simulation_time must be non-negative.")

    out = theoretical_approx_jump_count_vs_half_width_data(
        delta=delta,
        Gamma=Gamma,
        omega_over_omega_c=omega_over_omega_c,
        N=N,
        max_dN=max_dN,
        simulation_time=simulation_time,
    )
    coeffs = out["coeffs"]
    print(f"A: {coeffs['A']}")
    print(f"B: {coeffs['B']}")
    dN_values = out["dN_values"]
    n_sectors = out["n_sectors"]
    var_Nj = out["var_Nj"]
    jump_count = out["jump_count"]

    if ax is None:
        fig, ax = plt.subplots(figsize=(7, 4))
    else:
        fig = ax.figure

    if label is None:
        label = (
            rf"$\delta={delta}$, $\Gamma={Gamma}$, "
            rf"$\Omega/\Omega_c={omega_over_omega_c}$"
        )

    ax.plot(n_sectors, jump_count, marker="o", label=label)
    ax.set_xlabel(r"Number of included $N_J$ sectors, $2dN+1$")
    ax.set_ylabel(r"Approx. jump count")
    ax.set_title(r"Theoretical-approximation jump count vs sector number")
    ax.set_ylim(0, 1.1 * max(jump_count))
    ax.grid(True, alpha=0.3)
    ax.legend()

    out["fig"] = fig
    out["ax"] = ax
    return out
