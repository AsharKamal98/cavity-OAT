from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from common.utils import active_manifold_angles, phase_change_times, phase1_ss_angles_for_nj


def plot_trajectory_angles_and_excitation(
    result,
    phases,
    *,
    output_path=None,
    show_phase1_ss=False,
    show_spread: bool = False,
    axes=None,
    label=None
):
    obs = result.observables
    tlist = obs.t
    theta_mc = obs.theta
    phi_mc = obs.phi
    ne_mc = obs.N_e
    jump_rate = obs.jump_rate
    t_step1_end, t_step2_end = phase_change_times(phases)

    if axes is None:
        fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False)
        axes = np.asarray(axes)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    if label is None:
        label = "missing label"

    flat_axes = axes.ravel()
    flat_axes[0].plot(tlist, theta_mc, label=label, linewidth=1.8)
    flat_axes[1].plot(tlist, phi_mc, label=label, linewidth=1.8)
    flat_axes[2].plot(tlist, ne_mc, label=label, linewidth=1.8)
    flat_axes[3].plot(tlist, jump_rate, label=label, linewidth=1.8)

    if show_spread:
        spread_specs = [
            (flat_axes[0], None, theta_mc),
            (flat_axes[1], None, phi_mc),
            (flat_axes[2], getattr(obs, "N_e_std", None), ne_mc),
            (flat_axes[3], getattr(obs, "jump_rate_std", None), jump_rate),
        ]
        for ax, std, mean in spread_specs:
            if std is not None:
                ax.fill_between(tlist, mean - std, mean + std, alpha=0.2)

    for ax in flat_axes:
        ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.6)
        ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.6)
        ax.grid(alpha=0.3)

    if show_phase1_ss:
        Nj_ref = result.N // 2
        Omega1 = phases[0].omega
        theta_ss, phi_ss = phase1_ss_angles_for_nj(Nj_ref, Omega1, result.gamma)

        if np.isfinite(theta_ss):
            flat_axes[0].hlines(
                y=theta_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($N_J=N/2$)",
            )
            flat_axes[1].hlines(
                y=phi_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($\phi=\pi/2$)",
            )

    flat_axes[0].set_xlabel(r"$\Gamma t$")
    flat_axes[0].set_ylabel(r"Polar $\theta(t)$")
    flat_axes[0].legend()

    flat_axes[1].set_xlabel(r"$\Gamma t$")
    flat_axes[1].set_ylabel(r"Azimuthal $\phi(t)$")
    flat_axes[1].legend()

    flat_axes[2].set_xlabel(r"$\Gamma t$")
    flat_axes[2].set_ylabel(r"$\langle N_e(t)\rangle$")
    flat_axes[2].legend()

    flat_axes[3].set_xlabel(r"$\Gamma t$")
    flat_axes[3].set_ylabel(r"Jump rate $r(t)$")
    flat_axes[3].legend()

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
    return fig, axes


def plot_qutip_angles_and_excitation(
    qt_data,
    phases,
    *,
    N,
    output_path=None,
    show_phase1_ss=True,
    gamma=None,
):
    tlist = np.asarray(qt_data["t"], dtype=float)
    Jx = np.asarray(qt_data["Jx"], dtype=float)
    Jy = np.asarray(qt_data["Jy"], dtype=float)
    Jz = np.asarray(qt_data["Jz"], dtype=float)
    N_e = np.asarray(qt_data["N_e"], dtype=float)
    jump_rate = np.asarray(qt_data["jump_rate"], dtype=float)

    Nj = N // 2
    theta, phi, _, _, _, _ = active_manifold_angles(Jx, Jy, Jz, N_e)
    t_step1_end, t_step2_end = phase_change_times(phases)

    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=False)
    axes = np.asarray(axes)
    flat_axes = axes.ravel()

    flat_axes[0].plot(tlist, theta, label="qutip", linewidth=1.8)
    flat_axes[1].plot(tlist, phi, label="qutip", linewidth=1.8)
    flat_axes[2].plot(tlist, N_e, label="qutip", linewidth=1.8)
    flat_axes[3].plot(tlist, jump_rate, label="qutip", linewidth=1.8)

    for ax in flat_axes:
        ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.6)
        ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.6)
        ax.grid(alpha=0.3)

    if show_phase1_ss and gamma is not None:
        Omega1 = phases[0].omega
        theta_ss, phi_ss = phase1_ss_angles_for_nj(Nj, Omega1, gamma)

        if np.isfinite(theta_ss):
            flat_axes[0].hlines(
                y=theta_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($N_J=N/2$)",
            )
            flat_axes[1].hlines(
                y=phi_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss",
            )

    flat_axes[0].set_xlabel(r"$\Gamma t$")
    flat_axes[0].set_ylabel(r"Polar $\theta(t)$")
    flat_axes[0].legend()

    flat_axes[1].set_xlabel(r"$\Gamma t$")
    flat_axes[1].set_ylabel(r"Azimuthal $\phi(t)$")
    flat_axes[1].legend()

    flat_axes[2].set_xlabel(r"$\Gamma t$")
    flat_axes[2].set_ylabel(r"$\langle N_e(t)\rangle$")
    flat_axes[2].legend()

    flat_axes[3].set_xlabel(r"$\Gamma t$")
    flat_axes[3].set_ylabel(r"Jump rate $r(t)$")
    flat_axes[3].legend()

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes


def plot_mse_vs_time(
    mse_series_by_label,
    *,
    keys=("Jx", "Jy", "Jz", "N_e"),
    output_path=None,
):
    fig, axes = plt.subplots(len(keys), 1, figsize=(9, 3 * len(keys)), sharex=True)

    if len(keys) == 1:
        axes = [axes]

    for ax, key in zip(axes, keys):
        for label, mse_data in mse_series_by_label.items():
            ax.plot(
                np.asarray(mse_data[key]["t"], dtype=float),
                np.asarray(mse_data[key]["mse_t"], dtype=float),
                linewidth=1.8,
                label=label,
            )

        ax.set_ylabel(f"{key} MSE")
        ax.grid(alpha=0.3)
        ax.legend()

    axes[-1].set_xlabel(r"$\Gamma t$")
    fig.tight_layout()

    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes
