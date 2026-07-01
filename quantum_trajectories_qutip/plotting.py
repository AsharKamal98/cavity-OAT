from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MaxNLocator

from common.plotting_legacy import plot_qutip_angles_and_excitation
from common.utils import phase_change_times, phase1_ss_angles_for_nj
from quantum_trajectories_qutip.aggregator import qutip_group_observables


def plot_qutip_group_angles(
    sim_data,
    phases,
    *,
    axes=None,
    output_path=None,
    labels=None,
    show_phase1_ss=True,
    Gamma=None,
):
    """
    Plot subgroup polar and azimuthal angles from QuTiP output.

    The grid follows the custom inhomogeneous group-angle plot: one theta axis
    and one phi axis, with one curve per subgroup. Homogeneous input remains
    usable and simply produces one curve on each axis.
    """
    subgroup_results = qutip_group_observables(sim_data)
    t_step1_end, t_step2_end = phase_change_times(phases)

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure

    axes = np.asarray(axes).ravel()
    if axes.size != 2:
        raise ValueError("axes must contain exactly two axes for the 1x2 angle grid.")

    group_colors = ("tab:blue", "tab:orange")
    if labels is None:
        theta_labels = [rf"$\theta_{idx + 1}$" for idx in range(len(subgroup_results))]
        phi_labels = [rf"$\phi_{idx + 1}$" for idx in range(len(subgroup_results))]
    else:
        theta_labels = list(labels)
        phi_labels = list(labels)

    for idx, group_result in enumerate(subgroup_results):
        obs = group_result.observables
        color = group_colors[idx] if idx < len(group_colors) else None
        axes[0].plot(
            obs.t,
            obs.theta,
            linewidth=1.8,
            color=color,
            label=theta_labels[idx],
        )
        axes[1].plot(
            obs.t,
            obs.phi,
            linewidth=1.8,
            color=color,
            label=phi_labels[idx],
        )

    if show_phase1_ss and Gamma is not None:
        model = sim_data["model"]
        Nj_ref = int(getattr(model, "N_J", sim_data["N"] // 2))

        Omega_ref = float(phases[0].omega)
        theta_ss, phi_ss = phase1_ss_angles_for_nj(Nj_ref, Omega_ref, Gamma)

        if np.isfinite(theta_ss):
            axes[0].plot(
                [0.0, t_step1_end],
                [theta_ss, theta_ss],
                linestyle="--",
                linewidth=2.2,
                color="tab:green",
                alpha=0.95,
                label="phase 1 ss",
            )
        if np.isfinite(phi_ss):
            axes[1].plot(
                [0.0, t_step1_end],
                [phi_ss, phi_ss],
                linestyle="--",
                linewidth=2.2,
                color="tab:green",
                alpha=0.95,
                label="phase 1 ss",
            )

    angle_specs = [
        (axes[0], r"$\theta$"),
        (axes[1], r"$\phi$"),
    ]

    for ax, ylabel in angle_specs:
        ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.5)
        ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.5)
        ax.set_xlabel(r"$\Gamma t$")
        ax.set_ylabel(ylabel)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
        ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)
        ax.grid(alpha=0.3)
        ax.legend()

    fig.suptitle("Inhomogeneous group active-manifold angles")
    fig.tight_layout()
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes


__all__ = [
    "plot_qutip_angles_and_excitation",
    "plot_qutip_group_angles",
]
