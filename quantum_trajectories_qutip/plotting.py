import numpy as np
from utils import phase_change_times


def phase1_ss_angles_for_nj(Nj: int, Omega: float, Gamma: float):
    Omega_c = 0.5 * Nj * Gamma
    if Omega_c <= 0:
        raise ValueError("Omega_c must be positive.")
    ratio = Omega / Omega_c
    if abs(ratio) > 1.0:
        return np.nan, np.nan
    cos_theta = np.sqrt(1.0 - ratio**2)
    theta_ss = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    phi_ss = 0.5 * np.pi
    return theta_ss, phi_ss


def plot_qutip_angles_and_excitation(
    qt_data,
    phases,
    *,
    N,
    output_path=None,
    show_phase1_ss=True,
    gamma=None,
):
    import matplotlib.pyplot as plt

    tlist = np.asarray(qt_data["t"], dtype=float)
    Jx = np.asarray(qt_data["Jx"], dtype=float)
    Jy = np.asarray(qt_data["Jy"], dtype=float)
    Jz = np.asarray(qt_data["Jz"], dtype=float)
    N_e = np.asarray(qt_data["N_e"], dtype=float)

    Nj = N // 2

    N_active = 2.0 * (N_e - Jz)
    valid = N_active > 1e-12

    cos_theta = np.zeros_like(Jz)
    cos_theta[valid] = -2.0 * Jz[valid] / N_active[valid]
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    theta = np.zeros_like(Jz)
    theta[valid] = np.arccos(cos_theta[valid])

    phi = np.arctan2(Jy, Jx)
    r_perp = np.sqrt(Jx**2 + Jy**2)
    phi[r_perp < 1e-12] = 0.0

    t_step1_end, t_step2_end = phase_change_times(phases)

    fig, axes = plt.subplots(3, 1, figsize=(8, 10), sharex=False)

    axes[0].plot(tlist, theta, label="qutip", linewidth=1.8)
    axes[1].plot(tlist, phi, label="qutip", linewidth=1.8)
    axes[2].plot(tlist, N_e, label="qutip", linewidth=1.8)

    for ax in axes:
        ax.axvline(t_step1_end, linestyle="--", color="black", alpha=0.6)
        ax.axvline(t_step2_end, linestyle="--", color="black", alpha=0.6)
        ax.grid(alpha=0.3)

    if show_phase1_ss and gamma is not None:
        Omega1 = phases[0].omega
        theta_ss, phi_ss = phase1_ss_angles_for_nj(Nj, Omega1, gamma)

        if np.isfinite(theta_ss):
            axes[0].hlines(
                y=theta_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss ($N_J=N/2$)",
            )
            axes[1].hlines(
                y=phi_ss,
                xmin=0.0,
                xmax=t_step1_end,
                linestyle=":",
                alpha=0.9,
                label=r"phase-1 ss",
            )

    axes[0].set_xlabel(r"$\Gamma t$")
    axes[0].set_ylabel(r"Polar $\theta(t)$")
    axes[0].legend()

    axes[1].set_xlabel(r"$\Gamma t$")
    axes[1].set_ylabel(r"Azimuthal $\phi(t)$")
    axes[1].legend()

    axes[2].set_xlabel(r"$\Gamma t$")
    axes[2].set_ylabel(r"$\langle N_e(t)\rangle$")
    axes[2].legend()

    fig.tight_layout()
    if output_path is not None:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")

    return fig, axes
