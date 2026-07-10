from __future__ import annotations

from typing import Mapping, Sequence, Tuple

import numpy as np

def omega_G_from_weighted_average(
    omega_i: Sequence[float],
    N_i: Sequence[int],
    *,
    tol: float = 1e-12,
) -> float:
    """
    Choose the final group coupling so the atom-number weighted mean coupling is one.
    """
    N_G = float(N_i[-1])
    if abs(N_G) <= tol:
        return 1.0
    N_total = float(sum(N_i))
    weighted_partial_sum = float(
        sum(float(omega_g) * float(N_g) for omega_g, N_g in zip(omega_i, N_i[:-1]))
    )
    return float((N_total - weighted_partial_sum) / N_G)


def omega_c(N_J: int, Gamma: float) -> float:
    """Critical drive for the polarized-to-mixed transition at delta = 0."""
    return 0.5 * N_J * Gamma


def scaled_N_Gamma(factor: float, N: int, Gamma: float) -> float:
    """Scale a dimensionless factor by N * Gamma."""
    return float(factor) * float(N) * float(Gamma)


def inverse_scaled_N_Gamma(factor: float, N: int, Gamma: float) -> float:
    """Scale a dimensionless factor by 1 / (N * Gamma)."""
    return float(factor) / (float(N) * float(Gamma))


def mcwf_dt_from_scales(
    Omega0: float,
    delta0: float,
    N: int,
    Gamma: float,
    *,
    drive_factor: float = 0.05,
    decay_factor: float = 0.1,
) -> float:
    """Choose an MCWF timestep from drive, detuning, and collective decay scales."""
    scales = [decay_factor / (float(N) * float(Gamma))]
    if abs(Omega0) > 0.0:
        scales.append(drive_factor / abs(float(Omega0)))
    if abs(delta0) > 0.0:
        scales.append(drive_factor / abs(float(delta0)))
    return min(scales)


def delta0_from_N_Gamma(N: int, Gamma: float) -> float:
    """Protocol detuning used in the notebook scans: delta = 0.05 * N * Gamma."""
    return 0.05 * N * Gamma


def Omega0_from_N_Gamma(N: int, Gamma: float) -> float:
    """Protocol drive used in the notebook scans: Omega = 0.465 * N * Gamma."""
    return 0.465 * N * Gamma


def Omega_Gamma_from_cavity_parameters(
    epsilon: float,
    g_c: float,
    kappa: float,
    N_J: int,
    delta: float = 0.0,
    *,
    bad_cavity_factor: float = 10.0,
    round_digits: int = 6,
) -> Tuple[float, float]:
    """
    Convert driven-cavity parameters to the effective spin-model parameters.
    """
    if g_c <= 0.0:
        raise ValueError("g_c must be positive.")
    if kappa <= 0.0:
        raise ValueError("kappa must be positive.")
    if N_J < 0:
        raise ValueError("N_J must be non-negative.")
    if bad_cavity_factor <= 0.0:
        raise ValueError("bad_cavity_factor must be positive.")
    if round_digits < 0:
        raise ValueError("round_digits must be non-negative.")

    Omega = 4.0 * epsilon * g_c / kappa
    Gamma = 4.0 * g_c**2 / kappa

    collective_coupling_scale = np.sqrt(float(N_J)) * abs(g_c)
    detuning_scale = abs(delta)
    bad_cavity_scale = max(collective_coupling_scale, detuning_scale)
    required_kappa = bad_cavity_factor * bad_cavity_scale

    value_fmt = f".{round_digits}g"
    print(
        "Cavity-derived effective parameters: "
        f"Omega={Omega:{value_fmt}}, Gamma={Gamma:{value_fmt}}"
    )
    print(
        "Bad-cavity check: "
        f"kappa={kappa:{value_fmt}}, "
        f"sqrt(N_J)*g_c={collective_coupling_scale:{value_fmt}}, "
        f"|delta|={detuning_scale:{value_fmt}}, "
        f"required kappa>={required_kappa:{value_fmt}} "
        f"(factor={bad_cavity_factor:{value_fmt}})."
    )

    if kappa < required_kappa:
        print(
            "Bad-cavity limit not satisfied; exiting before using the "
            "adiabatically eliminated spin model."
        )
        raise SystemExit(1)

    print("Bad-cavity limit satisfied.")
    return float(Omega), float(Gamma)


def validated_mcwf_dt(
    dt: float,
    N: int,
    Gamma: float,
    *,
    safety_factor: float = 250.0,
) -> float:
    """
    Enforce the notebook MCWF timestep rule.
    """
    if dt <= 0.0:
        raise ValueError("dt must be positive.")
    if N <= 0:
        raise ValueError("N must be positive.")
    if Gamma <= 0.0:
        raise ValueError("Gamma must be positive.")
    if safety_factor <= 0.0:
        raise ValueError("safety_factor must be positive.")

    max_dt = 1.0 / (safety_factor * N * Gamma)
    if dt <= max_dt:
        return float(dt)

    print(
        "Warning: proposed dt was too large; using dt="
        f"{max_dt} instead of dt={dt} based on dt <= (N*Gamma)^(-1)/{safety_factor:.0f}."
    )
    return float(max_dt)


def check_initial_sector_omega_ratio(
    sector_coeffs: Mapping[object, complex],
    Omega: float,
    Gamma: float,
    *,
    ratio_limit: float = 1.0,
) -> dict:
    """
    Validate Omega / Omega_c against the smallest populated initial Nj sector.
    """
    if not sector_coeffs:
        raise ValueError("sector_coeffs must contain at least one populated sector.")

    min_nj = min(
        int(sum(sector_key)) if isinstance(sector_key, tuple) else int(sector_key)
        for sector_key in sector_coeffs
    )
    omega_crit = omega_c(min_nj, Gamma)

    if omega_crit <= 0.0:
        return {
            "is_valid": False,
            "min_nj": min_nj,
            "omega": float(Omega),
            "omega_c": float(omega_crit),
            "ratio": np.inf,
            "ratio_limit": float(ratio_limit),
        }

    ratio = float(Omega / omega_crit)
    return {
        "is_valid": bool(abs(ratio) < ratio_limit),
        "min_nj": min_nj,
        "omega": float(Omega),
        "omega_c": float(omega_crit),
        "ratio": ratio,
        "ratio_limit": float(ratio_limit),
    }
