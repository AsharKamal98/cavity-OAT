from __future__ import annotations

from typing import Iterable, List, Mapping, Sequence, Tuple

import numpy as np

from parser.common import Array, Phase


def phase_change_times(phases: Sequence[Phase]) -> Tuple[float, float]:
    if len(phases) < 2:
        raise ValueError("Need at least two phases to define change times.")
    t1 = phases[0].duration
    t2 = phases[0].duration + phases[1].duration
    return t1, t2


def phase_values_at_time(t: float, phases: Sequence[Phase]) -> Tuple[float, float]:
    """
    Return phase-local (Omega, delta) values for a piecewise-constant protocol.
    """
    if not phases:
        raise ValueError("Need at least one phase.")

    t_value = float(t)
    total_time = float(sum(phase.duration for phase in phases))
    if t_value < 0.0 or t_value > total_time:
        raise ValueError(f"t must lie in [0, {total_time}], got {t_value}.")

    phase_end = 0.0
    for index, phase in enumerate(phases):
        phase_end += phase.duration
        if t_value <= phase_end or index == len(phases) - 1:
            return float(phase.omega), float(phase.delta)

    # The loop always returns for non-empty phases, but keep type checkers happy.
    phase = phases[-1]
    return float(phase.omega), float(phase.delta)


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


def active_manifold_angles(
    Jx: Array,
    Jy: Array,
    Jz: Array,
    N_e: Array,
    *,
    tol: float = 1e-12,
) -> Tuple[Array, Array, Array, Array, Array, Array]:
    """
    Compute Bloch-sphere angles inside the active {|down>, |e>} manifold.

    The active-manifold population is
        N_active = <N_down + N_e> = 2 (N_e - J_z),
    since J_z = (N_e - N_down) / 2 in each fixed-Nj sector.

    We then normalize the collective spin components by N_active so that the
    returned (sx, sy, sz) correspond to the averaged single-particle Bloch
    vector within the active manifold, with sz using the same sign convention
    as J_z. The polar angle uses theta = arccos(-sz), so the |down> state
    points to the north pole in the J-Bloch convention.
    """
    Jx = np.asarray(Jx, dtype=float)
    Jy = np.asarray(Jy, dtype=float)
    Jz = np.asarray(Jz, dtype=float)
    N_e = np.asarray(N_e, dtype=float)

    N_active = 2.0 * (N_e - Jz)

    sx = np.zeros_like(Jx, dtype=float)
    sy = np.zeros_like(Jy, dtype=float)
    sz = np.zeros_like(Jz, dtype=float)

    valid = N_active > tol
    sx[valid] = 2.0 * Jx[valid] / N_active[valid]
    sy[valid] = 2.0 * Jy[valid] / N_active[valid]
    sz[valid] = 2.0 * Jz[valid] / N_active[valid]

    if np.any(sz < -1.0 - tol) or np.any(sz > 1.0 + tol):
        raise ValueError("sz values must lie in [-1, 1] to compute angles.")
    sz = np.clip(sz, -1.0, 1.0)

    theta, phi = angles_from_norm_spin_components(sx, sy, sz, valid=valid, tol=tol)
    return theta, phi, N_active, sx, sy, sz


def norm_spin_components(
    x: Array,
    y: Array,
    z: Array,
    *,
    tol: float = 1e-12,
) -> Tuple[Array, Array, Array, Array]:
    """
    Compute Euclidean length, normalized direction, and angles of a spin vector.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    len = np.sqrt(x**2 + y**2 + z**2)
    valid = len > tol

    sx = np.zeros_like(x, dtype=float)
    sy = np.zeros_like(y, dtype=float)
    sz = np.zeros_like(z, dtype=float)
    sx[valid] = x[valid] / len[valid]
    sy[valid] = y[valid] / len[valid]
    sz[valid] = z[valid] / len[valid]

    return len, sx, sy, sz


def angles_from_norm_spin_components(
    sx: Array,
    sy: Array,
    sz: Array,
    valid: Array,
    tol: float = 1e-12,
) -> Tuple[Array, Array]:
    sx = np.asarray(sx, dtype=float)
    sy = np.asarray(sy, dtype=float)
    sz = np.asarray(sz, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    theta = np.zeros_like(sz, dtype=float)
    theta[valid] = np.arccos(np.clip(-sz[valid], -1.0, 1.0))

    phi = np.arctan2(sy, sx)
    r_perp = np.sqrt(sx**2 + sy**2)
    phi[r_perp < tol] = 0.0

    return theta, phi


def observable_mse_by_time(
    candidate,
    reference,
    *,
    keys: Iterable[str] = ("Jx", "Jy", "Jz", "N_e"),
):
    """
    Compute per-timestep MSE between two observable containers.

    The reference series are linearly interpolated onto the candidate time grid,
    so the returned time array always matches candidate.observables.t.
    """
    o_c = candidate.observables
    o_r = reference.observables

    t_c = np.asarray(o_c.t, dtype=float)
    t_r = np.asarray(o_r.t, dtype=float)

    out = {}
    for key in keys:
        y_c = np.asarray(getattr(o_c, key), dtype=float)
        y_r = np.asarray(getattr(o_r, key), dtype=float)
        y_r_interp = np.interp(t_c, t_r, y_r)
        mse_t = (y_c - y_r_interp) ** 2

        out[key] = {
            "t": t_c,
            "mse_t": mse_t,
            "mean_mse": float(np.mean(mse_t)),
            "integrated_mse": float(np.trapezoid(mse_t, t_c) / (t_c[-1] - t_c[0])),
        }

    return out


def omega_c(N_J: int, Gamma: float) -> float:
    """Critical drive for the polarized-to-mixed transition at delta = 0."""

    return 0.5 * N_J * Gamma


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

    From the cavity-elimination notes,

        Omega = 4 * epsilon * g_c / kappa,
        Gamma = 4 * g_c**2 / kappa.

    The bad-cavity approximation assumes

        kappa >> sqrt(N_J) * g_c, delta.

    Since ``>>`` is not a sharp mathematical threshold, ``bad_cavity_factor``
    sets the requested scale separation. The default requires

        kappa >= bad_cavity_factor * max(sqrt(N_J) * |g_c|, |delta|).

    If the check fails, the function prints the offending values and exits the
    current run with ``SystemExit(1)`` so notebook cells can stop immediately.
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
    Enforce the notebook MCWF timestep rule

        dt <= (N * Gamma)^(-1) / safety_factor.

    If the proposed ``dt`` is already valid it is returned unchanged. If it is
    too large, a warning is printed and the largest allowed timestep is
    returned instead.
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

    This is useful for pre-validating an entire half-width sweep: if the check
    passes for the smallest Nj in the widest initial support, it will also pass
    for every narrower support centered at the same N/2.
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


def default_three_phase_protocol(
    T1: float,
    T2: float,
    T3: float,
    delta0: float,
    Omega0: float,
) -> List[Phase]:
    """Three-phase protocol."""
    return [
        Phase(duration=T1, omega=Omega0, delta=0.0, label="phase1"),
        Phase(duration=T2, omega=Omega0, delta=delta0, label="phase2"),
        Phase(duration=T3, omega=0.0, delta=0.0, label="phase3"),
    ]
