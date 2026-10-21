from __future__ import annotations

import numpy as np

from parser.common import Array
from parser.generator_moments import TrajectoryResolvedMomentSeries
from parser.harmonic_analysis import DampedOscillationSeries
from parser.moments import MomentSeries
from parser.short_time_oscillations import ShortTimeDarkModeAnalysis
from post_analysis.generator_moments import compute_mean_generator_moments
from post_analysis.harmonic_analysis import (
    compute_bootstrap_damped_oscillations,
    extract_family_phase_series,
)
from post_analysis.j_moments import compute_j_moments
from post_analysis.j_modes import compute_j_modes_from_j_moments


def compute_short_time_dark_mode_analysis(
    moments: MomentSeries,
    *,
    duration: float,
    reference_frequencies: Array,
    num_resamples: int = 200,
    confidence_level: float = 0.95,
    seed: int | None = None,
    verbose: bool = False,
) -> ShortTimeDarkModeAnalysis:
    """Fit early phase-2 dark oscillations using cached trajectory moments.

    Duration uses the simulation's time units and is capped at the end of
    phase 2. Positive reference frequencies set the fitting band to half to
    one-and-a-half times their median, capped below Nyquist. These are plain
    input values; this helper does not run windowed harmonic analysis.

    Resample complete trajectories and repeat the usual generator-mean,
    J-moment, and dark-mode calculations. Return plot-ready curves and SEMs
    without changing the simulation's stored moments or running a solver.
    """
    raw_moments = moments.generator_moments
    metadata = moments.metadata
    modes = moments.J_modes
    if raw_moments is None or raw_moments.num_trajectories < 2:
        raise ValueError("This analysis needs at least two cached trajectories.")
    if metadata is None or modes is None:
        raise ValueError("Simulation metadata and cached J modes are required.")
    if not np.isfinite(duration) or duration <= 0.0:
        raise ValueError("duration must be finite and positive.")
    if not np.array_equal(modes.t, raw_moments.t):
        raise ValueError("Cached J modes and generator moments must share a time grid.")

    dark_components = (modes.dark.x, modes.dark.y, modes.dark.z)
    phase_t, _ = extract_family_phase_series(
        raw_moments.t,
        dark_components,
        metadata.phase_protocol,
        family_phase_index=1,
    )
    if phase_t.size < 8:
        raise ValueError("Phase 2 needs at least eight cached samples.")
    interval_end = min(phase_t[0] + duration, phase_t[-1])
    indices = np.flatnonzero(
        (raw_moments.t >= phase_t[0]) & (raw_moments.t <= interval_end)
    )
    if indices.size < 8:
        raise ValueError("The early-time fit needs at least eight cached samples.")
    time_slice = slice(indices[0], indices[-1] + 1)
    selected_moments = TrajectoryResolvedMomentSeries(
        t=raw_moments.t[time_slice].copy(),
        first_moments=raw_moments.first_moments[:, time_slice].copy(),
        second_moments=raw_moments.second_moments[:, time_slice].copy(),
    )

    frequencies = np.asarray(reference_frequencies, dtype=float)
    valid_frequencies = frequencies[np.isfinite(frequencies) & (frequencies > 0.0)]
    if valid_frequencies.size == 0:
        raise ValueError("At least one positive reference frequency is required.")
    reference_frequency = float(np.median(valid_frequencies))
    dt = selected_moments.t[1] - selected_moments.t[0]
    frequency_bounds = (
        0.5 * reference_frequency,
        min(1.5 * reference_frequency, np.nextafter(0.5 / dt, 0.0)),
    )
    bootstrap_calls = 0

    def dark_curves_from_trajectory_indices(trajectory_indices):
        nonlocal bootstrap_calls
        mean_generators = compute_mean_generator_moments(
            selected_moments, trajectory_indices=trajectory_indices,
        )
        j_moments = compute_j_moments(mean_generators, metadata=metadata)
        resampled_modes = compute_j_modes_from_j_moments(
            j_moments, omega_groups=metadata.omega_groups,
        )
        bootstrap_calls += 1
        if verbose and bootstrap_calls > 1 and (bootstrap_calls - 1) % 25 == 0:
            print(f"Bootstrap: {bootstrap_calls - 1}/{num_resamples}")
        return resampled_modes.dark.x, resampled_modes.dark.y, resampled_modes.dark.z

    if verbose:
        print("Fitting a single exponentially damped sinusoid with a linear background.")
        print("Frequency band:", frequency_bounds, "; intervals are conditional on this model.")
    fits = compute_bootstrap_damped_oscillations(
        selected_moments.t,
        dark_curves_from_trajectory_indices,
        num_trajectories=selected_moments.num_trajectories,
        frequency_bounds=frequency_bounds,
        num_resamples=num_resamples,
        confidence_level=confidence_level,
        seed=seed,
    )
    if verbose:
        _print_fit_summary(fits, Gamma=metadata.Gamma)

    sem_fields = (modes.dark.x_sem, modes.dark.y_sem, modes.dark.z_sem)
    sem = None
    if all(values is not None for values in sem_fields):
        sem = tuple(values[time_slice].copy() for values in sem_fields)
    return ShortTimeDarkModeAnalysis(
        t=selected_moments.t,
        components=tuple(values[time_slice].copy() for values in dark_components),
        fits=fits,
        sem=sem,
    )


def _print_fit_summary(fits: DampedOscillationSeries, *, Gamma: float) -> None:
    """Print component fit quality and model-conditional decay intervals."""
    for component, fit in zip(("x", "y", "z"), fits.signals):
        if fit is None:
            print(f"Dark {component}: no oscillatory fit.")
            continue
        uncertainty = fit.uncertainty
        interval = uncertainty.decay_rate_interval
        scaled_interval = (
            None if interval is None else tuple(value / Gamma for value in interval)
        )
        print(
            f"Dark {component}: gamma/Gamma={fit.decay_rate / Gamma:.5g}, "
            f"f/Gamma={fit.frequency / Gamma:.5g}, R²={fit.fit_r_squared:.3f}"
        )
        print(
            f"  {100 * uncertainty.confidence_level:.0f}% bootstrap interval for "
            f"gamma/Gamma: {scaled_interval}; usable resamples: "
            f"{uncertainty.num_successful}/{uncertainty.num_resamples}; "
            f"parameter bound: {fit.at_parameter_bound}"
        )
        if interval is None:
            print("  The fit is not stable enough to report an interval.")
        elif interval[0] <= 0.0 <= interval[1]:
            print("  The interval includes zero damping.")
        elif interval[0] > 0.0:
            print("  Positive damping within the assumed single-oscillation model.")
        else:
            print("  The fit favors growth; inspect the model and frequency band.")


__all__ = ["compute_short_time_dark_mode_analysis"]
