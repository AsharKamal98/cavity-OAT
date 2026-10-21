from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from common.utils.moments import as_series_tuple
from common.utils.phases import phase_boundary_times
from parser.common import Array, PhaseProtocol
from parser.harmonic_analysis import (
    DampedOscillationAnalysis,
    DampedOscillationSeries,
    DampedOscillationUncertainty,
    HarmonicAnalysisSeries,
    HarmonicSignalAnalysis,
    WindowedHarmonicAnalysisSeries,
)

RELATIVE_AMPLITUDE_TOLERANCE = 1e-10


def _validate_time_grid(t: Array) -> tuple[Array, float]:
    """Return a finite, increasing, uniformly sampled time grid and its spacing."""
    t = np.asarray(t, dtype=float)
    if t.ndim != 1 or t.size < 4 or not np.all(np.isfinite(t)):
        raise ValueError("Harmonic analysis requires at least four finite time points.")

    time_steps = np.diff(t)
    if np.any(time_steps <= 0.0) or not np.allclose(
        time_steps,
        time_steps[0],
        rtol=1e-7,
        atol=1e-12,
    ):
        raise ValueError("Harmonic analysis requires a uniformly sampled time grid.")
    return t, float(time_steps[0])


def _harmonic_design_matrix(
    frequency: float,
    shifted_t: Array,
    harmonic_indices: Array,
) -> Array:
    """Build columns for an offset and integer sine/cosine harmonics."""
    columns = [np.ones_like(shifted_t)]
    for harmonic_index in harmonic_indices:
        angle = 2.0 * np.pi * harmonic_index * frequency * shifted_t
        columns.extend((np.cos(angle), np.sin(angle)))
    return np.column_stack(columns)


def _harmonic_fit_residual(
    frequency: float,
    shifted_t: Array,
    values: Array,
    harmonic_indices: Array,
) -> float:
    """Return the least-squares residual for one candidate fundamental."""
    design = _harmonic_design_matrix(frequency, shifted_t, harmonic_indices)
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ coefficients
    return float(np.dot(residual, residual))


def _estimate_fundamental_frequency(
    t: Array,
    values: Array,
    dt: float,
    max_harmonic: int,
) -> float:
    """Return zero for a constant series or estimate its fundamental frequency."""
    centered_values = values - np.mean(values)
    fft_coefficients = np.fft.rfft(centered_values)
    fft_frequencies = np.fft.rfftfreq(values.size, d=dt)

    scale = max(1.0, float(np.max(np.abs(values))))
    if np.max(np.abs(centered_values)) <= RELATIVE_AMPLITUDE_TOLERANCE * scale:
        return 0.0

    peak_bin = int(np.argmax(np.abs(fft_coefficients[1:])) + 1)
    if values.size % 2 == 0 and peak_bin == fft_frequencies.size - 1:
        raise ValueError(
            "The fundamental frequency must lie below the Nyquist frequency."
        )

    frequency_step = float(fft_frequencies[1])
    frequency_below_nyquist = np.nextafter(0.5 / dt, 0.0)
    frequency_guess = float(fft_frequencies[peak_bin])
    lower_bound = max(0.5 * frequency_step, frequency_guess - frequency_step)
    upper_bound = min(
        frequency_below_nyquist,
        frequency_guess + frequency_step,
    )
    fitted_harmonic_count = max(
        1,
        min(
            max_harmonic,
            int(np.floor(frequency_below_nyquist / upper_bound)),
            (values.size - 1) // 2,
        ),
    )
    harmonic_indices = np.arange(1, fitted_harmonic_count + 1, dtype=int)
    shifted_t = t - t[0]
    result = minimize_scalar(
        _harmonic_fit_residual,
        args=(shifted_t, values, harmonic_indices),
        bounds=(lower_bound, upper_bound),
        method="bounded",
        options={"xatol": frequency_step * 1e-10},
    )
    if not result.success:
        raise ValueError("Could not refine the fundamental frequency.")
    return float(result.x)


def _fit_harmonics(
    t: Array,
    values: Array,
    *,
    fundamental_frequency: float,
    dt: float,
    max_harmonic: int,
) -> HarmonicSignalAnalysis:
    """Fit an offset and integer harmonics at the refined fundamental frequency."""
    nyquist_frequency = 0.5 / dt
    frequency_below_nyquist = np.nextafter(nyquist_frequency, 0.0)
    supported_by_frequency = int(
        np.floor(frequency_below_nyquist / fundamental_frequency)
    )
    supported_by_samples = (values.size - 1) // 2
    harmonic_count = min(
        max_harmonic,
        supported_by_frequency,
        supported_by_samples,
    )
    if harmonic_count < 1:
        raise ValueError("No fitted harmonic lies below the Nyquist frequency.")

    harmonic_indices = np.arange(1, harmonic_count + 1, dtype=int)
    frequencies = fundamental_frequency * harmonic_indices
    shifted_t = t - t[0]

    design = _harmonic_design_matrix(
        fundamental_frequency,
        shifted_t,
        harmonic_indices,
    )
    coefficients, _, _, _ = np.linalg.lstsq(design, values, rcond=None)

    cosine_coefficients = coefficients[1::2]
    sine_coefficients = coefficients[2::2]
    amplitudes = np.hypot(cosine_coefficients, sine_coefficients)
    phase_offsets = np.arctan2(-sine_coefficients, cosine_coefficients)

    fundamental_amplitude = float(amplitudes[0])
    if fundamental_amplitude <= np.finfo(float).eps:
        raise ValueError("The fitted fundamental amplitude is zero.")
    total_harmonic_distortion = (
        float(np.linalg.norm(amplitudes[1:]) / fundamental_amplitude)
        if harmonic_count > 1 else np.nan
    )
    rms_amplitude = float(np.linalg.norm(amplitudes) / np.sqrt(2.0))
    residual = values - design @ coefficients
    centered_values = values - np.mean(values)
    fit_r_squared = 1.0 - float(np.dot(residual, residual)) / float(
        np.dot(centered_values, centered_values)
    )

    return HarmonicSignalAnalysis(
        offset=float(coefficients[0]),
        fundamental_frequency=fundamental_frequency,
        observed_cycles=float(fundamental_frequency * (t[-1] - t[0])),
        harmonic_indices=harmonic_indices,
        frequencies=frequencies,
        amplitudes=amplitudes,
        rms_amplitude=rms_amplitude,
        phase_offsets=phase_offsets,
        total_harmonic_distortion=total_harmonic_distortion,
        fit_r_squared=fit_r_squared,
    )


def extract_family_phase_series(
    t: Array,
    series: Sequence[Array] | Array,
    phase_protocol: PhaseProtocol,
    *,
    family_phase_index: int,
) -> tuple[Array, tuple[Array, ...]]:
    """Return the time grid and curves belonging to one family phase."""
    family_phases = phase_protocol.family_phases
    if not 0 <= family_phase_index < len(family_phases):
        raise ValueError("family_phase_index is outside the phase protocol.")

    t = np.asarray(t, dtype=float)
    curves = as_series_tuple(series)
    if t.ndim != 1 or any(np.asarray(curve).shape != t.shape for curve in curves):
        raise ValueError("Each input series must be one-dimensional and match t.")

    phase_ends = phase_boundary_times(family_phases)
    phase_start = (
        0.0
        if family_phase_index == 0
        else phase_ends[family_phase_index - 1]
    )
    phase_end = phase_ends[family_phase_index]
    phase_mask = (t >= phase_start) & (t <= phase_end)
    return t[phase_mask], tuple(np.asarray(curve)[phase_mask] for curve in curves)


def compute_harmonic_analysis(
    t: Array,
    series: Sequence[Array] | Array,
    *,
    max_harmonic: int = 5,
) -> HarmonicAnalysisSeries | None:
    """
    Analyze one or more real time series on a shared time grid.

    The FFT supplies an initial fundamental-frequency estimate. A subsequent
    harmonic fit provides the returned amplitudes and phase offsets, so the
    input need not contain an integer number of oscillations.
    """
    if max_harmonic < 1:
        raise ValueError("max_harmonic must be at least one.")

    t, dt = _validate_time_grid(t)
    signal_series = as_series_tuple(series)

    analyses = []
    for values in signal_series:
        values = np.asarray(values, dtype=float)
        if values.shape != t.shape or not np.all(np.isfinite(values)):
            raise ValueError("Each input series must be finite and match t.")

        fundamental_frequency = _estimate_fundamental_frequency(
            t,
            values,
            dt,
            max_harmonic,
        )
        if fundamental_frequency == 0.0:
            analyses.append(
                HarmonicSignalAnalysis(
                    offset=float(np.mean(values)),
                    fundamental_frequency=0.0,
                    observed_cycles=0.0,
                    harmonic_indices=np.array([], dtype=int),
                    frequencies=np.array([], dtype=float),
                    amplitudes=np.array([], dtype=float),
                    rms_amplitude=0.0,
                    phase_offsets=np.array([], dtype=float),
                    total_harmonic_distortion=np.nan,
                    fit_r_squared=1.0,
                )
            )
            continue

        observed_cycles = fundamental_frequency * (t[-1] - t[0])
        if observed_cycles < 1.0:
            return None
        analyses.append(
            _fit_harmonics(
                t,
                values,
                fundamental_frequency=fundamental_frequency,
                dt=dt,
                max_harmonic=max_harmonic,
            )
        )

    return HarmonicAnalysisSeries(signals=tuple(analyses))


def compute_windowed_harmonic_analysis(
    t: Array,
    series: Sequence[Array] | Array,
    *,
    window_duration: float,
    step_duration: float,
    max_harmonic: int = 5,
    min_fit_r_squared: float | None = None,
) -> WindowedHarmonicAnalysisSeries:
    """Fit each curve in successive complete windows using the existing fitter.

    Durations use the same units as t. Windows start at t[0], include their
    endpoints, and advance by step_duration; incomplete trailing windows are
    omitted. Metric arrays have shape (windows, curves), even for one curve.
    Fewer than four samples or one resolved oscillation leaves NaN metrics.
    Constant curves retain their offset and have zero amplitudes and frequency.
    An optional fit-quality threshold suppresses unresolved oscillation metrics
    while retaining the window mean as a background estimate.
    """
    t, dt = _validate_time_grid(t)
    curves = as_series_tuple(series)
    if any(
        values.shape != t.shape or not np.all(np.isfinite(values))
        for values in curves
    ):
        raise ValueError("Each input series must be finite and match t.")
    if not np.isfinite(window_duration) or not np.isfinite(step_duration):
        raise ValueError("Window and step durations must be finite.")
    if window_duration <= 0.0 or step_duration <= 0.0:
        raise ValueError("Window and step durations must be positive.")
    if max_harmonic < 1:
        raise ValueError("max_harmonic must be at least one.")
    if min_fit_r_squared is not None and (
        not np.isfinite(min_fit_r_squared) or not 0.0 <= min_fit_r_squared <= 1.0
    ):
        raise ValueError("min_fit_r_squared must lie between zero and one.")

    # Include grid-aligned endpoints despite floating-point roundoff.
    time_tolerance = dt * 1e-7
    if window_duration > t[-1] - t[0] + time_tolerance:
        raise ValueError("window_duration must not exceed the supplied time span.")
    if window_duration < 3.0 * dt - time_tolerance:
        raise ValueError("window_duration must span at least four saved samples.")
    last_start_offset = t[-1] - t[0] - window_duration
    window_starts = t[0] + np.arange(
        0.0, last_start_offset + time_tolerance, step_duration
    )
    shape = (window_starts.size, len(curves))
    offsets = np.full(shape, np.nan)
    fundamental_amplitudes = np.full(shape, np.nan)
    rms_amplitudes = np.full(shape, np.nan)
    fundamental_frequencies = np.full(shape, np.nan)
    total_harmonic_distortions = np.full(shape, np.nan)
    fit_r_squared = np.full(shape, np.nan)

    for window_index, start in enumerate(window_starts):
        window_mask = (t >= start - time_tolerance) & (
            t <= start + window_duration + time_tolerance
        )
        window_t = t[window_mask]
        if window_t.size < 4:
            continue
        for curve_index, values in enumerate(curves):
            window_values = values[window_mask]
            if min_fit_r_squared is not None:
                offsets[window_index, curve_index] = np.mean(window_values)
            try:
                analysis = compute_harmonic_analysis(
                    window_t, window_values, max_harmonic=max_harmonic
                )
            except ValueError as error:
                # A Nyquist-edge noise peak is unresolved, not a usable fit.
                if min_fit_r_squared is not None and "Nyquist" in str(error):
                    continue
                raise
            if analysis is None:
                continue
            signal = analysis.signals[0]
            fit_r_squared[window_index, curve_index] = signal.fit_r_squared
            if (
                min_fit_r_squared is not None
                and signal.fit_r_squared < min_fit_r_squared
            ):
                continue
            offsets[window_index, curve_index] = signal.offset
            fundamental_amplitudes[window_index, curve_index] = (
                signal.amplitudes[0] if signal.amplitudes.size else 0.0
            )
            rms_amplitudes[window_index, curve_index] = signal.rms_amplitude
            fundamental_frequencies[window_index, curve_index] = (
                signal.fundamental_frequency
            )
            total_harmonic_distortions[window_index, curve_index] = (
                signal.total_harmonic_distortion
            )

    return WindowedHarmonicAnalysisSeries(
        window_centers=window_starts + window_duration / 2.0,
        offsets=offsets,
        fundamental_amplitudes=fundamental_amplitudes,
        rms_amplitudes=rms_amplitudes,
        fundamental_frequencies=fundamental_frequencies,
        total_harmonic_distortions=total_harmonic_distortions,
        fit_r_squared=fit_r_squared,
    )


def compute_harmonic_quality(
    analysis: WindowedHarmonicAnalysisSeries,
    *,
    sample_spacing: float,
    window_duration: float,
    labels: Sequence[str],
    min_fit_r_squared: float = 0.9,
    minimum_cycles: float = 3.0,
    minimum_samples_per_period: float = 8.0,
    verbose: bool = False,
) -> Array:
    """Return accepted-window masks and optionally print why fits were flagged."""
    if any(
        not np.isfinite(value) or value <= 0.0
        for value in (sample_spacing, window_duration, minimum_cycles, minimum_samples_per_period)
    ):
        raise ValueError("Sampling, window, and cycle settings must be finite and positive.")
    if not np.isfinite(min_fit_r_squared) or not 0.0 <= min_fit_r_squared <= 1.0:
        raise ValueError("min_fit_r_squared must lie between zero and one.")
    frequency = analysis.fundamental_frequencies
    r_squared = analysis.fit_r_squared
    if r_squared is None or r_squared.shape != frequency.shape:
        raise ValueError("Fit R-squared values must match the frequency array.")
    if frequency.ndim != 2 or len(labels) != frequency.shape[1]:
        raise ValueError("Provide one label per harmonic curve.")

    candidate = np.isfinite(frequency) & (frequency > 0.0)
    low_quality = candidate & ~(r_squared >= min_fit_r_squared)
    few_cycles = candidate & (frequency * window_duration < minimum_cycles)
    too_coarse = candidate & (frequency * sample_spacing > 1.0 / minimum_samples_per_period)
    accepted = candidate & ~low_quality & ~few_cycles & ~too_coarse

    if verbose:
        print("Failure counts can overlap. A flagged fit is not a detection or a non-detection.")
        print(f"{'Curve':<17} {'Candidate':>10} {'No fit':>8} {'Low R²':>8} {'Few cycles':>11} {'Too coarse':>11} {'Accepted':>9}")
        for index, label in enumerate(labels):
            print(
                f"{label:<17} {candidate[:, index].sum():>10} "
                f"{(~candidate[:, index]).sum():>8} {low_quality[:, index].sum():>8} "
                f"{few_cycles[:, index].sum():>11} {too_coarse[:, index].sum():>11} "
                f"{accepted[:, index].sum():>9}"
            )
            quality = r_squared[:, index]
            finite_quality = quality[np.isfinite(quality)]
            if finite_quality.size:
                print(
                    "  R² min/median/max:",
                    ", ".join(f"{value:.3f}" for value in np.quantile(finite_quality, (0.0, 0.5, 1.0))),
                )
    return accepted


def compute_damped_oscillation(
    t: Array,
    series: Array,
    *,
    frequency_bounds: tuple[float, float],
) -> DampedOscillationAnalysis | None:
    """Fit c + b*u + exp(-gamma*u)*(a*cos(2*pi*f*u)+d*sin(2*pi*f*u)).

    Frequencies are cycles per time unit. Gamma is signed: zero and growth
    remain admissible, so positive damping is not imposed by the fit. The
    envelope may change by at most exp(10) over the supplied interval.
    """
    t, dt = _validate_time_grid(t)
    values = np.asarray(series)
    if np.iscomplexobj(values):
        raise ValueError("The input series must be real.")
    values = np.asarray(values, dtype=float)
    if values.shape != t.shape or not np.all(np.isfinite(values)):
        raise ValueError("The input series must be finite and match t.")
    if t.size < 8:
        raise ValueError("A damped-oscillation fit requires at least eight samples.")
    lower_frequency, upper_frequency = frequency_bounds
    if not (
        np.isfinite(lower_frequency)
        and np.isfinite(upper_frequency)
        and 0.0 < lower_frequency < upper_frequency < 0.5 / dt
    ):
        raise ValueError(
            "frequency_bounds must be positive, ordered, and below Nyquist."
        )

    duration = float(t[-1] - t[0])
    u = (t - t[0]) / duration
    background_design = np.column_stack((np.ones_like(u), u))
    background_coefficients = np.linalg.lstsq(
        background_design, values, rcond=None
    )[0]
    trend = background_design @ background_coefficients
    centered = values - trend
    scale = float(np.max(np.abs(centered)))
    signal_scale = max(1.0, float(np.max(np.abs(values))))
    if scale <= RELATIVE_AMPLITUDE_TOLERANCE * signal_scale:
        return None
    scaled_values = (values - values.mean()) / scale
    lower_cycles = lower_frequency * duration
    upper_cycles = upper_frequency * duration

    def solve_linear(parameters):
        decay, cycles = parameters
        envelope = np.exp(-decay * u)
        angle = 2.0 * np.pi * cycles * u
        design = np.column_stack(
            (np.ones_like(u), u, envelope * np.cos(angle), envelope * np.sin(angle))
        )
        coefficients = np.linalg.lstsq(design, scaled_values, rcond=None)[0]
        return design, coefficients

    def residual(parameters):
        design, coefficients = solve_linear(parameters)
        return design @ coefficients - scaled_values

    fft_frequencies = np.fft.rfftfreq(t.size, d=dt)
    fft_amplitudes = np.abs(np.fft.rfft(centered))
    in_band = (fft_frequencies >= lower_frequency) & (
        fft_frequencies <= upper_frequency
    )
    if np.any(in_band):
        candidates = np.flatnonzero(in_band)
        peak = candidates[np.argmax(fft_amplitudes[candidates])]
        frequency_guess = fft_frequencies[peak]
    else:
        frequency_guess = 0.5 * (lower_frequency + upper_frequency)

    initial_cycles = np.clip(
        frequency_guess * duration, lower_cycles, upper_cycles
    )
    best = None
    # A few starts distinguish a damped oscillation from a nearby noisy peak.
    for initial_decay in (0.0, 2.0, -2.0):
        fit = least_squares(
            residual,
            (initial_decay, initial_cycles),
            bounds=((-10.0, lower_cycles), (10.0, upper_cycles)),
            max_nfev=300,
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
        )
        if fit.success and (best is None or fit.cost < best.cost):
            best = fit
    if best is None:
        return None

    design, coefficients = solve_linear(best.x)
    coefficients *= scale
    coefficients[0] += values.mean()
    fitted_values = design @ coefficients
    background = coefficients[0] + coefficients[1] * u
    amplitude = float(np.hypot(coefficients[2], coefficients[3]))
    total_variation = np.sum((values - values.mean()) ** 2)
    r_squared = 1.0 - float(np.sum((values - fitted_values) ** 2) / total_variation)
    at_bound = bool(
        abs(best.x[0]) >= 10.0 - 1e-4
        or best.x[1] - lower_cycles <= 1e-4
        or upper_cycles - best.x[1] <= 1e-4
    )
    return DampedOscillationAnalysis(
        t=t.copy(),
        fitted_values=fitted_values,
        background=background,
        envelope=amplitude * np.exp(-best.x[0] * u),
        offset=float(coefficients[0]),
        background_slope=float(coefficients[1] / duration),
        amplitude=amplitude,
        decay_rate=float(best.x[0] / duration),
        frequency=float(best.x[1] / duration),
        phase_offset=float(np.arctan2(-coefficients[3], coefficients[2])),
        fit_r_squared=r_squared,
        at_parameter_bound=at_bound,
    )


def compute_bootstrap_damped_oscillations(
    t: Array,
    series_from_trajectory_indices: Callable[[Array], Sequence[Array] | Array],
    *,
    num_trajectories: int,
    frequency_bounds: tuple[float, float],
    num_resamples: int = 200,
    confidence_level: float = 0.95,
    seed: int | None = None,
) -> DampedOscillationSeries:
    """Refit observables after resampling complete independent trajectories.

    The callback receives trajectory indices, including repetitions, and must
    recompute ensemble observables in their usual averaging order. Percentile
    intervals are conditional on this single-decay model and frequency band.
    A bound-hitting fit is excluded from the intervals, and at least 80% of
    resamples must succeed before any interval is reported.
    """
    if num_trajectories < 2 or num_resamples < 2:
        raise ValueError(
            "Bootstrap fitting requires at least two trajectories and resamples."
        )
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must lie between zero and one.")
    t, _ = _validate_time_grid(t)
    curves = as_series_tuple(
        series_from_trajectory_indices(np.arange(num_trajectories))
    )
    fits = tuple(
        compute_damped_oscillation(t, curve, frequency_bounds=frequency_bounds)
        for curve in curves
    )
    parameters = np.full((num_resamples, len(fits), 3), np.nan)
    rng = np.random.default_rng(seed)
    for resample in range(num_resamples):
        indices = rng.integers(0, num_trajectories, size=num_trajectories)
        resampled_curves = as_series_tuple(series_from_trajectory_indices(indices))
        if len(resampled_curves) != len(fits):
            raise ValueError("The resampling callback must preserve the number of curves.")
        for curve_index, curve in enumerate(resampled_curves):
            fit = compute_damped_oscillation(
                t, curve, frequency_bounds=frequency_bounds
            )
            if fit is not None and not fit.at_parameter_bound:
                parameters[resample, curve_index] = (
                    fit.amplitude, fit.decay_rate, fit.frequency
                )

    tail_probability = 0.5 * (1.0 - confidence_level)
    for curve_index, fit in enumerate(fits):
        if fit is None:
            continue
        samples = parameters[:, curve_index]
        samples = samples[np.all(np.isfinite(samples), axis=1)]
        uncertainty = DampedOscillationUncertainty(
            confidence_level=confidence_level,
            num_resamples=num_resamples,
            num_successful=len(samples),
        )
        minimum_successful = max(2, int(np.ceil(0.8 * num_resamples)))
        if not fit.at_parameter_bound and len(samples) >= minimum_successful:
            intervals = np.quantile(
                samples, (tail_probability, 1.0 - tail_probability), axis=0
            )
            uncertainty.amplitude_interval = tuple(float(v) for v in intervals[:, 0])
            uncertainty.decay_rate_interval = tuple(float(v) for v in intervals[:, 1])
            uncertainty.frequency_interval = tuple(float(v) for v in intervals[:, 2])
        fit.uncertainty = uncertainty
    return DampedOscillationSeries(signals=fits)


__all__ = [
    "compute_bootstrap_damped_oscillations",
    "compute_damped_oscillation",
    "compute_harmonic_analysis",
    "compute_harmonic_quality",
    "compute_windowed_harmonic_analysis",
    "extract_family_phase_series",
]
