from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy.optimize import minimize_scalar

from common.utils.moments import as_series_tuple
from common.utils.phases import phase_boundary_times
from parser.common import Array, PhaseProtocol
from parser.harmonic_analysis import (
    HarmonicAnalysisSeries,
    HarmonicSignalAnalysis,
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
    if peak_bin == fft_frequencies.size - 1:
        raise ValueError(
            "The fundamental frequency must lie below the Nyquist frequency."
        )

    frequency_step = float(fft_frequencies[1])
    frequency_guess = float(fft_frequencies[peak_bin])
    lower_bound = max(0.5 * frequency_step, frequency_guess - frequency_step)
    upper_bound = min(
        float(fft_frequencies[-1]),
        frequency_guess + frequency_step,
    )
    fitted_harmonic_count = max(
        1,
        min(
            max_harmonic,
            int(np.floor(np.nextafter(fft_frequencies[-1], 0.0) / upper_bound)),
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
    total_harmonic_distortion = float(
        np.linalg.norm(amplitudes[1:]) / fundamental_amplitude
    )
    rms_amplitude = float(np.linalg.norm(amplitudes) / np.sqrt(2.0))

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


__all__ = ["compute_harmonic_analysis", "extract_family_phase_series"]
