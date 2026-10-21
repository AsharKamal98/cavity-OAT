from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from parser.common import Array


class HarmonicSignalAnalysis(BaseModel):
    """Harmonic content fitted to one real time series."""

    offset: float
    fundamental_frequency: float
    observed_cycles: float
    harmonic_indices: Array
    frequencies: Array
    amplitudes: Array
    rms_amplitude: float
    phase_offsets: Array
    total_harmonic_distortion: float
    fit_r_squared: float = float("nan")

    model_config = ConfigDict(arbitrary_types_allowed=True)


class HarmonicAnalysisSeries(BaseModel):
    """Harmonic analyses for one or more input curves."""

    signals: tuple[HarmonicSignalAnalysis, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)


class WindowedHarmonicAnalysisSeries(BaseModel):
    """Local fit metrics with shape (windows, curves), including for one curve."""

    window_centers: Array
    offsets: Array
    fundamental_amplitudes: Array
    rms_amplitudes: Array
    fundamental_frequencies: Array
    total_harmonic_distortions: Array
    fit_r_squared: Array | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)


class DampedOscillationUncertainty(BaseModel):
    """Whole-trajectory bootstrap intervals for an early-time decay fit."""

    confidence_level: float
    num_resamples: int
    num_successful: int
    amplitude_interval: tuple[float, float] | None = None
    decay_rate_interval: tuple[float, float] | None = None
    frequency_interval: tuple[float, float] | None = None


class DampedOscillationAnalysis(BaseModel):
    """A sinusoid with an exponential envelope and a linear background."""

    t: Array
    fitted_values: Array
    background: Array
    envelope: Array
    offset: float
    background_slope: float
    amplitude: float
    decay_rate: float
    frequency: float
    phase_offset: float
    fit_r_squared: float
    at_parameter_bound: bool
    uncertainty: DampedOscillationUncertainty | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)


class DampedOscillationSeries(BaseModel):
    """Early-time fits in input-curve order; an unavailable fit is None."""

    signals: tuple[DampedOscillationAnalysis | None, ...]


__all__ = [
    "DampedOscillationAnalysis",
    "DampedOscillationSeries",
    "DampedOscillationUncertainty",
    "HarmonicAnalysisSeries",
    "HarmonicSignalAnalysis",
    "WindowedHarmonicAnalysisSeries",
]
