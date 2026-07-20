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

    model_config = ConfigDict(arbitrary_types_allowed=True)


class HarmonicAnalysisSeries(BaseModel):
    """Harmonic analyses for one or more input curves."""

    signals: tuple[HarmonicSignalAnalysis, ...]

    model_config = ConfigDict(arbitrary_types_allowed=True)


__all__ = ["HarmonicAnalysisSeries", "HarmonicSignalAnalysis"]
