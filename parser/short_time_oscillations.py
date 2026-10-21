from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from parser.common import Array
from parser.harmonic_analysis import DampedOscillationSeries


class ShortTimeDarkModeAnalysis(BaseModel):
    """Early phase-2 dark x/y/z curves, optional SEMs, and bootstrap fits."""

    t: Array
    components: tuple[Array, Array, Array]
    fits: DampedOscillationSeries
    sem: tuple[Array, Array, Array] | None = None

    model_config = ConfigDict(arbitrary_types_allowed=True)


__all__ = ["ShortTimeDarkModeAnalysis"]
