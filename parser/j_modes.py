from __future__ import annotations

from pydantic import BaseModel, ConfigDict

from parser.bloch_vector import BlochVectorSeries
from parser.common import Array


class JModeSeries(BaseModel):
    """Common, contrast, bright, and dark two-group J-vector modes."""

    t: Array
    common: BlochVectorSeries
    contrast: BlochVectorSeries
    bright: BlochVectorSeries
    dark: BlochVectorSeries

    model_config = ConfigDict(arbitrary_types_allowed=True)
