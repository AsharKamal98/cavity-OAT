from __future__ import annotations

from pydantic import BaseModel

from parser.bloch_vector import BlochVectorSeries
from parser.common import Array


class JModeSeries(BaseModel):
    """Common, contrast, bright, and dark two-group J-vector modes."""

    t: Array
    common: BlochVectorSeries
    contrast: BlochVectorSeries
    bright: BlochVectorSeries
    dark: BlochVectorSeries

    class Config:
        arbitrary_types_allowed = True
