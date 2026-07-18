from __future__ import annotations

from pydantic import BaseModel

from parser.common import Array


class BlochVectorSeries(BaseModel):
    """Cartesian components and length of one derived Bloch vector."""

    x: Array
    y: Array
    z: Array
    length: Array

    class Config:
        arbitrary_types_allowed = True


class JModeSeries(BaseModel):
    """Common, contrast, bright, and dark two-group J-vector modes."""

    t: Array
    common: BlochVectorSeries
    contrast: BlochVectorSeries
    bright: BlochVectorSeries
    dark: BlochVectorSeries

    class Config:
        arbitrary_types_allowed = True
