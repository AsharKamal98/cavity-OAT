from __future__ import annotations

COLOUR_PALETTES = (
    ("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"),
    ("#800000", "#556b2f", "#2f4f4f", "#7a3e9d"),
    ("#56B4E9", "#E69F00", "#F0E442", "#8c564b"),
    ("#4d4d4d", "#7f7f7f", "#a6a6a6", "#c0c0c0"),
    ("#4d4d4d", "#4d4d4d", "#4d4d4d", "#4d4d4d"),
    ("#800000", "#800000", "#800000", "#800000"),
)

LINESTYLES = ("-", "--", ":")


def indexed_curve_color(colour_index: int = 0, curve_index: int = 0) -> str:
    palette = COLOUR_PALETTES[colour_index % len(COLOUR_PALETTES)]
    return palette[curve_index % len(palette)]


def validated_linestyle(linestyle: str = "-") -> str:
    if linestyle not in LINESTYLES:
        raise ValueError(f"linestyle must be one of {LINESTYLES}, got {linestyle!r}.")
    return linestyle
