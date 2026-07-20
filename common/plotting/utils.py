from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

import numpy as np
from matplotlib.markers import MarkerStyle

from common.utils.phases import phase_boundary_times
from parser.common import PhaseProtocol

FIGURE_FACE_COLOR = "white"
AXES_FACE_COLOR = "white"
GRID_COLOR = "#d7d7d7"
SPINE_COLOR = "#b8b8b8"
GRADIENT_COLOUR_PALETTE = (
    ("#93c5fd", "#60a5fa", "#3b82f6", "#2563eb", "#1d4ed8", "#1e3a8a"),
    ("#fed7aa", "#fdba74", "#fb923c", "#f97316", "#ea580c", "#9a3412"),
    ("#bbf7d0", "#86efac", "#4ade80", "#22c55e", "#16a34a", "#166534"),
    ("#fecaca", "#fca5a5", "#f87171", "#ef4444", "#dc2626", "#991b1b"),
    ("#e5e7eb", "#d1d5db", "#9ca3af", "#6b7280", "#4b5563", "#1f2937"),
    ("#fecdd3", "#fda4af", "#fb7185", "#e11d48", "#be123c", "#881337"),
)
LINESTYLES = ("-", "--", ":", "-.", "none")
SECTOR_CURVE_COLORS = (
    "#0072B2",
    "#D55E00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#E69F00",
    "#000000",
    "#F0E442",
)
PHASE_SHADE_COLORS = ("#efe7bd", "#dcecf2", "#f2ddd2")
PHASE_BOUNDARY_COLOR = "black"


def prepare_figure(fig) -> None:
    fig.patch.set_facecolor(FIGURE_FACE_COLOR)


def style_axis(ax) -> None:
    ax.set_facecolor(AXES_FACE_COLOR)
    ax.grid(color=GRID_COLOR, linewidth=0.8, alpha=0.85)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(SPINE_COLOR)
    ax.spines["bottom"].set_color(SPINE_COLOR)
    ax.tick_params(color=SPINE_COLOR)
    ax.margins(x=0.01)


def format_time_axis(ax) -> None:
    ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 0), useOffset=False)


def colour_palette(
    *,
    colour_family_index: Optional[int] = None,
    shade_index: Optional[int] = None,
) -> Union[tuple[str, ...], str]:
    """Return colors from the shared family-by-shade gradient palette."""
    if colour_family_index is None and shade_index is None:
        default_shade = len(GRADIENT_COLOUR_PALETTE[0]) // 2
        return tuple(family[default_shade] for family in GRADIENT_COLOUR_PALETTE)

    if colour_family_index is None:
        return tuple(
            family[shade_index % len(family)] for family in GRADIENT_COLOUR_PALETTE
        )

    family = GRADIENT_COLOUR_PALETTE[colour_family_index % len(GRADIENT_COLOUR_PALETTE)]
    if shade_index is None:
        return family

    return family[shade_index % len(family)]


def palette_curve_color(palette: Union[tuple[str, ...], str], curve_index: int) -> str:
    """Select one curve color from a shared palette."""
    if isinstance(palette, str):
        return palette
    return palette[curve_index % len(palette)]


def shade_axes_for_family_phase(axes, family_phase_index: int) -> None:
    """Shade axes backgrounds using one shared family-phase color."""
    if family_phase_index < 0:
        raise ValueError("family_phase_index must be non-negative.")
    color = PHASE_SHADE_COLORS[family_phase_index % len(PHASE_SHADE_COLORS)]
    for ax in np.asarray(axes).ravel():
        ax.set_facecolor(color)
        ax.patch.set_alpha(0.35)


def validated_linestyle(linestyle: str | None = "-") -> str:
    if linestyle is None:
        return "none"
    if linestyle not in LINESTYLES:
        raise ValueError(f"linestyle must be one of {LINESTYLES}, got {linestyle!r}.")
    return linestyle


def validated_marker(marker: str | None = None) -> str | None:
    """Return a Matplotlib marker after validating it."""
    if marker is None:
        return None
    try:
        MarkerStyle(marker)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Invalid Matplotlib marker {marker!r}.") from error
    return marker


def sector_curve_color(sector_index: int) -> str:
    return SECTOR_CURVE_COLORS[sector_index % len(SECTOR_CURVE_COLORS)]


def curve_label(base_label: str, *, label: Optional[str]) -> str:
    if label is None:
        return base_label
    return f"{label} {base_label}"


def get_axes(axes, *, n_axes: int, create_figure, error_message: str):
    if axes is None:
        fig, axes = create_figure()
    else:
        axes = np.asarray(axes)
        fig = axes.flat[0].figure
    prepare_figure(fig)

    axes = np.asarray(axes).ravel()
    if axes.size != n_axes:
        raise ValueError(error_message)
    return fig, axes


def set_bottom_figure_legend(
    fig,
    source_ax,
    *,
    max_columns: int = 4,
    row_major: bool = False,
) -> None:
    """Replace panel legends with one combined legend below the figure grid."""
    for ax in fig.axes:
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    for legend in tuple(fig.legends):
        legend.remove()

    handles, labels = source_ax.get_legend_handles_labels()
    if labels:
        column_count = min(len(labels), max_columns)
        if row_major and len(labels) > column_count:
            display_order = tuple(
                index
                for column_index in range(column_count)
                for index in range(column_index, len(labels), column_count)
            )
            handles = [handles[index] for index in display_order]
            labels = [labels[index] for index in display_order]
        fig.legend(
            handles,
            labels,
            loc="outside lower center",
            ncols=column_count,
        )


def add_phase_regions(axes, phase_protocol: PhaseProtocol | None) -> None:
    if phase_protocol is None:
        return

    ends = phase_boundary_times(phase_protocol.family_phases)
    starts = np.concatenate(([0.0], ends[:-1]))
    for ax in np.asarray(axes).ravel():
        has_shading = any(patch.get_gid() == "phase_shading" for patch in ax.patches)
        has_boundaries = any(line.get_gid() == "phase_boundary" for line in ax.lines)

        if not has_shading:
            for phase_index, (start, end) in enumerate(zip(starts, ends)):
                span = ax.axvspan(
                    start,
                    end,
                    color=PHASE_SHADE_COLORS[phase_index % len(PHASE_SHADE_COLORS)],
                    alpha=0.35,
                    zorder=0,
                )
                span.set_gid("phase_shading")

        if not has_boundaries:
            for boundary in starts[1:]:
                line = ax.axvline(
                    boundary,
                    color=PHASE_BOUNDARY_COLOR,
                    linewidth=1.2,
                    linestyle="--",
                    alpha=0.75,
                    zorder=10,
                )
                line.set_gid("phase_boundary")


def save_figure(fig, output_path: Optional[Union[str, Path]]) -> None:
    if output_path is None:
        return

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")


def finish_time_plot(
    fig,
    axes,
    *,
    phase_protocol: PhaseProtocol | None,
    title: str,
    output_path,
    title_y: float = 1.05,
) -> None:
    add_phase_regions(axes, phase_protocol)
    fig.supxlabel(r"$\Gamma t$")
    fig.suptitle(title, y=title_y, fontsize=14)
    save_figure(fig, output_path)
