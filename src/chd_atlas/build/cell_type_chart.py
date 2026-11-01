# src/chd_atlas/build/cell_type_chart.py
"""The gene page's cell-type block: one dot-and-whisker row per declared cell type.

Static SVG, no script. One shared log axis from `_AXIS_MIN` to `_AXIS_MAX` CPM
with a gridline at each decade and a shaded band below the dataset's detection
floor. A row's dot sits at the pooled median and its whisker spans the
quartiles across hearts. Its radius is `2.5 + 7.5 * sqrt(share)` (8 in place of
7.5 in the region grid), where share is the fraction of the cell type's cells
with at least one count: it grows with the share but is not proportional to
it, and a share of 0 still draws a dot of radius 2.5 (`_radius`). The tau line
follows the chart. A "By region" fold repeats the cell types against the
dataset's declared regions.

**Every declared cell type gets a row, in the dataset's own order.**
`cell_types.py` omits a cell type no heart reached 10 cells of, so this module
walks `Dataset.cell_types` and writes "no data" where `CellTypeDataset` has no
entry: an absent row would read as a measured zero. A region the payload has no
row for draws "–" for the same reason. A median at or below `_AXIS_MIN`
(including 0) is drawn at the axis minimum, and one above `_AXIS_MAX` at the axis
maximum; its number column still says what it is.

Every coordinate and radius goes through `charts.coordinate`, for the reason
`test_every_geometry_the_pages_publish_has_coordinate_s_fixed_shape` gives,
and every curated string -- label, compartment, token, region -- and the gene
label are escaped, as `charts.py` and `pages.py` escape theirs.
"""

from __future__ import annotations

import html
import math
from collections.abc import Mapping
from typing import Final

from chd_atlas.build.cell_types import CellTypeDataset, CellTypeEntry, CellTypeValues
from chd_atlas.build.charts import coordinate
from chd_atlas.models.dataset import Dataset

# What the fixed copy states of its dataset: `_HEADING`'s "(Farah 2024)", the
# two `pages._CELL_TYPE_*NOTICE` banners' "8 hearts" and `_READING_TEXT`'s "1 CPM".
_DESCRIBED: Final = (
    ("publication", "PMID:38480880"),
    ("n_samples", 8),
    ("detection_floor", 1.0),
)
_HEADING: Final = "Cell types (Farah 2024)"
_CAPTION: Final = "Log axis, 0.1 to 10,000 CPM; a dot's size is its share of cells."
_READING_SUMMARY: Final = "How to read this chart"
_READING_TEXT: Final = (
    "The dot is the median across hearts; the line spans the middle half of hearts (shown "
    "from three hearts). A heart counts only with 10 or more cells of that type. The shaded "
    'band is below this atlas\'s 1 CPM detection floor. "No data" means no heart reached 10 '
    "cells; it is not a low value. Values pool 9 to 15 weeks. "
)
_READING_LINK: Final = '<a href="../methods.html#expression">How to read cell-type figures</a>.'
_NO_DATA: Final = "no data"
_NO_REGION: Final = "–"
_NO_REGION_KEY: Final = "no heart had 10 or more cells there"
# The chart's dots, as an attribute: a CSS `fill-opacity` on `.ct-dot` would
# also override the region grid's per-decade attribute (final review,
# 2026-10-28, when it did).
_DOT_OPACITY: Final = "0.85"

_AXIS_MIN: Final = 0.1
_AXIS_MAX: Final = 10_000.0
_LOG_MIN: Final = math.log10(_AXIS_MIN)
_LOG_MAX: Final = math.log10(_AXIS_MAX)

_ROW_H: Final = 26
_TOP: Final = 40
_LEFT: Final = 190
_AXIS_W: Final = 430
_NUMBERS_W: Final = 230
_REGION_COL_W: Final = 62


def _round(value: float) -> int:
    """Half up, as the figures' reader expects; `round` would go to even."""
    return int(value + 0.5)


def format_median(median: float) -> str:
    """`0` reads `0`; from 100 up an integer with thousands separators; from 10 an
    integer; from 1 one decimal; below 1 two decimals."""
    if median == 0:
        return "0"
    if median >= 100:
        return f"{_round(median):,}"
    if median >= 10:
        return str(_round(median))
    if median >= 1:
        return f"{median:.1f}"
    return f"{median:.2f}"


def _hearts(values: CellTypeValues) -> str:
    hearts = values["n_donors"]
    return f"{hearts} heart" if hearts == 1 else f"{hearts} hearts"


def _numbers(values: CellTypeValues) -> str:
    share = _round(values["fraction_expressing"] * 100)
    return f"{format_median(values['median'])} CPM · {share}% · {_hearts(values)}"


def _title(label: str, values: CellTypeValues, where: str = "") -> str:
    """The dot's tooltip, escaped: `label` and `where` are curated text."""
    spread = (
        ""
        if values["q25"] is None or values["q75"] is None
        else f", quartiles {format_median(values['q25'])} to {format_median(values['q75'])}"
    )
    return html.escape(
        f"{label}{where}: median {format_median(values['median'])} CPM{spread}, "
        f"{_round(values['fraction_expressing'] * 100)}% of {values['n_cells']:,} cells, "
        f"{_hearts(values)}"
    )


def _x(value: float, left: float, width: float) -> float:
    """Log position of `value`, clamped to the axis at both ends.

    A value at or below `_AXIS_MIN` (including 0) draws at the axis minimum,
    and one above `_AXIS_MAX` (10,000 CPM) at the axis maximum. The chart
    marks neither clamp; the row's number column and the dot's title state
    the value itself.
    """
    clamped = min(max(value, _AXIS_MIN), _AXIS_MAX)
    return left + (math.log10(clamped) - _LOG_MIN) / (_LOG_MAX - _LOG_MIN) * width


def _radius(fraction: float, scale: float) -> float:
    return 2.5 + scale * math.sqrt(min(max(fraction, 0.0), 1.0))


def _text(x: float, y: float, css: str, body: str, extra: str = "") -> str:
    """One `<text>`; `body` is escaped here, so a caller passes it raw."""
    return (
        f'<text x="{coordinate(x)}" y="{coordinate(y)}" class="{css}"{extra}>'
        f"{html.escape(body)}</text>"
    )


def _tag(token: str, region: str | None = None) -> str:
    """The `data-` attributes naming a mark's cell type (and region), escaped."""
    tag = f' data-cell-type="{html.escape(token)}"'
    return tag if region is None else tag + f' data-region="{html.escape(region)}"'


def _rows(dataset: Dataset) -> list[tuple[str, str | None, str | None]]:
    """(compartment, token, label) per drawn row; a header row has no token."""
    rows: list[tuple[str, str | None, str | None]] = []
    previous: str | None = None
    for declared in dataset.cell_types:
        if declared.compartment != previous:
            rows.append((declared.compartment, None, None))
            previous = declared.compartment
        rows.append((declared.compartment, declared.token, declared.label))
    return rows


def _svg(width: float, height: float, label: str, body: str) -> str:
    w, h = coordinate(width), coordinate(height)
    return (
        f'<div class="ct-scroll"><svg class="ct-svg" viewBox="0 0 {w} {h}" '
        f'width="{w}" height="{h}" role="img" aria-label="{html.escape(label)}">'
        f"{body}</svg></div>"
    )


def _axis(top: float, bottom: float) -> str:
    parts = []
    y1, y2 = coordinate(top), coordinate(bottom)
    for power in range(int(_LOG_MIN), int(_LOG_MAX) + 1):
        value = 10.0**power
        x = _x(value, _LEFT, _AXIS_W)
        parts.append(
            f'<line x1="{coordinate(x)}" x2="{coordinate(x)}" y1="{y1}" y2="{y2}" class="ct-grid"/>'
        )
        label = "0.1" if value < 1 else f"{int(value):,}"
        parts.append(_text(x, top - 6, "ct-small ct-mid", label))
    parts.append(
        _text(_x(_AXIS_MAX, _LEFT, _AXIS_W), top - 20, "ct-small ct-end", "CPM, log scale")
    )
    return "".join(parts)


def _chart(
    entry: CellTypeDataset, dataset: Dataset, by_token: Mapping[str, CellTypeEntry], gene: str
) -> str:
    rows = _rows(dataset)
    bottom = _TOP + len(rows) * _ROW_H
    column = _LEFT + _AXIS_W + 40
    floor_x = _x(entry["detection_floor"], _LEFT, _AXIS_W)
    left_x = _x(_AXIS_MIN, _LEFT, _AXIS_W)
    top, low = coordinate(_TOP), coordinate(bottom)
    parts = [
        f'<rect x="{coordinate(left_x)}" y="{top}" width="{coordinate(floor_x - left_x)}" '
        f'height="{coordinate(bottom - _TOP)}" class="ct-below"/>',
        _axis(_TOP, bottom),
        f'<line x1="{coordinate(floor_x)}" x2="{coordinate(floor_x)}" y1="{top}" y2="{low}" '
        'class="ct-floor"/>',
    ]
    for index, (compartment, token, label) in enumerate(rows):
        y = _TOP + index * _ROW_H + _ROW_H / 2
        if token is None or label is None:
            parts.append(_text(0, y + 4, "ct-group", compartment))
            continue
        parts.append(_text(12, y + 4, "ct-label", label))
        found = by_token.get(token)
        if found is None:
            parts.append(_text(_LEFT + 4, y + 4, "ct-muted", _NO_DATA))
            continue
        values = found["pooled"]
        cy = coordinate(y)
        if values["q25"] is not None and values["q75"] is not None:
            parts.append(
                f'<line x1="{coordinate(_x(values["q25"], _LEFT, _AXIS_W))}" '
                f'x2="{coordinate(_x(values["q75"], _LEFT, _AXIS_W))}" y1="{cy}" y2="{cy}" '
                f'class="ct-whisker"{_tag(token)}/>'
            )
        parts.append(
            f'<circle cx="{coordinate(_x(values["median"], _LEFT, _AXIS_W))}" cy="{cy}" '
            f'r="{coordinate(_radius(values["fraction_expressing"], 7.5))}" class="ct-dot" '
            f'fill-opacity="{_DOT_OPACITY}"{_tag(token)}>'
            f"<title>{_title(label, values)}</title></circle>"
        )
        parts.append(_text(column, y + 4, "ct-num", _numbers(values), _tag(token)))
    parts.append(_text(column, _TOP - 6, "ct-small", "median · share · hearts"))
    height = bottom + 28
    return _svg(
        _LEFT + _AXIS_W + _NUMBERS_W,
        height,
        f"{gene} expression by cell type",
        "".join(parts),
    )


def _opacity(values: CellTypeValues) -> float:
    """Fill opacity stepped by CPM decade; a value below the floor is an open ring."""
    median = values["median"]
    for bound, opacity in ((10, 0.3), (100, 0.55), (1000, 0.8)):
        if median < bound:
            return opacity
    return 1.0


_LEGEND: Final = (
    ("below floor", None),
    ("under 10", 0.3),
    ("10–100", 0.55),
    ("100–1,000", 0.8),
    ("1,000+ CPM", 1.0),
)


def _region_grid(dataset: Dataset, by_token: Mapping[str, CellTypeEntry], gene: str) -> str:
    rows = _rows(dataset)
    top = 30
    parts = []
    for column, region in enumerate(dataset.regions):
        cx = _LEFT + column * _REGION_COL_W + _REGION_COL_W / 2
        parts.append(_text(cx, top - 10, "ct-small ct-mid", region))
    for index, (compartment, token, label) in enumerate(rows):
        y = top + index * _ROW_H + _ROW_H / 2
        if token is None or label is None:
            parts.append(_text(0, y + 4, "ct-group", compartment))
            continue
        parts.append(_text(12, y + 4, "ct-label", label))
        found = by_token.get(token)
        by_region = {} if found is None else {r["region"]: r for r in found["regions"]}
        for column, region in enumerate(dataset.regions):
            cx = _LEFT + column * _REGION_COL_W + _REGION_COL_W / 2
            cell = by_region.get(region)
            tag = _tag(token, region)
            if cell is None:
                parts.append(_text(cx, y + 4, "ct-muted ct-mid", _NO_REGION, tag))
                continue
            css = "ct-ring" if cell["below_floor"] else "ct-dot"
            opacity = "" if cell["below_floor"] else f' fill-opacity="{_opacity(cell)}"'
            radius = coordinate(_radius(cell["fraction_expressing"], 8))
            title = _title(label, cell, f" ({region})")
            parts.append(
                f'<circle cx="{coordinate(cx)}" cy="{coordinate(y)}" r="{radius}" '
                f'class="{css}"{opacity}{tag}><title>{title}</title></circle>'
            )
    legend_y = top + len(rows) * _ROW_H + 18
    for index, (name, level) in enumerate(_LEGEND):
        cx = 40 + index * 100
        mark = 'class="ct-ring"' if level is None else f'class="ct-dot" fill-opacity="{level}"'
        parts.append(
            f'<circle cx="{coordinate(cx)}" cy="{coordinate(legend_y)}" r="{coordinate(6)}" '
            f"{mark}/>"
        )
        parts.append(_text(cx + 10, legend_y + 4, "ct-small", name))
    # The dash's own key, on a line of its own: "no data" is defined in the
    # reading fold, and this mark is not that phrase.
    dash_y = legend_y + 20
    parts.append(_text(40, dash_y + 4, "ct-muted ct-mid", _NO_REGION))
    parts.append(_text(50, dash_y + 4, "ct-small", _NO_REGION_KEY))
    return _svg(
        _LEFT + len(dataset.regions) * _REGION_COL_W + 10,
        dash_y + 24,
        f"{gene} expression by cell type and region",
        "".join(parts),
    )


def _tau_line(entry: CellTypeDataset, dataset: Dataset) -> str:
    spec = entry["specificity"]
    if spec is None:
        return ""
    labels = {c.token: c.label for c in dataset.cell_types}
    head = f"τ = {spec['tau']:.2f} ({spec['n_cell_types']} cell types, log2(x+1)), "
    highest = spec["highest_in"]
    tail = (
        "with no single highest."
        if highest is None
        else f"highest in {html.escape(labels[highest])}."
    )
    return f'<p class="ct-tau">{head}{tail}</p>'


def cell_type_block(entry: CellTypeDataset, dataset: Dataset, gene: str) -> str:
    """One dataset's cell-type block: heading, caption, chart, tau line, two folds.

    `gene` is what the figures' accessible names call the gene -- its symbol,
    or its id where it has none.
    """
    # The copy here and in the page's banner states these values, so a dataset
    # that differs needs that copy generalised before it can render.
    for field, value in _DESCRIBED:
        if getattr(dataset, field) != value:
            raise ValueError(
                f"dataset {str(dataset.id)!r} has {field} {getattr(dataset, field)!r}; the "
                f"cell-type copy states {value!r}"
            )
    by_token = {item["cell_type"]: item for item in entry["cell_types"]}
    return (
        '<section class="cell-types">'
        f"<h3>{_HEADING}</h3>"
        f'<p class="ct-caption">{_CAPTION}</p>'
        + _chart(entry, dataset, by_token, gene)
        + _tau_line(entry, dataset)
        + "<details><summary>By region</summary>"
        + _region_grid(dataset, by_token, gene)
        + "</details>"
        + f"<details><summary>{_READING_SUMMARY}</summary>"
        + f"<p>{_READING_TEXT}{_READING_LINK}</p></details>"
        + "</section>"
    )
