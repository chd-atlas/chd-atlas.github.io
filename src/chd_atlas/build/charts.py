# src/chd_atlas/build/charts.py
"""Inline-SVG primitives. This module does not know what a gene is.

`render.py` turns text into a page and knows nothing about genes; this turns
numbers into geometry on the same terms. `pages.py` is the only module that
joins either to data, which is what keeps 2,000 lines of copy and a coordinate
system from growing into each other.

**No primitive here emits a colour.** Every stroke and fill comes from a class
defined once in `render.STYLESHEET`, against the custom properties that already
carry the page's `prefers-color-scheme: dark` palette. A hardcoded hex would
render a light chart on a dark page for every reader who has that preference
set, and nothing in the suite would fail.

Deterministic for the same reason every other builder here is: nothing reads a
clock, and every coordinate is formatted rather than `repr`'d, so two builds of
one commit are byte-identical.

**There is deliberately no primitive, and no `chart-diamond` class, for a
forest plot's summary diamond.** A diamond draws a *pooled* estimate, and D33
forbids one: `pages._POOLING_NOTICE` states the reason in full and
`build/burden.py::shared_cohorts` is what acts on it. Defining the class would
be the first half of that defect, and the forest plot's tests assert it never
appears -- so a caller wanting one has a data decision to make, not a missing
shape to add. Recorded here rather than in `render.STYLESHEET`, because that
string is inlined verbatim into every page: a comment there is published bytes,
and the token a page-wide check greps for must not be among them.
"""

from __future__ import annotations

import html
import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal


def coordinate(value: float) -> str:
    """One axis coordinate, fixed to a tenth of a pixel.

    A signed zero is normalised **after** formatting, not by adding zero
    before it, and the difference is load-bearing. Measured 2026-08-20:
    `-0.0 + 0.0` is `+0.0`, so adding zero handles an exact negative zero --
    but `-1e-9 + 0.0` is still `-1e-9`, and `f"{-1e-9:.1f}"` is `"-0.0"`. Any
    coordinate a hair below the axis origin formats to a signed zero, so the
    pre-format trick fixes the case nobody hits and misses the case a
    subtraction actually produces. `.1f` emits one decimal place, so `"-0.0"`
    is the only string form a negative zero can take here.

    It matters because `-0.0` and `0.0` are the same point: a subtraction
    landing on one rather than the other would change published bytes on the
    strength of operand order alone, and two builds of one commit must be
    byte-identical.
    """
    if not math.isfinite(value):
        raise ValueError(
            f"a coordinate must be finite, got {value!r}; `encode_json` refuses a "
            f"non-finite float for the JSON side and this is its counterpart for "
            f"the geometry, which that guard cannot reach"
        )
    text = f"{value:.1f}"
    return "0.0" if text == "-0.0" else text


@dataclass(frozen=True)
class LogScale:
    """Maps a positive value onto a pixel axis, base 10.

    Log rather than linear because both axes this serves span orders of
    magnitude: one gene's heart median runs 5 to 275 tpm across development,
    and published odds ratios run 0.0 to unbounded. On a linear axis the
    post-natal half of a trajectory collapses onto the floor.
    """

    low: float
    high: float
    left: float
    width: float

    def __post_init__(self) -> None:
        """Refuses an axis that cannot place anything, at construction.

        Each of the three costs something different, and one of them costs
        silence. Measured 2026-08-20 against the arithmetic in `x` below:

        - `low <= 0` raises `ValueError: math domain error` from `log10`, but
          only once a caller draws its first point, and the message names
          neither the axis nor the argument.
        - `high == low` gives a zero span and raises `ZeroDivisionError:
          float division by zero`, again at draw time and again opaquely.
          `high < low` raises nothing at all: the span goes negative and the
          axis silently inverts, so a rising trajectory is drawn falling.
        - `width <= 0` **raises nothing and returns `left` for every value.**
          Every measurement stacks on one pixel, which reads as "all these
          numbers are the same" rather than as a broken axis.
        - a **non-finite** bound. Added 2026-08-21 after a review measured it:
          `high=inf` gives an infinite span, `finite/inf` is `0.0`, and
          `x(1.0)`, `x(1e6)` and `x(1e300)` all returned `0.0` -- the zero-width
          failure exactly, through a guard that did not look for it. A `nan`
          bound passed every check below, because `nan <= x` is `False` for
          any x.

        Two of the four are silent (`high < low` and the non-finite bounds),
        not one; this paragraph said one until the same review.

        `raise`, never `assert` -- `-O` strips `assert`.
        """
        bounds = (self.low, self.high, self.left, self.width)
        if not all(math.isfinite(bound) for bound in bounds):
            raise ValueError(
                f"a log axis needs finite bounds, got low={self.low!r} high={self.high!r} "
                f"left={self.left!r} width={self.width!r}; an infinite span places every "
                f"value at the same pixel exactly as a zero width does, and a nan bound "
                f"passes every comparison below because nan <= x is always False"
            )
        if self.low <= 0:
            raise ValueError(f"a log axis needs a positive low bound, got {self.low!r}")
        if self.high <= self.low:
            raise ValueError(
                f"a log axis needs high above low, got low={self.low!r} high={self.high!r}"
            )
        if self.width <= 0:
            raise ValueError(
                f"a log axis needs a positive width, got {self.width!r}; a zero-width axis "
                f"places every value at the same pixel rather than failing"
            )

    def x(self, value: float) -> float:
        """Where `value` sits. Refuses a non-positive value rather than clamping.

        A clamped zero would land on the axis origin and read as "the lowest
        measurement here", which is a different claim from "below the detection
        floor" and from "not sampled". The caller owns that decision because
        only the caller knows which of the three it has.
        """
        if value <= 0:
            raise ValueError(f"a log axis needs a positive value, got {value!r}")
        span = math.log10(self.high) - math.log10(self.low)
        offset = (math.log10(value) - math.log10(self.low)) / span
        return self.left + self.width * offset


def svg_figure(*, width: int, height: int, title: str, body: str) -> str:
    """One `<svg>`, labelled for assistive technology.

    `role="img"` with a `<title>` rather than a bare `<svg>`: a screen reader
    otherwise announces the individual shapes or nothing at all. The `<details>`
    table each chart sits above is the full text alternative; this is the
    one-sentence summary that says whether opening it is worth it.
    """
    return (
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'preserveAspectRatio="xMinYMin meet">'
        f"<title>{html.escape(title)}</title>{body}</svg>"
    )


def polyline(points: Sequence[tuple[float, float]], *, css_class: str) -> str:
    """A connected run of measurements.

    Refuses fewer than two points. One point is not a line, and a renderer that
    quietly emitted a zero-length polyline for a gene measured at a single
    stage would draw a trend from one measurement. Live case, not
    hypothetical: measured 2026-08-20 on the committed corpus, FOXH1 and
    TFAP2B are each placed in heart at 1 of the 19 stages this dataset samples
    heart at, and 18 of the 581 (gene, organ) series with anything placed at
    all have exactly one placed point.
    """
    if len(points) < 2:
        raise ValueError(f"a polyline needs at least two points, got {len(points)}")
    drawn = " ".join(f"{coordinate(x)},{coordinate(y)}" for x, y in points)
    return f'<polyline class="{html.escape(css_class)}" points="{drawn}"/>'


def marker(x: float, y: float, *, css_class: str) -> str:
    """One measured point."""
    return (
        f'<circle class="{html.escape(css_class)}" '
        f'cx="{coordinate(x)}" cy="{coordinate(y)}" r="2.5"/>'
    )


Shape = Literal["circle", "square", "triangle_up", "triangle_down", "star", "ring"]

_GLYPH_CLASS: dict[Shape, str] = {
    "circle": "lit-mark",
    "square": "lit-mark",
    "triangle_up": "lit-mark",
    "triangle_down": "lit-mark",
    "star": "lit-mark",
    "ring": "lit-mark lit-mark-open",
}


def glyph(x: float, y: float, *, shape: Shape, title: str) -> str:
    """One titled mark for a single item, centred on `(x, y)` in a 10x10 box.

    Shape is the only channel: colour is the stylesheet's, and the open
    `ring` is the one hollow form, so the six read apart without colour. A
    `<title>` child carries the item's own text, so a hover or a screen reader
    names what the mark stands for.

    An unknown shape is refused, never drawn as a default. A vocabulary member
    added later must fail loudly rather than publish a mark the reader cannot
    tell from the others.

    The shape set deliberately has no diamond; see the module docstring.
    """
    if shape not in _GLYPH_CLASS:
        raise ValueError(f"unknown glyph shape {shape!r}; the shapes are {sorted(_GLYPH_CLASS)}")
    css = _GLYPH_CLASS[shape]
    label = f"<title>{html.escape(title)}</title>"
    cx, cy = coordinate(x), coordinate(y)
    if shape in ("circle", "ring"):
        r = "3.5" if shape == "ring" else "4"
        return f'<circle class="{css}" cx="{cx}" cy="{cy}" r="{r}">{label}</circle>'
    if shape == "square":
        return (
            f'<rect class="{css}" x="{coordinate(x - 4)}" y="{coordinate(y - 4)}" '
            f'width="8" height="8">{label}</rect>'
        )
    if shape in ("triangle_up", "triangle_down"):
        apex = -4.5 if shape == "triangle_up" else 4.5
        base = -apex
        points = [(x, y + apex), (x + 4.5, y + base), (x - 4.5, y + base)]
        return f'<polygon class="{css}" points="{_points(points)}">{label}</polygon>'
    if shape == "star":
        # Five outer and five inner vertices, alternating, first point up.
        vertices: list[tuple[float, float]] = []
        for i in range(10):
            angle = -math.pi / 2 + i * math.pi / 5
            reach = 5.0 if i % 2 == 0 else 2.2
            vertices.append((x + reach * math.cos(angle), y + reach * math.sin(angle)))
        return f'<polygon class="{css}" points="{_points(vertices)}">{label}</polygon>'
    raise ValueError(f"glyph shape {shape!r} has a class but no drawing")


def _points(points: Sequence[tuple[float, float]]) -> str:
    return " ".join(f"{coordinate(px)},{coordinate(py)}" for px, py in points)
