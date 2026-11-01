# tests/unit/test_build_cell_type_chart.py
"""The gene page's cell-type block: every figure on it is read off the SVG element it labels."""

from __future__ import annotations

import re

import pytest

from chd_atlas.build.cell_type_chart import cell_type_block, format_median
from chd_atlas.build.cell_types import (
    CellTypeDataset,
    CellTypeEntry,
    CellTypeSpecificity,
    CellTypeValues,
    RegionEntry,
)
from chd_atlas.build.render import STYLESHEET
from chd_atlas.models.dataset import Dataset
from tests.test_built_site_is_consumable import _FIXED_PRECISION, _GEOMETRY

GENE = "GATA4"
REGIONS = ("LA", "RA", "LV", "RV", "IVS")
AXIS_MIN_X = "190.0"  # _LEFT: the x of 0.1 CPM

_DECLARED = (
    ("aCM", "Atrial cardiomyocytes", "Cardiomyocyte"),
    ("vCM", "Ventricular cardiomyocytes", "Cardiomyocyte"),
    ("Fibro", "Fibroblasts", "Mesenchymal"),
    ("SMC", "Smooth muscle cells", "Mesenchymal"),
    ("WBC", "White blood cells", "Blood"),
)


def dataset(
    declared: tuple[tuple[str, str, str], ...] = _DECLARED,
    regions: tuple[str, ...] = REGIONS,
) -> Dataset:
    return Dataset(
        id="hoc",
        archive="ucsc_cell_browser",
        technology="scrnaseq",
        design="cell_type_profile",
        tissue="Heart",
        developmental_stage="fetal",
        organism="NCBITaxon:9606",
        n_samples=8,
        licence="CC BY 4.0",
        cell_types=[{"token": t, "label": n, "compartment": c} for t, n, c in declared],
        regions=list(regions),
        detection_floor=1.0,
        floor_source="source methods, section 2",
        quantile_estimator="linear",
        # The values the block's fixed copy states; `cell_type_block` refuses others.
        publication="PMID:38480880",
    )


def values(
    median: float,
    *,
    fraction: float,
    donors: int,
    cells: int,
    q25: float | None = None,
    q75: float | None = None,
    floor: float = 1.0,
) -> CellTypeValues:
    return CellTypeValues(
        median=median,
        q25=q25,
        q75=q75,
        n_donors=donors,
        n_cells=cells,
        fraction_expressing=fraction,
        below_floor=median < floor,
    )


def entry(token: str, pooled: CellTypeValues, **by_region: CellTypeValues) -> CellTypeEntry:
    label = {t: n for t, n, _ in _DECLARED}[token]
    return CellTypeEntry(
        cell_type=token,
        label=label,
        pooled=pooled,
        regions=[RegionEntry(region=r, **by_region[r]) for r in REGIONS if r in by_region],
    )


def specificity(highest: str | None = "aCM") -> CellTypeSpecificity:
    return CellTypeSpecificity(
        tau=0.8765,
        scale="log2(x+1)",
        method="m",
        cell_types=("aCM", "vCM", "Fibro", "SMC"),
        n_cell_types=4,
        highest_in=highest,
        medians={},
    )


def payload(
    *, spec: CellTypeSpecificity | None = None, extra: list[CellTypeEntry] | None = None
) -> CellTypeDataset:
    """Distinct figures per cell type, so no two rows can be mistaken for each other."""
    acm = values(1234.6, fraction=0.62, donors=7, cells=900, q25=800.0, q75=2000.0)
    vcm = values(56.7, fraction=0.38, donors=6, cells=1100, q25=30.0, q75=90.0)
    fibro = values(3.21, fraction=0.12, donors=8, cells=1300, q25=None, q75=None)
    return CellTypeDataset(
        dataset="hoc",
        unit="cpm",
        detection_floor=1.0,
        n_donors=8,
        specificity=spec,
        cell_types=[
            entry("aCM", acm, LA=acm, LV=vcm),
            entry("vCM", vcm),
            entry("Fibro", fibro, RA=fibro),
            *(extra or []),
        ],
    )


def element(block: str, tag: str, token: str, region: str | None = None) -> str:
    """The one `<tag ...>` element carrying `data-cell-type=token` (and `region`)."""
    attr = f'data-cell-type="{token}"' + (f' data-region="{region}"' if region else "")
    found = re.findall(rf"<{tag}\b[^>]*{re.escape(attr)}[^>]*>(?:.*?</{tag}>)?", block, re.S)
    chart_only = [f for f in found if region or "data-region" not in f]
    assert len(chart_only) == 1, (tag, token, region, found)
    return chart_only[0]


def attr(markup: str, name: str) -> float:
    return float(re.search(rf'\b{name}="([-0-9.]+)"', markup).group(1))  # type: ignore[union-attr]


def test_each_row_states_median_share_and_hearts_beside_its_dot() -> None:
    block = cell_type_block(payload(), dataset(), GENE)

    assert ">1,235 CPM · 62% · 7 hearts<" in element(block, "text", "aCM")
    assert ">57 CPM · 38% · 6 hearts<" in element(block, "text", "vCM")
    assert ">3.2 CPM · 12% · 8 hearts<" in element(block, "text", "Fibro")
    # The dot's own title carries the full numbers, and the row's position follows the median.
    title = re.search(r"<title>(.*?)</title>", element(block, "circle", "aCM"), re.S)
    assert title is not None
    assert title.group(1) == (
        "Atrial cardiomyocytes: median 1,235 CPM, quartiles 800 to 2,000, 62% of 900 cells, "
        "7 hearts"
    )
    assert attr(element(block, "circle", "aCM"), "cx") > attr(element(block, "circle", "vCM"), "cx")
    assert attr(element(block, "circle", "vCM"), "cx") > attr(
        element(block, "circle", "Fibro"), "cx"
    )
    # Area follows the share of cells: the 62% dot is larger than the 12% one.
    assert attr(element(block, "circle", "aCM"), "r") > attr(element(block, "circle", "Fibro"), "r")
    # A whisker spans the quartiles where there are any, and is absent where there are none.
    whisker = element(block, "line", "aCM")
    assert attr(whisker, "x1") < attr(element(block, "circle", "aCM"), "cx") < attr(whisker, "x2")
    assert 'data-cell-type="Fibro"' not in "".join(re.findall(r"<line\b[^>]*>", block))


def test_every_declared_cell_type_has_a_row_and_one_without_data_says_so() -> None:
    block = cell_type_block(payload(), dataset(), GENE)
    labels = re.findall(r'<text[^>]*class="ct-label"[^>]*>(.*?)</text>', block)

    # Declared order, each label once in the chart and once in the region grid.
    assert labels == [n for _, n, _ in _DECLARED] * 2
    groups = re.findall(r'<text[^>]*class="ct-group"[^>]*>(.*?)</text>', block)
    assert groups == ["Cardiomyocyte", "Mesenchymal", "Blood"] * 2
    chart = block.split("<details>")[0]
    assert re.findall(r'class="ct-muted"[^>]*>(.*?)</text>', chart) == ["no data", "no data"]
    assert 'data-cell-type="SMC"' not in chart


def test_a_zero_median_draws_at_the_axis_minimum_and_reads_0_cpm() -> None:
    zero = values(0.0, fraction=0.03, donors=5, cells=700, q25=0.0, q75=0.0)
    block = cell_type_block(payload(extra=[entry("SMC", zero)]), dataset(), GENE)

    assert ">0 CPM · 3% · 5 hearts<" in element(block, "text", "SMC")
    dot = element(block, "circle", "SMC")
    assert f'cx="{AXIS_MIN_X}"' in dot
    assert "median 0 CPM" in dot


def test_a_cell_type_with_no_regions_shows_a_dash_in_every_region_column() -> None:
    block = cell_type_block(payload(), dataset(), GENE)
    grid = block.split("<details>")[1]

    for region in REGIONS:
        assert re.search(r">–</text>", element(grid, "text", "vCM", region))
        assert "<circle" not in element(grid, "text", "vCM", region)
    # Where a region has data it is a dot, not a dash, and the others stay dashes.
    assert "<circle" in element(grid, "circle", "aCM", "LA")
    assert "<circle" in element(grid, "circle", "aCM", "LV")
    assert ">–</text>" in element(grid, "text", "aCM", "RA")
    # Columns follow the dataset's regions, in order.
    assert "<summary>By region</summary>" in block
    assert re.findall(r'class="ct-small ct-mid"[^>]*>(LA|RA|LV|RV|IVS)<', grid) == list(REGIONS)


def test_a_region_dot_is_a_ring_below_the_floor_and_steps_its_opacity_by_decade() -> None:
    low = values(0.4, fraction=0.21, donors=5, cells=300)
    high = values(5400.0, fraction=0.55, donors=6, cells=400)
    mid = values(42.0, fraction=0.33, donors=6, cells=500)
    pooled = values(10.0, fraction=0.4, donors=6, cells=500)
    block = cell_type_block(
        CellTypeDataset(
            dataset="hoc",
            unit="cpm",
            detection_floor=1.0,
            n_donors=8,
            specificity=None,
            cell_types=[entry("aCM", pooled, LA=low, RA=mid, LV=high)],
        ),
        dataset(),
        GENE,
    )
    grid = block.split("<details>")[1]

    assert 'class="ct-ring"' in element(grid, "circle", "aCM", "LA")
    assert 'fill-opacity="0.55"' in element(grid, "circle", "aCM", "RA")
    assert 'fill-opacity="1.0"' in element(grid, "circle", "aCM", "LV")


def test_the_tau_line_is_omitted_when_specificity_is_null() -> None:
    assert "τ =" not in cell_type_block(payload(spec=None), dataset(), GENE)


def test_the_tau_line_names_the_highest_cell_type_or_says_there_is_none() -> None:
    named = cell_type_block(payload(spec=specificity("aCM")), dataset(), GENE)
    tied = cell_type_block(payload(spec=specificity(None)), dataset(), GENE)

    assert (
        '<p class="ct-tau">τ = 0.88 (4 cell types, log2(x+1)), highest in '
        "Atrial cardiomyocytes.</p>"
    ) in named
    assert (
        '<p class="ct-tau">τ = 0.88 (4 cell types, log2(x+1)), with no single highest.</p>'
    ) in tied
    assert "highest in" not in tied


def test_the_heading_and_caption_are_the_pinned_copy() -> None:
    block = cell_type_block(payload(), dataset(), GENE)

    assert "<h3>Cell types (Farah 2024)</h3>" in block
    assert ("Log axis, 0.1 to 10,000 CPM; a dot's size is its share of cells.</p>") in block
    assert "<script" not in block


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("publication", "PMID:8988165"),
        ("publication", None),
        ("n_samples", 9),
        ("detection_floor", 0.5),
    ],
)
def test_a_dataset_the_fixed_copy_does_not_describe_is_refused(field: str, value: object) -> None:
    """The heading says "(Farah 2024)", both banners "8 hearts" and the fold
    "1 CPM", whatever the record says. A second dataset needs that copy
    generalised before it can render, so the block refuses it, naming the
    field, rather than print one study's facts above another's figures."""
    other = dataset().model_copy(update={field: value})
    with pytest.raises(ValueError, match=field):
        cell_type_block(payload(), other, GENE)


def test_the_how_to_read_fold_is_at_most_80_words_and_links_methods() -> None:
    block = cell_type_block(payload(), dataset(), GENE)
    fold = block.split("<summary>How to read this chart</summary>")[1].split("</details>")[0]
    text = re.sub(r"<[^>]+>", "", fold)

    assert text == (
        "The dot is the median across hearts; the line spans the middle half of hearts (shown "
        "from three hearts). A heart counts only with 10 or more cells of that type. The shaded "
        'band is below this atlas\'s 1 CPM detection floor. "No data" means no heart reached 10 '
        "cells; it is not a low value. Values pool 9 to 15 weeks. How to read cell-type figures."
    )
    assert len(text.split()) <= 80
    assert '<a href="../methods.html#expression">How to read cell-type figures</a>.' in fold


def test_the_shaded_band_ends_at_the_datasets_own_floor() -> None:
    block = cell_type_block(payload(), dataset(), GENE)

    band = re.search(r'<rect x="([0-9.]+)"[^>]*width="([0-9.]+)"[^>]*class="ct-below"', block)
    assert band is not None
    # 190 + (log10(1) - log10(0.1)) / 5 * 430 = 276
    assert float(band.group(1)) + float(band.group(2)) == pytest.approx(276.0)


@pytest.mark.parametrize(
    ("median", "text"),
    [
        (0.0, "0"),
        (0.004, "0.00"),
        (0.456, "0.46"),
        (1.0, "1.0"),
        (9.96, "10.0"),
        (10.0, "10"),
        (99.4, "99"),
        (100.0, "100"),
        (1234.5, "1,235"),
        (10000.0, "10,000"),
    ],
)
def test_the_median_format_steps_at_1_10_and_100(median: float, text: str) -> None:
    assert format_median(median) == text


def test_the_stylesheet_sets_no_fill_opacity_on_a_class_the_region_grid_draws() -> None:
    """The grid's decade shading is each dot's `fill-opacity` attribute, and a CSS
    rule setting that property on the same element wins over the attribute.
    `.ct-dot { fill-opacity: 0.85 }` drew every grid dot identically until the
    final review (2026-10-28); the chart's own dots carry their 0.85 as an
    attribute instead."""
    block = cell_type_block(payload(), dataset(), GENE)
    grid = block.split("<details>")[1]
    drawn = {name for value in re.findall(r'class="([^"]*)"', grid) for name in value.split()}
    assert {"ct-dot", "ct-ring"} <= drawn
    rules = re.findall(r"([^{}]+)\{([^{}]*)\}", re.sub(r"/\*.*?\*/", "", STYLESHEET, flags=re.S))
    styled = [
        (selector.strip(), body)
        for selector, body in rules
        if set(re.findall(r"\.([\w-]+)", selector)) & drawn or re.search(r"\bcircle\b", selector)
    ]
    assert any(selector == ".ct-dot" for selector, _ in styled), "no rule styles the grid's dots"
    assert [selector for selector, body in styled if "fill-opacity" in body] == []
    assert 'fill-opacity="0.85"' in element(block, "circle", "aCM")


def test_every_coordinate_the_block_draws_has_coordinate_s_fixed_shape() -> None:
    """The site-wide scan (`test_every_geometry_the_pages_publish_...`) reads the
    real build, which has no cell-type mirror yet, so it could not see this
    block's `r="6.97"`. The same check, over a fixture with whiskers, a dot in
    every grid state and the legend."""
    low = values(0.4, fraction=0.21, donors=5, cells=300)
    block = cell_type_block(
        payload(extra=[entry("SMC", low, LA=low, RA=low)], spec=specificity()), dataset(), GENE
    )
    tokens = [token for value in _GEOMETRY.findall(block) for token in value.split()]
    assert len(tokens) > 100, len(tokens)
    assert [token for token in tokens if not _FIXED_PRECISION.match(token)] == []


def test_one_heart_reads_singular_in_the_column_and_the_tooltip() -> None:
    one = values(12.0, fraction=0.5, donors=1, cells=40)
    block = cell_type_block(payload(extra=[entry("SMC", one, RV=one)]), dataset(), GENE)

    assert ">12 CPM · 50% · 1 heart<" in element(block, "text", "SMC")
    pooled = re.search(r"<title>(.*?)</title>", element(block, "circle", "SMC"))
    regional = re.search(r"<title>(.*?)</title>", element(block, "circle", "SMC", "RV"))
    assert pooled is not None and pooled.group(1).endswith(", 50% of 40 cells, 1 heart")
    assert regional is not None and regional.group(1).endswith(", 1 heart")
    assert "1 hearts" not in block


def test_curated_strings_are_escaped() -> None:
    """Labels, compartments and tokens are free text in a curated record; a
    `<` in one would open a tag on the page. A region is one of five tokens
    since review (2026-10-28) and is escaped all the same: the record's
    regions here bypass the model to show it."""
    declared = (("T<1", "Cells <10 um", "Blood & immune"),)
    odd = values(5.0, fraction=0.3, donors=4, cells=200)
    block = cell_type_block(
        CellTypeDataset(
            dataset="hoc",
            unit="cpm",
            detection_floor=1.0,
            n_donors=8,
            specificity=CellTypeSpecificity(
                tau=0.5,
                scale="log2(x+1)",
                method="m",
                cell_types=("T<1",),
                n_cell_types=2,
                highest_in="T<1",
                medians={},
            ),
            cell_types=[
                CellTypeEntry(
                    cell_type="T<1",
                    label="Cells <10 um",
                    pooled=odd,
                    regions=[RegionEntry(region="L<A", **odd)],
                )
            ],
        ),
        dataset(declared).model_copy(update={"regions": ("L<A",)}),
        "A<B",
    )

    assert "<10" not in block and "<1" not in block and "<A" not in block and "<B" not in block
    assert "& " not in block
    assert ">Cells &lt;10 um</text>" in block
    assert ">Blood &amp; immune</text>" in block
    assert 'data-cell-type="T&lt;1" data-region="L&lt;A"' in block
    assert ">L&lt;A</text>" in block
    assert "highest in Cells &lt;10 um." in block
    assert 'aria-label="A&lt;B expression by cell type"' in block


def test_the_tau_line_sits_below_the_chart_and_the_figures_name_the_gene() -> None:
    block = cell_type_block(payload(spec=specificity()), dataset(), GENE)
    chart_end = block.index("</svg>")

    assert chart_end < block.index('<p class="ct-tau">') < block.index("<details>")
    assert 'aria-label="GATA4 expression by cell type"' in block
    assert 'aria-label="GATA4 expression by cell type and region"' in block
    assert "hoc expression" not in block


def test_the_region_grid_s_key_defines_its_dash() -> None:
    grid = cell_type_block(payload(), dataset(), GENE).split("<details>")[1]
    key = grid.split('class="ct-small">1,000+ CPM</text>')[1]

    assert re.findall(r"<text[^>]*>([^<]*)</text>", key) == [
        "–",
        "no heart had 10 or more cells there",
    ]
