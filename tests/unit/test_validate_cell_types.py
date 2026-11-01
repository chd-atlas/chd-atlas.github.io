# tests/unit/test_validate_cell_types.py
"""CTP001-CTP011: a cell-type row checked against its dataset record.

Each case breaks exactly one rule on top of a clean base, so a failure names
the rule that stopped isolating. The base carries a row for every cell type and
region the dataset declares, because CTP008 reports a declaration no row backs.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import pytest

from chd_atlas.models.dataset import Dataset
from chd_atlas.validate.cell_types import validate_cell_type_profiles

_HEADER = (
    "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
    "n_donors\tn_cells\tfraction_expressing\tsource\n"
)


def _dataset() -> Dataset:
    return Dataset(
        id="hoc",
        archive="ucsc_cell_browser",
        technology="scrnaseq",
        design="cell_type_profile",
        tissue="Heart",
        developmental_stage="adult",
        organism="NCBITaxon:9606",
        n_samples=10,
        licence="CC BY 4.0",
        cell_types=[
            {"token": "Cardiomyocyte", "label": "Cardiomyocyte", "compartment": "myocardial"},
            {"token": "Fibroblast", "label": "Fibroblast", "compartment": "stromal"},
        ],
        regions=["LV", "RV"],
        detection_floor=0.1,
        floor_source="source methods, section 2",
        quantile_estimator="linear",
    )


def _row(
    *,
    dataset: str = "hoc",
    gene: str = "HGNC:1100",
    cell_type: str = "Cardiomyocyte",
    region: str = "all",
    median: str = "2.0",
    q25: str = "1.0",
    q75: str = "3.0",
    n_donors: str = "4",
    n_cells: str = "900",
) -> str:
    return (
        f"{dataset}\t{gene}\t{cell_type}\t{region}\t{median}\tcpm\t{q25}\t{q75}\t"
        f"{n_donors}\t{n_cells}\t0.4\ts\n"
    )


def _base(skip: tuple[str, ...] = (), **first: str) -> list[str]:
    """A row for every declared (cell type, region), plus `all`; `first`
    overrides the first row (Cardiomyocyte, `all`), and a cell type or region
    named in `skip` gets no row."""
    rows = []
    for cell_type in ("Cardiomyocyte", "Fibroblast"):
        for region in ("all", "LV", "RV"):
            if cell_type in skip or region in skip:
                continue
            fields = {"cell_type": cell_type, "region": region}
            if not rows:
                fields.update(first)
            rows.append(_row(**fields))
    return rows


def _run(tmp_path: Path, shards: Mapping[str, str] | None, genes: set[str] | None) -> list[str]:
    """Codes reported over `shards` (name -> file text); None writes no directory."""
    if shards is not None:
        directory = tmp_path / "mirrors" / "cell_type_profiles"
        directory.mkdir(parents=True)
        for name, text in shards.items():
            (directory / f"{name}.tsv").write_text(text, encoding="utf-8")
    issues = validate_cell_type_profiles(tmp_path, [_dataset()], genes)
    return [issue.code for issue in issues]


def _shard(rows: list[str]) -> str:
    return _HEADER + "".join(rows)


_GENES = {"HGNC:1100"}


@pytest.mark.parametrize(
    ("rows", "genes", "expected"),
    [
        pytest.param(_base(), _GENES, [], id="clean"),
        pytest.param([*_base(), _row(dataset="other")], _GENES, ["CTP001"], id="CTP001"),
        pytest.param([*_base(), _row(cell_type="Mast")], _GENES, ["CTP002"], id="CTP002"),
        pytest.param([*_base(), _row(region="RA")], _GENES, ["CTP003"], id="CTP003"),
        pytest.param(
            [r for r in _base() if "\tall\t" not in r or "Fibroblast" in r],
            _GENES,
            ["CTP004"],
            id="CTP004",
        ),
        pytest.param(_base(median="5.0"), _GENES, ["CTP005"], id="CTP005"),
        pytest.param(_base(n_donors="11"), _GENES, ["CTP006"], id="CTP006"),
        pytest.param(_base(), {"HGNC:9"}, ["CTP007"], id="CTP007"),
        pytest.param(_base(q25=""), _GENES, ["CTP005"], id="one-quartile-null"),
        pytest.param(_base(), None, [], id="registry-unloaded-skips-CTP007"),
        # A declared token no row carries: one issue per token.
        pytest.param(_base(skip=("RV",)), _GENES, ["CTP008"], id="CTP008-region"),
        pytest.param(_base(skip=("Fibroblast",)), _GENES, ["CTP008"], id="CTP008-cell-type"),
        pytest.param(
            _base(skip=("Fibroblast", "RV")), _GENES, ["CTP008", "CTP008"], id="CTP008-two"
        ),
        # One issue for the dataset, not one per declared token (four here).
        pytest.param([], _GENES, ["CTP008"], id="CTP008-no-rows"),
        # Quartiles from two hearts, and none from four.
        pytest.param(_base(n_donors="2", n_cells="20"), _GENES, ["CTP009"], id="CTP009-two"),
        pytest.param(_base(q25="", q75=""), _GENES, ["CTP009"], id="CTP009-withheld"),
        pytest.param(_base(n_donors="3", n_cells="30"), _GENES, [], id="CTP009-boundary"),
        # Four hearts of at least 10 cells each pool at least 40.
        pytest.param(_base(n_cells="39"), _GENES, ["CTP010"], id="CTP010"),
        pytest.param(_base(n_cells="40"), _GENES, [], id="CTP010-boundary"),
    ],
)
def test_each_cell_type_check_isolates_its_own_error(
    tmp_path: Path, rows: list[str], genes: set[str] | None, expected: list[str]
) -> None:
    assert _run(tmp_path, {"hoc": _shard(rows)}, genes) == expected


def test_a_key_repeated_across_shards_is_ctp011_and_within_one_is_left_to_tbl007(
    tmp_path: Path,
) -> None:
    """TBL007 refuses a repeated sort key within one shard, on the same
    (dataset, gene, cell_type, region) key, so CTP011 covers repeats across
    shards only. Measured before CTP011 (review, 2026-10-28): a rerun shard
    beside the first validated clean and the build published whichever row
    it read last."""
    rerun = _shard([_row(median="2.5")])
    assert _run(tmp_path, {"hoc": _shard(_base()), "hoc_rerun": rerun}, _GENES) == ["CTP011"]
    within = tmp_path / "within"
    within.mkdir()
    assert _run(within, {"hoc": _shard([*_base(), _row(median="2.5")])}, _GENES) == []


@pytest.mark.parametrize(
    "shards",
    [
        pytest.param(None, id="no-directory"),
        pytest.param({"other": _shard([_row(dataset="other")])}, id="rows-for-another-dataset"),
    ],
)
def test_a_declared_dataset_with_no_rows_is_reported_once(
    tmp_path: Path, shards: Mapping[str, str] | None
) -> None:
    """A missing `mirrors/cell_type_profiles/` is the case CTP008 exists for,
    not a reason to skip it."""
    codes = _run(tmp_path, shards, _GENES)
    assert codes.count("CTP008") == 1


@pytest.mark.parametrize(
    "broken",
    [
        pytest.param("", id="unreadable"),
        pytest.param(_HEADER.replace("\tn_cells", ""), id="missing-column"),
    ],
)
def test_ctp008_is_skipped_when_a_shard_could_not_be_read(tmp_path: Path, broken: str) -> None:
    """TBL000/TBL001 report the shard; the dataset's rows may be in it, so a
    "no data" error for every declaration would bury that cause."""
    assert _run(tmp_path, {"broken": broken, "hoc": _shard([])}, _GENES) == []
