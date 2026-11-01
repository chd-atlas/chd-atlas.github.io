# tests/unit/test_build_cell_types.py
"""A gene's cell-type expression, assembled from the `cell_type_profiles` mirror."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from chd_atlas.build.cell_types import (
    CELL_TYPE_TAU_METHOD,
    CellTypeExpression,
    cell_type_census,
    gene_cell_type_expression,
)
from chd_atlas.models.dataset import Dataset

_HEADER = (
    "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
    "n_donors\tn_cells\tfraction_expressing\tsource\n"
)
_GENE = "HGNC:1100"
_TOKENS = ("Cardiomyocyte", "Fibroblast", "Endothelial")


def _dataset(
    *,
    dataset_id: str = "hoc",
    tokens: tuple[str, ...] = _TOKENS,
    regions: tuple[str, ...] = ("LV", "RV"),
    floor: float = 0.5,
) -> Dataset:
    return Dataset(
        id=dataset_id,
        archive="ucsc_cell_browser",
        technology="scrnaseq",
        design="cell_type_profile",
        tissue="Heart",
        developmental_stage="adult",
        organism="NCBITaxon:9606",
        n_samples=10,
        licence="CC BY 4.0",
        cell_types=[{"token": t, "label": f"{t} label", "compartment": "x"} for t in tokens],
        regions=list(regions),
        detection_floor=floor,
        floor_source="source methods, section 2",
        quantile_estimator="linear",
    )


def _row(cell_type: str, region: str, median: float, *, gene: str = _GENE, ds: str = "hoc") -> str:
    return f"{ds}\t{gene}\t{cell_type}\t{region}\t{median}\tcpm\t0.1\t9.0\t4\t900\t0.4\ts\n"


def _write(root: Path, rows: list[str], name: str = "hoc") -> None:
    shard = root / "mirrors" / "cell_type_profiles"
    shard.mkdir(parents=True, exist_ok=True)
    (shard / f"{name}.tsv").write_text(_HEADER + "".join(rows), encoding="utf-8")


def _entry(root: Path, rows: list[str], dataset: Dataset | None = None) -> CellTypeExpression:
    _write(root, rows)
    return gene_cell_type_expression(root, [dataset or _dataset()])[_GENE]


def test_cell_types_and_regions_follow_the_declared_order(tmp_path: Path) -> None:
    script = """
import json, sys
from pathlib import Path

import pytest
from chd_atlas.build.cell_types import gene_cell_type_expression
from chd_atlas.models.dataset import Dataset

tokens = ["zebra", "apple", "mango", "kiwi"]
regions = ["RV", "IVS", "LA"]
dataset = Dataset(
    id="hoc", archive="ucsc_cell_browser", technology="scrnaseq",
    design="cell_type_profile", tissue="Heart", developmental_stage="adult",
    organism="NCBITaxon:9606", n_samples=10, licence="CC BY 4.0",
    cell_types=[{"token": t, "label": t, "compartment": "x"} for t in tokens],
    regions=regions, detection_floor=0.5, floor_source="s",
    quantile_estimator="linear",
)
entry = gene_cell_type_expression(Path(sys.argv[1]), [dataset])["HGNC:1100"]
ds = entry["datasets"][0]
print(json.dumps({
    "cell_types": [c["cell_type"] for c in ds["cell_types"]],
    "regions": [[r["region"] for r in c["regions"]] for c in ds["cell_types"]],
}))
"""
    rows = []
    # Rows in sorted-file order, which matches neither declared order.
    for token in sorted(["zebra", "apple", "mango", "kiwi"]):
        for region in ("IVS", "LA", "RV", "all"):
            rows.append(_row(token, region, 5.0))
    _write(tmp_path, rows)
    outputs = set()
    for seed in ("0", "1", "2", "3", "4"):
        env = {**os.environ, "PYTHONHASHSEED": seed}
        done = subprocess.run(
            [sys.executable, "-c", script, str(tmp_path)],
            env=env,
            capture_output=True,
            text=True,
            check=True,
        )
        outputs.add(done.stdout.strip())
    assert len(outputs) == 1
    result = json.loads(outputs.pop())
    assert result["cell_types"] == ["zebra", "apple", "mango", "kiwi"]
    assert result["regions"] == [["RV", "IVS", "LA"]] * 4


def test_an_absent_region_is_not_below_the_floor(tmp_path: Path) -> None:
    rows = [
        _row("Cardiomyocyte", "all", 5.0),
        _row("Cardiomyocyte", "LV", 0.4),
    ]
    dataset = _entry(tmp_path, rows)["datasets"][0]
    (cell_type,) = dataset["cell_types"]
    assert [r["region"] for r in cell_type["regions"]] == ["LV"]
    assert cell_type["regions"][0]["below_floor"] is True
    assert cell_type["regions"][0]["median"] == 0.4
    assert cell_type["pooled"]["below_floor"] is False
    assert dataset["detection_floor"] == 0.5
    assert dataset["n_donors"] == 10
    assert dataset["unit"] == "cpm"


def test_a_cell_type_without_a_pooled_row_is_omitted(tmp_path: Path) -> None:
    rows = [_row("Cardiomyocyte", "all", 5.0), _row("Fibroblast", "LV", 7.0)]
    dataset = _entry(tmp_path, rows)["datasets"][0]
    assert [c["cell_type"] for c in dataset["cell_types"]] == ["Cardiomyocyte"]
    assert dataset["cell_types"][0]["label"] == "Cardiomyocyte label"


def test_a_median_of_exactly_the_floor_is_not_below_it(tmp_path: Path) -> None:
    rows = [_row("Cardiomyocyte", "all", 0.5), _row("Cardiomyocyte", "LV", 0.5)]
    (cell_type,) = _entry(tmp_path, rows)["datasets"][0]["cell_types"]
    assert cell_type["pooled"]["below_floor"] is False
    assert cell_type["regions"][0]["below_floor"] is False


def test_tau_uses_pooled_medians_and_the_bulk_method(tmp_path: Path) -> None:
    rows = [
        _row("Cardiomyocyte", "all", 100.0),
        _row("Fibroblast", "all", 10.0),
        _row("Endothelial", "all", 1.0),
        # A regional row must not enter tau.
        _row("Endothelial", "LV", 900.0),
    ]
    # "Pericyte" is declared and has no rows, so tau is over the other three.
    declared = _dataset(tokens=("Cardiomyocyte", "Pericyte", "Fibroblast", "Endothelial"))
    specificity = _entry(tmp_path, rows, declared)["datasets"][0]["specificity"]
    assert specificity is not None
    assert specificity["tau"] == 0.6651179051469356
    assert specificity["highest_in"] == "Cardiomyocyte"
    assert specificity["n_cell_types"] == 3
    # The dataset's declared order, present types only; sorted would put
    # Endothelial second. `cell_types` beside it is in the same order.
    assert specificity["cell_types"] == ("Cardiomyocyte", "Fibroblast", "Endothelial")
    assert specificity["medians"] == {
        "Cardiomyocyte": 100.0,
        "Fibroblast": 10.0,
        "Endothelial": 1.0,
    }
    assert specificity["scale"] == "log2(x+1)"
    assert specificity["method"] == CELL_TYPE_TAU_METHOD


def test_tau_is_null_when_every_cell_type_is_below_the_floor(tmp_path: Path) -> None:
    rows = [_row("Cardiomyocyte", "all", 0.4), _row("Fibroblast", "all", 0.3)]
    dataset = _entry(tmp_path, rows)["datasets"][0]
    assert dataset["specificity"] is None
    assert all(c["pooled"]["below_floor"] for c in dataset["cell_types"])


def test_the_tau_method_names_cell_types_not_organs() -> None:
    assert "organ" not in CELL_TYPE_TAU_METHOD
    assert CELL_TYPE_TAU_METHOD == (
        "tau (Yanai et al. 2005): the sum over cell types of (1 - x_i/x_max) "
        "divided by (n - 1), where n is the number of cell types and "
        "x = log2(median+1); a negative median is clamped to 0 before the "
        "transform, and every cell type's raw median is used even below the "
        "dataset's detection floor"
    )


def test_the_census_counts_published_genes_only(tmp_path: Path) -> None:
    _write(
        tmp_path,
        [
            _row("Cardiomyocyte", "all", 5.0),
            _row("Cardiomyocyte", "all", 5.0, gene="HGNC:2200"),
        ],
    )
    expression = gene_cell_type_expression(tmp_path, [_dataset()])
    assert set(expression) == {_GENE, "HGNC:2200"}
    assert cell_type_census(expression, {_GENE, "HGNC:9999"}) == {"genes": 1, "datasets": 1}


def test_a_bypassed_gate_fails_rather_than_publishes(tmp_path: Path) -> None:
    """Two inputs only a bypassed validator lets through (CLAUDE.md section 9).

    A key repeated across shards (CTP011): the build kept whichever row it read
    last, and measured on the review's fixture (2026-10-28) published 230.0
    from a rerun shard over 220.67. A dataset with no detection floor, which
    the model refuses: an `assert` guarded it, and `-O` strips those.
    """
    _write(tmp_path, [_row("Cardiomyocyte", "all", 5.0)])
    _write(tmp_path, [_row("Cardiomyocyte", "all", 6.0)], name="hoc_rerun")
    with pytest.raises(ValueError, match=r"hoc.*HGNC:1100.*Cardiomyocyte.*all"):
        gene_cell_type_expression(tmp_path, [_dataset()])

    floorless = tmp_path / "floorless"
    _write(floorless, [_row("Cardiomyocyte", "all", 5.0)])
    dataset = _dataset().model_copy(update={"detection_floor": None})
    with pytest.raises(ValueError, match="detection_floor"):
        gene_cell_type_expression(floorless, [dataset])


def test_a_dataset_with_two_units_is_refused(tmp_path: Path) -> None:
    # Different genes, so no single gene's rows disagree: the check is per dataset.
    rows = [
        _row("Cardiomyocyte", "all", 5.0),
        _row("Cardiomyocyte", "all", 5.0, gene="HGNC:2200").replace("cpm", "tpm"),
    ]
    _write(tmp_path, rows)
    with pytest.raises(ValueError, match=r"hoc.*cpm.*tpm"):
        gene_cell_type_expression(tmp_path, [_dataset()])
