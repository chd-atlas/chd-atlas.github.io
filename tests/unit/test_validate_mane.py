"""MAN001-MAN004, each fired by its own corruption and by nothing else."""

from collections.abc import Callable
from pathlib import Path
from typing import Any

import polars as pl
import pytest

from chd_atlas.mane import ManeFrame, read_mane_frame
from chd_atlas.validate.mane import CODON_TABLE, validate_mane
from tests.unit.mane_fixtures import (
    GENE_MINUS,
    GENE_PLUS,
    genes_frame,
    make_gene,
    synthetic_rows,
    write_mane_mirrors,
)

Rows = tuple[list[dict[str, Any]], list[dict[str, Any]]]


def _frame(tmp_path: Path, rows: Rows) -> ManeFrame:
    write_mane_mirrors(tmp_path, *rows)
    frame, issues = read_mane_frame(tmp_path)
    assert frame is not None and issues == []
    return frame


def _flip_exon_base(rows: Rows) -> pl.DataFrame:
    # Flip a base inside exon 2 of the `+` gene, in its core rather than its flank.
    exon = rows[1][1]
    seq = exon["window_sequence"]
    at = 50 + 5
    exon["window_sequence"] = seq[:at] + ("A" if seq[at] != "A" else "C") + seq[at + 1 :]
    return genes_frame(rows[0])


def _gap_in_transcript_coordinates(rows: Rows) -> pl.DataFrame:
    # Exon 2 of the `-` gene starts one base late and ends one base late.
    exon = next(e for e in rows[1] if e["hgnc_id"] == GENE_MINUS and e["exon"] == 2)
    exon["transcript_start"] += 1
    exon["transcript_end"] += 1
    return genes_frame(rows[0])


def _truncated_window(rows: Rows) -> pl.DataFrame:
    exon = rows[1][0]
    exon["window_sequence"] = exon["window_sequence"][:-1]
    return genes_frame(rows[0])


def _exon_longer_in_genome_than_in_transcript(rows: Rows) -> pl.DataFrame:
    exon = rows[1][1]
    exon["genomic_end"] += 1
    return genes_frame(rows[0])


def _cds_not_a_multiple_of_three(rows: Rows) -> pl.DataFrame:
    rows[0][0]["cds_end"] += 1
    return genes_frame(rows[0])


def _cds_without_its_stop(rows: Rows) -> pl.DataFrame:
    rows[0][0]["cds_end"] -= 3
    return genes_frame(rows[0])


def _change_a_residue(rows: Rows) -> pl.DataFrame:
    t = rows[0][0]
    p = t["protein_sequence"]
    t["protein_sequence"] = p[:4] + ("W" if p[4] != "W" else "G") + p[5:]
    return genes_frame(rows[0])


def _version_mismatch(rows: Rows) -> pl.DataFrame:
    genes = genes_frame(rows[0])
    return genes.with_columns(
        pl.when(pl.col("hgnc_id") == GENE_PLUS)
        .then(pl.lit("NM_900001.3"))
        .otherwise(pl.col("mane_select"))
        .alias("mane_select")
    )


def _null_mane_select(rows: Rows) -> pl.DataFrame:
    return genes_frame(rows[0]).with_columns(
        pl.when(pl.col("hgnc_id") == GENE_MINUS)
        .then(None)
        .otherwise(pl.col("mane_select"))
        .alias("mane_select")
    )


def _registry_gene_without_a_row(rows: Rows) -> pl.DataFrame:
    return pl.concat(
        [genes_frame(rows[0]), pl.DataFrame({"hgnc_id": ["HGNC:90099"], "mane_select": ["NM_1.1"]})]
    )


@pytest.mark.parametrize(
    ("corrupt", "expected", "fragment"),
    [
        pytest.param(_flip_exon_base, {"MAN001"}, "differ from the transcript", id="flipped-base"),
        pytest.param(
            _gap_in_transcript_coordinates, {"MAN001"}, "expected 36", id="coordinate-gap"
        ),
        pytest.param(_truncated_window, {"MAN001"}, "window_sequence has", id="short-window"),
        pytest.param(
            _exon_longer_in_genome_than_in_transcript,
            {"MAN001"},
            "in the genome",
            id="span-mismatch",
        ),
        pytest.param(_change_a_residue, {"MAN002"}, "residue 5", id="changed-residue"),
        pytest.param(_cds_not_a_multiple_of_three, {"MAN002"}, "multiple of 3", id="cds-length"),
        pytest.param(_cds_without_its_stop, {"MAN002"}, "stop codon", id="no-stop"),
        pytest.param(_version_mismatch, {"MAN003"}, "NM_900001.3", id="version-mismatch"),
        pytest.param(_null_mane_select, {"MAN003"}, "None", id="null-mane-select"),
        pytest.param(_registry_gene_without_a_row, {"MAN004"}, "HGNC:90099", id="no-row"),
    ],
)
def test_each_mane_check_fires_on_its_own_corruption(
    tmp_path: Path,
    corrupt: Callable[[Rows], pl.DataFrame],
    expected: set[str],
    fragment: str,
) -> None:
    rows = synthetic_rows()
    genes = corrupt(rows)
    issues = validate_mane(_frame(tmp_path, rows), genes, tmp_path)
    assert {i.code for i in issues} == expected
    # One issue per corruption, and it says *which* fact failed: a different
    # check firing under the same code would not satisfy the fragment.
    assert len(issues) == 1 and fragment in issues[0].message
    assert issues == sorted(issues)
    assert all(i.severity.value == "error" for i in issues)
    assert "HGNC:" in issues[0].message


def test_a_clean_frame_adds_no_issue(tmp_path: Path) -> None:
    rows = synthetic_rows()
    frame = _frame(tmp_path, rows)
    assert validate_mane(frame, genes_frame(rows[0]), tmp_path) == []
    # Genes in the frame but absent from genes.tsv are not reported.
    assert validate_mane(frame, genes_frame(rows[0][:1]), tmp_path) == []
    # Without a registry only the frame's own checks run.
    assert validate_mane(frame, None, tmp_path) == []


def test_registry_checks_are_skipped_without_a_registry(tmp_path: Path) -> None:
    rows = synthetic_rows()
    rows[0][0]["protein_sequence"] = "G" + rows[0][0]["protein_sequence"][1:]
    issues = validate_mane(_frame(tmp_path, rows), None, tmp_path)
    assert {i.code for i in issues} == {"MAN002"}


def test_a_transcript_with_no_exon_rows_is_reported(tmp_path: Path) -> None:
    rows = synthetic_rows()
    no_exons = (rows[0], [e for e in rows[1] if e["hgnc_id"] != GENE_MINUS])
    issues = validate_mane(_frame(tmp_path, no_exons), None, tmp_path)
    assert [(i.code, GENE_MINUS in i.message) for i in issues] == [("MAN001", True)]
    assert "has no exon rows" in issues[0].message


def _one_gene(overrides: dict[int, str], protein: str) -> Rows:
    t, e = make_gene(
        "HGNC:90010", "NM_900010.1", protein, "1", "+", (30, 40, 26), (1000, 1200, 1500),
        seed=5, codon_overrides=overrides,
    )  # fmt: skip
    return [t], e


@pytest.mark.parametrize(
    ("overrides", "protein", "expected"),
    [
        pytest.param({0: "CTG"}, "MKTAYIAKQRQISFVKSHFSRQ", set(), id="non-ATG-start-reads-as-M"),
        pytest.param({4: "TGA"}, "MKTAUIAKQRQISFVKSHFSRQ", set(), id="TGA-reads-as-U"),
        pytest.param({4: "TGA"}, "MKTACIAKQRQISFVKSHFSRQ", {"MAN002"}, id="TGA-against-C"),
        pytest.param({3: "CTG"}, "MKTMYIAKQRQISFVKSHFSRQ", {"MAN002"}, id="M-allowance-only-at-1"),
        pytest.param({4: "TAA"}, "MKTAYIAKQRQISFVKSHFSRQ", {"MAN002"}, id="internal-stop"),
    ],
)
def test_the_two_translation_allowances_are_exactly_two(
    tmp_path: Path, overrides: dict[int, str], protein: str, expected: set[str]
) -> None:
    rows = _one_gene(overrides, protein)
    issues = validate_mane(_frame(tmp_path, rows), None, tmp_path)
    assert {i.code for i in issues} == expected


def test_the_codon_table_is_the_standard_genetic_code() -> None:
    known = {"ATG": "M", "TGG": "W", "TAA": "*", "TAG": "*", "TGA": "*", "TTT": "F", "GCC": "A"}
    assert all(CODON_TABLE[c] == a for c, a in known.items())
    assert len(CODON_TABLE) == 64 and sum(a == "*" for a in CODON_TABLE.values()) == 3


@pytest.mark.parametrize("missing", ["hgnc_id", "mane_select"])
def test_a_genes_frame_missing_a_needed_column_is_skipped_not_raised(
    tmp_path: Path, missing: str
) -> None:
    """TBL001 reports the column; MAN003/MAN004 must neither raise nor cascade."""
    rows = synthetic_rows()
    frame = _frame(tmp_path, rows)
    genes = genes_frame(rows[0]).drop(missing)
    assert validate_mane(frame, genes, tmp_path) == []
