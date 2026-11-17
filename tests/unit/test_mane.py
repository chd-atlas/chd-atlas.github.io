"""The MANE frame reader: ordering, absence, and frames it cannot build."""

from pathlib import Path
from typing import Any

import polars as pl

from chd_atlas.mane import read_mane_frame
from tests.unit.mane_fixtures import (
    GENE_MINUS,
    GENE_PLUS,
    make_gene,
    synthetic_rows,
    write_mane_mirrors,
)


def test_read_mane_frame_orders_exons_by_transcript_and_returns_none_when_absent(
    tmp_path: Path,
) -> None:
    transcripts, exons = synthetic_rows()
    # Absent: no files, then only one of the two.
    assert read_mane_frame(tmp_path) == (None, [])
    write_mane_mirrors(tmp_path, transcripts, exons)
    (tmp_path / "mirrors" / "mane_exons.tsv").unlink()
    assert read_mane_frame(tmp_path) == (None, [])

    # Written in descending exon order: the reader must sort by exon number, not
    # trust file order. On `-` exon 1 holds the highest genomic coordinates, so
    # exon-number order is descending genomic order there.
    write_mane_mirrors(tmp_path, transcripts, list(reversed(exons)))
    frame, issues = read_mane_frame(tmp_path)
    assert issues == []
    assert frame is not None and sorted(frame) == [GENE_PLUS, GENE_MINUS]
    assert [e.number for e in frame[GENE_PLUS].exons] == [1, 2, 3]
    assert [e.number for e in frame[GENE_MINUS].exons] == [1, 2, 3]
    plus = [e.genomic_start for e in frame[GENE_PLUS].exons]
    minus = [e.genomic_start for e in frame[GENE_MINUS].exons]
    assert plus == sorted(plus) and minus == sorted(minus, reverse=True)
    one = frame[GENE_PLUS]
    assert (one.transcript, one.strand, one.chrom) == ("NM_900001.1", "+", "1")
    first = next(e for e in exons if e["hgnc_id"] == GENE_PLUS and e["exon"] == 1)
    assert one.exons[0].window_sequence == first["window_sequence"]


def test_a_present_but_unreadable_file_returns_its_cause(tmp_path: Path) -> None:
    transcripts, exons = synthetic_rows()
    write_mane_mirrors(tmp_path, transcripts, exons)
    (tmp_path / "mirrors" / "mane_exons.tsv").write_bytes(b"")
    frame, issues = read_mane_frame(tmp_path)
    assert frame is None
    assert [i.code for i in issues] == ["TBL000"]


def test_a_frame_the_reader_cannot_build_returns_none_without_raising(tmp_path: Path) -> None:
    transcripts, exons = synthetic_rows()
    write_mane_mirrors(tmp_path, transcripts, exons)
    path = tmp_path / "mirrors" / "mane_exons.tsv"
    table = pl.read_csv(path, separator="\t", infer_schema_length=0)
    table.drop("window_end").write_csv(path, separator="\t")
    frame, issues = read_mane_frame(tmp_path)
    assert frame is None
    assert [i.code for i in issues] == ["TBL001"]

    # A null where a value is needed.
    write_mane_mirrors(tmp_path, transcripts, exons)
    table = pl.read_csv(path, separator="\t", infer_schema_length=0)
    table.with_columns(
        pl.when(pl.col("exon") == "2")
        .then(None)
        .otherwise(pl.col("genomic_start"))
        .alias("genomic_start")
    ).write_csv(path, separator="\t")
    frame, issues = read_mane_frame(tmp_path)
    assert frame is None
    assert [i.code for i in issues] == ["TBL001"]


def test_exon_rows_that_cannot_be_attached_are_reported_not_dropped(tmp_path: Path) -> None:
    transcripts, exons = synthetic_rows()
    _, stray = make_gene(
        "HGNC:90003", "NM_900003.1", "MKTAYIAKQRQISFVKSHFSRQ", "3", "+",
        (30, 40, 26), (1000, 1200, 1500), seed=3,
    )  # fmt: skip
    wrong_name: list[dict[str, Any]] = [dict(e) for e in exons]
    wrong_name[0]["refseq_transcript"] = "NM_900001.9"
    duplicate = dict(exons[1])
    write_mane_mirrors(tmp_path, transcripts, wrong_name + stray[:1] + [duplicate])
    frame, issues = read_mane_frame(tmp_path)
    assert frame is not None
    assert [i.code for i in issues] == ["MAN001"] * 3
    joined = " ".join(i.message for i in issues)
    assert "HGNC:90003" in joined and "NM_900001.9" in joined and "more than once" in joined
    # None of the three reached the frame: exon 1 of the first gene is gone, so
    # the validator will see a gap rather than a wrong exon.
    assert [e.number for e in frame[GENE_PLUS].exons] == [2, 3]


def test_a_gene_listed_twice_in_mane_transcripts_is_reported_and_the_first_row_wins(
    tmp_path: Path,
) -> None:
    transcripts, exons = synthetic_rows()
    second = dict(transcripts[0], refseq_transcript="NM_900001.9", refseq_protein="NP_555555.1")
    write_mane_mirrors(tmp_path, [transcripts[0], second, transcripts[1]], exons)
    frame, issues = read_mane_frame(tmp_path)
    assert frame is not None
    assert [(i.code, i.message) for i in issues] == [
        ("MAN001", "HGNC:90001 appears in mane_transcripts more than once")
    ]
    # The first row answers, for the frame and for the exons matched against it.
    assert (frame[GENE_PLUS].transcript, frame[GENE_PLUS].protein) == (
        "NM_900001.1",
        "NP_211112.9",
    )
    assert [e.number for e in frame[GENE_PLUS].exons] == [1, 2, 3]
