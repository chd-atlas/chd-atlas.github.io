"""The committed MANE mirror, read against the registry and against two anchors.

`mirrors/mane_transcripts.tsv` and `mirrors/mane_exons.tsv` are written only by
`chd-atlas ingest mane_select`, so these tests read the committed files. They
neither read `mirrors/sources.yaml` nor the package manifest, so they do not
catch a pin changed without its release (SRC005/SRC006 do). What they catch is
a mirror edited by hand: a gene outside the registry that MAN004 never looks
for, a transcript, strand, CDS start or genomic window that no longer sits where
an independent source puts it (CLAUDE.md section 4.23). That every registry gene
has a clean frame is `test_repository_validates.py`'s, which runs
`validate_repository` over this same tree: 0 errors (MAN001-MAN003, MAN005) and
exactly the baseline warnings, so a registry gene without a row (MAN004, a
warning) would raise its warning count.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from chd_atlas.mane import ManeFrame, read_mane_frame
from chd_atlas.tables import TABLE_SCHEMAS, read_table

REPO_ROOT = Path(__file__).resolve().parents[1]


def _frame() -> ManeFrame:
    frame, read_issues = read_mane_frame(REPO_ROOT)
    assert frame is not None, read_issues
    return frame


def _forward_bases(frame: ManeFrame, chrom: str, start: int, end: int) -> list[str]:
    """The forward-strand bases of `chrom:start-end`, one string per exon window
    that covers the whole span (an exon of any transcript, either strand)."""
    return [
        e.window_sequence[start - e.window_start : end - e.window_start + 1]
        for t in frame.values()
        if t.chrom == chrom
        for e in t.exons
        if e.window_start <= start and end <= e.window_end
    ]


def test_the_frame_holds_exactly_the_registry_genes() -> None:
    """MAN004 reports a registry gene with no row; nothing reports a row for a gene
    outside the registry, which is what a hand-edited mirror would add."""
    registry, issues = read_table(REPO_ROOT / "mirrors" / "genes.tsv", TABLE_SCHEMAS["genes"])
    assert registry is not None, issues
    assert isinstance(registry, pl.DataFrame)
    assert set(_frame()) == set(registry["hgnc_id"].to_list())


def test_tbx5_and_nr2f2_sit_where_the_spike_measured() -> None:
    """Two anchors, one per strand.

    TBX5 (`-`): the forward-strand base at 12:114366312 is `G`, read from the one
    exon window that contains it (TBX5 exon 8, window 114366115-114366441).

    NR2F2 (`+`): the start codon of `NM_021005.4` is chr15:96332106-96332108 and
    reads `ATG` (flank `GAT|ATG|GCA`). That coordinate was measured by the
    controller on 2026-11-18 from feeds' MANE v1.5 GTF `start_codon` feature and
    the GRCh38.p14 genome FASTA, independently of this mirror.
    """
    frame = _frame()
    by_name = {t.transcript: t for t in frame.values()}

    tbx5 = by_name["NM_181486.4"]
    assert (tbx5.strand, tbx5.cds_start, len(tbx5.protein_sequence)) == ("-", 556, 518)
    assert _forward_bases(frame, "12", 114_366_312, 114_366_312) == ["G"]

    nr2f2 = by_name["NM_021005.4"]
    assert (nr2f2.chrom, nr2f2.strand) == ("15", "+")
    assert (nr2f2.cds_start, len(nr2f2.protein_sequence)) == (1407, 414)
    assert _forward_bases(frame, "15", 96_332_106, 96_332_108) == ["ATG"]
