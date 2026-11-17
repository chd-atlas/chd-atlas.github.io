"""The committed MANE mirror against the registry and against the spike's measurements.

`mirrors/mane_transcripts.tsv` and `mirrors/mane_exons.tsv` are written only by
`chd-atlas ingest mane_select`, so these tests read the committed files: they
catch a mirror edited by hand or a pin changed without its release, not a recipe
defect (CLAUDE.md section 4.23).
"""

from __future__ import annotations

from pathlib import Path

import polars as pl

from chd_atlas.mane import read_mane_frame
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.mane import validate_mane

REPO_ROOT = Path(__file__).resolve().parents[1]


def _registry() -> pl.DataFrame:
    genes, issues = read_table(REPO_ROOT / "mirrors" / "genes.tsv", TABLE_SCHEMAS["genes"])
    assert genes is not None, issues
    return genes


def test_every_registry_gene_has_its_mane_select_and_the_frame_checks_clean() -> None:
    """Measured 2026-11-18 on the committed mirror: 186 registry genes, 186
    transcripts, 2,857 exons, every `genes.tsv` `mane_select` equal to the
    frame's, and one MAN002 allowance in use (VEGFA, NM_003376.6, whose protein
    starts at M but whose first codon is CTG); no selenoprotein `U` against `TGA`."""
    frame, read_issues = read_mane_frame(REPO_ROOT)
    assert frame is not None, read_issues
    registry = _registry()
    assert registry.height == 186
    assert len(frame) == 186
    assert set(frame) == set(registry["hgnc_id"].to_list())
    assert validate_mane(frame, registry, REPO_ROOT) == []


def test_tbx5_and_nr2f2_sit_where_the_spike_measured() -> None:
    frame, read_issues = read_mane_frame(REPO_ROOT)
    assert frame is not None, read_issues
    by_name = {t.transcript: t for t in frame.values()}

    tbx5 = by_name["NM_181486.4"]
    assert (tbx5.strand, tbx5.cds_start, len(tbx5.protein_sequence)) == ("-", 556, 518)

    nr2f2 = by_name["NM_021005.4"]
    assert (nr2f2.strand, nr2f2.cds_start, len(nr2f2.protein_sequence)) == ("+", 1407, 414)

    # The forward-strand base at 12:114366312, read from the one exon window that
    # contains it (TBX5 exon 8, window 114366115-114366441).
    position = 114_366_312
    bases = [
        e.window_sequence[position - e.window_start]
        for t in frame.values()
        if t.chrom == "12"
        for e in t.exons
        if e.window_start <= position <= e.window_end
    ]
    assert bases == ["G"]
