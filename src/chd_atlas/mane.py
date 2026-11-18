"""The pinned MANE frame: each gene's MANE Select transcript, its exons and genome windows.

Reads the two flat mirror tables `mane_transcripts` and `mane_exons` into
`ManeTranscript` records. This module only reads; `validate/mane.py` checks what
it returns.

Absent versus unreadable is a distinction the callers rely on. `read_mane_frame`
returns `(None, [])` for an absent file and `(None, issues)` for a present one
that cannot be read. Absence is not an error here, because a repository that has
not pinned `mane_select` carries no mirror; once `mirrors/sources.yaml` pins it,
absence is an error, reported by the runner as MAN005, not by this reader.

Exon rows that cannot be attached to a transcript are never silently dropped:
an exon whose gene has no transcript row, whose `refseq_transcript` differs from
its gene's, or whose number repeats is returned as a MAN001 issue (and left out
of the frame, so no later check reads it as that transcript's exon).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import polars as pl

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import TABLE_SCHEMAS, read_table

TRANSCRIPTS_FILE = "mane_transcripts.tsv"
EXONS_FILE = "mane_exons.tsv"


@dataclass(frozen=True)
class Exon:
    """One exon. All coordinates are 1-based and inclusive.

    `transcript_*` are positions in the transcript sequence, `genomic_*` on the
    forward strand of the chromosome. On `-`, `transcript_start` pairs with
    `genomic_end` (the exon's 5' base is its highest genomic coordinate) and
    `transcript_end` with `genomic_start`. `window_sequence` is the forward-strand
    genome from `window_start` to `window_end`: the exon plus its flank.
    """

    number: int
    transcript_start: int
    transcript_end: int
    genomic_start: int
    genomic_end: int
    window_start: int
    window_end: int
    window_sequence: str


@dataclass(frozen=True)
class ManeTranscript:
    """A MANE Select transcript. `exons` are in exon-number order, which MAN001
    checks is transcript order (on `-`, exon 1 then holds the highest genomic
    coordinates). `cds_start` is the first base of the
    start codon and `cds_end` the last base of the stop codon."""

    hgnc_id: str
    transcript: str
    protein: str
    chrom: str
    strand: str
    cds_start: int
    cds_end: int
    sequence: str
    protein_sequence: str
    exons: tuple[Exon, ...]


ManeFrame = Mapping[str, ManeTranscript]


def _unusable(path: Path, why: str) -> list[ValidationIssue]:
    return [
        ValidationIssue("TBL001", Severity.ERROR, str(path), f"unusable as the MANE frame: {why}")
    ]


def read_mane_frame(root: Path) -> tuple[ManeFrame | None, list[ValidationIssue]]:
    """Read the MANE frame under `root/mirrors/`.

    Returns `(None, [])` when either file is absent, `(None, issues)` when a
    present file cannot be read or lacks a column or a value the frame needs,
    and otherwise `(frame, issues)` where `issues` holds only MAN001 reports of
    exon rows that could not be attached (see the module docstring).
    """
    t_path = root / "mirrors" / TRANSCRIPTS_FILE
    e_path = root / "mirrors" / EXONS_FILE
    if not t_path.is_file() or not e_path.is_file():
        return None, []
    tables: list[pl.DataFrame] = []
    for path, name in ((t_path, "mane_transcripts"), (e_path, "mane_exons")):
        schema = TABLE_SCHEMAS[name]
        table, read_issues = read_table(path, schema)
        if table is None:
            return None, read_issues
        missing = [c for c in schema.column_names if c not in table.columns]
        if missing:
            return None, _unusable(path, f"missing column(s) {missing}")
        needed = [c for c in schema.column_names if c != "source"]
        if any(table[c].null_count() for c in needed):
            return None, _unusable(path, "a required value is empty")
        tables.append(table)
    transcripts, exon_rows = tables

    issues: list[ValidationIssue] = []
    exons_by_gene: dict[str, dict[int, Exon]] = {}
    # First row wins, as it does for `frame_out` below, so an exon is matched
    # against the transcript the frame will hold.
    transcript_name: dict[str, str] = {}
    for r in transcripts.iter_rows(named=True):
        transcript_name.setdefault(r["hgnc_id"], r["refseq_transcript"])
    for r in exon_rows.iter_rows(named=True):
        gene, number = r["hgnc_id"], r["exon"]
        if gene not in transcript_name:
            issues.append(
                ValidationIssue(
                    "MAN001",
                    Severity.ERROR,
                    str(e_path),
                    f"{gene} exon {number} has no row in mane_transcripts",
                )
            )
            continue
        if r["refseq_transcript"] != transcript_name[gene]:
            issues.append(
                ValidationIssue(
                    "MAN001",
                    Severity.ERROR,
                    str(e_path),
                    f"{gene} exon {number} names {r['refseq_transcript']} but the "
                    f"gene's transcript is {transcript_name[gene]}",
                )
            )
            continue
        known = exons_by_gene.setdefault(gene, {})
        if number in known:
            issues.append(
                ValidationIssue(
                    "MAN001",
                    Severity.ERROR,
                    str(e_path),
                    f"{gene} exon {number} appears more than once",
                )
            )
            continue
        known[number] = Exon(
            number=number,
            transcript_start=r["transcript_start"],
            transcript_end=r["transcript_end"],
            genomic_start=r["genomic_start"],
            genomic_end=r["genomic_end"],
            window_start=r["window_start"],
            window_end=r["window_end"],
            window_sequence=r["window_sequence"],
        )

    frame_out: dict[str, ManeTranscript] = {}
    for r in transcripts.iter_rows(named=True):
        gene = r["hgnc_id"]
        if gene in frame_out:
            issues.append(
                ValidationIssue(
                    "MAN001",
                    Severity.ERROR,
                    str(t_path),
                    f"{gene} appears in mane_transcripts more than once",
                )
            )
            continue
        exons = exons_by_gene.get(gene, {})
        frame_out[gene] = ManeTranscript(
            hgnc_id=gene,
            transcript=r["refseq_transcript"],
            protein=r["refseq_protein"],
            chrom=r["chrom"],
            strand=r["strand"],
            cds_start=r["cds_start"],
            cds_end=r["cds_end"],
            sequence=r["transcript_sequence"],
            protein_sequence=r["protein_sequence"],
            exons=tuple(exons[n] for n in sorted(exons)),
        )
    return frame_out, sorted(issues)
