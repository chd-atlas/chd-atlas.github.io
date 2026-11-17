"""MAN001-MAN004: the MANE frame against itself and against `mirrors/genes.tsv`.

Reports; never raises. MAN001 and MAN002 read the frame alone. MAN003 and MAN004
compare it with the gene registry and run only when one is supplied: an absent or
unreadable `genes.tsv` is reported once, elsewhere (TBL008), and must not cascade
into one MAN004 per gene.

Issues come back sorted, never in dict or set iteration order.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import polars as pl

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.mane import ManeFrame, ManeTranscript

_BASES = "TCAG"
_AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
# The standard genetic code, codons in TCAG order.
CODON_TABLE: dict[str, str] = {
    a + b + c: _AMINO[16 * i + 4 * j + k]
    for i, a in enumerate(_BASES)
    for j, b in enumerate(_BASES)
    for k, c in enumerate(_BASES)
}
_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _man(code: str, location: Path, message: str) -> ValidationIssue:
    return ValidationIssue(code, Severity.ERROR, str(location), message)


def _frame_fault(t: ManeTranscript) -> str | None:
    """The first way the transcript's exons contradict its sequence, or None."""
    if not t.exons:
        return "has no exon rows"
    expected_start = 1
    pieces: list[str] = []
    for exon in t.exons:
        label = f"exon {exon.number}"
        if len(exon.window_sequence) != exon.window_end - exon.window_start + 1:
            return (
                f"{label} window_sequence has {len(exon.window_sequence)} bases but "
                f"window {exon.window_start}-{exon.window_end} spans "
                f"{exon.window_end - exon.window_start + 1}"
            )
        if not (exon.window_start <= exon.genomic_start <= exon.genomic_end <= exon.window_end):
            return f"{label} lies outside its window"
        if exon.transcript_end - exon.transcript_start != exon.genomic_end - exon.genomic_start:
            return (
                f"{label} spans {exon.transcript_end - exon.transcript_start + 1} bases "
                f"in the transcript but {exon.genomic_end - exon.genomic_start + 1} in the genome"
            )
        if exon.transcript_start != expected_start:
            return (
                f"{label} starts at transcript position {exon.transcript_start}, "
                f"expected {expected_start}"
            )
        expected_start = exon.transcript_end + 1
        core = exon.window_sequence[
            exon.genomic_start - exon.window_start : exon.genomic_end - exon.window_start + 1
        ]
        pieces.append(reverse_complement(core) if t.strand == "-" else core)
    if expected_start - 1 != len(t.sequence):
        return (
            f"exons end at transcript position {expected_start - 1} but the "
            f"sequence has {len(t.sequence)} bases"
        )
    joined = "".join(pieces)
    if joined != t.sequence:
        first = next(
            (i + 1 for i, (a, b) in enumerate(zip(joined, t.sequence, strict=False)) if a != b),
            None,
        )
        return f"exon sequences differ from the transcript sequence at position {first}"
    return None


def _cds_fault(t: ManeTranscript) -> str | None:
    """The first way the coding sequence fails to translate to the protein, or None."""
    if not 1 <= t.cds_start <= t.cds_end <= len(t.sequence):
        return f"cds {t.cds_start}-{t.cds_end} is not inside the {len(t.sequence)}-base sequence"
    cds = t.sequence[t.cds_start - 1 : t.cds_end]
    if len(cds) % 3:
        return f"cds length {len(cds)} is not a multiple of 3"
    codons = [cds[i : i + 3] for i in range(0, len(cds), 3)]
    translated = "".join(CODON_TABLE.get(c, "?") for c in codons)
    if "?" in translated:
        return "cds contains a codon with an ambiguous base"
    if translated[-1] != "*":
        return "cds does not end in a stop codon"
    body = translated[:-1]
    if len(body) != len(t.protein_sequence):
        return f"cds encodes {len(body)} residues but the protein has {len(t.protein_sequence)}"
    for index, (got, want, codon) in enumerate(
        zip(body, t.protein_sequence, codons[:-1], strict=True)
    ):
        if got == want:
            continue
        if want == "U" and codon == "TGA":
            continue
        if index == 0 and want == "M" and codon != "ATG":
            continue
        return f"residue {index + 1} is {want} in the protein but {codon} encodes {got}"
    return None


def validate_mane(
    frame: ManeFrame, registry_genes: pl.DataFrame | None, root: Path = Path()
) -> list[ValidationIssue]:
    """MAN001-MAN004. `root` is the repository root, used only for issue locations."""
    mirrors = root / "mirrors"
    issues: list[ValidationIssue] = []
    for hgnc_id in sorted(frame):
        t = frame[hgnc_id]
        fault = _frame_fault(t)
        if fault is not None:
            issues.append(
                _man("MAN001", mirrors / "mane_exons.tsv", f"{hgnc_id} {t.transcript}: {fault}")
            )
        cds_fault = _cds_fault(t)
        if cds_fault is not None:
            issues.append(
                _man(
                    "MAN002",
                    mirrors / "mane_transcripts.tsv",
                    f"{hgnc_id} {t.transcript}: {cds_fault}",
                )
            )
    if registry_genes is not None:
        issues.extend(_against_registry(frame, registry_genes, mirrors / "genes.tsv"))
    return sorted(issues)


def _against_registry(
    frame: ManeFrame, genes: pl.DataFrame, location: Path
) -> Iterable[ValidationIssue]:
    for row in genes.sort("hgnc_id").iter_rows(named=True):
        hgnc_id = row["hgnc_id"]
        declared = row.get("mane_select")
        t = frame.get(hgnc_id)
        if t is None:
            yield _man("MAN004", location, f"{hgnc_id} has no row in mane_transcripts")
        elif declared != t.transcript:
            yield _man(
                "MAN003",
                location,
                f"{hgnc_id} mane_select is {declared!r} but the MANE frame has {t.transcript}",
            )
