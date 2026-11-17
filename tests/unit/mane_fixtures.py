"""A synthetic MANE frame, derived so every part of it agrees by construction.

`make_gene` first writes a transcript (5' UTR, a coding sequence back-translated
from a protein, 3' UTR) and then lays its exons into a made-up chromosome, on the
named strand, at the given genomic positions. Windows, exon cores and the
transcript are therefore consistent without any of them being typed twice.
`write_mane_mirrors` writes the two TSVs; `genes_frame` builds the matching
`genes.tsv` columns `validate_mane` reads.

Public surface (shared with later MANE tasks):
    synthetic_rows() -> tuple[list[dict], list[dict]]    # (transcripts, exons), two genes
    write_mane_mirrors(root, transcripts, exons) -> None
    genes_frame(transcripts) -> pl.DataFrame
    GENE_PLUS, GENE_MINUS                                  # the two HGNC ids
"""

from __future__ import annotations

import random
from pathlib import Path
from typing import Any

import polars as pl

from chd_atlas.tables import TABLE_SCHEMAS

GENE_PLUS = "HGNC:90001"
GENE_MINUS = "HGNC:90002"
FLANK = 50

# One codon per residue, written out here rather than inverted from the code
# under test.
_CODON = {
    "M": "ATG",
    "K": "AAA",
    "T": "ACC",
    "A": "GCC",
    "Y": "TAT",
    "I": "ATC",
    "Q": "CAG",
    "R": "CGG",
    "S": "TCC",
    "F": "TTC",
    "V": "GTG",
    "H": "CAC",
    "L": "CTG",
    "N": "AAC",
    "P": "CCC",
    "D": "GAC",
    "E": "GAG",
    "G": "GGC",
    "W": "TGG",
    "C": "TGC",
    "*": "TAA",
}
_COMP = str.maketrans("ACGT", "TGCA")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def _filler(rng: random.Random, n: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(n))


def make_gene(
    hgnc_id: str,
    transcript: str,
    protein: str,
    chrom: str,
    strand: str,
    exon_lengths: tuple[int, ...],
    exon_genomic_starts: tuple[int, ...],
    seed: int,
    utr5: int = 12,
    utr3: int = 15,
    codon_overrides: dict[int, str] | None = None,
    refseq_protein: str | None = None,
    flanks: tuple[tuple[int, int], ...] | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """One transcript row and its exon rows (exons in transcript order).

    `exon_genomic_starts[i]` is the forward-strand start of exon i+1's core;
    on `-` they must therefore descend. `exon_lengths` must sum to the
    transcript length (`utr5 + 3 * (len(protein) + 1) + utr3`).
    `codon_overrides` maps a 0-based residue index to the codon to use there
    instead of the default (a non-ATG start, or TGA for a selenocysteine).
    `refseq_protein` defaults to an accession whose number is unrelated to the
    transcript's, as in real MANE rows. `flanks[i]` is exon i+1's
    `(left, right)` flank in bases; the default is FLANK on both sides.
    """
    rng = random.Random(seed)
    overrides = codon_overrides or {}
    if refseq_protein is None:
        refseq_protein = f"NP_{(int(transcript[3:].split('.')[0]) + 111111) % 900000 + 100000}.9"
    sides = flanks or ((FLANK, FLANK),) * len(exon_lengths)
    cds = (
        "".join(overrides[i] if i in overrides else _CODON[r] for i, r in enumerate(protein))
        + _CODON["*"]
    )
    sequence = _filler(rng, utr5) + cds + _filler(rng, utr3)
    if sum(exon_lengths) != len(sequence):
        raise ValueError(f"exon lengths sum to {sum(exon_lengths)}, transcript is {len(sequence)}")
    genome_length = max(exon_genomic_starts) + max(exon_lengths) + 2 * FLANK + 100
    genome = list(_filler(rng, genome_length))  # 1-based positions: genome[p - 1]
    exons: list[dict[str, Any]] = []
    cursor = 1
    for number, (length, g_start, (left, right)) in enumerate(
        zip(exon_lengths, exon_genomic_starts, sides, strict=True), start=1
    ):
        piece = sequence[cursor - 1 : cursor - 1 + length]
        forward = _revcomp(piece) if strand == "-" else piece
        genome[g_start - 1 : g_start - 1 + length] = forward
        exons.append(
            {
                "hgnc_id": hgnc_id,
                "refseq_transcript": transcript,
                "exon": number,
                "transcript_start": cursor,
                "transcript_end": cursor + length - 1,
                "genomic_start": g_start,
                "genomic_end": g_start + length - 1,
                "window_start": g_start - left,
                "window_end": g_start + length - 1 + right,
                "source": "synthetic",
            }
        )
        cursor += length
    text = "".join(genome)
    for row in exons:
        row["window_sequence"] = text[row["window_start"] - 1 : row["window_end"]]
    transcript_row = {
        "hgnc_id": hgnc_id,
        "refseq_transcript": transcript,
        "refseq_protein": refseq_protein,
        "ensembl_transcript": f"ENST{hgnc_id[5:].zfill(11)}.1",
        "ensembl_protein": f"ENSP{hgnc_id[5:].zfill(11)}.1",
        "chrom": chrom,
        "strand": strand,
        "cds_start": utr5 + 1,
        "cds_end": utr5 + len(cds),
        "transcript_sequence": sequence,
        "protein_sequence": protein,
        "source": "synthetic",
    }
    return transcript_row, exons


def synthetic_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Two genes, three exons each: GENE_PLUS on chromosome 1 `+`, GENE_MINUS on 2 `-`.

    Fresh dicts on every call, sorted as the mirror's sort keys require.
    """
    t1, e1 = make_gene(
        GENE_PLUS,
        "NM_900001.1",
        "MKTAYIAKQRQISFVKSHFSRQ",
        "1",
        "+",
        (30, 40, 26),
        (1000, 1200, 1500),
        seed=1,
    )
    # `-`: exon 1 (transcript 5' end) has the highest genomic coordinates.
    t2, e2 = make_gene(
        GENE_MINUS,
        "NM_900002.3",
        "MSTNPKPQRKTKRNTNRRPQDVKFPGW",
        "2",
        "-",
        (35, 45, 31),
        (3000, 2400, 2000),
        seed=2,
    )
    return [t1, t2], e1 + e2


def write_mane_mirrors(
    root: Path, transcripts: list[dict[str, Any]], exons: list[dict[str, Any]]
) -> None:
    """Write `mirrors/mane_transcripts.tsv` and `mirrors/mane_exons.tsv` under `root`."""
    mirrors = root / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    for name, rows in (("mane_transcripts", transcripts), ("mane_exons", exons)):
        columns = TABLE_SCHEMAS[name].column_names
        frame = pl.DataFrame(
            {c: [r[c] for r in rows] for c in columns},
            schema={c: TABLE_SCHEMAS[name].polars_overrides().get(c, pl.String) for c in columns},
        )
        frame.write_csv(mirrors / f"{name}.tsv", separator="\t")


def genes_frame(transcripts: list[dict[str, Any]]) -> pl.DataFrame:
    """The `hgnc_id` and `mane_select` columns of a `genes.tsv` that agrees with the frame."""
    return pl.DataFrame(
        {
            "hgnc_id": [t["hgnc_id"] for t in transcripts],
            "mane_select": [t["refseq_transcript"] for t in transcripts],
        }
    )
