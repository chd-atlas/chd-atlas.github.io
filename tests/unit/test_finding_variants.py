"""A curated variant's names against the pinned MANE frame (LIT016-LIT021, LIT025).

The real cases read the committed mirror. Their values were measured by the C0
spike against MANE v1.5 and GRCh38 and are written here as literals, because the
repository may not read the spike's gitignored notes: TBX5 (HGNC:11604) is
NM_181486.4 / NP_852259.1 on `-`, NR2F2 (HGNC:7976) NM_021005.4 / NP_066285.1 on
`+`. The indel cases need a homopolymer, so they use synthetic genes built by
`mane_fixtures.make_gene`, one per strand.
"""

from pathlib import Path
from typing import Any

import pytest

from chd_atlas.hgvs import parse_c
from chd_atlas.mane import (
    Exon,
    ManeFrame,
    ManeTranscript,
    genome_position,
    genome_slice,
    read_mane_frame,
)
from chd_atlas.models.gene_literature import Variant
from chd_atlas.validate.finding_variants import reference_issues
from tests.unit.mane_fixtures import make_gene

REPO_ROOT = Path(__file__).parent.parent.parent
TBX5 = "HGNC:11604"
NR2F2 = "HGNC:7976"

# Porto 2010's variant as VariantValidator normalised it (spec 7.1).
TBX5_R279X: dict[str, Any] = {
    "hgvs_c": "NM_181486.4:c.835C>T",
    "hgvs_p": "NP_852259.1:p.(Arg279Ter)",
    "grch38": "12:114366312:G:A",
    "reported": {"c": "c.835C>T", "p": "p.R279X"},
}


@pytest.fixture(scope="module")
def frame() -> ManeFrame:
    loaded, issues = read_mane_frame(REPO_ROOT)
    assert loaded is not None and issues == []
    return loaded


def _variant(**fields: Any) -> Variant:
    entry: dict[str, Any] = {
        "id": 1,
        "as_reported": "as printed",
        "locator": {"kind": "in_text", "element": "Table 1"},
        "kind": "snv_indel",
        "placement": "mane",
    }
    entry.update(fields)
    return Variant.model_validate(entry)


def _entry(transcript: dict[str, Any], exons: list[dict[str, Any]]) -> ManeTranscript:
    return ManeTranscript(
        hgnc_id=transcript["hgnc_id"],
        transcript=transcript["refseq_transcript"],
        protein=transcript["refseq_protein"],
        chrom=transcript["chrom"],
        strand=transcript["strand"],
        cds_start=transcript["cds_start"],
        cds_end=transcript["cds_end"],
        sequence=transcript["transcript_sequence"],
        protein_sequence=transcript["protein_sequence"],
        exons=tuple(
            Exon(
                number=e["exon"],
                transcript_start=e["transcript_start"],
                transcript_end=e["transcript_end"],
                genomic_start=e["genomic_start"],
                genomic_end=e["genomic_end"],
                window_start=e["window_start"],
                window_end=e["window_end"],
                window_sequence=e["window_sequence"],
            )
            for e in exons
        ),
    )


# `ATG AAA AAA AAA GGC`: c.4-c.12 is a run of nine A between the G of c.3 and of
# c.13, all inside exon 1 (transcript 1-30, the coding sequence from 13).
_RUN_PROTEIN = "MKKKGAYIAKQRQISFVKSH"


def _run_gene(strand: str) -> ManeTranscript:
    starts = (1000, 1200, 1500) if strand == "+" else (3000, 2400, 2000)
    chrom = "1" if strand == "+" else "2"
    row, exons = make_gene(
        "HGNC:90009", "NM_900009.1", _RUN_PROTEIN, chrom, strand, (30, 40, 20), starts, seed=9
    )
    return _entry(row, exons)


def test_tbx5_c835_maps_to_the_genome_on_the_minus_strand(frame: ManeFrame) -> None:
    tbx5 = frame[TBX5]
    assert (tbx5.transcript, tbx5.strand, tbx5.cds_start) == ("NM_181486.4", "-", 556)
    at = genome_position(tbx5, parse_c("c.835C>T").start)
    assert at == 114366312
    assert genome_slice(tbx5, at, at) == "G"
    # c.510 is exon 5's last base (114398573 on `-`), so c.510+1 is the base below it.
    assert genome_position(tbx5, parse_c("c.510+1G>A").start) == 114398572
    # NR2F2 on `+`: the start codon, measured from the MANE GTF and the genome FASTA.
    nr2f2 = frame[NR2F2]
    assert genome_position(nr2f2, parse_c("c.1A>G").start) == 96332106
    assert genome_slice(nr2f2, 96332106, 96332108) == "ATG"


def test_the_spikes_three_mismatches_fire_and_its_matches_pass(frame: ManeFrame) -> None:
    tbx5, nr2f2 = frame[TBX5], frame[NR2F2]
    assert reference_issues(_variant(hgvs_c="NM_181486.4:c.998G>T"), tbx5) == [
        ("LIT017", "hgvs_c NM_181486.4:c.998G>T states G at c.998 but NM_181486.4 has C there")
    ]
    assert reference_issues(_variant(hgvs_c="NM_181486.4:c.*1101C>T"), tbx5) == [
        ("LIT017", "hgvs_c NM_181486.4:c.*1101C>T states C at c.*1101 but NM_181486.4 has A there")
    ]
    assert reference_issues(
        _variant(placement="protein", hgvs_p="NP_066285.1:p.Gln70dup"), nr2f2
    ) == [("LIT019", "hgvs_p NP_066285.1:p.Gln70dup names Q at 70 but NP_066285.1 has K")]

    assert reference_issues(_variant(**TBX5_R279X), tbx5) == []
    # The canonical GT donor: GRCh38 has C at 114398572 and A at 114398571, which
    # on TBX5's strand read G then T. A wrong strand reversal fails here.
    assert genome_slice(tbx5, 114398571, 114398572) == "AC"
    donor = _variant(hgvs_c="NM_181486.4:c.510+1G>A", grch38="12:114398572:C:T")
    assert reference_issues(donor, tbx5) == []


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    ("hgvs_c", "plus_key", "minus_key"),
    [
        # One A deleted from the run: HGVS names the 3'-most, VCF the leftmost.
        ("NM_900009.1:c.12del", "1:1014:GA:G", "2:3005:CT:C"),
        # The change does not depend on where in the run c. names it.
        ("NM_900009.1:c.4del", "1:1014:GA:G", "2:3005:CT:C"),
        ("NM_900009.1:c.12dup", "1:1014:G:GA", "2:3005:C:CT"),
        ("NM_900009.1:c.12_13insA", "1:1014:G:GA", "2:3005:C:CT"),
    ],
)
def test_an_indel_named_three_prime_shifted_matches_its_left_aligned_key(
    strand: str, hgvs_c: str, plus_key: str, minus_key: str
) -> None:
    gene = _run_gene(strand)
    key = plus_key if strand == "+" else minus_key
    assert reference_issues(_variant(hgvs_c=hgvs_c, grch38=key), gene) == []
    # And a key one A short of it is a different change.
    other = "1:1014:GAA:G" if strand == "+" else "2:3005:CTT:C"
    assert [code for code, _ in reference_issues(_variant(hgvs_c=hgvs_c, grch38=other), gene)] == [
        "LIT018"
    ]


@pytest.mark.parametrize(
    ("key", "reason"),
    [
        ("1:1015:AA:A", "not left-aligned: 1:1014:GA:G is the same change"),
        ("1:1015:A:AA", "not left-aligned: 1:1014:G:GA is the same change"),
        ("2:3006:TT:T", "not left-aligned: 2:3005:CT:C is the same change"),
        ("1:1014:GAA:GA", "not minimal: ref and alt end in one base"),
        ("1:1014:GA:GT", "not minimal: ref and alt share a first base beyond the VCF anchor"),
    ],
)
def test_a_key_not_left_aligned_is_refused(key: str, reason: str) -> None:
    gene = _run_gene("+" if key.startswith("1:") else "-")
    assert reference_issues(_variant(grch38=key), gene) == [("LIT018", f"grch38 {key} is {reason}")]


def test_a_paper_on_another_transcript_does_not_trip_lit025(frame: ManeFrame) -> None:
    # NM_000192.3 is another TBX5 transcript; its residue numbers may differ.
    elsewhere = {"transcript": "NM_000192.3", "c": "c.835C>T", "p": "p.R245X"}
    assert reference_issues(_variant(**{**TBX5_R279X, "reported": elsewhere}), frame[TBX5]) == []
    unstated = {"c": "c.835C>T", "p": "p.R245X"}
    assert [
        code
        for code, _ in reference_issues(
            _variant(**{**TBX5_R279X, "reported": unstated}), frame[TBX5]
        )
    ] == ["LIT025"]
    # The MANE transcript, stated with or without its version, is the same claim.
    for stated in ("NM_181486.4", "NM_181486"):
        mane = {**unstated, "transcript": stated}
        found = reference_issues(_variant(**{**TBX5_R279X, "reported": mane}), frame[TBX5])
        assert [code for code, _ in found] == ["LIT025"]


def test_a_variant_on_a_gene_outside_the_frame_is_lit016_until_it_is_unplaced(
    frame: ManeFrame,
) -> None:
    # Owner decision D2 (2026-11-18): a registry gene with no MANE row is a
    # MAN004 warning, and a variant on it cannot be placed or have its residue
    # checked. Listed `unplaced` with no names, it is listed without coordinates.
    absent = [
        (
            "LIT016",
            "the gene has no MANE Select in the pinned frame, so no variant on it can be "
            "placed or its residue checked; list it unplaced",
        )
    ]
    assert reference_issues(_variant(**TBX5_R279X), None) == absent
    assert reference_issues(_variant(placement="protein", hgvs_p="p.Arg279Ter"), None) == absent
    assert reference_issues(_variant(placement="unplaced", hgvs_p="p.Arg279Ter"), None) == absent
    assert reference_issues(_variant(placement="unplaced"), None) == []
    # An unplaced variant on a gene the frame holds keeps only its residue
    # checked: a mismatching paper p. is not LIT025 there.
    kept = _variant(
        placement="unplaced", hgvs_p="NP_852259.1:p.Gln279Ter", reported={"p": "p.R279Q"}
    )
    assert [code for code, _ in reference_issues(kept, frame[TBX5])] == ["LIT019"]
    for placement in ("structural", "outside_region"):
        assert reference_issues(_variant(placement=placement, **TBX5_R279X), None) == []


@pytest.mark.parametrize(
    ("hgvs_p", "paper_p", "codes"),
    [
        # D6: an hgvs_p naming no residue has nothing to compare with the paper.
        ("NP_852259.1:p.?", "p.(Arg170_Gln220del)", []),
        ("NP_852259.1:p.(?)", "p.(Arg170_Gln220del)", []),
        ("NP_852259.1:p.(=)", "p.R279X", []),
        ("NP_852259.1:p.0", "p.R279X", []),
        # An absent hgvs_p is not exempt when the paper names a residue.
        (None, "p.R279X", ["LIT025"]),
        (None, "p.?", []),
        # D5: the paper's short frameshift against VariantValidator's long form.
        ("NP_852259.1:p.(Arg279GlyfsTer5)", "p.R279fs", []),
        ("NP_852259.1:p.(Arg279GlyfsTer5)", "p.R280fs", ["LIT025"]),
    ],
)
def test_lit025_compares_only_what_both_names_state(
    frame: ManeFrame, hgvs_p: str | None, paper_p: str, codes: list[str]
) -> None:
    variant = _variant(
        hgvs_c="NM_181486.4:c.510+1G>A",
        grch38="12:114398572:C:T",
        hgvs_p=hgvs_p,
        reported={"c": "c.510+1G>A", "p": paper_p},
    )
    assert [code for code, _ in reference_issues(variant, frame[TBX5])] == codes


def _with(**changes: Any) -> dict[str, Any]:
    return {**TBX5_R279X, **changes}


@pytest.mark.parametrize(
    ("fields", "code"),
    [
        (_with(hgvs_c="NM_181486.3:c.835C>T"), "LIT016"),
        (_with(hgvs_p="NP_852259.2:p.(Arg279Ter)"), "LIT016"),
        (_with(hgvs_c="NM_181486.4:c.835A>T"), "LIT017"),
        # c.663 is exon 6's last base, so c.663-1 is anchored on no intron. Read
        # as if it were, it would land on an A (measured), so only the anchor
        # check can refuse it.
        ({"hgvs_c": "NM_181486.4:c.663-1A>G"}, "LIT017"),
        (_with(grch38="12:114366312:G:C"), "LIT018"),
        # The right alt on a wrong ref: the two edits agree, so only the ref check sees it.
        (_with(grch38="12:114366312:C:A"), "LIT018"),
        (_with(hgvs_p="NP_852259.1:p.(Gln279Ter)", reported={"c": "c.835C>T"}), "LIT019"),
        ({"placement": "protein", "hgvs_p": "NP_852259.1:p.Arg600Ter"}, "LIT019"),
        ({"hgvs_c": "NM_181486.4:c.510+60G>A"}, "LIT020"),
        ({"hgvs_c": "NM_181486.4:c.*1700C>T"}, "LIT020"),
        ({"hgvs_c": "NM_181486.4:c.-600G>A"}, "LIT020"),
        ({"hgvs_c": "NM_181486.4:c.1600G>A"}, "LIT020"),
        (_with(grch38="12:114370000:G:A"), "LIT020"),
        (_with(hgvs_c="NM_181486.4:c.835_836inv"), "LIT021"),
        (_with(hgvs_p="NP_852259.1:p.Ter519GlnextTer?"), "LIT021"),
        (_with(grch38="chr12:114366312:G:A"), "LIT021"),
        (_with(reported={"c": "c.835C>T", "p": "p.R279Q"}), "LIT025"),
        ({"hgvs_c": "NM_181486.4:c.510+1G>A", "reported": {"p": "p.R279X"}}, "LIT025"),
    ],
)
def test_each_reference_code_fires_alone(
    frame: ManeFrame, fields: dict[str, Any], code: str
) -> None:
    issues = reference_issues(_variant(**fields), frame[TBX5])
    assert issues and {found for found, _ in issues} == {code}, issues
