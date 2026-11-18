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
from chd_atlas.models.gene_literature import Finding, Variant
from chd_atlas.validate.finding_variants import reference_issues, shape_issues
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
        ("1:1018:AA:A", "not left-aligned: 1:1014:GA:G is the same change"),
        ("1:1015:A:AA", "not left-aligned: 1:1014:G:GA is the same change"),
        ("2:3006:TT:T", "not left-aligned: 2:3005:CT:C is the same change"),
        ("1:1014:GAA:GA", "not minimal: ref and alt end in one base"),
        ("1:1014:GA:GT", "not minimal: ref and alt share a first base beyond the VCF anchor"),
    ],
)
def test_a_key_not_left_aligned_is_refused(key: str, reason: str) -> None:
    gene = _run_gene("+" if key.startswith("1:") else "-")
    assert reference_issues(_variant(grch38=key), gene) == [("LIT018", f"grch38 {key} is {reason}")]


def test_a_key_not_left_aligned_is_refused_where_two_windows_meet(frame: ManeFrame) -> None:
    # NM_003737.4 (HGNC:13681, `-`): exon 9's window ends at 11:6631260 and exon
    # 8's begins at 6631261, so no one window holds the key and its shifted
    # twin. Ending in one base is the rule, not a comparison inside a window.
    key = "11:6631261:C:CC"
    assert reference_issues(_variant(grch38=key), frame["HGNC:13681"]) == [
        ("LIT018", f"grch38 {key} is not left-aligned: 11:6631260:G:GC is the same change")
    ]


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
    head = (
        "the gene has no MANE Select in the pinned frame, so no variant on it can be "
        "placed or its residue checked; "
    )
    absent = [("LIT016", head + "list it unplaced, with no hgvs_c or grch38 and no hgvs_p")]
    # An unplaced variant is already unplaced: its only fix is to drop the hgvs_p.
    dropped = [("LIT016", head + "drop its hgvs_p")]
    assert reference_issues(_variant(**TBX5_R279X), None) == absent
    assert reference_issues(_variant(placement="protein", hgvs_p="p.Arg279Ter"), None) == absent
    assert reference_issues(_variant(placement="unplaced", hgvs_p="p.Arg279Ter"), None) == dropped
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
        # D6': the paper's own residue-less p. has nothing to compare either,
        # whatever hgvs_p holds; a paper p. that does name a residue still fires.
        ("NP_852259.1:p.(Arg279Ter)", "p.?", []),
        ("NP_852259.1:p.(Arg279Ter)", "p.(?)", []),
        ("NP_852259.1:p.(Arg279Ter)", "p.=", []),
        ("NP_852259.1:p.(Arg279Ter)", "p.0", []),
        ("NP_852259.1:p.(Arg279Ter)", "p.R279Q", ["LIT025"]),
        # D5: the paper's short frameshift against VariantValidator's long form.
        ("NP_852259.1:p.(Arg279GlyfsTer5)", "p.R279fs", []),
        ("NP_852259.1:p.(Arg279GlyfsTer5)", "p.R280fs", ["LIT025"]),
        # D5': a counted frameshift naming no new residue, as a TBX5 pilot paper
        # wrote c.1333delC.
        ("NP_852259.1:p.(His445MetfsTer137)", "p.H445fsX136", []),
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
    ("gene", "fields"),
    [
        # `+`-strand intron: NR2F2 c.442 ends exon 1 (96332547); c.442+1 is
        # 96332548, where GRCh38 has G (and T at +2: the GT donor).
        (NR2F2, {"hgvs_c": "NM_021005.4:c.442+1G>A", "grch38": "15:96332548:G:A"}),
        # The 5' UTR: TBX5 c.-1 is transcript base 555, C (GRCh38 G at 114403899).
        (TBX5, {"hgvs_c": "NM_181486.4:c.-1C>A", "grch38": "12:114403899:G:T"}),
        # delins on `-`: c.835_836 read CG on both strands, so TT is AA forward.
        (TBX5, {"hgvs_c": "NM_181486.4:c.835_836delinsTT", "grch38": "12:114366311:CG:AA"}),
        # The stop, one past the 518 residues, named as itself and as a stop-loss.
        (TBX5, {"placement": "protein", "hgvs_p": "NP_852259.1:p.(Ter519=)"}),
        (TBX5, {"placement": "protein", "hgvs_p": "NP_852259.1:p.(Ter519GlnextTer?)"}),
        # D7: off the transcript's ends but inside the end exons' windows, on
        # each strand. TBX5 c.-560 is 5 bp upstream (GRCh38 T, A on TBX5's
        # strand) and c.*1622 1 bp downstream (GRCh38 A, T on its strand);
        # NR2F2 c.-1410 is 4 bp upstream (C) and c.*2640 4 bp downstream (A).
        (TBX5, {"hgvs_c": "NM_181486.4:c.-560A>G", "grch38": "12:114406149:T:C"}),
        (TBX5, {"hgvs_c": "NM_181486.4:c.*1622T>C", "grch38": "12:114353910:A:G"}),
        (NR2F2, {"hgvs_c": "NM_021005.4:c.-1410C>T", "grch38": "15:96330696:C:T"}),
        (NR2F2, {"hgvs_c": "NM_021005.4:c.*2640A>G", "grch38": "15:96340262:A:G"}),
    ],
)
def test_names_that_agree_with_the_frame_pass(
    frame: ManeFrame, gene: str, fields: dict[str, Any]
) -> None:
    assert reference_issues(_variant(**fields), frame[gene]) == []


@pytest.mark.parametrize(
    ("fields", "code"),
    [
        (_with(hgvs_c="NM_181486.3:c.835C>T"), "LIT016"),
        (_with(hgvs_p="NP_852259.2:p.(Arg279Ter)"), "LIT016"),
        # A bare name states no transcript or protein, so it is not the MANE one.
        (_with(hgvs_c="c.835C>T"), "LIT016"),
        ({"placement": "protein", "hgvs_p": "p.Arg279Ter"}, "LIT016"),
        (_with(hgvs_c="NM_181486.4:c.835A>T"), "LIT017"),
        # c.835 is C, so a dup naming A misstates it.
        ({"hgvs_c": "NM_181486.4:c.835dupA"}, "LIT017"),
        # Offsets with no intron to point into. Each names the base its position
        # would land on if read anyway (measured), so only the anchor check
        # refuses it: c.663 ends exon 6, read as c.663-1 it would land on an A;
        # c.*1621 ends the last exon (T below it on TBX5's strand) and c.-555
        # starts exon 1 (C above it); D7 places those as c.*1622 and c.-556.
        ({"hgvs_c": "NM_181486.4:c.663-1A>G"}, "LIT017"),
        ({"hgvs_c": "NM_181486.4:c.*1621+1T>G"}, "LIT017"),
        ({"hgvs_c": "NM_181486.4:c.-555-1C>T"}, "LIT017"),
        ({"hgvs_c": "NM_181486.4:c.-560+1A>G"}, "LIT017"),
        # D7: c.-560 is placed 5 bp upstream, so its base is checked (A, not G).
        ({"hgvs_c": "NM_181486.4:c.-560G>A"}, "LIT017"),
        (_with(grch38="12:114366312:G:C"), "LIT018"),
        # The right alt on a wrong ref: the two edits agree, so only the ref check sees it.
        (_with(grch38="12:114366312:C:A"), "LIT018"),
        # A correct key for another variant, in another window: the two names
        # lie in different places.
        (_with(grch38="12:114398572:C:T"), "LIT018"),
        (_with(hgvs_p="NP_852259.1:p.(Gln279Ter)", reported={"c": "c.835C>T"}), "LIT019"),
        ({"placement": "protein", "hgvs_p": "NP_852259.1:p.Arg600Ter"}, "LIT019"),
        # A stop-loss names the stop; residue 518 is not one.
        ({"placement": "protein", "hgvs_p": "NP_852259.1:p.(Ter518GlnextTer?)"}, "LIT019"),
        ({"hgvs_c": "NM_181486.4:c.510+60G>A"}, "LIT020"),
        ({"hgvs_c": "NM_181486.4:c.*1700C>T"}, "LIT020"),
        # 95 bp upstream: beyond exon 1's window, which ends 50 bp out.
        ({"hgvs_c": "NM_181486.4:c.-650G>A"}, "LIT020"),
        (_with(grch38="12:114370000:G:A"), "LIT020"),
        (_with(grch38="11:114366312:G:A"), "LIT020"),
        (_with(hgvs_c="NM_181486.4:c.835_836inv"), "LIT021"),
        # Past the 1,557-base coding sequence, where a position is written c.*N.
        ({"hgvs_c": "NM_181486.4:c.1600G>A"}, "LIT021"),
        ({"hgvs_c": "NM_181486.4:c.836_835del"}, "LIT021"),
        ({"hgvs_c": "NM_181486.4:c.834_836insA"}, "LIT021"),
        (_with(hgvs_p="NP_852259.1:p.Met1ext-5"), "LIT021"),
        (_with(grch38="chr12:114366312:G:A"), "LIT021"),
        (_with(reported={"c": "c.835C>T", "p": "p.R279Q"}), "LIT025"),
        ({"hgvs_c": "NM_181486.4:c.510+1G>A", "reported": {"p": "p.R279X"}}, "LIT025"),
        # A paper p. the checker cannot read still names something.
        (
            {
                "hgvs_c": "NM_181486.4:c.510+1G>A",
                "grch38": "12:114398572:C:T",
                "reported": {"p": "p.R279X;R280X"},
            },
            "LIT025",
        ),
    ],
)
def test_each_reference_code_fires_alone(
    frame: ManeFrame, fields: dict[str, Any], code: str
) -> None:
    issues = reference_issues(_variant(**fields), frame[TBX5])
    assert issues and {found for found, _ in issues} == {code}, issues


# ---- shape, linking and reconciliation (LIT022-LIT024) ----------------------------
#
# The finding carries {de_novo: 2, inherited: 1, origin_unknown: 3} (six in all) and
# the phenotypes HP:0001631 and HP:0001680. Every origin and the total differ from
# each other, so a check keyed on the wrong one cannot pass by coincidence (§4.30).

_MANE_FIELDS: dict[str, Any] = {
    "placement": "mane",
    "hgvs_c": "NM_181486.4:c.835C>T",
    "grch38": "12:114366312:G:A",
    "normalisation": 1,
}
_REVIEW: dict[str, Any] = {
    "review_status": "pending_review",
    "review_note": "numbering differs",
    "review_settled_by": "authors",
}


def _sv(vid: int, **fields: Any) -> dict[str, Any]:
    """A clean `mane` variant carried once, de novo; `fields` overrides."""
    entry: dict[str, Any] = {
        "id": vid,
        "as_reported": "as printed",
        "locator": {"kind": "in_text", "element": "Table 1"},
        "kind": "snv_indel",
        "carriers": {"de_novo": 1, "inherited": 0, "origin_unknown": 0},
        **_MANE_FIELDS,
        # one key per variant, or LIT022's one-name-per-finding rule fires
        "grch38": f"12:{114366300 + vid}:G:A",
    }
    entry.update(fields)
    return entry


def _finding(variants: list[dict[str, Any]], **changes: Any) -> Finding:
    raw: dict[str, Any] = {
        "id": "CHDA:LIT:0000001",
        "publication": "PMID:1001",
        "design": "cohort_screen",
        "methods": ["exome"],
        "count_unit": "probands",
        "screened": 7,
        "carriers": {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
        "variant_classes": ["snv_indel", "cnv"],
        "variants": variants,
        "phenotypes": ["HP:0001631", "HP:0001680"],
        "lesion_groups": ["septal"],
        "syndromic": "both",
        "shares_probands_with": "unknown",
        "locators": [{"kind": "in_text", "element": "Table 2"}],
        "summary": "Six of seven probands carry a variant.",
        "extracted_by": "agent:extract/test",
    }
    raw.update(changes)
    return Finding.model_validate(raw)


def _shape(variants: list[dict[str, Any]], **changes: Any) -> list[tuple[str, int | None, str]]:
    return shape_issues(_finding(variants, **changes), {1})


def _linked(vid: int, other: int, **fields: Any) -> dict[str, Any]:
    return _sv(vid, **{"zygosity": "heterozygous", "with": [other], "phase": "trans", **fields})


_COMMON: dict[str, Any] = {
    "design": "common_variant_association",
    "methods": ["genotyping_array"],
    "cases": 100,
    "controls": 200,
    "count_unit": None,
    "screened": None,
    "carriers": None,
}

# (the one defect, the finding's changes, the message fragment) -- each fires LIT022 alone.
_SHAPE_CASES: list[tuple[list[dict[str, Any]], dict[str, Any], str]] = [
    ([_sv(1, hgvs_c=None)], {}, "needs hgvs_c"),
    ([_sv(1, grch38=None)], {}, "needs grch38"),
    ([_sv(1, normalisation=None)], {}, "needs normalisation"),
    ([_sv(1, kind="cnv")], {}, "is kind snv_indel, not cnv"),
    (
        [
            _sv(
                1,
                placement="protein",
                hgvs_p="NP_852259.1:p.(Arg279Ter)",
                grch38=None,
                normalisation=None,
                hgvs_c=None,
            )
            | {"hgvs_p": None}
        ],
        {},
        "needs hgvs_p",
    ),
    (
        [
            _sv(
                1,
                placement="protein",
                hgvs_p="NP_852259.1:p.(Arg279Ter)",
                grch38=None,
                normalisation=None,
            )
        ],
        {},
        "must leave hgvs_c null",
    ),
    (
        [
            _sv(
                1,
                placement="protein",
                hgvs_p="NP_852259.1:p.(Arg279Ter)",
                hgvs_c=None,
                normalisation=None,
            )
        ],
        {},
        "must leave grch38 null",
    ),
    (
        [_sv(1, placement="protein", hgvs_p="NP_852259.1:p.(Arg279Ter)", hgvs_c=None, grch38=None)],
        {},
        "must leave normalisation null",
    ),
    (
        [_sv(1, placement="unplaced", hgvs_c=None, grch38=None, normalisation=None)],
        {},
        "needs review_status",
    ),
    (
        [_sv(1, placement="unplaced", grch38=None, normalisation=None, **_REVIEW)],
        {},
        "must leave hgvs_c null",
    ),
    (
        [_sv(1, placement="unplaced", hgvs_c=None, normalisation=None, **_REVIEW)],
        {},
        "must leave grch38 null",
    ),
    (
        [
            _sv(
                1,
                placement="structural",
                kind="snv_indel",
                hgvs_c=None,
                grch38=None,
                normalisation=None,
            )
        ],
        {},
        "is kind chromosomal or cnv, not snv_indel",
    ),
    (
        [_sv(1, placement="structural", kind="cnv", grch38=None, normalisation=None)],
        {},
        "must leave hgvs_c null",
    ),
    (
        [
            _sv(
                1,
                placement="structural",
                kind="cnv",
                hgvs_c=None,
                grch38=None,
                normalisation=None,
                hgvs_p="NP_852259.1:p.(Arg279Ter)",
            )
        ],
        {},
        "must leave hgvs_p null",
    ),
    (
        [_sv(1, placement="structural", kind="cnv", hgvs_c=None, normalisation=None)],
        {},
        "must leave grch38 null",
    ),
    (
        [
            _sv(
                1,
                placement="outside_region",
                kind="cnv",
                hgvs_c=None,
                grch38=None,
                normalisation=None,
            )
        ],
        {},
        "is kind snv_indel, not cnv",
    ),
    (
        [_sv(1, placement="outside_region", grch38=None, normalisation=None)],
        {},
        "must leave hgvs_c null",
    ),
    (
        [
            _sv(
                1,
                placement="outside_region",
                hgvs_c=None,
                grch38=None,
                normalisation=None,
                hgvs_p="NP_852259.1:p.(Arg279Ter)",
            )
        ],
        {},
        "must leave hgvs_p null",
    ),
    (
        [_sv(1, placement="outside_region", hgvs_c=None, normalisation=None)],
        {},
        "must leave grch38 null",
    ),
    (
        [_sv(1, placement="structural", kind="cnv", hgvs_c=None, grch38=None, normalisation=None)],
        {"variant_classes": ["snv_indel"]},
        "not among the finding's variant_classes",
    ),
    ([_sv(1, normalisation=9)], {}, "normalisation 9 is not in this file"),
    ([_sv(1), _sv(1)], {}, "variant id 1 is used more than once"),
    ([_sv(1, carriers=None)], {}, "carriers are required"),
    (
        [_sv(1, carriers={"de_novo": 1, "inherited": 0, "origin_unknown": 0})],
        _COMMON,
        "carriers on a common_variant_association",
    ),
    # `with` rules
    ([_linked(1, 9), _sv(2)], {}, "names variant 9, which this finding does not have"),
    ([_linked(1, 1)], {}, "names the variant itself"),
    (
        [_linked(1, 2, phase="not_stated"), _sv(2, zygosity="heterozygous")],
        {},
        "does not name 1 back",
    ),
    ([_linked(1, 2, phase="cis"), _linked(2, 1)], {}, "differ in phase"),
    ([_sv(1, phase="trans")], {}, "phase is only for a variant with a with link"),
    ([_linked(1, 2, zygosity="not_stated"), _linked(2, 1)], {}, "must be heterozygous"),
    ([_linked(1, 2, zygosity="homozygous"), _linked(2, 1)], {}, "has no partner"),
    ([_linked(1, 2, zygosity="hemizygous"), _linked(2, 1)], {}, "has no partner"),
    (
        [
            _linked(1, 2),
            _linked(2, 1, carriers={"de_novo": 1, "inherited": 1, "origin_unknown": 0}),
        ],
        {},
        "carry different totals",
    ),
    # D9: three variants cannot be pairwise in trans on a diploid genome.
    (
        [
            _linked(1, 2, **{"with": [2, 3]}),
            _linked(2, 1, **{"with": [1, 3]}),
            _linked(3, 1, **{"with": [1, 2]}),
        ],
        {},
        "cannot hold more than two variants pairwise in trans",
    ),
    # A chain 1-2, 2-3 with no 1-3: `with` lists every other entry the same carriers carry.
    (
        [
            _linked(1, 2, phase="cis"),
            _linked(2, 1, phase="cis", **{"with": [1, 3]}),
            _linked(3, 2, phase="cis"),
        ],
        {},
        "neither names the other",
    ),
]


@pytest.mark.parametrize(("variants", "changes", "fragment"), _SHAPE_CASES)
def test_each_shape_rule_fires_alone(
    variants: list[dict[str, Any]], changes: dict[str, Any], fragment: str
) -> None:
    issues = _shape(variants, **changes)
    assert [(code, fragment in message) for code, _, message in issues] == [("LIT022", True)], (
        issues
    )
    assert issues[0][1] is not None  # attributed to a variant, so the location can name it


def test_a_clean_finding_and_a_finding_with_no_variants_have_no_shape_issue() -> None:
    assert _shape([]) == []
    assert (
        _shape([_sv(1), _sv(2, carriers={"de_novo": 0, "inherited": 1, "origin_unknown": 0})]) == []
    )
    # A placement with its fields null is whole: nothing is inferred from nulls.
    assert (
        _shape(
            [
                _sv(
                    1,
                    placement="structural",
                    kind="cnv",
                    hgvs_c=None,
                    grch38=None,
                    normalisation=None,
                ),
                _sv(
                    2,
                    placement="outside_region",
                    hgvs_c=None,
                    grch38=None,
                    normalisation=None,
                    carriers={"de_novo": 0, "inherited": 1, "origin_unknown": 0},
                ),
                _sv(
                    3,
                    placement="unplaced",
                    hgvs_c=None,
                    grch38=None,
                    normalisation=None,
                    carriers={"de_novo": 0, "inherited": 0, "origin_unknown": 1},
                    **_REVIEW,
                ),
            ],
            variant_classes=["snv_indel", "cnv"],
        )
        == []
    )
    # A common-variant finding carries no per-variant counts.
    assert _shape([_sv(1, carriers=None)], **_COMMON) == []


def test_a_compound_heterozygote_with_mixed_origins_reconciles() -> None:
    """One person, two linked entries (de novo and inherited, one each).

    The finding names a single de novo carrier. Counted per entry, or per origin
    over linked entries, this overshoots; counted as one group it fits.
    """
    variants = [
        _linked(1, 2, carriers={"de_novo": 1, "inherited": 0, "origin_unknown": 0}),
        _linked(2, 1, carriers={"de_novo": 0, "inherited": 1, "origin_unknown": 0}),
    ]
    one_person = {"carriers": {"de_novo": 1, "inherited": 0, "origin_unknown": 0}}
    assert _shape(variants, **one_person) == []
    # The same pair beside a second, unlinked de novo person is two people: refused.
    issues = _shape([*variants, _sv(3)], **one_person)
    assert [code for code, _, _ in issues] == ["LIT023"]
    assert "name 2 carrier" in issues[0][2]


def test_a_finding_whose_variants_name_fewer_carriers_than_it_is_valid() -> None:
    """LIT:0000094's shape: 78 carriers and none of them named, or only some."""
    assert (
        _shape([], carriers={"de_novo": 30, "inherited": 20, "origin_unknown": 28}, screened=90)
        == []
    )
    assert (
        _shape(
            [_sv(1)], carriers={"de_novo": 30, "inherited": 20, "origin_unknown": 28}, screened=90
        )
        == []
    )


@pytest.mark.parametrize("origin", ["de_novo", "inherited"])
def test_variants_naming_more_carriers_than_the_finding_fire_lit023(origin: str) -> None:
    named = {"de_novo": 0, "inherited": 0, "origin_unknown": 0, origin: 1}
    other = "inherited" if origin == "de_novo" else "de_novo"
    inherited = {"de_novo": 0, "inherited": 1, "origin_unknown": 0}
    unknown = {"de_novo": 0, "inherited": 0, "origin_unknown": 1}
    # Each origin is checked on its own: three of this origin against the finding's
    # one, with a total (3) that is well inside the finding's six.
    finding = {"carriers": {"de_novo": 1, "inherited": 1, "origin_unknown": 4}}
    per_origin = _shape([_sv(i, carriers=named) for i in (1, 2, 3)], **finding)
    assert [(code, origin in message) for code, _, message in per_origin] == [("LIT023", True)]
    assert other not in per_origin[0][2]
    assert per_origin[0][1] is None  # about the finding's variants as a whole
    # The total is checked on its own: every origin fits, the sum does not. The
    # finding names {1, 0, 0}; one unlinked de novo plus one linked pair is two
    # people, and the pair (one inherited allele, one unknown) adds nothing to any
    # origin's lower bound.
    pair = [
        _linked(2, 3, carriers=inherited),
        _linked(3, 2, carriers=unknown),
    ]
    total = _shape(
        [_sv(1), *pair], carriers={"de_novo": 1, "inherited": 0, "origin_unknown": 0}, screened=1
    )
    assert [(code, "carrier(s)" in message) for code, _, message in total] == [("LIT023", True)]
    assert "name 2 carrier" in total[0][2] and "finding's 1" in total[0][2]
    # A linked pair is one person: the pair and two unlinked entries name three
    # people, which the finding's three de novo carriers cover.
    assert (
        _shape(
            [
                _linked(1, 2),
                _linked(2, 1),
                _sv(3),
                _sv(4),
            ],
            carriers={"de_novo": 3, "inherited": 0, "origin_unknown": 0},
            screened=3,
        )
        == []
    )


def test_a_variant_phenotype_outside_the_findings_fires_lit024() -> None:
    issues = _shape(
        [_sv(1, phenotypes=["HP:0001631"]), _sv(2, phenotypes=["HP:0001631", "HP:0009999"])]
    )
    assert [(code, number) for code, number, _ in issues] == [("LIT024", 2)]
    assert "HP:0009999" in issues[0][2]


_ONE = {"de_novo": 1, "inherited": 0, "origin_unknown": 0}
_INH = {"de_novo": 0, "inherited": 1, "origin_unknown": 0}
_UNK = {"de_novo": 0, "inherited": 0, "origin_unknown": 1}


@pytest.mark.parametrize("counted_as", [_ONE, _INH])
def test_a_mixed_origin_compound_heterozygote_passes_however_the_finding_counts_it(
    counted_as: dict[str, int],
) -> None:
    """One de novo allele and one inherited: the finding may call the person either."""
    pair = [_linked(1, 2, carriers=_ONE), _linked(2, 1, carriers=_INH)]
    assert _shape(pair, carriers=counted_as, screened=1) == []


def test_origin_unknown_claims_no_origin_and_draws_only_on_the_total() -> None:
    """D10, a departure from spec 7.2's literal "per origin".

    A paper reporting "1 de novo, 1 inherited" overall with two variants of unstated
    origin is consistent; naming a per-variant origin would be invented. Two such
    variants against a total of one are still too many. A linked pair of unknown
    origin carries no per-origin lower bound either.
    """
    both = {"carriers": {"de_novo": 1, "inherited": 1, "origin_unknown": 0}, "screened": 2}
    assert _shape([_sv(1, carriers=_UNK), _sv(2, carriers=_UNK)], **both) == []
    one = {"carriers": _ONE, "screened": 1}
    issues = _shape([_sv(1, carriers=_UNK), _sv(2, carriers=_UNK)], **one)
    assert [(code, "carrier(s)" in message) for code, _, message in issues] == [("LIT023", True)]
    three = {"carriers": {"de_novo": 3, "inherited": 0, "origin_unknown": 0}, "screened": 3}
    pair = [
        _linked(1, 2, carriers=_UNK | {"origin_unknown": 3}),
        _linked(2, 1, carriers=_UNK | {"origin_unknown": 3}),
    ]
    assert _shape(pair, **three) == []


def test_trans_needs_two_entries_but_cis_may_span_three() -> None:
    """D9, and a full three-way `with` is the only way to link three entries."""
    everyone = {1: [2, 3], 2: [1, 3], 3: [1, 2]}
    cis = [_linked(i, 0, phase="cis", **{"with": w}) for i, w in everyone.items()]
    assert _shape(cis) == []
    trans = [_linked(i, 0, phase="trans", **{"with": w}) for i, w in everyone.items()]
    assert [code for code, _, _ in _shape(trans)] == ["LIT022"]


def test_a_linked_entry_must_state_a_phase_even_if_it_is_not_stated() -> None:
    silent = [_linked(1, 2, phase=None), _linked(2, 1, phase=None)]
    issues = _shape(silent)
    assert [(code, number) for code, number, _ in issues] == [("LIT022", 1), ("LIT022", 2)]
    assert all("needs a phase" in message for _, _, message in issues)
    # One entry silent and its partner stating a phase also disagree on the group's phase.
    mixed = _shape([_linked(1, 2, phase=None), _linked(2, 1)])
    assert sorted(message[:30] for _, _, message in mixed) == [
        "a variant with a with link nee",
        "variants 1, 2 are linked but d",
    ]
    stated = [_linked(1, 2, phase="not_stated"), _linked(2, 1, phase="not_stated")]
    assert _shape(stated) == []


def test_a_repeated_id_builds_no_false_group_and_the_duplicate_refuses_the_file() -> None:
    """With two entries sharing an id a link names no one entry, so no group is built and
    LIT023 is skipped; the repeated id's own LIT022 still refuses the file."""
    issues = _shape([_sv(1), _sv(1), _sv(2)], carriers=_ONE, screened=1)
    assert [(code, number) for code, number, _ in issues] == [("LIT022", 1)]
    assert "more than once" in issues[0][2]


def _pair(first: dict[str, int], second: dict[str, int]) -> list[dict[str, Any]]:
    return [_linked(1, 2, carriers=first), _linked(2, 1, carriers=second)]


_DN2 = {"de_novo": 2, "inherited": 0, "origin_unknown": 0}
_IH2 = {"de_novo": 0, "inherited": 2, "origin_unknown": 0}


@pytest.mark.parametrize(
    ("variants", "finding_carriers", "fragment"),
    [
        # D8: {de novo 2} in each of two linked entries, two people: both are de novo
        # on every allele, and the finding names none.
        (_pair(_DN2, _DN2), _IH2, "at least 2 de_novo carrier(s)"),
        (_pair(_IH2, _IH2), _DN2, "at least 2 inherited carrier(s)"),
        # The pigeonhole bound is exact: {de novo 2} with {de novo 1, inherited 1}
        # forces one person de novo on both alleles.
        (
            _pair(_DN2, {"de_novo": 1, "inherited": 1, "origin_unknown": 0}),
            _IH2,
            "at least 1 de_novo",
        ),
    ],
)
def test_a_with_group_that_forces_an_origin_the_finding_lacks_fires_lit023(
    variants: list[dict[str, Any]], finding_carriers: dict[str, int], fragment: str
) -> None:
    issues = _shape(variants, carriers=finding_carriers, screened=2)
    assert [(code, fragment in message) for code, _, message in issues] == [("LIT023", True)], (
        issues
    )
    # The same group is clean once the finding counts those people with that origin.
    ok = {"de_novo": 2, "inherited": 0, "origin_unknown": 0}
    ok = ok if "de_novo" in fragment else {"de_novo": 0, "inherited": 2, "origin_unknown": 0}
    assert _shape(variants, carriers=ok, screened=2) == []
    # A group of unstated origin forces nothing (D10).
    unknown = {"de_novo": 0, "inherited": 0, "origin_unknown": 3}
    assert (
        _shape(
            _pair(unknown, unknown),
            carriers={"de_novo": 3, "inherited": 0, "origin_unknown": 0},
            screened=3,
        )
        == []
    )


def test_a_groups_surplus_in_one_origin_never_offsets_another_entrys_claim() -> None:
    """The lower bound is floored at zero per group.

    A pair of inherited entries has a de novo bound of 0 - 1 = -1, floored to 0:
    unfloored, that -1 would hide one of the two unlinked de novo carriers against
    a finding that names a single de novo carrier, with every total still fitting.
    """
    variants = [*_pair(_INH, _INH), _sv(3, carriers=_ONE), _sv(4, carriers=_ONE)]
    issues = _shape(
        variants, carriers={"de_novo": 1, "inherited": 2, "origin_unknown": 0}, screened=4
    )
    assert [(code, "at least 2 de_novo" in message) for code, _, message in issues] == [
        ("LIT023", True)
    ]


_PROTEIN = {"placement": "protein", "hgvs_c": None, "grch38": None, "normalisation": None}


@pytest.mark.parametrize("hgvs_p", ["NP_852259.1:p.?", "p.0?", "p.=", "p.0", "NP_852259.1:p.(=)"])
def test_a_protein_variant_naming_no_residue_is_lit022_and_the_builder_still_refuses_it(
    hgvs_p: str,
) -> None:
    from chd_atlas.build.gene_literature import variant_table

    issues = _shape([_sv(1, **_PROTEIN, hgvs_p=hgvs_p)])
    assert [(code, vid) for code, vid, _ in issues] == [("LIT022", 1)]
    assert "unplaced" in issues[0][2]
    published = _finding([_sv(1, **_PROTEIN, hgvs_p=hgvs_p)]).model_dump(mode="json", by_alias=True)
    with pytest.raises(ValueError):
        variant_table([published])


def test_an_unreadable_protein_name_is_not_also_lit022() -> None:
    assert _shape([_sv(1, **_PROTEIN, hgvs_p="not hgvs")]) == []


def test_two_variants_with_one_row_key_in_a_finding_are_lit022() -> None:
    mane = _shape([_sv(1), _sv(2, hgvs_c="NM_181486.4:c.836C>T", grch38="12:114366301:G:A")])
    assert [(c, v) for c, v, _ in mane] == [("LIT022", 2)]
    protein = _shape(
        [
            _sv(1, **_PROTEIN, hgvs_p="NP_852259.1:p.(Arg279Ter)"),
            _sv(2, **_PROTEIN, hgvs_p="p.Arg279Ter"),
        ]
    )
    assert [(c, v) for c, v, _ in protein] == [("LIT022", 2)]
    different = _shape([_sv(1), _sv(2, grch38="12:114366399:G:A")])
    assert different == []


def test_a_wrong_protein_accession_names_the_fix_for_its_placement(frame: ManeFrame) -> None:
    """A `mane` variant is re-normalised; a `protein` or `unplaced` one is never normalised."""
    wrong = "NP_852259.2:p.(Arg279Ter)"
    ((_, mane),) = reference_issues(_variant(**_with(hgvs_p=wrong)), frame[TBX5])
    assert "re-normalise" in mane
    for placement in ("protein", "unplaced"):
        issues = reference_issues(_variant(placement=placement, hgvs_p=wrong), frame[TBX5])
        assert [code for code, _ in issues] == ["LIT016"]
        assert "write hgvs_p on NP_852259.1" in issues[0][1]
        assert "re-normalise" not in issues[0][1]


def test_one_variant_alone_and_with_a_partner_is_two_entries_but_not_otherwise() -> None:
    """Spec 7.2 (ruling D16): the same name twice is allowed only in different genotype
    contexts, i.e. when the entries' `with` partner sets differ."""
    alone = _sv(1, grch38="12:114366312:G:A")
    partnered = _linked(2, 3, grch38="12:114366312:G:A")
    partner = _linked(3, 2, grch38="12:114366399:G:A")
    assert [i for i in _shape([alone, partnered, partner]) if i[0] == "LIT022"] == []
    # Both unlinked: refused.
    both_alone = _shape([alone, _sv(2, grch38="12:114366312:G:A")])
    assert [(c, v) for c, v, _ in both_alone] == [("LIT022", 2)]
    assert "genotype context" in both_alone[0][2]
    # Same partner set: refused, whatever order the partners are listed in.
    twice = _shape(
        [
            _linked(1, 3, grch38="12:114366312:G:A", **{"with": [3, 4]}),
            _linked(2, 3, grch38="12:114366312:G:A", **{"with": [4, 3]}),
            _linked(3, 1, grch38="12:114366399:G:A", **{"with": [1, 2, 4]}),
            _linked(4, 1, grch38="12:114366398:G:A", **{"with": [1, 2, 3]}),
        ]
    )
    assert [(c, v) for c, v, m in twice if "repeats" in m] == [("LIT022", 2)]
    # Different partners of the same count are still different contexts.
    apart = _shape(
        [
            _linked(1, 3, grch38="12:114366312:G:A"),
            _linked(2, 4, grch38="12:114366312:G:A"),
            _linked(3, 1, grch38="12:114366399:G:A"),
            _linked(4, 2, grch38="12:114366398:G:A"),
        ]
    )
    assert [i for i in apart if "repeats" in i[2]] == []


def test_a_variant_naming_zero_carriers_is_lit022() -> None:
    zero = {"de_novo": 0, "inherited": 0, "origin_unknown": 0}
    issues = _shape([_sv(1, carriers=zero)])
    assert [(c, v) for c, v, _ in issues] == [("LIT022", 1)]
    assert "0 carriers" in issues[0][2]
    one = {"de_novo": 0, "inherited": 0, "origin_unknown": 1}
    assert _shape([_sv(1, carriers=one)]) == []


def _repeats(issues: list[tuple[str, int | None, str]]) -> list[int | None]:
    return [vid for code, vid, message in issues if code == "LIT022" and "repeats" in message]


def test_a_repeated_name_needs_a_different_genotype_context() -> None:
    """Ruling D16': a context is the partners' row keys plus the phase, not their ids.

    A: three contexts (alone, with one partner, with another) pass. B: alone twice.
    C: one pair entered twice under fresh ids. D: same partners, other phase.
    E: the same change twice in one with group, where the contexts differ (so only
    the group rule refuses it). F: twice in one cis group.
    """
    k = "12:114366312:G:A"
    a = [
        _sv(1, grch38=k),
        _linked(2, 4, grch38=k),
        _linked(3, 5, grch38=k),
        _linked(4, 2, grch38="12:114366399:G:A"),
        _linked(5, 3, grch38="12:114366398:G:A"),
    ]
    assert _repeats(_shape(a)) == []
    assert _repeats(_shape([_sv(1, grch38=k), _sv(2, grch38=k)])) == [2]
    other = "12:114366399:G:A"
    split = [
        _linked(1, 2, grch38=k),
        _linked(2, 1, grch38=other),
        _linked(3, 4, grch38=k),
        _linked(4, 3, grch38=other),
    ]
    assert _repeats(_shape(split)) == [3, 4]
    phased = [
        _linked(1, 2, grch38=k),
        _linked(2, 1, grch38=other),
        _linked(3, 4, grch38=k, phase="cis"),
        _linked(4, 3, grch38=other, phase="cis"),
    ]
    assert _repeats(_shape(phased)) == []
    # Same group, different contexts (variant 1 has one partner, variant 2 has two).
    grouped = [
        _linked(1, 2, grch38=k),
        _linked(2, 1, grch38=k, **{"with": [1, 3]}),
        _linked(3, 2, grch38=other),
    ]
    assert _repeats(_shape(grouped)) == [2]
    cis = [
        _linked(1, 2, grch38=k, phase="cis", **{"with": [2, 3]}),
        _linked(2, 1, grch38=k, phase="cis", **{"with": [1, 3]}),
        _linked(3, 1, grch38=other, phase="cis", **{"with": [1, 2]}),
    ]
    assert _repeats(_shape(cis)) == [2]


def test_a_partner_with_no_row_key_is_told_apart_by_what_the_paper_printed() -> None:
    """Two unplaced partners have no row key; their `as_reported` names the context."""
    k = "12:114366312:G:A"
    unplaced = {"placement": "unplaced", "grch38": None, "hgvs_c": None, "normalisation": None}
    unplaced["review_status"] = "pending_review"
    unplaced["review_note"] = "numbering differs"
    unplaced["review_settled_by"] = "authors"
    pair = [
        _linked(1, 3, grch38=k),
        _linked(2, 4, grch38=k),
        _linked(3, 1, as_reported="partner X", **unplaced),
        _linked(4, 2, as_reported="partner Y", **unplaced),
    ]
    assert _repeats(_shape(pair)) == []
    same = [*pair[:3], _linked(4, 2, as_reported="partner X", **unplaced)]
    assert _repeats(_shape(same)) == [2]


_UNCLEAR_PHASE = {
    "reason": "the paper does not say whether the changes share an allele",
    "settled_by": "full_text",
}


def _unclear_trio(entries: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Three variants linked pairwise, each `phase: unclear` with its own entry."""
    everyone = {1: [2, 3], 2: [1, 3], 3: [1, 2]}
    return [
        _linked(i, 0, phase="unclear", unclear={"phase": entries[i - 1]}, **{"with": w})
        for i, w in everyone.items()
    ]


def test_a_with_group_sharing_an_unclear_phase_validates_when_its_entries_agree() -> None:
    assert _shape(_unclear_trio([_UNCLEAR_PHASE] * 3)) == []


@pytest.mark.parametrize(
    "different",
    [
        {**_UNCLEAR_PHASE, "reason": "the authors give no haplotype data"},
        {**_UNCLEAR_PHASE, "settled_by": "authors"},
    ],
    ids=["reason", "settler"],
)
def test_a_with_group_whose_unclear_phase_entries_differ_is_refused(
    different: dict[str, Any],
) -> None:
    issues = _shape(_unclear_trio([_UNCLEAR_PHASE, _UNCLEAR_PHASE, different]))
    assert [(code, number) for code, number, _ in issues] == [("LIT022", 1)]
    assert "variants 1, 2, 3 share an unclear phase" in issues[0][2]
    assert "one statement in the paper is one entry" in issues[0][2]


def test_a_linked_variant_cannot_be_zygosity_unclear() -> None:
    unclear = {"zygosity": {"reason": "not stated in the table", "settled_by": "full_text"}}
    pair = [
        _linked(1, 2, zygosity="unclear", unclear=unclear),
        _linked(2, 1),
    ]
    issues = _shape(pair)
    assert [(code, number) for code, number, _ in issues] == [("LIT022", 1)]
    assert "must be heterozygous" in issues[0][2]


def test_the_repeat_message_tells_the_curator_to_write_zygosity_unclear() -> None:
    twice = _shape([_sv(1, grch38="12:114366312:G:A"), _sv(2, grch38="12:114366312:G:A")])
    messages = [m for c, _, m in twice if c == "LIT022" and "repeats" in m]
    assert len(messages) == 1
    assert (
        "with zygosity unclear, and its reason, when its carriers differ in zygosity" in messages[0]
    )
    assert "not_stated when its carriers differ" not in messages[0]
