"""The HGVS subset the variant checker reads, and the forms it refuses (LIT021)."""

import pytest

from chd_atlas.hgvs import (
    CPosition,
    CVariant,
    GenomeKey,
    HgvsUnsupported,
    PVariant,
    Region,
    parse_c,
    parse_grch38,
    parse_p,
    same_protein_change,
)


def _at(region: Region, base: int, offset: int = 0) -> CPosition:
    return CPosition(region, base, offset)


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("c.835C>T", CVariant(None, _at("cds", 835), _at("cds", 835), "sub", "C", "T")),
        ("c.130del", CVariant(None, _at("cds", 130), _at("cds", 130), "del", None, None)),
        (
            "c.92_98delGCCCGCC",
            CVariant(None, _at("cds", 92), _at("cds", 98), "del", "GCCCGCC", None),
        ),
        ("c.220_222dup", CVariant(None, _at("cds", 220), _at("cds", 222), "dup", None, None)),
        (
            "c.680_681insCTGAGAATAAT",
            CVariant(None, _at("cds", 680), _at("cds", 681), "ins", None, "CTGAGAATAAT"),
        ),
        ("c.1_3delinsATG", CVariant(None, _at("cds", 1), _at("cds", 3), "delins", None, "ATG")),
        (
            "c.-45G>A",
            CVariant(None, _at("five_prime", -45), _at("five_prime", -45), "sub", "G", "A"),
        ),
        (
            "c.*1101C>T",
            CVariant(None, _at("three_prime", 1101), _at("three_prime", 1101), "sub", "C", "T"),
        ),
        ("c.510+1G>A", CVariant(None, _at("cds", 510, 1), _at("cds", 510, 1), "sub", "G", "A")),
        (
            "c.-39+1G>A",
            CVariant(None, _at("five_prime", -39, 1), _at("five_prime", -39, 1), "sub", "G", "A"),
        ),
        (
            "NM_181486.4:c.663-1G>A",
            CVariant("NM_181486.4", _at("cds", 663, -1), _at("cds", 663, -1), "sub", "G", "A"),
        ),
        ("c.(?_-1)_(*1_?)del", None),
        ("c.835C>T;c.836A>G", None),
        ("c.[835C>T;836A>G]", None),
        ("c.835C>C", None),
        ("c.835_836C>T", None),
        ("c.835insA", None),
        ("c.0A>G", None),
        ("c.510+0G>A", None),
        ("c.835_840inv", None),
        ("c.835_836delCGinsTT", None),
    ],
)
def test_parse_c_reads_each_supported_form_and_refuses_the_rest(
    text: str, expected: CVariant | None
) -> None:
    if expected is None:
        with pytest.raises(HgvsUnsupported):
            parse_c(text)
    else:
        assert parse_c(text) == expected


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("NP_852259.1:p.(Arg279Ter)", PVariant("NP_852259.1", (("R", 279),))),
        ("p.R279X", PVariant(None, (("R", 279),))),
        ("p.Ter519Gln", PVariant(None, (("*", 519),))),
        ("p.Gln70dup", PVariant(None, (("Q", 70),))),
        ("p.(Arg279GlyfsTer5)", PVariant(None, (("R", 279),))),
        ("p.Lys5_Gly8delinsTrp", PVariant(None, (("K", 5), ("G", 8)))),
        ("p.Lys2_Gly3insArg", PVariant(None, (("K", 2), ("G", 3)))),
        ("NP_852259.1:p.?", PVariant("NP_852259.1", ())),
        ("p.Lys2_Gly4insArg", None),
        ("p.Gly8_Lys5del", None),
        ("p.Arg279_Gly280Ter", None),
        ("p.Ter519GlnextTer?", None),
        ("p.Arg0Ter", None),
        ("R279X", None),
    ],
)
def test_parse_p_names_the_reference_residues_and_refuses_the_rest(
    text: str, expected: PVariant | None
) -> None:
    if expected is None:
        with pytest.raises(HgvsUnsupported):
            parse_p(text)
    else:
        assert parse_p(text) == expected


def test_parse_grch38_reads_a_vcf_key_and_refuses_the_rest() -> None:
    assert parse_grch38("12:114366312:G:A") == GenomeKey("12", 114366312, "G", "A")
    for text in ("chr12:114366312:G:A", "12:114366312:G:G", "12:0:G:A", "12:114366312:G:-"):
        with pytest.raises(HgvsUnsupported):
            parse_grch38(text)


def test_same_protein_change_equates_one_and_three_letter_forms() -> None:
    assert same_protein_change("p.R279X", "p.(Arg279Ter)")
    assert same_protein_change("p.R279X", "NP_852259.1:p.Arg279*")
    assert not same_protein_change("p.R279Q", "p.(Arg279Ter)")
    assert not same_protein_change("p.R279X", "p.(Arg280Ter)")
