# tests/unit/test_build_paths.py
from typing import Final, NewType

import pytest
from pydantic import TypeAdapter

from chd_atlas import identifiers
from chd_atlas.build.paths import gene_bundle_path, gene_page_path, slug
from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    CohortId,
    ContrastId,
    DiseaseId,
    Doi,
    FunctionalId,
    HgncId,
    LiteratureFindingId,
    ModelPhenotypeId,
    ModificationId,
    OrphaId,
    PhenotypeId,
    Pmcid,
    Pmid,
    SequenceOntologyId,
    TaxonId,
    UniprotAccession,
)

# At least one ASCII value per grammar in `identifiers.py` except `Doi`, which
# slug provably cannot separate and which has its own test below, plus the two
# values closest to colliding: `HGNC:11604` against a `ContrastId` spelling the
# same characters. Each is paired with its type so the injectivity claim is
# pinned against strings the grammars actually admit rather than against
# plausible-looking inventions. ASCII because injectivity holds only there —
# `test_slug_cannot_separate_non_ascii_digits_so_injectivity_is_ascii_only`.
SAMPLE: Final = [
    (HgncId, "HGNC:11604"),
    (Pmid, "PMID:8988165"),
    (Pmcid, "PMC1234567"),
    (PhenotypeId, "HP:0001631"),
    (PhenotypeId, "MONDO:0007259"),
    (DiseaseId, "MONDO:0007732"),
    (ModelPhenotypeId, "MP:0010402"),
    (ModelPhenotypeId, "ZP:0000706"),
    (SequenceOntologyId, "SO:0001587"),
    (ModificationId, "MOD:00046"),
    (TaxonId, "NCBITaxon:9606"),
    (AssertionId, "CHDA:AST:0000001"),
    (FunctionalId, "CHDA:FUN:0000042"),
    (LiteratureFindingId, "CHDA:LIT:0000007"),
    (UniprotAccession, "Q99593"),
    (UniprotAccession, "Q99593-2"),
    (AccessionId, "PXD012345"),
    (AccessionId, "GSE123456"),
    (AccessionId, "E-MTAB-1234"),
    (AccessionId, "EGAS00001000123"),
    (ContrastId, "tof_vs_control"),
    (ContrastId, "hgnc_11604"),
    (CohortId, "ddd"),
    (CohortId, "ukbb"),
    (OrphaId, "ORPHA:648"),
]


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        # One colon row, not one per CURIE grammar: the colon is the same
        # character in `PMID:8988165` and `HP:0001631`, and rows for those two
        # detected nothing this row does not. Per-grammar coverage is
        # `SAMPLE`'s job, checked by the injectivity test below.
        ("HGNC:11604", "HGNC_11604"),
        # Zygosity.NOT_APPLICABLE is the value "n/a"; an unslugged slash would
        # silently become a directory separator in any path built from it.
        ("n/a", "n_a"),
        ("PXD012345", "PXD012345"),
        # `ContrastId` is built from underscore-joined words, and slug leaves an
        # underscore as itself rather than rewriting it to something else.
        ("tof_vs_control", "tof_vs_control"),
        # The dot is safe, so a DOI keeps its "10.1000" and only the slash moves.
        # `Doi` is the only grammar containing a dot at all.
        ("10.1000/ab.cd", "10.1000_ab.cd"),
        # The hyphen is safe as well: a UniProt isoform suffix and an
        # ArrayExpress accession both survive whole. `slug`'s claim that it
        # rewrites only the colon and the underscore depends on this, so a
        # mutation dropping "-" from the safe set has to fail somewhere.
        ("Q99593-2", "Q99593-2"),
        ("E-MTAB-1234", "E-MTAB-1234"),
        # No grammar admits a space, but `slug` takes `str` and the "n/a" row
        # above establishes that vocabulary values reach it. A space surviving
        # into a published path is one every consumer must percent-encode.
        ("atrial septal defect", "atrial_septal_defect"),
    ],
)
def test_slug_replaces_every_character_unsafe_in_a_path(value: str, expected: str) -> None:
    assert slug(value) == expected


def test_the_injectivity_sample_is_drawn_from_the_real_grammars() -> None:
    """`slug` claims a property of `identifiers.py`, so the sample must come from it.

    Left unchecked, a sample drifts into strings no grammar admits — "HGNC_1"
    reads like an identifier and is not one — and the test below would then pin
    a property of values that cannot occur.
    """
    for type_, value in SAMPLE:
        assert TypeAdapter(type_).validate_python(value) == value


def test_slug_is_injective_over_the_ascii_values_of_every_grammar_but_doi() -> None:
    """Two distinct ASCII identifiers must not collapse to one filename.

    Within ASCII the only characters slug rewrites are the colon and the
    underscore, and no two values differ solely by swapping the two: the
    underscore appears only in `ContrastId` and `CohortId`, neither of which
    admits a capital, while every grammar bearing a colon starts with one — and
    the hyphen and dot, which the parametrised cases above pin, are left alone
    entirely.

    `CohortId` was added on 2026-08-04 and this sentence named `ContrastId`
    alone until then. The argument survived the addition because the two share a
    grammar shape, but it had to be re-checked rather than assumed: a grammar
    admitting both an underscore and a capital would break the claim outright.

    Three things this does not cover, each with its own test below: `Doi`, any
    non-ASCII value of any grammar, and two values a case-insensitive filesystem
    stores as one file. A fourth is beyond it: a grammar added to
    `identifiers.py` with no value in `SAMPLE`, which is why
    `test_the_sample_covers_every_grammar_that_exists` exists.
    """
    values = [value for _, value in SAMPLE]
    assert len({slug(value) for value in values}) == len(values)


def test_slug_cannot_separate_dois_so_no_path_may_be_built_from_one() -> None:
    r"""The one grammar excluded outright, pinned so the exclusion stays true.

    `Doi` is `^10\.\d{4,9}/\S+$`: the suffix admits every non-whitespace
    character, so slash, colon and underscore all arrive here and all leave as
    "_". This fails if something later relies on slug to name a DOI uniquely.
    """
    colliding = ["10.1000/ab/cd", "10.1000/ab:cd", "10.1000/ab_cd"]
    for value in colliding:
        assert TypeAdapter(Doi).validate_python(value) == value
    assert len(set(colliding)) == 3
    assert len({slug(value) for value in colliding}) == 1


def test_slug_cannot_separate_non_ascii_digits_so_injectivity_is_ascii_only() -> None:
    r"""The second exception, and the one that reaches `HgncId` itself.

    `\d` is Unicode-aware, so `^HGNC:\d+$` accepts far more than "HGNC:" and
    ASCII digits. The two values below are both accepted by the project's own
    validator, are genuinely distinct, and share one bundle path — which makes
    the injectivity above an ASCII property rather than a total one.

    This costs no gene page today: no curated source emits a non-ASCII digit,
    and `Emitter._write` raises on the second write to a path rather than
    replacing the first, so a collision would fail the build. That is a backstop
    elsewhere, not a property of `slug`, which is why it is pinned here.
    """
    # ARABIC-INDIC DIGIT ONE and TWO (U+0661, U+0662). They are near-identical to
    # ASCII in some fonts; they are not typos and must not be "corrected".
    colliding = ["HGNC:1160١", "HGNC:1160٢"]
    for value in colliding:
        assert TypeAdapter(HgncId).validate_python(value) == value
    assert len(set(colliding)) == 2
    assert len({slug(value) for value in colliding}) == 1
    assert len({gene_bundle_path(value) for value in colliding}) == 1


def test_slug_separates_identifiers_that_a_case_insensitive_filesystem_would_not() -> None:
    """The third exception, and the only one that leaves the string domain.

    `slug` is injective over strings, so these two get distinct slugs and the
    injectivity test above is satisfied. APFS and NTFS store them as one file:
    the second write replaces the first while the manifest advertises both, so a
    URL serves another record's bytes under a checksum that cannot verify.

    Note what the pair demonstrates about the proof: the only thing separating
    these two values is case, which is exactly what those filesystems discard.
    Nothing builds both paths today — `gene_bundle_path` takes an `HgncId` — and
    `Emitter._write` carries a casefolded guard for when something does.
    """
    gene, contrast = "HGNC:11604", "hgnc_11604"
    assert TypeAdapter(HgncId).validate_python(gene) == gene
    assert TypeAdapter(ContrastId).validate_python(contrast) == contrast

    assert slug(gene) != slug(contrast)
    assert slug(gene).casefold() == slug(contrast).casefold()


def test_the_sample_covers_every_grammar_that_exists() -> None:
    """`SAMPLE`'s comment claims full coverage; this is what makes that true.

    The injectivity argument is only as good as the sample it is checked
    against, and a grammar added to `identifiers.py` would otherwise be waved
    through silently — the hand-maintained weakness `slug`'s docstring admits to.
    `Doi` is the deliberate omission and is pinned by its own test.
    """
    declared = {name for name, obj in vars(identifiers).items() if isinstance(obj, NewType)}
    covered = {type_.__name__ for type_, _ in SAMPLE}
    assert declared, "no NewType grammars found; the introspection above has rotted"
    assert declared - covered == {"Doi"}


def test_gene_bundle_path_is_relative_and_forward_slashed() -> None:
    assert gene_bundle_path(HgncId("HGNC:11604")) == "genes/HGNC_11604.json"


def test_a_gene_page_sits_beside_its_bundle() -> None:
    """One slug, two extensions, so the pair cannot drift.

    `Emitter._write`'s case guard compares casefolded paths, and these two
    differ only after the dot, so they never collide.
    """
    assert gene_page_path(HgncId("HGNC:11604")) == "genes/HGNC_11604.html"
    assert gene_bundle_path(HgncId("HGNC:11604")) == "genes/HGNC_11604.json"
