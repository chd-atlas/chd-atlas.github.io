# tests/unit/test_build_validity.py
"""The three published validity states, and the discordance that travels with them."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.build.validity import (
    GeneValidity,
    ValidityRecord,
    _admitting_clingen,
    admission_provenance,
    admitting_grade,
    agreeing_submitters,
    gene_validity,
    published_genes,
    uncurated,
)
from chd_atlas.corpus import load_curation
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.vocab import Classification, ValiditySource, ValidityState

_CLINGEN_COLUMNS = (
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "sop",
    "classification",
    "classification_date",
    "gcep",
    "report_url",
)
_GENCC_COLUMNS = (
    "sgc_id",
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "classification",
    "submitter",
    "submitted_on",
    "report_url",
)

# Real HGNC/MONDO ids and labels, measured 2026-08-03 against the committed
# mirrors -- not synthesised, so a fixture cannot disagree with the mirror it is
# meant to stand in for.
_TBX5 = "HGNC:11604"
_GATA6 = "HGNC:4174"
_GATA6_CHD = "MONDO:0100540"
_GATA6_DCM = "MONDO:0005021"
_KIF26A = "HGNC:20226"


def _clingen_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "gene": _GATA6,
        "gene_symbol": "GATA6",
        "disease": _GATA6_CHD,
        "disease_label": "GATA6-related congenital heart disease",
        "moi": "AD",
        "sop": "SOP9",
        "classification": "Definitive",
        "classification_date": "2023-11-21",
        "gcep": "Congenital Heart Disease",
        "report_url": "https://search.clinicalgenome.org/a",
    }
    row.update(overrides)
    return row


def _gencc_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "sgc_id": "SGC-000001",
        "gene": _GATA6,
        "gene_symbol": "GATA6",
        "disease": _GATA6_CHD,
        "disease_label": "GATA6-related congenital heart disease",
        "moi": "Autosomal dominant",
        "classification": "Definitive",
        "submitter": "Ambry Genetics",
        "submitted_on": "2020-01-01",
        "report_url": None,
    }
    row.update(overrides)
    return row


@pytest.fixture(scope="module")
def repository_validity() -> dict[str, GeneValidity]:
    """`gene_validity()` over the committed mirrors and scope file.

    Module-scoped: the GenCC mirror is 30,376 rows (2026-11-14), and one call
    cost 0.22s on its 30,410 rows (mean of 5 warmed runs, measured 2026-08-04)
    -- more than a third of what this whole module took to run then.
    """
    root = Path(__file__).parent.parent.parent
    corpus, _ = load_curation(root)
    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    assert clingen is not None and gencc is not None
    return gene_validity(clingen, gencc, in_scope={str(entry.id) for entry in corpus.chd_scope})


def _clingen_frame(*rows: dict[str, object]) -> pl.DataFrame:
    return pl.DataFrame(list(rows) or None, schema={c: pl.String for c in _CLINGEN_COLUMNS})


def _gencc_frame(*rows: dict[str, object]) -> pl.DataFrame:
    return pl.DataFrame(list(rows) or None, schema={c: pl.String for c in _GENCC_COLUMNS})


def test_a_clingen_record_marks_expert_curated_and_gencc_records_still_publish() -> None:
    """GenCC rows for a gene ClinGen has also curated are not hidden.

    The point of publishing both mirrors is to show what the whole field says,
    not to hide agreement behind whichever source is "more authoritative".
    """
    clingen = _clingen_frame(_clingen_row(gene=_TBX5, disease="MONDO:0007732", moi="AD"))
    gencc = _gencc_frame(_gencc_row(gene=_TBX5, disease="MONDO:0007732", moi="Autosomal dominant"))

    result = gene_validity(clingen, gencc, in_scope={"MONDO:0007732"})

    tbx5 = result[_TBX5]
    assert tbx5.state is ValidityState.EXPERT_CURATED
    assert {record.source for record in tbx5.records} == {
        ValiditySource.CLINGEN,
        ValiditySource.GENCC,
    }


def test_a_gene_with_only_gencc_records_is_submitter_curated() -> None:
    result = gene_validity(_clingen_frame(), _gencc_frame(_gencc_row()), in_scope={_GATA6_CHD})

    assert result[_GATA6].state is ValidityState.SUBMITTER_CURATED


def test_out_of_scope_disease_rows_are_dropped_entirely() -> None:
    """GATA6 is Definitive for GATA6-related CHD and Limited for dilated cardiomyopathy.

    Both current, curated by different ClinGen expert panels (measured
    2026-08-03 against the committed mirror). Only the first is CHD; admitting
    the second would make the gene read as contested *about CHD* when the two
    verdicts concern different diseases.
    """
    clingen = _clingen_frame(
        _clingen_row(disease=_GATA6_CHD, classification="Definitive", gcep="CHD"),
        _clingen_row(
            disease=_GATA6_DCM,
            disease_label="dilated cardiomyopathy",
            classification="Limited",
            gcep="Dilated Cardiomyopathy",
        ),
    )

    result = gene_validity(clingen, _gencc_frame(), in_scope={_GATA6_CHD})

    gata6 = result[_GATA6]
    assert [record.disease for record in gata6.records] == [_GATA6_CHD]
    assert gata6.has_source_discordance is False


def test_records_and_genes_come_back_sorted_by_content_not_encounter_order() -> None:
    """Sorted against literals, not by building twice and comparing.

    `PYTHONHASHSEED` is fixed for the life of one interpreter, so two calls in
    the same process would iterate any unsorted dict or set identically either
    way -- a same-process comparison can never catch a missing sort (measured
    directly for this exact failure mode elsewhere in this project: dropping a
    sort and building under seeds 0 and 12345 gave byte-identical output over
    the one-gene committed corpus). Pinning the exact published order is the
    only check that can fail when the sort is removed.

    Both mirrors are handed rows out of the order this test expects back:
    ClinGen's two rows are inserted GATA6-CHD-then-dilated-cardiomyopathy, which
    sorts the other way round by `disease`; GenCC's two rows are inserted
    Zeta-then-Alpha, which sorts the other way round by `submitter`. TBX5 is
    included too, so the returned dict's own key order -- HGNC:11604 before
    HGNC:4174 -- is pinned in the same assertion.
    """
    clingen = _clingen_frame(
        _clingen_row(disease=_GATA6_CHD, moi="AD", classification="Definitive"),
        _clingen_row(
            disease=_GATA6_DCM,
            disease_label="dilated cardiomyopathy",
            moi="AD",
            classification="Limited",
        ),
    )
    gencc = _gencc_frame(
        _gencc_row(sgc_id="SGC-1", submitter="Zeta Labs", classification="Strong"),
        _gencc_row(sgc_id="SGC-2", submitter="Alpha Diagnostics", classification="Moderate"),
        _gencc_row(sgc_id="SGC-3", gene=_TBX5, disease="MONDO:0007732", submitter="Ambry"),
    )

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD, _GATA6_DCM, "MONDO:0007732"})

    assert list(result) == [_TBX5, _GATA6]
    assert [(r.source.value, r.disease, r.moi, r.submitter) for r in result[_GATA6].records] == [
        ("clingen", _GATA6_DCM, "AD", None),
        ("clingen", _GATA6_CHD, "AD", None),
        ("gencc", _GATA6_CHD, "Autosomal dominant", "Alpha Diagnostics"),
        ("gencc", _GATA6_CHD, "Autosomal dominant", "Zeta Labs"),
    ]


def test_supportive_publishes_the_verbatim_term_with_no_mapped_classification() -> None:
    """GenCC's `Supportive` is a mapping exception, not a rung on the ladder.

    5,274 rows in the committed GenCC mirror carry it, all from Orphanet
    (measured 2026-08-03) -- a submitter that asserts an association without
    grading its strength. Coercing it onto `Classification` would invent a
    confidence nobody stated, so the record must carry the verbatim term with
    `classification` left `None`.
    """
    result = gene_validity(
        _clingen_frame(),
        _gencc_frame(_gencc_row(classification="Supportive", submitter="Orphanet")),
        in_scope={_GATA6_CHD},
    )

    record = result[_GATA6].records[0]
    assert record.classification is None
    assert record.classification_term == "Supportive"


def test_an_unmapped_classification_raises_naming_the_term_and_the_gene() -> None:
    """A term either authority adds in a later release must reach a human.

    `raise`, never `assert`: `-O` strips `assert`, and this project's standing
    rule is that a bypassed gate must fail rather than publish a confidence
    nobody asserted.
    """
    clingen = _clingen_frame(_clingen_row(classification="Probably Fine"))

    with pytest.raises(ValueError, match="Probably Fine") as excinfo:
        gene_validity(clingen, _gencc_frame(), in_scope={_GATA6_CHD})
    assert _GATA6 in str(excinfo.value)


def test_discordance_fires_when_one_source_contests_and_another_supports() -> None:
    clingen = _clingen_frame(_clingen_row(classification="Disputed"))
    gencc = _gencc_frame(_gencc_row(classification="Strong", submitter="Ambry Genetics"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is True


def test_discordance_does_not_fire_for_differing_supportive_grades_across_sources() -> None:
    """Definitive and Strong differ in rank but agree the gene is supported.

    Worked example, measured 2026-08-03: ClinGen is Definitive for GATA6 /
    GATA6-related CHD, and PanelApp Australia (via GenCC) is Strong for the
    same gene-disease pair. `strongest()` already orders Definitive above
    Strong; this flag is not about rank, only about whether one source contests
    what another supports, and here neither does.
    """
    clingen = _clingen_frame(_clingen_row(classification="Definitive"))
    gencc = _gencc_frame(_gencc_row(classification="Strong", submitter="PanelApp Australia"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is False


def test_discordance_does_not_fire_when_supportive_is_the_only_record_on_one_side() -> None:
    """A `Supportive` record takes no side, so it cannot supply the "other" side.

    ClinGen contests the gene here; GenCC's only record is `Supportive`, which
    maps to no `Classification` at all. A submitter that declined to grade
    evidence has not thereby voted "supported", so this must not read as
    disagreement between the two sources -- only one of them has taken a
    position.
    """
    clingen = _clingen_frame(_clingen_row(classification="Disputed"))
    gencc = _gencc_frame(_gencc_row(classification="Supportive", submitter="Orphanet"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is False


def test_kif26a_multiple_gencc_submitters_disagreeing_on_grade_is_not_discordance() -> None:
    """Same-source multi-submitter disagreement on strength is not cross-source discordance.

    Worked example, measured 2026-08-03: KIF26A / MONDO:0859332 carries Ambry
    Genetics at both Limited and Strong, Baylor at Moderate, and Labcorp and
    PanelApp Australia at Strong -- five GenCC submissions, zero ClinGen
    records, and every classification supportive-level. `has_source_discordance`
    requires at least two distinct *sources*; GenCC disagreeing with itself
    five ways over is still one source.
    """
    disease = "MONDO:0859332"
    gencc = _gencc_frame(
        _gencc_row(
            sgc_id="SGC-102815",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Limited",
            submitter="Ambry Genetics",
        ),
        _gencc_row(
            sgc_id="SGC-104042",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Strong",
            submitter="Ambry Genetics",
        ),
        _gencc_row(
            sgc_id="SGC-125207",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Moderate",
            submitter="Baylor College of Medicine Research Center",
        ),
    )

    result = gene_validity(_clingen_frame(), gencc, in_scope={disease})

    kif26a = result[_KIF26A]
    assert kif26a.state is ValidityState.SUBMITTER_CURATED
    assert kif26a.has_source_discordance is False


def test_uncurated_is_the_shape_for_no_authority_and_a_fresh_instance_each_call() -> None:
    """A function, not a shared constant, so no caller can mutate a shared instance."""
    first = uncurated()
    second = uncurated()

    assert first == GeneValidity(
        state=ValidityState.UNCURATED, records=(), has_source_discordance=False
    )
    assert first is not second


def _clingen(classification: Classification) -> ValidityRecord:
    return ValidityRecord(
        source=ValiditySource.CLINGEN,
        classification=classification,
        classification_term=classification.value.title(),
        disease="MONDO:0007732",
        disease_label="Holt-Oram syndrome",
        moi="AD",
        report_url=None,
        gcep="Congenital Heart Disease",
    )


def _gencc(classification: Classification | None, submitter: str) -> ValidityRecord:
    """One GenCC submission. `None` is `Supportive`, which maps to no rung.

    `None` is accepted because it is an ordinary GenCC value, not an edge case:
    `vocab.GENCC_CLASSIFICATIONS` maps `Supportive` to it, six genes in the
    committed mirrors carry nothing else, and a helper that could not express it
    is why no test distinguished "asserted without a grade" from "not asserted".
    """
    return ValidityRecord(
        source=ValiditySource.GENCC,
        classification=classification,
        classification_term="Supportive"
        if classification is None
        else classification.value.title(),
        disease="MONDO:0007732",
        disease_label="Holt-Oram syndrome",
        moi="Autosomal dominant",
        report_url=None,
        submitter=submitter,
    )


def _entry(*records: ValidityRecord) -> GeneValidity:
    return GeneValidity(
        records=records,
        state=(
            ValidityState.EXPERT_CURATED
            if any(record.source is ValiditySource.CLINGEN for record in records)
            else ValidityState.SUBMITTER_CURATED
        ),
        has_source_discordance=False,
    )


def test_the_two_warrants_that_publish_a_gene_and_the_four_that_do_not() -> None:
    """D21 as it stands after the 2026-08-06 widening, one corpus per clause.

    Each gene below is the *only* one that fails if its clause is removed, which
    is why they share a corpus rather than being seven tests: the rule is one
    predicate and a per-clause test cannot see a clause interacting with its
    neighbour.

    - **1** ClinGen `Definitive` -- the old rule, still admitted.
    - **2** ClinGen `Limited` -- the widening. `PUBLICATION_FLOOR` is `LIMITED`,
      so this is the weakest ClinGen record that admits, and lowering the floor
      one rung would admit `DISPUTED` with it.
    - **3** ClinGen `Disputed` alone. A rank *floor* excludes it by construction;
      an enumeration of admissible rungs is what a later edit forgets.
    - **4** one GenCC submitter at `Definitive`. `SUBMITTER_AGREEMENT` is 2, and
      a single submitter is one laboratory's panel-inclusion decision however
      strongly it is graded -- so a grade cannot substitute for a second body.
    - **5** two GenCC submitters at `Limited`. The second warrant, at its
      weakest: two bodies agreeing admits where one at `Definitive` does not.
    - **6** ClinGen `Disputed` plus two agreeing submitters. **The veto.**
      Without it two laboratories overrule a chartered panel. LEFTY2 is the live
      case (HGNC:3122, measured 2026-08-06: CHD GCEP `Disputed`, G2P and PanelApp
      Australia `Limited`), and it is why `published_genes` is not simply the
      union of two independent predicates.
    - **7** ClinGen submitting to GenCC beside one other submitter. **ClinGen
      must not vote twice** -- it is the source the first warrant already trusts
      alone, and measured 2026-08-06 it is GenCC's largest in-scope submitter at
      111 rows over 109 genes. Counting it here would turn every one of those
      into a one-submitter admission wearing two hats.
    """
    published = published_genes(
        {
            "HGNC:1": _entry(_clingen(Classification.DEFINITIVE)),
            "HGNC:2": _entry(_clingen(Classification.LIMITED)),
            "HGNC:3": _entry(_clingen(Classification.DISPUTED)),
            "HGNC:4": _entry(_gencc(Classification.DEFINITIVE, "Ambry Genetics")),
            "HGNC:5": _entry(
                _gencc(Classification.LIMITED, "G2P"),
                _gencc(Classification.LIMITED, "PanelApp Australia"),
            ),
            "HGNC:6": _entry(
                _clingen(Classification.DISPUTED),
                _gencc(Classification.LIMITED, "G2P"),
                _gencc(Classification.LIMITED, "PanelApp Australia"),
            ),
            "HGNC:7": _entry(
                _gencc(Classification.DEFINITIVE, "ClinGen"),
                _gencc(Classification.DEFINITIVE, "Ambry Genetics"),
            ),
        }
    )
    assert published == {"HGNC:1", "HGNC:2", "HGNC:5"}


def test_no_known_association_is_not_a_veto() -> None:
    """Spec D34, applied to the gate rather than to `has_conflicting_evidence`.

    A stated absence of association is not a contest -- it takes neither side of
    the discordance test, and treating it as one here would give a single "we
    looked and found nothing" the force of a refutation by an expert panel.

    Unreachable on the committed mirrors, which is exactly why it is pinned:
    measured 2026-08-06, every one of the 16 genes admitted on GenCC agreement
    has no in-scope ClinGen record at all, so the corpus cannot distinguish this
    from `_clingen_contests` returning True for rank 0. `NO_KNOWN_ASSOCIATION`
    sits at rank 0, below `DISPUTED`, so a veto written as a rank test rather
    than as membership of `CONTESTED` would swallow it silently.
    """
    published = published_genes(
        {
            "HGNC:1": _entry(
                _clingen(Classification.NO_KNOWN_ASSOCIATION),
                _gencc(Classification.LIMITED, "G2P"),
                _gencc(Classification.LIMITED, "PanelApp Australia"),
            )
        }
    )
    assert published == {"HGNC:1"}


def test_a_submitter_below_the_floor_is_not_a_second_body_agreeing() -> None:
    """`agreeing_submitters` counts submitters *at or above* `PUBLICATION_FLOOR`.

    **Measured 2026-08-06: deleting the floor check from `agreeing_submitters`
    survived all 796 tests.** Every fixture that exercised the second warrant
    used two submitters at exactly `Limited`, so "a submitter" and "a submitter
    at or above the floor" were the same set, and a gate admitting on two
    submissions of `No Known Disease Relationship` was indistinguishable from
    the real one. That is CLAUDE.md section 4.30: two figures equal in the
    corpus are one figure to every test.

    The dissenting rungs are the point. Without the floor, a gene two
    laboratories looked at and reported *no known relationship* for is published
    as though they had agreed it is a CHD gene -- the exact inversion, and it
    would publish silently because the gene appears with a full validity table.

    `Supportive` maps to `None` (`vocab.GENCC_CLASSIFICATIONS`) and is checked
    here too: it asserts an association without grading it, so it cannot clear a
    floor expressed as a rank, and a submitter offering one is not a second body
    agreeing for the purpose of admission.
    """
    below = _entry(
        _gencc(Classification.NO_KNOWN_ASSOCIATION, "Illumina"),
        _gencc(Classification.LIMITED, "G2P"),
    )
    assert agreeing_submitters(below) == ("G2P",)

    refuted = _entry(
        _gencc(Classification.REFUTED, "Illumina"),
        _gencc(Classification.DISPUTED, "Ambry Genetics"),
        _gencc(Classification.LIMITED, "G2P"),
    )
    assert agreeing_submitters(refuted) == ("G2P",)

    # ... and the gate that reads it refuses the gene, which is the claim that
    # matters to a reader.
    assert published_genes({"HGNC:1": below, "HGNC:2": refuted}) == set()


def test_the_strongest_clingen_record_admits_when_a_gene_carries_several() -> None:
    """`_admitting_clingen` takes the strongest qualifying record, not any one.

    **Measured 2026-08-06: turning that `max` into `min` survived all 796
    tests.** No fixture and no gene in the committed mirrors carries two
    in-scope ClinGen records at different rungs -- measured, all 76 panel-graded
    genes have exactly one distinct label at their headline grade -- so the
    choice among several was unreachable and therefore unguarded.

    It is not hypothetical: `curation/chd_scope.yaml` is an editorial file a
    curator widens, and one added MONDO term can give a gene a second qualifying
    record without any code changing. Under `min` such a gene publishes the
    *weakest* thing a panel said about it, and `admitted_by.disease_label` names
    the disease of that weakest record -- so a gene ClinGen calls definitive for
    one CHD syndrome and limited for another would headline `limited`.
    """
    entry = _entry(_clingen(Classification.LIMITED), _clingen(Classification.DEFINITIVE))
    provenance = admission_provenance(entry)
    admitted = provenance["admitted_by"]
    assert isinstance(admitted, dict)
    assert admitted["classification"] == "definitive"
    assert admitting_grade(entry) is Classification.DEFINITIVE


def test_asserted_by_counts_institutions_not_records_and_never_a_dissent() -> None:
    """One name per institution, and only institutions that assert the gene.

    Two mutants this kills, both measured 2026-08-06 as surviving all 796 tests.

    **The dedup.** ClinGen submits to GenCC under its own name, so a gene with a
    panel record and a ClinGen GenCC submission is one institution twice.
    Counting `gcep` values instead overstates by exactly one per gene. The
    assertion that guarded this was `institutions.count("ClinGen") == 1` against
    a list built from a `set`, which cannot return anything but 0 or 1 -- a
    vacuous check that no mutation could fail.

    **The dissent.** `asserted_by` is a count of *assertions*, and
    `no_known_association`, `disputed` and `refuted` are not assertions.
    Illumina recorded `No Known Disease Relationship` for GDF1 and was counted
    among the institutions asserting it until 2026-08-06 -- on the very gene the
    two-submitter rule exists because of.

    `Supportive` (`None`) **is** counted: an ungraded assertion is still an
    assertion, and excluding it would silently drop Orphanet.
    """
    entry = _entry(
        _clingen(Classification.DEFINITIVE),
        _gencc(Classification.DEFINITIVE, "ClinGen"),
        _gencc(Classification.LIMITED, "Ambry Genetics"),
        _gencc(None, "Orphanet"),
        _gencc(Classification.NO_KNOWN_ASSOCIATION, "Illumina"),
        _gencc(Classification.REFUTED, "Natera"),
    )
    asserted = admission_provenance(entry)["asserted_by"]
    assert isinstance(asserted, dict)
    assert asserted["institutions"] == ["Ambry Genetics", "ClinGen", "Orphanet"]
    assert asserted["count"] == 3


def test_admission_provenance_names_the_one_warrant_that_admitted_the_gene(
    repository_validity: dict[str, GeneValidity],
) -> None:
    """`admitted_by` is an object, not an array, because exactly one warrant admits.

    A gene with both is published on ClinGen's word and says so: the difference
    between "a chartered expert panel graded this" and "two laboratories put it
    on a panel" is the whole reason the field exists.

    Three properties in one test because they share the fixture and the same
    mutant kills them together: every key is present on every gene (an object
    whose shape varies is a trap for a consumer reading a field off one gene and
    expecting it on the next), `submitters` is `[]` rather than absent under a
    ClinGen warrant, and `asserted_by.count` deduplicates ClinGen against its own
    GenCC submissions.

    TBX1 is asserted against by name, not by index. It is the 22q11.2 gene, it
    has no in-scope ClinGen record at all, and it is what a non-commercial
    requirement would have dropped -- so a mutant that admitted it under
    `"clingen"` would publish the most recognisable gene on the site under an
    authority that never graded it.
    """
    tbx1 = admission_provenance(repository_validity["HGNC:11592"])
    assert tbx1["admitted_by"] == {
        "authority": "gencc_agreement",
        "classification": None,
        "disease": None,
        "disease_label": None,
        "panel": None,
        "submitters": ["Ambry Genetics", "Labcorp Genetics (formerly Invitae)"],
    }

    tbx5 = admission_provenance(repository_validity["HGNC:11604"])
    admitted = tbx5["admitted_by"]
    assert isinstance(admitted, dict)
    assert admitted["authority"] == "clingen"
    assert admitted["classification"] == "definitive"
    assert admitted["submitters"] == []
    assert set(admitted) == set(tbx1["admitted_by"] or {})

    asserted = tbx5["asserted_by"]
    assert isinstance(asserted, dict)
    institutions = asserted["institutions"]
    assert isinstance(institutions, list)
    assert asserted["count"] == len(institutions)
    assert institutions.count("ClinGen") == 1, "ClinGen counted once, not once per source"


def test_the_committed_mirrors_publish_one_hundred_nineteen_genes_across_four_rungs(
    repository_validity: dict[str, GeneValidity],
) -> None:
    """A regression pin on the live corpus, not a restatement of the rule above.

    **Pinned as a composition, not as a total.** 119 alone is one number that a
    scope edit and a mirror refresh can move in opposite directions and leave
    unchanged, and it says nothing about which warrant admitted what. The split
    is what `docs/data-api.md`, `index.html`'s counts and every chip on the site
    are written against: measured 2026-11-14, 100 genes carry a ClinGen record
    at or above the floor and 19 are admitted on GenCC agreement alone -- the
    total held while SEMA3E and ZFPM2 moved warrant, when ClinGen's refresh
    graded both Limited (98 and 21 from 2026-10-25). (92 = 76 + 16 until the
    Orphanet scope path admitted Noonan, Costello, Adams-Oliver and the other
    syndromes where CHD is frequent but not definitional, and NOTCH1's own
    ClinGen term; issue #8.)

    The population was 23 until 2026-08-06 and every gene in it was ClinGen
    `Definitive`; those 23 are still here, which is what makes this a widening
    rather than a replacement.
    """
    published = published_genes(repository_validity)
    assert len(published) == 119

    by_warrant = Counter(
        admission_provenance(repository_validity[gene])["admitted_by"]["authority"]  # type: ignore[index]
        for gene in published
    )
    assert by_warrant == {"clingen": 100, "gencc_agreement": 19}

    by_grade = Counter(
        record.classification
        for gene in published
        if (record := _admitting_clingen(repository_validity[gene])) is not None
    )
    assert by_grade == {
        Classification.DEFINITIVE: 42,
        Classification.LIMITED: 47,
        Classification.MODERATE: 10,
        Classification.STRONG: 1,
    }


def test_a_gencc_record_carries_its_submission_date_and_id_and_a_clingen_one_neither() -> None:
    """D51: the date a GenCC row renders, and the id that keys it, reach the record.

    Both were read from the mirror and dropped here until schema 2.16, so every
    GenCC row of every page rendered a blank date -- 470 of 591 on 2026-10-25.
    `submitted_on` is carried verbatim, in the mirror's own shape; normalising
    it for display is `pages._date`'s job.
    """
    clingen = _clingen_frame(_clingen_row())
    gencc = _gencc_row(sgc_id="SGC-102815", submitted_on="2018-05-31 20:04:32")

    records = gene_validity(clingen, _gencc_frame(gencc), in_scope={_GATA6_CHD})[_GATA6].records

    panel = next(r for r in records if r.source is ValiditySource.CLINGEN)
    submission = next(r for r in records if r.source is ValiditySource.GENCC)
    assert submission.submitted_on == "2018-05-31 20:04:32"
    assert submission.sgc_id == "SGC-102815"
    assert panel.submitted_on is None
    assert panel.sgc_id is None


def test_two_submissions_tied_on_every_other_sort_field_order_by_sgc_id() -> None:
    """PKD1L1 (HGNC:18053) is the published case: two Ambry Genetics submissions
    for one disease and inheritance, which `(source, disease, moi, submitter)`
    cannot order. Fed in both orders, so the result cannot be encounter order."""
    rows = [
        _gencc_row(sgc_id="SGC-104080", classification="Definitive"),
        _gencc_row(sgc_id="SGC-101696", classification="Moderate"),
    ]
    for feed in (rows, rows[::-1]):
        validity = gene_validity(_clingen_frame(), _gencc_frame(*feed), in_scope={_GATA6_CHD})
        assert [record.sgc_id for record in validity[_GATA6].records] == [
            "SGC-101696",
            "SGC-104080",
        ]


_HCVD = "Hereditary Cardiovascular Disease Gene Curation Expert Panel"
_HTAAD = "Heritable Thoracic Aortic Aneurysm and Dissection Gene Curation Expert Panel"


def test_two_panels_on_one_triple_publish_both_in_panel_order() -> None:
    """Two ClinGen panels' records on one (gene, disease, moi) both publish, ordered by panel.

    FOXE3 / MONDO:0019625 / AD is the real case (ClinGen's 2026-11-13 file:
    SOP4 Moderate 2016 by the HTAAD GCEP, SOP12 Limited 2026 by the HCVD GCEP;
    F28). FOXE3 is not published, so the pair is placed on GATA6's in-scope
    triple here. `(source, disease, moi, submitter, sgc_id)` ties on both --
    neither carries a submitter or an sgc_id -- so without `gcep` in
    `_sort_key` the order is the mirror's. Fed in both orders.

    The headline keeps `strongest()` (plan Decision 2): the older Moderate
    admits over the newer Limited. On equal grades the first record in sort
    order admits, which is the panel that sorts first.
    """
    rows = [
        _clingen_row(gcep=_HTAAD, sop="SOP4", classification="Moderate"),
        _clingen_row(gcep=_HCVD, sop="SOP12", classification="Limited"),
    ]
    for feed in (rows, rows[::-1]):
        entry = gene_validity(_clingen_frame(*feed), _gencc_frame(), in_scope={_GATA6_CHD})[_GATA6]
        assert [record.gcep for record in entry.records] == [_HCVD, _HTAAD]
        assert admitting_grade(entry) is Classification.MODERATE

    tied = [{**row, "classification": "Limited"} for row in rows]
    for feed in (tied, tied[::-1]):
        entry = gene_validity(_clingen_frame(*feed), _gencc_frame(), in_scope={_GATA6_CHD})[_GATA6]
        admitting = _admitting_clingen(entry)
        assert admitting is not None
        assert admitting.gcep == _HCVD
