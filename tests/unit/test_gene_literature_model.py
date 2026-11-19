# tests/unit/test_gene_literature_model.py
from __future__ import annotations

from datetime import date
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.gene_literature import (
    Finding,
    GeneLiteratureFile,
    Search,
    Synthesis,
    TriageEntry,
)
from chd_atlas.models.review import ReviewStatus
from chd_atlas.vocab import StudyType


def _rare(**changes: Any) -> dict[str, Any]:
    """A valid cohort_screen finding. Counts differ from one another (§4.30):
    screened 7, de_novo 2, inherited 1, origin_unknown 3 (sum 6)."""
    base: dict[str, Any] = {
        "id": "CHDA:LIT:0000001",
        "publication": "PMID:8988165",
        "design": "cohort_screen",
        "methods": ["exome"],
        "count_unit": "probands",
        "screened": 7,
        "carriers": {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
        "variant_classes": ["snv_indel"],
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "syndromic": "both",
        "shares_probands_with": "unknown",
        "locators": [{"kind": "in_text", "element": "Table 2"}],
        "summary": "Six of seven probands carry a variant.",
        "extracted_by": "agent:extract/test",
    }
    base.update(changes)
    return base


def _common(**changes: Any) -> dict[str, Any]:
    base = _rare(
        id="CHDA:LIT:0000002",
        design="common_variant_association",
        methods=["genotyping_array"],
        cases=11,
        controls=13,
    )
    for name in ("count_unit", "screened", "carriers", "variant_classes"):
        del base[name]
    base.update(changes)
    return base


def _file(**changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "gene": "HGNC:11604",
        "searches": [
            {
                "id": 1,
                "source": "pubmed",
                "query": "TBX5[tiab]",
                "run_on": "2026-11-02",
                "hits": 5,
                "new_hits": 4,
            }
        ],
        "triage": [],
        "findings": [],
    }
    base.update(changes)
    return base


def test_valid_findings_load() -> None:
    rare = Finding.model_validate(_rare())
    assert rare.carriers is not None and rare.carriers.total == 6
    Finding.model_validate(_common())
    Finding.model_validate(
        _rare(
            design="case_report",
            screened=1,
            carriers={"de_novo": 1, "inherited": 0, "origin_unknown": 0},
        )
    )
    Finding.model_validate(_rare(design="rare_variant_association", comparator="control_cohort"))
    # A paper that names no detection method (owner, 2026-11-03): recorded as such,
    # never as the nearest method, which would publish a claim the paper never made.
    assert Finding.model_validate(_rare(methods=["not_stated"])).methods == ["not_stated"]
    Finding.model_validate(_common(methods=["not_stated"]))


@pytest.mark.parametrize(
    ("changes", "field"),
    [
        # 3 + 3 + 3 = 9 > screened 7, though every arm alone is under 7: only the sum exceeds it
        ({"carriers": {"de_novo": 3, "inherited": 3, "origin_unknown": 3}}, "carriers"),
        # 0 + 0 + 8 = 8 > screened 7: one arm alone exceeds it
        ({"carriers": {"de_novo": 0, "inherited": 0, "origin_unknown": 8}}, "carriers"),
        ({"design": "case_report"}, "screened"),  # screened is 7, not 1
        ({"count_unit": None}, "count_unit"),
        ({"screened": None}, "screened"),
        ({"carriers": None}, "carriers"),
        ({"variant_classes": []}, "variant_classes"),
        ({"design": "family_segregation", "count_unit": None}, "count_unit"),
        ({"design": "rare_variant_association"}, "comparator"),
        (
            {"design": "rare_variant_association", "comparator": "none"},
            "cohort_screen",
        ),
        ({"comparator": "control_cohort"}, "comparator"),  # on a cohort_screen
        ({"methods": ["exome", "genotyping_array"]}, "genotyping"),
        ({"methods": ["targeted_genotyping"]}, "genotyping"),
        ({"methods": []}, "methods"),
        # "not stated" beside a named method contradicts itself.
        ({"methods": ["exome", "not_stated"]}, "not_stated"),
        ({"locators": []}, "locators"),
        ({"phenotypes": []}, "phenotypes"),
        ({"lesion_groups": []}, "lesion_groups"),
        ({"summary": "  "}, "summary"),
        ({"extracted_by": " "}, "extracted_by"),
        ({"screened": 0}, "screened"),
        ({"p_value": 0.0}, "p_value"),
        ({"p_value": 1.5}, "p_value"),
        ({"odds_ratio": 0}, "odds_ratio"),
        # `.inf` and `.nan` are YAML floats; an infinite odds ratio validated and then
        # crashed the build (allow_nan=False) until allow_inf_nan=False (final review I5).
        ({"odds_ratio": float("inf")}, "odds_ratio"),
        ({"odds_ratio": float("nan")}, "odds_ratio"),
        ({"p_value": float("inf")}, "p_value"),
        ({"p_value": float("nan")}, "p_value"),
        ({"approved_by": "A. Curator"}, "approved_by"),  # no approved_on
        ({"approved_on": "2026-11-02"}, "approved_by"),  # no approved_by
        ({"review_status": "not_accessible"}, "review_note"),
        ({"review_note": "no full text"}, "review_status"),
        ({"surprise": 1}, "surprise"),
    ],
)
def test_a_finding_breaking_a_design_rule_is_refused(changes: dict[str, Any], field: str) -> None:
    with pytest.raises(ValidationError, match=field):
        Finding.model_validate(_rare(**changes))


@pytest.mark.parametrize(
    ("changes", "field"),
    [
        ({"cases": None}, "cases"),
        ({"controls": None}, "cases and controls"),
        ({"cases": 0}, "cases"),
        ({"screened": 7}, "carriers or screened"),
        (
            {"carriers": {"de_novo": 1, "inherited": 0, "origin_unknown": 0}},
            "carriers or screened",
        ),
        ({"comparator": "control_cohort"}, "comparator"),
    ],
)
def test_a_common_variant_finding_breaking_a_rule_is_refused(
    changes: dict[str, Any], field: str
) -> None:
    with pytest.raises(ValidationError, match=field):
        Finding.model_validate(_common(**changes))


def test_carriers_at_exactly_screened_are_accepted() -> None:
    Finding.model_validate(_rare(carriers={"de_novo": 2, "inherited": 2, "origin_unknown": 3}))


@pytest.mark.parametrize(
    ("changes", "valid"),
    [
        ({"decision": "excluded"}, False),  # no reason
        ({"decision": "excluded", "reason": "review_article"}, True),
        ({"decision": "included", "reason": "review_article"}, False),
        ({"decision": "pending", "reason": "review_article"}, False),
        ({"decision": "included"}, True),
        ({"decision": "pending"}, True),
        ({"decision": "excluded", "reason": "same_data_as"}, False),
        ({"decision": "excluded", "reason": "same_data_as", "same_as": "PMID:2"}, True),
        ({"decision": "excluded", "reason": "not_chd", "same_as": "PMID:2"}, False),
        ({"decision": "excluded", "reason": "other"}, False),
        ({"decision": "excluded", "reason": "other", "note": "  "}, False),
        ({"decision": "excluded", "reason": "other", "note": "a thesis"}, True),
        ({"decision": "excluded", "reason": "full_text_not_accessible"}, False),
        ({"decision": "excluded", "reason": "full_text_not_accessible", "note": "no count"}, True),
        ({"decision": "included", "search": 0}, False),
    ],
)
def test_a_triage_entry_states_its_reason_only_when_excluded(
    changes: dict[str, Any], valid: bool
) -> None:
    raw = {"pmid": "PMID:1", "search": 1, **changes}
    if valid:
        TriageEntry.model_validate(raw)
    else:
        with pytest.raises(ValidationError):
            TriageEntry.model_validate(raw)


def test_a_search_cannot_find_more_new_papers_than_papers() -> None:
    raw = _file()["searches"][0]
    Search.model_validate({**raw, "hits": 4, "new_hits": 4})
    with pytest.raises(ValidationError, match="new_hits"):
        Search.model_validate({**raw, "hits": 4, "new_hits": 5})
    with pytest.raises(ValidationError):
        Search.model_validate({**raw, "source": "europepmc"})
    with pytest.raises(ValidationError):
        Search.model_validate({**raw, "query": " "})
    with pytest.raises(ValidationError):
        Search.model_validate({**raw, "id": 0})


@pytest.mark.parametrize("who", ["agent:x", " Agent:x", "AGENT:skill/model", ""])
def test_an_agent_cannot_approve(who: str) -> None:
    with pytest.raises(ValidationError, match="approved_by"):
        Finding.model_validate(
            _rare(approved_by=who, approved_on="2026-11-02", approved_digest="sha256:x")
        )


def test_a_person_approves_and_an_agent_may_verify() -> None:
    finding = Finding.model_validate(
        _rare(
            approved_by="A. Curator",
            approved_on="2026-11-02",
            approved_digest="sha256:x",
            verified_by="agent:x",
        )
    )
    assert finding.approved_by == "A. Curator"
    assert finding.verified_by == "agent:x"


def test_an_unreadable_paper_is_kept_with_its_abstract() -> None:
    finding = Finding.model_validate(
        _rare(
            review_status="not_accessible",
            review_note="Full text not accessible; the abstract is all that was read.",
            locators=[{"kind": "in_text", "element": "Abstract"}],
        )
    )
    assert finding.review_status is not None
    assert finding.review_status.value == "not_accessible"
    assert finding.locators[0].element == "Abstract"  # type: ignore[union-attr]


@pytest.mark.parametrize(
    ("changes", "message"),
    [
        (
            {"findings": [_rare(), _rare()]},
            "duplicate finding ids: ['CHDA:LIT:0000001']",
        ),
        (
            {"searches": _file()["searches"] * 2},
            "duplicate search ids: [1]",
        ),
        ({"searches": []}, "searches"),
        ({"gene": "TBX5"}, "gene"),
    ],
)
def test_a_literature_file_refuses_what_it_cannot_resolve(
    changes: dict[str, Any], message: str
) -> None:
    with pytest.raises(ValidationError, match=message.replace("[", r"\[").replace("]", r"\]")):
        GeneLiteratureFile.model_validate(_file(**changes))


def test_a_literature_file_round_trips_its_parts() -> None:
    parsed = GeneLiteratureFile.model_validate(
        _file(
            triage=[{"pmid": "PMID:8988165", "search": 1, "decision": "included"}],
            findings=[_rare(), _common()],
        )
    )
    assert [f.id for f in parsed.findings] == ["CHDA:LIT:0000001", "CHDA:LIT:0000002"]


def _synthesis(**changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "sentences": [
            {"text": "TBX5 variants cause septal defects.", "cites": ["CHDA:LIT:0000001"]}
        ],
        "drafted_by": "agent:synthesis/test",
    }
    base.update(changes)
    return base


def test_a_literature_file_carries_one_optional_synthesis() -> None:
    assert GeneLiteratureFile.model_validate(_file()).synthesis is None
    loaded = GeneLiteratureFile.model_validate(_file(synthesis=_synthesis()))
    assert loaded.synthesis is not None
    assert loaded.synthesis.approved_by is None


_ONE = {"text": "A sentence.", "cites": ["CHDA:LIT:0000001"]}


@pytest.mark.parametrize(
    "changes",
    [
        {"sentences": []},
        {"sentences": [_ONE] * 6},
        {"sentences": [{"text": "A sentence.", "cites": []}]},
        {"sentences": [{"text": "A sentence.", "cites": ["CHDA:LIT:0000001", "CHDA:LIT:0000001"]}]},
        {"sentences": [{"text": "   ", "cites": ["CHDA:LIT:0000001"]}]},
        {"sentences": [{"text": "A sentence.", "cites": ["CHDA:AST:0000001", "PMID:1"]}]},
        {"approved_by": "A. Curator"},
        {"approved_on": "2026-11-05"},
        {"approved_digest": "sha256:abc"},
        {
            "approved_by": "agent:x",
            "approved_on": "2026-11-05",
            "approved_digest": "sha256:abc",
        },
        {"extra": 1},
    ],
)
def test_a_synthesis_breaking_a_rule_is_refused(changes: dict[str, Any]) -> None:
    with pytest.raises(ValidationError):
        Synthesis.model_validate(_synthesis(**changes))


def test_a_synthesis_may_cite_findings_and_assertions_and_have_five_sentences() -> None:
    five = [
        {"text": f"Sentence {i}.", "cites": ["CHDA:LIT:0000001", "CHDA:AST:0000001"]}
        for i in range(5)
    ]
    assert len(Synthesis.model_validate(_synthesis(sentences=five)).sentences) == 5


def test_a_synthesis_is_drafted_by_an_agent_and_a_person_is_refused() -> None:
    with pytest.raises(
        ValidationError,
        match=r"synthesis: drafted_by must name the drafting agent \(agent:<skill>/<model>\); "
        "a person approves, an agent drafts",
    ):
        Synthesis.model_validate(_synthesis(drafted_by="A. Curator"))
    assert Synthesis.model_validate(_synthesis(drafted_by=" Agent:gene-synthesis/x")).drafted_by
    assert Synthesis.model_validate(_synthesis(drafted_by="agent:gene-synthesis/x"))


def test_a_synthesis_approval_names_a_person_with_the_findings_messages() -> None:
    with pytest.raises(
        ValidationError,
        match="synthesis: approved_by, approved_on and approved_digest are set together",
    ):
        Synthesis.model_validate(_synthesis(approved_by="A. Curator"))
    with pytest.raises(
        ValidationError, match="synthesis: approved_by names an agent; a person must approve"
    ):
        Synthesis.model_validate(
            _synthesis(approved_by=" Agent:x", approved_on="2026-11-05", approved_digest="sha256:a")
        )


def _variant_entry(**changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": 2,
        "as_reported": "p.Arg9Cys",
        "locator": {"kind": "in_text", "element": "Table 1"},
        "kind": "snv_indel",
        "placement": "protein",
        "with": [3],
    }
    base.update(changes)
    return base


def test_a_variant_round_trips_with_its_with_key_and_defaults_omitted() -> None:
    from chd_atlas.models.gene_literature import Variant

    entry = _variant_entry()
    variant = Variant.model_validate(entry)
    assert variant.with_ == [3]
    # exclude_defaults also drops the locator's own default `kind`, as it does
    # for every finding locator in the digest.
    expected = {**entry, "locator": {"element": "Table 1"}}
    assert variant.model_dump(mode="json", by_alias=True, exclude_defaults=True) == expected


def test_a_variant_review_note_without_a_status_is_refused() -> None:
    with pytest.raises(ValidationError, match="variant 2: a review_note needs a review_status"):
        Finding.model_validate(_rare(variants=[_variant_entry(review_note="open")]))
    with pytest.raises(ValidationError, match="variant 2: review_status"):
        Finding.model_validate(_rare(variants=[_variant_entry(review_status="pending_review")]))


def test_duplicate_normalisation_ids_are_refused() -> None:
    run = {
        "tool": "variantvalidator",
        "version": "v",
        "hgvs_version": "h",
        "vvta": "t",
        "vvdb": "d",
        "seqrepo": "s",
        "run_on": "2026-11-17",
    }
    with pytest.raises(ValidationError, match="duplicate normalisation ids"):
        GeneLiteratureFile.model_validate(
            _file(normalisations=[{"id": 1, **run}, {"id": 1, **run}])
        )


def test_normalisations_sit_immediately_after_searches_in_the_file() -> None:
    names = list(GeneLiteratureFile.model_fields)
    assert names[names.index("searches") + 1] == "normalisations"


def test_with_findings_updated_places_a_whole_variants_list_after_variant_classes() -> None:
    from io import StringIO

    from ruamel.yaml import YAML

    from chd_atlas.gene_literature_io import dump_bytes, load_round_trip, with_findings_updated

    raw_file = _file(findings=[_rare()])
    data = load_round_trip(dump_bytes(raw_file))
    out = with_findings_updated(
        dump_bytes(data), {"CHDA:LIT:0000001": {"variants": [_variant_entry()]}}
    )
    keys = list(YAML().load(StringIO(out.decode()))["findings"][0])
    assert keys.index("variants") == keys.index("variant_classes") + 1


@pytest.mark.parametrize("record", ["finding", "variant"])
@pytest.mark.parametrize(
    ("mark", "message"),
    [
        ({"review_settled_by": "authors"}, "review_settled_by needs a review_status"),
        (
            {"review_status": ReviewStatus.PENDING_REVIEW, "review_note": "open"},
            "pending_review needs a review_settled_by",
        ),
        (
            {
                "review_status": "not_accessible",
                "review_note": "no full text",
                "review_settled_by": "full_text",
            },
            "the full text settles it",
        ),
    ],
)
def test_a_pending_review_needs_review_settled_by_and_not_accessible_refuses_it(
    record: str, mark: dict[str, Any], message: str
) -> None:
    if record == "finding":
        with pytest.raises(ValidationError, match=message):
            Finding.model_validate(_rare(**mark))
    else:
        with pytest.raises(ValidationError, match=message):
            Finding.model_validate(_rare(variants=[_variant_entry(**mark)]))
    # the paired marks load
    ok = (
        {
            "review_status": ReviewStatus.PENDING_REVIEW,
            "review_note": "open",
            "review_settled_by": "authors",
        }
        if record == "finding"
        else {}
    )
    if record == "finding":
        assert Finding.model_validate(_rare(**ok)).review_settled_by is not None
    else:
        paired = {
            "review_status": "not_accessible",
            "review_note": "no full text",
        }
        Finding.model_validate(_rare(variants=[_variant_entry(**paired)]))


@pytest.mark.parametrize("field", ["zygosity", "phase"])
def test_an_unclear_value_and_its_entry_come_together(field: str) -> None:
    entry = {"reason": "Table 2 and the text disagree.", "settled_by": "atlas_model"}
    both = _variant_entry(**{field: "unclear", "unclear": {field: entry}})
    loaded = Finding.model_validate(_rare(variants=[both])).variants[0]
    assert loaded.unclear[field].settled_by.value == "atlas_model"
    with pytest.raises(ValidationError, match=f"{field} is unclear but has no unclear entry"):
        Finding.model_validate(_rare(variants=[_variant_entry(**{field: "unclear"})]))
    with pytest.raises(ValidationError, match=f"unclear entry for {field}, which is not unclear"):
        Finding.model_validate(_rare(variants=[_variant_entry(unclear={field: entry})]))
    blank = {**entry, "reason": "  "}
    with pytest.raises(ValidationError):
        Finding.model_validate(
            _rare(variants=[_variant_entry(**{field: "unclear", "unclear": {field: blank}})])
        )


def _shown_pair() -> tuple[Finding, Any]:
    from chd_atlas.models.literature import Publication

    unclear = {"reason": "Table 2 and the text disagree.", "settled_by": "atlas_model"}
    finding = Finding.model_validate(
        _rare(
            variants=[_variant_entry(phase="unclear", unclear={"phase": unclear}, **{"with": []})]
        )
    )
    paper = Publication.model_validate(
        {
            "id": "PMID:8988165",
            "title": "T",
            "journal": "J",
            "year": 1997,
            "authors": ["A A"],
            "study_type": "case_report",
        }
    )
    return finding, paper


def _edited(edit: str) -> tuple[Finding, Any]:
    finding, paper = _shown_pair()
    variant = finding.variants[0]
    if edit == "finding review_note":
        finding = finding.model_copy(
            update={
                "review_status": ReviewStatus.PENDING_REVIEW,
                "review_note": "Text and table differ.",
            }
        )
    elif edit == "variant review_note":
        mark = {"review_status": ReviewStatus.PENDING_REVIEW, "review_note": "Name unchecked."}
        finding = finding.model_copy(update={"variants": [variant.model_copy(update=mark)]})
    elif edit == "variant unclear reason":
        entry = variant.unclear["phase"].model_copy(update={"reason": "Only the text says so."})
        changed = variant.model_copy(update={"unclear": {"phase": entry}})
        finding = finding.model_copy(update={"variants": [changed]})
    elif edit == "verified_by":
        finding = finding.model_copy(update={"verified_by": "agent:verify/blind"})
    elif edit == "publication study_type":
        paper = paper.model_copy(update={"study_type": StudyType.GWAS})
    elif edit == "publication own_lab":
        paper = paper.model_copy(update={"own_lab": True})
    elif edit == "summary":
        finding = finding.model_copy(update={"summary": "Another sentence."})
    return finding, paper


@pytest.mark.parametrize(
    "edit",
    [
        "finding review_note",
        "variant review_note",
        "variant unclear reason",
        "verified_by",
        "publication study_type",
        "publication own_lab",
        "summary",
    ],
)
def test_the_shown_digest_moves_with_every_field_the_page_shows(edit: str) -> None:
    from chd_atlas.models.gene_literature import shown_digest

    before = shown_digest(*_shown_pair())
    assert before.startswith("sha256:")
    assert shown_digest(*_edited(edit)) != before


def test_approving_does_not_move_the_shown_digest_and_a_missing_paper_does() -> None:
    from chd_atlas.models.gene_literature import shown_digest

    finding, paper = _shown_pair()
    approved = finding.model_copy(
        update={
            "approved_by": "E.A.M",
            "approved_on": date(2026, 11, 19),
            "approved_digest": "sha256:a",
        }
    )
    assert shown_digest(approved, paper) == shown_digest(finding, paper)
    assert shown_digest(finding, None) != shown_digest(finding, paper)
