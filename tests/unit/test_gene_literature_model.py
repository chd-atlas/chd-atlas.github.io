# tests/unit/test_gene_literature_model.py
from __future__ import annotations

from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.gene_literature import (
    Finding,
    GeneLiteratureFile,
    Search,
    TriageEntry,
)


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
