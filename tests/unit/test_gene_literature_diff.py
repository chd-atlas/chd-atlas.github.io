# tests/unit/test_gene_literature_diff.py
"""`diff_findings`: a blind second extraction against the first.

Every pair below differs in exactly one thing, so a test that names the field
names the only field the diff may report. Values within a pair are chosen to
differ from every other field's (§4.30): screened 7 vs 9, de_novo 2 vs 3.
"""

from __future__ import annotations

from typing import Any

import pytest

from chd_atlas.gene_literature_diff import COMPARED_FIELDS, FindingDiff, diff_findings
from chd_atlas.models.gene_literature import Finding


def _finding(**changes: Any) -> Finding:
    base: dict[str, Any] = {
        "id": "CHDA:LIT:0000001",
        "publication": "PMID:8988165",
        "design": "cohort_screen",
        "methods": ["exome", "genome"],
        "count_unit": "probands",
        "screened": 7,
        "carriers": {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
        "variant_classes": ["snv_indel", "cnv"],
        "phenotypes": ["HP:0001631", "HP:0001629"],
        "lesion_groups": ["septal", "avsd"],
        "syndromic": "both",
        "cohorts": ["alpha_cohort", "beta_cohort"],
        "shares_probands_with": "unknown",
        "locators": [{"kind": "in_text", "element": "Table 2"}],
        "summary": "Six of seven probands carry a variant.",
        "extracted_by": "agent:extract/test",
    }
    base.update(changes)
    return Finding.model_validate(base)


def _pair(first: dict[str, Any], second: dict[str, Any]) -> FindingDiff:
    return diff_findings([_finding(**first)], [_finding(id="CHDA:LIT:0000009", **second)])


def test_the_compared_set_is_pinned() -> None:
    assert COMPARED_FIELDS == (
        "publication",
        "design",
        "methods",
        "comparator",
        "burden_study",
        "count_unit",
        "screened",
        "carriers",
        "variant_classes",
        "cases",
        "controls",
        "p_value",
        "odds_ratio",
        "phenotypes",
        "lesion_groups",
        "syndromic",
        "cohorts",
        "shares_probands_with",
    )


@pytest.mark.parametrize(
    ("field", "first", "second"),
    [
        ("methods", ["exome"], ["gene_panel"]),
        ("burden_study", None, "PMID:42230622"),
        ("count_unit", "probands", "families"),
        ("screened", 7, 9),
        (
            "carriers",
            {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
            {"de_novo": 3, "inherited": 1, "origin_unknown": 3},
        ),
        ("variant_classes", ["snv_indel"], ["snv_indel", "chromosomal"]),
        ("cases", None, 11),
        ("controls", None, 13),
        ("p_value", None, 0.02),
        ("odds_ratio", 1.5, 2.5),
        ("phenotypes", ["HP:0001631"], ["HP:0001631", "HP:0001642"]),
        ("lesion_groups", ["septal"], ["lvoto"]),
        ("syndromic", "isolated", "syndromic"),
        ("shares_probands_with", "unknown", ["CHDA:LIT:0000004"]),
    ],
)
def test_a_pair_differing_in_one_field_reports_that_field_and_only_it(
    field: str, first: Any, second: Any
) -> None:
    diff = _pair({field: first}, {field: second})

    assert diff.only_first == () and diff.only_second == ()
    reported_first = first if not isinstance(first, list) else sorted(first)
    reported_second = second if not isinstance(second, list) else sorted(second)
    assert diff.fields == {"CHDA:LIT:0000001": [(field, reported_first, reported_second)]}
    assert diff.agreed == ()


def test_the_comparator_differs_inside_a_rare_variant_association() -> None:
    rva = {"design": "rare_variant_association"}
    diff = _pair({**rva, "comparator": "control_cohort"}, {**rva, "comparator": "mutation_model"})
    assert diff.fields == {"CHDA:LIT:0000001": [("comparator", "control_cohort", "mutation_model")]}


def test_fields_the_curator_reads_are_not_compared() -> None:
    diff = _pair(
        {},
        {
            "summary": "Worded differently.",
            "locators": [{"kind": "in_text", "element": "Figure 9"}],
            "extracted_by": "agent:verify/other",
            "verified_by": "agent:verify/other",
            "review_status": "pending_review",
            "review_note": "a note",
        },
    )
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)


def test_lists_that_carry_no_order_are_compared_as_sets() -> None:
    diff = _pair(
        {"shares_probands_with": ["CHDA:LIT:0000004", "CHDA:LIT:0000005"]},
        {
            "methods": ["genome", "exome"],
            "variant_classes": ["cnv", "snv_indel"],
            "phenotypes": ["HP:0001629", "HP:0001631"],
            "lesion_groups": ["avsd", "septal"],
            "cohorts": ["beta_cohort", "alpha_cohort"],
            "shares_probands_with": ["CHDA:LIT:0000005", "CHDA:LIT:0000004"],
        },
    )
    # Matched despite the cohort order, and agreeing in every field.
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)


def test_a_finding_in_one_pass_only_lands_in_that_side() -> None:
    first = [
        _finding(),
        _finding(id="CHDA:LIT:0000002", publication="PMID:111"),
    ]
    second = [
        _finding(id="CHDA:LIT:0000008"),
        _finding(id="CHDA:LIT:0000007", publication="PMID:222"),
        # same paper and design, a different cohort set: another study
        _finding(id="CHDA:LIT:0000006", cohorts=["alpha_cohort"]),
    ]
    diff = diff_findings(first, second)
    assert diff.only_first == ("CHDA:LIT:0000002",)
    assert diff.only_second == ("CHDA:LIT:0000006", "CHDA:LIT:0000007")
    assert diff.agreed == ("CHDA:LIT:0000001",)


@pytest.mark.parametrize("side", ["first", "second"])
def test_a_key_carried_twice_in_either_pass_is_never_paired(side: str) -> None:
    twin_a = _finding(id="CHDA:LIT:0000003", screened=8)
    twin_b = _finding(id="CHDA:LIT:0000002", screened=9)
    single = _finding(id="CHDA:LIT:0000001")
    other = _finding(id="CHDA:LIT:0000008", publication="PMID:111")
    other_again = _finding(id="CHDA:LIT:0000009", publication="PMID:111")
    if side == "first":
        diff = diff_findings([twin_a, twin_b, other], [single, other_again])
        assert diff.only_first == ("CHDA:LIT:0000002", "CHDA:LIT:0000003")
        assert diff.only_second == ("CHDA:LIT:0000001",)
    else:
        diff = diff_findings([single, other], [twin_a, twin_b, other_again])
        assert diff.only_first == ("CHDA:LIT:0000001",)
        assert diff.only_second == ("CHDA:LIT:0000002", "CHDA:LIT:0000003")
    # The unambiguous paper-111 pair is still paired and agrees.
    assert diff.agreed == ("CHDA:LIT:0000008",)
    assert diff.fields == {}


def test_the_report_is_sorted() -> None:
    first = [
        _finding(id="CHDA:LIT:0000030", publication="PMID:3", screened=8, odds_ratio=1.5),
        _finding(id="CHDA:LIT:0000010", publication="PMID:1", screened=8),
    ]
    second = [
        _finding(id="CHDA:LIT:0000041", publication="PMID:1", screened=9),
        _finding(id="CHDA:LIT:0000042", publication="PMID:3", screened=9, odds_ratio=2.5),
    ]
    diff = diff_findings(first, second)
    assert list(diff.fields) == ["CHDA:LIT:0000010", "CHDA:LIT:0000030"]
    assert [name for name, _, _ in diff.fields["CHDA:LIT:0000030"]] == [
        "odds_ratio",
        "screened",
    ]  # model field order puts screened first; the report is alphabetical
