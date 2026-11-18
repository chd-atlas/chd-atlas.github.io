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
    """The 18 fields a second pass is compared on.

    chd-atlas/curator keeps a twin of this pin, `test_the_compared_set_is_pinned`
    in its `tests/test_diff.py`, because `chd-curator diff` compares the same set
    from its own copy. The two change together: editing one without the other
    makes `literature review --against` and `chd-curator diff` disagree about
    what "agreed" means.
    """
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
        "variants",
        "cases",
        "controls",
        "p_value",
        "odds_ratio",
        "phenotypes",
        "lesion_groups",
        "syndromic",
        "cohorts",
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
            "review_settled_by": "authors",
        },
    )
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)


def test_a_shares_probands_with_difference_alone_leaves_the_pair_agreed() -> None:
    # R26: the blind pass's ids are local, so the field cannot be compared.
    diff = _pair(
        {"shares_probands_with": ["CHDA:LIT:0000004"]}, {"shares_probands_with": "unknown"}
    )
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)
    assert "shares_probands_with" not in COMPARED_FIELDS


def test_lists_that_carry_no_order_are_compared_as_sets() -> None:
    diff = _pair(
        {},
        {
            "methods": ["genome", "exome"],
            "variant_classes": ["cnv", "snv_indel"],
            "phenotypes": ["HP:0001629", "HP:0001631"],
            "lesion_groups": ["avsd", "septal"],
            "cohorts": ["beta_cohort", "alpha_cohort"],
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


def test_the_same_paper_and_cohorts_under_another_design_is_another_study() -> None:
    first = [_finding()]
    second = [
        _finding(
            id="CHDA:LIT:0000009",
            design="case_report",
            screened=1,
            carriers={"de_novo": 1, "inherited": 0, "origin_unknown": 0},
        )
    ]
    diff = diff_findings(first, second)
    assert diff.only_first == ("CHDA:LIT:0000001",)
    assert diff.only_second == ("CHDA:LIT:0000009",)
    assert diff.agreed == () and diff.fields == {}


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


def _variant(vid: int, as_reported: str, **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": vid,
        "as_reported": as_reported,
        "locator": {"kind": "in_text", "element": f"Table {vid}"},
        "kind": "snv_indel",
        "placement": "protein",
        "with": [],
    }
    base.update(changes)
    return base


def test_variant_order_and_local_ids_do_not_count_as_a_difference() -> None:
    """All bookkeeping keys differ between the passes: id, with, normalisation,
    locator, review_status, review_note and review_settled_by."""
    first = [
        _variant(1, "p.Arg1Trp"),
        _variant(2, "p.Gly2Ser", zygosity="homozygous", **{"with": [1]}),
    ]
    second = [
        _variant(
            7,
            "p.Gly2Ser",
            zygosity="homozygous",
            normalisation=3,
            review_status="pending_review",
            review_note="the paper's table disagrees with its text",
            review_settled_by="authors",
            **{"with": [5]},
        ),
        _variant(5, "p.Arg1Trp", locator={"kind": "in_text", "element": "Figure 9"}),
    ]
    diff = _pair({"variants": first}, {"variants": second})
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)

    changed = _pair(
        {"variants": first}, {"variants": [second[1], {**second[0], "zygosity": "hemizygous"}]}
    )
    assert [name for name, _, _ in changed.fields["CHDA:LIT:0000001"]] == ["variants"]


def test_two_passes_differing_only_in_settlers_or_unclear_reasons_agree() -> None:
    reason = {"reason": "Table 2 and the text disagree.", "settled_by": "atlas_model"}
    first = [_variant(1, "p.Arg1Trp", zygosity="unclear", unclear={"zygosity": reason})]
    second = [
        _variant(
            4,
            "p.Arg1Trp",
            zygosity="unclear",
            unclear={"zygosity": {"reason": "Worded again.", "settled_by": "authors"}},
            review_status="pending_review",
            review_note="open",
            review_settled_by="expert",
        )
    ]
    diff = _pair(
        {
            "variants": first,
            "review_status": "pending_review",
            "review_note": "n",
            "review_settled_by": "authors",
        },
        {"variants": second},
    )
    assert diff.fields == {} and diff.agreed == ("CHDA:LIT:0000001",)
    other = _pair(
        {"variants": first},
        {"variants": [_variant(1, "p.Arg1Trp", zygosity="not_stated")]},
    )
    assert [name for name, _, _ in other.fields["CHDA:LIT:0000001"]] == ["variants"]
