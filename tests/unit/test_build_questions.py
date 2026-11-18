"""The questions a gene's published findings leave open (`build/questions.py`).

Fixtures are built with the real `Finding` model and dumped the way the bundle
publishes them, so a renamed field breaks here and not silently in a consumer.
Expected lists are literals.
"""

from __future__ import annotations

from typing import Any

from chd_atlas.build.questions import (
    open_questions,
    question_anchor,
    question_rows,
)
from chd_atlas.models.gene_literature import Finding


def _finding(number: int, pmid: str, **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": f"CHDA:LIT:{number:07d}",
        "publication": pmid,
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
        "summary": f"Summary {number}.",
        "extracted_by": "agent:extract/test",
    }
    base.update(changes)
    return Finding.model_validate(base).model_dump(mode="json", by_alias=True)


def _variant(number: int, **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": number,
        "as_reported": f"as printed {number}",
        "locator": {"kind": "in_text", "element": "Table 1"},
        "kind": "snv_indel",
        "placement": "protein",
    }
    base.update(changes)
    return base


def _entry(reason: str, settled_by: str) -> dict[str, str]:
    return {"reason": reason, "settled_by": settled_by}


def _q(kind: str, reason: str, settled_by: str, *records: dict[str, Any]) -> dict[str, Any]:
    return {"kind": kind, "reason": reason, "settled_by": settled_by, "records": list(records)}


def _rec(
    finding: int, pmid: str, anchor: str, variant: int | None = None, name: str | None = None
) -> dict[str, Any]:
    return {
        "finding": f"CHDA:LIT:{finding:07d}",
        "publication": pmid,
        "variant": variant,
        "variant_name": name,
        "anchor": anchor,
    }


def test_a_findings_mark_is_one_question_with_its_note_and_settler() -> None:
    plain = _finding(1, "PMID:1")
    marked = _finding(
        2,
        "PMID:2",
        review_status="pending_review",
        review_note="Case counts differ between Table 1 and Table 2.",
        review_settled_by="authors",
    )
    assert open_questions([plain, marked]) == [
        _q(
            "pending_review",
            "Case counts differ between Table 1 and Table 2.",
            "authors",
            _rec(2, "PMID:2", "q-chda-lit-0000002"),
        )
    ]


def test_not_accessible_is_settled_by_the_full_text() -> None:
    abstract_only = _finding(
        3,
        "PMID:3",
        review_status="not_accessible",
        review_note="Only the abstract was read.",
    )
    assert open_questions([abstract_only]) == [
        _q(
            "not_accessible",
            "Only the abstract was read.",
            "full_text",
            _rec(3, "PMID:3", "q-chda-lit-0000003"),
        )
    ]


def test_a_variant_with_both_fields_unclear_is_two_questions() -> None:
    finding = _finding(
        4,
        "PMID:4",
        variants=[
            _variant(
                5,
                hgvs_p="p.Arg9Cys",
                hgvs_c="c.25C>T",
                zygosity="unclear",
                phase="unclear",
                review_status="pending_review",
                review_note="Table 2 and the text disagree on the allele.",
                review_settled_by="expert",
                unclear={
                    "phase": _entry("Cis in two carriers, trans in one.", "atlas_model"),
                    "zygosity": _entry("Het in two, hom in one.", "authors"),
                },
            ),
            _variant(6, as_reported="p.Gly1Ala"),
        ],
    )
    name = "c.25C>T"
    assert open_questions([finding]) == [
        _q(
            "pending_review",
            "Table 2 and the text disagree on the allele.",
            "expert",
            _rec(4, "PMID:4", "q-chda-lit-0000004-v5", 5, name),
        ),
        _q(
            "zygosity",
            "Het in two, hom in one.",
            "authors",
            _rec(4, "PMID:4", "q-chda-lit-0000004-v5-zygosity", 5, name),
        ),
        _q(
            "phase",
            "Cis in two carriers, trans in one.",
            "atlas_model",
            _rec(4, "PMID:4", "q-chda-lit-0000004-v5-phase", 5, name),
        ),
    ]


def test_a_variant_is_named_by_hgvs_c_then_hgvs_p_then_as_reported() -> None:
    finding = _finding(
        5,
        "PMID:5",
        variants=[
            _variant(
                1,
                hgvs_c="c.1A>G",
                hgvs_p="p.Met1Val",
                zygosity="unclear",
                unclear={"zygosity": _entry("r", "authors")},
            ),
            _variant(
                2,
                hgvs_p="p.Gly2Ala",
                zygosity="unclear",
                unclear={"zygosity": _entry("r", "authors")},
            ),
            _variant(3, zygosity="unclear", unclear={"zygosity": _entry("r", "authors")}),
        ],
    )
    names = [r["variant_name"] for q in open_questions([finding]) for r in q["records"]]
    assert names == ["c.1A>G", "p.Gly2Ala", "as printed 3"]


def test_a_with_group_sharing_an_unclear_phase_is_one_question_naming_every_member() -> None:
    shared = _entry("Two in cis, trans to the third.", "expert")

    def member(number: int, partners: list[int]) -> dict[str, Any]:
        return _variant(
            number,
            hgvs_c=f"c.{number}A>G",
            zygosity="heterozygous",
            phase="unclear",
            unclear={"phase": shared},
            **{"with": partners},
        )

    finding = _finding(
        6,
        "PMID:6",
        variants=[
            member(3, [1, 2]),
            _variant(
                4,
                hgvs_c="c.4A>G",
                zygosity="unclear",
                unclear={"zygosity": _entry("Het in two, hom in one.", "authors")},
            ),
            member(2, [1, 3]),
            member(1, [2, 3]),
        ],
    )
    expected = [
        _q(
            "phase",
            "Two in cis, trans to the third.",
            "expert",
            _rec(6, "PMID:6", "q-chda-lit-0000006-v1-phase", 1, "c.1A>G"),
            _rec(6, "PMID:6", "q-chda-lit-0000006-v2-phase", 2, "c.2A>G"),
            _rec(6, "PMID:6", "q-chda-lit-0000006-v3-phase", 3, "c.3A>G"),
        ),
        _q(
            "zygosity",
            "Het in two, hom in one.",
            "authors",
            _rec(6, "PMID:6", "q-chda-lit-0000006-v4-zygosity", 4, "c.4A>G"),
        ),
    ]
    assert open_questions([finding]) == expected


def test_identical_questions_share_a_row_and_different_ones_do_not() -> None:
    note = "Only the abstract was read."
    first = _finding(1, "PMID:1", review_status="not_accessible", review_note=note)
    second = _finding(2, "PMID:2", review_status="not_accessible", review_note=note)
    other_reason = _finding(
        3, "PMID:3", review_status="not_accessible", review_note="The PDF is scanned."
    )
    kind_differs = _finding(
        4,
        "PMID:4",
        review_status="pending_review",
        review_note=note,
        review_settled_by="full_text",
    )
    settler_differs = _finding(
        5,
        "PMID:5",
        review_status="pending_review",
        review_note=note,
        review_settled_by="authors",
    )
    rows = question_rows(
        open_questions([first, other_reason, second, kind_differs, settler_differs])
    )
    a1, a2, a3, a4, a5 = (f"q-chda-lit-000000{n}" for n in range(1, 6))
    assert rows == [
        _q("pending_review", note, "full_text", _rec(4, "PMID:4", a4)),
        _q(
            "not_accessible",
            note,
            "full_text",
            _rec(1, "PMID:1", a1),
            _rec(2, "PMID:2", a2),
        ),
        _q("not_accessible", "The PDF is scanned.", "full_text", _rec(3, "PMID:3", a3)),
        _q("pending_review", note, "authors", _rec(5, "PMID:5", a5)),
    ]


def test_rows_follow_full_text_authors_expert_then_a_change_to_the_atlas() -> None:
    def pending(number: int, settled_by: str) -> dict[str, Any]:
        return _finding(
            number,
            f"PMID:{number}",
            review_status="pending_review",
            review_note=f"Question {number}.",
            review_settled_by=settled_by,
        )

    variants = _finding(
        9,
        "PMID:9",
        variants=[
            _variant(
                1, hgvs_c="c.1A>G", phase="unclear", unclear={"phase": _entry("P.", "authors")}
            ),
            _variant(
                2,
                hgvs_c="c.2A>G",
                zygosity="unclear",
                unclear={"zygosity": _entry("Z.", "authors")},
            ),
        ],
    )
    findings = [
        pending(1, "atlas_model"),
        pending(2, "expert"),
        variants,
        pending(3, "authors"),
        _finding(4, "PMID:4", review_status="not_accessible", review_note="Abstract."),
        pending(5, "full_text"),
    ]
    order = [
        (row["settled_by"], row["kind"], row["reason"])
        for row in question_rows(open_questions(findings))
    ]
    assert order == [
        ("full_text", "pending_review", "Question 5."),
        ("full_text", "not_accessible", "Abstract."),
        ("authors", "zygosity", "Z."),
        ("authors", "phase", "P."),
        ("authors", "pending_review", "Question 3."),
        ("expert", "pending_review", "Question 2."),
        ("atlas_model", "pending_review", "Question 1."),
    ]


def test_anchors_are_shaped_and_unique_across_a_gene() -> None:
    assert question_anchor("CHDA:LIT:0000020") == "q-chda-lit-0000020"
    assert question_anchor("CHDA:LIT:0000020", 3) == "q-chda-lit-0000020-v3"
    assert question_anchor("CHDA:LIT:0000020", 3, "zygosity") == "q-chda-lit-0000020-v3-zygosity"
    finding = _finding(
        20,
        "PMID:20",
        review_status="pending_review",
        review_note="Whole.",
        review_settled_by="authors",
        variants=[
            _variant(
                3,
                zygosity="unclear",
                phase="unclear",
                review_status="pending_review",
                review_note="Own.",
                review_settled_by="authors",
                unclear={"zygosity": _entry("Z.", "authors"), "phase": _entry("P.", "authors")},
            )
        ],
    )
    anchors = [r["anchor"] for q in open_questions([finding]) for r in q["records"]]
    assert anchors == [
        "q-chda-lit-0000020",
        "q-chda-lit-0000020-v3",
        "q-chda-lit-0000020-v3-zygosity",
        "q-chda-lit-0000020-v3-phase",
    ]
    assert len(set(anchors)) == len(anchors)
