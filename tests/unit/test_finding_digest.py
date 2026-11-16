# tests/unit/test_finding_digest.py
"""An approval is bound to a digest of the finding it approved (LIT015).

Fixtures take their digest from `tests.unit._stamp.stamped`, which computes it
with `finding_digest`.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.gene_literature import (
    Finding,
    Synthesis,
    finding_digest,
    synthesis_digest,
)
from tests.unit._stamp import stamped
from tests.unit.test_runner import _codes_of, _literature_repo


def _raw(**changes: Any) -> dict[str, Any]:
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
        "approved_by": "A. Curator",
        "approved_on": "2026-11-05",
    }
    base.update(changes)
    return stamped(base)


def _digest(**changes: Any) -> str:
    return finding_digest(Finding.model_validate(_raw(**changes)))


def test_the_digest_covers_content_and_not_marks() -> None:
    base = _digest()
    assert base.startswith("sha256:") and len(base) == len("sha256:") + 64
    assert _digest(carriers={"de_novo": 3, "inherited": 1, "origin_unknown": 3}) != base
    assert _digest(summary="Another sentence.") != base
    assert _digest(id="CHDA:LIT:0000002") != base
    assert _digest(locators=[{"kind": "in_text", "element": "Table 3"}]) != base
    assert _digest(verified_by="agent:verify/x") == base
    assert _digest(review_status="pending_review", review_note="check Table 2") == base
    assert _digest(approved_by="B. Curator", approved_on="2026-12-01") == base


@pytest.mark.parametrize(
    "changes",
    [
        {"approved_digest": None},
        {"approved_by": None, "approved_on": None},
        {"approved_on": None},
    ],
)
def test_the_three_approval_fields_are_set_together(changes: dict[str, Any]) -> None:
    with pytest.raises(ValidationError, match="approved_digest"):
        Finding.model_validate({**_raw(), **changes})


def _edited(tmp_path: Path, **finding: Any) -> Path:
    """A repository whose only finding was stamped, then edited without re-approving."""
    root = _literature_repo(tmp_path)
    path = root / "curation" / "literature" / "TBX5.yaml"
    raw = json.loads(path.read_text())
    raw["findings"][0].update(finding)
    path.write_text(json.dumps(raw))
    return root


@pytest.mark.parametrize(
    ("changes", "refused"),
    [
        ({}, False),
        ({"summary": "Edited after approval."}, True),
        ({"screened": 8}, True),
        ({"locators": [{"kind": "in_text", "element": "Table 3"}]}, True),
        ({"review_status": "pending_review", "review_note": "look again"}, False),
        ({"verified_by": "agent:verify/later"}, False),
        ({"approved_digest": "sha256:" + "0" * 64}, True),
    ],
)
def test_editing_approved_content_voids_the_approval(
    tmp_path: Path, changes: dict[str, Any], refused: bool
) -> None:
    assert ("LIT015" in _codes_of(_edited(tmp_path, **changes))) is refused


def test_an_approved_finding_with_no_digest_never_validates_clean(tmp_path: Path) -> None:
    root = _literature_repo(tmp_path)
    path = root / "curation" / "literature" / "TBX5.yaml"
    raw = json.loads(path.read_text())
    del raw["findings"][0]["approved_digest"]
    path.write_text(json.dumps(raw))
    codes = _codes_of(root)
    assert "SCHEMA001" in codes  # the model refuses it before LIT015 can see it


def test_a_field_added_to_the_model_at_its_default_does_not_move_the_digest() -> None:
    """An additive schema change must not void every standing approval.

    The subclass stands for a later release gaining an optional field: a finding
    that leaves it at its default keeps its digest, one that sets it does not.
    """

    class Grown(Finding):
        new_optional: str | None = None

    raw = _raw()
    assert finding_digest(Grown.model_validate(raw)) == finding_digest(Finding.model_validate(raw))
    changed = finding_digest(Grown.model_validate({**raw, "new_optional": "x"}))
    assert changed != finding_digest(Finding.model_validate(raw))


def _assertion(**changes: Any) -> LesionAssertion:
    base: dict[str, Any] = {
        "id": "CHDA:AST:0000001",
        "gene": "HGNC:11604",
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "isolated",
        "evidence": [
            {
                "publication": "PMID:8988165",
                "evidence_class": "genetic_case",
                "locator": {"kind": "in_text", "element": "Abstract"},
                "summary": "A nonsense variant.",
                "strength": "strong",
            }
        ],
        "curator": "A. Curator",
        "curated_on": "2026-07-31",
        "last_reviewed": "2026-07-31",
    }
    base.update(changes)
    return LesionAssertion.model_validate(base)


def _synthesis() -> Synthesis:
    return Synthesis.model_validate(
        {
            "sentences": [{"text": "Septal.", "cites": ["CHDA:LIT:0000001", "CHDA:AST:0000001"]}],
            "drafted_by": "agent:synthesis/test",
        }
    )


def test_the_synthesis_digest_moves_with_an_approved_finding_and_an_assertion() -> None:
    first = Finding.model_validate(_raw())
    second = Finding.model_validate(
        _raw(id="CHDA:LIT:0000002", summary="Nine of eleven probands carry a variant.", screened=11)
    )
    unapproved = Finding.model_validate(
        _raw(id="CHDA:LIT:0000003", summary="Not yet approved.", approved_by=None, approved_on=None)
    )
    assert unapproved.approved_by is None
    syn = _synthesis()
    assertion = _assertion()
    base = synthesis_digest(syn, [first, second], [assertion])
    assert base.startswith("sha256:")
    # reordering distinct approved findings changes nothing
    assert synthesis_digest(syn, [second, first], [assertion]) == base
    # an unapproved finding is not covered
    assert synthesis_digest(syn, [first, second, unapproved], [assertion]) == base
    # approving a further finding moves it
    third = Finding.model_validate(_raw(id="CHDA:LIT:0000003", summary="Approved now."))
    assert synthesis_digest(syn, [first, second, third], [assertion]) != base
    # editing an approved finding moves it
    edited = Finding.model_validate(_raw(summary="Changed."))
    assert synthesis_digest(syn, [edited, second], [assertion]) != base
    # editing an assertion's mechanism moves it; so does adding one
    assert (
        synthesis_digest(syn, [first, second], [_assertion(mechanism="dominant_negative")]) != base
    )
    other = _assertion(id="CHDA:AST:0000002")
    assert synthesis_digest(syn, [first, second], [assertion, other]) != base
    # the text and the cites are covered
    edited_text = Synthesis.model_validate(
        {
            "sentences": [{"text": "Valvar.", "cites": ["CHDA:LIT:0000001", "CHDA:AST:0000001"]}],
            "drafted_by": "agent:synthesis/test",
        }
    )
    assert synthesis_digest(edited_text, [first, second], [assertion]) != base
    # an approved finding's review mark is covered by the synthesis digest, and
    # not by the finding's own digest (a mark never voids a finding's approval)
    marked_raw = _raw(review_status="pending_review", review_note="The paper's table disagrees.")
    marked = Finding.model_validate(marked_raw)
    assert finding_digest(marked) == finding_digest(first)
    with_mark = synthesis_digest(syn, [marked, second], [assertion])
    assert with_mark != base
    other_note = Finding.model_validate({**marked_raw, "review_note": "A different note."})
    assert synthesis_digest(syn, [other_note, second], [assertion]) not in {base, with_mark}
    not_accessible = Finding.model_validate(
        {
            **marked_raw,
            "review_status": "not_accessible",
            "review_note": "The paper's table disagrees.",
        }
    )
    assert synthesis_digest(syn, [not_accessible, second], [assertion]) not in {base, with_mark}
    # the assertion part is order-independent
    assert synthesis_digest(syn, [first, second], [other, assertion]) == synthesis_digest(
        syn, [first, second], [assertion, other]
    )
    # the synthesis's own approval is not covered; its drafter is (the card prints it)
    approved = Synthesis.model_validate(
        {
            **syn.model_dump(mode="json"),
            "approved_by": "A. Curator",
            "approved_on": "2026-11-05",
            "approved_digest": base,
        }
    )
    assert synthesis_digest(approved, [first, second], [assertion]) == base
    redrafted = Synthesis.model_validate(
        {**syn.model_dump(mode="json"), "drafted_by": "agent:other/model"}
    )
    assert synthesis_digest(redrafted, [first, second], [assertion]) != base


@pytest.mark.parametrize(
    "change",
    [
        {"approved_by": "C. Third"},
        {"approved_on": "2026-11-20"},
        {"approved_digest": "sha256:" + "f" * 64},
    ],
)
def test_re_approving_a_finding_moves_the_synthesis_digest_and_not_the_findings_own(
    change: dict[str, str],
) -> None:
    """GS5: a re-approved finding withholds the synthesis, even with unchanged content.

    Each approval field changes alone, so a digest that dropped any one of the
    three leaves its row equal to the base. The two findings carry different
    approvers and dates (section 4.30), so a digest reading one finding's
    approval for both cannot pass either.
    """
    first = Finding.model_validate(_raw())
    second = Finding.model_validate(
        _raw(
            id="CHDA:LIT:0000002",
            summary="Another.",
            approved_by="B. Other",
            approved_on="2026-11-07",
        )
    )
    syn = _synthesis()
    base = synthesis_digest(syn, [first, second], [])
    for target, untouched in ((first, second), (second, first)):
        again = Finding.model_validate({**target.model_dump(mode="json"), **change})
        assert finding_digest(again) == finding_digest(target)
        assert synthesis_digest(syn, [again, untouched], []) != base
