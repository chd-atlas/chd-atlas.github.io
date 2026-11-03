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

from chd_atlas.models.gene_literature import Finding, finding_digest
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
