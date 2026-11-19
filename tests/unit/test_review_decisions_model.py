"""The review decisions file: what a review page downloads and a converter
writes, and what `chd-atlas literature apply-review` reads."""

from __future__ import annotations

import copy
from datetime import UTC
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.review_decisions import ReviewDecisions
from chd_atlas.schema_export import EXPORTED_MODELS

_DIGEST = "sha256:" + "ab" * 32


def _sample() -> dict[str, Any]:
    return {
        "format": 1,
        "gene": "HGNC:11634",
        "symbol": "TBX5",
        "page": {
            "built_on": "2026-11-19T13:00:00.000Z",
            "atlas_commit": None,
            "literature_sha256": _DIGEST,
        },
        "reviewer": {"name": "E.A.M", "role": "owner"},
        "session": {
            "started_at": "2026-11-19T13:40:00.000Z",
            "finished_at": None,
            "active_ms": 1200,
        },
        "papers_seen": ["PMID:111", "PMID:222"],
        "findings": [
            {
                "id": "CHDA:LIT:0000001",
                "shown_digest": _DIGEST,
                "decision": "approve",
                "active_ms": 10,
                "first_opened_at": "2026-11-19T13:45:20.500Z",
                "decided_at": "2026-11-19T13:45:27.860Z",
            },
            {
                "id": "CHDA:LIT:0000002",
                "shown_digest": _DIGEST,
                "decision": "fix",
                "note": "wrong denominator",
                "active_ms": 0,
                "decided_at": "2026-11-19T13:46:00.000Z",
            },
        ],
        "triage": {
            "decision": "approve",
            "flips": [{"pmid": "PMID:333", "to": "included", "note": "reports a variant"}],
            "active_ms": 5,
        },
    }


def _edit(path: list[str | int], value: object) -> dict[str, Any]:
    data = copy.deepcopy(_sample())
    node: Any = data
    for key in path[:-1]:
        node = node[key]
    if value is _DELETE:
        del node[path[-1]]
    else:
        node[path[-1]] = value
    return data


_DELETE = object()


def test_a_valid_sample_parses_with_the_pages_datetime_shape() -> None:
    parsed = ReviewDecisions.model_validate(_sample())
    assert parsed.papers_seen == ["PMID:111", "PMID:222"]
    assert parsed.findings[0].decided_at.utcoffset() == UTC.utcoffset(None)
    assert parsed.findings[0].decided_at.microsecond == 860000
    assert parsed.findings[0].first_opened_at is not None
    assert parsed.findings[1].first_opened_at is None
    assert parsed.findings[0].note == ""
    assert ReviewDecisions.model_validate({**_sample(), "triage": None, "papers_seen": []})


def _dup_id() -> dict[str, Any]:
    return _edit(["findings", 1, "id"], "CHDA:LIT:0000001")


@pytest.mark.parametrize(
    "data",
    [
        pytest.param(_dup_id(), id="duplicate-finding-id"),
        pytest.param(_edit(["findings", 1, "note"], "  "), id="fix-blank-note"),
        pytest.param(_edit(["findings", 1, "note"], _DELETE), id="fix-no-note"),
        pytest.param(_edit(["extra"], 1), id="unknown-top-level-key"),
        pytest.param(_edit(["findings", 0, "bogus"], 1), id="unknown-nested-key"),
        pytest.param(_edit(["format"], 2), id="format-2"),
        pytest.param(_edit(["findings", 0, "shown_digest"], "sha256:xyz"), id="bad-digest"),
        pytest.param(
            _edit(["findings", 0, "shown_digest"], "sha256:" + "ab" * 31), id="short-digest"
        ),
        pytest.param(
            _edit(["findings", 0, "shown_digest"], "sha256:" + "ab" * 33), id="long-digest"
        ),
        pytest.param(_edit(["page", "literature_sha256"], "sha256:" + "AB" * 32), id="upper-hex"),
        pytest.param(_edit(["findings", 0, "id"], "CHDA:LIT:1"), id="bad-finding-id"),
        pytest.param(_edit(["gene"], "TBX5"), id="bad-gene"),
        pytest.param(_edit(["papers_seen"], ["PMID:1", "PMID:1"]), id="duplicate-paper-seen"),
        pytest.param(_edit(["papers_seen"], ["12"]), id="bad-paper-seen"),
        pytest.param(_edit(["findings", 0, "active_ms"], -1), id="negative-ms"),
        pytest.param(_edit(["triage", "flips", 0, "note"], " "), id="blank-flip-note"),
        pytest.param(_edit(["triage", "flips", 0, "to"], "maybe"), id="bad-flip-target"),
        pytest.param(_edit(["triage", "decision"], "fix"), id="triage-fix-no-note"),
        pytest.param(
            _edit(["triage"], {"decision": "later", "note": " ", "active_ms": 1}),
            id="triage-later-blank-note",
        ),
        pytest.param(_edit(["reviewer", "role"], "intern"), id="bad-role"),
        pytest.param(_edit(["reviewer", "name"], " "), id="blank-reviewer"),
    ],
)
def test_malformed_files_are_refused(data: dict[str, Any]) -> None:
    with pytest.raises(ValidationError):
        ReviewDecisions.model_validate(data)


def test_later_with_a_note_is_accepted() -> None:
    data = _edit(["findings", 1, "decision"], "later")
    assert ReviewDecisions.model_validate(data).findings[1].decision == "later"


def test_the_schema_is_registered_and_committed() -> None:
    assert EXPORTED_MODELS["review_decisions"] is ReviewDecisions
    path = Path(__file__).parents[2] / "schemas" / "review_decisions.schema.json"
    assert path.is_file()
