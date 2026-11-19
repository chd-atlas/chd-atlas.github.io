# src/chd_atlas/models/review_decisions.py
"""The review decisions file: what a review page downloads, or a converter
writes, and what `chd-atlas literature apply-review` reads.

A decision names the finding it is about by id and by the digest of what the
reviewer was shown, so a finding edited after the page was built can be
detected. This module defines the shape only; it records no approval.

JSON Schema cannot express the cross-field rules here (a `fix` or `later`
needs a note; finding ids and papers seen are unique), so they live in model
validators and the exported schema checks well-formedness only.
"""

from __future__ import annotations

from datetime import datetime
from enum import StrEnum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import HgncId, LiteratureFindingId, Pmid

_NonBlankStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]
_Sha256 = Annotated[str, StringConstraints(pattern=r"^sha256:[0-9a-f]{64}$")]


class ReviewerRole(StrEnum):
    OWNER = "owner"
    SECOND_OPINION = "second_opinion"


class Decision(StrEnum):
    APPROVE = "approve"
    FIX = "fix"
    LATER = "later"


class _Strict(BaseModel):
    model_config = ConfigDict(extra="forbid")


class Reviewer(_Strict):
    name: _NonBlankStr
    role: ReviewerRole


class PageProvenance(_Strict):
    built_on: datetime
    atlas_commit: str | None
    literature_sha256: _Sha256


def _require_note(decision: Decision, note: str) -> None:
    """A `fix` or `later` says why; only `approve` may go without."""
    if decision is not Decision.APPROVE and not note.strip():
        raise ValueError(f"a '{decision.value}' decision needs a non-blank note")


class FindingDecision(_Strict):
    id: LiteratureFindingId
    shown_digest: _Sha256
    decision: Decision
    note: str = ""
    active_ms: int = Field(ge=0)
    first_opened_at: datetime | None = None
    decided_at: datetime

    @model_validator(mode="after")
    def _a_fix_or_later_says_why(self) -> FindingDecision:
        _require_note(self.decision, self.note)
        return self


class TriageFlip(_Strict):
    pmid: Pmid
    to: Literal["included", "excluded"]
    note: _NonBlankStr


class TriageReview(_Strict):
    decision: Decision
    note: str = ""
    flips: list[TriageFlip] = []
    active_ms: int = Field(ge=0)

    @model_validator(mode="after")
    def _a_fix_or_later_says_why(self) -> TriageReview:
        _require_note(self.decision, self.note)
        return self


class ReviewSession(_Strict):
    started_at: datetime | None
    finished_at: datetime | None
    active_ms: int = Field(ge=0)


class ReviewDecisions(_Strict):
    format: Literal[1]
    gene: HgncId
    symbol: _NonBlankStr
    page: PageProvenance
    reviewer: Reviewer
    session: ReviewSession
    papers_seen: list[Pmid] = []
    findings: list[FindingDecision]
    triage: TriageReview | None = None

    @model_validator(mode="after")
    def _no_repeats(self) -> ReviewDecisions:
        if repeated := duplicates([str(f.id) for f in self.findings]):
            raise ValueError(f"finding ids repeated: {sorted(repeated)}")
        if repeated := duplicates([str(p) for p in self.papers_seen]):
            raise ValueError(f"papers_seen repeated: {sorted(repeated)}")
        return self
