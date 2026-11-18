# src/chd_atlas/models/review.py
"""An open question a curated record carries, stated rather than settled.

The owner's standing rule (answer of 2026-10-27, "from now on"): where a curated
record rests on a source that conflicts with the paper it cites, or on a source
that could not be read, the record says so. It carries a status naming which,
and a short note stating exactly what is unresolved. It never asserts either
side, never speculates about which is right, and is never dropped: dropping
the record would hide the question along with it.

**The vocabulary is wider than the rule that introduced it.** The same answer
marked two scope terms `pending_review` for a question that is not about
sources at all: whether each names a syndrome, left open for clinical review.
So the one definition, which `docs/data-api.md`, the manifest's changelog and
both curated files' headers repeat (the review of 2026-10-27 found the API
doc, the changelog and `chd_scope.yaml` defining it by a source conflict
alone, which is false for the two terms):

* `pending_review` -- a question about the record that this atlas leaves
  open. On the two cohorts marked today, the sources it rests on
  conflict; on the two scope terms, it is a clinical judgement.
* `not_accessible` -- a source the record rests on could not be read.

**One vocabulary for every record kind that carries one**, so "pending review"
means the same thing on a cohort as on a scope term, and a consumer reads one
token set. Today that is `Cohort`, `ScopeEntry` and the literature `Finding`
and `Variant`; each declares the two fields itself and calls `check_review`
from its own validator. A literature mark also names what would settle it
(`check_settled`, `gene_literature.py`): the authors, the full text, an expert
or a change to the atlas, so a literature question is not always one for an expert.

Absent means no open question has been recorded, not that the record was
checked and found clean.
"""

from __future__ import annotations

from enum import StrEnum


class ReviewStatus(StrEnum):
    """Why a curated record is marked as carrying an open question.

    `PENDING_REVIEW`: a question about the record that this atlas leaves open
    -- the sources it rests on conflict, as on the two cohorts
    the owner marked on 2026-10-27, or a clinical judgement is needed, as on
    the two scope terms marked the same day.
    `NOT_ACCESSIBLE`: a source it rests on could not be read, so what it says
    is unchecked. Neither is a verdict on the record.
    """

    PENDING_REVIEW = "pending_review"
    NOT_ACCESSIBLE = "not_accessible"


def check_review(record: str, status: ReviewStatus | None, note: str | None) -> None:
    """Refuse a status without its note, or a note without its status.

    A status alone tells a reader there is a question and not what it is, which
    invites the speculation the rule forbids; a note alone is an open question
    nothing marks as one, so it reads as part of the record. Raised from each
    model's validator, so the refusal arrives where the file loads (the build
    then refuses the corpus), and again from the two JSON emitters,
    `burden.build_cohorts` and `scope.build_scope`, so a record built without
    validation is refused before it publishes a status with a `null` note (a
    guard on a bypassed gate, like `render.review_status_note`'s for the
    pages). A blank note is refused before this, by the fields' own non-blank
    type.
    """
    if status is not None and note is None:
        raise ValueError(
            f"{record}: review_status {status.value!r} needs a review_note stating what is "
            f"unresolved"
        )
    if status is None and note is not None:
        raise ValueError(
            f"{record}: a review_note needs a review_status saying which open question it records"
        )
