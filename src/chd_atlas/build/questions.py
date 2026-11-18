"""A gene's open literature questions, gathered from its published findings.

Pure: no I/O, and nothing from `pages.py`. The input is the finding dicts as the
bundle publishes them (`Finding.model_dump(mode="json", by_alias=True)`), so the
`chd-atlas literature questions` command, the register page and the gene page's
"why" links all read one list and cannot disagree about what is open.

A question is one field of one record that the atlas leaves open, with the
curator's reason and what would settle it. Two shapes are not one-per-record:
a `with` group's shared unclear phase is one question (LIT022 makes every member
carry the identical entry), and `question_rows` merges equal questions so that,
for example, 22 abstract-only findings read as one row.

Every ordering here is over lists; nothing iterates a set (CLAUDE.md section 9).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Final, Literal, TypedDict

from chd_atlas.build.gene_literature import _partner_name

QuestionKind = Literal["zygosity", "phase", "pending_review", "not_accessible"]

QUESTION_LABELS: Final[dict[str, str]] = {
    "zygosity": "Zygosity unclear",
    "phase": "Phase unclear",
    "pending_review": "Pending review",
    "not_accessible": "Source not accessible",
}
SETTLED_LABELS: Final[dict[str, str]] = {
    "full_text": "the full text",
    "authors": "the authors",
    "expert": "an expert",
    "atlas_model": "a change to the atlas",
}
SETTLED_ORDER: Final = ("full_text", "authors", "expert", "atlas_model")
_KIND_ORDER: Final = ("zygosity", "phase", "pending_review", "not_accessible")


class QuestionRecord(TypedDict):
    finding: str
    publication: str
    variant: int | None
    variant_name: str | None
    anchor: str


class Question(TypedDict):
    kind: QuestionKind
    reason: str
    settled_by: str
    records: list[QuestionRecord]


def question_anchor(finding: str, variant: int | None = None, field: str | None = None) -> str:
    """`q-<finding id, lowercased, ':' as '-'>[-v<id>[-<field>]]`."""
    anchor = "q-" + finding.lower().replace(":", "-")
    if variant is not None:
        anchor += f"-v{variant}"
        if field is not None:
            anchor += f"-{field}"
    return anchor


def _mark(
    kind: QuestionKind,
    reason: str,
    settled_by: str,
    records: list[QuestionRecord],
) -> Question:
    return {"kind": kind, "reason": reason, "settled_by": settled_by, "records": records}


def _variant_record(
    finding: Mapping[str, Any], variant: Mapping[str, Any], field: str | None
) -> QuestionRecord:
    return {
        "finding": finding["id"],
        "publication": finding["publication"],
        "variant": variant["id"],
        "variant_name": _partner_name(variant),
        "anchor": question_anchor(finding["id"], variant["id"], field),
    }


def _review_mark(record: Mapping[str, Any], records: list[QuestionRecord]) -> Question | None:
    """The record's own mark as a question, or None. A status other than the two
    `ReviewStatus` values is refused rather than read as pending review, so a value
    added to the vocabulary later cannot publish under the wrong label."""
    status = record["review_status"]
    if status is None:
        return None
    if status == "not_accessible":
        # Settled by the full text by definition; the model refuses a settler here.
        return _mark("not_accessible", record["review_note"], "full_text", records)
    if status != "pending_review":
        raise ValueError(f"review_status {status!r} is neither pending_review nor not_accessible")
    return _mark("pending_review", record["review_note"], record["review_settled_by"], records)


def open_questions(findings: Sequence[Mapping[str, Any]]) -> list[Question]:
    """Every open question, in finding order.

    Within a finding: its own mark, then each variant by id (its own mark, then
    zygosity, then phase). A `with` group's unclear phase is emitted once, at the
    member with the lowest id, naming every member in id order.

    A member's group is read from its own `with` list, so this relies on LIT022's
    rule that every member names every other. In a chain 1-2, 2-3, which LIT022
    refuses, variant 3 reads its group as {2, 3} and leaves the question to 2, 2
    reads {1, 2, 3} and leaves it to 1, and 1 reads {1, 2}: the one question
    emitted names variants 1 and 2, and variant 3 appears nowhere (measured
    2026-11-18). The build refuses a corpus failing validation; the `literature
    questions` command parses the YAML through the models and runs no validate
    code, so on such a file its list can omit a record.
    """
    questions: list[Question] = []
    for finding in findings:
        fid: str = finding["id"]
        pmid: str = finding["publication"]
        own = _review_mark(
            finding,
            [
                {
                    "finding": fid,
                    "publication": pmid,
                    "variant": None,
                    "variant_name": None,
                    "anchor": question_anchor(fid),
                }
            ],
        )
        if own is not None:
            questions.append(own)
        variants = sorted(finding["variants"], key=lambda v: v["id"])
        by_id = {v["id"]: v for v in variants}

        for variant in variants:
            mark = _review_mark(variant, [_variant_record(finding, variant, None)])
            if mark is not None:
                questions.append(mark)
            unclear = variant["unclear"]
            if "zygosity" in unclear:
                entry = unclear["zygosity"]
                questions.append(
                    _mark(
                        "zygosity",
                        entry["reason"],
                        entry["settled_by"],
                        [_variant_record(finding, variant, "zygosity")],
                    )
                )
            if "phase" in unclear:
                entry = unclear["phase"]
                members = [
                    by_id[i]
                    for i in sorted({variant["id"], *variant["with"]})
                    if i in by_id and "phase" in by_id[i]["unclear"]
                ]
                if members[0]["id"] != variant["id"]:
                    continue
                questions.append(
                    _mark(
                        "phase",
                        entry["reason"],
                        entry["settled_by"],
                        [_variant_record(finding, m, "phase") for m in members],
                    )
                )
    return questions


def question_rows(questions: Sequence[Question]) -> list[Question]:
    """Merge equal `(kind, reason, settled_by)` questions into register rows.

    Records keep their first-appearance order. Rows run in `SETTLED_ORDER`, then
    kind (zygosity, phase, pending review, not accessible), then first appearance.
    """
    merged: dict[tuple[str, str, str], Question] = {}
    for question in questions:
        key = (question["kind"], question["reason"], question["settled_by"])
        row = merged.get(key)
        if row is None:
            merged[key] = _mark(
                question["kind"],
                question["reason"],
                question["settled_by"],
                list(question["records"]),
            )
        else:
            row["records"].extend(question["records"])
    rows = list(merged.values())  # insertion order is first appearance
    rows.sort(key=lambda r: (SETTLED_ORDER.index(r["settled_by"]), _KIND_ORDER.index(r["kind"])))
    return rows
