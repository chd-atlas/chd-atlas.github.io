# src/chd_atlas/gene_literature_diff.py
"""Comparing a blind second extraction against the first. Pure: no I/O.

Two findings are the same study when they share a key, `(publication, design,
sorted cohorts)`. A key carried by more than one finding *within either pass*
is ambiguous and is never paired (ruling R11): pairing it would let
`--mark-verified` verify a finding against an arbitrarily chosen twin, so every
finding carrying it is reported as present in one pass only and a curator
reads it.

Compared: every `Finding` field except `EXCLUDED_FIELDS`. The set is derived
from the model, so a field added later is compared by default instead of being
skipped without anyone noticing. Summaries and locators are free wording and
the curator reads them; the identifier, the extraction and approval
bookkeeping say who did what, not what the paper says. `shares_probands_with`
is excluded too (ruling R26): it names finding ids, the blind pass uses ids
local to its own file, so a list from the first pass could never agree with
anything; the curator reads it in the `review` listing.

Lists compared as sets, because their order carries no meaning: `methods`,
`variant_classes`, `phenotypes`, `lesion_groups` and `cohorts`. `carriers` is a mapping and is
compared by value. Reported values are the model's JSON-mode values, with the
set-like lists sorted.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any, Final

from chd_atlas.models.gene_literature import Finding

EXCLUDED_FIELDS: Final = (
    "id",
    "summary",
    "extracted_by",
    "verified_by",
    "approved_by",
    "approved_on",
    "approved_digest",
    "review_status",
    "review_note",
    "locators",
    "shares_probands_with",
)
COMPARED_FIELDS: Final[tuple[str, ...]] = tuple(
    name for name in Finding.model_fields if name not in EXCLUDED_FIELDS
)
_SET_LIKE: Final = frozenset(
    {
        "methods",
        "variant_classes",
        "phenotypes",
        "lesion_groups",
        "cohorts",
    }
)

FieldDifference = tuple[str, Any, Any]
_Key = tuple[str, str, tuple[str, ...]]


@dataclass(frozen=True)
class FindingDiff:
    """The result of `diff_findings`; every sequence is sorted.

    `only_first` / `only_second`: ids with no partner (including every finding
    of an ambiguous key). `agreed`: first-pass ids paired with a finding that
    has no field difference. `fields`: first-pass id -> its differences,
    sorted by field name, only for paired findings that differ.
    """

    only_first: tuple[str, ...]
    only_second: tuple[str, ...]
    agreed: tuple[str, ...]
    fields: dict[str, list[FieldDifference]]


def _key(finding: Finding) -> _Key:
    return (finding.publication, finding.design.value, tuple(sorted(finding.cohorts)))


def _values(finding: Finding) -> dict[str, Any]:
    dumped = finding.model_dump(mode="json")
    return {
        name: sorted(dumped[name])
        if name in _SET_LIKE and isinstance(dumped[name], list)
        else dumped[name]
        for name in COMPARED_FIELDS
    }


def diff_findings(first: Sequence[Finding], second: Sequence[Finding]) -> FindingDiff:
    first_counts = Counter(_key(f) for f in first)
    second_counts = Counter(_key(f) for f in second)
    pairable = {key for key, n in first_counts.items() if n == 1 and second_counts.get(key) == 1}
    partner = {_key(f): f for f in second if _key(f) in pairable}

    only_first: list[str] = []
    agreed: list[str] = []
    fields: dict[str, list[FieldDifference]] = {}
    for finding in sorted(first, key=lambda f: f.id):
        key = _key(finding)
        if key not in pairable:
            only_first.append(finding.id)
            continue
        mine, theirs = _values(finding), _values(partner[key])
        differences = [
            (n, mine[n], theirs[n]) for n in sorted(COMPARED_FIELDS) if mine[n] != theirs[n]
        ]
        if differences:
            fields[finding.id] = differences
        else:
            agreed.append(finding.id)
    only_second = sorted(f.id for f in second if _key(f) not in pairable)
    return FindingDiff(tuple(only_first), tuple(only_second), tuple(agreed), fields)
