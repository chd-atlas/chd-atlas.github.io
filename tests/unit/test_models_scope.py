# tests/unit/test_models_scope.py
"""The scope list is a selection among other people's claims, so it is shaped like one.

Every entry records which EXTERNAL authority treats the term as congenital heart
disease. That is not bookkeeping: "why is this gene here?" must have an answer,
and until 2026-08-06 the answer was a curator's name -- this file's own docstring
said inclusion "is the one thing about a gene the atlas asserts on its own
authority ... so it must have an answer with a name on it". It must have an
answer with an *authority* on it, and never an individual's: the atlas mirrors
every validity classification precisely because it cannot out-curate a chartered
expert panel, and deciding what counts as CHD is the same kind of judgement.

Adversarial review of #30 found this paragraph still asserting the old position
directly above a test named `..._and_names_no_individual`, in a file the change
had edited from line 12 down.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chd_atlas.models.scope import ChdScopeFile, ScopeAuthority

ENTRY = {
    "id": "MONDO:0007732",
    "label": "Holt-Oram syndrome",
    "reason": "cardiac septation defects are a defining feature",
    "admitted_by": "clingen_chd_panel",
    "attributed_to": "Congenital Heart Disease Gene Curation Expert Panel",
    "admitted_on": "2026-08-03",
}


def test_an_entry_requires_a_reason_an_admitter_and_a_date() -> None:
    assert ChdScopeFile.model_validate({"diseases": [ENTRY]}).diseases[0].id == "MONDO:0007732"
    for omitted in ("reason", "admitted_by", "attributed_to", "admitted_on"):
        incomplete = {k: v for k, v in ENTRY.items() if k != omitted}
        with pytest.raises(ValidationError):
            ChdScopeFile.model_validate({"diseases": [incomplete]})


def test_a_blank_reason_is_rejected() -> None:
    """An empty reason is worse than no field: it looks answered.

    `min_length` on a stripped value, because a reason of one space passes a
    bare `min_length` and renders as nothing -- the same defect `cell_text` in
    `tables.py` exists to catch on the mirror side.
    """
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [{**ENTRY, "reason": "   "}]})


def test_a_duplicate_term_is_rejected() -> None:
    """Two entries for one term means two reasons, and no rule picks between them."""
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [ENTRY, {**ENTRY, "reason": "again"}]})


def test_a_phenotype_term_is_not_a_disease_term() -> None:
    """`DiseaseId` is MONDO-only; an HP term names a feature, not an entity."""
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [{**ENTRY, "id": "HP:0001631"}]})


def test_the_committed_scope_list_loads_and_names_no_individual() -> None:
    """A smoke test over the real file, which no unit fixture can stand in for.

    The committed list is the atlas's scope *selection*. If it stops parsing,
    every gene falls out of scope at once.

    **The second assertion is the point of the 2026-08-06 change and is not a
    count.** Every entry named the project owner as the admitting authority,
    which made a non-clinician the authority on what counts as congenital heart
    disease. `ScopeAuthority` now makes a person's name unrepresentable in
    `admitted_by`, and this checks the other half: that no free-text field
    smuggles one back. Asserted against the raw bytes rather than the parsed
    model, because `reason` and `attributed_to` are both free text and either
    could carry a name the enum cannot.
    """
    from pathlib import Path

    from ruamel.yaml import YAML

    path = Path(__file__).parent.parent.parent / "curation" / "chd_scope.yaml"
    raw = path.read_text(encoding="utf-8")
    parsed = ChdScopeFile.model_validate(YAML(typ="safe").load(raw))

    # 68 until the Orphanet path landed (2026-10-25); 105 until MONDO:0016460 left (F29).
    assert len(parsed.diseases) == 104
    # Every authority is external, by construction of the enum.
    assert {entry.admitted_by for entry in parsed.diseases} <= set(ScopeAuthority)
    # And the file names no individual anywhere, comments included.
    assert "Audain" not in raw
    # A GenCC submitter must be named, never the aggregator: GenCC's largest
    # in-scope submitter is ClinGen itself, so "GenCC" hides who decided.
    for entry in parsed.diseases:
        if entry.admitted_by is ScopeAuthority.GENCC_SUBMITTER:
            assert entry.attributed_to != "GenCC", entry.id


# Noonan syndrome on Orphanet's frequency warrant (#8 part 1, 2026-10-25).
ORPHANET_ENTRY = {
    "id": "MONDO:0018997",
    "label": "Noonan syndrome",
    "reason": "Orphanet records abnormal pulmonary valve morphology at Frequent (79-30%).",
    "admitted_by": "orphanet_phenotype_frequency",
    "attributed_to": "Orphanet",
    "orphanet": "ORPHA:648",
    "admitted_on": "2026-10-25",
}


@pytest.mark.parametrize(
    "entry",
    [
        pytest.param({**ORPHANET_ENTRY, "orphanet": None}, id="frequency-warrant-without-a-code"),
        pytest.param({**ENTRY, "orphanet": "ORPHA:392"}, id="code-on-a-panel-warrant"),
        pytest.param({**ORPHANET_ENTRY, "attributed_to": "Orphanet experts"}, id="not-orphanet"),
        pytest.param({**ORPHANET_ENTRY, "orphanet": "Orphanet:648"}, id="mondo-spelling"),
    ],
)
def test_an_orphanet_code_belongs_to_the_orphanet_warrant_and_only_to_it(
    entry: dict[str, object],
) -> None:
    """The code names the disorder whose annotation is the warrant, so the two travel together.

    A frequency warrant with no code has nothing for SCP006 to check; a code on a
    panel's warrant would let a term look Orphanet-backed without being checked
    as one. `attributed_to` is pinned to "Orphanet" because the authority class
    already names the body, and a free-text variant would read as a different one.
    """
    assert ChdScopeFile.model_validate({"diseases": [ORPHANET_ENTRY]}).diseases[0].orphanet == (
        "ORPHA:648"
    )
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [entry]})
