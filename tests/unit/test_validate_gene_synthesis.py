"""SYN001-SYN004: a gene synthesis's citations, numbers, length and approval."""

from __future__ import annotations

import copy
from datetime import date
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.corpus import Corpus
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.gene_literature import GeneLiteratureFile, synthesis_digest
from chd_atlas.validate.gene_synthesis import validate_gene_synthesis
from tests.unit._stamp import stamped
from tests.unit.test_validate_gene_literature import F1, F2, OTHER, TBX5, _base

PATH = Path("/r/curation/literature/TBX5.yaml")
OWN_AST = "CHDA:AST:0000001"
OTHER_AST = "CHDA:AST:0000002"
F3 = "CHDA:LIT:0000003"


def _assertion(identifier: str, gene: str) -> LesionAssertion:
    return LesionAssertion.model_validate(
        {
            "id": identifier,
            "gene": gene,
            "phenotypes": ["HP:0001631"],
            "lesion_groups": ["septal"],
            "inheritance": ["AD"],
            "mechanism": "haploinsufficiency",
            "syndromic": "isolated",
            "evidence": [
                {
                    "publication": "PMID:8988165",
                    "evidence_class": "genetic_case",
                    "locator": {"kind": "in_text", "element": "Figure 1"},
                    "summary": "s",
                    "strength": "strong",
                }
            ],
            "curator": "c",
            "curated_on": date(2026, 7, 1),
            "last_reviewed": date(2026, 7, 1),
        }
    )


ASSERTIONS = (_assertion(OWN_AST, TBX5), _assertion(OTHER_AST, OTHER))


def _file(synthesis: dict[str, Any] | None, findings: list[dict[str, Any]] | None = None) -> Any:
    raw = copy.deepcopy(_base())
    if findings is not None:
        raw["findings"] = findings
    raw["findings"] = [stamped(f) for f in raw["findings"]]
    if synthesis is not None:
        raw["synthesis"] = synthesis
    return GeneLiteratureFile.model_validate(raw)


def _synthesis(*sentences: tuple[str, list[str]], **extra: Any) -> dict[str, Any]:
    return {
        "sentences": [{"text": t, "cites": c} for t, c in sentences],
        "drafted_by": "agent:gene-synthesis/test",
        **extra,
    }


def _run(
    synthesis: dict[str, Any],
    *others: tuple[Path, GeneLiteratureFile],
    findings: list[dict[str, Any]] | None = None,
) -> list[ValidationIssue]:
    corpus = Corpus(
        root=Path("/r"),
        assertions=ASSERTIONS,
        literature=((PATH, _file(synthesis, findings)), *others),
    )
    return validate_gene_synthesis(corpus)


def _codes(issues: list[ValidationIssue]) -> list[str]:
    return [i.code for i in issues]


def _approved(synthesis: dict[str, Any], file: Any = None) -> dict[str, Any]:
    from chd_atlas.models.gene_literature import Synthesis

    file = file or _file(None)
    own = [a for a in ASSERTIONS if a.gene == TBX5]
    digest = synthesis_digest(Synthesis.model_validate(synthesis), file.findings, own)
    return {
        **synthesis,
        "approved_by": "A. Curator",
        "approved_on": "2026-11-14",
        "approved_digest": digest,
    }


GOOD = _synthesis(("Septal defects are the usual presentation.", [F1, OWN_AST]))


def test_a_clean_unapproved_synthesis_and_a_file_without_one_raise_nothing() -> None:
    assert _run(GOOD) == []
    assert validate_gene_synthesis(Corpus(root=Path("/r"), literature=((PATH, _file(None)),))) == []


def test_an_approved_synthesis_with_a_current_digest_is_clean() -> None:
    assert _run(_approved(GOOD)) == []


@pytest.mark.parametrize(
    "cite",
    [
        OTHER_AST,  # an assertion of another gene
        F3,  # a finding id that exists in no file
    ],
)
def test_a_cite_that_is_not_an_approved_finding_or_own_assertion_is_syn001(cite: str) -> None:
    issues = _run(_synthesis(("Septal defects occur.", [F1, cite])))
    assert _codes(issues) == ["SYN001"]
    assert issues[0].severity is Severity.ERROR
    assert issues[0].location == f"{PATH}:synthesis.sentences[0]"
    assert cite in issues[0].message


def test_a_cite_of_an_unapproved_finding_of_the_same_file_is_syn001() -> None:
    raw = copy.deepcopy(_base())
    unapproved = {k: v for k, v in raw["findings"][1].items() if not k.startswith("approved")}
    issues = _run(
        _synthesis(("a", [F1]), ("Another sentence here.", [F2])),
        findings=[raw["findings"][0], unapproved],
    )
    assert _codes(issues) == ["SYN001"]
    assert issues[0].location == f"{PATH}:synthesis.sentences[1]"


def test_a_finding_that_exists_only_in_another_genes_file_is_syn001() -> None:
    other_raw = copy.deepcopy(_base())
    other_raw["gene"] = OTHER
    # Rename the other file's finding so it exists nowhere in TBX5's file.
    other_findings = [
        {**other_raw["findings"][0], "id": F3, "shares_probands_with": "unknown"},
    ]
    other_file = GeneLiteratureFile.model_validate(
        {**other_raw, "findings": [stamped(f) for f in other_findings]}
    )
    issues = _run(
        _synthesis(("Septal defects occur.", [F3])),
        (Path("/r/curation/literature/BRCA1.yaml"), other_file),
    )
    assert _codes(issues) == ["SYN001"]


@pytest.mark.parametrize(
    "text",
    ["a 3-generation family", "Two families", "in 2003", "most carriers", "one of twenty"],
)
def test_a_number_digit_or_word_is_syn002_naming_the_token(text: str) -> None:
    issues = _run(_synthesis((text, [F1])))
    assert _codes(issues) == ["SYN002"]
    assert issues[0].location == f"{PATH}:synthesis.sentences[0]"
    assert issues[0].severity is Severity.ERROR


@pytest.mark.parametrize(
    ("text", "token"),
    [("a 3-generation family", "3"), ("Two families", "Two"), ("most carriers", "most")],
)
def test_syn002_names_the_offending_token(text: str, token: str) -> None:
    assert f"'{token}'" in _run(_synthesis((text, [F1])))[0].message


@pytest.mark.parametrize(
    "text",
    [
        "NKX2-5 variants occur.",
        "A 22q11.2 deletion is syndromic.",
        "The variant p.Gly80Arg recurs.",
        "The variant c.2T>C recurs.",
        "Someone mostly agrees, and the honest answer is none.",
    ],
)
def test_mixed_letter_digit_tokens_and_embedded_words_are_not_numbers(text: str) -> None:
    assert _run(_synthesis((text, [F1]))) == []


def test_more_than_120_words_is_syn003_and_exactly_120_is_not() -> None:
    sentences = [("w" * 3 + " " + " ".join(["word"] * 39), [F1]) for _ in range(3)]  # 40 each
    assert _run(_synthesis(*sentences)) == []
    sentences[2] = (sentences[2][0] + " extra", [F1])
    issues = _run(_synthesis(*sentences))
    assert _codes(issues) == ["SYN003"]
    assert issues[0].location == f"{PATH}:synthesis"
    assert "121" in issues[0].message
    assert "120" in issues[0].message


def test_a_bare_punctuation_token_is_not_a_word() -> None:
    sentences = [("w " + " ".join(["word"] * 119) + " - - - -", [F1])]
    assert _run(_synthesis(*sentences)) == []


def test_a_stale_approval_is_a_syn004_warning_naming_the_review_command() -> None:
    stale = _approved(GOOD)
    stale["approved_digest"] = "sha256:" + "0" * 64
    issues = _run(stale)
    assert _codes(issues) == ["SYN004"]
    assert issues[0].severity is Severity.WARNING
    assert issues[0].location == f"{PATH}:synthesis"
    assert "chd-atlas literature review TBX5 --approve-synthesis" in issues[0].message


def test_an_edited_assertion_voids_the_approval() -> None:
    """The digest covers the gene's own assertions, so one that changed is stale."""
    approved = _approved(GOOD)
    changed = (_assertion(OWN_AST, TBX5), _assertion("CHDA:AST:0000003", TBX5), ASSERTIONS[1])
    corpus = Corpus(root=Path("/r"), assertions=changed, literature=((PATH, _file(approved)),))
    assert _codes(validate_gene_synthesis(corpus)) == ["SYN004"]


def test_another_genes_assertion_does_not_void_the_approval() -> None:
    approved = _approved(GOOD)
    more = (*ASSERTIONS, _assertion("CHDA:AST:0000004", OTHER))
    corpus = Corpus(root=Path("/r"), assertions=more, literature=((PATH, _file(approved)),))
    assert validate_gene_synthesis(corpus) == []
