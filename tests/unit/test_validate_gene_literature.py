# tests/unit/test_validate_gene_literature.py
"""LIT001-LIT014: one gene's literature file against itself and the registries.

The base file is clean and every count in it differs from its neighbours
(§4.30): search 1 found 3 new papers, search 2 found 1 new of 4 hits, so a
mutant that attributes by `hits`, by the latest search or by the file total
produces a different answer from the correct one.
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.corpus import Corpus
from chd_atlas.genes import GeneRegistry
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import PhenotypeTerm
from chd_atlas.validate.gene_literature import validate_gene_literature
from chd_atlas.vocab import LesionGroup
from tests.unit._stamp import stamped

TBX5 = "HGNC:11604"
OTHER = "HGNC:1100"
F1 = "CHDA:LIT:0000001"
F2 = "CHDA:LIT:0000002"
PHENOTYPES = (
    PhenotypeTerm(id="HP:0001631", label="atrial septal defect", lesion_group=LesionGroup.SEPTAL),
    PhenotypeTerm(id="HP:0001680", label="coarctation of aorta", lesion_group=LesionGroup.LVOTO),
)
# Gene -> studies with burden rows; PMID:1005 is the study the burden test links.
BURDEN = {TBX5: {"PMID:1005"}}
PUBS = {"PMID:1001", "PMID:1002", "PMID:1003", "PMID:1004", "PMID:1005"}

SYMBOLS = GeneRegistry(
    [
        {"hgnc_id": TBX5, "symbol": "TBX5", "aliases": "TBX5A", "prev_symbols": "HOS"},
        {"hgnc_id": OTHER, "symbol": "BRCA1", "aliases": "", "prev_symbols": ""},
    ]
)


def _finding(fid: str, pmid: str, **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": fid,
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
        "summary": "Six of seven probands carry a variant.",
        "extracted_by": "agent:extract/test",
        "approved_by": "A. Curator",
        "approved_on": "2026-11-05",
    }
    base.update(changes)
    return stamped(base)


def _base() -> dict[str, Any]:
    """Clean. Search 1 (run 2026-11-02): 3 hits, 3 new. Search 2 (run 2026-11-09):
    4 hits, 1 new. Triage 1001/1002/1003 -> search 1, 1004 -> search 2."""
    return {
        "gene": TBX5,
        "searches": [
            {
                "id": 1,
                "source": "pubmed",
                "query": "q1",
                "run_on": "2026-11-02",
                "hits": 3,
                "new_hits": 3,
            },
            {
                "id": 2,
                "source": "pubmed",
                "query": "q2",
                "run_on": "2026-11-09",
                "hits": 4,
                "new_hits": 1,
            },
        ],
        "triage": [
            {"pmid": "PMID:1001", "search": 1, "decision": "included"},
            {"pmid": "PMID:1002", "search": 1, "decision": "excluded", "reason": "not_chd"},
            {"pmid": "PMID:1003", "search": 1, "decision": "excluded", "reason": "review_article"},
            {"pmid": "PMID:1004", "search": 2, "decision": "included"},
        ],
        "findings": [
            _finding(F1, "PMID:1001"),
            _finding(F2, "PMID:1004", shares_probands_with=[F1], approved_on="2026-11-10"),
        ],
    }


def _corpus(*files: tuple[str, dict[str, Any]], phenotypes: Any = PHENOTYPES) -> Corpus:
    return Corpus(
        root=Path("/r"),
        phenotypes=tuple(phenotypes),
        literature=tuple(
            (
                Path("/r/curation/literature") / f"{stem}.yaml",
                # Tests edit findings after `_finding` built them; restamping keeps
                # LIT015 out of every assertion here (it has its own suite).
                GeneLiteratureFile.model_validate(
                    {**raw, "findings": [stamped(f) for f in raw["findings"]]}
                ),
            )
            for stem, raw in files
        ),
    )


def _run(raw: dict[str, Any] | None = None, stem: str = "TBX5", **kw: Any) -> list[ValidationIssue]:
    args: dict[str, Any] = {
        "symbols": SYMBOLS,
        "known_publications": PUBS,
        "known_cohorts": set(),
        "burden_studies": BURDEN,
        "published_genes": {TBX5},
    }
    args.update(kw)
    phenotypes = args.pop("phenotypes", PHENOTYPES)
    return validate_gene_literature(_corpus((stem, raw or _base()), phenotypes=phenotypes), **args)


def _codes(issues: list[ValidationIssue]) -> list[str]:
    return [i.code for i in issues]


def _edit(**paths: Any) -> dict[str, Any]:
    raw = copy.deepcopy(_base())
    for key, value in paths.items():
        raw[key] = value
    return raw


def test_the_base_file_is_clean() -> None:
    assert _run() == []


def test_a_rerun_search_counts_only_new_hits() -> None:
    """Attributed 3 and 1 is clean; all four on search 1 is LIT001 twice."""
    raw = _base()
    for entry in raw["triage"]:
        entry["search"] = 1
    issues = _run(raw)
    assert _codes(issues) == ["LIT001", "LIT001"]
    first, second = issues
    assert first.severity is Severity.ERROR
    assert first.location == "/r/curation/literature/TBX5.yaml"
    assert "search 1" in first.message and "4" in first.message and "3" in first.message
    assert "search 2" in second.message and "0" in second.message and "1" in second.message


def test_an_approval_earlier_than_its_search_is_refused_but_a_later_rerun_is_not() -> None:
    # F1's paper was found by search 1 (run 11-02); approved 11-05, which is
    # earlier than search 2's run 11-09. A re-run must not make that an error.
    assert _run() == []
    raw = _base()
    raw["findings"][0]["approved_on"] = "2026-11-01"  # before search 1
    raw["findings"][1]["approved_on"] = "2026-11-08"  # before search 2, which found 1004
    issues = _run(raw)
    assert _codes(issues) == ["LIT008", "LIT008"]
    assert [i.location for i in issues] == [
        f"/r/curation/literature/TBX5.yaml:finding {F1}",
        f"/r/curation/literature/TBX5.yaml:finding {F2}",
    ]


def _pending(raw: dict[str, Any]) -> None:
    raw["triage"][3] = {"pmid": "PMID:1004", "search": 2, "decision": "pending"}
    raw["findings"] = raw["findings"][:1]


def _three_pending(raw: dict[str, Any]) -> None:
    for entry in raw["triage"][1:]:
        entry.pop("reason", None)
        entry["decision"] = "pending"
    raw["findings"] = raw["findings"][:1]


def _drop_second_finding(raw: dict[str, Any]) -> None:
    raw["findings"] = raw["findings"][:1]


def _dup_triage(raw: dict[str, Any]) -> None:
    raw["triage"].append(
        {"pmid": "PMID:1002", "search": 2, "decision": "excluded", "reason": "not_chd"}
    )
    raw["searches"][1]["new_hits"] = 2


def _unincluded_finding(raw: dict[str, Any]) -> None:
    raw["findings"][1]["publication"] = "PMID:1002"


def _same_as_unincluded(raw: dict[str, Any]) -> None:
    raw["triage"][2] = {
        "pmid": "PMID:1003",
        "search": 1,
        "decision": "excluded",
        "reason": "same_data_as",
        "same_as": "PMID:1002",
    }


def _unknown_search(raw: dict[str, Any]) -> None:
    raw["triage"][3]["search"] = 7


def _dangling_share(raw: dict[str, Any]) -> None:
    raw["findings"][1]["shares_probands_with"] = ["CHDA:LIT:0000099"]


@pytest.mark.parametrize(
    ("mutate", "expected", "severity", "fragment"),
    [
        (_pending, ["LIT002"], Severity.WARNING, "1 "),
        (_three_pending, ["LIT002"], Severity.WARNING, "3 "),
        (_drop_second_finding, ["LIT003"], Severity.ERROR, "PMID:1004"),
        (_dup_triage, ["LIT005"], Severity.ERROR, "PMID:1002"),
        (_unincluded_finding, ["LIT004", "LIT003"], Severity.ERROR, "PMID:1002"),
        (_same_as_unincluded, ["LIT006"], Severity.ERROR, "PMID:1002"),
        (_unknown_search, ["LIT010", "LIT001"], Severity.ERROR, "7"),
        (_dangling_share, ["LIT007"], Severity.ERROR, "CHDA:LIT:0000099"),
    ],
)
def test_each_file_rule_fires_on_its_own_fault(
    mutate: Any, expected: list[str], severity: Severity, fragment: str
) -> None:
    raw = _base()
    mutate(raw)
    issues = _run(raw)
    assert sorted(_codes(issues)) == sorted(expected)
    assert issues[0].severity is severity or any(i.severity is severity for i in issues)
    assert any(fragment in i.message for i in issues)


def test_same_as_naming_an_included_paper_is_not_lit006_even_if_uncited() -> None:
    """1004 is included but has no finding: that is LIT003's fault, not LIT006's."""
    raw = _base()
    raw["triage"][2] = {
        "pmid": "PMID:1003",
        "search": 1,
        "decision": "excluded",
        "reason": "same_data_as",
        "same_as": "PMID:1004",
    }
    raw["findings"] = raw["findings"][:1]
    assert _codes(_run(raw)) == ["LIT003"]


def test_unknown_probands_are_not_lit007_and_a_pending_hit_is_only_a_warning() -> None:
    """The literal "unknown" is never LIT007; LIT002 (pending hits) is a warning, not an error."""
    assert "LIT007" not in _codes(_run())
    raw = _base()
    _pending(raw)
    issues = _run(raw)
    assert _codes(issues) == ["LIT002"]
    assert issues[0].severity is Severity.WARNING


def test_a_finding_cannot_share_probands_with_itself() -> None:
    raw = _base()
    raw["findings"][0]["shares_probands_with"] = [F1]
    issues = _run(raw)
    assert _codes(issues) == ["LIT007"]
    assert issues[0].location.endswith(f":finding {F1}") and "itself" in issues[0].message
    # Naming the other finding is not self-reference.
    raw["findings"][0]["shares_probands_with"] = [F2]
    assert _codes(_run(raw)) == []


@pytest.mark.parametrize(
    ("phenotypes", "groups", "expected"),
    [
        (["HP:0001631"], ["septal"], []),
        (["HP:0001631", "HP:0001680"], ["septal", "lvoto"], []),
        (["HP:0001631"], ["septal", "conotruncal"], ["REF009"]),  # a group nothing justifies
        (["HP:0001631", "HP:0001680"], ["septal"], ["REF007"]),  # lvoto is missing
        (["HP:0001631"], ["conotruncal"], ["REF007", "REF009"]),  # both directions at once
        # An unmapped term exempts the finding from both (no REF012 analogue, R29).
        (["HP:0001631", "HP:0009999"], ["septal", "conotruncal"], []),
        (["HP:0009999"], ["conotruncal"], []),
    ],
)
def test_a_findings_lesion_groups_must_agree_with_its_mapped_phenotypes(
    phenotypes: list[str], groups: list[str], expected: list[str]
) -> None:
    """Final review I3: the assertion rules REF007/REF009 now run on findings too.

    The cases differ in which of the two vocabularies is wrong, so a check keyed on
    the wrong one (or dropped) changes the answer (§4.30).
    """
    raw = _base()
    raw["findings"][0]["phenotypes"] = phenotypes
    raw["findings"][0]["lesion_groups"] = groups
    issues = _run(raw)
    assert _codes(issues) == expected
    assert all(i.location.endswith(f":finding {F1}") for i in issues)


@pytest.mark.parametrize(
    ("burden_studies", "expected"),
    [
        ({TBX5: {"PMID:1005"}}, []),
        ({TBX5: {"PMID:1001"}}, ["LIT014"]),  # a study, but not the one linked
        ({OTHER: {"PMID:1005"}}, ["LIT014"]),  # the right study, another gene's rows
        ({}, ["LIT014"]),  # no burden mirror at all
        (None, []),  # unreadable mirror: TBL000/TBL001 report it, so skip
    ],
)
def test_a_linked_burden_study_needs_a_burden_row_for_this_gene(
    burden_studies: dict[str, set[str]] | None, expected: list[str]
) -> None:
    raw = _base()
    raw["findings"][0]["burden_study"] = "PMID:1005"
    issues = _run(raw, burden_studies=burden_studies)
    assert _codes(issues) == expected
    assert all(i.location.endswith(f":finding {F1}") and "PMID:1005" in i.message for i in issues)


def test_included_pmids_without_a_registered_publication_are_refused() -> None:
    issues = _run(known_publications={"PMID:1001"})
    assert _codes(issues) == ["LIT011"]
    assert issues[0].location.endswith(f":finding {F2}") and "PMID:1004" in issues[0].message
    raw = _base()
    raw["findings"][0]["burden_study"] = "PMID:1005"
    assert _codes(_run(raw)) == []
    assert _codes(_run(raw, known_publications=PUBS - {"PMID:1005"})) == ["LIT011"]


def test_a_cohort_outside_the_registry_is_refused() -> None:
    raw = _base()
    raw["findings"][0]["cohorts"] = ["ddd", "pcgc"]
    issues = _run(raw, known_cohorts={"ddd"})
    assert _codes(issues) == ["LIT012"]
    assert "pcgc" in issues[0].message and "ddd" not in issues[0].message.replace("pcgc", "")
    assert _codes(_run(raw, known_cohorts={"ddd", "pcgc"})) == []


@pytest.mark.parametrize("kw", [{"known_publications": None}, {"known_cohorts": None}])
def test_a_missing_registry_skips_only_its_own_check(kw: dict[str, Any]) -> None:
    raw = _base()
    raw["findings"][0]["cohorts"] = ["pcgc"]
    full = _run(raw, known_publications={"PMID:1001"}, known_cohorts={"ddd"})
    assert sorted(_codes(full)) == ["LIT011", "LIT012"]
    key = next(iter(kw))
    kept = {"known_publications": {"PMID:1001"}, "known_cohorts": {"ddd"}}
    kept.update(kw)
    other = {"LIT011": "known_cohorts", "LIT012": "known_publications"}
    skipped = _run(raw, **kept)
    assert _codes(skipped) == [c for c, k in sorted(other.items()) if k == key]


@pytest.mark.parametrize(
    ("stem", "gene", "clean"),
    [
        ("TBX5", TBX5, True),
        ("TBX5", OTHER, False),  # another gene's id under TBX5's name
        ("HOS", TBX5, False),  # previous symbol of the right gene
        ("TBX5A", TBX5, False),  # alias of the right gene
        ("NOPE", TBX5, False),  # unknown
    ],
)
def test_the_file_stem_must_be_the_genes_approved_symbol(stem: str, gene: str, clean: bool) -> None:
    raw = _base()
    raw["gene"] = gene
    issues = _run(raw, stem=stem, published_genes={TBX5, OTHER})
    assert _codes(issues) == ([] if clean else ["LIT009"])
    if stem == "HOS":
        assert "TBX5" in issues[0].message and "approved" in issues[0].message
    if not clean:
        assert issues[0].location == f"/r/curation/literature/{stem}.yaml"


def test_no_symbol_registry_skips_lit009_only() -> None:
    assert _run(stem="HOS", symbols=None) == []


def test_a_gene_outside_the_publication_gate_is_a_warning_and_none_skips_it() -> None:
    issues = _run(published_genes={OTHER})
    assert _codes(issues) == ["LIT013"]
    assert issues[0].severity is Severity.WARNING
    assert issues[0].message == (
        f"curated literature for {TBX5} reaches no page: the gene is outside the publication gate"
    )
    assert _run(published_genes=None) == []


def test_findings_may_share_probands_across_files_and_order_is_deterministic() -> None:
    other = _base()
    other["gene"] = OTHER
    other["findings"] = [_finding("CHDA:LIT:0000003", "PMID:1001", shares_probands_with=[F1])]
    other["triage"] = [t for t in other["triage"] if t["pmid"] != "PMID:1004"]
    other["searches"][1]["hits"] = 1
    other["searches"][1]["new_hits"] = 0
    corpus = _corpus(("TBX5", _base()), ("BRCA1", other))
    issues = validate_gene_literature(
        corpus,
        symbols=SYMBOLS,
        known_publications=PUBS,
        known_cohorts=set(),
        burden_studies=BURDEN,
        published_genes={TBX5, OTHER},
    )
    # BRCA1's finding names a finding that lives in TBX5's file: resolved. Its
    # file leaves PMID:1002/1003 included-less, so only LIT-coverage is clean.
    assert "LIT007" not in _codes(issues)
