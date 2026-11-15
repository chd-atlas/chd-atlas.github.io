# tests/unit/test_build_gene_literature.py
"""What a gene's `literature` bundle key publishes, and the counts it feeds.

Every count in `_file` differs from the others (CLAUDE.md section 4.30):
6 included, 3 excluded as not_chd, 1 as review_article, 5 pending, 2 findings
awaiting approval, 4 approved. A mutant that publishes one count under another's
name, or counts the wrong decision, changes a number this file pins as a
literal.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.bundles import build_genes
from chd_atlas.build.emit import Emitter
from chd_atlas.build.gene_literature import (
    EMPTY_GENE_LITERATURE,
    LiteratureCensus,
    gene_literature,
    literature_census,
    literature_tallies,
)
from chd_atlas.corpus import Corpus
from chd_atlas.models.gene_literature import Finding, GeneLiteratureFile
from chd_atlas.models.literature import Publication
from tests.unit._stamp import stamped

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"
NKX2_5 = "HGNC:2488"


def _finding(number: int, pmid: str, **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": f"CHDA:LIT:{number:07d}",
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
        "summary": f"Summary {number}.",
        "extracted_by": "agent:extract/test",
        "approved_by": "A. Curator",
        "approved_on": "2026-11-05",
    }
    base.update(changes)
    return stamped(base)


def _file(gene: str = TBX5) -> GeneLiteratureFile:
    """6 included, 3 not_chd, 1 review_article, 5 pending; 4 approved, 2 not.

    Approved findings are listed out of order on purpose: PMID:9 and PMID:10 share
    a year (string order would put 10 first), PMID:2000 is earlier by year though
    its number is larger, and PMID:9 holds two findings listed 5 then 2.
    """
    triage: list[dict[str, Any]] = []
    number = 100

    def add(count: int, decision: str, reason: str | None = None) -> None:
        nonlocal number
        for _ in range(count):
            number += 1
            entry: dict[str, Any] = {"pmid": f"PMID:{number}", "search": 1, "decision": decision}
            if reason:
                entry["reason"] = reason
            triage.append(entry)

    add(6, "included")
    add(3, "excluded", "not_chd")
    add(1, "excluded", "review_article")
    add(5, "pending")
    return GeneLiteratureFile.model_validate(
        {
            "gene": gene,
            "searches": [
                {
                    "id": 2,
                    "source": "pubmed",
                    "query": "second",
                    "run_on": "2026-11-09",
                    "hits": 4,
                    "new_hits": 0,
                },
                {
                    "id": 1,
                    "source": "pubmed",
                    "query": "first",
                    "run_on": "2026-11-02",
                    "hits": 15,
                    "new_hits": 15,
                },
            ],
            "triage": triage,
            "findings": [
                _finding(5, "PMID:9"),
                _finding(7, "PMID:10"),
                _finding(2, "PMID:9"),
                _finding(4, "PMID:2000"),
                _finding(3, "PMID:10", approved_by=None, approved_on=None),
                _finding(6, "PMID:9", approved_by=None, approved_on=None),
            ],
        }
    )


def _publication(pmid: str, year: int) -> Publication:
    return Publication(
        id=pmid,  # type: ignore[arg-type]
        title=f"Title of {pmid}",
        journal="J",
        year=year,
        authors=["A B"],
        study_type="family_linkage",  # type: ignore[arg-type]
    )


def _corpus(*files: GeneLiteratureFile) -> Corpus:
    return Corpus(
        root=Path("."),
        publications=(
            _publication("PMID:9", 2005),
            _publication("PMID:10", 2005),
            _publication("PMID:2000", 1999),
        ),
        literature=tuple((Path(f"{f.gene}.yaml"), f) for f in files),
    )


def test_a_genes_literature_is_its_searches_its_tally_and_only_its_approved_findings() -> None:
    entry = gene_literature(_corpus(_file()))[TBX5]

    assert [s["query"] for s in entry["searches"]] == ["first", "second"]
    assert entry["searches"][0] == {
        "source": "pubmed",
        "query": "first",
        "run_on": "2026-11-02",
        "hits": 15,
        "new_hits": 15,
    }
    assert entry["triage"] == {
        "included": 6,
        "excluded": {"not_chd": 3, "review_article": 1},
        "pending": 5,
    }
    assert list(entry["triage"]["excluded"]) == ["not_chd", "review_article"]
    assert entry["pending_approval"] == 2
    # Approved only; sorted by (year, numeric PMID, id), not by PMID as a string.
    assert [(f["publication"], f["id"]) for f in entry["findings"]] == [
        ("PMID:2000", "CHDA:LIT:0000004"),
        ("PMID:9", "CHDA:LIT:0000002"),
        ("PMID:9", "CHDA:LIT:0000005"),
        ("PMID:10", "CHDA:LIT:0000007"),
    ]
    # JSON-ready: a date is its ISO string, not a `date`.
    assert entry["findings"][0]["approved_on"] == "2026-11-05"
    assert all("approved_digest" not in finding for finding in entry["findings"])
    assert json.dumps(entry)


def test_a_review_mark_is_published_with_an_approved_finding() -> None:
    raw = _file().model_dump(mode="json")
    raw["findings"][0].update(review_status="pending_review", review_note="Source disagrees.")
    entry = gene_literature(_corpus(GeneLiteratureFile.model_validate(raw)))[TBX5]
    marked = [f for f in entry["findings"] if f["review_status"] is not None]
    assert [(f["id"], f["review_note"]) for f in marked] == [
        ("CHDA:LIT:0000005", "Source disagrees.")
    ]


def test_a_file_whose_finding_cites_an_unregistered_paper_is_refused() -> None:
    corpus = Corpus(root=Path("."), publications=(), literature=((Path("TBX5.yaml"), _file()),))
    with pytest.raises(ValueError, match="not registered"):
        gene_literature(corpus)


def test_a_gene_with_no_file_publishes_the_empty_shape_and_one_with_a_file_its_own(
    tmp_path: Path,
) -> None:
    """Both outcomes in one build: GATA4 has no file, TBX5 has one, NKX2-5 is unpublished."""
    from tests.unit.test_build_bundles import _corpus as bundle_corpus  # noqa: PLC0415

    literature = gene_literature(_corpus(_file(), _file(NKX2_5)))
    emitter = Emitter(root=tmp_path)
    build_genes(
        bundle_corpus(),
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance={TBX5: {}, GATA4: {}},
        profiles={},
        cell_types={},
        literature=literature,
    )
    tbx5 = json.loads((tmp_path / "genes/HGNC_11604.json").read_text())
    gata4 = json.loads((tmp_path / "genes/HGNC_4173.json").read_text())
    assert tbx5["literature"] == json.loads(json.dumps(literature[TBX5]))
    assert gata4["literature"] == EMPTY_GENE_LITERATURE
    assert not (tmp_path / "genes/HGNC_2488.json").exists()


def test_the_empty_shape_is_literally_the_documented_one() -> None:
    assert EMPTY_GENE_LITERATURE == {
        "searches": [],
        "triage": {"included": 0, "excluded": {}, "pending": 0},
        "pending_approval": 0,
        "findings": [],
        "tallies": {
            "papers": 0,
            "first_year": None,
            "last_year": None,
            "designs": {},
            "lesion_groups": {},
            "de_novo_papers": 0,
            "inherited_papers": 0,
            "abstract_only_papers": 0,
            "pending_review_papers": 0,
        },
    }


def test_the_census_counts_published_genes_with_a_search_and_their_approved_findings() -> None:
    """GATA4 has a file with 4 approved findings but is unpublished: it counts for nothing.

    TBX5 -> 4 findings. A second published gene with one search and no findings
    makes `genes` (2) differ from `findings` (4) and from both gene ids' counts.
    """
    bare = GeneLiteratureFile.model_validate(
        {
            "gene": NKX2_5,
            "searches": [
                {
                    "id": 1,
                    "source": "pubmed",
                    "query": "q",
                    "run_on": "2026-11-02",
                    "hits": 0,
                    "new_hits": 0,
                }
            ],
        }
    )
    literature = gene_literature(_corpus(_file(), _file(GATA4), bare))
    census = literature_census(literature, {TBX5, NKX2_5})
    assert census == LiteratureCensus(genes=2, findings=4)
    assert literature_census(literature, {TBX5, NKX2_5, GATA4}) == LiteratureCensus(
        genes=3, findings=8
    )
    assert literature_census({}, {TBX5}) == LiteratureCensus(genes=0, findings=0)


def _f(number: int, pmid: str, **changes: Any) -> Finding:
    return Finding.model_validate(_finding(number, pmid, **changes))


_YEARS = {"PMID:9": 2003, "PMID:10": 2026, "PMID:11": 2010}


def test_tallies_count_papers_not_findings() -> None:
    """Four findings in two papers (two cohort screens in one): each category counts papers."""
    case = {
        "design": "case_report",
        "screened": 1,
        "carriers": {"de_novo": 1, "inherited": 0, "origin_unknown": 0},
    }
    tallies = literature_tallies(
        [
            _f(1, "PMID:9", lesion_groups=["septal", "avsd"]),
            _f(2, "PMID:9", **case, lesion_groups=["septal"]),
            _f(3, "PMID:10", lesion_groups=["septal"]),
            _f(4, "PMID:9", lesion_groups=["avsd"]),
        ],
        _YEARS,
    )
    assert tallies["papers"] == 2
    assert tallies["designs"] == {"case_report": 1, "cohort_screen": 2}
    assert tallies["lesion_groups"] == {"avsd": 1, "septal": 2}
    assert list(tallies["designs"]) == ["case_report", "cohort_screen"]


def test_tallies_read_null_carriers_as_none() -> None:
    common = {
        "design": "common_variant_association",
        "carriers": None,
        "screened": None,
        "count_unit": None,
        "cases": 40,
        "controls": 90,
        "variant_classes": [],
        "methods": ["exome"],
    }
    tallies = literature_tallies([_f(1, "PMID:9", **common)], _YEARS)
    assert tallies["papers"] == 1
    assert tallies["de_novo_papers"] == tallies["inherited_papers"] == 0


def test_tallies_count_carrier_origin_and_review_marks_by_paper() -> None:
    mark = {"review_note": "n"}
    tallies = literature_tallies(
        [
            _f(1, "PMID:9", review_status="not_accessible", **mark),
            _f(2, "PMID:9", review_status="not_accessible", **mark),
            _f(3, "PMID:10", review_status="pending_review", **mark),
            _f(
                4,
                "PMID:11",
                carriers={"de_novo": 0, "inherited": 4, "origin_unknown": 1},
                review_status="not_accessible",
                **mark,
            ),
        ],
        _YEARS,
    )
    assert tallies["abstract_only_papers"] == 2  # PMID:9 (twice) and PMID:11
    assert tallies["pending_review_papers"] == 1
    # PMID:9 and :10 carry de novo 2 each; PMID:11 has none but inherited 4.
    assert tallies["de_novo_papers"] == 2
    assert tallies["inherited_papers"] == 3


def test_tallies_span_the_publication_years() -> None:
    tallies = literature_tallies([_f(1, "PMID:10"), _f(2, "PMID:9"), _f(3, "PMID:11")], _YEARS)
    assert (tallies["first_year"], tallies["last_year"]) == (2003, 2026)


def test_a_gene_with_no_approved_finding_publishes_empty_tallies() -> None:
    raw = _file().model_dump(mode="json")
    for finding in raw["findings"]:
        finding.update(approved_by=None, approved_on=None, approved_digest=None)
    entry = gene_literature(_corpus(GeneLiteratureFile.model_validate(raw)))[TBX5]
    assert entry["findings"] == []
    assert entry["tallies"] == EMPTY_GENE_LITERATURE["tallies"]
    assert entry["tallies"]["first_year"] is None
    assert entry["tallies"]["last_year"] is None
