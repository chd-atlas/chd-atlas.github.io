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
from dataclasses import replace
from datetime import date
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
    variant_table,
)
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.gene_literature import (
    Finding,
    GeneLiteratureFile,
    synthesis_digest,
)
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


def _run(number: int) -> dict[str, Any]:
    return {
        "id": number,
        "tool": "variantvalidator",
        "version": f"v{number}",
        "hgvs_version": "2.2",
        "vvta": "vvta_1",
        "vvdb": "vvdb_1",
        "seqrepo": "seqrepo_1",
        "run_on": "2026-11-17",
    }


def _variant_file(*, second_approved: bool) -> GeneLiteratureFile:
    """Runs 3, 1 and 2 listed out of order; the approved finding cites 2 and 1, an
    unapproved one cites 3. Variant 2 carries a non-empty `with`."""
    raw = _file().model_dump(mode="json", by_alias=True)
    raw["normalisations"] = [_run(3), _run(1), _run(2)]
    by_id = {f["id"]: f for f in raw["findings"]}

    def variant(number: int, run: int, **changes: Any) -> dict[str, Any]:
        base = {
            "id": number,
            "as_reported": f"as printed {number}",
            "locator": {"kind": "in_text", "element": "Table 1"},
            "kind": "snv_indel",
            "placement": "mane",
            "grch38": f"12:{number}:G:A",
            "normalisation": run,
        }
        return {**base, **changes}

    by_id["CHDA:LIT:0000002"]["variants"] = [
        variant(1, 2),
        variant(2, 1, **{"with": [1]}),
    ]
    by_id["CHDA:LIT:0000003"]["variants"] = [variant(1, 3)]
    for finding in raw["findings"]:
        if finding["id"] == "CHDA:LIT:0000002":
            finding["approved_digest"] = "sha256:probe"
    if not second_approved:
        by_id["CHDA:LIT:0000002"].update(approved_by=None, approved_on=None, approved_digest=None)
    raw["findings"] = [
        stamped({k: v for k, v in f.items() if k != "approved_digest"}) for f in raw["findings"]
    ]
    return GeneLiteratureFile.model_validate(raw)


def test_an_approved_findings_variants_publish_with_their_with_key_and_cited_normalisations() -> (
    None
):
    entry = gene_literature(_corpus(_variant_file(second_approved=True)))[TBX5]

    finding = next(f for f in entry["findings"] if f["id"] == "CHDA:LIT:0000002")
    assert [v["id"] for v in finding["variants"]] == [1, 2]
    linked = finding["variants"][1]
    assert linked["with"] == [1]
    assert "with_" not in linked
    assert linked["reported"] == {"transcript": None, "c": None, "p": None, "rsid": None}
    # A finding with no variant publishes `[]`, not an absent key.
    assert next(f for f in entry["findings"] if f["id"] == "CHDA:LIT:0000004")["variants"] == []
    # Cited runs only (run 3 is cited by an unapproved finding), sorted by id, JSON-ready.
    assert [n["id"] for n in entry["normalisations"]] == [1, 2]
    assert entry["normalisations"][0]["run_on"] == "2026-11-17"
    assert json.dumps(entry)


def test_an_unapproved_findings_normalisation_does_not_publish() -> None:
    entry = gene_literature(_corpus(_variant_file(second_approved=False)))[TBX5]

    assert all(f["variants"] == [] for f in entry["findings"])
    assert entry["normalisations"] == []


def test_a_review_mark_is_published_with_an_approved_finding() -> None:
    raw = _file().model_dump(mode="json")
    raw["findings"][0].update(
        review_status="pending_review",
        review_note="Source disagrees.",
        review_settled_by="authors",
    )
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
        "normalisations": [],
        "synthesis": None,
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
            _f(3, "PMID:10", review_status="pending_review", review_settled_by="authors", **mark),
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


OWN_AST = "CHDA:AST:0000001"
OTHER_AST = "CHDA:AST:0000002"


def _assertion(
    identifier: str, gene: str, mechanism: str = "haploinsufficiency"
) -> LesionAssertion:
    return LesionAssertion.model_validate(
        {
            "id": identifier,
            "gene": gene,
            "phenotypes": ["HP:0001631"],
            "lesion_groups": ["septal"],
            "inheritance": ["AD"],
            "mechanism": mechanism,
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


def _synthesis_file(
    approved: bool = True, text: str = "Variants cause septal defects."
) -> tuple[GeneLiteratureFile, Corpus]:
    """TBX5's file with a two-sentence synthesis citing findings 5 and 2 and its assertion."""
    base = _file()
    raw = base.model_dump(mode="json")
    raw["synthesis"] = {
        "sentences": [
            {"text": text, "cites": ["CHDA:LIT:0000005", "CHDA:LIT:0000002"]},
            {"text": "Loss of function is the mechanism.", "cites": [OWN_AST]},
        ],
        "drafted_by": "agent:gene-synthesis/test",
    }
    file = GeneLiteratureFile.model_validate(raw)
    assert file.synthesis is not None
    if approved:
        digest = synthesis_digest(file.synthesis, file.findings, [_assertion(OWN_AST, TBX5)])
        raw["synthesis"].update(
            approved_by="B. Person", approved_on="2026-11-12", approved_digest=digest
        )
        file = GeneLiteratureFile.model_validate(raw)
    return file, _corpus(file)


def _with_assertions(corpus: Corpus, *assertions: LesionAssertion) -> Corpus:
    return replace(corpus, assertions=assertions)


def test_a_current_approved_synthesis_is_published_and_withheld_when_stale_or_unapproved() -> None:
    """Four outcomes of one fixture: current, unapproved, text edited, no synthesis."""
    own = _assertion(OWN_AST, TBX5)
    file, corpus = _synthesis_file()
    published = gene_literature(_with_assertions(corpus, own))[TBX5]["synthesis"]
    assert published == {
        "sentences": [
            {
                "text": "Variants cause septal defects.",
                "cites": ["CHDA:LIT:0000005", "CHDA:LIT:0000002"],
            },
            {"text": "Loss of function is the mechanism.", "cites": [OWN_AST]},
        ],
        "drafted_by": "agent:gene-synthesis/test",
        "approved_by": "B. Person",
        "approved_on": "2026-11-12",
    }

    _, unapproved = _synthesis_file(approved=False)
    assert gene_literature(_with_assertions(unapproved, own))[TBX5]["synthesis"] is None

    # An edit to the text after approval: re-use the old digest on new text.
    assert file.synthesis is not None
    edited_raw = file.model_dump(mode="json")
    edited_raw["synthesis"]["sentences"][0]["text"] = "Edited after approval."
    edited = GeneLiteratureFile.model_validate(edited_raw)
    assert edited.synthesis is not None
    assert edited.synthesis.approved_digest == file.synthesis.approved_digest
    assert gene_literature(_with_assertions(_corpus(edited), own))[TBX5]["synthesis"] is None

    assert gene_literature(_corpus(_file()))[TBX5]["synthesis"] is None


def test_a_synthesis_goes_stale_through_its_genes_assertions_only() -> None:
    """Editing the gene's own assertion withholds it; editing another gene's does not."""
    _, corpus = _synthesis_file()
    own = _assertion(OWN_AST, TBX5)
    other = _assertion(OTHER_AST, GATA4)
    assert gene_literature(_with_assertions(corpus, own, other))[TBX5]["synthesis"] is not None

    edited_own = _assertion(OWN_AST, TBX5, mechanism="dominant_negative")
    assert gene_literature(_with_assertions(corpus, edited_own, other))[TBX5]["synthesis"] is None

    edited_other = _assertion(OTHER_AST, GATA4, mechanism="dominant_negative")
    assert (
        gene_literature(_with_assertions(corpus, own, edited_other))[TBX5]["synthesis"] is not None
    )


# --- variant_table -----------------------------------------------------------


def _var(vid: int, placement: str = "mane", **changes: Any) -> dict[str, Any]:
    """A published variant dict, as ``Finding.model_dump(by_alias=True)`` writes it."""
    base: dict[str, Any] = {
        "id": vid,
        "as_reported": f"var{vid}",
        "placement": placement,
        "reported": {"transcript": None, "c": None, "p": None, "rsid": None},
        "hgvs_c": None,
        "hgvs_p": None,
        "grch38": None,
        "zygosity": "heterozygous",
        "with": [],
        "phase": None,
        "carriers": None,
        "phenotypes": [],
        "review_status": None,
        "review_note": None,
    }
    base.update(changes)
    return base


def _pub(fid: int, pmid: str, variants: list[dict[str, Any]], **changes: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": f"CHDA:LIT:{fid:07d}",
        "publication": pmid,
        "design": "case_report",
        "count_unit": "probands",
        "carriers": None,
        "variants": variants,
    }
    base.update(changes)
    return base


def _mane(vid: int, c: str, g: str, **changes: Any) -> dict[str, Any]:
    return _var(vid, "mane", hgvs_c=c, grch38=g, **changes)


def test_one_variant_in_two_papers_is_one_row_with_two_reports_and_no_summed_carriers() -> None:
    first = {"de_novo": 2, "inherited": 0, "origin_unknown": 0}
    second = {"de_novo": 0, "inherited": 3, "origin_unknown": 1}
    other = _mane(2, "c.700A>G", "12:1:A:G")
    table = variant_table(
        [
            # arrival order is (year, PMID, id): finding 2 arrives before finding 1.
            _pub(
                2,
                "PMID:5",
                [
                    _mane(1, "c.100C>T", "12:9:C:T", carriers=first, phenotypes=["HP:0001631"]),
                    _mane(3, "c.5C>T", "12:8:C:T", **{"with": [4]}, phase="in_trans"),
                    _mane(4, "c.6C>T", "12:7:C:T", hgvs_p="p.Arg2Ter"),
                ],
            ),
            _pub(1, "PMID:7", [_mane(1, "c.100C>T", "12:9:C:T", carriers=second), other]),
            _pub(3, "PMID:7", [_mane(1, "c.100C>T", "12:9:C:T")]),
        ]
    )
    row = next(r for r in table["rows"] if r["key"] == "g:12:9:C:T")
    assert set(row) == {"key", "hgvs_c", "hgvs_p", "papers", "reports"}
    assert row["papers"] == 2  # three reports, two distinct publications
    assert [r["finding"] for r in row["reports"]] == [
        "CHDA:LIT:0000002",
        "CHDA:LIT:0000001",
        "CHDA:LIT:0000003",
    ]
    assert [r["carriers"] for r in row["reports"]] == [first, second, None]
    assert row["reports"][0]["phenotypes"] == ["HP:0001631"]
    assert row["reports"][0]["count_unit"] == "probands"
    assert row["reports"][0]["design"] == "case_report"
    assert row["reports"][0]["publication"] == "PMID:5"
    # a partner is named by its hgvs_c, scoped to its own finding
    partnered = next(r for r in table["rows"] if r["key"] == "g:12:8:C:T")
    assert partnered["reports"][0]["with_names"] == ["c.6C>T"]
    assert partnered["reports"][0]["phase"] == "in_trans"
    assert row["reports"][1]["with_names"] == []
    # no total of any kind: the only counts in the table are per-report
    assert "carriers" not in row and "total" not in row


def test_one_finding_naming_a_variant_alone_and_with_a_partner_is_one_row_two_reports() -> None:
    """Ruling D16: the same change in two genotype contexts of one finding is two reports
    on one row, each with its own carriers and partner names."""
    alone = {"de_novo": 1, "inherited": 0, "origin_unknown": 0}
    paired = {"de_novo": 0, "inherited": 2, "origin_unknown": 0}
    table = variant_table(
        [
            _pub(
                1,
                "PMID:1",
                [
                    _mane(1, "c.10A>G", "1:1:A:G", carriers=alone),
                    _mane(2, "c.10A>G", "1:1:A:G", carriers=paired, **{"with": [3]}),
                    _mane(3, "c.20A>G", "1:2:A:G", **{"with": [2]}),
                ],
            )
        ]
    )
    row = next(r for r in table["rows"] if r["key"] == "g:1:1:A:G")
    assert row["papers"] == 1
    assert [(r["finding"], r["carriers"], r["with_names"]) for r in row["reports"]] == [
        ("CHDA:LIT:0000001", alone, []),
        ("CHDA:LIT:0000001", paired, ["c.20A>G"]),
    ]


def test_a_partner_is_named_by_hgvs_c_then_hgvs_p_then_as_reported() -> None:
    table = variant_table(
        [
            _pub(
                1,
                "PMID:1",
                [
                    _mane(1, "c.10A>G", "1:1:A:G", **{"with": [2, 3, 4]}),
                    _mane(2, "c.20A>G", "1:2:A:G", hgvs_p="p.Lys7Arg"),
                    _var(3, "protein", hgvs_p="p.Lys9Arg"),
                    _var(4, "unplaced"),
                ],
            )
        ]
    )
    row = next(r for r in table["rows"] if r["key"] == "g:1:1:A:G")
    assert row["reports"][0]["with_names"] == ["c.20A>G", "p.Lys9Arg", "var4"]


def test_the_same_deletion_named_at_two_shifted_positions_is_one_row() -> None:
    table = variant_table(
        [
            _pub(1, "PMID:1", [_mane(1, "c.100_101del", "5:50:AGT:A")]),
            _pub(2, "PMID:2", [_mane(1, "c.101_102del", "5:50:AGT:A")]),
        ]
    )
    assert [r["key"] for r in table["rows"]] == ["g:5:50:AGT:A"]
    assert table["rows"][0]["hgvs_c"] == "c.100_101del"
    assert table["rows"][0]["papers"] == 2


def test_a_protein_row_is_never_merged_into_a_mane_row_with_the_same_p() -> None:
    table = variant_table(
        [
            _pub(1, "PMID:1", [_mane(1, "c.7C>T", "1:7:C:T", hgvs_p="p.Arg3Cys")]),
            _pub(2, "PMID:2", [_var(1, "protein", hgvs_p="p.Arg3Cys")]),
        ]
    )
    assert [(r["key"], r["papers"]) for r in table["rows"]] == [
        ("g:1:7:C:T", 1),
        ("p:R3C", 1),
    ]


def test_rows_sort_by_transcript_position_with_protein_rows_at_their_codon() -> None:
    # Arrival order is scrambled; the literal expected order is below. Protein
    # Arg100 sits at c.298 (3*100-2): after c.-20, before c.300-3; Gly101 at
    # c.301, after c.300+5 (offset is compared after base, not before).
    table = variant_table(
        [
            _pub(
                1,
                "PMID:1",
                [
                    _mane(1, "c.*5A>G", "1:9:A:G"),
                    _mane(2, "c.300+5A>G", "1:5:A:G"),
                    _var(3, "protein", hgvs_p="p.Gly101Arg"),
                    _mane(4, "c.300A>G", "1:4:A:G"),
                    _mane(5, "c.300-3A>G", "1:3:A:G"),
                    _var(6, "protein", hgvs_p="p.Arg100Cys"),
                    _mane(7, "c.-20A>G", "1:1:A:G"),
                    _mane(8, "c.-5A>G", "1:2:A:G"),
                ],
            )
        ]
    )
    assert [r["key"] for r in table["rows"]] == [
        "g:1:1:A:G",  # c.-20
        "g:1:2:A:G",  # c.-5
        "p:R100C",  # c.298
        "g:1:3:A:G",  # c.300-3
        "g:1:4:A:G",  # c.300
        "g:1:5:A:G",  # c.300+5
        "p:G101R",  # c.301
        "g:1:9:A:G",  # c.*5
    ]


def test_the_lists_hold_outside_region_then_structural_then_unplaced_variants() -> None:
    table = variant_table(
        [
            _pub(
                1,
                "PMID:1",
                [
                    _var(
                        1,
                        "unplaced",
                        as_reported="u1",
                        review_status="pending_review",
                        review_note="n",
                        review_settled_by="authors",
                    ),
                    _var(2, "structural", as_reported="s1"),
                    _var(3, "outside_region", as_reported="o1", reported={"rsid": "rs12"}),
                ],
            ),
            _pub(
                2,
                "PMID:2",
                [_var(1, "unplaced", as_reported="u2"), _var(2, "structural", as_reported="s2")],
            ),
        ]
    )
    assert table["rows"] == []
    assert [v["as_reported"] for v in table["outside_region"]] == ["o1"]
    assert table["outside_region"][0] == {
        "as_reported": "o1",
        "rsid": "rs12",
        "publication": "PMID:1",
        "finding": "CHDA:LIT:0000001",
        "variant": 3,
        "zygosity": "heterozygous",
        "phase": None,
        "review_status": None,
        "review_note": None,
    }
    assert [v["as_reported"] for v in table["no_position"]] == ["s1", "s2", "u1", "u2"]
    assert table["no_position"][2]["review_status"] == "pending_review"
    assert table["no_position"][2]["review_note"] == "n"


def test_findings_naming_no_variant_counts_only_carrying_non_association_findings() -> None:
    carriers = {"de_novo": 0, "inherited": 0, "origin_unknown": 1}
    zero = {"de_novo": 0, "inherited": 0, "origin_unknown": 0}
    table = variant_table(
        [
            _pub(1, "PMID:1", [], carriers=carriers),  # counted
            _pub(2, "PMID:2", [], carriers=carriers, design="cohort_screen"),  # counted
            _pub(3, "PMID:3", [], carriers=zero),
            _pub(4, "PMID:4", [], carriers=None),
            _pub(5, "PMID:5", [], carriers=carriers, design="common_variant_association"),
            _pub(6, "PMID:6", [_var(1, "unplaced")], carriers=carriers),
        ]
    )
    assert table["findings_naming_no_variant"] == 2


def test_a_report_keeps_its_own_names_when_a_key_is_named_differently() -> None:
    table = variant_table(
        [
            _pub(1, "PMID:1", [_mane(1, "c.12del", "1:1014:GA:G", hgvs_p="p.Arg4fs")]),
            _pub(2, "PMID:2", [_mane(1, "c.4del", "1:1014:GA:G", hgvs_p="p.Arg4Profs")]),
        ]
    )
    (row,) = table["rows"]
    assert (row["hgvs_c"], row["hgvs_p"]) == ("c.12del", "p.Arg4fs")
    assert [(r["hgvs_c"], r["hgvs_p"]) for r in row["reports"]] == [
        ("c.12del", "p.Arg4fs"),
        ("c.4del", "p.Arg4Profs"),
    ]


def test_a_protein_key_is_the_canonical_reading_so_spellings_share_a_row() -> None:
    table = variant_table(
        [
            _pub(1, "PMID:1", [_var(1, "protein", hgvs_p="NP_852259.1:p.(Arg279Ter)")]),
            _pub(2, "PMID:2", [_var(1, "protein", hgvs_p="p.Arg279Ter")]),
            _pub(3, "PMID:3", [_var(1, "protein", hgvs_p="p.R279X")]),
        ]
    )
    assert [(r["key"], r["papers"]) for r in table["rows"]] == [("p:R279*", 3)]
    assert table["rows"][0]["hgvs_p"] == "NP_852259.1:p.(Arg279Ter)"


def test_the_outside_region_list_keeps_finding_then_variant_order() -> None:
    table = variant_table(
        [
            _pub(
                2,
                "PMID:9",
                [
                    _var(2, "outside_region", as_reported="b"),
                    _var(1, "outside_region", as_reported="a"),
                ],
            ),
            _pub(1, "PMID:10", [_var(1, "outside_region", as_reported="c")]),
        ]
    )
    assert [(v["as_reported"], v["finding"]) for v in table["outside_region"]] == [
        ("b", "CHDA:LIT:0000002"),
        ("a", "CHDA:LIT:0000002"),
        ("c", "CHDA:LIT:0000001"),
    ]


@pytest.mark.parametrize(
    "variants",
    [
        [_mane(1, "c.5C>T", "1:5:C:T", **{"with": [9]})],
        [_var(1, "mane", hgvs_c="c.5C>T")],
        [_var(1, "protein")],
        [_var(1, "protein", hgvs_p="NP_852259.1:p.?")],
        [_mane(1, "not hgvs", "1:5:C:T")],
        [_var(1, "protein", hgvs_p="not hgvs")],
    ],
    ids=["dangling-partner", "no-grch38", "no-hgvs_p", "no-residue", "bad-c", "bad-p"],
)
def test_a_variant_a_bypassed_gate_lets_through_raises_rather_than_dropping(
    variants: list[dict[str, Any]],
) -> None:
    with pytest.raises(ValueError):
        variant_table([_pub(1, "PMID:1", variants)])
