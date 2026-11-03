# tests/unit/test_cli_literature.py
"""`chd-atlas literature search | abstracts | add-publication`, with a fake PubMed.

The repository under test is a tmp directory: a four-gene registry, two
published genes (TBX5, GATA4 -- different diseases, so a label cannot pass by
coming from the wrong gene) and the real `curation/chd_scope.yaml`. Nothing here
reaches the network; `chd_atlas.cli._fetch` is the one transport hook and every
test replaces it.
"""

from __future__ import annotations

import json
import re
import urllib.error
from collections.abc import Callable
from datetime import date
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

import pytest
from ruamel.yaml import YAML
from typer.testing import CliRunner

from chd_atlas import cli
from chd_atlas.cli import app
from chd_atlas.models.literature import PublicationFile
from chd_atlas.pubmed import PubMedError, build_query
from chd_atlas.validate.runner import validate_repository

REPO = Path(__file__).resolve().parents[2]
runner = CliRunner()

GENES = (
    "hgnc_id\tsymbol\tname\taliases\tprev_symbols\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
    "HGNC:11604\tTBX5\tT-box 5\tTBX5A\tHOS\t\t\t\t\t\n"
    "HGNC:4173\tGATA4\tGATA binding protein 4\t\t\t\t\t\t\t\n"
    "HGNC:1100\tBRCA1\tBRCA1 DNA repair\t\t\t\t\t\t\t\n"
    "HGNC:6772\tSMAD6\tSMAD family member 6\t\tMADH7\t\t\t\t\t\n"
    "HGNC:6773\tSMAD7\tSMAD family member 7\t\tMADH7\t\t\t\t\t\n"
)
CLINGEN_HEADER = (
    "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
    "classification_date\tgcep\treport_url\n"
)
CLINGEN = (
    "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP11\tDefinitive\t"
    "2025-03-25T16:00:00.000Z\tSyndromic Disorders Gene Curation Expert Panel\thttps://x/1\n"
    "HGNC:4173\tGATA4\tMONDO:0100009\tstructural congenital heart disease, multiple types - GATA4\t"
    "AD\tSOP9\tDefinitive\t2023-02-07T17:00:00.000Z\t"
    "Congenital Heart Disease Gene Curation Expert Panel\thttps://x/2\n"
)
GENCC_HEADER = (
    "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\tsubmitter\t"
    "submitted_on\treport_url\n"
)

TBX5_QUERY = build_query("TBX5", ("HOS",), ("Holt-Oram syndrome",))

PUBLICATIONS = """\
# Head comment: every record is transcribed from PubMed.
publications:
  - id: PMID:111
    title: First paper.
    journal: Journal One
    year: 1999
    authors:
      - Aa A
    # A comment sitting above a neighbouring field.
    study_type: case_report
  # A comment sitting between two records.
  - id: PMID:222
    title: Second paper.
    journal: Journal Two
    year: 2001
    authors:
      - Bb B
    study_type: gwas
"""


@pytest.fixture
def root(tmp_path: Path) -> Path:
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "curation").mkdir()
    (tmp_path / "mirrors" / "genes.tsv").write_text(GENES)
    (tmp_path / "mirrors" / "clingen_gene_validity.tsv").write_text(CLINGEN_HEADER + CLINGEN)
    (tmp_path / "mirrors" / "gencc_submissions.tsv").write_text(GENCC_HEADER)
    (tmp_path / "curation" / "chd_scope.yaml").symlink_to(REPO / "curation" / "chd_scope.yaml")
    (tmp_path / "curation" / "publications.yaml").write_text(PUBLICATIONS)
    return tmp_path


class FakePubMed:
    """esearch / esummary / efetch over a dict of papers, every URL recorded."""

    def __init__(self, papers: dict[str, dict[str, object]]) -> None:
        self.papers = papers
        self.urls: list[str] = []

    def __call__(self, url: str) -> bytes:
        self.urls.append(url)
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        endpoint = parsed.path.rsplit("/", 1)[1]
        if endpoint == "esearch.fcgi":
            ids = list(self.papers)
            start, size = int(query["retstart"][0]), int(query["retmax"][0])
            return json.dumps(
                {
                    "esearchresult": {
                        "count": str(len(ids)),
                        "idlist": ids[start : start + size],
                    }
                }
            ).encode()
        wanted = query["id"][0].split(",")
        if endpoint == "esummary.fcgi":
            result: dict[str, object] = {"uids": [p for p in wanted if p in self.papers]}
            for pmid in wanted:
                if pmid in self.papers:
                    paper = self.papers[pmid]
                    result[pmid] = {
                        "uid": pmid,
                        "pubdate": f"{paper['year']} Jan",
                        "title": paper["title"],
                        "fulljournalname": "Journal of Tests",
                        "authors": [{"name": "Zed Z", "authtype": "Author"}],
                        "articleids": [{"idtype": "doi", "value": f"10.1000/t{pmid}"}],
                    }
            return json.dumps({"result": result}).encode()
        articles = ""
        for pmid in wanted:
            if pmid in self.papers:
                paper = self.papers[pmid]
                abstract = (
                    f"<Abstract><AbstractText>{paper['abstract']}</AbstractText></Abstract>"
                    if paper.get("abstract")
                    else ""
                )
                articles += (
                    f"<PubmedArticle><MedlineCitation><PMID>{pmid}</PMID><Article>"
                    f"<ArticleTitle>{paper['title']}</ArticleTitle>{abstract}"
                    f"</Article></MedlineCitation></PubmedArticle>"
                )
        return f"<PubmedArticleSet>{articles}</PubmedArticleSet>".encode()


def _papers(*pmids: str) -> dict[str, dict[str, object]]:
    return {
        p: {"title": f"Title of {p}", "year": 2000 + int(p) % 20, "abstract": f"Abstract of {p}."}
        for p in pmids
    }


def _install(monkeypatch: pytest.MonkeyPatch, fetch: Callable[[str], bytes]) -> None:
    monkeypatch.setattr(cli, "_fetch", fetch)
    monkeypatch.delenv("NCBI_API_KEY", raising=False)


def _search(root: Path, *args: str) -> object:
    return runner.invoke(app, ["literature", "search", *args, "--root", str(root)])


def _load(path: Path) -> dict[str, object]:
    loaded: dict[str, object] = YAML(typ="safe").load(path.read_text())
    return loaded


def test_a_first_search_creates_the_file_with_pending_entries(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Numeric order 9 < 25 < 100; string order would put "100" first.
    fake = FakePubMed(_papers("100", "9", "25"))
    _install(monkeypatch, fake)

    # "HOS" is TBX5's retired name (ruling R19): accepted, but the file is TBX5.yaml.
    result = _search(root, "HOS", "--run-on", "2026-11-02")

    assert result.exit_code == 0, result.stdout  # type: ignore[attr-defined]
    assert not (root / "curation" / "literature" / "HOS.yaml").exists()
    data = _load(root / "curation" / "literature" / "TBX5.yaml")
    assert data["gene"] == "HGNC:11604"
    assert data["searches"] == [
        {
            "id": 1,
            "source": "pubmed",
            "query": TBX5_QUERY,
            "run_on": date(2026, 11, 2),
            "hits": 3,
            "new_hits": 3,
        }
    ]
    assert data["triage"] == [
        {"pmid": "PMID:9", "search": 1, "decision": "pending"},
        {"pmid": "PMID:25", "search": 1, "decision": "pending"},
        {"pmid": "PMID:100", "search": 1, "decision": "pending"},
    ]
    # The query was built from the registry, not from what the curator typed.
    sent = parse_qs(urlparse(fake.urls[0]).query)["term"][0]
    assert sent == TBX5_QUERY
    # And the file leaves no LIT001 (every search's new_hits matches its entries).
    codes = {i.code for i in validate_repository(root).issues}
    assert "LIT001" not in codes
    assert "LIT002" in codes  # pending is a warning, by design


def test_a_rerun_appends_a_search_and_triages_only_the_new_papers(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _install(monkeypatch, FakePubMed(_papers("9", "25")))
    assert _search(root, "TBX5", "--run-on", "2026-11-02").exit_code == 0  # type: ignore[attr-defined]
    path = root / "curation" / "literature" / "TBX5.yaml"
    # A curator decides one paper and leaves comments.
    text = path.read_text()
    text = text.replace("triage:\n", "# curator: triage in progress\ntriage:\n", 1)
    text = text.replace("decision: pending", "decision: included # the TBX5 paper", 1)
    path.write_text(text)

    _install(monkeypatch, FakePubMed(_papers("9", "25", "7", "300")))
    result = _search(root, "TBX5", "--run-on", "2026-12-01")

    assert result.exit_code == 0, result.stdout  # type: ignore[attr-defined]
    after = path.read_text()
    assert "# curator: triage in progress" in after
    assert "decision: included # the TBX5 paper" in after
    data = _load(path)
    assert [(s["id"], s["hits"], s["new_hits"]) for s in data["searches"]] == [  # type: ignore[attr-defined]
        (1, 2, 2),
        (2, 4, 2),
    ]
    assert [(t["pmid"], t["search"], t["decision"]) for t in data["triage"]] == [  # type: ignore[attr-defined]
        ("PMID:9", 1, "included"),
        ("PMID:25", 1, "pending"),
        ("PMID:7", 2, "pending"),
        ("PMID:300", 2, "pending"),
    ]
    assert "LIT001" not in {i.code for i in validate_repository(root).issues}


class _Boom:
    def __init__(self, exc: Exception | None, body: bytes | None = None) -> None:
        self.exc, self.body = exc, body

    def __call__(self, url: str) -> bytes:
        if self.exc is not None:
            raise self.exc
        assert self.body is not None
        return self.body


def _count_changes(url: str) -> bytes:
    start = int(parse_qs(urlparse(url).query)["retstart"][0])
    count = 1200 if start == 0 else 1201
    ids = [str(n) for n in range(start, min(start + 500, count))]
    return json.dumps({"esearchresult": {"count": str(count), "idlist": ids}}).encode()


FAILURES = {
    "pubmed-error": _Boom(PubMedError("PubMed answered HTTP 500")),
    "url-error": _Boom(urllib.error.URLError("no route")),
    "http-429": _Boom(urllib.error.HTTPError("https://x", 429, "Too Many Requests", {}, None)),  # type: ignore[arg-type]
    "html-body": _Boom(None, b"<html>Service unavailable</html>"),
    "count-changes": _count_changes,
}


@pytest.mark.parametrize("name", sorted(FAILURES))
@pytest.mark.parametrize("existing", [False, True])
def test_a_pubmed_failure_exits_1_and_leaves_the_file_exactly_as_it_was(
    root: Path, monkeypatch: pytest.MonkeyPatch, name: str, existing: bool
) -> None:
    path = root / "curation" / "literature" / "TBX5.yaml"
    if existing:
        _install(monkeypatch, FakePubMed(_papers("9")))
        assert _search(root, "TBX5").exit_code == 0  # type: ignore[attr-defined]
        before = path.read_bytes()
    _install(monkeypatch, FAILURES[name])

    result = _search(root, "TBX5")

    assert result.exit_code == 1, result.stdout  # type: ignore[attr-defined]
    if existing:
        assert path.read_bytes() == before
    else:
        assert not path.exists()
        assert not path.parent.exists()  # not even the directory is created early
    # No stray temporary from the atomic write, and no directory litter.
    assert [p.name for p in (root / "curation").rglob("*.tmp")] == []


@pytest.mark.parametrize(
    "args",
    [
        ["NOSUCHGENE"],
        ["MADH7"],
        ["BRCA1"],
        ["TBX5", "--run-on", "yesterday"],
    ],
    ids=["unknown", "ambiguous", "outside-the-gate", "bad-date"],
)
def test_bad_arguments_exit_2_without_touching_pubmed_or_the_disk(
    root: Path, monkeypatch: pytest.MonkeyPatch, args: list[str]
) -> None:
    def never(url: str) -> bytes:
        raise AssertionError("the network was reached")

    _install(monkeypatch, never)

    result = _search(root, *args)

    assert result.exit_code == 2, result.stdout  # type: ignore[attr-defined]
    assert not (root / "curation" / "literature").exists()


def test_the_api_key_reaches_the_request_and_never_the_output(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seen: list[str] = []

    def leaky(url: str) -> bytes:
        seen.append(url)
        raise PubMedError(f"could not fetch {url}")

    monkeypatch.setattr(cli, "_fetch", leaky)
    monkeypatch.setenv("NCBI_API_KEY", "SEKRIT-KEY-123")

    result = _search(root, "TBX5")

    assert result.exit_code == 1  # type: ignore[attr-defined]
    assert "api_key=SEKRIT-KEY-123" in seen[0]
    assert "SEKRIT-KEY-123" not in result.stdout  # type: ignore[attr-defined]

    # A blank or whitespace-only variable is no key at all.
    seen.clear()
    monkeypatch.setenv("NCBI_API_KEY", "   ")
    _search(root, "TBX5")
    assert "api_key" not in seen[0]


def test_a_file_that_belongs_to_another_gene_is_refused_untouched(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(
        "gene: HGNC:4173\nsearches:\n"
        "  - {id: 1, source: pubmed, query: q, run_on: 2026-11-02, hits: 0, new_hits: 0}\n"
    )
    before = path.read_bytes()
    _install(monkeypatch, FakePubMed(_papers("9")))

    result = _search(root, "TBX5")

    assert result.exit_code == 1  # type: ignore[attr-defined]
    assert "HGNC:4173" in result.stdout  # type: ignore[attr-defined]
    assert path.read_bytes() == before


def _seed_triage(root: Path) -> Path:
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir(parents=True)
    path.write_text(
        "gene: HGNC:11604\n"
        "searches:\n"
        "  - {id: 1, source: pubmed, query: q, run_on: 2026-11-02, hits: 4, new_hits: 4}\n"
        "triage:\n"
        "  - {pmid: 'PMID:5', search: 1, decision: pending}\n"
        "  - {pmid: 'PMID:6', search: 1, decision: included}\n"
        "  - {pmid: 'PMID:7', search: 1, decision: pending}\n"
        "  - {pmid: 'PMID:8', search: 1, decision: excluded, reason: not_chd}\n"
    )
    return path


def test_abstracts_prints_pending_papers_and_marks_a_missing_abstract(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _seed_triage(root)
    papers = _papers("5", "6", "7", "8")
    papers["7"]["abstract"] = ""
    _install(monkeypatch, FakePubMed(papers))

    result = runner.invoke(app, ["literature", "abstracts", "TBX5", "--root", str(root)])

    assert result.exit_code == 0, result.stdout
    assert result.stdout.splitlines() == [
        "PMID:5 · 2005 · Title of 5",
        "Abstract of 5.",
        "",
        "PMID:7 · 2007 · Title of 7",
        "(no abstract in PubMed)",
        "",
    ]


def test_abstracts_filters_by_decision(root: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    _seed_triage(root)
    _install(monkeypatch, FakePubMed(_papers("5", "6", "7", "8")))

    included = runner.invoke(
        app, ["literature", "abstracts", "TBX5", "--decision", "included", "--root", str(root)]
    )
    bad = runner.invoke(
        app, ["literature", "abstracts", "TBX5", "--decision", "maybe", "--root", str(root)]
    )
    missing = runner.invoke(app, ["literature", "abstracts", "GATA4", "--root", str(root)])

    assert included.stdout.splitlines()[0] == "PMID:6 · 2006 · Title of 6"
    assert "PMID:5" not in included.stdout
    assert bad.exit_code == 2
    assert missing.exit_code == 1  # a known gene with no literature file yet


def test_add_publication_appends_and_keeps_every_existing_byte(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    papers = {
        "333": {"title": 'Variants in <TBX5> & "GATA4": a #study', "year": 2024, "abstract": "x"}
    }
    _install(monkeypatch, FakePubMed(papers))
    path = root / "curation" / "publications.yaml"
    old = path.read_text()

    result = runner.invoke(
        app,
        [
            "literature",
            "add-publication",
            "333",
            "--study-type",
            "cohort_sequencing",
            "--reading",
            "Methods: exome sequencing of 500 probands.\nNo controls are described.",
            "--root",
            str(root),
        ],
    )

    assert result.exit_code == 0, result.stdout
    new = path.read_text()
    assert new.startswith(old)
    assert "# A comment sitting between two records." in new
    assert "# A comment sitting above a neighbouring field." in new
    loaded = PublicationFile.model_validate(YAML(typ="safe").load(new))
    assert [p.id for p in loaded.publications] == ["PMID:111", "PMID:222", "PMID:333"]
    added = loaded.publications[-1]
    assert added.title == 'Variants in <TBX5> & "GATA4": a #study'
    assert (added.journal, added.year, added.authors, added.doi) == (
        "Journal of Tests",
        2024,
        ["Zed Z"],
        "10.1000/t333",
    )
    assert added.study_type.value == "cohort_sequencing"
    # The reading sits beside the value it justifies (publications.yaml's head comment).
    assert (
        new[len(old) :].count(
            "    # Methods: exome sequencing of 500 probands.\n"
            "    # No controls are described.\n"
            "    study_type: cohort_sequencing\n"
        )
        == 1
    )


def test_add_publication_to_a_file_without_a_final_newline(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "publications.yaml"
    path.write_text(PUBLICATIONS.rstrip("\n"))
    _install(monkeypatch, FakePubMed(_papers("333")))

    result = runner.invoke(
        app,
        [
            "literature",
            "add-publication",
            "PMID:333",
            "--study-type",
            "gwas",
            "--reading",
            "GWAS design in the abstract.",
            "--root",
            str(root),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert path.read_text().startswith(PUBLICATIONS.rstrip("\n") + "\n")
    assert (
        len(PublicationFile.model_validate(YAML(typ="safe").load(path.read_text())).publications)
        == 3
    )


@pytest.mark.parametrize(
    ("args", "code", "fetch"),
    [
        (["111", "--study-type", "gwas"], 1, FakePubMed(_papers("111"))),
        (["333", "--study-type", "not_a_type"], 2, FakePubMed(_papers("333"))),
        (["abc", "--study-type", "gwas"], 2, FakePubMed(_papers("333"))),
        (["333", "--study-type", "gwas"], 1, _Boom(urllib.error.URLError("down"))),
        (["999", "--study-type", "gwas"], 1, FakePubMed(_papers("333"))),
        # F8/R33: the reading is required, and a blank one is no reading.
        (["333", "--study-type", "gwas", "--reading", "  "], 2, FakePubMed(_papers("333"))),
        (["333", "--study-type", "gwas", "--reading", ""], 2, FakePubMed(_papers("333"))),
    ],
    ids=[
        "already-present",
        "bad-study-type",
        "bad-pmid",
        "network-down",
        "pmid-unknown-to-pubmed",
        "blank-reading",
        "empty-reading",
    ],
)
def test_add_publication_refusals_leave_the_file_untouched(
    root: Path,
    monkeypatch: pytest.MonkeyPatch,
    args: list[str],
    code: int,
    fetch: Callable[[str], bytes],
) -> None:
    path = root / "curation" / "publications.yaml"
    before = path.read_bytes()
    _install(monkeypatch, fetch)

    reading = [] if "--reading" in args else ["--reading", "r"]
    result = runner.invoke(
        app, ["literature", "add-publication", *args, *reading, "--root", str(root)]
    )

    assert result.exit_code == code, result.stdout
    assert path.read_bytes() == before
    if args[0] == "111":
        assert "already in publications.yaml" in result.stdout


def test_add_publication_without_a_reading_is_refused_by_the_command_line(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "publications.yaml"
    before = path.read_bytes()
    _install(monkeypatch, FakePubMed(_papers("333")))
    result = runner.invoke(
        app,
        ["literature", "add-publication", "333", "--study-type", "gwas", "--root", str(root)],
    )
    assert result.exit_code == 2
    assert path.read_bytes() == before


def test_a_malformed_esummary_exits_1_and_leaves_publications_untouched(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "publications.yaml"
    before = path.read_bytes()
    _install(monkeypatch, _Boom(None, b'{"result": {"uids": ["333"]}}'))

    result = runner.invoke(
        app,
        [
            "literature",
            "add-publication",
            "333",
            "--study-type",
            "gwas",
            "--reading",
            "r",
            "--root",
            str(root),
        ],
    )

    assert result.exit_code == 1, result.stdout
    assert result.stdout.startswith("error: esummary")
    assert path.read_bytes() == before


_CANONICAL = (
    "# top comment\n"
    "gene: HGNC:11604\n"
    "searches:\n"
    "  - id: 1\n"
    "    source: pubmed\n"
    "    query: q\n"
    "    run_on: 2026-11-02\n"
    "    hits: 1\n"
    "    new_hits: 1\n"
    "triage:\n"
    "  # about the first paper\n"
    "  - pmid: PMID:9\n"
    "    search: 1\n"
    "    decision: included\n"
)


def test_a_canonical_file_keeps_every_line_and_another_indent_is_normalised_once(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir(parents=True)
    _install(monkeypatch, FakePubMed(_papers("9", "25")))

    path.write_text(_CANONICAL)
    assert _search(root, "TBX5", "--run-on", "2026-12-01").exit_code == 0  # type: ignore[attr-defined]
    # Only additions: the old lines survive unchanged and in order.
    old_lines = _CANONICAL.splitlines()
    remaining = iter(path.read_text().splitlines())
    assert all(line in remaining for line in old_lines)

    # Column-0 dashes: not the tool's style.
    path.write_text(
        _CANONICAL.replace("  - ", "- ").replace("    ", "  ").replace("  # about", "# about")
    )
    old = _load(path)
    assert _search(root, "TBX5", "--run-on", "2026-12-01").exit_code == 0  # type: ignore[attr-defined]
    text = path.read_text()
    assert "# about the first paper" in text
    assert "# top comment" in text
    new = _load(path)
    assert new["searches"][0] == old["searches"][0]  # type: ignore[index]
    assert new["triage"][0] == old["triage"][0]  # type: ignore[index]
    assert "\n  - id: 1\n    source: pubmed\n" in text  # now canonical


# --------------------------------------------------------------------------
# diff / review (Task 4)
# --------------------------------------------------------------------------

LITERATURE = """\
# Curator's head comment -- must survive every write.
gene: HGNC:11604
searches:
  - id: 1
    source: pubmed
    query: TBX5[tiab]
    run_on: 2026-11-02
    hits: 3
    new_hits: 3
triage:
  - pmid: PMID:8988165
    search: 1
    decision: included
findings:
  - id: CHDA:LIT:0000001
    publication: PMID:8988165
    design: cohort_screen
    methods:
      - exome
    count_unit: probands
    screened: 7
    carriers:
      de_novo: 2
      inherited: 1
      origin_unknown: 3
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: both
    cohorts:
      - alpha_cohort
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Table 2 & <3>
    summary: Six of seven probands carry a variant.
    extracted_by: agent:extract/test
  # A comment between two findings.
  - id: CHDA:LIT:0000002
    publication: PMID:111
    design: common_variant_association
    methods:
      - genotyping_array
    cases: 11
    controls: 13
    p_value: 0.02
    odds_ratio: 1.5
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: supplementary
        file: s1.xlsx
        sheet: Sheet 1
        row: 4
    summary: A GWAS hit.
    extracted_by: agent:extract/test
    verified_by: agent:verify/earlier
  - id: CHDA:LIT:0000003
    publication: PMID:222
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - conotruncal
    syndromic: syndromic
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A proband.
    extracted_by: agent:extract/test
    approved_by: Dr. Ada Curator
    approved_on: 2026-10-01
  - id: CHDA:LIT:0000004
    publication: PMID:333
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 0
      inherited: 1
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Another proband.
    extracted_by: agent:extract/test
"""

VERIFIER = """\
findings:
  - id: CHDA:LIT:0000001
    publication: PMID:8988165
    design: cohort_screen
    methods:
      - exome
    count_unit: probands
    screened: 7
    carriers:
      de_novo: 2
      inherited: 1
      origin_unknown: 3
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: both
    cohorts:
      - alpha_cohort
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Table 2, worded differently
    summary: Different words, same finding.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000002
    publication: PMID:111
    design: common_variant_association
    methods:
      - genotyping_array
    cases: 11
    controls: 13
    p_value: 0.02
    odds_ratio: 1.5
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A GWAS hit.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000003
    publication: PMID:222
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - conotruncal
    syndromic: syndromic
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A proband.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000004
    publication: PMID:333
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Another proband, counted differently.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000005
    publication: PMID:999
    design: case_report
    methods:
      - exome
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Only the blind pass found this.
    extracted_by: agent:verify/blind
"""

FINDING_ID = re.compile(r"(?m)^  - id: (CHDA:LIT:\d{7})\n")


@pytest.fixture
def lit_path(root: Path) -> Path:
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir()
    path.write_text(LITERATURE)
    return path


def _blocks(text: str) -> dict[str, str]:
    """Each finding's exact text, by id; "head" is everything before the first."""
    marks = list(FINDING_ID.finditer(text))
    blocks = {"head": text[: marks[0].start()]}
    for mark, nxt in zip(marks, [*marks[1:], None], strict=True):
        blocks[mark.group(1)] = text[mark.start() : nxt.start() if nxt else len(text)]
    return blocks


def _second(root: Path, text: str = VERIFIER) -> Path:
    path = root / "local" / "second.yaml"
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return path


def _diff(root: Path, against: Path, *args: str) -> Any:
    return runner.invoke(
        app,
        ["literature", "diff", "TBX5", "--against", str(against), *args, "--root", str(root)],
    )


def _review(root: Path, *args: str) -> Any:
    return runner.invoke(app, ["literature", "review", "TBX5", *args, "--root", str(root)])


def test_the_diff_prints_what_differs_and_writes_nothing(root: Path, lit_path: Path) -> None:
    before = lit_path.read_bytes()
    result = _diff(root, _second(root))

    assert result.exit_code == 0, result.stdout  # differences are output, not an error
    assert result.stdout == (
        "agreed (3): CHDA:LIT:0000001, CHDA:LIT:0000002, CHDA:LIT:0000003\n"
        "only in the first pass (0): none\n"
        "only in the second pass (1): CHDA:LIT:0000005\n"
        "differ (1):\n"
        "  CHDA:LIT:0000004\n"
        '    carriers: first {"de_novo": 0, "inherited": 1, "origin_unknown": 0}'
        ' | second {"de_novo": 1, "inherited": 0, "origin_unknown": 0}\n'
    )
    assert lit_path.read_bytes() == before


def test_mark_verified_stamps_only_agreeing_unverified_findings(root: Path, lit_path: Path) -> None:
    before = _blocks(lit_path.read_text())
    result = _diff(root, _second(root), "--mark-verified", "agent:verify/blind")

    assert result.exit_code == 0, result.stdout
    assert result.stdout.endswith("marked verified (2): CHDA:LIT:0000001, CHDA:LIT:0000003\n")
    after = _blocks(lit_path.read_text())
    assert after["head"] == before["head"]
    # 1 and 3 agreed and were unverified: one line each, right after extracted_by.
    assert after["CHDA:LIT:0000001"] == before["CHDA:LIT:0000001"].replace(
        "    extracted_by: agent:extract/test\n",
        "    extracted_by: agent:extract/test\n    verified_by: agent:verify/blind\n",
    )
    assert after["CHDA:LIT:0000003"] == before["CHDA:LIT:0000003"].replace(
        "    extracted_by: agent:extract/test\n",
        "    extracted_by: agent:extract/test\n    verified_by: agent:verify/blind\n",
    )
    # 2 agreed but already carries a verifier: untouched. 4 differs: untouched.
    assert after["CHDA:LIT:0000002"] == before["CHDA:LIT:0000002"]
    assert after["CHDA:LIT:0000004"] == before["CHDA:LIT:0000004"]
    assert (
        "approved_by: Dr. Ada Curator\n    approved_on: 2026-10-01\n" in after["CHDA:LIT:0000003"]
    )


def test_a_second_pass_key_shared_by_two_findings_verifies_nothing(
    root: Path, lit_path: Path
) -> None:
    twin = (
        VERIFIER.split("  - id: CHDA:LIT:0000002")[0]
        .removeprefix("findings:\n")
        .replace("CHDA:LIT:0000001", "CHDA:LIT:0000006")
    )
    before = lit_path.read_bytes()
    second = _second(
        root, VERIFIER.replace("  - id: CHDA:LIT:0000002", twin + "  - id: CHDA:LIT:0000002", 1)
    )

    result = _diff(root, second, "--mark-verified", "agent:verify/blind")

    assert result.exit_code == 0, result.stdout
    assert "only in the first pass (1): CHDA:LIT:0000001" in result.stdout
    assert (
        "only in the second pass (3): CHDA:LIT:0000001, CHDA:LIT:0000005, CHDA:LIT:0000006"
        in result.stdout
    )
    assert result.stdout.endswith("marked verified (1): CHDA:LIT:0000003\n")
    assert (
        "verified_by: agent:verify/blind" not in _blocks(lit_path.read_text())["CHDA:LIT:0000001"]
    )
    assert lit_path.read_bytes() != before


@pytest.mark.parametrize(
    ("args", "needle"),
    [
        (["--mark-verified", "Dr. Ada"], "agent:"),
        (["--mark-verified", ""], "agent:"),
    ],
)
def test_diff_refuses_a_verifier_that_is_not_an_agent(
    root: Path, lit_path: Path, args: list[str], needle: str
) -> None:
    before = lit_path.read_bytes()
    result = _diff(root, _second(root), *args)
    assert result.exit_code == 2 and needle in result.stdout
    assert lit_path.read_bytes() == before


def test_diff_names_a_missing_file_a_missing_verifier_and_a_bad_verifier(
    root: Path, lit_path: Path
) -> None:
    missing = runner.invoke(
        app,
        ["literature", "diff", "GATA4", "--against", str(_second(root)), "--root", str(root)],
    )
    assert missing.exit_code == 2
    assert str(root / "curation" / "literature" / "GATA4.yaml") in missing.stdout

    absent = _diff(root, root / "local" / "nope.yaml")
    assert absent.exit_code == 2 and "nope.yaml" in absent.stdout

    broken = _diff(root, _second(root, "findings:\n  - id: CHDA:LIT:0000001\n"))
    assert broken.exit_code == 1 and "error:" in broken.stdout

    duplicate = _diff(root, _second(root, VERIFIER + VERIFIER.split("findings:\n", 1)[1]))
    assert duplicate.exit_code == 1 and "duplicate" in duplicate.stdout

    # a verifier pass is held to the same per-record rules (approved_by an agent)
    approved = _diff(
        root,
        _second(
            root,
            VERIFIER.replace(
                "    summary: A GWAS hit.\n",
                "    summary: A GWAS hit.\n    approved_by: agent:x\n    approved_on: 2026-11-02\n",
            ),
        ),
    )
    assert approved.exit_code == 1 and "approved_by" in approved.stdout


def test_the_review_listing_pins_one_finding(root: Path, lit_path: Path) -> None:
    result = _review(root)

    assert result.exit_code == 0, result.stdout
    blocks = result.stdout.split("\n\n")
    one = next(b for b in blocks if b.startswith("CHDA:LIT:0000001"))
    assert one == (
        "CHDA:LIT:0000001  https://pubmed.ncbi.nlm.nih.gov/8988165/  cohort_screen  UNVERIFIED\n"
        "  paper: not in curation/publications.yaml\n"
        "  methods: exome\n"
        "  comparator: none\n"
        "  burden study: none\n"
        "  variant classes: snv_indel\n"
        "  screened 7 probands; carriers 6 (de novo 2, inherited 1, origin unknown 3)\n"
        "  phenotypes: HP:0001631\n"
        "  lesion groups: septal\n"
        "  syndromic: both\n"
        "  cohorts: alpha_cohort\n"
        "  locators: Table 2 & <3>\n"
        "  shares probands with: unknown\n"
        "  extracted by: agent:extract/test\n"
        "  Six of seven probands carry a variant."
    )
    two = next(b for b in blocks if b.startswith("CHDA:LIT:0000002"))
    assert two == (
        "CHDA:LIT:0000002  https://pubmed.ncbi.nlm.nih.gov/111/  common_variant_association"
        "  agent:verify/earlier\n"
        "  paper: study_type case_report; own_lab no\n"
        "  methods: genotyping_array\n"
        "  comparator: none\n"
        "  burden study: none\n"
        "  variant classes: none\n"
        "  cases 11, controls 13, p 0.02, OR 1.5\n"
        "  phenotypes: HP:0001631\n"
        "  lesion groups: septal\n"
        "  syndromic: isolated\n"
        "  cohorts: none\n"
        "  locators: s1.xlsx / Sheet 1 / row 4\n"
        "  shares probands with: unknown\n"
        "  extracted by: agent:extract/test\n"
        "  A GWAS hit."
    )
    # CHDA:LIT:0000003 is approved: not listed. Order is by id.
    assert [b[:16] for b in blocks if b.strip()] == [
        "CHDA:LIT:0000001",
        "CHDA:LIT:0000002",
        "CHDA:LIT:0000004",
    ]


def test_the_review_listing_shows_floats_unrounded_and_the_papers_own_lab(
    root: Path, lit_path: Path
) -> None:
    """Final review I2: `:g` showed 1.23456789e-9 as 1.23457e-09; own_lab is from the registry."""
    lit_path.write_text(
        lit_path.read_text()
        .replace("p_value: 0.02", "p_value: 1.23456789e-9")
        .replace("odds_ratio: 1.5", "odds_ratio: 2.718281828")
    )
    (root / "curation" / "publications.yaml").write_text(
        PUBLICATIONS.replace("study_type: case_report", "study_type: gwas\n    own_lab: true")
    )
    two = next(b for b in _review(root).stdout.split("\n\n") if b.startswith("CHDA:LIT:0000002"))
    assert "\n  cases 11, controls 13, p 1.23456789e-09, OR 2.718281828\n" in two
    assert "\n  paper: study_type gwas; own_lab yes\n" in two


def test_the_review_listing_prints_the_verifiers_differences_under_each_finding(
    root: Path, lit_path: Path
) -> None:
    """R35: `--against` shows what the blind pass disagreed on, finding by finding."""
    # Only the first block, with `screened` changed: finding 1 is paired and differs;
    # findings 2 and 4 have no partner in this pass.
    first_only = VERIFIER.split("  - id: CHDA:LIT:0000002")[0].replace(
        "    screened: 7\n", "    screened: 8\n", 1
    )
    second = _second(root, first_only)
    result = _review(root, "--against", str(second))
    assert result.exit_code == 0, result.stdout
    blocks = result.stdout.split("\n\n")
    one = next(b for b in blocks if b.startswith("CHDA:LIT:0000001"))
    assert one.endswith("\n  differs: screened: first 7 | second 8")
    assert "no partner" not in one
    for other in ("CHDA:LIT:0000002", "CHDA:LIT:0000004"):
        block = next(b for b in blocks if b.startswith(other))
        assert block.endswith("\n  second pass: no partner finding")
        assert "differs:" not in block
    # A paired finding that agrees prints neither line.
    agreeing = _second(root, VERIFIER.split("  - id: CHDA:LIT:0000003")[0])
    agreed = _review(root, "--against", str(agreeing)).stdout.split("\n\n")
    two = next(b for b in agreed if b.startswith("CHDA:LIT:0000002"))
    assert "differs:" not in two and "no partner" not in two
    # Without the option nothing about the second pass is printed.
    assert "differs:" not in _review(root).stdout
    # Refused alongside an action; a missing file is a bad argument.
    assert (
        _review(root, "--against", str(second), "--pending-review", "CHDA:LIT:0000001").exit_code
        == 2
    )
    assert _review(root, "--against", str(root / "nope.yaml")).exit_code == 2


def test_the_review_listing_names_the_findings_that_share_probands(
    root: Path, lit_path: Path
) -> None:
    text = lit_path.read_text()
    assert text.count("    shares_probands_with: unknown\n") >= 1
    lit_path.write_text(
        text.replace(
            "    shares_probands_with: unknown\n",
            "    shares_probands_with:\n      - CHDA:LIT:0000004\n      - CHDA:LIT:0000002\n",
            1,
        )
    )
    one = next(b for b in _review(root).stdout.split("\n\n") if b.startswith("CHDA:LIT:0000001"))
    assert "\n  shares probands with: CHDA:LIT:0000004, CHDA:LIT:0000002\n" in one


def test_approve_sets_both_fields_and_changes_no_other_byte(root: Path, lit_path: Path) -> None:
    before = lit_path.read_text()
    result = _review(
        root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. Grace Hopper", "--on", "2026-11-03"
    )

    assert result.exit_code == 0, result.stdout
    assert result.stdout == "approved CHDA:LIT:0000002 (Dr. Grace Hopper, 2026-11-03)\n"
    old, new = _blocks(before), _blocks(lit_path.read_text())
    assert new["CHDA:LIT:0000002"] == old["CHDA:LIT:0000002"].replace(
        "    verified_by: agent:verify/earlier\n",
        "    verified_by: agent:verify/earlier\n"
        "    approved_by: Dr. Grace Hopper\n    approved_on: 2026-11-03\n",
    )
    for other in ("head", "CHDA:LIT:0000001", "CHDA:LIT:0000003", "CHDA:LIT:0000004"):
        assert new[other] == old[other]


def test_approve_defaults_the_date_to_today_and_takes_several_ids(
    root: Path, lit_path: Path
) -> None:
    result = _review(
        root,
        "--approve",
        "CHDA:LIT:0000004",
        "--approve",
        "CHDA:LIT:0000002",
        "--by",
        "Dr. Grace Hopper",
        "--allow-unverified",
    )
    assert result.exit_code == 0, result.stdout
    today = date.today().isoformat()
    assert result.stdout == (
        f"approved CHDA:LIT:0000002 (Dr. Grace Hopper, {today})\n"
        f"approved CHDA:LIT:0000004 (Dr. Grace Hopper, {today})\n"
    )
    new = _blocks(lit_path.read_text())
    assert new["CHDA:LIT:0000004"].endswith(
        f"    extracted_by: agent:extract/test\n    approved_by: Dr. Grace Hopper\n"
        f"    approved_on: {today}\n"
    )


@pytest.mark.parametrize(
    ("args", "needle"),
    [
        (["--approve", "CHDA:LIT:0000002", "--by", "agent:review/x"], "agent"),
        (["--approve", "CHDA:LIT:0000002", "--by", "  Agent:review/x"], "agent"),
        (["--approve", "CHDA:LIT:0000002", "--by", "   "], "--by"),
        (["--approve", "CHDA:LIT:0000002"], "--by"),
        (["--approve", "CHDA:LIT:0000002", "--by", "A", "--on", "2026-13-40"], "--on"),
        (["--approve", "CHDA:LIT:0000009", "--by", "A"], "CHDA:LIT:0000009"),
        # M1: a date LIT008 would refuse (search 1 ran 2026-11-02) or one not yet reached.
        (
            [
                "--approve",
                "CHDA:LIT:0000001",
                "--by",
                "A",
                "--allow-unverified",
                "--on",
                "2026-11-01",
            ],
            "LIT008",
        ),
        (["--approve", "CHDA:LIT:0000002", "--by", "A", "--on", "2999-01-01"], "future"),
        (["--approve", "CHDA:LIT:0000003", "--by", "Someone Else"], "already approved"),
        (["--approve", "CHDA:LIT:0000001", "--by", "A"], "unverified"),
        # one bad id refuses the whole request: nothing is written for the good one
        (
            ["--approve", "CHDA:LIT:0000002", "--approve", "CHDA:LIT:0000001", "--by", "A"],
            "unverified",
        ),
        (
            ["--approve", "CHDA:LIT:0000002", "--pending-review", "CHDA:LIT:0000001", "--by", "A"],
            "not both",
        ),
        (["--pending-review", "CHDA:LIT:0000001"], "--note"),
        (["--pending-review", "CHDA:LIT:0000001", "--note", "  "], "--note"),
        (["--pending-review", "CHDA:LIT:0000009", "--note", "x"], "CHDA:LIT:0000009"),
        (["--note", "x"], "--pending-review"),
        (["--by", "A"], "--approve"),
    ],
)
def test_every_refusal_exits_two_and_writes_nothing(
    root: Path, lit_path: Path, args: list[str], needle: str
) -> None:
    before = lit_path.read_bytes()
    result = _review(root, *args)
    assert result.exit_code == 2, result.stdout
    assert needle in result.stdout
    assert lit_path.read_bytes() == before


def test_approval_on_the_day_of_the_search_is_accepted(root: Path, lit_path: Path) -> None:
    """LIT008 refuses only a date *before* the search that found the paper."""
    result = _review(
        root,
        "--approve",
        "CHDA:LIT:0000001",
        "--by",
        "A",
        "--allow-unverified",
        "--on",
        "2026-11-02",
    )
    assert result.exit_code == 0, result.stdout


def test_pending_review_sets_the_pair_and_nothing_else(root: Path, lit_path: Path) -> None:
    before = _blocks(lit_path.read_text())
    result = _review(
        root, "--pending-review", "CHDA:LIT:0000003", "--note", "Table 2 says 7; the text says 8."
    )

    assert result.exit_code == 0, result.stdout
    after = _blocks(lit_path.read_text())
    assert after["CHDA:LIT:0000003"] == before["CHDA:LIT:0000003"] + (
        "    review_status: pending_review\n    review_note: Table 2 says 7; the text says 8.\n"
    )
    for other in ("head", "CHDA:LIT:0000001", "CHDA:LIT:0000002", "CHDA:LIT:0000004"):
        assert after[other] == before[other]
    # a second mark would overwrite the first note
    again = _review(root, "--pending-review", "CHDA:LIT:0000003", "--note", "other")
    assert again.exit_code == 2 and "already" in again.stdout


def test_review_without_a_file_names_the_expected_path(root: Path) -> None:
    result = _review(root)
    assert result.exit_code == 2
    assert str(root / "curation" / "literature" / "TBX5.yaml") in result.stdout


def test_with_publication_refuses_a_blank_reading() -> None:
    from chd_atlas.gene_literature_io import LiteratureFileError, with_publication
    from chd_atlas.models.literature import Publication

    publication = Publication(
        id="PMID:333",
        title="A paper",
        journal="Journal of Tests",
        year=2024,
        authors=["Zed Z"],
        study_type="gwas",
    )
    raw = PUBLICATIONS.encode()
    with pytest.raises(LiteratureFileError, match="reading of the study design must not be blank"):
        with_publication(raw, publication, reading="   ")


def test_a_write_the_model_would_refuse_is_never_made(lit_path: Path) -> None:
    from chd_atlas.gene_literature_io import LiteratureFileError, with_findings_updated

    raw = lit_path.read_bytes()
    with pytest.raises(LiteratureFileError):
        with_findings_updated(raw, {"CHDA:LIT:0000004": {"approved_by": "Dr. A"}})
    with pytest.raises(LiteratureFileError):
        with_findings_updated(
            raw,
            {"CHDA:LIT:0000004": {"approved_by": "agent:x", "approved_on": date(2026, 11, 2)}},
        )
    with pytest.raises(LiteratureFileError):
        with_findings_updated(raw, {"CHDA:LIT:0000099": {"verified_by": "agent:x"}})


def test_the_review_listing_is_ordered_by_id_not_by_file_order(root: Path, lit_path: Path) -> None:
    blocks = _blocks(lit_path.read_text())
    ids = ["CHDA:LIT:0000004", "CHDA:LIT:0000002", "CHDA:LIT:0000001"]
    lit_path.write_text(blocks["head"] + "".join(blocks[i] for i in ids))

    result = _review(root)

    assert result.exit_code == 0, result.stdout
    assert [line[:16] for line in result.stdout.splitlines() if line.startswith("CHDA:")] == [
        "CHDA:LIT:0000001",
        "CHDA:LIT:0000002",
        "CHDA:LIT:0000004",
    ]


# --------------------------------------------------------------------------
# what separates two findings must not end up inside one (fix round 1)
# --------------------------------------------------------------------------

_CHUNKS: dict[str, str] = {
    "publication": "    publication: PMID:{pmid}\n",
    "design": "    design: cohort_screen\n",
    "methods": "    methods:\n      - exome\n",
    "count_unit": "    count_unit: probands\n",
    "screened": "    screened: 7\n",
    "carriers": "    carriers:\n      de_novo: 2\n      inherited: 1\n      origin_unknown: 3\n",
    "variant_classes": "    variant_classes:\n      - snv_indel\n",
    "phenotypes": "    phenotypes:\n      - HP:0001631\n",
    "lesion_groups": "    lesion_groups:\n      - septal\n",
    "syndromic": "    syndromic: both\n",
    "shares_probands_with": "    shares_probands_with: unknown\n",
    "locators": "    locators:\n      - kind: in_text\n        element: Table 2\n",
    "summary": "    summary: Six of seven.\n",
    "extracted_by": "    extracted_by: agent:extract/test\n",
}
_HEADER = (
    "gene: HGNC:11604\nsearches:\n  - id: 1\n    source: pubmed\n    query: TBX5[tiab]\n"
    "    run_on: 2026-11-02\n    hits: 3\n    new_hits: 3\nfindings:\n"
)


def _entry(
    number: int,
    *,
    last: str | None = None,
    verified: bool,
    trailer: str = "",
    eol: str = "",
) -> str:
    """One finding; `last` moves that key to the end (hand-ordered files), `trailer`
    is whatever separates it from the next finding."""
    keys = list(_CHUNKS)
    if last is not None:
        keys.remove(last)
        keys.append(last)
    body = "".join(_CHUNKS[k].format(pmid=1000 + number) for k in keys)
    if eol:  # a same-line comment on the last key, which holds a collection
        body = body.replace(f"    {last}:\n", f"    {last}:  {eol}\n")
    if verified:
        body = body.replace(
            _CHUNKS["extracted_by"], _CHUNKS["extracted_by"] + "    verified_by: agent:v\n"
        )
    return f"  - id: CHDA:LIT:000000{number}\n" + body + trailer


LAYOUTS = {
    "blank line": (None, "\n", ""),
    "comment then blank": (None, "  # note about the next one\n\n", ""),
    "comment": (None, "  # note about the next one\n", ""),
    "list last then comment": ("locators", "  # note about the next one\n", ""),
    "list last then blank": ("locators", "\n", ""),
    "mapping last then comment and blank": ("carriers", "  # note\n\n", ""),
    "same-line comment on a mapping, then blank": ("carriers", "\n", "# from Table 2"),
    "same-line comment on a list, then blank": ("locators", "\n", "# from Table 2"),
}


def _assert_only_lines_added(before: str, after: str, edited: str, added: list[str]) -> None:
    old, new = _blocks(before), _blocks(after)
    for key in old:
        if key != edited:
            assert new[key] == old[key], key
    lines = old[edited].splitlines(keepends=True)
    content = [
        i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith("#")
    ]
    cut = content[-1] + 1
    assert new[edited] == "".join(lines[:cut] + added + lines[cut:])


def _layout_file(root: Path, layout: str, *, verified: bool) -> tuple[Path, str]:
    last, trailer, eol = LAYOUTS[layout]
    text = (
        _HEADER
        + _entry(1, verified=verified, trailer="\n")
        + _entry(2, last=last, verified=verified, trailer=trailer, eol=eol)
        + _entry(3, verified=verified)
    )
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return path, text


@pytest.mark.parametrize("layout", LAYOUTS)
def test_approve_adds_only_its_lines_whatever_separates_the_findings(
    root: Path, layout: str
) -> None:
    path, before = _layout_file(root, layout, verified=True)
    result = _review(root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. A", "--on", "2026-11-03")
    assert result.exit_code == 0, result.stdout
    _assert_only_lines_added(
        before,
        path.read_text(),
        "CHDA:LIT:0000002",
        ["    approved_by: Dr. A\n", "    approved_on: 2026-11-03\n"],
    )


@pytest.mark.parametrize("layout", LAYOUTS)
def test_pending_review_adds_only_its_lines_whatever_separates_the_findings(
    root: Path, layout: str
) -> None:
    path, before = _layout_file(root, layout, verified=False)
    result = _review(root, "--pending-review", "CHDA:LIT:0000002", "--note", "Table 2 says 7.")
    assert result.exit_code == 0, result.stdout
    _assert_only_lines_added(
        before,
        path.read_text(),
        "CHDA:LIT:0000002",
        ["    review_status: pending_review\n", "    review_note: Table 2 says 7.\n"],
    )


@pytest.mark.parametrize("layout", LAYOUTS)
def test_mark_verified_adds_only_its_line_whatever_separates_the_findings(
    root: Path, layout: str
) -> None:
    path, before = _layout_file(root, layout, verified=False)
    second = "findings:\n" + "".join(
        _entry(n + 6, verified=False).replace(f"PMID:{1000 + n + 6}", f"PMID:{1000 + n}")
        for n in (1, 2, 3)
    )
    result = _diff(root, _second(root, second), "--mark-verified", "agent:verify/blind")
    assert result.exit_code == 0, result.stdout
    assert result.stdout.endswith(
        "marked verified (3): CHDA:LIT:0000001, CHDA:LIT:0000002, CHDA:LIT:0000003\n"
    )
    after = path.read_text()
    for number in (1, 2, 3):
        # extracted_by is where the line lands unless a hand-ordered key follows it
        old, new = (
            _blocks(before)[f"CHDA:LIT:000000{number}"],
            _blocks(after)[f"CHDA:LIT:000000{number}"],
        )
        assert new.count("    verified_by: agent:verify/blind\n") == 1
        assert new.replace("    verified_by: agent:verify/blind\n", "") == old


def test_approving_a_finding_under_review_keeps_the_mark_and_the_note(
    root: Path, lit_path: Path
) -> None:
    _review(root, "--pending-review", "CHDA:LIT:0000002", "--note", "Table 2 says 7; text says 8.")
    result = _review(root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. A", "--on", "2026-11-03")

    assert result.exit_code == 0, result.stdout
    block = _blocks(lit_path.read_text())["CHDA:LIT:0000002"]
    assert block.endswith(
        "    verified_by: agent:verify/earlier\n"
        "    approved_by: Dr. A\n    approved_on: 2026-11-03\n"
        "    review_status: pending_review\n    review_note: Table 2 says 7; text says 8.\n"
    )
