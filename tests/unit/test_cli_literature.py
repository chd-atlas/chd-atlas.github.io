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
import urllib.error
from collections.abc import Callable
from datetime import date
from pathlib import Path
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


def test_add_publication_to_a_file_without_a_final_newline(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "publications.yaml"
    path.write_text(PUBLICATIONS.rstrip("\n"))
    _install(monkeypatch, FakePubMed(_papers("333")))

    result = runner.invoke(
        app,
        ["literature", "add-publication", "PMID:333", "--study-type", "gwas", "--root", str(root)],
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
    ],
    ids=["already-present", "bad-study-type", "bad-pmid", "network-down", "pmid-unknown-to-pubmed"],
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

    result = runner.invoke(app, ["literature", "add-publication", *args, "--root", str(root)])

    assert result.exit_code == code, result.stdout
    assert path.read_bytes() == before
    if args[0] == "111":
        assert "already in publications.yaml" in result.stdout


def test_a_malformed_esummary_exits_1_and_leaves_publications_untouched(
    root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = root / "curation" / "publications.yaml"
    before = path.read_bytes()
    _install(monkeypatch, _Boom(None, b'{"result": {"uids": ["333"]}}'))

    result = runner.invoke(
        app,
        ["literature", "add-publication", "333", "--study-type", "gwas", "--root", str(root)],
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
