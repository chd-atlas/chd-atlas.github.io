# tests/unit/test_pubmed.py
"""The PubMed client: query template, paging guards, parsers, request shape.

No test here reaches the network. `tests/fixtures/pubmed/` holds three
responses recorded from the live E-utilities on 2026-11-02 (one call each; the
commit body names the URLs); everything else is a fake `Fetch` that returns
bytes built in the test.
"""

from __future__ import annotations

import io
import json
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any
from urllib.parse import parse_qs, urlparse

import pytest

from chd_atlas import pubmed
from chd_atlas.models.literature import PublicationFile
from chd_atlas.pubmed import (
    CHD_TERMS,
    PaperMeta,
    PubMedError,
    build_query,
    efetch_abstracts,
    esearch_all,
    esummary,
    publication_from_meta,
    urllib_fetch,
)
from chd_atlas.vocab import StudyType

FIXTURES = Path(__file__).resolve().parents[1] / "fixtures" / "pubmed"
REPO = Path(__file__).resolve().parents[2]


def _esearch_page(count: int, ids: list[str]) -> bytes:
    return json.dumps(
        {
            "header": {"type": "esearch", "version": "0.3"},
            "esearchresult": {"count": str(count), "retmax": str(len(ids)), "idlist": ids},
        }
    ).encode()


class _Pages:
    """A fake esearch: `ids` served 500 at a time, each URL recorded."""

    def __init__(self, ids: list[str], counts: list[int] | None = None) -> None:
        self.ids = ids
        self.counts = counts
        self.urls: list[str] = []

    def __call__(self, url: str) -> bytes:
        self.urls.append(url)
        query = parse_qs(urlparse(url).query)
        start = int(query["retstart"][0])
        size = int(query["retmax"][0])
        count = self.counts[len(self.urls) - 1] if self.counts else len(self.ids)
        return _esearch_page(count, self.ids[start : start + size])


def test_chd_terms_are_pinned() -> None:
    """The spec's list with every "(s)" expanded singular-then-plural (ruling R6)."""
    assert CHD_TERMS == (
        "congenital heart",
        "heart defect",
        "heart defects",
        "cardiac defect",
        "cardiac defects",
        "heart malformation",
        "heart malformations",
        "cardiac malformation",
        "cardiac malformations",
        "septal defect",
        "septal defects",
        "tetralogy of Fallot",
        "transposition of the great arteries",
        "hypoplastic left heart",
        "heterotaxy",
        "conotruncal",
        "outflow tract",
        "atrioventricular septal defect",
        "coarctation",
    )


def test_build_query_is_the_template() -> None:
    query = build_query("TBX5", ("HOS",), ("Holt-Oram syndrome",))
    assert query == (
        '("TBX5"[tiab] OR "HOS"[tiab])'
        ' AND ("Heart Defects, Congenital"[mh] OR "congenital heart"[tiab]'
        ' OR "heart defect"[tiab] OR "heart defects"[tiab] OR "cardiac defect"[tiab]'
        ' OR "cardiac defects"[tiab] OR "heart malformation"[tiab]'
        ' OR "heart malformations"[tiab] OR "cardiac malformation"[tiab]'
        ' OR "cardiac malformations"[tiab] OR "septal defect"[tiab]'
        ' OR "septal defects"[tiab] OR "tetralogy of Fallot"[tiab]'
        ' OR "transposition of the great arteries"[tiab]'
        ' OR "hypoplastic left heart"[tiab] OR "heterotaxy"[tiab]'
        ' OR "conotruncal"[tiab] OR "outflow tract"[tiab]'
        ' OR "atrioventricular septal defect"[tiab] OR "coarctation"[tiab]'
        ' OR "Holt-Oram syndrome"[tiab])'
        " NOT preprint[pt] NOT review[pt]"
    )


def test_build_query_sorts_previous_symbols_and_diseases_and_drops_quotes() -> None:
    """Order of the arguments must not change the query; a quote cannot break it."""
    query = build_query("GATA4", ("ZZ1", "AA1"), ('beta "x" disease', "alpha disease"))
    assert query.startswith('("GATA4"[tiab] OR "AA1"[tiab] OR "ZZ1"[tiab]) AND (')
    assert query.endswith(
        ' OR "alpha disease"[tiab] OR "beta x disease"[tiab]) NOT preprint[pt] NOT review[pt]'
    )
    assert build_query("GATA4", (), ()).startswith('("GATA4"[tiab]) AND (')


def test_esearch_parses_the_recorded_page() -> None:
    body = (FIXTURES / "esearch_tbx5_holt_oram_1997.json").read_bytes()
    count, ids = esearch_all("q", lambda url: body, api_key=None)
    assert count == 4
    assert ids == ["9309177", "9207801", "8988165", "8988164"]


def test_esearch_pages_500_at_a_time_and_returns_every_id() -> None:
    ids = [str(n) for n in range(1, 1201)]
    pages = _Pages(ids)
    count, found = esearch_all("q", pages, api_key=None)
    assert (count, found) == (1200, ids)
    assert [parse_qs(urlparse(u).query)["retstart"][0] for u in pages.urls] == ["0", "500", "1000"]
    assert {parse_qs(urlparse(u).query)["retmax"][0] for u in pages.urls} == {"500"}


def test_a_count_that_changes_between_pages_refuses() -> None:
    pages = _Pages([str(n) for n in range(1, 1201)], counts=[1200, 1203, 1200])
    with pytest.raises(PubMedError, match="changed"):
        esearch_all("q", pages, api_key=None)


def test_a_count_past_the_esearch_ceiling_refuses_rather_than_truncating() -> None:
    """Measured 2026-11-02: retstart 9999 is refused by NCBI ("can only retrieve the
    first 9,999 records"). A silently truncated search would publish a coverage line
    over a subset of the hits, so the first page is enough to refuse."""
    pages = _Pages([str(n) for n in range(1, 501)], counts=[10000])
    with pytest.raises(PubMedError, match="9,999"):
        esearch_all("q", pages, api_key=None)
    assert len(pages.urls) == 1


@pytest.mark.parametrize(
    "body",
    [
        b"<html><body>Service Unavailable</body></html>",
        b'{"esearchresult": {"ERROR": "Search Backend failed"}}',
        b'{"esearchresult": {"count": "3", "idlist": ["1", "2"]}}',
        b'{"esearchresult": {"count": "2", "idlist": ["1", "1"]}}',
        b"[]",
    ],
    ids=["html", "error-key", "fewer-ids-than-count", "duplicate-ids", "not-an-object"],
)
def test_a_body_that_is_not_a_complete_result_refuses(body: bytes) -> None:
    with pytest.raises(PubMedError):
        esearch_all("q", lambda url: body, api_key=None)


@pytest.mark.parametrize(
    ("key", "expect_key"),
    [(None, False), ("", False), ("sekrit", True)],
)
def test_requests_carry_the_tool_and_the_key_only_when_there_is_one(
    key: str | None, expect_key: bool
) -> None:
    pages = _Pages(["1"])
    esearch_all("a b", pages, api_key=key)
    query = parse_qs(urlparse(pages.urls[0]).query)
    assert query["tool"] == ["chd-atlas"]
    assert query["term"] == ["a b"]
    assert query["db"] == ["pubmed"]
    assert ("api_key" in query) is expect_key
    assert "email" not in query
    assert "@" not in pages.urls[0]


def test_esummary_parses_the_recorded_record() -> None:
    body = (FIXTURES / "esummary_8988165.json").read_bytes()
    meta = esummary(["8988165"], lambda url: body)["8988165"]
    assert meta == PaperMeta(
        pmid="8988165",
        title="Mutations in human TBX5 [corrected] cause limb and cardiac malformation "
        "in Holt-Oram syndrome.",
        journal="Nature genetics",
        year=1997,
        authors=(
            "Basson CT",
            "Bachinsky DR",
            "Lin RC",
            "Levi T",
            "Elkins JA",
            "Soults J",
            "Grayzel D",
            "Kroumpouzou E",
            "Traill TA",
            "Leblanc-Straceski J",
            "Renault B",
            "Kucherlapati R",
            "Seidman JG",
            "Seidman CE",
        ),
        doi="10.1038/ng0197-30",
        pmcid=None,
    )


def test_the_built_publication_equals_the_committed_record() -> None:
    """Field by field against `curation/publications.yaml`'s own PMID:8988165."""
    body = (FIXTURES / "esummary_8988165.json").read_bytes()
    built = publication_from_meta(
        esummary(["8988165"], lambda url: body)["8988165"], StudyType.FAMILY_LINKAGE
    )
    from ruamel.yaml import YAML

    committed = {
        p.id: p
        for p in PublicationFile.model_validate(
            YAML(typ="safe").load((REPO / "curation" / "publications.yaml").read_text())
        ).publications
    }["PMID:8988165"]
    for field in ("id", "title", "journal", "year", "authors", "study_type", "doi", "pmcid"):
        assert getattr(built, field) == getattr(committed, field), field


def test_esummary_reports_a_pmid_pubmed_did_not_return_by_omitting_it() -> None:
    body = (FIXTURES / "esummary_8988165.json").read_bytes()
    assert set(esummary(["8988165", "1"], lambda url: body)) == {"8988165"}


def test_esummary_year_and_pmcid_and_collective_authors() -> None:
    record = {
        "uid": "5",
        "pubdate": "2023 Winter",
        "epubdate": "2022 Dec 30",
        "title": "A & B <i>x</i>.",
        "fulljournalname": "Some Journal",
        "authors": [
            {"name": "A B", "authtype": "Author"},
            {"name": "Group X", "authtype": "CollectiveName"},
        ],
        "articleids": [
            {"idtype": "pubmed", "value": "5"},
            {"idtype": "pmc", "value": "PMC123"},
        ],
    }
    body = json.dumps({"result": {"uids": ["5"], "5": record}}).encode()
    meta = esummary(["5"], lambda url: body)["5"]
    assert (meta.year, meta.pmcid, meta.doi, meta.authors) == (
        2023,
        "PMC123",
        None,
        ("A B", "Group X"),
    )
    assert meta.title == "A & B <i>x</i>."


def test_esummary_with_no_year_refuses() -> None:
    body = json.dumps(
        {"result": {"uids": ["5"], "5": {"title": "t", "fulljournalname": "J", "authors": []}}}
    ).encode()
    with pytest.raises(PubMedError):
        esummary(["5"], lambda url: body)


def test_efetch_parses_the_recorded_abstract() -> None:
    body = (FIXTURES / "efetch_8988165.xml").read_bytes()
    title, abstract = efetch_abstracts(["8988165"], lambda url: body)["8988165"]
    assert title.startswith("Mutations in human TBX5 [corrected] cause limb")
    assert abstract.startswith("Holt-Oram syndrome is characterized by upper limb")
    assert abstract.endswith("haploinsufficiency of TBX5 causes Holt-Oram syndrome.")


def test_efetch_labels_nested_markup_and_a_missing_abstract() -> None:
    """Synthetic XML (not recorded): structured labels, inline markup, no abstract."""
    xml = (
        b"<PubmedArticleSet>"
        b"<PubmedArticle><MedlineCitation><PMID>1</PMID><Article>"
        b"<ArticleTitle>T <i>one</i></ArticleTitle><Abstract>"
        b'<AbstractText Label="BACKGROUND">We saw <i>TBX5</i> p &lt; 0.05.</AbstractText>'
        b'<AbstractText Label="RESULTS">Four carriers.</AbstractText>'
        b"</Abstract></Article></MedlineCitation></PubmedArticle>"
        b"<PubmedArticle><MedlineCitation><PMID>2</PMID><Article>"
        b"<ArticleTitle>Two</ArticleTitle></Article></MedlineCitation></PubmedArticle>"
        b"</PubmedArticleSet>"
    )
    found = efetch_abstracts(["1", "2"], lambda url: xml)
    assert found["1"] == (
        "T one",
        "BACKGROUND: We saw TBX5 p < 0.05.\nRESULTS: Four carriers.",
    )
    assert found["2"] == ("Two", "")


def test_efetch_and_esummary_batch_long_id_lists() -> None:
    urls: list[str] = []

    def fetch(url: str) -> bytes:
        urls.append(url)
        return b"<PubmedArticleSet></PubmedArticleSet>"

    efetch_abstracts([str(n) for n in range(450)], fetch)
    assert len(urls) == 3
    assert all(len(parse_qs(urlparse(u).query)["id"][0].split(",")) <= 200 for u in urls)


def test_urllib_fetch_never_echoes_the_key_and_sets_a_timeout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, Any] = {}

    def refuse(request: Any, timeout: float | None = None) -> Any:
        seen["timeout"] = timeout
        raise urllib.error.HTTPError(
            "https://x/esearch.fcgi?api_key=sekrit&term=a",
            429,
            "Too Many Requests",
            {},
            None,  # type: ignore[arg-type]
        )

    monkeypatch.setattr(urllib.request, "urlopen", refuse)
    monkeypatch.setattr(pubmed, "_throttle", pubmed._Throttle(lambda: 0.0, lambda s: None))
    with pytest.raises(PubMedError) as caught:
        urllib_fetch("https://x/esearch.fcgi?api_key=sekrit&term=a")
    assert "sekrit" not in str(caught.value)
    assert "429" in str(caught.value)
    assert seen["timeout"] is not None and seen["timeout"] > 0


def test_urllib_fetch_returns_the_body(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(urllib.request, "urlopen", lambda request, timeout=None: io.BytesIO(b"ok"))
    monkeypatch.setattr(pubmed, "_throttle", pubmed._Throttle(lambda: 0.0, lambda s: None))
    assert urllib_fetch("https://x/y") == b"ok"


@pytest.mark.parametrize(
    ("url", "gap"), [("https://x/?a=1", 1 / 3), ("https://x/?api_key=k", 1 / 10)]
)
def test_requests_are_spaced_to_the_ncbi_limits(url: str, gap: float) -> None:
    """Three requests a second without a key, ten with one."""
    now = [100.0]
    slept: list[float] = []

    def sleep(seconds: float) -> None:
        slept.append(seconds)
        now[0] += seconds

    throttle = pubmed._Throttle(lambda: now[0], sleep)
    throttle.wait(url)
    assert slept == []
    throttle.wait(url)
    assert slept == [pytest.approx(gap)]
