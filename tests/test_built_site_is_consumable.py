# tests/test_built_site_is_consumable.py
import csv
import gzip
import html
import json
import math
import posixpath
import re
import shutil
import statistics
from collections import Counter
from html.parser import HTMLParser
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.burden import BurdenRow
from chd_atlas.build.concordance import survives_own_correction
from chd_atlas.build.render import RESEARCH_USE_NOTICE
from chd_atlas.build.runner import build_site
from chd_atlas.corpus import load_curation

REPO = Path(__file__).parent.parent


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One build of the real committed corpus, shared by every test here.

    The unit suites build from fixtures they control. This file exists to ask a
    different question — whether the data actually committed produces a site a
    frontend can consume — so it must build the real thing, and once.
    """
    base = tmp_path_factory.mktemp("consumable")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    out = base / "dist"
    build_site(root, out)
    return out


def _read(site: Path, relative: str) -> dict[str, Any]:
    payload: dict[str, Any] = json.loads(site.joinpath(*relative.split("/")).read_text())
    return payload


def _search_records(site: Path) -> list[dict[str, Any]]:
    raw = (site / "search" / "index.json.gz").read_bytes()
    records: list[dict[str, Any]] = json.loads(gzip.decompress(raw))["records"]
    return records


def _advertised(site: Path) -> list[tuple[str, str]]:
    """Every path any payload points a consumer at, paired with who advertised it."""
    references: list[tuple[str, str]] = []

    for relative in _read(site, "manifest.json")["files"]:
        references.append(("manifest.json", relative))

    index = _read(site, "genes/index.json")
    for entry in index["genes"]:
        references.append(("genes/index.json", entry["bundle"]))
        bundle = _read(site, entry["bundle"])
        for modality in bundle["omics"].values():
            for shard in modality["shards"]:
                references.append((entry["bundle"], shard))
        for dataset in bundle["expression_profile"]["datasets"]:
            # None on a corpus mid-curation, which advertises nothing.
            if dataset["quantile_shard"] is not None:
                references.append((entry["bundle"], dataset["quantile_shard"]))

    for shard in _read(site, "variants/index.json")["shards"]:
        references.append(("variants/index.json", shard))

    for record in _search_records(site):
        references.append(("search/index.json.gz", record["path"]))

    return references


def test_every_path_the_site_advertises_resolves(site: Path) -> None:
    """The dead-link sweep, across every payload that points at another.

    A dangling path here is the worst kind of failure this project has: the
    build is green, every checksum verifies, and a reader sees a gene with no
    data rather than a broken site. Nothing downstream can catch it —
    `Emitter.checksums` records what was written, never what a payload promised.

    This is the only test that checks the search index's `path` values against
    what was actually built. That matters because `search.py` names
    `publications.json` and `phenotypes.json` as string literals while
    `literature.py` is what writes them, with nothing tying the two together —
    amendment A36, whose structural fix is shared constants and is queued. Until
    then, this is the detection.

    The floors are asserted so the sweep cannot pass by finding nothing, which
    is how a loop over an empty collection reports success. They are floors
    rather than exact counts so that curating more data strengthens this test
    instead of breaking it. Two of the six sources contribute nothing today —
    no gene carries omics summaries and no variant shard exists in the committed
    corpus — so those loops are currently unexercised, and the floors are
    deliberately not raised to hide that. The bundles' floor is met by their
    `expression_profile` percentile grids alone, which reached no sweep of the
    JSON until 2026-10-26, although the gene pages link them.
    """
    references = _advertised(site)
    by_source: dict[str, int] = {}
    for source, _ in references:
        by_source[source] = by_source.get(source, 0) + 1

    assert by_source.get("manifest.json", 0) >= 8
    assert by_source.get("genes/index.json", 0) >= 1
    assert by_source.get("search/index.json.gz", 0) >= 3
    bundles = [source for source in by_source if source.startswith("genes/HGNC_")]
    assert sum(by_source[source] for source in bundles) >= 1

    for source, relative in references:
        assert site.joinpath(*relative.split("/")).is_file(), (
            f"{source} advertises {relative}, which was never written"
        )


def test_every_html_page_the_build_writes_carries_the_research_use_notice(site: Path) -> None:
    """Not one page kind. Every page, counted off the filesystem.

    `test_build_render.py` pins the sentence and pins that `document` emits it.
    This is the different question, and the one the defect actually was: *does
    every page the real build writes go through `document`?* Measured 2026-08-04
    on the build immediately before this change, the answer was no in effect --
    the notice lived in `landing.py`, so `grep -l` matched 1 of 25 files and the
    24 carrying gene-level classifications matched none.

    Enumerated with `rglob("*.html")` rather than against a list of expected
    paths, so a page kind added later is covered without anyone remembering to
    add it here -- which is the failure mode this replaces.

    **The count is asserted against the published population, not against a
    literal.** An `rglob` matching nothing would satisfy an all-pages loop
    vacuously, so a count is needed; a literal 25 was what the committed corpus
    published until the gate widened on 2026-08-06, and pinning the total meant
    re-editing this line for a change it has nothing to say about. Derived from
    the bundle count instead: one page per published gene, plus the landing page
    and the browse page. That fails on a gene page that was not written and on a
    stray page that was, which is what this guards, and it does not fail merely
    because the population moved.

    `+ 3` since 2026-10-26: `methods.html`, the site's reading notes (D50), is
    the third page that is not a gene's.
    """
    pages = sorted(site.rglob("*.html"))
    bundles = sorted((site / "genes").glob("HGNC_*.json"))
    expected = len(bundles) + 3

    assert len(bundles) > 1, "the build published no gene bundles; the fixture is broken"
    assert len(pages) == expected, (
        f"expected {expected} HTML pages (1 landing + 1 browse + 1 methods + "
        f"{len(bundles)} gene pages), found {len(pages)}"
    )
    for page in pages:
        assert RESEARCH_USE_NOTICE in page.read_text(encoding="utf-8"), (
            f"{page.relative_to(site)} carries no research-use notice"
        )


def test_no_advertised_path_would_leave_the_site_root(site: Path) -> None:
    """Every path in the output is a URL relative to the site root.

    An absolute path leaking from a build machine resolves on the machine that
    built it and 404s for every visitor, which is a failure only a reader ever
    sees. A `..` segment escapes the published tree entirely, and a backslash is
    a path separator on the filesystem a curator may be using but a literal
    character in a URL.

    Checked over what the payloads advertise rather than over the emitter's
    guard, because the two are different claims: `Emitter._write` refuses these
    shapes for files it writes, and this asserts that no payload *points* at one.
    """
    for source, relative in _advertised(site):
        assert not relative.startswith("/"), f"{source} advertises an absolute path {relative}"
        assert ".." not in relative.split("/"), f"{source} advertises an escaping path {relative}"
        assert "\\" not in relative, f"{source} advertises a backslash in {relative}"
        assert "://" not in relative, f"{source} advertises an off-site URL {relative}"


def test_the_seed_gene_and_its_manuscript_render_from_the_payloads_alone(site: Path) -> None:
    """What a frontend needs on first paint, from the data actually committed.

    Three separate promises, checked against the real corpus rather than a
    fixture: a browse row carries what it ranks on, a gene page is one fetch,
    and the landing page renders a featured manuscript without a second request.

    `featured[].publication` is asserted to be a resolved object, not a bare
    PMID. Amendment A24 records why: the specified fallback made it an object
    for resolvable entries and a string otherwise, so a consumer reading
    `entry.publication.title` got `undefined` with nothing raised — in the one
    payload the landing page renders without a further fetch.

    `headline_confidence` and `validity_state` read `"definitive"` and
    `"expert_curated"` here, not `null`/`"uncurated"`: the committed ClinGen
    mirror curates TBX5 Definitive for Holt-Oram syndrome
    (`MONDO:0007732`), `build_site` reads that mirror, and a gene the mirror
    plainly classifies must not publish as unassessed.
    """
    entry = next(
        item for item in _read(site, "genes/index.json")["genes"] if item["gene"] == "HGNC:11604"
    )
    assert entry["symbol"] == "TBX5"
    assert entry["headline_confidence"] == "definitive"
    assert entry["validity_state"] == "expert_curated"
    assert entry["has_conflicting_evidence"] is False

    bundle = _read(site, entry["bundle"])
    assert bundle["assertions"], "a gene page with no assertions is an empty page"
    assert bundle["publications"], "the page cites nothing"

    featured = _read(site, "featured.json")["featured"][0]
    assert isinstance(featured["publication"], dict), "a bare PMID needs a second fetch to render"
    assert featured["publication"]["title"]
    assert featured["blurb"]

    assert {record["kind"] for record in _search_records(site)} == {
        "gene",
        "publication",
        "phenotype",
    }


def test_tbx5_publishes_its_mirrored_validity_with_the_provenance_behind_it(
    site: Path,
) -> None:
    """The gene bundle's `validity` block, checked against the real mirrors.

    Not a fixture: `mirrors/clingen_gene_validity.tsv` curates TBX5 Definitive
    for Holt-Oram syndrome under SOP11, so the bundle's first record must carry
    that, and it must carry it *with its SOP version* — an attributed
    classification with no SOP is an unqualified claim, since the mirror spans
    SOP4 through SOP12 with no published crosswalk between them.

    `has_source_discordance` is asserted with `is`, not `==`: `{"f": 0} ==
    {"f": False}` is `True` in Python, so an equality alone would silently
    accept an `int` where the published contract is a `bool`.
    """
    bundle = _read(site, "genes/HGNC_11604.json")
    validity = bundle["validity"]

    assert validity["state"] == "expert_curated"
    assert validity["has_source_discordance"] is False

    clingen = validity["records"][0]
    assert clingen["source"] == "clingen"
    for field in ("source", "classification_term", "disease", "sop", "report_url"):
        assert clingen[field], f"{field} is empty on TBX5's ClinGen record"

    # R8 (owner decision 2026-10-26): the bundle keeps the URL GenCC gave,
    # PanelApp Australia's bare homepage, and the page links TBX5's own entry
    # there instead -- the one report URL the atlas builds.
    panelapp = [r for r in validity["records"] if r["submitter"] == "PanelApp Australia"]
    assert [r["report_url"] for r in panelapp] == ["https://panelapp-aus.org"]
    page = (site / "genes" / "HGNC_11604.html").read_text(encoding="utf-8")
    assert re.findall(r'href="(https://panelapp-aus\.org[^"]*)"', page) == [
        "https://panelapp-aus.org/panels/entities/TBX5"
    ]


def _burden_rows(gene: str, bundle: dict[str, Any]) -> list[BurdenRow]:
    """A bundle's published burden rows as the `BurdenRow`s the renderers read."""
    return [
        BurdenRow(
            **{
                **row,
                "gene": gene,
                "case_cohorts": tuple(row["case_cohorts"]),
                "control_cohorts": tuple(row["control_cohorts"]),
            }
        )
        for row in bundle["burden"]
    ]


def test_tbx5s_card_says_how_wide_each_surviving_result_is(site: Path) -> None:
    """R3 (owner decision 2026-10-26) on the real corpus, pinned on TBX5.

    PMID:40127276's q is a Benjamini-Hochberg correction over its 248 panel
    genes, and the study states its own genome-wide line at p < 2.6e-6. TBX5
    carries both kinds of survivor: its case-control damaging row (p 1.34e-4,
    q 0.00369) survives within the panel only, and its de novo damaging row (p
    1.07e-8) and loss-of-function row (p 1.35e-9) clear the genome-wide line
    too. Literals, never the constants that render them.

    **Every** surviving row below the line is marked on its card, the named
    ones as well as each line's own: derived from every bundle with the
    predicate the card reads, never pinned. Marking only each line's row left
    15 of the 25 such rows unmarked when measured (2026-10-26), and 6 of the 13
    genes carrying one, FLT4's case-control loss-of-function row (p 1.75e-15)
    among them. Both sides of that count must be non-zero for it to measure
    anything (CLAUDE.md section 4.30), so the line's own and the named marks are
    each required to occur.

    W1 is checked here too: PMID:40127276's `tests_reported` is null, because
    the 248 it held counted panel genes, not comparisons.
    """
    publications = {
        publication["id"]: publication
        for publication in _read(site, "publications.json")["publications"]
    }
    sierant = publications["PMID:40127276"]
    assert sierant["correction_genes"] == 248
    assert sierant["genome_wide_threshold"] == 2.6e-06
    assert sierant["tests_reported"] is None

    page = (site / "genes" / "HGNC_11604.html").read_text(encoding="utf-8")
    card = page[page.index('<div class="card-variants">') :]
    card = card[: card.index("</div>")]
    assert re.findall(r"<li>(.*?) · <a ", card) == [
        "Sierant MC et al. 2025 · case-control burden test — damaging (LOF + missense): "
        "22 alleles / 21,750 vs 103 alleles / 270,733 · q 0.00369 within its 248-gene panel "
        "· also surviving: loss-of-function",
        "Sierant MC et al. 2025 · de novo enrichment (trios vs mutation model) — damaging "
        "(LOF + missense): 5 de novo / 3,887 trios vs 0.0671 expected · q 5.31e-07 within its "
        "248-gene panel · genome-wide (p 1.07e-08) · also surviving: loss-of-function "
        "(genome-wide, p 1.35e-09)",
    ]

    thresholds = {
        study: publication["genome_wide_threshold"]
        for study, publication in publications.items()
        if publication["genome_wide_threshold"] is not None
    }
    expected = own = named = 0
    for entry in _read(site, "genes/index.json")["genes"]:
        rows = _burden_rows(entry["gene"], _read(site, entry["bundle"]))
        expected += sum(
            1
            for row in rows
            if row.study in thresholds
            and survives_own_correction(row)
            and row.pvalue is not None
            and row.pvalue < thresholds[row.study]
        )
        html_page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        rail = html_page[html_page.index('<aside class="rail">') : html_page.index("</aside>")]
        own += rail.count(" · genome-wide (p ")
        named += rail.count(" (genome-wide, p ")
    assert own > 0 and named > 0, (own, named)
    assert own + named == expected


def test_studies_on_collections_declared_overlapping_are_one_dataset_wherever_counted(
    site: Path,
) -> None:
    """Owner decision 2026-10-27 (I3), on the real corpus, read from the payloads.

    `cohorts.json` declares DDD and DECIPHER overlapping, so PMID:42230622
    (`ddd`) and PMID:34324492 (`decipher`) are one cohort family. Every
    consumer of families must agree (CLAUDE.md section 4.28): the manifest's
    and the front page's dataset count, every bundle's `families`, and the
    notice on every page whose rows cite both -- derived here from the rows
    and the declarations, never pinned per page. Measured on the build it was
    written against: 3 families became 2, and 72 pages gained the notice.
    """
    cohorts = _read(site, "cohorts.json")["cohorts"]
    declared = {(cohort["id"], other) for cohort in cohorts for other in cohort["overlaps_with"]}
    assert declared == {("ddd", "decipher"), ("decipher", "ddd")}
    families = _read(site, "manifest.json")["counts"]["cohort_families"]
    assert families == 2
    landing = (site / "index.html").read_text(encoding="utf-8")
    # The hero figure, with the qualifier owner decision R12 put beside it,
    # since the figure counts only what the registry records.
    assert re.search(
        rf'"figure-value">{families}</span>\s*<span class="figure-label">independent datasets'
        r'</span>\s*<span class="figure-note">\(as far as the cohort registry records\)</span>',
        landing,
    )

    noticed = 0
    for entry in _read(site, "genes/index.json")["genes"]:
        bundle = _read(site, entry["bundle"])
        cited: dict[str, set[str]] = {}
        for row in bundle["burden"]:
            cited.setdefault(row["study"], set()).update(
                row["case_cohorts"], row["control_cohorts"]
            )
        studies = sorted(cited)
        linked = [
            (left, right)
            for i, left in enumerate(studies)
            for right in studies[i + 1 :]
            if cited[left] & cited[right]
            or any((a, b) in declared for a in cited[left] for b in cited[right])
        ]
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        assert page.count("These two studies are not independent.") == len(linked), entry["gene"]
        noticed += bool(linked)
        assert len(bundle["independent_datasets"]["families"]) == families, entry["gene"]
        grouped = [set(family["studies"]) for family in bundle["independent_datasets"]["families"]]
        for left, right in linked:
            assert any({left, right} <= group for group in grouped), (entry["gene"], left, right)
    assert noticed == 72


# One card per form of the tested line that occurs on the real build, measured
# the day the qualifier reached the cards (owner answer 2026-10-27): 24, 4, 48,
# 7 and 3 cards, 86 in all.
_TESTED_LINES: dict[str, str] = {
    "HGNC_10252": "Tested in 1 independent dataset (as far as the cohort registry records) (a "
    "study by an atlas author).",
    "HGNC_13672": "Tested in 1 independent dataset (as far as the cohort registry records) (a "
    "study by an atlas author); none survives a correction its study published.",
    "HGNC_11291": "Tested in 1 independent dataset (as far as the cohort registry records) "
    "(each study by an atlas author); none survives a correction its study published.",
    "HGNC_6770": "Tested in 2 independent datasets (as far as the cohort registry records) (1 "
    "of the 2 studies by an atlas author); none survives a correction its study published.",
    "HGNC_7882": "Tested in 2 independent datasets (as far as the cohort registry records) (2 "
    "of the 3 studies by an atlas author); none survives a correction its study published.",
}


def test_every_dataset_count_on_a_card_says_what_the_registry_can_see(site: Path) -> None:
    """The front page's qualifier on the gene cards too (owner answer 2026-10-27).

    "Independent datasets 2 (as far as the cohort registry records)" on the
    front page, and on a card "Tested in N independent dataset(s)" with no
    qualifier, though both count the same cohort families: joined only through
    collections the registry names or declares overlapping. The qualifier
    follows the count and precedes the author declaration, which it must not
    read as qualifying. One card per form that occurs is pinned as a literal
    -- ROCK2, KLF13, SRF, SMAD4, NOTCH2 -- and every tested line on every card
    is held to the shape, so a form that appears later carries it too.
    """
    seen = {}
    for path in sorted((site / "genes").glob("HGNC_*.html")):
        page = path.read_text(encoding="utf-8")
        rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
        for line in re.findall(r"<p>(Tested in [^<]*)</p>", rail):
            assert re.fullmatch(
                r"Tested in (1 independent dataset|\d+ independent datasets) "
                r"\(as far as the cohort registry records\)( \([^)]*\))?"
                r"(; none survives a correction its study published)?\.",
                line,
            ), (path.stem, line)
            seen[path.stem] = line
    assert seen, "no card prints a dataset count; this checks nothing"
    assert {stem: seen.get(stem) for stem in _TESTED_LINES} == _TESTED_LINES


def test_every_published_bundle_carries_the_identity_its_registry_row_holds(site: Path) -> None:
    """D52 on the committed corpus. A runner that stops passing `identities` to
    `build_genes` publishes six nulls on every bundle -- which the API doc defines
    as "the registry holds no value" -- with every unit test green, because those
    pass `identities` themselves. This compares every published bundle against
    `mirrors/genes.tsv`. The count is derived, never pinned.

    Compared after stripping, as `tables.cell_text` strips: HGNC's own name for
    HGNC:2348 ends in a space ("CREB binding lysine acetyltransferase "), the one
    cell of 1,850 in the mirror with surrounding whitespace (measured 2026-10-25).
    """
    with (REPO / "mirrors" / "genes.tsv").open(encoding="utf-8") as handle:
        registry = {row["hgnc_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    entries = _read(site, "genes/index.json")["genes"]
    assert len(entries) > 1, "the build published no genes; the fixture is broken"
    for entry in entries:
        row = registry[entry["gene"]]
        identity = _read(site, entry["bundle"])["identity"]
        for key in ("name", "locus", "ensembl_gene", "uniprot", "mane_select"):
            assert identity[key] == (row[key].strip() or None), (entry["gene"], key)
        assert identity["ncbi_gene"] == int(row["ncbi_gene"]), entry["gene"]
        assert type(identity["ncbi_gene"]) is int, entry["gene"]


def _page_links(site: Path) -> list[tuple[str, str, str, str]]:
    """Every internal `href` on every built page, with the page that carries it.

    Returns `(page, href, resolved, fragment)` where `resolved` is the
    site-relative path the browser would request and `fragment` is what follows
    a `#`, or `""`. External schemes and pure fragments are skipped: this sweep
    is about links the build is responsible for, and a pure fragment is
    `test_every_fragment_link_on_a_gene_page_lands_on_one_element_of_that_page`'s.

    The fragment is split off before resolving, because a browser requests the
    file and scrolls to the fragment: `../methods.html#scope` is a request for
    `methods.html`. Resolved with the fragment attached it would name a file no
    build could write, and every short form linking `methods.html` (D50) would
    read as broken.
    """
    links: list[tuple[str, str, str, str]] = []
    for path in sorted(site.rglob("*.html")):
        page = path.relative_to(site).as_posix()
        for href in re.findall(r'href="([^"]+)"', path.read_text(encoding="utf-8")):
            if href.startswith(("http://", "https://", "mailto:", "#")):
                continue
            file, _, fragment = href.partition("#")
            target = posixpath.normpath(posixpath.join(posixpath.dirname(page), file))
            links.append((page, href, target, fragment))
    return links


def test_every_link_on_every_page_resolves_to_a_file_the_build_wrote(site: Path) -> None:
    """The HTML half of the sweep above, which had no HTML half at all.

    `_advertised` walks manifest entries, gene bundles, omics and variant shards
    and search `path`s -- every one of them JSON. It never opens an `.html`
    file. This release added 25 pages and, measured here, 100+ internal links,
    none of which anything checked.

    This project has already shipped exactly this failure once: `runner.py`
    records that the nav and the landing page both linked `genes/index.html`
    while no builder wrote one -- green build, verifying checksums, 404 for
    every visitor. It was fixed by hand and no sweep was added.

    Measured before this test existed: changing the browse row's href from
    `f"../{gene_page_path(...)}"` to `gene_page_path(...)` -- which points all
    23 rows at `genes/genes/HGNC_*.html` -- left the whole 624-test suite green.
    `test_every_browse_row_links_to_a_page_that_was_written` could not catch it
    because it asserts a bare substring (`name in page`) and a key in
    `emitter.checksums`, neither of which sees the `../` prefix.

    Resolution is `posixpath.normpath` against the page's own directory, so a
    link that escapes the site root normalises to a path starting `..` and fails
    the membership test rather than silently resolving somewhere on disk.

    **A fragment must land too (2026-10-26, D50).** Every gene page and the
    gene index link `../methods.html#scope`. A fragment that names no element
    opens the right file at the top -- a reader sent to "how a gene is
    admitted" lands on the page title and has to hunt, and nothing reports it.
    So each cross-page fragment must name exactly one `id` in the page it
    targets: absent is a link that scrolls nowhere, repeated is a link whose
    target the browser picks for you.
    """
    written = {path.relative_to(site).as_posix() for path in site.rglob("*") if path.is_file()}
    links = _page_links(site)

    assert links, "no page links found; the sweep would pass vacuously"

    def served(target: str) -> str | None:
        # A directory URL is served by its `index.html`, which is how
        # `genes/index.html` is reached as `genes/`.
        for candidate in (target, f"{target}/index.html"):
            if candidate in written:
                return candidate
        return None

    broken = [(page, href) for page, href, target, _ in links if served(target) is None]
    assert not broken, f"pages link to files the build never wrote: {broken}"

    ids: dict[str, list[str]] = {}
    fragments = 0
    for page, href, target, fragment in links:
        if not fragment:
            continue
        fragments += 1
        file = served(target)
        assert file is not None
        if file not in ids:
            ids[file] = re.findall(r' id="([^"]+)"', (site / file).read_text(encoding="utf-8"))
        found = ids[file].count(fragment)
        assert found == 1, f"{page} links {href}, and {file} carries that id {found} time(s)"
    assert fragments, "no page links a fragment on another page; the check measures nothing"


def test_every_cohort_a_published_burden_row_cites_resolves_with_its_caveat(site: Path) -> None:
    """A bare cohort id in a bundle must resolve to a description in `cohorts.json`.

    **This is the join the burden columns exist for.** A row names its sample
    collections by bare id -- `["taa_cases"]`, `["gnomad_controls"]` -- and the
    ids carry no meaning on their own. Until schema 2.9 nothing published mapped
    them, so a consumer had the statistics and none of the caveats: measured
    2026-08-06 over a real build, 13 distinct ids appeared across 915 published
    rows and `curation/cohorts.yaml` reached zero published bytes.

    Asserted over a *real* build rather than a fixture, because what could break
    it is not the emitter -- covered by
    `test_cohorts_json_publishes_the_registry_a_burden_row_resolves_against` --
    but the publication gate moving under it. A widened gate can admit a gene
    whose rows cite a collection nobody curated, and that is a dangling id in a
    published payload rather than a build failure. `BUR009` refuses it at the
    gate; this checks the guarantee survives to the bytes.

    `taa_cases` is named explicitly. It is 777 thoracic aortic aneurysm probands
    in a case set the study reports as CHD, and it is the id whose description
    most changes what a row means -- a reader who resolves it and a reader who
    does not are reading different data.

    **Its pinned phrase changed on 2026-10-26.** The description said "TAA IS
    NOT CONGENITAL HEART DISEASE"; it now says thoracic aortic aneurysm is
    outside this atlas's scope and that whether any of these probands also has
    a congenital heart lesion is not recorded here -- the study's abstract
    describes them as "probands with CHD including cases with sporadic
    thoracic aortic aneurysm", so the atlas cannot assert they have none. The
    caveat a reader needs is the scope one, and that is what is pinned.
    """
    registry = {
        cohort["id"]: cohort
        for cohort in json.loads((site / "cohorts.json").read_text(encoding="utf-8"))["cohorts"]
    }
    assert registry, "the build published no cohort registry"

    cited: set[str] = set()
    for bundle in sorted((site / "genes").glob("HGNC_*.json")):
        for row in json.loads(bundle.read_text(encoding="utf-8"))["burden"]:
            cited.update(row["case_cohorts"])
            cited.update(row["control_cohorts"])

    assert cited, "no published burden row names a cohort, so this proves nothing"
    dangling = sorted(cited - set(registry))
    assert not dangling, f"published rows cite cohorts that resolve to nothing: {dangling}"

    # Every resolved record carries the sentence a reader needs, not just a name.
    for identifier in sorted(cited):
        assert registry[identifier]["description"].strip(), f"{identifier} publishes no caveat"

    assert "outside this atlas's scope" in registry["taa_cases"]["description"]


def test_methods_html_lists_every_cohort_cohorts_json_publishes(site: Path) -> None:
    """The registry for a reader, `cohorts.json` for a program: the same records (P4).

    The unit test pins the markup on a fixture; this pins the wiring and the
    data -- that `runner.build_site` hands `build_methods` the real registry,
    and that every record reaches the page whole. Equality on the id list, so
    a cohort missing, listed twice or out of `cohorts.json`'s order (sorted by
    id) fails; then each name, description and url, escaped as the page
    escapes them, so a description cut short or a url dropped fails too. A
    record's open question (owner's standing rule, 2026-10-27) follows its
    description, before the link, where `cohorts.json` carries one.
    """
    cohorts = json.loads((site / "cohorts.json").read_text(encoding="utf-8"))["cohorts"]
    assert cohorts, "the build published no cohort registry; this measures nothing"
    page = (site / "methods.html").read_text(encoding="utf-8")
    start = page.index('<h2 id="cohorts">')
    ends = [i for i in (page.find("<h2", start + 1), page.find("</main>", start)) if i != -1]
    section = page[start : min(ends)]

    listed = re.findall(r"<dt>(.*?) \(<code>([^<]*)</code>\)</dt><dd>(.*?)</dd>", section, re.S)
    assert [identifier for _, identifier, _ in listed] == [cohort["id"] for cohort in cohorts]
    for (name, _, entry), cohort in zip(listed, cohorts, strict=True):
        assert name == html.escape(cohort["name"])
        expected = html.escape(cohort["description"])
        if cohort["review_status"] is not None:
            label = _REVIEW_LABELS[cohort["review_status"]]
            expected += f" <strong>{label}:</strong> {html.escape(cohort['review_note'])}"
        if cohort["url"] is not None:
            expected += f' <a href="{html.escape(cohort["url"])}">website</a>'
        assert entry == expected, f"{cohort['id']} is not listed whole"


# How a page labels each review status, as a literal (CLAUDE.md section 4.38).
_REVIEW_LABELS: dict[str, str] = {
    "pending_review": "Pending review",
    "not_accessible": "Source not accessible",
}

_SYNDROME_QUESTION = (
    "Whether this term names a syndrome -- which decides whether a gene page calls a grade "
    "for it a grade for a syndrome as a whole -- is left open for clinical review."
)


def test_exactly_the_four_records_the_owner_marked_carry_an_open_question(site: Path) -> None:
    """The owner's answer of 2026-10-27, on the real corpus, read from what is published.

    Four records carry `pending_review`, each note stating what is
    unresolved and taking no side. Two are S2 case series where PMID:34324492's
    S2 Table and the cited paper's PubMed abstract conflict (both read
    2026-10-27): S2's Affymetrix and Agilent arrays and 432 cases against
    PMID:22939634's Illumina 660W-Quad arrays in 2,256 individuals with CHD;
    S2's 178 individuals with CHD against PMID:22969434's "464 individuals
    with 174 affecteds (37.5%) in 59 multiplex families and 8 trios". Each
    note names the PMID whose abstract it quotes: "the paper's abstract", after
    a description naming PMID:34324492, read as that study's (review,
    2026-10-27). The other five series S2 identifies by PMID list no more
    individuals with CHD than their papers' abstracts report, on no array an
    abstract contradicts -- compared with the abstracts only, never the full
    texts -- which is no conflict, and carry none. Two are scope terms whose
    syndromic status the owner left open for clinical review, a clinical
    judgement rather than a conflict between sources; their `syndrome` mark
    stays unset.

    Pinned whole, as literals: every other record publishes both fields as
    `null`, never omits them. A cohort's question renders where its
    description does -- `methods.html`'s list, checked with every record by
    `test_methods_html_lists_every_cohort_cohorts_json_publishes`, and the
    cohort notes of every gene page whose rows cite it, here: the 77 pages
    carrying PMID:34324492 on the build this was written against (every one
    of its rows names all seven series), each note once, inside the fold and
    nowhere else on the page. A scope term's question is published in
    `scope.json` only, so its note is asserted absent from every page.
    """
    cohorts = {cohort["id"]: cohort for cohort in _read(site, "cohorts.json")["cohorts"]}
    terms = {term["id"]: term for term in _read(site, "scope.json")["terms"]}
    records = {**cohorts, **terms}
    for identifier, record in records.items():
        assert "review_status" in record and "review_note" in record, identifier
    marked = {
        identifier: (record["review_status"], record["review_note"])
        for identifier, record in records.items()
        if record["review_status"] is not None or record["review_note"] is not None
    }
    assert marked == {
        "pmid22939634_cases": (
            "pending_review",
            "S2 gives the platform as Affymetrix 6.0 array, Agilent CGH 44k/105K and lists 432 "
            "cases; PMID:22939634's abstract reports Illumina 660W-Quad SNP arrays in 2,256 "
            "individuals with CHD. How S2's 432 relate to the individuals that abstract "
            "reports, and which platform applies to them, is not stated.",
        ),
        "pmid22969434_cases": (
            "pending_review",
            "S2 lists 178 individuals with CHD from this source; PMID:22969434's abstract "
            'reports "464 individuals with 174 affecteds (37.5%) in 59 multiplex families and '
            "8 trios\". How S2's 178 relate to the individuals that abstract reports is not "
            "stated.",
        ),
        "MONDO:1010098": ("pending_review", _SYNDROME_QUESTION),
        "MONDO:1060150": ("pending_review", _SYNDROME_QUESTION),
    }
    assert terms["MONDO:1010098"]["syndrome"] is False
    assert terms["MONDO:1060150"]["syndrome"] is False

    # GATA4's page, sliced to the cohort-notes fold and then to the item.
    page = (site / "genes" / "HGNC_4173.html").read_text(encoding="utf-8")
    fold = page[page.index('<details class="cohort-notes">') :]
    fold = fold[: fold.index("</details>")]
    start = fold.index("<li><strong>CHD case series of PMID:22969434</strong>")
    assert fold[start : fold.index("</li>", start) + 5] == (
        "<li><strong>CHD case series of PMID:22969434</strong> &mdash; 178 individuals with "
        "congenital heart disease, listed in PMID:34324492&#x27;s S2 Table as a source of its "
        "CNV case set; S2 gives the platform as Affymetrix Array 6.0. <strong>Pending "
        "review:</strong> S2 lists 178 individuals with CHD from this source; "
        "PMID:22969434&#x27;s abstract reports &quot;464 individuals with 174 affecteds "
        "(37.5%) in 59 multiplex families and 8 trios&quot;. How S2&#x27;s 178 relate to the "
        "individuals that abstract reports is not stated.</li>"
    )

    carrying = 0
    for entry in _read(site, "genes/index.json")["genes"]:
        bundle = _read(site, entry["bundle"])
        cited = {
            c for row in bundle["burden"] for c in (*row["case_cohorts"], *row["control_cohorts"])
        }
        expected = sum(
            1 for cohort in ("pmid22939634_cases", "pmid22969434_cases") if cohort in cited
        )
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        notes = ""
        if '<details class="cohort-notes">' in page:
            notes = page[page.index('<details class="cohort-notes">') :]
            notes = notes[: notes.index("</details>")]
        assert notes.count("<strong>Pending review:</strong>") == expected, entry["gene"]
        assert page.count("Pending review") == expected, entry["gene"]
        carrying += bool(expected)
    assert carrying == 77
    for path in sorted(site.rglob("*.html")):
        text = path.read_text(encoding="utf-8")
        for question in (_SYNDROME_QUESTION, html.escape(_SYNDROME_QUESTION)):
            assert question not in text, path.name


def test_the_taa_description_states_the_denominator_its_rows_publish(site: Path) -> None:
    """The `taa_cases` caveat, measured against the rows it qualifies (2026-10-26).

    It said "Any row naming this cohort has a case set that is about 10%
    non-CHD (777 of 7,958), and the page must say so" -- a curator's
    instruction published as a caveat, beside rows whose case set is not
    7,958. S2 puts this series at 777 of the 7,958; every row naming the
    cohort is a PMID:34324492 deletion row (121 mirror rows), and those rows
    publish the 6,746 cases and 14,024 controls the article's methods say the
    permutation test was run on, after filtering. So the description gives
    both figures, and this asserts the second pair against every published row
    that cites the cohort: a re-curated mirror with other denominators fails
    here rather than leaving the caveat contradicting the table beside it.

    **The denominator changed the same day.** The rows published S1 Table's
    deletion column, 7,107 cases and 14,040 controls, and the description said
    a row "publishes the 7,107 cases the study's S1 Table gives for its
    deletion analysis" until the filtered set was settled as the one the
    statistic was computed on. That clause is asserted absent beside the other
    retired phrases, so it cannot come back beside the new one.

    The retired phrases are asserted absent from every published byte that
    carries the description -- `cohorts.json`, `methods.html` and each gene
    page's cohort notes -- both apostrophe spellings (section 4.35).
    """
    registry = {
        cohort["id"]: cohort
        for cohort in json.loads((site / "cohorts.json").read_text(encoding="utf-8"))["cohorts"]
    }
    description = registry["taa_cases"]["description"]
    denominators = {
        (row["n_cases"], row["n_controls"])
        for bundle in sorted((site / "genes").glob("HGNC_*.json"))
        for row in json.loads(bundle.read_text(encoding="utf-8"))["burden"]
        if "taa_cases" in row["case_cohorts"]
    }
    assert denominators, "no published row cites taa_cases; this measures nothing"
    for cases, controls in sorted(denominators):
        assert f"is computed on the {cases:,} cases and {controls:,} controls" in description, (
            f"rows citing taa_cases publish {cases:,} cases and {controls:,} controls; "
            "the caveat does not say so"
        )
    assert "777 of its 7,958 cases" in description

    published = [site / "cohorts.json", site / "methods.html", *sorted(site.glob("genes/*.html"))]
    for path in published:
        text = path.read_text(encoding="utf-8")
        for retired in (
            "must say so",
            "explicitly out of",
            "puts aortopathy outside",
            "S1 Table gives for its deletion analysis",
        ):
            assert retired not in text and retired not in html.unescape(text), (
                f"{retired!r} is published in {path.name}"
            )


def test_a_gene_an_authority_reported_no_association_for_says_so_in_both_payloads(
    site: Path,
) -> None:
    """Issue #13's third axis, on the one published gene that exercises it.

    GDF1 carries G2P `Definitive`, Labcorp `Strong` and Illumina `No Known
    Disease Relationship` in scope. Until schema 2.10 it published
    `has_conflicting_evidence: false` and nothing else, so a consumer was told
    the evidence did not conflict while two bodies disagreed about whether an
    association exists at all. `no_known_association` takes neither side of
    `CONTESTED` -- correctly, since a null result is not a refutation -- and the
    disagreement therefore reached no published byte.

    Asserted on the bundle **and** the browse row, because a flag a reader meets
    on one and not the other is missing exactly where they are standing, and
    asserted against a real build rather than a fixture because what breaks it is
    the gate moving: GDF1 is admitted on GenCC agreement, and it is the only
    published gene on this axis, so a widening or narrowing that drops it leaves
    the axis published on nothing while every unit test still passes.

    `has_conflicting_evidence` is asserted **false** in the same breath. The
    whole design decision is that these are two axes rather than one, and a test
    that checked only the new flag would pass with the two collapsed together.
    """
    bundle = json.loads((site / "genes" / "HGNC_4214.json").read_text(encoding="utf-8"))
    assert bundle["symbol"] == "GDF1", "this test is about GDF1"

    assert bundle["has_no_association_report"] is True
    assert bundle["no_association_reported_by"] == ["Illumina"]
    # The other axis stays untouched: a null result is not a contest.
    assert bundle["has_conflicting_evidence"] is False

    row = next(
        item
        for item in json.loads((site / "genes" / "index.json").read_text(encoding="utf-8"))["genes"]
        if item["gene"] == "HGNC:4214"
    )
    assert row["has_no_association_report"] is True
    assert row["no_association_reported_by"] == ["Illumina"]
    assert row["has_conflicting_evidence"] is False

    # The pair is consistent everywhere, not only on this gene: an authority
    # named while the flag reads false would say two things at once.
    for bundle_path in sorted((site / "genes").glob("HGNC_*.json")):
        payload = json.loads(bundle_path.read_text(encoding="utf-8"))
        assert payload["has_no_association_report"] == bool(payload["no_association_reported_by"])
        # And the card says it exactly where the bundle does (owner decision
        # 2026-10-26, R6): one line, naming every authority behind a no-known-
        # disease-relationship record, as the card names them -- a panel record
        # by its panel, a GenCC one by its submitter. **Read off the validity
        # records, not `no_association_reported_by`** (review, 2026-10-27): the
        # bundle names a panel record's GenCC copy "ClinGen", and the card folds
        # that copy into the panel's own clause, so the bundle's list fails here
        # the day a published gene carries a panel's record (RXRA has one, and
        # is not published). A copy folds where the panel holds a record for
        # the same disease; one that differs from it is listed in the table,
        # and the card may name it, which this does not require.
        rail = bundle_path.with_suffix(".html").read_text(encoding="utf-8")
        rail = rail[rail.index('<aside class="rail">') : rail.index("</aside>")]
        lines = re.findall(r'<p class="card-no-association">(.*?)</p>', rail)
        assert len(lines) == int(payload["has_no_association_report"]), payload["gene"]
        if not lines:
            continue
        silent = [
            record
            for record in payload["validity"]["records"]
            if record["classification"] == "no_known_association"
        ]
        panel = {record["disease"] for record in silent if record["source"] == "clingen"}
        authorities = {
            (f"the {record['gcep']}" if record["gcep"] else "ClinGen")
            if record["source"] == "clingen"
            else record["submitter"]
            for record in silent
            if not (record["submitter"] == "ClinGen" and record["disease"] in panel)
        }
        assert authorities, payload["gene"]
        for authority in sorted(authorities):
            clause = f"{html.escape(authority)} reports no known disease relationship"
            assert clause in lines[0] or clause[0].upper() + clause[1:] in lines[0], authority
    page = (site / "genes" / "HGNC_4214.html").read_text(encoding="utf-8")
    assert (
        '<p class="card-no-association">Illumina reports no known disease relationship for '
        "congenital heart defects, multiple types, 6 (2021); see the validity table.</p>"
    ) in page


# Every attribute the chart renderers position anything with. `points` and
# `viewBox` hold several numbers each and are split on whitespace and commas.
_GEOMETRY_ATTRIBUTES = (
    "data-scale-high",
    "viewBox",
    "height",
    "points",
    "width",
    "x1",
    "x2",
    "y1",
    "y2",
    "cx",
    "cy",
    "x",
    "y",
    "r",
)
# Longest-first. Measured: the order is not load-bearing -- the alternation
# backtracks, so `x` failing on the `1` of `x1="..."` lets `x1` match anyway --
# but written this way so the scan does not depend on that.
_GEOMETRY = re.compile(r"\b(?:" + "|".join(_GEOMETRY_ATTRIBUTES) + r')="([^"]*)"')

# What `charts.coordinate` produces: an integer, or an integer and exactly one
# decimal place. A raw `repr` reaches `0.30000000000000004` and `1e-05`.
_FIXED_PRECISION = re.compile(r"^-?\d+(\.\d)?$")


def test_every_geometry_the_pages_publish_has_coordinate_s_fixed_shape(site: Path) -> None:
    """`coordinate` is pinned; that it is *used* was not, until 2026-08-20.

    `test_build_charts.py::test_coordinates_are_fixed_precision_so_two_builds_
    agree` proves the function rounds. It cannot prove any renderer calls it,
    so a new chart interpolating a raw float -- the ordinary way to write an
    SVG attribute -- passes every test in this repository while publishing
    `y="43.33333333333333"`. That is the whole reason `coordinate` exists:
    binary floating point is not identical across platforms at the last bits,
    so a raw repr makes two builds of one commit differ by architecture while
    every checksum still verifies against itself.

    Scanned off the built pages rather than off the source, for the reason
    CLAUDE.md section 4.35 gives: the question is what a reader's browser
    receives, and `grep` over `src/` answers a different one -- these strings
    are assembled from adjacent literals and f-strings that no source-level
    search reconstructs.

    Measured 2026-08-20 over a real build: 94 pages, 41,833 geometry tokens,
    0 non-conforming.

    **The count is derived, never pinned to a literal**, for the reason
    `test_every_html_page_the_build_writes_carries_the_research_use_notice`
    derives its own: a scan matching nothing satisfies a loop vacuously, and a
    literal total is re-edited every time the population moves. Every gene page
    must carry geometry -- measured, all 92 do, though not all for the same
    reason: GDF1 draws three forest panels and no organ chart at all, and
    CFC1's 138 tokens are the fewest. The landing and browse pages carry none
    today and are scanned anyway, so a chart added to either is covered
    without anyone remembering this test.
    """
    pages = sorted(site.rglob("*.html"))
    gene_pages = sorted((site / "genes").glob("HGNC_*.html"))
    assert len(gene_pages) > 1, "the build published no gene pages; the fixture is broken"

    total = 0
    for page in pages:
        tokens = [
            token
            for value in _GEOMETRY.findall(page.read_text(encoding="utf-8"))
            for token in value.replace(",", " ").split()
        ]
        if page in gene_pages:
            assert tokens, f"{page.relative_to(site)} publishes no geometry at all"
        total += len(tokens)
        for token in tokens:
            assert _FIXED_PRECISION.match(token), (
                f"{page.relative_to(site)} publishes the geometry {token!r}, which did not "
                "come through charts.coordinate -- two builds on two platforms can differ"
            )

    assert total > len(gene_pages), f"only {total} geometry tokens over {len(pages)} pages"


# Each `<rect class="chart-band">` a trajectory draws, paired with the phase
# label its own `<title>` names. Anchored on the `<rect>` rather than on the
# bare class name, because `render.STYLESHEET` is inlined verbatim into every
# page and carries `.chart-band { ... }` — a `"chart-band" in page` check
# passes on all 94 pages with every band deleted (CLAUDE.md section 4.35).
_BAND = re.compile(r'<rect class="chart-band"[^>]*><title>([^<]*)</title>')


def test_the_curated_phase_vocabulary_reaches_the_pages_it_is_drawn_on(site: Path) -> None:
    """The wiring, not the renderer: `build_gene_pages(phases=...)` in `runner.py`.

    Every page test hands `phases=` a fixture explicitly, so the keyword
    `runner.py` actually passes was exercised by nothing. Measured 2026-08-20,
    changing that one argument to `None` passed all 1,053 tests while changing
    86 of 198 published files: **85 pages lost their trajectory entirely** and
    30 of them gained, in bold, "This gene reads *below this dataset's
    detection floor at every stage sampled* in whole heart (19 stages
    sampled)" — including HGNC:10249, whose own bundle places 18 of those 19
    heart stages against the percentile grid with the percentiles printed in
    the table directly below. A false claim, on 30 pages, in bold, under a
    green build. This repository's characteristic defect exactly.

    Three assertions, because the wiring can fail in three ways and only the
    first is the mutant above:

    * bands are drawn at all — `phases=None` publishes zero;
    * every label drawn is one a curator wrote in
      `curation/cardiac_phases.yaml` — a renderer inventing its own phase
      names would satisfy the first;
    * the phases banded are **exactly** those whose window a sampled stage
      falls inside — not a subset, not a superset. This is the only check that
      the renderer and `CardiacPhaseFile.phases_for` agree over the real
      corpus rather than over a fixture, and unlike a subset assertion it
      catches a band that should be drawn and is not.

    **The third assertion used to be "the one phase whose end the source never
    states is banded nowhere", and it retired on 2026-08-24 rather than being
    quietly rewritten.** `heart_looping` was that phase; it now carries a
    `derived` end from GO:0001947's own definition, so **no curated phase is
    unended** and the assertion had nothing left to range over. It failed
    loudly rather than passing vacuously, because whoever wrote it put
    `assert unended, "the third assertion below is vacuous"` above it — and
    that guard is the only reason this is a paragraph rather than a dead check
    nobody noticed. The rule itself is unchanged and still tested, at the layer
    where a fixture can build the state the real pipeline no longer produces:
    `test_build_pages.py::test_a_phase_the_source_never_ended_bands_nothing`.
    What is genuinely lost is real-pipeline confirmation of that rule, and it
    is lost because the corpus stopped containing an instance — not because
    anyone judged it less important.

    Measured 2026-08-24 on the committed corpus: 85 gene pages, 4 bands each,
    4 distinct labels out of 6 curated phases. **Two** phases close before this
    dataset's first sampled stage (4 wpc) and are therefore eligible and
    unbanded — `embryonic_heart_tube_morphogenesis` ending 3.29 and
    `heart_looping` ending 3.71. That is why the expected set is computed from
    the vocabulary and the dataset rather than pinned to a literal that a
    second dataset would falsify.
    """
    corpus = load_curation(REPO)[0]
    phases = corpus.cardiac_phases
    assert phases is not None, "the committed corpus has no phase vocabulary; fixture is broken"
    curated = {phase.label for phase in phases.phases}

    sampled = {
        stage.wpc
        for dataset in corpus.datasets
        for stage in dataset.stages
        if stage.wpc is not None
    }
    assert sampled, "no dataset declares a prenatal stage; the expected set would be empty"
    expected = {phase.label for wpc in sampled for phase in phases.phases_for(wpc)}
    assert expected, "no sampled stage falls in any curated phase; this test would assert nothing"

    banded: set[str] = set()
    pages_with_a_band = 0
    for page in sorted((site / "genes").glob("HGNC_*.html")):
        labels = _BAND.findall(page.read_text(encoding="utf-8"))
        if labels:
            pages_with_a_band += 1
        banded.update(labels)

    assert pages_with_a_band, "no page carries a phase band; the vocabulary reached no page"
    assert banded <= curated, f"a band names a phase no curator wrote: {sorted(banded - curated)}"
    assert banded == expected, (
        f"the bands drawn are not the phases the vocabulary places at a sampled stage; "
        f"drawn but not expected: {sorted(banded - expected)}; "
        f"expected but not drawn: {sorted(expected - banded)}"
    )


def _tau_denominator(method: str, n_organs: int) -> float:
    """The divisor the published `method` string *names*, read out of the words.

    Deliberately a parser and not a constant. D39(a)'s whole promise is that a
    consumer who reads `method` and never opens this repository reproduces the
    number; the only way to test that promise is to follow the sentence rather
    than the source. A constant `n - 1` here would agree with `specificity`
    for the same reason the code agrees with itself, and would have passed
    unchanged while the string said "mean over organs" — which is what
    happened, through four review passes.

    Raises rather than defaulting on an unrecognised wording: a rewrite this
    function cannot follow is a rewrite a consumer cannot follow either, and
    silently picking a denominator would restore exactly the hole this closes.
    """
    if "divided by (n - 1)" in method:
        return float(n_organs - 1)
    if "mean over organs" in method:
        return float(n_organs)
    raise AssertionError(
        f"the published tau method names no denominator this test can follow: {method!r}"
    )


def test_every_published_tau_is_reproducible_from_its_own_method_string(site: Path) -> None:
    """D39(a) made executable: `method` + `medians` must reproduce `tau`.

    `specificity` publishes the formula beside the number precisely so a
    consumer can re-derive it without reading `build/profiles.py`. Nothing
    checked that the sentence and the arithmetic agreed, and they did not:
    the string said "mean over organs of (1 - x_i/x_max)" — a divisor of *n* —
    while the code divides by *n − 1*, the standard Yanai τ. Measured
    2026-08-20 over the committed corpus: all 1817 published specificity
    blocks match `sum/(n − 1)` and 4 also match `sum/n`, those four being
    coincidences of ties rather than agreement.

    So this asserts equality *and* that the two candidate divisors are
    distinguishable on this corpus — per section 4.30, two figures that are
    equal in the fixture are one figure to the test. Without the second
    assertion a corpus of two-organ stages (where n − 1 == 1 and the
    difference is a factor of 2 on every block, but a degenerate one-organ
    corpus would collapse) could let a wrong string pass.

    `scale` is asserted alongside because `method` names the transform and the
    two are separately published; a consumer following one and reading the
    other must not get two different answers.
    """
    blocks = 0
    also_matching_plain_mean = 0
    for bundle in sorted((site / "genes").glob("HGNC_*.json")):
        payload = json.loads(bundle.read_text(encoding="utf-8"))
        for dataset in payload["expression_profile"]["datasets"]:
            for stage in dataset["stages"]:
                specificity = stage["specificity"]
                if specificity is None:
                    continue
                blocks += 1
                medians = specificity["medians"]
                method = specificity["method"]

                assert specificity["scale"] == "log2(x+1)"
                assert "x = log2(median+1)" in method
                assert "clamped to 0" in method

                logged = {
                    organ: math.log2(max(value, 0.0) + 1.0) for organ, value in medians.items()
                }
                x_max = max(logged.values())
                total = sum(1.0 - value / x_max for value in logged.values())

                derived = total / _tau_denominator(method, len(medians))
                assert derived == pytest.approx(specificity["tau"], abs=1e-12), (
                    f"{bundle.name}: tau {specificity['tau']} is not what the published "
                    f"method reproduces ({derived}) from the published medians"
                )
                if total / len(medians) == pytest.approx(specificity["tau"], abs=1e-12):
                    also_matching_plain_mean += 1

    assert blocks, "no specificity block was published; this test measured nothing"
    assert also_matching_plain_mean < blocks, (
        "every published tau matches sum/n as well as sum/(n-1); this corpus "
        "cannot distinguish the two divisors, so the assertion above is vacuous"
    )


def test_every_trajectory_mark_sits_where_the_payload_puts_it(site: Path) -> None:
    """Chart geometry, recomputed from the bundle and compared to the bytes.

    A mutation review measured ~30 survivors in this area: marker radius set
    to zero, `cx` and `cy` transposed, the stage pitch divided by `count`
    rather than `count - 1`, band width doubled, forest rows stacked on one
    line. Every one changes where a reader sees a measurement and none changed
    a test, because the suite asserted which classes appear and never where
    anything is.

    One test rather than thirty: it walks every published trajectory, derives
    each placed stage's x from its index and each below-floor stage's tick
    from the same pitch, and asserts the drawn marks are exactly that set.
    A transposition, a radius change, a pitch change or a dropped tick all
    fail it, and they fail it on the real corpus rather than on a fixture.
    """
    left, right, bottom = 46.0, 550.0, 124.0
    checked = 0

    for page in sorted((site / "genes").glob("HGNC_*.html")):
        bundle = json.loads((page.with_suffix(".json")).read_text())
        datasets = (bundle.get("expression_profile") or {}).get("datasets") or []
        if not datasets:
            continue
        markup = page.read_text()
        section = markup[markup.find("Developmental expression") :]
        figures = re.findall(r'<svg class="chart" viewBox="0 0 560 202".*?</svg>', section, re.S)
        if not figures:
            continue

        stages = datasets[0]["stages"]
        count = len(stages)
        pitch = (right - left) / (count - 1)
        # Every cardiac organ this dataset declares, in the order the page draws
        # them; one trajectory each, so figure i belongs to organ i.
        for figure in figures:
            title = re.search(r"<title>Median abundance in ([^:,]+)", figure)
            assert title is not None, "a trajectory names the organ it draws"
            organ = title.group(1).strip()

            placed_x, tick_x = [], []
            for index, stage in enumerate(stages):
                for cell in stage["tissues"]:
                    if cell["tissue"] != organ:
                        continue
                    x = round(left + pitch * index, 1)
                    if cell["placement"] is not None and (cell["median_abundance"] or 0) > 0:
                        placed_x.append(x)
                    else:
                        tick_x.append(x)

            drawn_ticks = [
                float(m) for m in re.findall(r'<circle class="chart-absent" cx="([^"]+)"', figure)
            ]
            assert drawn_ticks == pytest.approx(tick_x, abs=0.05), (
                f"{page.stem} {organ}: below-floor ticks are not at the stage pitch"
            )

            drawn_points = re.findall(
                r'<circle class="chart-point" cx="([^"]+)" cy="([^"]+)" r="([^"]+)"', figure
            )
            assert [float(x) for x, _, _ in drawn_points] == pytest.approx(placed_x, abs=0.05), (
                f"{page.stem} {organ}: placed marks are not at the stage pitch "
                f"(a cx/cy transposition or a wrong divisor fails here)"
            )
            for _, cy, radius in drawn_points:
                assert 0.0 < float(cy) <= bottom, f"{page.stem} {organ}: a mark is off the plot"
                assert float(radius) > 0.0, f"{page.stem} {organ}: an invisible mark"
            checked += 1

    assert checked >= 80, f"only {checked} trajectories checked; the sweep found too few"


def test_datasets_json_orders_its_stages_by_order_not_by_declaration(site: Path) -> None:
    """The published stage array holds chronology by construction.

    It agreed with `order` only because the curator typed it that way.
    Measured before this guard: moving the `elderly` block to the top of
    `curation/datasets/E-MTAB-6814.yaml` while leaving `order: 21` alone gave
    `validate` the exact documented baseline -- 0 errors, 4 warnings -- and
    published `datasets.json` with `elderly` first, with gene bundles and
    pages unchanged. PRF011 and PRF014 constrain the *values* of `order`,
    never the sequence a YAML file lists them in.

    `bundles._concordance_for`'s lesson (CLAUDE.md 4.28): the gene bundles
    already sorted their own stage array, and this second consumer of the same
    input did not.
    """
    payload = _read(site, "datasets.json")
    checked = 0
    for dataset in payload["datasets"]:
        stages = dataset.get("stages") or []
        if not stages:
            continue
        orders = [stage["order"] for stage in stages]
        assert orders == sorted(orders), (
            f"{dataset['id']}: stages are not in `order` sequence, so the array "
            f"carries declaration order rather than chronology"
        )
        assert orders == list(range(1, len(orders) + 1)), (
            f"{dataset['id']}: `order` is not a contiguous 1..N sequence"
        )
        checked += 1
    assert checked >= 1, "no dataset carried stages; this guard would pass vacuously"


def test_every_gene_page_links_gnomad_and_clinvar_from_its_registry_ids(site: Path) -> None:
    """D52, the page half: what the unit suites cannot see. They pass `identities`
    to `build_gene_pages` themselves or not at all, so a runner that stopped
    threading it would leave them green while every real page lost its link row.
    The count is derived, never pinned."""
    with (REPO / "mirrors" / "genes.tsv").open(encoding="utf-8") as handle:
        registry = {row["hgnc_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    entries = _read(site, "genes/index.json")["genes"]
    assert len(entries) > 1, "the build published no genes; the fixture is broken"
    for entry in entries:
        row = registry[entry["gene"]]
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        start = page.index('<p class="xrefs">')
        xrefs = page[start : page.index("</p>", start)]
        gnomad = f"https://gnomad.broadinstitute.org/gene/{row['ensembl_gene']}?dataset=gnomad_r4"
        clinvar = f"https://www.ncbi.nlm.nih.gov/clinvar/?term={row['ncbi_gene']}%5Bgeneid%5D"
        assert f'href="{gnomad}"' in xrefs, entry["gene"]
        assert f'href="{clinvar}"' in xrefs, entry["gene"]


def _s8_sentence(burden: list[dict[str, Any]], study: str) -> str:
    """What a card adds after `study`'s summary caveat, worded from the bundle (R2).

    One row per gene carries `n_chd_case_carriers` for PMID:34324492, the only
    study with one, so more than one is a fixture this cannot read.
    """
    rows = [
        row for row in burden if row["study"] == study and row["n_chd_case_carriers"] is not None
    ]
    assert len(rows) <= 1, study
    if not rows:
        return ""
    chd, carriers = rows[0]["n_chd_case_carriers"], rows[0]["n_case_carriers"]
    if chd == carriers:
        return (
            " For this gene, the study's S8 Table counts every carrier among its CHD cases, "
            "so none is an aneurysm proband."
        )
    return (
        f" For this gene, the study's S8 Table counts {chd:,} of the {carriers:,} carriers "
        "among its CHD cases."
    )


def test_every_card_lists_rare_variants_exactly_where_its_bundle_counts_a_corrected_dataset(
    site: Path,
) -> None:
    """D55/D56 against every real page: the card and the bundle cannot disagree.

    A gene's card lists a surviving rare-variant result iff its bundle's
    `independent_datasets.corrected` is above zero -- both read
    `concordance.survives_own_correction`, and this is what fails if a runner
    or a renderer stops agreeing with it on real data rather than on a
    fixture. Lines are counted inside the card's rare-variant block only, not
    the whole rail. A card listing none says why: the none-survives sentence
    where the bundle has burden rows, the no-study sentence where it has none.

    **No card may say a result failed a correction its study never made.**
    PMID:42230622 publishes no corrected p on any row, and SMAD4's card read
    "none survives its own study's correction" of it, while its block judged
    two of its rows against a scan that they clear. So every study that tested
    the gene and published no correction for it is named on its own line, and
    exactly those; the old sentence appears on no card. A study by an atlas
    author carries the declaration on every card line that names it.

    **A study's summary caveat follows its card lines**: once per study with a
    line on the card, in study-id order, with exactly the text
    `publications.json` publishes for it, and on no card without a line from
    that study -- so PMID:34324492's thoracic aortic aneurysm probands reach
    every card that shows its counts. The gene index, whose strips summarise
    every study's result for every gene, prints each caveat exactly once.
    After it, where the study's row carries `n_chd_case_carriers`, what S8
    Table says of the CHD cases among its carriers (owner decision R2,
    2026-10-27), worded from the bundle's two counts; GATA4's card is pinned
    as a literal, and the sentence must reach at least one card.

    And every card states once what the atlas itself curated, nothing in the
    main column repeats it, a card names who disputes the gene exactly where
    the bundle's `has_conflicting_evidence` is true, and every card ends on one
    link to `methods.html` and one to the gene's own JSON.

    Counts derived, never pinned. Both sides of each "iff" must be exercised,
    or the test is measuring one population twice (CLAUDE.md section 4.30).

    **`counts.burden_surviving_genes` is the number of cards listing a
    result** (owner decision 2026-10-26, R9): the census counts by the
    card's own predicate, and this is where the two would disagree on real
    data. 32 on the build it was written against, of 118 genes with a row.
    """
    entries = _read(site, "genes/index.json")["genes"]
    assert len(entries) > 1, "the build published no genes; the fixture is broken"
    publications = _read(site, "publications.json")["publications"]
    own_lab = {publication["id"]: publication["own_lab"] for publication in publications}
    caveats = {publication["id"]: publication["summary_caveat"] for publication in publications}
    listed = flagged = uncorrected_cards = noted_cards = split_cards = 0
    on_strips: set[str] = set()
    for entry in entries:
        bundle = _read(site, entry["bundle"])
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
        card = ""
        if '<div class="card-variants">' in rail:
            card = rail[rail.index('<div class="card-variants">') :]
            card = card[: card.index("</div>")]
        lines = card.count("<li>")
        corrected = bundle["independent_datasets"]["corrected"]
        assert (lines > 0) is (corrected > 0), (entry["gene"], lines, corrected)
        tested = {row["study"] for row in bundle["burden"]}
        published = {row["study"] for row in bundle["burden"] if row["pvalue_adjusted"] is not None}
        none_survives = "none survives a correction its study published." in card
        if lines:
            listed += 1
        elif bundle["burden"]:
            # "None survives" only where some study published a correction.
            assert none_survives is bool(published), entry["gene"]
            assert "<p>Tested in " in card, entry["gene"]
        else:
            assert "No rare-variant study in this atlas reports a result" in card, entry["gene"]
        assert "none survives its own study" not in page, entry["gene"]

        named = re.findall(r'<p class="card-uncorrected"><a href="#study-(\d+)">', card)
        assert sorted(f"PMID:{digits}" for digits in named) == sorted(tested - published), entry[
            "gene"
        ]
        uncorrected_cards += bool(named)
        for item in re.findall(r"<li>(.*?)</li>|<p class=\"card-uncorrected\">(.*?)</p>", card):
            text = item[0] or item[1]
            digits = re.search(r'href="#study-(\d+)"', text)
            assert digits, (entry["gene"], text)
            marked = "· study by an atlas author" in text
            assert marked is own_lab[f"PMID:{digits.group(1)}"], (entry["gene"], text)

        lined = {
            f"PMID:{digits}" for digits in re.findall(r'<li>[^<]*<a href="#study-(\d+)"', card)
        }
        notes = re.findall(
            r'<p class="card-note"><a href="#study-(\d+)">[^<]*</a>: (.*?)</p>', card
        )

        assert notes == [
            (
                study.removeprefix("PMID:"),
                html.escape(caveats[study] + _s8_sentence(bundle["burden"], study)),
            )
            for study in sorted(lined)
            if caveats[study]
        ], entry["gene"]
        assert card.count('class="card-note"') == len(notes), entry["gene"]
        noted_cards += bool(notes)
        split_cards += "S8 Table counts" in card
        if entry["gene"] == "HGNC:4173":  # GATA4: 29 carriers, all CHD cases in S8
            assert (
                card.count(
                    "For this gene, the study&#x27;s S8 Table counts every carrier among its CHD "
                    "cases, so none is an aneurysm proband.</p>"
                )
                == 1
            )
        on_strips.update(
            study
            for family in bundle["independent_datasets"]["families"]
            for study in family["studies"]
        )

        assert rail.count('<p class="card-curation">') == 1, entry["gene"]
        assert rail.count('<a href="../methods.html">How to read this page</a>') == 1, entry["gene"]
        assert rail.endswith(
            '<p><a href="../methods.html">How to read this page</a> &middot; '
            f'<a href="../{entry["bundle"]}">this gene as JSON</a></p>'
        ), entry["gene"]
        main = page[page.index("</aside>") :]
        assert "curated by the atlas" not in main.lower(), entry["gene"]
        disputed = rail.count('<p class="card-dispute">')
        assert disputed == int(bundle["has_conflicting_evidence"]), entry["gene"]
        flagged += disputed
    assert 0 < listed < len(entries), f"{listed} of {len(entries)} cards list a result"
    # R9: the census the front page and `manifest.json` state is these cards.
    surviving = _read(site, "manifest.json")["counts"]["burden_surviving_genes"]
    assert surviving == listed, (surviving, listed)
    assert 0 < flagged < len(entries), f"{flagged} of {len(entries)} cards flag a dispute"
    assert 0 < uncorrected_cards < len(entries), (
        f"{uncorrected_cards} of {len(entries)} cards name a study with no correction"
    )
    assert 0 < noted_cards < len(entries), f"{noted_cards} of {len(entries)} cards carry a note"
    assert 0 < split_cards <= noted_cards, f"{split_cards} cards carry the S8 Table sentence"

    browse = (site / "genes" / "index.html").read_text(encoding="utf-8")
    assert any(caveats.get(study) for study in on_strips), "no strip names a noted study"
    for study, caveat in caveats.items():
        if caveat:
            assert browse.count(html.escape(caveat)) == int(study in on_strips), study


def test_every_contested_gene_is_chipped_on_the_gene_index_and_no_other_is(site: Path) -> None:
    """The gene index never shows a contested gene as settled (2026-10-26).

    `docs/data-api.md`'s one consumer obligation is that a gene flagged
    `has_conflicting_evidence` is never displayed as settled, and the front page
    says the same. The site's own gene index broke it: its confidence column
    rendered BRAF, KRAS, MAP2K1, PTPN11, RAF1 and SOS1 as a bare `definitive`,
    while the `genes/index.json` behind the page flagged every one.

    The flagged set is read from that JSON, never written here, so this keeps
    testing the right genes as the population moves; it must be non-empty, or
    the positive half measures nothing. Every published gene must have exactly
    one row, so a gene the page dropped cannot pass as "unflagged". The cell is
    compared whole -- the bare grade (or the em dash a gene no panel graded
    shows) and the card's own chip -- and `data-confidence` must stay the bare
    grade, which is what the confidence facet matches on.
    """
    entries = _read(site, "genes/index.json")["genes"]
    flagged = {entry["gene"] for entry in entries if entry["has_conflicting_evidence"]}
    assert flagged, "no published gene is flagged has_conflicting_evidence; this measures nothing"
    browse = (site / "genes" / "index.html").read_text(encoding="utf-8")
    chip = '<span class="chip chip-warn">conflicting evidence</span>'
    for entry in entries:
        page = entry["gene"].replace(":", "_")
        rows = re.findall(
            rf"<tr[^>]*>(?:(?!</tr>).)*{page}\.html(?:(?!</tr>).)*</tr>", browse, re.S
        )
        assert len(rows) == 1, f"{entry['gene']} has {len(rows)} rows on the gene index"
        confidence = entry["headline_confidence"] or ""
        assert f'data-confidence="{confidence}"' in rows[0], entry["gene"]
        cell = re.findall(r"<td>(.*?)</td>", rows[0], re.S)[2]
        bare = confidence or "—"
        assert cell == (f"{bare} {chip}" if entry["gene"] in flagged else bare), entry["gene"]


def test_exactly_the_genes_graded_for_a_clingen_scoped_syndrome_say_so_on_their_card(
    site: Path,
) -> None:
    """R5 (owner decision 2026-10-26) on the real corpus.

    Ten genes are graded at their headline for a syndrome a ClinGen expert
    panel scoped, and their cards said nothing of it while the 21 Orphanet-
    scoped ones did. `scope.json` now publishes the curated `syndrome` mark,
    on exactly the nine terms those ten are graded for; ZIC3 and NODAL, graded
    for "congenital heart disease with heterotaxy syndrome", are CHD itself and
    stay unmarked. Two ways, both required to agree: the genes whose card
    carries the line are the literal ten, and they are the genes whose every
    headline disease `scope.json` marks and attributes to a ClinGen class --
    derived from the payloads, so a runner that stopped threading the registry
    to the pages fails here while every unit test passes. TBX5's line is pinned
    as a literal.

    So is EVC2's, the one line whose panel did not grade its gene (review,
    2026-10-27): it names the gene that panel graded for the term instead.
    Every line is held to that rule from the bundles' own `gcep`, and the
    lines it applies to must be exactly EVC2's, so a page that dropped the
    distinction, or drew it on another gene, fails here.
    """
    scope = {term["id"]: term for term in _read(site, "scope.json")["terms"]}
    assert sorted(mondo for mondo, term in scope.items() if term["syndrome"] is True) == [
        "MONDO:0007318",  # Alagille syndrome
        "MONDO:0007732",  # Holt-Oram syndrome
        "MONDO:0007843",  # Kabuki syndrome 1
        "MONDO:0008965",  # CHARGE syndrome
        "MONDO:0009162",  # Ellis-van Creveld syndrome
        "MONDO:0009341",  # Mowat-Wilson syndrome
        "MONDO:0010465",  # Kabuki syndrome 2
        "MONDO:0010711",  # TARP syndrome
        "MONDO:0800458",  # NR2F2-related multiple congenital anomalies/dysmorphic syndrome
    ]
    assert all(type(term["syndrome"]) is bool for term in scope.values())

    lined: set[str] = set()
    lined_by_another: set[str] = set()
    expected: set[str] = set()
    for entry in _read(site, "genes/index.json")["genes"]:
        bundle = _read(site, entry["bundle"])
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
        lines = re.findall(r'<p class="card-scope">([^<]*in scope on [^<]*)</p>', rail)
        if lines:
            lined.add(entry["symbol"])
        if entry["symbol"] == "TBX5":
            assert lines == [
                "Graded for Holt-Oram syndrome as a whole syndrome; in scope on its curation by "
                "ClinGen&#x27;s Syndromic Disorders Gene Curation Expert Panel."
            ]
        if entry["symbol"] == "EVC2":
            assert lines == [
                "Graded for Ellis-van Creveld syndrome as a whole syndrome; in scope on "
                "ClinGen&#x27;s Kidney Cystic and Ciliopathy Disorders Gene Curation Expert "
                "Panel&#x27;s curation of it for EVC."
            ]
        headline = bundle["headline_confidence"]
        graded = {
            record["disease"]
            for record in bundle["validity"]["records"]
            if record["source"] == "clingen" and record["classification"] == headline
        }
        # "its curation by <panel>" only where that panel graded this gene for
        # the term (review, 2026-10-27): the term's `attributed_to` is the panel
        # that scopes the term, and EVC2's grader is another. Read off the
        # bundle's own ClinGen records and `scope.json`, on every lined gene.
        for line in lines:
            [disease] = graded
            own = {
                record["gcep"]
                for record in bundle["validity"]["records"]
                if record["source"] == "clingen" and record["disease"] == disease
            }
            panel = scope[disease]["attributed_to"]
            if panel in own:
                assert line.endswith(f"in scope on its curation by ClinGen&#x27;s {panel}.")
            else:
                assert f"in scope on ClinGen&#x27;s {panel}&#x27;s curation of it for " in line
                lined_by_another.add(entry["symbol"])
        if (
            headline is not None
            and graded
            and all(
                scope[disease]["syndrome"]
                and scope[disease]["admitted_by"] in {"clingen_chd_panel", "clingen_expert_panel"}
                for disease in graded
            )
        ):
            expected.add(entry["symbol"])
    assert lined == expected
    assert lined_by_another == {"EVC2"}
    assert sorted(lined) == [
        "CHD7",
        "EVC",
        "EVC2",
        "KDM6A",
        "KMT2D",
        "NOTCH2",
        "NR2F2",
        "RBM10",
        "TBX5",
        "ZEB2",
    ]


def test_every_orphanet_warrant_scope_json_publishes_reaches_the_pages_that_list_its_disease(
    site: Path,
) -> None:
    """The Orphanet warrants reach the gene pages, not only `scope.json` (D55, issue #8).

    `build_site` derives the warrants once and hands them to `build_scope` and
    to `build_gene_pages`. The unit suites hand `build_gene_pages` warrants of
    their own or none, so a runner that stopped threading them left them green:
    measured 2026-10-26, `scope_warrants={}` in `build_site` passed all 1,297
    tests while the pages carrying "is in scope because Orphanet records" went
    from 27 to 0 and the card's scope line from 21 to 0, and `scope.json` still
    published all 36 warrants. PTPN11's Noonan grade then sits on a CHD atlas
    with nothing saying why Noonan counts, and reads as a grade for heart
    disease.

    Read from `scope.json` and each bundle, never pinned:

    * the validity section carries one note per distinct Orphanet-scoped
      disease among the bundle's records, in disease-id order, and nothing
      else on the page does -- a note for a disease the gene has no record
      for is a warrant stated about nothing the page lists. Each names its
      disease by the registry's label (the labels are threaded separately, and
      without them the note names a bare MONDO id) and states the code, the
      syndrome a subtype inherits from, and every phenotype at the band
      `scope.json` publishes, so the page and the registry state one warrant;
    * the card's scope line appears exactly where every disease the headline
      grade is for is Orphanet-scoped, names those diseases, and names the
      syndrome a subtype's frequency was recorded for.

    Both sides of each "exactly where" must be exercised (CLAUDE.md section
    4.30), and the subtype branch too, or it is measuring nothing.
    """
    orphanet = {term["id"]: term for term in _read(site, "scope.json")["terms"] if term["orphanet"]}
    assert orphanet, "scope.json publishes no Orphanet warrant; this measures nothing"
    entries = _read(site, "genes/index.json")["genes"]
    noted = carded = subtype_notes = subtype_cards = 0
    for entry in entries:
        bundle = _read(site, entry["bundle"])
        page = (site / entry["bundle"].replace(".json", ".html")).read_text(encoding="utf-8")
        start = page.index("<h2>Mirrored gene&ndash;disease validity</h2>")
        section = page[start : page.index("<h2", start + 1)]
        records = bundle["validity"]["records"]

        diseases = sorted({record["disease"] for record in records} & orphanet.keys())
        if diseases:
            # The scope table has exactly the disease rows, above the validity
            # groups; slice to it so a six-cell validity row cannot be confused
            # with a three-cell scope row by the regex below.
            scope_start = section.index("<h3>Why these diseases are in scope</h3>")
            scope_end = section.index('id="validity-table"')
            rows = re.findall(
                r"<tr><td>(.*?)</td><td>(.*?)</td><td>(.*?)</td></tr>",
                section[scope_start:scope_end],
            )
        else:
            assert "<h3>Why these diseases are in scope</h3>" not in section, entry["gene"]
            rows = []
        assert len(rows) == len(diseases), (entry["gene"], diseases)
        for disease, (label_cell, why_cell, malformations_cell) in zip(diseases, rows, strict=True):
            term = orphanet[disease]
            label = html.escape(term["label"])
            assert label_cell == label, (entry["gene"], label_cell)
            warrant = term["orphanet"]
            note = why_cell + malformations_cell
            assert warrant["code"] in note, (entry["gene"], disease)
            if warrant["via"] == "subtype":
                subtype_notes += 1
                assert f"Subtype of {html.escape(warrant['name'])}; " in why_cell, entry["gene"]
                assert warrant["omim"] in why_cell, entry["gene"]
                # The "its cases" in the why cell is the parent's, so the row
                # carries the code only in the subtype mapping, and the
                # malformations cell for a subtype carries no "Orphanet, ...".
                assert f"(Orphanet, {warrant['code']})" not in malformations_cell, entry["gene"]
            else:
                assert f"(Orphanet, {warrant['code']})" in malformations_cell, entry["gene"]
            for phenotype in warrant["phenotypes"]:
                stated = (
                    f"{html.escape(phenotype['label'])} at "
                    f"<em>{html.escape(phenotype['frequency'])}</em>"
                )
                assert stated in malformations_cell, (entry["gene"], disease, phenotype["hpo"])
        noted += bool(rows)

        headline = bundle["headline_confidence"]
        graded = [
            record
            for record in records
            if headline is not None
            and record["source"] == "clingen"
            and record["classification"] == headline
        ]
        scoped = bool(graded) and all(record["disease"] in orphanet for record in graded)
        rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
        every_line = re.findall(r'<p class="card-scope">(.*?)</p>', rail)
        assert page.count('class="card-scope"') == len(every_line), entry["gene"]
        # The Orphanet lines; a syndrome a ClinGen panel scoped has its own (R5).
        lines = [line for line in every_line if "in scope because Orphanet records" in line]
        assert len(lines) == int(scoped), (entry["gene"], lines)
        if scoped:
            carded += 1
            line = lines[0]
            for record in graded:
                assert html.escape(record["disease_label"]) in line, (entry["gene"], line)
            recorded_for = line[line.index("30% of cases") :]
            warrants = [orphanet[record["disease"]]["orphanet"] for record in graded]
            if all(warrant["via"] == "direct" for warrant in warrants):
                assert recorded_for == "30% of cases.", (entry["gene"], line)
            for warrant in warrants:
                if warrant["via"] == "subtype":
                    subtype_cards += 1
                    assert html.escape(warrant["name"]) in recorded_for, (entry["gene"], line)
    assert 0 < noted < len(entries), f"{noted} of {len(entries)} pages carry an Orphanet note"
    assert 0 < carded < len(entries), f"{carded} of {len(entries)} cards carry a scope line"
    assert subtype_notes, "no page notes a subtype's warrant; that branch measures nothing"
    assert subtype_cards, "no card names a subtype's syndrome; that branch measures nothing"


def test_every_fragment_link_on_a_gene_page_lands_on_one_element_of_that_page(
    site: Path,
) -> None:
    """The card's links into its own page (D55), which the link sweep skips.

    `test_every_link_on_every_page_resolves_to_a_file_the_build_wrote`
    deliberately ignores pure fragments, so a card line pointing at
    `#study-34324492` on a page whose study heading carried no such id would
    pass it -- a link that scrolls nowhere, beside the counts it was added to
    qualify. Every `href="#..."` on a gene page must name exactly one `id` on
    that page, and no `id` may repeat on a page at all. Counted, not pinned;
    the sweep must find at least one fragment link or it measures nothing.
    """
    fragments = 0
    for page_path in sorted((site / "genes").glob("HGNC_*.html")):
        page = page_path.read_text(encoding="utf-8")
        ids = re.findall(r' id="([^"]+)"', page)
        repeated = sorted({value for value in ids if ids.count(value) > 1})
        assert not repeated, f"{page_path.name} repeats ids {repeated}"
        for target in re.findall(r'href="#([^"]+)"', page):
            fragments += 1
            assert target in ids, f"{page_path.name} links #{target}, which no element carries"
    assert fragments, "no gene page carries a fragment link; the sweep would pass vacuously"


def test_the_front_page_declares_exactly_the_study_blocks_the_gene_pages_declare(
    site: Path,
) -> None:
    """The landing's own-lab sentence, against the blocks it vouches for (R7).

    "Of the 3 rare-variant studies whose results appear here, 2 are by an
    author of this atlas, and a gene page says so wherever it shows their
    results." Both numbers are read here off the gene pages -- the distinct
    study blocks they render, and which of those `publications.json` marks
    `own_lab` -- and every block of such a study must carry the declaration,
    while no other study's block may. Measured on a real build the day it was
    written: 190 own-lab blocks, all declared; 36 others, none.

    The declaration said "republished as published, no assessment added"
    until review (2026-10-27), on every one of those blocks, while
    `methods.html#rare-variants` states the rule the atlas applies to every
    study's rows. The phrase is asserted absent from every page built.
    """
    own = {
        entry["id"]
        for entry in json.loads((site / "publications.json").read_text(encoding="utf-8"))[
            "publications"
        ]
        if entry["own_lab"]
    }
    declaration = "Declaration: this study is by an author of this atlas"
    rendered: set[str] = set()
    for path in sorted((site / "genes").glob("HGNC_*.html")):
        page = path.read_text(encoding="utf-8")
        start = page.find('<h2 id="rare-variants">')
        if start == -1:
            continue
        stop = page.find("<h2", start + 1)
        section = page[start : stop if stop != -1 else len(page)]
        heads = list(re.finditer(r'<h3 id="study-(\d+)">', section))
        ends = [head.start() for head in heads[1:]] + [len(section)]
        for head, end in zip(heads, ends, strict=True):
            block = section[head.start() : end]
            study = f"PMID:{head.group(1)}"
            rendered.add(study)
            assert (declaration in block) == (study in own), (path.name, study)
    assert rendered, "no gene page renders a study block; the check would pass vacuously"
    mine = len(rendered & own)
    assert mine, "no rendered study is by an atlas author; the sentence's main form is unmeasured"
    landing = " ".join(
        html.unescape(re.sub(r"<[^>]+>", "", (site / "index.html").read_text())).split()
    )
    assert (
        f"Of the {len(rendered)} rare-variant studies whose results appear here, "
        f"{mine} {'is' if mine == 1 else 'are'} by an author of this atlas"
    ) in landing
    built = [*sorted((site / "genes").glob("*.html")), site / "index.html", site / "methods.html"]
    for path in built:
        assert "no assessment added" not in path.read_text(encoding="utf-8"), path.name


# Per source `sources.json` records as `permitted_with_attribution`, a phrase
# the front page must carry: its name, and its licence where that is CC BY.
_ATTRIBUTED_ON_THE_FRONT_PAGE: dict[str, str] = {
    "hpo": "Human Phenotype Ontology",
    "orphanet_product1": "mirrored under CC BY 4.0 — Orphadata Science",
    "orphanet_product4": "mirrored under CC BY 4.0 — Orphadata Science",
    "mondo_sssom": "Mondo's disease mappings are mirrored under CC BY 4.0",
    "expression_atlas_e_mtab_6814": "EMBL-EBI Expression Atlas's (E-MTAB-6814), "
    "under EMBL-EBI's Terms of Use",
    "audain2026_sd3": "Audain et al. 2026 and Sierant et al. 2025, both CC BY 4.0",
    "sierant2025_ds6": "Audain et al. 2026 and Sierant et al. 2025, both CC BY 4.0",
}
# A source the sentence names though `sources.json` requires no attribution for
# it, beside the two supplements that do.
_ALSO_ON_THE_FRONT_PAGE: dict[str, str] = {"audain2021_s6": "and of Audain et al. 2021, CC0."}
# The words the sentence states each `sources.json` licence in. HPO's licence is
# a URL the sentence does not state, so it is the one source held to its name
# alone.
_LICENCE_WORDS: dict[str, str] = {
    "CC-BY-4.0": "CC BY 4.0",
    "CC0-1.0": "CC0",
    "EMBL-EBI Terms of Use": "EMBL-EBI's Terms of Use",
}
_LICENCE_NOT_STATED: frozenset[str] = frozenset({"hpo"})


def test_the_front_page_names_every_source_that_requires_attribution(site: Path) -> None:
    """The attribution sentence is hand-written, so it is held to `sources.json` (R10).

    Owner decision 2026-10-26: "Where this data comes from" named ClinGen,
    GenCC, HPO and Orphanet, and not Mondo, the Expression Atlas or the CC BY
    rare-variant supplements, all of which `sources.json` records as
    `permitted_with_attribution`. One hand-written sentence then named them.
    The set of attribution-required sources must equal the keys pinned above,
    so a source added to `mirrors/sources.yaml` with that term fails here until
    the page names it and a phrase is pinned; each phrase must be on the page,
    with the licence where it is CC BY. Measured on a real build the day it was
    written: 7 such sources. Since the doorway page (2026-10-27) the section is
    "Sources", one paragraph of 123 words, and the sentence is three.

    **Every licence the sentence states is held to `sources.json`** (review,
    2026-10-27), not only CC BY: "Audain et al. 2021, CC0" and "EMBL-EBI's
    Terms of Use" were pinned by a unit literal alone, so a licence changed in
    `mirrors/sources.yaml` left the hand-written sentence stating the old one
    on a green build. Each phrase must carry the words for the licence its
    source records, and a licence the sentence has no words for fails here.
    """
    sources = json.loads((site / "sources.json").read_text(encoding="utf-8"))["sources"]
    required = {
        source["id"]
        for source in sources
        if source["redistribution"] == "permitted_with_attribution"
    }
    assert required == set(_ATTRIBUTED_ON_THE_FRONT_PAGE), (
        "sources.json's attribution-required sources changed; name each on index.html "
        f"and pin a phrase here: {sorted(required ^ set(_ATTRIBUTED_ON_THE_FRONT_PAGE))}"
    )
    licences = {source["id"]: source["licence"] for source in sources}
    page = (site / "index.html").read_text(encoding="utf-8")
    where = page[page.index("<h2>Sources</h2>") :]
    where = " ".join(html.unescape(re.sub(r"<[^>]+>", "", where[: where.index("<h2", 4)])).split())
    named = {**_ATTRIBUTED_ON_THE_FRONT_PAGE, **_ALSO_ON_THE_FRONT_PAGE}
    for source, phrase in sorted(named.items()):
        assert phrase in where, f"{source}: {phrase!r} is not on the front page"
        if source in _LICENCE_NOT_STATED:
            continue
        licence = licences[source]
        assert licence in _LICENCE_WORDS, f"{source}: the front page states no {licence!r}"
        assert _LICENCE_WORDS[licence] in phrase, f"{source}: {licence} is not what it states"

    # R11 (owner decision 2026-10-26): the Orphadata citation's "Data version"
    # is the release `sources.json` publishes for both products, read from the
    # files' own headers -- one release, one day, so stated once.
    releases = {
        (source["file_version"], source["file_date"].split(" ")[0])
        for source in sources
        if source["id"] in {"orphanet_product4", "orphanet_product1"}
    }
    assert releases == {("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23")}
    (version, day) = releases.pop()
    assert (
        f"http://sciences.orphadata.com/. Data version {version}, {day}. Developmental expression"
    ) in where


def test_no_page_says_a_sample_was_a_whole_organ(site: Path) -> None:
    """Owner decision 2026-10-26 (E2): "bulk tissue", never "whole".

    Every expression sentence said a sample was a whole organ -- "from whole,
    bulk tissue", "median abundance in whole heart", "below the detection floor
    in whole liver", methods.html's "A whole embryonic heart is ...", "diluted
    across the whole organ" -- while E-MTAB-6814's own curation record sets the
    Expression Atlas's organ tokens (forebrain, hindbrain) against the paper's
    (brain, cerebellum), so "whole" need not hold of every organ and stage. The
    dilution argument rests on bulk, mixed cell types, alone.

    Over the whole text of every gene page and `methods.html`, `<title>`s of
    the charts included: "whole" before "bulk", an organ, a tissue, a sample,
    "embryonic" or any organ token a published bundle names. Measured with
    this pattern on a build of 3b766ba: all 119 gene pages and methods.html
    matched; on a build of the change, none. The other meanings of "whole" --
    a whole scan, a syndrome as a whole, the whole gene set -- name none of
    these and are untouched.
    """
    tissues = {
        cell["tissue"]
        for path in sorted((site / "genes").glob("HGNC_*.json"))
        for dataset in json.loads(path.read_text(encoding="utf-8"))["expression_profile"][
            "datasets"
        ]
        for stage in dataset["stages"]
        for cell in stage["tissues"]
    }
    assert tissues, "no bundle names an organ; the check would pass vacuously"
    whole = re.compile(
        r"\bwhole[\s,-]+(?:bulk|organs?|tissues?|samples?|embryonic|"
        + "|".join(re.escape(tissue) for tissue in sorted(tissues, key=len, reverse=True))
        + r")\b",
        re.IGNORECASE,
    )
    pages = [*sorted((site / "genes").glob("HGNC_*.html")), site / "methods.html"]
    for path in pages:
        text = " ".join(html.unescape(path.read_text(encoding="utf-8")).split())
        found = sorted(set(whole.findall(text)))
        assert not found, f"{path.name} still says a sample was a whole {found}"


# --- The gene page's prose budget (D50, D57) ----------------------------------
#
# What a reader sees on load: text outside every `<details>`, or inside the
# `<summary>` of an outermost one. Text in `<svg>`, `<style>`, `<script>` and
# `<head>` is not prose. A block tag ends a run of text, so two paragraphs are
# never read as one sentence. Ported from the caveat inventory's parser
# (2026-10-25), with one difference measured below: `<br>` is not read as a
# space. Every `<br>` on a gene page sits in a table cell (measured on a build
# of 802d694), which neither figure counts.
_BLOCK_TAGS = frozenset(
    {
        "p", "li", "figcaption", "td", "th", "dt", "dd", "summary", "h1", "h2", "h3", "h4",
        "h5", "h6", "nav", "div", "ul", "ol", "table", "tr", "details", "figure", "aside",
        "main", "body", "dl", "thead", "tbody", "form", "select", "option",
    }
)  # fmt: skip
_SKIPPED = frozenset({"svg", "style", "script", "head"})
# Table and definition cells: data, not prose, and the one definition of a cell
# both figures share. These four are every cell and term tag a gene page uses
# (no `<caption>`), measured on a build of 802d694.
_CELLS = frozenset({"td", "th", "dt", "dd"})
# A word is a whitespace token carrying a letter, a digit, tau or infinity, so
# a bare dash or middle dot between two links is not one.
_WORD = re.compile(r"[A-Za-z0-9τ∞]")
_SENTENCE = re.compile(r"(?<=[.!?])\s+(?=[A-Z(\"'“])")
# (b)'s share: a sentence on at least this fraction of gene pages is boilerplate.
_BOILERPLATE_SHARE = 0.9

# Measured after the owner's answers of 2026-10-27 (the registry qualifier on
# the gene cards), with the headroom owner decision E4 keeps: (a) 964 x 1.03
# rounded up, (b) 21 + 1. See the test for the series and the rule. Lower each
# when a page sheds prose; raising one is a decision, and its commit says why.
_MEDIAN_VISIBLE_WORDS = 993
_BOILERPLATE_SENTENCES = 22


class _VisibleText(HTMLParser):
    """A gene page's text blocks as a reader sees them on load."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.skip = self.details = self.summary = self.cells = 0
        self.buffer: list[str] = []
        self.blocks: list[tuple[str, bool]] = []  # (text, inside a table/definition cell)

    def _flush(self) -> None:
        text = " ".join("".join(self.buffer).split())
        self.buffer = []
        if text and (self.details == 0 or self.summary):
            self.blocks.append((text, self.cells > 0))

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag in _SKIPPED:
            self.skip += 1
            return
        if tag in _BLOCK_TAGS:
            self._flush()
        if tag == "details":
            self.details += 1
        elif tag == "summary" and self.details == 1:
            self.summary += 1
        if tag in _CELLS:
            self.cells += 1

    def handle_endtag(self, tag: str) -> None:
        if tag in _SKIPPED:
            self.skip -= 1
            return
        if tag in _BLOCK_TAGS:
            self._flush()
        if tag == "details":
            self.details -= 1
        elif tag == "summary" and self.summary:
            self.summary -= 1
        if tag in _CELLS:
            self.cells -= 1

    def handle_data(self, data: str) -> None:
        if not self.skip:
            self.buffer.append(data)


def _visible_blocks(markup: str) -> list[tuple[str, bool]]:
    parser = _VisibleText()
    parser.feed(markup)
    parser.close()
    return parser.blocks


def _prose(markup: str) -> list[str]:
    """The visible blocks outside every cell: what both figures measure."""
    return [text for text, cell in _visible_blocks(markup) if not cell]


def _words(text: str) -> int:
    return sum(1 for token in text.split() if _WORD.search(token))


def _boilerplate(pages: list[list[str]]) -> list[str]:
    """(b): every sentence of six words or more on at least the share of `pages`
    the budget names, each once however many pages carry it."""
    carried: Counter[str] = Counter()
    for blocks in pages:
        carried.update(
            {
                sentence
                for text in blocks
                for sentence in _SENTENCE.split(text)
                if _words(sentence) >= 6
            }
        )
    return sorted(
        sentence for sentence, n in carried.items() if n >= _BOILERPLATE_SHARE * len(pages)
    )


def test_a_gene_page_carries_no_more_visible_prose_than_its_measured_budget(site: Path) -> None:
    """D50: site-wide explanation lives once, on `methods.html`; this keeps it there.

    Two figures over every `genes/HGNC_*.html` of the real build, each at most a
    constant set from a measurement, and both over **prose, not data**: the
    visible blocks outside every table and definition cell (`_CELLS`). (a) The
    median of the words visible on load per page; (b) the boilerplate, the
    visible sentences of six words or more that are identical on at least 90%
    of gene pages, each counted once however many pages carry it.

    **(a) counted cells until review (2026-10-26), and so measured data.** The
    burden, validity and stage tables a page shows are rows of mirrored data,
    and a curation-only change failed the suite: one more GenCC row for USP44
    (HGNC:20064), the median page at 957, moved the median to 958 -- and the
    deploy workflow runs this test. Re-run on a `git archive` copy of the
    branch with the row added (`SGC-999001`, PanelApp Australia, `Limited`
    for MONDO:0005453): (a) is 841 with and without it, this test passes on
    it, and USP44's own page falls from 862 to 844, because the sentence
    saying no other submitter is mirrored gives way to a table.

    Measured with this parser on builds of `git archive` copies
    (`PYTHONPATH=<copy>/src`), 119 gene pages each:

    =========================================  =====  =====  ======================
    build                                      (a)    cells  (b)
    =========================================  =====  =====  ======================
    cb326ca, `dev` before WS3 PR iii           1,327  1,456  32 sentences, 756 words
    0fc451a, after Tasks 1-5 (#60), #62, #63   1,066  1,194  25 sentences, 506 words
    452625c, after Task 6                        836    952  20 sentences, 386 words
    802d694, the gene card's link to methods     841    957  21 sentences, 395 words
    2f6300a, owner decisions 2026-10-26          895  1,016  22 sentences, 432 words
    9a8f14d, review fixes 2026-10-27             882  1,001  21 sentences, 404 words
    owner decisions 2026-10-26, batch 2 (R1-I3)  942  1,051  21 sentences, 415 words
    review fixes to batch 2, 2026-10-27          957  1,066  21 sentences, 430 words
    owner decisions batch 3 (R12, R2)            957  1,066  21 sentences, 430 words
    review fixes to batch 3, 2026-10-27          957  1,066  21 sentences, 430 words
    owner answers 2026-10-27 (A1, A2)            964  1,073  21 sentences, 430 words
    review fixes to owner answers, 2026-10-27    964  1,073  21 sentences, 430 words
    =========================================  =====  =====  ======================

    The 2f6300a row is the owner's wording decisions of 2026-10-26 (E1, E2, R4,
    R6, R7), measured on a build of the clean tree at each commit: 3b766ba
    read 841 and 21 / 395, as at 802d694. E1 lifted two particular caveats
    out of the chart folds (+12, to 853); E2 dropped "whole" (-2, to 851); R4's
    matrix caption sentence is on all 118 pages with a matrix (+38, to 889) and
    is the 22nd boilerplate sentence (38 words); R7's declaration on the
    card's tested line is on 86 cards (+6, to 895). R6's line is on one card
    and E3's section intro on one page, and neither moves either figure; R10,
    W2, W4 and W5 change no gene page's visible prose (W4's and W5's cohort
    descriptions render folded). Growth from owner-decided sentences is deliberate (E4), and
    none of them was folded or cut to fit.

    The row after it is the review of those decisions (2026-10-27), measured
    on a fresh build of the clean tree at that commit against one of ad03128
    (895, 1,016, 22 / 432, as at 2f6300a). R4's caption now renders only where
    the card lists a surviving row, so it left 86 pages (-38 each) and the
    boilerplate (-1 sentence, -38 words); the own-lab declaration's corrected
    wording is 10 words longer in each own-lab block (+10 or +20 on 118 pages,
    +10 boilerplate words). Per page: -28 on 27, -18 on 40, +16 on 16, and
    smaller moves elsewhere; the median fell by 13. Both constants were lowered
    to the new measurement, not left at their old headroom.

    The last row is the second batch of owner decisions (R1, R3, R5, R8, R9,
    R11, R13, I3; 2026-10-26/27), measured on fresh builds of the clean tree
    after each item against one of e0393f4 (882, 1,001, 21 / 404, as at
    9a8f14d). Two items moved the figures, both owner-decided sentences (E4):
    R5's card line naming the ClinGen panel that scoped a syndrome, on 10
    cards (+4, to 886); and I3, which merged PMID:34324492 and PMID:42230622
    into one cohort family -- both legends gained the clause saying a family
    also merges across collections recorded as overlapping, 11 words on the
    118 matrix legends (+11, to 897; boilerplate words 404 to 415, still 21
    sentences), and its notice naming the declared DDD-DECIPHER overlap is 45
    words on the 72 pages carrying both studies (+45, to 942; 897 is the
    median of the same build with that notice's words subtracted). R1, R3,
    R8, R9, R11 and R13 moved neither figure: R3's 35 card lines sit on 25
    pages, below the median's reach, and R8 changes only hrefs. The notice is
    particular to the page's data, so the fold rule keeps it visible; (a) was
    raised to 971, 942 x 1.03 rounded up, and (b) stays 22, one over 21.

    The last row is the review of that batch (2026-10-27), measured on a fresh
    build of the fixed tree against one of 5f7161b (942, 1,051, 21 / 415, as
    the row above). Two boilerplate sentences grew, both to stay true since
    I3: the pooling caption names the declared-overlap reason beside the
    shared-collection one (+13 words), and the matrix key's `corrected`
    caption no longer says "that study" of a family that may hold two (+2).
    So +15 on 117 of the 118 pages with a matrix, +17 on the 118th, EVC2's
    (its card's scope line now names the gene its panel graded, +2), and none
    on CRIPTO's, which has no matrix; the median rose by 15. (a) was
    reset to 986, 957 x 1.03 rounded up, and (b) stays 22, one over 21.

    The row after it is the third batch (R12, R2; decided 2026-10-26),
    measured on 2026-10-27 with this parser on fresh builds of the tree after
    each item against one of a12b897 (957, 1,066, 21 / 430). Neither item
    moved any figure. R12's seven newly registered case series are named only
    inside two folds on the 77 pages carrying PMID:34324492 (its block's "How
    this study counted" and the cohort notes). R2's card sentence,
    owner-decided (E4), is 14 words on 10 cards (+14 each, measured) and 0 on
    the other 109, below the median's reach. Both constants stay.

    The "review fixes to batch 3" row is the review of that batch
    (2026-10-27), measured with this
    parser on a fresh build of the fixed tree against one of 64ac83a (957,
    1,066, 21 / 430, as the row above). Only R2's card sentence moved: it
    gained "so none is an aneurysm proband", +6 words on each of the same 10
    cards and 0 on the other 109, below the median's reach. By the rule,
    957 x 1.03 rounded up is 986 and 21 + 1 is 22, so both constants stay.

    The "owner answers" row is the owner's answers of 2026-10-27, measured with this
    parser on fresh builds of the tree after each item against one of 3bf08fb
    (957, 1,066, 21 / 430, as the row above). A1's review statuses render
    only inside the cohort-notes fold (two notes on each of the 77 pages
    carrying PMID:34324492) and moved nothing. A2's qualifier on the card's
    "Tested in N independent dataset(s)" line, owner-decided (E4), is 7 words
    on 86 cards and 0 on the other 33 (measured per page), and the median
    rose by 7 to 964. (b) did not move: the line is on 86 of 119 pages, under
    the 90% share, and reads differently by count and declaration. By the
    rule (a) was raised to 993, 964 x 1.03 rounded up, and (b) stays 22, one
    over 21.

    The last row is the review of those answers (2026-10-27), measured with
    this parser on a fresh build of the fixed tree against one of ba2b2f9
    (964, 1,073, 21 / 430 on both). Its gene-page changes are the wording of
    the two cohort notes, inside the cohort-notes fold on the 77 pages
    carrying PMID:34324492, and they moved nothing; the sentence saying what
    "Pending review" marks is on methods.html, not on a gene page. By the
    rule, 964 x 1.03 rounded up is 993 and 21 + 1 is 22, so both constants
    stay.

    "cells" is the history: (a) with cells counted, the figure this test
    asserted until review, 957 at zero headroom. (b) never counted them, and
    the review's fixes to methods.html left every gene page's bytes as at
    802d694. The card's last line, "How to read this page · this gene as
    JSON", is one of the 21: nine words, on every page. The plan's baseline
    read 1,460 on cb326ca with the inventory's parser, which reads `<br>` as a
    space; under that reading the cells column is 4 higher on each build
    (1,460, 1,198, 956 and 961) and no other figure moves, since every `<br>`
    on a gene page sits in a table cell.

    **Headroom, because a curation-only change must not fail the suite.**
    (a) is held at 993, the last row's 964 x 1.03 rounded up, and (b) at 22,
    one more than 21 (owner decision E4 keeps both margins; they were 986 and
    22 over the three rows before it, 971 and 22 over the batch-2 row, 909 and
    22 over 9a8f14d's, 922 and 23 over 2f6300a's, and 867 and 22 over
    802d694's). **When this fails,
    re-measure before anything else.** If the growth is a new study's or
    dataset's own data-driven prose, raise the budget and put the measurement
    here. Never fold a caveat the fold rule
    requires visible to get under it: the budget exists to keep explanation
    off the page, not to take a warning out of sight.

    **Proved able to fail** (2026-10-26, each a build of a `git archive` copy
    of the branch, measured with this parser): the full scope rule
    (`methods.SCOPE_RULE`) rendered after the short form on every gene page
    reads 1,144 and 27 sentences, over both budgets; the grade key back on the
    page 986 and 27; the short form replaced by the full rule 1,112 and 25;
    the trajectory's toggle unfolded 889 and 22; the small multiples' 948 and
    24. Two now pass, and each is killed by a test of its own: the
    rare-variant reading notes re-added under a fold, 846 and 21 (a summary
    line is all that shows), and the pooling notice's long form, 863 and 22.
    That is the cost of headroom, and why this test is a ceiling on growth
    rather than a pin on any one sentence.

    The first assertions check the instrument on input this test writes, so
    a broken parser, word count, sentence split or share cannot pass the
    budget on any pages: with `_WORD` matching nothing every page reads 0
    words, and the short form replaced by the full scope rule passed this
    test (review, 2026-10-26). The real pages' median must also be positive.
    """
    assert _visible_blocks(
        "<html><head><title>no</title></head><body><p>One two.</p><details>"
        "<summary>Three</summary><p>four</p><details><summary>five</summary></details>"
        "</details><table><tr><td>six</td></tr></table><svg><text>seven</text></svg>"
        "<dl><dt>eight</dt><dd>nine</dd></dl></body></html>"
    ) == [("One two.", False), ("Three", False), ("six", True), ("eight", True), ("nine", True)]
    assert _prose("<p>One two.</p><table><tr><th>three</th><td>four</td></tr></table>") == [
        "One two."
    ]
    assert _words("τ = 0.71 — across ∞ organs · here") == 6
    assert _SENTENCE.split("First one here. Second (one) here? Third.") == [
        "First one here.",
        "Second (one) here?",
        "Third.",
    ]
    # (b)'s rule on ten pages written here: six words on 9 of 10 is boilerplate
    # and on 8 is not; five words on all 10 is not; a sentence in a cell is not.
    on_nine, on_eight = "Nine pages of ten say this.", "Eight pages of ten say this."
    synthetic = [
        _prose(
            f"<p>{on_nine if n >= 1 else ''} {on_eight if n >= 2 else ''} Every page says "
            "these five.</p><table><tr><td>This sentence sits in a table cell on all ten."
            "</td></tr></table>"
        )
        for n in range(10)
    ]
    assert _boilerplate(synthetic) == [on_nine]

    pages = [
        _prose(path.read_text(encoding="utf-8"))
        for path in sorted((site / "genes").glob("HGNC_*.html"))
    ]
    assert len(pages) > 1, "the build published no gene pages; the budget would pass vacuously"
    median = statistics.median([sum(_words(text) for text in blocks) for blocks in pages])
    boilerplate = _boilerplate(pages)
    # Both figures in either failure, so one run says how far over each is.
    report = (
        f"median visible words per gene page {median} (budget {_MEDIAN_VISIBLE_WORDS}); "
        f"{len(boilerplate)} boilerplate sentences (budget {_BOILERPLATE_SENTENCES}): "
        f"{boilerplate}"
    )
    assert median > 0, report
    assert median <= _MEDIAN_VISIBLE_WORDS, report
    assert len(boilerplate) <= _BOILERPLATE_SENTENCES, report
