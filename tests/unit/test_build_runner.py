# tests/unit/test_build_runner.py
import gzip
import hashlib
import json
import re
import shutil
import subprocess
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.burden import BurdenRow, load_burden
from chd_atlas.build.cell_types import CellTypeCensus
from chd_atlas.build.concordance import survives_own_correction
from chd_atlas.build.identity import GeneIdentity, read_gene_identities
from chd_atlas.build.profiles import ProfileCensus
from chd_atlas.build.runner import BuildRefused, _gene_registry, build_site
from chd_atlas.build.validity import gene_validity, published_genes
from chd_atlas.corpus import load_curation
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.runner import ValidationReport

REPO = Path(__file__).parent.parent.parent

# Every artifact the published API promises. Listed rather than derived, because
# this is the contract: a builder dropped from `build_site` would otherwise take
# its artifact out of the site and out of the assertion in the same edit.
PUBLISHED = (
    "manifest.json",
    "index.html",
    "genes/index.json",
    "publications.json",
    "featured.json",
    "phenotypes.json",
    "datasets.json",
    "sources.json",
    "variants/index.json",
    "search/index.json.gz",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A copy of the real committed corpus, which validates cleanly.

    `curation` and `mirrors` are copied because tests edit them; `ontologies` is
    symlinked because it is 11 MB of pinned OBO releases that nothing here
    changes, and copying it per test would cost more than the builds do.
    """
    source = tmp_path / "repo"
    source.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, source / name)
    (source / "ontologies").symlink_to(REPO / "ontologies")
    return source


def _break(repo: Path) -> None:
    """Make the corpus fail validation with one SCHEMA001.

    An empty `assertions` list, because `AssertionFile` requires at least one
    record — a file that parses as YAML and fails as a record, so the failure is
    the validator's rather than the loader's.
    """
    (repo / "curation" / "assertions" / "BROKEN.yaml").write_text("assertions: []\n")


def _search_records(out: Path) -> list[dict[str, Any]]:
    raw = (out / "search" / "index.json.gz").read_bytes()
    records: list[dict[str, Any]] = json.loads(gzip.decompress(raw))["records"]
    return records


def test_the_registry_reader_cleans_the_alias_cell_and_falls_a_blank_symbol_back_to_the_hgnc_id(
    tmp_path: Path,
) -> None:
    """Every way a mirror cell reaches this reader other than as a curator wrote it.

    Read directly rather than through `build_site`, because both inputs here are
    what the gate exists to prevent and one of them the gate does not catch.

    Whitespace and empty alias segments are reachable on a fully green build:
    nothing validates the contents of the `aliases` cell, so
    `T-box 5|  TBX5B  ||Chr12` is a legal mirror row, and unstripped it
    publishes a search term with leading spaces plus an empty term that every
    "contains" query matches.

    A blank `symbol` is not reachable behind `build_site` either way — `symbol`
    is non-nullable, so TBL003 refuses the build first — but a single space in
    the cell is: `read_table` maps only the true empty string to null, so
    `HGNC:9999\t \t...` is a legal mirror row TBL003 does not catch, and
    `HGNC:4173\t\t...` (truly empty) reaches here identically once `cell_text`
    strips it. This is the regression the fix corrects: the previous guard read
    `if gene is None or symbol is None: continue` and dropped both rows
    outright, discarding a validated `name` along with the blank symbol —
    watched failing here before the fix, with `registry` holding only
    `HGNC:11604` where it now holds all three keyed rows. The symbol now falls
    back to the HGNC id instead, the same fallback `build_genes` and
    `search.py` already apply to a gene missing from this registry entirely, so
    the row's `name` and `aliases` reach the site rather than being discarded
    with it.

    A row with no `hgnc_id` at all is the one case still, correctly, dropped:
    there is no gene to key it on and no row to correct it against, so `name`
    and `aliases` have nowhere to attach either.
    """
    mirror = tmp_path / "mirrors"
    mirror.mkdir()
    (mirror / "genes.tsv").write_text(
        "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
        "HGNC:11604\tTBX5\tT-box transcription factor 5\tT-box 5|  TBX5B  ||Chr12\t\t\t\t\t\n"
        "HGNC:4173\t\tGATA binding protein 4\t\t\t\t\t\t\n"
        # A single space, which TBL003 does not catch: `read_table` maps only the
        # empty string to null.
        "HGNC:9999\t \tSome other gene\t\t\t\t\t\t\n"
        "\tNONAME\tUnkeyed row\t\t\t\t\t\t\n"
    )

    registry = _gene_registry(tmp_path)

    assert list(registry) == ["HGNC:11604", "HGNC:4173", "HGNC:9999"]
    assert registry["HGNC:11604"].aliases == ("T-box 5", "TBX5B", "Chr12")
    assert registry["HGNC:11604"].name == "T-box transcription factor 5"
    assert registry["HGNC:11604"].symbol == "TBX5"

    # The regression this commit fixes: a blank symbol must not discard the
    # name beside it, whichever of the two blank shapes it arrives as.
    assert registry["HGNC:4173"].symbol == "HGNC:4173"
    assert registry["HGNC:4173"].name == "GATA binding protein 4"
    assert registry["HGNC:9999"].symbol == "HGNC:9999"
    assert registry["HGNC:9999"].name == "Some other gene"


def test_the_build_publishes_every_artifact_and_a_manifest_that_verifies(
    repo: Path, tmp_path: Path
) -> None:
    """The whole contract in one build, checked against the bytes on disk.

    Every checksum is verified rather than merely shaped, because that is the
    manifest's entire purpose and A43 records what the plan's suite missed: with
    every value replaced by 64 zeros, tests asserting only the `sha256:` prefix
    all passed. A manifest whose checksums are wrong is worse than none — a
    consumer verifying a download rejects a file that is correct.

    The file set is checked in both directions. A file on disk that the manifest
    omits is served with no checksum and invisible to any consumer that
    verifies; a manifest entry with no file is a published URL that 404s. Only
    `manifest.json` may be in one and not the other, since it cannot contain its
    own checksum.
    """
    out = tmp_path / "dist"

    written = build_site(repo, out)

    for relative in PUBLISHED:
        assert relative in written, f"{relative} is promised by the API and was not built"

    manifest = json.loads((out / "manifest.json").read_text())
    for relative, advertised in manifest["files"].items():
        actual = hashlib.sha256((out / relative).read_bytes()).hexdigest()
        assert advertised == f"sha256:{actual}", f"{relative} does not match its checksum"

    on_disk = {path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()}
    assert on_disk == set(manifest["files"]) | {"manifest.json"}

    # The return value accounts for the manifest as well as everything in it, so
    # a caller checking what was built does not have to know that one artifact
    # is missing from the manifest by construction.
    assert set(written) == on_disk


def test_a_build_from_a_checkout_stamps_the_manifest_with_that_commit(
    repo: Path, tmp_path: Path
) -> None:
    """The provenance branch that only ever runs in the deploy, and nowhere else.

    Every other test here builds from a `tmp_path` copy that is not a git
    repository, so `source_commit` returns `None` in all of them and the manifest
    records no commit. That is the opposite of what happens in CI, where the
    build always runs inside a checkout — so the wiring the published site
    actually depends on was reachable by no test, and a regression in how
    `build_site` obtains the commit would have shipped green.

    Made a checkout here rather than asserting against this repository's own
    HEAD, which would pass just as well against a `source_commit` that ignored
    `root` entirely.
    """
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "corpus"],
        cwd=repo,
        check=True,
    )
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True, check=True
    ).stdout.strip()
    out = tmp_path / "dist"

    build_site(repo, out)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["source_commit"] == head
    assert len(head) >= 40


def test_a_repository_that_fails_validation_is_refused_before_anything_is_written(
    repo: Path, tmp_path: Path
) -> None:
    """Publishing from a corpus that fails its own gate defeats having one.

    Two halves, and the second is the one that matters on a deploy: the failure
    must name the validation errors, so a curator who ran `build` knows what to
    fix without being told to go and run a second command; and nothing may be
    written, not even an empty directory, because a half-built `dist/` is
    something a deploy step will happily upload.
    """
    _break(repo)
    out = tmp_path / "dist"

    with pytest.raises(BuildRefused) as raised:
        build_site(repo, out)

    assert "SCHEMA001" in str(raised.value)
    assert not out.exists()


@pytest.mark.parametrize(
    ("code", "severity", "refuses"),
    [("SCHEMA001", Severity.ERROR, True), ("REF000", Severity.WARNING, False)],
)
def test_the_gate_refuses_on_an_error_and_builds_through_a_warning(
    repo: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    code: str,
    severity: Severity,
    refuses: bool,
) -> None:
    """Which condition the gate reads, decided here and not reimplemented.

    A stricter gate makes a repository CI calls green unpublishable; a looser one
    publishes a corpus the validator rejected. `ValidationReport.ok` ignores
    warnings by design, and that half cannot be provoked from a fixture: every
    warning this project emits is a "skipped checks" notice that fires *because*
    something else errored, so no corpus produces warnings alone. Leaving it
    untested is how a gate tightened to `warning_count == 0` reaches `main`.

    The report is therefore supplied rather than provoked, which makes this test
    about the condition alone. That the real validator's verdict is what reaches
    the gate is a different claim, pinned by the two tests above: a genuinely
    broken corpus refuses and names its SCHEMA001, and a genuinely clean one
    builds. Running the real validator again here cost 3 of the file's 11
    seconds and duplicated both.
    """
    report = ValidationReport(issues=[ValidationIssue(code, severity, str(repo), "constructed")])
    assert report.ok is not refuses
    monkeypatch.setattr("chd_atlas.build.runner.validate_repository", lambda _: report)
    out = tmp_path / "dist"

    if refuses:
        with pytest.raises(BuildRefused, match=code):
            build_site(repo, out)
        assert not out.exists()
    else:
        build_site(repo, out)
        assert (out / "manifest.json").is_file()


def test_the_censuses_reach_the_manifest_and_the_page_by_the_rules_handed_in(
    repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`profile_census`'s own published-vs-registry restriction is proven by a
    hand-built fixture in `test_build_profiles.py`; this proves the *wiring*
    `build_site` does with its answer, which no test reading a real build can.
    This said "the committed corpus has no `profiles` mirror, so both figures
    are 0 there" until 2026-08-20; `mirrors/profiles/E-MTAB-6814.tsv` landed on
    2026-08-19 and a real build now publishes `profile_genes: 92` and
    `profile_datasets: 1`. That does make a *transposition* visible on a real
    build, but nothing there asserts either value, and the mutant that matters
    is still invisible: 92 is also `counts.genes` and 1 is also
    `counts.datasets`, so an assignment reading the neighbouring count instead
    of `profile_stats` publishes exactly the same manifest and page as a
    correct build. Genes and datasets are given different values here for
    the same reason `BurdenCensus`'s three figures are all distinct in
    `test_the_burden_census_reaches_the_page_once_per_figure`: two figures
    equal to each other cannot catch a swap between their own keys.

    `profile_census` itself is monkeypatched to a fixed answer rather than
    given real mirror data -- the same technique
    `test_the_gate_refuses_on_an_error_and_builds_through_a_warning` already
    uses on `validate_repository` in this file -- so this test is cheap and
    about exactly one thing: `build_site` reads the return and publishes
    `genes` under `profile_genes` and `datasets` under `profile_datasets` in
    `manifest.json`, never the other way round, and `index.html` states the
    first.
    """
    monkeypatch.setattr(
        "chd_atlas.build.runner.profile_census",
        lambda profiles, published: ProfileCensus(genes=7, datasets=4),
    )

    # **And the surviving-genes census counts by the card's predicate** (review,
    # 2026-10-27). `build_site` hands `burden_census` `survives_own_correction`,
    # and swapping it for "a corrected p below 0.05" survived every test: no
    # committed row tells the two apart. KDM6A -- published, no surviving row --
    # is given two rows that do: its upward synonymous row and a depleted row,
    # each with a corrected p of 0.01. Neither supports the gene (a synonymous
    # row is the negative control; a depleted one points the wrong way), so the
    # census must not count it.
    kdm6a = "HGNC:12637"
    real_load = load_burden
    loaded: dict[str, list[BurdenRow]] = {}

    def with_decoys(root: Path) -> dict[str, list[BurdenRow]]:
        burden = {gene: list(rows) for gene, rows in real_load(root).items()}
        q = {"pvalue": 0.001, "pvalue_adjusted": 0.01, "pvalue_adjustment": "benjamini_hochberg"}
        burden[kdm6a] = [
            replace(row, effect=2.0, **q)  # type: ignore[arg-type]
            if row.consequence_class == "synonymous" and row.cohort_stratum == "all"
            else replace(row, **q)  # type: ignore[arg-type]
            if row.effect is not None and row.effect < 1 and row.pvalue_adjusted == 1.0
            else row
            for row in burden[kdm6a]
        ]
        loaded.update(burden)
        return burden

    monkeypatch.setattr("chd_atlas.build.runner.load_burden", with_decoys)
    out = tmp_path / "dist"

    build_site(repo, out)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["counts"]["profile_genes"] == 7
    assert manifest["counts"]["profile_datasets"] == 4

    page = (out / "index.html").read_text()
    assert re.search(r"<dt>Genes with developmental expression data</dt>\s*<dd>7</dd>", page)
    # The datasets card went with the doorway page (2026-10-27); the manifest
    # above states that count.
    assert "Developmental expression datasets" not in page

    published = [
        entry["gene"] for entry in json.loads((out / "genes/index.json").read_text())["genes"]
    ]
    surviving = sum(
        1 for gene in published if any(survives_own_correction(row) for row in loaded.get(gene, ()))
    )
    by_q_alone = sum(
        1
        for gene in published
        if any(
            row.pvalue_adjusted is not None and row.pvalue_adjusted < 0.05
            for row in loaded.get(gene, ())
        )
    )
    # The decoys separate the two predicates, or this measures one number twice.
    assert kdm6a in published and by_q_alone == surviving + 1
    assert manifest["counts"]["burden_surviving_genes"] == surviving


def test_every_gene_label_the_registry_holds_reaches_the_site(repo: Path, tmp_path: Path) -> None:
    """Symbol, approved name and aliases, from `mirrors/genes.tsv` to both payloads.

    TBX5's committed `aliases` cell is null, so the split is exercised only if
    this test supplies one.
    That matters more than it looks: the cell is pipe-separated and handing it
    over unsplit type-checks — `frame.to_dicts()` yields `dict[str, Any]`, and
    an `Any` satisfies `tuple[str, ...]` silently — so the mistake this asserts
    against is one mypy cannot catch. A reader that skipped the split publishes
    one term reading "T-box 5|Chr12q24.1", and a reader that passed the raw
    string publishes a term per character.

    The name is asserted because it is the column every valid mirror is required
    to have (TBL003) and the one A38 found unreachable: a visitor typing
    "T-box transcription factor" must find TBX5.

    `prev_symbols` is asserted beside `aliases` since schema 2.13, because it
    is the second pipe-separated cell the same reader splits and it fails in
    exactly the same two ways. The expected list pins the **order** as well as
    the membership: aliases sort among themselves and come first, retired names
    sort among themselves and follow. That is the precedence `genes.py`
    resolves by, so the index and the resolver cannot disagree about which tier
    a string belongs to, and a membership check would pass with the two
    interleaved.

    The row is located by id rather than taken as the first line, and the rest
    of the mirror is left intact. Rewriting the file down to its first row
    worked only while the registry held TBX5 alone; against the 154-gene mirror
    it deleted the asserted gene and the build refused with REF001.

    The search record is located by id for the same reason, and the reason is
    newer: D31 keyed the index on `published` rather than on the assertions, so
    it went from one gene record to 23 sorted by HGNC id, and `next(...)` on the
    first one started returning TBX20 (HGNC:11598). It was reading a real
    payload, just never the row this test names — the shape of failure that a
    positional lookup into a growing array always has.
    """
    mirror = repo / "mirrors" / "genes.tsv"
    header, *rows = mirror.read_text().splitlines()
    columns = header.split("\t")
    aliases, previous = columns.index("aliases"), columns.index("prev_symbols")

    def with_names(row: str) -> str:
        cells = row.split("\t")
        if cells[0] != "HGNC:11604":
            return row
        cells[aliases] = "T-box 5|Chr12q24.1"
        # Set here rather than left at the mirror's own value, so the expected
        # list below is a fact about this fixture and not about whatever HGNC
        # last withdrew from TBX5. Two retired names, both distinct from the two
        # aliases, so the assertion can tell the columns apart -- with one value
        # each, "read prev_symbols" and "read aliases twice" publish the same
        # list (CLAUDE.md section 4.36).
        #
        # **Written out of sorted order deliberately.** These arrived as
        # `HOS|HOS2`, which is what the mirror would plausibly hold and is
        # already sorted -- so dropping the `sorted()` in `search.py` published
        # the identical list and the mutant survived the whole matrix. A cell
        # whose own order differs from its sorted order is the only fixture
        # that can tell a sort from a pass-through.
        cells[previous] = "HOS2|HOS"
        return "\t".join(cells)

    edited = [with_names(row) for row in rows]
    assert edited != rows, "TBX5 left the registry; this test would assert nothing"
    mirror.write_text("\n".join((header, *edited)) + "\n")
    out = tmp_path / "dist"

    build_site(repo, out)

    index = json.loads((out / "genes" / "index.json").read_text())
    symbols = {entry["gene"]: entry["symbol"] for entry in index["genes"]}
    assert symbols["HGNC:11604"] == "TBX5"
    # Widened rather than dropped when D21 took the index from 1 row to 23: this
    # line used to read `== ["TBX5"]`, which covered the one published gene there
    # was. `bundles.py` falls back to the HGNC id for a gene absent from the
    # registry, so a row whose symbol *is* its id is a label the mirror was
    # supposed to hold and the site did not get — the same loss the assertion
    # above pins for one gene, checked over every gene the build publishes.
    assert [gene for gene, symbol in symbols.items() if symbol == gene] == []

    gene = next(
        record
        for record in _search_records(out)
        if record["kind"] == "gene" and record["id"] == "HGNC:11604"
    )
    assert gene["label"] == "TBX5"
    assert gene["terms"] == [
        "TBX5",
        "HGNC:11604",
        "T-box transcription factor 5",
        "Chr12q24.1",
        "T-box 5",
        "HOS",
        "HOS2",
    ]


def test_the_site_publishes_exactly_the_genes_the_gate_selects(repo: Path, tmp_path: Path) -> None:
    """`build_site` must publish `published_genes()`'s answer, not another set.

    Design decision D21 lives in `build/validity.py`, and
    `test_build_validity.py` pins what it returns over the committed mirrors --
    23 genes. Nothing pinned that `build_site` *uses* that return, and the gap
    was not academic: replacing `published = published_genes(validity)` with
    `published = set(genes)` -- the 154-gene registry -- published 154 rows and
    154 bundles, presenting 131 genes no ClinGen panel calls definitive inside
    the definitive browse set, and the whole 590-test suite passed.

    That is the charter's worse failure rather than its usual one. The
    characteristic defect here is curated work reaching no page; this is its
    inverse, a claim the sources do not support reaching every page, and "a
    wrong claim here is worse than a missing one".

    The expected set is recomputed from the same mirrors rather than written as
    23 literals, because a literal list would have to be re-typed whenever
    ClinGen curates another CHD gene and would fail as a mirror refresh rather
    than as a defect. The *count* is pinned in `test_build_validity.py`, which is
    where a moved number should be looked at; this test pins that the build and
    the gate agree, which is a different property and needs its own guard --
    `test_the_published_gene_count_agrees_with_a_real_build_of_genes_index_json`
    cannot cover it, because both of its figures derive from one object and move
    together under exactly this mutant.
    """
    out = tmp_path / "dist"

    build_site(repo, out)

    corpus, _ = load_curation(repo)
    clingen, _ = read_table(
        repo / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        repo / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    assert clingen is not None and gencc is not None
    expected = published_genes(
        gene_validity(clingen, gencc, in_scope={str(entry.id) for entry in corpus.chd_scope})
    )

    index = json.loads((out / "genes" / "index.json").read_text())
    assert {entry["gene"] for entry in index["genes"]} == expected
    # The bundles too, not only the index: the two are written by one loop today
    # and a reader of this test should not have to know that to trust it.
    assert {path.stem for path in (out / "genes").glob("*.json")} - {"index"} == {
        gene.replace(":", "_") for gene in expected
    }


def test_a_built_site_carries_a_page_for_every_published_gene(repo: Path, tmp_path: Path) -> None:
    """`build_site` is the only place the published population reaches the page builders.

    Nothing in `src/` imported `pages.py` until this call existed — only its own
    unit test did — while `build_landing` and the shared `<nav>` already linked
    to `genes/index.html`. That is a green build with every checksum verifying
    and every visitor who clicks "Genes" served a 404: the work reaching no
    page, which is the failure this project is written against.

    Both directions are asserted. A missing page is a dead link from a row the
    browse payload published; an extra page is a gene reachable by URL that the
    index does not list, which is how a de-published gene stays up. The expected
    names are derived from each row's own `bundle` rather than from a second
    `slug()` call, so the test cannot reproduce a page-naming bug and agree with
    it.

    The count is deliberately not re-pinned here — `test_build_validity.py` and
    `test_build_landing.py` already hold the figure of 23, and a third copy would
    fail as a mirror refresh rather than as a defect. What is guarded instead is
    that the index is non-empty, without which every assertion below is vacuous.

    Mutation matrix, measured 2026-08-04 against the full suite (616 tests), one
    fresh process per mutant. Each was killed by this test and by no other:

    * `build_gene_pages` dropped from `build_site` — 1 failed, 615 passed.
    * `build_gene_index_page` dropped from `build_site` — 1 failed, 615 passed.
    * `symbols={}` passed to `build_gene_index_page` — 1 failed, 615 passed,
      which is what earns the symbol loop its place beside the file-set
      assertion. `test_every_gene_label_the_registry_holds_reaches_the_site`
      does not cover it: that one reads `genes/index.json` and the search index,
      neither of which this builder writes.
    """
    out = tmp_path / "dist"

    build_site(repo, out)

    index = json.loads((out / "genes" / "index.json").read_text())["genes"]
    assert index, "an empty index would make every assertion below assert nothing"

    expected = {Path(row["bundle"]).name.replace(".json", ".html") for row in index}
    assert {path.name for path in (out / "genes").glob("*.html")} == expected | {"index.html"}

    browse = (out / "genes" / "index.html").read_text(encoding="utf-8")
    for row in index:
        assert row["symbol"] in browse, f"{row['symbol']} is published but is not browsable"


def test_the_site_carries_the_terms_of_everything_it_republishes(
    repo: Path, tmp_path: Path
) -> None:
    """Attribution has to travel with the data, not sit in a README.

    `mirrors/sources.yaml` records HPO as `permitted_with_attribution`, and the
    build republishes its term labels and synonyms verbatim in `phenotypes.json`
    and in the search index. A consumer fetching JSON never opens the README, and
    the repository's only `LICENSE` is Apache-2.0 — which a reader will
    reasonably take to cover the published files too.

    So the obligation is asserted where it is discharged: every source the
    registry knows appears in the published payload, with its licence and its
    redistribution terms intact.
    """
    out = tmp_path / "dist"

    build_site(repo, out)

    published = json.loads((out / "sources.json").read_text())["sources"]
    by_id = {source["id"]: source for source in published}
    assert "hpo" in by_id, "the ontology whose content is republished must be named"
    assert by_id["hpo"]["redistribution"] == "permitted_with_attribution"
    assert by_id["hpo"]["licence"], "a source with terms must carry them"
    for source in published:
        assert source["name"] and source["url"] and source["licence"]
    # Ordered, because this is a published array like every other.
    assert [source["id"] for source in published] == sorted(by_id)


def test_the_registry_reads_each_genes_identity_with_blank_cells_as_none(tmp_path: Path) -> None:
    """D52's ids come from one row of `mirrors/genes.tsv`, cleaned like the labels are.

    A whitespace-only cell is absent, by the same `cell_text` rule the symbol uses;
    a row with no `hgnc_id` has nothing to key on and is dropped. `ncbi_gene`
    must arrive as an `int` -- the bundle publishes it as one.

    `read_gene_identities` is the one reader the build and `chd-atlas validate`
    (GEN004) both use since review, 2026-10-27 -- validate held a copy of its
    own until then -- so a padded id is read the same way by both: stripped.
    A registry missing an id column is `None`, for TBL001 to report, never a
    mapping of genes with no ids for GEN004 to report one by one.
    """
    mirror = tmp_path / "mirrors"
    mirror.mkdir()
    (mirror / "genes.tsv").write_text(
        "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\t ENSG00000089225 \t6910\t12q24.21\t"
        "Q99593\tNM_181486.4\n"
        "HGNC:4173\tGATA4\tGATA binding protein 4\t\t\t\t \t\t\n"
        "\tNONAME\tUnkeyed row\t\t\t\t\t\t\n"
    )

    identities = read_gene_identities(tmp_path)

    assert identities is not None
    assert list(identities) == ["HGNC:11604", "HGNC:4173"]
    assert identities["HGNC:11604"] == GeneIdentity(
        name="T-box transcription factor 5",
        locus="12q24.21",
        ensembl_gene="ENSG00000089225",
        ncbi_gene=6910,
        uniprot="Q99593",
        mane_select="NM_181486.4",
    )
    assert type(identities["HGNC:11604"].ncbi_gene) is int
    assert identities["HGNC:4173"] == GeneIdentity(name="GATA binding protein 4")

    (mirror / "genes.tsv").write_text("hgnc_id\tsymbol\tname\tensembl_gene\nHGNC:1\tA\tB\tC\n")
    assert read_gene_identities(tmp_path) is None


def test_the_cell_type_census_reaches_the_manifest_by_key(
    repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`build_site` publishes `genes` under `cell_type_genes`, never swapped.

    Fixed, distinct answers (as `profile_census` is given above), because a
    real build has no cell-type mirror yet and publishes 0 for both; the empty
    case is asserted too, on the unpatched build.
    """
    out = tmp_path / "dist"
    build_site(repo, out)
    counts = json.loads((out / "manifest.json").read_text())["counts"]
    assert counts["cell_type_genes"] == 0
    assert counts["cell_type_datasets"] == 0

    monkeypatch.setattr(
        "chd_atlas.build.runner.cell_type_census",
        lambda expression, published: CellTypeCensus(genes=5, datasets=2),
    )
    patched = tmp_path / "patched"
    build_site(repo, patched)
    counts = json.loads((patched / "manifest.json").read_text())["counts"]
    assert counts["cell_type_genes"] == 5
    assert counts["cell_type_datasets"] == 2


# The cell-type layer on a real build. No committed dataset has `design:
# cell_type_profile` yet, so a build of the corpus alone publishes `{"datasets":
# []}` on every bundle and cannot see the layer: deleting `cell_types=` from
# either runner call passed the whole suite (review, 2026-10-28). One dataset
# over two published genes, so `cell_type_genes` (2) differs from both
# `cell_type_datasets` (1) and 0 (CLAUDE.md section 4.30). It carries the values
# the block's fixed copy states, which `cell_type_block` refuses to render
# without.
_HOC_DATASET = """\
id: hoc
archive: ucsc_cell_browser
technology: scrnaseq
design: cell_type_profile
tissue: heart
developmental_stage: 9 to 15 weeks post conception
organism: NCBITaxon:9606
n_samples: 8
licence: CC-BY-4.0
publication: PMID:38480880
cell_types:
  - {token: aCM, label: Atrial cardiomyocytes, compartment: Cardiomyocytes}
  - {token: Fibro, label: Fibroblasts, compartment: Stromal}
regions: [LA, LV]
detection_floor: 1.0
floor_source: test fixture
quantile_estimator: linear
"""
_HOC_PUBLICATION = """\
  - id: PMID:38480880
    title: Test fixture standing in for the cell-type dataset's publication.
    journal: Test fixture
    year: 2024
    authors:
      - Fixture A
    study_type: transcriptomics
"""
_HOC_SOURCE = """\
  - id: hoc_fixture
    name: Test fixture standing in for the cell-type pseudobulk
    version: "sha256:0000000000000000000000000000000000000000000000000000000000000000"
    retrieved_on: 2026-10-28
    url: https://cells.ucsc.edu/?ds=hoc
    licence: CC-BY-4.0
    redistribution: permitted_with_attribution
"""
# TBX5 and GATA4: each gene's pooled atrial-cardiomyocyte median, distinct so a
# bundle carrying the other gene's entry fails.
_HOC_ACM_MEDIANS = {"HGNC:11604": 220.67, "HGNC:4173": 41.48}


def _hoc_shard() -> str:
    rows = []
    for gene, acm in _HOC_ACM_MEDIANS.items():
        for cell_type, median in (("aCM", acm), ("Fibro", 11.16)):
            for region in ("LA", "LV", "all"):
                cells = 8000 if region == "all" else 4000
                rows.append(
                    (
                        gene,
                        cell_type,
                        region,
                        f"{median}\tcpm\t{median * 0.9}\t{median * 1.1}"
                        f"\t8\t{cells}\t0.5\thoc_fixture",
                    )
                )
    return (
        "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
        "n_donors\tn_cells\tfraction_expressing\tsource\n"
        + "".join(f"hoc\t{g}\t{c}\t{r}\t{rest}\n" for g, c, r, rest in sorted(rows))
    )


@pytest.fixture(scope="module")
def cell_type_build(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """The committed corpus plus the synthetic cell-type dataset above, built."""
    root = tmp_path_factory.mktemp("cell_type_build")
    repo = root / "repo"
    repo.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, repo / name)
    (repo / "ontologies").symlink_to(REPO / "ontologies")
    (repo / "curation" / "datasets" / "hoc.yaml").write_text(_HOC_DATASET, encoding="utf-8")
    for path, entry in (
        (repo / "curation" / "publications.yaml", _HOC_PUBLICATION),
        (repo / "mirrors" / "sources.yaml", _HOC_SOURCE),
    ):
        path.write_text(path.read_text(encoding="utf-8") + entry, encoding="utf-8")
    shard = repo / "mirrors" / "cell_type_profiles"
    shard.mkdir()
    (shard / "hoc.tsv").write_text(_hoc_shard(), encoding="utf-8")
    out = root / "dist"
    build_site(repo, out)
    return out


def test_a_cell_type_dataset_reaches_the_bundle_the_page_the_manifest_and_the_front_page(
    cell_type_build: Path,
) -> None:
    out = cell_type_build
    heading = '<section class="cell-types"><h3>Cell types (Farah 2024)</h3>'
    index = json.loads((out / "genes" / "index.json").read_text())
    published = [row["gene"] for row in index["genes"]]
    assert set(_HOC_ACM_MEDIANS) <= set(published)
    without = next(gene for gene in published if gene not in _HOC_ACM_MEDIANS)

    for gene, median in _HOC_ACM_MEDIANS.items():
        slug = gene.replace(":", "_")
        bundle = json.loads((out / "genes" / f"{slug}.json").read_text())
        (entry,) = bundle["cell_type_expression"]["datasets"]
        assert entry["dataset"] == "hoc"
        assert [c["cell_type"] for c in entry["cell_types"]] == ["aCM", "Fibro"]
        assert entry["cell_types"][0]["pooled"]["median"] == median
        page = (out / "genes" / f"{slug}.html").read_text()
        section = page.split("<h2>Developmental expression</h2>", 1)[1].split("<h2", 1)[0]
        assert heading in section

    other = without.replace(":", "_")
    assert json.loads((out / "genes" / f"{other}.json").read_text())["cell_type_expression"] == {
        "datasets": []
    }
    assert heading not in (out / "genes" / f"{other}.html").read_text()

    counts = json.loads((out / "manifest.json").read_text())["counts"]
    assert counts["cell_type_genes"] == len(_HOC_ACM_MEDIANS) == 2
    assert counts["cell_type_datasets"] == 1
    landing = (out / "index.html").read_text()
    assert "<dt>Genes with cell-type expression data</dt><dd>2</dd>" in landing
