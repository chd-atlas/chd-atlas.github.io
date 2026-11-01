# tests/test_docs_match_the_build.py
"""Pin the consumer documentation against what the build actually publishes.

`docs/data-api.md` opens by claiming every example is copied from a real build
rather than written by hand. Nothing checked that, and it drifted: the
`genes/index.json` example lost `conflicting_lesion_groups` — the field
`schema_version` 1.1 was minted for, and the one the document's own closing
section calls "the one consumer obligation this API imposes".

The consequence was not a missing line. A consumer typing their browse row off
an exhaustive, alphabetically-ordered example ships without the per-group
contested flag, and for a gene that is definitive-and-refuted on septal disease
while merely moderate on conotruncal, what they render inverts the science: the
definitive call reads as the solid one when it is precisely the disputed one.

So this is not a style check. It is the same guarantee the dead-link sweep gives
the payloads, applied to the document a consumer builds against.
"""

import csv
import gzip
import html
import inspect
import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.manifest import write_manifest
from chd_atlas.build.profiles import LOOKUP_RULE, TAU_METHOD, TAU_SCALE
from chd_atlas.build.runner import build_site
from chd_atlas.build.validity import GeneValidity, gene_validity
from chd_atlas.corpus import load_curation
from chd_atlas.models.dataset import Dataset, Stage
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.vocab import has_conflicting_evidence

REPO = Path(__file__).parent.parent
DOC = REPO / "docs" / "data-api.md"

# The doc writes small counts as words ("Six obligations"); a pin reads the
# number off the build and looks the word up here rather than hard-coding it.
_WORDS: dict[int, str] = dict(
    enumerate(
        (
            "zero",
            "one",
            "two",
            "three",
            "four",
            "five",
            "six",
            "seven",
            "eight",
            "nine",
            "ten",
            "eleven",
            "twelve",
        )
    )
)

# Each entry maps a documented example to the published object it claims to
# describe. Only *exhaustive* examples are listed — blocks written with an
# ellipsis are illustrations of shape, not of field sets, and are excluded here
# rather than silently half-checked.
#
# `heading` locates the fenced block; `select` pulls the object out of the built
# site that the block is an example of.
EXHAUSTIVE: dict[str, str] = {
    "## `manifest.json`": "manifest",
    "## `genes/index.json`": "index_row",
    "## `publications.json`": "publication",
    "## `featured.json`": "featured",
    "## `phenotypes.json`": "phenotype",
    "## `search/index.json.gz`": "search_record",
    "### The bundle's `validity` object: mirrored, attributed, never authored here": (
        "validity_record"
    ),
    "### The bundle's `burden_combined` array: one study's own combination of its arms": (
        "combined_row"
    ),
    "### The bundle's `burden` array: per study, never pooled": "burden_row",
    "## `scope.json`": "scope_term",
}


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    base = tmp_path_factory.mktemp("docs")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    out = base / "dist"
    build_site(root, out)
    return out


@pytest.fixture(scope="module")
def mirrored() -> dict[str, GeneValidity]:
    """Every gene the validity mirrors curate in scope, published or not.

    The doc states figures about the mirrors as well as the build -- how many
    genes set a flag, which genes carry only a `Supportive` record -- and the
    build publishes only the gate's side of them, so they are read here the
    way `build_site` reads them, before the gate.
    """
    corpus, _ = load_curation(REPO)
    clingen, _ = read_table(
        REPO / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        REPO / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    assert clingen is not None and gencc is not None
    return gene_validity(clingen, gencc, in_scope={str(entry.id) for entry in corpus.chd_scope})


def _bundles(site: Path) -> list[dict[str, Any]]:
    return [json.loads(path.read_text()) for path in sorted((site / "genes").glob("HGNC_*.json"))]


def _section(heading: str, end: str = "\n## ") -> str:
    """One section of the doc, whitespace collapsed: prose claims, not line breaks."""
    doc = DOC.read_text()
    start = doc.index(heading)
    stop = doc.find(end, start + len(heading))
    return " ".join(doc[start : stop if stop != -1 else len(doc)].split())


def _text_of(page: str) -> str:
    """One HTML page as its visible text, whitespace collapsed.

    Tags stripped rather than parsed: the assertion is about a sentence a
    reader sees, and `<strong>` inside it is exactly what a substring check on
    the raw markup would trip over.
    """
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", page))).strip()


def _first_block_quote(text: str) -> str:
    """The first `>` block in `text`, joined, stopping at the blank line after it."""
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            lines.append(line.lstrip("> ").strip())
        elif lines:
            break
    return " ".join(lines)


def _example(heading: str) -> dict[str, Any]:
    """The first fenced JSON block under one heading, parsed."""
    text = DOC.read_text()
    start = text.index(heading)
    nxt = text.find("\n## ", start + 1)
    section = text[start : nxt if nxt != -1 else len(text)]
    match = re.search(r"```json\n(.*?)```", section, re.S)
    assert match, f"no JSON example under {heading}"
    payload: dict[str, Any] = json.loads(match.group(1))
    return payload


def _published(site: Path, kind: str) -> dict[str, Any]:
    def load(rel: str) -> Any:
        raw = (site / rel).read_bytes()
        return json.loads(gzip.decompress(raw) if rel.endswith(".gz") else raw)

    if kind == "manifest":
        return dict(load("manifest.json"))
    if kind == "index_row":
        # TBX5 by id, not `genes[0]`. The example under that heading is TBX5's
        # row and the first row is now TBX20 (HGNC:11598) -- D21 took the index
        # from 1 gene to 23, sorted by HGNC id. Every row carries the same keys
        # by construction (`bundles.py::_headline` plus one literal dict), so
        # this test would still pass on the wrong row; selecting the documented
        # gene is what keeps it checking the example against what the example
        # claims to be.
        row = next(
            item for item in load("genes/index.json")["genes"] if item["gene"] == "HGNC:11604"
        )
        return dict(row)
    if kind == "publication":
        return dict(load("publications.json")["publications"][0])
    if kind == "featured":
        return dict(load("featured.json")["featured"][0])
    if kind == "phenotype":
        return dict(load("phenotypes.json")["phenotypes"][0])
    if kind == "search_record":
        return dict(load("search/index.json.gz")["records"][0])
    if kind == "validity_record":
        # TBX5 (HGNC:11604) is the committed corpus's one asserted gene, and its
        # first record is ClinGen's -- `build.validity._sort_key` sorts
        # `(source, disease, moi, submitter, sgc_id)`, and "clingen" precedes
        # "gencc" lexically.
        return dict(load("genes/HGNC_11604.json")["validity"]["records"][0])
    if kind == "scope_term":
        # The documented example is an Orphanet-path subtype, so the published
        # term compared against it is the same one: every key, warrant included.
        return dict(next(t for t in load("scope.json")["terms"] if t["id"] == "MONDO:0012547"))
    if kind == "burden_row":
        # TAB2, the documented row's gene. Its first row is PMID:34324492's,
        # not the documented one, which the key-set comparison does not mind:
        # every row carries every key (`burden.burden_payload`).
        return dict(load("genes/HGNC_17075.json")["burden"][0])
    if kind == "combined_row":
        # TBX5 again: it carries PMID:40127276's rows, so its first combined
        # object is that study's damaging-class combination.
        return dict(load("genes/HGNC_11604.json")["burden_combined"][0])
    raise AssertionError(f"unknown kind {kind}")


@pytest.mark.parametrize(("heading", "kind"), sorted(EXHAUSTIVE.items()))
def test_a_documented_example_carries_every_field_the_build_publishes(
    site: Path, heading: str, kind: str
) -> None:
    """Key sets, not values: the doc illustrates content, but promises shape.

    Values legitimately differ — the example may show a contested gene where the
    committed corpus has none — so comparing them would make the document a
    second copy of the corpus. What a consumer types off the block is the *field
    set*, and that is what has to be exhaustive.
    """
    documented = _example(heading)
    published = _published(site, kind)

    # The example may be a container (`{"genes": [ … ]}`) or the object itself.
    if len(documented) == 1:
        only = next(iter(documented.values()))
        if isinstance(only, list) and only and isinstance(only[0], dict):
            documented = only[0]

    missing = sorted(set(published) - set(documented))
    invented = sorted(set(documented) - set(published))
    assert not missing, f"{heading}: the build publishes {missing}, the doc does not show it"
    assert not invented, f"{heading}: the doc shows {invented}, the build does not publish it"


def test_the_documented_manifest_census_is_the_census_the_build_publishes(site: Path) -> None:
    """`counts` and `schema_version`, in the doc, measured against a real build.

    Values rather than key sets, unlike the shape test above, and deliberately:
    the `manifest.json` block is the one example in this document whose numbers a
    reader takes as a *census of the atlas* rather than as an illustration. It
    said `"publications": 1` for a corpus holding 4 — true when it was written,
    false for three releases after — and nothing could see it, because a wrong
    census stays internally consistent (CLAUDE.md section 4.24).

    `schema_version` is pinned here as well as in `test_build_manifest.py`. That
    file compares the published value against a literal; this one checks the
    *document* names the same version, which is the claim a consumer branches on
    after reading it. A release that bumps the constant and not the prose leaves
    the two disagreeing with nothing to notice.

    Checksums and the commit stay placeholders and are not asserted: both change
    with every commit, so pinning them would guarantee this document is wrong by
    the next one.

    **Parsed, not substring-matched, and the difference is not cosmetic.** The
    first version of this test asserted `f'"{key}": {value}' in section`, which
    is unanchored: with the build at `burden_rows: 29` and the doc still saying
    `290`, the needle `"burden_rows": 29` is a substring of `"burden_rows": 290`
    and the test passes. Any documented figure having the build's figure as a
    prefix slips through — exactly the stale-census failure this guard exists to
    catch. Found by review 2026-08-06. `_example` already parses the fenced
    block, so comparing dicts is both stronger and less code.

    The whole `counts` object is compared, not key by key, so a key the document
    invents or drops fails here too.
    """
    doc = DOC.read_text()
    start = doc.index("## `manifest.json`")
    section = doc[start : doc.index("\n## ", start)]
    documented = _example("## `manifest.json`")
    manifest = json.loads((site / "manifest.json").read_text())

    assert f'"schema_version": "{manifest["schema_version"]}"' in section
    assert documented["counts"] == manifest["counts"], (
        f"the build counts {manifest['counts']}; the documented example says {documented['counts']}"
    )

    # **Which keys count a curated file and which the build computes**, and the
    # number words in front of each list. The prose classified 9 of the 14 keys
    # until 2026-10-26: three releases had added keys to the object and none to
    # the sentence about it. The curated set is read off `write_manifest`'s own
    # `len(corpus.…)` literal, so a key moved into the wrong list fails here,
    # and a key the build adds without a sentence fails the partition.
    curated = set(re.findall(r'"(\w+)": len\(corpus\.', inspect.getsource(write_manifest)))
    computed = set(manifest["counts"]) - curated
    prose = " ".join(section.split())
    curated_clause = prose[prose.index("of its keys count a curated file") :]
    curated_clause = curated_clause[: curated_clause.index("A gene the atlas has not curated")]
    computed_block = section[section.index("keys are computed by the build") :]
    computed_block = computed_block[: computed_block.index("`index.html` states these")]
    # `(\w+)` cannot match the linked `cohorts.json`/`scope.json` file names.
    assert set(re.findall(r"`(\w+)`", curated_clause)) == curated
    assert set(re.findall(r"^  - `(\w+)`", computed_block, re.M)) == computed
    assert f"{_WORDS[len(curated)].capitalize()} of its keys count a curated file" in prose
    assert f"The other {_WORDS[len(computed)]} keys are computed by the build" in prose

    # The mirror-wider-than-the-gate figure: genes the burden mirror holds that
    # publish no page. It read "59 of its 150" for a mirror of 181 genes.
    with (REPO / "mirrors" / "burden.tsv").open() as handle:
        mirrored = {row["gene"] for row in csv.DictReader(handle, delimiter="\t")}
    published = {path.stem.replace("_", ":") for path in (site / "genes").glob("HGNC_*.json")}
    unpublished = len(mirrored - published)
    assert f"today {unpublished} of its {len(mirrored)} genes publish no page" in prose


def test_the_contested_example_is_labelled_with_the_payload_it_comes_from(site: Path) -> None:
    """`confidence_by_lesion_group` and its flag are browse-row fields only.

    The closing section shows the three contested fields in one object and tells
    a consumer to read them together. That object exists in `genes/index.json`
    and nowhere else — the gene bundle carries `has_conflicting_evidence` but
    neither per-group field — while the section above says the gene-level flag
    "appears in both the browse row and the bundle", which invites exactly the
    wrong inference.

    A detail-page author who takes the obligation to the bundle finds nothing to
    implement it with, and cannot derive it either: reconstructing per-group
    confidence means reimplementing `strongest()`'s rank and
    `has_conflicting_evidence`'s treatment of `no_known_association`, neither of
    which is published. So the document has to say where the fields live.
    """
    bundle = json.loads((site / "genes" / "HGNC_11604.json").read_text())
    # By id, for the reason `_published` gives: `genes[0]` is TBX20 since D21,
    # and the bundle this is compared against is TBX5's.
    index_row = next(
        item
        for item in json.loads((site / "genes" / "index.json").read_text())["genes"]
        if item["gene"] == "HGNC:11604"
    )

    assert "conflicting_lesion_groups" in index_row
    assert "conflicting_lesion_groups" not in bundle, "if this changes, the doc must too"
    assert "confidence_by_lesion_group" not in bundle

    text = DOC.read_text()
    section = text[text.index("## Contested genes") :]
    # Specific enough that the earlier sentence "it appears in both the browse
    # row and the bundle" — which is about the gene-level flag and is exactly
    # what misleads here — cannot satisfy it.
    assert "and **nowhere else**" in section, (
        "the section must say the per-group fields are browse-row only"
    )
    assert "conflicting_lesion_groups" in section


def _card_curation(page: str) -> str:
    """The visible text of a gene page's curation sentence, on its card (D55)."""
    rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
    found = re.findall(r'<p class="card-curation">(.*?)</p>', rail, re.S)
    assert len(found) == 1, f"expected one curation sentence on the card, found {len(found)}"
    return _text_of(found[0])


def test_the_documented_uncurated_sentences_are_the_ones_the_card_actually_shows(
    site: Path,
) -> None:
    """The doc quotes both sentences verbatim, so drift makes the doc a lie.

    Renamed 2026-10-26 from
    `test_the_documented_uncurated_notice_is_the_sentence_the_page_actually_shows`,
    when D55 moved the main-column notice onto the card as one sentence. What
    it guards is unchanged: the doc quotes what the page says, and both
    sentences are narrower than the obvious wording, deliberately -- "every
    *classification* on this page" rather than "nothing on this page", which
    would deny functional records the atlas curated for a gene with no lesion
    assertion; and "no ClinGen expert panel has graded this gene *for any of the
    diseases in scope here*", because unscoped it is false for five of the genes
    admitted on agreement. The doc was
    once written against the earlier, wider wording of the notice and kept it
    through the fix, reintroducing in prose the ambiguity the code had just
    removed. `EXHAUSTIVE` above pins documented *JSON* examples against the
    build and had nothing to say about a quoted sentence; this closes that.
    Raised by review on #14.

    Compared against real uncurated gene pages, sliced to the card's own
    curation sentence, rather than against the constants that compose it, so
    the assertion fails if the sentence stops reaching the card at all --
    comparing the doc to a constant would pass on a page that renders neither.
    """
    gata4 = json.loads((site / "genes" / "HGNC_4173.json").read_text())
    assert gata4["atlas_curation"] == "not_yet_curated", (
        "GATA4 must still be an uncurated gene for this to test"
    )
    page = (site / "genes" / "HGNC_4173.html").read_text()

    text = DOC.read_text()
    prose = re.sub(r"\s+", " ", text)
    # **Every count is read from the build, never written into the marker.**
    # The uncurated count was hard-coded as "the 22 genes" here and in the doc,
    # so the assertion pinned a stale number *in place*: the widening moved it
    # to 91 and this test went on passing, because it was comparing the doc to
    # itself.
    bundles = [
        json.loads(bundle.read_text())
        for bundle in sorted((site / "genes").glob("*.json"))
        if bundle.name != "index.json"
    ]
    uncurated = [b for b in bundles if b["atlas_curation"] == "not_yet_curated"]
    by_panel = sum(1 for b in uncurated if b["validity_state"] == "expert_curated")
    by_agreement = sum(1 for b in uncurated if b["validity_state"] == "submitter_curated")
    marker = (
        f"For the {len(uncurated)} genes published today with no curation here, the card's curation"
    )
    assert marker in text, f"the doc must state the {len(uncurated)} uncurated genes this build has"
    assert f"That wording renders on {by_panel} of those {len(uncurated)} pages" in prose
    assert f"The {by_agreement} admitted on submitter agreement" in prose

    # The two quoted blocks after the marker, each compared with a real card.
    # Joining every `>` line to the end of the file would weld them into a
    # sentence no page has ever rendered, so each is read on its own.
    after = text[text.index(marker) :]
    panel_quote = _first_block_quote(after)
    agreement_quote = _first_block_quote(after[after.index("say so instead:") :])
    assert panel_quote and agreement_quote, "the doc must still quote both sentences"
    assert panel_quote.replace("**", "") == _card_curation(page)

    # **The second sentence, on a gene no panel graded.** It is the half that
    # was wrong once: every one of the 16 submitter-admitted pages carried "the
    # classification above is an expert panel's" until 2026-08-06, describing a
    # commercial laboratory's submission as a chartered panel's classification.
    # Both that wording (guarding its return) and the panel-graded sentence are
    # asserted absent. A gene is chosen by reading `validity_state` off the
    # build rather than by naming an HGNC id, so this keeps testing the right
    # gene as the population moves.
    ungraded = next(
        bundle
        for bundle in sorted((site / "genes").glob("*.json"))
        if bundle.name != "index.json"
        and json.loads(bundle.read_text())["validity_state"] == "submitter_curated"
    )
    ungraded_html = ungraded.with_suffix(".html").read_text()
    assert agreement_quote.replace("**", "") == _card_curation(ungraded_html)
    ungraded_page = _text_of(ungraded_html)
    assert "no ClinGen expert panel has graded this gene" in ungraded_page
    assert "the classification above is an expert panel's" not in ungraded_page.lower()
    # A submitter-admitted page must not carry the panel-graded sentence either.
    assert "is a clingen expert panel's or a named curator's" not in ungraded_page.lower()


def test_the_burden_census_in_the_doc_is_the_census_the_build_publishes(site: Path) -> None:
    """Every count in the `burden` section, asserted against a real build.

    **This is the guard for the failure that actually happened.** Commit 770aa2a
    added a third study to the mirror -- 180 rows, 90 of them published, the
    first `mutation_model` comparator the atlas has ever carried -- and edited
    `docs/data-api.md` only to add the `count_unit` table. Every census figure in
    that section, and the sentence "Only `control_cohort` appears today", went on
    describing the two-study corpus. An adversarial review of the promotion diff
    measured them all false before the deploy.

    The numbers were internally consistent with each other, which is what made it
    survive reading: they were exactly the pre-Sierant subset. Only measurement
    catches that, so the doc's counts are now measured rather than read.

    The mirror figures are counted from the file rather than the build because
    the whole point of the paragraph is the gap between them -- the mirror is
    deliberately wider than the publication gate.
    """
    doc = DOC.read_text()
    start = doc.index("### The bundle's `burden` array")
    # Whitespace collapsed: these are claims in prose, and where a sentence
    # happens to wrap is not part of the claim. Asserting against the raw text
    # would make this test fail on a reflow and pass on a wrong number.
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    published = [
        row
        for path in sorted((site / "genes").glob("HGNC_*.json"))
        for row in json.loads(path.read_text()).get("burden", [])
    ]
    mirror = list(csv.DictReader((REPO / "mirrors" / "burden.tsv").open(), delimiter="\t"))

    studies = sorted({row["study"] for row in published})
    comparators = {
        c: sum(1 for r in published if r["comparator"] == c)
        for c in {r["comparator"] for r in published}
    }
    synonymous = sum(1 for r in published if r["consequence_class"] == "synonymous")
    lof = sum(1 for r in published if r["consequence_class"] == "lof")
    composite = sum(1 for r in published if r["consequence_class"] == "damaging")
    unbounded = sum(1 for r in published if r["effect_bound"] is not None)

    # Each claim, in the exact spelling the document uses.
    assert f"**{len(published)} rows reach the API**" in section
    assert f"{len(studies)} studies ({', '.join(studies)})" in section
    assert f"holds {len(mirror):,} rows for {len({r['gene'] for r in mirror})} genes" in section
    assert f"the other {len(mirror) - len(published):,}" in section

    # The sentence a consumer writes code against. It read "Only
    # `control_cohort` appears today" while 45 rows were `mutation_model`.
    assert f"`control_cohort` ({comparators['control_cohort']} rows)" in section
    assert f"`mutation_model` ({comparators['mutation_model']} rows)" in section

    assert f"{unbounded} published rows today" in section
    assert f"{synonymous} of the {len(published)} published rows are synonymous" in section
    assert f"loss-of-function is {lof} and synonymous {synonymous}" in section
    assert f"{composite} of the {len(published)} rows" in section

    # The per-study split, and the gene the widened gate admitted with no burden
    # evidence at all. Both are new in schema 2.8 and both are the kind of claim
    # section 4.24 says goes stale silently: the first is four numbers that stay
    # internally consistent while describing a corpus that no longer exists, and
    # the second is a sentence that was simply not true before the gate widened.
    per_study = ", ".join(
        str(sum(1 for r in published if r["study"] == study)) for study in studies
    )
    assert f"{per_study} rows" in section, f"the per-study split is not {per_study}"

    covered = {
        path.stem
        for path in sorted((site / "genes").glob("HGNC_*.json"))
        if json.loads(path.read_text()).get("burden")
    }
    total = len(sorted((site / "genes").glob("HGNC_*.json")))
    assert f"across {len(covered)} of the {total} published genes" in section
    assert (len(covered) < total) == ("carries no burden rows at all" in section), (
        "the doc and the build disagree about whether a published gene lacks burden rows"
    )

    # And the claim that the composite is a union rather than a partition, which
    # is what stops a consumer double-counting when it aggregates.
    assert "double-counts" in section

    # **An infinite odds ratio is a bound, not a strength.** Obligation 2 called
    # it "the strongest result in the data" until 2026-10-26, when 14 of the 29
    # published unbounded rows had a lower bound below 1. The count is pinned and
    # the retired phrase is asserted absent, for the reason section 4.35 gives.
    spans_one = sum(
        1 for r in published if r["effect_bound"] is not None and (r["ci_low"] or 0) < 1
    )
    assert f"{spans_one} of those {unbounded} rows have a `ci_low` below 1" in section
    assert "is the strongest result in the data" not in section

    # The obligations are counted, not announced: "Five obligations" sat over
    # six numbered items for a release.
    raw = doc[start : doc.index("\n## ", start)]
    block = raw[raw.index("obligations, each of which") : raw.index("`n_cases` and `n_controls`")]
    items = re.findall(r"^\d+\. \*\*", block, re.M)
    assert f"{_WORDS[len(items)].capitalize()} obligations, each of which" in section

    # Obligation 3's per-study rule, against the mirror it describes: no
    # PMID:42230622 row is empty on both sides, and the one row that is comes
    # from a study the doc names, for a gene the doc names.
    empty = [r for r in mirror if r["n_case_carriers"] == "0" and r["n_control_carriers"] == "0"]
    rows_42230622 = sum(1 for r in mirror if r["study"] == "PMID:42230622")
    assert not [r for r in empty if r["study"] == "PMID:42230622"]
    assert f"zero of its {rows_42230622:,} mirror rows have no case carrier" in section
    with (REPO / "mirrors" / "genes.tsv").open() as handle:
        symbols = {g["hgnc_id"]: g["symbol"] for g in csv.DictReader(handle, delimiter="\t")}
    assert len(empty) == 1, f"{len(empty)} rows are empty on both sides; the doc names one"
    assert f"`{empty[0]['study']}` does publish a row with zeros on both sides" in section
    # The row named whole: its gene, class and design (review, 2026-10-26).
    named = {"lof": "loss-of-function"}[empty[0]["consequence_class"]]
    design = {"control_cohort": "case-control"}[empty[0]["comparator"]]
    assert f"{symbols[empty[0]['gene']]}'s {named} {design} row" in section

    # PMID:34324492's denominators, which the doc names as the filtered set its
    # article's methods put the test on. Every row of that study carries the
    # same two, so the doc can state them once. Until 2026-10-26 the rows
    # carried S1 Table's deletion column instead and the doc called the choice
    # open; a mirror re-run from a converter reverted to S1 fails here.
    cnv = [r for r in mirror if r["study"] == "PMID:34324492"]
    [n_cases] = {int(r["n_cases"]) for r in cnv}
    [n_controls] = {int(r["n_controls"]) for r in cnv}
    assert (
        f"For `PMID:34324492` they are the {n_cases:,} cases and {n_controls:,} controls that "
        "the article's methods say remained after filtering" in section
    )


def test_the_independent_datasets_section_states_the_shape_the_build_publishes(
    site: Path,
) -> None:
    """The census guard, applied to the new section before it can go stale.

    The burden section drifted by a whole study because its counts were written
    once and never re-measured, and the numbers stayed internally consistent so
    reading them found nothing. This section is pinned from the day it lands
    rather than after the same thing happens to it.

    The three display obligations are asserted as *text*, not as numbers,
    because they are what stops the object being rendered as a verdict — and a
    later edit that softens them should fail here rather than ship.
    """
    doc = DOC.read_text()
    start = doc.index("### `independent_datasets`")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    published = [
        json.loads(path.read_text())["independent_datasets"]
        for path in sorted((site / "genes").glob("HGNC_*.json"))
    ]
    states = {entry["state"] for gene in published for entry in gene["families"]}

    # Every state the build can emit is documented; a fifth would fail here.
    assert states <= {"corrected", "nominal", "no_enrichment", "not_tested"}
    for state in states:
        assert f"| `{state}` |" in section, f"{state!r} is emitted and undocumented"

    # `tested` is the denominator, and it is never `len(families)` where any
    # family did not test the gene. Measured on the real corpus rather than
    # asserted, so the doc's claim is checked against what ships.
    assert any(gene["tested"] < len(gene["families"]) for gene in published), (
        "no published gene has an untested family, so the denominator claim is untested"
    )

    assert "never be rendered as one" in section
    assert "never `len(families)`" in section
    assert "KDM6A causes Kabuki syndrome" in section

    # **How many families, read off the build** (review, 2026-10-27). The
    # section said "three today" and "three `not_tested` entries" for a release
    # after I3 made them two: nothing here read the count, only the states, the
    # denominator and the obligations, so a sentence stating a number passed
    # beside a build publishing another.
    [count] = {len(gene["families"]) for gene in published}
    assert f"one entry per cohort family in the whole corpus** — {_WORDS[count]} today" in section
    assert f"A gene no study reported gets {_WORDS[count]} `not_tested` entries" in section

    # `members` (2.19): every study of the family, in its order, and the family's
    # state the strongest of theirs -- the claim the section makes, on what ships.
    rank = ["not_tested", "no_enrichment", "nominal", "corrected"]
    for gene in published:
        for family in gene["families"]:
            assert [member["study"] for member in family["members"]] == family["studies"]
            assert family["state"] == max(
                (member["state"] for member in family["members"]), key=rank.index
            )
    assert "Attribute a state to a study from `members`, never from the family" in section


def test_the_omics_section_documents_profiles_stratified_top_not_ranked_by_significance() -> None:
    """16a: `top` was documented as "ranked by significance" with no
    qualification, which Task 8 made false for `profiles` -- that table has
    no significance column at all, so its slice is stratified: the cardiac
    series leads on derived percentile, with slots reserved for the other
    tissues so a gene's own τ stays auditable from the same payload that
    publishes it. A consumer reading the old, unscoped sentence expects the 25
    highest-ranked rows for every modality and gets something deliberately
    different for this one.

    **The reservation is bounded and the doc must not overstate it.** An
    earlier wording promised "at least one row of every other tissue present
    in the shard", which `select_top` guarantees only while there are 12 or
    fewer such tissues -- `reserved` is `min(len(others), TOP_N // 2)`.
    Measured 2026-08-19 with a 19-row cardiac series: 6 and 12 other tissues
    all appear, 13 and 20 do not. True of this dataset, which has six, and
    stated as a general guarantee it is not. Raised by review on PR #39, and
    the negative assertion below is what stops the promise coming back.
    """
    doc = DOC.read_text()
    start = doc.index("## `omics/<modality>/<accession>.json`")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    # Scoped to the three modalities it is actually true of, not stated as a
    # blanket rule every modality's `top` follows.
    assert "For `expression`, `proteomics` and `phospho`" in section
    assert "ranked by significance" in section
    assert "profiles` ranks `top` differently" in section
    assert "stratified" in section
    assert "median_percentile" in section
    assert "reserved for the non-cardiac tissues" in section
    # The bound must be stated, not just the reservation.
    assert "12" in section
    # Negative: the unconditional promise must not return. A positive
    # assertion fails when a sentence is deleted; this fails when a deleted
    # one comes back, which is what nobody is watching for.
    assert "at least one row of every other tissue" not in section
    assert "reserved" in section
    assert "D39(b)" in section
    assert "ranking bug" in section


def test_the_datasets_field_list_names_every_field_the_model_declares() -> None:
    """16b: the field list for `datasets.json` must name every field `Dataset`
    publishes, so a field added to the model and not to this sentence fails a
    test instead of aging into a false sentence the way Task 3's six fields
    already did once, silently, before this guard existed (CLAUDE.md section
    4.24 -- true when written, false by the next release, and internally
    consistent throughout, so re-reading it finds nothing).

    Read from `Dataset.model_fields`, never from a real corpus record. That
    was originally because no dataset had ever been committed, so a real
    build's `datasets.json` was always `{"datasets": []}` and there was no
    record on disk to read -- **which stopped being true on 2026-08-19, when
    `E-MTAB-6814` landed**, and the sentence is corrected rather than deleted
    because the reason to keep reading the model is now a different and better
    one. `build/literature.py::_dump` calls `model_dump(mode="json")` with no
    field filtering, so the model's declared fields are exactly what any
    record publishes, present or future. Reading a committed record instead
    would only pin the fields *that* record happens to exercise, and would
    catch a field being removed from the doc while missing one being added to
    the model -- which is the direction this guard exists for.

    `Stage.model_fields` as well as `Dataset`'s, since schema 2.12. `stages`
    is the one field on this record that is a list of structured entries
    rather than a scalar, so it is the one whose sub-fields a consumer types
    against without the outer field list saying anything about them. `order`
    arrived there and not on `Dataset`, so the `Dataset` loop alone could not
    have seen it.
    """
    doc = DOC.read_text()
    start = doc.index("## `datasets.json`")
    section = doc[start : doc.index("\n## ", start)]

    for field in Dataset.model_fields:
        assert f"`{field}`" in section, (
            f"Dataset.{field} is not named in the datasets.json field list"
        )
    for field in Stage.model_fields:
        assert f"`{field}`" in section, (
            f"Stage.{field} is not named in the datasets.json stage description"
        )


def test_the_expression_profile_layer_is_documented_end_to_end() -> None:
    """16c: the two mirror tables, the `profile_quantiles` payload and the
    `expression_profile` bundle key were entirely undocumented before this.
    Checks the load-bearing claims a consumer needs, sharing one fixture
    because they are all facts about one new section:

    - the method travels with the number -- `LOOKUP_RULE`, `TAU_SCALE` and
      `TAU_METHOD` quoted from the constants that actually ship, not retyped,
      so a wording change to either constant must be echoed here or this
      fails;
    - why the quantile grid is published at all -- a one-level trade, stated
      as one, never as a proof;
    - percentiles are not comparable across organs, and -- more importantly,
      because the layer's headline question rides on it -- not comparable
      across developmental stages;
    - τ carries no adjective and no band (D39(c));
    - τ is blind to where a gene peaks. Verified directly against
      `profiles.specificity` before writing the doc, not reasoned from a
      number handed down: `specificity({"heart": 20.0, "liver": 200.0, "o1":
      5.0, "o2": 5.0, "o3": 5.0, "o4": 5.0, "o5": 5.0}, floor=0.0)` returns
      `tau == 0.6227723406454512` and `highest_in == "liver"`; the same inputs
      score `0.9625` computed by hand on the raw values with no log transform.
      Both round to the figures asserted below.
    """
    doc = DOC.read_text()

    assert "## `omics/profile_quantiles/<accession>.json`" in doc

    start = doc.index("### The bundle's `expression_profile` object")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    assert LOOKUP_RULE in section
    assert TAU_SCALE in section
    assert TAU_METHOD in section

    assert "one-level guarantee" in section
    assert "not an unbounded one" in section
    assert "re-derivable" in section

    assert "not comparable across organs" in section
    assert "not comparable across developmental stages" in section
    assert "matters more" in section

    assert "no adjective, no band" in section

    assert "0.963" in section
    assert "0.623" in section
    assert "liver" in section
    assert "highest_in" in section


def test_the_api_doc_states_the_expression_census_the_build_produces(site: Path) -> None:
    """Section 4.24: when a release adds a data source, every count in every
    prose artifact is suspect by default. Grep the numbers and re-derive them
    from a real build; do not re-read the sentences. Eight false claims
    survived multiple readings last time because they still agreed with each
    other.

    The figures moved on 2026-08-19, when E-MTAB-6814 landed: 0 genes across 0
    datasets became 92 across 1. That is exactly the event this test exists for
    -- the doc's sentence was internally consistent and false the moment a
    source arrived, which is how eight claims survived last time.

    Matched on the singular stem "dataset", not "datasets": the doc must read
    "across 1 dataset", and a pin that demanded "1 datasets" would force
    ungrammatical English into a public artifact to keep a test green.

    **Pinned to the sentence since 2026-10-25, not to the figure.** It asserted
    `f"{n} genes" in doc`, page-wide (4.19), so once the Orphanet widening took
    `profile_genes` to 119 it passed on "98 of the 119 genes" in the validity
    section while this sentence still said 92.
    """
    manifest = json.loads((site / "manifest.json").read_text())
    doc = " ".join(DOC.read_text().split())
    genes, datasets = manifest["counts"]["profile_genes"], manifest["counts"]["profile_datasets"]
    sentence = f"{genes} genes carry a developmental expression profile, across {datasets} dataset"

    assert sentence in doc


def test_the_doc_states_the_stage_order_datasets_json_actually_publishes(site: Path) -> None:
    """CLAUDE.md 4.24, caught in the artifact the section itself exists to guard.

    Through schema `2.11` this section told a consumer to **"Sort on `order`
    rather than trusting this array's own sequence"**, and gave the reason:
    `datasets.json` is serialised generically from the curated record, so
    `stages` arrived in declaration order and "nothing enforces it" to match
    `order`. Both halves were true when written. `build/literature.py` now
    sorts the array at emit, so the second half became false while the advice
    stayed harmless -- which is the shape that survives re-reading, because a
    consumer who follows it still gets the right answer and never learns the
    sentence is wrong.

    **No published byte distinguishes the two**, and that is what makes this a
    documentation guard rather than a build one: the curated file was already
    in `order`, so `datasets.json` is byte-identical before and after (measured
    2026-08-21, `cmp` over two real builds). `test_built_site_is_consumable.py
    ::test_datasets_json_orders_its_stages_by_order_not_by_declaration` pins
    the emitter's behaviour; nothing pinned the sentence a consumer reads about
    it, so reverting this paragraph passed the whole suite.

    Three assertions on one section, scoped to it rather than to the file
    (4.19 -- `order` and the word "sorted" both occur elsewhere in this
    document, so a document-wide check would pass on a section that says the
    opposite):

    - the retired instruction is **absent**. A negative assertion, for 4.35's
      reason: the positive one below fails when the sentence is deleted, this
      one fails when the deleted sentence comes back, and coming back is what
      nobody is watching for.
    - the section states the sort, against a **literal** rather than against
      any constant the code builds it from (4.38).
    - the claim is **true of a real build**, so the wording cannot be corrected
      into a different false sentence.
    """
    doc = DOC.read_text()
    start = doc.index("## `datasets.json`")
    section = doc[start : doc.index("\n## ", start)]

    assert "rather than trusting this array's own sequence" not in section
    assert "a curated record printed as written" not in section
    assert "This array is sorted by the build, on `order`" in section

    published = json.loads((site / "datasets.json").read_text())
    for record in published["datasets"]:
        orders = [stage["order"] for stage in record["stages"]]
        assert orders == sorted(orders), (
            f"{record['id']} publishes its stages out of `order`, which is what "
            f"the datasets.json section of the API doc tells a consumer it does not"
        )


def test_the_combined_census_in_the_doc_is_the_census_the_build_publishes(site: Path) -> None:
    """Every count in the `burden_combined` section, against a real build (§4.24).

    The published figure and the mirrored one differ on the committed corpus
    (108 and 123, measured 2026-10-26), so this distinguishes the two (§4.30): a
    section that quoted the mirror as what reaches the API would fail here.
    """
    doc = DOC.read_text()
    start = doc.index("### The bundle's `burden_combined` array")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n### ", start + 1)])
    bundles = [json.loads(p.read_text()) for p in sorted((site / "genes").glob("HGNC_*.json"))]
    published = [row for bundle in bundles for row in bundle["burden_combined"]]
    genes = sum(1 for bundle in bundles if bundle["burden_combined"])
    with (REPO / "mirrors" / "burden_combined.tsv").open() as handle:
        mirror = list(csv.DictReader(handle, delimiter="\t"))

    assert f"**{len(published)} combined results reach the API**" in section
    assert f"across {genes} published genes" in section
    assert f"holds {len(mirror)} rows for {len({r['gene'] for r in mirror})} genes" in section
    assert "is not a corrected p" in section

    # The example abbreviates the row's note with "…", so what it quotes must
    # still be the start of the note the build publishes -- batch D reworded
    # that note's Fisher sentence on 2026-10-26. And the prose beside it must not
    # keep the every-gene claim the note dropped: the meta p reproduces as
    # Fisher's method within 10% for 116, 133 and 144 of the 248 genes.
    quoted = _example("### The bundle's `burden_combined` array")["method_note"]
    note = _published(site, "combined_row")["method_note"]
    assert note.startswith(quoted.removesuffix("…").rstrip(" ,"))
    assert "It does not reproduce as Fisher's method" not in section
    assert "For many of the study's genes it does not reproduce as Fisher's method" in section


def test_the_scope_census_in_the_doc_is_the_registry_the_build_publishes(site: Path) -> None:
    """`scope.json`'s counts, and the Orphanet warrant's shape, against a real build."""
    doc = DOC.read_text()
    start = doc.index("## `scope.json`")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start + 1)])
    terms = json.loads((site / "scope.json").read_text())["terms"]
    orphanet = [term for term in terms if term["orphanet"] is not None]

    assert f"**{len(terms)} terms**" in section
    assert f"{len(terms) - len(orphanet)} are attributed to an authority that curates" in section
    assert f"**{len(orphanet)}** are admitted on Orphanet's phenotype frequencies" in section
    assert [term["id"] for term in terms] == sorted(term["id"] for term in terms)

    # The term-level vocabulary, which shipped in 2.15 with no table at all
    # (2026-10-26): every token the build emits has a row, and the doc names no
    # token the build cannot emit. The name it shares with the bundle's
    # `admitted_by` -- an object, with different values -- is called out.
    tokens = {term["admitted_by"] for term in terms}
    # The header row is a row too, and names the column rather than a token.
    rows = set(re.findall(r"\| `(\w+)` \| [^|]+ \| [^|]+ \|", section)) - {"admitted_by"}
    assert rows == tokens, f"the build emits {sorted(tokens)}; the table lists {sorted(rows)}"
    assert "is not a gene bundle's `admitted_by`" in section
    # `omim` on a direct term: the doc says null, every time.
    direct = [term for term in orphanet if term["orphanet"]["via"] == "direct"]
    assert direct and all(term["orphanet"]["omim"] is None for term in direct)
    assert '`omim` is `null` on every `"direct"` term' in section

    # **The attribution `scope.json` owes.** It republishes Orphanet's codes,
    # names and frequency bands, and `sources.json` records those as
    # attribution-required; the redistribution note covered HPO, ClinGen and
    # GenCC only until 2026-10-26. Every source the paragraph names must exist
    # and carry the terms the paragraph says it does.
    sources = {
        source["id"]: source
        for source in json.loads((site / "sources.json").read_text())["sources"]
    }
    note = _section("**`scope.json` republishes Orphanet content")
    named = re.findall(r"`(orphanet_\w+|mondo_\w+)`", note)
    assert {"orphanet_product1", "orphanet_product4", "mondo_sssom"} <= set(named)
    for source in named:
        assert sources[source]["licence"] == "CC-BY-4.0"
        assert sources[source]["redistribution"] == "permitted_with_attribution"


def test_the_doc_states_the_population_the_build_publishes(site: Path) -> None:
    """Four present-tense sentences that state how much the atlas publishes.

    Section 4.24, found 2026-10-25: the Orphanet widening took the build from 92
    genes and 94 pages to 119 and 121, the census refresh that went with it
    re-derived every pinned figure, and these four still said 92, 93 and 94
    because nothing pinned them. Pages and genes differ, so neither assertion
    can pass on the other's number (4.30).
    """
    doc = " ".join(DOC.read_text().split())
    genes = json.loads((site / "manifest.json").read_text())["counts"]["genes"]
    pages = len(list(site.rglob("*.html")))

    assert f"the entry point to the other {pages - 1}." in doc
    assert f"behind this document publishes {pages} of them" in doc
    assert f"never the {genes} this atlas publishes" in doc
    assert f"this atlas publishes {genes} of them" in doc

    # **Sub-populations, which went stale in the same widening and were caught
    # by review rather than by this test** (2026-10-26). Each pair below is two
    # figures that differ on the committed corpus, so neither can pass on the
    # other's number (4.30): uncurated genes with a panel grade and without
    # one; genes admitted by ClinGen and on agreement.
    bundles = _bundles(site)
    uncurated = [b for b in bundles if b["atlas_curation"] == "not_yet_curated"]
    graded = sum(1 for b in uncurated if b["headline_confidence"] is not None)
    ungraded = len(uncurated) - graded
    assert graded != ungraded
    # "presents 118 of the 119 ... as this atlas's assessment" sat beside "on
    # 21 it is null": only a grade can be mistaken for an assessment.
    assert f"presents {graded} of the {genes} published genes" in doc
    assert f"On the other {ungraded} genes not yet curated here it is `null`" in doc
    # The empty-map paragraph: the same figure twice, and it read 97 and 75.
    rows = json.loads((site / "genes" / "index.json").read_text())["genes"]
    empty = sum(
        1
        for r in rows
        if r["headline_confidence"] is not None and not r["confidence_by_lesion_group"]
    )
    assert f"{empty} of the {genes} published rows carry a non-null `headline_confidence`" in doc
    assert "understates 75" not in doc
    authority = [b["admitted_by"]["authority"] for b in bundles]
    listed = _section("**Which genes are listed.**", end="**No disease is in scope")
    assert f"for a disease in scope here ({authority.count('clingen')} genes)" in listed
    assert f"no ClinGen panel disputes it ({authority.count('gencc_agreement')} genes)" in listed
    # Withdrawn symbols in search: genes in the index with a `prev_symbols` entry.
    with (REPO / "mirrors" / "genes.tsv").open() as handle:
        previous = {
            g["hgnc_id"] for g in csv.DictReader(handle, delimiter="\t") if g["prev_symbols"]
        }
    indexed = {
        record["id"]
        for record in json.loads(gzip.decompress((site / "search" / "index.json.gz").read_bytes()))[
            "records"
        ]
        if record["kind"] == "gene"
    }
    assert f"{len(indexed & previous)} of the {len(indexed)} indexed genes have one" in doc


def test_the_report_url_census_in_the_doc_is_the_census_the_build_publishes(site: Path) -> None:
    """The `report_url` bullet's figures, measured against a real build (section 4.24).

    The sentence went stale once already -- it read "157 of the 266" for a build
    publishing 470 GenCC records -- and nothing could see it, because stale prose
    stays internally consistent. Each figure is anchored on the words either side
    of it, so a documented "189" cannot pass for a built "89" (the prefix trap
    `test_the_documented_manifest_census_is_the_census_the_build_publishes` records).
    """
    bundles = sorted((site / "genes").glob("HGNC_*.json"))
    gencc = [
        record
        for bundle in bundles
        for record in json.loads(bundle.read_text())["validity"]["records"]
        if record["source"] == "gencc"
    ]
    homepages = ("https://panelapp-aus.org", "http://panelapp-aus.org")
    null = sum(record["report_url"] is None for record in gencc)
    homepage = sum((record["report_url"] or "").rstrip("/") in homepages for record in gencc)
    doc = " ".join(DOC.read_text().split())

    assert f"Of the {len(gencc)} GenCC records across the {len(bundles)} gene bundles," in doc
    assert f", {null} publish `report_url: null`" in doc
    assert f"and {homepage} more point at" in doc
    assert f"so {null + homepage} of the {len(gencc)} offer" in doc


def test_the_cohorts_example_quotes_the_record_the_build_publishes(site: Path) -> None:
    """The `cohorts.json` example quotes `taa_cases`, and it went stale (2026-10-26).

    Its description was rewritten in `curation/cohorts.yaml` -- it published a
    curator's instruction ("and the page must say so") and a denominator no
    row uses -- and an example quoting the old text would teach a consumer the
    caveat that was withdrawn. The example abbreviates with a trailing "…", so
    what it quotes must be a prefix of the published description; the name
    and url are quoted whole.
    """
    documented = _example("## `cohorts.json`")["cohorts"][0]
    published = {
        cohort["id"]: cohort
        for cohort in json.loads((site / "cohorts.json").read_text())["cohorts"]
    }[documented["id"]]

    assert set(documented) == set(published)
    assert documented["name"] == published["name"]
    assert documented["url"] == published["url"]
    quoted = documented["description"].removesuffix("…").rstrip()
    assert len(quoted) > 100, "the example quotes too little of the description to check"
    assert published["description"].startswith(quoted)


@pytest.mark.parametrize(
    ("heading", "payload", "key", "pattern"),
    [
        ("## `cohorts.json`", "cohorts.json", "cohorts", r"`(pmid\w+)`"),
        ("## `scope.json`", "scope.json", "terms", r"(MONDO:\d+)"),
    ],
)
def test_the_review_status_bullets_name_the_records_the_build_marks(
    site: Path, heading: str, payload: str, key: str, pattern: str
) -> None:
    """Each section's `review_status` bullet names its records with an open question.

    Names, not a count, are what a reader acts on, and a count written beside
    them ("Two records carry one today") expires the day a curator marks a
    third. Read from the bullet itself, so a record marked or cleared in the
    curation fails here until the bullet is reworded. Measured on the build it
    was written against: two cohorts and two scope terms, each `pending_review`.

    **The bullet ends where its paragraph does**, at the next bullet or a
    blank line, found in the raw text and asserted found. It was sliced from
    the whitespace-collapsed section up to the next `" - "`, which `find`
    misses after a section's last bullet: `cohorts.json`'s ran to the end of
    the section, and a name moved out of the bullet into the paragraph after
    it passed the full suite (review, 2026-10-27, mutant M26).
    """
    doc = DOC.read_text()
    section = doc[doc.index(heading) : doc.index("\n## ", doc.index(heading) + len(heading))]
    start = section.index("- `review_status` and `review_note`")
    ends = [end for end in (section.find("\n- ", start), section.find("\n\n", start)) if end != -1]
    assert ends, f"the review_status bullet under {heading} has no end inside its section"
    bullet = " ".join(section[start : min(ends)].split())
    marked = sorted(
        record["id"]
        for record in json.loads((site / payload).read_text())[key]
        if record["review_status"] is not None
    )
    assert marked, f"no {payload} record carries a review status; this checks nothing"
    assert sorted(set(re.findall(pattern, bullet))) == marked
    assert f"{_WORDS[len(marked)].capitalize()} " in bullet
    # Review, 2026-10-27 (D1): the cohorts bullet defined `pending_review` as
    # "the sources disagree", false for the two scope terms it also marks, and
    # the words are the label of the `has_source_discordance` chip.
    assert "sources disagree" not in bullet.lower()


def test_the_contested_section_names_each_contested_gene_and_what_its_headline_is(
    site: Path, mirrored: dict[str, GeneValidity]
) -> None:
    """The closing section's census and its account of `headline_confidence`.

    **Both were false on 2026-10-26, and internally consistent.** The section
    said each of the six contested genes "is `Definitive` for Noonan syndrome
    (or a Noonan subtype)"; BRAF and MAP2K1 are admitted on, and Definitive
    for, cardiofaciocutaneous syndrome, and ClinGen grades them Moderate and
    Limited for Noonan. And it opened by defining `headline_confidence` as "the
    strongest classification the mirrored ClinGen and GenCC records assert",
    which measured false for 31 of the 119 published genes -- TBX1's submitters
    grade it Definitive and its headline is `null` -- and contradicted the
    browse-row section's own definition.

    So the table is derived here from the build rather than read: one row per
    contested bundle, its diseases at the headline grade (which must include
    the one `admitted_by` names) and its disputed ones, and no gene contested
    for a disease it is graded on -- the claim the section's heading makes. The
    mirror-wide figures come from the `mirrored` fixture, before the gate.
    """
    section = _section("## Contested genes")
    raw = DOC.read_text()[DOC.read_text().index("## Contested genes") :]
    bundles = _bundles(site)
    published = {b["gene"] for b in bundles}
    contested = [b for b in bundles if b["has_conflicting_evidence"]]
    assert contested, "no published gene is contested; the table measures nothing"
    assert {b["headline_confidence"] for b in contested} == {"definitive"}, (
        "the table's column is headed `Definitive`"
    )

    expected: set[str] = set()
    for bundle in contested:
        clingen = [r for r in bundle["validity"]["records"] if r["source"] == "clingen"]
        headline = bundle["headline_confidence"]
        graded = sorted({r["disease_label"] for r in clingen if r["classification"] == headline})
        against = [r for r in clingen if r["classification"] in ("disputed", "refuted")]
        disputed = sorted({r["disease_label"] for r in against})
        assert bundle["admitted_by"]["disease_label"] in graded
        assert not set(graded) & set(disputed), f"{bundle['symbol']} is contested where graded"
        for panel in {r["gcep"] for r in against}:
            assert panel in section, f"{panel} disputes {bundle['symbol']}, unnamed"
        expected.add(f"| {bundle['symbol']} | {'; '.join(graded)} | {'; '.join(disputed)} |")
    table = raw[raw.index("| gene | `Definitive` for |") :]
    table = table[: table.index("\n\n")]
    assert set(table.splitlines()[2:]) == expected
    assert "Each is `Definitive` for Noonan syndrome" not in section

    # The mirror-wide census beside the table.
    flagged = {
        gene
        for gene, entry in mirrored.items()
        if has_conflicting_evidence(r.classification for r in entry.records if r.classification)
    }
    discordant = {gene for gene, entry in mirrored.items() if entry.has_source_discordance}
    assert flagged == discordant, "the section says the same genes set both flags"
    assert (
        f"of the {len(mirrored)} mirrored genes, {len(flagged)} set `has_conflicting_evidence`, "
        f"and the same {len(discordant)} set `has_source_discordance`. "
        f"{len(flagged & published)} are published." in section
    )
    assert f"the {len(flagged)} genes that set one flag set both" in section
    for gene in sorted(flagged - published):
        assert f"({gene})" in section, f"{gene} is contested and unpublished; the doc names it not"

    # What `headline_confidence` is: ClinGen's grade or `null`, never a max over
    # GenCC submitters. Each clause the section states, measured.
    agreement = {b["gene"] for b in bundles if b["admitted_by"]["authority"] == "gencc_agreement"}
    assert {b["gene"] for b in bundles if b["headline_confidence"] is None} == agreement
    assert "which today is exactly the genes admitted on submitter agreement" in section
    tbx1 = json.loads((site / "genes" / "HGNC_11592.json").read_text())
    terms = {r["classification_term"] for r in tbx1["validity"]["records"]}
    assert tbx1["headline_confidence"] is None and {"Definitive", "Strong"} <= terms
    assert "TBX1's submitters grade it `Definitive` and `Strong`, and its headline is" in section
    assert "the strongest classification the mirrored ClinGen and GenCC records" not in section

    # `admitted_by.classification` and the headline: equal on every gene, by
    # measurement and not "by construction" -- they are two computations, and
    # would part on an agreement gene carrying a ClinGen NKA record (none is
    # published: every agreement gene has no ClinGen record in scope at all).
    assert all(b["headline_confidence"] == b["admitted_by"]["classification"] for b in bundles)
    assert {b["validity_state"] for b in bundles if b["gene"] in agreement} == {"submitter_curated"}
    provenance = _section("### `admitted_by` and `asserted_by`", end="\n### ")
    assert "by construction rather than by coincidence" not in provenance
    assert "and no such gene is published" in provenance

    # The genes whose only in-scope record is an ungraded `Supportive`, named by
    # id in the browse-row section. It listed NOTCH1 (HGNC:7881) for a release
    # after NOTCH1 was published on its own ClinGen term, and missed CBL.
    supportive = sorted(
        gene
        for gene, entry in mirrored.items()
        if all(record.classification is None for record in entry.records)
    )
    assert all(
        len(mirrored[gene].records) == 1 and mirrored[gene].records[0].submitter == "Orphanet"
        for gene in supportive
    )
    assert not set(supportive) & published
    doc = " ".join(DOC.read_text().split())
    word = _WORDS[len(supportive)]
    assert f"{word.capitalize()} genes in the committed mirrors ({', '.join(supportive)})" in doc
    assert f"{word} genes in the committed mirrors carry exactly one such Orphanet record" in doc


def test_the_phase_paragraphs_state_each_end_basis_and_the_phases_no_stage_reaches(
    site: Path,
) -> None:
    """`heart_looping` was documented as excluded at every `wpc` for two months
    after it stopped being excluded at all (found 2026-10-26).

    On 2026-08-24 its end became `derived` -- GO:0001947 defines looping as
    ending "preceding septation", and the source dates septation's start at
    CS12 -- and `phases_for` matches any phase with an end. The doc went on
    saying the source "never states when it ends", that `end_basis` is
    `not_stated`, and that the phase is excluded from every match. What is true
    is narrower and is a fact about the data: no published stage is early
    enough to fall inside it, nor inside `embryonic_heart_tube_morphogenesis`.

    So every clause is measured: the one `derived` phase and its boundary from
    the curation file, the absence of `not_stated` phases, and the unreached
    phases from the built bundles against the earliest published stage.
    """
    corpus, _ = load_curation(REPO)
    assert corpus.cardiac_phases is not None
    phases = corpus.cardiac_phases.phases
    section = _section("### The bundle's `expression_profile` object")
    six = _WORDS[len(phases)]
    assert f"curates {six} phases" in section

    derived = [phase for phase in phases if phase.end_basis.value == "derived"]
    assert len(derived) == 1, "the doc names one `derived` phase"
    phase = derived[0]
    assert f"`{phase.id}` is the one `derived` phase" in section
    assert (
        f"the source states when looping starts ({phase.start_carnegie_stage}) and never "
        "when it ends" in section
    )
    assert f"{phase.end_carnegie_stage}, {phase.end_wpc:g} wpc (`{phase.end_hsapdv_id}`)" in section
    assert not [p for p in phases if p.end_basis.value == "not_stated"]
    assert f"None of the {six} is `not_stated` today" in section

    reached = {
        phase_id
        for bundle in _bundles(site)
        for dataset in bundle["expression_profile"]["datasets"]
        for stage in dataset["stages"]
        for phase_id in (stage["phase"] or {}).get("phase_ids", [])
    }
    earliest = min(
        stage["wpc"]
        for dataset in json.loads((site / "datasets.json").read_text())["datasets"]
        for stage in dataset["stages"]
        if stage["wpc"] is not None
    )
    unreached = [p for p in phases if p.id not in reached]
    assert unreached == [p for p in phases if p.end_wpc is not None and p.end_wpc <= earliest]
    paragraph = section[section.index("never appear in a published `phase_ids` today") :]
    assert f"**{_WORDS[len(unreached)].capitalize()} of the {six} never appear" in section
    assert f"the earliest stage any published dataset samples is {earliest:g} wpc" in paragraph
    for p in unreached:
        assert f"`{p.id}`" in paragraph and f"{p.end_wpc:g}" in paragraph
    assert "never appears in `phase_ids` at all" not in section


def _abridged(block: str) -> Any:
    """A JSON example abridged with `…`, parsed, each `…` kept as a marker.

    An object's trailing `, … }` becomes a `"…"` key and any other `…` an array
    item `"…"`, so `_agrees` can tell an abridged value from a complete one.
    """
    keyed = re.sub(r",\s*…\s*\}", ', "…": "…" }', block)
    return json.loads(re.sub(r'(?<!")…(?!")', '"…"', keyed))


def _agrees(documented: Any, published: Any, at: str = "") -> None:
    """Every value an example shows is the published one, and only `…` omits.

    An object marked `…` may leave keys out and an array ending in `…` may
    leave trailing items out; an unmarked object or array must be complete.
    Leaves compare with their type as well as their value, because `0 ==
    False` (CLAUDE.md section 5).
    """
    if isinstance(documented, dict):
        assert isinstance(published, dict), at
        shown = {key: value for key, value in documented.items() if key != "…"}
        if "…" not in documented:
            assert set(shown) == set(published), f"{at}: keys differ and nothing is marked …"
        for key, value in shown.items():
            assert key in published, f"{at}.{key} is not published"
            _agrees(value, published[key], f"{at}.{key}")
    elif isinstance(documented, list):
        assert isinstance(published, list), at
        abridged = bool(documented) and documented[-1] == "…"
        shown_items = documented[:-1] if abridged else documented
        if abridged:
            assert len(shown_items) < len(published), f"{at}: marked … and omits nothing"
        else:
            assert len(shown_items) == len(published), (
                f"{at}: {len(published)} items published, {len(shown_items)} shown, no …"
            )
        for index, (item, actual) in enumerate(zip(shown_items, published, strict=False)):
            _agrees(item, actual, f"{at}[{index}]")
    else:
        assert (type(documented), documented) == (type(published), published), (
            f"{at}: the doc shows {documented!r}, the build publishes {published!r}"
        )


def test_the_examples_that_claim_to_be_real_are_the_ones_the_build_publishes(site: Path) -> None:
    """The opening line says every example is a real build's unless it says
    otherwise. Measured 2026-10-26, four said nothing and were not: an
    `expression` omics block for a dataset that does not exist, a populated
    variant index where the real one is empty, a search record missing a term
    the build publishes (TBX5's `HOS`), and a contested browse row no
    published gene has. The bundle example also dropped three of its 23 keys.

    Each real example is compared with the published object, value by value
    -- through `_agrees` where it is abridged with `…` -- and each
    illustrative one is asserted to be labelled so and to show a shape the
    build really has no instance of.

    **The abridged bundle and `omics` examples were compared by key set alone
    until review on 2026-10-26** (plus `admitted_by`'s values), so
    `asserted_by`'s `count` and `omics.profiles.count` could be edited to any
    number with every test green, and the `burden` array showed a
    PMID:42230622 row as its only item where TBX5's first of 14 is
    PMID:40127276's. `_agrees` reads each `…` as the one mark of an omission:
    every value shown must be published and every unmarked array and object
    complete. Two more unmarked abridgements in the bundle example were marked
    with it: `burden_combined` showed one of TBX5's 3 objects, and the
    assertion object 3 of its 13 keys.
    """
    doc = DOC.read_text()
    tbx5 = json.loads((site / "genes" / "HGNC_11604.json").read_text())

    def block(heading: str) -> str:
        text = doc[doc.index(heading) :]
        opening = text.index("```json\n") + len("```json\n")
        return text[opening : text.index("\n```", opening)]

    _agrees(_abridged(block("## `genes/<slug>.json`")), tbx5, "bundle")
    _agrees(_abridged(block("## `omics/<modality>/<accession>.json`")), {"omics": tbx5["omics"]})

    variants = json.loads((site / "variants" / "index.json").read_text())
    assert _example("## `variants/index.json`") == variants
    assert variants == {"shards": []}, "a variant is curated; the example need not be illustrative"

    records = {
        record["id"]: record
        for record in json.loads(gzip.decompress((site / "search" / "index.json.gz").read_bytes()))[
            "records"
        ]
    }
    for record in _example("## `search/index.json.gz`")["records"]:
        assert record == records[record["id"]], f"the search example misquotes {record['id']}"

    rows = json.loads((site / "genes" / "index.json").read_text())["genes"]
    assert not [r for r in rows if r["has_conflicting_evidence"] and r["lesion_groups"]]
    assert "*Illustrative* — no published gene today is both contested and curated" in (
        _section("## Contested genes")
    )
    assert "*Illustrative* — no proteomics or phospho mirror is committed" in _section(
        "## `omics/<modality>/<accession>.json`"
    )
    assert not list((REPO / "mirrors").glob("proteomics*")) and not list(
        (REPO / "mirrors").glob("phospho*")
    )


def _card_paragraph(page: str, css: str) -> str:
    """The visible text of one `<p class="...">` on a gene page's card, or ''."""
    rail = page[page.index('<aside class="rail">') : page.index("</aside>")]
    found = re.findall(rf'<p class="{css}">(.*?)</p>', rail, re.S)
    return _text_of(found[0]) if found else ""


def test_the_page_sentences_the_doc_reports_are_the_ones_the_pages_render(site: Path) -> None:
    """Sentences the doc quotes or restates from a page, against the page.

    Earlier batches of the 2026-10-26 review reworded several pages, and the
    doc kept the old wording of three: it quoted the card's Orphanet line as
    "Graded for Noonan syndrome as a syndrome" where every page says "as a
    whole syndrome"; it stated the scope rule as a sufficient condition
    ("for a disease that an external authority treats as congenital heart
    disease, or in which Orphanet records ...") after `methods.html#scope`
    had narrowed it to the terms the atlas lists; and it said the gene
    index's cells "are plain text, not chips" beside a `conflicting evidence`
    chip on every contested row. Each is pinned to a literal and to the page,
    and each retired wording is asserted absent from the doc. So are the gene
    page's short forms and the labels of its links to `methods.html`, which
    the doc quotes so a reader can find them on the page.
    """
    doc = " ".join(DOC.read_text().split())
    pages = {path.stem: path.read_text() for path in sorted((site / "genes").glob("HGNC_*.html"))}

    quoted = re.search(r'the reason: "(Graded for [^"]+)"', doc)
    assert quoted, "the doc must still quote the card's Orphanet line"
    # Owner decision R2: the doc quotes the card's S8 sentence; GATA4's card has it.
    s8 = (
        "For this gene, the study's S8 Table counts every carrier among its CHD cases, so "
        "none is an aneurysm proband."
    )
    assert f'shows: "{s8}"' in doc
    assert s8 in _card_paragraph(pages["HGNC_4173"], "card-note")
    scope_lines = {_card_paragraph(page, "card-scope") for page in pages.values()}
    assert quoted.group(1) in scope_lines
    assert "Graded for Noonan syndrome as a syndrome" not in doc

    rule = (
        "A disease is in scope when it is one of the terms this atlas lists, and a term is "
        "listed only if an external authority treats it as congenital heart disease, or if a "
        "curator nominates a syndrome and Orphanet records a structural heart malformation in "
        "at least 30% of its cases (for a subtype, of the syndrome it belongs to)"
    )
    assert rule in _text_of((site / "methods.html").read_text())
    assert rule in doc
    assert "or in which Orphanet records a structural heart malformation" not in doc

    index_html = (site / "genes" / "index.html").read_text()
    contested = {
        r["gene"]
        for r in json.loads((site / "genes" / "index.json").read_text())["genes"]
        if r["has_conflicting_evidence"]
    }
    assert contested
    chip = '<span class="chip chip-warn">conflicting evidence</span>'
    for search, grade, cells in re.findall(
        r'<tr data-search="([^"]+)"[^>]*data-confidence="([^"]*)"[^>]*>(.*?)</tr>', index_html, re.S
    ):
        gene = "HGNC:" + search.split()[0].removeprefix("hgnc:")
        confidence = re.findall(r"<td>(.*?)</td>", cells, re.S)[2]
        assert (chip in confidence) == (gene in contested), gene
        assert " " not in grade and "chip" not in grade
    assert "Cells are plain text, not chips" not in doc
    assert "a contested gene's confidence cell adds a `conflicting evidence` chip" in doc

    # Rendered strings the doc cites as examples of what a page shows.
    every_page = "".join(pages.values())
    for literal in (
        "OR ∞ (CI 0.786–∞)",
        "OR ∞ (95% CI 28.1–∞)",
        "About these cohorts, and the caveats recorded for them",
        "row is a union, and this study published no non-syndromic · loss-of-function row "
        "for this gene.",
        # The matrix label's mark for a family's study with no row for the
        # gene, and EVC2's card scope line, whose panel graded EVC rather than
        # EVC2 (both review, 2026-10-27).
        "(no row for this gene)",
        "in scope on ClinGen's Kidney Cystic and Ciliopathy Disorders Gene Curation Expert "
        "Panel's curation of it for EVC.",
        # The overlap notice the doc says a declared overlap renders (review,
        # 2026-10-27: the doc named only a shared collection as its trigger).
        "collections this atlas records as overlapping, so their results may describe partly "
        "the same people and must not be combined",
    ):
        assert literal in doc.replace("*", ""), literal
        assert literal in html.unescape(every_page), f"no page renders {literal!r}"
    # The short forms and toggles the doc names (2026-10-26), against the text a
    # reader sees: tags stripped with no space, because two carry a `<strong>`.
    # Each link label is pinned with the anchor the doc says it opens.
    shown = html.unescape(re.sub(r"<[^>]+>", "", every_page))
    for sentence in (
        "Every count and statistic in the study tables below is an upstream study's, "
        "republished as published; the atlas computes none of them.",
        "The atlas computes no pooled statistic across studies: pooling two that share a "
        "sample collection would count the same people twice, and pooling two that draw on "
        "collections this atlas records as overlapping could.",
        "A cell names the most significant of the rows behind its state; the card lists a "
        "study's surviving rows with a damaging union before its parts, so the two can lead "
        "with different rows.",
        "Declaration: this study is by an author of this atlas; its results are republished "
        "as published and read by the same rule as every other study's.",
        "Illumina reports no known disease relationship for congenital heart defects, "
        "multiple types, 6 (2021); see the validity table.",
        "Every measurement below is from bulk tissue: a low or absent figure is not "
        "evidence a gene is unimportant to heart development.",
    ):
        assert f'"{sentence}"' in doc, sentence
        assert sentence in shown, f"no page renders {sentence!r}"
    for anchor, label in (
        ("", "How to read this page"),
        ("#rare-variants", "How to read these results"),
        ("#rare-variants", "How to read rare-variant results"),
        ("#expression", "How to read these figures"),
        ("#expression", "How to read expression figures"),
        ("#curated-evidence", "How to read curated evidence"),
    ):
        assert f"*{label}*" in doc, label
        assert f'<a href="../methods.html{anchor}">{label}</a>' in every_page, label
    assert "*How to read this chart*" in doc
    assert "<summary>How to read this chart</summary>" in every_page
    methods = _text_of((site / "methods.html").read_text())
    assert "a flag on the gene, not on a disease" in methods
    assert "Both are flags on the gene, not on a disease" in doc
    # An unsplit study's `all cases` is not its whole case set: PMID:34324492's
    # rows count 6,746 of its 7,958 cases. The doc restated `#rare-variants`'
    # retired "whole case set" until review on 2026-10-26.
    assert "all cases is its only case group: every case the row counts, not a subgroup" in (
        methods
    )
    assert "only case group, not a subgroup, where it does not" in doc
    assert "whole case set" not in doc


def test_the_doc_documents_the_cell_type_expression_bundle_key(site: Path) -> None:
    """The section exists, its example's keys are the real shape, its claims are there.

    A real build has no cell-type mirror yet, so every bundle publishes the
    empty shape and the counts are 0; the section must not claim otherwise.
    """
    bundle = json.loads(next((site / "genes").glob("HGNC_*.json")).read_text())
    assert bundle["cell_type_expression"] == {"datasets": []}
    counts = json.loads((site / "manifest.json").read_text())["counts"]
    assert counts["cell_type_genes"] == counts["cell_type_datasets"] == 0

    doc = DOC.read_text()
    heading = (
        "### The bundle's `cell_type_expression` object: a cell-type profile, never a contrast"
    )
    assert heading in doc
    section = doc[doc.index(heading) + len(heading) :]
    section = section[: section.index("\n## ")]
    prose = " ".join(re.sub(r"```.*?```", "", section, flags=re.S).split())
    assert len(prose.split()) <= 150
    assert "not evidence against a gene" in prose
    assert "9 to 15 weeks" in prose
    assert "no heart had at least 10 cells" in prose
    example = json.loads(re.search(r"```json\n(.*?)```", section, re.S).group(1))  # type: ignore[union-attr]
    dataset = example["datasets"][0]
    assert set(dataset) == {
        "dataset", "unit", "detection_floor", "n_donors", "specificity", "cell_types",
    }  # fmt: skip
    assert "`2.20`" in " ".join(doc.split())
