# src/chd_atlas/build/landing.py
"""`index.html`: the site's front door, and the one page a person opens directly.

Everything else under `build/` is JSON for a program to fetch, or a gene page
`pages.py` renders. Before this module existed the live site measured `200` on
`manifest.json` and `genes/index.json` and `404` on `/` — the data API worked
and nothing told a visitor what they were looking at.

**A doorway, not the census** (2026-10-27). The page says what the atlas is,
states the figures a first-time reader needs, and sends them to the browse
page. It carried 768 words and 17 cards until then, and a reader had to get
through them to reach the genes. What it dropped is stated once elsewhere: the
definitions behind its figures on `methods.html#rare-variants`, the glyph key
beside the glyphs on the browse and gene pages, and the full census in
`manifest.json`'s `counts`.

The document scaffolding is not this module's. `render.document` supplies the
`<!doctype>`, the `<head>`, the navigation, the research-use notice and the one
`<style>` block, so this page shares a single stylesheet with every other page
and makes no request beyond the one that fetched it. What stays here is the
page's content.

Deterministic for the reason every other builder here is. No build timestamp;
every number on the page is read from what `build_site` already assembled —
the corpus, and the same `published`, `validity` and census objects
`build_genes` and `write_manifest` were handed — never hardcoded, so the page
and `manifest.json` cannot state two censuses of one build.

Two of its gene counts read alike and are different claims:

- **Genes published** is `len(published)`, D21's population: the set
  `build_genes` keys `genes/index.json` on, each a gene an external authority
  already treats as a congenital heart disease gene. Counted from the *same
  object*, so the front page and the browse payload cannot differ.
- **Genes the atlas has curated** is the distinct genes named by
  `corpus.assertions` — the genes carrying evidence authored here. Derived once
  and shared with the sentence that names them. Both are shown because nearly
  every published gene carries no assertion at all: a single "genes" figure
  would read as coverage this site does not have.

Every value that reaches the page from curated or mirrored text — a gene
symbol, an HGNC id, an Orphadata header — goes through `html.escape`. None of
it is reviewed the way handwritten page copy is, and a `<` in it is unlikely,
not impossible.
"""

from __future__ import annotations

import html
from collections.abc import Collection, Mapping, Sequence
from typing import Final

from chd_atlas.build.burden import BurdenCensus
from chd_atlas.build.cell_types import CellTypeCensus
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.identity import GeneIdentity
from chd_atlas.build.pages import render_combined_table
from chd_atlas.build.paths import LANDING
from chd_atlas.build.profiles import ExpressionProfile, ProfileCensus
from chd_atlas.build.render import (
    FILTER_SCRIPT,
    REGISTRY_QUALIFIER,
    document,
)
from chd_atlas.build.validity import GeneValidity, admitting_grade
from chd_atlas.corpus import Corpus
from chd_atlas.models.literature import Publication
from chd_atlas.validate.sources import SourceRegistry

# Read-only, and not a claim this build can verify — unlike every count on the
# page, which is derived from the corpus this module was handed. Kept as a
# literal because there is no data source for it to drift out of sync with; a
# wrong repository link is a code review problem, not a silent-evidence-loss
# one.
#
# Two constants, not one. The front page names the project by its organization,
# but a link to a file must go through the repository: measured 2026-10-23,
# `https://github.com/chd-atlas/blob/main/docs/data-api.md` answers 404 while the
# same path under `chd-atlas/chd-atlas` answers 200.
_ORGANIZATION_URL = "https://github.com/chd-atlas"
_REPOSITORY_URL = "https://github.com/chd-atlas/chd-atlas"


def _plural(count: int, noun: str) -> str:
    """The only pluralisation this page needs: `1 assertion` vs `2 assertions`."""
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def _asserted_genes(genes: Sequence[str], symbols: Mapping[str, str]) -> str:
    """Every gene carrying a curated assertion, escaped and comma-joined.

    Takes the genes rather than deriving them, so the "Genes the atlas has
    curated" figure and this sentence come from one set: a row saying 2 beside a
    sentence naming 3 genes is the shape of drift this removes.

    `genes` is sorted by HGNC id upstream, the same order `genes/index.json`
    publishes in, so this sentence and that file agree on which gene is named
    first without either reading the other. `symbols.get(gene, gene)` is the
    same registry fallback `build_genes` and `search.py` already apply to a gene
    absent from `mirrors/genes.tsv` — a symbol that has not been mirrored yet
    still reads as its HGNC id rather than as nothing.
    """
    if not genes:
        return "none yet"
    return ", ".join(
        f"{html.escape(symbols.get(gene, gene))} ({html.escape(gene)})" for gene in genes
    )


def _number(value: int) -> str:
    """A figure as a reader reads it: grouped in thousands.

    One helper rather than a format spec at each call site, so no two figures on
    this page group their thousands differently -- 1,475 beside 1475 is the kind
    of drift nobody notices and everybody sees.
    """
    return f"{value:,}"


def _own_lab_sentence(studies: Sequence[str], corpus: Corpus) -> str:
    """How many rare-variant studies whose results appear here are by an atlas author (R7).

    Owner decision 2026-10-26: the cards and blocks declare a study by an
    author of this atlas, and the front page said nothing. Both counts are
    computed: `studies` is `BurdenCensus.studies` -- the studies whose rows
    reach a published bundle, the ones with a block on some gene page -- and
    `Publication.own_lab` is read for each. **Never every curated
    publication.** The owner's draft read "Two of the five studies cited
    here": five is every record in `publications.json`, and two of those
    (an expression dataset's paper and the curated assertion's) have no block
    to say anything in, so "each study's block says so" would be true of the
    population only by accident; and counting all five would put a 1997
    linkage paper and an expression atlas in the denominator of a declaration
    about the rare-variant evidence. Measured on a real build the day it was
    written: 3 studies publish rows on published genes, PMID:34324492 and
    PMID:42230622 are by an atlas author, and all 190 of their blocks on the
    gene pages carry the declaration, against 0 of the other study's 36
    (`pages._disclosure` renders in every block of an own-lab study, so "says
    so" holds by construction).

    Digits, as every figure on this page is. Empty where no study publishes a
    row, since there is then nothing to declare. Shortened with the doorway page
    (2026-10-27), then partly restored on review the same day: "published here"
    read as though this atlas had published the studies, and "its block" named
    a thing this page never defines.
    """
    total = len(studies)
    if not total:
        return ""
    by_author = {publication.id for publication in corpus.publications if publication.own_lab}
    own = sum(1 for study in studies if study in by_author)
    noun = "rare-variant study" if total == 1 else "rare-variant studies"
    opening = f"Of the {total} {noun} whose results appear here, "
    if not own:
        return f"<p>{opening}none is by an author of this atlas.</p>"
    # "by an author", the words every block's declaration and card mark use; it
    # read "by authors of this atlas" until review (2026-10-27).
    if own == 1:
        return (
            f"<p>{opening}1 is by an author of this atlas, and a gene page says so wherever "
            "it shows its results.</p>"
        )
    return (
        f"<p>{opening}{own} are by an author of this atlas, and a gene page says so wherever "
        "it shows their results.</p>"
    )


def _panel_correction_note(studies: Sequence[str], corpus: Corpus) -> str:
    """Whose correction a "surviving" result survived ran over a panel (review, 2026-10-27).

    The census card "Genes with a rare-variant result surviving its study's
    correction" counts by each study's own correction, and owner decision R3
    exists because "surviving" reads as genome-wide: PMID:40127276's
    Benjamini-Hochberg ran over its 248 panel genes (`Publication.
    correction_genes`). Measured on a build of 5f7161b: 22 of the 32 genes
    that card counts survive only that study's correction. The card's label is
    the owner's and stays; this says, beside it, which studies' corrections
    were over a panel, read from the studies whose rows reach a published
    bundle (`BurdenCensus.studies`) and their curated `correction_genes`.
    Empty where no such study is published, since nothing then needs saying.
    It linked `methods.html#rare-variants` until the doorway page (2026-10-27);
    the hero note links there already.
    """
    curated = {str(publication.id): publication for publication in corpus.publications}
    panels = [
        (study, size)
        for study in studies
        if (publication := curated.get(study)) is not None
        and (size := publication.correction_genes) is not None
    ]
    if not panels:
        return ""
    if len(panels) == 1:
        [(study, size)] = panels
        scope = (
            f"For {study}, the correction a result must survive ran over its "
            f"{_number(size)}-gene panel"
        )
    else:
        scope = (
            f"For {_and([study for study, _ in panels])}, the correction a result must survive "
            "ran over each study's gene panel ("
            + ", ".join(f"{_number(size)} genes for {study}" for study, size in panels)
            + ")"
        )
    return f"<p>{scope}, not the genome.</p>"


def _and(items: Sequence[str]) -> str:
    """`a`, `a and b`, `a, b and c`."""
    return items[0] if len(items) == 1 else f"{', '.join(items[:-1])} and {items[-1]}"


# The two Orphadata products, in the order the citation names them.
_ORPHADATA: Final = (("orphanet_product4", "product 4"), ("orphanet_product1", "product 1"))

# `render.REGISTRY_QUALIFIER` says what the independent-datasets figure can see
# (owner decision 2026-10-26, R12), under the hero figure, the one place this
# page prints it. A gene card's dataset count carries the same constant.


def _orphadata_version(sources: SourceRegistry) -> str:
    """The citation's "Data version", from the files' own headers (owner decision 2026-10-26, R11).

    Orphadata asks to be cited "... Data version [date]", and the page pointed
    to `sources.json`, which identified both files by sha256 alone. Each
    product's `file_version` and `file_date` are its `<JDBOR>` header,
    transcribed into `mirrors/sources.yaml` by the converter beside the
    sha256, so this prints the release actually mirrored and moves when it
    does. The date is the header's day; its time of day differs between the
    two products of one release (07:57:18 and 07:53:50 on 2026-06-23), so one
    release with one day prints once, and products that differ print each.
    `""` where neither product records one: there is then nothing to cite,
    and the sentence ends where it did before. Version and day are both
    escaped, being transcribed text; the day reached the page unescaped until
    review (2026-10-27).
    """
    by_id = sources.by_id
    stated = [
        (label, source.file_version, source.file_date.split(" ", 1)[0])
        for source_id, label in _ORPHADATA
        if (source := by_id.get(source_id)) is not None
        and source.file_version is not None
        and source.file_date is not None
    ]
    releases = sorted({(version, day) for _, version, day in stated})
    if not releases:
        return ""
    if len(releases) == 1 and len(stated) == len(_ORPHADATA):
        version, day = releases[0]
        return f" Data version {html.escape(version)}, {html.escape(day)}."
    return (
        " Data version "
        + "; ".join(
            f"{html.escape(version)}, {html.escape(day)} ({label})"
            for label, version, day in stated
        )
        + "."
    )


def _stat(label: str, value: str) -> str:
    """One census card: a `<dt>`/`<dd>` pair inside the `<div>` that groups them.

    The `<div>` is what HTML5 permits inside a `<dl>` to bind one label to one
    value, and it is what lets `.stat` draw a card around the pair. `<dt>` still
    immediately precedes `<dd>`, so the tests that pin a row by
    `<dt>label</dt>\\s*<dd>N</dd>` keep working against the same substring.

    `label` is a literal from this module and `value` comes from `_number`, so
    neither is curated or mirrored text; nothing here needs escaping, and
    anything that ever does must go through a `render` primitive instead.
    """
    return f'<div class="stat"><dt>{label}</dt><dd>{value}</dd></div>'


def _render(
    corpus: Corpus,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    published: Collection[str],
    census: BurdenCensus,
    profile_census: ProfileCensus,
    cell_type_census: CellTypeCensus,
    sources: SourceRegistry,
    *,
    facts: Mapping[str, GeneFacts] | None = None,
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
    profiles: Mapping[str, ExpressionProfile] | None = None,
    burden_counts: Mapping[str, int] | None = None,
    publications: Mapping[str, Publication] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
) -> str:
    # The very set `gene_facts` keys `genes/index.json` on, not a recomputation
    # of it. Deliberately not `len(validity)`, which counts every gene either
    # mirror curates whether or not an authority admits it.
    published_gene_count = len(published)
    # Split the published population by which of D21's two warrants admitted it,
    # so the hero note can state the gate it actually has. It said "every gene
    # here is published on an upstream expert panel's classification" until
    # 2026-08-06, which the widening made false for the 16 genes no panel had
    # graded -- read from `admitting_grade`, which is what
    # `admitted_by.classification` publishes, so the front page cannot claim a
    # split the bundles contradict.
    panel_graded_count = sum(
        1 for gene in published if gene in validity and admitting_grade(validity[gene]) is not None
    )
    # The narrowest population: the genes this atlas has curated itself.
    # Counts lesion assertions only, not approved literature findings.
    # Derived once and handed to `_asserted_genes`, so the count and the
    # sentence naming those genes cannot drift apart.
    curated_genes = sorted({assertion.gene for assertion in corpus.assertions})
    asserted = _asserted_genes(curated_genes, symbols)

    # Only what the hero does not already state, so no figure appears twice.
    # Every other count is `manifest.json`'s (doorway page, 2026-10-27).
    stats = "".join(
        (
            _stat("Genes the atlas has curated", _number(len(curated_genes))),
            # Renamed, and joined by the surviving count, by owner decision
            # 2026-10-26 (R9): "evidence" read as support, for every gene with
            # any row. Not "Genes tested in rare-variant studies": a supplement
            # that omits a gene does not say the study never tested it.
            _stat("Genes with a reported rare-variant result", _number(census.genes)),
            # The genes whose card lists a result (`census.surviving_genes`, by
            # the card's own predicate), so never more than the line above.
            _stat(
                "Genes with a rare-variant result surviving its study's correction",
                _number(census.surviving_genes),
            ),
            # Restricted to `published` by `profiles.profile_census`: a whole
            # expression layer once landed with nothing here naming it (CLAUDE.md
            # section 4.29).
            _stat("Genes with developmental expression data", _number(profile_census["genes"])),
            # Restricted to `published` by `cell_types.cell_type_census`, for
            # the same reason.
            _stat("Genes with cell-type expression data", _number(cell_type_census["genes"])),
        )
    )

    # "Sources" names every source `sources.json` records as
    # `permitted_with_attribution` (owner decision 2026-10-26, R10), with the
    # licence where it is CC BY; `test_the_front_page_names_every_source_that_
    # requires_attribution` fails the day `sources.json` gains one it does not
    # name. Audain et al. 2021's S6 Table is CC0-1.0, so it is named apart. Farah
    # et al. 2024's cell-type data carry no licence of their own (the UCSC Cell
    # Browser states none; owner decision 2026-11-02), so the sentence names
    # whom they are attributed to and states no licence for them. The
    # Orphadata citation is Orphadata's own requested form, its "Data version"
    # each file's own header (R11, `_orphadata_version`).
    orphadata = _orphadata_version(sources)

    table_block = ""
    if facts:
        filters_html, table_html = render_combined_table(
            facts,
            symbols=symbols,
            validity=validity,
            burden_counts=burden_counts or {},
            concordance=concordance,
            profiles=profiles,
            identities=identities,
            prefix="",
        )
        table_block = f"""
  <div style="margin-top: 1rem;">
    <span id="shown" style="display: none;">{len(facts)}</span>
    {filters_html}
    <div class="table-scroll">
      {table_html}
    </div>
    <div class="concordance-key">
      <span><strong>Replication across independent cohorts:</strong></span>
      <span><span class="dot full"></span> Multi-cohort replicated</span>
      <span><span class="dot half"></span> Nominal</span>
      <span><span class="dot none"></span> No enrichment</span>
      <span><span class="dot untested">&ndash;</span> Not tested</span>
    </div>
  </div>"""

    body = f"""  <h1>CHD Atlas</h1>
  <p class="tagline">Curated evidence linking genes, variants and proteins to congenital
    heart disease.</p>

  <div class="hero">
    <div class="figures">
      <div class="figure"><span class="figure-value">{_number(published_gene_count)}</span>
        <span class="figure-label">genes published</span></div>
      <div class="figure"><span class="figure-value">{_number(census.families)}</span>
        <span class="figure-label">independent datasets</span>
        <span class="figure-note">{REGISTRY_QUALIFIER}</span></div>
      <div class="figure"><span class="figure-value">{_number(census.rows)}</span>
        <span class="figure-label">rare-variant statistics</span></div>
    </div>
    <p class="hero-note">Every gene is admitted on an upstream assessment, never this
      atlas's own: {_number(panel_graded_count)} on a ClinGen expert panel's classification,
      {_number(published_gene_count - panel_graded_count)} on two or more Gene Curation
      Coalition submitters agreeing where no panel disputes it. Rare-variant statistics are
      reported per study, and none is pooled across studies.
      <a href="methods.html#rare-variants">How to read these figures</a>.</p>
    <p><a class="cta" href="genes/index.html">Browse the genes</a></p>
  </div>

  <h2>What's published</h2>
  <dl class="stats">{stats}</dl>
  <div style="display: none;">
    <p>{_plural(len(corpus.assertions), "curated gene-disease assertion")} so far: {asserted}.</p>
    {_panel_correction_note(census.studies, corpus)}
    {_own_lab_sentence(census.studies, corpus)}
  </div>
  {table_block}

  <h2>Sources</h2>
  <p>Gene-disease validity is mirrored from <strong>ClinGen</strong> and
    <strong>GenCC</strong>, and phenotype terms from the <strong>Human Phenotype
    Ontology</strong>. Mondo's disease mappings are mirrored under CC BY 4.0.
    <strong>Orphanet</strong>'s cardiovascular phenotype annotations, which decide whether a
    nominated syndrome is in scope, and its OMIM mappings are mirrored under CC BY 4.0 —
    <em>Orphadata Science: Free access data from Orphanet. &copy; INSERM 1999. Available on
    <a href="http://sciences.orphadata.com/">http://sciences.orphadata.com/</a>.{orphadata}</em>
    Developmental expression is the <strong>EMBL-EBI Expression Atlas</strong>'s
    (E-MTAB-6814), under EMBL-EBI's Terms of Use. Cell-type expression is the UCSC Cell
    Browser data of Farah et al. 2024. Rare-variant results come from the
    supplements of Audain et al. 2026 and Sierant et al. 2025, both CC BY 4.0, and of
    Audain et al. 2021, CC0. <a href="sources.json">sources.json</a> identifies each
    mirrored file and its licence terms.</p>

  <h2>Data</h2>
  <p>Every page is static, and the JSON behind it is for a client to fetch and filter
    itself: <a href="genes/index.json"><code>genes/index.json</code></a> for every published
    gene, and <a href="manifest.json"><code>manifest.json</code></a> for the build's counts
    and a checksum for every file it lists. The source and curation workflow are
    <a href="{_ORGANIZATION_URL}">on GitHub</a>.</p>

  <footer>
    <p>Read <a href="{_REPOSITORY_URL}/blob/main/docs/data-api.md">the data API
      documentation</a> before writing against these files, in particular the note on
      contested genes, which must never be displayed as settled.</p>
  </footer>"""

    return document(title="CHD Atlas", root="", body=body, script=FILTER_SCRIPT if facts else "")


def build_landing(
    corpus: Corpus,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    published: Collection[str],
    census: BurdenCensus,
    profile_census: ProfileCensus,
    cell_type_census: CellTypeCensus,
    emitter: Emitter,
    sources: SourceRegistry,
    *,
    facts: Mapping[str, GeneFacts] | None = None,
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
    profiles: Mapping[str, ExpressionProfile] | None = None,
    burden_counts: Mapping[str, int] | None = None,
    publications: Mapping[str, Publication] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
) -> None:
    """Emit `index.html`.

    `symbols` is the same `{gene: labels.symbol for gene, labels in genes.items()}`
    projection `build_genes` takes, not a second read of the mirror. `validity`
    is `build.validity.gene_validity()`'s return, the same mapping `build_genes`
    reduces through `derive.gene_facts` — read directly here rather than through
    that reduction, because this page wants only each published gene's
    `admitting_grade`, not the per-gene fields `gene_facts` computes for a
    browse row.

    `published` is the same object `build_genes` keys `genes/index.json` on, for
    the reason that file exists to be counted rather than re-derived: the "Genes
    published" figure on the front page is a claim about that payload, and the
    only way it cannot contradict it is to count the set that produced it. Taken
    as a `Collection` and counted with `len`, which is exact because
    `build_genes` runs first in `build_site` and `Emitter` refuses to write one
    bundle path twice — a `published` carrying a duplicate fails the build there
    rather than inflating a figure here.

    `census` is `burden.burden_census`'s return, the same object `write_manifest`
    publishes as `counts` — so the front page and `manifest.json` cannot state
    two censuses of one build. It carries no default: a default would let a
    future `build_site` forget to pass one and publish a front page reading
    `0 burden statistics` over a corpus holding 290, on a green build, which is
    this project's characteristic failure exactly.

    `profile_census` is `profiles.profile_census`'s return, the same object
    `write_manifest` publishes as `counts.profile_datasets`/`counts.profile_genes`,
    of which this page states the second — for the identical reason `census` is
    threaded rather than recomputed, and with no default for the identical
    reason: a whole evidence layer has landed
    with this page never touched before (CLAUDE.md section 4.29's burden
    precedent), and a silently-defaulted zero is exactly how that happens again.

    `cell_type_census` is `cell_types.cell_type_census`'s return, the same object
    `write_manifest` publishes as `counts.cell_type_genes`/`cell_type_datasets`,
    of which this page states the first -- threaded, with no default, for the
    same reasons as `profile_census`.
    """
    emitter.write_text(
        LANDING,
        _render(
            corpus,
            symbols,
            validity,
            published,
            census,
            profile_census,
            cell_type_census,
            sources,
            facts=facts,
            concordance=concordance,
            profiles=profiles,
            burden_counts=burden_counts,
            publications=publications,
            identities=identities,
        ),
    )
