# src/chd_atlas/build/runner.py
"""Orchestration: gate on validation, then emit every artifact in order.

This is the only module that decides *whether* to build. Everything under
`build/` publishes what it is given; the judgement that the corpus is fit to
publish is made once, here, and made before the first byte is written.

The gate is `validate_repository` itself rather than a rule of this module's
own. A stricter rule here would make a repository that CI calls green
unpublishable; a looser one would publish a corpus the validator rejected. Both
are drift, and the only way to not drift is to ask the same question the same
way — hence `report.ok`, which ignores warnings by design, rather than a
condition assembled here.

Order matters in exactly one place: the manifest is last, because
`emitter.checksums` is complete only once every other builder has run. That is
no longer a convention to remember — `write_manifest` seals the emitter, so a
builder added below it fails the build rather than publishing a file the
manifest does not list.
"""

from __future__ import annotations

from pathlib import Path

from chd_atlas.build.bundles import build_genes
from chd_atlas.build.burden import (
    build_cohorts,
    burden_census,
    cohort_registry,
    combined_census,
    load_burden,
    load_combined,
)
from chd_atlas.build.cell_types import cell_type_census, gene_cell_type_expression
from chd_atlas.build.concordance import (
    cohort_families,
    evidence_axes,
    gene_concordance,
    survives_own_correction,
)
from chd_atlas.build.emit import Emitter
from chd_atlas.build.gene_literature import gene_literature, literature_census
from chd_atlas.build.identity import read_gene_identities
from chd_atlas.build.landing import build_landing
from chd_atlas.build.literature import build_literature, build_sources
from chd_atlas.build.manifest import source_commit, write_manifest
from chd_atlas.build.methods import build_methods
from chd_atlas.build.omics import build_omics
from chd_atlas.build.pages import build_gene_index_page, build_gene_pages
from chd_atlas.build.profiles import (
    build_profile_quantiles,
    gene_expression_profiles,
    percentile_annotations,
    profile_census,
)
from chd_atlas.build.scope import build_scope, scope_warrants
from chd_atlas.build.search import GeneLabels, build_search
from chd_atlas.build.validity import gene_validity, published_genes
from chd_atlas.build.variants import build_variants
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.orphanet_scope import load_orphanet_tables, pinned_structural_phenotypes
from chd_atlas.tables import TABLE_SCHEMAS, cell_text, read_table
from chd_atlas.validate.runner import validate_repository
from chd_atlas.validate.sources import load_sources


class BuildRefused(Exception):
    """The repository did not pass validation, so nothing was written.

    Publishing from a corpus that fails its own gate would defeat the point of
    having one: the site would carry claims the validator rejected, and nobody
    would find out until a reader did. The message carries the rendered report,
    because a curator who ran `build` needs to know what to fix without being
    told to go and run a second command.
    """


def _gene_registry(root: Path) -> dict[str, GeneLabels]:
    """What `mirrors/genes.tsv` says about each gene, keyed on HGNC id.

    One `GeneLabels` per row rather than parallel mappings, so a gene cannot end
    up with a symbol and no name. `build_genes` needs only the symbol and is
    handed a projection of this rather than a second read of the mirror.

    The aliases cell is split here, and this is the only place that knows the
    separator. The mirror stores them pipe-separated because the schema keeps
    every TSV one value per cell rather than inventing a nested format, so
    unpacking that is the reader's job. Handing the raw cell to `GeneLabels`
    would type-check — `frame.to_dicts()` yields `dict[str, Any]`, and an `Any`
    satisfies `tuple[str, ...]` silently — and would publish a search index of
    single characters; the guard in `GeneLabels.__post_init__` is what catches
    that, and it exists because this line is where the mistake is made.

    Returning empty for a missing or unreadable mirror is unreachable behind
    `build_site`, which refuses on any error and so never reaches here without
    TBL008 having fired. It is not a fallback anyone should rely on: it is what
    keeps this function total for a direct caller, in the same way
    `search.py`'s symbol fallback is.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns:
        return {}

    registry: dict[str, GeneLabels] = {}
    for row in frame.to_dicts():
        gene = cell_text(row.get("hgnc_id"))
        # Both `hgnc_id` and `symbol` are non-nullable in the schema, so a blank
        # here is a bypassed TBL003. The two blanks are not the same failure,
        # though. A row with no id has nothing to key the registry on and no
        # row to correct it against, so it is dropped — there is no gene here to
        # attach `name` or `aliases` to. A row with no symbol is still a real,
        # validated gene: `name` and `aliases` are perfectly good and reached
        # this reader having passed every check that applies to them, so
        # dropping the row along with the blank symbol would discard both for
        # no reason a curator chose. The symbol falls back to the HGNC id
        # instead — the same fallback `build_genes` and `search.py` already
        # apply to a gene missing from this registry entirely — so the row, its
        # name and its aliases all still reach `genes/index.json` and the
        # search index, rather than the gene silently losing everything but its
        # id. Found by adversarial review after release: the previous version
        # of this guard dropped the row and, with it, a name and aliases that
        # had already been validated and were fit to publish.
        if gene is None:
            continue
        registry[gene] = GeneLabels(
            symbol=cell_text(row.get("symbol")) or gene,
            name=cell_text(row.get("name")),
            aliases=_pipes(row.get("aliases")),
            prev_symbols=_pipes(row.get("prev_symbols")),
        )
    return registry


def _pipes(raw: object) -> tuple[str, ...]:
    """One pipe-separated mirror cell as its terms, blanks dropped.

    Shared by `aliases` and `prev_symbols` rather than written twice. The two
    columns have identical shape and identical failure mode -- handing the raw
    `str` to `GeneLabels` publishes one search term per *character* -- and
    `GeneLabels.__post_init__` refuses that for both. A second inline copy is
    how the two drift, and the one that drifts is the one nobody was looking at.

    Returns `()` for a null, an empty cell, or a cell of only separators, so a
    caller cannot tell "no aliases" from "the column is absent". Neither should
    reach a search term, and `TBL001`/`TBL003` are what report the difference.
    """
    if not raw:
        return ()
    return tuple(part.strip() for part in str(raw).split("|") if part.strip())


def _assertions_by_gene(corpus: Corpus) -> dict[str, list[LesionAssertion]]:
    """Every curated assertion, grouped by the gene it is about, ordered by id.

    `pages.py` renders one `<h3>` block per assertion in the order it is handed
    them, so this order is part of a page's bytes and therefore its checksum.
    Ordered by id rather than left in corpus order for the reason
    `bundles._records_by_gene` gives: corpus order is `load_curation`'s — record
    files sorted by filename, then each file's own order within it — which makes
    it an artefact of how the curation is filed. Without this, moving an
    assertion from `TBX5.yaml` into a second file would reorder a page that is
    otherwise unchanged.

    Not shared with `bundles._records_by_gene`, which groups the same records
    the same way and then dumps each one to JSON: these builders need the models
    themselves, and a helper returning both shapes would be doing two jobs.

    Genes outside the published population are grouped here too and never looked
    up — `build_gene_pages` iterates `facts`, so what it renders is decided by
    the population, not by this mapping.
    """
    grouped: dict[str, list[LesionAssertion]] = {}
    for assertion in corpus.assertions:
        grouped.setdefault(assertion.gene, []).append(assertion)
    return {gene: sorted(items, key=lambda record: record.id) for gene, items in grouped.items()}


def build_site(root: Path, out: Path) -> dict[str, str]:
    """Build the published API into `out`. Returns path-to-checksum per file.

    Raises `BuildRefused` before writing anything if validation reports errors.
    Nothing is written on *that* refusal — not a partial `dist/`, not an empty
    directory — because a half-built site a deploy step then uploads is worse
    than no build at all.

    That is not a guarantee about every failure, and saying so once read as one.
    The builders below carry their own guards, and `validate/` does not yet check
    everything they do — a variant shard named `chr12.tsv` validates clean and
    raises here — so a `ValueError` from one of them can arrive after earlier
    builders have written. `out` may then hold a partial site with no manifest.
    The manifest is the signal: it is written last and seals the emitter, so a
    `dist/` without one is a build that did not finish. `cli.build` says so in
    as many words rather than leaving a caller to infer it.

    The returned mapping covers every artifact including `manifest.json`, which
    the manifest itself cannot list. `dict()` is hygiene and nothing more: the
    emitter is local and sealed, so it is unreachable once this returns and no
    test can tell a copy from the original. Measured, not assumed — the mutant
    returning `emitter.checksums` directly survives the whole suite, which is
    what an equivalent mutant looks like rather than a gap in the tests.
    """
    report = validate_repository(root)
    if not report.ok:
        raise BuildRefused(
            f"refusing to build: {report.error_count} validation error(s)\n{report.render()}"
        )

    # Issues are discarded rather than re-reported: `validate_repository` above
    # loaded the same corpus and its report is what the gate just passed, so a
    # second copy here could only ever say the same thing or disagree with the
    # decision already made.
    corpus, _ = load_curation(root)
    genes = _gene_registry(root)
    # The reader `chd-atlas validate` reports GEN004 through, so the two
    # cannot read a cell two ways (review, 2026-10-27). `None` is an
    # unreadable registry, which TBL008 refuses at the gate above; `{}` is
    # then the fallback `_gene_registry` has, unreachable in the same way.
    identities = read_gene_identities(root) or {}

    # Unreachable behind the gate above in the same way `_gene_registry`'s
    # empty-mirror fallback is: TBL012 is an error, so `validate_repository`
    # already refused a repository whose validity mirrors do not read, and this
    # line never runs on one that reaches here. Read again anyway and `raise`
    # rather than assume it: `-O` strips `assert`, and this project keeps a
    # guard on every bypassed gate (`Emitter.seal`, `encode_json`'s
    # `allow_nan=False`) rather than trusting "the gate already checked this" to
    # stay true as the two modules evolve apart.
    scope_terms = {str(entry.id) for entry in corpus.chd_scope}
    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    if clingen is None or gencc is None:
        raise ValueError("a validity mirror could not be read; the gate should have refused first")
    validity = gene_validity(clingen, gencc, in_scope=scope_terms)

    # Computed once and reused by both callers below: `build_genes` needs the
    # symbol alone, and so does `build_landing`, which is not a second read of
    # the mirror either — both project the same `genes` registry the same way.
    symbols = {gene: labels.symbol for gene, labels in genes.items()}

    # D21's population, computed once for the same reason `symbols` is, and with
    # more riding on it: `build_genes` publishes one index row and one bundle per
    # member, `build_search` indexes one record per member so that every one of
    # those bundles is reachable from the search box, and `build_landing` prints
    # how many there are. Deriving it three times would let the front page state a
    # figure the browse payload contradicts, or let the search index advertise a
    # `genes/<id>.json` no builder wrote — the shape the "154 genes published"
    # defect already had once.
    published = published_genes(validity)

    # Read once and handed to both builders below, for the same reason
    # `published` and `symbols` are: a gene page and the bundle it links to
    # rendering different burden rows is the failure a single source removes.
    burden = load_burden(root)
    # Loaded once beside `burden` for the same reason: the bundle, the page and
    # the manifest census all read this one mapping.
    combined = load_combined(root)

    emitter = Emitter(root=out)
    # `build_omics` skips this table outright -- it is keyed on `_GENE_COLUMN`,
    # and a quantile grid has no gene column -- so without this call
    # `mirrors/profile_quantiles/` is mirrored, schema-validated, sort-checked
    # and sha256'd, and reaches no published byte. That is the one file whose
    # entire purpose is letting a consumer check a published percentile
    # (D39(b)), so it would be the one thing nothing could fetch, on a build
    # every other check reports clean. Must run before `write_manifest`,
    # which seals the emitter and would refuse a write placed after it.
    #
    # Moved ahead of `build_omics` (Task 8b, which closes what a comment here
    # used to record as open wiring): `select_top` ranks a `profiles` row on
    # the percentile `gene_expression_profiles` derives below, so that figure
    # has to exist before `build_omics` runs, not be computed from the rows
    # it has already ranked.
    #
    # Returns `{accession: shard path}`, consumed immediately below by
    # `gene_expression_profiles` -- the same way `ModalitySummary.shards`
    # links a gene bundle to an omics shard today. Reusing the exact path
    # this call wrote, rather than reconstructing it with a second call to
    # `slug`, is what keeps the two from drifting apart.
    quantile_shards = build_profile_quantiles(root, emitter)
    # Pure derivation, no `emitter`: reads `mirrors/profiles/*.tsv` and
    # `mirrors/profile_quantiles/*.tsv` directly and returns one
    # `ExpressionProfile` per gene those mirrors mention, published or not --
    # `build_genes` is what restricts the result to `published` genes, the
    # same restriction it already applies to `omics` and `variants`.
    # `corpus.cardiac_phases` is `None` only when `curation/cardiac_phases.yaml`
    # is absent or unparsable -- the committed corpus's own file loads to a
    # real `CardiacPhaseFile` with zero phases *declared* instead (Task 4's
    # placeholder, boundaries not yet transcribed from a verified source).
    # `assign_phase`/`_phase_entry` treat both states identically, as "outside
    # the curated window" rather than raising, so this call is correct either
    # way without this function needing to tell the two apart.
    expression_profiles = gene_expression_profiles(
        root, corpus.datasets, corpus.cardiac_phases, quantile_shards
    )
    # One derivation, two consumers -- `build_landing`'s new card and
    # `write_manifest`'s `counts` below -- the same discipline `census`
    # (burden) already keeps: computed once here so the front page and the
    # manifest cannot state two different censuses of one build. Restricted to
    # `published`, not to what `mirrors/profiles/*.tsv` happens to mention --
    # see `profile_census`'s own docstring for why neither figure has the
    # `families`-style exception `burden_census` carries.
    profile_stats = profile_census(expression_profiles, published)
    # Derived once, like the profiles above: `build_genes` restricts it to
    # `published`, and `cell_type_census` counts the same objects the bundles
    # carry, so the manifest cannot disagree with them. The gene pages and the
    # front page are handed the same mapping and census below.
    cell_type_expression = gene_cell_type_expression(root, corpus.datasets)
    cell_type_stats = cell_type_census(cell_type_expression, published)
    # Derived once, like the layers above, and handed to the bundles and the
    # pages alike; `literature_census` counts the published genes' entries, so
    # the manifest states what the bundles carry (ruling R16).
    literature = gene_literature(corpus)
    literature_stats = literature_census(literature, published)
    # The one flat read of the `Placement`s just computed above -- never a
    # second derivation of a percentile. `percentile_annotations` only walks
    # the nested `ExpressionProfile`s and copies `median_percentile` back out;
    # see its own docstring for why that keeps the two uses below from
    # drifting apart. `build_omics` is the sole consumer of this return; the
    # bundle itself gets `expression_profiles` directly, a few lines below,
    # from this same call's result.
    percentile_lookup = percentile_annotations(expression_profiles)
    # Keyed by accession (as `str`, not the `AccessionId` newtype, to match
    # `profiles.dataset` cell-for-cell) so `build_omics` can resolve, for
    # whichever datasets contributed rows to one gene's profiles slice, which
    # of *that* dataset's own tissue tokens is the heart. `frozenset()` for a
    # contrast dataset is the correct answer, not a gap: that design never
    # reaches the profiles path at all.
    cardiac_tissues = {
        str(dataset.id): frozenset(dataset.cardiac_tissues) for dataset in corpus.datasets
    }
    # The whole curated record per dataset, not only its `cardiac_tissues`
    # projection above -- `build_gene_pages`'s expression section also needs
    # `detection_floor` to name the value a below-floor caveat requires, the
    # same "publish the whole record, not a projection" rule `cohort_registry`
    # follows for `Cohort`. A second dict rather than widening `cardiac_tissues`
    # itself: that mapping is `build_omics`'s own parameter today, and a second
    # consumer reading a wider shape through the same name is a needless
    # coupling between the two.
    dataset_registry = {str(dataset.id): dataset for dataset in corpus.datasets}
    # `select_top` ranks the cardiac series on `percentile_lookup` above, so
    # the rank a reader's bundle preview is chosen by and the percentile the
    # bundle itself publishes (via `expression_profiles`, handed to
    # `build_genes` below) are the same number by construction -- there is
    # exactly one computation of any gene's percentile in this build.
    #
    # This call is not safe to move back above `quantile_shards` /
    # `expression_profiles` / `percentile_lookup`: `percentile_lookup` would
    # not exist yet at that point in the function, and Python raises
    # `UnboundLocalError` rather than silently falling back to
    # `percentiles=None` -- measured directly, by making that exact edit and
    # running the suite, in preference to assuming it (every test that calls
    # `build_site` fails on it, since the name is unbound regardless of
    # whether the corpus being built has any profiles data at all).
    omics = build_omics(root, emitter, cardiac=cardiac_tissues, percentiles=percentile_lookup)
    variants = build_variants(root, emitter)
    # `facts` rather than a second `gene_facts` call below: the pages and the
    # bundles render from one derivation, so a page cannot state a confidence the
    # bundle it links to contradicts. See `build_genes`' docstring.
    # Derived once, here, and handed to all three consumers -- the bundles that
    # publish it, the gene pages that render it as a matrix, and the browse page
    # that renders it as a strip. Deriving it in each would type-check and leave
    # three call sites free to drift, which is the reason `facts` and `validity`
    # are threaded down rather than recomputed.
    #
    # A cohort family is a property of the *corpus*: which studies share sample
    # collections. So it is computed over every burden row, not per gene, or two
    # genes would show different numbers of family slots and an untested dataset
    # would stop looking like one.
    all_burden_rows = [row for rows in burden.values() for row in rows]
    families = cohort_families(all_burden_rows, cohort_registry(corpus.cohorts))
    axes = evidence_axes(all_burden_rows)
    # One census, two consumers: the front page states it in prose and the
    # manifest publishes it as `counts`. Derived here rather than in either, so
    # neither can drift from the other.
    census = burden_census(burden, published, families, survives=survives_own_correction)
    concordance = {
        gene: gene_concordance(burden.get(gene, ()), families) for gene in sorted(published)
    }

    facts = build_genes(
        corpus,
        emitter,
        symbols=symbols,
        identities=identities,
        omics=omics,
        variants=variants,
        validity=validity,
        published=published,
        burden=burden,
        concordance=concordance,
        profiles=expression_profiles,
        combined=combined,
        cell_types=cell_type_expression,
        literature=literature,
    )
    build_literature(corpus, emitter)
    # The resolution table for the bare cohort ids every burden row carries.
    # Emitted from `corpus.cohorts` -- the whole curated registry, not the subset
    # the published rows happen to cite -- so that every id in every row resolves
    # regardless of where the publication gate sits. See `build_cohorts`.
    build_cohorts(corpus.cohorts, emitter)
    # Read again rather than threaded down from the gate: `validate_repository`
    # discards its registry, and re-reading one small YAML costs less than
    # widening the gate's return type to carry a value only the build wants.
    # It cannot fail here — SRC001 is an error, so the gate refused already.
    registry, _ = load_sources(root)
    build_sources(registry, emitter)
    # The Orphanet warrants, derived once from the mirrors by the function SCP006
    # refused on, and handed to both `scope.json` and the gene pages, so the
    # registry and the note under a gene's validity table state one warrant.
    warrants = scope_warrants(
        corpus.chd_scope,
        load_orphanet_tables(root),
        pinned_structural_phenotypes(root, registry.ontology_files()),
    )
    build_scope(corpus.chd_scope, warrants, emitter)
    # The site's reading notes (D50): the full forms of what the gene pages'
    # short forms and chart toggles link to. Every page's nav links it, so a
    # build that skipped this would publish a dead link on all of them.
    build_methods(emitter, cohorts=cohort_registry(corpus.cohorts))
    build_search(corpus, emitter, genes=genes, published=published)
    # One integer, for the manifest below; the front page stopped stating it
    # with the doorway page (2026-10-27).
    combined_results = combined_census(combined, published)
    build_landing(
        corpus,
        symbols=symbols,
        validity=validity,
        published=published,
        census=census,
        profile_census=profile_stats,
        cell_type_census=cell_type_stats,
        emitter=emitter,
        # For the Orphadata citation's data version (R11), read off the files'
        # own headers into `mirrors/sources.yaml` by the converter.
        sources=registry,
        facts=facts,
        concordance=concordance,
        profiles=expression_profiles,
        burden_counts={gene: len(rows) for gene, rows in burden.items()},
        publications={publication.id: publication for publication in corpus.publications},
        identities=identities,
    )
    # The HTML over everything above. Wired here and nowhere else: until this
    # call existed, `pages.py` was imported by its unit test alone, `build_landing`
    # and the shared `<nav>` both linked to `genes/index.html`, and no build wrote
    # one — a green build, every checksum verifying, and every visitor clicking
    # "Genes" served a 404.
    #
    # Keyed by id rather than by position, because both builders look records up
    # by identifier: `publications` resolves the PMID an assertion's evidence
    # cites, the same lookup `literature._featured` performs over the same tuple.
    build_gene_pages(
        facts,
        emitter,
        symbols=symbols,
        identities=identities,
        validity=validity,
        assertions=_assertions_by_gene(corpus),
        publications={publication.id: publication for publication in corpus.publications},
        burden=burden,
        cohorts=cohort_registry(corpus.cohorts),
        families=families,
        axes=axes,
        # The same `expression_profiles` the bundles were built from, a few
        # lines above -- so a gene page cannot show developmental expression
        # data its own bundle does not carry, or vice versa.
        profiles=expression_profiles,
        datasets=dataset_registry,
        # The same vocabulary `gene_expression_profiles` assigned each stage
        # its phases from, a hundred lines above -- so a band drawn on a
        # trajectory and the phase sentence under the stage it covers cannot
        # come from two different readings of `curation/cardiac_phases.yaml`.
        phases=corpus.cardiac_phases,
        # The same mapping the bundles were built from, so a page cannot state a
        # combination its own bundle does not carry.
        combined=combined,
        cell_types=cell_type_expression,
        # The mapping the bundles carry, so a page cannot state a search or a
        # finding its own bundle does not.
        literature=literature,
        scope_warrants=warrants,
        scope_labels={str(e.id): e.label for e in corpus.chd_scope},
        # The registry itself, for the card's line on a syndrome a ClinGen
        # panel scoped (R5): it reads each term's `syndrome` mark and authority.
        scope_terms={str(e.id): e for e in corpus.chd_scope},
    )
    # `validity` again, and the same object `build_gene_pages` was handed: the
    # browse row's `definitive for` cell and the gene page's `definitive for`
    # rail row must name one disease, and they do because both read this mapping
    # rather than either re-deriving it from the mirrors.
    build_gene_index_page(
        facts,
        emitter,
        symbols=symbols,
        validity=validity,
        # The same count `build_genes` publishes as `burden_row_count`, from
        # the same mapping, so the browse page and the payload behind it
        # cannot disagree about how much evidence a gene has.
        burden_counts={gene: len(rows) for gene, rows in burden.items()},
        # The same objects the bundles publish, so the strip a reader scans
        # and the JSON a program reads cannot disagree.
        concordance=concordance,
        # For the summary caveats of the studies the strips name, which the
        # gene pages' cards print from the same records.
        publications={publication.id: publication for publication in corpus.publications},
        profiles=expression_profiles,
        identities=identities,
    )
    # Last, and enforced as last: this seals the emitter.
    #
    # The build counts come from `census` and `profile_stats`, the same objects
    # `build_landing` was handed a moment ago, so `manifest.json` and
    # `index.html` cannot publish two censuses of one build. `genes` is
    # `len(published)` rather than `census.genes`: they answer different
    # questions — how many genes the site publishes, and how many of those
    # carry burden evidence — and they differ by every published gene no
    # rare-variant study covered. They were equal at 23 of 23 when this was
    # written and have since parted, which is the day this was written for.
    write_manifest(
        corpus,
        emitter,
        commit=source_commit(root),
        counts={
            "genes": len(published),
            "burden_rows": census.rows,
            # Combined results reaching a published bundle, never the mirror's
            # row count -- the same convention as `burden_rows`.
            "burden_combined": combined_results,
            "cohort_families": census.families,
            # The genes whose card lists a surviving result (R9), by the
            # card's own predicate through `census`.
            "burden_surviving_genes": census.surviving_genes,
            "profile_genes": profile_stats["genes"],
            "profile_datasets": profile_stats["datasets"],
            "cell_type_genes": cell_type_stats["genes"],
            "cell_type_datasets": cell_type_stats["datasets"],
            "literature_genes": literature_stats.genes,
            "literature_findings": literature_stats.findings,
        },
    )
    return dict(emitter.checksums)
