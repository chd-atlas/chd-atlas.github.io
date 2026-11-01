# src/chd_atlas/build/bundles.py
"""The gene index and the per-gene bundles.

A bundle is what replaces a backend: a gene detail page is one fetch of one
small file, joined here so that no client ever has to assemble a gene's
assertions, functional records, variants and omics summaries itself.

`genes/index.json` is the other half: the browse payload, downloaded by every
visitor before they have picked a gene. It carries what the browser ranks and
filters on plus the path to fetch the rest, and none of the evidence itself —
that is what keeps it small as curation grows, and why the bundle is a separate
file rather than a section of this one.

Two rules decide what goes where.

A contested gene must never be displayed as settled. `strongest()` ranks
classifications on one linear scale where `definitive` outranks `refuted`, while
ClinGen treats disputed and refuted as a separate axis rather than weaker rungs
of the same ladder, so `headline_confidence` on its own hides a refutation
behind a stronger assertion about the same gene. `has_conflicting_evidence` is
the other half of that pair; `_headline` below is the single place either is
written, so a payload cannot come to carry the confidence without the flag.

Every path the index advertises must be a file that was written. `bundle` is the
URL a detail page is fetched from, and nothing downstream can catch a wrong one:
`Emitter.checksums` records what the build wrote, never what a payload promised,
so a dead link is published with every checksum verifying and every gate green.
One `gene_bundle_path` call per gene supplies both the index entry and the
write, which is what makes the two agree by construction rather than by review.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping, Sequence
from typing import Any, cast

from chd_atlas.build.burden import BurdenRow, CombinedRow, burden_payload, combined_payload
from chd_atlas.build.cell_types import EMPTY_CELL_TYPE_EXPRESSION, CellTypeExpression
from chd_atlas.build.derive import GeneFacts, gene_facts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.identity import EMPTY_IDENTITY, GeneIdentity, identity_payload
from chd_atlas.build.omics import ModalitySummary
from chd_atlas.build.paths import gene_bundle_path
from chd_atlas.build.profiles import EMPTY_EXPRESSION_PROFILE, ExpressionProfile
from chd_atlas.build.validity import (
    GeneValidity,
    ValidityRecord,
    admission_provenance,
    uncurated,
)
from chd_atlas.corpus import Corpus
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.functional import FunctionalEvidence


def _records_by_gene[RecordT: (LesionAssertion, FunctionalEvidence)](
    records: Iterable[RecordT],
) -> dict[str, list[dict[str, Json]]]:
    """Group one kind of curated record by gene, ordered by id, JSON-ready.

    Ordered by id rather than left in corpus order, because these are published
    as JSON arrays and `encode_json`'s `sort_keys` orders dict keys and nothing
    else: whatever order they reach the emitter in is the order a consumer
    downloads. Corpus order is `load_curation`'s: record files sorted by
    filename — `_record_files` sorts the glob — and then each file's own order
    within it, which makes it an artefact of how the curation is filed rather
    than of the records. Without this, moving a record from `TBX5.yaml` into a
    second file would reorder a bundle that is otherwise unchanged.

    `mode="json"` so every value is a JSON primitive by construction.
    `LesionAssertion` carries two `date` fields, which `json.dumps` refuses
    outright — without this the build dies part-way through, with `dist/` already
    half written. `FunctionalEvidence` carries none today and would serialise
    either way, which is exactly why the conversion is applied to both rather
    than to the one that needs it.

    The type parameter is constrained to the two record types rather than bound
    to `BaseModel`, because the body reads `.gene` and `.id`, which `BaseModel`
    does not declare. Constrained means mypy checks this body once per record
    type, so a third record type has to be added here deliberately.
    """
    grouped: dict[str, list[RecordT]] = {}
    for record in records:
        grouped.setdefault(record.gene, []).append(record)
    return {
        gene: [record.model_dump(mode="json") for record in sorted(items, key=lambda r: r.id)]
        for gene, items in grouped.items()
    }


def _concordance_for(concordance: Mapping[str, Json], gene: str) -> Json:
    """One gene's `independent_datasets`, or a refusal naming the gap.

    `raise`, never a silent default. A gene present in `published` and absent
    from this mapping means the caller derived it over a different population,
    and the failure a default would produce -- a bundle publishing
    `tested: 0` -- is indistinguishable from the true statement that no study
    reported that gene. `raise` rather than `assert`: `-O` strips `assert`.
    """
    if gene not in concordance:
        raise KeyError(
            f"no concordance derived for published gene {gene}; the mapping was built over "
            f"a different population, and defaulting it would publish 'no dataset tested "
            f"this' as if it had been measured"
        )
    return concordance[gene]


def _summaries(modalities: Mapping[str, ModalitySummary]) -> Mapping[str, Json]:
    """One gene's omics summaries, retyped rather than rebuilt.

    `ModalitySummary` is a `TypedDict`, and a `TypedDict` is assignable only to
    `Mapping[str, object]` — never to `Mapping[str, Json]`, whatever its fields
    hold. That is a structural fact about `TypedDict`, not a judgement about
    these values, so `emit.Json` cannot accept one and the cast is the whole fix.

    The alternative is to rebuild the summary from the three fields it declares
    today. That type-checks with nothing cast, and it is worse: a field added to
    `ModalitySummary` would be dropped here silently, and the failure would be a
    gene page missing data that the shard it links to still contains. Passing the
    value through instead makes the same change either work or fail loudly, since
    `encode_json` raises on anything `json.dumps` cannot serialise.
    """
    return cast(Mapping[str, Json], modalities)


def _expression_profile(profiles: Mapping[str, ExpressionProfile], gene: str) -> Json:
    """One gene's developmental expression profile, always present.

    `profiles.get(gene, EMPTY_EXPRESSION_PROFILE)` -- never a `KeyError` guard
    like `_concordance_for`'s, because a gene absent from `profiles` is not a
    caller error the way a gene absent from `concordance` is: most genes have
    no profile dataset covering them at all (there is no profile mirror
    committed today), and `EMPTY_EXPRESSION_PROFILE` is exactly the shape
    `gene_expression_profiles` itself never emits an entry for, so a consumer
    reads one shape rather than guarding for a missing key.

    `ExpressionProfile` is a `TypedDict`, assignable only to
    `Mapping[str, object]` and never to `Mapping[str, Json]` whatever its
    fields hold -- the same structural fact `_summaries` documents for
    `ModalitySummary`, and the cast is the whole difference from rebuilding
    the payload key by key, which would silently drop a field a later task
    adds to `ExpressionProfile`.
    """
    return cast(Json, profiles.get(gene, EMPTY_EXPRESSION_PROFILE))


def _cell_type_expression(cell_types: Mapping[str, CellTypeExpression], gene: str) -> Json:
    """One gene's single-cell expression by cell type, always present.

    `EMPTY_CELL_TYPE_EXPRESSION` for a gene no cell-type dataset covers -- the
    rule `_expression_profile` keeps -- so a consumer reads one shape. Cast
    for the reason `_expression_profile` is: a `TypedDict` is never a `Json`.
    """
    return cast(Json, cell_types.get(gene, EMPTY_CELL_TYPE_EXPRESSION))


def _validity_record(record: ValidityRecord) -> dict[str, Json]:
    """One mirrored classification, JSON-ready.

    Every record in the published array carries the same thirteen keys, whichever
    source it came from. `sop`, `classification_date` and `gcep` are ClinGen-only
    and `submitter`, `submitted_on` and `sgc_id` are GenCC-only on `ValidityRecord` itself (see its
    docstring), and each is `null` on a record from the other source rather than
    left out of the object. The alternative — omitting a source's absent fields
    — would publish an array whose objects differ in shape, which is a trap for
    any consumer that reads one field off `records[0]` and expects the same key
    to exist on `records[1]`.

    `classification` publishes `null` for GenCC's `Supportive`, which
    `ValidityRecord.classification` already represents as `None` — the mapping
    exception documented on that field — while `classification_term` always
    carries the authority's own word, mapped or not.
    """
    return {
        "source": record.source.value,
        "classification_term": record.classification_term,
        "classification": record.classification.value if record.classification else None,
        "disease": record.disease,
        "disease_label": record.disease_label,
        "moi": record.moi,
        "sop": record.sop,
        "classification_date": record.classification_date,
        "gcep": record.gcep,
        "report_url": record.report_url,
        "submitter": record.submitter,
        "submitted_on": record.submitted_on,
        "sgc_id": record.sgc_id,
    }


def _validity(gene_validity: GeneValidity) -> dict[str, Json]:
    """The bundle's `validity` block: how curated the gene is, and by whom.

    `records` is serialised in the order `build.validity.gene_validity` already
    sorted it into — content order, not encounter order, so two builds of one
    commit agree on it. Never re-sorted or re-derived here: that module already
    did the work of deciding the order, and doing it again risks disagreeing
    with it instead of reusing it.

    `state` and `has_source_discordance` duplicate the top-level
    `validity_state` and `has_source_discordance` that `_headline` already
    writes into this same bundle. That is deliberate rather than an oversight:
    the top-level fields are what the browse row and the bundle publish
    identically (`_headline`'s whole reason to exist), while this block is the
    self-contained provenance object — a consumer who only wants "who curated
    this gene, and what did each of them say" reads `validity` alone and never
    has to cross-reference the flat fields beside it.
    """
    return {
        "state": gene_validity.state.value,
        "has_source_discordance": gene_validity.has_source_discordance,
        "records": [_validity_record(record) for record in gene_validity.records],
    }


def _headline(gene: str, symbol: str, fact: GeneFacts) -> dict[str, Json]:
    """What the index row and the bundle must say identically about one gene.

    One function rather than two literals so the browse row and the page it
    opens cannot disagree, and — the reason it exists — so `headline_confidence`,
    `validity_state`, `atlas_curation`, `has_conflicting_evidence`,
    `has_source_discordance`, `has_no_association_report` and
    `no_association_reported_by` are written together. Dropping any one of them
    from one payload alone is not an edit that can be made here by accident.

    `.value` publishes the vocabulary's string rather than leaning on a `StrEnum`
    member *being* a `str`. Measured, not assumed: mypy accepts either spelling
    inside a payload this size, and `json.dumps` writes the same bytes for both
    today. What differs is a vocabulary that stops being a `StrEnum` — the bare
    member would then fail type checking here and be refused by `json.dumps`,
    while `.value` keeps publishing the same JSON. Neither is a trap; `.value` is
    the one that states at the call site what the file will contain.

    `headline_confidence` publishes `null` for a gene no authority has assessed
    — `fact.headline_confidence` is `None` in that case, per `derive.gene_facts`
    — rather than coercing it to a string, which would invent a confidence
    nobody stated.
    """
    return {
        "gene": gene,
        "symbol": symbol,
        "headline_confidence": (
            fact.headline_confidence.value if fact.headline_confidence else None
        ),
        "validity_state": fact.validity_state.value,
        # Written inside `_headline` for the same reason `has_conflicting_evidence`
        # is: the browse row and the page it opens must not be able to disagree
        # about whether the atlas has curated a gene, and nearly every published
        # gene has not (118 of 119, measured 2026-10-26).
        "atlas_curation": fact.atlas_curation.value,
        "has_conflicting_evidence": fact.has_conflicting_evidence,
        "has_source_discordance": fact.has_source_discordance,
        # The third axis, written here for exactly the reason the two above are:
        # a reader who meets the flag on the browse row and not on the page, or
        # the other way round, meets it exactly where it is missing. The list is
        # empty whenever the flag is false, so the two cannot say different
        # things about the same gene. Issue #13.
        "has_no_association_report": fact.has_no_association_report,
        "no_association_reported_by": list(fact.no_association_reported_by),
        "lesion_groups": [group.value for group in fact.lesion_groups],
    }


def build_genes(
    corpus: Corpus,
    emitter: Emitter,
    symbols: Mapping[str, str],
    omics: Mapping[str, Mapping[str, ModalitySummary]],
    variants: Mapping[str, list[dict[str, Any]]],
    validity: dict[str, GeneValidity],
    published: Collection[str],
    burden: Mapping[str, Sequence[BurdenRow]],
    concordance: Mapping[str, Json],
    profiles: Mapping[str, ExpressionProfile],
    combined: Mapping[str, Sequence[CombinedRow]] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
    *,
    cell_types: Mapping[str, CellTypeExpression],
) -> dict[str, GeneFacts]:
    """Emit `genes/index.json` and one bundle per published gene, and return the facts.

    The return is `gene_facts`' mapping, handed back rather than discarded so
    that `build_gene_pages` and `build_gene_index_page` render from the same
    derivation these bundles were written from. Deriving it a second time in the
    runner would type-check, cost another pass over the corpus, and leave two
    call sites free to drift apart — a page and the bundle it links to
    disagreeing about a gene's confidence is exactly the failure `_headline`
    exists to make unrepresentable within one payload, and it would reappear
    between payloads. Returning it does not widen what this function promises:
    the mapping is `gene_facts`' own return, keyed on the same population, and
    this function does not mutate it.

    `published` is D21's population: `build.validity.published_genes()`'s
    return, the genes a ClinGen expert panel classifies `Limited` or better for
    an in-scope disease, plus the genes two or more GenCC submitters agree on
    that no ClinGen panel contests. It is passed to `gene_facts`, which keys on
    it, so this function writes exactly one bundle per member and the index lists
    exactly the same set. Nearly all carry no curated assertion at all
    (118 of 119, measured 2026-10-26); `atlas_curation` is what says so, in the
    index row and the bundle alike.

    `symbols` comes from `mirrors/genes.tsv`, keyed on HGNC id. A gene absent
    from it keeps its HGNC id as its display label: the browse row and the page
    heading both render this value, so blank reads as a broken site while
    "HGNC:4173" reads as a symbol that has not been mirrored yet and is still
    something a reader can search for.

    `validity` is passed straight to `gene_facts` — see its docstring for what
    a gene absent from it publishes. The bundle additionally embeds the whole
    `GeneValidity` as a `validity` object (`_validity`), which `gene_facts`
    does not carry forward because the browse layer never needs the individual
    mirrored records, only the fields it already reduces `GeneValidity` to.

    `omics` and `variants` are what `build_omics` and `build_variants` returned,
    taken as `Mapping` because nothing here mutates them. Genes they carry that
    are outside `published` are ignored, which is the same rule `gene_facts`
    applies: no external authority has admitted them under D21, so the atlas
    publishes no page for them. Their rows are still published in the shards
    those two modules wrote — the gene index simply does not link to them.

    `combined` is `burden.load_combined`'s return. It defaults to `None` so a
    caller that predates it publishes an empty `burden_combined` array rather
    than failing.

    `identities` is `identity.read_gene_identities`' return: `mirrors/genes.tsv`'s
    identity columns (D52). A gene absent from it publishes an `identity`
    object of six nulls, the shape every bundle carries -- the same "keep the
    record, mark what is missing" fallback `symbols` applies. Defaulted, so a
    caller that forgets it is not refused here: it publishes six nulls on every
    bundle, which the API doc defines as "the registry holds no value", and
    every unit test that passes `identities` itself stays green.
    `test_every_published_bundle_carries_the_identity_its_registry_row_holds`
    is what catches it, by comparing every bundle of a real build with
    `mirrors/genes.tsv`.

    `cell_types` is `gene_cell_type_expression`'s return, published or not, and
    restricted here the way `profiles` is below. Required, with no default,
    unlike `combined` and `identities`: no committed dataset carries cell-type
    data yet, so a real build cannot tell a dropped mapping from an empty one,
    and with a default, deleting it from the runner's call passed the whole
    suite while every bundle published `{"datasets": []}` (review, 2026-10-28).
    A caller with no cell-type data passes `{}`.

    `profiles` is `gene_expression_profiles`' return (`build/profiles.py`):
    the derived percentile placement, tau and phase assignment for every gene
    `mirrors/profiles/*.tsv` mentions, published or not — that function takes
    no `published` argument at all, so restricting it is entirely this
    function's job, by the same iteration over `published` that already
    restricts `omics` and `variants`. The same asymmetry as those two
    applies: a merely-registered gene's *derived* facts reach no file (there
    is no `expression_profile`-only shard, any more than there is a
    `ModalitySummary`-only one), but its *raw* `profiles.tsv` measurements are
    not lost — `build_omics` already treats `profiles` as one of its
    dataset-linked tables and shards it exactly like `expression`,
    `proteomics` and `phospho`, independent of gene publication.

    That last paragraph, and the one before it, cover `omics`, `variants` and
    `profiles` and nothing else. A
    `FunctionalEvidence` record about a gene outside `published` is worse off: no
    other build module reads `corpus.functional` — `derive.py` reads it only to
    count records per gene, and this is its only writer — so it reaches no
    published file at all rather than an unlinked one. That is evidence loss,
    not a missing link. It still does not raise here, because such a record is
    legal under every validator in the project and `gene_facts` is what decides
    which genes exist; the cost is pinned instead, by
    `test_a_gene_with_evidence_but_outside_the_published_set_is_not_published`.
    """
    facts = gene_facts(corpus, validity, published=published)
    assertions = _records_by_gene(corpus.assertions)
    functional = _records_by_gene(corpus.functional)

    index: list[Json] = []
    # Sorted here as well as in `gene_facts`, which returns its genes in this
    # same order — pinned by `test_derive.py`'s
    # `test_assertions_are_counted_per_gene_and_the_genes_come_back_sorted`.
    # Redundant today, and kept anyway: this array is what consumers download,
    # so its order is the published contract and should not rest on a guarantee
    # another module is free to withdraw.
    for gene in sorted(facts):
        fact = facts[gene]
        # `HgncId` rather than `str` at this boundary: `gene_facts` keys on a
        # plain `str` by choice, but every key it returns is a member of
        # `published`, which `build_site` derives from the `gene` column of
        # `mirrors/clingen_gene_validity.tsv` — checked against `HGNC_PATTERN` by
        # that table's schema, and TBL-level failures are errors the gate refuses
        # on. The NewType is being recovered, not asserted. `gene_bundle_path` requires
        # it so that only HGNC ids can name a file in `genes/`, where a
        # `ContrastId` would slug to a name that collides with one on a
        # case-insensitive filesystem.
        bundle = gene_bundle_path(HgncId(gene))
        headline = _headline(gene, symbols.get(gene, gene), fact)
        rows = variants.get(gene, [])
        # Looked up again rather than threaded through `GeneFacts`: `gene_facts`
        # already reduced `GeneValidity` to the scalar fields the browse layer
        # needs (`headline_confidence`, `validity_state`, ...) and does not carry
        # the record list forward, so the bundle -- the one payload that needs
        # the records themselves -- reads the same `validity` mapping this
        # function was handed, with the same `uncurated()` fallback `gene_facts`
        # documents for a gene neither mirror curates.
        gene_validity = validity.get(gene, uncurated())

        index.append(
            {
                **headline,
                "confidence_by_lesion_group": {
                    group.value: classification.value
                    for group, classification in fact.confidence_by_lesion_group.items()
                },
                # Written immediately beside the map it qualifies, for the reason
                # `_headline` exists: a confidence a reader can act on must never
                # be published without the flag that says it is contested. At
                # gene level one function guarantees that; here it is adjacency
                # plus the test that reads both, because the two are different
                # shapes and cannot share a writer.
                "conflicting_lesion_groups": [
                    group.value for group in fact.conflicting_lesion_groups
                ],
                "evidence_counts": {
                    evidence_class.value: count
                    for evidence_class, count in fact.evidence_counts.items()
                },
                "assertion_count": fact.assertion_count,
                "functional_count": fact.functional_count,
                # Counted from the same list the bundle embeds, so the browse
                # promise and the page cannot differ.
                "variant_count": len(rows),
                # Counted from the same list the bundle embeds, for the reason
                # `variant_count` is: a browse row promising burden evidence a
                # page does not show is the browse layer lying about the page.
                "burden_row_count": len(burden.get(gene, ())),
                # Computed once and published on both the browse row and the
                # bundle below, so the two layers cannot disagree -- the rule
                # `burden_row_count` follows. It reaches JSON at all because the
                # 2026-08-05 review found `curation/cohorts.yaml` rendering only
                # in HTML and reaching zero published bytes; a fact that lives
                # only on a page is a fact no consumer has.
                # `concordance` is required rather than defaulted: a default
                # would let a caller forget it and publish zeros indistinguishable
                # from a measured "no dataset found anything", which is this
                # project's characteristic failure. A missing gene is a caller
                # that built the mapping over the wrong population, so it fails
                # loudly here rather than writing a bundle that quietly claims
                # nothing was tested.
                "independent_datasets": _concordance_for(concordance, gene),
                "bundle": bundle,
            }
        )

        emitter.write_json(
            bundle,
            {
                **headline,
                # Who the gene is, from the registry row the symbol came from:
                # name, locus and the ids its page links out by (D52).
                "identity": identity_payload((identities or {}).get(gene, EMPTY_IDENTITY)),
                # The full mirrored picture, alongside the flat fields
                # `headline` already carries -- see `_validity`'s docstring for
                # why both are published rather than one replacing the other.
                "validity": _validity(gene_validity),
                # Why this gene is published at all, as a payload rather than a
                # rule a reader has to trust. `methods.SCOPE_RULE` tells a reader
                # that no disease is in scope on this atlas's own judgement;
                # this is what lets them check which authority admitted the gene
                # in front of them. Kept beside `validity` rather than inside it
                # because `validity.records` is the verbatim mirror and this is
                # a derivation over it -- and its `asserted_by.count` dedupes
                # ClinGen's own GenCC submissions, which the mirror does not.
                **admission_provenance(gene_validity),
                "publications": list(fact.publications),
                "assertions": assertions.get(gene, []),
                "functional": functional.get(gene, []),
                # Embedded rather than linked. Not because variants are bounded
                # per gene — nothing bounds them, and TTN alone carries tens of
                # thousands of ClinVar records — but because this atlas curates
                # variants by hand as CHD evidence, so the count is a curation
                # policy rather than a property of the data. If that policy
                # changes, this has to become a link the way omics rows are.
                "variants": rows,
                # Always present, empty when the gene has none, so a consumer
                # reads one shape rather than guarding for a missing key.
                "omics": _summaries(omics.get(gene, {})),
                # Always present, empty when the gene has none -- same rule as
                # `omics` above. An absent key would make "no study reported
                # this gene" indistinguishable from "the build dropped it".
                "burden": burden_payload(burden.get(gene, ())),
                # Always present, empty when the gene has none, for the reason
                # `burden` is. A study's own combination of its arms -- never a
                # pooling across studies -- whose `local_fdr` is not a corrected p.
                "burden_combined": combined_payload((combined or {}).get(gene, ())),
                # Always present, empty ({"datasets": []}) when no profile
                # dataset's mirror covers this gene -- same rule as `omics`
                # and `burden` above, for the same reason: an absent key
                # would make "this gene was not covered by that source
                # matrix" indistinguishable from "the build dropped it".
                "expression_profile": _expression_profile(profiles, gene),
                # Always present, empty ({"datasets": []}) when no cell-type
                # dataset covers this gene, for the reason `omics` is. A
                # profile of cell types, never a contrast: a low or absent
                # value is not evidence against the gene.
                "cell_type_expression": _cell_type_expression(cell_types, gene),
                # The same object the browse row carries -- literally the
                # same mapping, handed in by the runner -- so a consumer that
                # filtered on the index and then fetched the bundle reads one
                # answer rather than two.
                "independent_datasets": _concordance_for(concordance, gene),
            },
        )

    # Written last, so a gene whose bundle could not be written — two genes
    # normalising to one filename is what `Emitter` refuses — fails the build
    # with no index at all, rather than leaving one that advertises a file
    # nothing wrote.
    #
    # Emitted even when there is nothing to list: a consumer that has to tell
    # "no genes curated yet" from "wrong URL" by reading a 404 will get it wrong,
    # and this is the first thing the browse page fetches.
    emitter.write_json("genes/index.json", {"genes": index})
    return facts
