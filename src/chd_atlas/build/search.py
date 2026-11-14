# src/chd_atlas/build/search.py
"""The search index: one flat record per gene, publication and phenotype.

Deliberately a list of records rather than an inverted index. At this corpus size
a client filters the whole array in a fraction of a frame, while an inverted
index is a structure to maintain — kept in step with every payload it summarises,
and debugged when a query misses — before there is the query volume that would
need one.

Variants are excluded because the variant space is the one that grows without
bound, so indexing it would let the least likely query dominate the size of the
file every visitor downloads before they have typed anything.

Datasets are excluded for a different reason, which that sentence used to be read
as covering. Size is not the argument: a dataset is bounded by hand curation like
everything `literature.py` writes. What a dataset lacks is a name — the model
carries an accession, an archive, a technology, a tissue, a stage, an organism, a
sample count, a licence and its contrasts, and no title or description of its
own; the only free text describing it belongs to a contrast, whose `description`
names a comparison rather than the dataset — so the only string that would name
one record rather than the class is the accession. An accession already
resolves where a reader meets it: every omics row a gene bundle embeds names its
dataset in the `dataset` column, and `datasets.json` is what that resolves
against, exactly as `publications.json` resolves an assertion's PMID. A dataset
whose measurements were reported somewhere also reaches a reader through
`Dataset.publication`, and that paper is indexed here by title and author.
Revisit this when a dataset gains a title, or when the site gains a page for one
to land on; today an indexed dataset would be a row labelled with an accession.

The gene loop iterates `published` -- `build.validity.published_genes()`'s
return, the same population `build_genes` writes a bundle for. Keyed on the
asserted genes instead, as it was, the index would hold one gene while the
site published 23 bundles: a visitor typing GATA4 gets nothing while
`genes/HGNC_4173.json` is served. (The plan's wording here said
`genes/HGNC_4173.html`; no builder writes an HTML gene page yet, so the
extension is corrected rather than copied.) The dead-link sweep in
`tests/test_built_site_is_consumable.py` does not cover that direction --
it walks each published `path` and checks the file exists, so a record that
was never emitted breaks no link and leaves the sweep green.

Keyed on `genes` -- the whole 186-gene registry -- it would be worse: 67
results whose `path` names a bundle no builder wrote (2026-11-14). That direction the
sweep does catch. `published` is the only population that is right in both.

Amendment A29 records the sharper case, which this module does not change: a
`FunctionalEvidence` record about a gene outside `published` reaches no
published file at all, so nothing here can make it findable.
"""

from __future__ import annotations

from collections.abc import Collection, Iterable, Mapping
from dataclasses import dataclass

from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import PHENOTYPES, PUBLICATIONS, gene_bundle_path
from chd_atlas.corpus import Corpus
from chd_atlas.identifiers import HgncId


@dataclass(frozen=True)
class GeneLabels:
    """What `mirrors/genes.tsv` says about one gene that a visitor might type.

    One record rather than three `Mapping`s keyed on HGNC id. All three fields
    come from one row of one file, so parallel mappings could only ever disagree:
    a gene present in `symbols` and absent from `names` is a half-labelled record
    nothing would report, and every caller would carry three lookups and three
    fallbacks where one will do. `build_genes` is deliberately untouched — it
    needs the symbol alone and keeps taking `Mapping[str, str]`, which the build
    runner, once it is written, can derive from these in one comprehension rather
    than reading the mirror a second time. Nothing calls either function yet.

    `aliases` is `tuple[str, ...]` because `str` satisfies `Sequence[str]`,
    `Collection[str]` and `Iterable[str]`. Measured under `mypy --strict` rather
    than assumed: all three accept the raw pipe-separated cell, and a `str` handed
    to any of them explodes into single characters, publishing a search index of
    punctuation with every checksum verifying. `tuple`, `list` and `frozenset` all
    reject it; `tuple` is the immutable one, which this frozen record needs, and
    it keeps the mirror's order until the index sorts it.

    The guard closes what the annotation cannot. The registry reader that will
    build these reads the mirror through polars, whose `frame.to_dicts()` yields
    `dict[str, Any]`, and an `Any` satisfies `tuple[str, ...]` silently:
    constructing this straight from a row cell type-checks under `--strict`, and
    `reveal_type` still reports `tuple[str, ...]` while the value is a string.
    That is the footing `variants.py`'s shard guards and `encode_json`'s
    `allow_nan=False` already occupy — reaching it means a caller skipped the
    split, and a bypassed gate must fail rather than publish.

    `name` is optional here although the schema requires it on every row, because
    the same reasoning applies: a null reaching this means TBL003 was bypassed,
    and a gene searchable by symbol and id beats a build that dies on one blank
    cell. Callers holding only a symbol construct `GeneLabels(symbol=...)`.
    """

    symbol: str
    name: str | None = None
    aliases: tuple[str, ...] = ()
    prev_symbols: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        # Both sequence fields, not just the one that had the bug. They are
        # built from the same pipe-separated mirror cells by the same loop, so
        # a caller that can forget to split one can forget to split the other,
        # and a guard on one layer is not a guard (CLAUDE.md section 4.28).
        for field, value in (("aliases", self.aliases), ("prev_symbols", self.prev_symbols)):
            if isinstance(value, str):
                raise ValueError(
                    f"{field} for {self.symbol} arrived as one string, {value!r}; a `str` "
                    f"is a sequence of characters, so this would publish "
                    f"{len(set(value))} single-character search terms. Split the "
                    f"pipe-separated mirror cell before constructing this."
                )


def _terms(values: Iterable[str]) -> list[str]:
    """The strings a visitor might type, each kept once, in the order given.

    `dict.fromkeys` rather than a `set`: this list is published, so its order is
    part of the file that every visitor downloads, and set order is neither
    stable across interpreters nor meaningful to a reader. Insertion order is a
    language guarantee, so keeping the first occurrence of each term costs no
    determinism.

    Duplicates arrive on a green gate rather than only through a bypassed one.
    `mirrors/genes.tsv` is free to list a symbol among that gene's own aliases,
    an ontology export routinely repeats a term's label as one of its synonyms,
    and nothing forbids a publication naming one author twice. None of those
    changes a query's outcome — a record matching a term twice is one hit — so
    this buys bytes and a file that reads as intended, not correctness.
    """
    return list(dict.fromkeys(values))


def build_search(
    corpus: Corpus,
    emitter: Emitter,
    genes: Mapping[str, GeneLabels],
    published: Collection[str],
) -> None:
    """Emit `search/index.json.gz`, the whole searchable surface in one file.

    Each record carries `terms`, the strings a visitor might type, and `path`,
    the payload that answers the query once they pick the result. `label` is what
    the result row displays; `id` is what identifies the thing it names.

    Every record's `terms` lead with its label, then its own identifier, then
    whatever else names it. The identifier is there for all three kinds rather
    than only for genes: a PMID is the string `publications.yaml` is keyed on and
    the string every assertion's `evidence.publication` cites, so it is among the
    likeliest things typed into a search box, and an HPO id is the same for a
    reader who has one from a lesion facet or a clinical note. `_terms`
    deduplicates, so a record whose label is its identifier — an unmirrored gene
    — still publishes it once.

    The approved name is a term and not the label. `mirrors/genes.tsv` requires
    `name` on every row and leaves `aliases` nullable, and measured on the mirror
    as committed, `name` is non-null in 186/186 rows while `aliases` is null in
    41 (re-measured 2026-11-14; 185/185 and 41 on 2026-10-25, 154/154 and 38
    before the registry widened):
    the one piece of gene text a reader is guaranteed to be able to search for
    is the one that was not indexed, and the one that was is absent for more
    than a fifth of the registry. This paragraph read "1/1 rows" and "null in 1/1"
    until the registry was mirrored from HGNC in full; the earlier figures were
    measured against a two-line file holding TBX5 alone, and the conclusion they
    supported -- that `aliases` "contributes nothing today" -- was an artefact of
    that one row rather than a property of the column. It stays out of `label`
    because `label` is the result row and
    `bundles.py` publishes the symbol as both the browse row's label and the page
    heading; a search hit reading "T-box transcription factor 5" for the gene the
    rest of the site calls TBX5 would read as two different things.

    `genes` comes from `mirrors/genes.tsv`, keyed on HGNC id, and is taken as a
    `Mapping` because nothing here mutates it — the concrete type is invariant,
    so spelling it `dict` would reject legitimate arguments for no benefit.
    `build_genes` takes `symbols: Mapping[str, str]` the same way for the same
    reason, and keeps taking exactly that: it needs the symbol alone, so
    `build_site` passes `{gene: labels.symbol for gene, labels in genes.items()}`
    there and this whole record here.

    `published` is D21's population and is what the gene loop iterates — the
    module docstring says why that and not `corpus.assertions` or `genes`. Taken
    as a `Collection[str]` rather than a `set` for the reason `genes` is a
    `Mapping`, and the same annotation `build_genes` and `build_landing` take it
    under; `build_site` derives it once and hands the same object to all three,
    so the search index, the browse payload and the front page's count cannot
    name three different populations.

    A gene absent from `genes` labels itself with its HGNC id, so a result row
    reads "HGNC:4173" rather than blank. That fallback is reachable behind a
    green gate, which it was not while this loop drew genes from the assertions:
    REF001 makes an *asserted* gene missing from the registry an error, and no
    rule pairs `mirrors/clingen_gene_validity.tsv` with `mirrors/genes.tsv`, so a
    a gene an authority has classified and HGNC has not been mirrored for reaches
    here labelled with its id. Measured 2026-08-06 against the committed mirrors:
    0 of the 92 published genes are missing from the registry. `build_genes` and
    `build_landing` apply the identical `symbols.get(gene, gene)` fallback, so
    such a gene reads "HGNC:4173" consistently across the browse row, the page
    heading and the search hit rather than being unlabelled in one of them.

    Aliases are sorted; authors and synonyms are not. An alias cell is a
    pipe-separated dump whose order is whatever the upstream release emitted, so
    sorting it makes two mirrors listing the same aliases publish identical
    bytes. An author list is data — first author, last author — and a synonym
    list is the curator's, so reordering either would publish something the
    record does not say.

    The two literal filenames below are written by `literature.py`, and this is
    the only place that dependency exists; renaming either payload there means
    renaming it here, and nothing checks the pair. Gene paths carry no such risk
    because `gene_bundle_path` is what `build_genes` writes with.
    """
    records: list[Json] = []

    # Sorted so that the published array does not follow set iteration order,
    # which varies with `PYTHONHASHSEED` between two builds of one commit and
    # would move every checksum downstream of it. `published_genes` returns a
    # `set`, so that hazard is exactly as live as it was when this line built a
    # set from the assertions.
    for gene in sorted(published):
        labels = genes.get(gene, GeneLabels(symbol=gene))
        # Filtered rather than published as a term: a null `name` is a bypassed
        # TBL003 (see `GeneLabels`), and an empty string in `terms` is a term
        # that every "contains" query matches, so one blank cell would make every
        # gene it reached a hit for anything a visitor types.
        name = [labels.name] if labels.name else []
        records.append(
            {
                "kind": "gene",
                "id": gene,
                "label": labels.symbol,
                # Retired symbols last, after the aliases. `_terms` keeps the
                # first occurrence of each string, so a name that is both an
                # alias and a retired name of this gene is published once, in
                # the alias position -- the same precedence `genes.py` resolves
                # by, so the search index and the resolver cannot disagree
                # about which tier a string belongs to.
                "terms": _terms(
                    [
                        labels.symbol,
                        gene,
                        *name,
                        *sorted(labels.aliases),
                        *sorted(labels.prev_symbols),
                    ]
                ),
                # The NewType is recovered rather than asserted, exactly as
                # `bundles.py` recovers it over the same population: `published`
                # is annotated `Collection[str]`, but every member `build_site`
                # puts in it comes from the `gene` column of
                # `mirrors/clingen_gene_validity.tsv`, checked against
                # `HGNC_PATTERN` by that table's schema, and a TBL failure is an
                # error the gate refuses on. `gene_bundle_path` requires `HgncId`
                # so that only HGNC ids can name a file under `genes/`, where a
                # `ContrastId` would slug to a name that collides with one on a
                # case-insensitive filesystem.
                "path": gene_bundle_path(HgncId(gene)),
            }
        )

    # These two follow their corpus order, which is one YAML file's line order
    # and so is already stable between two builds of one commit — unlike the gene
    # set above, this is not what keeps the build byte-identical. Sorted anyway,
    # for two other reasons: it keeps the published order independent of where a
    # curator inserts a record, and it matches the order `publications.json` and
    # `phenotypes.json` already use. It ranks nothing; a search result is ranked
    # by whatever the client matched on rather than by the order it read them in.
    for publication in sorted(corpus.publications, key=lambda item: item.id):
        records.append(
            {
                "kind": "publication",
                "id": publication.id,
                "label": publication.title,
                "terms": _terms(
                    [publication.title, publication.id, publication.journal, *publication.authors]
                ),
                "path": PUBLICATIONS,
            }
        )

    for term in sorted(corpus.phenotypes, key=lambda item: item.id):
        records.append(
            {
                "kind": "phenotype",
                "id": term.id,
                "label": term.label,
                "terms": _terms([term.label, term.id, *term.synonyms]),
                "path": PHENOTYPES,
            }
        )

    # Gzipped because this is the one payload downloaded in full by a visitor who
    # has not asked for anything yet, and it is almost entirely repeated English.
    # Emitted even when empty, so a search box fetching it before the first
    # assertion is curated reads an empty array rather than having to tell "no
    # records" from "wrong URL" by interpreting a 404.
    emitter.write_json_gz("search/index.json.gz", {"records": records})
