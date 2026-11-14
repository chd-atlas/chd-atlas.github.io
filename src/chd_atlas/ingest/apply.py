"""Apply a verified package to the repository: `mirrors/`, its pin, and the cache.

`ingest_package` is the one writer of a single package (`ingest_all` also rewrites
mirrors to prune, below). It verifies the package (`load_package`),
resolves each gene column to HGNC ids, filters to the registry, merges the
package's rows into the current mirrors, commits the manifest and its pin, and
runs `validate_repository`. **If anything fails from the first write onwards,
every file it touched is restored and every directory it created removed**, so
a refused ingest leaves the checkout as it found it. A restore that itself
fails raises `OSError` naming each path it could not put back.

Which rows a package owns (plan Decision 6, as the controller composed it for a
sharded table with a ``source`` column):

- a flat table with a ``source`` column: the rows whose ``source`` is the
  package id. Every other source's rows are kept, byte for byte where the
  mirror's header is in the schema's column order (every committed mirror's
  is); a permuted header is rewritten in schema order, every cell kept;
- a flat table without one: the whole file;
- a sharded table: the shard file ``<shard>.tsv``, and if the table has a
  ``source`` column (``cell_type_profiles``, ``variants``) only that shard's
  rows whose ``source`` is the package id.

Only the tables the manifest names are replaced, and two refusals guard what
that cannot see. Both are checked before anything is written.

- ``ING013``: a package omits a (table, shard) that the source's committed
  manifest (``mirrors/packages/<source_id>.json``, when one exists) carried.
  Keeping those rows would leave rows no pinned manifest lists, and deleting
  them would remove evidence nobody chose to remove, so neither happens.
- ``ING014``: a package writes a table without a ``source`` column that
  another source owns, which `_merge` would replace whole. A flat one is owned
  by the source `FLAT_TABLE_OWNERS` names; a shard is owned by the source whose
  committed manifest lists it. The converter-written shards no manifest lists
  yet, ``profiles/E-MTAB-6814.tsv`` and ``profile_quantiles/E-MTAB-6814.tsv``,
  stay unprotected until phase 3 pins them.

- ``ING016``: a package carries the version the source is already pinned at but
  a manifest whose sha256 differs from the pin's. A version names exactly one
  set of bytes (spec F13), so changed bytes must change the version. The same
  version with the same bytes is allowed, because re-ingesting a pinned package
  from its cache is how a scope change re-filters it: `chd-atlas ingest <id>`
  without ``--package`` does that for one source, `ingest_all` for every one.

`ingest_all` re-ingests every pinned package from its cache as ONE transaction,
in `ingest_order`: the pinned owners of the registry's input tables (ClinGen,
GenCC), then ``hgnc``, whose filter they derive, then every other source by id.
Given ``repin`` (``ingest --all --repin SOURCE=DIR``), it ingests each named
source from that directory instead, in the same order. That is the one route
for a ClinGen or GenCC package, whatever it does to the registry. Ingested
alone it is refused only when it brings a gene into the registry (validation
reports the registry gene `genes.tsv` lacks, GEN006, before ``hgnc`` has
re-derived it), which is why a single ingest's validation failure ends with
`REPIN_HINT`. A package that removes a gene, or adds a disease the pinned Mondo
package carries, is accepted alone and leaves the pruning and `mondo_xrefs`
re-filtering undone (reviewer's measurement on copies of the committed corpus,
2026-11-12: GenCC v2 without HGNC:10727's in-scope rows, 185 -> 184 genes, 0
errors; HGNC:10727, SEMA3E, has had an in-scope ClinGen record since 2026-11-14,
so the same edit no longer removes it).
No step validates. Validating each step would deadlock a scope shrink: ``hgnc``
alone drops a gene that other rows still name (BUR011), and any other package
alone keeps its rows, because its filter is the `genes.tsv` only ``hgnc``
rewrites. After the steps it prunes (the rule is `ingest_all`'s docstring),
then `validate_repository` runs once. On any failure, every file under
`mirrors/` is put back to its bytes before the call and everything the call
created there is removed.

``ING017``: a pinned package is not in its cache, ``local/packages/<id>/<version>/``:
the directory is missing, its ``manifest.json`` cannot be read or is not the
bytes the pin's sha256 names, or those bytes name another package or version
than the pin (a pin SRC005 reports). The message names the `gh release download`
command that fetches it; the atlas never downloads (plan Decision 1).
`ingest_all` checks every pin but a repinned source's before ingesting anything.
So every other step reads a cache directory that holds the pinned manifest and
names the pinned package and version, and its cache write (`_cache_writes`)
rewrites that same directory with the bytes it read, changing nothing under
`local/`. That holds only because of the ING017 checks; `ingest_package` called
directly with another directory writes a new cache directory, as it should, and
so does a repinned step. So before any step `ingest_all` records each file a
repinned step's cache write will touch (`_cache_paths` of the manifest in the
repinned directory), and a failure puts those back with `mirrors/`: a cache
directory the call created is removed, and one that already held the package
(where ``gh release download`` put it) keeps its bytes. Nothing else under
`local/` is snapshotted.

Cells are written as polars parsed them, which Decision 8 calls verbatim: a
package table is read as quoted TSV and a mirror is written unquoted, so a cell
that opens with a quote loses its quoting (measured on polars 1.43.1: ``"x"``
is written ``x``, ``"p""q"`` is written ``p"q``). A quote inside a cell is
kept, never repaired, and what follows depends on the quote count of its row.
Measured on 2026-11-11 (polars 1.43.1) by ingesting a six-row package into a
copy of the corpus, one ``method_note`` cell changed per run: ``m"id"dle`` (an
even count) on the first row, or on the last, is accepted and written
verbatim; ``mid"dle`` (an odd count) on the first or a middle row is refused as
ING011 (``ComputeError``, CSV malformed); ``mid"dle`` on the last row parses,
is written, and `validate` then cannot read the mirror (TBL000), so the ingest
is refused as ING010. Not measured over a large table. Only the gene column is
rewritten, to the HGNC id it resolved to.
A current mirror is read as text, never through `read_table`'s dtypes, so the
rows another source owns are rewritten as the bytes they were (measured over
every committed mirror; see
`test_every_committed_mirror_reads_and_rewrites_byte_for_byte`). Rows are
sorted by the schema's ``sort_key``, each column compared as its schema dtype
with an empty cell smallest: the order `validate.sort_order._precedes` checks.

`ING010` is every refusal about the repository rather than the package: a
current mirror, ``sources.yaml`` or committed package manifest the ingest
cannot read before writing, a registry or table filter the ingest cannot use,
and any validation error once the package is applied. An ``hgnc`` ingest's validation
failure ends with one sentence pointing to ``chd-atlas ingest --all`` (`ALL_HINT`),
and a ClinGen or GenCC ingest's with one naming ``chd-atlas ingest --all --repin``
and the package's directory (`REPIN_HINT`).

The registry filter is ``mirrors/genes.tsv``'s genes for every package but
``hgnc``. The ``hgnc`` package writes ``genes.tsv`` itself, so filtering it to
the file it replaces would keep a stale registry stale: its filter is the
registry `chd_atlas.registry` derives from the scope and the two validity
mirrors. So an ``hgnc`` ingest needs no ``genes.tsv`` at all, and refuses when
the registry cannot be derived (`RegistryUnavailable`, its reason quoted) or
derives empty; any other package refuses when ``genes.tsv`` is absent or lists
no gene. Either empty registry would drop every row.
Two tables are also filtered to the atlas's scope (spec F23), after resolution
and the registry filter; neither has a gene column, so this is their only
filter. ``orphanet_phenotypes`` keeps the rows under HP:0001626 in the pinned
HPO release and ``mondo_xrefs`` the rows whose term a validity mirror names
(the rules, and when a filter cannot be built, are `ingest.filters`'s
docstring). A row the filter does not keep is counted in
`IngestResult.rows_outside_filter`; a row whose term the release does not
define is counted in `IngestResult.rows_unknown_term` instead and its term
named in `IngestResult.unknown_terms`, never counted as outside. A filter that
cannot be built (`FilterUnavailable`) is ING010, before anything is written.

Validation warnings never refuse an ingest (`ValidationReport.ok`); they are
returned for the caller to print: on `IngestResult.warnings`, which is empty for
an ingest run with ``validate=False``; `ingest_all` returns its one
validation's warnings on `IngestAllResult.warnings`.
"""

from __future__ import annotations

import hashlib
import shlex
from collections import Counter
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
from io import StringIO
from pathlib import Path
from typing import Any, Final, Literal

import polars as pl
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.error import CommentMark, YAMLError
from ruamel.yaml.scalarstring import FoldedScalarString, LiteralScalarString
from ruamel.yaml.tokens import CommentToken

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.genes import Resolved
from chd_atlas.ingest.filters import FilterUnavailable, table_filter
from chd_atlas.ingest.package import (
    BLANK,
    HGNC_SOURCE,
    IngestRefused,
    LoadedPackage,
    LoadedTable,
    load_package,
)
from chd_atlas.ingest.resolve import GeneIdResolver, package_cache, resolve_frame
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.package import (
    DroppedRows,
    ManifestError,
    PackageManifest,
    PackageSource,
    RegistrySlice,
    load_manifest,
)
from chd_atlas.registry import SCOPE_FILE, RegistryUnavailable, derive_registry_from_root
from chd_atlas.tables import (
    GENES,
    PACKAGES_DIR,
    TABLE_SCHEMAS,
    TableSchema,
    cell_text,
    mirror_paths,
    table_layout,
)
from chd_atlas.validate.runner import ValidationReport, validate_repository
from chd_atlas.validate.sources import PackagePin, load_sources

# Written whole, never filtered to the registry: spec §4.4 derives the registry
# from these two, so filtering them would discard the rows the next derivation
# needs. Measured 2026-11-14 against the 186-gene `mirrors/genes.tsv`: 3,501 of
# `clingen_gene_validity.tsv`'s 3,733 rows and 29,222 of
# `gencc_submissions.tsv`'s 30,376 rows name a gene outside it (3,429 of 3,653
# and 29,257 of 30,410 against 185 genes on 2026-11-10), and the committed
# mirrors carry them all. Their gene column is still resolved.
UNFILTERED_TABLES: Final[frozenset[str]] = frozenset({"clingen_validity", "gencc_submissions"})

# The source that owns each flat table without a ``source`` column (ING014).
# `_merge` replaces such a table whole, so any other source's package would
# overwrite the owner's rows. Ids as `mirrors/sources.yaml` records them
# (measured 2026-11-10). `genes`, the third such table, belongs to `hgnc`, and
# `load_package` refuses it from any other package (ING009) before this runs.
FLAT_TABLE_OWNERS: Final[dict[str, str]] = {
    "clingen_validity": "clingen",
    "gencc_submissions": "gencc",
}

# The sources whose tables the registry is derived from (`chd_atlas.registry`):
# `ingest_order` puts them before hgnc, and their single ingest's validation
# failure ends with `REPIN_HINT`.
REGISTRY_INPUTS: Final[frozenset[str]] = frozenset(FLAT_TABLE_OWNERS.values())

# Where `chd-atlas/feeds` publishes packages, as GitHub releases tagged
# ``<id>/<version>`` (spec F14). Only ever named in a message (ING017).
FEEDS_REPOSITORY: Final = "chd-atlas/feeds"

# Appended to an hgnc ingest's validation failure (see the module docstring).
ALL_HINT: Final = (
    "If curation/chd_scope.yaml dropped a gene that other mirror rows still name, "
    "`chd-atlas ingest --all` re-ingests every pinned package, removes the rows of genes "
    "that left the registry, and validates once."
)

# Appended to a ClinGen or GenCC ingest's validation failure, formatted with the
# source id and the package directory as given (shell-quoted).
REPIN_HINT: Final = (
    "If this package moved the gene registry, `chd-atlas ingest --all --repin "
    "{source}={directory}` ingests it with every other pinned package, hgnc after it, "
    "and validates once."
)


@dataclass(frozen=True)
class IngestResult:
    source_id: str
    version: str
    # rows_in == rows_written + rows_outside_registry + rows_unresolved
    #   + rows_outside_filter + rows_unknown_term,
    # which `ingest_package` checks before writing.
    rows_in: int
    rows_written: int
    rows_outside_registry: int
    # Rows whose gene id is one of `unresolved`'s, and so were not written.
    rows_unresolved: int
    unresolved: tuple[Resolved, ...]
    # Rows a table filter (`ingest.filters`) did not keep although the
    # filter's release defines their term, and the reason of each filter the
    # package's tables went through, once each in table order (empty when none
    # did, as for every table but orphanet_phenotypes and mondo_xrefs).
    rows_outside_filter: int
    filter_reasons: tuple[str, ...]
    # Rows whose term the filter's release does not define, not written and
    # never counted as outside the filter; those terms, sorted.
    rows_unknown_term: int
    unknown_terms: tuple[str, ...]
    # Registry genes with no resolved row in any gene-keyed table of the
    # package, sorted; empty for a package with no gene-keyed table.
    registry_genes_absent: tuple[str, ...]
    # How many genes the registry filter held: those `mirrors/genes.tsv` lists,
    # or for `hgnc` those the scope and the validity mirrors derive.
    registry_genes: int
    coverage: Literal["every_gene"] | RegistrySlice
    dropped: tuple[DroppedRows, ...]
    # The warnings `validate_repository` reported with the package applied; an
    # error would have refused the ingest. Empty when the ingest did not
    # validate (``validate=False``), as no step of `ingest_all` does.
    warnings: tuple[ValidationIssue, ...]


@dataclass(frozen=True)
class PrunedRows:
    """Rows `ingest_all` removed from one mirror file because their gene left the registry."""

    table: str
    path: str  # relative to the repository root, forward-slashed
    gene: str
    rows: int


@dataclass(frozen=True)
class IngestAllResult:
    # One per pinned or repinned source, in `ingest_order`; each step's
    # `warnings` is empty.
    results: tuple[IngestResult, ...]
    # Sorted by (path, gene); empty unless the hgnc step removed genes.
    pruned: tuple[PrunedRows, ...]
    # The one validation's warnings; an error would have refused the call.
    warnings: tuple[ValidationIssue, ...]


Row = tuple[str, ...]


def ingest_package(
    root: Path, source_id: str, package_dir: Path, *, validate: bool = True
) -> IngestResult:
    """Merge ``package_dir`` into ``root`` as source ``source_id``, or refuse and restore.

    With ``validate=False`` the repository is not validated and the result has
    no warnings: `ingest_all` validates once, after every step.
    """
    package = load_package(package_dir, source_id)
    _refuse_omitted_or_owned_tables(root, source_id, package)
    resolver = GeneIdResolver.from_root(root, package.manifest.ids)
    registry = _registry(root, source_id)

    _refuse_changed_bytes_under_a_pinned_version(root, source_id, package)

    writes: dict[Path, bytes] = {}
    unresolved: dict[str, Resolved] = {}
    present: set[str] = set()
    gene_keyed = False
    rows_in = rows_written = rows_outside = rows_unresolved = 0
    rows_outside_filter = rows_unknown_term = 0
    filter_reasons: list[str] = []
    unknown_terms: set[str] = set()
    for table in package.tables:
        schema = table.schema
        frame = table.frame
        rows_in += frame.height
        gene_column = schema.gene_column
        if gene_column is not None:
            gene_keyed = True
            raw_ids = frame[gene_column].to_list()
            frame, missed = resolve_frame(frame, gene_column, resolver)
            unresolved.update((item.symbol, item) for item in missed)
            # Counted by the names reported, not by the rows `resolve_frame`
            # kept, so a row dropped without being named breaks the identity
            # checked below rather than passing as unresolved.
            named = {item.symbol for item in missed}
            rows_unresolved += sum((cell_text(raw) or BLANK) in named for raw in raw_ids)
            present.update(str(gene) for gene in frame[gene_column].to_list())
            if schema.name not in UNFILTERED_TABLES:
                kept = frame.filter(pl.col(gene_column).is_in(sorted(registry)))
                rows_outside += frame.height - kept.height
                frame = kept
        try:
            row_filter = table_filter(root, schema.name)
        except FilterUnavailable as exc:
            raise IngestRefused(
                "ING010", f"the {schema.name} filter cannot be built: {exc}"
            ) from exc
        if row_filter is not None:
            filtered = row_filter.apply(frame)
            frame = filtered.frame
            rows_outside_filter += filtered.outside
            rows_unknown_term += filtered.unknown_rows
            unknown_terms.update(filtered.unknown_terms)
            if row_filter.reason not in filter_reasons:
                filter_reasons.append(row_filter.reason)
        rows = [
            tuple("" if cell is None else str(cell) for cell in row)
            for row in frame.select(schema.column_names).iter_rows()
        ]
        rows_written += len(rows)
        path = _mirror_path(root, table)
        writes[path] = _render(schema, _merge(_read_mirror(path, schema), rows, schema, source_id))
    not_written = rows_outside + rows_unresolved + rows_outside_filter + rows_unknown_term
    if rows_in != rows_written + not_written:
        # Never reached through data: every row is written, counted outside the
        # registry or a table filter, or named as unresolved or as an unknown
        # term. A raise, not an assert (`-O`).
        raise RuntimeError(
            f"ingest row accounting does not add up: {rows_in} in, {rows_written} written, "
            f"{rows_outside} outside the registry, {rows_unresolved} unresolved, "
            f"{rows_outside_filter} outside the filter, {rows_unknown_term} naming an "
            "unknown term"
        )

    manifest = package.manifest
    digest = hashlib.sha256(package.manifest_bytes).hexdigest()
    pin = PackagePin(version=manifest.version, manifest_sha256=digest)
    writes[root / "mirrors" / PACKAGES_DIR / f"{source_id}.json"] = package.manifest_bytes
    sources_path = root / "mirrors" / "sources.yaml"
    try:
        sources = sources_path.read_bytes()
    except OSError as exc:
        raise IngestRefused("ING010", f"mirrors/sources.yaml cannot be read: {exc}") from exc
    writes[sources_path] = with_source_entry(sources, source_id, manifest.source, pin)
    writes.update(_cache_writes(root, package))

    warnings = _apply(root, writes, _hint(source_id, package_dir), validate=validate)

    coverage = manifest.coverage
    return IngestResult(
        source_id=source_id,
        version=manifest.version,
        rows_in=rows_in,
        rows_written=rows_written,
        rows_outside_registry=rows_outside,
        rows_unresolved=rows_unresolved,
        unresolved=tuple(unresolved[symbol] for symbol in sorted(unresolved)),
        rows_outside_filter=rows_outside_filter,
        filter_reasons=tuple(filter_reasons),
        rows_unknown_term=rows_unknown_term,
        unknown_terms=tuple(sorted(unknown_terms)),
        registry_genes_absent=tuple(sorted(registry - present)) if gene_keyed else (),
        registry_genes=len(registry),
        coverage=coverage if coverage == "every_gene" else coverage.registry_slice,
        dropped=tuple(manifest.dropped),
        warnings=warnings,
    )


def ingest_order(pinned: Iterable[str]) -> tuple[str, ...]:
    """The order `ingest_all` applies pinned and repinned sources in (see the module docstring).

    ``hgnc``'s filter is derived from the ClinGen and GenCC mirrors, so a
    re-pinned one must land before it; a gene-keyed package's filter is the
    `genes.tsv` that ``hgnc`` writes. Two others filter on something else:
    ``orphanet_product4`` on the pinned HPO release's cardiovascular branch and
    ``mondo_sssom`` on the terms the validity mirrors name.
    """
    sources = set(pinned)
    inputs = sorted(sources & REGISTRY_INPUTS)
    hgnc = [HGNC_SOURCE] if HGNC_SOURCE in sources else []
    rest = sorted(sources - set(inputs) - {HGNC_SOURCE})
    return (*inputs, *hgnc, *rest)


def package_pins(root: Path) -> dict[str, PackagePin]:
    """Each pinned source's pin in ``mirrors/sources.yaml``; ING010 if the file cannot be read."""
    registry, issues = load_sources(root)
    if any(issue.severity is Severity.ERROR for issue in issues):
        raise IngestRefused(
            "ING010",
            "mirrors/sources.yaml cannot be read: " + "; ".join(i.message for i in issues),
        )
    return {source.id: source.package for source in registry.sources if source.package}


def pinned_package_dir(root: Path, source_id: str, pin: PackagePin) -> Path:
    """The cache directory holding the pinned package; ING017 when it does not."""
    cache = package_cache(root, source_id, pin.version)
    where = cache.relative_to(root).as_posix()
    download = f"gh release download {source_id}/{pin.version} -R {FEEDS_REPOSITORY} -D {where}"
    if not cache.is_dir():
        raise IngestRefused(
            "ING017",
            f"{source_id} {pin.version} is pinned in mirrors/sources.yaml but not cached at "
            f"{where}; download it: {download}",
        )
    try:
        data = (cache / "manifest.json").read_bytes()
    except OSError as exc:
        raise IngestRefused(
            "ING017",
            f"{where}/manifest.json cannot be read ({exc.strerror or exc}); remove the "
            f"directory and download it again: {download}",
        ) from exc
    if hashlib.sha256(data).hexdigest() != pin.manifest_sha256:
        raise IngestRefused(
            "ING017",
            f"{where}/manifest.json is not the manifest {source_id} {pin.version} is pinned at "
            f"(sha256:{pin.manifest_sha256[:12]}\u2026); remove the directory and download it "
            f"again: {download}",
        )
    try:
        manifest = load_manifest(cache / "manifest.json")
    except ManifestError:
        # The pinned bytes are not a manifest: `load_package` refuses them
        # (ING004) before that step writes anything.
        return cache
    if (manifest.package, manifest.version) != (source_id, pin.version):
        raise IngestRefused(
            "ING017",
            f"{where}/manifest.json has the pinned bytes but is {manifest.package} "
            f"{manifest.version}, so the pin disagrees with its own manifest (validate reports "
            f"it as SRC005); ingesting it would move the pin",
        )
    return cache


def ingest_all(root: Path, repin: Mapping[str, Path] | None = None) -> IngestAllResult:
    """Re-ingest every pinned package from its cache as one transaction, then validate once.

    ``repin`` maps a source id to a package directory. That source's step
    reads the directory instead of the pinned cache, which is therefore not
    checked (ING017); its cache and pin are written as a single ingest writes
    them, and ING016 applies to it as to every step, since every step is
    `ingest_package`. A repinned source need not be pinned or even listed in
    ``mirrors/sources.yaml``: its step then appends the entry, as a single
    ingest with ``--package`` does.

    The steps run in `ingest_order`, each with ``validate=False``. Then, **only
    when the ``hgnc`` step ran in this call**, the call prunes: the genes it
    removed are the HGNC ids `mirrors/genes.tsv` held before the call and does
    not hold after it, and every row naming one of them is removed from every
    mirror file (each flat file, each shard) of every table with a gene column
    except ``genes`` and `UNFILTERED_TABLES`, whichever package or converter
    wrote the row. That is the atlas applying its own registry (spec F4), and
    what re-running each converter against the shrunk registry would do
    (reasoned by reading the five converters, not measured). The files are read
    and rewritten as text (`_read_mirror`, `_render`), so every other row keeps
    its bytes provided the header is in schema order (as `_merge` documents).
    Each removal is returned in `IngestAllResult.pruned`, never silent. Then
    `validate_repository` runs once, also when nothing is pinned; an error
    refuses the call (ING010). On any exception everything under
    `mirrors/` is restored, pruned rows included, and so is every file a
    repinned step's cache write touched (see the module docstring).
    `curation/` is never edited: a curated record naming a removed gene is left
    for validation to report. A single `ingest_package` never prunes.
    """
    repinned = dict(repin or {})
    pins = package_pins(root)
    # Every pinned cache a step will read, before any step, so a missing one
    # refuses with nothing written.
    packages = {
        source_id: (
            repinned[source_id]
            if source_id in repinned
            else pinned_package_dir(root, source_id, pins[source_id])
        )
        for source_id in ingest_order({*pins, *repinned})
    }
    cache_before, cache_created = _journal(
        path for directory in repinned.values() for path in _repinned_cache_paths(root, directory)
    )
    hgnc_runs = HGNC_SOURCE in packages
    genes_before = _registry_ids(root) if hgnc_runs else frozenset()
    mirrors = root / "mirrors"
    files, directories = _snapshot(mirrors)
    try:
        results = tuple(
            ingest_package(root, source_id, directory, validate=False)
            for source_id, directory in packages.items()
        )
        pruned = _prune(root, genes_before - _registry_ids(root)) if hgnc_runs else ()
        report = validate_repository(root)
        _refuse_invalid(
            report,
            "with every pinned package re-ingested" if packages else "with no package pinned",
        )
    # BaseException, and always re-raised, as in `_apply`.
    except BaseException:
        before, created = _tree_journal(mirrors, files, directories)
        _restore(
            {**before, **cache_before},
            sorted([*created, *cache_created], key=lambda path: len(path.parts), reverse=True),
        )
        raise
    return IngestAllResult(results, pruned, _warnings(report))


def _repinned_cache_paths(root: Path, directory: Path) -> list[Path]:
    """`_cache_paths` of the package in ``directory``; none if its manifest cannot be read.

    A manifest that cannot be read is refused (ING004) before its step writes.
    """
    try:
        return _cache_paths(root, load_manifest(directory / "manifest.json"))
    except ManifestError:
        return []


def _registry_ids(root: Path) -> frozenset[str]:
    """The HGNC ids `mirrors/genes.tsv` lists (none when it is absent)."""
    return frozenset(_column(_read_mirror(root / "mirrors" / "genes.tsv", GENES), GENES, "hgnc_id"))


def _prune(root: Path, genes: frozenset[str]) -> tuple[PrunedRows, ...]:
    """Remove every row naming one of ``genes`` from the registry-filtered mirrors."""
    if not genes:
        return ()
    pruned: list[PrunedRows] = []
    for path, name in mirror_paths(root):
        schema = TABLE_SCHEMAS[name]
        column = schema.gene_column
        if column is None or name == GENES.name or name in UNFILTERED_TABLES:
            continue
        rows = _read_mirror(path, schema)
        index = schema.column_names.index(column)
        removed = Counter(row[index] for row in rows if row[index] in genes)
        if not removed:
            continue
        kept = [row for row in rows if row[index] not in genes]
        write_bytes_atomically(path, _render(schema, kept))
        relative = path.relative_to(root).as_posix()
        pruned += [PrunedRows(name, relative, gene, count) for gene, count in removed.items()]
    return tuple(sorted(pruned, key=lambda item: (item.path, item.gene)))


def _snapshot(directory: Path) -> tuple[dict[Path, bytes], frozenset[Path]]:
    """Every file's bytes and every directory under ``directory``."""
    files: dict[Path, bytes] = {}
    directories: set[Path] = set()
    for path in directory.rglob("*"):
        if path.is_file():
            files[path] = path.read_bytes()
        elif path.is_dir():
            directories.add(path)
    return files, frozenset(directories)


def _tree_journal(
    directory: Path, files: dict[Path, bytes], directories: frozenset[Path]
) -> tuple[dict[Path, bytes | None], list[Path]]:
    """What `_restore` needs to put ``directory`` back to a `_snapshot`.

    Every file's bytes in the snapshot (None for a file created since), and
    every directory created since, deepest first.
    """
    now = list(directory.rglob("*"))
    before: dict[Path, bytes | None] = {
        path: None for path in now if path.is_file() and path not in files
    }
    before.update(files)
    created = sorted(
        (path for path in now if path.is_dir() and path not in directories),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    return before, created


def _journal(paths: Iterable[Path]) -> tuple[dict[Path, bytes | None], list[Path]]:
    """Each path's bytes now (None when absent) and the directories writing it would create.

    The directories are deepest first, the order `_restore` removes them in.
    """
    before = {path: path.read_bytes() if path.is_file() else None for path in paths}
    created = sorted(
        {parent for path in before for parent in path.parents if not parent.exists()},
        key=lambda directory: len(directory.parts),
        reverse=True,
    )
    return before, created


def _warnings(report: ValidationReport) -> tuple[ValidationIssue, ...]:
    return tuple(issue for issue in report.issues if issue.severity is Severity.WARNING)


def _refuse_invalid(report: ValidationReport, applied: str, hint: str | None = None) -> None:
    """ING010 naming every validation error, and ``hint`` on its own line after them."""
    if report.ok:
        return
    errors = [issue.format() for issue in report.issues if issue.severity is Severity.ERROR]
    raise IngestRefused(
        "ING010",
        f"the repository fails validation {applied} ({len(errors)} error(s)):\n"
        + "\n".join(errors)
        + ("" if hint is None else f"\n{hint}"),
    )


def _registry(root: Path, source_id: str) -> frozenset[str]:
    """The genes the filter keeps (see the module docstring); refuses an empty one.

    Not left to validation: with no registry the filter drops every row, and a
    report saying so would read as a fact about the package.
    """
    if source_id == HGNC_SOURCE:
        try:
            derived = derive_registry_from_root(root)
        except RegistryUnavailable as exc:
            raise IngestRefused("ING010", f"the gene registry cannot be derived: {exc}") from exc
        if not derived:
            raise IngestRefused(
                "ING010",
                "the derived registry is empty: no validity mirror row names a gene against "
                f"a disease in {SCOPE_FILE}, so the filter would drop every row",
            )
        return derived
    current = frozenset(
        _column(_read_mirror(root / "mirrors" / "genes.tsv", GENES), GENES, "hgnc_id")
    )
    if not current:
        raise IngestRefused(
            "ING010",
            "mirrors/genes.tsv is absent or lists no gene, so the registry filter "
            "would drop every row",
        )
    return current


def _label(name: str, shard: str | None) -> str:
    return name if shard is None else f"{name} shard {shard}"


def _committed_manifests(root: Path) -> dict[str, PackageManifest]:
    """Every committed package manifest, by the source id its file is named for."""
    found: dict[str, PackageManifest] = {}
    for path in sorted((root / "mirrors" / PACKAGES_DIR).glob("*.json")):
        try:
            found[path.stem] = load_manifest(path)
        except ManifestError as exc:
            raise IngestRefused(
                "ING010",
                f"a committed package manifest cannot be read, so which tables it "
                f"carried is unknown: {exc}",
            ) from exc
    return found


def _refuse_changed_bytes_under_a_pinned_version(
    root: Path, source_id: str, package: LoadedPackage
) -> None:
    """ING016 (see the module docstring); reads, never writes."""
    pin = package_pins(root).get(source_id)
    if pin is None or pin.version != package.manifest.version:
        return
    digest = hashlib.sha256(package.manifest_bytes).hexdigest()
    if pin.manifest_sha256 != digest:
        raise IngestRefused(
            "ING016",
            f"{source_id} {pin.version} is pinned at sha256:{pin.manifest_sha256[:12]}\u2026; "
            f"this package's manifest is sha256:{digest[:12]}\u2026; "
            f"a package that changes its bytes must change its version",
        )


def _refuse_omitted_or_owned_tables(root: Path, source_id: str, package: LoadedPackage) -> None:
    """ING013 and ING014 (see the module docstring); reads, never writes."""
    committed = _committed_manifests(root)
    carried = {(table.spec.name, table.spec.shard) for table in package.tables}
    previous = committed.get(source_id)
    if previous is not None:
        omitted = sorted(
            {(table.name, table.shard) for table in previous.tables} - carried,
            key=lambda pair: (pair[0], pair[1] or ""),
        )
        if omitted:
            raise IngestRefused(
                "ING013",
                f"{source_id} {package.manifest.version} omits what its committed manifest "
                f"({previous.version}) carried: {', '.join(_label(*pair) for pair in omitted)}; "
                f"those rows would be left under no pinned manifest",
            )
    for table in package.tables:
        name, shard = table.spec.name, table.spec.shard
        if table.schema.source_column is not None:
            continue
        if table_layout(name)[0] == "flat":
            owner = FLAT_TABLE_OWNERS.get(name)
            if owner is not None and owner != source_id:
                raise IngestRefused(
                    "ING014",
                    f"table '{name}' has no source column and belongs to source '{owner}'; "
                    f"package '{source_id}' would replace all of it",
                )
            continue
        claimants = sorted(
            other
            for other, manifest in committed.items()
            if other != source_id
            and any((t.name, t.shard) == (name, shard) for t in manifest.tables)
        )
        if claimants:
            raise IngestRefused(
                "ING014",
                f"{_label(name, shard)} has no source column and the committed manifest of "
                f"{', '.join(claimants)} carries it; package '{source_id}' would replace it",
            )


def _hint(source_id: str, package_dir: Path) -> str | None:
    """The sentence a single ingest's validation failure ends with (see the module docstring)."""
    if source_id == HGNC_SOURCE:
        return ALL_HINT
    if source_id in REGISTRY_INPUTS:
        return REPIN_HINT.format(source=source_id, directory=shlex.quote(str(package_dir)))
    return None


def _apply(
    root: Path, writes: dict[Path, bytes], hint: str | None, *, validate: bool
) -> tuple[ValidationIssue, ...]:
    """Write every payload, then validate; on any exception restore all of it and re-raise.

    A validation failure ends with ``hint``. Returns the validation's warnings,
    none when ``validate`` is false.
    """
    before, created = _journal(writes)
    try:
        for path, payload in writes.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            write_bytes_atomically(path, payload)
        if not validate:
            return ()
        report = validate_repository(root)
        _refuse_invalid(report, "with this package applied", hint)
    # BaseException, and always re-raised: a Ctrl-C mid-write must restore too.
    except BaseException:
        _restore(before, created)
        raise
    return _warnings(report)


def _restore(before: dict[Path, bytes | None], created: list[Path]) -> None:
    """Put back every path that differs from its snapshot; remove created directories.

    A path already as it was is left alone. A write that failed (a full disk, a
    read-only ``sources.yaml``) left its target untouched, and rewriting it
    would fail the same way: the restore would then report a file nobody
    changed as one it could not put back, and the real error would be buried
    under that one. So "could not restore" names only paths that still differ.
    """
    failures: list[str] = []
    for path, data in before.items():
        try:
            if (path.read_bytes() if path.is_file() else None) == data:
                continue
            if data is None:
                path.unlink()
            else:
                write_bytes_atomically(path, data)
        except OSError as exc:
            failures.append(f"{path}: {exc}")
    for directory in created:
        try:
            if directory.is_dir():
                directory.rmdir()
        except OSError as exc:
            failures.append(f"{directory}: {exc}")
    if failures:
        raise OSError("the ingest failed and could not restore: " + "; ".join(failures))


def _mirror_path(root: Path, table: LoadedTable) -> Path:
    layout, where = table_layout(table.schema.name)
    if layout == "flat":
        return root / "mirrors" / where
    return root / "mirrors" / where / f"{table.spec.shard}.tsv"


def _cache_writes(root: Path, package: LoadedPackage) -> dict[Path, bytes]:
    """The manifest and the table files it names, keyed by their cache path.

    Nothing else in the package directory is cached. The manifest is the bytes
    `load_package` verified; the tables are read once more and not re-hashed.
    A cache is read by `resolve._hgnc_rows`, `pinned_package_dir` and every
    re-ingest from it. That is safe by reasoning from the code, not by a
    measurement of this function: `load_package` re-hashes every table against
    the manifest, and ING017 and ING016 check the manifest against the pin. A
    package directory that *is* the cache rewrites its own bytes, which changes
    nothing.
    """
    manifest_path, *table_paths = _cache_paths(root, package.manifest)
    writes = {manifest_path: package.manifest_bytes}
    try:
        for path, table in zip(table_paths, package.tables, strict=True):
            writes[path] = (package.directory / table.spec.file).read_bytes()
    except OSError as exc:
        raise IngestRefused("ING001", f"a package table cannot be read again: {exc}") from exc
    return writes


def _cache_paths(root: Path, manifest: PackageManifest) -> list[Path]:
    """Where `_cache_writes` puts a package: its manifest, then each table it lists, in order."""
    cache = package_cache(root, manifest.package, manifest.version)
    return [cache / "manifest.json", *(cache / table.file for table in manifest.tables)]


def _read_mirror(path: Path, schema: TableSchema) -> list[Row]:
    """A current mirror's rows as text, in the schema's column order; none if absent.

    Split on tabs and newlines exactly, so writing a row back reproduces its
    bytes. A file the ingest cannot read, or whose columns are not the
    schema's, is refused rather than rewritten: rewriting would drop a column.
    """
    if not path.exists():
        return []
    try:
        text = path.read_bytes().decode("utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        raise IngestRefused("ING010", f"{path} cannot be read: {exc}") from exc
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    if not lines:
        raise IngestRefused("ING010", f"{path} is empty; a mirror starts with its header")
    header = lines[0].split("\t")
    if sorted(header) != sorted(schema.column_names):
        raise IngestRefused(
            "ING010",
            f"{path} has columns {header}; "
            f"the {schema.name} schema has {list(schema.column_names)}",
        )
    order = [header.index(name) for name in schema.column_names]
    rows: list[Row] = []
    for number, line in enumerate(lines[1:], start=2):
        cells = line.split("\t")
        if len(cells) != len(header):
            raise IngestRefused(
                "ING010",
                f"{path} row {number} has {len(cells)} fields; the header has {len(header)}",
            )
        rows.append(tuple(cells[index] for index in order))
    return rows


def _column(rows: Iterable[Row], schema: TableSchema, name: str) -> list[str]:
    index = schema.column_names.index(name)
    return [row[index] for row in rows]


def _merge(
    current: list[Row], package: list[Row], schema: TableSchema, source_id: str
) -> list[Row]:
    """The current rows the package does not own, then the package's rows (unsorted)."""
    if schema.source_column is None:
        return list(package)
    index = schema.column_names.index(schema.source_column)
    return [row for row in current if row[index] != source_id] + list(package)


def _render(schema: TableSchema, rows: list[Row]) -> bytes:
    """Header and rows sorted by the schema's key, tab-separated, newline-terminated."""
    ordered = sorted(rows, key=_sort_key(schema))
    lines = ["\t".join(schema.column_names), *("\t".join(row) for row in ordered)]
    return ("\n".join(lines) + "\n").encode("utf-8")


def _sort_key(schema: TableSchema) -> Callable[[Row], tuple[tuple[Any, ...], ...]]:
    """Each key column as its schema dtype; an empty cell sorts smallest, as `_precedes`.

    A cell that does not parse as its dtype sorts as empty, and validation then
    refuses the file. Measured: a burden ``maf_max`` of ``abc`` is reported as
    TBL000 "could not parse `abc` as dtype `f64` at column 'maf_max'", which
    names the cell as a raise here would not.
    Python's sort is stable, so rows with equal keys keep their order (TBL007
    refuses such a pair anyway).
    """
    dtypes = {column.name: column.dtype for column in schema.columns}
    converters: list[tuple[int, Callable[[str], Any]]] = []
    for name in schema.sort_key:
        dtype = dtypes[name]
        if dtype is pl.String:
            convert: Callable[[str], Any] = str
        elif dtype is pl.Int64:
            convert = int
        elif dtype is pl.Float64:
            convert = float
        else:
            raise TypeError(f"{schema.name}.{name} is {dtype}, which no sort order is defined for")
        converters.append((schema.column_names.index(name), convert))

    def key(row: Row) -> tuple[tuple[Any, ...], ...]:
        parts: list[tuple[Any, ...]] = []
        for index, convert in converters:
            cell = row[index]
            try:
                parts.append((0,) if cell == "" else (1, convert(cell)))
            except ValueError:
                parts.append((0,))
        return tuple(parts)

    return key


# --- mirrors/sources.yaml ---------------------------------------------------


def _yaml() -> YAML:
    """Round-trip settings that reproduce the committed file byte for byte.

    Measured by the controller 2026-11-10: load and dump of mirrors/sources.yaml
    (14,583 bytes) is identical with these settings; with ruamel's default
    width it is not, because long comment-free lines are re-wrapped.
    """
    yaml = YAML(typ="rt")
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.width = 4096
    return yaml


def with_source_entry(raw: bytes, source_id: str, source: PackageSource, pin: PackagePin) -> bytes:
    """``sources.yaml`` with the entry for ``source_id`` set from a manifest, or appended.

    An existing entry is updated in place: a field whose value is unchanged is
    not touched, ``file_version``/``file_date`` are deleted when the manifest
    gives none, and ``package`` is set. Other entries and every comment are
    left where they are (see `_add_key` and `_delete_key` for the two edits
    that would otherwise move one).
    """
    try:
        data = _yaml().load(raw.decode("utf-8"))
    except (YAMLError, UnicodeDecodeError, ValueError) as exc:
        raise IngestRefused("ING010", f"mirrors/sources.yaml cannot be parsed: {exc}") from exc
    entries = data.get("sources") if isinstance(data, CommentedMap) else None
    if not isinstance(entries, CommentedSeq):
        raise IngestRefused("ING010", "mirrors/sources.yaml has no `sources:` list")

    # Every `PackageSource` field, in its declaration order, so a field added to
    # `Source` (and, by the drift test, to `PackageSource`) is written on the
    # next ingest. An enum is written as its value, a date as a YAML date.
    fields: dict[str, object] = {
        name: value.value if isinstance(value, Enum) else value
        for name, value in source.model_dump().items()
    }
    entry = next(
        (
            item
            for item in entries
            if isinstance(item, CommentedMap) and item.get("id") == source_id
        ),
        None,
    )
    if entry is None:
        entry = CommentedMap([("id", source_id)])
        entry.update((key, value) for key, value in fields.items() if value is not None)
        entry["package"] = CommentedMap(
            [("version", pin.version), ("manifest_sha256", pin.manifest_sha256)]
        )
        entries.append(entry)
    else:
        for key, value in fields.items():
            if value is None:
                if key in entry:
                    _delete_key(entry, key)
            else:
                _set(entry, key, value)
        current = entry.get("package")
        if isinstance(current, CommentedMap):
            _set(current, "version", pin.version)
            _set(current, "manifest_sha256", pin.manifest_sha256)
        else:
            _set(
                entry,
                "package",
                CommentedMap([("version", pin.version), ("manifest_sha256", pin.manifest_sha256)]),
            )
    buffer = StringIO()
    _yaml().dump(data, buffer)
    return buffer.getvalue().encode("utf-8")


def _set(mapping: CommentedMap, key: str, value: object) -> None:
    if key not in mapping:
        _add_key(mapping, key, value)
    elif mapping[key] != value:
        mapping[key] = value


# ruamel attaches the comment lines that follow a value to that value's key
# (measured on ruamel.yaml 0.19.1): a comment above the next entry belongs to
# this entry's last key, and one between two keys to the first of them. So a
# key appended after the last one would land *below* a comment written above
# the next entry, and deleting a key would delete the comment that follows it.
# These two helpers move that block of comment lines to where it still reads
# in the same place; a comment on the end of a deleted key's own line goes
# with the line.


def _add_key(mapping: CommentedMap, key: str, value: object) -> None:
    last = next(reversed(mapping))
    block = _take_comment_block(*_deepest(mapping, last))
    mapping[key] = value
    if block:
        _put_comment_block(*_deepest(mapping, key), block)


def _delete_key(mapping: CommentedMap, key: str) -> None:
    keys = list(mapping)
    # `keys[index - 1]` wraps to the last key when the deleted key opens the
    # entry, keeping its comment (at the entry's end) rather than losing it.
    # Every committed entry opens with `id`, which is never deleted.
    previous = keys[keys.index(key) - 1]
    block = _take_comment_block(*_deepest(mapping, key))
    del mapping[key]
    if block:
        _put_comment_block(*_deepest(mapping, previous), block)


def _deepest(mapping: CommentedMap, key: str) -> tuple[CommentedMap, str]:
    """Where the comments after ``mapping[key]`` live: a nested map's last key, else here."""
    value = mapping[key]
    while isinstance(value, CommentedMap) and value:
        mapping, key = value, next(reversed(value))
        value = mapping[key]
    return mapping, key


def _is_block_scalar(value: object) -> bool:
    # A folded or literal scalar ends its own line, so the comment token after
    # it holds whole comment lines with no end-of-line part.
    return isinstance(value, FoldedScalarString | LiteralScalarString)


def _take_comment_block(mapping: CommentedMap, key: str) -> str | None:
    """Detach the whole comment lines after ``mapping[key]``; keep its end-of-line comment."""
    slot = mapping.ca.items.get(key)
    if not slot or len(slot) < 3 or slot[2] is None:
        return None
    token = slot[2]
    if _is_block_scalar(mapping[key]):
        end_of_line, block = "", token.value
    else:
        head, newline, block = token.value.partition("\n")
        end_of_line = head + newline
    if end_of_line.strip():
        token.value = end_of_line
    else:
        slot[2] = None
    return str(block) or None


def _put_comment_block(mapping: CommentedMap, key: str, block: str) -> None:
    slot = mapping.ca.items.setdefault(key, [None, None, None, None])
    if slot[2] is None:
        prefix = "" if _is_block_scalar(mapping[key]) else "\n"
        slot[2] = CommentToken(prefix + block, CommentMark(0), None)
    else:
        slot[2].value += block
