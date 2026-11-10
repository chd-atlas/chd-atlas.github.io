"""Apply a verified package to the repository: `mirrors/`, its pin, and the cache.

`ingest_package` is the one writer. It verifies the package (`load_package`),
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

Cells are written as polars parsed them, which Decision 8 calls verbatim: a
package table is read as quoted TSV and a mirror is written unquoted, so a cell
that opens with a quote loses its quoting (measured on polars 1.43.1: ``"x"``
is written ``x``, ``"p""q"`` is written ``p"q``, while ``mid"dle`` is written
as it came). Only the gene column is rewritten, to the HGNC id it resolved to.
A current mirror is read as text, never through `read_table`'s dtypes, so the
rows another source owns are rewritten as the bytes they were (measured over
every committed mirror; see
`test_every_committed_mirror_reads_and_rewrites_byte_for_byte`). Rows are
sorted by the schema's ``sort_key``, each column compared as its schema dtype
with an empty cell smallest: the order `validate.sort_order._precedes` checks.

`ING010` is every refusal about the repository rather than the package: a
current mirror, ``sources.yaml`` or committed package manifest the ingest
cannot read before writing, a ``mirrors/genes.tsv`` that is absent or lists no
gene, and any validation error once the package is applied.
Validation warnings never refuse an ingest (`ValidationReport.ok`); they are
returned for the caller to print.
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable, Iterable
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
from chd_atlas.ingest.package import (
    BLANK,
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
from chd_atlas.tables import GENES, PACKAGES_DIR, TableSchema, cell_text, table_layout
from chd_atlas.validate.runner import ValidationReport, validate_repository
from chd_atlas.validate.sources import PackagePin

# Written whole, never filtered to the registry: spec §4.4 derives the registry
# from these two, so filtering them would discard the rows the next derivation
# needs. Measured 2026-11-10 against the 185-gene `mirrors/genes.tsv`: 3,429 of
# `clingen_gene_validity.tsv`'s 3,653 rows and 29,257 of
# `gencc_submissions.tsv`'s 30,410 rows name a gene outside it, and the
# committed mirrors carry them all. Their gene column is still resolved.
UNFILTERED_TABLES: Final[frozenset[str]] = frozenset({"clingen_validity", "gencc_submissions"})

# The source that owns each flat table without a ``source`` column (ING014).
# `_merge` replaces such a table whole, so any other source's package would
# overwrite the owner's rows. Ids as `mirrors/sources.yaml` records them
# (measured 2026-11-10). `genes`, the third such table, is refused as ING009.
FLAT_TABLE_OWNERS: Final[dict[str, str]] = {
    "clingen_validity": "clingen",
    "gencc_submissions": "gencc",
}


@dataclass(frozen=True)
class IngestResult:
    source_id: str
    version: str
    # rows_in == rows_written + rows_outside_registry + rows_unresolved, which
    # `ingest_package` checks before writing.
    rows_in: int
    rows_written: int
    rows_outside_registry: int
    # Rows whose gene id is one of `unresolved`'s, and so were not written.
    rows_unresolved: int
    unresolved: tuple[Resolved, ...]
    # Registry genes with no resolved row in any gene-keyed table of the
    # package, sorted; empty for a package with no gene-keyed table.
    registry_genes_absent: tuple[str, ...]
    # How many genes `mirrors/genes.tsv` lists.
    registry_genes: int
    coverage: Literal["every_gene"] | RegistrySlice
    dropped: tuple[DroppedRows, ...]
    # The warnings `validate_repository` reported with the package applied; an
    # error would have refused the ingest.
    warnings: tuple[ValidationIssue, ...]


Row = tuple[str, ...]


def ingest_package(root: Path, source_id: str, package_dir: Path) -> IngestResult:
    """Merge ``package_dir`` into ``root`` as source ``source_id``, or refuse and restore."""
    package = load_package(package_dir, source_id)
    _refuse_omitted_or_owned_tables(root, source_id, package)
    resolver = GeneIdResolver.from_root(root, package.manifest.ids)
    registry = frozenset(
        _column(_read_mirror(root / "mirrors" / "genes.tsv", GENES), GENES, "hgnc_id")
    )
    if not registry:
        # Not left to validation: with no registry the filter drops every row,
        # and a report saying so would read as a fact about the package.
        raise IngestRefused(
            "ING010",
            "mirrors/genes.tsv is absent or lists no gene, so the registry filter "
            "would drop every row",
        )

    writes: dict[Path, bytes] = {}
    unresolved: dict[str, Resolved] = {}
    present: set[str] = set()
    gene_keyed = False
    rows_in = rows_written = rows_outside = rows_unresolved = 0
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
        rows = [
            tuple("" if cell is None else str(cell) for cell in row)
            for row in frame.select(schema.column_names).iter_rows()
        ]
        rows_written += len(rows)
        path = _mirror_path(root, table)
        writes[path] = _render(schema, _merge(_read_mirror(path, schema), rows, schema, source_id))
    if rows_in != rows_written + rows_outside + rows_unresolved:
        # Never reached through data: every row is written, counted outside the
        # registry or named as unresolved. A raise, not an assert (`-O`).
        raise RuntimeError(
            f"ingest row accounting does not add up: {rows_in} in, {rows_written} written, "
            f"{rows_outside} outside the registry, {rows_unresolved} unresolved"
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

    report = _apply(root, writes)

    coverage = manifest.coverage
    return IngestResult(
        source_id=source_id,
        version=manifest.version,
        rows_in=rows_in,
        rows_written=rows_written,
        rows_outside_registry=rows_outside,
        rows_unresolved=rows_unresolved,
        unresolved=tuple(unresolved[symbol] for symbol in sorted(unresolved)),
        registry_genes_absent=tuple(sorted(registry - present)) if gene_keyed else (),
        registry_genes=len(registry),
        coverage=coverage if coverage == "every_gene" else coverage.registry_slice,
        dropped=tuple(manifest.dropped),
        warnings=tuple(issue for issue in report.issues if issue.severity is Severity.WARNING),
    )


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


def _apply(root: Path, writes: dict[Path, bytes]) -> ValidationReport:
    """Write every payload, then validate; on any exception restore all of it and re-raise."""
    before = {path: path.read_bytes() if path.is_file() else None for path in writes}
    created = sorted(
        {parent for path in writes for parent in path.parents if not parent.exists()},
        key=lambda directory: len(directory.parts),
        reverse=True,
    )
    try:
        for path, payload in writes.items():
            path.parent.mkdir(parents=True, exist_ok=True)
            write_bytes_atomically(path, payload)
        report = validate_repository(root)
        if not report.ok:
            errors = [issue.format() for issue in report.issues if issue.severity is Severity.ERROR]
            raise IngestRefused(
                "ING010",
                f"the repository fails validation with this package applied "
                f"({len(errors)} error(s)):\n" + "\n".join(errors),
            )
    # BaseException, and always re-raised: a Ctrl-C mid-write must restore too.
    except BaseException:
        _restore(before, created)
        raise
    return report


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
    The one reader of a cache today, `resolve._hgnc_rows`, verifies what it
    reads again: the manifest against the pin, and the one table it reads,
    ``genes``, against the manifest. A package directory that *is* the cache
    rewrites its own bytes, which changes nothing.
    """
    cache = package_cache(root, package.manifest.package, package.manifest.version)
    writes = {cache / "manifest.json": package.manifest_bytes}
    try:
        for table in package.tables:
            writes[cache / table.spec.file] = (package.directory / table.spec.file).read_bytes()
    except OSError as exc:
        raise IngestRefused("ING001", f"a package table cannot be read again: {exc}") from exc
    return writes


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
