"""Load a source package from disk and verify it before anything reads its rows.

A package is a directory holding ``manifest.json`` and the table files it lists.
`load_package` is the gate: every way a package can be wrong is its own
`IngestRefused` code, so a curator reads which promise broke. Nothing but
`IngestRefused` escapes it; the unreadable-body cases (ING011) are caught by
exception type, never by a bare ``except``.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

from chd_atlas.models.package import (
    ManifestError,
    PackageManifest,
    PackageTable,
    load_manifest,
)
from chd_atlas.tables import TABLE_SCHEMAS, TableSchema, cell_text, table_layout
from chd_atlas.vocab import GeneIdKind

# How a blank or whitespace-only cell is named in a refusal or a report.
BLANK: Final = "<blank>"

# The one source whose package ships `genes`, the registry table (ING009).
HGNC_SOURCE: Final = "hgnc"


class IngestRefused(Exception):
    """A package that must not be ingested; ``code`` says which promise broke."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code
        self.message = message


@dataclass(frozen=True)
class LoadedTable:
    spec: PackageTable
    schema: TableSchema
    frame: pl.DataFrame  # every column `pl.String`, as read


@dataclass(frozen=True)
class LoadedPackage:
    directory: Path
    manifest: PackageManifest
    manifest_bytes: bytes
    tables: tuple[LoadedTable, ...]


def load_package(directory: Path, source_id: str) -> LoadedPackage:
    """Verify ``directory`` as the package for ``source_id``; refuse, never partly load."""
    manifest_path = directory / "manifest.json"
    try:
        manifest_bytes = manifest_path.read_bytes()
        manifest = load_manifest(manifest_path)
    except (OSError, ManifestError) as exc:
        raise IngestRefused("ING004", f"manifest.json is missing or invalid: {exc}") from exc
    if manifest.package != source_id:
        raise IngestRefused(
            "ING006", f"the package is '{manifest.package}', not source '{source_id}'"
        )
    _refuse_a_misplaced_registry(manifest)
    tables = tuple(_load_table(directory, spec, source_id) for spec in manifest.tables)
    return LoadedPackage(directory, manifest, manifest_bytes, tables)


def _refuse_a_misplaced_registry(manifest: PackageManifest) -> None:
    """ING009: `genes` ships in the ``hgnc`` package alone, which always ships it by HGNC id.

    `genes` has no ``source`` column, so the package that writes it replaces
    the whole registry (`apply._merge`); only HGNC's may. Checked before any
    table is read, so the message names the rule rather than a later symptom.
    """
    ships_genes = any(table.name == "genes" for table in manifest.tables)
    if manifest.package != HGNC_SOURCE:
        if ships_genes:
            raise IngestRefused(
                "ING009",
                f"only the hgnc package may ship the genes table; "
                f"package '{manifest.package}' does",
            )
        return
    if not ships_genes:
        raise IngestRefused(
            "ING009", "the hgnc package must ship the genes table; this one does not"
        )
    if manifest.ids is not GeneIdKind.HGNC_ID:
        raise IngestRefused(
            "ING009",
            f"the hgnc package's gene ids must be hgnc_id; this one declares "
            f"'{manifest.ids.value}'",
        )


def _load_table(directory: Path, spec: PackageTable, source_id: str) -> LoadedTable:
    schema = _schema_for(spec)
    if spec.format_version != schema.version:
        raise IngestRefused(
            "ING002",
            f"table '{spec.name}' is format version {spec.format_version}; "
            f"this atlas reads version {schema.version}",
        )
    on_disk = _read_bytes(directory, spec)
    text = _decompress(on_disk, spec) if spec.file.endswith(".gz") else on_disk
    if b"\r" in text:
        raise IngestRefused("ING011", f"table '{spec.name}' contains a carriage return")
    frame = _parse(text, spec)
    missing = [name for name in schema.column_names if name not in frame.columns]
    extra = [name for name in frame.columns if name not in schema.column_names]
    if missing or extra:
        raise IngestRefused(
            "ING012",
            f"table '{spec.name}' columns differ from the schema; "
            f"missing: {missing or 'none'}; extra: {extra or 'none'}",
        )
    opening = {
        name: frame[name].str.starts_with('"').fill_null(False) for name in schema.column_names
    }
    total = sum(int(flags.sum()) for flags in opening.values())
    if total:
        first = next(name for name, flags in opening.items() if flags.any())
        row = int(opening[first].arg_max() or 0) + 1
        raise IngestRefused(
            "ING015",
            f"table '{spec.name}' has {total} {'cell' if total == 1 else 'cells'} opening "
            f"with a double quote; first column '{first}', first at data row {row} "
            "(1 = the first row after the header)",
        )
    if frame.height != spec.rows:
        raise IngestRefused(
            "ING003", f"table '{spec.name}' has {frame.height} rows; the manifest says {spec.rows}"
        )
    gene_column = schema.gene_column
    if gene_column is not None:
        genes = {cell_text(value) for value in frame[gene_column].to_list()} - {None}
        if len(genes) != spec.genes:
            raise IngestRefused(
                "ING003",
                f"table '{spec.name}' has {len(genes)} distinct genes; "
                f"the manifest says {spec.genes}",
            )
    source_column = schema.source_column
    if source_column is not None:
        foreign = sorted(
            {
                BLANK if cell_text(v) is None else str(v)
                for v in frame[source_column].to_list()
                if v != source_id
            }
        )
        if foreign:
            raise IngestRefused(
                "ING007",
                f"table '{spec.name}' has rows whose source is not '{source_id}': {foreign}",
            )
    return LoadedTable(spec, schema, frame)


def _schema_for(spec: PackageTable) -> TableSchema:
    schema = TABLE_SCHEMAS.get(spec.name)
    if schema is None:
        raise IngestRefused("ING005", f"'{spec.name}' is not a table this atlas knows")
    sharded = table_layout(spec.name)[0] == "sharded"
    if sharded and spec.shard is None:
        raise IngestRefused(
            "ING005", f"table '{spec.name}' is sharded and the manifest gives no shard"
        )
    if not sharded and spec.shard is not None:
        raise IngestRefused("ING005", f"table '{spec.name}' is flat and the manifest gives a shard")
    if spec.genes is None and schema.gene_column is not None:
        raise IngestRefused(
            "ING005",
            f"table '{spec.name}' has a gene column but the manifest gives no gene count",
        )
    if spec.genes is not None and schema.gene_column is None:
        raise IngestRefused(
            "ING005",
            f"table '{spec.name}' has no gene column but the manifest gives a gene count",
        )
    return schema


def _read_bytes(directory: Path, spec: PackageTable) -> bytes:
    try:
        data = (directory / spec.file).read_bytes()
    except OSError as exc:
        raise IngestRefused("ING001", f"table file '{spec.file}' cannot be read: {exc}") from exc
    if hashlib.sha256(data).hexdigest() != spec.sha256:
        raise IngestRefused("ING001", f"table file '{spec.file}' does not match its sha256")
    return data


def _decompress(data: bytes, spec: PackageTable) -> bytes:
    try:
        return gzip.decompress(data)
    except (OSError, EOFError, zlib.error) as exc:
        raise IngestRefused(
            "ING011", f"table file '{spec.file}' is not readable gzip: {exc}"
        ) from exc


def _parse(text: bytes, spec: PackageTable) -> pl.DataFrame:
    """Read every cell as text, refusing a body polars would quietly repair.

    Measured on polars 1.43.1 by reading small byte strings with ``pl.read_csv``:
    an invalid byte in a data cell raises ``ComputeError``, but one in the header
    row is decoded lossily (the column name becomes ``a\ufffd``), which would
    surface as a misleading ING012; and a row with too few fields is padded with
    nulls, which would pass as a row. The explicit decode and the field-count
    check close those two. The field-count check is line-based, so a quoted
    cell spanning lines is refused too; whether any upstream table carries one
    is unmeasured.
    """
    try:
        decoded = text.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise IngestRefused("ING011", f"table '{spec.name}' is not UTF-8: {exc}") from exc
    lines = decoded.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    widths = {line.count("\t") for line in lines}
    if len(widths) > 1:
        raise IngestRefused(
            "ING011", f"table '{spec.name}' has rows with differing numbers of fields"
        )
    try:
        return pl.read_csv(
            io.BytesIO(text), separator="\t", infer_schema_length=0, null_values=[""]
        )
    except pl.exceptions.PolarsError as exc:
        raise IngestRefused("ING011", f"table '{spec.name}' cannot be parsed: {exc}") from exc
