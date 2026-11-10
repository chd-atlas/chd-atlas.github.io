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

import polars as pl

from chd_atlas.models.package import (
    ManifestError,
    PackageManifest,
    PackageTable,
    load_manifest,
)
from chd_atlas.tables import TABLE_SCHEMAS, TableSchema, cell_text, table_layout


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
    tables = tuple(_load_table(directory, spec, source_id) for spec in manifest.tables)
    return LoadedPackage(directory, manifest, manifest_bytes, tables)


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
        foreign = sorted({str(v) for v in frame[source_column].to_list() if v != source_id})
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
    if (spec.genes is None) != (schema.gene_column is None):
        raise IngestRefused(
            "ING005",
            f"table '{spec.name}' must give a gene count exactly when it has a gene column",
        )
    if spec.name == "genes":
        raise IngestRefused(
            "ING009", "the genes table ships in the HGNC package, which is not ingested yet"
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

    Measured with polars 1.43: ``pl.read_csv`` decodes invalid UTF-8 lossily
    rather than raising, and pads a row with too few fields with nulls. Both
    would turn a damaged file into plausible rows, so they are checked here. The
    field-count check is line-based, so a quoted cell spanning lines is refused
    too; whether any upstream table carries one is unmeasured.
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
