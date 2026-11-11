"""Deterministic schema export, written to the committed ``schemas/`` directory.

JSON Schema for the YAML record models and for a source package's
``manifest.json`` (`EXPORTED_MODELS`), and a Frictionless-style Table Schema
per mirror table (`export_table_schemas`), which is also the format of a
package's tables.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Final

import polars as pl
from pydantic import BaseModel

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.models.assertion import AssertionFile
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalFile
from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import FeaturedFile, PhenotypeFile, PublicationFile
from chd_atlas.models.package import PackageManifest
from chd_atlas.tables import TABLE_SCHEMAS, Column, TableSchema, table_layout

EXPORTED_MODELS: Final[dict[str, type[BaseModel]]] = {
    "assertion_file": AssertionFile,
    "functional_file": FunctionalFile,
    "gene_literature_file": GeneLiteratureFile,
    "publication_file": PublicationFile,
    "featured_file": FeaturedFile,
    "phenotype_file": PhenotypeFile,
    "dataset": Dataset,
    "package_manifest": PackageManifest,
}


def export_schemas(target: Path) -> list[Path]:
    """Write one ``*.schema.json`` per model. Output is stable across runs."""
    target.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, model in sorted(EXPORTED_MODELS.items()):
        schema = model.model_json_schema()
        schema["title"] = model.__name__
        path = target / f"{name}.schema.json"
        payload = json.dumps(schema, indent=2, sort_keys=True) + "\n"
        write_bytes_atomically(path, payload.encode("utf-8"))
        written.append(path)
    return written


TABLE_SCHEMA_DIR: Final = "tables"

_TABLE_COMMENT: Final = (
    "Exported by `chd-atlas schemas export` from src/chd_atlas/tables.py; do not edit."
)

_FIELD_TYPES: Final[dict[Any, str]] = {
    pl.String: "string",
    pl.Int64: "integer",
    pl.Float64: "number",
    pl.Boolean: "boolean",
}


# Plan Decision 5, stated where a feeds author reads the table's contract.
_GENE_COLUMN_NOTE: Final = (
    "In a source package this column holds the ids the manifest's `ids` names. "
    "When that is not `hgnc_id`, the cells are the source's own ids (symbols or "
    "Ensembl gene ids) and this field's `pattern` does not apply to them: "
    "`chd-atlas ingest` rewrites each to the HGNC id it resolves to, and names "
    "every id that does not resolve and writes none of its rows."
)

# The `genes` table is the registry itself, so the generic note would mislead:
# it ships only in the `hgnc` package and its gene column is already `hgnc_id`.
_GENES_TABLE_NOTE: Final = (
    "The `genes` table ships only in the `hgnc` package, whose manifest `ids` is "
    "`hgnc_id`, so this column holds HGNC ids as they stand. `chd-atlas ingest` "
    "writes only the rows of genes in the registry the atlas derives from its scope "
    "and the validity mirrors, not the whole package."
)


def _field(column: Column, *, gene: bool = False, table: str = "") -> dict[str, Any]:
    field_type = _FIELD_TYPES.get(column.dtype)
    if field_type is None:
        raise ValueError(f"column {column.name!r} has unexportable dtype {column.dtype!r}")
    constraints: dict[str, Any] = {"required": not column.nullable}
    if column.allowed is not None:
        constraints["enum"] = sorted(column.allowed)
    if column.pattern is not None:
        constraints["pattern"] = column.pattern
    if column.minimum is not None:
        constraints["minimum"] = column.minimum
    if column.maximum is not None:
        constraints["maximum"] = column.maximum
    field: dict[str, Any] = {"name": column.name, "type": field_type, "constraints": constraints}
    if gene:
        field["description"] = _GENES_TABLE_NOTE if table == "genes" else _GENE_COLUMN_NOTE
    return field


def _table_document(schema: TableSchema) -> dict[str, Any]:
    kind, location = table_layout(schema.name)
    return {
        "$comment": _TABLE_COMMENT,
        "name": schema.name,
        "version": schema.version,
        "layout": {kind: location},
        "fields": [
            _field(column, gene=column.name == schema.gene_column, table=schema.name)
            for column in schema.columns
        ],
        "missingValues": [""],
        # TBL007 refuses a duplicate sort key, so the sort key is a primary key.
        "primaryKey": list(schema.sort_key),
        "sortKey": list(schema.sort_key),
        "geneColumn": schema.gene_column,
        "sourceColumn": schema.source_column,
    }


def export_table_schemas(target: Path) -> list[Path]:
    """Write ``<name>.v<version>.json`` per mirror table; remove superseded versions."""
    target.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, schema in sorted(TABLE_SCHEMAS.items()):
        path = target / f"{name}.v{schema.version}.json"
        payload = json.dumps(_table_document(schema), indent=2, sort_keys=True) + "\n"
        write_bytes_atomically(path, payload.encode("utf-8"))
        written.append(path)
    for stale in sorted(target.glob("*.v*.json")):
        if stale not in written:
            stale.unlink()
    return written
