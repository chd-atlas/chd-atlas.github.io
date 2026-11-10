"""Deterministic JSON Schema emission for the YAML record models."""

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


def _field(column: Column) -> dict[str, Any]:
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
    return {"name": column.name, "type": field_type, "constraints": constraints}


def _table_document(schema: TableSchema) -> dict[str, Any]:
    kind, location = table_layout(schema.name)
    return {
        "$comment": _TABLE_COMMENT,
        "name": schema.name,
        "version": schema.version,
        "layout": {kind: location},
        "fields": [_field(column) for column in schema.columns],
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
