"""Deterministic JSON Schema emission for the YAML record models."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Final

from pydantic import BaseModel

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.models.assertion import AssertionFile
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalFile
from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import FeaturedFile, PhenotypeFile, PublicationFile

EXPORTED_MODELS: Final[dict[str, type[BaseModel]]] = {
    "assertion_file": AssertionFile,
    "functional_file": FunctionalFile,
    "gene_literature_file": GeneLiteratureFile,
    "publication_file": PublicationFile,
    "featured_file": FeaturedFile,
    "phenotype_file": PhenotypeFile,
    "dataset": Dataset,
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
