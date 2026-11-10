"""The ``manifest.json`` an outside source tool ships with each data package.

The atlas does not build packages; it verifies and ingests them. This model is
the contract, exported as ``schemas/package_manifest.schema.json`` so a feeds
tool can validate its own output before delivering it.
"""

from __future__ import annotations

import json
from datetime import date
from pathlib import Path, PurePosixPath
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator, model_validator

from chd_atlas.vocab import GeneIdKind, Redistribution

_SHA256 = r"^[0-9a-f]{64}$"
_COMMIT = r"^[0-9a-f]{40}$"


class ManifestError(Exception):
    """A manifest that cannot be read or does not satisfy `PackageManifest`."""


class PackageTable(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1)
    format_version: int = Field(ge=1)
    file: str = Field(min_length=1)
    sha256: str = Field(pattern=_SHA256)
    rows: int = Field(ge=0)
    genes: int | None = Field(default=None, ge=0)
    shard: str | None = Field(default=None, pattern=r"^[A-Za-z0-9._-]+$")

    @field_validator("file")
    @classmethod
    def file_stays_inside_the_package(cls, value: str) -> str:
        """Resolved as ``package_dir / file``; an absolute path or ``..`` would escape it."""
        path = PurePosixPath(value)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(f"file '{value}' must be a relative path inside the package")
        return value


class PackageInput(BaseModel):
    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1)
    url: str = Field(min_length=1)
    sha256: str = Field(pattern=_SHA256)
    retrieved_on: date
    upstream_release: str | None = None


class InputCounts(BaseModel):
    model_config = ConfigDict(extra="forbid")

    rows: int = Field(ge=0)
    genes: int | None = Field(default=None, ge=0)


class DroppedRows(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reason: str = Field(min_length=1)
    count: int = Field(ge=1)
    examples: list[str] = Field(default_factory=list, max_length=5)


class RegistrySlice(BaseModel):
    model_config = ConfigDict(extra="forbid")

    atlas_commit: str = Field(pattern=_COMMIT)
    genes: int = Field(ge=0)


class SliceCoverage(BaseModel):
    """The spec's nested form: ``{"registry_slice": {"atlas_commit": ..., "genes": ...}}``."""

    model_config = ConfigDict(extra="forbid")

    registry_slice: RegistrySlice


class BuiltBy(BaseModel):
    model_config = ConfigDict(extra="forbid")

    recipe_commit: str = Field(pattern=_COMMIT)
    hvantk_release: str | None = None
    image_digest: str | None = Field(default=None, pattern=r"^sha256:[0-9a-f]{64}$")
    host: str = Field(min_length=1)
    job: str | None = None


class PackageSource(BaseModel):
    """`Source`'s fields minus ``id``, ``ontology_*`` and ``package``.

    A test pins every name and type here against `Source`.
    """

    model_config = ConfigDict(extra="forbid")

    name: str = Field(min_length=1)
    version: str = Field(min_length=1)
    retrieved_on: date
    url: str = Field(min_length=1)
    licence: str = Field(min_length=1)
    redistribution: Redistribution
    file_version: str | None = None
    file_date: str | None = None


class PackageManifest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    package: str = Field(pattern=r"^[a-z][a-z0-9_]*$")
    version: str = Field(pattern=r"^[A-Za-z0-9._-]+$")
    ids: GeneIdKind
    coverage: Literal["every_gene"] | SliceCoverage
    tables: list[PackageTable] = Field(min_length=1)
    inputs: list[PackageInput] = Field(min_length=1)
    input_counts: InputCounts
    dropped: list[DroppedRows]
    built_by: BuiltBy
    source: PackageSource
    licence_evidence: str = Field(min_length=1)

    @field_validator("licence_evidence")
    @classmethod
    def licence_evidence_is_not_blank(cls, value: str) -> str:
        if not value.strip():
            raise ValueError("licence_evidence must not be blank")
        return value

    @model_validator(mode="after")
    def table_names_are_unique(self) -> PackageManifest:
        seen: set[str] = set()
        repeated: set[str] = set()
        for table in self.tables:
            if table.name in seen:
                repeated.add(table.name)
            seen.add(table.name)
        if repeated:
            raise ValueError(f"table names repeated: {', '.join(sorted(repeated))}")
        return self


def load_manifest(path: Path) -> PackageManifest:
    """Read and validate a manifest; every failure is a `ManifestError` naming the file."""
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise ManifestError(f"{path}: cannot read manifest: {exc}") from exc
    try:
        return PackageManifest.model_validate(raw)
    except ValidationError as exc:
        problems = "; ".join(
            f"{'.'.join(str(part) for part in error['loc']) or '<root>'}: {error['msg']}"
            for error in exc.errors()
        )
        raise ManifestError(f"{path}: {problems}") from exc
