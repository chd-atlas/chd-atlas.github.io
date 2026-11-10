from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.package import ManifestError, PackageManifest, PackageSource, load_manifest
from chd_atlas.schema_export import EXPORTED_MODELS
from chd_atlas.validate.sources import Source
from tests.unit.package_fixtures import manifest_dict


def test_a_complete_manifest_round_trips() -> None:
    raw = manifest_dict()
    manifest = PackageManifest.model_validate(raw)
    assert manifest.model_dump(mode="json") == raw
    assert manifest_dict() is not manifest_dict()


def test_a_manifest_without_a_recipe_commit_is_refused() -> None:
    raw = manifest_dict()
    del raw["built_by"]["recipe_commit"]
    with pytest.raises(ValidationError, match="recipe_commit"):
        PackageManifest.model_validate(raw)


@pytest.mark.parametrize("file", ["../x.tsv", "/etc/x.tsv"])
def test_a_table_file_outside_the_package_is_refused(file: str) -> None:
    raw = manifest_dict()
    raw["tables"][0]["file"] = file
    with pytest.raises(ValidationError, match="relative path inside the package"):
        PackageManifest.model_validate(raw)


def test_two_tables_with_one_name_are_refused() -> None:
    raw = manifest_dict()
    raw["tables"].append(dict(raw["tables"][0], file="burden2.tsv"))
    with pytest.raises(ValidationError, match="burden"):
        PackageManifest.model_validate(raw)


def test_package_source_fields_are_sources_fields_with_the_same_types() -> None:
    for name, field in PackageSource.model_fields.items():
        assert name in Source.model_fields, name
        assert field.annotation == Source.model_fields[name].annotation, name
    assert "package" not in PackageSource.model_fields


def test_the_manifest_schema_is_exported() -> None:
    assert EXPORTED_MODELS["package_manifest"] is PackageManifest


def test_load_manifest_names_the_failing_field(tmp_path: Path) -> None:
    raw: dict[str, Any] = manifest_dict()
    raw["tables"][0]["sha256"] = "nothex"
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(raw))
    with pytest.raises(ManifestError, match=r"tables\.0\.sha256"):
        load_manifest(path)
    path.write_text("{not json")
    with pytest.raises(ManifestError, match="manifest.json"):
        load_manifest(path)
    with pytest.raises(ManifestError, match="absent.json"):
        load_manifest(tmp_path / "absent.json")
