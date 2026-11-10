from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.package import ManifestError, PackageManifest, PackageSource, load_manifest
from chd_atlas.schema_export import EXPORTED_MODELS
from chd_atlas.validate.sources import PackagePin, Source
from tests.unit.package_fixtures import manifest_dict


@pytest.mark.parametrize(
    "coverage",
    ["every_gene", {"registry_slice": {"atlas_commit": "b" * 40, "genes": 181}}],
)
def test_a_complete_manifest_round_trips(coverage: Any) -> None:
    raw = manifest_dict(coverage=coverage)
    manifest = PackageManifest.model_validate(raw)
    assert manifest.model_dump(mode="json") == raw
    assert manifest_dict() is not manifest_dict()


def test_a_bare_registry_slice_without_its_wrapper_is_refused() -> None:
    raw = manifest_dict(coverage={"atlas_commit": "b" * 40, "genes": 181})
    with pytest.raises(ValidationError, match="coverage"):
        PackageManifest.model_validate(raw)


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


@pytest.mark.parametrize("file", ["a\x00.tsv", "a\n.tsv", "a\t.tsv", "a\x7f.tsv", "a\x85.tsv"])
def test_a_table_file_with_a_control_character_is_refused(file: str) -> None:
    """A NUL reached `Path.read_bytes` as a ValueError, past `load_package`'s refusals."""
    raw = manifest_dict()
    raw["tables"][0]["file"] = file
    with pytest.raises(ValidationError, match=r"tables\.0\.file"):
        PackageManifest.model_validate(raw)


@pytest.mark.parametrize("name", [".", ".."])
@pytest.mark.parametrize("where", ["version", "shard", "pin"])
def test_a_dot_or_dot_dot_version_or_shard_is_refused(name: str, where: str) -> None:
    """Each becomes a path segment: `..` moved the package cache to `local/packages/`."""
    raw = manifest_dict()
    with pytest.raises(ValidationError):
        if where == "version":
            PackageManifest.model_validate({**raw, "version": name})
        elif where == "shard":
            raw["tables"][0]["shard"] = name
            PackageManifest.model_validate(raw)
        else:
            PackagePin.model_validate({"version": name, "manifest_sha256": "a" * 64})


@pytest.mark.parametrize("name", ["2026-11-10.1", "...", ".a", "..a", "a..", "v_1-2"])
def test_the_rest_of_the_version_and_shard_alphabet_is_accepted(name: str) -> None:
    raw = manifest_dict(version=name)
    raw["tables"][0]["shard"] = name
    assert PackageManifest.model_validate(raw).version == name
    assert PackagePin.model_validate({"version": name, "manifest_sha256": "a" * 64}).version == name


def test_two_tables_with_one_name_are_refused() -> None:
    raw = manifest_dict()
    raw["tables"].append(dict(raw["tables"][0], file="burden2.tsv"))
    with pytest.raises(ValidationError, match="table names repeated: burden"):
        PackageManifest.model_validate(raw)


def test_package_source_fields_are_sources_fields_with_the_same_types() -> None:
    """Both ways: a field added to `Source` must reach the manifest's `source` block."""
    copied = set(Source.model_fields) - {"id", "ontology_prefix", "ontology_file", "package"}
    assert set(PackageSource.model_fields) == copied
    for name, field in PackageSource.model_fields.items():
        assert field.annotation == Source.model_fields[name].annotation, name


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
