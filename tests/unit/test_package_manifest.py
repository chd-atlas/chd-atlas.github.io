from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from pydantic import ValidationError

from chd_atlas.models.package import ManifestError, PackageManifest, PackageSource, load_manifest
from chd_atlas.validate.sources import PackagePin, Source
from tests.unit.package_fixtures import manifest_dict


def test_coverage_is_every_gene_or_the_registry_slice_wrapper() -> None:
    """Spec §4.3's shape: a slice is `{"registry_slice": {...}}`, never a bare slice."""
    for coverage in ["every_gene", {"registry_slice": {"atlas_commit": "b" * 40, "genes": 181}}]:
        raw = manifest_dict(coverage=coverage)
        assert PackageManifest.model_validate(raw).model_dump(mode="json") == raw
    assert manifest_dict() is not manifest_dict()

    bare = manifest_dict(coverage={"atlas_commit": "b" * 40, "genes": 181})
    with pytest.raises(ValidationError, match="coverage"):
        PackageManifest.model_validate(bare)


@pytest.mark.parametrize(
    "file", ["../x.tsv", "/etc/x.tsv", "a\x00.tsv", "a\n.tsv", "a\t.tsv", "a\x7f.tsv", "a\x85.tsv"]
)
def test_a_table_file_that_is_not_a_plain_path_inside_the_package_is_refused(file: str) -> None:
    """A NUL reached `Path.read_bytes` as a ValueError, past `load_package`'s refusals."""
    raw = manifest_dict()
    raw["tables"][0]["file"] = file
    with pytest.raises(ValidationError, match=r"tables\.0\.file"):
        PackageManifest.model_validate(raw)


@pytest.mark.parametrize(
    ("name", "accepted"),
    [(".", False), ("..", False)]
    + [(name, True) for name in ["2026-11-10.1", "...", ".a", "..a", "a..", "v_1-2"]],
)
def test_a_version_or_shard_is_a_path_segment_never_dot_or_dot_dot(
    name: str, accepted: bool
) -> None:
    """Each becomes a path segment: `..` moved the package cache to `local/packages/`."""
    sharded = manifest_dict()
    sharded["tables"][0]["shard"] = name
    for check in (
        lambda: PackageManifest.model_validate(manifest_dict(version=name)),
        lambda: PackageManifest.model_validate(sharded),
        lambda: PackagePin.model_validate({"version": name, "manifest_sha256": "a" * 64}),
    ):
        if accepted:
            check()
        else:
            with pytest.raises(ValidationError):
                check()


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
