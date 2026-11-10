import json
from pathlib import Path

import pytest

from chd_atlas.schema_export import (
    EXPORTED_MODELS,
    TABLE_SCHEMA_DIR,
    export_schemas,
    export_table_schemas,
)
from chd_atlas.tables import BURDEN, TABLE_SCHEMAS


def test_writes_one_schema_file_per_model(tmp_path: Path) -> None:
    written = export_schemas(tmp_path)
    assert {path.name for path in written} == {f"{name}.schema.json" for name in EXPORTED_MODELS}
    for path in written:
        assert path.is_file()


def test_emitted_schema_is_valid_json_with_a_title(tmp_path: Path) -> None:
    export_schemas(tmp_path)
    payload = json.loads((tmp_path / "assertion_file.schema.json").read_text(encoding="utf-8"))
    assert payload["title"] == "AssertionFile"


def test_written_schemas_are_world_readable(tmp_path: Path) -> None:
    """These are committed artifacts; mkstemp's 0600 would leave them owner-only."""
    for path in export_schemas(tmp_path):
        assert path.stat().st_mode & 0o777 == 0o644, path


def test_export_is_byte_identical_across_runs(tmp_path: Path) -> None:
    export_schemas(tmp_path)
    first = (tmp_path / "assertion_file.schema.json").read_bytes()
    export_schemas(tmp_path)
    assert (tmp_path / "assertion_file.schema.json").read_bytes() == first


def test_a_failed_write_leaves_the_previous_schema_intact(tmp_path: Path) -> None:
    """A truncate-then-fail would leave a schema the drift test reports as stale
    rather than corrupt.

    The containing directory rather than the file itself is made read-only:
    `os.replace` renames over a destination regardless of that destination's own
    mode, so a read-only file would be overwritten and nothing would raise.
    Removing write permission from the directory is what stops the write.
    `assertion_file.schema.json` sorts first, so it is the file the failing
    export would otherwise have written over.
    """
    export_schemas(tmp_path)
    target = tmp_path / "assertion_file.schema.json"
    before = target.read_bytes()

    tmp_path.chmod(0o555)
    try:
        with pytest.raises(OSError):
            export_schemas(tmp_path)
        assert target.read_bytes() == before
    finally:
        tmp_path.chmod(0o755)


def test_committed_schemas_match_the_models() -> None:
    """Fails when a model changes without re-running `chd-atlas schemas export`."""
    committed = Path("schemas")
    if not committed.is_dir():
        return
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        fresh = Path(tmp)
        export_schemas(fresh)
        for path in sorted(fresh.glob("*.schema.json")):
            target = committed / path.name
            assert target.is_file(), f"{path.name} is not committed"
            assert target.read_text(encoding="utf-8") == path.read_text(encoding="utf-8"), (
                f"{path.name} is stale; run `chd-atlas schemas export`"
            )
        fresh_tables = fresh / TABLE_SCHEMA_DIR
        export_table_schemas(fresh_tables)
        committed_tables = committed / TABLE_SCHEMA_DIR
        assert {p.name for p in committed_tables.glob("*.json")} == {
            p.name for p in fresh_tables.glob("*.json")
        }
        for path in sorted(fresh_tables.glob("*.json")):
            assert (committed_tables / path.name).read_bytes() == path.read_bytes(), (
                f"{path.name} is stale; run `chd-atlas schemas export`"
            )


def test_every_table_schema_is_exported_once_under_its_version(tmp_path: Path) -> None:
    written = export_table_schemas(tmp_path)
    assert len(written) == 16
    assert {path.name for path in written} == {f"{name}.v1.json" for name in TABLE_SCHEMAS}


def test_an_exported_table_states_what_validate_table_checks(tmp_path: Path) -> None:
    export_table_schemas(tmp_path)

    def load(name: str) -> dict:  # type: ignore[type-arg]
        return json.loads((tmp_path / f"{name}.v1.json").read_text(encoding="utf-8"))

    burden = load("burden")
    fields = {f["name"]: f for f in burden["fields"]}
    for column in BURDEN.columns:
        assert fields[column.name]["constraints"]["required"] is (not column.nullable)
    assert fields["count_unit"]["constraints"]["enum"] == sorted(
        ["alleles", "de_novo_mutations", "individuals"]
    )
    assert "pattern" in fields["study"]["constraints"]
    assert "enum" not in fields["study"]["constraints"]
    assert burden["primaryKey"] == list(BURDEN.sort_key)
    assert burden["sortKey"] == list(BURDEN.sort_key)
    assert burden["geneColumn"] == "gene"
    assert burden["sourceColumn"] == "source"
    assert burden["layout"] == {"flat": "burden.tsv"}
    assert burden["missingValues"] == [""]

    assert load("genes")["geneColumn"] == "hgnc_id"
    assert load("profile_quantiles")["geneColumn"] is None
    assert load("profile_quantiles")["sourceColumn"] is None
    assert load("profiles")["layout"] == {"sharded": "profiles"}


def test_a_superseded_table_version_is_removed(tmp_path: Path) -> None:
    stray = tmp_path / "burden.v0.json"
    stray.write_text("{}\n", encoding="utf-8")
    export_table_schemas(tmp_path)
    assert not stray.exists()
