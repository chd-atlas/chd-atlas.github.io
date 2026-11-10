from pathlib import Path

from typer.testing import CliRunner

from chd_atlas.cli import app

runner = CliRunner()


def test_validate_exits_nonzero_on_an_invalid_repository(tmp_path: Path) -> None:
    result = runner.invoke(app, ["validate", "--root", str(tmp_path)])
    assert result.exit_code == 1
    assert "CORPUS001" in result.stdout
    assert "error(s)" in result.stdout


def test_validate_rejects_a_root_that_does_not_exist(tmp_path: Path) -> None:
    """Exit 2, not 1: CI must be able to tell a typo'd --root from a broken repo.

    Without the guard this reports the same content errors as a real but empty
    repository, so a mistyped path looks exactly like a genuine failure.
    """
    missing = tmp_path / "nope"

    result = runner.invoke(app, ["validate", "--root", str(missing)])

    assert result.exit_code == 2
    assert str(missing) in result.stdout
    assert "not a directory" in result.stdout
    assert "CORPUS001" not in result.stdout


def test_validate_rejects_a_root_that_is_not_a_directory(tmp_path: Path) -> None:
    file_root = tmp_path / "README.md"
    file_root.write_text("# not a repository\n")

    result = runner.invoke(app, ["validate", "--root", str(file_root)])

    assert result.exit_code == 2
    assert "not a directory" in result.stdout
    assert "CORPUS001" not in result.stdout


def test_schemas_export_reports_a_target_that_is_not_a_directory(tmp_path: Path) -> None:
    """A bad --target is a bad argument, not a traceback for the user to read."""
    target = tmp_path / "README.md"
    target.write_text("# not a directory\n")

    result = runner.invoke(app, ["schemas", "export", "--target", str(target)])

    assert result.exit_code == 2
    assert not isinstance(result.exception, OSError)
    assert str(target) in result.stdout


def test_schemas_export_reports_an_unwritable_target(tmp_path: Path) -> None:
    target = tmp_path / "readonly" / "schemas"
    target.parent.mkdir()
    target.parent.chmod(0o500)
    try:
        result = runner.invoke(app, ["schemas", "export", "--target", str(target)])
    finally:
        target.parent.chmod(0o700)

    assert result.exit_code == 2
    assert not isinstance(result.exception, OSError)
    assert str(target) in result.stdout


def test_schemas_export_writes_files(tmp_path: Path) -> None:
    result = runner.invoke(app, ["schemas", "export", "--target", str(tmp_path)])
    assert result.exit_code == 0
    assert (tmp_path / "assertion_file.schema.json").is_file()
    assert (tmp_path / "package_manifest.schema.json").is_file()
    assert "wrote 8 schema(s)" in result.stdout
    assert "wrote 16 table schema(s)" in result.stdout


def test_help_lists_both_commands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "validate" in result.stdout
    assert "schemas" in result.stdout
