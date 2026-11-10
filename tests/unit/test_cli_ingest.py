"""`chd-atlas ingest`: exit 0 prints the result, 1 a refusal, 2 a bad argument or OSError."""

from __future__ import annotations

from pathlib import Path
from typing import Never

import pytest
from typer.testing import CliRunner

from chd_atlas.cli import app
from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest.apply import IngestResult
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.package import DroppedRows, RegistrySlice
from tests.unit.package_fixtures import write_package

runner = CliRunner()


def _never_called(*_: object) -> Never:
    raise AssertionError("ingest_package must not run when an argument is not a directory")


@pytest.mark.parametrize("bad", ["--root", "--package"])
def test_an_argument_that_is_not_a_directory_exits_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)
    missing = tmp_path / "nope"
    paths = {"--root": tmp_path, "--package": tmp_path, bad: missing}

    result = runner.invoke(
        app,
        [
            "ingest",
            "toy_burden",
            "--root",
            str(paths["--root"]),
            "--package",
            str(paths["--package"]),
        ],
    )

    assert result.exit_code == 2
    assert f"{bad} {missing} is not a directory" in result.stdout


def test_a_refused_package_exits_1_and_prints_its_code(tmp_path: Path) -> None:
    package = write_package(tmp_path)  # package `toy_burden`, ingested as another source

    result = runner.invoke(
        app, ["ingest", "other_source", "--root", str(tmp_path), "--package", str(package)]
    )

    assert result.exit_code == 1
    assert result.stdout.startswith(
        "ING006: the package is 'toy_burden', not source 'other_source'"
    )


SLICE = RegistrySlice(atlas_commit="a" * 40, genes=183)


@pytest.mark.parametrize(
    ("coverage", "absent"),
    [
        ("every_gene", "2 registry gene(s) have no resolved row in the package: HGNC:1, HGNC:2"),
        (
            SLICE,
            f"2 registry gene(s) have no resolved row in this registry slice (cut at {'a' * 40} "
            "over 183 genes; the registry now lists 185): HGNC:1, HGNC:2",
        ),
    ],
    ids=["every_gene", "registry_slice"],
)
def test_an_ingested_package_exits_0_and_names_everything_it_did_not_write(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, coverage: object, absent: str
) -> None:
    """The absent-genes line states the slice and the registry, never a cause.

    The committed Farah shard lacks 4 of the 185 registry genes (POPDC1, CRIPTO,
    CIROZ, CIROP; measured 2026-11-10), and CRIPTO at least is a symbol gap (the
    source names it TDGF1), not a registry that grew after the slice was cut.
    The ingest itself is `test_ingest_apply`'s (on the real corpus); this pins
    what the command prints of its result.
    """
    calls: list[tuple[Path, str, Path]] = []
    result_ = IngestResult(
        source_id="farah2024_hoc",
        version="v2",
        rows_in=6,
        rows_written=3,
        rows_outside_registry=1,
        rows_unresolved=2,
        unresolved=(
            Resolved("MADH7", Resolution.AMBIGUOUS, None, ("HGNC:6772", "HGNC:6773")),
            Resolved("NOSUCH", Resolution.UNKNOWN),
        ),
        registry_genes_absent=("HGNC:1", "HGNC:2"),
        registry_genes=185,
        coverage=coverage,  # type: ignore[arg-type]
        dropped=(
            DroppedRows(reason="too few cells", count=4, examples=[]),
            DroppedRows(reason="symbol not in HGNC", count=2, examples=["FOO1", "BAR2"]),
        ),
        warnings=(ValidationIssue("SCP003", Severity.WARNING, "curation/chd_scope.yaml", "msg"),),
    )

    def fake(root: Path, source_id: str, package: Path) -> IngestResult:
        calls.append((root, source_id, package))
        return result_

    monkeypatch.setattr("chd_atlas.cli.ingest_package", fake)

    result = runner.invoke(
        app, ["ingest", "farah2024_hoc", "--root", str(tmp_path), "--package", str(tmp_path)]
    )

    assert result.exit_code == 0
    assert calls == [(tmp_path, "farah2024_hoc", tmp_path)]
    assert result.stdout.splitlines() == [
        "ingested farah2024_hoc v2: 6 row(s) in the package, 3 written",
        "1 row(s) name a gene outside the registry and were not written",
        "2 row(s) name a gene id that did not resolve and were not written",
        "2 gene id(s) did not resolve and their rows were not written: "
        "MADH7 (ambiguous: HGNC:6772, HGNC:6773); NOSUCH (unknown)",
        absent,
        "dropped before packaging: too few cells: 4",
        "dropped before packaging: symbol not in HGNC: 2 (e.g. FOO1, BAR2)",
        "validate: 0 error(s), 1 warning(s)",
        "WARNING [SCP003] curation/chd_scope.yaml: msg",
    ]
    assert "grew" not in result.stdout


def test_a_filesystem_error_exits_2_with_its_message(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Exit 1 means refused with the checkout restored; an OSError is exit 2, as in `build`.

    Whether the write failed and was restored or the restore failed too, the
    command prints the error as it came (`apply._restore` names the paths).
    """
    message = "the ingest failed and could not restore: /r/mirrors/burden.tsv: injected"

    def fail(*_: object) -> Never:
        raise OSError(message)

    monkeypatch.setattr("chd_atlas.cli.ingest_package", fail)

    result = runner.invoke(
        app, ["ingest", "toy_burden", "--root", str(tmp_path), "--package", str(tmp_path)]
    )

    assert result.exit_code == 2
    assert result.stdout == f"error: {message}\n"
