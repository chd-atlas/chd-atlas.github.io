"""`chd-atlas ingest`: exit 0 prints the result, 1 a refusal, 2 a bad argument or OSError."""

from __future__ import annotations

import hashlib
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
    raise AssertionError("no ingest may run when an argument is unusable")


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


def _result(source_id: str, warnings: tuple[ValidationIssue, ...] = ()) -> IngestResult:
    return IngestResult(
        source_id=source_id,
        version="v1",
        rows_in=2,
        rows_written=2,
        rows_outside_registry=0,
        rows_unresolved=0,
        unresolved=(),
        registry_genes_absent=(),
        registry_genes=185,
        coverage="every_gene",
        dropped=(),
        warnings=warnings,
    )


def test_ingest_all_prints_one_summary_per_source_and_the_validation_once(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """`ingest_all` returns the one validation's warnings on its last result."""
    warning = ValidationIssue("SCP003", Severity.WARNING, "curation/chd_scope.yaml", "msg")
    calls: list[Path] = []

    def fake(root: Path) -> tuple[IngestResult, ...]:
        calls.append(root)
        return (_result("hgnc"), _result("toy_burden", (warning,)))

    monkeypatch.setattr("chd_atlas.cli.ingest_all", fake)

    result = runner.invoke(app, ["ingest", "--all", "--root", str(tmp_path)])

    assert result.exit_code == 0
    assert calls == [tmp_path]
    assert result.stdout.splitlines() == [
        "ingested hgnc v1: 2 row(s) in the package, 2 written",
        "0 row(s) name a gene outside the registry and were not written",
        "0 row(s) name a gene id that did not resolve and were not written",
        "ingested toy_burden v1: 2 row(s) in the package, 2 written",
        "0 row(s) name a gene outside the registry and were not written",
        "0 row(s) name a gene id that did not resolve and were not written",
        "validate: 0 error(s), 1 warning(s)",
        "WARNING [SCP003] curation/chd_scope.yaml: msg",
    ]


def test_ingest_all_with_nothing_pinned_says_so(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_all", lambda _: ())

    result = runner.invoke(app, ["ingest", "--all", "--root", str(tmp_path)])

    assert result.exit_code == 0
    assert result.stdout == (
        "no source in mirrors/sources.yaml is pinned to a package; nothing was ingested\n"
    )


@pytest.mark.parametrize(
    "args",
    [["--all", "toy_burden"], ["--all", "--package", "."], []],
    ids=["all-and-a-source", "all-and-a-package", "neither"],
)
def test_all_takes_no_source_or_package_and_one_of_them_is_required(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, args: list[str]
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_all", _never_called)
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)

    result = runner.invoke(app, ["ingest", *args, "--root", str(tmp_path)])

    assert result.exit_code == 2
    assert result.stdout == "error: give a SOURCE_ID, or --all with no SOURCE_ID or --package\n"


TOY_PIN = (
    "sources:\n"
    "  - id: toy_burden\n"
    "    name: Toy burden study\n"
    '    version: "v1"\n'
    "    retrieved_on: 2026-11-01\n"
    "    url: https://example.org/toy\n"
    "    licence: CC-BY-4.0\n"
    "    redistribution: permitted\n"
)


def test_a_pinned_source_with_no_cache_exits_1_with_ing017(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "sources.yaml").write_text(
        TOY_PIN + f'    package: {{version: "v1", manifest_sha256: "{"0" * 64}"}}\n',
        encoding="utf-8",
    )

    result = runner.invoke(app, ["ingest", "toy_burden", "--root", str(tmp_path)])

    assert result.exit_code == 1
    assert result.stdout.startswith("ING017: ")
    assert (
        "gh release download toy_burden/v1 -R chd-atlas/feeds -D local/packages/toy_burden/v1"
        in result.stdout
    )


def test_a_source_with_no_pin_and_no_package_exits_2_naming_the_flag(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "sources.yaml").write_text(TOY_PIN, encoding="utf-8")

    result = runner.invoke(app, ["ingest", "toy_burden", "--root", str(tmp_path)])

    assert result.exit_code == 2
    assert result.stdout == (
        "error: source 'toy_burden' has no package pin in mirrors/sources.yaml; "
        "give the package directory with --package\n"
    )


def test_a_pinned_and_cached_source_ingests_from_its_cache(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Without --package, the pinned version's cache is the package directory."""
    cache = tmp_path / "local" / "packages" / "toy_burden" / "v1"
    cache.mkdir(parents=True)
    (cache / "manifest.json").write_bytes(b"{}")
    digest = hashlib.sha256(b"{}").hexdigest()
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "sources.yaml").write_text(
        TOY_PIN + f'    package: {{version: "v1", manifest_sha256: "{digest}"}}\n',
        encoding="utf-8",
    )
    calls: list[tuple[Path, str, Path]] = []

    def fake(root: Path, source_id: str, package: Path) -> IngestResult:
        calls.append((root, source_id, package))
        return _result(source_id)

    monkeypatch.setattr("chd_atlas.cli.ingest_package", fake)

    result = runner.invoke(app, ["ingest", "toy_burden", "--root", str(tmp_path)])

    assert result.exit_code == 0
    assert calls == [(tmp_path, "toy_burden", cache)]
