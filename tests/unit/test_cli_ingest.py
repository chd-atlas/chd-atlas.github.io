"""`chd-atlas ingest`: exit 0 prints the result, 1 a refusal, 2 a bad argument."""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Never

import pytest
from typer.testing import CliRunner

from chd_atlas.cli import app
from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest.apply import IngestResult
from chd_atlas.models.package import DroppedRows, RegistrySlice
from tests.unit.package_fixtures import write_package

REPO_ROOT = Path(__file__).parent.parent.parent
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


def test_an_ingested_package_exits_0_and_names_what_the_source_lacks(tmp_path: Path) -> None:
    """End to end on a copy of the real repository, with an `every_gene` package."""
    repo = tmp_path / "repo"
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, repo / item)
    header, *lines = (REPO_ROOT / "mirrors" / "burden.tsv").read_text().splitlines()
    columns = header.split("\t")
    rows = [
        {**dict(zip(columns, line.split("\t"), strict=True)), "source": "toy_burden"}
        | {"maf_max": "0.0001"}
        for line in lines
        if line.startswith("PMID:40127276\tHGNC:10023\t")
    ]
    (tmp_path / "pkg").mkdir()
    package = write_package(tmp_path / "pkg", rows, ids="hgnc_id")
    registry = len((repo / "mirrors" / "genes.tsv").read_text().splitlines()) - 1

    result = runner.invoke(
        app, ["ingest", "toy_burden", "--root", str(repo), "--package", str(package)]
    )

    assert result.exit_code == 0, result.stdout
    out = result.stdout.splitlines()
    assert out[0] == "ingested toy_burden 2026-11-10.1: 6 row(s) in the package, 6 written"
    assert "0 row(s) name a gene outside the registry and were not written" in out
    absent = f"the source reports no row for {registry - 1} registry gene(s): HGNC:"
    assert any(line.startswith(absent) for line in out)
    assert "dropped before packaging: symbol not in HGNC: 2 (e.g. FOO1)" in out


def test_a_registry_slice_says_the_registry_grew(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls: list[tuple[Path, str, Path]] = []
    result_ = IngestResult(
        source_id="farah2024_hoc",
        version="v2",
        rows_in=5,
        rows_written=3,
        rows_outside_registry=1,
        unresolved=(
            Resolved("MADH7", Resolution.AMBIGUOUS, None, ("HGNC:6772", "HGNC:6773")),
            Resolved("NOSUCH", Resolution.UNKNOWN),
        ),
        registry_genes_absent=("HGNC:1", "HGNC:2"),
        coverage=RegistrySlice(atlas_commit="a" * 40, genes=183),
        dropped=(DroppedRows(reason="too few cells", count=4, examples=[]),),
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
        "ingested farah2024_hoc v2: 5 row(s) in the package, 3 written",
        "1 row(s) name a gene outside the registry and were not written",
        "2 gene id(s) did not resolve and their rows were not written: "
        "MADH7 (ambiguous: HGNC:6772, HGNC:6773); NOSUCH (unknown)",
        f"2 registry gene(s) absent because the registry grew after the slice was cut at "
        f"{'a' * 40}: HGNC:1, HGNC:2",
        "dropped before packaging: too few cells: 4",
    ]
