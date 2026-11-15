"""`chd-atlas ingest`: exit 0 prints the result, 1 a refusal, 2 a bad argument or OSError."""

from __future__ import annotations

import hashlib
from dataclasses import replace
from pathlib import Path
from typing import Never

import pytest
from typer.testing import CliRunner

from chd_atlas.cli import app
from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest.apply import IngestAllResult, IngestResult, PrunedRows
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

    The committed Farah shard lacks 5 of the 186 registry genes: POPDC1, CRIPTO,
    CIROZ and CIROP (measured 2026-11-10 against 185), and PLAGL1, which entered
    the registry on 2026-11-14 after the slice was cut. CRIPTO at least is a
    symbol gap (the source names it TDGF1), not a registry that grew after the
    slice was cut; PLAGL1 is the opposite case.
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
        rows_outside_filter=0,
        filter_reasons=(),
        rows_unknown_term=0,
        unknown_terms=(),
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


def test_a_filtered_ingest_names_the_rows_outside_the_filter_and_the_unknown_terms(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Two separate lines: a term the release does not know is not "outside the filter".

    The filter line prints only for a package one of `ingest.filters` applied
    to (`test_ingest_all_prints_...` holds that it does not print otherwise),
    the unknown-term line only when there is a term to name.
    """
    reason = "hpo not under HP:0001626 in ontologies/hp-2026-06-23.obo"
    result_ = IngestResult(
        source_id="orphanet_product4",
        version="4f44e8a6",
        rows_in=10,
        rows_written=5,
        rows_outside_registry=0,
        rows_unresolved=0,
        unresolved=(),
        rows_outside_filter=3,
        filter_reasons=(reason,),
        rows_unknown_term=2,
        unknown_terms=("HP:7777777",),
        registry_genes_absent=(),
        registry_genes=185,
        coverage="every_gene",
        dropped=(),
        warnings=(),
    )
    monkeypatch.setattr("chd_atlas.cli.ingest_package", lambda *_: result_)

    result = runner.invoke(
        app, ["ingest", "orphanet_product4", "--root", str(tmp_path), "--package", str(tmp_path)]
    )

    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        "ingested orphanet_product4 4f44e8a6: 10 row(s) in the package, 5 written",
        "0 row(s) name a gene outside the registry and were not written",
        "0 row(s) name a gene id that did not resolve and were not written",
        f"3 row(s) are outside the filter ({reason})",
        "2 row(s) name a term the pinned release does not know: HP:7777777",
        "validate: 0 error(s), 0 warning(s)",
    ]


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
        rows_outside_filter=0,
        filter_reasons=(),
        rows_unknown_term=0,
        unknown_terms=(),
        registry_genes_absent=(),
        registry_genes=185,
        coverage="every_gene",
        dropped=(),
        warnings=warnings,
    )


def test_ingest_all_prints_each_summary_every_pruned_row_and_the_validation_once(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """One summary per source, the pruned rows by file with their total, one validation."""
    warning = ValidationIssue("SCP003", Severity.WARNING, "curation/chd_scope.yaml", "msg")
    calls: list[tuple[Path, dict[str, Path]]] = []
    pruned = (
        PrunedRows("burden", "mirrors/burden.tsv", "HGNC:1", 6),
        PrunedRows("burden", "mirrors/burden.tsv", "HGNC:2", 2),
        PrunedRows("profiles", "mirrors/profiles/E-MTAB-6814.tsv", "HGNC:1", 119),
    )

    def fake(root: Path, repin: dict[str, Path]) -> IngestAllResult:
        calls.append((root, repin))
        return IngestAllResult((_result("hgnc"), _result("toy_burden")), pruned, (warning,))

    monkeypatch.setattr("chd_atlas.cli.ingest_all", fake)

    result = runner.invoke(app, ["ingest", "--all", "--root", str(tmp_path)])

    assert result.exit_code == 0
    assert calls == [(tmp_path, {})]
    assert result.stdout.splitlines() == [
        "ingested hgnc v1: 2 row(s) in the package, 2 written",
        "0 row(s) name a gene outside the registry and were not written",
        "0 row(s) name a gene id that did not resolve and were not written",
        "ingested toy_burden v1: 2 row(s) in the package, 2 written",
        "0 row(s) name a gene outside the registry and were not written",
        "0 row(s) name a gene id that did not resolve and were not written",
        "removed 127 row(s) naming 2 gene(s) that left the registry: HGNC:1, HGNC:2",
        "removed from mirrors/burden.tsv: HGNC:1 (6), HGNC:2 (2)",
        "removed from mirrors/profiles/E-MTAB-6814.tsv: HGNC:1 (119)",
        "validate: 0 error(s), 1 warning(s)",
        "WARNING [SCP003] curation/chd_scope.yaml: msg",
    ]


def test_ingest_all_with_nothing_pinned_says_so_and_prints_the_validation(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    warning = ValidationIssue("SCP003", Severity.WARNING, "curation/chd_scope.yaml", "msg")
    monkeypatch.setattr("chd_atlas.cli.ingest_all", lambda *_: IngestAllResult((), (), (warning,)))

    result = runner.invoke(app, ["ingest", "--all", "--root", str(tmp_path)])

    assert result.exit_code == 0
    assert result.stdout.splitlines() == [
        "no source in mirrors/sources.yaml is pinned to a package; nothing was ingested",
        "validate: 0 error(s), 1 warning(s)",
        "WARNING [SCP003] curation/chd_scope.yaml: msg",
    ]


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


DOWNLOAD = "gh release download toy_burden/v1 -R chd-atlas/feeds -D local/packages/toy_burden/v1"


@pytest.mark.parametrize(
    ("directory", "message"),
    [
        (
            False,
            "toy_burden v1 is pinned in mirrors/sources.yaml but not cached at "
            f"local/packages/toy_burden/v1; download it: {DOWNLOAD}",
        ),
        (
            True,
            "local/packages/toy_burden/v1/manifest.json cannot be read (No such file or "
            f"directory); remove the directory and download it again: {DOWNLOAD}",
        ),
    ],
    ids=["no-directory", "directory-without-manifest"],
)
def test_a_pinned_source_with_no_cache_exits_1_with_ing017(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, directory: bool, message: str
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "sources.yaml").write_text(
        TOY_PIN + f'    package: {{version: "v1", manifest_sha256: "{"0" * 64}"}}\n',
        encoding="utf-8",
    )
    if directory:
        (tmp_path / "local" / "packages" / "toy_burden" / "v1").mkdir(parents=True)

    result = runner.invoke(app, ["ingest", "toy_burden", "--root", str(tmp_path)])

    assert result.exit_code == 1
    assert result.stdout == f"ING017: {message}\n"


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


@pytest.mark.parametrize(
    ("args", "message"),
    [
        (
            ["toy_burden", "--repin", "clingen={dir}"],
            "--repin needs --all (a single ingest takes its package with --package)",
        ),
        (["--all", "--repin", "clingen"], "--repin clingen is not SOURCE=DIR"),
        (["--all", "--repin", "=x"], "--repin =x is not SOURCE=DIR"),
        (
            ["--all", "--repin", "clingen={dir}", "--repin", "clingen={dir}"],
            "--repin names clingen more than once",
        ),
        (
            ["--all", "--repin", "clingen={dir}/nope"],
            "--repin clingen={dir}/nope: {dir}/nope is not a directory",
        ),
    ],
    ids=["without-all", "no-equals", "no-source", "twice", "not-a-directory"],
)
def test_an_unusable_repin_exits_2(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, args: list[str], message: str
) -> None:
    monkeypatch.setattr("chd_atlas.cli.ingest_all", _never_called)
    monkeypatch.setattr("chd_atlas.cli.ingest_package", _never_called)

    result = runner.invoke(
        app, ["ingest", *(arg.format(dir=tmp_path) for arg in args), "--root", str(tmp_path)]
    )

    assert result.exit_code == 2
    assert result.stdout == f"error: {message.format(dir=tmp_path)}\n"


def test_ingest_all_repins_a_source_from_its_directory_and_prints_the_new_version(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Each `--repin SOURCE=DIR` reaches `ingest_all`; the summary names the version it pinned."""
    clingen, gencc = tmp_path / "clingen", tmp_path / "gencc"
    clingen.mkdir()
    gencc.mkdir()
    calls: list[tuple[Path, dict[str, Path]]] = []

    def fake(root: Path, repin: dict[str, Path]) -> IngestAllResult:
        calls.append((root, repin))
        repinned = replace(_result("clingen"), version="2026-11-12")
        return IngestAllResult((repinned, _result("hgnc")), (), ())

    monkeypatch.setattr("chd_atlas.cli.ingest_all", fake)

    result = runner.invoke(
        app,
        [
            "ingest",
            "--all",
            "--repin",
            f"clingen={clingen}",
            "--repin",
            f"gencc={gencc}",
            "--root",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0
    assert calls == [(tmp_path, {"clingen": clingen, "gencc": gencc})]
    assert result.stdout.splitlines()[0] == (
        "ingested clingen 2026-11-12: 2 row(s) in the package, 2 written"
    )


def test_the_summary_shows_25_unknown_ids_and_writes_every_one(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """30 unknown ids and 1 ambiguous: the ambiguous in full, 25 unknowns, all 31 in a file."""
    unknown = [Resolved(f"ENSG{n:03d}", Resolution.UNKNOWN) for n in range(30)]
    ambiguous = Resolved("ENSGAMB", Resolution.AMBIGUOUS, None, ("HGNC:6772", "HGNC:6773"))
    # Out of order on purpose: the summary and the file sort them.
    unresolved = (ambiguous, *reversed(unknown))
    result_ = replace(
        _result("e_mtab_6814"),
        version="v9",
        rows_in=33,
        rows_written=2,
        rows_unresolved=31,
        unresolved=unresolved,
    )
    monkeypatch.setattr("chd_atlas.cli.ingest_package", lambda *_: result_)

    result = runner.invoke(
        app, ["ingest", "e_mtab_6814", "--root", str(tmp_path), "--package", str(tmp_path)]
    )

    path = "local/ingest/e_mtab_6814-v9.unresolved.tsv"
    shown = "; ".join(f"ENSG{n:03d} (unknown)" for n in range(25))
    assert result.exit_code == 0
    assert (
        "31 gene id(s) did not resolve and their rows were not written: "
        f"ENSGAMB (ambiguous: HGNC:6772, HGNC:6773); {shown}; "
        f"… and 5 more (all 31 in {path})"
    ) in result.stdout.splitlines()
    assert "ENSG025" not in result.stdout
    assert (tmp_path / path).is_file(), f"{path} was not written"
    lines = (tmp_path / path).read_text(encoding="utf-8").splitlines()
    assert lines[0] == "id\tresolution\tcandidates"
    assert lines[1:] == [
        "ENSG000\tunknown\t",
        *[f"ENSG{n:03d}\tunknown\t" for n in range(1, 30)],
        "ENSGAMB\tambiguous\tHGNC:6772|HGNC:6773",
    ]


@pytest.mark.parametrize("every", [False, True], ids=["single", "all"])
def test_an_unwritable_unresolved_list_is_named_and_the_summary_still_prints(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, every: bool
) -> None:
    """The ingest has committed by the time the list is written, so a failed write is a line.

    `local/ingest` is a file here, so the mkdir raises OSError. The exit code stays 0, every
    summary line of every source is printed, and one line names the path and the error.
    """
    (tmp_path / "local").mkdir()
    (tmp_path / "local" / "ingest").write_text("in the way", encoding="utf-8")
    first = replace(
        _result("alpha"),
        rows_unresolved=1,
        unresolved=(Resolved("NOSUCH", Resolution.UNKNOWN),),
    )
    second = replace(_result("beta"), rows_unresolved=1, unresolved=first.unresolved)
    if every:
        monkeypatch.setattr(
            "chd_atlas.cli.ingest_all", lambda *_: IngestAllResult((first, second), (), ())
        )
        args = ["ingest", "--all", "--root", str(tmp_path)]
    else:
        monkeypatch.setattr("chd_atlas.cli.ingest_package", lambda *_: first)
        args = ["ingest", "alpha", "--root", str(tmp_path), "--package", str(tmp_path)]

    result = runner.invoke(app, args)

    assert result.exit_code == 0
    lines = result.stdout.splitlines()
    for source in ("alpha", "beta") if every else ("alpha",):
        assert f"ingested {source} v1: 2 row(s) in the package, 2 written" in lines
        assert (
            "1 gene id(s) did not resolve and their rows were not written: NOSUCH (unknown)"
        ) in lines
        path = f"local/ingest/{source}-v1.unresolved.tsv"
        failed = [line for line in lines if line.startswith("could not write") and path in line]
        assert len(failed) == 1
        assert "Errno" in failed[0] or "exist" in failed[0].lower()


_ALL_NOTE = (
    "note: a single ingest writes only this mirror; `chd-atlas ingest --all` (with "
    "`--repin <id>=<dir>` for a new version) re-filters the registry, the "
    "registry-filtered mirrors and mondo_xrefs"
)


@pytest.mark.parametrize("count", [25, 26], ids=["exactly-25", "26"])
def test_the_summary_caps_unknown_ids_at_25_and_names_the_file_only_past_it(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, count: int
) -> None:
    """The boundary: 25 unknowns are all shown with no suffix and no path; 26 add `… and 1 more`.

    Kills `more >= 0` (a suffix at exactly 25) and printing the path unconditionally.
    The file is written either way.
    """
    unknown = [Resolved(f"ENSG{n:03d}", Resolution.UNKNOWN) for n in range(count)]
    result_ = replace(
        _result("e_mtab_6814"),
        version="v9",
        rows_unresolved=count,
        unresolved=tuple(reversed(unknown)),
    )
    monkeypatch.setattr("chd_atlas.cli.ingest_package", lambda *_: result_)

    result = runner.invoke(
        app, ["ingest", "e_mtab_6814", "--root", str(tmp_path), "--package", str(tmp_path)]
    )

    path = "local/ingest/e_mtab_6814-v9.unresolved.tsv"
    shown = "; ".join(f"ENSG{n:03d} (unknown)" for n in range(25))
    head = f"{count} gene id(s) did not resolve and their rows were not written: {shown}"
    expected = head if count == 25 else f"{head}; … and 1 more (all 26 in {path})"
    assert result.exit_code == 0
    assert expected in result.stdout.splitlines()
    assert (path in result.stdout) is (count == 26)
    assert "… and" not in result.stdout.replace(expected, "")
    assert len((tmp_path / path).read_text(encoding="utf-8").splitlines()) == count + 1


@pytest.mark.parametrize(
    ("source", "noted"),
    [("clingen", True), ("gencc", True), ("hgnc", False), ("alpha", False)],
)
def test_a_single_validity_ingest_says_only_ingest_all_refilters_the_registry(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, source: str, noted: bool
) -> None:
    """The note is pinned against a literal; ingest --all and other sources print none."""
    monkeypatch.setattr("chd_atlas.cli.ingest_package", lambda *_: _result(source))

    single = runner.invoke(
        app, ["ingest", source, "--root", str(tmp_path), "--package", str(tmp_path)]
    )
    assert single.exit_code == 0
    assert (_ALL_NOTE in single.stdout.splitlines()) is noted

    if noted:
        monkeypatch.setattr(
            "chd_atlas.cli.ingest_all",
            lambda *_: IngestAllResult((_result(source),), (), ()),
        )
        everything = runner.invoke(app, ["ingest", "--all", "--root", str(tmp_path)])
        assert everything.exit_code == 0
        assert _ALL_NOTE not in everything.stdout.splitlines()
