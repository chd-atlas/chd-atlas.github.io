"""`ingest_package`: a verified package merged into a copy of the real repository.

Every test that writes works on a temporary copy of `curation/`, `mirrors/` and
`ontologies/` (as `tests/test_repository_validates.py` does), never on the
checkout. `validate_repository` over the real corpus costs seconds, so the
happy path is ingested once per module, twice (the second time to prove a
re-ingest changes nothing), and its tests read the snapshots the fixture took.
The second ingest runs with validation stubbed: its bytes are asserted to be
the first's, which validation already passed, so a real run would re-answer it.

Package rows are `sierant2025_ds6`'s own rows for HGNC:10023 with `source`
rewritten to `toy_burden` and `maf_max` set to `0.0001`, a value no committed
burden row uses (measured: sierant2025_ds6 uses only null and 1e-05,
audain2026_sd3 0.001, audain2021_s6 0.01). Copying real rows keeps every
study, cohort and publication reference valid, and the new `maf_max` keeps the
sort key unique (TBL007). It is also the value that makes the sort observable:
as text "0.0001" sorts *before* "1e-05", as a number it sorts after, and
`validate_sort_order` compares numbers.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest import apply
from chd_atlas.ingest.apply import IngestResult, ingest_all, ingest_package, with_source_entry
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.ingest.resolve import package_cache
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.package import PackageSource, RegistrySlice
from chd_atlas.registry import RegistryUnavailable, derive_registry_from_root
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, table_layout
from chd_atlas.validate.runner import ValidationReport, validate_repository
from chd_atlas.validate.sort_order import validate_sort_order
from chd_atlas.validate.sources import PackagePin, load_sources
from tests.unit.package_fixtures import manifest_dict, tsv_bytes, write_hgnc_root, write_package

REPO_ROOT = Path(__file__).parent.parent.parent
SOURCE = "toy_burden"
VERSION = "2026-11-10.1"
SLICE_COMMIT = "0123456789abcdef0123456789abcdef01234567"
SENTINEL = "  # SENTINEL: a comment above clingen, which an ingest must not move\n"
TOY_MAF = "0.0001"
OUTSIDE_GENE = "HGNC:5"  # A1BG: a valid HGNC id, not in mirrors/genes.tsv


def copy_repo(target: Path) -> Path:
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, target / item)
    return target


def burden_lines(repo: Path) -> list[str]:
    return (repo / "mirrors" / "burden.tsv").read_text(encoding="utf-8").splitlines()


def sierant_rows(gene: str) -> list[dict[str, str]]:
    """`sierant2025_ds6`'s committed rows for one gene, cells verbatim."""
    header, *lines = burden_lines(REPO_ROOT)
    columns = header.split("\t")
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in lines]
    return [row for row in rows if row["gene"] == gene and row["source"] == "sierant2025_ds6"]


def toy_rows() -> list[dict[str, str]]:
    """Six registry rows, one outside the registry, one id that does not resolve."""
    rows = [{**row, "source": SOURCE, "maf_max": TOY_MAF} for row in sierant_rows("HGNC:10023")]
    assert len(rows) == 6, "the fixture is sized on HGNC:10023's six sierant2025_ds6 rows"
    rows.append({**rows[1], "gene": OUTSIDE_GENE})
    # A symbol in an `hgnc_id` package: the identity resolver names it UNKNOWN.
    rows.append({**rows[1], "gene": "TBX5"})
    return rows


def tree(repo: Path) -> dict[str, bytes | None]:
    """Every file (bytes) and directory (None) under mirrors/ and local/."""
    found: dict[str, bytes | None] = {}
    for top in ("mirrors", "local"):
        base = repo / top
        if not base.exists():
            continue
        found[top] = None
        for path in sorted(base.rglob("*")):
            found[path.relative_to(repo).as_posix()] = path.read_bytes() if path.is_file() else None
    return found


def registry(repo: Path) -> set[str]:
    header, *lines = (repo / "mirrors" / "genes.tsv").read_text(encoding="utf-8").splitlines()
    column = header.split("\t").index("hgnc_id")
    return {line.split("\t")[column] for line in lines}


class Ingested:
    def __init__(self, base: Path) -> None:
        self.repo = copy_repo(base / "repo")
        sources = self.repo / "mirrors" / "sources.yaml"
        text = sources.read_text(encoding="utf-8")
        assert text.count("  - id: clingen\n") == 1
        sources.write_text(text.replace("  - id: clingen\n", SENTINEL + "  - id: clingen\n"))
        (base / "pkg").mkdir()
        self.package = write_package(
            base / "pkg",
            toy_rows(),
            ids="hgnc_id",
            coverage={"registry_slice": {"atlas_commit": SLICE_COMMIT, "genes": 185}},
        )
        self.before = tree(self.repo)
        self.burden_before = burden_lines(self.repo)
        self.first = ingest_package(self.repo, SOURCE, self.package)
        self.after = tree(self.repo)
        self.burden_after = burden_lines(self.repo)
        with pytest.MonkeyPatch.context() as patch:
            report = ValidationReport(list(self.first.warnings))
            patch.setattr(apply, "validate_repository", lambda _: report)
            self.second = ingest_package(self.repo, SOURCE, self.package)
        self.again = tree(self.repo)


@pytest.fixture(scope="module")
def ingested(tmp_path_factory: pytest.TempPathFactory) -> Ingested:
    return Ingested(tmp_path_factory.mktemp("ingested"))


def test_ingest_replaces_only_the_packages_own_burden_rows(ingested: Ingested) -> None:
    others = [line for line in ingested.burden_after if not line.endswith(f"\t{SOURCE}")]
    assert others == ingested.burden_before, "another source's row moved or changed"

    columns = ingested.burden_after[0].split("\t")

    def cell(line: str, name: str) -> str:
        return line.split("\t")[columns.index(name)]

    block = [
        (cell(line, "study"), cell(line, "consequence_class"), cell(line, "origin"))
        + (cell(line, "maf_max"), cell(line, "source"))
        for line in ingested.burden_after
        if cell(line, "gene") == "HGNC:10023" and cell(line, "study") == "PMID:40127276"
    ]
    # Null first, then 1e-05 before 0.0001 (numbers, not text), ties broken by
    # the key and never by source: the order `validate_sort_order` checks.
    sierant, toy = "sierant2025_ds6", SOURCE
    assert block == [
        ("PMID:40127276", "damaging", "de_novo", "", sierant),
        ("PMID:40127276", "damaging", "de_novo", TOY_MAF, toy),
        ("PMID:40127276", "damaging", "transmitted_or_unphased", "1e-05", sierant),
        ("PMID:40127276", "damaging", "transmitted_or_unphased", TOY_MAF, toy),
        ("PMID:40127276", "lof", "de_novo", "", sierant),
        ("PMID:40127276", "lof", "de_novo", TOY_MAF, toy),
        ("PMID:40127276", "lof", "transmitted_or_unphased", "1e-05", sierant),
        ("PMID:40127276", "lof", "transmitted_or_unphased", TOY_MAF, toy),
        ("PMID:40127276", "missense_damaging", "de_novo", "", sierant),
        ("PMID:40127276", "missense_damaging", "de_novo", TOY_MAF, toy),
        ("PMID:40127276", "missense_damaging", "transmitted_or_unphased", "1e-05", sierant),
        ("PMID:40127276", "missense_damaging", "transmitted_or_unphased", TOY_MAF, toy),
    ]
    path = ingested.repo / "mirrors" / "burden.tsv"
    assert validate_sort_order(path, TABLE_SCHEMAS["burden"]) == []


def test_the_result_counts_and_names_what_was_not_written(ingested: Ingested) -> None:
    result = ingested.first
    assert (
        result.rows_in,
        result.rows_written,
        result.rows_outside_registry,
        result.rows_unresolved,
    ) == (8, 6, 1, 1)
    assert result.unresolved == (Resolved("TBX5", Resolution.UNKNOWN),)
    genes = {line.split("\t")[1] for line in ingested.burden_after[1:]}
    assert OUTSIDE_GENE not in genes
    assert "TBX5" not in genes
    # The warnings validation reported with the package applied (the corpus's
    # own SCP003/SCP007), which `chd-atlas ingest` prints after its summary.
    assert result.warnings
    assert all(issue.severity is Severity.WARNING for issue in result.warnings)

    assert result.coverage == RegistrySlice(atlas_commit=SLICE_COMMIT, genes=185)
    assert result.registry_genes == len(registry(ingested.repo))
    expected = registry(ingested.repo) - {"HGNC:10023"}
    assert len(expected) == len(registry(ingested.repo)) - 1
    # HGNC:11604 (TBX5) is in the registry and was named only by symbol, so it
    # is absent; HGNC:5 is outside the registry, so it is not "absent" either.
    assert "HGNC:11604" in result.registry_genes_absent
    assert result.registry_genes_absent == tuple(sorted(expected))
    assert result.source_id == SOURCE
    assert result.version == VERSION
    assert [entry.reason for entry in result.dropped] == ["symbol not in HGNC"]


def test_sources_yaml_comments_survive_and_a_new_entry_is_appended(ingested: Ingested) -> None:
    before = ingested.before["mirrors/sources.yaml"]
    after = ingested.after["mirrors/sources.yaml"]
    assert before is not None and after is not None
    assert SENTINEL.encode() in before
    assert after.startswith(before), "a byte before the appended entry changed"
    manifest = ingested.after[f"mirrors/packages/{SOURCE}.json"]
    assert manifest == (ingested.package / "manifest.json").read_bytes()
    digest = hashlib.sha256(manifest).hexdigest()
    assert after[len(before) :].decode() == (
        f"  - id: {SOURCE}\n"
        "    name: Toy burden study\n"
        f"    version: {VERSION}\n"
        "    retrieved_on: 2026-11-01\n"
        "    url: https://example.org/toy\n"
        "    licence: CC-BY-4.0\n"
        "    redistribution: permitted\n"
        "    package:\n"
        f"      version: {VERSION}\n"
        f"      manifest_sha256: {digest}\n"
    )
    registry_, issues = load_sources(ingested.repo)
    assert issues == []
    assert registry_.sources[-1].id == SOURCE
    assert registry_.sources[-1].package == PackagePin(version=VERSION, manifest_sha256=digest)


def test_reingesting_the_same_version_changes_no_byte(ingested: Ingested) -> None:
    assert ingested.again == ingested.after
    assert ingested.second == ingested.first
    cache = f"local/packages/{SOURCE}/{VERSION}"
    for name in ("manifest.json", "burden.tsv"):
        assert ingested.after[f"{cache}/{name}"] == (ingested.package / name).read_bytes()
    # Only what the ingest wrote changed: the burden mirror, the committed
    # manifest, sources.yaml and the cache.
    assert changed_paths(ingested.before, ingested.after) == {
        "mirrors/burden.tsv",
        "mirrors/sources.yaml",
        "mirrors/packages",
        f"mirrors/packages/{SOURCE}.json",
        "local",
        "local/packages",
        f"local/packages/{SOURCE}",
        cache,
        f"{cache}/manifest.json",
        f"{cache}/burden.tsv",
    }


def recording_writes(
    monkeypatch: pytest.MonkeyPatch,
    *,
    refuse: str | None = None,
    refuse_after_first: str | None = None,
) -> list[Path]:
    """Record every path `apply` writes, and fail writes by file name.

    Every write to a file named ``refuse`` raises OSError and leaves it as it
    was, as a read-only file or a full disk would. Writes to a file named
    ``refuse_after_first`` succeed once and raise from then on, so the
    restore's write of it fails.
    """
    written: list[Path] = []

    def write(path: Path, payload: bytes) -> None:
        earlier = written.count(path)
        written.append(path)
        if path.name == refuse or (path.name == refuse_after_first and earlier):
            raise OSError(f"injected: cannot write {path}")
        write_bytes_atomically(path, payload)

    monkeypatch.setattr(apply, "write_bytes_atomically", write)
    return written


def test_a_package_failing_validation_restores_every_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = copy_repo(tmp_path / "repo")
    rows = toy_rows()[:6]
    # BUR019: a de novo mutation count on a row whose variant set is not de novo.
    rows[1] = {**rows[1], "count_unit": "de_novo_mutations"}
    assert rows[1]["origin"] == "transmitted_or_unphased"
    (tmp_path / "pkg").mkdir()
    package = write_package(tmp_path / "pkg", rows, ids="hgnc_id")
    before = tree(repo)
    written = recording_writes(monkeypatch)

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, SOURCE, package)

    assert excinfo.value.code == "ING010"
    assert "BUR019" in excinfo.value.message
    # The refusal came after the writes, so the restore had work to undo.
    mirrors = repo / "mirrors"
    cache = repo / "local" / "packages" / SOURCE / VERSION
    for path in (
        mirrors / "burden.tsv",
        mirrors / "packages" / f"{SOURCE}.json",
        mirrors / "sources.yaml",
        cache / "manifest.json",
    ):
        assert path in written
    assert tree(repo) == before
    assert not (repo / "local").exists()
    assert not (mirrors / "packages").exists()


def test_a_write_that_raises_restores_every_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = copy_repo(tmp_path / "repo")
    (tmp_path / "pkg").mkdir()
    package = write_package(tmp_path / "pkg", toy_rows()[:6], ids="hgnc_id")
    before = tree(repo)
    written = recording_writes(monkeypatch, refuse="sources.yaml")
    validated: list[Path] = []
    monkeypatch.setattr(apply, "validate_repository", validated.append)

    # The write's own error, not "could not restore": sources.yaml was never
    # changed, so the restore does not try to rewrite it (and fail again).
    with pytest.raises(OSError, match=r"^injected: cannot write \S+/mirrors/sources\.yaml$"):
        ingest_package(repo, SOURCE, package)

    assert repo / "mirrors" / "burden.tsv" in written
    assert repo / "mirrors" / "packages" / f"{SOURCE}.json" in written
    assert validated == []
    assert tree(repo) == before


def test_a_restore_that_cannot_complete_names_only_what_it_left_changed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """sources.yaml refuses every write; burden.tsv takes the ingest's write, not the restore's.

    The message names burden.tsv alone, the one path left changed, and the
    write's own error survives as its context.
    """
    repo = copy_repo(tmp_path / "repo")
    (tmp_path / "pkg").mkdir()
    package = write_package(tmp_path / "pkg", toy_rows()[:6], ids="hgnc_id")
    before = tree(repo)
    recording_writes(monkeypatch, refuse="sources.yaml", refuse_after_first="burden.tsv")
    monkeypatch.setattr(apply, "validate_repository", _never_validated)

    with pytest.raises(OSError) as excinfo:
        ingest_package(repo, SOURCE, package)

    burden = repo / "mirrors" / "burden.tsv"
    assert str(excinfo.value) == (
        f"the ingest failed and could not restore: {burden}: injected: cannot write {burden}"
    )
    original = excinfo.value.__context__
    assert str(original) == f"injected: cannot write {repo / 'mirrors' / 'sources.yaml'}"
    assert changed_paths(before, tree(repo)) == {"mirrors/burden.tsv"}, (
        "the message must name exactly what was left changed"
    )


def changed_paths(before: dict[str, bytes | None], after: dict[str, bytes | None]) -> set[str]:
    """Paths added, removed, or (files) with different bytes."""
    return {
        key
        for key in before.keys() | after.keys()
        if (key in before) != (key in after) or before.get(key) != after.get(key)
    }


def _never_validated(_: Path) -> ValidationReport:
    raise AssertionError("validation must not run after a write fails")


def test_a_validity_table_is_replaced_whole_and_not_filtered(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`clingen_validity` has no source column, so a package owns the whole file.

    And it is one of `UNFILTERED_TABLES`: the registry is derived from it, so a
    row for a gene outside the registry (HGNC:10004, RGS9) is written, not
    counted out. Validation is stubbed to pass: a two-row ClinGen mirror fails
    it for reasons (scope, the publication gate) that are not this test's.
    """
    repo = tmp_path / "repo"
    (repo / "mirrors").mkdir(parents=True)
    for name in ("genes.tsv", "sources.yaml", "clingen_gene_validity.tsv"):
        shutil.copy(REPO_ROOT / "mirrors" / name, repo / "mirrors" / name)
    header, *lines = (REPO_ROOT / "mirrors" / "clingen_gene_validity.tsv").read_text().splitlines()
    chosen = [
        next(line for line in lines if line.startswith("HGNC:10023\t")),
        next(line for line in lines if line.startswith("HGNC:10004\t")),
    ]
    assert "HGNC:10004" not in registry(repo) and "HGNC:10023" in registry(repo)
    columns = header.split("\t")
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in chosen]
    data = tsv_bytes(rows, columns)
    package = tmp_path / "pkg"
    package.mkdir()
    (package / "clingen.tsv").write_bytes(data)
    table = {
        "name": "clingen_validity",
        "format_version": 1,
        "file": "clingen.tsv",
        "sha256": hashlib.sha256(data).hexdigest(),
        "rows": 2,
        "genes": 2,
        "shard": None,
    }
    manifest = manifest_dict(package="clingen", ids="hgnc_id", tables=[table])
    (package / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport())

    result = ingest_package(repo, "clingen", package)

    assert (result.rows_written, result.rows_outside_registry) == (2, 0)
    mirror = (repo / "mirrors" / "clingen_gene_validity.tsv").read_text().splitlines()
    assert mirror == [header, chosen[1], chosen[0]]


def test_a_shard_keeps_other_sources_rows_and_every_other_shard(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A `cell_type_profiles` package owns its own rows of shard ``hoc`` and nothing else.

    `hoc.tsv` starts with one `farah2024_hoc` row (foreign, kept byte for byte)
    and one stale `toy_cells` row (owned, replaced). A second shard holding a
    `toy_cells` row is another shard, so it is kept. Rows are four real
    committed hoc rows, chosen so the foreign one sorts between the package's
    two. Validation is stubbed, as in the validity-table test.
    """
    repo = tmp_path / "repo"
    shards = repo / "mirrors" / "cell_type_profiles"
    shards.mkdir(parents=True)
    for name in ("genes.tsv", "sources.yaml"):
        shutil.copy(REPO_ROOT / "mirrors" / name, repo / "mirrors" / name)
    header, *lines = (REPO_ROOT / "mirrors/cell_type_profiles/hoc.tsv").read_text().splitlines()
    ivs, la, lv, fourth = (line.rsplit("\t", 1)[0] for line in lines[:4])
    assert lines[1].endswith("\tfarah2024_hoc")
    foreign = f"{la}\tfarah2024_hoc"
    (shards / "hoc.tsv").write_text(f"{header}\n{foreign}\n{fourth}\ttoy_cells\n")
    (shards / "other.tsv").write_text(f"{header}\n{fourth}\ttoy_cells\n")
    columns = header.split("\t")
    package_lines = [f"{ivs}\ttoy_cells", f"{lv}\ttoy_cells"]
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in package_lines]
    data = tsv_bytes(rows, columns)
    package = tmp_path / "pkg"
    package.mkdir()
    (package / "cells.tsv").write_bytes(data)
    table = {
        "name": "cell_type_profiles",
        "format_version": 1,
        "file": "cells.tsv",
        "sha256": hashlib.sha256(data).hexdigest(),
        "rows": 2,
        "genes": 1,
        "shard": "hoc",
    }
    manifest = manifest_dict(package="toy_cells", ids="hgnc_id", tables=[table])
    (package / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport())
    before = tree(repo)

    ingest_package(repo, "toy_cells", package)

    after = tree(repo)
    assert (shards / "hoc.tsv").read_text().splitlines() == [
        header,
        package_lines[0],
        foreign,
        package_lines[1],
    ]
    cache = "local/packages/toy_cells/2026-11-10.1"
    assert changed_paths(before, after) == {
        "mirrors/cell_type_profiles/hoc.tsv",
        "mirrors/sources.yaml",
        "mirrors/packages",
        "mirrors/packages/toy_cells.json",
        "local",
        "local/packages",
        "local/packages/toy_cells",
        cache,
        f"{cache}/manifest.json",
        f"{cache}/cells.tsv",
    }


@pytest.mark.parametrize(
    ("extra_header", "extra_cell", "message"),
    [("\tnote", "\tx", "has columns"), ("", "\tx", "row 2 has 30 fields; the header has 29")],
    ids=["a column the schema lacks", "a row longer than its header"],
)
def test_a_mirror_the_ingest_would_have_to_trim_is_refused(
    tmp_path: Path, extra_header: str, extra_cell: str, message: str
) -> None:
    """Rewriting either in the schema's columns would drop a cell nobody chose to drop."""
    header, first, *_ = burden_lines(REPO_ROOT)
    path = tmp_path / "burden.tsv"
    path.write_text(f"{header}{extra_header}\n{first}{extra_cell}\n", encoding="utf-8")

    with pytest.raises(IngestRefused) as excinfo:
        apply._read_mirror(path, TABLE_SCHEMAS["burden"])

    assert excinfo.value.code == "ING010"
    assert message in excinfo.value.message


@pytest.mark.parametrize("genes_tsv", [None, "header only"])
def test_a_missing_or_empty_registry_is_refused_before_any_write(
    tmp_path: Path, genes_tsv: str | None
) -> None:
    """With no registry the filter would drop every row and call it a fact about the package."""
    mirrors = tmp_path / "repo" / "mirrors"
    mirrors.mkdir(parents=True)
    if genes_tsv is not None:
        (mirrors / "genes.tsv").write_text("\t".join(TABLE_SCHEMAS["genes"].column_names) + "\n")
    (tmp_path / "pkg").mkdir()
    package = write_package(tmp_path / "pkg", toy_rows()[:6], ids="hgnc_id")
    before = tree(tmp_path / "repo")

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(tmp_path / "repo", SOURCE, package)

    assert excinfo.value.code == "ING010"
    assert "mirrors/genes.tsv is absent or lists no gene" in excinfo.value.message
    assert tree(tmp_path / "repo") == before


@pytest.mark.parametrize(
    ("path", "table"),
    [(path.relative_to(REPO_ROOT), table) for path, table in mirror_paths(REPO_ROOT)],
    ids=lambda value: str(value),
)
def test_every_committed_mirror_reads_and_rewrites_byte_for_byte(path: Path, table: str) -> None:
    """Rows another source owns survive an ingest as the bytes they were.

    The ingest reads a mirror as text, re-sorts it with the package's rows and
    writes it back; with no package rows that must reproduce the committed file
    exactly. Measured 2026-11-10 over all 11 committed mirrors (burden.tsv's
    1,805 rows, gencc_submissions.tsv's 30,410, ...): byte-identical, every
    header already in the schema's column order, no file containing a quote
    character. A typed read would not be: `read_table` then polars `write_csv`
    rewrites 171 lines of burden.tsv (measured the same day), `1e-05` becoming
    `0.00001`.
    """
    schema = TABLE_SCHEMAS[table]
    data = (REPO_ROOT / path).read_bytes()
    assert apply._render(schema, apply._read_mirror(REPO_ROOT / path, schema)) == data


# --- the sources.yaml entry, as a pure function (no ingest, no validation) ---

ORPHANET_ENTRY_TAIL = "    redistribution: permitted_with_attribution\n  - id: orphanet_product1\n"
OLD_FILE_LINES = (
    '    version: "sha256:4f44e8a61201399911aa1ba44a293c0ccaa5ce11272c47d40862255cb72f6b32"\n'
    '    file_version: "1.3.42 / 4.1.8 [2025-03-03]"\n'
    '    file_date: "2026-06-23 07:57:18"\n'
)


def _replace_once(text: str, old: str, new: str) -> str:
    assert text.count(old) == 1, old
    return text.replace(old, new)


def _source_as_package(source_id: str, **changes: Any) -> PackageSource:
    registry_, issues = load_sources(REPO_ROOT)
    assert issues == []
    fields = registry_.by_id[source_id].model_dump(include=set(PackageSource.model_fields))
    return PackageSource.model_validate({**fields, **changes})


def test_an_entry_updated_in_place_changes_no_line_outside_it() -> None:
    """orphanet_product4, in the middle of the file, loses its file header fields.

    Two comments a careless edit would move: one after ``file_date`` (which
    ruamel attaches to the deleted key) and one above the next entry (which it
    attaches to this entry's last key, so a key appended after it would land
    below the comment). The end-of-line comment on ``file_date`` goes with its
    line. Everything else is byte-identical, and a second identical update is a
    no-op.
    """
    text = (REPO_ROOT / "mirrors" / "sources.yaml").read_text(encoding="utf-8")
    text = _replace_once(
        text,
        OLD_FILE_LINES,
        OLD_FILE_LINES.replace('07:57:18"\n', '07:57:18"  # eol on file_date\n')
        + "    # SENTINEL after file_date\n",
    )
    text = _replace_once(
        text,
        ORPHANET_ENTRY_TAIL,
        ORPHANET_ENTRY_TAIL.replace("  - id:", "  # SENTINEL above orphanet_product1\n  - id:"),
    )
    source = _source_as_package(
        "orphanet_product4", version=VERSION, file_version=None, file_date=None
    )
    pin = PackagePin(version=VERSION, manifest_sha256="ab" * 32)

    after = with_source_entry(text.encode(), "orphanet_product4", source, pin)

    expected = _replace_once(
        text,
        OLD_FILE_LINES.replace('07:57:18"\n', '07:57:18"  # eol on file_date\n'),
        f'    version: "{VERSION}"\n',
    )
    expected = _replace_once(
        expected,
        "    redistribution: permitted_with_attribution\n  # SENTINEL above orphanet_product1\n",
        "    redistribution: permitted_with_attribution\n"
        "    package:\n"
        f"      version: {VERSION}\n"
        f"      manifest_sha256: {'ab' * 32}\n"
        "  # SENTINEL above orphanet_product1\n",
    )
    assert after.decode() == expected
    assert with_source_entry(after, "orphanet_product4", source, pin) == after


def test_the_committed_sources_yaml_round_trips_byte_for_byte() -> None:
    """Re-applying an entry's own values changes nothing (the controller's measurement).

    And a field added to the source block is written with no writer change: the
    entry's fields come from the model, not a hand-written list (the drift test
    makes a new `Source` field a `PackageSource` field).
    """
    raw = (REPO_ROOT / "mirrors" / "sources.yaml").read_bytes()
    pin = PackagePin(version="v1", manifest_sha256="cd" * 32)
    tail = f"    package:\n      version: v1\n      manifest_sha256: {'cd' * 32}\n"

    pinned = with_source_entry(raw, "farah2024_hoc", _source_as_package("farah2024_hoc"), pin)

    assert pinned == raw + tail.encode()

    class Wider(PackageSource):
        citation: str | None = None

    source = Wider.model_validate(
        {**_source_as_package("farah2024_hoc").model_dump(), "citation": "Farah 2024"}
    )

    after = with_source_entry(raw, "farah2024_hoc", source, pin)

    assert after == raw + ("    citation: Farah 2024\n" + tail).encode()


# --- what a package may not do: ING013, ING014, and the row accounting ---------


def small_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """`genes.tsv` and `sources.yaml` only, validation stubbed to pass.

    Each test below refuses before anything is written, so validation never
    runs on the refused package; the stub lets the setup ingest succeed, and
    makes a removed guard show as a package that was written.
    """
    repo = tmp_path / "repo"
    (repo / "mirrors").mkdir(parents=True)
    for name in ("genes.tsv", "sources.yaml"):
        shutil.copy(REPO_ROOT / "mirrors" / name, repo / "mirrors" / name)
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport())
    return repo


def mirror_rows(relative: str, gene: str) -> list[dict[str, str]]:
    """A committed mirror's rows for one gene, cells verbatim."""
    header, *lines = (REPO_ROOT / relative).read_text(encoding="utf-8").splitlines()
    columns = header.split("\t")
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in lines]
    return [row for row in rows if row["gene"] == gene]


def write_tables(
    directory: Path, tables: list[tuple[str, str | None, list[dict[str, str]]]], **manifest: Any
) -> Path:
    """A package of ``(table, shard, rows)``, each file in its schema's column order."""
    directory.mkdir(parents=True)
    specs = []
    for name, shard, rows in tables:
        schema = TABLE_SCHEMAS[name]
        data = tsv_bytes(rows, schema.column_names)
        file = f"{name}.tsv" if shard is None else f"{name}-{shard}.tsv"
        (directory / file).write_bytes(data)
        gene = schema.gene_column
        specs.append(
            {
                "name": name,
                "format_version": schema.version,
                "file": file,
                "sha256": hashlib.sha256(data).hexdigest(),
                "rows": len(rows),
                "genes": None if gene is None else len({row[gene] for row in rows}),
                "shard": shard,
            }
        )
    payload = manifest_dict(**{"ids": "hgnc_id", **manifest, "tables": specs})
    (directory / "manifest.json").write_text(json.dumps(payload), encoding="utf-8")
    return directory


PROFILE = "mirrors/profiles/E-MTAB-6814.tsv"


def test_a_version_omitting_a_table_its_pinned_version_carried_is_ing013(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = small_repo(tmp_path, monkeypatch)
    profile = mirror_rows(PROFILE, "HGNC:10023")[:1]
    burden = toy_rows()[:6]
    v1 = write_tables(
        tmp_path / "v1", [("burden", None, burden), ("profiles", "toy_a", profile)], version="v1"
    )
    ingest_package(repo, SOURCE, v1)
    assert (repo / "mirrors" / "profiles" / "toy_a.tsv").is_file()
    before = tree(repo)
    v2 = write_tables(tmp_path / "v2", [("burden", None, burden)], version="v2")

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, SOURCE, v2)

    assert excinfo.value.code == "ING013"
    assert excinfo.value.message == (
        "toy_burden v2 omits what its committed manifest (v1) carried: profiles shard toy_a; "
        "those rows would be left under no pinned manifest"
    )
    assert tree(repo) == before


def test_a_flat_table_without_a_source_column_is_written_only_by_its_owner(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`_merge` replaces such a table whole; `toy_burden` would erase every ClinGen row.

    `test_a_validity_table_is_replaced_whole_and_not_filtered` is the owner's
    own package, which is written.
    """
    repo = small_repo(tmp_path, monkeypatch)
    shutil.copy(REPO_ROOT / "mirrors" / "clingen_gene_validity.tsv", repo / "mirrors")
    rows = mirror_rows("mirrors/clingen_gene_validity.tsv", "HGNC:10023")
    package = write_tables(tmp_path / "pkg", [("clingen_validity", None, rows)])
    before = tree(repo)

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, SOURCE, package)

    assert excinfo.value.code == "ING014"
    assert excinfo.value.message == (
        "table 'clingen_validity' has no source column and belongs to source 'clingen'; "
        "package 'toy_burden' would replace all of it"
    )
    assert tree(repo) == before
    assert apply.FLAT_TABLE_OWNERS == {"clingen_validity": "clingen", "gencc_submissions": "gencc"}
    # Every flat table without a source column has an owner but `genes`, which is ING009.
    unowned = {
        name
        for name, schema in TABLE_SCHEMAS.items()
        if schema.source_column is None and table_layout(name)[0] == "flat"
    } - set(apply.FLAT_TABLE_OWNERS)
    assert unowned == {"genes"}


def test_a_shard_without_a_source_column_another_source_pinned_is_ing014(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = small_repo(tmp_path, monkeypatch)
    profile = mirror_rows(PROFILE, "HGNC:10023")[:1]
    owner = write_tables(tmp_path / "owner", [("profiles", "toy_a", profile)], package="toy_cells")
    ingest_package(repo, "toy_cells", owner)
    before = tree(repo)
    package = write_tables(tmp_path / "pkg", [("profiles", "toy_a", profile)])

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, SOURCE, package)

    assert excinfo.value.code == "ING014"
    assert excinfo.value.message == (
        "profiles shard toy_a has no source column and the committed manifest of toy_cells "
        "carries it; package 'toy_burden' would replace it"
    )
    assert tree(repo) == before


def test_the_registry_filter_drops_a_row_outside_the_registry_and_counts_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The filter alone, with validation stubbed (mutation M10).

    Without the filter, the module-scoped ingest fixture fails in fixture setup because
    validation (BUR011) refuses the stray row; that detects the mutant by an
    error, not by an assertion about the filter. Validation is stubbed here, so
    the only thing keeping HGNC:5 out of `burden.tsv` is the filter.
    """
    repo = small_repo(tmp_path, monkeypatch)
    inside = toy_rows()[:6]
    outside = {**inside[1], "gene": OUTSIDE_GENE}
    assert OUTSIDE_GENE not in registry(repo)
    package = write_tables(tmp_path / "pkg", [("burden", None, [*inside, outside])])

    result = ingest_package(repo, SOURCE, package)

    assert (result.rows_in, result.rows_written, result.rows_outside_registry) == (7, 6, 1)
    lines = (repo / "mirrors" / "burden.tsv").read_text(encoding="utf-8").splitlines()
    genes = [line.split("\t")[1] for line in lines[1:]]
    assert OUTSIDE_GENE not in genes
    assert genes == ["HGNC:10023"] * 6


def test_a_row_dropped_without_being_named_stops_the_ingest_before_any_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The row identity is a guard on `resolve_frame`'s promise: a dropped id is named."""
    repo = small_repo(tmp_path, monkeypatch)
    resolve = apply.resolve_frame

    def lossy(*args: Any) -> Any:
        kept, missed = resolve(*args)
        return kept.head(kept.height - 1), missed

    monkeypatch.setattr(apply, "resolve_frame", lossy)
    package = write_tables(tmp_path / "pkg", [("burden", None, toy_rows()[:6])])
    before = tree(repo)

    with pytest.raises(RuntimeError) as excinfo:
        ingest_package(repo, SOURCE, package)

    assert str(excinfo.value) == (
        "ingest row accounting does not add up: 6 in, 5 written, "
        "0 outside the registry, 0 unresolved"
    )
    assert tree(repo) == before


def test_a_repin_that_keeps_its_version_must_keep_its_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = small_repo(tmp_path, monkeypatch)
    rows = toy_rows()[:6]
    v1 = write_tables(tmp_path / "v1", [("burden", None, rows)], version="v1")
    ingest_package(repo, SOURCE, v1)
    pinned = load_sources(repo)[0]
    pin = next(s.package for s in pinned.sources if s.id == SOURCE)
    assert pin is not None
    before = tree(repo)

    # Same version, different bytes: one row fewer, hashes and counts honest.
    edited = write_tables(tmp_path / "edited", [("burden", None, rows[:5])], version="v1")
    digest = hashlib.sha256((edited / "manifest.json").read_bytes()).hexdigest()
    assert digest != pin.manifest_sha256
    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, SOURCE, edited)

    assert excinfo.value.code == "ING016"
    assert excinfo.value.message == (
        f"toy_burden v1 is pinned at sha256:{pin.manifest_sha256[:12]}\u2026; "
        f"this package's manifest is sha256:{digest[:12]}\u2026; "
        f"a package that changes its bytes must change its version"
    )
    assert tree(repo) == before

    # The identical package is allowed and changes no byte.
    ingest_package(repo, SOURCE, v1)
    assert tree(repo) == before


# --- the hgnc package: `genes.tsv` is the registry the atlas derives ----------

HGNC_VERSION = "2026-10-16"
TBX5, TBX1, RGS9 = "HGNC:11604", "HGNC:11592", "HGNC:10004"
A1BG = {
    "hgnc_id": "HGNC:5",
    "symbol": "A1BG",
    "name": "alpha-1-B glycoprotein",
    "ensembl_gene": "ENSG00000121410",
    "ncbi_gene": "1",
    "locus": "19q13.43",
    "uniprot": "P04217",
    "mane_select": "NM_130786.4",
}


def committed_lines(relative: str, *starts: str) -> list[str]:
    """The header and, in order, the first committed line starting with each prefix."""
    header, *lines = (REPO_ROOT / relative).read_text(encoding="utf-8").splitlines()
    return [header, *(next(line for line in lines if line.startswith(s)) for s in starts)]


def hgnc_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genes: str | None) -> Path:
    """The committed scope, and validity mirrors that derive {TBX5, TBX1}.

    ClinGen names TBX5 for Holt-Oram syndrome (in scope) and RGS9 for
    bradyopsia (not); GenCC names TBX1 for conotruncal heart malformations
    (SGC-102133, in scope) and RGS9 for bradyopsia (SGC-101840). Every line is
    copied from the committed mirrors. ``genes.tsv`` lists ``genes`` alone, or
    is absent when it is None. Validation is stubbed, as in `small_repo`.
    """
    repo = tmp_path / "repo"
    (repo / "curation").mkdir(parents=True)
    (repo / "mirrors").mkdir()
    shutil.copy(REPO_ROOT / "curation" / "chd_scope.yaml", repo / "curation")
    shutil.copy(REPO_ROOT / "mirrors" / "sources.yaml", repo / "mirrors")
    mirrors = {
        "clingen_gene_validity.tsv": (f"{RGS9}\t", f"{TBX5}\tTBX5\tMONDO:0007732\t"),
        "gencc_submissions.tsv": ("SGC-101840\t", "SGC-102133\t"),
    }
    for name, starts in mirrors.items():
        lines = committed_lines(f"mirrors/{name}", *starts)
        (repo / "mirrors" / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    if genes is not None:
        lines = committed_lines("mirrors/genes.tsv", f"{genes}\t")
        (repo / "mirrors" / "genes.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport())
    return repo


def hgnc_package(tmp_path: Path) -> Path:
    """TBX5's and TBX1's committed registry rows, and A1BG, which the registry lacks."""
    header, *lines = committed_lines("mirrors/genes.tsv", f"{TBX5}\t", f"{TBX1}\t")
    columns = header.split("\t")
    rows = [dict(zip(columns, line.split("\t"), strict=True)) for line in lines] + [A1BG]
    return write_tables(
        tmp_path / "pkg", [("genes", None, rows)], package="hgnc", version=HGNC_VERSION
    )


@pytest.mark.parametrize("current", [TBX5, None], ids=["genes-lists-tbx5-only", "no-genes-tsv"])
def test_hgnc_ingest_writes_the_derived_registry_not_the_current_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, current: str | None
) -> None:
    """The filter is the registry derived from the scope, not the `genes.tsv` it replaces.

    TBX1 is in the derived registry and not in the current file, so a filter
    on the current file would drop it (or, with no file, refuse as ING010).
    """
    repo = hgnc_repo(tmp_path, monkeypatch, current)
    package = hgnc_package(tmp_path)

    result = ingest_package(repo, "hgnc", package)

    expected = committed_lines("mirrors/genes.tsv", f"{TBX1}\t", f"{TBX5}\t")
    assert (repo / "mirrors" / "genes.tsv").read_text(encoding="utf-8") == "\n".join(
        expected
    ) + "\n"
    counts = (result.rows_in, result.rows_written, result.rows_outside_registry)
    assert counts == (3, 2, 1)
    assert (result.registry_genes, result.registry_genes_absent) == (2, ())
    cache = repo / "local" / "packages" / "hgnc" / HGNC_VERSION
    assert (cache / "manifest.json").read_bytes() == (package / "manifest.json").read_bytes()
    assert (cache / "genes.tsv").read_bytes() == (package / "genes.tsv").read_bytes()


def _keep_only_rgs9(repo: Path) -> None:
    for name, start in (
        ("clingen_gene_validity.tsv", f"{RGS9}\t"),
        ("gencc_submissions.tsv", "SGC-101840\t"),
    ):
        lines = committed_lines(f"mirrors/{name}", start)
        (repo / "mirrors" / name).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _symbol_for_an_id(repo: Path) -> None:
    path = repo / "mirrors" / "clingen_gene_validity.tsv"
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace(f"{TBX5}\t", "TBX5\t", 1), encoding="utf-8")


def _no_validity_mirrors(repo: Path) -> None:
    for name in ("clingen_gene_validity.tsv", "gencc_submissions.tsv"):
        (repo / "mirrors" / name).unlink()


def _scope_fails_its_model(repo: Path) -> None:
    """An authority outside `ScopeAuthority`: the file parses as YAML, not as `ChdScopeFile`."""
    path = repo / "curation" / "chd_scope.yaml"
    text = path.read_text(encoding="utf-8")
    old = "    admitted_by: clingen_chd_panel\n"
    assert old in text
    path.write_text(text.replace(old, "    admitted_by: the_atlas_itself\n", 1), encoding="utf-8")


def every_file(repo: Path) -> dict[str, bytes | None]:
    """Every file (bytes) and directory (None) under the repository, curation included."""
    return {
        path.relative_to(repo).as_posix(): path.read_bytes() if path.is_file() else None
        for path in sorted(repo.rglob("*"))
    }


_UNDERIVABLE = "the gene registry cannot be derived: "
_EMPTY = (
    "the derived registry is empty: no validity mirror row names a gene against a "
    "disease in curation/chd_scope.yaml, so the filter would drop every row"
)


@pytest.mark.parametrize(
    ("damage", "fragment"),
    [
        (_keep_only_rgs9, _EMPTY),
        # An absent mirror has no rows, so both absent derive nothing.
        (_no_validity_mirrors, _EMPTY),
        (
            lambda repo: (repo / "curation" / "chd_scope.yaml").unlink(),
            _UNDERIVABLE + "curation/chd_scope.yaml is missing",
        ),
        (
            _scope_fails_its_model,
            _UNDERIVABLE + "curation/chd_scope.yaml cannot be read: SCHEMA001 ",
        ),
        (
            _symbol_for_an_id,
            _UNDERIVABLE + "mirrors/clingen_gene_validity.tsv fails its schema",
        ),
    ],
    ids=[
        "derives-nothing",
        "no-validity-mirrors",
        "no-scope",
        "scope-fails-its-model",
        "clingen-fails-its-schema",
    ],
)
def test_an_hgnc_ingest_refuses_a_registry_it_cannot_derive_before_any_write(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    damage: Callable[[Path], None],
    fragment: str,
) -> None:
    """ING010, with every file of the repository byte-unchanged.

    Where the registry cannot be derived, `derive_registry_from_root` itself
    raises `RegistryUnavailable` with the reason the refusal quotes.
    """
    repo = hgnc_repo(tmp_path, monkeypatch, TBX5)
    damage(repo)
    package = hgnc_package(tmp_path)
    before = every_file(repo)

    if fragment.startswith(_UNDERIVABLE):
        with pytest.raises(RegistryUnavailable) as unavailable:
            derive_registry_from_root(repo)
        assert str(unavailable.value).startswith(fragment.removeprefix(_UNDERIVABLE))
    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, "hgnc", package)

    assert excinfo.value.code == "ING010"
    assert excinfo.value.message.startswith(fragment)
    assert every_file(repo) == before


@pytest.mark.parametrize(
    ("absent", "kept"),
    [("gencc_submissions.tsv", TBX5), ("clingen_gene_validity.tsv", TBX1)],
    ids=["gencc-absent", "clingen-absent"],
)
def test_an_hgnc_ingest_with_one_validity_mirror_absent_derives_from_the_other(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, absent: str, kept: str
) -> None:
    """The rule GEN006 applies: an absent mirror has no rows (controller ruling, 2026-11-11).

    So the registry written is the one validation will then check against, and
    a later ingest of the absent mirror that widens it is reported as GEN006.
    """
    repo = hgnc_repo(tmp_path, monkeypatch, TBX5)
    (repo / "mirrors" / absent).unlink()

    result = ingest_package(repo, "hgnc", hgnc_package(tmp_path))

    written = (repo / "mirrors" / "genes.tsv").read_text(encoding="utf-8")
    assert written == "\n".join(committed_lines("mirrors/genes.tsv", f"{kept}\t")) + "\n"
    assert (result.rows_written, result.rows_outside_registry, result.registry_genes) == (1, 2, 1)


# --- ingest --all: every pinned package from its cache, validated once --------

SCOPE = "curation/chd_scope.yaml"
ENO1 = "HGNC:3350"
# The scope term these tests add, then remove. ENO1's one validity row is GenCC's
# PanelApp Australia submission for polymicrogyria (SGC-178789); no other row
# of either validity mirror names MONDO:0000087, and no other committed mirror
# names ENO1 (measured 2026-11-11 with polars over both validity mirrors and
# grep over mirrors/). So this term brings ENO1 alone into the registry, and
# removing it drops ENO1 alone, leaving it with rows only where the two
# ingests below wrote them. Not congenital heart disease: a fixture.
POLYMICROGYRIA = (
    "  - id: MONDO:0000087\n"
    "    label: polymicrogyria\n"
    "    reason: >-\n"
    "      A test fixture's term, not congenital heart disease.\n"
    "    admitted_by: gencc_submitter\n"
    "    attributed_to: PanelApp Australia\n"
    "    admitted_on: 2026-11-11\n"
)
# Only `hgnc_id` and `symbol` matter here; the other cells have HGNC's shapes
# and were not checked against HGNC.
ENO1_ROW = {
    "hgnc_id": ENO1,
    "symbol": "ENO1",
    "name": "enolase 1",
    "ensembl_gene": "ENSG00000074800",
    "ncbi_gene": "2023",
    "locus": "1p36.23",
    "uniprot": "P06733",
    "mane_select": "NM_001428.5",
}
# The committed `hgnc` entry's fields, so the ingest rewrites it as itself.
HGNC_FIELDS = {
    "name": "HUGO Gene Nomenclature Committee",
    "version": HGNC_VERSION,
    "retrieved_on": "2026-10-25",
    "url": "https://www.genenames.org/",
    "licence": "CC0-1.0",
    "redistribution": "permitted",
    "file_version": None,
    "file_date": None,
}


def committed_rows(relative: str) -> list[dict[str, str]]:
    header, *lines = (REPO_ROOT / relative).read_text(encoding="utf-8").splitlines()
    columns = header.split("\t")
    return [dict(zip(columns, line.split("\t"), strict=True)) for line in lines]


class Pinned:
    """The repository with ENO1's disease in scope and `hgnc` and `toy_burden` pinned and cached.

    The hgnc package is every committed `genes.tsv` row plus ENO1's; the burden
    package is `toy_rows()`'s six HGNC:10023 rows and the same six for ENO1.
    Both are ingested with validation stubbed, to save two runs: this state
    validates with 0 errors and the corpus's 3 warnings (SCP003 twice, SCP007;
    measured 2026-11-11, 7.3 s per run under load). Each test copies `repo`.
    """

    def __init__(self, base: Path) -> None:
        self.repo = copy_repo(base / "repo")
        scope = self.repo / SCOPE
        text = scope.read_text(encoding="utf-8")
        assert text.endswith("\n") and POLYMICROGYRIA not in text
        scope.write_text(text + POLYMICROGYRIA, encoding="utf-8")
        genes = [*committed_rows("mirrors/genes.tsv"), ENO1_ROW]
        hgnc = write_tables(
            base / "hgnc",
            [("genes", None, genes)],
            package="hgnc",
            version=HGNC_VERSION,
            source=HGNC_FIELDS,
        )
        burden = toy_rows()[:6]
        toy = write_tables(
            base / "toy", [("burden", None, burden + [{**row, "gene": ENO1} for row in burden])]
        )
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(apply, "validate_repository", lambda _: ValidationReport())
            ingest_package(self.repo, "hgnc", hgnc)
            ingest_package(self.repo, SOURCE, toy)
        assert ENO1 in registry(self.repo)
        assert sum(f"\t{ENO1}\t" in line for line in burden_lines(self.repo)) == 6

    def copy(self, tmp_path: Path) -> Path:
        return Path(shutil.copytree(self.repo, tmp_path / "repo"))


@pytest.fixture(scope="module")
def pinned(tmp_path_factory: pytest.TempPathFactory) -> Pinned:
    return Pinned(tmp_path_factory.mktemp("pinned"))


def shrink(repo: Path) -> None:
    """Remove ENO1's disease from the scope, leaving the committed file."""
    path = repo / SCOPE
    text = path.read_text(encoding="utf-8")
    assert text.count(POLYMICROGYRIA) == 1
    path.write_text(text.replace(POLYMICROGYRIA, ""), encoding="utf-8")
    assert path.read_bytes() == (REPO_ROOT / SCOPE).read_bytes()


def digests(repo: Path, top: str = "mirrors") -> dict[str, str | None]:
    """The sha256 of every file under ``top`` (None for a directory), by relative path."""
    return {
        path.relative_to(repo).as_posix(): (
            hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
        )
        for path in sorted((repo / top).rglob("*"))
    }


def test_ingest_all_survives_a_scope_shrink(pinned: Pinned, tmp_path: Path) -> None:
    """hgnc drops ENO1, the burden package is re-filtered to match, and one validation passes.

    `genes.tsv` comes back byte for byte as committed: the package holds every
    committed row and the shrunk scope derives exactly the committed registry.
    The final validation's warnings are returned once, on the last result.
    """
    repo = pinned.copy(tmp_path)
    shrink(repo)

    results = ingest_all(repo)

    assert [(r.source_id, r.rows_written, r.rows_outside_registry) for r in results] == [
        ("hgnc", 185, 1),
        (SOURCE, 6, 6),
    ]
    genes = repo / "mirrors" / "genes.tsv"
    assert genes.read_bytes() == (REPO_ROOT / "mirrors" / "genes.tsv").read_bytes()
    assert ENO1 not in registry(repo)
    lines = burden_lines(repo)
    assert not [line for line in lines if f"\t{ENO1}\t" in line]
    toy = [line for line in lines if line.endswith(f"\t{SOURCE}")]
    assert len(toy) == 6 and all("\tHGNC:10023\t" in line for line in toy)
    report = validate_repository(repo)
    assert [issue.format() for issue in report.issues if issue.severity is Severity.ERROR] == []
    warnings = tuple(issue for issue in report.issues if issue.severity is Severity.WARNING)
    assert warnings
    assert [r.warnings for r in results] == [(), warnings]


ALL_HINT = (
    "If curation/chd_scope.yaml dropped a gene that another pinned package still has rows "
    "for, `chd-atlas ingest --all` re-filters every pinned package and validates once."
)


def test_a_single_hgnc_ingest_after_a_shrink_refuses_and_names_the_all_flag(
    pinned: Pinned, tmp_path: Path
) -> None:
    """The deadlock `ingest --all` exists for: hgnc alone drops ENO1 under the burden rows."""
    repo = pinned.copy(tmp_path)
    shrink(repo)
    before = digests(repo)

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, "hgnc", package_cache(repo, "hgnc", HGNC_VERSION))

    assert excinfo.value.code == "ING010"
    assert excinfo.value.message == (
        "the repository fails validation with this package applied (1 error(s)):\n"
        f"ERROR [BUR011] mirrors/burden.tsv: gene {ENO1} is not in mirrors/genes.tsv\n" + ALL_HINT
    )
    assert digests(repo) == before


def _corrupt_the_burden_cache(repo: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    path = package_cache(repo, SOURCE, VERSION) / "burden.tsv"
    data = bytearray(path.read_bytes())
    data[-2] ^= 1  # one byte of the last row's `source` cell
    path.write_bytes(bytes(data))


ERROR = ValidationIssue("BUR011", Severity.ERROR, "mirrors/burden.tsv", "injected")


def _fail_the_final_validation(repo: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport([ERROR]))


def _fail_it_with_the_manifests_directory_missing(
    repo: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """hgnc's step then creates `mirrors/packages/`, which the restore must remove."""
    shutil.rmtree(repo / "mirrors" / "packages")
    _fail_the_final_validation(repo, monkeypatch)


def _interrupt_the_final_validation(repo: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def interrupt(_: Path) -> ValidationReport:
        raise KeyboardInterrupt

    monkeypatch.setattr(apply, "validate_repository", interrupt)


@pytest.mark.parametrize(
    ("fail", "raised", "code"),
    [
        (_corrupt_the_burden_cache, IngestRefused, "ING001"),
        (_fail_the_final_validation, IngestRefused, "ING010"),
        (_fail_it_with_the_manifests_directory_missing, IngestRefused, "ING010"),
        (_interrupt_the_final_validation, KeyboardInterrupt, None),
    ],
    ids=["second-cache-corrupted", "final-validation-fails", "a-created-directory", "interrupted"],
)
def test_ingest_all_restores_mirrors_when_any_step_fails(
    pinned: Pinned,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fail: Callable[[Path, pytest.MonkeyPatch], None],
    raised: type[BaseException],
    code: str | None,
) -> None:
    """After a shrink, hgnc has already rewritten `genes.tsv` when the failure lands.

    So every byte of `mirrors/` the restore puts back is one an earlier step
    changed, which `_apply`'s per-step restore cannot see; with the manifests
    directory missing, the steps also create a directory and two files there.
    """
    repo = pinned.copy(tmp_path)
    shrink(repo)
    fail(repo, monkeypatch)
    before = digests(repo)
    cache = digests(repo, "local")
    written = recording_writes(monkeypatch)

    with pytest.raises(raised) as excinfo:
        ingest_all(repo)

    if code is not None:
        assert isinstance(excinfo.value, IngestRefused)
        assert excinfo.value.code == code
    if code == "ING010":
        assert excinfo.value.message == (
            "the repository fails validation with every pinned package re-ingested "
            "(1 error(s)):\n" + ERROR.format()
        )
    assert written[0] == repo / "mirrors" / "genes.tsv"
    assert digests(repo) == before
    assert digests(repo, "local") == cache


@pytest.mark.parametrize("source_id", ["hgnc", SOURCE])
def test_reingesting_from_the_cache_directory_changes_nothing(
    pinned: Pinned, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, source_id: str
) -> None:
    """`--package local/packages/<id>/<version>` is the directory ingest also writes to.

    Validation is stubbed: it writes nothing, and what is asserted is that no
    byte under `mirrors/` or `local/` moved.
    """
    repo = pinned.copy(tmp_path)
    before = tree(repo)
    monkeypatch.setattr(apply, "validate_repository", lambda _: ValidationReport())
    version = HGNC_VERSION if source_id == "hgnc" else VERSION

    ingest_package(repo, source_id, package_cache(repo, source_id, version))

    assert tree(repo) == before


def add_pin(
    repo: Path,
    source_id: str,
    *,
    version: str = "v1",
    manifest: bytes | None = b"{}",
    sha256: str | None = None,
) -> None:
    """Append ``source_id`` pinned at ``version`` to `sources.yaml`.

    ``manifest`` is cached as the package's `manifest.json` unless None; the
    pin is its sha256 unless ``sha256`` is given.
    """
    if manifest is not None:
        cache = package_cache(repo, source_id, version)
        cache.mkdir(parents=True)
        (cache / "manifest.json").write_bytes(manifest)
    digest = sha256 or hashlib.sha256(manifest or b"").hexdigest()
    with (repo / "mirrors" / "sources.yaml").open("a", encoding="utf-8") as handle:
        handle.write(
            f"  - id: {source_id}\n"
            f"    name: {source_id}\n"
            f'    version: "{version}"\n'
            "    retrieved_on: 2026-11-01\n"
            "    url: https://example.org/\n"
            "    licence: CC-BY-4.0\n"
            "    redistribution: permitted\n"
            f'    package: {{version: "{version}", manifest_sha256: "{digest}"}}\n'
        )


def _step_result(source_id: str) -> IngestResult:
    return IngestResult(
        source_id=source_id,
        version="v1",
        rows_in=0,
        rows_written=0,
        rows_outside_registry=0,
        rows_unresolved=0,
        unresolved=(),
        registry_genes_absent=(),
        registry_genes=0,
        coverage="every_gene",
        dropped=(),
        warnings=(),
    )


def test_ingest_all_takes_the_registry_inputs_then_hgnc_then_the_rest_validating_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """ClinGen and GenCC first (hgnc's filter is derived from them), then hgnc, then by id.

    `aaa_source` sorts before `hgnc` and `clingen`, so neither "sorted" nor
    "hgnc first" passes. Every step runs unvalidated, from its pinned cache.
    """
    repo = write_hgnc_root(tmp_path, [A1BG], version=HGNC_VERSION)
    for source_id in (SOURCE, "gencc", "aaa_source", "clingen"):
        add_pin(repo, source_id)
    calls: list[tuple[str, Path, bool]] = []

    def step(root: Path, source_id: str, package_dir: Path, *, validate: bool = True) -> Any:
        calls.append((source_id, package_dir, validate))
        return _step_result(source_id)

    validated: list[Path] = []
    warning = ValidationIssue("SCP003", Severity.WARNING, SCOPE, "msg")

    def validate(root: Path) -> ValidationReport:
        validated.append(root)
        return ValidationReport([warning])

    monkeypatch.setattr(apply, "ingest_package", step)
    monkeypatch.setattr(apply, "validate_repository", validate)

    results = ingest_all(repo)

    order = [
        ("clingen", "v1"),
        ("gencc", "v1"),
        ("hgnc", HGNC_VERSION),
        ("aaa_source", "v1"),
        (SOURCE, "v1"),
    ]
    assert calls == [(s, package_cache(repo, s, v), False) for s, v in order]
    assert validated == [repo]
    assert [(r.source_id, r.warnings) for r in results] == [
        ("clingen", ()),
        ("gencc", ()),
        ("hgnc", ()),
        ("aaa_source", ()),
        (SOURCE, (warning,)),
    ]


def _never_stepped(*_: object, **__: object) -> Any:
    raise AssertionError("no package may be ingested while a pinned cache is unusable")


@pytest.mark.parametrize(
    ("cache", "detail"),
    [
        (
            {"manifest": None},
            "toy_burden v1 is pinned in mirrors/sources.yaml but not cached at "
            "local/packages/toy_burden/v1; download it: ",
        ),
        (
            {"sha256": "f" * 64},
            "local/packages/toy_burden/v1/manifest.json is not the manifest toy_burden v1 "
            "is pinned at (sha256:ffffffffffff…); remove the directory and download it "
            "again: ",
        ),
    ],
    ids=["missing", "not-the-pinned-manifest"],
)
def test_a_missing_cache_names_the_download_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cache: dict[str, Any], detail: str
) -> None:
    """ING017, checked for every pin before any package is ingested.

    hgnc's cache is present and sorts first, so a check made step by step
    would ingest it before reaching `toy_burden`.
    """
    repo = write_hgnc_root(tmp_path, [A1BG], version=HGNC_VERSION)
    add_pin(repo, SOURCE, **cache)
    monkeypatch.setattr(apply, "ingest_package", _never_stepped)
    before = every_file(repo)

    with pytest.raises(IngestRefused) as excinfo:
        ingest_all(repo)

    assert excinfo.value.code == "ING017"
    assert excinfo.value.message == detail + (
        "gh release download toy_burden/v1 -R chd-atlas/feeds -D local/packages/toy_burden/v1"
    )
    assert every_file(repo) == before


def test_ingest_all_with_nothing_pinned_ingests_and_validates_nothing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = write_hgnc_root(tmp_path, [A1BG], pinned=False)
    monkeypatch.setattr(apply, "ingest_package", _never_stepped)
    monkeypatch.setattr(apply, "validate_repository", _never_validated)

    assert ingest_all(repo) == ()
