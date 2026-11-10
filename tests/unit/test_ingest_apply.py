"""`ingest_package`: a verified package merged into a copy of the real repository.

Every test that writes works on a temporary copy of `curation/`, `mirrors/` and
`ontologies/` (as `tests/test_repository_validates.py` does), never on the
checkout. `validate_repository` over the real corpus costs seconds, so the
happy path is ingested once per module, twice (the second time to prove a
re-ingest changes nothing), and its tests read the snapshots the fixture took.

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
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest import apply
from chd_atlas.ingest.apply import ingest_package, with_source_entry
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.models.package import PackageSource, RegistrySlice
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, table_layout
from chd_atlas.validate.runner import ValidationReport
from chd_atlas.validate.sort_order import validate_sort_order
from chd_atlas.validate.sources import PackagePin, load_sources
from tests.unit.package_fixtures import manifest_dict, tsv_bytes, write_package

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


def test_rows_outside_the_registry_are_counted_not_written(ingested: Ingested) -> None:
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


def test_a_registry_slice_names_the_registry_genes_it_lacks(ingested: Ingested) -> None:
    result = ingested.first
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
    """Re-applying an entry's own values changes nothing (the controller's measurement)."""
    raw = (REPO_ROOT / "mirrors" / "sources.yaml").read_bytes()
    pinned = with_source_entry(
        raw,
        "farah2024_hoc",
        _source_as_package("farah2024_hoc"),
        PackagePin(version="v1", manifest_sha256="cd" * 32),
    )
    tail = f"    package:\n      version: v1\n      manifest_sha256: {'cd' * 32}\n"
    assert pinned == raw + tail.encode()


def test_a_field_added_to_the_source_block_is_written_with_no_writer_change() -> None:
    """The entry's fields come from the model, not a hand-written list (the drift
    test makes a new `Source` field a `PackageSource` field)."""

    class Wider(PackageSource):
        citation: str | None = None

    source = Wider.model_validate(
        {**_source_as_package("farah2024_hoc").model_dump(), "citation": "Farah 2024"}
    )
    raw = (REPO_ROOT / "mirrors" / "sources.yaml").read_bytes()
    pin = PackagePin(version="v1", manifest_sha256="cd" * 32)

    after = with_source_entry(raw, "farah2024_hoc", source, pin)

    tail = (
        "    citation: Farah 2024\n"
        f"    package:\n      version: v1\n      manifest_sha256: {'cd' * 32}\n"
    )
    assert after == raw + tail.encode()


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
                "format_version": 1,
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
