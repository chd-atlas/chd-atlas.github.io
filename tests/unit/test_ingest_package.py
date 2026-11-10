from __future__ import annotations

import gzip
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.ingest.package import IngestRefused, load_package
from chd_atlas.tables import TABLE_SCHEMAS
from tests.unit.package_fixtures import (
    burden_row,
    default_burden_rows,
    edit_table,
    replace_table_bytes,
    tsv_bytes,
    write_package,
)

SOURCE = "toy_burden"
BURDEN_COLUMNS = TABLE_SCHEMAS["burden"].column_names


def refusal(directory: Path, source_id: str = SOURCE) -> IngestRefused:
    with pytest.raises(IngestRefused) as excinfo:
        load_package(directory, source_id)
    return excinfo.value


def test_a_valid_package_loads_every_cell_as_text(tmp_path: Path) -> None:
    """Written with its header reversed: ING012 compares column sets, not order (Task 5)."""
    rows = default_burden_rows()
    rows[2]["maf_max"] = "0.001"
    columns = list(reversed(BURDEN_COLUMNS))
    package = load_package(write_package(tmp_path, rows, columns=columns), SOURCE)
    (table,) = package.tables
    assert table.spec.name == "burden"
    assert table.schema.name == "burden"
    assert table.frame.height == 3
    assert set(table.frame.dtypes) == {pl.String}
    assert table.frame["maf_max"].to_list() == [None, None, "0.001"]
    assert table.frame["n_cases"].to_list() == ["3887", "3887", "22558"]
    assert package.manifest.package == SOURCE
    assert package.manifest_bytes == (package.directory / "manifest.json").read_bytes()


def test_a_gzipped_table_is_checked_compressed_and_counted_decompressed(tmp_path: Path) -> None:
    directory = write_package(tmp_path, gz=True)
    package = load_package(directory, SOURCE)
    assert package.tables[0].frame.height == 3
    on_disk = (directory / "burden.tsv.gz").read_bytes()
    assert on_disk[:2] == b"\x1f\x8b"
    # The recorded sha256 is of the compressed bytes: re-compressing the same
    # rows with a different level changes the bytes on disk and must refuse.
    plain = gzip.decompress(on_disk)
    (directory / "burden.tsv.gz").write_bytes(gzip.compress(plain, compresslevel=1, mtime=7))
    assert refusal(directory).code == "ING001"


def test_ING004_manifest_missing_or_invalid(tmp_path: Path) -> None:
    """A NUL in a table's file escaped `load_package` as `Path.read_bytes`' ValueError;
    the manifest refuses it, and the refusal names the field."""
    directory = write_package(tmp_path)
    edit_table(directory, file="burden\x00.tsv")
    refused = refusal(directory)
    assert refused.code == "ING004"
    assert "tables.0.file" in refused.message
    (directory / "manifest.json").write_text("{not json", encoding="utf-8")
    assert refusal(directory).code == "ING004"
    (directory / "manifest.json").unlink()
    error = refusal(directory)
    assert error.code == "ING004"
    assert str(error).startswith("ING004: ")


def test_ING006_package_id_differs_from_the_source(tmp_path: Path) -> None:
    assert refusal(write_package(tmp_path), "another_source").code == "ING006"


@pytest.mark.parametrize(
    ("fields", "fragment"),
    [
        ({"name": "no_such_table"}, "not a table"),
        ({"name": "cell_type_profiles", "shard": None}, "sharded"),
        ({"shard": "hoc"}, "flat"),
        ({"name": "orphanet_omim", "genes": 2}, "no gene column"),
        ({"genes": None}, "has a gene column"),
    ],
    ids=[
        "unknown",
        "sharded-without-shard",
        "flat-with-shard",
        "count-on-geneless",
        "gene-keyed-no-count",
    ],
)
def test_ING005_each_shape_condition_reaches_its_own_refusal(
    tmp_path: Path, fields: dict[str, object], fragment: str
) -> None:
    directory = write_package(tmp_path)
    edit_table(directory, **fields)
    error = refusal(directory)
    assert error.code == "ING005"
    assert fragment in error.message


@pytest.mark.parametrize(
    ("fields", "code"),
    [
        # The genes table is refused before its version or its bytes are read.
        ({"name": "genes", "genes": 2, "format_version": 9, "sha256": "0" * 64}, "ING009"),
        ({"format_version": 2}, "ING002"),
        ({"rows": 99}, "ING003"),
        ({"genes": 99}, "ING003"),
    ],
    ids=["ING009-before-version-and-bytes", "ING002-format-version", "ING003-rows", "ING003-genes"],
)
def test_a_manifest_entry_the_table_does_not_honour_is_its_own_code(
    tmp_path: Path, fields: dict[str, object], code: str
) -> None:
    directory = write_package(tmp_path)
    edit_table(directory, **fields)
    assert refusal(directory).code == code


def test_ING001_table_file_absent_from_the_package(tmp_path: Path) -> None:
    directory = write_package(tmp_path)
    (directory / "burden.tsv").unlink()
    assert refusal(directory).code == "ING001"


@pytest.mark.parametrize("gz", [False, True])
def test_ING011_carriage_return_is_checked_on_decompressed_bytes(tmp_path: Path, gz: bool) -> None:
    directory = write_package(tmp_path, gz=gz)
    raw = tsv_bytes(default_burden_rows(), BURDEN_COLUMNS).replace(b"\n", b"\r\n")
    replace_table_bytes(directory, gzip.compress(raw, mtime=0) if gz else raw)
    assert refusal(directory).code == "ING011"


@pytest.mark.parametrize(
    ("gz", "body"),
    [
        (True, b"not gzip at all"),
        (True, gzip.compress(b"a\tb\n")[:-8]),
        (False, "x".encode("utf-16")),
        (False, tsv_bytes([burden_row()], BURDEN_COLUMNS) + b"HGNC:1\textra\n"),
        (False, b""),
        (False, b"stud\xffy\t" + tsv_bytes([burden_row()], BURDEN_COLUMNS)[6:]),
    ],
    ids=["bad-gzip", "truncated-gzip", "not-utf8", "ragged-row", "empty", "bad-byte-in-header"],
)
def test_ING011_every_unreadable_body_is_refused_not_raised(
    tmp_path: Path, gz: bool, body: bytes
) -> None:
    directory = write_package(tmp_path, gz=gz)
    replace_table_bytes(directory, body)
    assert refusal(directory).code == "ING011"


def test_ING012_names_the_missing_and_the_extra_columns(tmp_path: Path) -> None:
    columns = [c for c in BURDEN_COLUMNS if c != "effect"] + ["surprise"]
    error = refusal(write_package(tmp_path, columns=columns))
    assert error.code == "ING012"
    assert "effect" in error.message
    assert "surprise" in error.message


def test_ING003_a_blank_gene_cell_is_not_a_gene(tmp_path: Path) -> None:
    rows = [*default_burden_rows(), burden_row(gene="  ", consequence_class="lof")]
    directory = write_package(tmp_path, rows)
    # write_package counted two genes (blank excluded) and four rows.
    assert load_package(directory, SOURCE).tables[0].frame.height == 4
    edit_table(directory, genes=3)
    assert refusal(directory).code == "ING003"


@pytest.mark.parametrize(
    ("source", "named"), [("someone_else", "'someone_else'"), ("", "'<blank>'"), (" ", "'<blank>'")]
)
def test_ING007_a_row_from_another_source(tmp_path: Path, source: str, named: str) -> None:
    """A blank source is named `<blank>`, as an unresolved blank gene id is, never 'None'."""
    rows = default_burden_rows()
    rows[1]["source"] = source
    refused = refusal(write_package(tmp_path, rows))
    assert refused.code == "ING007"
    assert refused.message == f"table 'burden' has rows whose source is not 'toy_burden': [{named}]"


def test_a_package_cell_that_opens_with_a_quote_is_refused(tmp_path: Path) -> None:
    rows = default_burden_rows()
    rows[1]["method_note"] = '"""open x"'  # CSV-quoted: parses to "open x
    refused = refusal(write_package(tmp_path, rows))
    assert refused.code == "ING015"
    assert "burden" in refused.message
    assert "method_note" in refused.message
    assert "1 cell " in refused.message
    assert "data row 2" in refused.message
