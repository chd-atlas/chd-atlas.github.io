from __future__ import annotations

import json
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.genes import Resolution, Resolved
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.ingest.resolve import GeneIdResolver, package_cache, resolve_frame
from chd_atlas.vocab import GeneIdKind
from tests.unit.package_fixtures import write_hgnc_root

# Synthetic HGNC rows. The ids and the ODD / MADH7 / BVES facts are the real ones
# (CLAUDE.md section 10, #33); the Ensembl ids are made up.
HGNC_ROWS = [
    {"hgnc_id": "HGNC:4274", "symbol": "GJA1", "aliases": "ODD|CX43", "ensembl_gene": "ENSG1"},
    {"hgnc_id": "HGNC:8111", "symbol": "OSR1", "prev_symbols": "ODD", "ensembl_gene": "ENSG2"},
    {"hgnc_id": "HGNC:6772", "symbol": "SMAD6", "prev_symbols": "MADH7", "ensembl_gene": "ENSG3"},
    {"hgnc_id": "HGNC:6773", "symbol": "SMAD7", "prev_symbols": "MADH7", "ensembl_gene": "ENSG3"},
    {"hgnc_id": "HGNC:11604", "symbol": "TBX5", "aliases": "HOS"},
    {"hgnc_id": "HGNC:1152", "symbol": "POPDC1", "prev_symbols": "BVES"},
]


def resolver(
    tmp_path: Path, kind: GeneIdKind = GeneIdKind.SYMBOL, **kwargs: bool
) -> GeneIdResolver:
    write_hgnc_root(tmp_path, HGNC_ROWS, **kwargs)
    return GeneIdResolver.from_root(tmp_path, kind)


def frame(genes: list[str | None]) -> pl.DataFrame:
    return pl.DataFrame(
        {"gene": genes, "n": list(range(len(genes)))}, schema_overrides={"gene": pl.String}
    )


def test_the_cache_path_is_under_local_packages(tmp_path: Path) -> None:
    assert (
        package_cache(tmp_path, "hgnc", "2026-10-16") == tmp_path / "local/packages/hgnc/2026-10-16"
    )


@pytest.mark.parametrize("gz", [False, True])
def test_symbols_resolve_with_the_registry_precedence(tmp_path: Path, gz: bool) -> None:
    # ODD is GJA1's live alias and OSR1's retired name: the alias wins. BVES is
    # only ever a retired name. A gz-compressed genes table reads the same.
    rows, unresolved = resolve_frame(
        frame(["ODD", "BVES", "TBX5", "HOS"]), "gene", resolver(tmp_path, gz=gz)
    )
    assert rows["gene"].to_list() == ["HGNC:4274", "HGNC:1152", "HGNC:11604", "HGNC:11604"]
    assert unresolved == ()


def test_an_ambiguous_symbol_drops_its_rows_and_is_named(tmp_path: Path) -> None:
    rows, unresolved = resolve_frame(
        frame(["TBX5", "MADH7", "MADH7", "NOSUCH"]), "gene", resolver(tmp_path)
    )
    assert rows["gene"].to_list() == ["HGNC:11604"]
    assert unresolved == (
        Resolved("MADH7", Resolution.AMBIGUOUS, None, ("HGNC:6772", "HGNC:6773")),
        Resolved("NOSUCH", Resolution.UNKNOWN),
    )


def test_an_ensembl_id_two_genes_claim_is_ambiguous(tmp_path: Path) -> None:
    ensembl = resolver(tmp_path, GeneIdKind.ENSEMBL_GENE_ID)
    assert ensembl.resolve("ENSG3") == Resolved(
        "ENSG3", Resolution.AMBIGUOUS, None, ("HGNC:6772", "HGNC:6773")
    )
    # TBX5 has a blank ensembl_gene and must not be indexed under "";
    # a version suffix is not stripped, so ENSG1.4 is not ENSG1.
    rows, unresolved = resolve_frame(
        frame(["ENSG1", "ENSG3", "ENSG1.4", " ENSG2 "]), "gene", ensembl
    )
    assert rows["gene"].to_list() == ["HGNC:4274", "HGNC:8111"]
    assert [item.symbol for item in unresolved] == ["ENSG1.4", "ENSG3"]
    assert ensembl.resolve("").resolution is Resolution.UNKNOWN


def test_blank_and_padded_ids(tmp_path: Path) -> None:
    # Review Focus 5: padded ids resolve, blanks are listed (once) as <blank>,
    # a null cell is a blank, and no row is silently lost: kept + dropped = all.
    rows, unresolved = resolve_frame(
        frame([" TBX5 ", "", "  ", None, "TBX5"]), "gene", resolver(tmp_path)
    )
    assert rows["gene"].to_list() == ["HGNC:11604", "HGNC:11604"]
    assert rows["n"].to_list() == [0, 4]
    assert unresolved == (Resolved("<blank>", Resolution.UNKNOWN),)


def test_hgnc_ids_pass_through_and_a_malformed_one_is_unknown(tmp_path: Path) -> None:
    # No HGNC package is needed (or read) for an HGNC-keyed package, and an id
    # is not checked against the registry here: Task 7 filters to it next.
    identity = GeneIdResolver.from_root(tmp_path, GeneIdKind.HGNC_ID)
    rows, unresolved = resolve_frame(
        frame(["HGNC:99999999", " HGNC:11604 ", "HGNC:", "hgnc:1", "TBX5"]), "gene", identity
    )
    assert rows["gene"].to_list() == ["HGNC:99999999", "HGNC:11604"]
    assert [item.symbol for item in unresolved] == ["HGNC:", "TBX5", "hgnc:1"]
    assert {item.resolution for item in unresolved} == {Resolution.UNKNOWN}


@pytest.mark.parametrize(
    "damage", ["no_pin", "no_cache", "no_manifest", "no_genes_table", "bad_registry"]
)
def test_an_unusable_hgnc_package_is_ing008(tmp_path: Path, damage: str) -> None:
    # A committed mirrors/genes.tsv exists and must NOT be used as a fallback.
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "genes.tsv").write_text(
        "hgnc_id\tsymbol\nHGNC:11604\tTBX5\n", encoding="utf-8"
    )
    write_hgnc_root(tmp_path, HGNC_ROWS, pinned=damage != "no_pin", cached=damage != "no_cache")
    cache = package_cache(tmp_path, "hgnc", "2026-10-16")
    if damage == "no_manifest":
        (cache / "manifest.json").write_text("{not json", encoding="utf-8")
    if damage == "no_genes_table":
        manifest = json.loads((cache / "manifest.json").read_text(encoding="utf-8"))
        manifest["tables"][0]["name"] = "burden"
        (cache / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    if damage == "bad_registry":
        (tmp_path / "mirrors" / "sources.yaml").write_text("sources: []\n", encoding="utf-8")
    with pytest.raises(IngestRefused) as excinfo:
        GeneIdResolver.from_root(tmp_path, GeneIdKind.SYMBOL)
    assert excinfo.value.code == "ING008"
    if damage in {"no_pin", "no_cache"}:
        assert "HGNC package is ingested first" in excinfo.value.message
