"""`chd_atlas.registry`: the gene registry the atlas derives from its scope and mirrors."""

from __future__ import annotations

from pathlib import Path

import polars as pl

from chd_atlas.registry import derive_registry, derive_registry_from_root
from chd_atlas.tables import TABLE_SCHEMAS, read_table

REPO_ROOT = Path(__file__).parent.parent.parent


def test_the_derived_registry_is_the_committed_one() -> None:
    """The derivation reproduces `mirrors/genes.tsv`'s 185 genes exactly.

    Measured 2026-11-10 on the committed corpus (105 scope terms): 185 genes,
    equal to the registry `scripts/convert_hgnc.py` wrote. The ClinGen mirror
    names 134 of them against an in-scope disease and the GenCC mirror all 185,
    so 51 come from GenCC alone and none from ClinGen alone: on this corpus a
    derivation that ignored ClinGen would still pass here, which is why the
    synthetic test below keeps a ClinGen-only gene.
    """
    genes, _ = read_table(REPO_ROOT / "mirrors" / "genes.tsv", TABLE_SCHEMAS["genes"])
    assert genes is not None
    committed = frozenset(genes["hgnc_id"].to_list())

    derived = derive_registry_from_root(REPO_ROOT)

    assert len(committed) == 185
    assert derived == committed


def _validity(rows: list[tuple[str | None, str]]) -> pl.DataFrame:
    return pl.DataFrame(
        {"gene": [gene for gene, _ in rows], "disease": [disease for _, disease in rows]},
        schema={"gene": pl.String, "disease": pl.String},
    )


def test_the_registry_counts_only_in_scope_diseases_from_both_mirrors() -> None:
    """A is ClinGen in scope, B is GenCC in scope, C is named in both, never in scope.

    A row with no gene names nothing, and the scope is read once from an iterator.
    """
    in_scope, out_of_scope = "MONDO:0007732", "MONDO:0012033"
    clingen = _validity([("HGNC:A", in_scope), ("HGNC:C", out_of_scope), (None, in_scope)])
    gencc = _validity([("HGNC:B", in_scope), ("HGNC:C", out_of_scope)])

    derived = derive_registry(iter([in_scope, "MONDO:0000001"]), clingen, gencc)

    assert derived == frozenset({"HGNC:A", "HGNC:B"})
    assert isinstance(derived, frozenset)
