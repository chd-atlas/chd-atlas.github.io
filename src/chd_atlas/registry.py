"""The atlas's gene registry, derived from its scope and the two validity mirrors.

The registry is every gene that a row of `mirrors/clingen_gene_validity.tsv`
or `mirrors/gencc_submissions.tsv` names against a disease listed in
`curation/chd_scope.yaml`, at any classification. It is deliberately wider than
the genes the site publishes: a lookup table, so that widening the publication
gate does not require re-mirroring.

`scripts/convert_hgnc.py::_in_scope_genes` held this rule until feeds phase 2,
reading the scope ids with a regular expression over the YAML text. Here they
are `ChdScopeFile.diseases[].id`; measured 2026-11-10, the file's 105 ``- id:``
lines are all MONDO ids under ``diseases``, so the two readings agree, and the
derivation over the committed corpus is the 185 genes of `mirrors/genes.tsv`
(`tests/unit/test_registry.py` pins both).

Two readers use it: `chd-atlas ingest hgnc` writes `mirrors/genes.tsv` filtered
to it, and `validate_repository` reports a gene it holds that `genes.tsv`
lacks (GEN006).
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Final

import polars as pl

from chd_atlas.corpus import _Accumulator, _parse
from chd_atlas.issues import Severity
from chd_atlas.models.scope import ChdScopeFile
from chd_atlas.tables import CLINGEN_VALIDITY, GENCC_SUBMISSIONS, read_table, validate_table

SCOPE_FILE: Final = "curation/chd_scope.yaml"
_MIRRORS: Final = (
    ("mirrors/clingen_gene_validity.tsv", CLINGEN_VALIDITY),
    ("mirrors/gencc_submissions.tsv", GENCC_SUBMISSIONS),
)


class RegistryUnavailable(Exception):
    """The registry cannot be derived; the message says which input failed and how."""


def derive_registry(
    scope_ids: Iterable[str], clingen: pl.DataFrame, gencc: pl.DataFrame
) -> frozenset[str]:
    """Every gene either mirror names against a disease in ``scope_ids``.

    The frames are the two validity mirrors as `read_table` returns them. A row
    with no gene names nothing. A frame without a ``gene`` or ``disease``
    column raises `RegistryUnavailable` rather than a polars error, so a caller
    that reports instead of raising can catch one exception.
    """
    scope = sorted(set(scope_ids))
    genes: set[str] = set()
    for name, frame in (("clingen_validity", clingen), ("gencc_submissions", gencc)):
        missing = [column for column in ("gene", "disease") if column not in frame.columns]
        if missing:
            raise RegistryUnavailable(f"the {name} mirror has no {' or '.join(missing)} column")
        named = frame.filter(pl.col("disease").is_in(scope))["gene"].to_list()
        genes.update(str(gene) for gene in named if gene is not None)
    return frozenset(genes)


def derive_registry_from_root(root: Path) -> frozenset[str]:
    """`derive_registry` over ``root``'s scope file and validity mirrors.

    Raises `RegistryUnavailable` when the scope file is missing, unreadable or
    fails `ChdScopeFile`, or when either mirror is missing or reports an error
    against its table schema (`validate_table`, which also covers a file
    polars cannot read). Measured 2026-11-10: `validate_table` over both
    committed mirrors takes about 0.16 s and reports nothing.
    """
    scope_path = root / SCOPE_FILE
    if not scope_path.is_file():
        raise RegistryUnavailable(f"{SCOPE_FILE} is missing")
    acc = _Accumulator()
    scope = _parse(ChdScopeFile, scope_path, acc)
    if scope is None:
        problems = "; ".join(f"{issue.code} {issue.message}" for issue in acc.issues[:3])
        raise RegistryUnavailable(f"{SCOPE_FILE} cannot be read: {problems}")
    frames: list[pl.DataFrame] = []
    for relative, schema in _MIRRORS:
        path = root / relative
        if not path.is_file():
            raise RegistryUnavailable(f"{relative} is missing")
        errors = [i for i in validate_table(path, schema) if i.severity is Severity.ERROR]
        if errors:
            raise RegistryUnavailable(
                f"{relative} fails its schema ({len(errors)} error(s); first: "
                f"{errors[0].code} {errors[0].message})"
            )
        frame, _ = read_table(path, schema)
        if frame is None:
            # Not reached: `validate_table` reads the file with `read_table` and
            # reports a failed read as TBL000. Kept for the type, as a raise.
            raise RegistryUnavailable(f"{relative} cannot be read")
        frames.append(frame)
    clingen, gencc = frames
    return derive_registry((str(entry.id) for entry in scope.diseases), clingen, gencc)
