"""The tables `chd-atlas ingest` keeps only part of, and what it keeps (spec F23).

Two packages carry their whole source where the atlas mirrors a slice:

- ``orphanet_phenotypes`` (Orphanet product 4, every phenotype annotation):
  ingest keeps a row whose ``hpo`` is under HP:0001626 (``CARDIOVASCULAR_ROOT``)
  by ``is_a``, alt ids included, in the HPO release ``mirrors/sources.yaml``
  pins for ``HP``, the release SCP006 and the build read. A re-pinned release
  therefore reaches the mirror at the next ingest of the package, which
  ``ingest --all`` does from its cache; until then the mirror was filtered
  through the release pinned when it was last ingested.
  A row whose ``hpo`` the release does not define at all (neither an ``id``
  nor an ``alt_id`` of a ``[Term]``) is neither kept nor counted as outside the
  filter: whether it is cardiovascular is unknown, so its term is named
  (`IngestResult.unknown_terms`) and its rows are not written. A blank cell is
  named ``<blank>``.
- ``mondo_xrefs`` (Mondo SSSOM, every exact match to Orphanet or OMIM): ingest
  keeps a row whose ``mondo`` is the ``disease`` of any row of either validity
  mirror, as they stand when the filter is built. ``ingest --all`` applies
  ``mondo_sssom`` after ClinGen and GenCC (`apply.ingest_order`), so a validity
  refresh re-filters Mondo in the same transaction.

Both are the rules the retired `scripts/convert_orphanet.py` applied
(`_cardiovascular`, `_mirrored_terms`; removed once the chd-atlas/feeds
packages replaced it), so a package of the whole source ingests to the
committed mirrors. Every committed ``orphanet_phenotypes`` row
is under the pinned branch (`test_every_committed_orphanet_phenotype_is_in_
the_pinned_cardiovascular_branch`); every committed ``mondo_xrefs`` row names
a disease of a committed validity mirror (measured 2026-11-12 with a csv read
of the three files: 10,011 of 10,011 rows, 7,418 distinct terms; no test holds
it).

`FilterUnavailable` is raised, and `chd-atlas ingest` refuses with ING010,
when a filter cannot be built: ``sources.yaml`` cannot be read or pins no HPO
release, the pinned release cannot be read or does not define HP:0001626, a
present validity mirror fails its schema, or neither validity mirror names a
disease. An absent validity mirror counts as one with no rows, as it does for
the registry (`chd_atlas.registry`). A release without HP:0001626 and no
disease in either mirror are refused for the reason `apply._registry` refuses
an empty registry: a filter keeping nothing would drop every row, and a count
of rows "outside the filter" would read as a fact about the package.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

from chd_atlas.ingest.package import BLANK
from chd_atlas.issues import Severity
from chd_atlas.orphanet_scope import CARDIOVASCULAR_ROOT, pinned_hpo_release
from chd_atlas.registry import VALIDITY_MIRRORS, RegistryUnavailable, read_validity_mirrors
from chd_atlas.tables import MONDO_XREFS, ORPHANET_PHENOTYPES
from chd_atlas.validate.sources import load_sources

HPO_PREFIX: Final = "HP"


class FilterUnavailable(Exception):
    """A table filter cannot be built; the message says which input failed and how."""


@dataclass(frozen=True)
class Filtered:
    """One table's rows after a `TableFilter`, and what it did not keep."""

    frame: pl.DataFrame
    # Rows whose value is a known term outside ``keep``.
    outside: int
    # Rows whose value is not in ``known``, and the distinct values, sorted.
    unknown_rows: int
    unknown_terms: tuple[str, ...]


@dataclass(frozen=True)
class TableFilter:
    """Keep a row whose ``column`` is in ``keep``; name one whose value ``known`` lacks."""

    column: str
    keep: frozenset[str]
    # None when every value is a term the filter can judge (Mondo's).
    known: frozenset[str] | None
    # Why a row is outside the filter, as the ingest summary prints it.
    reason: str

    def apply(self, frame: pl.DataFrame) -> Filtered:
        """Split ``frame``; rows keep their order. A null cell is the term ``<blank>``."""
        cells = [BLANK if cell is None else str(cell) for cell in frame[self.column].to_list()]
        known = self.known
        unknown = sorted({cell for cell in cells if known is not None and cell not in known})
        named = set(unknown)
        kept_flags = [cell in self.keep and cell not in named for cell in cells]
        # Each count is taken from the cells, not as what the kept frame left
        # over, so `ingest_package`'s row identity checks the split.
        return Filtered(
            frame=frame.filter(pl.Series("keep", kept_flags, dtype=pl.Boolean)),
            outside=sum(cell not in self.keep and cell not in named for cell in cells),
            unknown_rows=sum(cell in named for cell in cells),
            unknown_terms=tuple(unknown),
        )


def table_filter(root: Path, table: str) -> TableFilter | None:
    """The filter ingest applies to ``table`` (see the module docstring); None for most tables."""
    if table == ORPHANET_PHENOTYPES.name:
        return _orphanet_phenotypes(root)
    if table == MONDO_XREFS.name:
        return _mondo_xrefs(root)
    return None


def _orphanet_phenotypes(root: Path) -> TableFilter:
    registry, issues = load_sources(root)
    errors = [issue.message for issue in issues if issue.severity is Severity.ERROR]
    if errors:
        raise FilterUnavailable(
            "mirrors/sources.yaml cannot be read, so the pinned HPO release is unknown: "
            + "; ".join(errors)
        )
    pinned = registry.ontology_files()
    relative = pinned.get(HPO_PREFIX)
    if relative is None:
        raise FilterUnavailable("mirrors/sources.yaml pins no HPO release (ontology_prefix: HP)")
    release = pinned_hpo_release(root, pinned)
    if release is None:
        raise FilterUnavailable(f"the pinned HPO release {relative} cannot be read")
    if CARDIOVASCULAR_ROOT not in release.known:
        raise FilterUnavailable(
            f"the pinned HPO release {relative} has no term {CARDIOVASCULAR_ROOT}, so the "
            "filter would drop every row"
        )
    return TableFilter(
        column="hpo",
        keep=release.cardiovascular,
        known=release.known,
        reason=f"hpo not under {CARDIOVASCULAR_ROOT} in {relative}",
    )


def _mondo_xrefs(root: Path) -> TableFilter:
    clingen_path, gencc_path = (relative for relative, _ in VALIDITY_MIRRORS)
    try:
        clingen, gencc = read_validity_mirrors(root)
    except RegistryUnavailable as exc:
        raise FilterUnavailable(str(exc)) from exc
    keep = frozenset(
        str(disease)
        for frame in (clingen, gencc)
        for disease in frame["disease"].to_list()
        if disease is not None
    )
    if not keep:
        raise FilterUnavailable(
            f"neither {clingen_path} nor {gencc_path} names a disease, so the filter would "
            "drop every row"
        )
    return TableFilter(
        column="mondo",
        keep=keep,
        known=None,
        reason=f"mondo named by neither {clingen_path} nor {gencc_path}",
    )
