"""Rewrite a package's gene ids to HGNC ids against the pinned HGNC package.

A package carries gene ids as its source gives them (`GeneIdKind`). The atlas
keys every join on HGNC id, so `resolve_frame` rewrites the column and returns,
beside the rows that resolved, one `Resolved` per distinct id that did not.
**An unresolved id is named, never dropped silently**: that list is what the
ingest report and the refusal threshold read.

The registry is the pinned HGNC package's `genes` table read from the local
cache (`local/packages/hgnc/<version>/`), never `mirrors/genes.tsv`: the two
differ by construction (the mirror is a filtered registry of the atlas's own
genes) and falling back would resolve against the wrong population. A missing
pin or cache is `ING008`. Symbols go through `GeneRegistry.resolve`, so the
precedence (approved, alias, previous; ambiguous refused) is the registry's own.

Two deliberate limits:

- An `ENSEMBL_GENE_ID` is matched exactly. A version suffix (`ENSG....4`) is
  **not** stripped: an id with a suffix the HGNC table does not carry is
  `UNKNOWN`, which is reported, where guessing at a stripped match is not.
- An `HGNC_ID` that has the right shape resolves to itself without being looked
  up. Whether the registry contains it is the next step's filter, which counts
  the rest.
"""

from __future__ import annotations

import hashlib
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import polars as pl

from chd_atlas.genes import GeneRegistry, Resolution, Resolved
from chd_atlas.identifiers import HGNC_PATTERN
from chd_atlas.ingest.package import IngestRefused, _decompress, _parse, _read_bytes
from chd_atlas.models.package import ManifestError, load_manifest
from chd_atlas.tables import cell_text
from chd_atlas.validate.sources import load_sources
from chd_atlas.vocab import GeneIdKind

_HGNC = re.compile(HGNC_PATTERN)
_BLANK = "<blank>"
_HGNC_SOURCE = "hgnc"


def package_cache(root: Path, source_id: str, version: str) -> Path:
    """Where a curator keeps a source package's files; ``local/`` is gitignored."""
    return root / "local" / "packages" / source_id / version


@dataclass(frozen=True)
class GeneIdResolver:
    kind: GeneIdKind
    _registry: GeneRegistry | None = field(default=None, repr=False)
    _ensembl: Mapping[str, tuple[str, ...]] = field(default_factory=dict, repr=False)

    def resolve(self, raw: object) -> Resolved:
        """One raw id, stripped, resolved or reported. A blank is ``<blank>``, UNKNOWN."""
        name = cell_text(raw)
        if name is None:
            return Resolved(_BLANK, Resolution.UNKNOWN)
        if self.kind is GeneIdKind.HGNC_ID:
            if _HGNC.match(name):
                return Resolved(name, Resolution.APPROVED, name, (name,))
            return Resolved(name, Resolution.UNKNOWN)
        if self.kind is GeneIdKind.ENSEMBL_GENE_ID:
            claimed = self._ensembl.get(name, ())
            if not claimed:
                return Resolved(name, Resolution.UNKNOWN)
            if len(claimed) > 1:
                return Resolved(name, Resolution.AMBIGUOUS, None, claimed)
            return Resolved(name, Resolution.APPROVED, claimed[0], claimed)
        if self._registry is None:
            raise IngestRefused("ING008", "no HGNC registry was loaded for a symbol resolver")
        return self._registry.resolve(name)

    @classmethod
    def from_root(cls, root: Path, kind: GeneIdKind) -> GeneIdResolver:
        if kind is GeneIdKind.HGNC_ID:
            return cls(kind)
        rows = _hgnc_rows(root)
        if kind is GeneIdKind.SYMBOL:
            return cls(kind, GeneRegistry(rows))
        claims: dict[str, set[str]] = {}
        for row in rows:
            ensembl = cell_text(row.get("ensembl_gene"))
            gene = cell_text(row.get("hgnc_id"))
            if ensembl is not None and gene is not None:
                claims.setdefault(ensembl, set()).add(gene)
        return cls(kind, None, {key: tuple(sorted(value)) for key, value in claims.items()})


def _hgnc_rows(root: Path) -> list[dict[str, object]]:
    registry, issues = load_sources(root)
    if issues:
        raise IngestRefused(
            "ING008",
            "mirrors/sources.yaml could not be loaded, so the HGNC pin cannot be "
            f"read: {'; '.join(issue.message for issue in issues)}",
        )
    source = next((item for item in registry.sources if item.id == _HGNC_SOURCE), None)
    if source is None or source.package is None:
        raise IngestRefused(
            "ING008",
            "source 'hgnc' has no package pin; the HGNC package is ingested first, in phase 2",
        )
    directory = package_cache(root, _HGNC_SOURCE, source.package.version)
    if not directory.is_dir():
        raise IngestRefused(
            "ING008",
            f"the pinned HGNC package {source.package.version} is not cached at {directory}; "
            "the HGNC package is ingested first, in phase 2",
        )
    try:
        manifest_bytes = (directory / "manifest.json").read_bytes()
        digest = hashlib.sha256(manifest_bytes).hexdigest()
        if digest != source.package.manifest_sha256:
            raise IngestRefused(
                "ING008",
                f"the cached HGNC manifest's sha256 is {digest}, but the pin says "
                f"{source.package.manifest_sha256}; the cache is not the pinned package",
            )
        manifest = load_manifest(directory / "manifest.json")
        if manifest.package != _HGNC_SOURCE or manifest.version != source.package.version:
            raise IngestRefused(
                "ING008",
                f"the cached manifest is package '{manifest.package}' version "
                f"'{manifest.version}', not the pinned hgnc {source.package.version}",
            )
        spec = next((table for table in manifest.tables if table.name == "genes"), None)
        if spec is None:
            raise IngestRefused("ING008", f"the cached HGNC package {directory} has no genes table")
        data = _read_bytes(directory, spec)
        frame = _parse(_decompress(data, spec) if spec.file.endswith(".gz") else data, spec)
    except OSError as exc:
        raise IngestRefused("ING008", f"the cached HGNC manifest cannot be read: {exc}") from exc
    except ManifestError as exc:
        raise IngestRefused("ING008", f"the cached HGNC package is unreadable: {exc}") from exc
    except IngestRefused as exc:
        if exc.code == "ING008":
            raise
        raise IngestRefused("ING008", f"the cached HGNC genes table is unreadable: {exc}") from exc
    if "hgnc_id" not in frame.columns or "symbol" not in frame.columns:
        raise IngestRefused("ING008", "the cached HGNC genes table lacks hgnc_id or symbol")
    return frame.to_dicts()


def resolve_frame(
    frame: pl.DataFrame, column: str, resolver: GeneIdResolver
) -> tuple[pl.DataFrame, tuple[Resolved, ...]]:
    """Keep the rows whose id resolved, rewriting ``column`` to HGNC ids.

    Rows keep their order. The second value has one `Resolved` per distinct
    unresolved id (after stripping), sorted by symbol; ambiguous ids carry
    their sorted candidates.
    """
    raw = frame[column].to_list()
    genes: list[str | None] = []
    unresolved: dict[str, Resolved] = {}
    cache: dict[str, Resolved] = {}
    for value in raw:
        key = cell_text(value) or ""
        result = cache.get(key)
        if result is None:
            result = cache[key] = resolver.resolve(value)
        genes.append(result.gene)
        if result.gene is None:
            unresolved[result.symbol] = result
    keep = pl.Series("keep", [gene is not None for gene in genes])
    rewritten = frame.with_columns(pl.Series(column, genes, dtype=pl.String)).filter(keep)
    return rewritten, tuple(sorted(unresolved.values(), key=lambda item: item.symbol))
