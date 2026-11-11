# src/chd_atlas/literature_query.py
"""What the atlas answers about a gene, for the literature tooling.

`literature_target` resolves a typed symbol, says whether the gene is inside the
publication gate, and builds the recorded PubMed query. `chd-atlas literature
target` prints it as one JSON line so a separate program can ask the atlas
without importing it. The query code lives here, moved verbatim from `pubmed.py`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from chd_atlas.build.validity import gene_validity
from chd_atlas.corpus import load_curation
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.runner import _gate_published_genes, _symbol_registry

TARGET_FORMAT: Final = 1


class TargetError(Exception):
    """The atlas cannot answer about this gene; the message names the cause."""


# The spec's term list (docs/superpowers/specs/2026-11-02-literature-timeline-
# design.md, section 4, item 1) with every "(s)" expanded into the singular and
# then the plural as two terms -- "(s)" is not PubMed syntax (ruling R6). The
# spec lets the owner edit this list before the pilot; a test pins the tuple so
# a change to it is a visible diff.
CHD_TERMS: Final[tuple[str, ...]] = (
    "congenital heart",
    "heart defect",
    "heart defects",
    "cardiac defect",
    "cardiac defects",
    "heart malformation",
    "heart malformations",
    "cardiac malformation",
    "cardiac malformations",
    "septal defect",
    "septal defects",
    "tetralogy of Fallot",
    "transposition of the great arteries",
    "hypoplastic left heart",
    "heterotaxy",
    "conotruncal",
    "outflow tract",
    "atrioventricular septal defect",
    "coarctation",
)


def _phrase(text: str) -> str:
    """One quoted `[tiab]` phrase. PubMed has no escape for `"`, so it is dropped."""
    return '"' + " ".join(text.replace('"', " ").split()) + '"[tiab]'


def build_query(symbol: str, previous: Sequence[str], diseases: Sequence[str]) -> str:
    """The recorded search for one gene.

    The symbol comes first, then the previous symbols and the disease labels in
    sorted order; the CHD terms keep `CHD_TERMS` order. Nothing is filtered: the
    caller decides which previous symbols and diseases belong.
    """
    gene_block = " OR ".join([_phrase(symbol), *(_phrase(p) for p in sorted(previous))])
    chd_block = " OR ".join(
        ['"Heart Defects, Congenital"[mh]']
        + [_phrase(t) for t in CHD_TERMS]
        + [_phrase(d) for d in sorted(diseases)]
    )
    return f"({gene_block}) AND ({chd_block}) NOT preprint[pt] NOT review[pt]"


@dataclass(frozen=True)
class LiteratureTarget:
    """The atlas's answer about one gene."""

    gene: str
    symbol: str
    published: bool
    query: str

    def as_json(self) -> dict[str, object]:
        return {
            "format": TARGET_FORMAT,
            "gene": self.gene,
            "symbol": self.symbol,
            "published": self.published,
            "query": self.query,
        }


def resolve_symbol(root: Path, symbol: str) -> tuple[str, str, tuple[str, ...]]:
    """`(HGNC id, approved symbol, previous symbols)` for a typed symbol.

    Any unambiguous resolution is accepted -- approved, alias or previous
    (ruling R19) -- because curators type old names; the caller writes the
    approved one.
    """
    registry = _symbol_registry(root)
    if registry is None:
        raise TargetError("mirrors/genes.tsv is missing or unreadable")
    resolved = registry.resolve(symbol)
    if resolved.gene is None:
        detail = (
            f"{symbol!r} is claimed by {', '.join(resolved.candidates)}"
            if resolved.candidates
            else f"{symbol!r} is not a known gene symbol"
        )
        raise TargetError(detail)
    approved = registry.approved_symbol(resolved.gene)
    if approved is None:
        raise TargetError(f"{resolved.gene} has no approved symbol in mirrors/genes.tsv")
    return resolved.gene, approved, tuple(registry.previous_symbols(resolved.gene))


def literature_target(root: Path, symbol: str) -> LiteratureTarget:
    """Resolve `symbol`, test the publication gate and build the gene's query.

    The validity mirrors are read before the gate because the gate answers an
    empty set when either is unreadable; reading them first keeps "unreadable"
    from reporting as "empty". A gene with no in-scope validity record gets a
    query with no disease names.
    """
    gene, approved, previous = resolve_symbol(root, symbol)
    corpus, _ = load_curation(root)
    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    if clingen is None or gencc is None:
        raise TargetError("the validity mirrors are unreadable")
    published_genes = _gate_published_genes(root, corpus)
    if not published_genes:
        raise TargetError("the publication gate is empty")
    scope = {str(entry.id) for entry in corpus.chd_scope}
    validity = gene_validity(clingen, gencc, in_scope=scope).get(gene)
    diseases = sorted({r.disease_label for r in validity.records}) if validity else []
    return LiteratureTarget(
        gene=gene,
        symbol=approved,
        published=gene in published_genes,
        query=build_query(approved, previous, diseases),
    )
