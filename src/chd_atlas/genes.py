# src/chd_atlas/genes.py
"""The one place a gene symbol becomes an HGNC id. Issue #33.

**Every cross-source join in this atlas keys on HGNC id, never on symbol**, and
this module exists so that a source which publishes symbols alone can be brought
onto that key in exactly one place, under exactly one set of rules.

The owner's instruction, 2026-08-06:

> there will be a lot of gene mapping across sources. you need to ensure this
> occur at the gene id level and not at the gene symbol name, since the later is
> more prone to mismatch due to aliases.

The hazard is measured, not hypothetical. Across the two committed validity
mirrors (2026-08-14) **six HGNC ids carry two different symbols**, because
ClinGen and GenCC disagree about the approved name for one gene: HGNC:1152
(`BVES`/`POPDC1`), HGNC:17272 (`CENPJ`/`CPAP`), HGNC:17861
(`IFT54`/`TRAF3IP1`), HGNC:26742 (`ASPNAT`/`NAT8L`), HGNC:28178
(`CCDC115`/`VMA22`) and HGNC:7687 (`COXFA4`/`NDUFA4`). A symbol-keyed map built
first-seen-wins picks whichever the read order surfaced, which is silent and not
even stable across a change in read order.

**An unresolved symbol is reported, never dropped.** That is the load-bearing
rule and the reason `resolve` returns a `Resolution` rather than `str | None`: a
caller that cannot tell "no such symbol" from "this symbol names two genes"
cannot report either, and a gene that fails to map and simply vanishes is
precisely the silent evidence loss this project treats as its characteristic
failure -- the build stays green and the gene is not there.

**Ambiguity is refused, never guessed.** The retired Audain 2026 converter declined to
match aliases at all, and its stated reason was exactly right: "an alias that is
another gene's current symbol would silently attribute one gene's burden to
another". Its feeds recipe (`audain2026_sd3`) now ships symbols verbatim and
`chd-atlas ingest` resolves them here (approved, alias, previous), which is how
CRIPTO, POPDC1 and CIROZ were recovered. Measured against the committed
registry, that collision does not occur
today -- 186 approved symbols, 384 alias symbols, **zero** overlap between the
two sets and **zero** aliases claimed by more than one gene (re-measured
2026-11-14; 185 and 382 on 2026-10-25, after the registry widened from 154
genes) -- but the guard is
structural rather than a bet on that staying true, because the registry is
regenerated from upstream and the number that matters is the one after the next
refresh.

**An approved symbol first; then an alias or previous symbol, answered
together.** `prev_symbols` arrived with the mirror column of the same name;
before it this module could not resolve a retired name at all, which is what
left `BVES` unresolvable while ClinGen used it (the standing `GEN003`).

**A name the alias and previous tiers give to more than one gene between them
is refused** (owner decision F27, 2026-11-13). Until then a live alias outranked
a retired name, on the argument that HGNC leaves current synonyms in
`alias_symbol`. Measured on the HGNC release the pins use (2026-10-06, 45,233
genes): 356 names sit in both tiers and 352 of them name different genes, and
no order is right for all of them. The old name means the previous gene for `MLL2`
(KMT2D's withdrawn symbol, KMT2B's alias) and `HRAS1` (HRAS's, RASSF7's), so
alias-first lands a source's rows on another gene, reported only as "outside
the registry"; it means the alias gene, or nobody, for `CBP` (CREBBP's and
PAG1's alias, OPN1LW's withdrawn symbol) and `SIP1` (ZEB2, SCAF11; GEMIN2), so
previous-first answers a colour-vision gene for CREBBP. Refused, each is named
in the ingest report with every claimant. In the pinned packages that is SD3's
TAZ, LOR, SARS, MPP6 and DUSP27 and S6's DEC1, MPP6 and SEPT2 (42 rows), each
one gene's withdrawn symbol and another's alias, which alias-first had answered
with the alias gene; none is a registry gene, and no mirror or published byte
changed (three real `ingest --all` and builds, one per order).

The committed registry (186 genes, re-measured 2026-11-14) has one such name,
**`ODD`**: GJA1's (HGNC:4274) alias and OSR1's (HGNC:8111) withdrawn symbol. It
resolved to GJA1 until F27 and is refused now. **`MADH7`** is claimed twice
*within* the previous tier -- SMAD6 (HGNC:6772) and SMAD7 (HGNC:6773) both
retired it -- and was always refused. A name one gene holds as both an alias
and a withdrawn symbol (`HIGM1`, CD40LG's) resolves, labelled an alias.

**The tiers stay separate columns because the distinction is still read.**
`previous_symbols` lists only names that resolve as `PREVIOUS`, and the search
index publishes a gene's aliases ahead of its withdrawn names.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum


class Resolution(StrEnum):
    """Why `GeneRegistry.resolve` returned what it did."""

    APPROVED = "approved"
    """The symbol is some gene's current HGNC-approved symbol."""

    ALIAS = "alias"
    """The symbol is a current alias of exactly one gene, and no other gene
    claims it as an alias or a previous symbol.

    Also the label when that one gene holds the name as both its alias and its
    withdrawn symbol (`HIGM1`, CD40LG's).
    """

    PREVIOUS = "previous"
    """The symbol is the retired HGNC symbol of exactly one gene, and no gene
    claims it as an approved symbol or a current alias.

    A name that is one gene's retired symbol and another's alias is
    `AMBIGUOUS`, not this (F27; the module docstring has the measurement).
    """

    UNKNOWN = "unknown"
    """No gene in the registry claims this symbol. Report it; do not drop the row."""

    AMBIGUOUS = "ambiguous"
    """More than one gene claims it. Report it; never pick one."""


@dataclass(frozen=True)
class Resolved:
    """One symbol lookup: the outcome, and the id if there is exactly one.

    `gene` is `None` for `UNKNOWN` and for `AMBIGUOUS`, and `candidates` carries
    every id that claimed the symbol so a report can name them. The two are kept
    apart deliberately -- a caller that unpacked a bare `str | None` would treat
    "nobody claims this" and "two genes claim this" identically, and those need
    opposite responses from a curator.
    """

    symbol: str
    resolution: Resolution
    gene: str | None = None
    candidates: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.gene is not None


class GeneRegistry:
    """`mirrors/genes.tsv` as a symbol -> HGNC id lookup.

    Built from the rows of the registry rather than reading the file, so this
    stays a pure structure over already-loaded data -- the same shape
    `validate_scope_terms` and `validate_source_references` have, and what lets
    it be tested without a filesystem.

    **An approved symbol wins over every alias and previous symbol, and that
    precedence is not a tie-break to be tidied away.** HGNC retires a symbol by
    making it a *previous* symbol of the gene that kept the name, so a string
    that is one gene's approved symbol and another's alias names the first gene
    today. Resolving to the alias would attribute a row to a gene the source did
    not mean. Below the approved tier there is no precedence: an alias and a
    previous symbol are answered together, and a name they give to two genes is
    refused (F27; the module docstring has the measurement).
    """

    def __init__(self, rows: Iterable[Mapping[str, object]]) -> None:
        approved: dict[str, str] = {}
        aliases: dict[str, set[str]] = {}
        previous: dict[str, set[str]] = {}
        for row in rows:
            gene = _text(row.get("hgnc_id"))
            if gene is None:
                # No id is nothing to key on. `_gene_registry` in
                # `build/runner.py` drops these too, and TBL003 reports the
                # blank, so this is a bypassed-gate path rather than a silent
                # policy of its own.
                continue
            symbol = _text(row.get("symbol"))
            if symbol is not None:
                approved[symbol] = gene
            for column, index in (("aliases", aliases), ("prev_symbols", previous)):
                raw = _text(row.get(column))
                for name in (part.strip() for part in (raw or "").split("|")):
                    if name:
                        index.setdefault(name, set()).add(gene)
        self._approved = approved
        self._aliases = aliases
        self._previous = previous

    def __len__(self) -> int:
        return len(self._approved)

    def approved_symbols(self) -> tuple[str, ...]:
        """Every current HGNC-approved symbol the registry carries, sorted.

        Public so that a caller needing the registry's population does not reach
        into `_approved`. Sorted because a caller iterating it may put the result
        in a payload, and `encode_json`'s `sort_keys` orders dict keys only.
        """
        return tuple(sorted(self._approved))

    def approved_symbol(self, gene: str) -> str | None:
        """The current approved symbol of an HGNC id, or None when the id is unknown."""
        for symbol, hgnc_id in self._approved.items():
            if hgnc_id == gene:
                return symbol
        return None

    def previous_symbols(self, gene: str) -> tuple[str, ...]:
        """The retired names that `resolve` answers with `gene`, sorted.

        Defined through `resolve` rather than the raw index so the answer is
        exactly the set of names a search may attach to this gene: a name two
        genes retired (`MADH7`), one another gene holds as an alias (`ODD`) and
        one this gene also holds as an alias (`HIGM1`) all resolve to something
        other than `Resolution.PREVIOUS` for this gene, and are left out.
        """
        names = []
        for name in self._previous:
            resolved = self.resolve(name)
            if resolved.resolution is Resolution.PREVIOUS and resolved.gene == gene:
                names.append(name)
        return tuple(sorted(names))

    def resolve(self, symbol: str) -> Resolved:
        """One symbol, resolved or reported. Never guesses, never raises.

        An approved symbol answers outright. Otherwise every gene that claims
        the name as an alias **or** a previous symbol is a candidate, and one
        candidate is an answer while two are a refusal naming both. Taking the
        alias tier alone when it speaks would answer `ODD` with GJA1 and drop
        OSR1 from the report; taking the previous tier first would answer `CBP`
        with OPN1LW (F27). The empty check and the ambiguity check are separate
        branches, not one `len(claimed) == 1`.
        """
        name = symbol.strip()
        if not name:
            return Resolved(symbol=symbol, resolution=Resolution.UNKNOWN)
        gene = self._approved.get(name)
        if gene is not None:
            return Resolved(
                symbol=name, resolution=Resolution.APPROVED, gene=gene, candidates=(gene,)
            )
        aliased = self._aliases.get(name, set())
        claimed = aliased | self._previous.get(name, set())
        if not claimed:
            return Resolved(symbol=name, resolution=Resolution.UNKNOWN)
        if len(claimed) > 1:
            return Resolved(
                symbol=name,
                resolution=Resolution.AMBIGUOUS,
                candidates=tuple(sorted(claimed)),
            )
        only = next(iter(claimed))
        resolution = Resolution.ALIAS if aliased else Resolution.PREVIOUS
        return Resolved(symbol=name, resolution=resolution, gene=only, candidates=(only,))


def _text(value: object) -> str | None:
    """One registry cell as text, with blank and whitespace-only treated as absent.

    The same rule `tables.cell_text` applies, and for the same reason: a
    single space survives a `min_length` check and renders as nothing, so it
    would enter the index as a symbol matching a query nobody typed.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text or None
