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

**Ambiguity is refused, never guessed.** `scripts/convert_burden.py` declines to
match aliases at all, and its stated reason is exactly right: "an alias that is
another gene's current symbol would silently attribute one gene's burden to
another". Measured against the committed registry, that collision does not occur
today -- 185 approved symbols, 382 alias symbols, **zero** overlap between the
two sets and **zero** aliases claimed by more than one gene (re-measured
2026-10-25, after the registry widened from 154 genes) -- but the guard is
structural rather than a bet on that staying true, because the registry is
regenerated from upstream and the number that matters is the one after the next
refresh.

**Three tiers, not two: approved, then alias, then previous.** `prev_symbols`
arrived with the mirror column of the same name; before it this module could
not resolve a retired name at all, which is what left `BVES` unresolvable while
ClinGen used it (the standing `GEN003`). Measured over the regenerated registry
on 2026-08-24 -- 154 approved, 309 alias, **69 distinct** previous symbols
across 56 genes:

- `approved & alias` and `approved & previous` are both **empty**.
- `alias & previous` holds exactly one string, **`ODD`**: a live alias of GJA1
  (HGNC:4274) and the retired name of OSR1 (HGNC:8111). It is the whole reason
  the previous symbols are a *third* tier rather than more entries in the alias
  map. Merged, `ODD` is claimed twice and resolves to nothing; ranked, the live
  alias wins, which is what a source writing `ODD` today means. HGNC retires a
  name into `prev_symbol` and leaves the live synonyms in `alias_symbol`, so
  "currently used" beating "formerly official" is the upstream's own semantics
  rather than a preference invented here.
- **`MADH7`** is claimed by two genes *within* the previous tier -- SMAD6
  (HGNC:6772) and SMAD7 (HGNC:6773) both retired it -- and is refused, exactly
  as an ambiguous alias would be. Precedence resolves *between* tiers; it never
  breaks a tie inside one.

**A tier is consulted only when every tier above it is silent, and an ambiguous
tier stops the search.** A symbol two genes claim as an alias must not be
quietly answered by some third gene's retired name: that would turn a reportable
collision into a confident wrong id, which is the failure this module exists to
prevent.
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
    """The symbol is a current alias of exactly one gene.

    **Aliases only.** This docstring read "an alias or previous symbol" while
    the registry carried no previous symbols at all, so it described a tier
    that did not exist and would have been the natural place to merge them
    when one did -- see `PREVIOUS` for why they are ranked instead.
    """

    PREVIOUS = "previous"
    """The symbol is the retired HGNC symbol of exactly one gene, and no gene
    claims it as an approved symbol or a current alias.

    Ranked below `ALIAS` rather than merged into it, on the module docstring's
    measurement: `ODD` is GJA1's live alias and OSR1's retired name, and only
    an ordering answers it.
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

    **Approved symbols win over aliases, which win over previous symbols, and
    the precedence is not a tie-break to be tidied away.** HGNC retires a symbol
    by making it a *previous* symbol of the gene that kept the name, so a string
    that is one gene's approved symbol and another's alias names the first gene
    today. Resolving to the alias would attribute a row to a gene the source did
    not mean. The same argument one rung down puts a live alias ahead of a
    retired name -- see the module docstring, where the one string this actually
    decides (`ODD`) is named and measured.
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
        genes retired (`MADH7`) and one a live alias outranks (`ODD`) both
        resolve to something other than `Resolution.PREVIOUS` for this gene, and
        are left out.
        """
        names = []
        for name in self._previous:
            resolved = self.resolve(name)
            if resolved.resolution is Resolution.PREVIOUS and resolved.gene == gene:
                names.append(name)
        return tuple(sorted(names))

    def resolve(self, symbol: str) -> Resolved:
        """One symbol, resolved or reported. Never guesses, never raises.

        Walks the three tiers in order and **stops at the first that claims the
        name at all** -- including when that tier claims it twice. Falling
        through an ambiguous tier to a quieter one below would answer `MADH7`
        (SMAD6 and SMAD7 both retired it) with whichever gene happened to hold
        it further down, turning a reportable collision into a confident wrong
        id. The empty-tier check and the ambiguity check are therefore separate
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
        for index, resolution in (
            (self._aliases, Resolution.ALIAS),
            (self._previous, Resolution.PREVIOUS),
        ):
            claimed = index.get(name)
            if not claimed:
                continue
            if len(claimed) > 1:
                return Resolved(
                    symbol=name,
                    resolution=Resolution.AMBIGUOUS,
                    candidates=tuple(sorted(claimed)),
                )
            only = next(iter(claimed))
            return Resolved(symbol=name, resolution=resolution, gene=only, candidates=(only,))
        return Resolved(symbol=name, resolution=Resolution.UNKNOWN)


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
