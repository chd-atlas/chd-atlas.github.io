# src/chd_atlas/validate/genes.py
"""Symbol-to-id resolution, checked. Issue #33: report, never drop.

`GEN000` follows `REF000`/`SRC000`/`ONT000`/`SCP000`: when the registry cannot
be read, say "checks skipped" once instead of reporting every symbol in the
corpus as unresolvable. A missing reference table must never present as hundreds
of dangling references -- that buries the real error under cascading false
positives, and distinguishing "the data is wrong" from "the reference data never
loaded" is one of this project's two load-bearing guarantees.

`GEN004` (owner decision 2026-10-26, R13) is the registry's other job: the
Ensembl and NCBI Gene ids every published gene page links by, held to the rule
the build refuses a page on (`build/identity.py::identity_id_problems`).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

import polars as pl

from chd_atlas.build.identity import EMPTY_IDENTITY, GeneIdentity, identity_id_problems
from chd_atlas.genes import GeneRegistry, Resolution
from chd_atlas.issues import Severity, ValidationIssue


def validate_symbol_resolution(
    symbols: Iterable[tuple[str, str]],
    registry: GeneRegistry | None,
    location: str,
) -> list[ValidationIssue]:
    """Every symbol a symbol-keyed source names must resolve to exactly one gene.

    `symbols` is `(symbol, where)` pairs -- the caller says where each came from,
    because "SYMBOL column, row 412" is what a curator needs and this function
    has no way to know it.

    **GEN001 (unmapped) and GEN002 (ambiguous) are the whole point of the
    issue.** A source that publishes symbols alone -- CHDgene is the likely
    first, see #31 -- gets joined onto HGNC ids somewhere, and the failure mode
    is not a crash: it is a row quietly not appearing. `scripts/convert_burden.py`
    already demonstrates the discipline, naming all 25 Excel date-corrupted
    symbols (`2001-03-01` = MARCH1/MARC1 collapsed) on every run, which is the
    only reason anyone knows about them.

    Warnings rather than errors, deliberately. An unresolvable symbol in an
    upstream file is not a defect in *this* repository's curation, and refusing
    to build over one would make a third party's naming choice able to take the
    site down. The gate is 0 errors; what this must do is be impossible to miss.
    """
    wanted = sorted(set(symbols))
    if registry is None:
        # **No skip warning when there was nothing to skip.** A skip means "a
        # check that would have run did not", and with no symbols to resolve no
        # check would have run. Reporting one anyway would put a lone warning in
        # a report that has no error to explain it -- and `ValidationReport.ok`
        # ignores warnings, so a warnings-only report exits 0. Pinned by
        # `test_every_skip_warning_arrives_with_the_error_that_caused_it`, which
        # caught exactly this.
        return _skipped(location, "symbol resolution") if wanted else []

    issues: list[ValidationIssue] = []
    for symbol, where in wanted:
        found = registry.resolve(symbol)
        if found.resolution is Resolution.UNKNOWN:
            issues.append(
                ValidationIssue(
                    code="GEN001",
                    severity=Severity.WARNING,
                    location=where,
                    message=(
                        f"gene symbol {symbol!r} resolves to no gene in mirrors/genes.tsv; "
                        "the row cannot be keyed on an HGNC id and would be dropped"
                    ),
                )
            )
        elif found.resolution is Resolution.AMBIGUOUS:
            issues.append(
                ValidationIssue(
                    code="GEN002",
                    severity=Severity.WARNING,
                    location=where,
                    message=(
                        f"gene symbol {symbol!r} is claimed by "
                        f"{', '.join(found.candidates)}; refusing to guess which"
                    ),
                )
            )
    return issues


def validate_mirror_symbols(
    labelled: Iterable[tuple[str, str, str]],
    registry: GeneRegistry | None,
    location: str,
) -> list[ValidationIssue]:
    """Every symbol our own mirrors print must resolve back to the id beside it.

    `labelled` is `(hgnc_id, symbol, where)` -- what a validity mirror already
    carries in both columns, so this is a *self-consistency* check on data that
    is already keyed correctly rather than a join.

    **GEN003 exists because the registry goes stale and nothing else notices,
    and its first finding was not staleness at all.** From 2026-08-06 to
    2026-08-24 it stood on HGNC:1152: `mirrors/genes.tsv` recorded `POPDC1`
    with aliases `POP1|HBVES` while ClinGen labelled the same gene `BVES`, and
    `BVES` was not among those aliases (`HBVES` is a different string). The
    mirror was **current**. `BVES` is HGNC:1152's *previous* symbol, the
    converter read `alias_symbol` and nothing else, and the mirror had no
    column to hold a retired name -- so the recorded remedy, "regenerate the
    registry", reproduced the same row byte-for-byte and the warning survived
    every refresh. It cleared when `prev_symbols` and this module's third
    resolution tier landed, not when anyone re-ran the converter.

    Kept, because the failure mode it was written for is real and now
    reachable: a mirror that genuinely lags upstream. A symbol-keyed source
    naming a symbol the registry cannot resolve would drop the gene silently,
    which is the exact failure #33 exists to prevent. Design decision D6 keeps
    regeneration a human action rather than something a build does, so this
    reports and stops.

    A symbol resolving to a *different* gene is the serious case and is reported
    the same way rather than more loudly, because it is not reachable from a
    correct mirror: both columns come from the same upstream row. It is here so
    that a hand-edited mirror cannot introduce it undetected.
    """
    rows = sorted(set(labelled))
    if registry is None:
        # Same rule as above: a skip is only a skip when work existed.
        return _skipped(location, "mirror symbol") if rows else []

    known = {gene for gene, _, _ in rows} & _registry_ids(registry)
    issues: list[ValidationIssue] = []
    for gene, symbol, where in rows:
        # Only genes the registry actually carries. A mirror row for a gene
        # outside the 185-gene registry is not a naming disagreement -- the
        # registry is deliberately narrower than the mirrors, and reporting
        # those would drown the one real finding in 3,000 non-findings.
        if gene not in known:
            continue
        found = registry.resolve(symbol)
        if found.gene == gene:
            continue
        detail = (
            f"resolves to {found.gene}"
            if found.gene is not None
            else f"resolves to nothing ({found.resolution.value})"
        )
        issues.append(
            ValidationIssue(
                code="GEN003",
                severity=Severity.WARNING,
                location=where,
                message=(
                    f"{gene} is labelled {symbol!r} here, but that symbol {detail} "
                    "in mirrors/genes.tsv; regenerate the registry with "
                    "scripts/convert_hgnc.py so a symbol-keyed source can resolve it"
                ),
            )
        )
    return issues


def validate_identity_ids(
    identities: Mapping[str, GeneIdentity] | None,
    published: Iterable[str],
    location: str,
) -> list[ValidationIssue]:
    """GEN004: every published gene's identity ids, by the rule the build refuses a page on.

    Owner decision 2026-10-26 (R13). `build/identity.py::external_links`
    refuses a published gene whose Ensembl or NCBI Gene id is missing or
    malformed, and nothing in `chd-atlas validate` reported either: measured
    2026-10-26 on a copy of the committed corpus, TBX5's Ensembl id cut to
    `ENSG0089225` validated with 0 errors and the build refused the repository
    after writing had begun (a missing NCBI Gene id behaved the same). The
    rule is `identity.identity_id_problems`, shared, so the two cannot
    disagree; the population is the build's too, the published genes, so a
    registry row for a gene the gate does not admit is not this check's
    business. A published gene with no registry row at all lacks both ids.

    ERROR, because the build refuses on it. A validator reports and never
    raises. `identities` is `build.identity.read_gene_identities`' return, the
    reader the build uses. `None` is a registry that is absent, unreadable or
    missing an id column, which TBL008 or TBL001 reports as an error, so this
    returns nothing rather than one GEN004 per published gene.

    **An empty mapping is not `None`, and is one GEN004** (review,
    2026-10-27). A header-only `mirrors/genes.tsv` reads, as no rows, and
    neither TBL008 nor TBL001 fires on it. The caller passed `or None` to
    spare the report one error per published gene, which switched the check
    off: measured on a copy with TBX5's assertion moved out (so no REF001
    either), validate reported 0 errors and the build raised from
    `external_links` after writing had begun. So an empty registry with
    published genes is one error naming the file, and the per-gene loop does
    not run.
    """
    if identities is None:
        return []
    genes = sorted(published)
    if not identities and genes:
        return [
            ValidationIssue(
                code="GEN004",
                severity=Severity.ERROR,
                location=location,
                message=(
                    f"the gene registry holds no gene rows, and {len(genes)} published genes' "
                    "pages link gnomAD and ClinVar by the ids it carries; the build refuses "
                    "every such page"
                ),
            )
        ]
    issues: list[ValidationIssue] = []
    for gene in genes:
        missing, malformed = identity_id_problems(identities.get(gene, EMPTY_IDENTITY))
        problems = ([f"no {' or '.join(missing)}"] if missing else []) + malformed
        if problems:
            issues.append(
                ValidationIssue(
                    code="GEN004",
                    severity=Severity.ERROR,
                    location=location,
                    message=(
                        f"{gene} is published and its gene page links gnomAD and ClinVar by "
                        f"these ids: {'; '.join(problems)}; the build refuses such a page"
                    ),
                )
            )
    return issues


def _skipped(location: str, what: str) -> list[ValidationIssue]:
    """The one GEN000, in one place, so both entry points word it identically."""
    return [
        ValidationIssue(
            code="GEN000",
            severity=Severity.WARNING,
            location=location,
            message=f"mirrors/genes.tsv could not be read; {what} checks skipped",
        )
    ]


def _registry_ids(registry: GeneRegistry) -> set[str]:
    """Every HGNC id the registry carries, via its own public lookup.

    Reaches through `resolve` on the approved symbols rather than touching a
    private attribute, so this module has no back channel into `GeneRegistry`'s
    internals and the two can be changed apart.
    """
    return {
        found.gene
        for found in (registry.resolve(symbol) for symbol in registry.approved_symbols())
        if found.gene is not None
    }


def validate_one_uniprot_per_gene(genes: pl.DataFrame) -> list[ValidationIssue]:
    """GEN005: a row of `mirrors/genes.tsv` names at most one UniProt accession.

    The `genes` table is format version 2, whose `uniprot` column accepts
    HGNC's `|`-separated list (HGNC lists several for 57 genes), so that an
    every-gene HGNC package conforms. Every row of `mirrors/genes.tsv` is a
    registry gene, though, and a registry gene must name exactly one protein --
    RBM10 (HGNC:9896) lists `P0DW28|P98175` and is pinned to `P98175`.

    An ERROR, not a warning, because consumers of `genes.uniprot` assume one
    accession (measured by reading the call sites): `build/identity.py` builds a
    UniProt link from the cell, `build/omics.py` maps accession to gene, and
    `validate/referential.py` reads it too. A list reaching the build would
    publish a broken link. A null or blank cell yields nothing.
    """
    if "uniprot" not in genes.columns or "hgnc_id" not in genes.columns:
        return []
    issues: list[ValidationIssue] = []
    for gene, cell in zip(genes["hgnc_id"].to_list(), genes["uniprot"].to_list(), strict=True):
        if cell is None or "|" not in str(cell):
            continue
        accessions = [part for part in str(cell).split("|") if part]
        issues.append(
            ValidationIssue(
                code="GEN005",
                severity=Severity.ERROR,
                location="mirrors/genes.tsv",
                message=(
                    f"{gene} lists {len(accessions)} UniProt accessions "
                    f"({', '.join(accessions)}); a registry gene names exactly one "
                    "protein, so pin the gene's own protein in the HGNC recipe"
                ),
            )
        )
    return sorted(issues, key=lambda issue: issue.message)


# How many missing ids GEN006 names before it gives the total instead.
_GEN006_SHOWN = 10


def validate_registry_is_current(
    derived: frozenset[str], genes: pl.DataFrame
) -> list[ValidationIssue]:
    """GEN006: every gene of the derived registry has a row in `mirrors/genes.tsv`.

    `derived` is `registry.derive_registry` over the scope and the two validity
    mirrors: the genes `chd-atlas ingest hgnc` writes. A gene it holds that the
    file lacks has no symbol, identity ids or UniProt accession in the atlas,
    and an ingest filtered to the file drops its rows as outside the registry.
    A scope or validity-mirror change after the last HGNC ingest leaves the
    file in that state, as does a hand edit. One ERROR naming
    the missing ids (sorted, the first ten and then the total), never one per
    gene. A gene the file holds and the registry does not is not reported here.

    A frame without an ``hgnc_id`` column yields nothing: `validate_table`
    reports that as TBL001.
    """
    if "hgnc_id" not in genes.columns:
        return []
    present = {str(gene) for gene in genes["hgnc_id"].to_list() if gene is not None}
    missing = sorted(derived - present)
    if not missing:
        return []
    shown = ", ".join(missing[:_GEN006_SHOWN])
    more = f" (first {_GEN006_SHOWN} of {len(missing)})" if len(missing) > _GEN006_SHOWN else ""
    noun = "gene" if len(missing) == 1 else "genes"
    return [
        ValidationIssue(
            code="GEN006",
            severity=Severity.ERROR,
            location="mirrors/genes.tsv",
            message=(
                f"mirrors/genes.tsv lacks {len(missing)} {noun} the registry derives from "
                f"curation/chd_scope.yaml and the validity mirrors{more}: {shown}; "
                "re-ingest hgnc or run `chd-atlas ingest --all` (before hgnc is pinned, "
                "regenerate mirrors/genes.tsv with scripts/convert_hgnc.py); if the pinned "
                "HGNC package lacks the gene, pin a newer HGNC release"
            ),
        )
    ]
