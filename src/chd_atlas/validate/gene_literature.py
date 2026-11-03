# src/chd_atlas/validate/gene_literature.py
"""LIT001-LIT015: a gene's literature file against itself and the registries.

The single-record rules live in the model (`models/gene_literature.py`); these
are the ones a model cannot see because they compare a triage list with a
search list, a finding with a triage entry, a file name with a gene, or a
finding with another file.

Location format, mirroring `burden.py`'s ``path:row N``: the file's path, with
``:finding <id>`` or ``:pmid <PMID>`` appended where one record is at fault.
File-level issues carry the bare path.

New-hit counting (plan P4): a PMID is triaged once, attributed to the first
search that returned it, so a search's ``new_hits`` must equal the number of
triage entries attributed to it (LIT001). A later re-run that re-finds earlier
papers therefore adds nothing to the triage and changes no earlier approval.

Pure functions over loaded data, like `referential.py`. A registry passed as
`None` skips only the checks that need it; the caller pairs each skip with the
error that caused it.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

from chd_atlas.corpus import Corpus
from chd_atlas.genes import GeneRegistry, Resolution
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.gene_literature import GeneLiteratureFile, finding_digest
from chd_atlas.validate.referential import lesion_group_faults
from chd_atlas.vocab import LesionGroup, TriageDecision


def _approved_symbol_of(symbols: GeneRegistry, gene: str) -> str | None:
    for symbol in symbols.approved_symbols():
        if symbols.resolve(symbol).gene == gene:
            return symbol
    return None


def validate_gene_literature(
    corpus: Corpus,
    *,
    symbols: GeneRegistry | None,
    known_publications: set[str] | None,
    known_cohorts: set[str] | None,
    burden_studies: dict[str, set[str]] | None,
    published_genes: set[str] | None = None,
) -> list[ValidationIssue]:
    """Check every literature file; `None` for a registry skips only its checks.

    - `symbols` None skips LIT009 (TBL008 reports an unreadable registry).
    - `known_publications` None skips LIT011; `known_cohorts` None skips LIT012.
    - `published_genes` None skips LIT013.
    - `burden_studies` (gene -> the studies `mirrors/burden.tsv` has rows for) None
      skips LIT014; it is None only when that table is unreadable, which TBL000 or
      TBL001 reports. A missing table is `{}`: a finding linking a burden study then
      has nothing to link to.

    Lesion-group consistency (REF007/REF009, the assertion rules) runs on every
    finding against `corpus.phenotypes`; an unmapped phenotype is exempt.
    """
    issues: list[ValidationIssue] = []
    all_finding_ids: set[str] = {
        finding.id for _, literature in corpus.literature for finding in literature.findings
    }

    lesion_group_of: dict[str, LesionGroup | None] = {
        term.id: term.lesion_group for term in corpus.phenotypes
    }

    for path, literature in corpus.literature:
        _check_file(
            issues,
            path,
            literature,
            all_finding_ids,
            lesion_group_of,
            symbols,
            known_publications,
            known_cohorts,
            burden_studies,
            published_genes,
        )
    return issues


def _check_file(
    issues: list[ValidationIssue],
    path: Path,
    literature: GeneLiteratureFile,
    all_finding_ids: set[str],
    lesion_group_of: dict[str, LesionGroup | None],
    symbols: GeneRegistry | None,
    known_publications: set[str] | None,
    known_cohorts: set[str] | None,
    burden_studies: dict[str, set[str]] | None,
    published_genes: set[str] | None,
) -> None:
    def add(code: str, severity: Severity, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, severity, location, message))

    def error(code: str, location: str, message: str) -> None:
        add(code, Severity.ERROR, location, message)

    here = str(path)
    searches = {search.id: search for search in literature.searches}

    # LIT009: one literature file per gene, named by the approved symbol.
    if symbols is not None:
        resolved = symbols.resolve(path.stem)
        if resolved.resolution is not Resolution.APPROVED or resolved.gene != literature.gene:
            approved = _approved_symbol_of(symbols, literature.gene)
            named = f"; the file for {literature.gene} must be {approved}.yaml" if approved else ""
            error(
                "LIT009",
                here,
                f"file name {path.stem!r} is not the approved symbol of {literature.gene}"
                f" (an alias or previous symbol does not count){named}",
            )

    if published_genes is not None and literature.gene not in published_genes:
        add(
            "LIT013",
            Severity.WARNING,
            here,
            f"curated literature for {literature.gene} reaches no page: "
            "the gene is outside the publication gate",
        )

    # Triage: one entry per PMID (LIT005), attributed to a real search (LIT010).
    for pmid, count in sorted(Counter(entry.pmid for entry in literature.triage).items()):
        if count > 1:
            error("LIT005", f"{here}:pmid {pmid}", f"{pmid} is triaged {count} times; once only")
    for entry in literature.triage:
        if entry.search not in searches:
            error(
                "LIT010",
                f"{here}:pmid {entry.pmid}",
                f"triage names search {entry.search}, which this file does not record",
            )

    # LIT001: the arithmetic of each search.
    attributed = Counter(entry.search for entry in literature.triage)
    for search in sorted(literature.searches, key=lambda s: s.id):
        if attributed[search.id] != search.new_hits:
            error(
                "LIT001",
                here,
                f"search {search.id} found {search.new_hits} new paper(s) but "
                f"{attributed[search.id]} triage entries are attributed to it",
            )

    pending = sum(1 for entry in literature.triage if entry.decision is TriageDecision.PENDING)
    if pending:
        add(
            "LIT002",
            Severity.WARNING,
            here,
            f"{pending} hit(s) still pending triage",
        )

    included = {e.pmid for e in literature.triage if e.decision is TriageDecision.INCLUDED}
    cited = {finding.publication for finding in literature.findings}
    for pmid in sorted(included - cited):
        error("LIT003", f"{here}:pmid {pmid}", f"{pmid} is included but no finding cites it")
    for entry in literature.triage:
        if entry.same_as is not None and entry.same_as not in included:
            error(
                "LIT006",
                f"{here}:pmid {entry.pmid}",
                f"same_as names {entry.same_as}, which is not included in this file",
            )

    first_search: dict[str, int] = {}
    for entry in literature.triage:
        first_search.setdefault(entry.pmid, entry.search)

    for finding in literature.findings:
        where = f"{here}:finding {finding.id}"
        if finding.publication not in included:
            error(
                "LIT004",
                where,
                f"cites {finding.publication}, which is not included in this file's triage",
            )
        if isinstance(finding.shares_probands_with, list):
            for other in finding.shares_probands_with:
                if other not in all_finding_ids:
                    error("LIT007", where, f"shares_probands_with names unknown finding {other}")
                elif other == finding.id:
                    error("LIT007", where, "shares_probands_with names the finding itself")
        for code, message in lesion_group_faults(
            finding.phenotypes, finding.lesion_groups, lesion_group_of
        ):
            error(code, f"{here}:finding {finding.id}", message)
        if (
            burden_studies is not None
            and finding.burden_study is not None
            and finding.burden_study not in burden_studies.get(literature.gene, set())
        ):
            error(
                "LIT014",
                where,
                f"burden_study {finding.burden_study} has no row for {literature.gene} "
                f"in mirrors/burden.tsv, so there is nothing for the page to link",
            )
        if finding.approved_by is not None and finding.approved_digest != finding_digest(finding):
            error(
                "LIT015",
                where,
                "approved content changed since approval; re-approve",
            )
        found_by = searches.get(first_search.get(finding.publication, -1))
        if (
            finding.approved_on is not None
            and found_by is not None
            and finding.approved_on < found_by.run_on
        ):
            error(
                "LIT008",
                where,
                f"approved on {finding.approved_on}, before search {found_by.id} "
                f"({found_by.run_on}) found {finding.publication}",
            )
        if known_publications is not None:
            for pmid in sorted(
                {finding.publication, *([finding.burden_study] if finding.burden_study else [])}
            ):
                if pmid not in known_publications:
                    error("LIT011", where, f"{pmid} is not in curation/publications.yaml")
        if known_cohorts is not None:
            for cohort in sorted(set(finding.cohorts) - known_cohorts):
                error("LIT012", where, f"cohort {cohort} is not in curation/cohorts.yaml")
