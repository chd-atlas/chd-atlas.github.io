# src/chd_atlas/build/gene_literature.py
"""A gene's curated literature record, as the bundle publishes it.

Named ``gene_literature`` because ``build/literature.py`` already builds the
publications, featured and phenotype registries (plan decision P1).

Three things are published for a gene with a ``curation/literature`` file: the
searches that were run (``source``, ``query``, ``run_on``, ``hits``,
``new_hits``); the triage *tally* -- how many hits were included, excluded by
reason, or still pending, never the per-paper decisions; and the findings a
named person approved. A finding that has not been approved is counted in
``pending_approval`` and its content reaches no byte: publishing it would put
an agent's unreviewed extraction in front of a reader as curated evidence. The
count is what keeps that omission from being silent.

A ``Publication`` that a curator added for a finding's paper
(``literature add-publication``) carries the ``study_type`` the extracting agent
read, so ``publications.json`` withholds it until a finding of that paper is
approved, unless something else the build publishes cites it
(``build/literature.py``; owner decision 3, superseding ruling R32).

A finding that carries a ``review_status`` is published all the same, with its
mark and note (ruling R24): the mark is a published qualifier, never a
withdrawal.

Findings are ordered by (publication year, numeric PMID, finding id) -- numeric,
because as strings ``PMID:9`` sorts after ``PMID:10`` (ruling R14).
"""

from __future__ import annotations

from collections.abc import Collection, Mapping
from typing import Any, Final, NamedTuple, TypedDict

from chd_atlas.corpus import Corpus
from chd_atlas.models.gene_literature import Finding
from chd_atlas.vocab import TriageDecision


class SearchRecord(TypedDict):
    """One recorded search, as published."""

    source: str
    query: str
    run_on: str
    hits: int
    new_hits: int


class TriageTally(TypedDict):
    """What was decided about the hits, counted. ``excluded`` is keyed by reason."""

    included: int
    excluded: dict[str, int]
    pending: int


class GeneLiterature(TypedDict):
    """The ``literature`` key of a gene bundle."""

    searches: list[SearchRecord]
    triage: TriageTally
    pending_approval: int
    findings: list[dict[str, Any]]


class LiteratureCensus(NamedTuple):
    """What the literature layer this site publishes amounts to (ruling R16)."""

    genes: int
    findings: int


EMPTY_GENE_LITERATURE: Final[GeneLiterature] = GeneLiterature(
    searches=[],
    triage=TriageTally(included=0, excluded={}, pending=0),
    pending_approval=0,
    findings=[],
)


def _pmid_number(pmid: str) -> int:
    return int(pmid.removeprefix("PMID:"))


def gene_literature(corpus: Corpus) -> dict[str, GeneLiterature]:
    """Every gene's literature record, keyed by HGNC id, published or not.

    Restricting to the published genes is the callers' job (``build_genes``
    looks a gene up, ``literature_census`` takes ``published``): a file for an
    unpublished gene reaches no bundle, and validation warns (LIT013).

    Raises ``ValueError`` on a repeated gene or a finding whose publication is
    not registered. Both are refused by validation first (LIT009, LIT011); the
    guard is for a bypassed gate, and ``raise`` rather than ``assert`` because
    ``-O`` strips it.
    """
    years = {str(publication.id): publication.year for publication in corpus.publications}
    result: dict[str, GeneLiterature] = {}
    for path, file in corpus.literature:
        gene = str(file.gene)
        if gene in result:
            raise ValueError(f"{path}: a second literature file for {gene}")
        excluded: dict[str, int] = {}
        included = pending = 0
        for entry in file.triage:
            if entry.decision is TriageDecision.INCLUDED:
                included += 1
            elif entry.decision is TriageDecision.PENDING:
                pending += 1
            elif entry.reason is not None:
                excluded[entry.reason.value] = excluded.get(entry.reason.value, 0) + 1
            else:
                raise ValueError(f"{path}: {entry.pmid} is excluded with no reason")
        approved: list[Finding] = []
        unapproved = 0
        for finding in file.findings:
            if finding.approved_by is None:
                unapproved += 1
                continue
            if str(finding.publication) not in years:
                raise ValueError(
                    f"{path}: {finding.id} cites {finding.publication}, which is not registered"
                )
            approved.append(finding)
        approved.sort(
            key=lambda f: (years[str(f.publication)], _pmid_number(str(f.publication)), str(f.id))
        )
        result[gene] = GeneLiterature(
            searches=[
                SearchRecord(
                    source=search.source,
                    query=search.query,
                    run_on=search.run_on.isoformat(),
                    hits=search.hits,
                    new_hits=search.new_hits,
                )
                for search in sorted(file.searches, key=lambda s: s.id)
            ],
            triage=TriageTally(
                included=included,
                excluded={reason: excluded[reason] for reason in sorted(excluded)},
                pending=pending,
            ),
            pending_approval=unapproved,
            findings=[finding.model_dump(mode="json") for finding in approved],
        )
    return result


def literature_census(
    literature: Mapping[str, GeneLiterature], published: Collection[str]
) -> LiteratureCensus:
    """Genes with a recorded search and approved findings, over the published genes only.

    ``genes`` counts published genes whose record holds at least one search;
    ``findings`` the approved findings in those same bundles. Counting the
    mapping would advertise literature behind genes this site publishes no
    page for.
    """
    genes = findings = 0
    for gene in published:
        entry = literature.get(gene)
        if entry is None or not entry["searches"]:
            continue
        genes += 1
        findings += len(entry["findings"])
    return LiteratureCensus(genes=genes, findings=findings)
