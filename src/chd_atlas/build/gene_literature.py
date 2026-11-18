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
(``chd-curator add-publication``) carries the ``study_type`` the extracting agent
read, so ``publications.json`` withholds it while literature findings cite it
and none of them is approved, unless something else the build publishes cites
it; a record no finding cites publishes as before
(``build/literature.withheld_publications``; owner decision 3, superseding
ruling R32).

A finding that carries a ``review_status`` is published all the same, with its
mark and note (ruling R24): the mark is a published qualifier, never a
withdrawal.

``tallies`` counts the approved findings by *paper* (design, lesion group,
carrier origin, review mark). It never sums carriers across papers: two papers
may report the same probands, so a total would count people twice.

Findings are ordered by (publication year, numeric PMID, finding id) -- numeric,
because as strings ``PMID:9`` sorts after ``PMID:10`` (ruling R14).
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from typing import Any, Final, NamedTuple, TypedDict

from chd_atlas.corpus import Corpus
from chd_atlas.hgvs import parse_c, parse_p, protein_key
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.gene_literature import Finding, GeneLiteratureFile, synthesis_digest
from chd_atlas.models.review import ReviewStatus
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


class LiteratureTallies(TypedDict):
    """Papers (not findings) behind a gene's approved findings, by category.

    ``designs`` and ``lesion_groups`` list only categories with at least one paper.
    A paper is counted once per category however many findings it holds.
    """

    papers: int
    first_year: int | None
    last_year: int | None
    designs: dict[str, int]
    lesion_groups: dict[str, int]
    de_novo_papers: int
    inherited_papers: int
    abstract_only_papers: int
    pending_review_papers: int


class PublishedSynthesis(TypedDict):
    """A gene's approved "In short" narrative as published.

    ``sentences`` items are ``{"text": ..., "cites": [...]}`` with the cites in the
    curator's order; ``approved_on`` is ISO ``YYYY-MM-DD``. The approval digest is
    not published (as for findings), and no figure is computed here.
    """

    sentences: list[dict[str, Any]]
    drafted_by: str
    approved_by: str
    approved_on: str


class GeneLiterature(TypedDict):
    """The ``literature`` key of a gene bundle.

    ``synthesis`` is ``None`` when the gene has no synthesis, its synthesis is
    unapproved, or its approval no longer matches ``synthesis_digest``.
    """

    searches: list[SearchRecord]
    triage: TriageTally
    pending_approval: int
    findings: list[dict[str, Any]]
    normalisations: list[dict[str, Any]]
    tallies: LiteratureTallies
    synthesis: PublishedSynthesis | None


class LiteratureCensus(NamedTuple):
    """What the literature layer this site publishes amounts to (ruling R16)."""

    genes: int
    findings: int


def literature_tallies(findings: Sequence[Finding], years: Mapping[str, int]) -> LiteratureTallies:
    """Count ``findings`` by paper. ``years`` maps ``str(publication id)`` to its year."""
    papers: set[str] = set()
    designs: dict[str, set[str]] = {}
    lesions: dict[str, set[str]] = {}
    de_novo: set[str] = set()
    inherited: set[str] = set()
    abstract_only: set[str] = set()
    pending_review: set[str] = set()
    for finding in findings:
        paper = str(finding.publication)
        papers.add(paper)
        designs.setdefault(finding.design.value, set()).add(paper)
        for group in finding.lesion_groups:
            lesions.setdefault(group.value, set()).add(paper)
        if finding.carriers is not None:
            if finding.carriers.de_novo > 0:
                de_novo.add(paper)
            if finding.carriers.inherited > 0:
                inherited.add(paper)
        if finding.review_status is ReviewStatus.NOT_ACCESSIBLE:
            abstract_only.add(paper)
        elif finding.review_status is ReviewStatus.PENDING_REVIEW:
            pending_review.add(paper)
    paper_years = [years[paper] for paper in papers]
    return LiteratureTallies(
        papers=len(papers),
        first_year=min(paper_years) if paper_years else None,
        last_year=max(paper_years) if paper_years else None,
        designs={key: len(designs[key]) for key in sorted(designs)},
        lesion_groups={key: len(lesions[key]) for key in sorted(lesions)},
        de_novo_papers=len(de_novo),
        inherited_papers=len(inherited),
        abstract_only_papers=len(abstract_only),
        pending_review_papers=len(pending_review),
    )


EMPTY_GENE_LITERATURE: Final[GeneLiterature] = GeneLiterature(
    searches=[],
    triage=TriageTally(included=0, excluded={}, pending=0),
    pending_approval=0,
    findings=[],
    normalisations=[],
    tallies=literature_tallies([], {}),
    synthesis=None,
)


def pmid_number(pmid: str) -> int:
    return int(pmid.removeprefix("PMID:"))


def _published_synthesis(
    file: GeneLiteratureFile, assertions: Sequence[LesionAssertion]
) -> PublishedSynthesis | None:
    """The file's synthesis if it is approved and its approval is current, else ``None``.

    Withheld when there is none, when ``approved_by`` is unset, or when
    ``approved_digest`` differs from ``synthesis_digest`` over the synthesis, the
    file's findings and ``assertions`` (the gene's own). The last case is validation's SYN004, a
    warning, so the build can reach it and must not publish an edited text under an
    old approval.
    """
    synthesis = file.synthesis
    if (
        synthesis is None
        or synthesis.approved_by is None
        or synthesis.approved_on is None
        or synthesis.approved_digest != synthesis_digest(synthesis, file.findings, assertions)
    ):
        return None
    return PublishedSynthesis(
        sentences=[
            {"text": s.text, "cites": [str(c) for c in s.cites]} for s in synthesis.sentences
        ],
        drafted_by=synthesis.drafted_by,
        approved_by=synthesis.approved_by,
        approved_on=synthesis.approved_on.isoformat(),
    )


def gene_literature(corpus: Corpus) -> dict[str, GeneLiterature]:
    """Every gene's literature record, keyed by HGNC id, published or not.

    Restricting to the published genes is the callers' job (``build_genes``
    looks a gene up, ``literature_census`` takes ``published``): a file for an
    unpublished gene reaches no bundle, and validation warns (LIT013).

    Raises ``ValueError`` on a repeated gene or a finding whose publication is
    not registered. Both are refused by validation first (LIT009, LIT011); the
    guard is for a bypassed gate, and ``raise`` rather than ``assert`` because
    ``-O`` strips it.

    ``synthesis`` is published only for an approved synthesis whose digest still
    matches the synthesis, the file's findings and the gene's own assertions
    (``_published_synthesis``); otherwise it is ``None``.
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
            key=lambda f: (years[str(f.publication)], pmid_number(str(f.publication)), str(f.id))
        )
        # A run publishes only when a published variant cites it: an unapproved
        # finding's variants, and the runs only they cite, reach no file.
        cited = {
            variant.normalisation
            for finding in approved
            for variant in finding.variants
            if variant.normalisation is not None
        }
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
            findings=[
                finding.model_dump(mode="json", by_alias=True, exclude={"approved_digest"})
                for finding in approved
            ],
            normalisations=[
                run.model_dump(mode="json")
                for run in sorted(file.normalisations, key=lambda n: n.id)
                if run.id in cited
            ],
            tallies=literature_tallies(approved, years),
            synthesis=_published_synthesis(
                file, [a for a in corpus.assertions if str(a.gene) == gene]
            ),
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


class VariantReport(TypedDict):
    """One finding's report of a variant. ``carriers`` is that finding's own count.

    ``hgvs_c`` and ``hgvs_p`` are the names this report's variant carries. A row
    shows its first report's; a paper may name the same key differently (LIT018
    accepts shifted names), so the page renders a report's name where it differs
    from the row's.

    No count is ever summed across reports: two papers may report the same
    probands (D33), so a row carries reports and no total.
    """

    finding: str
    variant: int
    publication: str
    design: str
    hgvs_c: str | None
    hgvs_p: str | None
    count_unit: str | None
    carriers: dict[str, int] | None
    zygosity: str
    with_names: list[str]
    phase: str | None
    phenotypes: list[str]
    review_status: str | None
    review_note: str | None


class VariantRow(TypedDict):
    """One distinct variant: ``papers`` counts distinct publications among ``reports``."""

    key: str
    hgvs_c: str | None
    hgvs_p: str | None
    papers: int
    reports: list[VariantReport]


class ListedVariant(TypedDict):
    """A variant listed beside the table, never grouped."""

    as_reported: str
    rsid: str | None
    publication: str
    finding: str
    variant: int
    zygosity: str
    phase: str | None
    review_status: str | None
    review_note: str | None


class VariantTable(TypedDict):
    """A gene's approved variants for its page: rows, two lists and a count.

    ``no_position`` holds the ``structural`` variants, then the ``unplaced`` ones.
    ``findings_naming_no_variant`` counts findings of any design but a
    common-variant association that report carriers and list no variant.
    """

    rows: list[VariantRow]
    outside_region: list[ListedVariant]
    no_position: list[ListedVariant]
    findings_naming_no_variant: int


_REGION_RANK: Final[dict[str, int]] = {"five_prime": 0, "cds": 1, "three_prime": 2}


def _partner_name(variant: Mapping[str, Any]) -> str:
    return variant["hgvs_c"] or variant["hgvs_p"] or str(variant["as_reported"])


def _row_order(key: str, hgvs_c: str | None, hgvs_p: str | None) -> tuple[int, int, int, str]:
    """``(region rank, base, offset, key)``; a protein row sits at its codon's first base."""
    if key.startswith("g:"):
        if hgvs_c is None:
            raise ValueError(f"{key}: a mane row with no hgvs_c has no position to sort by")
        start = parse_c(hgvs_c).start
        return (_REGION_RANK[start.region], start.base, start.offset, key)
    if hgvs_p is None:
        raise ValueError(f"{key}: a protein row with no hgvs_p has no position to sort by")
    residues = parse_p(hgvs_p).residues
    if not residues:
        raise ValueError(f"{key}: {hgvs_p} names no residue to sort by")
    return (_REGION_RANK["cds"], 3 * residues[0][1] - 2, 0, key)


def variant_table(findings: Sequence[Mapping[str, Any]]) -> VariantTable:
    """Group approved published findings' variants into rows, one per distinct variant.

    ``findings`` are the bundle's ``literature.findings`` (already in (year, PMID,
    finding id) order, which is the order each row's reports keep). A ``mane``
    variant is keyed ``g:<grch38>`` and a ``protein`` one ``p:<hgvs_p>``, so a
    protein row is never merged into a mane row that shares its p. A protein key
    is the canonical reading of its ``hgvs_p`` (``protein_key``): accession and
    parentheses dropped, ``Ter`` and ``X`` as ``*``, no frameshift equivalence. Variant and
    ``with`` ids are scoped per finding. Raises ``ValueError`` on a variant a
    bypassed gate lets through: no key, an unreadable HGVS, a protein ``hgvs_p`` naming
    no residue (LIT022 refuses it) or a dangling partner.
    """
    rows: dict[str, VariantRow] = {}
    papers: dict[str, set[str]] = {}
    outside: list[ListedVariant] = []
    structural: list[ListedVariant] = []
    unplaced: list[ListedVariant] = []
    naming_none = 0
    for finding in findings:
        variants: Sequence[Mapping[str, Any]] = finding["variants"]
        carriers = finding.get("carriers")
        if (
            not variants
            and finding["design"] != "common_variant_association"
            and carriers is not None
            and sum(carriers.values()) > 0
        ):
            naming_none += 1
        by_id = {v["id"]: v for v in variants}
        for variant in variants:
            placement = variant["placement"]
            if placement in ("mane", "protein"):
                if placement == "mane":
                    if variant["grch38"] is None:
                        raise ValueError(
                            f"{finding['id']}: mane variant {variant['id']} has no grch38"
                        )
                    key = "g:" + variant["grch38"]
                else:
                    if variant["hgvs_p"] is None:
                        raise ValueError(
                            f"{finding['id']}: protein variant {variant['id']} has no hgvs_p"
                        )
                    key = "p:" + protein_key(variant["hgvs_p"])
                partners = []
                for partner in variant["with"]:
                    if partner not in by_id:
                        raise ValueError(
                            f"{finding['id']}: variant {variant['id']} is with {partner}, "
                            "which is not in the finding"
                        )
                    partners.append(_partner_name(by_id[partner]))
                row = rows.get(key)
                if row is None:
                    row = rows[key] = VariantRow(
                        key=key,
                        hgvs_c=variant["hgvs_c"],
                        hgvs_p=variant["hgvs_p"],
                        papers=0,
                        reports=[],
                    )
                    papers[key] = set()
                row["hgvs_c"] = row["hgvs_c"] or variant["hgvs_c"]
                row["hgvs_p"] = row["hgvs_p"] or variant["hgvs_p"]
                papers[key].add(finding["publication"])
                row["reports"].append(
                    VariantReport(
                        finding=finding["id"],
                        variant=variant["id"],
                        publication=finding["publication"],
                        design=finding["design"],
                        hgvs_c=variant["hgvs_c"],
                        hgvs_p=variant["hgvs_p"],
                        count_unit=finding.get("count_unit"),
                        carriers=None if variant["carriers"] is None else dict(variant["carriers"]),
                        zygosity=variant["zygosity"],
                        with_names=partners,
                        phase=variant["phase"],
                        phenotypes=list(variant["phenotypes"]),
                        review_status=variant["review_status"],
                        review_note=variant["review_note"],
                    )
                )
                continue
            listed = ListedVariant(
                as_reported=variant["as_reported"],
                rsid=variant["reported"]["rsid"],
                publication=finding["publication"],
                finding=finding["id"],
                variant=variant["id"],
                zygosity=variant["zygosity"],
                phase=variant["phase"],
                review_status=variant["review_status"],
                review_note=variant["review_note"],
            )
            if placement == "outside_region":
                outside.append(listed)
            elif placement == "structural":
                structural.append(listed)
            else:
                unplaced.append(listed)
    for key, row in rows.items():
        row["papers"] = len(papers[key])
    ordered = sorted(rows.values(), key=lambda r: _row_order(r["key"], r["hgvs_c"], r["hgvs_p"]))
    return VariantTable(
        rows=ordered,
        outside_region=outside,
        no_position=structural + unplaced,
        findings_naming_no_variant=naming_none,
    )
