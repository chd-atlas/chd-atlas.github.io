# src/chd_atlas/build/validity.py
"""Gene-disease validity, mirrored and attributed -- never authored here.

Design decision D12: the atlas publishes no validity classification of its own. One curator cannot
out-curate a chartered expert panel, and the resource that tried is the
cautionary case: CHDgene publishes a documented 3-tier ranking under which all
189 genes carry 5/5 stars.

What this module does instead is join the two mirrors -- `mirrors/clingen_gene_
validity.tsv` and `mirrors/gencc_submissions.tsv` -- to `curation/chd_scope.yaml`
and say, per gene, which of three things is true: an expert panel has assessed
it (`EXPERT_CURATED`), only submitters have (`SUBMITTER_CURATED`), or nobody has
(`UNCURATED`). The third is published as a fact rather than left as an absent
key, because a consumer cannot tell "no authority has assessed this gene" from
"the build dropped it" by looking at a missing field -- and that distinction is
exactly the silent evidence loss `CLAUDE.md` calls this project's characteristic
failure.

`has_source_discordance` is decided and returned alongside `records` in the same
`GeneValidity`, for the reason `bundles.py::_headline` writes
`headline_confidence` and `has_conflicting_evidence` together: a payload must
never be able to carry the records without the flag that says whether they
disagree.

Every published array here is sorted by the caller. `encode_json` (`build/
emit.py`) passes `sort_keys=True`, which orders dict *keys* and nothing else --
an array built by iterating a `dict` populated one mirror row at a time checksums
differently between processes unless this module sorts it before returning.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Final

import polars as pl

from chd_atlas.build.emit import Json
from chd_atlas.vocab import (
    CLASSIFICATION_RANK,
    CLINGEN_CLASSIFICATIONS,
    CONTESTED,
    GENCC_CLASSIFICATIONS,
    Classification,
    ValiditySource,
    ValidityState,
)


@dataclass(frozen=True)
class ValidityRecord:
    """One authority's verdict on one gene-disease pair.

    `classification_term` is what the authority published, verbatim.
    `classification` is that term mapped onto the atlas's `Classification`
    ladder, or `None` where the term is not a rung on it at all -- GenCC's
    `Supportive`, the mapping exception for a submitter that asserts an
    association without grading its strength. Publishing both lets a consumer
    render the authority's own word and still filter on a common scale, without
    the atlas inventing a confidence nobody stated.

    The last six fields are authority-specific: `sop`, `classification_date`
    and `gcep` come from ClinGen only; `submitter`, `submitted_on` and `sgc_id`
    from GenCC only. Each is `None` on a record from the other source rather
    than the two sources getting separate record types, because every consumer
    of `records` wants one homogeneous list to sort, filter and render, not a
    union type to match on first.

    **`submitted_on` and `sgc_id` reached no published byte until schema
    2.16.** They were a recorded scope decision rather than a dropped join, and
    the cost was measured before it was paid off: on 2026-10-25, 470 of the 591
    validity records the 119 published bundles carry were GenCC, every one
    with a non-null `submitted_on` and `sgc_id` in the mirror and a blank date
    on its page. `submitted_on` is carried verbatim -- the mirror holds three
    timestamp shapes (`2018-05-31 20:04:32`, `2025-01-17T00:00:00.000000Z` and
    a millisecond variant) -- because normalising a date for display is the
    page's job (`pages._date`), and a field that rewrote the upstream value
    would stop being a mirror of it.
    """

    source: ValiditySource
    classification_term: str
    classification: Classification | None
    disease: str
    disease_label: str
    moi: str
    report_url: str | None
    sop: str | None = None
    classification_date: str | None = None
    gcep: str | None = None
    submitter: str | None = None
    submitted_on: str | None = None
    sgc_id: str | None = None


@dataclass(frozen=True)
class GeneValidity:
    """One gene's whole validity picture.

    How curated it is (`state`), by whom and under what verdict (`records`),
    and whether the mirrors agree it is settled (`has_source_discordance`).
    """

    state: ValidityState
    records: tuple[ValidityRecord, ...]
    has_source_discordance: bool


def uncurated() -> GeneValidity:
    """The published shape for a gene no authority has assessed at all.

    Neither an expert panel nor a submitter.

    A function, not a module-level constant: a shared `GeneValidity` instance
    handed to every uncurated gene would still be one object underneath
    "frozen", and nothing stops a future caller mutating through it or relying
    on identity. Returning a fresh instance per call costs nothing and
    forecloses the question.
    """
    return GeneValidity(
        state=ValidityState.UNCURATED,
        records=(),
        has_source_discordance=False,
    )


def _sort_key(record: ValidityRecord) -> tuple[str, str, str, str, str, str]:
    """`(source, disease, moi, submitter or "", sgc_id or "", gcep or "")`.

    So two builds of one commit agree on array order regardless of which
    mirror row was read first. `submitter` is `None` on every ClinGen record;
    `or ""` sorts those first, which never has to arbitrate a real tie because
    a null-vs-populated `submitter` split never happens within one source --
    only GenCC rows carry it, and every GenCC row has one. The fifth element
    breaks the one published tie, PKD1L1's two Ambry Genetics submissions for
    MONDO:0014967, whose order otherwise followed mirror row order; measured
    2026-10-25, adding it reorders zero published records -- `gene_validity`'s
    output for all 185 registry genes sorted under both keys reorders none, and
    a byte diff of two real builds differed only in `manifest.json`'s
    `source_commit`. The sixth orders two ClinGen panels' records on one
    (gene, disease, moi), which carry neither a submitter nor an sgc_id
    (`clingen_validity` v2, F28: FOXE3 / MONDO:0019625 / AD in ClinGen's
    2026-11-13 file). No published gene carries such a pair (measured
    2026-11-13), so it reorders no published record.
    """
    return (
        record.source.value,
        record.disease,
        record.moi,
        record.submitter or "",
        record.sgc_id or "",
        record.gcep or "",
    )


def _clingen_records(frame: pl.DataFrame, in_scope: set[str]) -> dict[str, list[ValidityRecord]]:
    by_gene: dict[str, list[ValidityRecord]] = {}
    for row in frame.to_dicts():
        if row["disease"] not in in_scope:
            continue
        term = row["classification"]
        if term not in CLINGEN_CLASSIFICATIONS:
            raise ValueError(
                f"unmapped ClinGen classification {term!r} for gene {row['gene']}; "
                f"add it to vocab.CLINGEN_CLASSIFICATIONS rather than coercing it "
                f"to an adjacent rung"
            )
        by_gene.setdefault(row["gene"], []).append(
            ValidityRecord(
                source=ValiditySource.CLINGEN,
                classification_term=term,
                classification=CLINGEN_CLASSIFICATIONS[term],
                disease=row["disease"],
                disease_label=row["disease_label"],
                moi=row["moi"],
                report_url=row["report_url"],
                sop=row["sop"],
                classification_date=row["classification_date"],
                gcep=row["gcep"],
            )
        )
    return by_gene


def _gencc_records(frame: pl.DataFrame, in_scope: set[str]) -> dict[str, list[ValidityRecord]]:
    by_gene: dict[str, list[ValidityRecord]] = {}
    for row in frame.to_dicts():
        if row["disease"] not in in_scope:
            continue
        term = row["classification"]
        if term not in GENCC_CLASSIFICATIONS:
            raise ValueError(
                f"unmapped GenCC classification {term!r} for gene {row['gene']}; "
                f"add it to vocab.GENCC_CLASSIFICATIONS rather than coercing it "
                f"to an adjacent rung"
            )
        by_gene.setdefault(row["gene"], []).append(
            ValidityRecord(
                source=ValiditySource.GENCC,
                classification_term=term,
                classification=GENCC_CLASSIFICATIONS[term],
                disease=row["disease"],
                disease_label=row["disease_label"],
                moi=row["moi"],
                report_url=row["report_url"],
                submitter=row["submitter"],
                submitted_on=row["submitted_on"],
                sgc_id=row["sgc_id"],
            )
        )
    return by_gene


def _has_source_discordance(records: tuple[ValidityRecord, ...]) -> bool:
    """True iff some source contests the gene while a *different* source supports it.

    Not "the classifications differ": Definitive and Strong differ in rank but
    both sit outside `CONTESTED`, so two sources reporting one each still agree
    the gene is not contested. Not "one source alone carries both a supportive
    and a contesting call" either -- 90 genes in the committed ClinGen mirror do
    exactly that across two diseases or two GCEPs (measured 2026-08-03), and
    that is ClinGen being internally split, not ClinGen disagreeing with
    anyone. (None of those 90 splits survives the CHD-scope filter this
    module applies before a gene is even keyed -- 13 of the 90 have some
    in-scope ClinGen record, but never both sides of the split within scope,
    so this is a definitional justification, not a case seen in published
    data; measured 2026-08-04 against `gene_validity`'s real output, where
    exactly one gene sets `has_conflicting_evidence` -- LEFTY2, HGNC:3122, via
    a *cross-source* split -- and it also sets `has_source_discordance`.) Only
    a record with a mapped `classification` takes a side: GenCC's
    `Supportive` maps to `None` and is a submitter declining to grade evidence,
    not evidence *for* the gene, and `NO_KNOWN_ASSOCIATION` is a stated absence
    of association, not support for one.

    Implemented as the literal reading of "one authority contesting a gene
    another supports": true iff some `c` in the contesting sources and some `s`
    in the supporting sources are different sources. Both sets are subsets of
    `{ValiditySource.CLINGEN, ValiditySource.GENCC}` here, so on every case this
    module's own tests exercise this agrees with the simpler `contested !=
    supporting`; the pairwise form is kept because it also stays correct in the
    one case the two diverge -- both sources internally split identically
    (`contested == supporting == {CLINGEN, GENCC}`), where each source still
    contests something a *different* source supports, and the simpler
    set-inequality check would miss it.
    """
    contesting = {
        record.source
        for record in records
        if record.classification is not None and record.classification in CONTESTED
    }
    supporting = {
        record.source
        for record in records
        if record.classification is not None
        and record.classification not in CONTESTED
        and record.classification is not Classification.NO_KNOWN_ASSOCIATION
    }
    return any(c != s for c in contesting for s in supporting)


def gene_validity(
    clingen: pl.DataFrame,
    gencc: pl.DataFrame,
    in_scope: set[str],
) -> dict[str, GeneValidity]:
    """One `GeneValidity` per gene either mirror curates within CHD scope.

    `in_scope` is the set of MONDO disease terms from `curation/chd_scope.yaml`
    (`Corpus.chd_scope`) -- the atlas's *selection* among disease terms external
    authorities already treat as congenital heart disease, never a claim authored
    here (`models/scope.py`; every term names its authority and SCP005 checks
    that name against the mirrors). A row naming any other disease
    is dropped here, before a gene is even keyed, rather than filtered
    downstream: GATA6 is Definitive for GATA6-related congenital heart disease
    (`MONDO:0100540`) and Limited for dilated cardiomyopathy (`MONDO:0005021`,
    not CHD), both current and both curated by different ClinGen expert panels
    (measured 2026-08-03 against the committed mirror). Admitting the second row
    would make the gene read as contested *about CHD* when the two verdicts
    concern different diseases entirely.

    A gene reached through any in-scope ClinGen record is `EXPERT_CURATED`, even
    when GenCC has also curated it -- an expert panel's own verdict does not
    change because submitters have opinions too, and this module publishes
    every in-scope GenCC record regardless, so nothing about the field's wider
    view is hidden by the state a gene is filed under. A gene with only in-scope
    GenCC records is `SUBMITTER_CURATED`. A gene reached through neither mirror
    at all is not a key in the returned dict -- callers wanting the third state
    use `uncurated()` for a gene they otherwise know about (e.g. every gene in
    `mirrors/genes.tsv`), the same way `derive.gene_facts` omits a gene with no
    assertion rather than emitting an empty one.

    An unmapped classification term raises `ValueError` naming the term and the
    gene, rather than falling back to any adjacent rung: a term either mirror
    adds in a later release must reach a human, per `CLAUDE.md`'s standing rule
    that a bypassed gate must fail rather than publish a confidence nobody
    asserted.
    """
    from_clingen = _clingen_records(clingen, in_scope)
    from_gencc = _gencc_records(gencc, in_scope)

    result: dict[str, GeneValidity] = {}
    # Sorted so the dict's own key order -- and so anything a caller builds by
    # iterating it into a JSON array or object -- does not depend on which
    # mirror happened to name a gene first.
    for gene in sorted(set(from_clingen) | set(from_gencc)):
        records = tuple(
            sorted(from_clingen.get(gene, []) + from_gencc.get(gene, []), key=_sort_key)
        )
        result[gene] = GeneValidity(
            state=(
                ValidityState.EXPERT_CURATED
                if gene in from_clingen
                else ValidityState.SUBMITTER_CURATED
            ),
            records=records,
            has_source_discordance=_has_source_discordance(records),
        )
    return result


# The weakest ClinGen rung that admits a gene on ClinGen's word alone.
#
# **A rank floor, never a set of admissible rungs.** `DISPUTED` (rank 2),
# `REFUTED` (1) and `NO_KNOWN_ASSOCIATION` (0) sit below `LIMITED` (3), so they
# are excluded by construction rather than by an enumeration a later edit can
# forget. ClinGen treats contested as a *separate axis* rather than a weaker
# rung, and a linear rank is exactly where that distinction is easy to lose:
# lowering this one step admits 24 genes an expert panel actively disputes and
# publishes them as supported.
PUBLICATION_FLOOR: Final = Classification.LIMITED

# How many GenCC submitters must independently assert a gene for it to be
# admitted without ClinGen. Two, not one: a single submitter is one laboratory's
# panel-inclusion decision, and `Supportive` maps to `None` so an ungraded
# assertion counts for nobody.
SUBMITTER_AGREEMENT: Final = 2


def _asserts(record: ValidityRecord) -> bool:
    """True iff this record is its author *asserting* a gene-disease relationship.

    **A dissent is not an assertion, and `asserted_by` is a count of assertions.**
    `NO_KNOWN_ASSOCIATION` is an authority saying it looked and found no
    relationship; `DISPUTED` and `REFUTED` are it pushing back on one. Counting
    any of the three inflates `asserted_by.count` with the institutions that
    disagree, which is the opposite of what the field's name claims. Measured
    2026-08-06 before this filter existed: GDF1 published `"count": 6` including
    Illumina's `No Known Disease Relationship` -- and GDF1 is the gene every
    docstring here cites as *why* two submitters are required.

    `None` **does** assert. GenCC's `Supportive` maps to `None`
    (`vocab.GENCC_CLASSIFICATIONS`) and means the submitter asserted an
    association without grading its evidence -- an assertion, merely ungraded.
    Excluding it would drop a submitter from every gene its only in-scope record
    is a `Supportive` one for, and leave a gene whose sole in-scope record is
    such a submission asserted by nobody.
    """
    return record.classification not in CONTESTED and (
        record.classification is not Classification.NO_KNOWN_ASSOCIATION
    )


def _admitting_clingen(entry: GeneValidity) -> ValidityRecord | None:
    """The strongest ClinGen record at or above the floor, or None."""
    qualifying = [
        record
        for record in entry.records
        if record.source is ValiditySource.CLINGEN
        and record.classification is not None
        and CLASSIFICATION_RANK[record.classification] >= CLASSIFICATION_RANK[PUBLICATION_FLOOR]
    ]
    if not qualifying:
        return None
    # `records` is already sorted; `max` is stable, so ties resolve to the
    # earliest in that order and two builds agree. The `or` is unreachable --
    # every member of `qualifying` was filtered on a non-None classification --
    # and exists because the narrowing does not survive into the lambda.
    return max(
        qualifying,
        key=lambda record: CLASSIFICATION_RANK[record.classification or PUBLICATION_FLOOR],
    )


def admitting_grade(entry: GeneValidity) -> Classification | None:
    """The grade of the ClinGen record that admitted this gene, or `None`.

    `admitted_by.classification` publishes the same record's grade, and the
    landing page counts panel-graded genes with this.

    **`derive.gene_facts` does not headline this**, and this docstring said it
    did until 2026-10-26. The headline is `strongest()` over every in-scope
    ClinGen record with a mapped grade; for a gene ClinGen admitted that is
    this record, so the two agree on every published gene (measured 2026-10-26,
    0 of 119 differ). They part on an input the gate already admits:
    `NO_KNOWN_ASSOCIATION` is deliberately not a veto (`_clingen_contests`), so
    a gene ClinGen recorded it for, with two GenCC submitters agreeing, is
    published with `None` here and `no_known_association` in the headline.
    Keying the headline on this function was tried and reverted for exactly
    that gene: a `null` headline renders "no ClinGen expert panel has graded
    this gene", false for a gene a panel graded -- see the comment above
    `panel_graded` in `derive.gene_facts`. No such gene is in the committed
    mirrors.
    """
    record = _admitting_clingen(entry)
    return record.classification if record is not None else None


def agreeing_submitters(entry: GeneValidity) -> tuple[str, ...]:
    """GenCC submitters asserting this gene at or above the floor, sorted.

    **ClinGen's own GenCC submissions are excluded, and that exclusion is the
    point.** GenCC aggregates rather than adjudicates, and measured 2026-08-06
    its largest in-scope submitter is ClinGen itself -- 111 rows over 109 genes.
    Counting those would let one body vote twice: once as the source this gate
    already trusts alone, and again as corroboration of itself. That is the same
    error `build/concordance.py::cohort_families` exists to prevent for studies
    sharing a sample collection, one layer up.

    Sorted, because it reaches published JSON and `encode_json`'s `sort_keys`
    orders dict keys only.
    """
    return tuple(
        sorted({record.submitter for record in agreeing_records(entry) if record.submitter})
    )


def agreeing_records(entry: GeneValidity) -> tuple[ValidityRecord, ...]:
    """The GenCC records whose submitters `agreeing_submitters` counts, in `records` order.

    The rule lives here once, because two consumers read it: the gate, which
    counts the submitters, and the gene page's card, which states the mode of
    inheritance those submitters' qualifying records carry (D55). A submitter's
    *other* records -- a `Supportive` one, or one for a second disease below the
    floor -- are not what admitted the gene and are not returned.
    """
    return tuple(
        record
        for record in entry.records
        if record.source is ValiditySource.GENCC
        and record.submitter
        and record.submitter != "ClinGen"
        and record.classification is not None
        and CLASSIFICATION_RANK[record.classification] >= CLASSIFICATION_RANK[PUBLICATION_FLOOR]
    )


def _clingen_contests(entry: GeneValidity) -> bool:
    """True iff a ClinGen expert panel disputed or refuted this gene in scope.

    **A ClinGen contested record vetoes admission on GenCC agreement**, and
    without this the gate is asymmetric in the one direction that matters:
    ClinGen is trusted to admit a gene alone, so it must also be trusted to
    refuse one. Otherwise two submitters overrule a chartered panel, and the
    atlas publishes a gene the panel looked at and pushed back on.

    LEFTY2 is the case. Measured 2026-08-06: ClinGen's Congenital Heart Disease
    GCEP records it `Disputed`; G2P and PanelApp Australia both record `Limited`.
    Without this veto it is admitted on their agreement -- and before the
    headline became panel-only it published as `limited`, burying the dispute
    entirely.

    `NO_KNOWN_ASSOCIATION` is deliberately NOT a veto. It belongs to neither side
    of `has_conflicting_evidence`'s test (spec D34) because a stated absence of
    association is not a contest, and treating it as one here would give a
    single "we found nothing" submission the force of a refutation.
    """
    return any(
        record.source is ValiditySource.CLINGEN
        and record.classification is not None
        and record.classification in CONTESTED
        for record in entry.records
    )


def published_genes(validity: Mapping[str, GeneValidity]) -> set[str]:
    """The genes the atlas publishes a page for. Design decision D21.

    A gene qualifies on **either** of two independent warrants:

    - a ClinGen record at or above `PUBLICATION_FLOOR`, or
    - `SUBMITTER_AGREEMENT` distinct GenCC submitters, ClinGen excluded.

    Every record reaching here already names a disease in
    `curation/chd_scope.yaml` -- `gene_validity` filters on that -- so this does
    not re-check scope. Which authority admitted a gene is published per gene by
    `admission_provenance` below; the gate is not a claim a reader has to take
    on trust.

    **The rule widened on 2026-08-06**, from ClinGen `Definitive` alone, on the
    owner's decision that the atlas should not be the authority on which genes
    are CHD genes and should admit on agreement with external sources. Measured
    against the committed mirrors: 23 genes -> **92**, and burden rows reaching a
    page 290 -> **915**. Those are this function's own return, after the veto
    below; the two warrants alone select 93 genes and 916 rows, and this
    docstring recorded that larger pair until 2026-08-06, describing the gate
    without the refusal the same function applies. Only ClinGen and GenCC
    participate, because they are the
    only two sources whose terms permit republishing (CC0-1.0 both; see issues
    #31 and #32 for the others).

    **ClinGen alone suffices and GenCC alone does not**, and the asymmetry is
    deliberate: ClinGen is the only source here with chartered expert panels and
    published SOPs, while GenCC harmonises submissions rather than adjudicating
    between them. GDF1 is the case that fixes the shape of the problem -- its
    in-scope submissions run from G2P's `Definitive` to Illumina's `No Known
    Disease Relationship`, and `has_conflicting_evidence` reports `False` for
    that pair because `no_known_association` takes neither side (spec D34). A
    one-submitter gate would publish that as settled; two submitters is the
    weakest bar that requires a second body to agree.

    Requiring one of those two to be non-commercial was measured and rejected:
    it drops **TBX1** -- the 22q11.2 gene, which ClinGen has no in-scope record
    for at all -- along with CFC1, FLNA and SEMA3E. Publishing the provenance is
    what makes the looser bar safe: a reader sees that TBX1 rests on two
    clinical laboratories and judges it themselves.

    **Not** "headline confidence clears the floor and validity state is expert
    curated." `state` records only that ClinGen has *a* row for the gene, never
    what that row says, so that rule admits a gene ClinGen graded below the floor
    on the strength of a submitter grading it higher. `test_build_validity.py`
    constructs one.
    """
    return {
        gene
        for gene, entry in validity.items()
        if _admitting_clingen(entry) is not None
        or (len(agreeing_submitters(entry)) >= SUBMITTER_AGREEMENT and not _clingen_contests(entry))
    }


def admission_provenance(entry: GeneValidity) -> dict[str, Json]:
    """Why this gene is published, as a payload. `admitted_by` + `asserted_by`.

    **The gate is a claim about external authorities, and this is what lets a
    consumer check it rather than trust it.** `methods.SCOPE_RULE` tells a
    reader that no disease is in scope on the atlas's own judgement; without
    this field the same reader has no way to see which authority actually
    admitted the gene in front of them.

    `admitted_by` is the single warrant that cleared the gate:

    - `{"authority": "clingen", "classification": ..., "disease": ...,
       "disease_label": ..., "panel": ...}` where a ClinGen record at or above
      `PUBLICATION_FLOOR` exists -- the strongest one, since a gene may carry
      several.
    - `{"authority": "gencc_agreement", "submitters": [...]}` otherwise.

    An object rather than an array because **exactly one warrant admits a
    gene**, and ClinGen is checked first: a gene with both is published on
    ClinGen's word, and saying so is the difference between "an expert panel
    graded this" and "two laboratories put it on a panel".

    `asserted_by` is every distinct institution asserting the gene in scope,
    **deduped by institution and filtered to actual assertions** (`_asserts`),
    with a `count`. Measured 2026-08-06 over the 23 genes published before the
    widening: counting `gcep` and `submitter` values naively gives **134**,
    deduped gives **111** -- an overstatement of exactly one per gene, because
    ClinGen submits to GenCC under its own name. (This docstring said 132 and
    109 until 2026-08-06; the difference of one per gene was right, both
    absolutes were not.) A consumer computing agreement from `validity.records`
    gets that wrong on every gene, which is the reason this is derived here
    rather than left to them.

    **`count` is not a score, and must never be rendered as one.** Eight
    authorities asserting a gene is not evidence it is eight times better
    supported than one with a single warrant; it frequently means the gene sits
    on more commercial test panels. D12 applies: the atlas publishes no validity
    call of its own, and a rank derived from this count would be exactly that.

    Every key is present on every gene, and `submitters` is `[]` rather than
    absent where ClinGen admitted the gene -- the rule `burden_payload` and
    `gene_concordance` already follow, because an object whose shape varies is a
    trap for a consumer reading a field off one gene and expecting it on the
    next.
    """
    clingen = _admitting_clingen(entry)
    if clingen is not None:
        admitted: dict[str, Json] = {
            "authority": "clingen",
            "classification": clingen.classification.value if clingen.classification else None,
            "disease": clingen.disease,
            "disease_label": clingen.disease_label,
            "panel": clingen.gcep,
            "submitters": [],
        }
    else:
        admitted = {
            "authority": "gencc_agreement",
            "classification": None,
            "disease": None,
            "disease_label": None,
            "panel": None,
            "submitters": list(agreeing_submitters(entry)),
        }

    institutions = asserting_institutions(entry)
    return {
        "admitted_by": admitted,
        "asserted_by": {"count": len(institutions), "institutions": list(institutions)},
    }


def asserting_institutions(entry: GeneValidity) -> tuple[str, ...]:
    """`asserted_by.institutions`: every institution asserting the gene in scope, sorted.

    Its own function because two consumers read it -- `admission_provenance`,
    which publishes it, and the gene page's card, which counts the curators
    beside the one that admitted the gene (D55) -- and the card must count
    exactly the list the bundle publishes.

    One name per institution, and only institutions that actually assert the
    gene -- `_asserts` drops the dissenting rungs, which would otherwise be
    counted as agreement. A ClinGen expert panel is ClinGen; its GenCC
    submissions are the same body and must not count again.
    """
    institutions = {
        record.submitter
        for record in entry.records
        if record.source is ValiditySource.GENCC and record.submitter and _asserts(record)
    }
    if any(
        record.source is ValiditySource.CLINGEN and _asserts(record) for record in entry.records
    ):
        institutions.add("ClinGen")
    return tuple(sorted(institutions))
