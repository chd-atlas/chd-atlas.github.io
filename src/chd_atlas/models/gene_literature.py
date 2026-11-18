# src/chd_atlas/models/gene_literature.py
"""The per-gene literature record: the searches run, what was decided about
every hit, and the findings extracted from the included papers, each approved
by a human or not yet.

Named ``gene_literature`` because ``models/literature.py`` already holds the
publications, featured and phenotype registries.

A finding states what a paper's *rare-variant or association evidence* said
about one gene, with the denominators it was measured against. The shape rules
below are the ones a number cannot be read without: a count of carriers means
nothing without the count of people screened and what was counted.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Sequence
from datetime import date
from typing import Annotated, Any, Final, Literal

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import (
    AssertionId,
    CohortId,
    HgncId,
    LiteratureFindingId,
    PhenotypeId,
    Pmid,
)
from chd_atlas.models.assertion import LesionAssertion, Locator
from chd_atlas.models.review import ReviewStatus, check_review
from chd_atlas.vocab import (
    BurdenComparator,
    DetectionMethod,
    ExclusionReason,
    FindingDesign,
    FindingVariantClass,
    LesionGroup,
    Phase,
    ScreenUnit,
    SyndromicStatus,
    TriageDecision,
    VariantPlacement,
    VariantZygosity,
)

_NonBlankStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]

# Designs whose evidence is a count of carriers among people screened.
_RARE_VARIANT_DESIGNS: Final[frozenset[FindingDesign]] = frozenset(
    {
        FindingDesign.CASE_REPORT,
        FindingDesign.FAMILY_SEGREGATION,
        FindingDesign.COHORT_SCREEN,
        FindingDesign.RARE_VARIANT_ASSOCIATION,
    }
)

# Methods that genotype known sites; they cannot find a rare variant.
_GENOTYPING_METHODS: Final[frozenset[DetectionMethod]] = frozenset(
    {DetectionMethod.GENOTYPING_ARRAY, DetectionMethod.TARGETED_GENOTYPING}
)

# A person approves; an agent only proposes. Compared after strip and lower().
_AGENT_PREFIX: Final = "agent:"


# Fields outside the digest: the approval itself, and marks that do not change
# what the finding says (a review mark does not withdraw an approval; a verifier
# mark is bookkeeping).
_DIGEST_EXCLUDED: Final[frozenset[str]] = frozenset(
    {
        "approved_by",
        "approved_on",
        "approved_digest",
        "verified_by",
        "review_status",
        "review_note",
    }
)

# A variant's review mark is not an edit to what the finding says, as the
# finding's own is not; the synthesis covers it (`_covered_finding`).
_VARIANT_DIGEST_EXCLUDED: Final[frozenset[str]] = frozenset({"review_status", "review_note"})


class Search(BaseModel):
    """One recorded PubMed search and its arithmetic."""

    model_config = ConfigDict(extra="forbid")

    id: int = Field(ge=1)
    source: Literal["pubmed"]
    query: _NonBlankStr
    run_on: date
    hits: int = Field(ge=0)
    new_hits: int = Field(ge=0)

    @model_validator(mode="after")
    def new_hits_are_among_the_hits(self) -> Search:
        if self.new_hits > self.hits:
            raise ValueError(
                f"search {self.id}: new_hits {self.new_hits} exceeds hits {self.hits}; "
                f"a search cannot find more new papers than papers"
            )
        return self


class TriageEntry(BaseModel):
    """What was decided about one PubMed hit, attributed to the search that first found it."""

    model_config = ConfigDict(extra="forbid")

    pmid: Pmid
    search: int = Field(ge=1)
    decision: TriageDecision
    reason: ExclusionReason | None = None
    same_as: Pmid | None = None
    note: _NonBlankStr | None = None

    @model_validator(mode="after")
    def reason_matches_decision(self) -> TriageEntry:
        if self.decision is TriageDecision.EXCLUDED and self.reason is None:
            raise ValueError(f"{self.pmid}: an excluded paper needs a reason")
        if self.decision is not TriageDecision.EXCLUDED and self.reason is not None:
            raise ValueError(
                f"{self.pmid}: reason is only for an excluded paper, not {self.decision.value}"
            )
        if self.reason is ExclusionReason.SAME_DATA_AS and self.same_as is None:
            raise ValueError(f"{self.pmid}: same_data_as needs same_as naming the other paper")
        if self.reason is not ExclusionReason.SAME_DATA_AS and self.same_as is not None:
            raise ValueError(f"{self.pmid}: same_as is only for the reason same_data_as")
        if self.reason is ExclusionReason.OTHER and self.note is None:
            raise ValueError(f"{self.pmid}: the reason 'other' needs a note")
        if self.reason is ExclusionReason.FULL_TEXT_NOT_ACCESSIBLE and self.note is None:
            raise ValueError(
                f"{self.pmid}: the reason 'full_text_not_accessible' needs a note naming "
                "what the abstract lacks"
            )
        return self


class Carriers(BaseModel):
    """Carriers of a reported variant, split by what is known of its origin."""

    model_config = ConfigDict(extra="forbid")

    de_novo: int = Field(ge=0)
    inherited: int = Field(ge=0)
    origin_unknown: int = Field(ge=0)

    @property
    def total(self) -> int:
        return self.de_novo + self.inherited + self.origin_unknown


class ReportedNames(BaseModel):
    """The names a paper itself gave a variant, as printed, none of them checked."""

    model_config = ConfigDict(extra="forbid")

    transcript: _NonBlankStr | None = None
    c: _NonBlankStr | None = None
    p: _NonBlankStr | None = None
    rsid: _NonBlankStr | None = None


class Normalisation(BaseModel):
    """One run of the normaliser, with the versions it ran against.

    Variants point at it by `id`, so the versions are stated once per run.
    """

    model_config = ConfigDict(extra="forbid")

    id: int = Field(ge=1)
    tool: Literal["variantvalidator"]
    version: _NonBlankStr
    hgvs_version: _NonBlankStr
    vvta: _NonBlankStr
    vvdb: _NonBlankStr
    seqrepo: _NonBlankStr
    run_on: date


class Variant(BaseModel):
    """One variant a finding names, placed on the MANE Select frame.

    The YAML key `with` is a Python keyword, so the field is `with_` and is
    read and dumped through its alias. Only the review pairing is checked
    here; every other rule is a validate code, so one bad variant costs one
    issue and not the file.
    """

    model_config = ConfigDict(extra="forbid", populate_by_name=True)

    id: int = Field(ge=1)
    as_reported: _NonBlankStr
    locator: Locator
    kind: FindingVariantClass
    placement: VariantPlacement
    reported: ReportedNames = ReportedNames()
    hgvs_c: _NonBlankStr | None = None
    hgvs_p: _NonBlankStr | None = None
    grch38: _NonBlankStr | None = None
    normalisation: int | None = None
    zygosity: VariantZygosity = VariantZygosity.NOT_STATED
    with_: list[int] = Field(default_factory=list, alias="with")
    phase: Phase | None = None
    carriers: Carriers | None = None
    phenotypes: list[PhenotypeId] = Field(default_factory=list)
    review_status: ReviewStatus | None = None
    review_note: _NonBlankStr | None = None

    @model_validator(mode="after")
    def review_pairing(self) -> Variant:
        check_review(f"variant {self.id}", self.review_status, self.review_note)
        return self


class Finding(BaseModel):
    """What one paper reports about one gene, approved by a named person."""

    model_config = ConfigDict(extra="forbid")

    id: LiteratureFindingId
    publication: Pmid
    design: FindingDesign
    methods: list[DetectionMethod] = Field(min_length=1)
    comparator: BurdenComparator | None = None
    burden_study: Pmid | None = None
    count_unit: ScreenUnit | None = None
    screened: int | None = Field(default=None, ge=1)
    carriers: Carriers | None = None
    variant_classes: list[FindingVariantClass] = Field(default_factory=list)
    variants: list[Variant] = Field(default_factory=list)
    cases: int | None = Field(default=None, ge=1)
    controls: int | None = Field(default=None, ge=1)
    p_value: float | None = Field(default=None, gt=0, le=1, allow_inf_nan=False)
    odds_ratio: float | None = Field(default=None, gt=0, allow_inf_nan=False)
    phenotypes: list[PhenotypeId] = Field(min_length=1)
    lesion_groups: list[LesionGroup] = Field(min_length=1)
    syndromic: SyndromicStatus
    cohorts: list[CohortId] = Field(default_factory=list)
    shares_probands_with: list[LiteratureFindingId] | Literal["unknown"]
    locators: list[Locator] = Field(min_length=1)
    summary: _NonBlankStr
    extracted_by: _NonBlankStr
    verified_by: _NonBlankStr | None = None
    approved_by: _NonBlankStr | None = None
    approved_on: date | None = None
    approved_digest: _NonBlankStr | None = None
    review_status: ReviewStatus | None = None
    review_note: _NonBlankStr | None = None

    @model_validator(mode="after")
    def design_rules(self) -> Finding:
        design = self.design
        if design in _RARE_VARIANT_DESIGNS:
            for name in ("count_unit", "screened", "carriers"):
                if getattr(self, name) is None:
                    raise ValueError(f"{self.id}: {design.value} needs {name}")
            if not self.variant_classes:
                raise ValueError(f"{self.id}: {design.value} needs at least one variant_classes")
            if (
                self.carriers is not None
                and self.screened is not None
                and self.carriers.total > self.screened
            ):
                raise ValueError(
                    f"{self.id}: carriers sum to {self.carriers.total}, more than "
                    f"screened {self.screened}"
                )
        if design is FindingDesign.CASE_REPORT and self.screened != 1:
            raise ValueError(f"{self.id}: a case_report screens exactly one; screened is not 1")
        if design is FindingDesign.RARE_VARIANT_ASSOCIATION:
            if self.comparator is None:
                raise ValueError(f"{self.id}: rare_variant_association needs a comparator")
            if self.comparator is BurdenComparator.NONE:
                raise ValueError(
                    f"{self.id}: comparator 'none' is a case series; "
                    f"that is a cohort_screen, not a rare_variant_association"
                )
        elif self.comparator is not None:
            raise ValueError(f"{self.id}: comparator is only for rare_variant_association")
        if design is FindingDesign.COMMON_VARIANT_ASSOCIATION:
            if self.cases is None or self.controls is None:
                raise ValueError(f"{self.id}: common_variant_association needs cases and controls")
            if self.carriers is not None or self.screened is not None:
                raise ValueError(
                    f"{self.id}: common_variant_association counts cases and controls, "
                    f"not carriers or screened"
                )
        elif _GENOTYPING_METHODS.intersection(self.methods):
            raise ValueError(
                f"{self.id}: genotyping_array and targeted_genotyping appear only with "
                f"common_variant_association"
            )
        if DetectionMethod.NOT_STATED in self.methods and len(self.methods) > 1:
            raise ValueError(f"{self.id}: not_stated stands alone; the paper names no method")
        return self

    @model_validator(mode="after")
    def approval_is_a_named_person(self) -> Finding:
        approval = (self.approved_by, self.approved_on, self.approved_digest)
        if any(part is None for part in approval) and any(part is not None for part in approval):
            raise ValueError(
                f"{self.id}: approved_by, approved_on and approved_digest are set together "
                f"or not at all"
            )
        if self.approved_by is not None and self.approved_by.strip().lower().startswith(
            _AGENT_PREFIX
        ):
            raise ValueError(f"{self.id}: approved_by names an agent; a person must approve")
        check_review(self.id, self.review_status, self.review_note)
        return self


def finding_digest(finding: Finding) -> str:
    """`sha256:` and the hex digest of the finding's content, approval excluded.

    The canonical JSON (sorted keys, no ASCII escaping, compact separators) of
    `model_dump(mode="json", exclude_defaults=True, by_alias=True)` without
    `_DIGEST_EXCLUDED`, and without each variant's review mark
    (`_VARIANT_DIGEST_EXCLUDED`). `by_alias` keys a variant's links `with`, as
    in the YAML. The file's `normalisations` are not on the finding and enter
    no digest.
    An approval records this value, so a later edit to anything it covers no
    longer matches (LIT015). A field left at its default is not hashed, so a
    release that adds an optional field does not void standing approvals; the
    cost is that setting a field explicitly to its default value is not an edit.
    """
    content = {
        k: v
        for k, v in finding.model_dump(mode="json", exclude_defaults=True, by_alias=True).items()
        if k not in _DIGEST_EXCLUDED
    }
    if "variants" in content:
        content["variants"] = [
            {k: v for k, v in variant.items() if k not in _VARIANT_DIGEST_EXCLUDED}
            for variant in content["variants"]
        ]
    canonical = json.dumps(content, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()


class SynthesisSentence(BaseModel):
    """One sentence of a gene's synthesis and the records it rests on."""

    model_config = ConfigDict(extra="forbid")

    text: _NonBlankStr
    cites: list[LiteratureFindingId | AssertionId] = Field(min_length=1)

    @model_validator(mode="after")
    def cites_are_distinct(self) -> SynthesisSentence:
        repeated = duplicates(self.cites)
        if repeated:
            raise ValueError(f"a sentence cites the same record twice: {repeated}")
        return self


class Synthesis(BaseModel):
    """A gene's "In short" narrative: one to five sentences, each citing records.

    Drafted by an agent and approved by a named person; the approval is bound to
    `synthesis_digest`, so it is void after a change to the text or the cites, to
    `drafted_by`, to any approved finding of the file (its content, review mark or
    approval, so re-approving one too), to which findings are approved, or to any
    lesion assertion of the gene.
    """

    model_config = ConfigDict(extra="forbid")

    sentences: list[SynthesisSentence] = Field(min_length=1, max_length=5)
    drafted_by: _NonBlankStr
    approved_by: _NonBlankStr | None = None
    approved_on: date | None = None
    approved_digest: _NonBlankStr | None = None

    @model_validator(mode="after")
    def approval_is_a_named_person(self) -> Synthesis:
        approval = (self.approved_by, self.approved_on, self.approved_digest)
        if any(part is None for part in approval) and any(part is not None for part in approval):
            raise ValueError(
                "synthesis: approved_by, approved_on and approved_digest are set together "
                "or not at all"
            )
        if self.approved_by is not None and self.approved_by.strip().lower().startswith(
            _AGENT_PREFIX
        ):
            raise ValueError("synthesis: approved_by names an agent; a person must approve")
        if not self.drafted_by.strip().lower().startswith(_AGENT_PREFIX):
            raise ValueError(
                "synthesis: drafted_by must name the drafting agent (agent:<skill>/<model>); "
                "a person approves, an agent drafts"
            )
        return self


def _canonical_digest(content: object) -> str:
    canonical = json.dumps(content, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _covered_finding(finding: Finding) -> str:
    """A finding's digest plus what `finding_digest` leaves out and a synthesis rests on.

    The review mark, and the approval itself (`approved_by`, `approved_on` as an
    ISO date, `approved_digest`): a re-approved finding, even with unchanged
    content, withholds the synthesis until it is approved again (GS5).
    `verified_by` stays out, as it does of the finding's own digest.
    """
    covered: dict[str, Any] = {
        "finding": finding_digest(finding),
        "review_status": finding.review_status.value if finding.review_status else None,
        "review_note": finding.review_note or None,
        "approved_by": finding.approved_by,
        "approved_on": finding.approved_on.isoformat() if finding.approved_on else None,
        "approved_digest": finding.approved_digest,
    }
    if finding.variants:
        # Only when there are variants, so a finding with none keeps its digest.
        covered["variant_reviews"] = [
            [v.id, v.review_status.value if v.review_status else None, v.review_note or None]
            for v in sorted(finding.variants, key=lambda v: v.id)
        ]
    return _canonical_digest(covered)


def synthesis_digest(
    synthesis: Synthesis,
    findings: Sequence[Finding],
    assertions: Sequence[LesionAssertion],
) -> str:
    """`sha256:` and the hex digest of what a synthesis says and rests on.

    Canonical JSON (as `finding_digest`) of the sentences' `text` and `cites`
    (in order), `drafted_by` (the card prints it as the drafter), the sorted
    `_covered_finding` digest of every finding in `findings` with `approved_by`
    set (its `finding_digest`, its review mark and its approval), and the sorted
    digests of `assertions`, each the canonical JSON of
    `model_dump(mode="json", exclude_defaults=True)`. Only the synthesis's own
    approval fields are not covered. Findings are filtered here, so the caller
    passes all of the file's; assertions are NOT filtered by gene, so the caller
    passes the gene's own lesion assertions. Input order does not matter.
    Editing the text, the cites or `drafted_by`, editing an approved finding,
    marking it for review, re-approving it, approving another finding, or
    editing, adding or removing an assertion voids a standing approval.
    """
    content = {
        "sentences": [{"text": s.text, "cites": list(s.cites)} for s in synthesis.sentences],
        "drafted_by": synthesis.drafted_by,
        "findings": sorted(_covered_finding(f) for f in findings if f.approved_by is not None),
        "assertions": sorted(
            _canonical_digest(a.model_dump(mode="json", exclude_defaults=True)) for a in assertions
        ),
    }
    return _canonical_digest(content)


class GeneLiteratureFile(BaseModel):
    """``curation/literature/<SYMBOL>.yaml`` -- one gene's searches, triage and findings."""

    model_config = ConfigDict(extra="forbid")

    gene: HgncId
    searches: list[Search] = Field(min_length=1)
    normalisations: list[Normalisation] = Field(default_factory=list)
    triage: list[TriageEntry] = Field(default_factory=list)
    findings: list[Finding] = Field(default_factory=list)
    synthesis: Synthesis | None = None

    @model_validator(mode="after")
    def ids_are_unique(self) -> GeneLiteratureFile:
        duplicate_findings = duplicates(f.id for f in self.findings)
        if duplicate_findings:
            raise ValueError(f"duplicate finding ids: {duplicate_findings}")
        duplicate_searches = duplicates(s.id for s in self.searches)
        if duplicate_searches:
            raise ValueError(f"duplicate search ids: {duplicate_searches}")
        duplicate_normalisations = duplicates(n.id for n in self.normalisations)
        if duplicate_normalisations:
            raise ValueError(f"duplicate normalisation ids: {duplicate_normalisations}")
        return self
