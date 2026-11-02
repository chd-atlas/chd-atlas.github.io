# src/chd_atlas/models/gene_literature.py
"""The per-gene literature record: the searches run, what was decided about
every hit, and the findings a human approved.

Named ``gene_literature`` because ``models/literature.py`` already holds the
publications, featured and phenotype registries.

A finding states what a paper's *rare-variant or association evidence* said
about one gene, with the denominators it was measured against. The shape rules
below are the ones a number cannot be read without: a count of carriers means
nothing without the count of people screened and what was counted.
"""

from __future__ import annotations

from datetime import date
from typing import Annotated, Final, Literal

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import (
    CohortId,
    HgncId,
    LiteratureFindingId,
    PhenotypeId,
    Pmid,
)
from chd_atlas.models.assertion import Locator
from chd_atlas.models.review import ReviewStatus, check_review
from chd_atlas.vocab import (
    BurdenComparator,
    DetectionMethod,
    ExclusionReason,
    FindingDesign,
    FindingVariantClass,
    LesionGroup,
    ScreenUnit,
    SyndromicStatus,
    TriageDecision,
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
    cases: int | None = Field(default=None, ge=1)
    controls: int | None = Field(default=None, ge=1)
    p_value: float | None = Field(default=None, gt=0, le=1)
    odds_ratio: float | None = Field(default=None, gt=0)
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
        return self

    @model_validator(mode="after")
    def approval_is_a_named_person(self) -> Finding:
        if (self.approved_by is None) != (self.approved_on is None):
            raise ValueError(
                f"{self.id}: approved_by and approved_on are set together or not at all"
            )
        if self.approved_by is not None and self.approved_by.strip().lower().startswith(
            _AGENT_PREFIX
        ):
            raise ValueError(f"{self.id}: approved_by names an agent; a person must approve")
        check_review(self.id, self.review_status, self.review_note)
        return self


class GeneLiteratureFile(BaseModel):
    """``curation/literature/<SYMBOL>.yaml`` -- one gene's searches, triage and findings."""

    model_config = ConfigDict(extra="forbid")

    gene: HgncId
    searches: list[Search] = Field(min_length=1)
    triage: list[TriageEntry] = Field(default_factory=list)
    findings: list[Finding] = Field(default_factory=list)

    @model_validator(mode="after")
    def ids_are_unique(self) -> GeneLiteratureFile:
        duplicate_findings = duplicates(f.id for f in self.findings)
        if duplicate_findings:
            raise ValueError(f"duplicate finding ids: {duplicate_findings}")
        duplicate_searches = duplicates(s.id for s in self.searches)
        if duplicate_searches:
            raise ValueError(f"duplicate search ids: {duplicate_searches}")
        return self
