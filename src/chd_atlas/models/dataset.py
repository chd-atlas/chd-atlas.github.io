# src/chd_atlas/models/dataset.py
from __future__ import annotations

import re
from typing import Final

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationInfo,
    field_validator,
    model_validator,
)

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import AccessionId, ContrastId, Pmid, TaxonId
from chd_atlas.tables import CELL_TYPE_REGIONS
from chd_atlas.vocab import Archive, DatasetDesign, Technology

_ACCESSION_PATTERN: Final[dict[Archive, re.Pattern[str]]] = {
    Archive.PRIDE: re.compile(r"^PXD"),
    Archive.GEO: re.compile(r"^GSE"),
    Archive.ARRAYEXPRESS: re.compile(r"^E-"),
    Archive.EGA: re.compile(r"^EGA[SD]"),
    # A slug, not an accession: the Cell Browser has no accession scheme.
    Archive.UCSC_CELL_BROWSER: re.compile(r"^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$"),
}

_MASS_SPEC: Final[frozenset[Technology]] = frozenset(
    {Technology.TMT_MS, Technology.LFQ_MS, Technology.PHOSPHOPROTEOMICS}
)
_SEQUENCING: Final[frozenset[Technology]] = frozenset(
    {Technology.BULK_RNASEQ, Technology.SCRNASEQ, Technology.MICROARRAY}
)

# PRIDE is the mass-spectrometry repository; GEO, ArrayExpress and EGA hold
# sequencing and array data. A dataset filed under the wrong one is a curator
# error, not an unusual deposit.
_ARCHIVE_TECHNOLOGIES: Final[dict[Archive, frozenset[Technology]]] = {
    Archive.PRIDE: _MASS_SPEC,
    Archive.GEO: _SEQUENCING,
    Archive.ARRAYEXPRESS: _SEQUENCING,
    Archive.EGA: _SEQUENCING,
    Archive.UCSC_CELL_BROWSER: frozenset({Technology.SCRNASEQ}),
}


class Contrast(BaseModel):
    """A named comparison within a dataset.

    A fold change is uninterpretable without knowing what was compared to what,
    so contrasts are modelled explicitly rather than left as a free-text column.
    """

    model_config = ConfigDict(extra="forbid")

    id: ContrastId
    description: str = Field(min_length=1)
    case_group: str = Field(min_length=1)
    control_group: str = Field(min_length=1)
    statistical_method: str = Field(min_length=1)
    software: str = Field(min_length=1)
    covariates: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def groups_differ(self) -> Contrast:
        if self.case_group.strip().casefold() == self.control_group.strip().casefold():
            raise ValueError(
                f"case_group and control_group must differ, both are '{self.case_group}'"
            )
        return self


class Stage(BaseModel):
    """One of a dataset's own stage tokens, with its position in development.

    `wpc` is null post-natally, which is not a gap: a post-natal stage is
    outside every cardiac morphogenetic phase by definition, and the page says
    so rather than dropping the row.

    `order` is the curated chronological position, 1-based and ascending, and
    is what every renderer sorts on. It is **required and not defaulted**,
    because `wpc` cannot do this job: it is null for every post-natal stage,
    so a dataset's post-natal tokens have no numeric order at all. Sorting
    those tokens as strings is what published `elderly` second of eight and
    `4 week post conception` after `19 week post conception` -- deterministic,
    reproducible and wrong, under schema 2.11.

    **2.11 and no other version**, measured rather than estimated: the
    profiles mirror was first committed on 2026-08-19, and before it a gene
    bundle's `expression_profile` was `{"datasets": []}` with no stage array
    published at all. So the dictionary order reached published bytes under
    exactly one schema version, promoted to `main` that same day and corrected
    the next.

    Chosen over sorting by position in the YAML file, which would make
    declaration order silently load-bearing: a curator tidying the block would
    change published output with nothing flagging it.
    """

    model_config = ConfigDict(extra="forbid")

    token: str = Field(min_length=1)
    wpc: float | None = Field(default=None, gt=0)
    order: int = Field(ge=1)


class CellTypeDecl(BaseModel):
    """One of a cell-type dataset's own cell-type tokens.

    `token` is the string the source's matrix uses; `label` is what a page
    shows; `compartment` groups tokens for display. All three are required
    and must not be blank, because a blank one renders as an unlabelled row.
    """

    model_config = ConfigDict(extra="forbid")

    token: str = Field(min_length=1)
    label: str = Field(min_length=1)
    compartment: str = Field(min_length=1)

    # `min_length=1` alone accepted whitespace, `'   '` (review, 2026-10-28).
    @field_validator("token", "label", "compartment")
    @classmethod
    def is_not_blank(cls, value: str, info: ValidationInfo) -> str:
        if not value.strip():
            raise ValueError(f"{info.field_name} must not be blank")
        return value


class Dataset(BaseModel):
    """Top level of ``curation/datasets/<ACCESSION>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    id: AccessionId
    archive: Archive
    technology: Technology
    tissue: str = Field(min_length=1)
    developmental_stage: str = Field(min_length=1)
    organism: TaxonId
    # A dataset with contrasts declares at least a case and a control sample.
    # The same `ge=2` also binds a profile dataset, where the rationale is
    # unstated (see D40/D41 for what a profile dataset must declare instead).
    n_samples: int = Field(ge=2)
    licence: str = Field(min_length=1)
    contrasts: list[Contrast] = Field(default_factory=list)
    design: DatasetDesign
    # D40. A set, not a scalar: a dataset with two cardiac tokens would file the
    # second as a comparison organ, and tau would then measure specificity
    # against the heart's own other chamber -- driving it down for exactly the
    # most cardiac genes.
    cardiac_tissues: tuple[str, ...] = ()
    stages: tuple[Stage, ...] = ()
    # Only a `cell_type_profile` dataset declares these; see
    # `a_cell_type_profile_dataset_is_fully_declared`.
    cell_types: tuple[CellTypeDecl, ...] = ()
    regions: tuple[str, ...] = ()
    # D41. Per dataset, never a constant: a floor of 1.0 is a defensible RPKM
    # convention, is meaningless in `nx`, and is actively wrong in `lfq`, where
    # intensities are log-scale and 0 is a middling abundance.
    # `gt=0` because a floor of zero is not a floor. Without it a curated
    # `0.0` places a median of `0.0` -- the gate is `median < floor`, and
    # `0.0 < 0.0` is False -- and that cell then reaches `pages._Median` with
    # a placement AND no gap reason, which is the one pairing that type's
    # docstring says cannot happen. Measured: 1,309 of the 18,326 rows in
    # `mirrors/profiles/E-MTAB-6814.tsv` carry a median of exactly 0.0, so the
    # only thing that stood between the corpus and a page contradicting its
    # own table was the curated value being 1.0.
    detection_floor: float | None = Field(default=None, gt=0)
    floor_source: str | None = None
    quantile_estimator: str | None = None
    publication: Pmid | None = None

    @field_validator("licence")
    @classmethod
    def licence_is_not_blank(cls, value: str) -> str:
        """Same rule, same reason, as `Source.licence`.

        `min_length=1` accepts "   ", which records no provenance at all. A
        dataset redistributes third-party measurements exactly as a source does,
        so the two must not disagree about what counts as a recorded licence.
        """
        if not value.strip():
            raise ValueError("licence must not be blank")
        return value

    @field_validator("cardiac_tissues")
    @classmethod
    def cardiac_tissues_are_not_blank(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        """Same rule as `licence_is_not_blank`: a blank entry names no tissue at all."""
        if any(not tissue.strip() for tissue in value):
            raise ValueError("cardiac_tissues must not contain a blank entry")
        return value

    @model_validator(mode="after")
    def accession_matches_archive(self) -> Dataset:
        if not _ACCESSION_PATTERN[self.archive].match(self.id):
            raise ValueError(f"accession {self.id} does not match archive '{self.archive.value}'")
        return self

    @model_validator(mode="after")
    def contrast_ids_are_unique(self) -> Dataset:
        found = duplicates(contrast.id for contrast in self.contrasts)
        if found:
            raise ValueError(f"duplicate contrast ids: {found}")
        return self

    @model_validator(mode="after")
    def technology_matches_archive(self) -> Dataset:
        if self.technology not in _ARCHIVE_TECHNOLOGIES[self.archive]:
            raise ValueError(
                f"technology '{self.technology.value}' is not deposited in "
                f"archive '{self.archive.value}'"
            )
        return self

    @model_validator(mode="after")
    def contrasts_match_design(self) -> Dataset:
        if self.design is DatasetDesign.CONTRAST and not self.contrasts:
            raise ValueError("design 'contrast' requires at least one contrast")
        if self.design is DatasetDesign.PROFILE and self.contrasts:
            raise ValueError("design 'profile' declares no contrasts")
        if self.design is DatasetDesign.CELL_TYPE_PROFILE and self.contrasts:
            raise ValueError("design 'cell_type_profile' declares no contrasts")
        return self

    @model_validator(mode="after")
    def cell_types_and_regions_belong_to_their_design(self) -> Dataset:
        if self.design is not DatasetDesign.CELL_TYPE_PROFILE and (self.cell_types or self.regions):
            raise ValueError(f"design '{self.design.value}' declares no cell_types or regions")
        return self

    @model_validator(mode="after")
    def a_profile_dataset_is_fully_declared(self) -> Dataset:
        """Everything the derivation needs, required where it is used.

        Left optional on a contrast dataset because none of it applies there;
        required here because each field is the input to a published figure,
        and a null would publish the figure computed against a default nobody
        curated.
        """
        if self.design is not DatasetDesign.PROFILE:
            return self
        if not self.cardiac_tissues:
            raise ValueError("cardiac_tissues must not be empty for a profile dataset")
        if not self.stages:
            raise ValueError("stages must not be empty for a profile dataset")
        if self.detection_floor is None or not self.floor_source:
            raise ValueError("a profile dataset declares detection_floor and floor_source")
        if not self.quantile_estimator:
            raise ValueError("a profile dataset declares quantile_estimator")
        return self

    @model_validator(mode="after")
    def a_cell_type_profile_dataset_is_fully_declared(self) -> Dataset:
        """The cell-type counterpart of `a_profile_dataset_is_fully_declared`.

        Needs no `stages` or `cardiac_tissues`: it is indexed by cell type and
        region, so those two are what it declares instead. The floor and the
        estimator are required for the reason stated there.
        """
        if self.design is not DatasetDesign.CELL_TYPE_PROFILE:
            return self
        if not self.cell_types:
            raise ValueError("cell_types must not be empty for a cell_type_profile dataset")
        if not self.regions:
            raise ValueError("regions must not be empty for a cell_type_profile dataset")
        if self.detection_floor is None or not self.floor_source:
            raise ValueError(
                "a cell_type_profile dataset declares detection_floor and floor_source"
            )
        if not self.quantile_estimator:
            raise ValueError("a cell_type_profile dataset declares quantile_estimator")
        return self

    @model_validator(mode="after")
    def cell_type_tokens_are_unique(self) -> Dataset:
        found = duplicates(cell_type.token for cell_type in self.cell_types)
        if found:
            raise ValueError(f"duplicate cell_type tokens: {found}")
        return self

    @field_validator("regions")
    @classmethod
    def regions_are_the_mirror_s_own(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        """Unique, and each a region the mirror can carry, `all` excepted.

        The page draws one column per declared region. A region the mirror's
        `region` column cannot hold (`tables.CELL_TYPE_REGIONS`) would draw a
        column of dashes captioned "no heart had 10 or more cells there", which
        no heart could have; `all` is the pooled row, not a region. The model
        accepted `['LV', 'LV', 'all', 'Apex', ' ']` until review (2026-10-28).
        """
        found = duplicates(value)
        if found:
            raise ValueError(f"duplicate regions: {found}")
        unknown = [region for region in value if region not in CELL_TYPE_REGIONS]
        if unknown:
            raise ValueError(
                f"regions {unknown} are not among {sorted(CELL_TYPE_REGIONS)}; the pooled "
                "row 'all' is not declared"
            )
        return value

    @model_validator(mode="after")
    def stage_tokens_are_unique(self) -> Dataset:
        found = duplicates(stage.token for stage in self.stages)
        if found:
            raise ValueError(f"duplicate stage tokens: {found}")
        return self

    @model_validator(mode="after")
    def cardiac_tissues_are_unique(self) -> Dataset:
        found = duplicates(self.cardiac_tissues)
        if found:
            raise ValueError(f"duplicate cardiac_tissues: {found}")
        return self
