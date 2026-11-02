# src/chd_atlas/models/phases.py
"""The cardiac morphogenetic phases a developmental stage can fall in.

This is the first developmental-biology claim in the atlas, and it is
transcribed rather than authored: every phase carries the source that defines
its boundaries. `chd_scope.yaml` is the precedent (`cohorts.yaml` is not --
measured 2026-08-14, it carries no attribution at all).

Intervals are half-open, `[start_wpc, end_wpc)`, so a stage sitting exactly on
a boundary falls in exactly one phase's own edge case rather than a coin flip.

**A stage may carry SEVERAL phases at once, and this is a fact about the
biology, not a modelling looseness.** An earlier version of this file enforced
a disjoint partition -- one phase per stage -- and that was wrong. Human
cardiac morphogenetic processes run concurrently; a disjoint partition would
force picking one label per stage and silently discarding the others, and the
discarded one may be exactly the lesion-relevant process a reader came for.

**This table was corrected once and the correction is recorded here on
purpose.** It originally read "4 | looping, atrial septation, ventricular
septation" and so on, transcribed from the source's *ordinal*-week prose --
which is exactly the off-by-one-week error this file's own next section warns
against, and it produced a wrong row (it placed outflow tract septation at 7
wpc; the numeric boundaries below place it only through 6.86 wpc). The table
below is derived from the boundaries `curation/cardiac_phases.yaml` actually
curates, verified against a real built page rather than trusted by
arithmetic. `heart_looping` is absent from every row because its end, 3.71
wpc (`EndBasis.DERIVED` below), falls before the first of these stages. This
paragraph said `phases_for` excluded it for having no stated end, which
stopped being true on 2026-08-24, when that end was derived:

| elapsed wpc | phases matched |
|---|---|
| 4 | atrial septation, ventricular septation |
| 6 | atrial septation, ventricular septation, outflow tract septation |
| 7 | ventricular septation, valve formation |

At 6 wpc an AVSD gene needs "atrial septation" and a TGA gene needs "outflow
tract septation", and both are true at once -- a disjoint model can represent
only one of them. So `CardiacPhaseFile` enforces unique phase **ids** only; it
does not, and must not, enforce that phases are pairwise disjoint.

**Boundaries are computed from HsapDv day-post-fertilization ranges (day /
7.0), never transcribed from a source's ordinal-week prose.** A source stating
"the sixth week of development" means elapsed wpc **5 to 6** (ordinal week N
spans elapsed (N-1) to N), not 6 to 7 -- HsapDv gives Carnegie stage 16 as day
37-42 (5.29-6.00 elapsed wpc), and Buijtendijk et al. 2020 calls CS16 "the
sixth week", which is consistent with the ordinal reading and NOT with reading
"sixth week" as elapsed wpc 6-7. Transcribing the prose directly would shift
every boundary a full Carnegie stage too late, silently and consistently
enough that nothing downstream would catch it -- exactly the kind of error
this file's own docstring exists to name rather than let recur, and exactly
the error the table above was corrected from. See `curation/cardiac_
phases.yaml`'s header comment for the same warning at the point a curator is
most likely to re-introduce it.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from chd_atlas.duplicates import duplicates


class EndBasis(StrEnum):
    """Where a phase's `end_wpc` comes from: the source, an ontology, or nowhere.

    **A retired member first, because the third one below has to be read
    against it.** An earlier revision had `CURATOR_CAPPED`, for a phase the
    source names no end for -- "capped" at the last Carnegie stage the source
    discusses at all, on the reasoning that the heart "already resembles the
    postnatal configuration" by then. That statement turned out to describe
    arterial-pole remodelling, not any of the three processes it was used to
    cap, and was this atlas inventing a boundary the source never gave -- once
    for a genuine gap (`heart_looping`) and twice for phases whose real
    completion had simply not been extracted yet
    (`ventricular_septum_morphogenesis`, `heart_valve_morphogenesis` -- both
    now `STATED`, at the source's own words). It is retired rather than kept
    unused: a vocabulary member no curated phase exercises is a guard nobody
    has seen fail.

    - `STATED` -- the source names this boundary directly. Every end field
      (`end_wpc`, `end_carnegie_stage`, `end_hsapdv_id`) is populated.
    - `DERIVED` -- the boundary is not in the source's prose, but follows from
      the phase's **own ontology term** plus a boundary this file already
      carries as `STATED`. Every end field is populated, exactly as for
      `STATED`, and the two are kept apart so a reader can tell which
      sentences came from the cited paper and which this atlas concluded.
    - `NOT_STATED` -- nothing gives an end. Every end field is `None`.

    `CardiacPhase.end_fields_match_end_basis` enforces the pairing in both
    directions, so the basis and the fields can never disagree.

    **`DERIVED` is not `CURATOR_CAPPED` under a new name, and the difference is
    the whole reason it is allowed.** `CURATOR_CAPPED` took a statement about a
    *different* process and applied it to this one; there was no relation
    between the sentence and the boundary beyond both appearing in the same
    paper. `DERIVED` requires an actual entailment, and the live case is the
    only kind that qualifies: GO:0001947 -- the term `heart_looping` already
    carries in its own `go_id` -- *defines* the process as ending "preceding
    septation", and this file already dates the start of septation at CS12
    from the source, `end_basis: stated`. The end is not a guess at what the
    paper might have meant; it is what the phase's own definition says, read
    against a boundary already curated here.

    A `DERIVED` end is still weaker than a `STATED` one and must stay
    visible as such: it inherits every uncertainty of the boundary it is
    anchored to, and it would move if that anchor moved. A future phase whose
    end follows from nothing at all is `NOT_STATED`, with `end_wpc` left null
    rather than a boundary invented to fill it -- that is what this vocabulary
    exists to make expressible.
    """

    STATED = "stated"
    DERIVED = "derived"
    NOT_STATED = "not_stated"


class CardiacPhase(BaseModel):
    """One named developmental process and the window it occupies.

    `id` is this atlas's own slug, rendered on a gene page by replacing
    underscores with spaces (`build/pages.py::_phases_cell`) -- distinct
    from `go_id`, which is the transcribed Gene Ontology accession and is
    never rendered directly. Every boundary field below exists so a reader can
    check this atlas's transcription against the source without taking this
    file's word for it: `start_carnegie_stage`/`end_carnegie_stage` name the
    Carnegie stages the source itself cites, `start_hsapdv_id`/`end_hsapdv_id`
    name the HsapDv ontology terms those stages resolve to.

    **The start is always known; the end may not be.** Every phase this
    atlas curates has a source-stated start, so `start_wpc`/
    `start_carnegie_stage`/`start_hsapdv_id` are always required. The three
    end fields are `None` together exactly when `end_basis` is `NOT_STATED`
    (`end_fields_match_end_basis` enforces this) -- a phase with no stated
    end is recorded with its real start and an honest gap where its end
    would be, never a curator's guess at one.
    """

    model_config = ConfigDict(extra="forbid")

    id: str = Field(min_length=1)
    go_id: str = Field(min_length=1)
    label: str = Field(min_length=1)
    start_wpc: float = Field(gt=0)
    end_wpc: float | None = Field(default=None, gt=0)
    start_carnegie_stage: str = Field(min_length=1)
    end_carnegie_stage: str | None = Field(default=None, min_length=1)
    start_hsapdv_id: str = Field(min_length=1)
    end_hsapdv_id: str | None = Field(default=None, min_length=1)
    end_basis: EndBasis

    @model_validator(mode="after")
    def interval_is_non_empty(self) -> CardiacPhase:
        """`end_wpc`, when present, must exceed `start_wpc`.

        Skipped entirely when `end_wpc is None` -- there is no interval to
        check, and `end_fields_match_end_basis` below is what enforces that a
        null `end_wpc` is itself a deliberate, `NOT_STATED` fact rather than
        an oversight.
        """
        if self.end_wpc is not None and self.end_wpc <= self.start_wpc:
            raise ValueError(f"phase {self.id}: end_wpc must exceed start_wpc")
        return self

    @model_validator(mode="after")
    def end_fields_match_end_basis(self) -> CardiacPhase:
        """`end_basis` and the three end fields must never disagree.

        `STATED` and `DERIVED` require all three (`end_wpc`,
        `end_carnegie_stage`, `end_hsapdv_id`) to be present; `NOT_STATED`
        requires all three to be absent. Without this check a curator could
        write `end_basis: stated` with `end_wpc: null` (a claimed boundary
        with no number behind it) or `end_basis: not_stated` with a real
        `end_wpc` (a value this atlas would then treat as unsourced everywhere
        it matters -- `phases_for` and PRF006's coverage-span merge both key
        on `end_wpc is None`, not on `end_basis`, so a mismatch here would
        silently reintroduce exactly the invented-boundary defect
        `NOT_STATED` exists to rule out).

        **The "has an end" test is written as `is not NOT_STATED`, not as a
        list of the two members that do.** A fourth member added later
        inherits the requirement by default, so the failure mode of forgetting
        to extend this check is a phase that *cannot* be curated rather than
        one that publishes a boundary with nothing behind it. Only
        `NOT_STATED` may leave the fields empty, and that is the property
        worth stating positively.
        """
        end_fields = (self.end_wpc, self.end_carnegie_stage, self.end_hsapdv_id)
        if self.end_basis is not EndBasis.NOT_STATED and any(field is None for field in end_fields):
            raise ValueError(
                f"phase {self.id}: end_basis is {self.end_basis.value} but an end field is missing"
            )
        if self.end_basis is EndBasis.NOT_STATED and any(field is not None for field in end_fields):
            raise ValueError(f"phase {self.id}: end_basis is not_stated but an end field is set")
        return self


class CardiacPhaseFile(BaseModel):
    """Top level of ``curation/cardiac_phases.yaml``."""

    model_config = ConfigDict(extra="forbid")

    attributed_to: str = Field(min_length=1)
    citation: str = Field(min_length=1)
    phases: list[CardiacPhase] = Field(default_factory=list)

    @field_validator("attributed_to")
    @classmethod
    def attributed_to_is_not_blank(cls, value: str) -> str:
        """Same rule and reason as `Dataset.licence_is_not_blank`.

        `min_length=1` accepts "   ", which records no provenance at all --
        and provenance is the entire justification for a developmental-biology
        claim reaching this atlas. The atlas transcribes embryology; it does
        not author it, so a blank here is worse than an absent file: it looks
        attributed and is not.
        """
        if not value.strip():
            raise ValueError("attributed_to must not be blank")
        return value

    @field_validator("citation")
    @classmethod
    def citation_is_not_blank(cls, value: str) -> str:
        """Same rule as `attributed_to_is_not_blank`: a blank citation names no source."""
        if not value.strip():
            raise ValueError("citation must not be blank")
        return value

    @model_validator(mode="after")
    def phase_ids_are_unique(self) -> CardiacPhaseFile:
        """Ids must be unique; **intervals need not be disjoint** (see the module docstring).

        An earlier version of this validator also rejected overlapping
        intervals. That was a modelling error, not a data-quality check: human
        cardiac morphogenesis runs several processes at once, so a stage
        legitimately belongs to more than one phase, and rejecting that at
        load time would make it impossible to curate the real vocabulary at
        all rather than catch a curator's mistake.
        """
        found = duplicates(phase.id for phase in self.phases)
        if found:
            raise ValueError(f"duplicate phase ids: {found}")
        return self

    def phases_for(self, wpc: float | None) -> tuple[CardiacPhase, ...]:
        """Every phase containing `wpc`, sorted by `(start_wpc, id)` -- possibly several.

        Empty, never `None`, for three different reasons, and a caller must
        not conflate them: a null `wpc` (post-natal, outside every cardiac
        morphogenetic phase by definition), a real `wpc` matched by no
        curated phase, and a real `wpc` that falls only inside a phase with
        no stated end. That third case is deliberate and unconditional: a
        phase with `end_wpc is None` (`EndBasis.NOT_STATED`) is excluded from
        this match at every `wpc`, not merely at wpc values "too far" past
        its start -- there is no source-given point at which this atlas could
        say it has stopped, so it must never be asserted as still running
        rather than picking an arbitrary cutoff nobody stated. Returning an
        empty tuple rather than raising in any of the three cases is what
        lets a page say "outside the curated window" or "post-natal" instead
        of dropping the row. Sorted deterministically so two callers -- and
        two builds of one commit -- see phases in the same order regardless
        of this list's own YAML declaration order or of `PYTHONHASHSEED`.
        """
        if wpc is None:
            return ()
        matches = [
            phase
            for phase in self.phases
            if phase.end_wpc is not None and phase.start_wpc <= wpc < phase.end_wpc
        ]
        return tuple(sorted(matches, key=lambda phase: (phase.start_wpc, phase.id)))
