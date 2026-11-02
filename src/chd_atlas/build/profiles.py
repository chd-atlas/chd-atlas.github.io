# src/chd_atlas/build/profiles.py
"""The figures this atlas derives from a developmental transcriptome.

Neither figure is a mirror column. The build derives both from the two
published tables, so a consumer's re-derivation and the build's are the same
arithmetic over the same bytes -- the pattern `omics.py` uses for `genes`,
where `count` is built from the attribution the rows publish so the two cannot
drift. Emitting them as columns instead would mean two computations of each
number with no way to tell which is right.

D39(b) is a one-level guarantee, not an unbounded one: a gene's percentile is
re-derivable from its median and the published breakpoints; the breakpoints
are not, because D32 forbids publishing the matrix they came from. Stated as a
trade rather than claimed as a proof.

This module carries the percentile band (Task 9), tau (Task 10), phase
assignment (Task 11), quantile shard emission (Task 12), the per-gene
assembly that reads both mirrors and turns them into one `ExpressionProfile`
per gene (Task 13), and `percentile_annotations` (Task 8b). Task 12 is the
other half of D39(b)'s bargain: Task 9's percentile is re-derivable only if
the breakpoints it was read against are themselves fetchable, and
`build_omics` never emits this table -- it skips every schema absent from its
own `_GENE_COLUMN`, and a quantile grid has no gene column at all.

`percentile_annotations` is the other half of a different bargain: `omics.py`
ranks a `profiles` row's stratified slice on the percentile the row is never
told at read time, so this module hands the figure back out, flattened, from
the exact `Placement`s Task 13 already computed -- read once here, never
derived a second time. See that function's own docstring.

Task 13's join is deliberately a LEFT join on (dataset, tissue, stage), never
an inner one. `profiles.stage` is nullable and `profile_quantiles.stage` is
not (see `tables.py`'s own comment on that column), so a null-stage row can
never have a quantile partner *by construction* -- an inner join drops it
silently (measured: two rows in, one out, no error), while a left join keeps
it and states why it carries no percentile. `ProfileGap` is that stated
reason, published as a slug rather than a rendered sentence: the exact wording
a reader sees is `build/pages.py`'s job (a later task), and baking English
prose into this payload would make the atlas author phrasing here that a page
might need to phrase two different ways.
"""

from __future__ import annotations

import bisect
import math
import sys
from collections.abc import Collection, Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any, Final, TypedDict

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import slug
from chd_atlas.models.dataset import Dataset, Stage
from chd_atlas.models.phases import CardiacPhaseFile
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table

# D39(a): the method travels with the number. Published beside every
# percentile so a consumer can re-derive it without reading this module's
# source to learn which side of a tied breakpoint was chosen.
LOOKUP_RULE: Final = "lowest percentile of a tied breakpoint (bisect_left)"


def percentile_of(value: float, breakpoints: Sequence[float]) -> int:
    """Where `value` falls in a 101-point grid, as a percentile 0-100.

    `bisect_left`, so a value tied with a run of identical breakpoints reports
    the *lowest* percentile in that run. The alternative reads the top of the
    tie, and roughly half a bulk library ties at zero -- which would publish an
    unexpressed gene above the median. Understating is the safe direction here,
    and D41 keeps those genes out of the published figure regardless (see
    `placement` below).

    Clamped at 100 rather than allowed to return 101 for a value above every
    breakpoint: the grid's top is the maximum, so exceeding it is a rounding
    artefact, not a 101st percentile. A value that merely *equals* the top
    breakpoint already reaches 100 through `bisect_left` alone -- the clamp is
    only ever exercised by a value that exceeds every breakpoint, which is a
    genuinely different input, not the same case asserted twice.
    """
    return min(bisect.bisect_left(breakpoints, value), 100)


def band(
    median: float,
    q25: float | None,
    q75: float | None,
    breakpoints: Sequence[float],
) -> tuple[int | None, int, int | None]:
    """The percentile band for one gene: (q25's percentile, median's, q75's).

    Published as a band rather than a point. At `n_samples` = 2 a point
    estimate of "top 4%" has an honest range of roughly top 2% to top 15%;
    publishing it costs one extra lookup per side and no sampling model,
    because `q25`/`q75` are already columns in `mirrors/profiles/`.

    `q25`/`q75` are `None` exactly when the mirror's are -- the schema nulls
    them below `n_samples` = 3, where a quartile of two points is not a
    quartile -- and `None` passes straight through rather than falling back to
    the median's own percentile, which would publish a band with no reported
    width as though it had been measured.
    """
    q25_percentile = percentile_of(q25, breakpoints) if q25 is not None else None
    q75_percentile = percentile_of(q75, breakpoints) if q75 is not None else None
    return (q25_percentile, percentile_of(median, breakpoints), q75_percentile)


class Placement(TypedDict):
    """One gene's percentile band in one (dataset, tissue, stage) cell.

    Every field a page needs to render the figure honestly travels with it,
    the same discipline `count_unit` enforces for every burden count:

    - `n_samples` -- the schema permits `n_samples` = 1, so "median" can be a
      single observation, and a page must be able to say so.
    - `n_genes` -- "top 4%" names no population without its denominator
      (D39(b)); this is the source matrix's whole gene universe, not the 92
      genes this atlas publishes.
    - `unit` -- a median with no unit is not a quantity.
    - `method` -- `LOOKUP_RULE` (D39(a)), so a consumer re-deriving this
      figure resolves a tied breakpoint the same way this module did.
    """

    q25_percentile: int | None
    median_percentile: int
    q75_percentile: int | None
    median_abundance: float
    unit: str
    n_samples: int
    n_genes: int
    method: str


def placement(
    median: float,
    q25: float | None,
    q75: float | None,
    breakpoints: Sequence[float],
    floor: float,
    unit: str,
    n_samples: int,
    n_genes: int,
) -> Placement | None:
    """The published figure for one gene in one cell, or `None` below the floor.

    D41: below the detection floor, a percentile is an artefact of the tied
    breakpoints at the bottom of the grid, not a measurement -- so `None` is
    returned rather than a number a page would render as though it meant
    something.

    **Judgement call: the comparison is `median < floor`, strictly.** A
    median exactly equal to the declared floor counts as detected, not below
    it. The floor is the dataset record's declared threshold for "present at
    all" (`Dataset.detection_floor`): the source's own where it declares one,
    and curated by this atlas where it declares none -- as for E-MTAB-6814,
    whose `floor_source` reads "Curated, not mirrored: the Expression Atlas
    declares no threshold for this experiment". Whoever drew it, the record
    names that value as the line, so the boundary value itself is a detected
    value under the record's own convention; `<=` would exclude the one point
    the record chose to draw the line at, which is stricter than the rule it
    states. Task 10's
    `specificity` gates tau's peak the same way (`peak_raw < floor`), so one
    gene cannot read "detected" under one figure and "below the floor" under
    the other from the same value. Pinned at the exact boundary by
    `test_the_floor_gate_is_strict_less_than`.
    """
    if median < floor:
        return None
    q25_percentile, median_percentile, q75_percentile = band(median, q25, q75, breakpoints)
    return Placement(
        q25_percentile=q25_percentile,
        median_percentile=median_percentile,
        q75_percentile=q75_percentile,
        median_abundance=median,
        unit=unit,
        n_samples=n_samples,
        n_genes=n_genes,
        method=LOOKUP_RULE,
    )


# D39(a): the scale is part of the method, not an implementation detail --
# see `specificity`'s docstring for the measurement that makes this load-
# bearing (a 10-fold enriched gene reads tissue-specific on linear and broad
# on log2).
TAU_SCALE: Final = "log2(x+1)"

# D39(a): named beside every tau so a consumer can re-derive the number from
# `medians` without reading this module's source. The dataset's own detection
# floor is deliberately absent from this description -- it gates whether tau
# is published at all (see `specificity`), but no longer changes what any
# organ contributes once it is, so re-deriving tau needs nothing but
# `medians` itself. Same one-level guarantee `LOOKUP_RULE` documents for the
# percentile side (see the module docstring).
#
# **It said "mean over organs" until 2026-08-21, and that is a different
# number.** A mean over organs divides by *n*; Yanai's tau -- what
# `specificity` below has always computed -- divides by *n - 1*. Measured over
# the committed corpus: all 1817 published specificity blocks match
# `sum/(n - 1)`, and the 4 that also match `sum/n` do so by coincidence. The
# code was right and the sentence was wrong, which is the worse way round for
# the one artifact whose entire purpose is letting a consumer re-derive the
# figure without reading this module. `docs/data-api.md` quoted the wrong
# string beside its own correct prose ("tau's denominator is n - 1"), and that
# self-consistency is why four review passes read past it.
#
# The denominator is now spelled out as an operation rather than named by a
# statistic ("divided by (n - 1)", not "mean"), so following the sentence and
# running the code cannot diverge again without the words visibly changing.
# `test_every_published_tau_is_reproducible_from_its_own_method_string` reads
# the divisor back out of this string and re-derives every published tau from
# `medians`; it fails on the old wording.
TAU_METHOD: Final = (
    "tau (Yanai et al. 2005): the sum over organs of (1 - x_i/x_max) "
    "divided by (n - 1), where n is the number of organs and "
    "x = log2(median+1); a negative median is clamped to 0 before the "
    "transform, and every organ's raw median is used even below the "
    "dataset's detection floor"
)


class Specificity(TypedDict):
    """Tau with everything needed to read it, and nothing that reads as a verdict.

    No adjective and no band: tau 0.71 is "intermediate" under the common
    banding (<=0.5 broad, >=0.8-0.9 specific), and choosing a threshold is a
    classification the atlas would author, which D39(c) forbids. `highest_in`
    is what makes a page's wording safe -- tau measures concentration, not
    location, so an organ-naming gloss keyed on tau alone can state the
    opposite of the truth for a gene concentrated elsewhere. (The page said
    "heart-preferential" until 2026-08-21; `pages._highest_cell`
    and `_tau_cell` record why the adjective went and why this gate stayed.)

    - `scale` -- names the transform (D39(a)); the same input reads
      tissue-specific on linear and broad on log2, so this is not decoration.
    - `method` -- the formula and citation (D39(a)), `LOOKUP_RULE`'s
      counterpart for tau.
    - `tissues` -- the ordered organ panel tau was computed over. A count
      cannot say *which* organs contributed; a real panel can contain a
      correlated pair (cerebrum/cerebellum), so this matters.
    - `n_tissues` -- organs *sampled*, not organs *detected*: an organ below
      the detection floor is retained here, never dropped.
    - `highest_in` -- the argmax organ, or `None` when two or more organs tie
      at the peak, where "highest in X" would be an arbitrary choice.
    - `medians` -- tau's actual inputs (D39(b)), because the bundle's own
      slice of a gene's organs is truncated elsewhere (`omics.select_top`).
      Published exactly as measured, and tau's own arithmetic uses these
      same raw values -- no internal flooring -- so recomputing tau from
      `medians` alone reproduces `tau` exactly. See `specificity`'s
      docstring for why an earlier version of this function broke that.
    """

    tau: float
    scale: str
    method: str
    tissues: tuple[str, ...]
    n_tissues: int
    highest_in: str | None
    medians: dict[str, float]


def specificity(medians: Mapping[str, float], floor: float) -> Specificity | None:
    """Yanai tau for one gene's per-organ medians at one developmental stage.

    tau = sum(1 - x_i/x_max) / (n - 1), computed on `x = log2(median + 1)`.
    Measured: a 10-fold enriched gene (heart 100, six organs at 10) scores
    0.900 on linear RPKM -- "tissue-specific" -- and 0.480 on log2 --
    "broadly expressed". Yanai et al. 2005 and the Kryuchkova-Mostacci &
    Robinson-Rechavi benchmark both compute tau on log expression; computing
    on linear would publish a number that reads tissue-specific for ordinary
    genes across the board.

    Returns `None`, never `0.0`, when tau is undefined:

    - fewer than two organs were sampled at this stage (the denominator is
      `n - 1`), or
    - the gene's peak across every sampled organ is below `floor` (the
      normaliser is `x_max`, and computing one from an unreliable peak would
      publish a ratio measured against noise).

    A genuinely ubiquitous gene -- every organ equal -- measures *exactly*
    0.000, so a guard that defaults an undefined tau to `0.0` (`n < 2`, or
    `except ZeroDivisionError`) publishes "expressed identically everywhere"
    for a gene measured in a single organ -- the opposite claim, and it lands
    on the earliest stages, where the organ panel is smallest.

    The peak's floor comparison is `peak_raw < floor`, strictly -- the same
    comparison `placement` uses (D41), so one gene cannot read "detected"
    under one figure and "below the floor" under the other from the same
    value; pinned at the boundary by `test_the_floor_gate_on_tau_is_strict_
    less_than`, matching Task 9's own boundary test.

    **The floor gates whether tau is published at all. It does not change
    what any organ contributes once tau is computed.** A non-peak organ
    below the floor is retained in `tissues`/`n_tissues` (so those keep
    meaning "organs sampled", not "organs detected" -- dropping is
    self-defeating, since the most heart-exclusive gene in the atlas,
    detected in heart alone, would fall to n=1 and lose tau entirely), *and
    its raw median is used in the sum exactly as measured* -- not floored to
    zero. An earlier version of this function floored it, reasoning that a
    sub-floor value is unreliable noise. Reverted on review, for two reasons:

    1. D39(b) requires tau to be re-derivable from the `medians` this same
       payload publishes, and `medians` is (and was always) published raw.
       Flooring internally while publishing raw inputs makes those two
       requirements contradict each other: a consumer who recomputes tau
       from the published `medians` gets a different number than the one
       published beside them (measured: 0.696 re-derived vs. 0.740
       published, same fixture). `test_tau_is_re_derivable_from_its_own_
       published_medians` makes this a checked invariant rather than a claim
       nobody watches.
    2. Flooring is not a neutral rounding choice: it measurably *raises*
       tau (0.740 vs. 0.696 on that same fixture -- `test_a_below_floor_
       nonzero_value_uses_its_raw_median_not_a_floored_zero`), pushing
       toward a stronger specificity claim than the data supports, in the
       same direction tau's own single-organ normaliser is already biased
       toward. Where two readings of an under-specified rule differ, this
       atlas takes the one that claims less.

    "Below the detection floor" means the dataset's record does not vouch
    for the value as a detection -- a floor this atlas curated where the
    source declares none (see `placement`) -- not that the source measured
    zero, and substituting zero for it is this atlas authoring a number the
    source did not report (the act D12/D39 exist to forbid). Gating *publication* on
    the floor and *rewriting a value* because of it are different acts; this
    function does the first and, deliberately, not the second.

    The one rewrite tau's arithmetic does make is a negative median, clamped
    to zero before the log transform: a negative abundance is not a value,
    and `log2` of one is not a number. `test_a_negative_median_is_clamped_
    to_zero_before_log2` pins this, including that it does not raise.

    `medians` is published alongside because the bundle's own slice of a
    gene's organs is truncated elsewhere (`omics.select_top`), and D39(b)
    requires tau's actual inputs to be reachable from the same payload that
    carries tau.
    """
    tissues = tuple(sorted(medians))
    if len(tissues) < 2:
        return None

    peak_raw = max(medians[tissue] for tissue in tissues)
    if peak_raw < floor:
        return None

    def _log2p1(value: float) -> float:
        # The only rewrite applied to any organ's raw median: a negative
        # abundance is not a value, and log2 of one is not a number. A
        # value below the detection floor is NOT rewritten here -- see the
        # docstring's "gate publication, don't rewrite the value" argument.
        return math.log2(max(value, 0.0) + 1.0)

    logged = {tissue: _log2p1(value) for tissue, value in medians.items()}
    x_max = _log2p1(peak_raw)
    if x_max <= 0:
        return None

    tau = sum(1.0 - value / x_max for value in logged.values()) / (len(tissues) - 1)

    leaders = [tissue for tissue in tissues if medians[tissue] == peak_raw]
    return Specificity(
        tau=tau,
        scale=TAU_SCALE,
        method=TAU_METHOD,
        tissues=tissues,
        n_tissues=len(tissues),
        highest_in=leaders[0] if len(leaders) == 1 else None,
        medians=dict(medians),
    )


class PhaseOutcome(StrEnum):
    """Why `assign_phase` returned the `phase_ids` (or lack of them) that it did.

    Four outcomes, never a bare emptiness check: a page that only checked
    whether `phase_ids` was empty would render one sentence for a post-natal
    stage, a stage past the curated window, and a stage the dataset's own
    record never declared -- three different facts, and `UNDECLARED` is not a
    biological statement at all, it is a curation gap PRF005 also reports from
    the other direction. `Resolution` in `chd_atlas/genes.py` is the
    precedent: a caller branches on the named outcome, never on which fields
    happen to be empty or on matching the display text in `reason`.
    """

    MATCHED = "matched"
    """`wpc` falls inside one or more declared cardiac phases; `phase_ids`
    names all of them. Human cardiac morphogenesis runs several processes
    concurrently (see `models/phases.py`'s module docstring for the measured
    concurrency table), so MATCHED means "one or more phases", never "exactly
    one" -- a caller that assumes a single match would silently drop every
    phase but the first for any stage inside an overlap."""

    OUTSIDE_WINDOW = "outside_window"
    """A real `wpc` exists, but no curated phase covers it -- before the
    first, after the last, or in an interior gap no phase's interval reaches
    (PRF006 flags the gap shape as a probable transcription slip; this is the
    runtime fact a page renders regardless of *why* the gap exists)."""

    POST_NATAL = "post_natal"
    """The stage's own `wpc` is null. Every `CardiacPhase` is prenatal by
    construction (`start_wpc`/`end_wpc` are both `gt=0`), so this is a
    different fact than OUTSIDE_WINDOW: there is no developmental age to
    place at all, and a wider curated window could never change the
    answer."""

    UNDECLARED = "undeclared"
    """`token` is not in the dataset's own `stages`, so there is no `Stage`
    to read a `wpc` from in the first place. PRF005 reports the reverse
    direction (a `profiles` row naming a stage the dataset record does not
    declare); this is what a caller sees asking about a token with nothing
    behind it, and `assign_phase` must not guess a biological placement for
    a stage it was never told about."""


@dataclass(frozen=True)
class PhaseAssignment:
    """Where one of a dataset's own stage tokens falls in the curated phase vocabulary.

    `phase_ids` is non-empty only for `PhaseOutcome.MATCHED`, and may name
    several phases at once -- see `CardiacPhaseFile.phases_for`, which this
    function reads from, and `PhaseOutcome.MATCHED`'s own docstring for why a
    single-match assumption would be wrong. `reason` is a short, renderable
    phrase explaining an empty `phase_ids`, and is `None` exactly when
    `phase_ids` is not empty -- a page renders one or the other, never both,
    and never neither.
    """

    token: str
    outcome: PhaseOutcome
    phase_ids: tuple[str, ...] = ()
    reason: str | None = None


def assign_phase(
    token: str, stages: Sequence[Stage], phases: CardiacPhaseFile | None
) -> PhaseAssignment:
    """Place one of a dataset's own stage tokens in the curated phase vocabulary.

    Checked in this order, each a stronger claim than the last is absent:

    1. `token` is not in `stages` at all -- `UNDECLARED`. Nothing about
       development can be said for a token with no `Stage` behind it, so this
       is decided before either of the wpc-based questions below are even
       asked.
    2. The matching `Stage.wpc` is `None` -- `POST_NATAL`. Every `CardiacPhase`
       is prenatal by construction, so a null wpc can never match one
       regardless of what `phases` declares.
    3. `phases` is absent entirely, or declares zero phases -- `OUTSIDE_WINDOW`,
       reason "outside the curated window". Checked before indexing into
       `phases.phases` for the boundary comparison in step 4, which would
       raise `IndexError` on an empty sequence otherwise -- and this is not a
       hypothetical: it is the committed corpus's own state today.
    4. Otherwise `phases.phases_for(wpc)` is tried. One or more hits is
       `MATCHED`, naming every matching phase -- never just the first, since
       overlapping phases are the ordinary case (see `models/phases.py`'s
       module docstring). No hits is `OUTSIDE_WINDOW` again, with a reason
       naming *where*: "before the curated window" (earlier than every
       phase's `start_wpc` -- every phase has one, `CardiacPhase.start_wpc`
       is never null), "after the curated window" (at or past the latest
       *known* `end_wpc` across every phase that has one -- the common case,
       since a real developmental series runs well past any morphogenetic
       window, not an edge case), or the same generic phrase as step 3 for
       the rarer case of an interior gap no phase's interval reaches, which
       PRF006 already names as a probable transcription slip. A phase with
       no stated end (`EndBasis.NOT_STATED`) contributes its `start_wpc` to
       the "before" comparison -- a real, sourced fact -- but nothing to the
       "after" one: there is no known end for it to be the latest of, the
       same exclusion `phases_for` itself applies. If literally every phase
       has an unstated end, there is no "after" boundary to compare against
       at all, and every wpc past the last start reads as the generic
       "outside the curated window" rather than crashing on an empty `max()`.
    """
    stage = next((candidate for candidate in stages if candidate.token == token), None)
    if stage is None:
        return PhaseAssignment(
            token=token,
            outcome=PhaseOutcome.UNDECLARED,
            reason="stage not declared by this dataset",
        )
    if stage.wpc is None:
        return PhaseAssignment(token=token, outcome=PhaseOutcome.POST_NATAL, reason="post-natal")
    if phases is None or not phases.phases:
        return PhaseAssignment(
            token=token,
            outcome=PhaseOutcome.OUTSIDE_WINDOW,
            reason="outside the curated window",
        )
    found = phases.phases_for(stage.wpc)
    if found:
        return PhaseAssignment(
            token=token,
            outcome=PhaseOutcome.MATCHED,
            phase_ids=tuple(phase.id for phase in found),
        )

    starts = [phase.start_wpc for phase in phases.phases]
    ends = [phase.end_wpc for phase in phases.phases if phase.end_wpc is not None]
    if stage.wpc < min(starts):
        reason = "before the curated window"
    elif ends and stage.wpc >= max(ends):
        reason = "after the curated window"
    else:
        reason = "outside the curated window"
    return PhaseAssignment(token=token, outcome=PhaseOutcome.OUTSIDE_WINDOW, reason=reason)


def build_profile_quantiles(root: Path, emitter: Emitter) -> dict[str, str]:
    """Publish every quantile grid, so a percentile this atlas derives can be checked.

    D39(b) is unmet without this. A gene page states a percentile whose input
    is the whole source transcriptome (tens of thousands of genes), and D32
    forbids re-hosting that matrix -- so `mirrors/profile_quantiles/<accession>.tsv`
    exists solely to publish the 101 breakpoints a consumer needs to place any
    abundance and reproduce the number themselves. `build_omics` never emits
    it: that function skips every schema absent from its own `_GENE_COLUMN`
    (`omics.py`'s `if schema_name not in _GENE_COLUMN: continue`), and a
    quantile grid has no gene column -- there is no gene to attribute a
    breakpoint row to. Without this function the table is mirrored,
    schema-validated, sort-checked and sha256'd, and reaches no published
    byte: the one file whose entire purpose is auditability would be the one
    thing a consumer could not fetch, and nothing would catch it -- the
    validators check the mirror against itself, `percentile_of`'s arithmetic
    is tested against unpublished input either way, and a `diff -rq` between
    two builds shows a file that appears, never one that should have.

    Emits `omics/profile_quantiles/<accession>.json` as
    `{"table": "profile_quantiles", "rows": [...]}` -- the same envelope
    `build_omics` writes for its own shards, so a consumer reads one shape
    regardless of which table it fetched. The accession is the shard's own
    filename stem, exactly as `build_omics` reads a dataset accession from a
    `profiles`/`expression`/`proteomics`/`phospho` shard's name; it is put
    through `paths.slug` for the same reason every shard stem is -- the stem
    becomes a URL, and a space or a colon in one would need escaping before it
    could be fetched.

    Rows are sorted by the table's own canonical `sort_key`
    (`dataset, tissue, stage, percentile`) before being written, never trusted
    to already be in that order: `read_table` preserves a TSV's file order
    verbatim (`pl.read_csv`, no grouping or sort of its own), so an unsorted
    row order on disk would otherwise publish unsorted, silently, on a build
    that never re-validates a mirror it has already read. Sorting on the full
    key rather than on `percentile` alone matters the moment one shard holds
    more than one `(tissue, stage)` grid, which is the ordinary shape for a
    real dataset (many organs, many stages) rather than the single-grid shape
    of the test fixture -- sorting on `percentile` alone would interleave
    every grid in the shard by breakpoint value instead of keeping each one
    contiguous and internally ordered.

    Shards are emitted in `mirror_paths`' own order, which is sorted by
    filename -- so two builds of one commit write these files in the same
    order, though (unlike row order within one file) that has no bearing on
    any file's *content*.

    Returns `{dataset accession: relative shard path}`. Read by
    `gene_expression_profiles` below, which threads a gene's
    `expression_profile` bundle key to the grid its percentile came from, the
    way `ModalitySummary.shards` lets a gene bundle reach an omics shard
    today. Returning the exact path this function wrote -- rather than
    leaving the caller to reconstruct it via a second call to `slug` -- is
    what keeps that link from becoming a second computation that could drift
    from the first, the same discipline `omics.py`'s own module docstring
    states for `count`.
    """
    shards: dict[str, str] = {}
    for path, schema_name in mirror_paths(root):
        if schema_name != "profile_quantiles":
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None:
            # Unreadable is `validate_table`'s to report against this same
            # path; failing the build here would say it a second time with
            # less context. Unreachable behind `build_site`'s validation
            # gate in the same way `build_omics`'s equivalent branch is --
            # kept as a guard on that bypassed gate rather than trusted to
            # stay true as the two modules evolve apart.
            continue
        rows = frame.sort(list(TABLE_SCHEMAS[schema_name].sort_key)).to_dicts()
        accession = path.stem
        relative = f"omics/profile_quantiles/{slug(accession)}.json"
        emitter.write_json(relative, {"table": schema_name, "rows": rows})
        shards[accession] = relative
    return shards


class ProfileGap(StrEnum):
    """Why a derived expression figure (a percentile placement or tau) is absent.

    One vocabulary shared between `TissueProfileEntry.not_placed_reason` and
    `StageProfileEntry.specificity_unavailable_reason`, rather than two,
    because the first two members mean exactly the same fact in both places
    -- the dataset itself carries no floor to gate on -- and a reader
    comparing the two fields on one page should not have to learn two
    spellings of one fact. Published as a slug rather than a sentence; see
    the module docstring for why the wording itself is not this module's job.
    """

    DATASET_NOT_REGISTERED = "dataset_not_registered"
    """The `profiles.tsv` row names a dataset accession with no `Dataset`
    record in `corpus.datasets` at all. Reachable on a validated repository,
    not only a bypassed gate: `validate_profile_references`'s own
    `_prf004_issues` docstring states this is "out of scope for this check"
    -- there is no PRF code for it."""

    FLOOR_UNDECLARED = "detection_floor_undeclared"
    """The dataset record exists but `detection_floor` is `None`. Unreachable
    for a `design="profile"` dataset in a corpus that passed pydantic
    validation (`Dataset.a_profile_dataset_is_fully_declared` refuses to load
    one with no floor) -- reachable for a `design="contrast"` record, whose
    floor is optional and defaults to `None`, referenced by a `profiles.tsv`
    row naming the wrong accession. Guarded anyway, per this project's stance
    against trusting "the gate already checked this" to stay true forever."""

    NO_QUANTILE_GRID = "no_quantile_grid"
    """No complete 101-point breakpoint grid exists for this
    (dataset, tissue, stage) cell -- includes every null-stage row, which can
    never have one by construction, and a grid PRF003 would have flagged as
    non-monotone had validation run over one that is merely incomplete
    (PRF003 checks direction, not completeness). Not to be confused with
    `BELOW_DETECTION_FLOOR` below, which requires a grid to exist at all."""

    BELOW_DETECTION_FLOOR = "below_detection_floor"
    """A grid and a floor both exist; `placement()` itself refused (D41):
    the median is below the dataset's own detection floor."""

    ONE_ORGAN_SAMPLED = "one_organ_sampled"
    """Tau is undefined at n<2 organs (the denominator is n-1)."""

    PEAK_BELOW_DETECTION_FLOOR = "peak_below_detection_floor"
    """Tau is undefined when every sampled organ's median is below the floor
    (the normaliser is the peak)."""

    NOT_ON_A_LOG_AXIS = "not_on_a_log_axis"
    """The atlas placed this cell, and a figure still cannot draw it: the
    median is `<= 0`, which a logarithmic axis has no position for.

    Distinct from `BELOW_DETECTION_FLOOR`, which says the atlas declined to
    place the measurement at all. This one says it placed it and the picture
    cannot show it, and the two owe a reader different sentences. Unreachable
    while `Dataset.detection_floor` is `gt=0` -- a positive floor cannot admit
    a non-positive median -- and kept because the alternative is
    `pages._Median` carrying neither a value nor a reason, which publishes
    "no reason was recorded for this gap" about a cell whose percentile is
    printed three lines below."""

    UNDEFINED = "undefined"
    """`specificity()` refused for a reason neither of the two checks above
    predicts -- reachable only for a degenerate floor `<= 0` whose peak is
    exactly `0.0` (`specificity`'s own docstring and
    `test_a_degenerate_zero_floor_refuses_rather_than_dividing_by_zero`),
    where `peak_raw < floor` does not fire but the internal log-transform
    normaliser does. Kept as a residual rather than asserting a specific
    cause this module has not actually distinguished."""


# Reason text for a stage with no token at all -- distinct from
# `PhaseOutcome.UNDECLARED`, which names a real token the dataset's own record
# does not declare. `assign_phase` requires a `str` token; a null-stage row
# has none to offer it, so this is stated directly rather than manufactured by
# calling `assign_phase("")` and hoping no dataset ever declares that token.
_NO_STAGE_PHASE_REASON: Final = "no developmental stage recorded for this measurement"


class TissueProfileEntry(TypedDict):
    """One gene's raw measurement in one (dataset, tissue, stage) cell, and its
    derived placement if one could be computed.

    `median_abundance`/`unit`/`n_samples` are published here even though
    `Placement` (when present) repeats them -- so a consumer reads one shape
    for "what was measured" whether or not a percentile could be derived from
    it. `not_placed_reason` is `None` exactly when `placement` is not,
    matching the discipline `PhaseAssignment.reason` already keeps.
    """

    tissue: str
    median_abundance: float
    unit: str
    n_samples: int
    placement: Placement | None
    # `ProfileGap`, not `str`: a `StrEnum` member serialises to the same
    # bytes, so the published payload is unchanged, and mypy starts
    # checking the vocabulary the enum was written to constrain. Bare
    # `str` meant nothing stopped a specificity-side member landing here,
    # where `_PLACEMENT_GAP_CLAUSE` would miss it and the page would print
    # "no reason was recorded for this gap" about a gap whose reason was.
    not_placed_reason: ProfileGap | None


class PhaseInfo(TypedDict):
    """Where one stage falls in the curated phase vocabulary, or why not.

    `outcome` carries `PhaseOutcome.value` for a real, non-null stage token;
    `None` exactly for a stage with no token at all (a null `profiles.stage`),
    which is not one of `PhaseOutcome`'s four members because no token means
    there was nothing to ask `assign_phase` about in the first place.

    `phase_ids` may name more than one phase -- human cardiac morphogenesis
    runs several processes concurrently, so a stage legitimately matches every
    phase whose interval contains its `wpc` (see `models/phases.py`'s module
    docstring for the measured concurrency this exists to preserve). Empty
    exactly when `reason` is not `None`.
    """

    outcome: str | None
    phase_ids: tuple[str, ...]
    reason: str | None


class StageProfileEntry(TypedDict):
    """One (dataset, stage) cross-section: tau across organs, the phase, and
    each organ's own placement.

    `specificity` sits here, never inside a `TissueProfileEntry` -- tau is
    per (gene, dataset, stage) across every organ sampled there, not a
    per-organ figure. `specificity_unavailable_reason` is `None` exactly when
    `specificity` is not, the same pairing `not_placed_reason` keeps for
    `placement`.
    """

    stage: str | None
    phase: PhaseInfo
    specificity: Specificity | None
    specificity_unavailable_reason: ProfileGap | None
    tissues: list[TissueProfileEntry]


class DatasetProfileEntry(TypedDict):
    """One dataset's whole contribution to one gene: every stage it covers,
    and the shard a consumer re-derives its percentiles against.

    `quantile_shard` is `None` when `quantile_shards` (this dataset's own
    `profiles.dataset` column value) has no entry -- a corpus mid-curation,
    the gate bypassed, or `build_profile_quantiles`'s returned key naming a
    shard's *filename* rather than the accession inside it, if the two ever
    disagree (nothing validates that they must) -- and never omitted,
    matching every other "always present" key in this atlas. The percentile
    itself is unaffected either way: `_read_grids` below keys on the `dataset`
    *column*, not the filename, so a mismatched shard name costs only this
    link, never the figure it would have linked to.
    """

    dataset: str
    quantile_shard: str | None
    stages: list[StageProfileEntry]


class ExpressionProfile(TypedDict):
    """One gene's whole developmental expression profile.

    `datasets` is empty for a gene no profile dataset's mirror covers -- the
    always-present shape `bundles.py` publishes for exactly that gene,
    `EMPTY_EXPRESSION_PROFILE` below.
    """

    datasets: list[DatasetProfileEntry]


EMPTY_EXPRESSION_PROFILE: Final[ExpressionProfile] = ExpressionProfile(datasets=[])

# One cell's key throughout this module. `str | None` for `stage` because a
# lookup may legitimately ask about a null-stage row; a real grid never has
# one (`profile_quantiles.stage` is not nullable), so no key ever gets
# *inserted* with a `None` third element, but the type has to admit the
# lookup side too.
_Cell = tuple[str, str, str | None]


@dataclass(frozen=True)
class _ProfileRow:
    """One `mirrors/profiles/*.tsv` row, cast out of polars' `Any`-typed dict."""

    dataset: str
    tissue: str
    stage: str | None
    median_abundance: float
    unit: str
    q25: float | None
    q75: float | None
    n_samples: int

    @staticmethod
    def from_record(record: Mapping[str, Any]) -> _ProfileRow:
        return _ProfileRow(
            dataset=str(record["dataset"]),
            tissue=str(record["tissue"]),
            stage=str(record["stage"]) if record["stage"] is not None else None,
            median_abundance=float(record["median_abundance"]),
            unit=str(record["unit"]),
            q25=float(record["q25"]) if record["q25"] is not None else None,
            q75=float(record["q75"]) if record["q75"] is not None else None,
            n_samples=int(record["n_samples"]),
        )


def _read_grids(root: Path) -> dict[_Cell, dict[int, tuple[float, int]]]:
    """Every complete-or-not breakpoint grid, keyed by (dataset, tissue, stage).

    `{percentile: (value, n_genes)}` per cell, deliberately not yet reduced to
    a 101-element list: `_breakpoints` below is what decides whether a cell's
    grid is complete enough to place a gene against, and it needs the
    percentile keys to check that, not just the values in whatever order this
    function happened to see them in.
    """
    grids: dict[_Cell, dict[int, tuple[float, int]]] = {}
    for path, schema_name in mirror_paths(root):
        if schema_name != "profile_quantiles":
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None:
            # Unreadable is `validate_table`'s to report against this same
            # path; see `build_profile_quantiles`'s identical guard.
            continue
        for record in frame.iter_rows(named=True):
            dataset, tissue, stage = record["dataset"], record["tissue"], record["stage"]
            percentile, value, n_genes = record["percentile"], record["value"], record["n_genes"]
            if dataset is None or tissue is None or stage is None:
                continue
            if percentile is None or value is None or n_genes is None:
                continue
            cell: _Cell = (str(dataset), str(tissue), str(stage))
            grids.setdefault(cell, {})[int(percentile)] = (float(value), int(n_genes))
    return grids


def _breakpoints(grid: Mapping[int, tuple[float, int]] | None) -> tuple[list[float], int] | None:
    """The 101-point lookup table for one cell, or `None` if it is incomplete.

    `percentile_of` (`bisect.bisect_left`) assumes position i holds
    percentile i's own value -- this module's own docstring for `Placement`
    says so -- so a grid missing even one of 0..100 would silently place a
    gene against a *misaligned* array rather than merely cost one point of
    precision. No PRF check pins grid completeness today (PRF001-009 check
    units, cell presence and monotonicity, never the count of rows in one
    cell), so this function refuses rather than build one.

    `n_genes` is read from percentile 100's own row -- guaranteed present
    once completeness is confirmed -- rather than an arbitrary row, so a
    grid whose `n_genes` genuinely varies row to row (nothing enforces that it
    does not) gives a reproducible answer rather than one that depends on
    dict iteration order.
    """
    if grid is None or set(grid) != set(range(101)):
        return None
    return [grid[percentile][0] for percentile in range(101)], grid[100][1]


def _dataset_gap(dataset: Dataset | None) -> ProfileGap | None:
    """Why a dataset cannot gate a percentile or tau, or `None` when it can.

    Two facts kept distinct rather than folded into one "no floor" catch-all
    -- see `ProfileGap.DATASET_NOT_REGISTERED` and `.FLOOR_UNDECLARED` for why
    each is independently reachable.
    """
    if dataset is None:
        return ProfileGap.DATASET_NOT_REGISTERED
    if dataset.detection_floor is None:
        return ProfileGap.FLOOR_UNDECLARED
    return None


def _tissue_entry(
    row: _ProfileRow,
    floor: float | None,
    dataset_gap: ProfileGap | None,
    grid: Mapping[int, tuple[float, int]] | None,
) -> TissueProfileEntry:
    """One tissue's raw measurement, placed if a floor and a grid both exist.

    `floor` and `dataset_gap` are threaded in rather than re-derived from a
    `Dataset` here, because the caller already computed both once per
    dataset; the `elif floor is None` branch below is therefore a guard on an
    invariant the caller establishes, not a path this module's own tests
    reach independently of it -- kept explicit (never `assert`, which `-O`
    strips) rather than trusted to hold silently.
    """
    placed: Placement | None = None
    reason: ProfileGap | None = None
    if dataset_gap is not None:
        reason = dataset_gap
    elif floor is None:
        reason = ProfileGap.FLOOR_UNDECLARED
    else:
        complete = _breakpoints(grid)
        if complete is None:
            reason = ProfileGap.NO_QUANTILE_GRID
        else:
            breakpoints, n_genes = complete
            placed = placement(
                median=row.median_abundance,
                q25=row.q25,
                q75=row.q75,
                breakpoints=breakpoints,
                floor=floor,
                unit=row.unit,
                n_samples=row.n_samples,
                n_genes=n_genes,
            )
            if placed is None:
                reason = ProfileGap.BELOW_DETECTION_FLOOR
    return TissueProfileEntry(
        tissue=row.tissue,
        median_abundance=row.median_abundance,
        unit=row.unit,
        n_samples=row.n_samples,
        placement=placed,
        not_placed_reason=reason,
    )


def _specificity_entry(
    medians: Mapping[str, float], floor: float | None, dataset_gap: ProfileGap | None
) -> tuple[Specificity | None, ProfileGap | None]:
    """Tau over one stage's per-organ medians, or the reason it is absent.

    Checked in the same order `specificity()` itself would refuse, so the
    reason named is the first one that actually applies rather than a
    generic catch-all: no floor, then too few organs, then a peak below the
    floor, then (a residual this module cannot itself predict) whatever
    `specificity()` still refused for.
    """
    if dataset_gap is not None:
        return None, dataset_gap
    if floor is None:
        return None, ProfileGap.FLOOR_UNDECLARED
    if len(medians) < 2:
        return None, ProfileGap.ONE_ORGAN_SAMPLED
    if max(medians.values()) < floor:
        return None, ProfileGap.PEAK_BELOW_DETECTION_FLOOR
    result = specificity(medians, floor)
    if result is None:
        return None, ProfileGap.UNDEFINED
    return result, None


def _phase_entry(
    stage_token: str | None, dataset: Dataset | None, phases: CardiacPhaseFile | None
) -> PhaseInfo:
    """Where one stage token falls in the curated phase vocabulary.

    A null `stage_token` has nothing to hand `assign_phase` -- it requires a
    `str` -- so this states the fact directly rather than manufacturing a
    call with an empty string and hoping no dataset ever declares that as a
    real token.
    """
    if stage_token is None:
        return PhaseInfo(outcome=None, phase_ids=(), reason=_NO_STAGE_PHASE_REASON)
    stages = dataset.stages if dataset is not None else ()
    assignment = assign_phase(stage_token, stages, phases)
    return PhaseInfo(
        outcome=assignment.outcome.value, phase_ids=assignment.phase_ids, reason=assignment.reason
    )


def _stage_entry(
    dataset_id: str,
    stage_token: str | None,
    rows: Sequence[_ProfileRow],
    dataset: Dataset | None,
    dataset_gap: ProfileGap | None,
    phases: CardiacPhaseFile | None,
    grids: Mapping[_Cell, dict[int, tuple[float, int]]],
) -> StageProfileEntry:
    """One (dataset, stage) cross-section, tissues sorted for determinism."""
    floor = dataset.detection_floor if dataset is not None else None
    ordered = sorted(rows, key=lambda row: row.tissue)
    tissues: list[TissueProfileEntry] = []
    medians: dict[str, float] = {}
    for row in ordered:
        grid = grids.get((dataset_id, row.tissue, stage_token))
        tissues.append(_tissue_entry(row, floor, dataset_gap, grid))
        medians[row.tissue] = row.median_abundance
    specificity_result, specificity_reason = _specificity_entry(medians, floor, dataset_gap)
    return StageProfileEntry(
        stage=stage_token,
        phase=_phase_entry(stage_token, dataset, phases),
        specificity=specificity_result,
        specificity_unavailable_reason=specificity_reason,
        tissues=tissues,
    )


# Where a token the dataset never declared sorts: after every declared stage,
# before the null-stage bucket.
#
# Not `len(order)`. This comment said that would collide "whenever a curator
# numbers from 0 or leaves a gap", and both halves were wrong -- `Stage.order`
# is `ge=1`, so 0-based numbering cannot happen, and a gap is the one shape
# that is safe. Measured: `{a:1, b:2, c:3}` gives `len(order) == 3`, which is
# `c`'s own position, while `{a:1, b:2, c:5}` gives 3 and collides with
# nothing. On the committed corpus `len(order)` is 21, which is `elderly`.
#
# So the sentinel collides in the **ordinary** case -- contiguous 1..N, where
# `len(order)` is always the last declared position -- and an undeclared token
# would tie with the final stage and fall back to the token tie-break, which
# is the alphabetical ordering this whole field exists to remove.
_UNDECLARED_STAGE: Final = sys.maxsize


def _stage_sort_key(token: str | None, order: Mapping[str, int]) -> tuple[bool, int, str]:
    """Curated chronological order, undeclared tokens after it, null last.

    **This sorted alphabetically until 2026-08-20**, and its own docstring
    said so approvingly. Alphabetical is deterministic and reproducible and
    wrong for a time series: it published `4 week post conception` eighth,
    after `19 week post conception`, and `elderly` second of the eight
    post-natal stages -- in the HTML page and in the published JSON alike.
    The defect survived an adversarial review, a mutation matrix and a
    promotion to `main`, because sentences cannot show a trajectory and every
    block was individually correct.

    `order` is the dataset's own `{token: Stage.order}` map, so the curator's
    curated chronology is what sorts the axis. `Stage.wpc` cannot do this job:
    it is null for every post-natal stage, so the eight tokens that most need
    ordering have no number at all.

    A `bool` first component rather than interleaving `None` with ints --
    Python raises `TypeError` comparing `None` to `int` directly, so a plain
    key would crash the moment one dataset carries a null-stage row alongside
    a declared one. `token` is the final component so two undeclared tokens
    still order deterministically rather than by dict insertion.
    """
    if token is None:
        return (True, 0, "")
    return (False, order.get(token, _UNDECLARED_STAGE), token)


def _dataset_entry(
    dataset_id: str,
    rows: Sequence[_ProfileRow],
    dataset_by_id: Mapping[str, Dataset],
    phases: CardiacPhaseFile | None,
    grids: Mapping[_Cell, dict[int, tuple[float, int]]],
    quantile_shards: Mapping[str, str],
) -> DatasetProfileEntry:
    """One dataset's whole contribution to one gene, stages sorted for determinism."""
    dataset = dataset_by_id.get(dataset_id)
    dataset_gap = _dataset_gap(dataset)
    by_stage: dict[str | None, list[_ProfileRow]] = {}
    for row in rows:
        by_stage.setdefault(row.stage, []).append(row)
    # `dataset` is `None` for an accession no curated record names
    # (`ProfileGap.DATASET_NOT_REGISTERED`), which degrades every token to
    # undeclared and therefore to token order. That is the honest fallback: an
    # unregistered dataset has no curated chronology to sort by, and the gap is
    # already stated on every cell it affects.
    order = {stage.token: stage.order for stage in dataset.stages} if dataset is not None else {}
    stages = [
        _stage_entry(dataset_id, token, by_stage[token], dataset, dataset_gap, phases, grids)
        for token in sorted(by_stage, key=lambda token: _stage_sort_key(token, order))
    ]
    return DatasetProfileEntry(
        dataset=dataset_id,
        quantile_shard=quantile_shards.get(dataset_id),
        stages=stages,
    )


def gene_expression_profiles(
    root: Path,
    datasets: Sequence[Dataset],
    phases: CardiacPhaseFile | None,
    quantile_shards: Mapping[str, str],
) -> dict[str, ExpressionProfile]:
    """Read both profile mirrors and assemble one `ExpressionProfile` per gene.

    Returns only the genes `mirrors/profiles/*.tsv` actually mentions -- a
    gene absent from the return is a gene this function never saw, not a
    published claim of "no data". `bundles.py` is what turns that absence
    into the always-present empty shape a gene bundle publishes; deliberately
    not done here, so this function's own contract stays "report what the
    mirrors say" rather than "decide what a bundle looks like".

    Not gated on gene publication (`build.validity.published_genes()`) at
    all -- this function does not take that population as an argument. A
    profiles row can name a gene that is registered but never published (the
    same 154-vs-92 asymmetry `mirrors/genes.tsv` already has), and computing
    its facts anyway is harmless: `build_genes` is the single place that
    restricts what reaches a bundle, by iterating `published` alone, and
    restricting *here* as well would risk the two gates disagreeing about
    which genes exist.

    `datasets` is `corpus.datasets`; a dataset id a `profiles.tsv` row names
    with no matching record here is not an error this function raises over
    -- see `ProfileGap.DATASET_NOT_REGISTERED`. `quantile_shards` is
    `build_profile_quantiles`'s own return, reused rather than recomputed so
    the two cannot name different files for one accession.
    """
    # `str(dataset.id)`, not the bare `AccessionId`: `profiles.dataset` is a
    # plain string column, and `Mapping`'s key type is invariant, so a
    # `dict[AccessionId, Dataset]` does not satisfy `Mapping[str, Dataset]`
    # below even though every `AccessionId` is itself a `str` at runtime --
    # the same normalisation `runner.py` already applies for `cardiac_tissues`.
    dataset_by_id = {str(dataset.id): dataset for dataset in datasets}
    grids = _read_grids(root)

    rows_by_gene: dict[str, list[_ProfileRow]] = {}
    for path, schema_name in mirror_paths(root):
        if schema_name != "profiles":
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None:
            # Unreadable is `validate_table`'s to report against this same
            # path; see `build_profile_quantiles`'s identical guard.
            continue
        for record in frame.iter_rows(named=True):
            gene = record["gene"]
            if gene is None:
                continue
            rows_by_gene.setdefault(str(gene), []).append(_ProfileRow.from_record(record))

    result: dict[str, ExpressionProfile] = {}
    for gene, rows in rows_by_gene.items():
        by_dataset: dict[str, list[_ProfileRow]] = {}
        for row in rows:
            by_dataset.setdefault(row.dataset, []).append(row)
        datasets_entries = [
            _dataset_entry(
                dataset_id, by_dataset[dataset_id], dataset_by_id, phases, grids, quantile_shards
            )
            for dataset_id in sorted(by_dataset)
        ]
        result[gene] = ExpressionProfile(datasets=datasets_entries)
    return result


def percentile_annotations(
    profiles: Mapping[str, ExpressionProfile],
) -> dict[tuple[str, str, str, str | None], int]:
    """Every published placement's percentile, flattened to one lookup.

    Keyed `(gene, dataset, tissue, stage)`. Not a second computation: every
    value here is a `Placement["median_percentile"]` read back out of the
    exact `ExpressionProfile`s the caller is about to hand `build_genes` for
    publication -- there is exactly one call to `placement()` in the whole
    build, inside `gene_expression_profiles` above, and this function only
    ever reads its answer back rather than re-deriving it from a median and a
    grid a second way.

    `build_omics` is the only consumer (`omics._profile_percentile`), which
    annotates a `profiles` row with this figure before `select_top` ranks the
    cardiac series (`omics._by_percentile_then_stage`). So the rank a reader's
    bundle preview is chosen by and the percentile the bundle itself publishes
    cannot drift apart by construction -- not merely by agreeing today, the
    way two independent computations of the same number would.

    A cell with no placement (`not_placed_reason` set -- below the detection
    floor, no quantile grid, the dataset itself unregistered, and so on)
    contributes no entry. That is not a loss: `_by_percentile_then_stage`
    already treats a missing key as "sorts last", the same "not evidence of
    high expression" reading `placement`'s own `None` return gets everywhere
    else this atlas renders it.
    """
    annotations: dict[tuple[str, str, str, str | None], int] = {}
    for gene, profile in profiles.items():
        for dataset_entry in profile["datasets"]:
            for stage_entry in dataset_entry["stages"]:
                for tissue_entry in stage_entry["tissues"]:
                    placed = tissue_entry["placement"]
                    if placed is None:
                        continue
                    key = (
                        gene,
                        dataset_entry["dataset"],
                        tissue_entry["tissue"],
                        stage_entry["stage"],
                    )
                    annotations[key] = placed["median_percentile"]
    return annotations


class ProfileCensus(TypedDict):
    """The two developmental-expression figures this atlas derives, once, for
    two consumers: `index.html`, whose card states the genes (its datasets card
    went with the doorway page, 2026-10-27), and `manifest.json`'s
    `counts.profile_datasets`/`counts.profile_genes`.

    A `TypedDict`, not a scalar pair, for the same reason every other derived
    shape in this module is one: a caller unpacking `profile_census(...)["genes"]`
    at two call sites cannot transpose the two fields the way two positional
    return values could.

    Lives beside `gene_expression_profiles` rather than in `landing.py` or
    `manifest.py`, the same placement `burden.py` gives `BurdenCensus` and for
    the same reason: a figure with two consumers in different layers belongs to
    neither of them.
    """

    genes: int
    datasets: int


def profile_census(
    profiles: Mapping[str, ExpressionProfile], published: Collection[str]
) -> ProfileCensus:
    """What the developmental-expression layer *this site publishes* amounts to.

    **Both figures are restricted to `published`, with no exception.** Unlike
    `burden_census`'s `families` -- deliberately left unrestricted because a
    cohort family is a property of the corpus's sample collections, and
    dropping one whose studies tested no published gene can split it into
    what then reads as two independent datasets (D33) -- a profile dataset
    carries no such connectedness to protect. There is no reason for either
    figure here to read wider than what a consumer can actually reach.

    `genes` counts published genes whose `expression_profile` is non-empty --
    never `len(profiles)`. `gene_expression_profiles` returns an entry for
    every gene `mirrors/profiles/*.tsv` mentions "published or not" (its own
    docstring), the same 154-vs-92 asymmetry `mirrors/genes.tsv` already has
    for burden, so counting the mapping directly would advertise
    developmental-expression evidence behind genes this API publishes no
    bundle for.

    `datasets` counts the distinct dataset accessions named by those same
    genes' own `expression_profile.datasets` entries -- never
    `len(corpus.datasets)`, which already means something else: *every*
    registered omics dataset, `profile`-design and `contrast`-design alike,
    with no restriction to `published` at all (that count predates this
    layer and is deliberately left as-is; see `manifest.py`'s schema-history
    comment for 2.11). Two keys both named "datasets" and both meaning
    something different is exactly the confusion `cohort_families` was kept
    apart from `independent_datasets` to avoid, so this counts only
    accessions actually reachable from a published gene's own bundle.

    This paragraph read "both figures are 0 on every corpus this atlas has
    built so far -- no `profiles` mirror has ever been committed" until
    2026-08-20. `mirrors/profiles/E-MTAB-6814.tsv` landed on 2026-08-19 and
    both figures moved: measured on the committed corpus, `genes` is 92 and
    `datasets` is 1.

    **The fixture is still what proves the restriction, for a narrower
    reason.** One wrong implementation is now visible on a real build and one
    is not. Counting the mirror gives 154, because
    `mirrors/profiles/E-MTAB-6814.tsv` covers every registered gene; reading
    `len(published)` and `len(corpus.datasets)` gives 92 and 1, which is
    exactly what a correct implementation publishes. Two figures that are
    equal today are one figure to every test (CLAUDE.md section 4.30), so the
    case that separates them is hand-built: a gene outside `published`, and a
    dataset only that gene cites. See `tests/unit/test_build_profiles.py`.
    """
    genes = 0
    datasets: set[str] = set()
    for gene in published:
        # `.get`, not `profiles[gene]` with a `KeyError` guard: a published
        # gene absent from `profiles` is not a caller error the way a gene
        # absent from `concordance` is (`bundles._concordance_for`'s own
        # docstring draws the same distinction) -- most published genes have
        # no profile dataset covering them at all.
        profile = profiles.get(gene)
        if profile is None or not profile["datasets"]:
            continue
        genes += 1
        datasets.update(entry["dataset"] for entry in profile["datasets"])
    return ProfileCensus(genes=genes, datasets=len(datasets))
