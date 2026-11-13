# src/chd_atlas/validate/burden.py
"""Cross-column rules for ``mirrors/burden.tsv``.

`validate_table` checks each column against its own declaration and cannot see
that `comparator` decides which *other* columns must be present. That dependency
is the whole point of the burden schema: the four kinds of published burden
evidence this atlas harmonises -- case-control, de novo enrichment against a
mutation model, de novo against unaffected siblings, and a bare case series --
differ from one another in exactly one place, which is what "expected" was.

A row whose statistic does not match its comparator is not a formatting problem.
It is an odds ratio presented as a de novo enrichment, or a p-value on a
comparison that had nothing to compare against -- a claim the published study
did not make, which the charter ranks worse than publishing nothing.

Every check here reports rather than raises, per the project's split: validators
collect issues, and `build_site` is what refuses.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import BURDEN, BURDEN_COMBINED, CHD_CASE_SPLITS, read_table
from chd_atlas.vocab import (
    EFFECT_MEASURES,
    BurdenComparator,
    CountUnit,
    EffectBound,
    VariantOrigin,
)

# What each comparator requires and what it forbids. Read as: a case-control row
# must carry both control columns and must not carry a modelled expectation,
# because it did not model one -- it observed it.
_REQUIRED: dict[BurdenComparator, tuple[str, ...]] = {
    BurdenComparator.CONTROL_COHORT: ("n_control_carriers", "n_controls", "control_cohorts"),
    BurdenComparator.MUTATION_MODEL: ("expected_count",),
    BurdenComparator.NONE: (),
}
_FORBIDDEN: dict[BurdenComparator, tuple[str, ...]] = {
    BurdenComparator.CONTROL_COHORT: ("expected_count",),
    BurdenComparator.MUTATION_MODEL: ("n_control_carriers", "n_controls", "control_cohorts"),
    BurdenComparator.NONE: (
        "n_control_carriers",
        "n_controls",
        "control_cohorts",
        "expected_count",
    ),
}

# A case series compares the count with nothing, so it can support no statistic
# at all. Listed separately from `_FORBIDDEN` because the defect is a different
# one -- not a mismatched comparator but a claim with no basis -- and a curator
# fixes it differently: by deleting the statistic, or by naming the comparator
# that actually produced it.
# Hoisted out of the per-row loop below. Rebuilt per row it cost 0.92 ms over
# the committed 1,192-row mirror against 0.03 ms hoisted -- a saving of under a
# millisecond, so this is a readability change and not a performance one.
_COMPARATORS: frozenset[str] = frozenset(item.value for item in BurdenComparator)

_STATISTIC_COLUMNS = (
    "effect",
    "effect_measure",
    "effect_bound",
    "ci_low",
    "ci_high",
    "pvalue",
    "pvalue_test",
)


def _split(value: Any) -> list[str]:
    """A `;`-joined cohort cell as a list, empty for a null or blank cell."""
    if value is None:
        return []
    return [part for part in str(value).split(";") if part.strip()]


def _blank(value: Any) -> bool:
    """True for a null cell, or a string that is empty once stripped.

    `read_table` maps `""` to null, so the string branch only fires on a value
    that is whitespace -- which renders as nothing on a page while passing every
    presence check, the same defect `_NonBlankStr` exists to catch on the
    curation side.
    """
    return value is None or (isinstance(value, str) and not value.strip())


def validate_burden(root: Path) -> list[ValidationIssue]:
    """Check every burden row's statistic against the comparator it declares.

    Returns an empty list when the table is absent -- there is no burden mirror
    before any is curated -- and also when it is unreadable or missing a column,
    both of which `validate_table` reports as TBL000/TBL001 at the same path. No
    skip warning is emitted for those two, unlike REF000/SRC000/ONT000: those
    exist because a failed load would bury one real error under hundreds of
    cascading false ones, whereas here the cause is already reported as an
    ERROR against the same file, so the build refuses either way.
    """
    path = root / "mirrors" / "burden.tsv"
    if not path.is_file():
        return []
    frame, _ = read_table(path, BURDEN)
    if frame is None:
        return []
    if any(column not in frame.columns for column in BURDEN.column_names):
        return []

    issues: list[ValidationIssue] = []

    def error(code: str, row: int, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, f"{path}:row {row}", message))

    for index, row in enumerate(frame.iter_rows(named=True)):
        # +2 for the header line and 0-based index, matching every other
        # row-level location in this project.
        line = index + 2
        raw = row["comparator"]
        # An unknown comparator is TBL004's to report; carrying on here would
        # KeyError on the tables below and abort validation of every later row.
        if raw not in _COMPARATORS:
            continue
        comparator = BurdenComparator(raw)

        for column in _REQUIRED[comparator]:
            if _blank(row[column]):
                error(
                    "BUR001",
                    line,
                    f"comparator '{comparator}' requires '{column}', which is empty",
                )
        for column in _FORBIDDEN[comparator]:
            if not _blank(row[column]):
                error(
                    "BUR001",
                    line,
                    f"comparator '{comparator}' cannot carry '{column}', which holds "
                    f"{row[column]!r}",
                )

        if comparator is BurdenComparator.NONE:
            for column in _STATISTIC_COLUMNS:
                if not _blank(row[column]):
                    error(
                        "BUR002",
                        line,
                        f"a case series compares its count with nothing, so it cannot "
                        f"report '{column}' ({row[column]!r})",
                    )

        issues.extend(_effect_issues(row, comparator, line, path))

        # A carrier count above its denominator is arithmetically impossible and
        # means a column was transcribed into the wrong slot -- the failure mode
        # a per-column range check cannot see, because both numbers are
        # individually valid non-negative integers. The third pair is the same
        # rule one level in (owner decision R2, 2026-10-26): the CHD cases
        # among a row's case carriers are a subset of them, so more is
        # impossible. It is a subset invariant and nothing more. It cannot see
        # the wrong band of S8 Table read in place of "only CHD": every band's
        # `case` is within it (measured 2026-10-27, "only CHD" plus "only TAA"
        # is "CHD+TAA", which is S6's `case`, for all 528 genes, a blank read
        # as none). What catches a misread band is
        # the chd-atlas/feeds `audain2021_s6` recipe, which locates the column from
        # S8's header rows and refuses a sheet where it cannot.
        for carriers, total in (
            ("n_case_carriers", "n_cases"),
            ("n_control_carriers", "n_controls"),
            ("n_chd_case_carriers", "n_case_carriers"),
        ):
            a, b = row[carriers], row[total]
            if a is not None and b is not None and a > b:
                error(
                    "BUR003",
                    line,
                    f"'{carriers}' ({a}) exceeds '{total}' ({b})",
                )

        # BUR024 -- a CHD-only count on a row of a study that gives none (owner
        # decision R2). `CHD_CASE_SPLITS` is the one record of which studies
        # split their case carriers and in which table; a gene page's card
        # names that table, so a count on another study's row would publish
        # under a table that study never had.
        chd = row["n_chd_case_carriers"]
        if chd is not None and row["study"] not in CHD_CASE_SPLITS:
            error(
                "BUR024",
                line,
                f"'n_chd_case_carriers' is {chd} on a row of {row['study']}, which "
                f"CHD_CASE_SPLITS records no CHD-only count for (only "
                f"{', '.join(sorted(CHD_CASE_SPLITS))})",
            )

        low, high = row["ci_low"], row["ci_high"]
        if low is not None and high is not None and low > high:
            error("BUR004", line, f"confidence interval is inverted: [{low}, {high}]")

        # BUR012 -- the interval must contain the estimate it qualifies.
        #
        # This check was deliberately left out of the first version, whose test
        # docstring recorded why: "a bracketing rule was considered and left out
        # because it was not measured against the real data". It has now been
        # measured -- 1,158 rows of `mirrors/burden.tsv` carry a finite effect
        # and both bounds, and **zero** violate bracketing -- so the reason not
        # to have it is gone.
        #
        # It is pure arithmetic with no statistical assumption behind it, and it
        # is the exact signature of a column transcribed one cell out of place:
        # `effect=5.0` with `[0.1, 0.2]` says a five-fold risk and an interval
        # entirely below 1 in the same cell, and a reader cannot tell which half
        # is the typo. `BUR004` catches only an inverted interval, which is the
        # easier half of the same defect.
        effect_value = row["effect"]
        if effect_value is not None:
            if low is not None and effect_value < low:
                error(
                    "BUR012",
                    line,
                    f"effect {effect_value} lies below its own confidence interval [{low}, {high}]",
                )
            elif high is not None and effect_value > high:
                error(
                    "BUR012",
                    line,
                    f"effect {effect_value} lies above its own confidence interval [{low}, {high}]",
                )

        # BUR013 -- no test returns a p-value of exactly zero. It is what an
        # underflowed transcription of a very small published p looks like
        # (1e-400 reads back as 0.0 through float64), and a page printing it
        # claims certainty no study can have. `p == 1` is left alone: Fisher
        # returns exactly 1 routinely, on 383 of the 1,192 committed rows.
        if row["pvalue"] == 0:
            error(
                "BUR013",
                line,
                "p-value is exactly 0; no test returns that, so this is an "
                "underflowed or mistranscribed value",
            )

        # BUR016 -- a collection cannot be its own control. `shared_cohorts`
        # unions the two columns per study to surface reuse *between* studies, so
        # reuse *within one row* is invisible to the one mechanism built to catch
        # it. Zero rows violate this today.
        overlap = set(_split(row["case_cohorts"])) & set(_split(row["control_cohorts"]))
        if overlap:
            error(
                "BUR016",
                line,
                f"cohort(s) {sorted(overlap)} appear as both cases and controls in one comparison",
            )

        # BUR018 -- a corrected p-value must name its correction, and must have
        # a raw p-value to correct. Same pairing rule as BUR005 one line below,
        # and for a sharper reason: `0.99` means "not significant after
        # family-wise correction over 11,515 permutation tests" and `0.99` with
        # no method named means nothing at all. A corrected p is also never
        # smaller than its raw p, which is arithmetic, not an assumption.
        adjusted, method = row["pvalue_adjusted"], row["pvalue_adjustment"]
        if _blank(adjusted) != _blank(method):
            error(
                "BUR018",
                line,
                f"'pvalue_adjusted' ({adjusted!r}) and 'pvalue_adjustment' "
                f"({method!r}) must be given together",
            )
        if adjusted is not None and _blank(row["pvalue"]):
            error(
                "BUR018",
                line,
                f"'pvalue_adjusted' is {adjusted} with no 'pvalue' to correct",
            )
        if adjusted is not None and row["pvalue"] is not None and adjusted < row["pvalue"]:
            error(
                "BUR018",
                line,
                f"corrected p {adjusted} is smaller than the raw p {row['pvalue']} it corrects",
            )

        # BUR019 -- a row counting de novo mutations must be counting a de novo
        # variant set. This is the *only* implication `count_unit` supports, and
        # the narrowness is deliberate: a case-control study may legitimately
        # count people (PMID:42230622) or alleles (PMID:40127276), and a de novo
        # analysis may legitimately count carriers rather than mutations, so
        # every richer rule available here would encode a convention this
        # project has not measured -- the reason `StatisticalTest` declines to
        # constrain itself against `comparator`. What cannot happen is a de novo
        # mutation count over a set that was never restricted to de novo
        # variants: the numerator would be counting something the row's own
        # `origin` says was not selected for.
        if row["count_unit"] == CountUnit.DE_NOVO_MUTATIONS and row["origin"] != (
            VariantOrigin.DE_NOVO
        ):
            error(
                "BUR019",
                line,
                f"'count_unit' is {CountUnit.DE_NOVO_MUTATIONS.value!r} but 'origin' is "
                f"{row['origin']!r}; de novo mutations cannot be counted in a variant set "
                f"that was not restricted to them",
            )

        # Neither direction alone is publishable: a p-value whose test is unnamed
        # cannot be interpreted, and a named test with no p-value is a column
        # that lost its number in transcription.
        if _blank(row["pvalue"]) != _blank(row["pvalue_test"]):
            error(
                "BUR005",
                line,
                f"'pvalue' ({row['pvalue']!r}) and 'pvalue_test' "
                f"({row['pvalue_test']!r}) must be given together",
            )

    return issues


def burden_studies_by_gene(root: Path) -> dict[str, set[str]] | None:
    """Gene -> the studies `mirrors/burden.tsv` has at least one row for.

    `{}` when the table is absent (there is no burden mirror before one is
    curated); `None` when it is unreadable or lacks `study` or `gene`, which
    `validate_table` reports as TBL000/TBL001. Read by LIT014: a literature
    finding that links a burden study must find that study's rows for its gene.
    """
    path = root / "mirrors" / "burden.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, BURDEN)
    if frame is None or not {"study", "gene"} <= set(frame.columns):
        return None
    studies: dict[str, set[str]] = {}
    for study, gene in zip(frame["study"].to_list(), frame["gene"].to_list(), strict=True):
        if study is not None and gene is not None:
            studies.setdefault(str(gene), set()).add(str(study))
    return studies


def validate_burden_references(
    root: Path,
    *,
    known_cohorts: set[str] | None,
    known_genes: set[str] | None,
    known_studies: set[str] | None,
) -> list[ValidationIssue]:
    """Check that every burden row points at something that exists.

    Kept separate from `validate_burden` because the two need different inputs
    and must fail independently: the comparator rules read nothing but the table
    and so must still run when the corpus fails to load, while these need
    registries that a failed load empties. Merging them would make one bad
    publication file silently take the comparator checks with it.

    Each registry is `None` when it could not be read, which is deliberately
    distinct from empty -- the same rule `_known_genes` follows in the runner.
    Checking thousands of rows against an empty registry would report one issue
    per row and bury the single error that caused it.

    Distinct values are checked once rather than per row, so a table with one
    unknown cohort yields one issue and not one per row that names it.
    """
    path = root / "mirrors" / "burden.tsv"
    if not path.is_file():
        return []
    frame, _ = read_table(path, BURDEN)
    if frame is None:
        return []

    issues: list[ValidationIssue] = []

    def report(code: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, str(path), message))

    if known_cohorts is not None and {"case_cohorts", "control_cohorts"} <= set(frame.columns):
        cited: set[str] = set()
        # `_split`, not a second inline split. The two disagreed on a
        # whitespace-only token: `"ddd; "` yielded `['ddd', ' ']` here and
        # `['ddd']` in `validate_burden`, so one validator reported
        # `BUR009 cohort ' ' is not in curation/cohorts.yaml` while the other
        # silently ignored it. Raised by review on #17.
        for column in ("case_cohorts", "control_cohorts"):
            for value in frame[column].to_list():
                cited.update(_split(value))
        for cohort in sorted(cited - known_cohorts):
            report("BUR009", f"cohort '{cohort}' is not in curation/cohorts.yaml")

    if known_studies is not None and "study" in frame.columns:
        studies = {str(v) for v in frame["study"].to_list() if v is not None}
        for study in sorted(studies - known_studies):
            report("BUR010", f"study {study} is not in curation/publications.yaml")

    # A burden row whose gene `mirrors/genes.tsv` does not carry is unreachable
    # by any widening of the publication gate, since every downstream population
    # derives from that registry.
    #
    # This comment claimed "`build_genes` iterates the registry" until
    # 2026-08-05. It does not -- it iterates `gene_facts(..., published=...)`,
    # the 23-gene published set -- so the stated rationale was false, and it
    # pointed away from the real gap: 127 of the 150 genes in the mirror ARE in
    # `genes.tsv` and still reach no published byte. BUR017 below reports that;
    # this rule guards only the outer boundary.
    if known_genes is not None and "gene" in frame.columns:
        genes = {str(v) for v in frame["gene"].to_list() if v is not None}
        for gene in sorted(genes - known_genes):
            report("BUR011", f"gene {gene} is not in mirrors/genes.tsv")

    # There is deliberately no rule here reporting that a burden row's gene is
    # not *published*. 1,005 of the 1,192 committed rows, covering 122 of 145
    # genes, reach no bundle and no page -- the mirror is wider than the gate so
    # that widening D21 later needs no re-mirroring -- and two review lenses
    # rightly flagged that nothing said so. But saying it here would mean
    # deriving the published set inside `validate/`, which needs
    # `build.validity.published_genes` and inverts the layering: the validator
    # would depend on the builder it exists to gate.
    #
    # The gap is recorded where each audience meets it instead -- in
    # `docs/data-api.md` for a consumer, and in CLAUDE.md's open queue for a
    # curator. What was actually wrong was the *documentation*: BUR011's comment
    # above claimed `build_genes` iterates the registry, which pointed away from
    # the gap entirely. That claim is fixed; a warning on every build for a
    # deliberate design choice is not the remedy.
    return issues


def validate_burden_combined(root: Path) -> list[ValidationIssue]:
    """Check each combined row against the arms it claims to combine.

    A combined row is a statement about two rows of `burden.tsv`, so it is only
    publishable over both of them:

    * **BUR020** -- every comparator in `combines` has a burden row with the same
      study, gene, stratum, lesion group and consequence class. Without this a
      combined p could reach a page over arms that never did.
    * **BUR021** -- `combines` names at least two distinct burden comparators.
    * **BUR023** -- neither statistic is exactly 0, the underflow BUR013 refuses
      for an arm's p.

    Reads both mirrors and nothing else, so the runner calls it outside the
    corpus branch, like `validate_burden`. An unreadable `burden.tsv` skips the
    arm check (TBL000/TBL001 already report it as an error); an *absent* one is
    one BUR020, not one per combined row.
    """
    path = root / "mirrors" / "burden_combined.tsv"
    if not path.is_file():
        return []
    frame, _ = read_table(path, BURDEN_COMBINED)
    if frame is None or any(c not in frame.columns for c in BURDEN_COMBINED.column_names):
        return []

    issues: list[ValidationIssue] = []

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    burden_path = root / "mirrors" / "burden.tsv"
    arms: set[tuple[Any, ...]] | None = None
    if not burden_path.is_file():
        if frame.height:
            error(
                "BUR020",
                str(path),
                f"{frame.height} combined row(s), but mirrors/burden.tsv is absent, so "
                f"none of the arms they combine exists",
            )
    else:
        burden, _ = read_table(burden_path, BURDEN)
        if burden is not None and all(c in burden.columns for c in BURDEN.column_names):
            arms = {
                (
                    row["study"],
                    row["gene"],
                    row["cohort_stratum"],
                    row["lesion_group"],
                    row["consequence_class"],
                    row["comparator"],
                )
                for row in burden.iter_rows(named=True)
            }

    for index, row in enumerate(frame.iter_rows(named=True)):
        location = f"{path}:row {index + 2}"
        distinct = sorted(set(_split(row["combines"])))
        unknown = [comparator for comparator in distinct if comparator not in _COMPARATORS]
        if unknown:
            error("BUR021", location, f"'combines' names {unknown}, which no burden comparator is")
        elif len(distinct) < 2:
            error(
                "BUR021",
                location,
                f"'combines' names {distinct}; a combination needs at least two distinct "
                f"arms, and one arm under a new name is a single test",
            )
        for column in ("pvalue", "local_fdr"):
            if row[column] == 0:
                error(
                    "BUR023",
                    location,
                    f"'{column}' is exactly 0; no method returns that, so this is an "
                    f"underflowed or mistranscribed value",
                )
        if arms is None:
            continue
        key = (
            row["study"],
            row["gene"],
            row["cohort_stratum"],
            row["lesion_group"],
            row["consequence_class"],
        )
        for comparator in distinct:
            if comparator in _COMPARATORS and (*key, comparator) not in arms:
                error(
                    "BUR020",
                    location,
                    f"combines a '{comparator}' arm, but mirrors/burden.tsv has no such row "
                    f"for {row['study']} {row['gene']} "
                    f"({row['cohort_stratum']}, {row['consequence_class']})",
                )
    return issues


def validate_burden_combined_references(
    root: Path,
    *,
    known_genes: set[str] | None,
    known_studies: set[str] | None,
) -> list[ValidationIssue]:
    """BUR022: a combined row citing a study or gene no registry defines.

    Separate from `validate_burden_combined` for the reason
    `validate_burden_references` is separate from `validate_burden`: these need
    registries a failed corpus load empties, and `None` means "could not be
    read", which skips the check rather than reporting every row. `source` is
    not checked here; `runner._used_sources` already feeds it to SRC003.
    """
    path = root / "mirrors" / "burden_combined.tsv"
    if not path.is_file():
        return []
    frame, _ = read_table(path, BURDEN_COMBINED)
    if frame is None:
        return []

    issues: list[ValidationIssue] = []
    for column, known, registry in (
        ("study", known_studies, "curation/publications.yaml"),
        ("gene", known_genes, "mirrors/genes.tsv"),
    ):
        if known is None or column not in frame.columns:
            continue
        values = {str(value) for value in frame[column].to_list() if value is not None}
        for value in sorted(values - known):
            issues.append(
                ValidationIssue(
                    "BUR022", Severity.ERROR, str(path), f"{column} {value} is not in {registry}"
                )
            )
    return issues


def _effect_issues(
    row: dict[str, Any], comparator: BurdenComparator, line: int, path: Path
) -> list[ValidationIssue]:
    """Rules tying `effect`, `effect_measure` and `effect_bound` together.

    Split out because these three columns constrain each other in both
    directions and inlining them would bury the comparator contract above.

    The measure is mandatory whenever there is an effect *or* a bound, because
    "unbounded" is as meaningless without it as "3.1" is: the reader has to know
    whether an unbounded odds ratio or an unbounded enrichment is being claimed.
    """
    issues: list[ValidationIssue] = []

    def error(code: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, f"{path}:row {line}", message))

    effect, measure, bound = row["effect"], row["effect_measure"], row["effect_bound"]
    has_effect = not _blank(effect)
    has_bound = not _blank(bound)

    if (has_effect or has_bound) and _blank(measure):
        error(
            "BUR006",
            "an effect size must name its measure: an odds ratio of 3.1 and a de "
            "novo enrichment of 3.1 are different claims",
        )
    if not has_effect and not has_bound and not _blank(measure):
        error("BUR006", f"'effect_measure' is {measure!r} but no effect is reported")

    if not _blank(measure) and measure not in {item.value for item in EFFECT_MEASURES[comparator]}:
        allowed = sorted(item.value for item in EFFECT_MEASURES[comparator])
        error(
            "BUR007",
            f"comparator '{comparator}' cannot produce effect_measure "
            f"'{measure}'; expected one of {allowed}",
        )

    # BUR014 -- half an interval. `ci_low` and `ci_high` arrive together or not
    # at all, for the reason BUR005 pairs `pvalue` with `pvalue_test`: one
    # without the other is a column that lost its number in transcription, and
    # `pages._effect` drops the survivor silently rather than rendering half an
    # interval. The one legitimate exception is an unbounded effect, whose
    # `ci_high` is absent by construction -- all 34 such rows in the committed
    # mirror -- so the rule is suspended exactly there.
    low, high = row["ci_low"], row["ci_high"]
    if not has_bound and (low is None) != (high is None):
        issues.append(
            ValidationIssue(
                "BUR014",
                Severity.ERROR,
                f"{path}:row {line}",
                f"a confidence interval needs both bounds; got ci_low={low!r}, "
                f"ci_high={high!r} with no effect_bound to explain the missing one",
            )
        )

    # BUR015 -- an interval belonging to nothing. BUR006 above requires a
    # measure whenever there is an effect *or* a bound; a bare interval
    # satisfies neither, so it slipped through both. A consumer reading `ci_low`
    # off the bundle would get a number with no unit and no quantity, and the
    # page renders an em dash and throws it away.
    if not has_effect and not has_bound and (low is not None or high is not None):
        issues.append(
            ValidationIssue(
                "BUR015",
                Severity.ERROR,
                f"{path}:row {line}",
                f"confidence interval [{low}, {high}] is reported with no effect and "
                f"no effect_bound, so it qualifies no quantity",
            )
        )

    # An unbounded row that also carries a number contradicts itself, and which
    # of the two a reader believes decides whether the gene looks significant.
    if has_bound and bound == EffectBound.UNBOUNDED_ABOVE.value:
        if has_effect:
            error(
                "BUR008",
                f"'effect_bound' is 'unbounded_above' but 'effect' holds {effect!r}; "
                f"an unbounded effect has no finite value",
            )
        if not _blank(row["ci_high"]):
            error(
                "BUR008",
                f"'effect_bound' is 'unbounded_above' but 'ci_high' holds {row['ci_high']!r}",
            )
    return issues
