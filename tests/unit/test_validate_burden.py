# tests/unit/test_validate_burden.py
"""The rules that stop a burden row claiming more than its study measured.

Every case below is built by overriding one field of `_ROW`, which is a real
row: GATA6 loss-of-function in the syndromic stratum of PMID 42230622,
transcribed from Supplementary Data 3. Using a real row rather than a synthetic
one matters here because several of these checks (the `Infinity` odds ratio, the
synonymous negative control) exist only because the real data has that shape.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from chd_atlas.tables import BURDEN, BURDEN_COMBINED, validate_table
from chd_atlas.validate.burden import (
    validate_burden,
    validate_burden_combined,
    validate_burden_combined_references,
    validate_burden_references,
)

_ROW: dict[str, str] = {
    "study": "PMID:42230622",
    "gene": "HGNC:4174",
    "cohort_stratum": "syndromic",
    "lesion_group": "",
    "variant_class": "snv_indel",
    "consequence_class": "lof",
    "origin": "any",
    "maf_max": "0.001",
    "count_unit": "individuals",
    "n_case_carriers": "4",
    "n_cases": "1471",
    "n_chd_case_carriers": "",
    "comparator": "control_cohort",
    "n_control_carriers": "0",
    "n_controls": "45082",
    "expected_count": "",
    "effect": "",
    "effect_measure": "odds_ratio",
    "effect_bound": "unbounded_above",
    "ci_low": "3.72",
    "ci_high": "",
    "pvalue": "1.05e-06",
    "pvalue_test": "fisher_exact",
    "pvalue_adjusted": "",
    "pvalue_adjustment": "",
    "case_cohorts": "cnchd;ddd;nottingham",
    "control_cohorts": "ukbb",
    "method_note": "",
    "source": "audain2026_sd3",
}

# A de novo row, whose comparator is a mutation model rather than a control
# cohort. Kept beside `_ROW` because half the point of this schema is that both
# shapes live in one table, and a rule that only ever sees case-control rows
# would not be tested against the branch that admits neither control column.
_DE_NOVO: dict[str, str] = {
    **_ROW,
    "origin": "de_novo",
    "comparator": "mutation_model",
    "n_control_carriers": "",
    "n_controls": "",
    "control_cohorts": "",
    "expected_count": "0.42",
    "effect": "9.5",
    "effect_measure": "enrichment_ratio",
    "effect_bound": "",
    "ci_high": "24.1",
    "pvalue_test": "poisson",
}


# A case series: carriers, a denominator, and nothing to compare them with.
# Module-level rather than built inside one test because two separate rules
# constrain it -- it may carry no comparator column (BUR001) and no statistic
# (BUR002) -- and a mutation dropping either was measured to survive while this
# lived in only one of them.
_SERIES: dict[str, str] = {
    **_ROW,
    "comparator": "none",
    "n_control_carriers": "",
    "n_controls": "",
    "control_cohorts": "",
    "effect": "",
    "effect_measure": "",
    "effect_bound": "",
    "ci_low": "",
    "ci_high": "",
    "pvalue": "",
    "pvalue_test": "",
}


def _write(tmp_path: Path, *rows: dict[str, str]) -> Path:
    """Write `rows` as `mirrors/burden.tsv` under a fresh root."""
    mirrors = tmp_path / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    header = "\t".join(BURDEN.column_names)
    body = ["\t".join(row[name] for name in BURDEN.column_names) for row in rows]
    (mirrors / "burden.tsv").write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    return tmp_path


def _codes(tmp_path: Path, *rows: dict[str, str]) -> list[str]:
    return sorted(issue.code for issue in validate_burden(_write(tmp_path, *rows)))


def test_a_real_case_control_row_and_a_real_de_novo_row_both_pass(tmp_path: Path) -> None:
    """The smoke test every check below rests on.

    Without it, a rule that fired unconditionally would still make each negative
    case go green, and the whole suite would agree the validator worked.
    """
    assert validate_burden(_write(tmp_path, _ROW, _DE_NOVO)) == []


@pytest.mark.parametrize(
    ("row", "field", "value"),
    [
        # Case-control: both control columns and the control cohort are required.
        pytest.param(_ROW, "n_control_carriers", "", id="case_control-needs-carriers"),
        pytest.param(_ROW, "n_controls", "", id="case_control-needs-denominator"),
        pytest.param(_ROW, "control_cohorts", "", id="case_control-needs-cohort"),
        # ...and a modelled expectation is not something it has.
        pytest.param(_ROW, "expected_count", "0.42", id="case_control-cannot-model"),
        # De novo: the mirror image. An expectation is required, and a control
        # count is a column it cannot have observed.
        pytest.param(_DE_NOVO, "expected_count", "", id="mutation_model-needs-expected"),
        pytest.param(_DE_NOVO, "n_controls", "45082", id="mutation_model-has-no-controls"),
        pytest.param(
            _DE_NOVO, "control_cohorts", "ukbb", id="mutation_model-has-no-controls-named"
        ),
        # A case series: none of the three. `expected_count` is the one a
        # mutation was measured to slip past when this case lived only in
        # `test_a_case_series_cannot_report_any_statistic` -- a row claiming to
        # compare against nothing while carrying a modelled expectation
        # validated clean, which is a comparator the study did use, unrecorded.
        pytest.param(_SERIES, "expected_count", "0.42", id="case_series-cannot-model"),
        pytest.param(_SERIES, "n_controls", "45082", id="case_series-has-no-controls"),
    ],
)
def test_a_row_must_carry_exactly_the_columns_its_comparator_implies(
    tmp_path: Path, row: dict[str, str], field: str, value: str
) -> None:
    """BUR001, in both directions: a missing requirement and a stray field.

    Both are the same defect seen from either side -- the row's declared
    comparator does not match the numbers it carries -- so both report one code.
    A curator fixes either by deciding which of the two is the truth.
    """
    assert "BUR001" in _codes(tmp_path, {**row, field: value})


@pytest.mark.parametrize(
    "field",
    ["effect", "effect_measure", "effect_bound", "ci_low", "ci_high", "pvalue", "pvalue_test"],
)
def test_a_case_series_cannot_report_any_statistic(tmp_path: Path, field: str) -> None:
    """BUR002. A case series counted carriers and compared them with nothing.

    Every one of these columns is therefore a claim with no basis, and the
    parametrisation is exhaustive over `_STATISTIC_COLUMNS` on purpose: dropping
    one from that tuple is a mutation this test must kill, and it can only do
    that if it names them all.
    """
    assert validate_burden(_write(tmp_path, _SERIES)) == []
    assert "BUR002" in _codes(tmp_path, {**_SERIES, field: "1.5"})


def test_a_mutation_model_cannot_report_an_odds_ratio(tmp_path: Path) -> None:
    """BUR007 -- the guard the whole `effect`/`effect_measure` pair exists for.

    A single `effect` column holding both odds ratios and de novo enrichments is
    what makes this schema extend to a fifth study without a migration, and it
    is also the one place two incomparable quantities could silently merge. A
    mutation model has no control odds, so it cannot yield an odds ratio; a
    control cohort observes its expectation rather than modelling one, so it
    cannot yield an enrichment over expectation.
    """
    assert "BUR007" in _codes(tmp_path, {**_DE_NOVO, "effect_measure": "odds_ratio"})
    assert "BUR007" in _codes(tmp_path, {**_ROW, "effect_measure": "enrichment_ratio"})


def test_an_effect_and_its_measure_are_useless_without_each_other(tmp_path: Path) -> None:
    """BUR006, both directions.

    A bare 3.1 does not say whether it is an odds ratio or an enrichment, and a
    measure with nothing to measure is a number lost in transcription. The
    unbounded case is asserted too, because "unbounded" needs its measure for
    exactly the same reason a finite value does -- and it is the branch a naive
    `if effect is None: skip` would miss.
    """
    assert "BUR006" in _codes(tmp_path, {**_DE_NOVO, "effect_measure": ""})
    assert "BUR006" in _codes(tmp_path, {**_ROW, "effect_measure": ""})
    assert "BUR006" in _codes(tmp_path, {**_DE_NOVO, "effect": "", "ci_low": "", "ci_high": ""})


def test_an_unbounded_effect_may_not_also_carry_a_number(tmp_path: Path) -> None:
    """BUR008. The flag and a finite value contradict each other.

    Which one a reader believes decides whether the gene looks significant, so
    this cannot be left to whichever the renderer happens to reach for.
    """
    assert "BUR008" in _codes(tmp_path, {**_ROW, "effect": "12.4"})
    assert "BUR008" in _codes(tmp_path, {**_ROW, "ci_high": "99.0"})


def test_a_carrier_count_above_its_denominator_is_reported(tmp_path: Path) -> None:
    """BUR003. Both numbers are individually valid non-negative integers, so no
    per-column range check can see this -- it is what a column transcribed into
    the wrong slot looks like, which is the single likeliest converter bug.

    The CHD cases among a row's case carriers are a subset of them (owner
    decision R2, 2026-10-26), so more is the same defect one level in, and as
    many is the ordinary case: 11 of the 12 rows carrying the column when it
    was transcribed (2026-10-27; ANKRD11's 12 of 13 is the twelfth). Fewer
    than none is the schema's to refuse, as TBL006. These are PMID:34324492's
    rows, the one study whose rows carry the column.
    """
    split = {**_ROW, "study": "PMID:34324492"}
    assert "BUR003" in _codes(tmp_path, {**_ROW, "n_case_carriers": "1472"})
    assert "BUR003" in _codes(tmp_path, {**_ROW, "n_control_carriers": "45083"})
    assert "BUR003" in _codes(tmp_path, {**split, "n_chd_case_carriers": "5"})
    assert _codes(tmp_path, {**split, "n_chd_case_carriers": "4"}) == []
    negative = _write(tmp_path, {**split, "n_chd_case_carriers": "-1"}) / "mirrors" / "burden.tsv"
    assert [issue.code for issue in validate_table(negative, BURDEN)] == ["TBL006"]


def test_a_chd_only_count_on_a_study_that_gives_none_is_reported(tmp_path: Path) -> None:
    """BUR024. A gene page's card names the table the count came from, read
    from `CHD_CASE_SPLITS` by the row's study. Before this code, a count on a
    PMID:40127276 row validated clean and its card published "Sierant MC et al.
    2025: ... the study's S8 Table counts ..." (review, 2026-10-27), a table
    that study does not have. Reported once per row, like every row rule here;
    the study `CHD_CASE_SPLITS` lists passes (the BUR003 test above).
    """
    other = {**_ROW, "n_chd_case_carriers": "4"}
    assert _codes(tmp_path, other, {**other, "gene": "HGNC:4173"}) == ["BUR024", "BUR024"]


def test_an_inverted_confidence_interval_is_reported(tmp_path: Path) -> None:
    """BUR004. Measured 2026-08-04 against all three strata of Supplementary
    Data 3 (138,609 rows): zero rows have `ci_lower > ci_upper`, so this fires
    on a transcription error rather than on anything the study published.

    It deliberately does *not* check that the interval brackets the point
    estimate. That sentence is load-bearing: a bracketing rule was considered
    and left out because it was not measured against the real data, and
    asserting an invariant nobody verified is how this project's docstrings have
    been wrong before.

    The inversion below is 0.5 wide, not the 5.9 it was first written with. An
    off-by-one mutant (`low > high + 1`) survived the wider gap: the check still
    fired, so the test still passed, while a genuinely inverted [24.6, 24.1]
    would have validated clean. The boundary is what needs pinning, not the
    direction.
    """
    assert "BUR004" in _codes(tmp_path, {**_DE_NOVO, "ci_low": "24.6"})


def test_a_pvalue_and_the_test_that_produced_it_must_arrive_together(tmp_path: Path) -> None:
    """BUR005. A p-value whose test is unnamed cannot be interpreted, and a
    named test with no p-value lost its number somewhere.
    """
    assert "BUR005" in _codes(tmp_path, {**_ROW, "pvalue": ""})
    assert "BUR005" in _codes(tmp_path, {**_ROW, "pvalue_test": ""})


def test_a_whitespace_only_cell_counts_as_empty(tmp_path: Path) -> None:
    """`read_table` maps `""` to null but leaves `" "` a string, which passes
    every presence check while rendering as nothing on a page. Same defect
    `_NonBlankStr` exists to catch on the curation side.
    """
    assert "BUR001" in _codes(tmp_path, {**_ROW, "control_cohorts": "   "})


def test_an_absent_or_unreadable_table_reports_nothing_and_does_not_raise(
    tmp_path: Path,
) -> None:
    """Three distinct absences, none of which is this validator's to report.

    No burden mirror exists before any is curated. An unreadable file is
    TBL000's and a missing column is TBL001's, both ERRORs against the same
    path, so the build refuses either way and a second report would only say it
    twice. Raising instead would abort validation of every other table.
    """
    assert validate_burden(tmp_path) == []

    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "burden.tsv").write_bytes(b"\xff\xfe not utf-8")
    assert validate_burden(tmp_path) == []

    (tmp_path / "mirrors" / "burden.tsv").write_text("gene\tstudy\nHGNC:4174\tPMID:1\n")
    assert validate_burden(tmp_path) == []


_KNOWN = {
    "known_cohorts": {"cnchd", "ddd", "nottingham", "ukbb"},
    "known_genes": {"HGNC:4174", "HGNC:11603"},
    "known_studies": {"PMID:42230622"},
}


@pytest.mark.parametrize(
    ("field", "value", "code"),
    [
        pytest.param("case_cohorts", "cnchd;pcgc", "BUR009", id="unknown-case-cohort"),
        pytest.param("control_cohorts", "gnomad", "BUR009", id="unknown-control-cohort"),
        pytest.param("study", "PMID:99999999", "BUR010", id="unknown-study"),
        pytest.param("gene", "HGNC:99999", "BUR011", id="unknown-gene"),
    ],
)
def test_a_row_pointing_at_something_no_registry_defines_is_reported(
    tmp_path: Path, field: str, value: str, code: str
) -> None:
    """Each dangling reference has its own code because each is fixed elsewhere.

    The gene case is the one on this project's characteristic failure: a burden
    row for a gene `mirrors/genes.tsv` does not carry can never reach a page,
    however the publication gate later widens, because `build_genes` iterates the
    registry. It costs nothing to hold and reaches no reader.

    The `case_cohorts` case also proves the `;` split: `cnchd` resolves and
    `pcgc` does not, so a check that compared the whole cell against the registry
    would report both or neither.
    """
    root = _write(tmp_path, {**_ROW, field: value})
    issues = validate_burden_references(root, **_KNOWN)  # type: ignore[arg-type]

    assert [issue.code for issue in issues] == [code]


def test_a_registry_that_could_not_be_read_skips_its_own_check_only(tmp_path: Path) -> None:
    """`None` means "could not be read" and is deliberately not the empty set.

    Against an empty registry every row reports a dangling reference -- hundreds
    of issues naming the symptom, none naming the cause, which is the cascade
    `_known_genes` and `_mirrored_validity` return `None` to prevent. The second
    assertion is what stops a lazy `if not known_genes: return []`: an unreadable
    gene registry must not also silence the cohort and study checks.
    """
    root = _write(tmp_path, {**_ROW, "gene": "HGNC:99999", "study": "PMID:99999999"})

    assert (
        validate_burden_references(root, known_cohorts=None, known_genes=None, known_studies=None)
        == []
    )

    codes = [
        issue.code
        for issue in validate_burden_references(
            root,
            known_cohorts=_KNOWN["known_cohorts"],  # type: ignore[arg-type]
            known_genes=None,
            known_studies=_KNOWN["known_studies"],  # type: ignore[arg-type]
        )
    ]
    assert codes == ["BUR010"]


def test_one_unknown_cohort_is_reported_once_however_many_rows_cite_it(
    tmp_path: Path,
) -> None:
    """The mirror will carry thousands of rows drawn from one supplement, so a
    per-row report would bury every other issue in the file under one repeated
    line. Same reasoning as `validate_mirror_references`' distinct-pair check.
    """
    rows = [{**_ROW, "gene": gene, "case_cohorts": "pcgc"} for gene in ("HGNC:4174", "HGNC:11603")]
    issues = validate_burden_references(_write(tmp_path, *rows), **_KNOWN)  # type: ignore[arg-type]

    assert len(issues) == 1
    assert "pcgc" in issues[0].message


def test_an_unknown_comparator_is_left_to_the_column_check(tmp_path: Path) -> None:
    """TBL004 reports a value outside the allowed set. Reporting it again here
    would say the same thing twice; worse, indexing `_REQUIRED` with it would
    raise `KeyError` and abort validation of every row below it in the file.

    The second row proves the loop keeps going: without the `continue`, the
    clean row after the bad one is never reached.
    """
    codes = _codes(tmp_path, {**_ROW, "comparator": "sibling_rate"}, {**_ROW, "gene": "HGNC:11603"})
    assert [code for code in codes if code.startswith("BUR")] == []


@pytest.mark.parametrize(
    ("mutation", "code"),
    [
        # BUR012 -- the interval must contain the estimate it qualifies. The
        # signature of a column transcribed one cell out of place: "OR 5.0" and
        # an interval entirely below 1 say opposite clinical things in one cell,
        # and a reader cannot tell which half is the typo. BUR004 catches only
        # the inverted interval, which is the easier half of the same defect.
        pytest.param(
            {"effect": "5.0", "effect_bound": "", "ci_low": "0.1", "ci_high": "0.2"},
            "BUR012",
            id="effect-above-its-interval",
        ),
        pytest.param(
            {"effect": "0.05", "effect_bound": "", "ci_low": "0.1", "ci_high": "0.2"},
            "BUR012",
            id="effect-below-its-interval",
        ),
        # BUR013 -- no test returns exactly 0. This is what 1e-400 looks like
        # after float64 underflow, and a page printing it claims certainty no
        # study can have.
        pytest.param({"pvalue": "0"}, "BUR013", id="p-exactly-zero"),
        # BUR014 -- half an interval, with no effect_bound to explain it. Same
        # rule and rationale as BUR005's pvalue/pvalue_test pairing.
        pytest.param(
            {"effect": "5.0", "effect_bound": "", "ci_low": "0.5", "ci_high": ""},
            "BUR014",
            id="lower-bound-only",
        ),
        pytest.param(
            {"effect": "5.0", "effect_bound": "", "ci_low": "", "ci_high": "9.9"},
            "BUR014",
            id="upper-bound-only",
        ),
        # BUR015 -- an interval qualifying nothing. BUR006 requires a measure
        # when there is an effect *or* a bound; a bare interval satisfies
        # neither, so it slipped past both.
        pytest.param(
            {
                "effect": "",
                "effect_measure": "",
                "effect_bound": "",
                "ci_low": "0.5",
                "ci_high": "9.9",
            },
            "BUR015",
            id="interval-with-no-quantity",
        ),
        # BUR016 -- a collection cannot be its own control. `shared_cohorts`
        # unions the two columns to surface reuse *between* studies, so reuse
        # *within one row* is invisible to the one mechanism built to catch it.
        pytest.param({"case_cohorts": "cnchd;ukbb"}, "BUR016", id="cohort-on-both-sides"),
    ],
)
def test_a_statistic_that_contradicts_itself_is_refused(
    tmp_path: Path, mutation: dict[str, str], code: str
) -> None:
    """Six rules added after adversarial review, each measured before adding.

    None fires on the committed mirror: 0 of 1,158 rows with a finite effect and
    both bounds violate bracketing, 0 of 1,192 report p exactly 0 or carry half
    an interval outside the 34 unbounded rows, 0 carry an interval with no
    quantity, and 0 name a cohort on both sides. So each is a guard on data that
    is clean today, added because the defect it catches is silent and the fix a
    curator would apply is different in each case.
    """
    assert code in _codes(tmp_path, {**_ROW, **mutation})


def test_an_unbounded_row_is_the_one_interval_allowed_to_be_half_present(
    tmp_path: Path,
) -> None:
    """BUR014's single exception, asserted so it cannot be tightened away.

    All 34 unbounded rows in the committed mirror have `ci_low` and no
    `ci_high`. A rule requiring both bounds unconditionally would reject every
    one of them -- which is to say, the strongest results in the study.
    """
    assert validate_burden(_write(tmp_path, _ROW)) == []
    assert _ROW["effect_bound"] == "unbounded_above"
    assert _ROW["ci_high"] == ""


def test_both_validators_split_a_cohort_cell_the_same_way(tmp_path: Path) -> None:
    """One column, two parsers, and they disagreed on a whitespace token.

    Raised by review on #17. `validate_burden_references` split with
    `[p for p in value.split(";") if p]` while `validate_burden` used `_split`,
    which also strips. So `"cnchd; "` yielded `['cnchd', ' ']` in one and
    `['cnchd']` in the other, and a curator got

        BUR009 cohort ' ' is not in curation/cohorts.yaml

    from the reference check while BUR016 saw nothing there at all -- an error
    naming a cohort that does not exist, pointing at a registry that could never
    contain it.

    A guard on a bypassed gate: `COHORT_LIST_PATTERN` admits no space, so
    TBL005 refuses this cell first. It is asserted anyway because the defect is
    two parsers of one column drifting apart, which no gate prevents.
    """
    issues = validate_burden_references(
        _write(tmp_path, {**_ROW, "case_cohorts": "cnchd; "}),
        **_KNOWN,  # type: ignore[arg-type]
    )

    assert [issue.message for issue in issues] == []


@pytest.mark.parametrize(
    "mutation",
    [
        pytest.param({"pvalue_adjusted": "0.9"}, id="corrected-p-with-no-method"),
        pytest.param({"pvalue_adjustment": "bonferroni"}, id="method-with-no-corrected-p"),
        pytest.param(
            {
                "pvalue": "",
                "pvalue_test": "",
                "pvalue_adjusted": "0.9",
                "pvalue_adjustment": "bonferroni",
            },
            id="corrected-p-with-nothing-to-correct",
        ),
        pytest.param(
            {"pvalue": "0.5", "pvalue_adjusted": "0.01", "pvalue_adjustment": "bonferroni"},
            id="correction-made-the-p-smaller",
        ),
    ],
)
def test_a_corrected_p_value_must_name_its_method_and_exceed_its_raw_p(
    tmp_path: Path, mutation: dict[str, str]
) -> None:
    """BUR018. `0.99` with no method named means nothing at all.

    The pairing rule is BUR005's, applied one column over, and for a sharper
    reason: a family-wise permutation correction over 11,515 tests and a
    Bonferroni factor are different claims, and the number alone does not say
    which. A correction is also never *smaller* than the p it corrects, which is
    arithmetic rather than an assumption -- the signature of the two columns
    being transposed.

    Zero of the 1,295 committed rows violate any of these.
    """
    assert "BUR018" in _codes(tmp_path, {**_ROW, **mutation})


def test_de_novo_mutations_cannot_be_counted_in_a_set_that_is_not_de_novo(
    tmp_path: Path,
) -> None:
    """BUR019, and the *only* rule `count_unit` supports.

    A numerator counting de novo mutations over a variant set the row's own
    `origin` says was never restricted to de novo variants is counting something
    that was not selected for. Every other pairing this column suggests is a
    convention rather than an implication -- a case-control study may count
    people (PMID:42230622) or alleles (PMID:40127276), and a de novo analysis may
    count carriers rather than mutations -- so asserting more here would encode a
    rule this project has not measured, the reason `StatisticalTest` declines to
    constrain itself against `comparator`.

    The negative half is what makes that narrowness real: the same unit with a
    de novo origin is clean, and so is every other unit against every origin.
    """
    offending = {**_ROW, "count_unit": "de_novo_mutations", "origin": "any"}
    assert "BUR019" in _codes(tmp_path, offending)

    assert "BUR019" not in _codes(tmp_path, {**offending, "origin": "de_novo"})
    assert "BUR019" not in _codes(tmp_path, {**_ROW, "count_unit": "alleles"})
    assert "BUR019" not in _codes(tmp_path, {**_ROW, "count_unit": "individuals"})


# PMID:40127276's two arms for one gene and consequence, and the study's own
# combination of them. Only the fields BUR020 matches on differ from `_ROW`.
_ARM_KEY = {
    "study": "PMID:40127276",
    "gene": "HGNC:7133",
    "cohort_stratum": "all",
    "consequence_class": "lof",
}
_CASE_CONTROL_ARM = {
    **_ROW,
    **_ARM_KEY,
    "case_cohorts": "pcgc",
    "control_cohorts": "gnomad_controls",
}
_DE_NOVO_ARM = {**_DE_NOVO, **_ARM_KEY, "case_cohorts": "pcgc"}
_COMBINED: dict[str, str] = {
    **_ARM_KEY,
    "lesion_group": "",
    "combines": "mutation_model;control_cohort",
    "method": "joint_local_fdr",
    "pvalue": "1.6582e-38",
    "local_fdr": "2.6973e-34",
    "method_note": "",
    "source": "sierant2025_ds6",
}


def _write_combined(root: Path, *rows: dict[str, str]) -> Path:
    mirrors = root / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    header = "\t".join(BURDEN_COMBINED.column_names)
    body = ["\t".join(row[name] for name in BURDEN_COMBINED.column_names) for row in rows]
    path = mirrors / "burden_combined.tsv"
    path.write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    return root


def _combined_codes(
    tmp_path: Path, arms: list[dict[str, str]], *combined: dict[str, str]
) -> list[str]:
    root = _write_combined(_write(tmp_path, *arms) if arms else tmp_path, *combined)
    return sorted(issue.code for issue in validate_burden_combined(root))


def test_a_combination_over_both_of_its_arms_passes(tmp_path: Path) -> None:
    """The smoke test the cases below rest on (see the first test in this file)."""
    assert _combined_codes(tmp_path, [_DE_NOVO_ARM, _CASE_CONTROL_ARM], _COMBINED) == []


@pytest.mark.parametrize(
    "arms",
    [
        pytest.param([_CASE_CONTROL_ARM], id="de-novo-arm-missing"),
        pytest.param([_DE_NOVO_ARM], id="case-control-arm-missing"),
        # The arm exists for the same gene and stratum but another consequence:
        # a check that matched on (study, gene, stratum) alone would pass this.
        pytest.param(
            [{**_DE_NOVO_ARM, "consequence_class": "damaging"}, _CASE_CONTROL_ARM],
            id="arm-only-for-another-consequence",
        ),
    ],
)
def test_a_combination_whose_arm_is_not_published_is_refused(
    tmp_path: Path, arms: list[dict[str, str]]
) -> None:
    """BUR020: a combined p over an arm no reader can see is a claim with no rows under it."""
    assert _combined_codes(tmp_path, arms, _COMBINED) == ["BUR020"]


def test_a_combination_with_no_burden_mirror_at_all_is_reported_once(tmp_path: Path) -> None:
    """Absent `burden.tsv`: one issue naming the cause, not one per combined row."""
    codes = _combined_codes(tmp_path, [], _COMBINED, {**_COMBINED, "consequence_class": "damaging"})
    assert codes == ["BUR020"]


@pytest.mark.parametrize(
    "combines",
    [
        pytest.param("mutation_model", id="one-arm"),
        pytest.param("mutation_model;mutation_model", id="one-arm-twice"),
        pytest.param("mutation_model;case_control", id="unknown-comparator"),
    ],
)
def test_a_combination_must_name_at_least_two_distinct_known_arms(
    tmp_path: Path, combines: str
) -> None:
    """BUR021: one arm relabelled as a combination is a single test with a new name."""
    arms = [_DE_NOVO_ARM, _CASE_CONTROL_ARM]
    assert _combined_codes(tmp_path, arms, {**_COMBINED, "combines": combines}) == ["BUR021"]


@pytest.mark.parametrize("column", ["pvalue", "local_fdr"])
def test_a_combined_statistic_of_exactly_zero_is_refused(tmp_path: Path, column: str) -> None:
    """BUR023: an underflowed transcription, as BUR013 is for an arm's p."""
    arms = [_DE_NOVO_ARM, _CASE_CONTROL_ARM]
    assert _combined_codes(tmp_path, arms, {**_COMBINED, column: "0"}) == ["BUR023"]


@pytest.mark.parametrize(
    ("field", "value"),
    [
        pytest.param("study", "PMID:99999999", id="unknown-study"),
        pytest.param("gene", "HGNC:99999", id="unknown-gene"),
    ],
)
def test_a_combined_row_citing_nothing_a_registry_defines_is_reported(
    tmp_path: Path, field: str, value: str
) -> None:
    """BUR022, with the `None`-means-unreadable rule the arm references follow."""
    root = _write_combined(tmp_path, {**_COMBINED, field: value})
    issues = validate_burden_combined_references(
        root, known_genes={"HGNC:7133"}, known_studies={"PMID:40127276"}
    )
    assert [issue.code for issue in issues] == ["BUR022"]
    assert validate_burden_combined_references(root, known_genes=None, known_studies=None) == []


# --- committed-mirror invariants the deleted Sierant converter suite carried ---


def _mirror_rows(name: str) -> list[dict[str, str]]:
    path = Path(__file__).parent.parent.parent / "mirrors" / f"{name}.tsv"
    with path.open(encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _composite_gaps(rows: list[dict[str, str]]) -> list[tuple[str, str, int, int]]:
    """(gene, comparator, damaging, lof + missense_damaging) wherever the sum fails."""
    counts = {
        (r["gene"], r["consequence_class"], r["comparator"]): int(r["n_case_carriers"])
        for r in rows
        if r["study"] == "PMID:40127276"
    }
    return [
        (
            gene,
            comparator,
            n,
            counts[(gene, "lof", comparator)] + counts[(gene, "missense_damaging", comparator)],
        )
        for (gene, cls, comparator), n in sorted(counts.items())
        if cls == "damaging"
        and n != counts[(gene, "lof", comparator)] + counts[(gene, "missense_damaging", comparator)]
    ]


def _pair_gaps(arms: list[dict[str, str]], combined: list[dict[str, str]]) -> list[str]:
    """Problems with "one combined row per arm pair" for PMID:40127276, in both directions."""
    pairs = sorted(
        {
            (r["gene"], r["cohort_stratum"], r["lesion_group"], r["consequence_class"])
            for r in arms
            if r["study"] == "PMID:40127276"
        }
    )
    mine = [
        (r["gene"], r["cohort_stratum"], r["lesion_group"], r["consequence_class"])
        for r in combined
        if r["study"] == "PMID:40127276"
    ]
    problems = [f"no combined row for {p}" for p in pairs if p not in mine]
    problems += [f"combined row over no arm pair {p}" for p in sorted(set(mine) - set(pairs))]
    problems += [
        f"{p} combined {mine.count(p)} times" for p in sorted(set(mine)) if mine.count(p) > 1
    ]
    return problems


def test_the_committed_sierant_composite_is_the_sum_of_its_components() -> None:
    """`damaging` is `lof` + `missense_damaging` on case carriers, for each gene and arm.

    Catches a stale or mis-transcribed mirror; the arithmetic holds only if the
    three blocks were read from the bands of the sheet they were meant to be.
    Measured 2026-11-11: no gap at f094c28 (41 genes) and after the Sierant ingest (42).
    """
    rows = _mirror_rows("burden")
    assert len({r["gene"] for r in rows if r["study"] == "PMID:40127276"}) > 1
    assert _composite_gaps(rows) == []


def test_the_committed_combined_mirror_has_one_row_per_arm_pair() -> None:
    """Every Sierant arm pair has exactly one combined row and no combined row is orphaned.

    BUR020 checks only that each combined row has its arms; this is the other
    direction, plus no repeat. The arm pair is the arm rows' key without the
    comparator, which the two arms of a pair do not share.
    """
    arms = _mirror_rows("burden")
    combined = _mirror_rows("burden_combined")
    assert [r for r in combined if r["study"] == "PMID:40127276"]
    assert _pair_gaps(arms, combined) == []
