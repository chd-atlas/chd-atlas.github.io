"""Invariants of the committed mirrors that no ingest filter states.

Holds three kinds of guard: the population mirrors (Orphanet routes, the
recorded Orphadata release), the evidence mirror `burden.tsv` (the Audain 2026
and Audain 2021 rows), and the real-corpus run of SRC007/SRC008.

Since the converters left for chd-atlas/feeds, the mirrors are written only by
`chd-atlas ingest`, so nothing in the atlas can show a recipe defect; these
tests read the committed files and the pinned values, which is what catches a
mirror edited by hand or a pin changed without its release (CLAUDE.md section
4.23).
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

import polars as pl

from chd_atlas.orphanet_scope import load_orphanet_tables, reach
from chd_atlas.tables import BURDEN, read_table
from chd_atlas.validate.sources import load_sources

REPO_ROOT = Path(__file__).resolve().parents[1]

_RECORDED_RELEASE: Final = {
    "orphanet_product4": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:57:18"),
    "orphanet_product1": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:53:50"),
}


def test_the_committed_mirrors_carry_noonan_4s_route_to_its_parent() -> None:
    """MONDO:0012547 -> OMIM:610733 -> ORPHA:648 (BTNT, Validated), measured 2026-10-24.

    Guards the committed rows the Noonan 4S route needs. It does not guard the
    ingest filters that choose the rows (tests/unit/test_ingest_filters.py does).
    """
    tables = load_orphanet_tables(REPO_ROOT)
    assert tables is not None
    assert reach("MONDO:0012547", "ORPHA:648", tables) == ("subtype", "OMIM:610733")
    assert reach("MONDO:0018997", "ORPHA:648", tables) == ("direct", None)


def test_the_recorded_orphadata_release_is_pinned() -> None:
    """The `<JDBOR>` header values `sources.json` publishes and the landing page cites (R11).

    The feeds recipes read them from the two files' `<JDBOR>` headers into each
    package's manifest, and ingest copies them into `mirrors/sources.yaml`
    (2026-06-23 release, read 2026-10-26). A re-pin that moves them, or a hand
    edit of either field, fails here: it is a change to a published value. It
    does not check the values against the pinned manifests; ingest writes them
    from the manifest.
    """
    registry, _ = load_sources(REPO_ROOT)
    recorded = {
        source_id: (registry.by_id[source_id].file_version, registry.by_id[source_id].file_date)
        for source_id in ("orphanet_product4", "orphanet_product1")
    }
    assert recorded == _RECORDED_RELEASE


def _burden_rows(study: str) -> pl.DataFrame:
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    rows = frame.filter(pl.col("study") == study)
    assert rows.height, f"no rows for {study}; did an ingest drop them?"
    return rows


def test_each_audain_2026_odds_ratio_matches_its_own_stratum() -> None:
    """PMID:42230622's odds ratio, recomputed from the four counts on its own row.

    The study's sheets give every row all four carrier counts and differ only in
    which comparison their `fet.*` columns describe, so a recipe pairing the
    wrong numerator with a stratum publishes a plausible number that is the
    syndromic-only odds ratio against the whole cohort. Guards the committed
    rows: each odds ratio equals (a / (n - a)) / (c / (m - c)) within
    max(0.02, 0.02 * expected) -- 2% relative with an absolute floor of 0.02.
    Measured on the committed rows, over those with a case carrier, a control
    carrier and a finite effect: the largest relative deviation is 0.47% and the
    largest absolute 0.173. The largest relative one is HGNC:16053 syndromic
    synonymous, published as 1.04 against a computed 1.0449, which fits rounding
    to about three significant figures; the tolerance is therefore roughly four
    times the publication precision, not equal to it. Rows with 0
    or unbounded ratios carry no information about the pairing and are skipped.
    It does not guard how the feeds recipe pairs the columns, only what ingest
    wrote.
    """
    checked = 0
    for row in _burden_rows("PMID:42230622").iter_rows(named=True):
        a, n = row["n_case_carriers"], row["n_cases"]
        c, m = row["n_control_carriers"], row["n_controls"]
        if row["effect"] is None or a == 0 or c == 0:
            continue
        expected = (a / (n - a)) / (c / (m - c))
        assert abs(expected - row["effect"]) <= max(0.02, 0.02 * expected), row
        checked += 1
    # 1,041 until PLAGL1 entered the registry (2026-11-14) with three checkable rows.
    assert checked == 1044


def test_no_audain_2026_row_has_no_carrier_on_either_side() -> None:
    """PMID:42230622 has no row with 0 case and 0 control carriers.

    Guards the committed rows against such a row, which keeps the sentence in
    `docs/data-api.md` (burden obligation 3) true, and pins the two one-sided
    counts: 387 rows without a case carrier and 37 without a control carrier
    (384 until PLAGL1 entered the registry on 2026-11-14 with three).
    It does NOT establish what an absent (stratum, consequence) cell means in
    general. The gene page and `methods.html` say an absent cell is not a null
    result and does not always mean no such variant was seen; the stronger
    reading, as a claim about every study, was withdrawn on 2026-08-05
    (`build/pages.py`, the bullet beginning "A consequence class with no row").
    Only for PMID:42230622 does `docs/data-api.md` add "no carrier on either
    side". Against the workbook, 417 of
    the 417 checkable absences were consistent with "no carrier either side"
    when measured over the study's earlier 1,192 rows (2026-08-05); that was not
    re-measured over the 1,468 rows now committed, and the rest cannot be
    checked.
    """
    rows = _burden_rows("PMID:42230622")
    assert rows.filter(
        (pl.col("n_case_carriers") == 0) & (pl.col("n_control_carriers") == 0)
    ).is_empty()
    assert rows.filter(pl.col("n_case_carriers") == 0).height == 387
    assert rows.filter(pl.col("n_control_carriers") == 0).height == 37


def test_audain_2026_denominators_are_the_ones_its_statistics_used() -> None:
    """3,876 cases (1,471 syndromic + 2,405 non-syndromic) against 45,082 controls.

    The abstract states 4,747 cases and 52,881 controls. The paper never states
    3,876, 1,471 or 2,405 and does not reconcile the gap, so the reason for the
    difference is unstated and is not asserted here (the retired converter's
    docstring says the same: `git show 8e6710a:scripts/convert_burden.py`).
    What is certain is which numbers the published odds ratios were computed
    from, and the mirror carries those; a page printing the abstract's figures
    beside them would show a denominator no statistic used. Guards the
    committed denominators. It does not check them against the workbook.
    """
    rows = _burden_rows("PMID:42230622")
    by_stratum = {
        stratum: sorted(
            rows.filter(pl.col("cohort_stratum") == stratum)["n_cases"].unique().to_list()
        )
        for stratum in ("all", "syndromic", "nonsyndromic")
    }
    assert by_stratum == {"all": [3876], "syndromic": [1471], "nonsyndromic": [2405]}
    assert rows["n_controls"].unique().to_list() == [45082]


_AUDAIN_2021_NOTE: Final = (
    "Cases combine congenital heart disease with a series of 777 sporadic "
    "thoracic aortic aneurysm probands among the study's 7,958 cases; the study "
    "does not report how many of the 777 are in this row's case set. The study "
    "reports the CHD-only contribution separately for significant genes. A "
    "carrier is an individual with a rare deletion (MAF below 0.01) overlapping "
    "the gene."
)


def test_audain_2021_rows_carry_their_note_denominators_and_s8_split() -> None:
    """PMID:34324492's note, its post-filtering denominators and its S8 CHD-only counts.

    The note renders unfolded under every carrier count of this study, so it
    must say the 777 aneurysm probands are among the cases without a count of
    them in the row's set. The denominators are the 6,746 cases and 14,024
    controls the article's methods put the permutation test on (S1's deletion
    column, 7,107 and 14,040, is the cohort as assembled). The S8 split
    (`n_chd_case_carriers`) is on the 12 registered genes S8 lists and only on
    this study's rows. Guards the committed rows; it does not guard the feeds
    recipe that locates S8's column.
    """
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    mine = _burden_rows("PMID:34324492")
    assert mine["method_note"].unique().to_list() == [_AUDAIN_2021_NOTE]
    assert mine["n_cases"].unique().to_list() == [6746]
    assert mine["n_controls"].unique().to_list() == [14024]

    split = frame.filter(pl.col("n_chd_case_carriers").is_not_null())
    assert split["study"].unique().to_list() == ["PMID:34324492"]
    assert {
        gene: (carriers, chd)
        for gene, carriers, chd in split.select(
            "gene", "n_case_carriers", "n_chd_case_carriers"
        ).iter_rows()
    } == {
        "HGNC:11592": (120, 120),  # TBX1
        "HGNC:17075": (12, 12),  # TAB2
        "HGNC:17985": (10, 10),  # ROBO4
        "HGNC:18196": (21, 21),  # SOX7
        "HGNC:21316": (13, 12),  # ANKRD11, unpublished
        "HGNC:2488": (10, 10),  # NKX2-5
        "HGNC:3227": (12, 12),  # EFNB2
        "HGNC:3327": (30, 30),  # ELN
        "HGNC:3488": (16, 16),  # ETS1
        "HGNC:4173": (29, 29),  # GATA4
        "HGNC:7976": (10, 10),  # NR2F2
        "HGNC:8569": (10, 10),  # PCSK6, unpublished
    }


def test_every_committed_mirror_file_is_listed_by_a_pinned_manifest() -> None:
    """The real repository, with no `sources.yaml` gate (unlike the validator).

    SRC007 and SRC008 only run when some source pins a package, so a repository
    that lost every pin would validate clean; this holds the real corpus
    directly. Measured at dbaa36f: 11 manifests list all 11 mirror files, and
    every source a mirror relies on is pinned except `audain2021_s8`, which no
    row's `source` names: `_used_sources` derives it from `n_chd_case_carriers`,
    a column of `audain2021_s6`'s rows (Decision 8 of the phase 3 plan). Every
    pinned manifest must load, so a corrupted one cannot make this pass by
    having SRC007 skipped.
    """
    from chd_atlas.models.package import load_manifest
    from chd_atlas.validate.runner import _used_sources
    from chd_atlas.validate.sources import validate_pins_cover_mirrors

    registry, load_issues = load_sources(REPO_ROOT)
    assert load_issues == []
    pinned = [source for source in registry.sources if source.package is not None]
    assert pinned
    for source in pinned:
        load_manifest(REPO_ROOT / "mirrors" / "packages" / f"{source.id}.json")

    issues = validate_pins_cover_mirrors(REPO_ROOT, registry, _used_sources(REPO_ROOT))

    assert [i.format() for i in issues] == []
