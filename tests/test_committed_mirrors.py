"""Invariants of the committed population mirrors that no ingest filter states.

Since the ClinGen, GenCC and Orphanet converters left for chd-atlas/feeds, the
mirrors are written only by `chd-atlas ingest`, so nothing in the atlas can
show a converter defect; these tests read the committed files and the pinned
values, which is what catches a mirror edited by hand or a pin changed without
its release (CLAUDE.md section 4.23).
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
    rows: each odds ratio equals (a / (n - a)) / (c / (m - c)) within 2% (the
    precision the paper published at), over the rows with a case carrier, a
    control carrier and a finite effect. Rows with 0 or unbounded ratios carry
    no information about the pairing and are skipped. It does not guard how the
    feeds recipe pairs the columns, only what ingest wrote.
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
    assert checked == 1041


def test_no_audain_2026_row_has_no_carrier_on_either_side() -> None:
    """An absent cell means "no qualifying variant in either group", which the gene page says.

    A 2x2 table of zeros supports no test, so the study emitted no row for it;
    that is what lets the renderer read a blank cell as "no carrier either
    side" rather than "not tested". Guards the committed rows against a row
    with 0 case and 0 control carriers, and pins the two one-sided counts: 384
    rows without a case carrier and 37 without a control carrier. It does not
    prove the study's emission rule (417 of the 417 checkable absences agreed
    with it when measured against the workbook; the rest cannot be checked).
    """
    rows = _burden_rows("PMID:42230622")
    assert rows.filter(
        (pl.col("n_case_carriers") == 0) & (pl.col("n_control_carriers") == 0)
    ).is_empty()
    assert rows.filter(pl.col("n_case_carriers") == 0).height == 384
    assert rows.filter(pl.col("n_control_carriers") == 0).height == 37


def test_audain_2026_denominators_are_its_post_qc_set() -> None:
    """3,876 cases (1,471 syndromic + 2,405 non-syndromic) against 45,082 controls.

    The abstract's 4,747 cases and 52,881 controls are recruitment figures;
    every statistic in the mirror was computed on the post-QC numbers, so a page
    printing the abstract's beside them would show a denominator no statistic
    used. Guards the committed denominators and that the two strata partition
    the cases. It does not check them against the workbook.
    """
    rows = _burden_rows("PMID:42230622")
    by_stratum = {
        stratum: sorted(
            rows.filter(pl.col("cohort_stratum") == stratum)["n_cases"].unique().to_list()
        )
        for stratum in ("all", "syndromic", "nonsyndromic")
    }
    assert by_stratum == {"all": [3876], "syndromic": [1471], "nonsyndromic": [2405]}
    assert by_stratum["syndromic"][0] + by_stratum["nonsyndromic"][0] == by_stratum["all"][0]
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
