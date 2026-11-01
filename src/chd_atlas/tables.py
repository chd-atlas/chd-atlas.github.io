# src/chd_atlas/tables.py
"""Schema definition and validation for the sorted TSV tables under ``mirrors/``.

Polars rather than pandas: pandas coerces a nullable integer column to float,
which would silently turn a genomic coordinate of 12345 into 12345.0.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

# Imported from polars.datatypes, NOT as pl.DataTypeClass: the top-level alias was
# deprecated in polars 1.0 and survives only as a runtime shim, so `pl.DataTypeClass`
# fails `mypy --strict` with "Name is not defined" on every version we support.
from polars.datatypes import DataTypeClass
from polars.exceptions import PolarsError

from chd_atlas.duplicates import duplicates
from chd_atlas.fs import list_dir
from chd_atlas.identifiers import (
    COHORT_LIST_PATTERN,
    HGNC_PATTERN,
    HPO_PATTERN,
    MODIFICATION_PATTERN,
    MONDO_PATTERN,
    OMIM_PATTERN,
    ORPHA_PATTERN,
    PMID_PATTERN,
    SEQUENCE_ONTOLOGY_PATTERN,
    UNIPROT_PATTERN,
    XREF_TARGET_PATTERN,
)
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.vocab import (
    BurdenComparator,
    CohortStratum,
    CombinationMethod,
    ConsequenceClass,
    CountUnit,
    EffectBound,
    EffectMeasure,
    LesionGroup,
    PvalueAdjustment,
    SourceTier,
    StatisticalTest,
    VariantClass,
    VariantOrigin,
)


@dataclass(frozen=True)
class Column:
    name: str
    # DataTypeClass, not DataType: we always pass the class (pl.Int64), never an instance.
    dtype: DataTypeClass
    nullable: bool = False
    allowed: frozenset[str] | None = None
    pattern: str | None = None
    minimum: float | None = None
    maximum: float | None = None


@dataclass(frozen=True)
class TableSchema:
    name: str
    columns: tuple[Column, ...]
    sort_key: tuple[str, ...]

    @property
    def column_names(self) -> tuple[str, ...]:
        return tuple(column.name for column in self.columns)

    def polars_overrides(self) -> dict[str, DataTypeClass]:
        return {column.name: column.dtype for column in self.columns}


def cell_text(value: object) -> str | None:
    """One mirror cell as text, with blank and whitespace-only treated as absent.

    `read_table` maps the empty string to null, so TBL003 catches `symbol` left
    empty — and catches nothing when the cell holds a single space. A generated
    mirror produces exactly that, and the space survived every check to be
    published as the gene's browse label and its search result: a row rendering
    as nothing, on a build reporting zero errors, against a `docs/data-api.md`
    that promises `symbol` is "always a non-empty string you can render and
    search". Found by adversarial review of #2.

    Stripped rather than reported, because the fallback the caller already has —
    keep the HGNC id — is better than either publishing the blank or refusing to
    build over one stray character.

    Here rather than in `build/runner.py`, where it was `_cell` until review
    (2026-10-27): `build.identity.read_gene_identities` reads the registry's
    identity columns for `chd-atlas validate` and for the build alike, and the
    validator cannot import the build's runner, which imports it.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def read_table(
    path: Path, schema: TableSchema
) -> tuple[pl.DataFrame | None, list[ValidationIssue]]:
    """Read a TSV with declared dtypes, treating the empty string as null.

    Returns ``(None, issues)`` when the file cannot be read at all — a non-UTF-8
    byte, a zero-length file, or a missing path. Polars raises for the first two
    and a missing path raises ``FileNotFoundError``; an uncaught raise here would
    abort validation of every other shard in the repository rather than
    reporting one unreadable file.
    """
    try:
        # `scan_csv(...).collect_schema()` rather than `read_csv(path, n_rows=0)`:
        # both discover only the header in the success case (measured identical
        # column lists), but a lazy scan does not validate the file's body, so it
        # is ~2.5x faster to fail past on a large file (measured on the 30,410-row
        # GenCC mirror: 17ms vs 45ms, averaged over 20 warmed calls).
        #
        # This does change *which* call raises on a body-only defect (a bad byte
        # past the header, valid UTF-8 in the header itself): the lazy scan does
        # not see it, so the second, real `read_csv` below raises instead. Measured
        # identical for all three unreadable-file cases this module is tested
        # against -- empty file, non-UTF-8 byte, missing path -- same exception
        # type and same `str(exc)` in both the old eager-first-read order and this
        # one, so the `TBL000` message is byte-identical either way.
        present = pl.scan_csv(path, separator="\t").collect_schema().names()
        overrides = {
            name: dtype for name, dtype in schema.polars_overrides().items() if name in present
        }
        frame = pl.read_csv(
            path,
            separator="\t",
            schema_overrides=overrides,
            null_values=[""],
            infer_schema_length=0 if not overrides else None,
        )
    except (PolarsError, OSError) as exc:
        issue = ValidationIssue("TBL000", Severity.ERROR, str(path), f"could not read TSV: {exc}")
        return None, [issue]
    return frame, []


# Past this many offending rows in one column, report a count instead. A column
# renamed upstream makes every row fail, and mirrors are reviewed by digest.
_MAX_ROW_ISSUES: Final = 20


def _row_issues(
    code: str, path: Path, offenders: list[tuple[int, str]], summary: str
) -> list[ValidationIssue]:
    """One issue per offending row, collapsed to a count past the threshold."""
    if not offenders:
        return []
    if len(offenders) > _MAX_ROW_ISSUES:
        first_row, first_detail = offenders[0]
        return [
            ValidationIssue(
                code,
                Severity.ERROR,
                str(path),
                f"{len(offenders)} rows: {summary} (first at row {first_row}: {first_detail})",
            )
        ]
    return [
        ValidationIssue(code, Severity.ERROR, f"{path}:row {row}", f"{summary}: {detail}")
        for row, detail in offenders
    ]


def validate_table(path: Path, schema: TableSchema) -> list[ValidationIssue]:
    """Check one TSV file against its schema, returning every problem found."""
    issues: list[ValidationIssue] = []

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    frame, read_issues = read_table(path, schema)
    if frame is None:
        return read_issues
    present = set(frame.columns)
    expected = set(schema.column_names)

    for missing in sorted(expected - present):
        error("TBL001", str(path), f"missing required column '{missing}'")
    for unexpected in sorted(present - expected):
        error("TBL002", str(path), f"unexpected column '{unexpected}'")

    for column in schema.columns:
        if column.name not in present:
            continue
        series = frame[column.name]

        if not column.nullable:
            offenders_null = [(row, "empty") for row in _rows_where(series.is_null())]
            issues.extend(
                _row_issues(
                    "TBL003", path, offenders_null, f"column '{column.name}' must not be empty"
                )
            )

        # Keyed on the dtype rather than on the bounds, because a column with no
        # natural range (a log2 fold change) needs this exactly as much as a
        # bounded one. Null and NaN are distinct states in polars: the non-null
        # check above tests `value is None` and so never sees a NaN, while every
        # comparison against NaN is False by IEEE-754 and so the bounds check
        # below cannot fire for one either. Without this, a p-value written as
        # `NaN` — what R's `write.table` and numpy's `savetxt` emit for a missing
        # statistic — validates clean and ships as a real measurement.
        #
        # `Series.is_finite()` returns null for a null entry rather than False, so
        # `.not_()` on it is also null there (Kleene NOT), and `_value_offenders`'s
        # `mask.any()` / `.arg_true()` both treat null as "not an offender" — a
        # null cell is skipped here exactly as `_numeric_values` used to skip it
        # by testing `value is not None` before `math.isfinite`.
        if column.dtype is pl.Float64:
            offenders_finite = _value_offenders(series, series.is_finite().not_())
            issues.extend(
                _row_issues(
                    "TBL010",
                    path,
                    offenders_finite,
                    f"column '{column.name}' has non-finite values",
                )
            )

        if column.allowed is not None:
            allowed_sorted = sorted(column.allowed)
            # `is_in` also returns null for a null cell, so `.not_()` on it stays
            # null and `_value_offenders` skips it — a null never reported here,
            # matching `_string_values`'s `value is not None` filter.
            offenders_allowed = _value_offenders(series, series.is_in(list(column.allowed)).not_())
            issues.extend(
                _row_issues(
                    "TBL004",
                    path,
                    offenders_allowed,
                    f"column '{column.name}' has values outside {allowed_sorted}",
                )
            )

        if column.pattern is not None:
            # Vectorized with polars' `str.contains` (the Rust `regex` crate)
            # rather than Python `re.fullmatch`. The two engines are not the same
            # implementation, so this is only safe because it has been measured,
            # not reasoned about: every value of every pattern-bearing column in
            # every committed mirror (212,113 non-null values, as pinned by the
            # test below on 2026-10-25) plus an adversarial
            # synthetic set (trailing/embedded/leading newlines, case, leading
            # zeros, full-width and Arabic-Indic Unicode digits, stray `^`/`$`
            # characters, a 1,000-digit value) produced byte-identical offender
            # sets under both engines for every pattern used in this file —
            # `HGNC_PATTERN`, `MONDO_PATTERN`, `SEQUENCE_ONTOLOGY_PATTERN`,
            # `MODIFICATION_PATTERN`, `UNIPROT_PATTERN`, the PMID, cohort-list
            # and Orphanet scope patterns (`ORPHA:`, `OMIM:`, `HP:` and the
            # `ORPHA|OMIM` cross-reference target), and the inline
            # `gencc_submissions.sgc_id` and `submitted_on` (ISO date prefix)
            # patterns. See
            # `tests/unit/test_tables.py::test_polars_pattern_engine_agrees_with_python_re`.
            # `re.fullmatch` anchors both ends implicitly; `str.contains` does
            # not, but every pattern above already carries `^` and `$`, so no
            # pattern needed a Python fallback.
            offenders_pattern = _value_offenders(series, series.str.contains(column.pattern).not_())
            issues.extend(
                _row_issues(
                    "TBL005",
                    path,
                    offenders_pattern,
                    f"column '{column.name}' has values not matching {column.pattern}",
                )
            )

        if column.minimum is not None or column.maximum is not None:
            low = column.minimum if column.minimum is not None else float("-inf")
            high = column.maximum if column.maximum is not None else float("inf")
            # Non-finite values are excluded so exactly one code fires per
            # problem: an infinity is out of every finite range, and reporting it
            # as both TBL010 and TBL006 would say the same thing twice while
            # TBL010 names the cause more precisely. `is_finite()` is False for
            # NaN/inf and null for a null cell; combined with the bounds test via
            # `&`, a null cell yields null (Kleene AND) and is skipped by
            # `_value_offenders`, same as `_numeric_values`'s `value is not None`.
            #
            # `series.cast(pl.Float64)` for the gathered values only, not for the
            # mask: this check also runs on Int64 columns (`pos`, `n_case`, ...),
            # and `_numeric_values` used to `float()` every value before `repr`,
            # so an offending int reported "5.0", not "5" — the cast reproduces
            # that here rather than silently changing the message.
            mask_bounds = series.is_finite() & ((series < low) | (series > high))
            offenders_bounds = _value_offenders(series.cast(pl.Float64), mask_bounds)
            issues.extend(
                _row_issues(
                    "TBL006",
                    path,
                    offenders_bounds,
                    f"column '{column.name}' has values outside [{low}, {high}]",
                )
            )

    if all(column in present for column in schema.sort_key):
        keys = frame.select(schema.sort_key).rows()
        repeated = set(duplicates(keys))
        offenders_dupe: list[tuple[int, str]] = []
        seen: set[tuple[object, ...]] = set()
        for index, key in enumerate(keys):
            if key not in repeated:
                continue
            if key in seen:
                offenders_dupe.append((index + 2, repr(key)))
            else:
                seen.add(key)
        issues.extend(
            _row_issues("TBL007", path, offenders_dupe, f"duplicate {list(schema.sort_key)} key")
        )

    return issues


def _rows_where(mask: pl.Series) -> list[int]:
    """Row numbers (1-based, header counted as row 1) where `mask` is true.

    `mask.any()` runs entirely inside polars and treats a null entry as
    not-true (its default `ignore_nulls=True`), so a clean column — the normal
    case for every check in this module — costs one boolean reduction and
    never calls `.to_list()`. `arg_true()` likewise selects exact `True`
    values only, so a null in `mask` (every check below produces one for a
    null cell, via `is_null`, `is_finite`, `is_in` or `str.contains` on a null
    input) is never counted as an offender.
    """
    if not mask.any():
        return []
    return [int(position) + 2 for position in mask.arg_true().to_list()]


def _value_offenders(series: pl.Series, mask: pl.Series) -> list[tuple[int, str]]:
    """(row, repr(value)) pairs for the positions where `mask` is true.

    Only the offending positions are gathered out of `series` — never the
    whole column — so this stays cheap even when `series` holds tens of
    thousands of rows and the offender count is the normal zero.
    """
    if not mask.any():
        return []
    positions = mask.arg_true()
    values = series.gather(positions).to_list()
    return [
        (int(position) + 2, repr(value))
        for position, value in zip(positions.to_list(), values, strict=True)
    ]


_SOURCE_TIERS = frozenset(tier.value for tier in SourceTier)
_DIRECTIONS = frozenset({"up", "down", "unchanged"})

GENES = TableSchema(
    name="genes",
    columns=(
        Column("hgnc_id", pl.String, pattern=HGNC_PATTERN),
        Column("symbol", pl.String),
        Column("name", pl.String),
        Column("aliases", pl.String, nullable=True),
        # Separate from `aliases`, never merged into it. HGNC retires a symbol
        # into `prev_symbol` rather than into `alias_symbol`, and `genes.py`
        # resolves an alias *ahead* of a previous symbol -- so a string that is
        # one gene's live alias and another's retired name resolves to the
        # first. Merging the two columns here would erase the distinction that
        # precedence turns on and make `ODD` (alias of GJA1, previous of OSR1)
        # ambiguous instead of answerable.
        Column("prev_symbols", pl.String, nullable=True),
        Column("ensembl_gene", pl.String, nullable=True),
        Column("ncbi_gene", pl.Int64, nullable=True),
        Column("locus", pl.String, nullable=True),
        Column("uniprot", pl.String, nullable=True, pattern=UNIPROT_PATTERN),
        Column("mane_select", pl.String, nullable=True),
    ),
    sort_key=("hgnc_id",),
)

PTM_SITES = TableSchema(
    name="ptm_sites",
    columns=(
        Column("site_id", pl.String),
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y", "K", "R", "C", "N", "Q"})),
        Column("position", pl.Int64, minimum=1),
        Column("mod_type", pl.String, pattern=MODIFICATION_PATTERN),
        Column("flanking_sequence", pl.String, nullable=True),
        Column("known_kinases", pl.String, nullable=True),
        Column("source", pl.String),
    ),
    sort_key=("protein", "position", "mod_type"),
)

# ClinGen's own vocabulary, stored verbatim. Mapping onto `Classification`
# happens once in `vocab.CLINGEN_CLASSIFICATIONS`, so the mirror stays a
# faithful copy of what the authority published and a term ClinGen adds later
# fails loudly here rather than being silently coerced downstream.
_CLINGEN_CLASSIFICATIONS: Final[frozenset[str]] = frozenset(
    {
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed",
        "Refuted",
        "No Known Disease Relationship",
    }
)
# Measured 2026-08-03 against the ClinGen bulk CSV's own MOI column: {'AD',
# 'AR', 'MT', 'SD', 'UD', 'XL'}. ClinGen publishes "UD", not "Undetermined
# MOI" -- the latter never appears in the file.
_CLINGEN_MOI: Final[frozenset[str]] = frozenset({"AD", "AR", "XL", "MT", "SD", "UD"})

CLINGEN_VALIDITY = TableSchema(
    name="clingen_validity",
    columns=(
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("gene_symbol", pl.String),
        Column("disease", pl.String, pattern=MONDO_PATTERN),
        Column("disease_label", pl.String),
        Column("moi", pl.String, allowed=_CLINGEN_MOI),
        Column("sop", pl.String),
        Column("classification", pl.String, allowed=_CLINGEN_CLASSIFICATIONS),
        Column("classification_date", pl.String),
        Column("gcep", pl.String),
        Column("report_url", pl.String),
    ),
    # The triple, not the pair. Measured 2026-08-03 against the ClinGen bulk
    # CSV (3,653 rows): (gene, disease, moi) is unique at 3,653/3,653;
    # (gene, disease) at 3,594 -- 59 pairs are curated twice under different
    # MOI with different classifications.
    sort_key=("gene", "disease", "moi"),
)

# GenCC's own vocabulary (the harmonised `classification_title` column of its
# submissions-export TSV), stored verbatim -- same rationale as
# `_CLINGEN_CLASSIFICATIONS`. Measured 2026-08-03 against the real export
# (30,410 rows): {'Definitive', 'Disputed Evidence', 'Limited', 'Moderate',
# 'No Known Disease Relationship', 'Refuted Evidence', 'Strong', 'Supportive'}.
# 'Supportive' is a real GenCC term, not a plan guess -- it is a mapping
# exception for submitters (e.g. Orphanet) that do not grade evidence on
# ClinGen's ladder; `vocab.GENCC_CLASSIFICATIONS` maps it to `None`.
_GENCC_CLASSIFICATIONS: Final[frozenset[str]] = frozenset(
    {
        "Definitive",
        "Strong",
        "Moderate",
        "Limited",
        "Disputed Evidence",
        "Refuted Evidence",
        "No Known Disease Relationship",
        "Supportive",
    }
)

GENCC_SUBMISSIONS = TableSchema(
    name="gencc_submissions",
    columns=(
        # GenCC's own row id, e.g. "SGC-102815". Measured 2026-08-03 against
        # the real export: matches `^SGC-\d+$` and is unique on 30,410/30,410
        # rows -- it is GenCC's primary key, not a value this atlas invents.
        Column("sgc_id", pl.String, pattern=r"^SGC-\d+$"),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("gene_symbol", pl.String),
        Column("disease", pl.String, pattern=MONDO_PATTERN),
        Column("disease_label", pl.String),
        Column("moi", pl.String),
        Column("classification", pl.String, allowed=_GENCC_CLASSIFICATIONS),
        Column("submitter", pl.String),
        # Every GenCC value begins with an ISO calendar date (30,410 of 30,410
        # measured); what follows varies. The column stays nullable on purpose: a
        # third party's missing date must not fail the whole build. The pattern
        # skips nulls, so `docs/data-api.md` promises the prefix for any
        # non-null value only, and the gene page truncates to it. The trailing `.*$` is not
        # decoration: the engine-agreement test compares `re.fullmatch` with
        # `str.contains`, and a bare prefix pattern would disagree under it.
        Column("submitted_on", pl.String, nullable=True, pattern=r"^\d{4}-\d{2}-\d{2}.*$"),
        Column("report_url", pl.String, nullable=True),
    ),
    # The submitter, and then its own row id -- see
    # `scripts/convert_gencc.py` for why: GenCC publishes every submitter's
    # verdict with none adjudicated, so two submitters disagreeing about one
    # gene-disease-moi is the normal case, and `submitter` alone is still not
    # enough to key on.
    #
    # Measured 2026-08-03 against the real export (30,410 rows, ?format=new):
    # (gene, disease, moi, submitter) collides on 133 groups (134 extra
    # rows), and 70 of those 133 groups carry more than one distinct
    # `classification_title` -- e.g. Ambry Genetics submitting both `Limited`
    # (SGC-102815) and `Strong` (SGC-104042) for HGNC:20226 / MONDO:0859332 /
    # Autosomal recessive. `version_number` is identical within all 133
    # groups, so GenCC marks neither row as superseding the other, and the
    # `submitted_as_date` in the disagreeing cases differs by as little as
    # one second -- a batch-load artifact, not a real submission gap. There
    # is no field in this export that says which of two same-submitter,
    # same-gene-disease-moi rows is the "current" one. Picking a winner would
    # synthesise a verdict GenCC itself does not assert -- the plan's design
    # decision D12 ("the atlas publishes no validity classification of its
    # own") is written about the derive layer built on top of this mirror,
    # not this table directly, but the same reasoning applies one layer
    # earlier: a submitter disagreeing with itself is still discordance, and
    # GenCC publishes both rows, so this table does too. `sgc_id` is GenCC's
    # own row id and is unique across the whole export, so appending it to
    # the key is the only column that actually makes every real GenCC row
    # distinguishable without dropping any of them.
    sort_key=("gene", "disease", "moi", "submitter", "sgc_id"),
)

# The chromosomes, in karyotype order. Ordered because the published variant
# index drives a chromosome picker, and sorting shard names as strings gives
# 1, 10, 11, … 2, 20 with MT ahead of X — jumbled in both the numeric and the
# non-numeric part.
#
# It lives here, beside the schema, rather than in `build/variants.py` where it
# began. The `allowed` set below is derived from it, so the vocabulary and the
# order cannot drift — which retired the test that used to check they had not.
# It also lets `unexpected_mirror_entries` report a shard named for no
# chromosome without a validator importing from the build layer.
CHROMOSOMES: Final[tuple[str, ...]] = (*(str(n) for n in range(1, 23)), "X", "Y", "MT")
_CHROMOSOME_SET: Final[frozenset[str]] = frozenset(CHROMOSOMES)

VARIANTS = TableSchema(
    name="variants",
    columns=(
        Column("vrs_id", pl.String),
        Column("assembly", pl.String, allowed=frozenset({"GRCh38"})),
        Column("chrom", pl.String, allowed=_CHROMOSOME_SET),
        Column("pos", pl.Int64, minimum=1),
        Column("ref", pl.String),
        Column("alt", pl.String),
        Column("hgvs_g", pl.String, nullable=True),
        Column("hgvs_c", pl.String, nullable=True),
        Column("hgvs_p", pl.String, nullable=True),
        Column("mane_transcript", pl.String, nullable=True),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("consequence", pl.String, pattern=SEQUENCE_ONTOLOGY_PATTERN),
        Column("clinvar_vcv", pl.String, nullable=True),
        Column("clinvar_significance", pl.String, nullable=True),
        Column("clinvar_review_status", pl.String, nullable=True),
        Column("gnomad_af_popmax", pl.Float64, nullable=True),
        Column("gnomad_version", pl.String, nullable=True),
        Column("acmg_criteria", pl.String, nullable=True),
        Column("inheritance", pl.String, nullable=True),
        Column("de_novo", pl.Boolean, nullable=True),
        Column("zygosity", pl.String, nullable=True),
        Column("phenotypes", pl.String, nullable=True),
        Column("publication", pl.String, nullable=True),
        Column("locator", pl.String, nullable=True),
        Column("source_tier", pl.String, allowed=_SOURCE_TIERS),
        Column("source", pl.String),
        # Non-coding extension (design spec D10); null for coding variants.
        Column("regulatory_element", pl.String, nullable=True),
        Column("target_gene", pl.String, nullable=True),
        Column("assay_evidence", pl.String, nullable=True),
    ),
    sort_key=("chrom", "pos", "ref", "alt", "vrs_id"),
)

EXPRESSION = TableSchema(
    name="expression",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_case", pl.Int64, nullable=True, minimum=0),
        Column("n_control", pl.Int64, nullable=True, minimum=0),
        Column("tissue", pl.String),
        Column("stage", pl.String, nullable=True),
    ),
    sort_key=("contrast", "gene"),
)

# The unit vocabulary is shared with PROFILES rather than repeated: PRF001
# compares a dataset's two tables cell for cell, and two literals that drifted
# apart would make that check compare different alphabets.
_ABUNDANCE_UNITS: Final = frozenset({"tpm", "nx", "cpm", "lfq", "rpkm"})

PROFILES = TableSchema(
    name="profiles",
    columns=(
        Column("dataset", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("tissue", pl.String),
        Column("stage", pl.String, nullable=True),
        Column("median_abundance", pl.Float64),
        # `rpkm` added 2026-08-14 with the bulk developmental layer. The source
        # publishes RPKM; converting to TPM would be the atlas authoring a
        # number, and this column exists precisely so two incomparable
        # quantities cannot merge under one header.
        Column("unit", pl.String, allowed=_ABUNDANCE_UNITS),
        Column("q25", pl.Float64, nullable=True),
        Column("q75", pl.Float64, nullable=True),
        Column("n_samples", pl.Int64, minimum=1),
    ),
    sort_key=("dataset", "gene", "tissue", "stage"),
)

# Four chambers and the interventricular septum: the regions a cell-type
# dataset may declare (`Dataset.regions`). A row's `region` is one of them or
# `POOLED_REGION`, the row that pools them; one column, so a reader filters one
# field rather than two. Defined once: the model, the converter, the validator
# and the build import these rather than retype them.
CELL_TYPE_REGIONS: Final = frozenset({"LA", "RA", "LV", "RV", "IVS"})
POOLED_REGION: Final = "all"
# The two rules the cell-type block prints of every value, defined once. A
# heart counts towards a (cell type, region) only with at least this many cells
# of it there: hvantk drops a smaller group (`--min-cells`, 10 in the cluster
# job), `convert_farah2024.py` refuses one that arrives anyway, and CTP010
# checks what the mirror can show of it.
MIN_CELLS_PER_HEART: Final = 10
# Quartiles are published from this many hearts; below it they are null
# (`convert_farah2024.py` writes them so, and CTP009 holds the mirror to it).
MIN_HEARTS_FOR_QUARTILES: Final = 3

CELL_TYPE_PROFILES = TableSchema(
    name="cell_type_profiles",
    columns=(
        Column("dataset", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("cell_type", pl.String),
        Column("region", pl.String, allowed=CELL_TYPE_REGIONS | {POOLED_REGION}),
        Column("median_abundance", pl.Float64, minimum=0),
        # `cpm` only, not `_ABUNDANCE_UNITS`: the gene page writes "CPM" on its
        # axis and beside every value whatever this column says, so a second
        # unit would publish under the wrong name.
        Column("unit", pl.String, allowed=frozenset({"cpm"})),
        Column("q25", pl.Float64, nullable=True, minimum=0),
        Column("q75", pl.Float64, nullable=True, minimum=0),
        Column("n_donors", pl.Int64, minimum=1),
        Column("n_cells", pl.Int64, minimum=1),
        Column("fraction_expressing", pl.Float64, minimum=0, maximum=1),
        # Names the file the numbers were taken from; its licence and sha256 live
        # in sources.yaml, and carrying the column is what holds each value to
        # SRC003 (`_used_sources` collects it from any table that has one).
        Column("source", pl.String),
    ),
    sort_key=("dataset", "gene", "cell_type", "region"),
)

PROFILE_QUANTILES = TableSchema(
    name="profile_quantiles",
    columns=(
        Column("dataset", pl.String),
        Column("tissue", pl.String),
        # Not nullable, unlike `profiles.stage`: a quantile grid always names
        # the stage it summarizes. The asymmetry is a hazard for whoever joins
        # the two tables -- a null-stage `profiles` row can never have a
        # quantile partner by construction, so a naive inner join on
        # (dataset, tissue, stage) silently drops it (measured: 2 rows in, 1
        # row out, no error). The join must be a left join, so the null-stage
        # row still publishes, with a null percentile.
        Column("stage", pl.String),
        Column("percentile", pl.Int64, minimum=0, maximum=100),
        Column("value", pl.Float64),
        Column("unit", pl.String, allowed=_ABUNDANCE_UNITS),
        Column("n_genes", pl.Int64, minimum=1),
    ),
    sort_key=("dataset", "tissue", "stage", "percentile"),
)

PROTEOMICS = TableSchema(
    name="proteomics",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("gene", pl.String, pattern=HGNC_PATTERN, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_peptides", pl.Int64, nullable=True, minimum=1),
        Column("sequence_coverage", pl.Float64, nullable=True, minimum=0, maximum=100),
        Column("quant_method", pl.String),
    ),
    sort_key=("contrast", "protein"),
)

PHOSPHO = TableSchema(
    name="phospho",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("site_id", pl.String),
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y"})),
        Column("position", pl.Int64, minimum=1),
        Column("mod_type", pl.String, pattern=MODIFICATION_PATTERN),
        Column("flanking_sequence", pl.String, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
        # Mandatory: whether the change was corrected for protein abundance decides
        # whether the site is genuinely regulated or merely tracking total protein.
        Column("protein_normalized", pl.Boolean),
        Column("occupancy", pl.Float64, nullable=True),
        Column("known_kinases", pl.String, nullable=True),
    ),
    sort_key=("contrast", "protein", "position", "mod_type"),
)

# Derived from the enums so the mirror's vocabulary and the code's cannot drift.
_COHORT_STRATA = frozenset(item.value for item in CohortStratum)
_VARIANT_CLASSES = frozenset(item.value for item in VariantClass)
_CONSEQUENCE_CLASSES = frozenset(item.value for item in ConsequenceClass)
_VARIANT_ORIGINS = frozenset(item.value for item in VariantOrigin)
_COMPARATORS = frozenset(item.value for item in BurdenComparator)
_EFFECT_MEASURES = frozenset(item.value for item in EffectMeasure)
_EFFECT_BOUNDS = frozenset(item.value for item in EffectBound)
_STATISTICAL_TESTS = frozenset(item.value for item in StatisticalTest)
_COUNT_UNITS = frozenset(item.value for item in CountUnit)
_PVALUE_ADJUSTMENTS = frozenset(item.value for item in PvalueAdjustment)
_LESION_GROUPS = frozenset(item.value for item in LesionGroup)
_COMBINATION_METHODS = frozenset(item.value for item in CombinationMethod)

BURDEN = TableSchema(
    name="burden",
    columns=(
        # Identity. `study` cites the paper; `source` names the file the numbers
        # were transcribed out of, whose licence and sha256 live in sources.yaml.
        # Both are needed and neither implies the other: one supplement can carry
        # several studies' results, and one study can publish several supplements.
        Column("study", pl.String, pattern=PMID_PATTERN),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        # The partition. Every column here is a degree of freedom that changes
        # what the counts mean, and the set exists to stop two incomparable rows
        # from *looking* comparable -- not to make them easier to compare. Drop
        # `maf_max` and a MAF<1e-4 test sits beside a MAF<1e-2 test unlabelled;
        # drop `origin` and a trio's de novo count reads like a case-control
        # rare-variant count.
        Column("cohort_stratum", pl.String, allowed=_COHORT_STRATA),
        # Null unless the study analysed one lesion group separately. Reuses the
        # curated `LesionGroup` vocabulary rather than inventing a parallel one,
        # so a lesion-stratified burden row and a curated assertion name the same
        # thing with the same token.
        Column("lesion_group", pl.String, nullable=True, allowed=_LESION_GROUPS),
        Column("variant_class", pl.String, allowed=_VARIANT_CLASSES),
        Column("consequence_class", pl.String, allowed=_CONSEQUENCE_CLASSES),
        Column("origin", pl.String, allowed=_VARIANT_ORIGINS),
        # Null where the analysis applied no frequency filter, which is normal
        # for CNV burden. Not defaulted to 1.0: "no filter" and "filtered at
        # fixation" are different claims and only one of them was made.
        Column("maf_max", pl.Float64, nullable=True, minimum=0, maximum=1),
        # The observation. Both mandatory: a burden row without a denominator is
        # a count nobody can interpret. `count_unit` is mandatory for the same
        # reason one step further in -- a numerator and a denominator nobody can
        # interpret *the units of* is no better. It qualifies all four count
        # columns below, cases and controls alike, because a study that counts
        # case alleles counts control alleles too.
        Column("count_unit", pl.String, allowed=_COUNT_UNITS),
        Column("n_case_carriers", pl.Int64, minimum=0),
        Column("n_cases", pl.Int64, minimum=1),
        # Of `n_case_carriers`, how many are congenital heart disease cases,
        # where the study's case set also holds cases of a condition outside
        # this atlas's scope and the study counts its CHD cases' carriers apart
        # (owner decision R2, 2026-10-26). Empty where it does not. Only the
        # studies in `CHD_CASE_SPLITS` below may carry it (BUR024), and only
        # PMID:34324492 is one: its case set combines CHD with a series of 777
        # sporadic thoracic aortic aneurysm probands, and its S8 Table gives
        # the "only CHD" case count for 521 of the 528 genes it lists, leaving
        # 7 blank, which `convert_audain2021.py` transcribes. `source` still
        # names the file the row's statistic came from. Same `count_unit` as
        # `n_case_carriers`, and never above it (BUR003).
        Column("n_chd_case_carriers", pl.Int64, nullable=True, minimum=0),
        # The comparator, and the fields each kind requires. `validate_burden`
        # enforces which of the next three must be present and which must be
        # empty; the schema alone cannot express that dependency.
        Column("comparator", pl.String, allowed=_COMPARATORS),
        Column("n_control_carriers", pl.Int64, nullable=True, minimum=0),
        Column("n_controls", pl.Int64, nullable=True, minimum=1),
        Column("expected_count", pl.Float64, nullable=True, minimum=0),
        # The statistic. `effect` is null both when none was published and when
        # the published one was non-finite -- `effect_bound` distinguishes those,
        # and `validate_burden` requires `effect_measure` for either.
        Column("effect", pl.Float64, nullable=True, minimum=0),
        Column("effect_measure", pl.String, nullable=True, allowed=_EFFECT_MEASURES),
        Column("effect_bound", pl.String, nullable=True, allowed=_EFFECT_BOUNDS),
        Column("ci_low", pl.Float64, nullable=True, minimum=0),
        Column("ci_high", pl.Float64, nullable=True, minimum=0),
        Column("pvalue", pl.Float64, nullable=True, minimum=0, maximum=1),
        Column("pvalue_test", pl.String, nullable=True, allowed=_STATISTICAL_TESTS),
        # A corrected p-value the *study* published, never one the atlas
        # computed. Added 2026-08-05: the review found gene pages showing
        # uncorrected p-values with nothing to judge them against, and
        # PMID:34324492 publishes a correction beside every raw p. Dropping it
        # would discard the one number that answers that.
        Column("pvalue_adjusted", pl.Float64, nullable=True, minimum=0, maximum=1),
        Column("pvalue_adjustment", pl.String, nullable=True, allowed=_PVALUE_ADJUSTMENTS),
        # Provenance. Cohort membership is what makes overlap between studies
        # visible: DDD contributes cases to more than one paper here, so two rows
        # a reader would otherwise mentally meta-analyse can be shown to share
        # samples. `;`-joined, same convention as `genes.aliases`.
        Column("case_cohorts", pl.String, pattern=COHORT_LIST_PATTERN),
        Column("control_cohorts", pl.String, nullable=True, pattern=COHORT_LIST_PATTERN),
        # Free text, rendered verbatim as a footnote under the study's table by
        # `pages._footnotes`. That sentence was false from this column's
        # introduction until 2026-08-05: nothing read it, so a curator stating a
        # CNV carrier rule would have got a green build, a correct bundle and no
        # page text. Deliberately not a
        # controlled vocabulary yet: what "carrier" means for a CNV differs
        # between papers (any overlap, exonic overlap, whole gene), and inventing
        # the axis from a single example is how a vocabulary ends up wrong.
        Column("method_note", pl.String, nullable=True),
        Column("source", pl.String),
    ),
    # The full partition, which is what makes a row unique. Gene first so the
    # file reads gene by gene, which is also how a curator reviews a diff.
    sort_key=(
        "gene",
        "study",
        "cohort_stratum",
        "lesion_group",
        "variant_class",
        "consequence_class",
        "origin",
        "maf_max",
    ),
)


@dataclass(frozen=True)
class CaseSplit:
    """Where one study's `n_chd_case_carriers` comes from, and what it excludes.

    `source` is the `mirrors/sources.yaml` id of the file the column is
    transcribed from -- not the row's own `source`, which names the file its
    statistic came from -- and `table` is that file's name as a gene page's
    card prints it ("the study's S8 Table"). `other_case` is one case of the
    condition outside this atlas's scope that the study's case set also holds,
    with its article, as the card's sentence ends where every carrier is a CHD
    case ("so none is an aneurysm proband").
    """

    source: str
    table: str
    other_case: str


# The studies whose burden rows may carry `n_chd_case_carriers`, keyed by study
# (owner decision R2). One definition for its three readers: `validate_burden`
# refuses the column on any other study's row (BUR024); the validation runner
# holds `source` to the registry wherever a row carries the column, as it holds
# every row's own `source` (SRC003, SRC004); and `pages._chd_only` takes the
# table and the excluded condition from here, so a card never names one study's
# table, or its other condition, for another study's row.
CHD_CASE_SPLITS: Final[dict[str, CaseSplit]] = {
    "PMID:34324492": CaseSplit(
        source="audain2021_s8", table="S8 Table", other_case="an aneurysm proband"
    ),
}

# One study's own combination of two of its burden arms. Added 2026-10-24 for
# PMID:40127276's JL-FDR meta-analysis, which `convert_sierant2025.py` declined to
# put in `pvalue_adjusted` -- correctly, because it is not a corrected p.
BURDEN_COMBINED = TableSchema(
    name="burden_combined",
    columns=(
        Column("study", pl.String, pattern=PMID_PATTERN),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        # The burden partition minus what the arms differ in (`origin`,
        # `maf_max`, `variant_class`, and the comparator itself) -- exactly what
        # a combination spans. `lesion_group` stays so the paper's
        # per-phenotype meta-analysis can land without a migration.
        Column("cohort_stratum", pl.String, allowed=_COHORT_STRATA),
        Column("lesion_group", pl.String, nullable=True, allowed=_LESION_GROUPS),
        Column("consequence_class", pl.String, allowed=_CONSEQUENCE_CLASSES),
        # The comparators of the arms combined, `;`-joined, in the order the page
        # names them. Token *syntax* reuses `COHORT_LIST_PATTERN` rather than
        # adding a pattern to the engine-agreement test; vocabulary and arity
        # (at least two distinct comparators) are BUR021's.
        Column("combines", pl.String, pattern=COHORT_LIST_PATTERN),
        Column("method", pl.String, allowed=_COMBINATION_METHODS),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        # A *local* false discovery rate, never a corrected p: measured over DS6
        # it is smaller than its own p for 118, 135 and 129 of 248 genes by
        # consequence class. Named for what it is so no renderer can file it
        # under "corrected p".
        Column("local_fdr", pl.Float64, minimum=0, maximum=1),
        Column("method_note", pl.String, nullable=True),
        Column("source", pl.String),
    ),
    # Gene first, then study, for the reason `BURDEN`'s key is: the file reads
    # gene by gene, and `burden_mirror.merge_study` sorts both the same way.
    sort_key=("gene", "study", "cohort_stratum", "lesion_group", "consequence_class"),
)

# The three mirrors the Orphanet scope rule reads (`orphanet_scope.py`), added
# 2026-10-25 for issue #8. Each is wider than the rule needs, so changing the
# rule never needs a re-mirror: product 4 is kept for the whole cardiovascular
# branch, not only the structural subset the rule counts.
#
# Orphanet's frequency band names, verbatim from product 4. Measured 2026-10-25:
# these six and no others over its 116,664 annotations.
ORPHANET_FREQUENCIES: Final = frozenset(
    {
        "Obligate (100%)",
        "Very frequent (99-80%)",
        "Frequent (79-30%)",
        "Occasional (29-5%)",
        "Very rare (<4-1%)",
        "Excluded (0%)",
    }
)
# Product 1's mapping-relation codes, the token before each gloss. Measured: E
# 3,781, BTNT 4,138, NTBT 797, ND 29 over 8,745 Orphanet-to-OMIM references.
ORPHANET_RELATIONS: Final = frozenset({"E", "BTNT", "NTBT", "ND"})
ORPHANET_VALIDATION: Final = frozenset({"Validated", "Not yet validated"})

ORPHANET_PHENOTYPES = TableSchema(
    name="orphanet_phenotypes",
    columns=(
        Column("orpha", pl.String, pattern=ORPHA_PATTERN),
        Column("disorder_name", pl.String),
        Column("hpo", pl.String, pattern=HPO_PATTERN),
        Column("hpo_label", pl.String),
        # Null when Orphanet gave no band; such an annotation never clears the
        # rule's floor.
        Column("frequency", pl.String, nullable=True, allowed=ORPHANET_FREQUENCIES),
        Column("source", pl.String),
    ),
    # `frequency` is in the key because product 4 annotates 38 (orpha, hpo)
    # pairs twice; both are transcribed and neither is chosen.
    sort_key=("orpha", "hpo", "frequency"),
)
ORPHANET_OMIM = TableSchema(
    name="orphanet_omim",
    columns=(
        Column("orpha", pl.String, pattern=ORPHA_PATTERN),
        Column("omim", pl.String, pattern=OMIM_PATTERN),
        Column("relation", pl.String, allowed=ORPHANET_RELATIONS),
        Column("validation", pl.String, allowed=ORPHANET_VALIDATION),
        Column("source", pl.String),
    ),
    sort_key=("orpha", "omim"),
)
MONDO_XREFS = TableSchema(
    name="mondo_xrefs",
    columns=(
        Column("mondo", pl.String, pattern=MONDO_PATTERN),
        # `skos:exactMatch` only; the converter drops every other predicate.
        Column("target", pl.String, pattern=XREF_TARGET_PATTERN),
        Column("source", pl.String),
    ),
    sort_key=("mondo", "target"),
)

TABLE_SCHEMAS: Final[dict[str, TableSchema]] = {
    schema.name: schema
    for schema in (
        GENES,
        PTM_SITES,
        VARIANTS,
        EXPRESSION,
        PROFILES,
        CELL_TYPE_PROFILES,
        PROFILE_QUANTILES,
        PROTEOMICS,
        PHOSPHO,
        CLINGEN_VALIDITY,
        GENCC_SUBMISSIONS,
        BURDEN,
        BURDEN_COMBINED,
        ORPHANET_PHENOTYPES,
        ORPHANET_OMIM,
        MONDO_XREFS,
    )
}

# Which directory layout each schema lives under, relative to ``mirrors/``.
SHARDED_TABLES: Final[dict[str, str]] = {
    "variants": "variants",
    "expression": "expression",
    "profiles": "profiles",
    "cell_type_profiles": "cell_type_profiles",
    "profile_quantiles": "profile_quantiles",
    "proteomics": "proteomics",
    "phospho": "phospho",
}
FLAT_TABLES: Final[dict[str, str]] = {
    "genes": "genes.tsv",
    "ptm_sites": "ptm_sites.tsv",
    "clingen_validity": "clingen_gene_validity.tsv",
    "gencc_submissions": "gencc_submissions.tsv",
    # Flat rather than sharded: one study contributes a few thousand in-scope
    # rows, and the whole table is read at once to group a gene's rows by
    # (study, stratum, maf) — a per-gene shard would make that the expensive
    # access pattern rather than the cheap one.
    "burden": "burden.tsv",
    # Beside `burden`, not inside it: a combination is a statement about two
    # rows already in that table, and folding it in would need a fake arm.
    "burden_combined": "burden_combined.tsv",
    "orphanet_phenotypes": "orphanet_phenotypes.tsv",
    "orphanet_omim": "orphanet_omim.tsv",
    "mondo_xrefs": "mondo_xrefs.tsv",
}


def unexpected_mirror_entries(root: Path) -> list[ValidationIssue]:
    """Anything under ``mirrors/`` that no schema claims.

    ``Path.glob`` yields nothing for a missing or non-directory path, so a shard
    directory lost to a typo — ``mirrors/varaints`` — simply vanishes from
    validation and the repository reports clean. Naming the stray entry is what
    catches that. An *absent* shard directory stays legitimate: before any
    variants are curated there is no ``mirrors/variants/``.
    """
    mirrors = root / "mirrors"
    if not mirrors.is_dir():
        return []

    expected_dirs = set(SHARDED_TABLES.values())
    expected_files = set(FLAT_TABLES.values()) | {"sources.yaml"}
    issues: list[ValidationIssue] = []

    entries, listing_issues = list_dir(mirrors, "TBL009")
    issues.extend(listing_issues)

    for entry in entries:
        if entry.is_dir():
            if entry.name not in expected_dirs:
                issues.append(
                    ValidationIssue(
                        "TBL009",
                        Severity.ERROR,
                        str(entry),
                        f"unexpected directory under mirrors/; expected one of "
                        f"{sorted(expected_dirs)}",
                    )
                )
                continue
            # mirror_paths globs *.tsv, so a shard given the wrong extension is
            # invisible to every table check while the gate still passes.
            shards, shard_listing_issues = list_dir(entry, "TBL009")
            issues.extend(shard_listing_issues)
            for shard in shards:
                if not shard.is_file():
                    issues.append(
                        ValidationIssue(
                            "TBL009",
                            Severity.ERROR,
                            str(shard),
                            f"'{shard.name}' should be a .tsv shard, not a directory",
                        )
                    )
                elif shard.suffix != ".tsv":
                    issues.append(
                        ValidationIssue(
                            "TBL009",
                            Severity.ERROR,
                            str(shard),
                            f"shard files must end .tsv; '{shard.name}' does not",
                        )
                    )
                # A variant shard's filename is its whole addressing scheme: it
                # names the published `variants/<chrom>.json.gz` and a consumer
                # turns a chromosome into a URL by that rule alone. Nothing else
                # checked it — `validate_table` checks the `chrom` column, and
                # the sweep above checks only the extension — so `chr12.tsv`
                # validated at 0 errors and then raised inside the build, giving
                # a curator a traceback where a report belongs. The vocabulary is
                # the schema's own, derived from `CHROMOSOMES` above, so this
                # cannot disagree with the column it addresses.
                elif entry.name == SHARDED_TABLES["variants"] and shard.stem not in _CHROMOSOME_SET:
                    issues.append(
                        ValidationIssue(
                            "TBL011",
                            Severity.ERROR,
                            str(shard),
                            f"a variant shard is named for the chromosome it holds, and "
                            f"'{shard.stem}' is not one of {list(CHROMOSOMES)}",
                        )
                    )
        # A shard directory replaced by a regular file is checked before the
        # unexpected-file case: the two name sets are disjoint, so testing
        # membership of `expected_files` first would report `mirrors/variants`
        # as merely unnamed rather than as the wrong kind of entry.
        elif entry.name in expected_dirs:
            issues.append(
                ValidationIssue(
                    "TBL009",
                    Severity.ERROR,
                    str(entry),
                    f"'{entry.name}' should be a directory of shards, not a file",
                )
            )
        elif entry.name not in expected_files:
            issues.append(
                ValidationIssue(
                    "TBL009",
                    Severity.ERROR,
                    str(entry),
                    f"unexpected file under mirrors/; expected one of {sorted(expected_files)}",
                )
            )
    return issues


def mirror_paths(root: Path) -> list[tuple[Path, str]]:
    """Every mirror TSV under ``root/mirrors``, paired with its schema name."""
    paths: list[tuple[Path, str]] = []
    mirrors = root / "mirrors"
    for schema_name, filename in sorted(FLAT_TABLES.items()):
        path = mirrors / filename
        if path.is_file():
            paths.append((path, schema_name))
    for schema_name, directory in sorted(SHARDED_TABLES.items()):
        for path in sorted((mirrors / directory).glob("*.tsv")):
            paths.append((path, schema_name))
    return paths
