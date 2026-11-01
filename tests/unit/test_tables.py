# tests/unit/test_tables.py
import math
import re
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.duplicates import duplicates
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import (
    BURDEN,
    BURDEN_COMBINED,
    CHROMOSOMES,
    CLINGEN_VALIDITY,
    FLAT_TABLES,
    GENCC_SUBMISSIONS,
    SHARDED_TABLES,
    TABLE_SCHEMAS,
    Column,
    TableSchema,
    _row_issues,
    mirror_paths,
    read_table,
    unexpected_mirror_entries,
    validate_table,
)

REPO_ROOT = Path(__file__).parent.parent.parent

SCHEMA = TableSchema(
    name="demo",
    columns=(
        Column("gene", pl.String, pattern=r"^HGNC:\d+$"),
        Column("pos", pl.Int64),
        Column("tier", pl.String, allowed=frozenset({"a", "b"})),
        Column("note", pl.String, nullable=True),
    ),
    sort_key=("gene", "pos"),
)


def _write(path: Path, rows: str) -> Path:
    path.write_text(rows)
    return path


def test_accepts_a_conforming_table(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:2\t200\tb\thello\n",
    )
    assert validate_table(path, SCHEMA) == []


def test_reports_a_missing_column(tmp_path: Path) -> None:
    path = _write(tmp_path / "demo.tsv", "gene\tpos\ttier\nHGNC:1\t100\ta\n")
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL001"]
    assert "note" in issues[0].message


def test_reports_an_unexpected_column(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\textra\nHGNC:1\t100\ta\t\tx\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL002"]


def test_reports_a_null_in_a_non_nullable_column(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t\ta\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL003"]
    assert "row 2" in issues[0].location


def test_reports_a_value_outside_the_allowed_set(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\tz\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL004"]


def test_reports_a_value_failing_its_pattern(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nENSG001\t100\ta\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL005"]


@pytest.mark.parametrize(
    ("submitted_on", "codes"),
    [
        ("2018-05-31 20:04:32", []),
        ("2025-01-17T00:00:00.000000Z", []),
        ("2021-09-08T16:00:00.000Z", []),
        ("", []),  # null is allowed
        ("31/05/2018", ["TBL005"]),
        ("May 2018", ["TBL005"]),
        ("20180531", ["TBL005"]),
    ],
)
def test_gencc_submitted_on_must_begin_with_a_calendar_date(
    tmp_path: Path, submitted_on: str, codes: list[str]
) -> None:
    """The three real shapes pass; a value the date prefix would not truncate
    cleanly is reported, because `docs/data-api.md` promises the prefix."""
    row = "\t".join(
        [
            "SGC-1",
            "HGNC:1",
            "G",
            "MONDO:0000001",
            "d",
            "Unknown",
            "Strong",
            "S",
            submitted_on,
            "",
        ]
    )
    header = "\t".join(c.name for c in GENCC_SUBMISSIONS.columns)
    path = _write(tmp_path / "gencc_submissions.tsv", f"{header}\n{row}\n")
    issues = validate_table(path, GENCC_SUBMISSIONS)
    assert [i.code for i in issues] == codes


def test_integer_columns_stay_integers(tmp_path: Path) -> None:
    """Guards the reason polars was chosen over pandas: nullable ints must not become floats."""
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:2\t200\tb\t\n",
    )
    frame, issues = read_table(path, SCHEMA)
    assert issues == []
    assert frame is not None
    assert frame["pos"].dtype == pl.Int64
    assert frame["pos"].to_list() == [100, 200]


def test_registry_covers_every_mirror_table() -> None:
    assert set(TABLE_SCHEMAS) == {
        "genes",
        "ptm_sites",
        "variants",
        "expression",
        "profiles",
        "cell_type_profiles",
        "profile_quantiles",
        "proteomics",
        "phospho",
        "clingen_validity",
        "gencc_submissions",
        "burden",
        "burden_combined",
        "orphanet_phenotypes",
        "orphanet_omim",
        "mondo_xrefs",
    }


def test_a_combined_row_is_keyed_on_the_burden_partition_minus_what_its_arms_span() -> None:
    """The combined table's key is a contract, not a convenience.

    A combination spans its arms, so it cannot carry the columns the arms
    differ in: `origin` (de novo vs transmitted), `maf_max` (none vs 1e-5) and
    `variant_class`. Adding any of them back would make one combination two rows.
    Dropping `lesion_group` would leave the paper's per-phenotype meta-analysis
    needing a migration. And `pvalue_adjusted` must never appear: the JL-FDR is
    a local FDR, smaller than its own p for about half the panel, which BUR018
    exists to refuse under that name.
    """
    assert BURDEN_COMBINED.sort_key == (
        "gene",
        "study",
        "cohort_stratum",
        "lesion_group",
        "consequence_class",
    )
    assert set(BURDEN.sort_key) - set(BURDEN_COMBINED.sort_key) == {
        "variant_class",
        "origin",
        "maf_max",
    }
    assert "local_fdr" in BURDEN_COMBINED.column_names
    assert "pvalue_adjusted" not in BURDEN_COMBINED.column_names


def test_clingen_validity_is_keyed_on_the_triple_not_the_pair() -> None:
    """(gene, disease, moi) is the key; (gene, disease) is not.

    Measured 2026-08-03 against the ClinGen bulk CSV: the triple is unique at
    3,653/3,653 rows, while (gene, disease) collapses 59 of them -- pairs curated
    twice under different modes of inheritance with different classifications.
    Sorting on the pair would make those rows compare equal and let a real
    ordering violation pass; more importantly it documents which columns a
    consumer must key on to avoid dropping one of each duplicate.
    """
    assert CLINGEN_VALIDITY.sort_key == ("gene", "disease", "moi")
    assert "sop" in CLINGEN_VALIDITY.column_names
    sop = next(c for c in CLINGEN_VALIDITY.columns if c.name == "sop")
    # Not nullable: 22% of the file predates SOP8 and ClinGen publishes no
    # crosswalk between framework versions, so a classification without the
    # framework that produced it is an unqualified claim.
    assert sop.nullable is False


def test_gencc_keys_on_the_submitter_because_it_publishes_no_verdict() -> None:
    """GenCC aggregates; it does not adjudicate.

    Its own terms: "The GenCC does not independently verify the submitted
    information." Two submitters may disagree about the same gene-disease pair
    and GenCC publishes both. Omitting `submitter` from the key would make those
    rows compare equal and silently keep whichever sorted last.

    `submitter` alone is still not enough. Measured 2026-08-03 against the real
    export (30,410 rows, ?format=new): (gene, disease, moi, submitter) collides
    on 133 groups (134 extra rows), and 70 of those 133 carry more than one
    distinct `classification_title` -- a submitter disagreeing with *itself*,
    e.g. Ambry Genetics submitting both `Limited` and `Strong` for HGNC:20226 /
    MONDO:0859332 / Autosomal recessive (see
    `test_a_submitter_can_disagree_with_itself_and_the_mirror_keeps_both`
    below). `version_number` is identical within every one of the 133 groups,
    so GenCC marks neither row as superseding the other -- there is no field in
    the export that says which is current. Picking one would synthesise a
    verdict GenCC itself does not assert, so the mirror keys on GenCC's own row
    id (`sgc_id`, unique across all 30,410 rows) instead of resolving the
    collision.

    `Supportive` is in the allowed set but is NOT a rung on the evidence ladder:
    it is a mapping-exception bucket for submitters that do not grade evidence.
    `vocab.GENCC_CLASSIFICATIONS` (Task 5) maps it to None.
    """
    assert GENCC_SUBMISSIONS.sort_key == ("gene", "disease", "moi", "submitter", "sgc_id")
    classification = next(c for c in GENCC_SUBMISSIONS.columns if c.name == "classification")
    assert classification.allowed is not None
    assert "Supportive" in classification.allowed


def test_a_submitter_can_disagree_with_itself_and_the_mirror_keeps_both() -> None:
    """The guard against "tidying up" the 133 real duplicate-looking groups.

    Ambry Genetics submitted both `Limited` (SGC-102815) and `Strong`
    (SGC-104042) for KIF26A / cortical dysplasia (HGNC:20226 / MONDO:0859332)
    under the same mode of inheritance -- measured directly in the committed
    snapshot, not synthesised for this test. A future edit that deduplicates
    `mirrors/gencc_submissions.tsv` on (gene, disease, moi, submitter) -- the
    key GenCC's own `version_number` does not distinguish within -- would
    collapse this to one classification and silently drop the disagreement
    this table exists to surface.
    """
    path = Path(__file__).parent.parent.parent / "mirrors" / "gencc_submissions.tsv"
    frame = pl.read_csv(path, separator="\t")
    rows = frame.filter(
        (pl.col("gene") == "HGNC:20226") & (pl.col("disease") == "MONDO:0859332")
    ).filter(pl.col("submitter") == "Ambry Genetics")

    assert set(rows["sgc_id"].to_list()) == {"SGC-102815", "SGC-104042"}
    assert set(rows["classification"].to_list()) == {"Limited", "Strong"}


def test_phospho_protein_normalized_is_mandatory_and_non_nullable() -> None:
    column = next(c for c in TABLE_SCHEMAS["phospho"].columns if c.name == "protein_normalized")
    assert column.nullable is False
    assert column.dtype == pl.Boolean


def test_unreadable_file_is_reported_not_raised(tmp_path: Path) -> None:
    """One bad byte must not abort validation of every other shard."""
    path = tmp_path / "demo.tsv"
    path.write_bytes("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\tFran\xe7ois\n".encode("latin-1"))

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_empty_file_is_reported_not_raised(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_bytes(b"")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_column_wide_failure_collapses_to_one_issue(tmp_path: Path) -> None:
    rows = "".join(f"HGNC:{n}\t{n}\tz\t\n" for n in range(1, 51))
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\n" + rows)

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL004"]
    assert "50 rows" in issues[0].message


def test_few_failures_are_reported_per_row(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t1\tz\t\nHGNC:2\t2\tz\t\n")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL004", "TBL004"]
    assert all("row" in i.location for i in issues)


def test_reports_a_value_outside_its_numeric_bounds(tmp_path: Path) -> None:
    schema = TableSchema(
        name="bounded",
        columns=(Column("p", pl.Float64, minimum=0, maximum=1),),
        sort_key=("p",),
    )
    path = tmp_path / "bounded.tsv"
    path.write_text("p\n5.0\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL006"]


@pytest.mark.parametrize("literal", ["NaN", "nan", "inf", "-inf", "Infinity"])
def test_reports_a_non_finite_value_in_a_bounded_column(tmp_path: Path, literal: str) -> None:
    """A missing statistic written as NaN must not validate as a present one.

    R's `write.table` emits `NaN` and numpy's `savetxt` emits `nan`; both parse
    into a Float64 column without a read error. Every comparison against NaN is
    False, so the bounds check cannot fire for it, and polars models NaN and null
    as distinct states, so the non-null check does not see it either.
    """
    schema = TableSchema(
        name="bounded",
        columns=(Column("p", pl.Float64, minimum=0, maximum=1),),
        sort_key=("p",),
    )
    path = tmp_path / "bounded.tsv"
    path.write_text(f"p\n{literal}\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]
    assert "row 2" in issues[0].location


def test_reports_a_non_finite_value_in_an_unbounded_column(tmp_path: Path) -> None:
    """The check belongs to the dtype, not to the bounds.

    A log2 fold change has no natural range, so it declares no minimum or
    maximum. It is still a statistic that must be present to be reported.
    """
    schema = TableSchema(
        name="unbounded",
        columns=(Column("log2fc", pl.Float64),),
        sort_key=("log2fc",),
    )
    path = tmp_path / "unbounded.tsv"
    path.write_text("log2fc\nNaN\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]


def test_a_nullable_float_column_still_rejects_nan(tmp_path: Path) -> None:
    """Nullable means "may be absent", which the empty string already says.

    NaN is a third state that reads as present, so it stays an error even where
    a null is allowed.
    """
    schema = TableSchema(
        name="optional",
        columns=(Column("score", pl.Float64, nullable=True),),
        sort_key=("score",),
    )
    path = tmp_path / "optional.tsv"
    path.write_text("score\n\nNaN\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]
    assert "row 3" in issues[0].location


def test_variants_chrom_rejects_a_chr_prefix() -> None:
    column = next(c for c in TABLE_SCHEMAS["variants"].columns if c.name == "chrom")
    assert column.allowed is not None
    assert "1" in column.allowed
    assert "chr1" not in column.allowed


def test_profiles_admits_rpkm_and_still_refuses_an_unknown_unit() -> None:
    """`unit` admits `rpkm` so a source publishing it need not be converted.

    **This layer does not publish RPKM, and the distinction matters.** The
    paper (Cardoso-Moreira et al. 2019) reports RPKM; the mirror this atlas
    actually carries is the EMBL-EBI Expression Atlas *reprocessing* of the
    same experiment, which requantified everything and publishes TPM and FPKM.
    `mirrors/profiles/E-MTAB-6814.tsv` is 100% `tpm` — measured, 18,326 of
    18,326 rows — and `scripts/convert_cardoso_moreira.py` sets `UNIT = "tpm"`.
    An earlier revision of this docstring said "the source publishes RPKM",
    which was true of the paper and false of the mirror, on a test shipped in
    the same PR as the TPM converter. Raised by review on PR #39.

    `rpkm` stays in the vocabulary regardless: `unit` exists so two
    incomparable quantities cannot merge, and the next source to publish RPKM
    must be able to say so rather than be converted. The negative half is the
    half that matters: a typo must still fail.
    """
    column = next(c for c in TABLE_SCHEMAS["profiles"].columns if c.name == "unit")
    allowed = column.allowed
    assert allowed is not None
    assert "rpkm" in allowed
    assert allowed == frozenset({"tpm", "nx", "cpm", "lfq", "rpkm"})
    assert "rpkms" not in allowed


def test_profile_quantiles_is_registered_and_sharded_per_dataset() -> None:
    """The table that makes the percentile auditable must itself be reachable.

    The two registries fail differently, not equivalently. Forget
    `SHARDED_TABLES` and the file is both unread (`mirror_paths` only globs
    directories it names) and reported as a stray (`unexpected_mirror_entries`
    reports TBL009 for a directory no schema claims). Forget `TABLE_SCHEMAS`
    instead and neither function notices -- both read `SHARDED_TABLES` only --
    so the file is found and accepted, and validation crashes instead: the
    first unguarded `TABLE_SCHEMAS[schema_name]` lookup downstream (the
    per-mirror loop in `validate_repository` reaches one first) raises an
    uncaught `KeyError` rather than reporting an issue.
    """
    schema = TABLE_SCHEMAS["profile_quantiles"]
    assert schema.column_names == (
        "dataset",
        "tissue",
        "stage",
        "percentile",
        "value",
        "unit",
        "n_genes",
    )
    assert schema.sort_key == ("dataset", "tissue", "stage", "percentile")
    assert SHARDED_TABLES["profile_quantiles"] == "profile_quantiles"
    assert "profile_quantiles" not in FLAT_TABLES

    percentile = next(c for c in schema.columns if c.name == "percentile")
    assert (percentile.minimum, percentile.maximum) == (0, 100)
    # The unit vocabulary must match `profiles` exactly, or PRF001 is comparing
    # two different alphabets rather than two values. Identity, not just
    # equality: a fresh, value-equal frozenset literal would pass an `==`
    # check today and still leave the two tables free to drift apart later.
    unit = next(c for c in schema.columns if c.name == "unit")
    profiles_unit = next(c for c in TABLE_SCHEMAS["profiles"].columns if c.name == "unit")
    assert unit.allowed == profiles_unit.allowed
    assert unit.allowed is profiles_unit.allowed


def test_mirror_paths_finds_flat_and_sharded_tables(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "variants" / "12.tsv").write_text("vrs_id\n")

    found = mirror_paths(tmp_path)

    assert [(p.name, schema) for p, schema in found] == [
        ("genes.tsv", "genes"),
        ("12.tsv", "variants"),
    ]


@pytest.mark.parametrize("name", sorted(TABLE_SCHEMAS))
def test_every_real_schema_accepts_a_header_only_file(name: str, tmp_path: Path) -> None:
    """A header-only shard must validate cleanly, and pins each schema's column list."""
    schema = TABLE_SCHEMAS[name]
    path = tmp_path / f"{name}.tsv"
    path.write_text("\t".join(schema.column_names) + "\n")

    assert validate_table(path, schema) == []


def test_reports_a_duplicate_sort_key(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:1\t100\tb\t\nHGNC:2\t200\ta\t\n")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL007"]
    assert "row 3" in issues[0].location


def test_distinct_sort_keys_are_not_reported_as_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:1\t200\ta\t\n")

    assert validate_table(path, SCHEMA) == []


def test_missing_file_is_reported_not_raised(tmp_path: Path) -> None:
    issues = validate_table(tmp_path / "absent.tsv", SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_unexpected_mirror_entries_reports_a_directory_no_schema_claims(tmp_path: Path) -> None:
    """`mirrors/variants` renamed to `mirrors/varaints` glob-matches nothing at all."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "unexpected directory" in issues[0].message
    assert issues[0].location.endswith("varaints")


def test_unexpected_mirror_entries_reports_a_shard_directory_that_is_a_file(
    tmp_path: Path,
) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants").write_text("oops\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "should be a directory of shards" in issues[0].message


def test_unexpected_mirror_entries_reports_a_file_no_schema_claims(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "notes.txt").write_text("stray\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "unexpected file" in issues[0].message


def test_unexpected_mirror_entries_accepts_the_declared_layout(tmp_path: Path) -> None:
    """An absent shard directory is legitimate: nothing curated yet."""
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "ptm_sites.tsv").write_text("site_id\n")
    (tmp_path / "mirrors" / "sources.yaml").write_text("sources: []\n")

    assert unexpected_mirror_entries(tmp_path) == []


def test_unexpected_mirror_entries_ignores_a_missing_mirrors_directory(tmp_path: Path) -> None:
    assert unexpected_mirror_entries(tmp_path) == []


def test_a_shard_with_the_wrong_extension_is_reported(tmp_path: Path) -> None:
    """mirror_paths globs *.tsv, so 12.txt is invisible to every table check."""
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants" / "12.txt").write_text("vrs_id\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "12.txt" in issues[0].location


def test_a_correctly_named_shard_is_not_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants" / "12.tsv").write_text("vrs_id\n")

    assert unexpected_mirror_entries(tmp_path) == []


@pytest.mark.parametrize(
    "schema_name,column_name",
    [
        ("genes", "uniprot"),
        ("ptm_sites", "protein"),
        ("proteomics", "protein"),
        ("phospho", "protein"),
    ],
)
def test_every_protein_column_constrains_the_accession(schema_name: str, column_name: str) -> None:
    """A protein identifier column must not accept free text.

    `UniprotAccession` constrains the curated YAML side, but the mirror tables
    are where the bulk protein data lives, and a column with no pattern took
    anything at all — a gene symbol, a RefSeq accession, an empty placeholder.
    """
    column = next(c for c in TABLE_SCHEMAS[schema_name].columns if c.name == column_name)
    assert column.pattern is not None
    compiled = re.compile(column.pattern)
    assert compiled.fullmatch("Q99593")
    assert compiled.fullmatch("P12345-2")
    assert not compiled.fullmatch("TBX5")
    assert not compiled.fullmatch("NP_852259.1")


def test_the_chromosome_vocabulary_is_free_of_duplicates() -> None:
    """What deriving the schema from this tuple does not make impossible.

    `VARIANTS`'s `chrom` column takes its `allowed` set from `CHROMOSOMES`, so
    the two can no longer disagree about membership — the test that used to
    check that was retired rather than kept as decoration.

    A duplicate entry survives that refactor and fails quietly in two directions
    at once: `frozenset` silently collapses it, so the schema looks right, while
    `build/variants.py`'s `_ORDER` keeps the *last* index, so one chromosome
    sorts where its twin sat. Neither surfaces until a curator files a shard for
    that chromosome.
    """
    assert len(CHROMOSOMES) == len(set(CHROMOSOMES))
    # Karyotype order is the published contract: the variant index drives a
    # chromosome picker, and it is the one property no other test would notice
    # losing, since every shard-emitting test uses a handful of chromosomes.
    assert CHROMOSOMES[:3] == ("1", "2", "3")
    assert CHROMOSOMES[-3:] == ("X", "Y", "MT")


def test_a_variant_shard_named_for_no_chromosome_is_reported(tmp_path: Path) -> None:
    """The filename is the whole addressing scheme, and nothing else checked it.

    `mirrors/variants/<chrom>.tsv` becomes `variants/<chrom>.json.gz`, and a
    consumer turns a chromosome into a URL by that rule alone. `validate_table`
    checks the `chrom` column; the sweep here checked only the extension. So
    `chr12.tsv` validated at 0 errors and 0 warnings and then raised inside the
    build, giving a curator a traceback where a report belongs — measured before
    this rule existed.

    The vocabulary is the schema's own, derived from the same `CHROMOSOMES`
    tuple as the `chrom` column's `allowed` set, so this cannot come to disagree
    with the column it addresses.

    A correctly named shard is asserted alongside, because a rule that fired on
    every shard would be worse than none.
    """
    shards = tmp_path / "mirrors" / "variants"
    shards.mkdir(parents=True)
    (shards / "12.tsv").write_text("")
    (shards / "MT.tsv").write_text("")
    (shards / "chr12.tsv").write_text("")
    (shards / "12_part2.tsv").write_text("")

    issues = unexpected_mirror_entries(tmp_path)

    reported = sorted(Path(issue.location).name for issue in issues if issue.code == "TBL011")
    assert reported == ["12_part2.tsv", "chr12.tsv"]
    assert all("is not one of" in issue.message for issue in issues if issue.code == "TBL011")


# --- Equivalence guard -------------------------------------------------------
#
# `validate_table` used to materialize every column into a Python list on every
# call, even when nothing was wrong -- the normal case. `_slow_*` below is a
# line-for-line copy of that pre-vectorization code, kept here (never in
# `src/`) purely so the tests that follow can assert the fast, polars-expression
# version reports byte-identical `ValidationIssue` lists: same codes, same
# severities, same locations, same messages, same order. `_row_issues` is
# imported rather than copied, because its collapse-past-20 behaviour was not
# touched by the vectorization and is already covered by its own tests
# (`test_column_wide_failure_collapses_to_one_issue`); copying it here would
# test nothing extra while risking the copy drifting from the real one.


def _slow_null_row_indices(series: pl.Series) -> list[int]:
    return [index + 2 for index, value in enumerate(series.to_list()) if value is None]


def _slow_string_values(series: pl.Series) -> list[tuple[int, str]]:
    return [
        (index + 2, str(value)) for index, value in enumerate(series.to_list()) if value is not None
    ]


def _slow_numeric_values(series: pl.Series) -> list[tuple[int, float]]:
    return [
        (index + 2, float(value))
        for index, value in enumerate(series.to_list())
        if value is not None
    ]


def _slow_validate_table(path: Path, schema: TableSchema) -> list[ValidationIssue]:
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
            offenders_null = [(row, "empty") for row in _slow_null_row_indices(series)]
            issues.extend(
                _row_issues(
                    "TBL003", path, offenders_null, f"column '{column.name}' must not be empty"
                )
            )

        if column.dtype is pl.Float64:
            offenders_finite = [
                (row, repr(value))
                for row, value in _slow_numeric_values(series)
                if not math.isfinite(value)
            ]
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
            offenders_allowed = [
                (row, repr(value))
                for row, value in _slow_string_values(series)
                if value not in column.allowed
            ]
            issues.extend(
                _row_issues(
                    "TBL004",
                    path,
                    offenders_allowed,
                    f"column '{column.name}' has values outside {allowed_sorted}",
                )
            )

        if column.pattern is not None:
            compiled = re.compile(column.pattern)
            offenders_pattern = [
                (row, repr(value))
                for row, value in _slow_string_values(series)
                if not compiled.fullmatch(value)
            ]
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
            offenders_bounds = [
                (row, repr(value))
                for row, value in _slow_numeric_values(series)
                if math.isfinite(value) and (value < low or value > high)
            ]
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


@pytest.mark.parametrize(
    "path,schema_name",
    mirror_paths(REPO_ROOT),
    ids=[f"{schema}:{path.name}" for path, schema in mirror_paths(REPO_ROOT)],
)
def test_vectorized_checks_match_the_original_on_every_committed_mirror(
    path: Path, schema_name: str
) -> None:
    """Real curated/mirrored data holds zero offenders in every checked column --
    the case this vectorization optimizes for, and the case a hand-written
    synthetic fixture would not think to cover. `clingen_gene_validity.tsv`
    (3,653 rows) and `gencc_submissions.tsv` (30,410 rows) are the tables that
    exposed the regression this change fixes.
    """
    schema = TABLE_SCHEMAS[schema_name]
    assert validate_table(path, schema) == _slow_validate_table(path, schema)


_OFFENDER_SCHEMA = TableSchema(
    name="offenders",
    columns=(
        Column("gene", pl.String, pattern=r"^HGNC:\d+$"),
        Column("tier", pl.String, allowed=frozenset({"a", "b"})),
        Column("score", pl.Float64),
        Column("n", pl.Int64, minimum=0, maximum=10),
    ),
    sort_key=("gene",),
)


def test_vectorized_checks_report_one_offender_of_each_kind_at_distinct_rows(
    tmp_path: Path,
) -> None:
    """One planted offender per check this task vectorized (TBL003, TBL004,
    TBL010, TBL005) plus TBL006 (also converted, though not required by the
    task, since it shared the same per-row materialization this change
    removes), each on a different, non-first data row.

    An off-by-one in a mask-based rewrite -- e.g. always reporting position 0,
    or reporting an `arg_true()` index instead of a row position -- would
    either misreport the row number or drop the offender entirely, and a
    fixture where every offender sat on row 2 would leave exactly that bug
    invisible. Row order here is asserted too, and on two different axes that
    must not be conflated:

    - *Within* one column, `validate_table` runs whichever of its checks that
      column's schema entry declares, always in the order TBL003, then
      TBL010, then TBL004, then TBL005, then TBL006. No column in this
      fixture reaches more than two of them (`gene`: TBL003, TBL005; `tier`:
      TBL004; `score`: TBL010; `n`: TBL006), so this axis alone does not
      explain the asserted order below -- it only ever arbitrates two checks
      landing on the *same* column.
    - *Across* columns, the outer loop visits `schema.columns` in schema
      order -- gene, tier, score, n -- before row numbers enter into it at
      all, so an offender is grouped with its column rather than sorted by
      which row it sits on. That is the axis this fixture actually exercises:
      `gene`'s TBL005 (row 6) is reported ahead of `tier`'s TBL004 (row 4) and
      `score`'s TBL010 (row 5) -- both on numerically *earlier* rows -- solely
      because `gene` is the first column in the schema.

    Measured by running this exact fixture, not inferred: guessing this order
    from the source got it wrong once already while drafting this test.

    Also the equivalence proof for the offender case: `validate_table`
    (vectorized) and `_slow_validate_table` (the original per-row Python
    loops) must produce the identical list, in the identical order.
    """
    path = tmp_path / "offenders.tsv"
    path.write_text(
        "gene\ttier\tscore\tn\n"
        "HGNC:1\ta\t0.5\t5\n"  # row 2: clean
        "\tb\t0.25\t3\n"  # row 3: TBL003, empty gene
        "HGNC:2\tz\t0.75\t2\n"  # row 4: TBL004, tier outside {a, b}
        "HGNC:3\ta\tNaN\t1\n"  # row 5: TBL010, non-finite score
        "BAD_ID\ta\t1.0\t4\n"  # row 6: TBL005, gene fails the pattern
        "HGNC:4\ta\t0.1\t999\n"  # row 7: TBL006, n outside [0, 10]
    )

    issues = validate_table(path, _OFFENDER_SCHEMA)

    assert [(i.code, i.location) for i in issues] == [
        ("TBL003", f"{path}:row 3"),
        ("TBL005", f"{path}:row 6"),
        ("TBL004", f"{path}:row 4"),
        ("TBL010", f"{path}:row 5"),
        ("TBL006", f"{path}:row 7"),
    ]
    # An Int64 offender must still report as a float, matching the original
    # `_numeric_values`'s `float(value)` conversion -- "999", not "999.0",
    # would be a silent message change a cast-free rewrite could introduce.
    assert "999.0" in issues[4].message
    assert issues == _slow_validate_table(path, _OFFENDER_SCHEMA)


_COLLAPSE_SCHEMA = TableSchema(
    name="collapse",
    columns=(
        Column("gene", pl.String, pattern=r"^HGNC:\d+$"),
        Column("tier", pl.String, allowed=frozenset({"a", "b"})),
    ),
    sort_key=("gene",),
)


def test_the_20_offender_collapse_still_names_the_correct_first_row(tmp_path: Path) -> None:
    """`_MAX_ROW_ISSUES` collapses past 20 offenders into one issue naming the
    first offending row -- proven here with the first offender NOT on row 2,
    so a mask-based rewrite that always reported position 0 (right for every
    existing collapse test, all of which start their offenders on row 2)
    would be caught.
    """
    path = tmp_path / "demo.tsv"
    rows = "HGNC:1\ta\n" + "".join(f"HGNC:{n}\tz\n" for n in range(2, 27))
    path.write_text("gene\ttier\n" + rows)

    issues = validate_table(path, _COLLAPSE_SCHEMA)

    assert [i.code for i in issues] == ["TBL004"]
    assert issues[0].message == (
        "25 rows: column 'tier' has values outside ['a', 'b'] (first at row 3: 'z')"
    )
    assert issues == _slow_validate_table(path, _COLLAPSE_SCHEMA)


def test_polars_pattern_engine_agrees_with_python_re(tmp_path: Path) -> None:
    """The dangerous vectorization: polars' `str.contains` runs the Rust `regex`
    crate, not Python's `re` -- a different implementation that could disagree
    on some input `re.fullmatch` would accept or reject differently.

    Proven two ways, both required before `tables.py` was allowed to rely on
    `str.contains`:

    1. Every non-null value of every pattern-bearing column in every committed
       mirror table -- 104,024 values total, counted per table beside the
       assertion below -- produces
       the identical offender set (here, the empty set: curated data is clean)
       under both engines, for every pattern actually used in a schema:
       `HGNC_PATTERN`, `MONDO_PATTERN`, `UNIPROT_PATTERN` and the inline
       `gencc_submissions.sgc_id` pattern are exercised by the committed data
       directly -- `UNIPROT_PATTERN` since `mirrors/genes.tsv` was widened from
       TBX5 alone to all 154 in-scope genes, every one of which carries an
       accession. `SEQUENCE_ONTOLOGY_PATTERN` and `MODIFICATION_PATTERN` still
       have no non-null values in the committed corpus (no variants or PTM
       sites are curated yet), so they are checked only against the adversarial
       set below.
    2. An adversarial synthetic set -- trailing/leading/embedded newlines
       (`re.fullmatch`'s `$` allows a trailing newline that `str.contains`'s
       `$` does not; `fullmatch` requires consuming the whole string regardless,
       so this was measured to make no difference), case, leading zeros,
       full-width and Arabic-Indic Unicode digits (both engines' `\\d` are
       Unicode-aware), stray `^`/`$` characters, and a 1,000-digit value --
       run against every pattern. `PMID_PATTERN` and `COHORT_LIST_PATTERN`
       joined the set on 2026-08-04 with `mirrors/burden.tsv`, and brought the
       `;` separator cases with them: `COHORT_LIST_PATTERN` is the first
       pattern here with alternation under repetition, so the empty token a
       leading, trailing or doubled `;` produces is a disagreement this set
       would otherwise never have looked for.

    Every pattern here already carries `^` and `$` (`re.fullmatch` anchors both
    ends implicitly; `str.contains` does not), so no pattern needed the Python
    fallback the task allowed for one that disagreed.
    """
    patterns = {
        column.pattern
        for schema in TABLE_SCHEMAS.values()
        for column in schema.columns
        if column.pattern is not None
    }
    assert patterns == {
        r"^HGNC:\d+$",
        r"^MONDO:\d{7}$",
        r"^SO:\d{7}$",
        r"^MOD:\d{5}$",
        r"^SGC-\d+$",
        r"^\d{4}-\d{2}-\d{2}.*$",
        r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$",
        r"^PMID:\d+$",
        r"^[a-z0-9]+(_[a-z0-9]+)*(;[a-z0-9]+(_[a-z0-9]+)*)*$",
        # The Orphanet scope mirrors, 2026-10-25.
        r"^ORPHA:\d+$",
        r"^OMIM:\d+$",
        r"^HP:\d{7}$",
        r"^(ORPHA|OMIM):\d+$",
    }

    adversarial = [
        "HGNC:123",
        "HGNC:123\n",
        "HGNC:123\n\n",
        "\nHGNC:123",
        "HGNC:12\n3",
        "hgnc:123",
        "HGNC:123 ",
        " HGNC:123",
        "HGNC:0123",
        "HGNC:",
        "HGNC:123x",
        "MONDO:1234567",
        "MONDO:1234567\n",
        "SO:1234567",
        "SO:1234567\n",
        "MOD:12345",
        "MOD:12345\n",
        "Q99593",
        "Q99593\n",
        "Q99593-2",
        "Q99593-2\n",
        "A0A123",
        "SGC-102815",
        "SGC-102815\n",
        "PMID:8988165",
        "PMID:8988165\n",
        "pmid:8988165",
        "PMID:",
        # `COHORT_LIST_PATTERN` is the first pattern here with alternation under
        # repetition, which is where two regex engines are likeliest to differ.
        # The separator cases matter most: a trailing, leading or doubled `;`
        # each produces an empty token, and an engine that let one through would
        # publish a cohort id no registry can resolve.
        "ddd",
        "ddd;ukbb",
        "cnchd;ddd;nottingham",
        "ukbb_500k",
        "ddd;",
        ";ddd",
        "ddd;;ukbb",
        "DDD",
        "ddd ukbb",
        "ddd;ukbb\n",
        "_ddd",
        "ddd_",
        "",
        "\n",
        "HGNC:123$",
        "HGNC:123^",
        "HGNC:1٣٣",  # Arabic-Indic digits mixed into an ASCII-digit run
        "HGNC:１２３",  # full-width digits
        "HGNC:" + "9" * 1000,
        # The Orphanet scope patterns. `Orphanet:648` is MONDO's spelling, which
        # the converter must normalise; `XREF_TARGET_PATTERN` is the second
        # pattern here with alternation, so its two arms are probed separately.
        "ORPHA:648",
        "ORPHA:648\n",
        "ORPHA:",
        "Orphanet:648",
        "OMIM:163950",
        "OMIM:163950\n",
        "HP:0001627",
        "HP:001627",
        "HP:0001627\n",
        "ORPHAOMIM:1",
        # `gencc_submissions.submitted_on`: a date prefix, then anything.
        "2018-05-31 20:04:32",
        "2025-01-17T00:00:00.000000Z",
        "2018-05-31",
        "2018-5-31",
        " 2018-05-31",
        "x2018-05-31",
        "2018-05-31\n",
        "2018-05-31\nz",
        "２０１８-05-31",
    ]
    adversarial_series = pl.Series("v", adversarial, dtype=pl.String)

    total_compared = 0
    for path, schema_name in mirror_paths(REPO_ROOT):
        schema = TABLE_SCHEMAS[schema_name]
        frame, read_issues = read_table(path, schema)
        assert frame is not None, read_issues
        for column in schema.columns:
            if column.pattern is None or column.name not in frame.columns:
                continue
            series = frame[column.name]
            compiled = re.compile(column.pattern)
            values = series.to_list()
            py_offenders = {
                i for i, v in enumerate(values) if v is not None and not compiled.fullmatch(v)
            }
            pl_mask = series.str.contains(column.pattern).not_()
            pl_offenders = set(pl_mask.arg_true().to_list()) if pl_mask.any() else set()
            assert py_offenders == pl_offenders, (path, column.name, column.pattern)
            total_compared += sum(1 for v in values if v is not None)

    # 3,653 x 2 in `clingen_gene_validity.tsv`, 30,410 x 3 in
    # `gencc_submissions.tsv`, 154 x 2 in `genes.tsv`, and 4,425 + 1,385 in
    # `burden.tsv`: 1,475 rows x 3 always-present columns (`study`, `gene`,
    # `case_cohorts`) plus `control_cohorts`, which is non-null on only 1,385.
    #
    # **`control_cohorts` stopped being a full column on 2026-08-05**, and the
    # comment here asserted that it was one until then: "every burden row across
    # both studies is case-control, so `control_cohorts` is non-null on all of
    # them". PMID:40127276 contributes the mirror's first `mutation_model` rows,
    # 90 of them, whose comparator is a mutability model and which therefore
    # carry no control cohort at all -- BUR001 refuses one that does. The
    # arithmetic is what caught the stale sentence: 1,475 x 4 is 5,900, and the
    # measured figure is 5,810.
    #
    # The burden figure has now risen twice: from 1,192 x 4 when PMID:34324492
    # added 103 rows, and again here.
    #
    # **104,654 -> 122,980 on 2026-08-19**, when the first profiles mirror
    # landed: E-MTAB-6814 contributes 18,326 rows whose `gene` column is
    # pattern-checked, plus 12,019 quantile rows whose `unit` column is. This
    # figure is a census of the committed corpus, not a property of the code,
    # so it moves whenever a mirror does -- which is the point of pinning it:
    # a silent change in what gets pattern-checked is exactly what it catches.
    #
    # **122,980 -> 123,250 on 2026-10-24**, when `burden_combined.tsv` landed:
    # 90 rows x 3 pattern-checked columns (`study`, `gene`, `combines`) = 270.
    #
    # **123,250 -> 176,566 on 2026-10-25**, when the three Orphanet scope
    # mirrors landed, each with two pattern-checked columns: (7,902 phenotype
    # + 8,745 OMIM + 10,011 xref rows) x 2 = 53,316.
    #
    # **176,566 -> 178,014 the same day**, when the registry widened from 154
    # genes to 185 and the burden converters re-ran against it: 31 genes x 2
    # (`hgnc_id`, `uniprot`) + 330 burden rows x 3 + the 297 of them with a
    # control cohort + 33 combined rows x 3 = 1,448.
    #
    # **178,014 -> 181,703 the same day**, when `mirrors/profiles` was
    # regenerated for the 185-gene registry: 31 genes x 119 assay groups = 3,689
    # rows, each with its `gene` column pattern-checked. The quantile mirror is
    # computed over the whole transcriptome, not the registry, and came out
    # byte-identical.
    #
    # **181,703 -> 212,113**, when `gencc_submissions.submitted_on` gained its
    # date-prefix pattern: all 30,410 rows carry a value (awk count).
    assert total_compared == 212_113

    for pattern in patterns:
        compiled = re.compile(pattern)
        py_offenders = {i for i, v in enumerate(adversarial) if not compiled.fullmatch(v)}
        pl_mask = adversarial_series.str.contains(pattern).not_()
        pl_offenders = set(pl_mask.arg_true().to_list()) if pl_mask.any() else set()
        assert py_offenders == pl_offenders, pattern


def test_cell_type_profiles_refuses_an_unknown_region_a_unit_but_cpm_and_a_fraction_above_one(
    tmp_path: Path,
) -> None:
    """One row, three defects: region `LX` is outside the vocabulary, `tpm` is
    not `cpm`, and 1.2 is outside 0..1.

    `cpm` alone because the gene page writes "CPM" beside every value and on
    its axis whatever the row says; a `tpm` row would publish under the wrong
    unit. Registered and sharded per dataset, like `profiles`; `source` is what
    makes the SRC003 gate apply to it, so it must stay a column.
    """
    schema = TABLE_SCHEMAS["cell_type_profiles"]
    assert SHARDED_TABLES["cell_type_profiles"] == "cell_type_profiles"
    assert schema.column_names == (
        "dataset",
        "gene",
        "cell_type",
        "region",
        "median_abundance",
        "unit",
        "q25",
        "q75",
        "n_donors",
        "n_cells",
        "fraction_expressing",
        "source",
    )
    assert schema.sort_key == ("dataset", "gene", "cell_type", "region")
    header = "\t".join(schema.column_names)
    good = "hoc\tHGNC:1100\tcardiomyocyte\tLV\t2.5\tcpm\t1.0\t3.0\t4\t900\t0.4\tsrc"
    bad = "hoc\tHGNC:1100\tcardiomyocyte\tLX\t2.5\ttpm\t1.0\t3.0\t4\t900\t1.2\tsrc"
    path = _write(tmp_path / "hoc.tsv", f"{header}\n{good}\n")
    assert validate_table(path, schema) == []
    path = _write(tmp_path / "hoc.tsv", f"{header}\n{bad}\n")
    issues = validate_table(path, schema)
    named = sorted(
        (i.code, name)
        for i in issues
        for name in ("region", "unit", "fraction_expressing")
        if f"column '{name}'" in i.message
    )
    assert named == [("TBL004", "region"), ("TBL004", "unit"), ("TBL006", "fraction_expressing")]
    assert len(issues) == 3
