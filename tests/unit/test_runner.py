import json
import shutil
from collections.abc import Callable
from pathlib import Path

import pytest

from chd_atlas.build.runner import BuildRefused, build_site
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.validate.runner import ValidationReport, validate_repository
from tests.unit._stamp import stamped

REPO_ROOT = Path(__file__).parent.parent.parent

ASSERTION_YAML = """\
assertions:
  - id: CHDA:AST:0000001
    gene: HGNC:11604
    phenotypes: [HP:0001631]
    lesion_groups: [septal]
    inheritance: [AD]
    mechanism: haploinsufficiency
    syndromic: both
    curator: 0000-0002-1825-0097
    curated_on: 2026-07-01
    last_reviewed: 2026-07-15
    evidence:
      - publication: PMID:8988165
        evidence_class: genetic_case
        locator: {kind: supplementary, file: TableS3.xlsx, sheet: S3, row: 42}
        summary: De novo nonsense variant in a proband with ASD.
        strength: strong
"""

VALID_SOURCES_YAML = """sources:
  - id: clinvar
    name: ClinVar
    version: "2026-07"
    retrieved_on: 2026-07-15
    url: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
    licence: public-domain
    redistribution: permitted
"""


def _corpus_load_fails(root: Path) -> None:
    """No curation directory, but a source registry that loads cleanly."""
    (root / "mirrors").mkdir()
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)


def _source_registry_load_fails(root: Path) -> None:
    """A curation directory that loads cleanly, but no source registry."""
    (root / "curation").mkdir()


def _both_loads_fail(root: Path) -> None:
    """Neither the corpus nor the source registry is present."""


def _validity_mirrors_missing(root: Path) -> None:
    """A corpus, gene registry, ID registry and source registry that all load
    cleanly, but neither validity mirror is present.

    Everything else is made to load cleanly on purpose, so that TBL012 is the
    *only* error in the report — unlike the "source-registry-failed" case
    below, where SCP000 rides alongside SRC001 by fixture coincidence, not
    because TBL012 caused it. This isolates SCP000's own causing error.
    """
    (root / "curation").mkdir()
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)


def _orphanet_mirrors_missing(root: Path) -> None:
    """One Orphanet-scoped term, and none of the three mirrors its warrant is read from.

    Everything else loads cleanly -- the term's label agrees with the one GenCC
    row naming it, and the gene registry resolves that row's symbol -- so SCP006
    is the *only* error. Measured when written: the report is exactly SCP006
    (error) and SCP008 (warning). Without the error, SCP008 would arrive alone in
    a report `ok` calls clean, and every gene the term admits would publish on a
    warrant nobody checked.
    """
    (root / "curation").mkdir()
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "curation" / "chd_scope.yaml").write_text(
        "diseases:\n"
        "  - id: MONDO:0018997\n"
        "    label: Noonan syndrome\n"
        "    reason: Orphanet records abnormal pulmonary valve morphology at Frequent.\n"
        "    admitted_by: orphanet_phenotype_frequency\n"
        "    attributed_to: Orphanet\n"
        "    orphanet: ORPHA:648\n"
        "    admitted_on: 2026-10-25\n"
    )
    (root / "mirrors").mkdir()
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "genes.tsv").write_text(
        "hgnc_id\tsymbol\tname\taliases\tprev_symbols\tensembl_gene\tncbi_gene\tlocus\t"
        "uniprot\tmane_select\n"
        "HGNC:9644\tPTPN11\tprotein tyrosine phosphatase non-receptor type 11\t\t\t"
        "ENSG00000179295\t5781\t12q24.13\tQ06124\tNM_002834.5\n"
    )
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(
        "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\t"
        "submitter\tsubmitted_on\treport_url\n"
        "SGC-1\tHGNC:9644\tPTPN11\tMONDO:0018997\tNoonan syndrome\tAD\tDefinitive\tG2P\t"
        "2023-01-01\thttps://x\n"
    )


def _orphanet_mirror_missing_a_column(root: Path) -> None:
    """The three Orphanet mirrors present, and one of them missing a column the rule reads.

    `load_orphanet_tables` must read that as unreadable -- None, so SCP006 and
    SCP008 name it once -- rather than index a column that is not there.
    Measured 2026-10-26: this fixture reports exactly SCP006 and TBL001
    (errors) and SCP008 (warning); with the column check removed,
    `validate_repository` raises `KeyError: 'hpo_label'` instead of reporting,
    and validators report, never raise (CLAUDE.md section 1). The other two
    mirrors are whole, so nothing but the column check can return None; and
    their sources are registered, so TBL001 is the only error besides the
    SCP006 it causes.
    """
    _orphanet_mirrors_missing(root)
    with (root / "mirrors" / "sources.yaml").open("a") as registry:
        for source in ("orphanet_product4", "orphanet_product1", "mondo_sssom"):
            registry.write(
                f"  - id: {source}\n"
                f"    name: {source}\n"
                '    version: "2026-10"\n'
                "    retrieved_on: 2026-10-25\n"
                "    url: https://example.org/\n"
                "    licence: CC-BY-4.0\n"
                "    redistribution: permitted\n"
            )
    (root / "mirrors" / "orphanet_phenotypes.tsv").write_text(
        "orpha\tdisorder_name\thpo\tfrequency\tsource\n"
        "ORPHA:648\tNoonan syndrome\tHP:0001631\tFrequent (79-30%)\torphanet_product4\n"
    )
    (root / "mirrors" / "orphanet_omim.tsv").write_text(
        "orpha\tomim\trelation\tvalidation\tsource\n"
        "ORPHA:648\tOMIM:163950\tE\tValidated\torphanet_product1\n"
    )
    (root / "mirrors" / "mondo_xrefs.tsv").write_text(
        "mondo\ttarget\tsource\nMONDO:0018997\tORPHA:648\tmondo_sssom\n"
    )


def _gene_registry_missing(root: Path) -> None:
    """Validity mirrors present, gene registry absent. Isolates GEN000's cause.

    The inverse of `_validity_mirrors_missing`: here there *is* symbol work to
    do -- both validity mirrors print a `gene_symbol` for every row -- and
    nothing to resolve it against, which is the only state in which GEN000 is
    supposed to fire. Both mirrors carry a header alone, which is enough for
    `_mirror_gene_labels` to read the columns and find no rows... so a single
    real row is written into ClinGen's, because a fixture with nothing to check
    would take the "no work, no skip" branch and pin the opposite of what this
    case is for.
    """
    (root / "curation").mkdir()
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\tDefinitive\t"
        "2023-01-01\tCHD\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(
        "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\t"
        "submitter\tsubmitted_on\treport_url\n"
    )


_GENES_TSV_HEADER = (
    "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
)
_PROFILES_TSV_HEADER = "dataset\tgene\ttissue\tstage\tmedian_abundance\tunit\tq25\tq75\tn_samples\n"
_EMPTY_CLINGEN_TSV = (
    "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
    "classification_date\tgcep\treport_url\n"
)
_EMPTY_GENCC_TSV = (
    "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\t"
    "submitter\tsubmitted_on\treport_url\n"
)
# A minimal, valid profile-design `Dataset` YAML, transcribed from
# `test_validate_profiles.py::_profile_dataset`'s field values so a profiles
# row naming this accession, tissue and stage does not also trip PRF004/005
# or REF005 -- noise this file's fixtures otherwise isolate away from.
_PROFILE_DATASET_YAML = """\
id: GSE999999
archive: geo
technology: bulk_rnaseq
design: profile
tissue: Heart
developmental_stage: embryonic
organism: NCBITaxon:9606
n_samples: 12
licence: CC BY 4.0
contrasts: []
cardiac_tissues: [Heart]
detection_floor: 1.0
floor_source: test fixture
quantile_estimator: linear
stages:
  - token: 7wpc
    wpc: 7.0
    order: 1
"""


_QUANTILES_TSV_HEADER = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"

# A minimal, valid phase vocabulary covering the fixture dataset's one 7 wpc
# stage. Present in every fixture that is *not* about the vocabulary being
# missing, so PRF013/PRF000 stay isolated to the fixtures that mean to trip
# them -- the same isolation `_PROFILE_DATASET_YAML` buys against PRF004/005.
# One phase, so `_prf006_issues` finds no interior gap to report either.
_CARDIAC_PHASES_YAML = """\
attributed_to: test fixture
citation: PMID:32048790
phases:
  - id: ventricular_septum_morphogenesis
    go_id: "GO:0060412"
    label: "Ventricular septum morphogenesis"
    start_wpc: 3.71
    end_wpc: 7.29
    start_carnegie_stage: "CS12"
    end_carnegie_stage: "CS19"
    start_hsapdv_id: "HsapDv:0000019"
    end_hsapdv_id: "HsapDv:0000026"
    end_basis: stated
"""

# The same dataset with its one stage post-natal. Phase matching is a no-op
# for a null `wpc` by construction (`Stage.wpc`), so this is a repository with
# profiles rows and *no work for the phase vocabulary to do*: section 4.41's
# case, where a skip must not be emitted.
_POSTNATAL_DATASET_YAML = _PROFILE_DATASET_YAML.replace(
    "  - token: 7wpc\n    wpc: 7.0\n    order: 1\n",
    "  - token: adult\n    wpc: null\n    order: 1\n",
)


def _repo_with_profiles_but_no_quantiles(root: Path) -> Path:
    """Profiles committed, quantiles not yet -- the ordinary state mid-curation.

    Absent is not unreadable, and only unreadable raises TBL000 (`mirror_paths`
    in `tables.py` yields a shard only when `path.is_file()`). Without a
    partner error, PRF000 would be a lone warning; `ValidationReport.ok`
    ignores warnings; the build would publish every gene with its percentile
    band silently missing and nothing on the site saying a figure is missing.

    Everything else here is made to load cleanly -- the gene registry (which
    also carries the one gene the profiles row cites, so PRF007 does not
    additionally fire), the curated dataset record (so neither PRF004/005 nor
    REF005 fires), the id registry, the source registry and both (empty)
    validity mirrors -- following `_gene_registry_missing`'s and
    `_validity_mirrors_missing`'s precedent, so PRF000/PRF010 are isolated to
    their own causing pair rather than riding along with unrelated noise.

    The cardiac phase vocabulary is here for that same isolation and is not
    incidental: this fixture's profiles row names a 7 wpc stage, so without a
    vocabulary it would trip PRF013 as well and stop isolating PRF010.
    """
    (root / "curation" / "datasets").mkdir(parents=True)
    (root / "curation" / "datasets" / "GSE999999.yaml").write_text(_PROFILE_DATASET_YAML)
    (root / "curation" / "cardiac_phases.yaml").write_text(_CARDIAC_PHASES_YAML)
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box transcription factor 5", "", "", "", "", "", ""])
        + "\n"
    )
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(_EMPTY_CLINGEN_TSV)
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)
    (root / "mirrors" / "profiles").mkdir()
    (root / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:11604", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
    )
    return root


def _quantiles_tsv(dataset: str, tissue: str, stage: str, unit: str = "tpm") -> str:
    """A complete 101-point breakpoint grid for one cell, non-decreasing.

    Written out in full rather than truncated: a short grid is a different
    defect with its own codes, and a fixture carrying one would leave the
    PRF013 cases below reporting two things at once.
    """
    rows = "".join(
        "\t".join([dataset, tissue, stage, str(percentile), str(float(percentile)), unit, "20"])
        + "\n"
        for percentile in range(101)
    )
    return _QUANTILES_TSV_HEADER + rows


def _repo_with_profiles_and_quantiles(root: Path, dataset_yaml: str, stage: str) -> Path:
    """Everything a profiles corpus needs except the cardiac phase vocabulary.

    Deliberately complete on every other axis -- gene registry, id registry,
    source registry, both validity mirrors, the curated dataset record, and a
    full quantile grid for the one cell the profiles row occupies -- so the
    phase-vocabulary codes below arrive alone rather than beside PRF010,
    PRF004/005 or TBL012. Same isolation discipline as
    `_repo_with_profiles_but_no_quantiles`, from which this is otherwise
    copied.
    """
    (root / "curation" / "datasets").mkdir(parents=True)
    (root / "curation" / "datasets" / "GSE999999.yaml").write_text(dataset_yaml)
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box transcription factor 5", "", "", "", "", "", ""])
        + "\n"
    )
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(_EMPTY_CLINGEN_TSV)
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)
    (root / "mirrors" / "profiles").mkdir()
    (root / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:11604", "Heart", stage, "10.0", "tpm", "", "", "5"])
        + "\n"
    )
    (root / "mirrors" / "profile_quantiles").mkdir()
    (root / "mirrors" / "profile_quantiles" / "GSE999999.tsv").write_text(
        _quantiles_tsv("GSE999999", "Heart", stage)
    )
    return root


def _repo_with_no_phase_vocabulary(root: Path) -> Path:
    """Prenatal profiles rows and no `curation/cardiac_phases.yaml` at all.

    The absent case. `corpus.cardiac_phases` is `None`, `_prf006_issues`
    returns `[]` on sight of it, and every check that reads the vocabulary is
    skipped in silence. Measured on the committed corpus before PRF013:
    deleting the file took `validate` from 4 warnings to **3**, and the build
    stayed green while 460 stage entries flipped from `matched` to
    `outside_window` and all 85 charted pages lost their phase bands.
    """
    return _repo_with_profiles_and_quantiles(root, _PROFILE_DATASET_YAML, "7wpc")


def _repo_with_an_empty_phase_vocabulary(root: Path) -> Path:
    """A well-formed vocabulary file declaring no phases -- equally silent.

    Not the same input as a *malformed* file: `phases: []` on its own fails
    `CardiacPhaseFile`'s required `attributed_to`/`citation` and is already an
    error (SCHEMA001), which is why the brief's shorthand does not reproduce.
    A file carrying both and an empty list parses cleanly, so `_prf006_issues`
    finds no gap in nothing and reports nothing, exactly as an absent file
    does. Both must reach PRF013; a check written as `phases is None` catches
    only the first.
    """
    root = _repo_with_profiles_and_quantiles(root, _PROFILE_DATASET_YAML, "7wpc")
    (root / "curation" / "cardiac_phases.yaml").write_text(
        "attributed_to: test fixture\ncitation: PMID:32048790\nphases: []\n"
    )
    return root


def _repo_with_only_post_natal_profiles(root: Path) -> Path:
    """Profiles rows, no vocabulary, and nothing for a vocabulary to do.

    Section 4.41: emit a skip only when there was work to skip. A stage with a
    null `wpc` is post-natal *by construction* (`Stage.wpc`) and falls outside
    every cardiac morphogenetic phase whether or not one is curated, so no
    phase assignment changes here — and PRF013 firing would be the lone-skip
    defect that section exists to prevent, in reverse.
    """
    return _repo_with_profiles_and_quantiles(root, _POSTNATAL_DATASET_YAML, "adult")


def _repo_with_no_profiles(root: Path) -> Path:
    """No `mirrors/profiles/` at all -- nothing here for PRF000 to check.

    Not simply an empty `tmp_path`: PRF000 must be provably absent when a
    `mirrors/` directory exists for other tables too, not merely when the
    whole repository happens to be empty by coincidence.
    """
    (root / "mirrors").mkdir(parents=True)
    return root


def test_report_counts_errors_and_warnings_separately() -> None:
    report = ValidationReport(
        issues=[
            ValidationIssue("A001", Severity.ERROR, "x", "boom"),
            ValidationIssue("B001", Severity.WARNING, "y", "hmm"),
        ]
    )
    assert report.error_count == 1
    assert report.warning_count == 1
    assert report.ok is False


def test_empty_report_is_ok() -> None:
    assert ValidationReport(issues=[]).ok is True


def test_issues_are_sorted_and_deduplicated() -> None:
    duplicate = ValidationIssue("A001", Severity.ERROR, "x", "boom")
    report = ValidationReport(
        issues=[
            ValidationIssue("B001", Severity.ERROR, "y", "later"),
            duplicate,
            duplicate,
        ]
    )
    assert [i.code for i in report.issues] == ["A001", "B001"]


def test_validate_repository_reports_a_missing_corpus(tmp_path: Path) -> None:
    report = validate_repository(tmp_path)
    codes = {issue.code for issue in report.issues}
    assert "CORPUS001" in codes
    assert "SRC001" in codes
    assert report.ok is False


def test_a_misnamed_mirror_directory_is_reported(tmp_path: Path) -> None:
    """A typo'd rename must not make the repository look clean."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert any(issue.code == "TBL009" for issue in report.issues)


def test_a_shard_directory_replaced_by_a_file_is_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants").write_text("oops\n")

    report = validate_repository(tmp_path)

    assert any(issue.code == "TBL009" for issue in report.issues)


def test_expected_mirror_entries_are_not_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "sources.yaml").write_text("sources: []\n")

    report = validate_repository(tmp_path)

    assert not any(issue.code == "TBL009" for issue in report.issues)


def test_a_missing_shard_directory_is_not_reported(tmp_path: Path) -> None:
    """No variants curated yet is a legitimate state, not an error."""
    (tmp_path / "mirrors").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert not any(issue.code == "TBL009" for issue in report.issues)


def test_ontology_checks_are_skipped_when_the_registry_fails(tmp_path: Path) -> None:
    """41 ONT003s would otherwise bury the one SRC001 that caused them."""
    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert "ONT000" in codes
    assert "ONT003" not in codes


def test_a_failed_registry_does_not_cascade_one_error_per_term(tmp_path: Path) -> None:
    """With terms to check, the unguarded run reports ONT003 once per term."""
    (tmp_path / "curation").mkdir()
    (tmp_path / "curation" / "phenotypes.yaml").write_text(
        "phenotypes:\n"
        "  - id: HP:0001631\n"
        "    label: Atrial septal defect\n"
        "    lesion_group: septal\n"
        "  - id: HP:0001629\n"
        "    label: Ventricular septal defect\n"
        "    lesion_group: septal\n"
    )

    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert codes.count("ONT000") == 1
    assert "ONT003" not in codes


def test_a_gap_warning_is_reported_without_blocking_the_build() -> None:
    """The second kind of warning, and why `ok` may ignore one that arrives alone.

    A *skip* warning means a check did not run, and is ignored safely only
    because an error always accompanies it — the three parametrised cases below
    pin that. REF013 is a *gap* warning: the check ran and found curated evidence
    that will reach no reader, and it deliberately arrives alone.

    That is a real weakening of the gate, not an oversight, and it is asserted
    here so that it is a decision on record. The assertion is sound and the site
    is publishable; what is missing is one cell in a mirror, and refusing a whole
    deploy over it is the wrong trade for a corpus curated incrementally. What
    makes it acceptable is the alternative: before REF013 the same gap was
    reported nowhere at all, on a build that said 0 errors and 0 warnings.

    The report is constructed rather than provoked because the claim is about
    `ok`'s semantics, not about when the validator fires — that is
    `test_referential.py`'s subject, where the rule lives.
    """
    gap = ValidationReport(
        issues=[
            ValidationIssue("REF013", Severity.WARNING, "mirrors/genes.tsv", "no uniprot accession")
        ]
    )

    assert gap.warning_count == 1
    assert gap.error_count == 0
    assert gap.ok is True, "a gap must not block a deploy; it must be visible in the report"


@pytest.mark.parametrize(
    ("setup", "expected_warnings", "causing_errors"),
    [
        pytest.param(_corpus_load_fails, {"REF000"}, {"CORPUS001"}, id="corpus-failed"),
        pytest.param(
            _source_registry_load_fails,
            {"SRC000", "ONT000", "SCP000", "SCP005"},
            {"SRC001"},
            id="source-registry-failed",
        ),
        pytest.param(
            _both_loads_fail,
            {"REF000", "SRC000", "ONT000"},
            {"CORPUS001", "SRC001"},
            id="both-failed",
        ),
        pytest.param(
            _validity_mirrors_missing,
            {"SCP000", "SCP005"},
            {"TBL012"},
            id="validity-mirrors-missing",
        ),
        # GEN000's own cause, isolated the same way TBL012 isolates SCP000's.
        # A missing gene registry is TBL008, and without that error this skip
        # would arrive alone in a report `ok` calls clean -- a repository whose
        # every symbol was unresolvable, validating green. Issue #33.
        pytest.param(
            _gene_registry_missing,
            {"GEN000"},
            {"TBL008"},
            id="gene-registry-missing",
        ),
        # PRF000's own cause, isolated the same way TBL008 isolates GEN000's.
        # Profiles committed, quantiles not yet -- the ordinary state
        # mid-curation -- and without the partner error PRF010, this skip
        # would arrive alone in a report `ok` calls clean: every gene loses
        # its percentile band, silently, on a green build.
        pytest.param(
            _repo_with_profiles_but_no_quantiles,
            {"PRF000"},
            {"PRF010"},
            id="profiles-quantiles-absent",
        ),
        # PRF013's own causing error, isolated: profiles rows that would have
        # been phase-matched and no vocabulary to match them against. Without
        # the error, PRF000 arrives alone in a report `ok` calls clean and the
        # build publishes 85 pages saying the *dataset's* stages fall in no
        # curated phase, when what actually happened is that the atlas's own
        # vocabulary never loaded.
        pytest.param(
            _repo_with_no_phase_vocabulary,
            {"PRF000"},
            {"PRF013"},
            id="phase-vocabulary-absent",
        ),
        # And the same for a file that parses cleanly and declares nothing.
        pytest.param(
            _repo_with_an_empty_phase_vocabulary,
            {"PRF000"},
            {"PRF013"},
            id="phase-vocabulary-empty",
        ),
        # SCP008's own cause: an Orphanet-scoped term with no Orphanet mirrors.
        pytest.param(
            _orphanet_mirrors_missing,
            {"SCP008"},
            {"SCP006"},
            id="orphanet-mirrors-missing",
        ),
        # And for a mirror that is present but missing a column: the same skip,
        # the same cause, and the table error naming the column.
        pytest.param(
            _orphanet_mirror_missing_a_column,
            {"SCP008"},
            {"SCP006", "TBL001"},
            id="orphanet-mirror-missing-a-column",
        ),
    ],
)
def test_every_skip_warning_arrives_with_the_error_that_caused_it(
    tmp_path: Path,
    setup: Callable[[Path], None],
    expected_warnings: set[str],
    causing_errors: set[str],
) -> None:
    """`ValidationReport.ok` ignores warnings, so a warnings-only report exits 0.

    That is only safe while every warning means "a check was skipped" and each
    arrives with the error that caused the skip. This pins both halves for all
    three skip paths: the exact set of warnings each state produces, and the
    error that explains them.

    SRC000 and ONT000 share one guard — a source registry that did not load —
    so no repository state triggers either alone; they are pinned as the pair
    they are. REF000 is isolable, and the third state covers all three at once.

    SCP000 also appears in the "source-registry-failed" case, but not because
    it shares SRC000/ONT000's guard: `_source_registry_load_fails` builds a
    `curation/` directory with no `mirrors/` at all, so the validity-mirror
    read (`_mirrored_validity`) independently returns None. It is a fixture
    coincidence, not a shared cause — `causing_errors` still names only SRC001
    there, and that case's `report.ok is False` holds because of SRC001
    regardless of SCP000.

    "validity-mirrors-missing" is what actually isolates SCP000's own cause:
    every other check is made to load cleanly (gene registry, ID registry,
    source registry all present and valid) and only the two validity mirrors
    are absent, so TBL012 is the *only* error in that report. Before TBL012
    existed, this exact state reported SCP000 as its only warning and
    `report.ok` was True — a corpus whose scope could not be checked at all,
    validating clean. Measured directly: with the TBL012 emission removed from
    `validate_repository`, this case's `report.ok is False` assertion fails
    (`assert True is False`).
    """
    setup(tmp_path)

    report = validate_repository(tmp_path)

    warnings = {issue.code for issue in report.issues if issue.severity is Severity.WARNING}
    errors = {issue.code for issue in report.issues if issue.severity is Severity.ERROR}
    assert warnings == expected_warnings
    assert causing_errors <= errors
    assert report.error_count > 0
    assert report.ok is False


def test_render_separates_the_summary_from_the_issues(tmp_path: Path) -> None:
    """Without a blank line the count runs straight on from the last issue."""
    report = validate_repository(tmp_path)

    lines = report.render().splitlines()
    assert lines[-1] == f"{report.error_count} error(s), {report.warning_count} warning(s)"
    assert lines[-2] == ""


def test_locations_are_repo_relative(tmp_path: Path) -> None:
    """Absolute paths would make CI output differ from a curator's machine."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert report.issues
    for issue in report.issues:
        assert not issue.location.startswith(str(tmp_path)), issue.location


def _seeded(tmp_path: Path) -> Path:
    """A corpus with one assertion, enough to make the guards below reachable."""
    (tmp_path / "curation" / "assertions").mkdir(parents=True)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(ASSERTION_YAML)
    return tmp_path


def test_a_missing_phenotype_vocabulary_is_an_error(tmp_path: Path) -> None:
    """Without it REF007/REF009/REF010 silently no-op. An error, not a warning,
    because ValidationReport.ok ignores warnings."""
    report = validate_repository(_seeded(tmp_path))

    assert "CUR002" in [issue.code for issue in report.issues]
    assert report.ok is False


def test_a_missing_publication_registry_is_an_error(tmp_path: Path) -> None:
    """Otherwise every citation reports REF002 and nothing names the cause."""
    report = validate_repository(_seeded(tmp_path))

    assert "CUR003" in [issue.code for issue in report.issues]


def test_ontology_checks_are_skipped_when_a_pinned_release_fails(tmp_path: Path) -> None:
    """ONT003 says the prefix is unpinned, which is false when the file is
    merely unreadable - it would send a curator to add a pin that exists."""
    root = _seeded(tmp_path)
    (root / "mirrors").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "sources.yaml").write_text(
        "sources:\n"
        "  - id: hpo\n"
        "    name: Human Phenotype Ontology\n"
        '    version: "2026-06-23"\n'
        "    retrieved_on: 2026-07-31\n"
        "    url: https://hpo.jax.org/\n"
        "    licence: custom\n"
        "    redistribution: permitted_with_attribution\n"
        "    ontology_prefix: HP\n"
        "    ontology_file: ontologies/broken.obo\n"
    )
    (root / "ontologies").mkdir()
    (root / "ontologies" / "broken.obo").write_text("not an ontology\n")

    report = validate_repository(root)

    codes = [issue.code for issue in report.issues]
    assert "ONT004" in codes
    assert "ONT000" in codes
    assert "ONT003" not in codes


def test_scope_checks_run_against_the_real_repository() -> None:
    """The committed repository is in scope and reports no scope error.

    A smoke test in the sense CLAUDE.md permits: it is the only thing that
    exercises the wiring from `validate_repository` down through the mirror
    read into `validate_scope_terms`, `scope_candidates` and
    `validate_curation_is_in_scope`.

    Asserting only `report.error_count == 0` is not enough here -- that
    passes identically whether the scope checks ran or were never called.
    Measured directly: with the scope-check call site in `validate_repository`
    reverted entirely, `.venv/bin/python -m pytest` on a naive version of this
    test asserting only `report.error_count == 0` still passed (1 passed) --
    it cannot tell a live wire from a dead one. So this test also asserts that
    a scope check actually produced output.

    SCP003 (a mirrored cardiac-looking disease term that is not on the scope
    list) is emitted whenever such a term exists, and 268 (net 1, label
    keyword) plus 223 (net 2, gene bridge) were measured on 2026-08-03 against
    the committed 68-term `curation/chd_scope.yaml` and the committed
    ClinGen/GenCC mirrors, so it fires today. If the scope list ever grows to
    cover every candidate, replace the `"SCP003" in codes` assertion with one
    that the mirror-reading helper (`validate/runner.py::_mirrored_validity`)
    returned a non-empty result -- do not delete the guard, or this test
    degrades back into the worthless version described above.

    SCP007 is asserted the same way and for the same reason. Its digest logic is
    unit-tested, but the 2026-10-25 mutation matrix deleted its one call site in
    `validate_repository` and the whole suite stayed green: nothing watched it
    fire from the real validator. On the committed mirrors it lists 248 terms
    that would qualify on Orphanet's frequencies if nominated (measured
    2026-10-26; 239 before R1 widened the structural set); Fabry disease is one
    (valve findings), and Noonan syndrome, which is nominated, is not. Noonan
    syndrome-like disorder with loose anagen hair qualifies only through
    pulmonic stenosis, one of the lesions R1 added, so it is the real-data
    witness that the widening reaches the validator.
    """
    report = validate_repository(REPO_ROOT)

    codes = {issue.code for issue in report.issues}
    assert "SCP003" in codes, "scope checks did not run"
    assert report.error_count == 0, report.render()

    digests = [issue for issue in report.issues if issue.code == "SCP007"]
    assert len(digests) == 1, "the Orphanet candidate digest did not run"
    assert "MONDO:0010526 (Fabry disease)" in digests[0].message
    assert "MONDO:0018997 (Noonan syndrome)" not in digests[0].message
    assert (
        "MONDO:0011899 (Noonan syndrome-like disorder with loose anagen hair) via ORPHA:2701"
        in digests[0].message
    )


def test_a_validity_mirror_missing_a_selected_column_is_reported_not_raised(
    tmp_path: Path,
) -> None:
    """Reported by Copilot on the promotion PR, reproduced here.

    `_mirrored_validity` guarded only `"disease" not in frame.columns`, then
    unconditionally did `frame.select(["disease", "disease_label", "gene"])`.
    An upstream rename dropping `disease_label` or `gene` -- not `disease`
    itself -- left `"disease" not in frame.columns` false, so the guard passed
    and `frame.select` raised `polars.exceptions.ColumnNotFoundError`,
    unhandled, straight out of `validate_repository`. That breaks the one
    guarantee this module exists to keep (`CLAUDE.md`: "validators report; the
    build refuses") -- one malformed mirror must cost one issue, not the
    whole run, and `validate_table` already reports exactly this as TBL001
    over the same file, so nothing is lost by skipping the mirror here rather
    than crashing on it.

    Watched failing first: with the fix in `_mirrored_validity` reverted to
    checking only `"disease" not in frame.columns`, this test raised
    `polars.exceptions.ColumnNotFoundError: unable to find column
    "disease_label"` instead of returning a report -- confirmed both via this
    test and via `chd-atlas validate` against a scratch copy of the real
    repository with the same column dropped.

    A copy of the real, valid committed corpus rather than a minimal fixture,
    so `disease_label` is the *only* thing wrong with it and every other
    check -- including the second mirror, GenCC, read successfully -- still
    runs, isolating this guard's own effect from `TBL012`'s ("both mirrors
    unreadable") different failure mode.
    """
    source = tmp_path / "repo"
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO_ROOT / name, source / name)
    (source / "ontologies").symlink_to(REPO_ROOT / "ontologies")

    clingen = source / "mirrors" / "clingen_gene_validity.tsv"
    header, *rows = clingen.read_text(encoding="utf-8").splitlines()
    columns = header.split("\t")
    label_index = columns.index("disease_label")
    del columns[label_index]
    lines = ["\t".join(columns)]
    for row in rows:
        cells = row.split("\t")
        del cells[label_index]
        lines.append("\t".join(cells))
    clingen.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = validate_repository(source)

    codes = [issue.code for issue in report.issues]
    assert "TBL001" in codes
    assert "TBL012" not in codes, "GenCC alone should still populate the mirrored-validity map"
    assert report.ok is False


def _real_repo(tmp_path: Path) -> Path:
    """A copy of the committed repository, which validates at the baseline.

    Copied rather than fixtured because the defect below is about *scale* -- one
    error naming a cause versus one error per row naming symptoms -- and a
    two-row fixture cannot show the difference.

    The copy drops `curation/literature/`: the committed corpus now carries
    literature files, and a fixture that writes its own literature file would
    otherwise hit `FileExistsError` on `mkdir`, or compare against a baseline that
    already holds the real file's issues. Every fixture that builds on this copy
    gets the same empty literature directory.
    """
    for name in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / name, tmp_path / name)
    shutil.rmtree(tmp_path / "curation" / "literature", ignore_errors=True)
    return tmp_path


def test_an_absent_cohort_registry_reports_its_cause_once_not_every_row(
    tmp_path: Path,
) -> None:
    """CUR004, and the guarantee in CLAUDE.md §1 it exists to keep.

    "A missing registry must not report as hundreds of dangling references."
    The first version of `validate_burden_references` passed
    `{str(c.id) for c in corpus.cohorts}` unconditionally, so an absent
    `cohorts.yaml` -- which `load_curation` treats as optional and reports
    nothing for -- produced four BUR009 errors naming four cohorts and nothing
    naming the file. Measured 2026-08-05 before the fix: 4 errors, all BUR009.
    """
    root = _real_repo(tmp_path)
    (root / "curation" / "cohorts.yaml").unlink()

    report = validate_repository(root)
    codes = [issue.code for issue in report.issues if issue.severity is Severity.ERROR]

    assert codes == ["CUR004"]
    assert report.ok is False


def test_an_unreadable_gene_registry_does_not_multiply_by_the_burden_mirror(
    tmp_path: Path,
) -> None:
    """The same cascade from the other side, and the one that grew 146x.

    A header-only `mirrors/genes.tsv` reads fine and carries the `hgnc_id`
    column, so `_known_genes` returns an empty set rather than None and TBL008
    does not fire. Before `mirrors/burden.tsv` existed that input produced
    **1** error (REF001, the one curated assertion's gene). Measured 2026-08-05
    with the burden table added and the registries passed unconditionally:
    **146** -- 1 REF001 plus one BUR011 per distinct gene in the mirror -- none
    of which named the registry.

    Passing `known_genes or None` restores it to 1. The count is asserted, not
    just the code, because the defect was never a wrong code: BUR011 was right
    about all 145 genes. It was the burying that was wrong.

    **Plus one GEN004 naming the empty registry** (review, 2026-10-27). GEN004
    was handed the registry `or None` for the same reason, which switched it
    off, and then the only error left was REF001 -- the curated assertion's.
    With that assertion moved out, as here second, the report was **0
    errors** and `build_site` raised from `external_links` after writing had
    begun (measured on a copy): the failure R13 exists to report from
    validate. Now one GEN004, and the build refuses before writing anything.
    """
    root = _real_repo(tmp_path / "repo")
    genes = root / "mirrors" / "genes.tsv"
    genes.write_text(genes.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")

    report = validate_repository(root)
    codes = [issue.code for issue in report.issues if issue.severity is Severity.ERROR]

    # One REF001 per curated assertion (every committed assertion's gene is now
    # unknown under a header-only registry) plus one GEN004 naming the empty
    # registry. Counted, not a literal list, so adding another committed
    # assertion (as NR2F2/CHDA:AST:0000002 did) does not stale this test.
    n_assertions = sum(1 for _ in (root / "curation" / "assertions").glob("*.yaml"))
    assert codes == ["REF001"] * n_assertions + ["GEN004"]

    for assertion in (root / "curation" / "assertions").iterdir():
        assertion.unlink()
    report = validate_repository(root)
    errors = [issue for issue in report.issues if issue.severity is Severity.ERROR]

    assert [issue.code for issue in errors] == ["GEN004"]
    assert "the gene registry holds no gene rows" in errors[0].message
    out = tmp_path / "dist"
    with pytest.raises(BuildRefused, match="GEN004"):
        build_site(root, out)
    assert not out.exists()


def test_a_false_scope_attribution_is_caught_on_the_branch_that_actually_runs(
    tmp_path: Path,
) -> None:
    """SCP005 through `validate_repository`, on the mirror-READABLE branch.

    **The branch that runs on every real build was killed by 0 of 789 tests.**
    `validate_repository` calls `validate_scope_attribution` twice -- once when
    the mirrors are unreadable, once when they are -- and only the first was
    defended. Deleting the second call left the whole suite green while
    `chd-atlas validate` reported 0 errors on a corpus with an invented
    attribution, and the atlas would publish 24 pages asserting an external
    warrant nothing had checked. Found by adversarial review of #30; it is the
    §4.28 shape, a guard tested on the layer it was not needed.

    The forgery used here is the one the check exists for: a term whose only
    warrant is a commercial laboratory's GenCC submission, relabelled as
    ClinGen's chartered CHD panel. The unit test proves the function refuses it;
    this proves the function is reached.
    """
    root = _real_repo(tmp_path)
    scope = root / "curation" / "chd_scope.yaml"
    text = scope.read_text(encoding="utf-8")
    forged = text.replace(
        "    admitted_by: gencc_submitter\n    attributed_to: Ambry Genetics\n",
        "    admitted_by: clingen_chd_panel\n    attributed_to: Ambry Genetics\n",
        1,
    )
    assert forged != text, "no commercial-lab attribution left to forge; update this fixture"
    scope.write_text(forged, encoding="utf-8")

    report = validate_repository(root)

    assert "SCP005" in {issue.code for issue in report.issues}
    assert report.ok is False


def test_an_absent_quantile_mirror_is_an_error_not_a_lone_warning(tmp_path: Path) -> None:
    """Absent is not unreadable, and only unreadable raises TBL000.

    Profiles committed, quantiles not yet -- the ordinary state mid-curation.
    Without a partner error PRF000 is a lone warning, `ok` ignores warnings,
    the build publishes, and every gene loses its percentile band with nothing
    on the site saying a figure is missing.
    """
    root = _repo_with_profiles_but_no_quantiles(tmp_path)
    report = validate_repository(root)
    codes = {issue.code for issue in report.issues}
    assert "PRF000" in codes
    assert "PRF010" in codes  # the partner error
    assert not report.ok


def test_no_profiles_means_no_skip_warning_at_all(tmp_path: Path) -> None:
    """Emit a skip only when there was work to skip."""
    report = validate_repository(_repo_with_no_profiles(tmp_path))
    assert "PRF000" not in {issue.code for issue in report.issues}


def test_post_natal_only_profiles_need_no_phase_vocabulary(tmp_path: Path) -> None:
    """Section 4.41, for PRF013: a skip only when there was work to skip.

    A null `wpc` is post-natal by construction and matches no cardiac phase
    whether or not one is curated, so an absent vocabulary costs this corpus
    nothing. Written as its own case rather than trusted to the parametrised
    table above, which can only assert that codes *do* fire: a PRF013 keyed on
    "profiles rows exist" rather than "a stage that would be matched exists"
    passes every case in that table and fails only here.
    """
    report = validate_repository(_repo_with_only_post_natal_profiles(tmp_path))

    codes = {issue.code for issue in report.issues}
    assert "PRF013" not in codes
    assert "PRF000" not in codes


@pytest.mark.parametrize(
    ("include_validity_mirrors", "expect_tbl012"),
    [
        pytest.param(True, False, id="mirrors-readable"),
        pytest.param(False, True, id="mirrors-missing"),
    ],
)
def test_prf_references_run_on_both_scope_branches(
    tmp_path: Path, include_validity_mirrors: bool, expect_tbl012: bool
) -> None:
    """`validate_profile_references` has one call site, and both branches of
    the `mirrored` if/else below it must still reach it.

    Section 4.34: a validator called on two branches is tested on the branch
    you were thinking about. `validate_scope_attribution` is called once on
    each side of this exact if/else, and only the mirror-unreadable side had a
    test -- deleting the mirror-readable call passed all 789 tests (see
    `test_a_false_scope_attribution_is_caught_on_the_branch_that_actually_runs`
    above). `validate_profile_references`'s inputs (`datasets`, `known_genes`,
    `phases`) do not depend on `mirrored` at all, so this module places one
    call site *before* that split rather than duplicating it inside both
    branches -- but "one call site, unconditional" is a claim about the
    source, and this is the measurement: a profiles row citing an
    unregistered gene must still report PRF007 whichever way the split
    resolves, proving a future refactor that "helpfully" moves the call inside
    just one branch -- the shape of the #30 bug -- would be caught.
    """
    (tmp_path / "curation").mkdir()
    (tmp_path / "mirrors" / "profiles").mkdir(parents=True)
    (tmp_path / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:99999", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
    )
    # A registry with a real row, NOT a bare header. An empty-but-readable
    # registry reaches `known_genes or None` as None, which skips the gene
    # check entirely -- so a header-only fixture would make PRF007 fire for the
    # wrong reason (the registry "did not load") and would keep passing even if
    # the unregistered-gene rule were deleted. The row below is a gene the
    # profiles mirror does NOT cite, so PRF007 fires on HGNC:99999 specifically.
    (tmp_path / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box 5", "", "ENSG00000089225", "6910", "", "", ""])
        + "\n"
    )
    if include_validity_mirrors:
        (tmp_path / "mirrors" / "clingen_gene_validity.tsv").write_text(_EMPTY_CLINGEN_TSV)
        (tmp_path / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)

    report = validate_repository(tmp_path)

    codes = {issue.code for issue in report.issues}
    assert ("TBL012" in codes) is expect_tbl012, "fixture landed on the wrong scope branch"
    assert "PRF007" in codes, "validate_profile_references did not run on this branch"


def test_validate_profiles_internal_consistency_runs_through_the_repository(
    tmp_path: Path,
) -> None:
    """`validate_profiles`'s own checks (PRF001/002/003/008) are unit-tested
    directly in `test_validate_profiles.py`; this proves the *wiring* reaches
    them too, the same distinction `test_scope_checks_run_against_the_real_
    repository` draws for the scope validators.

    Called unconditionally, right beside `validate_burden`'s own unconditional
    call -- deleting either is a defect of the same shape, and PRF008 is the
    cheapest of the four to provoke: one dataset's `profiles` rows reporting
    two different units needs no `profile_quantiles` shard at all.
    """
    (tmp_path / "mirrors" / "profiles").mkdir(parents=True)
    (tmp_path / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:11604", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
        + "\t".join(["GSE999999", "HGNC:11604", "Liver", "7wpc", "10.0", "rpkm", "", "", "5"])
        + "\n"
    )

    report = validate_repository(tmp_path)

    assert "PRF008" in {issue.code for issue in report.issues}


def test_prf009_fires_against_a_real_published_gene_population(tmp_path: Path) -> None:
    """`_gate_published_genes` computes an actual population, not a permanent no-op.

    PRF009 (a published gene missing a `profiles` row in a cell that was
    assayed) can only ever fire when `published_genes` is non-empty -- an
    implementation that always handed `validate_profile_references` an empty
    set would make this check permanently dead, and every other test in this
    file is deliberately built so PRF009 never fires, which would hide exactly
    that. This builds the smallest repository where `build.validity.
    published_genes` admits a real gene -- one ClinGen `Definitive` record for
    a disease `chd_scope.yaml` places in scope -- and proves that population
    reaches `validate_profile_references` through `validate_repository`, not
    merely through the private helper in isolation.
    """
    root = tmp_path
    (root / "curation").mkdir()
    (root / "curation" / "chd_scope.yaml").write_text(
        "diseases:\n"
        "  - id: MONDO:0007732\n"
        "    label: Holt-Oram syndrome\n"
        "    reason: ClinGen classifies this gene-disease pair as CHD.\n"
        "    admitted_by: clingen_chd_panel\n"
        "    attributed_to: Congenital Heart Disease\n"
        "    admitted_on: 2026-01-01\n"
    )
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box transcription factor 5", "", "", "", "", "", ""])
        + "\n"
    )
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\tDefinitive\t"
        "2023-01-01\tCongenital Heart Disease\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)
    (root / "mirrors" / "profile_quantiles").mkdir(parents=True)
    (root / "mirrors" / "profile_quantiles" / "GSE999999.tsv").write_text(
        "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
        "GSE999999\tHeart\t7wpc\t50\t10.0\ttpm\t20000\n"
    )

    report = validate_repository(root)

    issue = next(i for i in report.issues if i.code == "PRF009")
    assert "HGNC:11604" in issue.message


def test_an_unmapped_classification_does_not_crash_validate_repository(tmp_path: Path) -> None:
    """`_gate_published_genes` must never let `build.validity.gene_validity` raise
    out of `validate_repository`.

    `gene_validity` is written for post-gate data and raises `ValueError` by
    design when a mirror carries a classification term the atlas vocabulary
    has not mapped yet -- reproduced directly against a synthetic frame while
    writing `_gate_published_genes` (see its docstring). That is the right
    behaviour for `build_site`, which must refuse rather than publish a
    confidence nobody vetted, but this call sits inside `validate_repository`
    to compute one WARNING-level check's population, and "validators report,
    they do not raise" (CLAUDE.md section 1) must hold regardless of what a
    third-party mirror contains -- including mid-curation, before
    `validate_table` has even finished reporting the same value as TBL004.
    """
    root = tmp_path
    (root / "curation").mkdir()
    (root / "curation" / "chd_scope.yaml").write_text(
        "diseases:\n"
        "  - id: MONDO:0007732\n"
        "    label: Holt-Oram syndrome\n"
        "    reason: ClinGen classifies this gene-disease pair as CHD.\n"
        "    admitted_by: clingen_chd_panel\n"
        "    attributed_to: Congenital Heart Disease\n"
        "    admitted_on: 2026-01-01\n"
    )
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(_GENES_TSV_HEADER)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\t"
        "Not A Real Classification\t2023-01-01\tCongenital Heart Disease\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)

    report = validate_repository(root)  # must return, never raise

    assert isinstance(report, ValidationReport)
    # The malformed value is still visible in the report, via the schema
    # check that runs independently of `_gate_published_genes` -- proving the
    # defensive fallback silences the crash, not the evidence of the defect.
    assert any(issue.code == "TBL004" for issue in report.issues)


def test_a_programmer_error_in_the_gate_propagates_rather_than_silencing_prf009(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sibling of the test above, and the case it does NOT cover.

    `_gate_published_genes` absorbs `KeyError` and `ValueError` because those
    are the *data* failures its docstring names, and each already arrives with
    its own ERROR from `validate_table` -- so degrading to an empty published
    set cannot hide anything.

    A blanket `except Exception` also absorbed programmer errors, and those
    have no accompanying validation error to make them visible. Measured
    2026-08-19: injecting a `TypeError` into `gene_validity` -- the shape of an
    ordinary refactor regression -- produced `0 error(s), 4 warning(s)` with
    `ok=True` while PRF009 checked nothing at all. The report was green and one
    validator was dead.

    So this asserts the opposite of its sibling on purpose: a `ValueError` is
    absorbed and a `TypeError` is not. Raised by review on PR #39; the
    docstring's own guarantee had covered data failures and not code ones.
    """

    def boom(*args: object, **kwargs: object) -> dict[str, object]:
        raise TypeError("simulated refactor regression, not a data problem")

    monkeypatch.setattr("chd_atlas.validate.runner.gene_validity", boom)
    with pytest.raises(TypeError, match="simulated refactor regression"):
        validate_repository(Path("."))


@pytest.mark.parametrize(
    ("column", "value", "expected"),
    [
        # No burden row has comparator `none` for this study: BUR020 only.
        pytest.param("combines", "mutation_model;none", ["BUR020"], id="arm-missing"),
        # An unknown study: BUR022 from the references call, and one BUR020 per
        # arm, because neither arm carries that study.
        pytest.param("study", "PMID:99999999", ["BUR020", "BUR020", "BUR022"], id="unknown-study"),
    ],
)
def test_a_combined_row_is_checked_by_the_real_validator(
    tmp_path: Path, column: str, value: str, expected: list[str]
) -> None:
    """Both combined validators are reached from `validate_repository` (§4.34).

    The *last* row is edited so the file stays sorted: the key is gene, then
    study, and the last row is its gene's last.
    """
    root = _real_repo(tmp_path)
    path = root / "mirrors" / "burden_combined.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    last = lines[-1].split("\t")
    last[header.index(column)] = value
    path.write_text("\n".join([*lines[:-1], "\t".join(last)]) + "\n", encoding="utf-8")

    report = validate_repository(root)
    codes = sorted(issue.code for issue in report.issues if issue.severity is Severity.ERROR)
    assert codes == expected


@pytest.mark.parametrize("validity_mirrors", ["readable", "missing"])
def test_an_orphanet_warrant_is_checked_on_both_validity_branches(
    tmp_path: Path, validity_mirrors: str
) -> None:
    """SCP006 is one call site, reached whether or not the validity mirrors read (§4.34).

    KBG syndrome, which nobody has nominated, put forward on Noonan syndrome's
    Orphanet code: MONDO:0007846 does not reach ORPHA:648, so the warrant fails
    on either branch. (This used Noonan itself until Noonan was nominated for
    real, when the appended entry became a duplicate id and the corpus stopped
    loading: the fixture must pick a term outside the committed scope.)
    """
    root = _real_repo(tmp_path)
    scope = root / "curation" / "chd_scope.yaml"
    scope.write_text(
        scope.read_text(encoding="utf-8") + "  - id: MONDO:0007846\n"
        "    label: KBG syndrome\n"
        "    reason: A nomination on the wrong Orphanet code.\n"
        "    admitted_by: orphanet_phenotype_frequency\n"
        "    attributed_to: Orphanet\n"
        "    orphanet: ORPHA:648\n"
        "    admitted_on: 2026-10-25\n",
        encoding="utf-8",
    )
    if validity_mirrors == "missing":
        (root / "mirrors" / "clingen_gene_validity.tsv").unlink()
        (root / "mirrors" / "gencc_submissions.tsv").unlink()

    report = validate_repository(root)
    refusals = [i for i in report.issues if i.code == "SCP006"]
    assert len(refusals) == 1
    assert "does not reach ORPHA:648" in refusals[0].message


def test_a_malformed_identity_id_on_a_published_gene_fails_validate_as_the_build_refuses_it(
    tmp_path: Path,
) -> None:
    """GEN004 (owner decision 2026-10-26, R13), from the real validator on a real corpus.

    Measured that day on a copy of the committed repository with TBX5's
    Ensembl id cut to `ENSG0089225`: `chd-atlas validate` reported 0 errors and
    `chd-atlas build` refused the repository after writing had begun,
    `external_links` raising on the page. Validate now reports exactly that
    one error, so its verdict is the build's, and the build refuses at the
    gate, before anything is written -- the refusal still holds, earlier.
    """
    repo = _real_repo(tmp_path / "repo")
    registry = repo / "mirrors" / "genes.tsv"
    text = registry.read_text(encoding="utf-8")
    assert text.count("\tENSG00000089225\t") == 1, "TBX5's registry row moved"
    registry.write_text(text.replace("\tENSG00000089225\t", "\tENSG0089225\t"), encoding="utf-8")

    report = validate_repository(repo)

    errors = [issue for issue in report.issues if issue.severity is Severity.ERROR]
    assert [issue.code for issue in errors] == ["GEN004"]
    assert "HGNC:11604" in errors[0].message and "'ENSG0089225'" in errors[0].message
    assert report.ok is False
    out = tmp_path / "dist"
    with pytest.raises(BuildRefused, match="GEN004"):
        build_site(repo, out)
    assert not out.exists()


def test_the_table_a_chd_only_count_comes_from_is_held_to_the_source_registry(
    tmp_path: Path,
) -> None:
    """`n_chd_case_carriers` is transcribed from a file no row's `source` names.

    S8 Table, `CHD_CASE_SPLITS`'s source for PMID:34324492, whose rows name S6
    as their `source`. Nothing read its registry entry until the review
    (2026-10-27) renamed it, and marked it `redistribution: prohibited`, and
    the whole suite passed both times. `_used_sources` now counts it as a
    source the burden mirror relies on wherever a row carries the column, so
    it is held to SRC003 and SRC004 like every row's own: renamed here, it is
    SRC003, the one error. Both edits to the committed registry fail
    `test_the_committed_repository_validates_cleanly`; this is the test that
    fails if `_used_sources` stops collecting it.
    """
    root = _real_repo(tmp_path)
    path = root / "mirrors" / "sources.yaml"
    text = path.read_text(encoding="utf-8")
    assert text.count("id: audain2021_s8\n") == 1
    path.write_text(text.replace("id: audain2021_s8\n", "id: audain2021_s8_renamed\n"))

    report = validate_repository(root)
    errors = [issue for issue in report.issues if issue.severity is Severity.ERROR]
    assert [(issue.code, "'audain2021_s8'" in issue.message) for issue in errors] == [
        ("SRC003", True)
    ]


def test_a_cell_type_mirror_with_an_unregistered_source_fails_validation(
    tmp_path: Path,
) -> None:
    """`cell_type_profiles` carries `source`, so `_used_sources` holds it to SRC003.

    Other checks (no `hoc` dataset record, an unregistered gene) may flag the
    row as well; this test is about the source gate, so it asserts SRC003 only.
    """
    root = _real_repo(tmp_path)
    shard = root / "mirrors" / "cell_type_profiles"
    shard.mkdir(parents=True, exist_ok=True)
    (shard / "hoc.tsv").write_text(
        "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
        "n_donors\tn_cells\tfraction_expressing\tsource\n"
        "hoc\tHGNC:1100\tcardiomyocyte\tLV\t2.5\tcpm\t1.0\t3.0\t4\t900\t0.4\tnot_in_sources\n",
        encoding="utf-8",
    )

    report = validate_repository(root)

    assert "SRC003" in [issue.code for issue in report.issues]
    assert any("'not_in_sources'" in i.message for i in report.issues if i.code == "SRC003")


def test_a_header_only_gene_registry_skips_only_the_cell_type_gene_check(tmp_path: Path) -> None:
    """CTP007 is skipped for an empty registry; REF001 names the cause instead.

    The other cell-type checks still run: a second shard's row names a dataset
    with no record, and CTP001 says so. That assertion is the one that fails if
    `validate_repository` stops calling `validate_cell_type_profiles` at all
    (review, 2026-10-28: with the call deleted, the suite passed and a build
    skipped every row naming an undeclared cell type or region in silence).

    **The committed `hoc.tsv` stays beside it since 2026-11-02**, when Farah et
    al. 2024's 12,308 rows over 181 genes landed. The test wrote its row over
    that file until then, and against the real record the row's undeclared cell
    type raised CTP002 and the emptied dataset CTP008 in CTP001's place. Kept,
    the real rows make the skip a measurement at scale: every one of those 181
    genes would be a CTP007 against the empty registry, and none may be. They
    are clean, so CTP001 must be the only cell-type code, which also fails if
    the skip turns into a flood of anything else. Measured 2026-11-02: the
    old fixture on the real corpus raised CTP002 once and CTP008 17 times; the
    committed rows against an empty registry raise CTP007 181 times; deleting
    the call, or passing the registry bare rather than `or None`, fails here.
    """
    root = _real_repo(tmp_path)
    genes = root / "mirrors" / "genes.tsv"
    genes.write_text(genes.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")
    shard = root / "mirrors" / "cell_type_profiles"
    assert (shard / "hoc.tsv").is_file(), "the committed cell-type mirror is gone; re-read this"
    (shard / "unrecorded.tsv").write_text(
        "dataset\tgene\tcell_type\tregion\tmedian_abundance\tunit\tq25\tq75\t"
        "n_donors\tn_cells\tfraction_expressing\tsource\n"
        "unrecorded\tHGNC:1100\tcardiomyocyte\tall\t2.5\tcpm\t1.0\t3.0\t4\t900\t0.4\t"
        "farah2024_hoc\n",
        encoding="utf-8",
    )

    codes = [issue.code for issue in validate_repository(root).issues]

    assert "REF001" in codes
    assert sorted({code for code in codes if code.startswith("CTP")}) == ["CTP001"]


# ---- literature records (LIT001-LIT014) reached through validate_repository ----


def _literature_repo(tmp_path: Path, *, stem: str = "TBX5", gene: str = "HGNC:11604", **tweak):  # type: ignore[no-untyped-def]
    """A copy of the committed repository plus one clean literature file.

    PMID:8988165 is the real TBX5 paper already in `publications.yaml`. Counts
    differ from one another: one search, one new hit, one finding.
    """
    root = _real_repo(tmp_path)
    (root / "curation" / ".id_registry.yaml").write_text("prefixes:\n  AST: 2\n  LIT: 2\n")
    finding = {
        "id": "CHDA:LIT:0000001",
        "publication": "PMID:8988165",
        "design": "cohort_screen",
        "methods": ["exome"],
        "count_unit": "probands",
        "screened": 7,
        "carriers": {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
        "variant_classes": ["snv_indel"],
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "syndromic": "both",
        "shares_probands_with": "unknown",
        "locators": [{"kind": "in_text", "element": "Table 2"}],
        "summary": "Six of seven probands carry a variant.",
        "extracted_by": "agent:extract/test",
        "approved_by": "A. Curator",
        "approved_on": "2026-11-05",
    }
    finding.update(tweak.pop("finding", {}))
    finding = stamped(finding)
    raw = {
        "gene": gene,
        "searches": [
            {
                "id": 1,
                "source": "pubmed",
                "query": "q",
                "run_on": "2026-11-02",
                "hits": 1,
                "new_hits": 1,
            }
        ],
        "triage": [{"pmid": finding["publication"], "search": 1, "decision": "included"}],
        "findings": [finding],
    }
    raw.update(tweak)
    directory = root / "curation" / "literature"
    directory.mkdir()
    (directory / f"{stem}.yaml").write_text(json.dumps(raw))  # JSON is YAML
    return root


def _codes_of(root: Path) -> set[str]:
    return {issue.code for issue in validate_repository(root).issues}


def test_a_clean_literature_file_adds_no_issue(tmp_path: Path) -> None:
    baseline = validate_repository(_real_repo(tmp_path / "a")).issues
    report = validate_repository(_literature_repo(tmp_path / "b"))
    assert report.issues == baseline


@pytest.mark.parametrize(
    ("kwargs", "code"),
    [
        (
            {
                "searches": [
                    {
                        "id": 1,
                        "source": "pubmed",
                        "query": "q",
                        "run_on": "2026-11-02",
                        "hits": 2,
                        "new_hits": 2,
                    }
                ]
            },
            "LIT001",
        ),
        ({"stem": "HOS"}, "LIT009"),
        ({"gene": "HGNC:10477", "stem": "RXRA"}, "LIT013"),
        ({"finding": {"publication": "PMID:1004"}}, "LIT011"),
        ({"finding": {"cohorts": ["no_such_cohort"]}}, "LIT012"),
        ({"finding": {"phenotypes": ["HP:9999999"]}}, "ONT001"),
        ({"finding": {"id": "CHDA:LIT:0000009"}}, "ID002"),
        # PMID:8988165 has no burden rows (final review I6); septal-only HP:0001631 (I3).
        ({"finding": {"burden_study": "PMID:8988165"}}, "LIT014"),
        ({"finding": {"lesion_groups": ["septal", "conotruncal"]}}, "REF009"),
    ],
)
def test_each_literature_check_is_reached_through_the_repository_validator(
    tmp_path: Path, kwargs: dict[str, object], code: str
) -> None:
    root = _literature_repo(tmp_path, **kwargs)  # type: ignore[arg-type]
    assert code in _codes_of(root)


@pytest.mark.parametrize(
    ("removed", "cause", "silenced"),
    [("publications.yaml", "CUR003", "LIT011"), ("cohorts.yaml", "CUR004", "LIT012")],
)
def test_a_missing_registry_a_finding_needs_is_reported_once_not_skipped(
    tmp_path: Path, removed: str, cause: str, silenced: str
) -> None:
    """R4: the registry's absence skips LIT011/LIT012, so the cause is an error.

    Everything else that would raise the same code is removed, so only the
    literature file can be the reason it appears.
    """
    root = _literature_repo(
        tmp_path, finding={"cohorts": ["ddd"]} if removed == "cohorts.yaml" else {}
    )
    (root / "mirrors" / "burden.tsv").unlink()
    (root / "mirrors" / "burden_combined.tsv").unlink(missing_ok=True)
    shutil.rmtree(root / "curation" / "assertions")
    (root / "curation" / "assertions").mkdir()
    (root / "curation" / removed).unlink()

    report = validate_repository(root)
    codes = {issue.code for issue in report.issues if issue.severity is Severity.ERROR}

    assert cause in codes
    assert silenced not in codes
    assert report.ok is False


def test_a_degraded_publication_gate_does_not_report_every_literature_file(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`_gate_published_genes` degrades to an empty set; LIT013 must skip, not
    fire for every file. The empty set is passed as None (`published or None`)."""

    def boom(*args: object, **kwargs: object) -> None:
        raise ValueError("unmapped classification")

    monkeypatch.setattr("chd_atlas.validate.runner.gene_validity", boom)

    codes = _codes_of(_literature_repo(tmp_path))

    assert "LIT013" not in codes


def test_a_header_only_gene_registry_skips_lit009_and_still_names_its_cause(
    tmp_path: Path,
) -> None:
    """An empty registry is not None; without `or None` every literature file
    reported LIT009 and buried the registry. The cause stays an error."""
    root = _literature_repo(tmp_path)
    header = (root / "mirrors" / "genes.tsv").read_text().splitlines()[0]
    (root / "mirrors" / "genes.tsv").write_text(header + "\n")

    report = validate_repository(root)
    errors = {i.code for i in report.issues if i.severity is Severity.ERROR}

    assert "LIT009" not in errors
    assert errors & {"REF001", "GEN004"}
    assert report.ok is False
