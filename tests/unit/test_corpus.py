# tests/unit/test_corpus.py
from pathlib import Path

import pytest

from chd_atlas.corpus import load_curation, unexpected_curation_entries

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

PUBLICATIONS_YAML = """\
publications:
  - id: PMID:8988165
    title: Mutations in human TBX5 cause limb and cardiac malformation
    journal: Nature Genetics
    year: 1997
    authors: [Basson CT]
    study_type: case_control
"""


def _write_minimal_corpus(root: Path) -> None:
    (root / "curation" / "assertions").mkdir(parents=True)
    (root / "curation" / "assertions" / "TBX5.yaml").write_text(ASSERTION_YAML)
    (root / "curation" / "publications.yaml").write_text(PUBLICATIONS_YAML)


def test_loads_assertions_and_publications(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert len(corpus.assertions) == 1
    assert corpus.assertions[0].gene == "HGNC:11604"
    assert len(corpus.publications) == 1


def test_schema_violation_becomes_an_issue_rather_than_an_exception(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    broken = ASSERTION_YAML.replace("mechanism: haploinsufficiency", "mechanism: very_sure")
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(broken)

    corpus, issues = load_curation(tmp_path)

    assert corpus.assertions == ()
    assert len(issues) == 1
    assert issues[0].code == "SCHEMA001"
    assert "mechanism" in issues[0].location


def test_reports_every_schema_violation_not_just_the_first(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    broken = ASSERTION_YAML.replace("mechanism: haploinsufficiency", "mechanism: magic")
    broken = broken.replace("syndromic: both", "syndromic: sometimes")
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(broken)

    _, issues = load_curation(tmp_path)

    assert len(issues) == 2


def test_malformed_yaml_becomes_a_parse_issue(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text("assertions: [oops\n")

    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["YAML001"]


def test_out_of_range_date_is_reported_rather_than_raised(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(
        ASSERTION_YAML.replace("2026-07-01", "2026-13-45")
    )

    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["YAML001"]


def test_missing_curation_directory_is_reported(tmp_path: Path) -> None:
    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["CORPUS001"]


def test_non_utf8_bytes_are_reported_rather_than_raised(tmp_path: Path) -> None:
    """A curator name saved as Latin-1 must not crash the whole validation run."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_bytes(
        ASSERTION_YAML.replace("0000-0002-1825-0097", "François").encode("latin-1")
    )

    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["YAML001"]


def test_empty_yaml_file_is_reported(tmp_path: Path) -> None:
    """yaml.load("") returns None; that must not read as 'no records, no problem'."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "publications.yaml").write_text("")

    corpus, issues = load_curation(tmp_path)

    assert corpus.publications == ()
    assert [issue.code for issue in issues] == ["SCHEMA001"]


def test_assertions_from_multiple_files_are_concatenated(tmp_path: Path) -> None:
    """A regression to assignment instead of extend would pass every other test."""
    _write_minimal_corpus(tmp_path)
    second = ASSERTION_YAML.replace("CHDA:AST:0000001", "CHDA:AST:0000002").replace(
        "HGNC:11604", "HGNC:4173"
    )
    (tmp_path / "curation" / "assertions" / "GATA4.yaml").write_text(second)

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert {a.id for a in corpus.assertions} == {"CHDA:AST:0000001", "CHDA:AST:0000002"}


def test_loads_phenotypes_and_featured(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "phenotypes.yaml").write_text(
        "phenotypes:\n"
        "  - id: HP:0001631\n"
        "    label: Atrial septal defect\n"
        "    lesion_group: septal\n"
    )
    (tmp_path / "curation" / "featured.yaml").write_text(
        "featured:\n"
        "  - publication: PMID:8988165\n"
        "    order: 1\n"
        "    blurb: Founding TBX5 paper.\n"
        "    topic: genomics\n"
    )

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert len(corpus.phenotypes) == 1
    assert len(corpus.featured) == 1


def test_an_absent_cardiac_phases_file_is_legitimate_and_loads_as_none(tmp_path: Path) -> None:
    """No curator has verified a phase source yet, so the file may not exist at all.

    Distinct from the model-level tests in test_models_phases.py: this exercises
    `load_curation` itself, which is the one place a wiring mistake -- forgetting
    to pass the parsed value into the `Corpus(...)` call -- would otherwise reach
    no test. Nothing under `build/` reads this field yet, so a silently dropped
    value would be invisible everywhere else in the suite.
    """
    _write_minimal_corpus(tmp_path)

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert corpus.cardiac_phases is None


def test_loads_cardiac_phases(tmp_path: Path) -> None:
    """A real (non-placeholder) file must reach `corpus.cardiac_phases` intact."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "cardiac_phases.yaml").write_text(
        "attributed_to: O'Rahilly & Muller 1987\n"
        "citation: ISBN:0872796248\n"
        "phases:\n"
        "  - id: looping\n"
        "    go_id: GO:0001947\n"
        "    label: Cardiac looping\n"
        "    start_wpc: 3.0\n"
        "    end_wpc: 5.0\n"
        "    start_carnegie_stage: CS10\n"
        "    end_carnegie_stage: CS12\n"
        "    start_hsapdv_id: HsapDv:0000017\n"
        "    end_hsapdv_id: HsapDv:0000019\n"
        "    end_basis: stated\n"
    )

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert corpus.cardiac_phases is not None
    assert [phase.id for phase in corpus.cardiac_phases.phases] == ["looping"]


def test_a_misnamed_assertion_file_is_reported(tmp_path: Path) -> None:
    """Renaming TBX5.yaml to TBX5.yml made the assertion invisible to every
    check while the gate still exited 0."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").rename(
        tmp_path / "curation" / "assertions" / "TBX5.yml"
    )

    issues = unexpected_curation_entries(tmp_path)

    assert [i.code for i in issues] == ["CUR001"]
    assert "TBX5.yml" in issues[0].location


def test_a_misspelled_record_directory_is_reported(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions").rename(tmp_path / "curation" / "assertion")

    issues = unexpected_curation_entries(tmp_path)

    assert [i.code for i in issues] == ["CUR001"]


def test_a_directory_where_a_record_belongs_is_reported(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").unlink()
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").mkdir()

    issues = unexpected_curation_entries(tmp_path)

    assert [i.code for i in issues] == ["CUR001"]


def test_an_absent_record_directory_is_legitimate(tmp_path: Path) -> None:
    """There is no functional/ before any functional evidence is curated."""
    _write_minimal_corpus(tmp_path)

    assert unexpected_curation_entries(tmp_path) == []


def test_an_unreadable_file_is_reported_rather_than_raised(tmp_path: Path) -> None:
    """OSError escaped the guard: an unreadable mode, a dangling symlink, or a
    directory where a file was expected each aborted the whole run."""
    _write_minimal_corpus(tmp_path)
    target = tmp_path / "curation" / "assertions" / "TBX5.yaml"
    target.chmod(0o000)
    try:
        _, issues = load_curation(tmp_path)
        assert [issue.code for issue in issues] == ["YAML001"]
    finally:
        target.chmod(0o644)


LITERATURE_YAML = """\
gene: HGNC:11604
searches:
  - {id: 1, source: pubmed, query: "TBX5[tiab]", run_on: 2026-11-02, hits: 2, new_hits: 2}
triage:
  - {pmid: "PMID:8988165", search: 1, decision: pending}
findings: []
"""


def test_literature_files_load_in_filename_order_with_their_paths(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    directory = tmp_path / "curation" / "literature"
    directory.mkdir()
    (directory / "TBX5.yaml").write_text(LITERATURE_YAML)
    (directory / "GATA4.yaml").write_text(LITERATURE_YAML.replace("11604", "4172"))

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert [(p.name, f.gene) for p, f in corpus.literature] == [
        ("GATA4.yaml", "HGNC:4172"),
        ("TBX5.yaml", "HGNC:11604"),
    ]
    assert unexpected_curation_entries(tmp_path) == []


def test_a_corpus_without_literature_has_none(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    corpus, _ = load_curation(tmp_path)
    assert corpus.literature == ()


def test_a_bad_literature_file_is_an_issue_and_a_stray_file_is_reported(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    directory = tmp_path / "curation" / "literature"
    directory.mkdir()
    (directory / "TBX5.yaml").write_text(LITERATURE_YAML.replace("hits: 2, new", "hits: 1, new"))
    (directory / "x.txt").write_text("stray")

    corpus, issues = load_curation(tmp_path)

    assert corpus.literature == ()
    assert [i.code for i in issues] == ["SCHEMA001"]
    stray = unexpected_curation_entries(tmp_path)
    assert [i.code for i in stray] == ["CUR001"]
    assert "x.txt" in stray[0].location


_SHORTHAND_BODIES = [
    ("anchor &s", "x: &s 1\n"),
    ("anchor &s", "x: &s 1\ny: *s\n"),  # an alias needs an anchor, which is reported first
    ("merge key <<", "x: {a: 1}\ny:\n  <<: {a: 2}\n"),
]


@pytest.mark.parametrize(("named", "extra"), _SHORTHAND_BODIES)
def test_a_literature_file_with_yaml_shorthand_is_refused_and_publishes_nothing(
    tmp_path: Path, named: str, extra: str
) -> None:
    """A merge key copies one finding's approval onto another, so the file is refused (R34).

    Measured on the pre-fix loader: `- <<: *a` with a new id loaded as a valid,
    approved finding. Anchors and aliases alone are refused too, because the
    loader cannot tell which of them carries an approval. The extra keys sit
    outside the model, so `extra="forbid"` would also refuse them; the test
    therefore asserts the code is YAML002 and that nothing else is reported.
    """
    _write_minimal_corpus(tmp_path)
    directory = tmp_path / "curation" / "literature"
    directory.mkdir()
    (directory / "TBX5.yaml").write_text(LITERATURE_YAML + extra)

    corpus, issues = load_curation(tmp_path)

    assert corpus.literature == ()
    assert [i.code for i in issues] == ["YAML002"]
    assert issues[0].message.startswith(named)
    assert issues[0].location.endswith("TBX5.yaml")
