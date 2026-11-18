# tests/unit/test_cli_literature.py
"""`chd-atlas literature review | target | frame`.

The repository under test is a tmp directory: a four-gene registry, two
published genes (TBX5, GATA4 -- different diseases, so a label cannot pass by
coming from the wrong gene) and the real `curation/chd_scope.yaml`. Nothing here
reaches the network. `search`, `abstracts`, `add-publication` and `diff` moved
to chd-atlas/curator, which tests them.
"""

from __future__ import annotations

import json
import re
from collections.abc import Callable
from datetime import date
from pathlib import Path
from typing import Any

import pytest
from ruamel.yaml import YAML
from typer.testing import CliRunner

from chd_atlas.cli import app
from chd_atlas.models.gene_literature import Finding, finding_digest
from tests.unit._stamp import stamped_yaml
from tests.unit.mane_fixtures import synthetic_rows, write_mane_mirrors

REPO = Path(__file__).resolve().parents[2]
runner = CliRunner()

GENES = (
    "hgnc_id\tsymbol\tname\taliases\tprev_symbols\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
    "HGNC:11604\tTBX5\tT-box 5\tTBX5A\tHOS\t\t\t\t\t\n"
    "HGNC:4173\tGATA4\tGATA binding protein 4\t\t\t\t\t\t\t\n"
    "HGNC:1100\tBRCA1\tBRCA1 DNA repair\t\t\t\t\t\t\t\n"
    "HGNC:6772\tSMAD6\tSMAD family member 6\t\tMADH7\t\t\t\t\t\n"
    "HGNC:6773\tSMAD7\tSMAD family member 7\t\tMADH7\t\t\t\t\t\n"
)
CLINGEN_HEADER = (
    "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
    "classification_date\tgcep\treport_url\n"
)
CLINGEN = (
    "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP11\tDefinitive\t"
    "2025-03-25T16:00:00.000Z\tSyndromic Disorders Gene Curation Expert Panel\thttps://x/1\n"
    "HGNC:4173\tGATA4\tMONDO:0100009\tstructural congenital heart disease, multiple types - GATA4\t"
    "AD\tSOP9\tDefinitive\t2023-02-07T17:00:00.000Z\t"
    "Congenital Heart Disease Gene Curation Expert Panel\thttps://x/2\n"
)
GENCC_HEADER = (
    "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\tsubmitter\t"
    "submitted_on\treport_url\n"
)

PUBLICATIONS = """\
# Head comment: every record is transcribed from PubMed.
publications:
  - id: PMID:111
    title: First paper.
    journal: Journal One
    year: 1999
    authors:
      - Aa A
    # A comment sitting above a neighbouring field.
    study_type: case_report
  # A comment sitting between two records.
  - id: PMID:222
    title: Second paper.
    journal: Journal Two
    year: 2001
    authors:
      - Bb B
    study_type: gwas
"""


@pytest.fixture
def root(tmp_path: Path) -> Path:
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "curation").mkdir()
    (tmp_path / "mirrors" / "genes.tsv").write_text(GENES)
    (tmp_path / "mirrors" / "clingen_gene_validity.tsv").write_text(CLINGEN_HEADER + CLINGEN)
    (tmp_path / "mirrors" / "gencc_submissions.tsv").write_text(GENCC_HEADER)
    (tmp_path / "curation" / "chd_scope.yaml").symlink_to(REPO / "curation" / "chd_scope.yaml")
    (tmp_path / "curation" / "publications.yaml").write_text(PUBLICATIONS)
    return tmp_path


# --------------------------------------------------------------------------
# review (Task 4)
# --------------------------------------------------------------------------

LITERATURE = """\
# Curator's head comment -- must survive every write.
gene: HGNC:11604
searches:
  - id: 1
    source: pubmed
    query: TBX5[tiab]
    run_on: 2026-11-02
    hits: 3
    new_hits: 3
triage:
  - pmid: PMID:8988165
    search: 1
    decision: included
findings:
  - id: CHDA:LIT:0000001
    publication: PMID:8988165
    design: cohort_screen
    methods:
      - exome
    count_unit: probands
    screened: 7
    carriers:
      de_novo: 2
      inherited: 1
      origin_unknown: 3
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: both
    cohorts:
      - alpha_cohort
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Table 2 & <3>
    summary: Six of seven probands carry a variant.
    extracted_by: agent:extract/test
  # A comment between two findings.
  - id: CHDA:LIT:0000002
    publication: PMID:111
    design: common_variant_association
    methods:
      - genotyping_array
    cases: 11
    controls: 13
    p_value: 0.02
    odds_ratio: 1.5
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: supplementary
        file: s1.xlsx
        sheet: Sheet 1
        row: 4
    summary: A GWAS hit.
    extracted_by: agent:extract/test
    verified_by: agent:verify/earlier
  - id: CHDA:LIT:0000003
    publication: PMID:222
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - conotruncal
    syndromic: syndromic
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A proband.
    extracted_by: agent:extract/test
    approved_by: Dr. Ada Curator
    approved_on: 2026-10-01
  - id: CHDA:LIT:0000004
    publication: PMID:333
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 0
      inherited: 1
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Another proband.
    extracted_by: agent:extract/test
"""

VERIFIER = """\
findings:
  - id: CHDA:LIT:0000001
    publication: PMID:8988165
    design: cohort_screen
    methods:
      - exome
    count_unit: probands
    screened: 7
    carriers:
      de_novo: 2
      inherited: 1
      origin_unknown: 3
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: both
    cohorts:
      - alpha_cohort
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Table 2, worded differently
    summary: Different words, same finding.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000002
    publication: PMID:111
    design: common_variant_association
    methods:
      - genotyping_array
    cases: 11
    controls: 13
    p_value: 0.02
    odds_ratio: 1.5
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A GWAS hit.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000003
    publication: PMID:222
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - conotruncal
    syndromic: syndromic
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: A proband.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000004
    publication: PMID:333
    design: case_report
    methods:
      - targeted_sequencing
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Another proband, counted differently.
    extracted_by: agent:verify/blind
  - id: CHDA:LIT:0000005
    publication: PMID:999
    design: case_report
    methods:
      - exome
    count_unit: probands
    screened: 1
    carriers:
      de_novo: 1
      inherited: 0
      origin_unknown: 0
    variant_classes:
      - snv_indel
    phenotypes:
      - HP:0001631
    lesion_groups:
      - septal
    syndromic: isolated
    shares_probands_with: unknown
    locators:
      - kind: in_text
        element: Abstract
    summary: Only the blind pass found this.
    extracted_by: agent:verify/blind
"""

FINDING_ID = re.compile(r"(?m)^  - id: (CHDA:LIT:\d{7})\n")


@pytest.fixture
def lit_path(root: Path) -> Path:
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir()
    path.write_text(stamped_yaml(LITERATURE))
    return path


def _blocks(text: str) -> dict[str, str]:
    """Each finding's exact text, by id; "head" is everything before the first."""
    marks = list(FINDING_ID.finditer(text))
    blocks = {"head": text[: marks[0].start()]}
    for mark, nxt in zip(marks, [*marks[1:], None], strict=True):
        blocks[mark.group(1)] = text[mark.start() : nxt.start() if nxt else len(text)]
    return blocks


def _second(root: Path, text: str = VERIFIER) -> Path:
    path = root / "local" / "second.yaml"
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return path


def _digest_in(text: str, finding_id: str) -> str:
    """The digest `finding_digest` gives the named finding as `text` has it."""
    (raw,) = (f for f in YAML(typ="safe").load(text)["findings"] if f["id"] == finding_id)
    return finding_digest(
        Finding.model_validate(
            {**raw, "approved_by": "p", "approved_on": "2026-11-03", "approved_digest": "x"}
        )
    )


def _review(root: Path, *args: str) -> Any:
    return runner.invoke(app, ["literature", "review", "TBX5", *args, "--root", str(root)])


def test_the_review_listing_pins_one_finding(root: Path, lit_path: Path) -> None:
    result = _review(root)

    assert result.exit_code == 0, result.stdout
    blocks = result.stdout.split("\n\n")
    one = next(b for b in blocks if b.startswith("CHDA:LIT:0000001"))
    assert one == (
        "CHDA:LIT:0000001  https://pubmed.ncbi.nlm.nih.gov/8988165/  cohort_screen  UNVERIFIED\n"
        "  paper: not in curation/publications.yaml\n"
        "  methods: exome\n"
        "  comparator: none\n"
        "  burden study: none\n"
        "  variant classes: snv_indel\n"
        "  screened 7 probands; carriers 6 (de novo 2, inherited 1, origin unknown 3)\n"
        "  phenotypes: HP:0001631\n"
        "  lesion groups: septal\n"
        "  syndromic: both\n"
        "  cohorts: alpha_cohort\n"
        "  locators: Table 2 & <3>\n"
        "  shares probands with: unknown\n"
        "  extracted by: agent:extract/test\n"
        "  Six of seven probands carry a variant."
    )
    two = next(b for b in blocks if b.startswith("CHDA:LIT:0000002"))
    assert two == (
        "CHDA:LIT:0000002  https://pubmed.ncbi.nlm.nih.gov/111/  common_variant_association"
        "  agent:verify/earlier\n"
        "  paper: study_type case_report; own_lab no\n"
        "  methods: genotyping_array\n"
        "  comparator: none\n"
        "  burden study: none\n"
        "  variant classes: none\n"
        "  cases 11, controls 13, p 0.02, OR 1.5\n"
        "  phenotypes: HP:0001631\n"
        "  lesion groups: septal\n"
        "  syndromic: isolated\n"
        "  cohorts: none\n"
        "  locators: s1.xlsx / Sheet 1 / row 4\n"
        "  shares probands with: unknown\n"
        "  extracted by: agent:extract/test\n"
        "  A GWAS hit."
    )
    # CHDA:LIT:0000003 is approved: not listed. Order is by id.
    assert [b[:16] for b in blocks if b.strip()] == [
        "CHDA:LIT:0000001",
        "CHDA:LIT:0000002",
        "CHDA:LIT:0000004",
    ]


def test_the_review_listing_shows_floats_unrounded_and_the_papers_own_lab(
    root: Path, lit_path: Path
) -> None:
    """Final review I2: `:g` showed 1.23456789e-9 as 1.23457e-09; own_lab is from the registry."""
    lit_path.write_text(
        lit_path.read_text()
        .replace("p_value: 0.02", "p_value: 1.23456789e-9")
        .replace("odds_ratio: 1.5", "odds_ratio: 2.718281828")
    )
    (root / "curation" / "publications.yaml").write_text(
        PUBLICATIONS.replace("study_type: case_report", "study_type: gwas\n    own_lab: true")
    )
    two = next(b for b in _review(root).stdout.split("\n\n") if b.startswith("CHDA:LIT:0000002"))
    assert "\n  cases 11, controls 13, p 1.23456789e-09, OR 2.718281828\n" in two
    assert "\n  paper: study_type gwas; own_lab yes\n" in two


def test_the_review_listing_prints_the_verifiers_differences_under_each_finding(
    root: Path, lit_path: Path
) -> None:
    """R35: `--against` shows what the blind pass disagreed on, finding by finding."""
    # Only the first block, with `screened` changed: finding 1 is paired and differs;
    # findings 2 and 4 have no partner in this pass.
    first_only = VERIFIER.split("  - id: CHDA:LIT:0000002")[0].replace(
        "    screened: 7\n", "    screened: 8\n", 1
    )
    second = _second(root, first_only)
    result = _review(root, "--against", str(second))
    assert result.exit_code == 0, result.stdout
    blocks = result.stdout.split("\n\n")
    one = next(b for b in blocks if b.startswith("CHDA:LIT:0000001"))
    assert one.endswith("\n  differs: screened: first 7 | second 8")
    assert "no partner" not in one
    for other in ("CHDA:LIT:0000002", "CHDA:LIT:0000004"):
        block = next(b for b in blocks if b.startswith(other))
        assert block.endswith("\n  second pass: no partner finding")
        assert "differs:" not in block
    # A paired finding that agrees prints neither line.
    agreeing = _second(root, VERIFIER.split("  - id: CHDA:LIT:0000003")[0])
    agreed = _review(root, "--against", str(agreeing)).stdout.split("\n\n")
    two = next(b for b in agreed if b.startswith("CHDA:LIT:0000002"))
    assert "differs:" not in two and "no partner" not in two
    # Without the option nothing about the second pass is printed.
    assert "differs:" not in _review(root).stdout
    # Refused alongside an action; a missing file is a bad argument.
    assert (
        _review(root, "--against", str(second), "--pending-review", "CHDA:LIT:0000001").exit_code
        == 2
    )
    assert _review(root, "--against", str(root / "nope.yaml")).exit_code == 2


def test_the_review_listing_names_the_findings_that_share_probands(
    root: Path, lit_path: Path
) -> None:
    text = lit_path.read_text()
    assert text.count("    shares_probands_with: unknown\n") >= 1
    lit_path.write_text(
        text.replace(
            "    shares_probands_with: unknown\n",
            "    shares_probands_with:\n      - CHDA:LIT:0000004\n      - CHDA:LIT:0000002\n",
            1,
        )
    )
    one = next(b for b in _review(root).stdout.split("\n\n") if b.startswith("CHDA:LIT:0000001"))
    assert "\n  shares probands with: CHDA:LIT:0000004, CHDA:LIT:0000002\n" in one


def test_approve_sets_both_fields_and_changes_no_other_byte(root: Path, lit_path: Path) -> None:
    before = lit_path.read_text()
    result = _review(
        root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. Grace Hopper", "--on", "2026-11-03"
    )

    assert result.exit_code == 0, result.stdout
    assert result.stdout == "approved CHDA:LIT:0000002 (Dr. Grace Hopper, 2026-11-03)\n"
    old, new = _blocks(before), _blocks(lit_path.read_text())
    assert new["CHDA:LIT:0000002"] == old["CHDA:LIT:0000002"].replace(
        "    verified_by: agent:verify/earlier\n",
        "    verified_by: agent:verify/earlier\n"
        "    approved_by: Dr. Grace Hopper\n    approved_on: 2026-11-03\n"
        f"    approved_digest: {_digest_in(before, 'CHDA:LIT:0000002')}\n",
    )
    for other in ("head", "CHDA:LIT:0000001", "CHDA:LIT:0000003", "CHDA:LIT:0000004"):
        assert new[other] == old[other]


def test_the_digest_review_writes_is_the_one_the_validator_recomputes(
    root: Path, lit_path: Path
) -> None:
    assert _review(root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. A").exit_code == 0
    data = YAML(typ="safe").load(lit_path.read_text())
    (written,) = (f for f in data["findings"] if f["id"] == "CHDA:LIT:0000002")
    assert written["approved_digest"] == _digest_in(lit_path.read_text(), "CHDA:LIT:0000002")


def test_approve_defaults_the_date_to_today_and_takes_several_ids(
    root: Path, lit_path: Path
) -> None:
    before = lit_path.read_text()
    result = _review(
        root,
        "--approve",
        "CHDA:LIT:0000004",
        "--approve",
        "CHDA:LIT:0000002",
        "--by",
        "Dr. Grace Hopper",
        "--allow-unverified",
    )
    assert result.exit_code == 0, result.stdout
    today = date.today().isoformat()
    assert result.stdout == (
        f"approved CHDA:LIT:0000002 (Dr. Grace Hopper, {today})\n"
        f"approved CHDA:LIT:0000004 (Dr. Grace Hopper, {today})\n"
    )
    new = _blocks(lit_path.read_text())
    assert new["CHDA:LIT:0000004"].endswith(
        f"    extracted_by: agent:extract/test\n    approved_by: Dr. Grace Hopper\n"
        f"    approved_on: {today}\n"
        f"    approved_digest: {_digest_in(before, 'CHDA:LIT:0000004')}\n"
    )


@pytest.mark.parametrize(
    ("args", "needle"),
    [
        (["--approve", "CHDA:LIT:0000002", "--by", "agent:review/x"], "agent"),
        (["--approve", "CHDA:LIT:0000002", "--by", "  Agent:review/x"], "agent"),
        (["--approve", "CHDA:LIT:0000002", "--by", "   "], "--by"),
        (["--approve", "CHDA:LIT:0000002"], "--by"),
        (["--approve", "CHDA:LIT:0000002", "--by", "A", "--on", "2026-13-40"], "--on"),
        (["--approve", "CHDA:LIT:0000009", "--by", "A"], "CHDA:LIT:0000009"),
        # M1: a date LIT008 would refuse (search 1 ran 2026-11-02) or one not yet reached.
        (
            [
                "--approve",
                "CHDA:LIT:0000001",
                "--by",
                "A",
                "--allow-unverified",
                "--on",
                "2026-11-01",
            ],
            "LIT008",
        ),
        (["--approve", "CHDA:LIT:0000002", "--by", "A", "--on", "2999-01-01"], "future"),
        (["--approve", "CHDA:LIT:0000003", "--by", "Someone Else"], "already approved"),
        (["--approve", "CHDA:LIT:0000001", "--by", "A"], "unverified"),
        # one bad id refuses the whole request: nothing is written for the good one
        (
            ["--approve", "CHDA:LIT:0000002", "--approve", "CHDA:LIT:0000001", "--by", "A"],
            "unverified",
        ),
        (
            ["--approve", "CHDA:LIT:0000002", "--pending-review", "CHDA:LIT:0000001", "--by", "A"],
            "separate actions; use only one of them",
        ),
        (["--pending-review", "CHDA:LIT:0000001"], "--note"),
        (["--pending-review", "CHDA:LIT:0000001", "--note", "  "], "--note"),
        (["--pending-review", "CHDA:LIT:0000009", "--note", "x"], "CHDA:LIT:0000009"),
        (["--note", "x"], "--pending-review"),
        (["--by", "A"], "--approve"),
    ],
)
def test_every_refusal_exits_two_and_writes_nothing(
    root: Path, lit_path: Path, args: list[str], needle: str
) -> None:
    before = lit_path.read_bytes()
    result = _review(root, *args)
    assert result.exit_code == 2, result.stdout
    assert needle in result.stdout
    assert lit_path.read_bytes() == before


def test_approval_on_the_day_of_the_search_is_accepted(root: Path, lit_path: Path) -> None:
    """LIT008 refuses only a date *before* the search that found the paper."""
    result = _review(
        root,
        "--approve",
        "CHDA:LIT:0000001",
        "--by",
        "A",
        "--allow-unverified",
        "--on",
        "2026-11-02",
    )
    assert result.exit_code == 0, result.stdout


def test_pending_review_sets_the_pair_and_nothing_else(root: Path, lit_path: Path) -> None:
    before = _blocks(lit_path.read_text())
    result = _review(
        root, "--pending-review", "CHDA:LIT:0000003", "--note", "Table 2 says 7; the text says 8."
    )

    assert result.exit_code == 0, result.stdout
    after = _blocks(lit_path.read_text())
    assert after["CHDA:LIT:0000003"] == before["CHDA:LIT:0000003"] + (
        "    review_status: pending_review\n    review_note: Table 2 says 7; the text says 8.\n"
    )
    for other in ("head", "CHDA:LIT:0000001", "CHDA:LIT:0000002", "CHDA:LIT:0000004"):
        assert after[other] == before[other]
    # a second mark would overwrite the first note
    again = _review(root, "--pending-review", "CHDA:LIT:0000003", "--note", "other")
    assert again.exit_code == 2 and "already" in again.stdout


def test_review_without_a_file_names_the_expected_path(root: Path) -> None:
    result = _review(root)
    assert result.exit_code == 2
    assert str(root / "curation" / "literature" / "TBX5.yaml") in result.stdout


def test_a_write_the_model_would_refuse_is_never_made(lit_path: Path) -> None:
    from chd_atlas.gene_literature_io import LiteratureFileError, with_findings_updated

    raw = lit_path.read_bytes()
    with pytest.raises(LiteratureFileError):
        with_findings_updated(raw, {"CHDA:LIT:0000004": {"approved_by": "Dr. A"}})
    with pytest.raises(LiteratureFileError):
        with_findings_updated(
            raw,
            {"CHDA:LIT:0000004": {"approved_by": "agent:x", "approved_on": date(2026, 11, 2)}},
        )
    with pytest.raises(LiteratureFileError):
        with_findings_updated(raw, {"CHDA:LIT:0000099": {"verified_by": "agent:x"}})


def test_the_review_listing_is_ordered_by_id_not_by_file_order(root: Path, lit_path: Path) -> None:
    blocks = _blocks(lit_path.read_text())
    ids = ["CHDA:LIT:0000004", "CHDA:LIT:0000002", "CHDA:LIT:0000001"]
    lit_path.write_text(blocks["head"] + "".join(blocks[i] for i in ids))

    result = _review(root)

    assert result.exit_code == 0, result.stdout
    assert [line[:16] for line in result.stdout.splitlines() if line.startswith("CHDA:")] == [
        "CHDA:LIT:0000001",
        "CHDA:LIT:0000002",
        "CHDA:LIT:0000004",
    ]


# --------------------------------------------------------------------------
# what separates two findings must not end up inside one (fix round 1)
# --------------------------------------------------------------------------

_CHUNKS: dict[str, str] = {
    "publication": "    publication: PMID:{pmid}\n",
    "design": "    design: cohort_screen\n",
    "methods": "    methods:\n      - exome\n",
    "count_unit": "    count_unit: probands\n",
    "screened": "    screened: 7\n",
    "carriers": "    carriers:\n      de_novo: 2\n      inherited: 1\n      origin_unknown: 3\n",
    "variant_classes": "    variant_classes:\n      - snv_indel\n",
    "phenotypes": "    phenotypes:\n      - HP:0001631\n",
    "lesion_groups": "    lesion_groups:\n      - septal\n",
    "syndromic": "    syndromic: both\n",
    "shares_probands_with": "    shares_probands_with: unknown\n",
    "locators": "    locators:\n      - kind: in_text\n        element: Table 2\n",
    "summary": "    summary: Six of seven.\n",
    "extracted_by": "    extracted_by: agent:extract/test\n",
}
_HEADER = (
    "gene: HGNC:11604\nsearches:\n  - id: 1\n    source: pubmed\n    query: TBX5[tiab]\n"
    "    run_on: 2026-11-02\n    hits: 3\n    new_hits: 3\nfindings:\n"
)


def _entry(
    number: int,
    *,
    last: str | None = None,
    verified: bool,
    trailer: str = "",
    eol: str = "",
) -> str:
    """One finding; `last` moves that key to the end (hand-ordered files), `trailer`
    is whatever separates it from the next finding."""
    keys = list(_CHUNKS)
    if last is not None:
        keys.remove(last)
        keys.append(last)
    body = "".join(_CHUNKS[k].format(pmid=1000 + number) for k in keys)
    if eol:  # a same-line comment on the last key, which holds a collection
        body = body.replace(f"    {last}:\n", f"    {last}:  {eol}\n")
    if verified:
        body = body.replace(
            _CHUNKS["extracted_by"], _CHUNKS["extracted_by"] + "    verified_by: agent:v\n"
        )
    return f"  - id: CHDA:LIT:000000{number}\n" + body + trailer


LAYOUTS = {
    "blank line": (None, "\n", ""),
    "comment then blank": (None, "  # note about the next one\n\n", ""),
    "comment": (None, "  # note about the next one\n", ""),
    "list last then comment": ("locators", "  # note about the next one\n", ""),
    "list last then blank": ("locators", "\n", ""),
    "mapping last then comment and blank": ("carriers", "  # note\n\n", ""),
    "same-line comment on a mapping, then blank": ("carriers", "\n", "# from Table 2"),
    "same-line comment on a list, then blank": ("locators", "\n", "# from Table 2"),
}


def _assert_only_lines_added(before: str, after: str, edited: str, added: list[str]) -> None:
    old, new = _blocks(before), _blocks(after)
    for key in old:
        if key != edited:
            assert new[key] == old[key], key
    lines = old[edited].splitlines(keepends=True)
    content = [
        i for i, line in enumerate(lines) if line.strip() and not line.lstrip().startswith("#")
    ]
    cut = content[-1] + 1
    assert new[edited] == "".join(lines[:cut] + added + lines[cut:])


def _layout_file(root: Path, layout: str, *, verified: bool) -> tuple[Path, str]:
    last, trailer, eol = LAYOUTS[layout]
    text = (
        _HEADER
        + _entry(1, verified=verified, trailer="\n")
        + _entry(2, last=last, verified=verified, trailer=trailer, eol=eol)
        + _entry(3, verified=verified)
    )
    path = root / "curation" / "literature" / "TBX5.yaml"
    path.parent.mkdir(exist_ok=True)
    path.write_text(text)
    return path, text


@pytest.mark.parametrize("layout", LAYOUTS)
def test_approve_adds_only_its_lines_whatever_separates_the_findings(
    root: Path, layout: str
) -> None:
    path, before = _layout_file(root, layout, verified=True)
    result = _review(root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. A", "--on", "2026-11-03")
    assert result.exit_code == 0, result.stdout
    _assert_only_lines_added(
        before,
        path.read_text(),
        "CHDA:LIT:0000002",
        [
            "    approved_by: Dr. A\n",
            "    approved_on: 2026-11-03\n",
            f"    approved_digest: {_digest_in(before, 'CHDA:LIT:0000002')}\n",
        ],
    )


@pytest.mark.parametrize("layout", LAYOUTS)
def test_pending_review_adds_only_its_lines_whatever_separates_the_findings(
    root: Path, layout: str
) -> None:
    path, before = _layout_file(root, layout, verified=False)
    result = _review(root, "--pending-review", "CHDA:LIT:0000002", "--note", "Table 2 says 7.")
    assert result.exit_code == 0, result.stdout
    _assert_only_lines_added(
        before,
        path.read_text(),
        "CHDA:LIT:0000002",
        ["    review_status: pending_review\n", "    review_note: Table 2 says 7.\n"],
    )


def test_approving_a_finding_under_review_keeps_the_mark_and_the_note(
    root: Path, lit_path: Path
) -> None:
    before = lit_path.read_text()
    _review(root, "--pending-review", "CHDA:LIT:0000002", "--note", "Table 2 says 7; text says 8.")
    result = _review(root, "--approve", "CHDA:LIT:0000002", "--by", "Dr. A", "--on", "2026-11-03")

    assert result.exit_code == 0, result.stdout
    block = _blocks(lit_path.read_text())["CHDA:LIT:0000002"]
    assert block.endswith(
        "    verified_by: agent:verify/earlier\n"
        "    approved_by: Dr. A\n    approved_on: 2026-11-03\n"
        f"    approved_digest: {_digest_in(before, 'CHDA:LIT:0000002')}\n"
        "    review_status: pending_review\n    review_note: Table 2 says 7; text says 8.\n"
    )


TBX5_LITERAL = (
    '("TBX5"[tiab] OR "HOS"[tiab]) AND ("Heart Defects, Congenital"[mh]'
    ' OR "congenital heart"[tiab] OR "heart defect"[tiab]'
    ' OR "heart defects"[tiab] OR "cardiac defect"[tiab]'
    ' OR "cardiac defects"[tiab] OR "heart malformation"[tiab]'
    ' OR "heart malformations"[tiab] OR "cardiac malformation"[tiab]'
    ' OR "cardiac malformations"[tiab] OR "septal defect"[tiab]'
    ' OR "septal defects"[tiab] OR "tetralogy of Fallot"[tiab]'
    ' OR "transposition of the great arteries"[tiab]'
    ' OR "hypoplastic left heart"[tiab] OR "heterotaxy"[tiab]'
    ' OR "conotruncal"[tiab] OR "outflow tract"[tiab]'
    ' OR "atrioventricular septal defect"[tiab] OR "coarctation"[tiab]'
    ' OR "Holt-Oram syndrome"[tiab]) NOT preprint[pt] NOT review[pt]'
)
BRCA1_LITERAL = (
    '("BRCA1"[tiab]) AND ("Heart Defects, Congenital"[mh]'
    ' OR "congenital heart"[tiab] OR "heart defect"[tiab]'
    ' OR "heart defects"[tiab] OR "cardiac defect"[tiab]'
    ' OR "cardiac defects"[tiab] OR "heart malformation"[tiab]'
    ' OR "heart malformations"[tiab] OR "cardiac malformation"[tiab]'
    ' OR "cardiac malformations"[tiab] OR "septal defect"[tiab]'
    ' OR "septal defects"[tiab] OR "tetralogy of Fallot"[tiab]'
    ' OR "transposition of the great arteries"[tiab]'
    ' OR "hypoplastic left heart"[tiab] OR "heterotaxy"[tiab]'
    ' OR "conotruncal"[tiab] OR "outflow tract"[tiab]'
    ' OR "atrioventricular septal defect"[tiab]'
    ' OR "coarctation"[tiab]) NOT preprint[pt] NOT review[pt]'
)


def _target(root: Path, *args: str) -> Any:
    return runner.invoke(app, ["literature", "target", *args, "--root", str(root)])


def test_target_prints_the_gene_the_gate_and_the_query(root: Path) -> None:
    """The query is a literal (4.38); a former symbol resolves to the approved one."""
    result = _target(root, "HOS")
    assert result.exit_code == 0
    lines = result.output.splitlines()
    assert len(lines) == 1
    printed = json.loads(lines[0])
    assert printed == {
        "format": 1,
        "gene": "HGNC:11604",
        "symbol": "TBX5",
        "published": True,
        "query": TBX5_LITERAL,
    }
    assert printed["published"] is True
    assert printed["format"] == 1 and printed["format"] is not True

    unpublished = json.loads(_target(root, "BRCA1").output)
    assert unpublished["published"] is False
    assert unpublished["gene"] == "HGNC:1100"
    assert unpublished["query"] == BRCA1_LITERAL


def test_target_takes_only_in_scope_disease_names(root: Path) -> None:
    scope = (root / "curation" / "chd_scope.yaml").read_text()
    assert "MONDO:0005021" not in scope
    with (root / "mirrors" / "clingen_gene_validity.tsv").open("a") as handle:
        handle.write(
            "HGNC:11604\tTBX5\tMONDO:0005021\tdilated cardiomyopathy\tAD\tSOP11\tDefinitive\t"
            "2025-03-25T16:00:00.000Z\tSome Expert Panel\thttps://x/3\n"
        )
    result = _target(root, "TBX5")
    assert result.exit_code == 0
    assert "dilated cardiomyopathy" not in result.output
    assert json.loads(result.output)["query"] == TBX5_LITERAL


def _delete(name: str) -> Callable[[Path], None]:
    return lambda root: (root / "mirrors" / name).unlink()


def _header_only(root: Path) -> None:
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(CLINGEN_HEADER)
    (root / "mirrors" / "gencc_submissions.tsv").write_text(GENCC_HEADER)


@pytest.mark.parametrize(
    ("symbol", "damage", "message"),
    [
        ("NOPE", None, "is not a known gene symbol"),
        ("MADH7", None, "HGNC:6772, HGNC:6773"),
        ("TBX5", _delete("genes.tsv"), "mirrors/genes.tsv is missing or unreadable"),
        ("TBX5", _delete("clingen_gene_validity.tsv"), "the validity mirrors are unreadable"),
        ("TBX5", _header_only, "the publication gate is empty"),
    ],
)
def test_target_refusals_exit_2_and_name_the_cause(
    root: Path, symbol: str, damage: Callable[[Path], None] | None, message: str
) -> None:
    if damage is not None:
        damage(root)
    result = _target(root, symbol)
    assert result.exit_code == 2
    assert message in result.output


# --------------------------------------------------------------------------
# `literature frame`
# --------------------------------------------------------------------------

FRAME_GENES = (
    "HGNC:90001\tFRAMEPLUS\tsynthetic plus\tFPALIAS\t\t\t\t\t\tNM_900001.1\n"
    "HGNC:90002\tFRAMEMINUS\tsynthetic minus\t\t\t\t\t\t\tNM_900002.3\n"
    "HGNC:90003\tFRAMEABSENT\tin genes, not in the frame\t\t\t\t\t\t\t\n"
)


MISSING_MIRROR = (
    "the MANE mirror is missing or incomplete in mirrors/ "
    "(needs mane_transcripts.tsv and mane_exons.tsv)"
)


def _frame(root: Path, *args: str) -> Any:
    return runner.invoke(app, ["literature", "frame", *args, "--root", str(root)])


@pytest.fixture
def frame_root(root: Path) -> Path:
    with (root / "mirrors" / "genes.tsv").open("a") as handle:
        handle.write(FRAME_GENES)
    transcripts, exons = synthetic_rows()
    write_mane_mirrors(root, transcripts, exons)
    return root


def test_frame_prints_the_pinned_transcript_and_protein(frame_root: Path) -> None:
    result = _frame(frame_root, "FRAMEMINUS")
    assert result.exit_code == 0
    assert result.output.splitlines() == [
        '{"chrom": "2", "format": 1, "gene": "HGNC:90002", "protein": "NP_211113.9", '
        '"strand": "-", "symbol": "FRAMEMINUS", "transcript": "NM_900002.3"}'
    ]
    plus = json.loads(_frame(frame_root, "FRAMEPLUS").output)
    assert plus["strand"] == "+" and plus["chrom"] == "1"
    assert plus["transcript"] == "NM_900001.1" and plus["protein"] == "NP_211112.9"
    assert (
        json.loads(_frame(frame_root, "FPALIAS").output) == plus
    )  # an alias answers as the approved symbol
    assert plus["format"] == 1 and plus["format"] is not True


@pytest.mark.parametrize(
    ("symbol", "damage", "message"),
    [
        ("NOPE", None, "is not a known gene symbol"),
        ("MADH7", None, "HGNC:6772, HGNC:6773"),
        ("FRAMEPLUS", "no_exons", MISSING_MIRROR),
        ("FRAMEPLUS", "no_transcripts", MISSING_MIRROR),
        ("FRAMEPLUS", "unreadable", "the MANE mirror is unreadable: "),
        ("FRAMEABSENT", None, "HGNC:90003 (FRAMEABSENT) has no row in the MANE mirror"),
    ],
)
def test_frame_refusals_exit_2_and_name_the_cause(
    frame_root: Path, symbol: str, damage: str | None, message: str
) -> None:
    if damage == "no_exons":
        (frame_root / "mirrors" / "mane_exons.tsv").unlink()
    elif damage == "no_transcripts":
        (frame_root / "mirrors" / "mane_transcripts.tsv").unlink()
    elif damage == "unreadable":
        (frame_root / "mirrors" / "mane_transcripts.tsv").write_text("not\ta\tmirror\n")
    result = _frame(frame_root, symbol)
    assert result.exit_code == 2
    assert message in result.output
    if damage == "unreadable":
        assert "mane_transcripts.tsv" in result.output  # the file is named
    assert "{" not in result.output


# --------------------------------------------------------------------------
# --approve-synthesis and the synthesis listing (gene synthesis, Task 6)
# --------------------------------------------------------------------------

ASSERTION_FILE = """\
assertions:
  - id: CHDA:AST:0000001
    gene: HGNC:11604
    phenotypes:
      - HP:0001631
      - HP:0001629
    lesion_groups:
      - septal
      - conotruncal
    inheritance:
      - AD
    mechanism: haploinsufficiency
    syndromic: isolated
    evidence:
      - publication: PMID:222
        evidence_class: genetic_case
        locator:
          kind: in_text
          element: Abstract
        summary: A nonsense variant.
        strength: strong
      - publication: PMID:111
        evidence_class: genetic_case
        locator:
          kind: in_text
          element: Figure 2
        summary: A segregating missense variant.
        strength: moderate
    curator: Dr. Curator
    curated_on: 2026-07-01
    last_reviewed: 2026-07-01
  - id: CHDA:AST:0000002
    gene: HGNC:4173
    phenotypes:
      - HP:0001629
    lesion_groups:
      - conotruncal
    inheritance:
      - AD
    mechanism: dominant_negative
    syndromic: isolated
    evidence:
      - publication: PMID:111
        evidence_class: genetic_case
        locator:
          kind: in_text
          element: Abstract
        summary: A missense variant in another gene.
        strength: moderate
    curator: Dr. Other
    curated_on: 2026-07-02
    last_reviewed: 2026-07-02
"""

SYNTHESIS = """\
synthesis:
  sentences:
    - text: Septal defects are the usual presentation.
      cites:
        - CHDA:LIT:0000003
    - text: A lesion record backs the inheritance pattern.
      cites:
        - CHDA:AST:0000001
  drafted_by: agent:gene-synthesis/test
"""


def _with_synthesis(root: Path, lit_path: Path, block: str = SYNTHESIS) -> str:
    """Write the assertions (TBX5's and GATA4's) and append `block` to the literature file."""
    (root / "curation" / "assertions").mkdir(exist_ok=True)
    (root / "curation" / "assertions" / "TBX5.yaml").write_text(ASSERTION_FILE)
    lit_path.write_text(lit_path.read_text() + block)
    return lit_path.read_text()


def _current_digest(root: Path, lit_path: Path) -> str:
    from chd_atlas.corpus import load_curation
    from chd_atlas.gene_literature_io import parse_literature
    from chd_atlas.models.gene_literature import synthesis_digest

    record = parse_literature(lit_path.read_bytes())
    assert record.synthesis is not None
    own = [a for a in load_curation(root)[0].assertions if str(a.gene) == str(record.gene)]
    return synthesis_digest(record.synthesis, record.findings, own)


def _listed(root: Path, lit_path: Path) -> str:
    """The `--digest` value the listing prints: the current digest's first 12 hex characters."""
    return _current_digest(root, lit_path).removeprefix("sha256:")[:12]


def _approve_synthesis(root: Path, lit_path: Path, *args: str) -> Any:
    """`--approve-synthesis` with the digest the listing would print now."""
    return _review(root, "--approve-synthesis", *args, "--digest", _listed(root, lit_path))


def _approve_in_text(text: str, finding_id: str, by: str, on: str) -> str:
    """`text` with the named finding approved, its digest computed from `text`, not typed."""
    block = _blocks(text)[finding_id]  # the last finding's block runs to the end of the file
    anchor = "    extracted_by: agent:extract/test\n"
    approval = (
        f"    approved_by: {by}\n    approved_on: {on}\n"
        f"    approved_digest: {_digest_in(text, finding_id)}\n"
    )
    approved = block.replace(anchor, anchor + approval, 1)
    assert approved != block and text.count(block) == 1
    return text.replace(block, approved)


def test_approve_synthesis_stamps_three_fields_and_changes_no_other_byte(
    root: Path, lit_path: Path
) -> None:
    before = _with_synthesis(root, lit_path)
    result = _approve_synthesis(root, lit_path, "--by", "Dr. Grace Hopper", "--on", "2026-11-03")

    assert result.exit_code == 0, result.stdout
    assert result.stdout == "approved synthesis (Dr. Grace Hopper, 2026-11-03)\n"
    digest = _current_digest(root, lit_path)
    assert lit_path.read_text() == before + (
        f"  approved_by: Dr. Grace Hopper\n  approved_on: 2026-11-03\n  approved_digest: {digest}\n"
    )
    # the digest is the one the validator recomputes: no SYN004 afterwards
    from chd_atlas.validate.runner import validate_repository

    assert not [i for i in validate_repository(root).issues if i.code.startswith("SYN")]


def test_a_stale_approval_is_re_approved_and_a_current_one_is_never_overwritten(
    root: Path, lit_path: Path
) -> None:
    _with_synthesis(
        root,
        lit_path,
        SYNTHESIS
        + "  approved_by: Dr. First Approver\n  approved_on: 2026-10-05\n"
        + "  approved_digest: sha256:stale\n",
    )
    stale = _approve_synthesis(root, lit_path, "--by", "Dr. Second", "--on", "2026-11-04")
    assert stale.exit_code == 0, stale.stdout
    after = lit_path.read_text()
    assert "approved_by: Dr. Second\n  approved_on: 2026-11-04\n" in after
    assert "Dr. First Approver" not in after
    assert f"approved_digest: {_current_digest(root, lit_path)}\n" in after
    from chd_atlas.validate.runner import validate_repository

    assert not [i for i in validate_repository(root).issues if i.code.startswith("SYN")]

    again = _approve_synthesis(root, lit_path, "--by", "Dr. Third")
    assert again.exit_code == 2
    assert "already approved by Dr. Second; approving again would erase who approved it" in (
        again.stdout
    )
    assert lit_path.read_text() == after


@pytest.mark.parametrize(
    ("block", "args", "needle"),
    [
        (None, ["--by", "A"], "no synthesis"),
        (SYNTHESIS, [], "--by"),
        (SYNTHESIS, ["--by", "agent:x"], "agent"),
        (SYNTHESIS, ["--by", "A", "--on", "2999-01-01"], "future"),
        (SYNTHESIS, ["--by", "A", "--allow-unverified"], "--allow-unverified"),
        (
            SYNTHESIS,
            ["--by", "A", "--approve", "CHDA:LIT:0000002"],
            "--approve, --pending-review and --approve-synthesis are separate actions; "
            "use only one of them",
        ),
        (SYNTHESIS, ["--by", "A", "--pending-review", "CHDA:LIT:0000001"], "separate actions"),
        (
            SYNTHESIS,
            ["--by", "A", "--against", "x.yaml"],
            "--against only applies to the listing, not to an action",
        ),
        (
            SYNTHESIS,
            ["--by", "A", "--on", "2026-09-30"],
            "--on 2026-09-30 is before 2026-10-01, when a finding this synthesis covers was "
            "approved; an approval cannot predate what it covers",
        ),
        (SYNTHESIS.replace("usual presentation", "usual in 12 families"), ["--by", "A"], "SYN002"),
        (
            SYNTHESIS.replace("usual presentation", "usual " + "word " * 120),
            ["--by", "A"],
            "SYN003",
        ),
        (SYNTHESIS.replace("CHDA:LIT:0000003", "CHDA:LIT:0000004"), ["--by", "A"], "SYN001"),
    ],
)
def test_every_synthesis_refusal_exits_two_and_writes_nothing(
    root: Path, lit_path: Path, block: str | None, args: list[str], needle: str
) -> None:
    """Each row passes the digest the listing would print, so the refusal is its own."""
    digest = "0" * 12
    if block is not None:
        _with_synthesis(root, lit_path, block)
        digest = _listed(root, lit_path)
    before = lit_path.read_bytes()
    result = _review(root, "--approve-synthesis", *args, "--digest", digest)
    assert result.exit_code == 2, result.stdout
    assert needle in result.stdout
    assert lit_path.read_bytes() == before


def test_approve_synthesis_requires_the_digest_the_listing_printed(
    root: Path, lit_path: Path
) -> None:
    """K-3: the person approves the text they listed, not whatever is on disk by then."""
    _with_synthesis(root, lit_path)
    listed = _current_digest(root, lit_path).removeprefix("sha256:")
    before = lit_path.read_bytes()

    def refused(*args: str) -> str:
        result = _review(root, "--approve-synthesis", "--by", "Dr. A", *args)
        assert result.exit_code == 2, result.stdout
        assert lit_path.read_bytes() == before
        return str(result.stdout)

    assert refused() == (
        "error: --digest is required with --approve-synthesis: run "
        "`chd-atlas literature review TBX5` and copy the command it prints\n"
    )
    not_listed = (
        f"error: the synthesis is not the one listed (its digest is {listed[:12]}); "
        "list it again before approving\n"
    )
    assert refused("--digest", listed[:11]) == not_listed  # a true prefix, but too short
    assert refused("--digest", "f" * 12) == not_listed
    assert refused("--digest", "sha256:" + listed) == not_listed
    # an agent redrafts between the listing and the approval: the listed digest no longer fits
    lit_path.write_text(lit_path.read_text().replace("usual presentation", "usual finding"))
    before = lit_path.read_bytes()
    moved = _current_digest(root, lit_path).removeprefix("sha256:")
    assert refused("--digest", listed[:12]) == (
        f"error: the synthesis is not the one listed (its digest is {moved[:12]}); "
        "list it again before approving\n"
    )
    # `--digest` is an option of the action, like `--note`
    plain = _review(root, "--digest", moved[:12])
    assert plain.exit_code == 2
    assert plain.stdout == "error: --digest only applies with --approve-synthesis\n"
    # the whole digest, or any prefix of at least 12, is accepted
    whole = _review(root, "--approve-synthesis", "--by", "Dr. A", "--digest", moved)
    assert whole.exit_code == 0, whole.stdout


def test_approve_synthesis_refuses_while_an_approved_finding_fails_lit015(
    root: Path, lit_path: Path
) -> None:
    """X-1a: an approval over a finding whose own approval is void is refused, and marked.

    Finding 3 is cited and finding 4 is not: the guard covers every approved
    finding of the file, as the synthesis digest does. Their approval dates
    differ (2026-10-01, 2026-10-15), so the listing's mark is per finding.
    """
    _with_synthesis(
        root,
        lit_path,
        SYNTHESIS.replace(
            "        - CHDA:LIT:0000003\n",
            "        - CHDA:LIT:0000003\n        - CHDA:LIT:0000002\n",
        ),
    )
    text = _approve_in_text(lit_path.read_text(), "CHDA:LIT:0000004", "Dr. Bee", "2026-10-15")
    text = _approve_in_text(text, "CHDA:LIT:0000002", "Dr. Cee", "2026-10-16")
    lit_path.write_text(text)
    clean = _approve_synthesis(root, lit_path, "--by", "Dr. A", "--on", "2026-11-03")
    assert clean.exit_code == 0, clean.stdout  # the fixture approves cleanly before the edits
    lit_path.write_text(text)

    edited = text.replace("summary: A proband.\n", "summary: An edited proband.\n").replace(
        "summary: Another proband.\n", "summary: Another, edited.\n"
    )
    lit_path.write_text(edited)
    result = _approve_synthesis(root, lit_path, "--by", "Dr. A", "--on", "2026-11-03")
    assert result.exit_code == 2, result.stdout
    assert result.stdout == (
        "error: cannot approve the synthesis: CHDA:LIT:0000003, CHDA:LIT:0000004 fail(s) "
        "LIT015 (approved content changed since approval); re-approve or correct them first\n"
    )
    assert lit_path.read_text() == edited

    # the listing marks the void approval, ahead of a review mark; a current one shows none
    marked = edited.replace(
        "summary: An edited proband.\n",
        "summary: An edited proband.\n    review_status: pending_review\n"
        "    review_note: The table disagrees.\n",
    )
    lit_path.write_text(marked)
    out = _review(root).stdout
    assert "     CHDA:LIT:0000003: An edited proband. [approval void: LIT015]\n" in out
    assert "     CHDA:LIT:0000002: A GWAS hit.\n" in out


def test_a_synthesis_approval_may_not_predate_the_latest_finding_approval(
    root: Path, lit_path: Path
) -> None:
    """Two approval dates are needed: with one, the earliest and latest are the same date,
    so `min` and `max` in the guard cannot be told apart."""
    from chd_atlas.gene_literature_io import parse_literature

    _with_synthesis(root, lit_path)
    lit_path.write_text(
        _approve_in_text(lit_path.read_text(), "CHDA:LIT:0000004", "Dr. Later", "2026-10-20")
    )
    assert parse_literature(lit_path.read_bytes()).findings  # the edit is still a valid file

    before = lit_path.read_bytes()
    early = _approve_synthesis(root, lit_path, "--by", "Dr. Grace Hopper", "--on", "2026-10-10")
    assert early.exit_code == 2, early.stdout
    assert (
        "--on 2026-10-10 is before 2026-10-20, when a finding this synthesis covers was "
        "approved; an approval cannot predate what it covers"
    ) in early.stdout
    assert lit_path.read_bytes() == before

    equal = _approve_synthesis(root, lit_path, "--by", "Dr. Grace Hopper", "--on", "2026-10-20")
    assert equal.exit_code == 0, equal.stdout
    assert equal.stdout == "approved synthesis (Dr. Grace Hopper, 2026-10-20)\n"


@pytest.mark.parametrize(
    ("args", "builder"),
    [
        (["--approve", "CHDA:LIT:0000002", "--by", "Dr. A"], "with_findings_updated"),
        (["--pending-review", "CHDA:LIT:0000003", "--note", "Unclear."], "with_findings_updated"),
        (["--approve-synthesis", "--by", "Dr. A"], "with_synthesis_approval"),
    ],
)
def test_a_write_racing_another_writer_is_refused_and_keeps_its_bytes(
    root: Path, lit_path: Path, monkeypatch: pytest.MonkeyPatch, args: list[str], builder: str
) -> None:
    """K-1: the file changes after it was read, while the new content is built.

    The other writer lands inside the call that builds the bytes to write, the
    last step before the write, so a check made any earlier cannot see it.
    """
    import chd_atlas.cli as cli

    _with_synthesis(root, lit_path)
    if "--approve-synthesis" in args:
        args = [*args, "--digest", _listed(root, lit_path)]
    theirs = lit_path.read_bytes() + b"# another writer's line\n"
    original = getattr(cli, builder)

    def racing(*a: Any, **k: Any) -> bytes:
        content: bytes = original(*a, **k)
        lit_path.write_bytes(theirs)
        return content

    monkeypatch.setattr(cli, builder, racing)
    result = _review(root, *args)
    assert result.exit_code == 2, result.stdout
    assert result.stdout == f"error: {lit_path} changed since it was read; nothing was written\n"
    assert lit_path.read_bytes() == theirs


def test_the_listing_prints_the_synthesis_state_sentences_and_what_each_cites(
    root: Path, lit_path: Path
) -> None:
    plain = _review(root)
    assert "synthesis" not in plain.stdout  # a file with none prints nothing about one

    _with_synthesis(
        root,
        lit_path,
        SYNTHESIS.replace(
            "        - CHDA:LIT:0000003\n",
            "        - CHDA:LIT:0000003\n        - CHDA:LIT:0000009\n",
        ),
    )
    tail = _review(root).stdout.split("synthesis: ", 1)[1]
    assert tail == (
        "unapproved\n"
        "  approve with: chd-atlas literature review TBX5 --approve-synthesis --by <name> "
        f"--digest {_listed(root, lit_path)}\n"
        "  1. Septal defects are the usual presentation.\n"
        "     CHDA:LIT:0000003: A proband.\n"
        "     CHDA:LIT:0000009: (not found)\n"
        "  2. A lesion record backs the inheritance pattern.\n"
        "     CHDA:AST:0000001: lesion assertion for HGNC:11604\n"
        "       lesion groups: septal, conotruncal\n"
        "       phenotypes: HP:0001631, HP:0001629\n"
        "       evidence: PMID:222 — A nonsense variant.\n"
        "       evidence: PMID:111 — A segregating missense variant.\n"
    )
    # the approved symbol, whichever name the person typed
    alias = runner.invoke(app, ["literature", "review", "HOS", "--root", str(root)])
    assert "  approve with: chd-atlas literature review TBX5 --approve-synthesis" in alias.stdout
    # an unapproved finding's summary is shown too, so the approver sees what is cited
    lit_path.write_text(lit_path.read_text().replace("- CHDA:LIT:0000009", "- CHDA:LIT:0000004"))
    assert "     CHDA:LIT:0000004: Another proband. [unapproved]\n" in _review(root).stdout

    both = "        - CHDA:LIT:0000003\n        - CHDA:LIT:0000004\n"
    lit_path.write_text(lit_path.read_text().replace(both, "        - CHDA:LIT:0000003\n"))
    done = _approve_synthesis(root, lit_path, "--by", "Dr. A", "--on", "2026-11-03")
    assert done.exit_code == 0, done.stdout
    current = _review(root).stdout
    assert "synthesis: approved by Dr. A on 2026-11-03\n  1. Septal" in current
    assert "approve with" not in current  # nothing to approve
    # a changed covered record voids the approval: STALE, with who and when
    lit_path.write_text(lit_path.read_text().replace("A proband.", "A different proband."))
    assert (
        "synthesis: STALE (approved by Dr. A on 2026-11-03; "
        "the synthesis or the records it covers changed since)\n"
        "  approve with: chd-atlas literature review TBX5 --approve-synthesis --by <name> "
        f"--digest {_listed(root, lit_path)}\n" in _review(root).stdout
    )


def test_a_synthesis_write_the_model_would_refuse_is_never_made(lit_path: Path) -> None:
    from chd_atlas.gene_literature_io import LiteratureFileError, with_synthesis_approval

    approval = {"approved_by": "Dr. A", "approved_on": date(2026, 11, 2), "approved_digest": "d"}
    with pytest.raises(LiteratureFileError, match="no synthesis block"):
        with_synthesis_approval(lit_path.read_bytes(), approval)
    with_block = lit_path.read_bytes() + SYNTHESIS.encode()
    with pytest.raises(LiteratureFileError):
        with_synthesis_approval(with_block, {**approval, "approved_by": "agent:x"})
    with pytest.raises(LiteratureFileError):
        with_synthesis_approval(with_block, {"approved_by": "Dr. A"})


def test_the_listing_marks_a_cited_findings_review_state(root: Path, lit_path: Path) -> None:
    """A cited finding shows `[pending_review]` or `[not_accessible]`; a plain one shows nothing."""
    block = SYNTHESIS.replace(
        "        - CHDA:LIT:0000003\n",
        "        - CHDA:LIT:0000003\n        - CHDA:LIT:0000004\n",
    )
    _with_synthesis(root, lit_path, block)
    marked = lit_path.read_text().replace(
        "    summary: A proband.\n",
        "    summary: A proband.\n    review_status: pending_review\n"
        "    review_note: The table disagrees.\n",
    )
    lit_path.write_text(marked)
    out = _review(root).stdout
    assert "     CHDA:LIT:0000003: A proband. [pending_review]\n" in out
    assert "     CHDA:LIT:0000004: Another proband. [unapproved]\n" in out
    lit_path.write_text(
        marked.replace("review_status: pending_review", "review_status: not_accessible")
    )
    assert "     CHDA:LIT:0000003: A proband. [not_accessible]\n" in _review(root).stdout


def test_approval_flags_without_an_approval_are_refused(root: Path, lit_path: Path) -> None:
    result = _review(root, "--by", "A")
    assert result.exit_code == 2
    assert (
        "--by and --on only apply with --approve or --approve-synthesis, "
        "and --allow-unverified only with --approve"
    ) in result.stdout


# --------------------------------------------------------------------------
# variants in the review listing, and re-approval of a void approval (V5)
# --------------------------------------------------------------------------

_TWO_VARIANTS = """\
    variants:
      - id: 1
        as_reported: c.100A>G
        locator: {kind: in_text, element: Abstract}
        kind: snv_indel
        placement: mane
        hgvs_c: c.100A>G
        hgvs_p: p.Lys34Glu
        grch38: 12:114000000:A:G
        zygosity: heterozygous
        carriers: {de_novo: 1, inherited: 0, origin_unknown: 0}
        phenotypes:
          - HP:0001631
      - id: 2
        as_reported: p.Arg5*
        locator: {kind: in_text, element: Abstract}
        kind: snv_indel
        placement: unplaced
        zygosity: heterozygous
        with: [3, 4]
        phase: trans
        carriers: {de_novo: 0, inherited: 2, origin_unknown: 3}
        review_status: pending_review
        review_note: the paper's c. name disagrees
"""


def _give_variants(lit_path: Path, finding_id: str, summary: str) -> None:
    text = lit_path.read_text()
    marker = f"    summary: {summary}\n"
    assert text.count(marker) == 1
    lit_path.write_text(text.replace(marker, _TWO_VARIANTS + marker))
    assert finding_id in lit_path.read_text()


def test_a_listing_shows_one_line_per_variant(root: Path, lit_path: Path) -> None:
    _give_variants(lit_path, "CHDA:LIT:0000004", "Another proband.")

    result = _review(root)

    assert result.exit_code == 0, result.stdout
    four = next(b for b in result.stdout.split("\n\n") if b.startswith("CHDA:LIT:0000004"))
    assert four.splitlines()[5:9] == [
        "  variant classes: snv_indel",
        "  variant 1: c.100A>G | mane | c.100A>G p.Lys34Glu 12:114000000:A:G"
        " | heterozygous | carriers 1/0/0 | HP:0001631"
        " | kind snv_indel | reported none | locator in_text Abstract",
        "  variant 2: p.Arg5* | unplaced | - - - | heterozygous with 3, 4 (trans)"
        " | carriers 0/2/3 | no phenotypes | kind snv_indel | reported none"
        " | locator in_text Abstract"
        " | review: pending_review: the paper's c. name disagrees",
        "  screened 1 probands; carriers 1 (de novo 0, inherited 1, origin unknown 0)",
    ]
    one = next(b for b in result.stdout.split("\n\n") if b.startswith("CHDA:LIT:0000001"))
    assert "variant 1:" not in one


def test_the_variant_line_pins_missing_carriers_and_a_with_clause_without_a_phase(
    root: Path, lit_path: Path
) -> None:
    text = lit_path.read_text()
    marker = "    summary: Another proband.\n"
    bare = """\
    variants:
      - id: 1
        as_reported: c.5del
        locator: {kind: in_text, element: Abstract}
        kind: snv_indel
        placement: unplaced
        with: [2]
"""
    lit_path.write_text(text.replace(marker, bare + marker))

    result = _review(root)

    assert result.exit_code == 0, result.stdout
    assert (
        "  variant 1: c.5del | unplaced | - - - | not_stated with 2"
        " | carriers not stated | no phenotypes | kind snv_indel | reported none"
        " | locator in_text Abstract"
    ) in result.stdout.splitlines()


def test_the_variant_line_shows_kind_reported_names_and_locator(root: Path, lit_path: Path) -> None:
    """The digest covers `kind`, `reported` and `locator`, so the listing must show them."""
    text = lit_path.read_text()
    marker = "    summary: Another proband.\n"
    one = """\
    variants:
      - id: 1
        as_reported: Arg5X
        locator: {kind: supplementary, file: s1.xlsx, sheet: Variants, row: 7}
        kind: cnv
        placement: unplaced
        reported: {transcript: NM_000001.1, p: p.R5X, rsid: rs123}
"""
    lit_path.write_text(text.replace(marker, one + marker))

    result = _review(root)

    assert result.exit_code == 0, result.stdout
    assert (
        "  variant 1: Arg5X | unplaced | - - - | not_stated | carriers not stated"
        " | no phenotypes | kind cnv | reported transcript NM_000001.1, p p.R5X, rsid rs123"
        " | locator supplementary s1.xlsx / Variants / row 7"
    ) in result.stdout.splitlines()


def test_a_void_finding_is_listed_marked_and_can_be_approved_again(
    root: Path, lit_path: Path
) -> None:
    _give_variants(lit_path, "CHDA:LIT:0000003", "A proband.")
    void_header = (
        "CHDA:LIT:0000003  https://pubmed.ncbi.nlm.nih.gov/222/  case_report  UNVERIFIED"
        "  APPROVAL VOID (LIT015)"
    )

    listed = _review(root)
    assert listed.exit_code == 0, listed.stdout
    assert void_header in listed.stdout.splitlines()
    assert "APPROVAL VOID" not in "".join(
        b for b in listed.stdout.split("\n\n") if not b.startswith("CHDA:LIT:0000003")
    )

    done = _review(
        root, "--approve", "CHDA:LIT:0000003", "--by", "Dr. New", "--on", "2026-11-04",
        "--allow-unverified",
    )  # fmt: skip
    assert done.exit_code == 0, done.stdout
    after = lit_path.read_text()
    assert "Dr. Ada Curator" not in after
    block = _blocks(after)["CHDA:LIT:0000003"]
    assert "approved_by: Dr. New\n" in block
    assert "approved_on: 2026-11-04\n" in block
    assert done.stdout == (
        "re-approved CHDA:LIT:0000003 (Dr. New, 2026-11-04; was Dr. Ada Curator on 2026-10-01)\n"
    )
    current = _finding_on_disk(lit_path, "CHDA:LIT:0000003")
    assert current.approved_digest == finding_digest(current)
    assert "APPROVAL VOID" not in _review(root).stdout
    headers = [line for line in _review(root).stdout.splitlines() if line.startswith("CHDA:LIT:")]
    assert [h[:16] for h in headers] == ["CHDA:LIT:0000001", "CHDA:LIT:0000002", "CHDA:LIT:0000004"]


def test_re_approving_a_void_finding_is_refused_if_the_file_changed_since_it_was_read(
    root: Path, lit_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import chd_atlas.cli as cli

    _give_variants(lit_path, "CHDA:LIT:0000003", "A proband.")
    theirs = lit_path.read_bytes() + b"# another writer's line\n"
    original = cli.with_findings_updated

    def racing(*a: Any, **k: Any) -> bytes:
        content: bytes = original(*a, **k)
        lit_path.write_bytes(theirs)
        return content

    monkeypatch.setattr(cli, "with_findings_updated", racing)
    result = _review(root, "--approve", "CHDA:LIT:0000003", "--by", "Dr. New", "--allow-unverified")
    assert result.exit_code == 2, result.stdout
    assert result.stdout == f"error: {lit_path} changed since it was read; nothing was written\n"
    assert lit_path.read_bytes() == theirs


def _finding_on_disk(path: Path, finding_id: str) -> Finding:
    data = YAML(typ="safe").load(path.read_text())
    return next(Finding.model_validate(f) for f in data["findings"] if f["id"] == finding_id)


def test_a_current_approval_is_still_refused(root: Path, lit_path: Path) -> None:
    before = lit_path.read_text()

    result = _review(root, "--approve", "CHDA:LIT:0000003", "--by", "Someone Else")

    assert result.exit_code == 2
    assert (
        "CHDA:LIT:0000003 is already approved by Dr. Ada Curator; "
        "approving again would erase who approved it"
    ) in result.output
    assert lit_path.read_text() == before
