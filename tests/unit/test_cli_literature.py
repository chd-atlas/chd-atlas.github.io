# tests/unit/test_cli_literature.py
"""`chd-atlas literature review | target`.

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
            "not both",
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
