"""The curator plugin's sentences, held to the code they describe (CLAUDE.md 4.20, 4.25).

`curator/` is plain Markdown that tells an agent which commands to run, which
options they take, which validator codes catch which mistake and which
vocabulary values a record may hold. Nothing executes it, so a misspelt option
or a code that does not exist would pass every other test and mislead the next
agent. Every name below is read out of the code at test time -- the typer app,
the validator source, the enums, the models -- rather than from a list kept
here, with two deliberate exceptions: the design -> StudyType table and the
list of skills are literals, because they are the contract under test.
"""

from __future__ import annotations

import json
import re
import textwrap
from collections.abc import Iterator
from enum import StrEnum
from pathlib import Path
from typing import Any

import pytest
import typer
from ruamel.yaml import YAML

from chd_atlas import vocab
from chd_atlas.cli import app
from chd_atlas.models.assertion import InTextLocator, SupplementaryLocator
from chd_atlas.models.gene_literature import Carriers, Finding, Search, TriageEntry
from chd_atlas.models.literature import Publication
from chd_atlas.models.review import ReviewStatus

REPO = Path(__file__).parent.parent
CURATOR = REPO / "curator"
SKILLS = ("literature-search", "literature-triage", "literature-extract", "literature-verify")
SKILL_FILES = [CURATOR / "skills" / name / "SKILL.md" for name in SKILLS]
ALL_FILES = [*SKILL_FILES, CURATOR / "curator-rules.md", CURATOR / "README.md"]

# The contract: the mapping `literature-extract` tells an agent to apply.
DESIGN_TO_STUDY_TYPE = [
    ("case_report", "case_report"),
    ("family_segregation", "family_linkage"),
    ("cohort_screen", "cohort_sequencing"),
    ("rare_variant_association with control_cohort", "case_control"),
    ("rare_variant_association with mutation_model and exome", "trio_wes"),
    ("rare_variant_association with mutation_model and genome", "trio_wgs"),
    (
        "rare_variant_association with mutation_model and any other method",
        "ask the curator; no StudyType is guessed",
    ),
    ("common_variant_association with genotyping_array", "gwas"),
    ("common_variant_association otherwise", "candidate_variant"),
]

# Options a file may name that belong to no literature command.
NON_LITERATURE_OPTIONS = {"--plugin-dir"}


def _text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def _literature_commands() -> dict[str, Any]:
    # typer 0.27 vendors click, so the command objects are duck-typed here.
    group = typer.main.get_command(app)
    return dict(group.commands["literature"].commands)  # type: ignore[attr-defined]


def _options(command: Any) -> set[str]:
    return {opt for param in command.params for opt in getattr(param, "opts", [])}


def _enum_values() -> dict[str, set[str]]:
    found = {
        name: {m.value for m in member}
        for name, member in vars(vocab).items()
        if isinstance(member, type) and issubclass(member, StrEnum) and member is not StrEnum
    }
    return {**found, "ReviewStatus": {m.value for m in ReviewStatus}}


def _all_enum_values() -> set[str]:
    return set().union(*_enum_values().values())


def _model_fields() -> set[str]:
    return (
        set(Finding.model_fields)
        | set(Search.model_fields)
        | set(TriageEntry.model_fields)
        | set(Carriers.model_fields)
        | set(Publication.model_fields)
        | set(InTextLocator.model_fields)
        | set(SupplementaryLocator.model_fields)
        | {"findings", "triage", "searches", "gene"}
    )


def _command_lines(text: str) -> Iterator[str]:
    """Lines naming a literature command; backslash continuations are joined."""
    joined = re.sub(r"\\\n\s*", " ", text)
    for line in joined.splitlines():
        if "chd-atlas literature " in line:
            yield line


def _fenced_yaml(text: str) -> list[Any]:
    # Fences inside a numbered step are indented; dedent before parsing.
    blocks = re.findall(r"^[ ]*```yaml\n(.*?)^[ ]*```", text, flags=re.S | re.M)
    return [YAML(typ="safe").load(textwrap.dedent(block)) for block in blocks]


def _walk(node: Any) -> Iterator[tuple[str, Any]]:
    if isinstance(node, dict):
        for key, value in node.items():
            yield str(key), value
            yield from _walk(value)
    elif isinstance(node, list):
        for item in node:
            yield from _walk(item)


def test_the_plugin_manifest_names_the_plugin() -> None:
    manifest = json.loads(_text(CURATOR / ".claude-plugin" / "plugin.json"))
    assert manifest["name"] == "chd-atlas-curator"
    assert manifest["version"] == "0.1.0"
    assert manifest["description"].strip()
    assert not (CURATOR / ".claude-plugin" / "marketplace.json").exists()


def test_the_plugin_ships_exactly_these_skills() -> None:
    found = sorted(p.parent.name for p in (CURATOR / "skills").glob("*/SKILL.md"))
    assert found == sorted(SKILLS)
    assert all(path.is_file() for path in ALL_FILES)


@pytest.mark.parametrize("path", SKILL_FILES, ids=SKILLS)
def test_a_skill_has_its_frontmatter_its_sections_and_its_hand_offs(path: Path) -> None:
    text = _text(path)
    match = re.match(r"---\nname: (\S+)\ndescription: (.+)\n---\n", text)
    assert match is not None, "frontmatter must be exactly `name` then `description`"
    assert match.group(1) == path.parent.name
    assert match.group(2).startswith("Use when ")
    for heading in ("## Inputs", "## Steps", "## Done when", "## Never"):
        assert len(re.findall(rf"^{heading}$", text, flags=re.M)) == 1, heading
    assert len(text.splitlines()) <= 120
    assert "curator/curator-rules.md" in text
    # Every skill ends by validating and handing the work to the human.
    tail = text[text.rindex("## Done when") :]
    assert ".venv/bin/chd-atlas validate" in tail
    assert ".venv/bin/chd-atlas literature review <SYMBOL>" in tail


def test_every_command_and_option_the_files_name_exists_in_the_cli() -> None:
    commands = _literature_commands()
    every_option = set().union(*(_options(c) for c in commands.values()))
    named_commands: set[str] = set()
    for path in ALL_FILES:
        text = _text(path)
        for line in _command_lines(text):
            for form in re.finditer(r"chd-atlas literature (\S+)", line):
                name = form.group(1).strip("`.,;:)")
                assert name in commands, f"{path.name}: no command `literature {name}`"
                named_commands.add(name)
                rest = line.split(form.group(0), 1)[1]
                for option in re.findall(r"(?<![\w-])(--[a-z][a-z-]*)", rest):
                    if option in NON_LITERATURE_OPTIONS:
                        continue
                    assert option in _options(commands[name]), f"{path.name}: `{name}` {option}"
        # Prose mentions of `literature <cmd>` and bare `--option`s.
        for name in re.findall(r"`literature ([a-z][a-z-]*)", text):
            assert name in commands, f"{path.name}: no command `literature {name}`"
        for option in re.findall(r"(?<![\w-])(--[a-z][a-z-]+)", text):
            assert option in every_option | NON_LITERATURE_OPTIONS, f"{path.name}: {option}"
    # Each command a curator runs is taught somewhere, so a removed one is noticed.
    assert named_commands == set(commands)


def test_every_code_the_files_name_is_raised_by_the_code() -> None:
    source = "".join(
        p.read_text(encoding="utf-8") for p in (REPO / "src" / "chd_atlas").rglob("*.py")
    )
    validator = (REPO / "src/chd_atlas/validate/gene_literature.py").read_text(encoding="utf-8")
    literature_codes = set(re.findall(r'"(LIT\d{3})"', validator))
    for path in ALL_FILES:
        for code in re.findall(r"\b([A-Z]{2,6}\d{3})\b", _text(path)):
            assert f'"{code}"' in source, f"{path.name}: {code} is raised nowhere"
            if code.startswith("LIT"):
                assert code in literature_codes, f"{path.name}: {code}"
    # The rules cite every literature code, so a validator added later is noticed.
    rules = set(re.findall(r"\b(LIT\d{3})\b", _text(CURATOR / "curator-rules.md")))
    assert rules == literature_codes


def test_the_yaml_examples_hold_real_vocabulary_and_validate_as_records() -> None:
    enums = _enum_values()
    by_field = {
        "design": enums["FindingDesign"],
        "methods": enums["DetectionMethod"],
        "count_unit": enums["ScreenUnit"],
        "variant_classes": enums["FindingVariantClass"],
        "comparator": enums["BurdenComparator"],
        "decision": enums["TriageDecision"],
        "reason": enums["ExclusionReason"],
        "syndromic": enums["SyndromicStatus"],
        "lesion_groups": enums["LesionGroup"],
        "review_status": enums["ReviewStatus"],
    }
    seen_fields: set[str] = set()
    for path in ALL_FILES:
        for block in _fenced_yaml(_text(path)):
            for key, value in _walk(block):
                assert key in _model_fields(), f"{path.name}: unknown field {key}"
                seen_fields.add(key)
                if key in by_field:
                    values = value if isinstance(value, list) else [value]
                    assert set(values) <= by_field[key], f"{path.name}: {key}: {values}"
            if isinstance(block, dict) and "findings" in block:
                for record in block["findings"]:
                    Finding.model_validate(record)
            if isinstance(block, list):
                for record in block:
                    if "pmid" in record:
                        TriageEntry.model_validate(record)
                    if "publication" in record:
                        Finding.model_validate(record)
    # The examples cover both record kinds, so they are not vacuous.
    assert {"findings", "decision", "design", "locators"} <= seen_fields


def test_backticked_names_are_real_values_fields_or_options() -> None:
    parameters = {
        param.name
        for command in _literature_commands().values()
        for param in command.params
        if param.name
    }
    known = _all_enum_values() | _model_fields() | parameters
    for path in ALL_FILES:
        # Fenced blocks are checked structurally above.
        prose = re.sub(r"^[ ]*```.*?^[ ]*```", "", _text(path), flags=re.S | re.M)
        for token in re.findall(r"`([a-z]+(?:_[a-z]+)+)`", prose):
            assert token in known, f"{path.name}: `{token}` is not a value, field or option"


def test_every_study_type_the_files_name_exists_and_the_table_is_the_contract() -> None:
    study_types = {m.value for m in vocab.StudyType}
    text = _text(CURATOR / "skills" / "literature-extract" / "SKILL.md")
    # A study cell is a backticked StudyType, or the one plain-text row that tells the
    # agent to stop and ask instead of guessing one.
    rows = re.findall(
        r"^\s*\| (`[^|]+?`(?: with `[^|]+?`(?: and [^|`]+?)?| otherwise)?) \| "
        r"(?:`([^|]+?)`|([^|`]+?)) \|$",
        text,
        re.M,
    )
    table = [(design.replace("`", ""), study or ask) for design, study, ask in rows]
    assert table == DESIGN_TO_STUDY_TYPE
    assert {study for _, study in table} <= study_types | {
        "ask the curator; no StudyType is guessed"
    }
    for path in ALL_FILES:
        for study in re.findall(r"--study-type (\S+)", _text(path)):
            assert study.strip("`") in study_types | {"<StudyType>"}, path.name


def test_the_files_state_the_verified_rulings_literally() -> None:
    extract = _text(CURATOR / "skills" / "literature-extract" / "SKILL.md")
    verify = _text(CURATOR / "skills" / "literature-verify" / "SKILL.md")
    rules = _text(CURATOR / "curator-rules.md")
    assert "local/literature/<SYMBOL>.verify.yaml" in verify
    assert "--mark-verified agent:literature-verify/<model>" in verify
    assert "curation/.id_registry.yaml" in extract
    assert "agent:literature-extract/<model>" in extract
    assert "agent:<skill>/<model>" in rules
    # R26: the verifier writes `unknown`; the curator reads the field in the listing.
    assert "always set `shares_probands_with: unknown`" in verify
    assert "reading this field in the `chd-atlas literature review` listing" in rules
    assert "field in the `review` listing" in extract
    for path in ALL_FILES:
        text = _text(path)
        assert not re.search(r"NCBI_API_KEY\s*=\s*\S", text), path.name
        assert not re.search(r"api_key\s*[=:]\s*\S", text), path.name
