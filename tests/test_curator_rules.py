"""The curator rules file, held to the code it describes (CLAUDE.md 4.20, 4.25).

`curator/curator-rules.md` is plain Markdown read by the agents of
chd-atlas/curator. Nothing executes it, so a code that does not exist or an
option that was renamed would pass every other test and mislead the next agent.
Every name below is read out of the code at test time -- the typer app, the
validator source, the enums, the models -- rather than from a list kept here.

The `chd-curator` commands the file names belong to another repository; its
tests (`test_every_curator_command_and_option_the_files_name_exists`) read this
file from `dev` and check them. This file checks only `chd-atlas` commands.
"""

from __future__ import annotations

import re
import textwrap
from collections.abc import Iterator
from enum import StrEnum
from pathlib import Path
from typing import Any

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
RULES = CURATOR / "curator-rules.md"


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


def _fenced_yaml(text: str) -> list[Any]:
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


def test_curator_holds_only_the_rules_file() -> None:
    # F20: the plugin, its skills and its README moved to chd-atlas/curator.
    assert sorted(p.name for p in CURATOR.iterdir()) == ["curator-rules.md"]


def test_the_literature_commands_are_review_and_target() -> None:
    assert set(_literature_commands()) == {"review", "target"}


def test_every_atlas_command_and_option_the_rules_name_exists() -> None:
    commands = _literature_commands()
    every_option = set().union(*(_options(c) for c in commands.values()))
    text = _text(RULES)
    joined = re.sub(r"\\\n\s*", " ", text)
    named: set[str] = set()
    for line in joined.splitlines():
        if "chd-curator" in line:
            continue  # checked by chd-atlas/curator's own tests
        for form in re.finditer(r"chd-atlas literature (\S+)", line):
            name = form.group(1).strip("`.,;:)")
            assert name in commands, f"no command `literature {name}`"
            named.add(name)
            rest = line.split(form.group(0), 1)[1]
            for option in re.findall(r"(?<![\w-])(--[a-z][a-z-]*)", rest):
                assert option in _options(commands[name]), f"`{name}` {option}"
    assert "review" in named
    for name in re.findall(r"`literature ([a-z][a-z-]*)", text):
        assert name in commands, f"no command `literature {name}`"
    # A bare `--option` is an atlas option unless the line is a curator command.
    for line in joined.splitlines():
        if "chd-curator" in line:
            continue
        for option in re.findall(r"(?<![\w-])(--[a-z][a-z-]+)", line):
            assert option in every_option, option


def test_every_code_the_rules_name_is_raised_and_they_cite_every_literature_code() -> None:
    source = "".join(
        p.read_text(encoding="utf-8") for p in (REPO / "src" / "chd_atlas").rglob("*.py")
    )
    validator = (REPO / "src/chd_atlas/validate/gene_literature.py").read_text(encoding="utf-8")
    literature_codes = set(re.findall(r'"(LIT\d{3})"', validator))
    text = _text(RULES)
    for code in re.findall(r"\b([A-Z]{2,6}\d{3})\b", text):
        assert f'"{code}"' in source, f"{code} is raised nowhere"
        if code.startswith("LIT"):
            assert code in literature_codes, code
    assert set(re.findall(r"\b(LIT\d{3})\b", text)) == literature_codes


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
    for block in _fenced_yaml(_text(RULES)):
        for key, value in _walk(block):
            assert key in _model_fields(), f"unknown field {key}"
            if key in by_field:
                values = value if isinstance(value, list) else [value]
                assert set(values) <= by_field[key], f"{key}: {values}"


def test_backticked_names_are_real_values_fields_or_options() -> None:
    parameters = {
        param.name
        for command in _literature_commands().values()
        for param in command.params
        if param.name
    }
    known = set().union(*_enum_values().values()) | _model_fields() | parameters
    prose = re.sub(r"^[ ]*```.*?^[ ]*```", "", _text(RULES), flags=re.S | re.M)
    for token in re.findall(r"`([a-z]+(?:_[a-z]+)+)`", prose):
        assert token in known, f"`{token}` is not a value, field or option"


def test_the_rules_state_the_verified_rulings_literally() -> None:
    rules = _text(RULES)
    assert "agent:<skill>/<model>" in rules
    assert "leaves the approval\n  stale: the validator reports `LIT015`" in rules
    assert "reading this field in the `chd-atlas literature review` listing" in rules
    assert "`chd-curator search`" in rules
    assert "`chd-curator abstracts`" in rules
    assert "`chd-curator diff`" in rules
    assert "`chd-curator add-publication`" in rules
    assert not re.search(r"NCBI_API_KEY\s*=\s*\S", rules)
    assert not re.search(r"api_key\s*[=:]\s*\S", rules)


def test_the_readme_points_at_curator() -> None:
    readme = _text(REPO / "README.md")
    assert "https://github.com/chd-atlas/curator" in readme
    assert "chd-atlas literature review" in readme
