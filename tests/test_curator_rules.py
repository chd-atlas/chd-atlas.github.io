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

import pytest
import typer
from ruamel.yaml import YAML

from chd_atlas import vocab
from chd_atlas.cli import app
from chd_atlas.models.assertion import InTextLocator, SupplementaryLocator
from chd_atlas.models.gene_literature import (
    Carriers,
    Finding,
    ReportedNames,
    Search,
    Synthesis,
    SynthesisSentence,
    TriageEntry,
    Unclear,
    Variant,
)
from chd_atlas.models.literature import Publication
from chd_atlas.models.review import ReviewStatus
from chd_atlas.validate.finding_variants import shape_issues
from chd_atlas.validate.gene_synthesis import _first_number

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
        | set(Synthesis.model_fields)
        | set(SynthesisSentence.model_fields)
        | set(Variant.model_fields)
        | set(ReportedNames.model_fields)
        | set(Unclear.model_fields)
        | {"with"}  # Variant.with_, read and written through its alias
        | {"findings", "triage", "searches", "gene", "synthesis"}
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


def test_the_literature_commands_are_the_six_the_rules_may_name() -> None:
    assert set(_literature_commands()) == {
        "review",
        "questions",
        "target",
        "frame",
        "digests",
        "apply-review",
    }


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
    validators = "".join(
        (REPO / "src/chd_atlas/validate" / name).read_text(encoding="utf-8")
        for name in ("gene_literature.py", "finding_variants.py")
    )
    literature_codes = set(re.findall(r'"(LIT\d{3})"', validators))
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
        "placement": enums["VariantPlacement"],
        "zygosity": enums["VariantZygosity"],
        "phase": enums["Phase"],
        "settled_by": enums["SettledBy"],
        "review_settled_by": enums["SettledBy"],
    }
    for block in _fenced_yaml(_text(RULES)):
        # `reason` is an ExclusionReason on a triage entry and free text on an
        # `unclear` entry; the second is read through the model instead.
        free = {
            Unclear.model_validate(entry).reason
            for key, value in _walk(block)
            if key == "unclear"
            for entry in value.values()
        }
        for key, value in _walk(block):
            assert key in _model_fields(), f"unknown field {key}"
            if key == "reason" and value in free:
                continue
            if key in by_field and not isinstance(value, dict):  # `unclear.zygosity` is an entry
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


def _variant_examples() -> list[list[Variant]]:
    blocks = [b["variants"] for b in _fenced_yaml(_text(RULES)) if "variants" in b]
    return [[Variant.model_validate(v) for v in block] for block in blocks]


def test_the_variant_examples_raise_only_the_before_normalise_messages() -> None:
    """Each `variants` example, set in a finding that admits it, raises what the rules say.

    Each `mane` entry (one in the first block, two in the `with` block) is as an
    agent writes it, before `normalise`: it must raise exactly the three LIT022
    messages the rules quote and nothing else (the LIT025 message needs the MANE
    frame, so `validate` over a copy of the repository measures it, not this test).
    Every other entry raises nothing. The first block holds one entry per
    placement state, so deleting a field a state requires from its example
    fails here.
    """
    blocks = _variant_examples()
    assert len(blocks) == 3
    assert [v.placement.value for v in blocks[0]] == [
        "mane", "protein", "structural", "outside_region", "unplaced",
    ]  # fmt: skip
    issues = []
    per_block = []
    for variants in blocks:
        finding = Finding.model_validate(
            {
                "id": "CHDA:LIT:0000001",
                "publication": "PMID:10000001",
                "design": "cohort_screen",
                "methods": ["targeted_sequencing"],
                "count_unit": "probands",
                "screened": 100,
                "carriers": {"de_novo": 4, "inherited": 4, "origin_unknown": 20},
                "variant_classes": ["snv_indel", "cnv"],
                "variants": [v.model_dump(by_alias=True, mode="json") for v in variants],
                "phenotypes": ["HP:0001631"],
                "lesion_groups": ["septal"],
                "syndromic": "isolated",
                "shares_probands_with": "unknown",
                "locators": [{"kind": "in_text", "element": "Table 1"}],
                "summary": "Illustrative.",
                "extracted_by": "agent:literature-extract/model",
            }
        )
        found = shape_issues(finding, {1})
        per_block.append([(c, n) for c, n, _ in found])
        issues += found
    rules = " ".join(_text(RULES).split())
    assert per_block == [
        [("LIT022", 1)] * 3,
        [("LIT022", 1)] * 3 + [("LIT022", 2)] * 3,
        [],
    ]  # one `mane` entry in the first block, two in the second
    for _, _, message in issues:
        assert f'"{message}"' in rules, message


def test_the_rules_state_the_variant_rulings_literally() -> None:
    rules = " ".join(_text(RULES).split())
    for phrase in (
        "One entry per distinct variant per genotype context.",
        "`not_stated` only when the paper is silent about it.",
        "Never choose between the paper's own names",
        "Every `pending_review` names `review_settled_by`, and `not_accessible` takes none.",
        "a group of more than two is never `trans`",
        "`phase` is set exactly when `with` is",
        "carries one identical `unclear.phase` entry on every member",
    ):
        assert phrase in rules, phrase


def test_the_rules_number_examples_hold_against_syn002() -> None:
    """X-4: the rules prescribed "a sibling sharing one parent", which SYN002 refuses.

    Each quoted example the no-numbers bullet says SYN002 refuses is refused,
    and each rewording it offers, and each name it says passes, passes: measured
    with the validator's own `_first_number`, against literals the rules must
    carry.
    """
    rules = " ".join(_text(RULES).split())
    refused = [
        "twofold", "3rd", "1990s", "mid-2000s", "3x", "½", "³",
        "half-sibling", "trisomy 21", "45,X", "46,XX", "Noonan syndrome 1",
        "one of the reports", "p.Gly35Argfs*75", "NM_181486.4", "NKX2–5",
    ]  # fmt: skip
    for term in refused:
        assert f'"{term}"' in rules, term
        assert _first_number(term) is not None, term
    rewordings = [
        "a sibling sharing only the mother", "only the father", "Down syndrome",
        "Turner syndrome", "an XX karyotype",
    ]  # fmt: skip
    for term in rewordings:
        assert f'"{term}"' in rules, term
        assert _first_number(term) is None, term
    for name in ("NKX2-5", "22q11.2", "HP:0001631", "p.Gly80Arg"):
        assert f" {name}" in rules, name
        assert _first_number(name) is None, name


def test_the_readme_points_at_curator() -> None:
    readme = _text(REPO / "README.md")
    assert "https://github.com/chd-atlas/curator" in readme
    flat = " ".join(readme.split())
    assert (
        "The searches it records are published as they are, but no finding is published "
        "until a named person approves it here, with `chd-atlas literature review` or with "
        "`chd-atlas literature apply-review` on the decisions file of a review page that "
        "person has finished."
    ) in flat
    assert "Nothing it writes is published" not in flat


def _entry(**fields: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "id": 1,
        "as_reported": "x",
        "locator": {"kind": "in_text", "element": "Table 1"},
        "kind": "snv_indel",
        "placement": "mane",
        "zygosity": "heterozygous",
        "carriers": {"de_novo": 0, "inherited": 0, "origin_unknown": 1},
    }
    return {**base, **fields}


def _verdict(*entries: dict[str, Any]) -> list[str]:
    finding = Finding.model_validate(
        {
            "id": "CHDA:LIT:0000001",
            "publication": "PMID:10000001",
            "design": "cohort_screen",
            "methods": ["targeted_sequencing"],
            "count_unit": "probands",
            "screened": 100,
            "carriers": {"de_novo": 0, "inherited": 0, "origin_unknown": 50},
            "variant_classes": ["snv_indel", "cnv"],
            "variants": list(entries),
            "phenotypes": ["HP:0001631"],
            "lesion_groups": ["septal"],
            "syndromic": "isolated",
            "shares_probands_with": "unknown",
            "locators": [{"kind": "in_text", "element": "Table 1"}],
            "summary": "Illustrative.",
            "extracted_by": "agent:literature-extract/model",
        }
    )
    return [m for _, _, m in shape_issues(finding, {1})]


_FILLED = {"hgvs_c": "c", "grch38": "g", "normalisation": 1, "hgvs_p": "p"}
_MARK = {
    "review_status": "pending_review",
    "review_note": "n",
    "review_settled_by": "atlas_model",
}


@pytest.mark.parametrize(
    ("sentence", "entry", "needle"),
    [
        ("`mane` holds `hgvs_c`, `grch38` and `normalisation`", _entry(), "needs hgvs_c"),
        (
            "`protein` holds `hgvs_p` and leaves `hgvs_c`, `grch38` and `normalisation` null",
            _entry(placement="protein", hgvs_p="p", hgvs_c="c"),
            "must leave hgvs_c null",
        ),
        (
            "`protein` holds `hgvs_p` and leaves `hgvs_c`, `grch38` and `normalisation` null",
            _entry(placement="protein", hgvs_p="p", normalisation=1),
            "must leave normalisation null",
        ),
        (
            "`structural` and `outside_region` leave `hgvs_c`, `hgvs_p` and `grch38` null",
            _entry(placement="structural", kind="cnv", grch38="g"),
            "must leave grch38 null",
        ),
        (
            "`structural` and `outside_region` leave `hgvs_c`, `hgvs_p` and `grch38` null",
            _entry(placement="outside_region", hgvs_p="p"),
            "must leave hgvs_p null",
        ),
        (
            "`unplaced` holds `review_status` and leaves `hgvs_c` and `grch38` null",
            _entry(placement="unplaced", grch38="g", **_MARK),
            "must leave grch38 null",
        ),
        (
            "`unplaced` holds `review_status` and leaves `hgvs_c` and `grch38` null",
            _entry(placement="unplaced"),
            "needs review_status",
        ),
        (
            "Links are mutual",
            _entry(id=1, with_=[2], phase="trans"),
            "does not name 1 back",
        ),
        (
            "a group of more than two is never `trans`",
            [_entry(id=i, **{"with": [j for j in (1, 2, 3) if j != i]}, phase="trans")
             for i in (1, 2, 3)],
            "more than two variants pairwise in trans",
        ),
    ],
)  # fmt: skip
def test_each_stated_hold_null_and_link_rule_is_the_validators(
    sentence: str, entry: Any, needle: str
) -> None:
    assert sentence in " ".join(_text(RULES).split()), sentence
    entries = entry if isinstance(entry, list) else [entry]
    entries = [{("with" if k == "with_" else k): v for k, v in e.items()} for e in entries]
    if not isinstance(entry, list) and entry.get("with_"):
        entries.append(_entry(id=2, placement="mane", **_FILLED))
    assert any(needle in m for m in _verdict(*entries)), _verdict(*entries)


def test_a_group_of_two_may_be_trans_and_a_mane_entry_with_every_name_is_clean() -> None:
    pair = [
        _entry(id=1, **_FILLED, **{"with": [2]}, phase="trans"),
        _entry(id=2, **{**_FILLED, "grch38": "h"}, **{"with": [1]}, phase="trans"),
    ]
    assert _verdict(*pair) == []


@pytest.mark.parametrize(
    ("entry", "settled", "ok"),
    [
        ({"placement": "unplaced"}, "atlas_model", True),  # rsID alone, or a `g.` name alone
        ({"placement": "unplaced"}, "expert", True),  # `normalise` could not place it
        ({"placement": "unplaced"}, None, False),
    ],
)
def test_pending_review_needs_a_settler_and_not_accessible_takes_none(
    entry: dict[str, Any], settled: str | None, ok: bool
) -> None:
    mark = {"review_status": "pending_review", "review_note": "n"}
    if settled:
        mark["review_settled_by"] = settled
    try:
        Variant.model_validate(_entry(**entry, **mark))
    except ValueError:
        assert not ok
    else:
        assert ok
    with pytest.raises(ValueError, match="not_accessible takes no review_settled_by"):
        Variant.model_validate(
            _entry(
                review_status="not_accessible",
                review_note="n",
                review_settled_by="full_text",
            )
        )


def test_the_settler_sentence_and_the_with_rule_are_stated_literally() -> None:
    rules = " ".join(_text(RULES).split())
    for phrase in (
        "An `unplaced` variant named by an rsID or a `g.` name alone that the paper does "
        "not place is `atlas_model`; a `protein` variant whose "
        "residue is refused is `authors`, or `atlas_model` when its `p.` names no residue; "
        "a `mane` variant `normalise` cannot place becomes `unplaced` with `expert`.",
        "a `protein` variant's `hgvs_p` is the agent's, and is checked at the first `validate`",
        "needs them on every design but `common_variant_association`",
        "a `protein` variant whose residue validate then refuses",
        "whatever names it: a `c.`, an rsID, a `g.` position or a description",
        "Position beats naming, and a position is never looked up in dbSNP, a genome "
        "browser or the agent's own knowledge",
        'A coding-DNA name printed without its `c.` prefix ("456delC", "998G>T") is still '
        "a `c.` name: write it `mane` with `reported.c` exactly as printed; `normalise` "
        "decides, and an agent never marks it `unplaced` itself",
        "for a variant named by an rsID or a `g.` name alone that the paper does not place",
    ):
        assert phrase in rules, phrase


def test_an_outside_region_variant_named_by_a_g_position_alone_validates() -> None:
    """R-F2: position beats naming, so a `g.`-only `reported` needs no mark."""
    entry = _entry(
        placement="outside_region",
        as_reported="g.12345678A>G",
        reported={"rsid": None, "c": None},
    )
    assert _verdict(entry) == []
    assert _verdict(_entry(placement="outside_region", reported={"rsid": "rs1"})) == []
