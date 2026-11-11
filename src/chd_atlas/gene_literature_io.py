# src/chd_atlas/gene_literature_io.py
"""Reading and writing the curator-edited literature files without losing a value
or a comment of what the curator wrote.

``curation/literature/<SYMBOL>.yaml`` is rewritten through ruamel's
round-trip loader with ONE canonical indent (mappings 2, sequences 4, dashes
at offset 2; ruling R20, no indent detection). Values and comments are kept.
A file already in that style changes only in the lines a command adds
(measured 2026-11-03: marking a finding verified on such a file changed one line, the
added ``verified_by``). A wrapped plain scalar is reflowed onto one line, and
a file in another indent style is re-indented, both once, on its first write
(measured: a summary hand-wrapped over two lines became one, and a file at
mapping indent 4 changed in 19 lines for one added line).

Of the writers that once lived here, ``chd-atlas/curator`` now owns the ones
that append searches and publications; this module keeps the reader and
``with_findings_updated`` that ``literature review`` uses for approvals.

`with_findings_updated` sets fields on named findings under the same rule: on a
file in the canonical style the lines it adds are the only change, and a blank
line or comment that followed the entry (even after a nested list or mapping)
moves to the new last line with it.

It returns **bytes** and writes nothing: the caller serialises, validates and
only then writes (`fs.write_bytes_atomically`), so a failure anywhere leaves the
file as it was.
"""

from __future__ import annotations

import copy
from datetime import date
from io import BytesIO
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ConfigDict, Field
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

from chd_atlas.duplicates import duplicates
from chd_atlas.models.gene_literature import Finding, GeneLiteratureFile


class LiteratureFileError(Exception):
    """A literature or publications file cannot be safely extended."""


def literature_path(root: Path, symbol: str) -> Path:
    return root / "curation" / "literature" / f"{symbol}.yaml"


def _round_trip_yaml() -> YAML:
    yaml = YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.preserve_quotes = True
    # Wide enough that a long recorded query is never folded, and an existing
    # line is never re-wrapped on the way through.
    yaml.width = 4096
    return yaml


def dump_bytes(data: Any) -> bytes:
    buffer = BytesIO()
    _round_trip_yaml().dump(data, buffer)
    return buffer.getvalue()


def load_round_trip(raw: bytes) -> CommentedMap:
    try:
        data = _round_trip_yaml().load(raw)
    except Exception as exc:  # ruamel raises several YAMLError subclasses
        raise LiteratureFileError(f"not valid YAML: {exc}") from exc
    if not isinstance(data, CommentedMap):
        raise LiteratureFileError("the top level is not a mapping")
    return data


def parse_literature(raw: bytes) -> GeneLiteratureFile:
    """The file as the model sees it, or `LiteratureFileError` naming why not."""
    try:
        data = YAML(typ="safe").load(raw)
        return GeneLiteratureFile.model_validate(data)
    except Exception as exc:  # YAMLError or pydantic.ValidationError
        raise LiteratureFileError(str(exc)) from exc


class _FindingsFile(BaseModel):
    """A verifier's output: `findings:` and nothing else (ruling R10)."""

    model_config = ConfigDict(extra="forbid")

    findings: list[Finding] = Field(default_factory=list)


def parse_findings_file(raw: bytes) -> list[Finding]:
    """The findings of a blind second pass, or `LiteratureFileError` naming why not.

    Each is held to the same per-record rules as a curated finding. Ids are
    local to the file and are not checked against the id registry; duplicates
    within the file are refused.
    """
    try:
        loaded = _FindingsFile.model_validate(YAML(typ="safe").load(raw))
    except Exception as exc:  # YAMLError or pydantic.ValidationError
        raise LiteratureFileError(str(exc)) from exc
    repeated = duplicates(f.id for f in loaded.findings)
    if repeated:
        raise LiteratureFileError(f"duplicate finding ids: {repeated}")
    return loaded.findings


def _trailing_slot(node: Any, key: Any) -> tuple[Any, Any, int] | None:
    """`(container, key, slot)` of the deepest comment trailing `node[key]`.

    ruamel hangs whatever follows a block (a `# note`, or a blank line before
    the next list item) on the last element inside it. At each level this
    checks the element's own comment slot -- slot 2 for a mapping key, slot 0
    for a sequence index -- and then walks into the last value if it is a
    non-empty collection. The deepest slot holding a comment wins, so a
    same-line comment on a collection's own key (`carriers:  # why`) does not
    shadow the blank line below it. Returns None when no level holds one.
    """
    found: tuple[Any, Any, int] | None = None
    while True:
        slot = 2 if isinstance(node, CommentedMap) else 0
        row = node.ca.items.get(key)
        if row and row[slot] is not None:
            found = (node, key, slot)
        child = node[key]
        if not isinstance(child, CommentedMap | CommentedSeq) or not len(child):
            return found
        node = child
        key = list(child)[-1] if isinstance(child, CommentedMap) else len(child) - 1


def _carry_trailing_comment(entry: CommentedMap, new_last: str) -> None:
    """Move what followed an entry's last element to its new last key.

    Appending a key after the old last one would otherwise leave a blank line
    or a note that separated two findings stranded between two fields of the
    same finding. A same-line comment (`key: v  # why`) stays where it was.
    Call it after the new key is inserted at the end; `new_last` is that key.
    """
    found = _trailing_slot(entry, list(entry)[-2])
    if found is None:
        return
    container, key, slot = found
    row = container.ca.items[key]
    token = row[slot]
    first, newline, rest = token.value.partition("\n")
    moved = copy.copy(token)
    if first.strip():  # a same-line comment stays; only what follows moves
        if not rest:
            return
        moved.value = "\n" + rest
        token.value = first + newline
    else:
        row[slot] = None
        if all(item is None for item in row):
            container.ca.items.pop(key)
    entry.ca.items[new_last] = [None, None, moved, None]


def with_findings_updated(existing: bytes, updates: dict[str, dict[str, Any]]) -> bytes:
    """The file's new content with fields set on the named findings.

    `updates` maps a finding id to `{field: value}`. A field already in the
    entry keeps its place; a new one goes right after the last field that
    precedes it in `Finding`'s own order, so the entry reads in model order.
    Everything else is round-tripped as `load_round_trip` keeps it. A date is
    copied per entry, never shared. The result
    is checked through `GeneLiteratureFile` before it is returned, so a change
    the model would refuse is never produced; an unknown id is refused too.
    """
    data = load_round_trip(existing)
    entries = {str(f.get("id")): f for f in data.get("findings") or []}
    missing = sorted(set(updates) - set(entries))
    if missing:
        raise LiteratureFileError(f"no such finding: {', '.join(missing)}")
    order = list(Finding.model_fields)
    for finding_id, fields in updates.items():
        entry = entries[finding_id]
        for name, value in fields.items():
            if isinstance(value, date):
                # One shared date object would be dumped once as a YAML anchor
                # and aliased (`*id001`) in every other finding it was set on.
                value = date(value.year, value.month, value.day)
            if name in entry:
                entry[name] = value
                continue
            before = [k for k in entry if k in order and order.index(k) < order.index(name)]
            position = max((list(entry).index(k) for k in before), default=-1) + 1
            last = list(entry)[-1] if position == len(entry) else None
            entry.insert(position, name, value)
            if last is not None:
                _carry_trailing_comment(entry, name)
    content = dump_bytes(data)
    parse_literature(content)
    return content
