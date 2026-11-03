# tests/unit/_stamp.py
"""Give a fixture finding a correct `approved_digest`.

The digest is computed with `finding_digest`, never typed: a hand-typed hash
would pass only until the fixture changed.
"""

from __future__ import annotations

import re
from typing import Any

from ruamel.yaml import YAML

from chd_atlas.models.gene_literature import Finding, finding_digest

_APPROVED_ON = re.compile(r"^( +)approved_on: .*\n", re.MULTILINE)


def stamped(raw: dict[str, Any]) -> dict[str, Any]:
    """`raw` with `approved_digest` set to the digest of the rest of it.

    A finding with no `approved_by` is returned without a digest.
    """
    without = {k: v for k, v in raw.items() if k != "approved_digest"}
    if without.get("approved_by") is None:
        return without
    probe = {**without, "approved_digest": "sha256:probe"}
    return {**without, "approved_digest": finding_digest(Finding.model_validate(probe))}


def stamped_yaml(text: str) -> str:
    """A literature file's text with an `approved_digest` line after each `approved_on`.

    The digest is that of the finding the line sits in, as the file parses. An
    unapproved finding is left as it is. Nothing else in the text moves.
    """
    findings = YAML(typ="safe").load(text).get("findings") or []
    digests = [
        finding_digest(Finding.model_validate({**f, "approved_digest": "sha256:probe"}))
        for f in findings
        if f.get("approved_by") is not None
    ]
    lines = _APPROVED_ON.findall(text)
    if len(lines) != len(digests):
        raise ValueError(f"{len(lines)} approved_on lines for {len(digests)} approved findings")
    parts = iter(digests)
    return _APPROVED_ON.sub(
        lambda m: f"{m.group(0)}{m.group(1)}approved_digest: {next(parts)}\n", text
    )
