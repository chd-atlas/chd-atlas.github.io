# tests/unit/_svg.py
"""Cut one `<g>` group out of drawn markup, nested groups and all.

Slicing at the first `</g>` after the opening tag stops inside the first
nested group, so any mark drawn after it is never seen and a count of marks
comes out low without failing. This walks to the group's own closing tag.
"""

from __future__ import annotations

import re

# An opening `<g ...>` (group 1 is "/" when it is self-closing) or a `</g>`.
# `\b` keeps tags such as `<glyph>` out.
_G_TAG = re.compile(r"<g\b[^>]*?(/?)>|</g>")


def g_part(markup: str, css: str) -> str:
    """`<g class="css">` and its content, up to but not including its own `</g>`.

    Raises `ValueError` when the group is missing or never closed.
    """
    start = markup.index(f'<g class="{css}">')
    depth = 0
    for tag in _G_TAG.finditer(markup, start):
        if tag.group(0) == "</g>":
            depth -= 1
            if depth == 0:
                return markup[start : tag.start()]
        elif not tag.group(1):
            depth += 1
    raise ValueError(f'<g class="{css}"> is never closed')
