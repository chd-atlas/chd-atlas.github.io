# tests/unit/test_build_render.py
"""The shell every page shares, and the escaping every page depends on."""

import re

import pytest

from chd_atlas.build.render import (
    RESEARCH_USE_NOTICE,
    STYLESHEET,
    Link,
    Row,
    chip,
    data_table,
    definition_list,
    document,
    review_status_note,
)
from chd_atlas.models.review import ReviewStatus


def test_the_research_use_notice_is_the_owner_s_sentence_character_for_character() -> None:
    """The one string in this build that may not be reworded, shortened or split.

    Asserted against a literal written out here rather than against
    `RESEARCH_USE_NOTICE` itself, because a test that compares the constant to
    itself passes on every possible wording -- including a "shortened for the
    gene pages" variant, which is exactly the edit this forbids. The `<strong>`
    is part of what was specified and is compared with the rest of it.
    """
    assert RESEARCH_USE_NOTICE == (
        "This atlas is under active development. It is "
        "<strong>not a clinical decision-support tool</strong> and must not be used to make "
        "or guide a diagnostic, treatment or any other clinical decision."
    )


def test_every_page_carries_the_research_use_notice_before_its_own_content() -> None:
    """It reached 1 of the 25 pages the build writes, and it was not a gene page.

    Measured 2026-08-04 against a real build: `grep -l "not a clinical
    decision-support tool"` matched `index.html` alone. The 24 pages that carry a
    gene-level classification -- 23 gene pages and `genes/index.html` -- carried
    no research-use statement at all, so a clinician arriving at
    `genes/HGNC_4173.html` from a search engine met `GATA4`, an HGNC id and a
    green `definitive` chip and nothing else above the fold.

    `document` is where it lives now, which is what makes it structurally
    impossible to add a page kind without it: `landing.py` and `pages.py` both
    render through this function and neither mentions the notice.

    **Position is asserted, not merely presence.** A notice in the footer is a
    notice the reader meets after the classification it qualifies, and on a gene
    page the classification is the first thing in `<main>`. `body` here stands in
    for a gene page's rail, and the assertion is that the sentence precedes it in
    DOM order -- which is reading order, screen-reader order and, on a viewport
    under 46rem where `.layout` collapses to one column, painting order too.
    """
    page = document(title="GATA4 — CHD Atlas", root="../", body="<h1>GATA4</h1>")

    assert page.count(RESEARCH_USE_NOTICE) == 1
    assert f'<div class="notice" role="note"><p>{RESEARCH_USE_NOTICE}</p></div>' in page
    assert page.index(RESEARCH_USE_NOTICE) < page.index("<main>")
    assert page.index(RESEARCH_USE_NOTICE) < page.index("<h1>GATA4</h1>")


def test_a_document_is_self_contained_and_carries_no_external_reference() -> None:
    """No `<link>`, `<script src>` or `@import`: the page renders with no
    network request beyond the one that fetched it, so an offline reader and a
    reviewer unpacking `dist/` see the same site.

    Not a strict-CSP claim, which this docstring used to make. `default-src
    'self'` blocks an inline `<style>` and an inline `<script>` exactly as it
    blocks a remote one, so under such a policy the styling and the browse
    filter die. Self-contained and CSP-safe are different properties and only
    the first is asserted here.
    """
    page = document(title="CHD Atlas", root="", body="<p>hi</p>")

    assert page.startswith("<!doctype html>")
    assert "<link" not in page
    assert "<script src" not in page
    assert "@import" not in page
    assert "<style>" in page


def test_the_title_is_escaped() -> None:
    assert "<script>" not in document(title="<script>x</script>", root="", body="")


def test_a_nested_page_reaches_the_root_through_the_prefix() -> None:
    """A gene page lives in `genes/`, so its links back need `../`. Passed in
    rather than computed here: `render.py` knows nothing about where a caller
    put its page.
    """
    assert 'href="../index.html"' in document(title="TBX5", root="../", body="")
    assert 'href="index.html"' in document(title="CHD Atlas", root="", body="")


def test_every_cell_of_a_table_is_escaped_including_a_link() -> None:
    """`data_table` takes text and `Link`, never HTML.

    A curator-supplied or mirror-supplied string reaching a page as markup is
    the one way this build turns data into script, and `bundles.py` and
    `search.py` already treat both as data for the same reason. Taking
    pre-rendered HTML here would move that decision to every call site.
    """
    html = data_table(
        headers=["source", "report"],
        rows=[Row(cells=("<b>clingen</b>", Link(text="open <x>", href="https://e.org/?a=1&b=2")))],
    )

    assert "<b>" not in html
    assert "&lt;b&gt;clingen&lt;/b&gt;" in html
    assert "&lt;x&gt;" in html
    assert "https://e.org/?a=1&amp;b=2" in html


def test_a_quote_in_an_href_cannot_close_the_attribute_early() -> None:
    """The escape that stops a URL from becoming a second attribute.

    `html.escape` rewrites `"` only under `quote=True`, its default; passing
    `quote=False` anywhere in `render.py` would still escape every angle bracket
    — so `test_every_cell_of_a_table_is_escaped_including_a_link` above stays
    green — while letting any href close its own attribute and open an event
    handler. Asserted against the whole anchor, because a bare `"` is legal in
    the surrounding markup and only its position inside the value is the defect.
    """
    href = 'https://e.org/?a=" onmouseover="alert(1)'
    html = data_table(headers=["report"], rows=[Row(cells=(Link(text="open", href=href),))])

    assert '<a href="https://e.org/?a=&quot; onmouseover=&quot;alert(1)">open</a>' in html
    assert ' onmouseover="' not in html


def test_a_row_attribute_value_is_escaped() -> None:
    """`data-*` carries curated text too — lesion group names reach it verbatim."""
    html = data_table(
        headers=["gene"],
        rows=[Row(cells=("TBX5",), attributes=(("lesion", '<b> "septal"'),))],
    )

    assert "<b>" not in html
    assert 'data-lesion="&lt;b&gt; &quot;septal&quot;"' in html


def test_an_attribute_name_that_is_not_a_bare_word_is_refused() -> None:
    """Escaping cannot make an attribute *name* safe, so it is refused instead.

    This test replaces one whose docstring was measured false on 2026-08-04. It
    passed `'lesion" onclick="x'` and asserted `' onclick="' not in html`,
    reasoning that a `"` "closes the attribute and starts another one the row
    never declared". Both halves were wrong. The bytes did contain
    `onclick=&quot;`, so the assertion passed — while `html.parser` read the
    result as `[('data-lesion"', None), ('onclick', '"x="v"')]`, meaning a
    second attribute *was* declared. And the `"` is not what does it: `<`, `>`,
    `&`, `"` and `'` are all escaped, so the character that opens a second
    attribute is the **space**, which `html.escape` does not touch and which
    that payload smuggled in unnoticed.

    The payload below carries no quote and no angle bracket at all — nothing
    `html.escape` would alter — and without the guard it renders
    `<tr data-x onmouseover=alert(1) y="v">`, a live event handler.

    Every name the build passes is a literal in `pages.py`, so this is a guard
    on a bypassed gate: reaching it means a caller derived a name from data.
    """
    with pytest.raises(ValueError, match="escaping cannot"):
        data_table(
            headers=["gene"],
            rows=[Row(cells=("TBX5",), attributes=(("x onmouseover=alert(1) y", "v"),))],
        )


def test_a_review_status_without_its_note_is_refused_rather_than_rendered_bare() -> None:
    """A guard on a bypassed gate (owner's standing rule, 2026-10-27).

    `models/review.py` refuses a status with no note when the file loads, so
    reaching this means a record was built some other way. Rendering it would
    print "Pending review:" with nothing after it: a question and no
    statement of what it is, which is what the rule exists to prevent. With no
    status there is nothing to say.
    """
    with pytest.raises(ValueError, match="no review_note to render"):
        review_status_note("pmid22969434_cases", ReviewStatus.PENDING_REVIEW, None)
    assert review_status_note("cnchd", None, None) == ""


def test_a_table_with_no_rows_still_renders_its_header() -> None:
    """A gene with no variants and a table that was never written must not read
    the same to a reader. The header names the columns, so an empty body is an
    answer; a missing table is a defect that looks like one.
    """
    html = data_table(headers=["variant", "consequence"], rows=[])

    assert "<thead>" in html
    assert '<th scope="col">variant</th><th scope="col">consequence</th>' in html
    assert "<tbody></tbody>" in html


def test_a_chip_and_a_definition_list_escape_their_values() -> None:
    assert "<i>" not in chip(text="<i>definitive</i>", kind="definitive")
    assert "<i>" not in definition_list([("<i>k</i>", "<i>v</i>")])


def test_the_page_column_cannot_grow_past_the_viewport() -> None:
    """`.layout`'s content track is `minmax(0, 1fr)`, never a bare `1fr`.

    A bare `1fr` floors at its widest child's min-content, so one wide table
    stretched the whole column past the window: measured 2026-11-02 at a
    1,280px window, 47 of 119 gene pages scrolled sideways (RIT1's 1,411px
    wide), and 0 with this track. Wide tables keep their own `.scroll`.
    """
    tracks = re.findall(r"\.layout\s*\{[^}]*grid-template-columns:\s*([^;]+);", STYLESHEET)

    assert len(tracks) == 2, tracks  # the two-column rule and the narrow-screen one
    assert all(track.endswith("minmax(0, 1fr)") for track in tracks), tracks
