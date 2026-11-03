# src/chd_atlas/build/render.py
"""The document shell and the HTML primitives every page is built from.

This module knows nothing about genes, assertions or validity. It knows how to
turn text into a page. It holds every *primitive* that escapes -- so escaping is
decided here rather than at each call site -- but it is **not** the only place
that produces markup: `pages.py` and `landing.py` both assemble their own
elements around these primitives, and this paragraph claimed otherwise until
review on #15 measured it. The claim mattered, because it is what made an
unescaped `<option value>` in `pages.py` look impossible. What is true is
narrower and is the rule to hold: **a value derived from curated or mirrored
data must reach markup only through a primitive here.**

Self-contained by construction, like every artifact this build writes: no
`<link>`, no `<script src>`, no `@import`, no `url(` in the stylesheet, no
external font, no `<img>` -- only an inline `<style>` and, on the browse page,
an inline `<script>`, both literals in this file. A page renders with no network
request beyond the one that fetched it, from a local checkout and from the
`dist/` directory a reviewer unpacks.

That is not the same as surviving a strict Content-Security-Policy, which this
paragraph also used to claim. Measured on a real build: `default-src 'self'`
blocks an inline `<style>` and an inline `<script>` exactly as it blocks a
remote one, so under such a policy the browse filter and all styling die. GitHub
Pages sets no CSP, so nothing breaks on deploy; a consumer self-hosting these
pages behind one needs a hash or nonce allowance.

Deterministic for the reason every other builder here is: nothing reads a clock,
and every value interpolated comes from the caller.

`data_table`, `chip` and `definition_list` take **text**, never HTML. A curated
gene symbol or a mirrored disease label reaching a page as markup is the one way
this build turns data into script -- `bundles.py` and `search.py` already treat
curator- and mirror-supplied strings as data for the same reason -- and taking
pre-rendered HTML would move that decision to every caller. `Link` exists so a
cell can carry a URL without any caller needing to write an `<a>` tag.

`document` renders `RESEARCH_USE_NOTICE` itself, so the statement is on the page
before any caller's body is. It used to live in `landing.py`, which meant it
reached `index.html` and nothing else: measured 2026-08-04 against a real build,
`grep -l "not a clinical decision-support tool" ` matched 1 of the 25 HTML files
the build writes, and the 24 that carry a gene-level classification -- 23 gene
pages and the browse page -- carried no research-use statement at all. A
clinician arriving at `genes/HGNC_4173.html` from a search engine met `GATA4`, an
HGNC id and a green `definitive` chip, first in DOM order and filling the first
screen under 46rem. Putting it here is what makes it structurally impossible to
add a page kind without it.
"""

from __future__ import annotations

import html
import re
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Final, assert_never

from chd_atlas.models.review import ReviewStatus

# What `data_table` accepts as a `data-*` attribute name. Deliberately narrower
# than HTML permits: every name this build passes is a literal, so the only way
# to reach a rejection is for a caller to have derived one from data.
_ATTRIBUTE_NAME: Final = re.compile(r"[a-z][a-z0-9-]*")

# Moved here from `landing.py`, which now renders through `document` below.
# One stylesheet rather than one per page kind: three `<style>` blocks are three
# things that drift, and the landing page's palette is the one this extends.
STYLESHEET: Final = """
  :root {
    color-scheme: light dark;
    --fg: #0f172a; --bg: #f8fafc; --muted: #64748b; --border: #e2e8f0;
    --notice-bg: #fffbeb; --notice-border: #f59e0b; --link: #0284c7;
    --chip-bg: rgba(100, 116, 139, 0.08); --definitive: #16a34a; --warn: #dc2626;
    /* The validity ladder, one hue at four intensities -- the same idiom as
       --evidence below and for the same reason. ClinGen's four supportive rungs
       are ordinal, so a ramp is what they mean; a second hue would invent an
       axis. `ungraded` is drawn as absence (outline, muted) rather than as the
       bottom of the ramp: no panel graded the gene, which is not a weaker grade
       but a different question, and every gene admitted on submitter
       agreement is in it. */
    --grade-4: #16a34a; --grade-4-fg: #ffffff;
    --grade-3: #0d9488; --grade-3-fg: #ffffff;
    --grade-2: #0284c7; --grade-2-fg: #ffffff;
    --grade-1: #64748b; --grade-1-fg: #ffffff;
    /* The evidence hue, deliberately NOT the green of --definitive and not
       a red/green pair. `significant` is not `true` on this site, and a
       green tick beside a nominal-only result would undo every caveat
       around it. One hue at three intensities encodes how much support
       there is; absence is drawn as absence rather than as a third colour. */
    --evidence: #0284c7; --evidence-soft: #e0f2fe; --on-evidence: #ffffff;
    --surface: #ffffff;
    /* The organ grid's four fixed bins, one hue at four steps (an ordered
       quantity). The dark steps are their own values, not the light ones
       flipped: on a dark page the low bin must sit near the background and
       the high bin must be the brightest. */
    --heat-1: #bae6fd; --heat-2: #38bdf8; --heat-3: #0284c7; --heat-4: #075985;
  }
  @media (prefers-color-scheme: dark) {
    :root {
      --fg: #f1f5f9; --bg: #0b0f17; --muted: #94a3b8; --border: #1e293b;
      --notice-bg: #451a03; --notice-border: #f59e0b; --link: #38bdf8;
      --definitive: #22c55e; --warn: #ef4444;
      --grade-4: #22c55e; --grade-4-fg: #052e16;
      --grade-3: #14b8a6; --grade-3-fg: #042f2e;
      --grade-2: #38bdf8; --grade-2-fg: #082f49;
      --grade-1: #94a3b8; --grade-1-fg: #0f172a;
      --evidence: #38bdf8; --evidence-soft: #082f49; --on-evidence: #0b0f17;
      --surface: #131b28;
      --heat-1: #0c4a6e; --heat-2: #0369a1; --heat-3: #0ea5e9; --heat-4: #7dd3fc;
    }
  }
  * { box-sizing: border-box; }
  body {
    margin: 0; padding: 1.5rem 1.5rem 5rem; background: var(--bg); color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica,
      Arial, sans-serif;
    line-height: 1.5;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
  }
  main { max-width: 76rem; margin: 0 auto; }
  nav {
    max-width: 76rem; margin: 0 auto 1.5rem; font-size: 0.92rem;
    display: flex; align-items: center; gap: 0.5rem;
    padding: 0.75rem 1rem; border: 1px solid var(--border);
    border-radius: 10px; background: var(--surface);
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
  }
  nav a {
    text-decoration: none; font-weight: 600; color: var(--muted);
    padding: 0.35rem 0.65rem; border-radius: 6px; transition: all 0.15s ease;
  }
  nav a:first-child {
    font-size: 1.05rem; font-weight: 750; color: var(--fg); margin-right: auto;
    letter-spacing: -0.02em;
  }
  nav a:hover { color: var(--link); background: var(--bg); }
  h1 { margin-bottom: 0.25rem; font-size: 2.2rem; letter-spacing: -0.02em; }
  .tagline { color: var(--muted); margin-top: 0; font-size: 1.05rem; }
  h2 {
    margin-top: 2.5rem; border-bottom: 1px solid var(--border); padding-bottom: 0.4rem;
    font-size: 1.35rem; letter-spacing: -0.01em;
  }
  .notice {
    display: none;
  }
  .notice p { margin: 0; }
  .layout .scope-rule { display: none; }
  .scope-rule { color: var(--muted); font-size: 0.9rem; }
  /* On methods.html the rule is the section's content, not a gloss beside a
     gene: body colour and size there, and only there. The gene pages' and the
     gene index's short form keeps the muted style above. */
  h2#scope + .scope-rule { color: var(--fg); font-size: 1rem; }
  /* The burden section's four paragraph kinds. `.notice-inline` and
     `.provenance` were emitted by `pages.py` with no rule at all, so the single
     most important safety message this site can print -- "these two studies are
     not independent ... must not be combined" -- rendered as ordinary body text
     while the research-use notice got a bordered amber box. Found by review
     2026-08-05. */
  .notice-inline {
    background: var(--surface); border: 1px solid var(--border);
    border-left: 4px solid var(--evidence); border-radius: 6px;
    padding: 0.75rem 1rem; margin: 1rem 0; font-size: 0.86rem; line-height: 1.5;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.02);
  }
  .method, .provenance {
    color: var(--muted); font-size: 0.85rem; margin: 0.35rem 0; line-height: 1.4;
  }
  .disclosure {
    color: var(--muted); font-size: 0.82rem; margin: 0.5rem 0; font-style: italic;
    background: var(--surface); padding: 0.5rem 0.75rem; border-radius: 6px;
    border: 1px solid var(--border);
  }
  .footnotes { color: var(--muted); font-size: 0.85rem; margin: 0.5rem 0 0; }
  /* The concordance strip and matrix. One evidence hue at three intensities
     rather than a red/green pair: "significant" is not "true" on this site, and
     a green tick beside a nominal-only result would undo every caveat around
     it. `not tested` is drawn as an absence -- a dash, a dashed border -- and
     never as a weak positive, because it is the state that would otherwise read
     as "this dataset looked and found nothing". */
  .strip { display: inline-flex; gap: 0.3rem; align-items: center; }
  .dot {
    width: 15px; height: 15px; border-radius: 50%; box-sizing: border-box;
    border: 1.5px solid var(--evidence); display: inline-grid; place-items: center;
  }
  .dot.full { background: var(--evidence); }
  .dot.half {
    background: linear-gradient(90deg, var(--evidence) 50%, transparent 50%);
  }
  .dot.none { border-color: var(--muted); }
  .dot.untested {
    border: 0; color: var(--muted); font-size: 15px; line-height: 1;
  }
  .strip-tally {
    font-size: 0.78rem; color: var(--muted); margin-left: 0.35rem;
    font-variant-numeric: tabular-nums; white-space: nowrap;
  }
  .strip-legend { color: var(--muted); font-size: 0.85rem; margin: 0.5rem 0 0; }
  .glance { margin: 0.25rem 0 0.5rem; }
  /* One legend entry -- a glyph and the sentence it captions -- as a single
     unbreakable unit. Without this the line can break immediately after the
     glyph, orphaning "not tested by that dataset" from any mark at ordinary
     desktop widths. `&nbsp;` does not fix it: `.dot` is `display: inline-grid`,
     an atomic inline box, and the break opportunity after it survives. */
  .legend-item { display: inline-block; }
  /* The matrix key: a swatch per state, drawn from the same rules as the cells
     so the key cannot drift from what it describes. */
  .cell-key {
    display: inline-block; width: 13px; height: 13px; border-radius: 3px;
    vertical-align: -2px; border: 1px solid transparent;
  }
  .cell-key.corrected { background: var(--evidence); }
  .cell-key.nominal { background: var(--evidence-soft); border-color: var(--evidence); }
  .cell-key.no-enrichment { background: var(--surface); border-color: var(--border); }
  .cell-key.not-tested { border: 1px dashed var(--border); }

  /* The matrix is the widest thing on the page. Without this the *body*
     scrolls sideways on a narrow screen, which moves the whole layout rather
     than the table. `.scroll` was emitted by `_evidence_matrix` with no rule
     anywhere -- found by review before the deploy. */
  .scroll { overflow-x: auto; margin: 0.75rem 0; }
  table.matrix {
    border-collapse: separate; border-spacing: 6px; margin: 0.5rem 0;
    background: var(--surface); padding: 0.75rem 1rem; border-radius: 10px;
    border: 1px solid var(--border); box-shadow: 0 1px 3px rgba(0, 0, 0, 0.03);
  }
  table.matrix th {
    font-size: 0.74rem; font-weight: 650; text-transform: uppercase;
    letter-spacing: 0.04em; color: var(--muted);
  }
  table.matrix th[scope="row"] {
    text-align: right; white-space: nowrap; font-size: 0.82rem;
    font-weight: 600; color: var(--fg); padding-right: 0.75rem;
  }
  .cell {
    display: flex; flex-direction: column; align-items: center; justify-content: center;
    min-width: 6.5rem; min-height: 2.8rem;
    padding: 0.4rem 0.65rem; border-radius: 6px; text-align: center;
    font-variant-numeric: tabular-nums; font-size: 0.84rem; line-height: 1.25;
    border: 1px solid transparent; white-space: nowrap; transition: transform 0.15s ease;
  }
  .cell:hover { transform: scale(1.02); }
  .cell.corrected {
    background: var(--evidence); color: var(--on-evidence); font-weight: 650;
    box-shadow: 0 1px 3px rgba(2, 132, 199, 0.25);
  }
  .cell.nominal {
    background: var(--evidence-soft); border-color: var(--evidence);
    color: var(--fg); font-weight: 600;
  }
  .cell.no-enrichment {
    background: var(--surface); color: var(--muted); border: 1px solid var(--border);
  }
  .cell.not-tested {
    border: 1px dashed var(--border); color: var(--muted); font-size: 0.75rem;
  }
  .cell .sub { font-size: 0.72rem; opacity: 0.9; margin-top: 0.1rem; }

  /* The front page's hero band, census cards and call to action. Every rule
     below is written against a token declared above, so the dark theme is
     covered by construction rather than by a second media query that someone
     has to remember to update. */
  .hero {
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 12px; padding: 1.6rem 2rem; margin: 1.5rem 0 0;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
  }
  .figures { display: flex; flex-wrap: wrap; gap: 1.5rem 3.5rem; margin: 0 0 1.25rem; }
  .figure { display: flex; flex-direction: column; }
  .figure-value {
    font-size: 2.4rem; line-height: 1; font-weight: 700;
    font-variant-numeric: tabular-nums; color: var(--evidence);
  }
  .figure-label {
    font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.04em;
    font-weight: 600; color: var(--muted); margin-top: 0.35rem;
  }
  /* A qualifier under a hero figure (`render.REGISTRY_QUALIFIER`): read with
     the number, so smaller and muted rather than in the value's size and
     weight. */
  .figure-note {
    font-size: 0.78rem; font-weight: 400; line-height: 1.3; color: var(--muted);
  }
  .hero-note { margin: 0 0 1.25rem; font-size: 0.92rem; color: var(--muted); line-height: 1.6; }
  .cta {
    display: inline-block; background: var(--evidence); color: var(--on-evidence);
    text-decoration: none; border-radius: 6px; padding: 0.55rem 1.15rem; font-weight: 550;
    box-shadow: 0 1px 2px rgba(0,0,0,0.08); transition: background 0.15s ease;
  }
  .cta:hover { background: #0369a1; text-decoration: none; }
  .cta:focus-visible { outline: 2px solid var(--fg); outline-offset: 2px; }
  /* A `<dl>` because it is one -- label and value -- restyled from the base
     two-column list into cards. `dl.stats` is (0,1,1) and the base `dl` rule
     below is (0,0,1), so the two do not fight over `grid-template-columns`.
     Each pair is wrapped in a `.stat` div, which HTML5 permits inside a `<dl>`
     and which is what lets one label and its value be one card. */
  dl.stats {
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(13rem, 1fr));
    gap: 0.75rem; margin: 1.25rem 0;
  }
  .stat {
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 8px; padding: 0.85rem 1rem;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.02);
  }
  .stat dt { font-size: 0.82rem; line-height: 1.3; color: var(--muted); }
  .stat dd {
    font-size: 1.4rem; font-weight: 700; font-variant-numeric: tabular-nums;
    margin-top: 0.25rem; color: var(--fg);
  }

  .reading-notes { margin: 1rem 0; }
  .reading-notes summary { cursor: pointer; color: var(--muted); }
  .reading-notes summary:focus-visible, .cohort-notes summary:focus-visible,
  .forest summary:focus-visible, .chart-notes summary:focus-visible {
    outline: 2px solid var(--evidence); outline-offset: 2px;
  }
  .cohort-notes { margin: 1rem 0; font-size: 0.9rem; }
  .cohort-notes summary { cursor: pointer; color: var(--muted); }
  .cohort-notes li { margin: 0.5rem 0; }
  /* `methods.html#cohorts`: each name above its description rather than beside
     it, since the base `dl` grid sizes its label column to the widest label.
     (0,1,1) against the base rule's (0,0,1). */
  dl.cohort-list { display: block; }
  .cohort-list dd { margin: 0.15rem 0 0.9rem; }
  dl { display: grid; grid-template-columns: auto 1fr; gap: 0.35rem 1rem; margin: 0; }
  dt { color: var(--muted); }
  dd { margin: 0; overflow-wrap: break-word; }
  a { color: var(--link); }
  code { background: var(--chip-bg); padding: 0.1rem 0.3rem; border-radius: 3px; }
  ul.links { list-style: none; padding: 0; margin: 0; }
  ul.links li { margin: 0.35rem 0; }
  footer { margin-top: 3rem; color: var(--muted); font-size: 0.9rem; }
  table { border-collapse: separate; border-spacing: 0; width: 100%; font-size: 0.88rem; }
  th, td { text-align: left; padding: 0.65rem 0.85rem; border-bottom: 1px solid var(--border); }
  th {
    color: var(--muted); font-weight: 600; font-size: 0.74rem; text-transform: uppercase;
    letter-spacing: 0.04em; background: var(--surface);
  }
  tr:last-child td { border-bottom: none; }
  tr:hover td { background-color: rgba(2, 132, 199, 0.03); }
  .chip {
    display: inline-flex; align-items: center; background: var(--chip-bg); border-radius: 9999px;
    padding: 0.15rem 0.55rem; font-size: 0.76rem; font-weight: 600; text-transform: uppercase;
    letter-spacing: 0.04em; margin: 0 0.25rem 0.25rem 0; line-height: 1.2;
  }
  .chip-definitive { background: var(--grade-4); color: var(--grade-4-fg); }
  .chip-strong { background: var(--grade-3); color: var(--grade-3-fg); }
  .chip-moderate { background: var(--grade-2); color: var(--grade-2-fg); }
  .chip-limited { background: var(--grade-1); color: var(--grade-1-fg); }
  .chip-ungraded {
    background: transparent; color: var(--muted); border: 1px dashed var(--border);
  }
  .chip-warn { background: var(--warn); color: #fff; }
  .grade-key { margin: 1rem 0; font-size: 0.9rem; }
  /* A row gap rather than a margin on `dt`: the margin pushed each chip below
     the first line of its own gloss, since `dd` had none. */
  .grade-key dl { row-gap: 0.6rem; }
  /* `minmax(0, 1fr)`: a bare `1fr` grows to its widest table, which pushed
     47 of 119 gene pages past a 1,280px window. */
  .layout { display: grid; grid-template-columns: 18rem minmax(0, 1fr); gap: 2.5rem; }
  @media (max-width: 52rem) { .layout { grid-template-columns: minmax(0, 1fr); } }
  .rail {
    font-size: 0.88rem; background: var(--surface); padding: 1.5rem;
    border-radius: 12px; border: 1px solid var(--border); height: fit-content;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.05);
  }
  .rail h1 { font-size: 2rem; margin: 0 0 0.35rem 0; letter-spacing: -0.02em; }
  .rail p { margin: 0.35rem 0; }
  .rail dl { margin: 0.75rem 0; }
  .rail .xrefs { margin: 0.75rem 0; display: flex; flex-wrap: wrap; gap: 0.35rem; }
  .rail .xrefs a {
    display: inline-block; padding: 0.2rem 0.5rem; background: var(--chip-bg);
    border-radius: 4px; border: 1px solid var(--border); text-decoration: none;
    font-size: 0.78rem; font-weight: 550; color: var(--fg);
  }
  .rail .card-variants {
    margin-top: 1.25rem; padding-top: 1rem; border-top: 1px solid var(--border);
  }
  .rail .card-variants ul { list-style: none; padding: 0; margin: 0.5rem 0; }
  .rail .card-variants li {
    background: var(--bg); padding: 0.5rem 0.65rem; border-radius: 6px;
    border-left: 3px solid var(--evidence); margin-bottom: 0.5rem;
    font-size: 0.8rem; line-height: 1.4;
  }
  .rail .card-variants li a { font-weight: 550; }
  .rail .card-gloss, .rail .card-scope, .rail .card-curation, .rail .card-uncorrected {
    font-size: 0.78rem; color: var(--muted); line-height: 1.4;
    background: var(--bg); padding: 0.5rem 0.65rem; border-radius: 6px;
    border: 1px solid var(--border); margin: 0.5rem 0;
  }
  .rail a[href*="methods.html"] {
    display: inline-block; font-size: 0.8rem; font-weight: 550; margin-top: 0.5rem;
  }

  /* Main column cards and tables on gene page */
  .layout table {
    border: 1px solid var(--border); border-radius: 8px; overflow: hidden;
    background: var(--surface); box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03);
    margin: 0.75rem 0 1.5rem;
  }
  .layout table th {
    background: #f8fafc; font-size: 0.72rem; font-weight: 650;
    text-transform: uppercase; letter-spacing: 0.04em; padding: 0.65rem 0.85rem;
    border-bottom: 2px solid var(--border); color: #475569;
  }
  @media (prefers-color-scheme: dark) {
    .layout table th { background: #1e293b; color: #94a3b8; }
  }
  .layout table td {
    padding: 0.65rem 0.85rem; font-size: 0.86rem;
    border-bottom: 1px solid var(--border);
  }
  .layout table tr:last-child td { border-bottom: none; }
  .layout table tr:hover td { background-color: rgba(2, 132, 199, 0.03); }

  /* The ClinGen expert panel's rows are the gold standard; GenCC rows stay plain. */
  .layout tr[data-authority="clingen"] td { background: #f0fdf4; }
  .layout tr[data-authority="clingen"] td:first-child { box-shadow: inset 3px 0 0 #16a34a; }
  #validity-table td:nth-child(5) { white-space: nowrap; }
  #validity-table td:first-child { min-width: 11rem; }
  #validity-table + p { font-size: 0.8rem; color: var(--muted); margin: 0.5rem 0 0; }
  .gold-standard {
    display: inline-block; margin-left: 0.4rem; padding: 0.1rem 0.45rem; border-radius: 9999px;
    font-size: 0.68rem; font-weight: 650; text-transform: uppercase; letter-spacing: 0.04em;
    background: #dcfce7; color: #15803d; border: 1px solid #bbf7d0; white-space: nowrap;
  }

  .layout h2 + p {
    background: var(--surface); border: 1px solid var(--border); border-radius: 8px;
    padding: 0.65rem 0.95rem 0.65rem 2.4rem; font-size: 0.84rem; color: var(--muted);
    position: relative; margin: 0.75rem 0 1rem; line-height: 1.5;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.02);
  }
  .layout h2 + p::before {
    content: "ⓘ"; position: absolute; left: 0.75rem; top: 0.65rem;
    display: inline-flex; align-items: center; justify-content: center;
    width: 1.25rem; height: 1.25rem; border-radius: 50%; background: #e0f2fe; color: #0284c7;
    font-weight: 700; font-size: 0.85rem;
  }

  /* Evidence assertion box styling */
  .assertion-block {
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 8px; padding: 1.2rem; margin: 1rem 0;
    box-shadow: 0 1px 2px rgba(0, 0, 0, 0.02);
  }
  .assertion-block dl { margin: 0.5rem 0; }
  .filters {
    display: flex; flex-wrap: wrap; gap: 0.6rem; margin: 1.25rem 0;
    background: var(--surface); padding: 0.85rem 1rem; border: 1px solid var(--border);
    border-radius: 8px; box-shadow: 0 1px 2px rgba(0, 0, 0, 0.03); align-items: center;
  }
  .filters input, .filters select {
    padding: 0.4rem 0.65rem; border: 1px solid var(--border); border-radius: 6px;
    background: var(--bg); color: var(--fg); font-size: 0.85rem;
  }
  .filters input[type="search"] { min-width: 16rem; }
  .filters input:focus, .filters select:focus {
    outline: 2px solid var(--evidence); outline-offset: 1px;
  }

  /* Sortable table headers and gene registry table styling */
  th.sortable {
    cursor: pointer; user-select: none; position: relative; padding-right: 1.4rem;
    transition: background 0.15s ease, color 0.15s ease;
  }
  th.sortable:hover { color: var(--fg); background: var(--border); }
  th.sortable::after {
    content: " ⇅"; opacity: 0.35; font-size: 0.72rem; position: absolute;
    right: 0.4rem; top: 50%; transform: translateY(-50%);
  }
  th.sortable.sort-asc::after {
    content: " ▲"; opacity: 1; color: var(--evidence);
  }
  th.sortable.sort-desc::after {
    content: " ▼"; opacity: 1; color: var(--evidence);
  }
  .table-scroll {
    overflow-x: auto; margin: 1.25rem 0; border: 1px solid var(--border);
    border-radius: 10px; background: var(--surface); box-shadow: 0 1px 3px rgba(0, 0, 0, 0.04);
  }
  .table-scroll table { margin: 0; }
  .table-scroll th { background: var(--bg); position: sticky; top: 0; z-index: 1; }
  .gene-desc { font-size: 0.78rem; color: var(--muted); margin-top: 0.15rem; line-height: 1.25; }
  .expression-cell {
    font-variant-numeric: tabular-nums; white-space: nowrap; font-weight: 600; color: var(--fg);
  }
  .expression-pct {
    display: inline-block; font-size: 0.72rem; font-weight: 550; color: #0284c7;
    background: #e0f2fe; padding: 0.1rem 0.4rem; border-radius: 4px; margin-left: 0.25rem;
  }
  .concordance-key {
    display: flex; flex-wrap: wrap; gap: 1rem; font-size: 0.82rem;
    color: var(--muted); margin: 0.5rem 0 1rem; align-items: center;
  }

  .chip-curated {
    display: inline-block; background: #e0f2fe; color: #0369a1; border: 1px solid #bae6fd;
    font-size: 0.68rem; font-weight: 650; text-transform: uppercase; letter-spacing: 0.04em;
    padding: 0.12rem 0.4rem; border-radius: 4px; vertical-align: 1px; margin-left: 0.25rem;
  }
  .source-tag { font-size: 0.72rem; color: var(--muted); margin-top: 0.2rem; font-weight: 500; }
  .gene-id {
    font-size: 0.76rem; color: var(--muted); background: var(--bg); border: 1px solid var(--border);
    padding: 0.1rem 0.35rem; border-radius: 4px;
    font-family: ui-monospace, SFMono-Regular, monospace;
    font-weight: 500;
  }

  .matrix-card {
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 12px; padding: 1.5rem; margin: 1rem 0;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.03);
  }
  .study-card {
    background: var(--surface); border: 1px solid var(--border);
    border-radius: 12px; padding: 1.5rem; margin: 1.5rem 0;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.03);
  }
  .study-figures {
    margin: 1rem 0; background: var(--bg); border: 1px solid var(--border);
    border-radius: 8px; padding: 0.75rem 1rem;
  }
  .study-figures summary {
    font-weight: 600; font-size: 0.85rem; color: var(--fg); cursor: pointer;
  }

  h3[id^="study-"] {
    margin-top: 0.5rem; margin-bottom: 0.4rem; font-size: 1.25rem; font-weight: 700;
  }
  h3[id^="study-"] a { text-decoration: none; color: var(--fg); }
  h3[id^="study-"] a:hover { color: var(--link); text-decoration: underline; }
  .study-analyses {
    font-size: 0.78rem; font-weight: 600; color: var(--evidence);
    background: var(--evidence-soft); padding: 0.2rem 0.6rem; border-radius: 4px;
    display: inline-block; margin: 0.2rem 0 0.75rem;
  }

  /* The finding callout that leads each study card: which consequence classes
     survive the study's own correction, visible at a glance instead of buried
     in repeated forest-plot rows. Emitted by `_study_finding_summary`. */
  .finding-card {
    border-radius: 10px; padding: 0.9rem 1.05rem; margin: 0.75rem 0 0.5rem;
    border: 1px solid var(--border); border-left-width: 4px;
    background: var(--bg); font-size: 0.9rem; line-height: 1.45;
  }
  .finding-badge {
    font-weight: 700; font-size: 0.82rem; text-transform: uppercase;
    letter-spacing: 0.03em; margin-bottom: 0.3rem;
  }
  .finding-summary { color: var(--fg); }
  .finding-missense { border-left-color: #16a34a; background: #f0fdf4; }
  .finding-missense .finding-badge { color: #15803d; }
  .finding-lof { border-left-color: #0891b2; background: #ecfeff; }
  .finding-lof .finding-badge { color: #0e7490; }
  .finding-both { border-left-color: #7c3aed; background: #f5f3ff; }
  .finding-both .finding-badge { color: #6d28d9; }
  .finding-union { border-left-color: #0284c7; background: #e0f2fe; }
  .finding-union .finding-badge { color: #0369a1; }
  .finding-nominal { border-left-color: var(--notice-border); background: var(--notice-bg); }
  .finding-nominal .finding-badge { color: var(--muted); }

  /* Table cells styling */
  #gene-table { width: 100%; table-layout: auto; }
  #gene-table th, #gene-table td { vertical-align: top; }
  #gene-table th:nth-child(1) { width: 12%; min-width: 7.5rem; }
  #gene-table th:nth-child(2) { width: 20%; min-width: 11rem; }
  #gene-table th:nth-child(3) { width: 16%; min-width: 10rem; }
  #gene-table th:nth-child(4) { width: 24%; min-width: 13rem; }
  #gene-table th:nth-child(5) { width: 14%; min-width: 9.5rem; }
  #gene-table th:nth-child(6) { width: 14%; min-width: 8.5rem; }

  #gene-table td:nth-child(1) a {
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    font-size: 0.8rem; color: var(--muted); text-decoration: none;
  }
  #gene-table td:nth-child(1) a:hover { color: var(--link); text-decoration: underline; }
  #gene-table td:nth-child(2) a {
    font-weight: 700; font-size: 0.95rem; text-decoration: none; color: var(--link);
  }
  #gene-table td:nth-child(2) a:hover { text-decoration: underline; }
  #gene-table td:nth-child(4) {
    overflow-wrap: break-word; word-break: normal; line-height: 1.4;
  }
  .cohort-tests {
    font-size: 0.75rem; color: var(--muted); margin-top: 0.25rem; font-weight: 500;
  }
  /* Charts (build/charts.py). Every colour resolves through a custom property,
     so the dark palette above governs charts with no second definition to
     drift. `charts.py` emits no colour of its own -- a hex baked into a shape
     would render a light chart on a dark page, and it is the class names below
     that make that constraint enforceable in one place. */
  /* The cell-type block (`cell_type_chart`): the same hues as the charts above,
     one rule per mark, light and dark through the variables. */
  .cell-types { margin: 1rem 0 1.5rem; }
  .ct-caption, .ct-tau { color: var(--muted); font-size: 0.9rem; margin: 0.35rem 0; }
  .ct-scroll { overflow-x: auto; }
  .ct-svg { display: block; }
  .ct-svg text { fill: var(--fg); font-size: 12px; }
  .ct-svg .ct-small, .ct-svg .ct-group, .ct-svg .ct-muted { fill: var(--muted); font-size: 10.5px; }
  .ct-svg .ct-group { letter-spacing: 0.06em; text-transform: uppercase; }
  .ct-svg .ct-mid { text-anchor: middle; }
  .ct-svg .ct-end { text-anchor: end; }
  .ct-svg .ct-num { font-variant-numeric: tabular-nums; }
  .ct-grid { stroke: var(--border); stroke-dasharray: 2 3; }
  .ct-floor { stroke: var(--notice-border); stroke-width: 1.2; stroke-dasharray: 4 3; }
  .ct-below { fill: var(--notice-bg); }
  .ct-whisker { stroke: var(--evidence); stroke-width: 2; }
  /* No fill-opacity here: the region grid steps each dot's by decade in an
     attribute, and a CSS property would override it on every dot. */
  .ct-dot { fill: var(--evidence); }
  .ct-ring { fill: none; stroke: var(--muted); stroke-width: 1.2; }
  .chart { width: 100%; height: auto; max-width: 40rem; overflow: visible; }
  .chart-line { fill: none; stroke: var(--link); stroke-width: 2; }
  .chart-hit { fill: none; stroke: none; pointer-events: all; }
  .chart-point { fill: var(--link); }
  .chart-point-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }
  .lit-mark { fill: var(--link); }
  .lit-mark-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }
  .chart-axis { stroke: var(--border); }
  .chart-band { fill: var(--notice-bg); }
  .chart-strip { fill: var(--muted); fill-opacity: 0.55; }
  .chart-absent { fill: var(--border); }
  .chart-control { stroke: var(--muted); fill: none; }
  .chart-control-point { fill: var(--muted); }
  .chart-nointerval { stroke: var(--border); stroke-width: 4; stroke-dasharray: 2 3; }
  .chart-arrow { fill: var(--link); }
  .chart-arrow-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }
  .chart-control-arrow { fill: var(--muted); }
  .chart-union { fill: var(--muted); font-size: 0.55rem; }
  /* SVG `<text>` takes no `color`, only `fill`, so a label with neither
     renders black on a dark page. Every `<text>` the charts emit carries this
     class; the `<tspan>` rule above is the other one that sets a text fill. */
  .chart-label { fill: var(--fg); font-size: 0.6rem; }
  .chart-tick-label { fill: var(--muted); font-size: 0.5rem; }
  .chart-tick-label.cardiac { fill: var(--fg); font-weight: 700; }
  .chart-tick-label.unsampled { fill: var(--muted); fill-opacity: 0.5; }
  .chart-axis-title { fill: var(--muted); font-size: 0.55rem; }
  .chart-grid { stroke: var(--border); stroke-dasharray: 2 3; fill: none; }
  /* One marker a forest plot might want is deliberately absent from the list
     above; `build/charts.py` names it and says why, and says it there rather
     than here on purpose. **This stylesheet is inlined verbatim into every
     page, so a sentence in this comment is a sentence in the published
     bytes.** This block stated the reason in full for one revision and thereby
     reintroduced a phrase that
     `test_the_pooling_notice_is_the_matrix_caption_and_no_longer_conditional`
     asserts is absent from a one-study page -- and naming the class here would
     put the very token a page-wide check greps for onto every page. Reason
     and class name both live in the module that emits the geometry. */
  /* The bulk figure's own scroller: the svg keeps a floor width, and the figure
     scrolls inside this box instead of the page scrolling sideways. Measured in
     headless Chrome in a 390 px iframe (it will not go below 500 px itself):
     at 560 px the stage labels rendered at 8.96 px and the smallest text at
     8.80 px; at 640 px they render at 10.24 px and the smallest at 10.06 px. */
  .chart-scroll { overflow-x: auto; margin: 0.5rem 0; }
  .chart-scroll .chart { min-width: 640px; max-width: none; }
  /* The organ grid: a fixed bin per shade (--heat-N), then the two unplaced
     states drawn as absence rather than as the bottom of the ramp. */
  .heat-1 { fill: var(--heat-1); }
  .heat-2 { fill: var(--heat-2); }
  .heat-3 { fill: var(--heat-3); }
  .heat-4 { fill: var(--heat-4); }
  .heat-floor { fill: var(--surface); stroke: var(--border); }
  .heat-gap { fill: var(--surface); stroke: var(--muted); stroke-dasharray: 2 2; }
  .cell-dot { fill: var(--muted); }
  .heat-key .legend-item { margin-right: 0.9rem; }
  .heat-swatch {
    display: inline-block; width: 0.8em; height: 0.8em; border-radius: 2px;
    vertical-align: -0.08em; border: 1px solid var(--border);
  }
  .heat-swatch.heat-1 { background: var(--heat-1); }
  .heat-swatch.heat-2 { background: var(--heat-2); }
  .heat-swatch.heat-3 { background: var(--heat-3); }
  .heat-swatch.heat-4 { background: var(--heat-4); }
  .heat-swatch.heat-floor { background: var(--surface); position: relative; }
  .heat-swatch.heat-floor::after {
    content: ""; position: absolute; left: 50%; top: 50%; width: 3px; height: 3px;
    margin: -1.5px 0 0 -1.5px; border-radius: 50%; background: var(--muted);
  }
  .heat-swatch.heat-gap { background: var(--surface); border: 1px dashed var(--muted); }
  .heat-swatch.heat-blank { background: transparent; }
  .forest {
    margin: 1rem 0; border: 1px solid var(--border); border-radius: 8px;
    background: var(--surface); padding: 1.25rem;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.03);
  }
  .forest figcaption { font-size: 0.82rem; color: var(--muted); margin-top: 0.5rem; }
  .forest figcaption p { margin: 0.25rem 0; }
  /* The forest fills its card instead of stopping at 40rem, and its marks scale with it. */
  .forest .chart { max-width: none; }
  svg.chart .chart-tick-label { font-size: 0.56rem; }
  .forest .chart-point, .forest .chart-point-open { r: 4.5; }
  .forest .chart-line, .forest .chart-control { stroke-width: 2.6; }
  .forest .chart-label { font-size: 0.64rem; }
  .forest .chart-union { font-size: 0.56rem; }
  .forest .chart-axis { stroke: var(--muted); stroke-opacity: 0.45; }
  .forest summary, .chart-notes summary { cursor: pointer; color: var(--muted); font-weight: 600; }
  .denovo {
    margin: 1rem 0; border: 1px solid var(--border); border-radius: 8px;
    background: var(--surface); padding: 1.25rem;
    box-shadow: 0 1px 3px rgba(0, 0, 0, 0.03);
  }
  .denovo .table-scroll { margin: 0; }
  .denovo table { font-size: 0.88rem; }
  .denovo th, .denovo td { padding: 0.5rem 0.75rem; }
  .denovo td:first-child { font-weight: 600; }
  .denovo figcaption { font-size: 0.82rem; color: var(--muted); margin-top: 0.5rem; }
  .denovo figcaption p { margin: 0.25rem 0; }
  .denovo-chip {
    display: inline-block; padding: 0.12rem 0.55rem; border-radius: 999px;
    font-size: 0.74rem; font-weight: 650; text-transform: uppercase; letter-spacing: 0.02em;
    white-space: nowrap;
  }
  .denovo-chip.denovo-sig { background: #dcfce7; color: #15803d; border: 1px solid #bbf7d0; }
  .denovo-chip.denovo-nominal { background: #fef3c7; color: #b45309; border: 1px solid #fde68a; }
  .denovo-chip.denovo-none {
    background: var(--surface); color: var(--muted); border: 1px solid var(--border);
  }
"""

# The one piece of client-side behaviour on the site, and a literal so it does
# not affect determinism. It only *hides* rows that the build already rendered,
# so the browse page is complete with JavaScript disabled -- which is why the
# table is server-rendered rather than fetched (D29).
FILTER_SCRIPT: Final = """
(function () {
  var form = document.getElementById('filters');
  var table = document.getElementById('gene-table');
  if (!table) return;
  var rows = Array.prototype.slice.call(
    table.querySelectorAll('tbody tr')
  );
  var count = document.getElementById('shown');

  function apply() {
    var q = form && form.querySelector('[name=q]');
    var text = (q ? q.value || '' : '').toLowerCase();
    var facets = form ? Array.prototype.slice.call(form.querySelectorAll('select')) : [];
    var shown = 0;
    rows.forEach(function (row) {
      var ok = (row.getAttribute('data-search') || '').indexOf(text) !== -1;
      facets.forEach(function (select) {
        var want = select.value;
        if (ok && want) {
          var attr = row.getAttribute('data-' + select.name) || '';
          ok = attr.split(' ').indexOf(want) !== -1;
        }
      });
      row.hidden = !ok;
      if (ok) shown += 1;
    });
    if (count) count.textContent = String(shown);
  }

  if (form) {
    form.addEventListener('input', apply);
    form.addEventListener('change', apply);
  }

  var headers = Array.prototype.slice.call(table.querySelectorAll('thead th'));
  var currentSortCol = -1;
  var currentSortAsc = true;

  var rankMap = {
    'definitive': 6,
    'strong': 5,
    'moderate': 4,
    'limited': 3,
    'disputed': 2,
    'refuted': 1,
    'expert_curated': 2,
    'submitter_curated': 1
  };

  function parseNumber(s) {
    if (!s || s === '—' || s.indexOf('<') === 0) return -1;
    var m = s.match(/([0-9.]+)/);
    return m ? parseFloat(m[1]) : -1;
  }

  function getCellValue(row, idx) {
    var td = row.children[idx];
    if (!td) return '';
    var custom = td.getAttribute('data-sort-val');
    if (custom !== null) return custom;
    return td.textContent.trim();
  }

  headers.forEach(function (th, colIdx) {
    th.classList.add('sortable');
    th.setAttribute('title', 'Click to sort');
    th.addEventListener('click', function () {
      if (currentSortCol === colIdx) {
        currentSortAsc = !currentSortAsc;
      } else {
        currentSortCol = colIdx;
        currentSortAsc = true;
      }

      headers.forEach(function (h, i) {
        h.classList.remove('sort-asc', 'sort-desc');
        if (i === colIdx) {
          h.classList.add(currentSortAsc ? 'sort-asc' : 'sort-desc');
        }
      });

      rows.sort(function (a, b) {
        var va = getCellValue(a, colIdx);
        var vb = getCellValue(b, colIdx);

        var na = parseNumber(va);
        var nb = parseNumber(vb);
        var isNumA = (na !== -1);
        var isNumB = (nb !== -1);

        var diff = 0;
        if (isNumA && isNumB) {
          diff = na - nb;
        } else if (isNumA && !isNumB) {
          diff = 1;
        } else if (!isNumA && isNumB) {
          diff = -1;
        } else {
          var rA = rankMap[va.toLowerCase()] !== undefined ? rankMap[va.toLowerCase()] : -1;
          var rB = rankMap[vb.toLowerCase()] !== undefined ? rankMap[vb.toLowerCase()] : -1;
          if (rA !== -1 || rB !== -1) {
            diff = rA - rB;
          } else {
            diff = va.localeCompare(vb, undefined, { numeric: true, sensitivity: 'base' });
          }
        }
        return currentSortAsc ? diff : -diff;
      });

      var tbody = table.querySelector('tbody');
      rows.forEach(function (row) {
        tbody.appendChild(row);
      });
      apply();
    });
  });

  apply();
})();
"""


# The project owner's own words, and the one string in this build that must
# reach a reader character for character. Kept as a single line with single
# spaces so a check for it is a plain substring check on the page's bytes rather
# than one that has to normalise a source wrap away first; `<strong>` is markup
# from this module, in the same class as `FILTER_SCRIPT` above, never a value
# from a caller.
#
# Do not reword, shorten or split this sentence. It is the site's only statement
# about what the atlas may and may not be used for, and every page now carries
# it, so a change here changes 25 pages at once.
RESEARCH_USE_NOTICE: Final = (
    "This atlas is under active development. It is "
    "<strong>not a clinical decision-support tool</strong> and must not be used to make "
    "or guide a diagnostic, treatment or any other clinical decision."
)


# What each of the four `FamilyState`s means, keyed by the state's own name and
# ordered as `concordance` ranks them: strongest support first, absence last.
#
# **One definition, two legends.** The browse page's dot strip and the gene
# page's evidence matrix both caption the same four states (the front page
# carried a third key until the doorway page, 2026-10-27). They draw different
# markup -- `.dot` glyphs, `.cell-key` swatches -- but must describe the states
# in the same words, and this project's single most repeated defect is
# duplicated prose going stale in one copy while staying internally consistent
# (CLAUDE.md section 4.24).
#
# **Keyed by state, never positional, and that is load-bearing.** The first
# version of this constant was a 4-tuple zipped against a separate glyph tuple
# at each call site, which made the caption-to-glyph pairing a relation between
# four independently reorderable sequences. Measured 2026-08-06 by mutation:
# swapping two entries survives all 787 tests, and the build then publishes a
# gene page whose key says the `not-tested` swatch means "tested, no enrichment
# detected" while that same page's cells say "this dataset did not test this
# gene". That is exactly the `not_tested`-reads-as-`no_enrichment` collapse this
# project forbids, on a green build. `zip(strict=True)` does not help: it checks
# length, not order.
#
# Keying by state makes the mispairing unrepresentable rather than merely
# detectable, which is the same reason `pages._STATE_TITLE` and
# `pages._STATE_GLYPH` are dicts. Reordering this mapping now changes only the
# order states are listed in; a caption cannot leave its state.
#
# The apostrophe is `&#x27;` rather than a literal, because these strings are
# interpolated into HTML by more than one call site and one of them may one day
# put one inside an attribute. `render` is where escaping decisions live.
#
# `corrected` read "survives that study's own correction" until review
# (2026-10-27): a glyph is a cohort family's, and since owner decision
# 2026-10-27 (I3) one family holds two studies, so "that study" had no single
# referent. A family is `corrected` when a row of one of its studies survives
# that row's own study's correction (`concordance.family_state`).
EVIDENCE_STATE_LABELS: Final[dict[str, str]] = {
    "corrected": "enriched, and a result survives its own study&#x27;s correction",
    "nominal": "enriched nominally, or no correction published",
    "no_enrichment": "tested, no enrichment detected",
    "not_tested": "not tested by that dataset",
}

# The glyph each state draws, per legend kind, keyed by the same state names.
# A state missing from either mapping raises `KeyError` in `evidence_legend`
# rather than rendering a legend one row short.
_DOT_GLYPHS: Final[dict[str, tuple[str, str]]] = {
    "corrected": ("dot full", ""),
    "nominal": ("dot half", ""),
    "no_enrichment": ("dot none", ""),
    "not_tested": ("dot untested", "&ndash;"),
}
_SWATCH_GLYPHS: Final[dict[str, tuple[str, str]]] = {
    "corrected": ("cell-key corrected", ""),
    "nominal": ("cell-key nominal", ""),
    "no_enrichment": ("cell-key no-enrichment", ""),
    "not_tested": ("cell-key not-tested", ""),
}

# The sentence that stops an empty matrix or a hollow strip reading as a verdict.
#
# Byte-identical in `_MATRIX_LEGEND` and `_STRIP_LEGEND` before this constant
# existed, which is the duplication the paragraph above says to avoid. It is the
# most load-bearing sentence the concordance layer publishes: KDM6A is ClinGen
# definitive, causes Kabuki syndrome, and shows nothing in either dataset that
# tested it.
EVIDENCE_POWER_CAVEAT: Final = (
    "<strong>No enrichment at these cohort sizes is not evidence against a gene</strong> "
    "&mdash; rare-variant tests routinely detect nothing for genes with overwhelming family and "
    "functional evidence."
)


# What an independent-datasets count can see, said wherever a page prints one:
# the front page's figure (owner decision 2026-10-26, R12) and each gene card's
# "Tested in N independent dataset(s)" line (owner answer 2026-10-27). One
# definition, so the two cannot come to say it differently.
# `concordance.cohort_families` groups studies through the collections their
# rows name and the overlaps `curation/cohorts.yaml` declares, and through
# nothing else: two studies whose case sets share people the registry does not
# record -- PMID:34324492's S2 Table lists case sources no row names
# (`methods.html#cohorts` counts them) -- count as two. Measured 2026-10-27: 2,
# of three studies, two of them joined by the one overlap the registry declares
# (DDD and DECIPHER) and by no shared collection.
REGISTRY_QUALIFIER: Final = "(as far as the cohort registry records)"


def review_status_label(status: ReviewStatus) -> str:
    """The words a page prints for a review status, before its note.

    One place, so the label after a description and the methods page's
    sentence saying what each label marks (`methods._review_key`) cannot come
    to name a status differently. Exhaustive: a status added to the vocabulary
    and not here fails `mypy --strict` rather than rendering as nothing.
    """
    match status:
        case ReviewStatus.PENDING_REVIEW:
            return "Pending review"
        case ReviewStatus.NOT_ACCESSIBLE:
            return "Source not accessible"
        case _:
            assert_never(status)


def review_status_note(record: str, status: ReviewStatus | None, note: str | None) -> str:
    """A curated record's open question, for after its description; `""` where none.

    The owner's standing rule (2026-10-27, `models/review.py`): a record
    carrying a question this atlas leaves open says so where it is read. One
    function, so the methods page's cohort list and a gene page's cohort notes
    cannot come to word it differently. The label says which kind of question
    (`review_status_label`), the note what is unresolved; the note is
    curator's free text and is escaped. A status arriving without its note
    means the model's refusal was bypassed, and that raises rather than
    printing a label with nothing after it.
    """
    if status is None:
        return ""
    if note is None:
        raise ValueError(f"{record}: review_status {status.value!r} has no review_note to render")
    return f" <strong>{review_status_label(status)}:</strong> {html.escape(note)}"


def evidence_legend(*, swatches: bool) -> str:
    """The four-state key, as one `<p>`. `swatches` picks square keys over dots.

    Every item is wrapped in a `.legend-item`, which the stylesheet makes
    `inline-block` so a line break cannot fall between a glyph and its own
    caption. `&nbsp;` does not achieve this and was measured not to: `.dot` is
    `display: inline-grid`, an atomic inline box, and a browser takes a break
    opportunity after it regardless of the whitespace that follows. Without the
    wrapper the front page orphaned "not tested by that dataset" from any glyph
    at ordinary desktop widths.
    """
    glyphs = _SWATCH_GLYPHS if swatches else _DOT_GLYPHS
    items = " &nbsp; ".join(
        f'<span class="legend-item"><span class="{glyphs[state][0]}">{glyphs[state][1]}</span>'
        f" {label}</span>"
        for state, label in EVIDENCE_STATE_LABELS.items()
    )
    return f'<p class="strip-legend">{items}</p>'


# The chip class each mirrored grade renders in, keyed by the grade's own value.
#
# **Keyed, never zipped against an ordered ladder.** `EVIDENCE_STATE_LABELS`
# above was a positional zip until a mutant swapping two of its entries survived
# all 787 tests and published a page that contradicted itself; this is the same
# shape of mapping and gets the same treatment.
#
# The three contested rungs map to `warn` and are here even though
# `validity.published_genes` admits no gene headlined by `disputed` or `refuted`
# -- the floor excludes them, and `_clingen_contests` vetoes admission on GenCC
# agreement where a panel contests. That makes this a guard on a bypassed gate,
# the same idiom as `encode_json`'s `allow_nan=False`: if the gate ever changes,
# a contested gene renders in the caution colour rather than raising `KeyError`
# mid-build or, worse, borrowing the green of `definitive`.
# `no_known_association` is *not* on a bypassed path -- it sits at rank 0 but
# outside `CONTESTED`, so a gene ClinGen recorded it for can still be admitted on
# GenCC agreement -- which is exactly why it must not fall through to a default.
GRADE_CHIP_KIND: Final[dict[str, str]] = {
    "definitive": "definitive",
    "strong": "strong",
    "moderate": "moderate",
    "limited": "limited",
    "disputed": "warn",
    "refuted": "warn",
    "no_known_association": "warn",
}

# What each rung means, in the words of the authority that defined it.
#
# **These gloss ClinGen's Gene-Disease Validity SOP; the atlas authors none of
# them** (D12). Before 2026-08-06 the site needed no key: every published gene
# was `definitive` and the only chip was green. The widened gate publishes four
# rungs, and a reader who meets `limited` on a clinical genetics resource with no
# statement of what a panel meant by it will supply their own meaning.
#
# `ungraded` is the state, not a rung: no expert panel has graded the gene for a
# disease in scope, and it is published on submitter agreement instead. Drawn as
# an outline rather than as the bottom of the ramp, because "nobody graded it"
# and "graded weakest" are different facts and the second is not implied by the
# first.
#
# Its gloss said the gene was ungraded "for a disease an external authority
# treats as congenital heart disease, or in which Orphanet records a structural
# heart malformation in at least 30% of cases", and that two or more submitters
# "independently assert it", until 2026-10-26. The first reads as a property
# of a disease, and a disease is in scope only once a curator lists it
# (`pages._CURATION_BY_AGREEMENT` records the gene it was false for); the
# second is weaker than the gate, which counts a submitter only at `Limited` or
# better (`validity.agreeing_records`) -- a `Supportive` submission asserts the
# gene and counts for nothing.
GRADE_GLOSS: Final[dict[str, str]] = {
    "definitive": (
        "repeatedly demonstrated in research and clinical settings and upheld over time, "
        "with no convincing evidence to the contrary"
    ),
    "strong": (
        "independently demonstrated by several studies, without the years of replication "
        "definitive requires"
    ),
    "moderate": "several unrelated probands and some experimental support",
    "limited": (
        "few probands, or evidence that is suggestive but not compelling. "
        "<strong>Limited is not a weak yes</strong> &mdash; it is a panel saying the case is "
        "not yet made"
    ),
    "no_known_association": (
        "a panel looked and found no evidence of a relationship to this disease"
    ),
    "ungraded": (
        "no ClinGen expert panel has graded this gene for any of the diseases in scope here. "
        "It is published because two or more Gene Curation Coalition submitters classify it "
        "Limited or better, and the gene page names them"
    ),
}

_GRADE_KEY_ORDER: Final = (
    "definitive",
    "strong",
    "moderate",
    "limited",
    "no_known_association",
    "ungraded",
)


def grade_key() -> str:
    """The validity-chip key, rendered open, for `methods.html#scope`.

    It was `grade_legend()` until 2026-10-26: a collapsed `<details>` rendered on
    the browse page and on every gene page, from one constant, for the reason
    `EVIDENCE_POWER_CAVEAT` is one constant: a key that exists on the page a
    reader arrives at and not on the page they land on is a key that is absent
    exactly when it is needed. Measured 2026-08-05 in the concordance review,
    where the browse page carried the evidence legend and the gene page did not.

    **Moved to `methods.html#scope` on 2026-10-26 (D50).** Its dispute sentence
    said "the record is in the table below", which was false on the gene index
    (no validity table) and is false here; it now points at the gene's validity
    table. "This page is not a complete account" became "a gene page" for the
    same reason. Both page kinds still reach it: each carries the scope short
    form, whose link names the key, and a gene page's card glosses its own
    headline grade in the key's words (`GRADE_GLOSS`). Rendered open, because on
    the page that exists to explain it the key is not a digression to fold.

    **It must render after `methods.SCOPE_RULE`.** The dispute sentence's
    condition read "Where that disease is one an external authority treats as
    congenital heart disease", which named one path into scope of two: Costello
    and cardiofaciocutaneous syndrome, the diseases the six disputed genes are
    disputed for, are in scope only on Orphanet's phenotype frequencies. It now
    reads "in scope here, by either path above" (2026-10-26), which points at
    the rule rather than restating it -- so it is true only below the rule, and
    `test_the_key_follows_the_rule_its_dispute_sentence_points_back_to` pins
    that order.

    `disputed` and `refuted` are deliberately not listed. They are in
    `GRADE_CHIP_KIND` as a guard, but no gene headlined by either is published,
    so the key says that instead of glossing a state a reader cannot meet.

    **What that sentence may claim is bounded by scope, and it once claimed more
    than that.** It read "a gene an expert panel disputed or refuted is not
    published here at all" until 2026-08-06. In scope that is exactly true --
    measured over the 92 published genes, 0 carry an in-scope ClinGen record in
    `CONTESTED`. Unqualified it is false: 5 published genes carry a ClinGen
    `Disputed` record for a disease outside this atlas's remit (ANKRD1 and MYH6
    for hypertrophic cardiomyopathy, MNS1 for primary ciliary dyskinesia, GJA1
    for nonsyndromic hearing loss, SMAD1 for pulmonary arterial hypertension),
    and `gene_validity` filters those rows out, so the dispute appears nowhere on
    the page. ANKRD1's page contained the word "disputed" exactly once -- in the
    sentence denying it. Hence the second half, which tells the reader the table
    is scope-filtered rather than complete.

    **And then the in-scope half went false too, on 2026-10-25.** The sentence
    read "No gene here was disputed or refuted by an expert panel for a disease an
    external authority treats as congenital heart disease". The Orphanet scope
    path brought Costello and cardiofaciocutaneous syndrome into scope, and
    ClinGen's RASopathy panel disputes six published genes for them -- BRAF,
    KRAS, MAP2K1, PTPN11, RAF1 and SOS1, each Definitive for the syndrome it
    is admitted on: Noonan syndrome or Noonan syndrome with multiple
    lentigines, or cardiofaciocutaneous syndrome itself for BRAF and MAP2K1,
    which the same panel disputes only for Costello syndrome.
    PTPN11's page carried the denial beside its own `conflicting evidence` chip
    and two `Disputed` rows. What stays true on every page is narrower and is
    what the sentence now says: no gene is *headlined* by a contested grade (the
    headline is the admitting record's, at `Limited` or better), an in-scope
    dispute is listed in the table and sets `has_conflicting_evidence` (every
    published gene has a supportive record, so the flag's other half always
    holds), and an out-of-scope one is not shown.

    The qualifier deliberately avoids the phrase "this atlas's CHD scope", which
    `test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas`
    asserts is absent from both page kinds, and `test_build_methods.py` from
    `methods.html`: it told 24 pages that this atlas decides what counts as
    congenital heart disease, and it does not.
    """
    rows = "".join(
        f'<dt><span class="chip chip-{GRADE_CHIP_KIND.get(grade, "ungraded")}">'
        f"{html.escape(grade.replace('_', ' ') if grade != 'ungraded' else 'not classified')}"
        f"</span></dt><dd>{GRADE_GLOSS[grade]}</dd>"
        for grade in _GRADE_KEY_ORDER
    )
    return (
        '<div class="grade-key">'
        "<p>Every classification on this site is mirrored from the authority that made it. "
        "The rung names are "
        '<a href="https://clinicalgenome.org/curation-activities/gene-disease-validity/">'
        "ClinGen&#x27;s</a>, glossed here from their Gene-Disease Validity SOP; Gene "
        "Curation Coalition submitters grade on the same ladder, and also use "
        "<em>Supportive</em>, which asserts an association without grading its evidence "
        "and which this atlas maps to no rung. "
        "<strong>No gene here is published on a disputed or refuted grade.</strong> A panel "
        "may still have disputed the same gene for a different disease. Where that disease "
        "is in scope here, by either path above, the record is in that gene&#x27;s validity "
        "table, marked <em>Disputed</em> or <em>Refuted</em>, and the gene carries a "
        "<em>conflicting evidence</em> flag. Where it is not, the record "
        "is not shown: a gene page lists only records for diseases in this atlas&#x27;s "
        "remit, so it is not a complete account of what ClinGen has said about the gene.</p>"
        f"<dl>{rows}</dl></div>"
    )


@dataclass(frozen=True)
class Link:
    """A table cell that is a hyperlink. `text` and `href` are both escaped."""

    text: str
    href: str


@dataclass(frozen=True)
class Markup:
    """A table cell whose HTML the *caller* has already made safe.

    **The escaping obligation moves to the caller, and this type exists so that
    is visible at the call site rather than implied.** Every other `Cell` is
    escaped by `_cell`; this one is not, so a `Markup` built from curated or
    mirrored text is an injection. `data_table`'s own docstring records why that
    matters here: a browse row's attributes are refused rather than escaped
    because escaping cannot make a name safe, and the same reasoning applies to
    a cell that carries markup.

    As a table cell, two callers today, both in `build_gene_index_page` (the
    forest plot's SVG text uses the type for the same obligation, outside any
    table): `pages._dot_strip`,
    whose HTML is a fixed set of `<span>` elements with class names from a
    literal map and `title` text escaped where it is built; and a contested
    gene's confidence cell, the grade (a vocabulary member, escaped anyway)
    followed by a `chip`, which escapes its own label. Nothing either
    interpolates comes from a mirror or from curation unescaped.
    """

    html: str


Cell = str | Link | Markup


@dataclass(frozen=True)
class Row:
    """One table row: its cells, and any `data-*` attributes it carries.

    `attributes` exists for the browse page, whose inline filter reads
    `data-lesion`, `data-confidence` and the rest off each `<tr>` rather than
    re-deriving them from the cell text. Names and values are escaped here, so
    no caller writes an attribute itself.
    """

    cells: tuple[Cell, ...]
    attributes: tuple[tuple[str, str], ...] = ()


def _cell(value: Cell) -> str:
    if isinstance(value, Link):
        return f'<a href="{html.escape(value.href)}">{html.escape(value.text)}</a>'
    # Deliberately unescaped; see `Markup`, which exists to make that visible.
    if isinstance(value, Markup):
        return value.html
    return html.escape(value)


def chip(text: str, kind: str = "", title: str = "") -> str:
    """A small labelled pill. `kind` selects a `.chip-*` modifier class.

    `title`, when given, is the pill's tooltip, escaped like the label. For a
    label too short to say what it measures: `sources disagree` is the one
    caller today (`pages._SOURCE_DISCORDANCE_TITLE`).
    """
    suffix = f" chip-{html.escape(kind)}" if kind else ""
    tooltip = f' title="{html.escape(title)}"' if title else ""
    return f'<span class="chip{suffix}"{tooltip}>{html.escape(text)}</span>'


def definition_list(items: Sequence[tuple[str, str]]) -> str:
    """A `<dl>` of label/value pairs, both escaped."""
    body = "".join(
        f"<dt>{html.escape(label)}</dt><dd>{html.escape(value)}</dd>" for label, value in items
    )
    return f"<dl>{body}</dl>"


def data_table(headers: Sequence[str], rows: Sequence[Row], table_id: str = "") -> str:
    """A `<table>` whose every cell and attribute is escaped by this function.

    An empty `rows` still renders the header, because a table with no rows and
    a table that was never written read the same to a caller otherwise.

    An attribute *name* is refused rather than escaped, because escaping cannot
    make it safe. `html.escape` rewrites `<`, `>`, `&`, `"` and `'` and leaves
    the space and the `=` alone, and a name sits outside any quote delimiter --
    so `("x onmouseover=alert(1) y", "v")` renders
    `<tr data-x onmouseover=alert(1) y="v">`, which the HTML tokenizer reads as
    three attributes, the second a live event handler, from a payload carrying
    no quote and no angle bracket. Measured with `html.parser` on 2026-08-04:
    `[('data-x', None), ('onmouseover', 'alert(1)'), ('y', 'v')]`.

    Every name the build passes is a literal in `pages.py`, so this is a guard
    on a bypassed gate in the same idiom as `encode_json`'s `allow_nan=False`:
    reaching it means a caller derived an attribute name from data, and that
    must fail rather than publish. `raise`, never `assert` -- `-O` strips
    `assert`.
    """
    ident = f' id="{html.escape(table_id)}"' if table_id else ""
    head = "".join(f'<th scope="col">{html.escape(header)}</th>' for header in headers)
    body = ""
    for row in rows:
        for name, _ in row.attributes:
            if not _ATTRIBUTE_NAME.fullmatch(name):
                raise ValueError(
                    f"data attribute name {name!r} is not [a-z][a-z0-9-]*; escaping cannot "
                    f"make a name safe, because a space in one opens a second attribute"
                )
        attributes = "".join(
            f' data-{name}="{html.escape(value)}"' for name, value in row.attributes
        )
        cells = "".join(f"<td>{_cell(cell)}</td>" for cell in row.cells)
        body += f"<tr{attributes}>{cells}</tr>"
    return f"<table{ident}><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"


def document(*, title: str, root: str, body: str, script: str = "") -> str:
    """The shell every page shares.

    `root` is the relative prefix from this page back to the site root -- "" for
    a page at the top level and "../" for one in `genes/`. Passed in rather than
    derived from a path, because this module has no opinion about where a caller
    put its output.

    `body` is HTML the caller has already assembled from the primitives above;
    it is the one argument not escaped here, which is why every primitive exists.
    `script` is likewise inlined verbatim and must be a literal from this module.

    `RESEARCH_USE_NOTICE` is emitted **before** `<main>`, so it precedes whatever
    the caller assembled in DOM order, in reading order and in the order a screen
    reader reaches it -- on a gene page that means it comes before the symbol,
    the HGNC id and the `definitive` chip rather than after them. It carries
    `role="note"` and no heading of its own: a gene page's `<h1>` lives inside
    `pages._rail`, and a heading here would put an `<h2>` above every `<h1>` on
    the site.
    """
    tail = f"<script>{script}</script>" if script else ""
    # Assembled ahead of the template only to keep the line under 100 characters;
    # the markup is one line in the output, exactly as it reads here.
    home = html.escape(root)
    nav = (
        f'<nav><a href="{home}index.html">CHD Atlas</a> &middot; '
        f'<a href="{home}genes/index.html">Genes</a> &middot; '
        f'<a href="{home}methods.html">How to read this atlas</a></nav>'
    )
    notice = f'<div class="notice" role="note"><p>{RESEARCH_USE_NOTICE}</p></div>'
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(title)}</title>
<style>{STYLESHEET}</style>
</head>
<body>
{nav}
{notice}
<main>
{body}
</main>
{tail}
</body>
</html>
"""
