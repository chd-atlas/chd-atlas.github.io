# src/chd_atlas/build/paths.py
"""Turning identifiers into path segments, in exactly one place.

The index and the bundle writer must agree on every filename or the site links to
files that were never written, so neither constructs a path itself.
"""

from __future__ import annotations

import re
from typing import Final

from chd_atlas.identifiers import HgncId

# Colon is illegal in a Windows filename and reserved in a URL authority; a slash
# would silently become a directory separator rather than fail, and
# `Zygosity.NOT_APPLICABLE` is the literal value "n/a". The dot is kept: it is
# legal in a filename, and `Emitter._write` is what rejects a path segment that
# is exactly "." or "..". The hyphen is kept because two grammars carry one and
# would otherwise be published in a form that matches nothing upstream:
# `UniprotAccession`'s isoform suffix ("Q99593-2") and `AccessionId`'s
# ArrayExpress form ("E-MTAB-1234"). Keeping both also holds `slug`'s injectivity
# argument to reasoning about two rewritten characters rather than four.
_UNSAFE: Final = re.compile(r"[^A-Za-z0-9.-]")


def slug(value: str) -> str:
    r"""A path-safe form of one identifier.

    Injective *as a function on strings*, over the ASCII values of every grammar
    in `identifiers.py` but one. Within ASCII the only characters this rewrites
    are the colon and the underscore, and no two values differ solely by swapping
    the two: the underscore appears only in `ContrastId`, which admits no
    capital, while every grammar bearing a colon starts with one. The hyphen and
    the dot are left alone, which is what keeps that list of two complete.

    Distinct strings are not distinct filenames, and the gap is not academic:
    the sole thing separating that worst-case pair is case, which is exactly what
    a case-insensitive filesystem discards. `HGNC:11604` (an `HgncId`) and
    `hgnc_11604` (a `ContrastId`) slug to "HGNC_11604" and "hgnc_11604", two
    files on the Linux CI and one on a curator's APFS. `Emitter._write` does not
    backstop the string-level guarantee here — its duplicate check compares
    Python strings — so it carries a second, casefolded check for this. Nothing
    can collide today regardless, because `gene_bundle_path` is the only path
    constructor in this module and only `HgncId` reaches it.

    The excepted grammar is `Doi` (`^10\.\d{4,9}/\S+$`), whose suffix admits
    every non-whitespace character, so "10.1000/ab/cd", "10.1000/ab:cd" and
    "10.1000/ab_cd" all reduce to one slug. Nothing builds a path from a DOI, and
    anything that starts to needs a rule that encodes rather than substitutes.

    "Within ASCII" is load-bearing rather than pedantic. `\d` is Unicode-aware in
    both Python's `re` and the rust engine pydantic validates with, so
    "HGNC:1160" suffixed with an Arabic-indic one and the same with an
    Arabic-indic two are two distinct values that `HgncId` accepts and that share
    the slug "HGNC_1160_". `\d` appears in every grammar but `ContrastId` —
    alone in spelling its alphabet out as `[a-z0-9]` — so this is not a corner of
    the schema: it reaches `HgncId`, the one grammar `gene_bundle_path` exists to
    serve. Two things keep it from quietly costing a gene page, neither of them
    this function: no curated file contains a non-ASCII digit today, and
    `Emitter._write` raises on the second write to a path rather than replacing
    the first, so a collision fails the build instead of publishing one page
    where two were expected. Note that the first of those is an observation
    about the data, not an invariant anything enforces.

    Three things the string-level claim does not cover, each with its own test:
    `Doi`, any non-ASCII value of any grammar, and the case collision above. A
    grammar added to `identifiers.py` and not to the sample fails
    `test_the_sample_covers_every_grammar_that_exists`, which is the intended
    outcome: only a reader can decide whether the argument above survives it.
    """
    return _UNSAFE.sub("_", value)


def gene_bundle_path(gene: HgncId) -> str:
    """Relative URL of one gene's bundle. Forward-slashed: consumers see a URL.

    `HgncId` rather than `str` so that the one sentence above holding the case
    collision at bay — only `HgncId` reaches this — is enforced by mypy instead
    of remembered: mypy rejects both the colliding `ContrastId` and a bare `str`.
    `slug` itself stays `str`, because vocabulary values such as `Zygosity`'s
    "n/a" belong to no grammar and must still be able to reach it.
    """
    return f"genes/{slug(gene)}.json"


def gene_page_path(gene: HgncId) -> str:
    """Relative URL of one gene's HTML page, beside its JSON bundle.

    Built from the same `slug` as `gene_bundle_path` so the pair cannot drift:
    a page whose name is derived differently from its bundle is a page linking
    to a file that was never written. `HgncId` rather than `str` for the reason
    `gene_bundle_path` takes one — only HGNC ids may name a file in `genes/`.
    """
    return f"genes/{slug(gene)}.html"


def register_page_path(gene: HgncId) -> str:
    """Relative URL of one gene's register (back) page, beside `genes/`.

    A folder of its own, not `genes/<slug>-register.html`, so the suite's and
    the consumers' `genes/HGNC_*.html` globs never read it as a gene page.
    """
    return f"register/{slug(gene)}.html"


# The whole-corpus payloads, named here because two modules have to agree on each
# and neither owns it: `literature.py` writes these files and `search.py`
# advertises two of them as the `path` a search result resolves to. Held apart as
# string literals in both, a rename in one published dead links from the other,
# with nothing in code tying the pair — `Emitter.checksums` records what was
# written, never what a payload promised. That is amendment A36; the dead-link
# sweep in `tests/test_built_site_is_consumable.py` detects the drift, and these
# constants are what make it unrepresentable.
#
# Gene bundles never needed this: `gene_bundle_path` above is already the single
# place that names one, which is why they carry no such risk.
PUBLICATIONS: Final = "publications.json"
FEATURED: Final = "featured.json"
PHENOTYPES: Final = "phenotypes.json"
DATASETS: Final = "datasets.json"
# Not a curated payload but a licensing obligation. `mirrors/sources.yaml`
# records HPO as `permitted_with_attribution`, and the site republishes its term
# labels and synonyms verbatim — so the attribution has to travel with the data,
# not sit only in a README a JSON consumer never opens.
SOURCES: Final = "sources.json"
# The sample collections every burden row names by bare id. Published because a
# consumer reading `case_cohorts: ["taa_cases"]` out of a bundle otherwise gets
# the numbers with none of the caveats that qualify them -- and `taa_cases` is
# 777 thoracic aortic aneurysm probands, a condition outside this atlas's scope.
COHORTS: Final = "cohorts.json"
# The CHD scope registry, and for each term the authority it is attributed to.
# Published because a bundle names a disease by bare MONDO id, and since the
# Orphanet scope path (2026-10-25) the reason a syndrome counts as CHD can be a
# frequency Orphanet records rather than a curation -- which a reader of the
# JSON has no other way to learn.
SCOPE: Final = "scope.json"

# The site's front door, not identifier-based like `gene_bundle_path` — there is
# exactly one, so nothing here builds it from a value. Named as a constant
# anyway, for the same reason as the four above: `build/landing.py` writes it
# and `tests/unit/test_build_landing.py` has to name the same file without
# either one hardcoding a string the other could drift from.
LANDING: Final = "index.html"

# The browse page. `genes/index.json` is its data and this is the page; Pages
# serves this at `/genes/`. Named here beside the other whole-site payloads
# because `pages.py` writes it and `landing.py` links to it, and nothing in code
# ties a written path to a path a payload promised.
GENE_INDEX_PAGE: Final = "genes/index.html"

# The site's reading notes: what every gene page assumes, written once (D50).
# `methods.py` writes it here. The shared nav and the gene pages' short forms
# link to it in literal markup, as the nav links `genes/index.html`, so what
# ties those links to this file is the real-build link sweep in
# `tests/test_built_site_is_consumable.py`, fragments included.
METHODS_PAGE: Final = "methods.html"
