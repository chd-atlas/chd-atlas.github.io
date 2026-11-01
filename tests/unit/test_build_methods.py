# tests/unit/test_build_methods.py
"""`methods.html`: what every gene page assumes, written once (D50).

Every sentence on this page was moved off a gene page, or states in full what a
gene page now says in one line, and the move is where it could go false: a
sentence's truth conditions travel with its position (CLAUDE.md section 4.27).
So each one is pinned here against a **literal**, sliced to its own section,
never against the constant that produced it (section 4.38) -- and the phrases
the move retired are asserted absent from the whole page, because a pasted old
constant is how they would come back.
"""

from __future__ import annotations

import html
import re
from pathlib import Path
from typing import Final

import pytest

from chd_atlas.build.burden import BurdenRow, cohort_registry, load_burden
from chd_atlas.build.concordance import ALPHA
from chd_atlas.build.emit import Emitter
from chd_atlas.build.methods import build_methods
from chd_atlas.build.paths import METHODS_PAGE
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.review import ReviewStatus
from chd_atlas.vocab import DatasetDesign, Technology

REPO_ROOT = Path(__file__).resolve().parents[2]


@pytest.fixture(scope="module")
def page(tmp_path_factory: pytest.TempPathFactory) -> str:
    """The page as built, with no cohorts: only `#cohorts` reads them, and
    `cohort_page` below builds it with a registry."""
    root = tmp_path_factory.mktemp("methods")
    build_methods(Emitter(root=root), cohorts={})
    return (root / METHODS_PAGE).read_text(encoding="utf-8")


# Inserted out of id order, because `burden.cohort_registry` returns the YAML's
# order and the section must sort. Insertion order is not seed-dependent, so a
# dropped sort fails here under every PYTHONHASHSEED (CLAUDE.md section 4.42),
# and the names sort differently from the ids (D, G, U against c, d, u), so a
# sort on the name fails too. Every curated field a reader sees carries a
# character that must be escaped. One cohort carries an open question (owner's
# standing rule, 2026-10-27) and two do not, so a status rendered on every
# cohort fails as surely as one rendered on none; it is the one with a url, so
# the note's place between the description and the link is pinned too.
_REGISTRY: Final = {
    "ukbb": Cohort(
        id="ukbb",
        name="UK Biobank & partners",
        description="Adults recruited at 40-69.",
        url="https://example.org/ukbb?a=1&b=2",
        review_status=ReviewStatus.PENDING_REVIEW,
        review_note="S2 lists <432> & the abstract 2,256.",
    ),
    "ddd": Cohort(
        id="ddd",
        name="Deciphering Developmental Disorders (DDD)",
        description="Enriched for <b>syndromic</b> CHD & more.",
    ),
    "cnchd": Cohort(id="cnchd", name="German Competence Network", description="German registry."),
}


@pytest.fixture(scope="module")
def cohort_page(tmp_path_factory: pytest.TempPathFactory) -> str:
    """The page as built from a three-cohort registry."""
    root = tmp_path_factory.mktemp("methods-cohorts")
    build_methods(Emitter(root=root), cohorts=_REGISTRY)
    return (root / METHODS_PAGE).read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def corpus() -> Corpus:
    """The committed curation, loaded the way `runner.build_site` loads it."""
    loaded, issues = load_curation(REPO_ROOT)
    assert issues == [], issues
    return loaded


@pytest.fixture(scope="module")
def registry(corpus: Corpus) -> dict[str, Cohort]:
    """The committed registry, loaded the way `runner.build_site` loads it."""
    return cohort_registry(corpus.cohorts)


@pytest.fixture(scope="module")
def real_page(tmp_path_factory: pytest.TempPathFactory, registry: dict[str, Cohort]) -> str:
    """The page as built from the committed registry: `#cohorts` as published.

    The fixtures above build it with no cohorts or three invented ones, so
    nothing a curator wrote reached any check on this page until 2026-10-26 --
    when the published page was found carrying a retired phrase ("the same
    children") from `pcgc`'s description.
    """
    root = tmp_path_factory.mktemp("methods-real")
    build_methods(Emitter(root=root), cohorts=registry)
    return (root / METHODS_PAGE).read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def mirrored() -> list[BurdenRow]:
    """Every row of the committed burden mirror."""
    rows = [row for rows in load_burden(REPO_ROOT).values() for row in rows]
    assert rows, "no burden rows loaded; every check on them would pass vacuously"
    return rows


def _section(page: str, anchor: str) -> str:
    """One section: from its `<h2 id=...>` to the next `<h2`, or to `</main>`.

    `</main>` rather than the end of the page (2026-10-26), so the last section
    can be pinned whole without the shell's closing tags -- or the newline the
    shell puts before them -- in its literal.
    """
    start = page.index(f'<h2 id="{anchor}">')
    ends = [i for i in (page.find("<h2", start + 1), page.find("</main>", start)) if i != -1]
    return (page[start : min(ends)] if ends else page[start:]).rstrip("\n")


def test_the_page_sits_at_the_site_root_under_its_own_title(page: str) -> None:
    """P1: one page at the root, named for what a reader comes to it for.

    `root=""` is what makes its own nav resolve from the site root; a page
    rendered with `"../"` would point every visitor's "Genes" link above the
    deploy root.
    """
    assert METHODS_PAGE == "methods.html"
    assert "<title>How to read this atlas — CHD Atlas</title>" in page
    assert "<h1>How to read this atlas</h1>" in page
    assert '<a href="genes/index.html">Genes</a>' in page
    assert '<a href="methods.html">How to read this atlas</a>' in page
    # The sections in the order the introduction names them and a gene page
    # renders them, and nothing else at that rank (written once all five
    # existed, 2026-10-26; `#curated-evidence` made six the same day, E3, in
    # gene-page order: a curated gene's evidence follows its validity table).
    assert re.findall(r'<h2 id="([^"]+)"', page) == [
        "scope",
        "validity",
        "curated-evidence",
        "rare-variants",
        "cohorts",
        "expression",
    ]
    assert page.count("<h2") == 6
    # The page's own account of what it holds (CLAUDE.md section 4.29: a
    # summary that omits what the page publishes passes every other check).
    # The cohort list joined it on 2026-10-26, and the same day it stopped
    # claiming the list is every collection a study "draws on": it is what the
    # rows name (`test_the_cohort_list_says_what_its_rows_name_and_what_they_do_not`).
    # The expression section joined it the same day, once the page held one,
    # and the curated evidence with its own section (E3).
    assert (
        "<p>What every gene page assumes, written once: which genes are here and what each "
        "classification means, how the validity records are listed, what the atlas's own "
        "curated evidence is, how to read the rare-variant results, every sample collection "
        "the rare-variant results name, and how to read the developmental expression "
        "figures.</p>"
    ) in page


def test_scope_states_the_rule_that_admits_a_gene_in_full(page: str) -> None:
    """The rule both page kinds stated in full until 2026-10-26, and now link to (D50).

    Moved from `pages._SCOPE_RULE`, which was pinned by
    `test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas` and,
    for its Orphanet clause, by
    `test_an_orphanet_scoped_syndrome_says_why_it_counts_as_chd_above_both_groups`.
    Those tests' recorded decisions travel with it:

    * **Both warrants, or the sentence describes a gate the site does not
      have.** Stating only the ClinGen half would tell a reader that every gene
      here was graded by an expert panel; 21 of the 119 were not (measured
      2026-10-26), and they are the ones whose provenance most needs stating.
    * **The negative is the load-bearing half.** The rule used to say a gene is
      definitive "for a disease in this atlas's CHD scope", which told 24
      published pages that this atlas decides what counts as congenital heart
      disease. It does not; scope is a selection among terms external
      authorities already use. Without the negative a future edit could
      restore the old sentence beside the new one and satisfy every positive.
    * **Sliced to the section, because page-wide the fragments could not
      fail.** The grade key's `ungraded` gloss carries "two or more Gene
      Curation Coalition submitters classify it Limited or better" in the same
      section's `<dl>`, so the slice is the rule's own `<p>`: rewriting the
      rule to state a **one**-submitter gate survived all 796 tests when the
      check was page-wide (measured 2026-08-06, CLAUDE.md section 4.19).

    Its last clause said the disease "is shown beside the gene", true on a gene
    page and on the gene index and unanchored here, so it now names both
    places -- and only for a gene a ClinGen panel graded: the 21 genes admitted
    on submitter agreement have no panel-named disease to show (2026-10-26).
    The full literal pins the wording; the fragments beside it are the
    claims any rewording must keep, so a deliberate edit to the literal is not
    also a silent edit to the gate it describes.

    **Its first part names a list, never a property of a disease**
    (2026-10-26). It read "for a disease that an external authority treats as
    congenital heart disease (or in which Orphanet records a structural heart
    malformation in at least 30% of cases)", which reads as a sufficient
    condition the build does not apply: SCP007 listed 239 unnominated terms
    that would qualify on Orphanet's frequencies, FLNA's periventricular
    nodular heterotopia among them, and Orphanet's own GenCC submissions name
    Ebstein anomaly, which is not a scope term. So the fragments now include
    what a term needs to be listed and that a qualifying disease is not in
    scope until it is; the retired clause is asserted absent from the whole
    page by `test_no_sentence_the_move_retired_comes_back_here`. Its ontology
    clauses are held to the set the build computes by
    `test_orphanet_scope.py::test_the_structural_set_the_build_applies_is_the_one_methods_html_states`.
    """
    section = _section(page, "scope")
    rule = re.search(r'<p class="scope-rule">.*?</p>', section, re.S)
    assert rule, "#scope carries no scope-rule element"
    sentence = rule.group(0)
    assert sentence == (
        '<p class="scope-rule">A gene is published here when an external authority already '
        "treats it as a congenital heart disease gene: either a <strong>ClinGen expert "
        "panel</strong> classifies it <strong>Limited or better</strong> for a disease in scope "
        "here, or <strong>two or more Gene Curation Coalition submitters</strong> classify it "
        "Limited or better for one and no ClinGen panel disputes it. A disease is in scope when "
        "it is one of the terms this atlas lists, and a term is listed only if an external "
        "authority treats it as congenital heart disease, or if a curator nominates a syndrome "
        "and <strong>Orphanet records a structural heart malformation in at least 30% of its "
        "cases</strong> (for a subtype, of the syndrome it belongs to); a disease that would "
        "qualify is not in scope until it is listed. Here a structural heart malformation is a "
        "finding under the Human Phenotype Ontology's <em>Abnormal heart morphology</em> "
        "(HP:0001627), less branches of it excluded as not malformations (cardiomyopathy, "
        "tumours, the pericardium and the like), plus eight lesions the ontology files "
        "elsewhere, with any subtypes: coarctation of the aorta, a patent ductus "
        "arteriosus, pulmonic stenosis, aortic valve stenosis, supravalvular aortic stenosis, "
        "truncus arteriosus, an interrupted aortic arch and anomalous pulmonary venous return. "
        "The findings counted still include some that are not malformations, such as a "
        "calcified or prolapsing valve, and aortic valve stenosis can be acquired; other "
        "lesions the ontology files elsewhere, among them peripheral pulmonary artery "
        "stenosis and an aortopulmonary window, are not counted. "
        "<strong>No disease is in scope on this atlas's own judgement</strong>, and neither "
        "is any gene: every scope term names the authority it came from "
        "(<code>scope.json</code>, from <code>curation/chd_scope.yaml</code>), and every gene "
        "page names the authority that admitted the gene. A classification is not the same "
        "as definitive for congenital heart disease: for many of these genes the disease "
        "graded is a syndrome of which a cardiac lesion is one feature, and for a gene a "
        "ClinGen panel graded, the disease the panel named is shown beside the gene, on "
        "its page and in the gene index.</p>"
    )
    for claim in (
        "Limited or better",
        "two or more Gene Curation Coalition submitters</strong> classify it Limited or better",
        "no ClinGen panel disputes it",
        "a term is listed only if an external authority treats it as congenital heart disease",
        "if a curator nominates a syndrome and <strong>Orphanet records",
        "a disease that would qualify is not in scope until it is listed",
        "No disease is in scope on this atlas's own judgement",
        "not the same as definitive for congenital heart disease",
        "(<code>scope.json</code>, from <code>curation/chd_scope.yaml</code>)",
    ):
        assert claim in sentence, f"the scope rule no longer says {claim!r}"


def test_scope_glosses_every_rung_as_one_label_chip_and_gloss(page: str) -> None:
    """The grade key, open, under the rule (P3). Moved from both page kinds (D50).

    Moved from `test_both_page_kinds_gloss_every_grade_they_can_render`, with its
    recorded decisions: a colour ramp with no key is a ranking the reader has to
    guess at, and since 2026-08-06 the site renders four rungs and an ungraded
    state -- `limited`, a panel saying the case is not yet made, is the
    plurality. Each rung's gloss is checked by a distinctive phrase rather than
    in full, so this fails when a rung stops being explained and not when a
    comma moves.

    **The label, the chip class and the gloss are asserted as one triple, and
    that is the whole point of this test.** It checked only that seven phrases
    appeared *somewhere on the page* until 2026-08-06, which says nothing about
    which rung each is attached to. Measured that day: painting `moderate` with
    `strong`'s chip -- inverting the colour ramp mid-ladder -- **survived all 796
    tests**, as did swapping two glosses. Scoped to the key's `<dl>`, because
    these words also occur in the scope rule above it (section 4.19).

    **Rendered open, not folded** (P3). On a gene page the key was a collapsed
    `<details>` because it was not what the page was for; here it is. Asserted
    on the whole page since 2026-10-26: `_section` began stopping at the next
    `<h2`, and a slice of `#scope` alone would let a fold appear anywhere below
    it -- nothing on this page is meant to be folded.
    """
    section = _section(page, "scope")
    assert "<details" not in page, "a fold on the page that exists to be read"
    key = re.search(r'<div class="grade-key">.*?</div>', section, re.S)
    assert key, "#scope carries no grade key"
    legend = re.search(r"<dl>(.*?)</dl>", key.group(0), re.S)
    assert legend, "the grade key carries no rows"

    # (label a reader sees, chip class that colours it, phrase unique to its gloss)
    rungs = (
        ("definitive", "chip-definitive", "upheld over time"),
        ("strong", "chip-strong", "without the years of replication"),
        ("moderate", "chip-moderate", "several unrelated probands"),
        ("limited", "chip-limited", "not yet made"),
        ("no known association", "chip-warn", "found no evidence of a relationship"),
        (
            "not classified",
            "chip-ungraded",
            "two or more Gene Curation Coalition submitters classify it Limited or better",
        ),
    )
    pairs = re.findall(r"<dt>(.*?)</dt><dd>(.*?)</dd>", legend.group(1), re.S)
    assert len(pairs) == len(rungs), f"the key glosses {len(pairs)} rungs, expected {len(rungs)}"
    for (label, chip, gloss), (dt, dd) in zip(rungs, pairs, strict=True):
        assert chip in dt, f"{label} is not chipped {chip} -- got {dt}"
        assert f">{label}<" in dt, f"expected the label {label}, got {dt}"
        assert gloss in dd, f"{label} is not glossed by {gloss!r}"


def test_scope_states_where_a_dispute_is_listed_without_pointing_below(page: str) -> None:
    """The dispute paragraph, reworded for a page with no validity table on it (P3).

    Moved from `test_the_grade_key_no_longer_denies_an_in_scope_dispute_it_now_displays`.
    The sentence went false once already, when Costello and CFC syndrome entered
    scope on 2026-10-25: it read "No gene here was disputed or refuted by an
    expert panel for a disease an external authority treats as congenital heart
    disease", and ClinGen's RASopathy panel disputes six published genes for
    those two syndromes -- PTPN11's page carried the denial beside two
    `Disputed` rows.

    **Reworded again on the move (2026-10-26, D50).** It said "the record is in
    the table below", which was false on the gene index (no validity table) and
    is false here; it now names the gene's validity table. "This page is not a
    complete account" became "a gene page", for the same reason.

    **And its condition named one path of two.** "Where that disease is one an
    external authority treats as congenital heart disease" excluded the
    Orphanet path -- and Costello and cardiofaciocutaneous syndrome, the
    diseases all six disputed genes (BRAF, KRAS, MAP2K1, PTPN11, RAF1, SOS1) are
    disputed for, are in scope only through it (`admitted_by:
    orphanet_phenotype_frequency`). It now reads "in scope here, by either path
    above", which `test_the_key_follows_the_rule_its_dispute_sentence_points_back_to`
    keeps true.
    """
    section = _section(page, "scope")
    key = re.search(r'<div class="grade-key"><p>(.*?)</p>', section, re.S)
    assert key, "#scope carries no grade-key paragraph"
    assert key.group(1) == (
        "Every classification on this site is mirrored from the authority that made it. "
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
        "remit, so it is not a complete account of what ClinGen has said about the gene."
    )


def test_validity_states_how_a_gene_page_lists_its_records(page: str) -> None:
    """How the validity table is listed, in full; each gene page links it (D50).

    A gene page's validity intro says only records for in-scope diseases are
    listed, and links here for the rest: the two groups, which of ClinGen's own
    GenCC copies are listed, and what is not shown at all.

    The copy sentence describes `pages._repeats_the_panel`, read 2026-10-26: a
    copy is hidden only when a panel record matches it on disease, mapped grade,
    inheritance (in the page's vocabulary) and calendar date, and a copy with no
    date matches nothing -- so it is listed where any of the four differs, where
    it has no date, or where there is no panel record to match.

    **Two clauses differ from the plan's text, each forced by a measurement**
    (2026-10-26, on the committed mirrors and a real build):

    * "Each gene's JSON keeps every record" read, beside the next paragraph's
      "that record is not shown", as though the JSON held the out-of-scope
      records. It holds none: 591 bundle records, 0 for a disease outside
      scope (`validity.gene_validity` drops them before a gene is keyed). It
      does hold every copy -- all 121 -- so it now says that.
    * "so a gene page is not every classification the gene has" was false for
      43 of the 119 published genes, which have no ClinGen or GenCC record for
      any disease outside scope. It says "may not be", as the gene page does.

    **Reworded after review, 2026-10-26.** Its opening sentence ("Every
    classification in this atlas is an upstream panel's or submitter's,
    republished with its provenance intact") repeated the grade key's first
    sentence two paragraphs up, so it is written once, in the key. "the
    ClinGen expert panel's own records" became "ClinGen expert panels'", and
    "A panel may have classified" became "A panel or submitter": the scope
    filter drops GenCC rows too, and most out-of-scope records are GenCC's (76
    published genes have one, against 27 with an out-of-scope ClinGen record).

    **The third paragraph glosses the card's two warning chips (2026-10-26).**
    `sources disagree` rendered on six pages with no title and was mentioned 0
    times on this page, and on all six ClinGen and GenCC agree on every
    gene-disease pair they both classify: the flag is per gene, by the owner's
    decision, and fires because one source disputes the gene for one disease
    while the other supports it for another. Each clause is read off the code
    (`vocab.has_conflicting_evidence`, `validity._has_source_discordance`), and
    the paragraph states what the per-gene flag permits rather than narrowing
    the flag.
    """
    section = _section(page, "validity")
    assert section == (
        '<h2 id="validity">Mirrored validity records</h2>'
        "<p>A gene page lists its classifications in two groups: ClinGen expert panels' own "
        "records first, then other curators' submissions to the Gene Curation Coalition "
        "(GenCC). ClinGen also submits its own classifications to GenCC; a gene page lists "
        "such a copy only where its mapped grade, date, inheritance or disease differs from "
        "the panel's own record, or where it carries no date or has no panel record to "
        "match. Each gene's JSON keeps every copy, listed or not.</p>"
        "<p>Only records for diseases in this atlas's remit are listed. A panel or submitter "
        "may have classified the same gene for a disease outside it, and that record is not "
        "shown, so a gene page may not be every classification the gene has.</p>"
        "<p>A gene page's card may carry two flags beside its grade. <em>conflicting "
        "evidence</em> means that among the gene's records for diseases in scope, at least one "
        "is <em>Disputed</em> or <em>Refuted</em> and at least one classifies the gene Limited "
        "or better, whoever made each and whatever disease each is for; the gene index marks "
        "the same genes beside their confidence. <em>sources disagree</em> is narrower: a "
        "record that disputes or refutes the gene and one that classifies it Limited or better "
        "come from different sources, one a ClinGen expert panel's own record and the other a "
        "submission to GenCC (ClinGen's own copies count as GenCC's here). It too is a flag on "
        "the gene, not on a disease: the two records may be for different diseases, and the "
        "two sources may agree on every disease they both classify. Neither flag counts a "
        "<em>No Known Disease Relationship</em> record on either side, nor GenCC's "
        "<em>Supportive</em>, which grades no evidence.</p>"
    )
    # Written once on this page, in the grade key.
    assert page.count("Every classification on this site is mirrored") == 1
    assert "Every classification in this atlas is an upstream" not in page


def test_curated_evidence_says_what_an_assertion_is_and_that_it_changes_nothing_mirrored(
    page: str,
) -> None:
    """The two sections only a curated gene's page carries (owner decision 2026-10-26, E3).

    TBX5's page has "Curated evidence" and "Publications" and this page
    explained every other section but those, while the card's "How to read
    this page" promised it. Each clause read off the code that day:
    `LesionAssertion` (lesion groups, inheritance, mechanism, syndromic status
    -- the four rows the page's definition list renders), `Evidence` (one
    required publication, a class, a strength of four values, the curator's
    summary -- the table's four columns), `validity.published_genes` and
    `derive.gene_facts` (the gate and the headline read the validity records
    alone, D37), and `pages._evidence_section` (the publications cited, with
    title and year).

    **One clause of the owner's draft was false.** "Each claim tied to a
    publication and quoted evidence": an evidence item's `summary` is the
    curator's summary, not a quotation, and the record locates the passage
    instead (`locator`), which no page renders. It says "summary".

    Pinned whole, as a literal, and it names no count: a second curated gene
    must not make it stale.
    """
    assert _section(page, "curated-evidence") == (
        '<h2 id="curated-evidence">Curated evidence</h2>'
        "<p>Where this atlas has curated a gene itself, its page carries a <strong>Curated "
        "evidence</strong> section. Each entry is a <strong>lesion assertion</strong>: a "
        "curator's record that the gene belongs to one or more congenital heart disease lesion "
        "groups, with its inheritance, its mechanism and whether it is isolated, syndromic or "
        "both. Each item of evidence behind an assertion cites one publication, and gives the "
        "class of the evidence, its strength as the curator grades it (strong, moderate, "
        "supporting or contradictory) and the curator's summary of what the publication "
        "shows.</p>"
        "<p>An assertion is this atlas's own work, and it <strong>changes nothing "
        "mirrored</strong>: a gene is published, and its grade is headlined, on the upstream "
        "validity records alone, so no assertion admits a gene or moves its grade. A curated "
        "gene's page also has a <strong>Publications</strong> section, which lists every "
        "publication its evidence cites, with its title and year.</p>"
    )


def test_rare_variants_states_what_every_gene_page_folded_or_shortened(page: str) -> None:
    """The rare-variant reading notes, moved off every gene page (2026-10-26, D50).

    Nine paragraphs, each pinned whole. The first is the front page's two
    definitions (below). The next two are the full forms of the two notices a
    gene page keeps in one line (`pages._MIRRORED_NOTICE`, which links here,
    and `pages._POOLING_NOTICE`); the two after those were the page's folded
    "How to read these numbers" notes; the next was the second sentence of its
    composite note; the next says what `all cases` is; the last two are below.
    Each was re-read for a page with no table on it: "In each study's table"
    anchors what the units sentence describes, and "the rows above it" became
    "beside it".

    Recorded decisions that moved with them:

    * **The absent-cell sentence must not say why a row is absent.** Moved from
      `test_the_section_says_what_a_missing_row_means_without_asserting_why`.
      Until 2026-08-05 it asserted "no variant of that class was seen in either
      group" -- measured over PMID:42230622 alone, and false for the other two
      studies: PMID:34324492 tests one consequence class by construction, and
      PMID:40127276 observed 14,364 synonymous variants (its Dataset S4) and
      published no synonymous row. That phrase is asserted absent from this page
      by `test_no_sentence_the_move_retired_comes_back_here`, and from the gene
      page by `test_the_general_reading_notes_are_on_methods_html_and_not_on_the_gene_page`.
    * **The negative control carries its own condition.** Moved from
      `test_the_synonymous_negative_control_is_explained_only_where_one_exists`.
      On a gene page it rendered only where a synonymous row did, because 2 of
      the 3 studies publish none (every synonymous row is PMID:42230622's); here
      it is stated once, with that condition in its own words ("Only a study
      that publishes one has it").
    * **The pooling paragraph is policy, in the present tense about nothing.**
      "these cohorts overlap" sent a one-study page's reader hunting for a
      second study; it is asserted absent here too.

    **Three clauses are narrower than the plan's text, each measured 2026-10-26
    on a real build of 4233eca** (`methods._RARE_VARIANTS` records the detail):
    "in a gene page's rare-variant section" became "in a gene page's study
    tables", because the evidence matrix in that section shows tallies the atlas
    derives ("n = N" on 86 pages, "N rows" on 48); "republished exactly as
    published" lost "exactly", because 376 of the 1,210 published raw p-values
    render rounded to three significant figures; and the units sentence says
    what `controls` shows on a row tested against a mutation model, because 108
    such cells on 36 pages read "N expected", not observed / total.

    **Reworded after review, the same day.** "children" became "people": a
    shared collection may be an adult control cohort. The overlap clause says
    what the page can know -- collections its studies share *as this atlas has
    recorded them* -- because `shared_cohorts` sees only shared registered ids,
    and two studies may share individuals no record names. The synonymous
    paragraph lost "the same uncorrected scale" (true of one study, stated in
    general) and "a scan of that size" (no antecedent here). The units
    parenthesis keys on what `pages._controls` keys on: the comparator, not the
    variant's origin.

    **Pinned whole, not paragraph by paragraph.** A `re.findall(r"<p>...")`
    over the section skips any `<p class=...>` or `<details>` appended to it:
    measured in review, a composite note, a disclosure and a "How to read these
    numbers" fold appended to this section all passed.

    **Two paragraphs were added with the forest key's fold (2026-10-26,
    D57)**: the key each forest plot now folds behind "How to read this
    chart", stated once in general, and what a block with no effect measure
    shows instead of a plot. Every mark clause was measured against a real
    build of 9f24190 before it was written in general: all 29 right arrows
    are odds ratios with 0 control carriers and a published lower bound; all
    257 left arrows have 0 case carriers; every plotted study either corrects
    every row or none. One clause is narrower than the plan's text: "contains
    other rows in the same panel" was false for 10 union-tagged rows, on 10
    panels, that contain exactly one other drawn row -- "at least one other
    row", as each panel's own key says.

    **Three changes after review (2026-10-26).** The composite paragraph
    explained `damaging` and nothing said what the `all cases` case group is,
    though a forest key says it in one clause; the `all cases` paragraph does, in
    general, and `test_an_all_cases_row_is_the_union_of_the_case_groups_its_
    study_splits` holds it to the mirror. "these have none" had no antecedent
    nearer than "a forest plot", so it names this atlas's. And the muted row was
    singular, while every panel that draws one draws three -- one synonymous row
    per case group, 113 panels on a real build -- so it is "one of a study's
    synonymous rows".

    **One sentence added the same day**, when a marker's fill began reading
    `concordance.survives_own_correction` for every row but a synonymous one:
    a depleted row that cleared its correction is drawn hollow, and "a hollow
    marker is a row that did not" survive is false of it. Latent: 0 published
    rows are in that state.

    **Two clauses narrowed after a later review (2026-10-26).** That sentence
    said "a row", and a synonymous row that survived its correction is drawn
    filled and muted, not hollow, so it says "a non-synonymous row". And the
    case-group paragraph called an unsplit study's `all cases` its "whole case
    set", a denominator equivalence no study's rows bear out (PMID:34324492's
    count 6,746 of its 7,958 cases); it now says only that the group is not a
    subgroup, and the old phrase is asserted absent.

    **The paragraph after the definitions states the one rule the atlas applies** (owner
    decision 2026-10-26, W2). It said the atlas "adds no assessment" of the
    studies' numbers, while the card, the matrix and the forest fill each
    decide what survives by `concordance.survives_own_correction`: not
    synonymous, a raw p below 0.05 and a result pointing to enrichment in
    cases, then a study-published corrected p below 0.05, both strictly. The
    owner's draft gave the rule as "a row whose study-published corrected p is
    below 0.05", which is the second half alone: a depleted row, or a
    synonymous one, clears that and survives nothing.

    **0.05 is the atlas's, and the page renders it from `concordance.ALPHA`**
    (review, 2026-10-27). It read "at the 0.05 the studies themselves use",
    which is false of a raw p: PMID:42230622's block tells a reader to judge
    raw p against its whole scan. The literal below pins the published number
    and the assertion beside it pins the rule's, so a change to either fails
    here until the other is reconsidered.

    **It names how wide a survival is** (owner decision 2026-10-26, R3): a
    correction over a gene panel is not a genome-wide one, and the card says
    which, beside the q. PMID:40127276's per-arm q is Benjamini-Hochberg over
    its 248 panel genes, reproduced from Dataset S6 on every arm; the card
    pins on TBX5 are in `test_built_site_is_consumable.py`.

    **The section opens with the front page's two definitions** (doorway page,
    2026-10-27), which `index.html`'s hero note now links here instead of
    stating: what one rare-variant statistic is (D49: all three analyses, and
    "copy-number variants" so a future non-deletion CNV analysis does not make
    it false) and what one independent dataset is.
    """
    assert ALPHA == 0.05, "the rule's threshold moved; re-read #rare-variants' rule paragraph"
    expected = [
        "A <strong>rare-variant statistic</strong> is one study's rare-variant test of one "
        "gene under one design: a case-control burden test of rare SNVs and indels or of "
        "copy-number variants, or a de novo enrichment test, which are different analyses and "
        "are named as such on each gene's page. An <strong>independent dataset</strong> is a "
        "cohort family rather than a paper: two studies drawing on the same sample "
        "collection, or on collections this atlas records as overlapping, can describe the "
        "same people and count once. How many datasets agree about a gene is shown on every "
        "row of the browse page, and it is a count of evidence &mdash; never a verdict on the "
        "gene.",
        "Every count and statistic in a gene page's study tables is an upstream study's, "
        "republished as published. The atlas computes none of them. It applies one rule of "
        "its own, at 0.05: a row that is not synonymous "
        "<strong>supports</strong> the gene when its result points to enrichment in cases and "
        "its p is below 0.05, and <strong>survives</strong> its study's correction when the "
        "corrected p its study published is below 0.05 as well; a row whose study published "
        "no corrected p survives none. A study's correction covers only the genes it was run "
        "over: where those were a gene panel rather than the genome, the card says so beside "
        "the corrected p, and marks any surviving row whose raw p is also below the threshold "
        "the study itself states for genome-wide significance. A gene page's card lists the "
        "rows that survive, the "
        "evidence matrix marks a dataset whose rows survive apart from one whose rows only "
        "support, and a forest plot draws a surviving row filled. Where a study is by an "
        "author of this atlas, its block says so.",
        "The atlas computes <strong>no pooled statistic across studies</strong>. Where two "
        "studies draw on the same sample collection they describe partly the same people, and "
        "where they draw on collections this atlas records as overlapping they may, so "
        "combining their results could count those people twice; a gene page names any "
        "sample collection its studies share, and any overlap between theirs, as this atlas "
        "has recorded them.",
        "In each study's table, <code>cases</code> and <code>controls</code> each show "
        "<em>observed / total</em> (for a row tested against a mutation model, "
        "<code>controls</code> is instead the count that model expected), and <strong>every "
        "cell names what it counted</strong> &mdash; carriers, alleles, or de novo mutations. "
        "The three are not interchangeable: a study counting alleles counts a person with two "
        "qualifying variants twice, and its denominator is roughly twice the number of people "
        "sequenced. A consequence class or case group with no row is one the study "
        "<strong>did not report</strong>, which is not a null result and does not always "
        "mean no such variant was seen &mdash; studies differ in which classes their "
        "gene-level tables cover.",
        "The <strong>synonymous</strong> row is a study's own negative control &mdash; "
        "synonymous variants should show no enrichment. Read it on the same scale as the rows "
        "beside it: a <em>strongly</em> enriched synonymous row is a warning about that "
        "gene's comparison, while a nominally significant one is what a scan over thousands "
        "of genes produces by chance. <strong>Only a study that publishes one has "
        "it</strong>; a table with no synonymous row offers no such check.",
        "A <strong>damaging (LOF + missense)</strong> row is the <strong>union</strong> of "
        "the study's loss-of-function and damaging-missense rows, not a third independent "
        "result &mdash; the same variants are counted in both. It is shown because a study "
        "that defines its findings by the composite publishes a p-value for it that is not "
        "derivable from the other two.",
        "Where a study splits its cases into <strong>syndromic</strong> and "
        "<strong>non-syndromic</strong> groups, its <strong>all cases</strong> group is "
        "likewise the union of the two, and an <strong>all cases</strong> row is not a result "
        "independent of theirs. Where a study does not split its cases, <strong>all "
        "cases</strong> is its only case group: every case the row counts, not a subgroup.",
        "A forest plot shows each row's effect and interval on a logarithmic axis. A "
        "<strong>filled</strong> marker is a row that survived its own study's correction for "
        "multiple testing; a <strong>hollow</strong> marker is a row that did not, or whose "
        "study published no correction at all &mdash; the atlas computes none of its own. A "
        "non-synonymous row that survived its correction without pointing to enrichment in "
        "cases is drawn hollow too: it is not evidence for the gene. A forest plot usually "
        "ends in a summary diamond; this atlas's have none, because no pooled estimate is "
        "computed anywhere on this site. A <strong>hatched</strong> bar "
        "means no interval was published for that row, not a wide one. An arrow to the "
        "<strong>right</strong> is an effect unbounded above &mdash; no control carried the "
        "variant &mdash; so the study's lower bound is the whole finding; an arrow to the "
        "<strong>left</strong> is an effect of zero &mdash; no case carried the variant "
        "&mdash; which a logarithmic axis has no position for. A row tagged "
        "<strong>union</strong> contains at least one other row in the same panel and is not "
        "an independent finding. A <strong>muted</strong> row is one of a study's synonymous "
        "rows, drawn apart from the results it calibrates.",
        "A study that publishes <strong>no effect measure</strong> for a gene &mdash; no "
        "odds ratio, enrichment or rate ratio &mdash; has nothing to place on an effect axis, "
        "so its block shows its counts and p-values in a table instead of a plot.",
    ]
    assert _section(page, "rare-variants") == (
        '<h2 id="rare-variants">Rare-variant association</h2>'
        + "".join(f"<p>{paragraph}</p>" for paragraph in expected)
    )


def test_an_all_cases_row_is_the_union_of_the_case_groups_its_study_splits(
    mirrored: list[BurdenRow],
) -> None:
    """`#rare-variants`'s case-group sentence, held to the mirror it generalises.

    It says, of every study, that `all cases` is the union of `syndromic` and
    `non-syndromic` where the study splits its cases, and its only case group --
    not a subgroup -- where it does not. Measured on the committed mirror
    (2026-10-26): all three studies publish `all`; only PMID:42230622 splits,
    into 1,471 syndromic and 2,405 non-syndromic cases; and on all 470 cells
    where it publishes all three strata, the `all` row's cases and carriers are
    the sum of the other two's and its controls are theirs. A study that arrives
    with a third case group, or an `all` that is not the sum, makes the sentence
    false, and fails here rather than on a reader.
    """
    strata: dict[str, set[str]] = {}
    cells: dict[tuple[object, ...], dict[str, BurdenRow]] = {}
    for row in mirrored:
        strata.setdefault(row.study, set()).add(row.cohort_stratum)
        cell = (row.study, row.gene, row.variant_class, row.consequence_class, row.origin)
        key = (*cell, row.lesion_group, row.count_unit, row.maf_max, row.comparator)
        assert row.cohort_stratum not in cells.setdefault(key, {}), f"two rows share {key}"
        cells[key][row.cohort_stratum] = row
    for study, found in strata.items():
        assert found in ({"all"}, {"all", "syndromic", "nonsyndromic"}), (study, found)
    split = [cell for cell in cells.values() if cell.keys() == {"all", "syndromic", "nonsyndromic"}]
    assert split, "no study splits its cases; the sentence's first half is untested"
    for cell in split:
        whole, parts = cell["all"], (cell["syndromic"], cell["nonsyndromic"])
        assert whole.n_cases == sum(part.n_cases for part in parts), whole
        assert whole.n_case_carriers == sum(part.n_case_carriers for part in parts), whole
        assert all(part.n_controls == whole.n_controls for part in parts), whole


def test_cohorts_lists_the_whole_registry_sorted_and_escaped(cohort_page: str) -> None:
    """Every cohort, with its curated description, once, for a reader (P4, 2026-10-26).

    A gene page folds the notes of only the collections its own rows cite
    (`pages._cohort_notes`) and ends the fold on a link here; this is the
    whole registry, the records `cohorts.json` publishes for a program. The
    fold stays on the gene page: CLAUDE.md section 9 forbids narrowing where a
    cohort caveat reaches, and this list adds a place rather than moving one.

    Pinned whole, as one literal. **Sorted by id**, because the registry
    arrives in the YAML's order and a published list from a mapping is the
    caller's to sort (CLAUDE.md section 9); the fixture's insertion order and
    its names both disagree with id order, so a missing sort and a sort on the
    wrong key both fail. **Escaped, every field**: names, descriptions and
    urls are free text a curator types. The url is the cohort's public page,
    so it is a link inside its own description rather than a bare string, and
    a cohort without one gets nothing rather than an empty link.

    The class on the `<dl>` is layout: the site's base `dl` is a two-column
    grid sized to its widest label, which suits the grade key's chips and
    would give a cohort's full name most of the line here. The link reads
    "website" (2026-10-26): it leads to the collection's home page, and
    "source" beside a description read as the citation for it.

    **A record's open question follows its description** (owner's standing
    rule, 2026-10-27), before the link, escaped like every curated field, and
    only on the record carrying one; the sentence saying what its label marks
    precedes the list (review, 2026-10-27, M8; the next test holds when it
    renders).

    The intro's numbers are pinned against the mirror and the registry by
    `test_the_cohort_list_says_what_its_rows_name_and_what_they_do_not`; here
    they are only the literal the fixture renders, since the intro is a
    constant and the fixture's three cohorts are not the ones it counts.
    """
    assert _section(cohort_page, "cohorts") == (
        '<h2 id="cohorts">Cohorts</h2>'
        "<p>Every sample collection a rare-variant row in this atlas names, with the caveats "
        "its curator recorded. A row may not name every collection its study drew on: the "
        "rows for PMID:34324492 name 11 of the 189 case sources its S2 Table lists, and those "
        "11 hold 6,211 of its 7,958 cases. Two studies that name the same collection describe "
        "partly the same people; two that name collections this registry records as "
        "overlapping may; and two that name neither may still share people this registry "
        "does not record. The same records are published as "
        '<a href="cohorts.json"><code>cohorts.json</code></a>.</p>'
        "<p>A record marked <strong>Pending review</strong> carries a question about it that "
        "this atlas leaves open for an expert: its note states what is unresolved, and the "
        "atlas takes no side.</p>"
        '<dl class="cohort-list">'
        "<dt>German Competence Network (<code>cnchd</code>)</dt><dd>German registry.</dd>"
        "<dt>Deciphering Developmental Disorders (DDD) (<code>ddd</code>)</dt>"
        "<dd>Enriched for &lt;b&gt;syndromic&lt;/b&gt; CHD &amp; more.</dd>"
        "<dt>UK Biobank &amp; partners (<code>ukbb</code>)</dt>"
        "<dd>Adults recruited at 40-69. <strong>Pending review:</strong> S2 lists "
        '&lt;432&gt; &amp; the abstract 2,256. <a href="https://example.org/ukbb?a=1&amp;b=2">'
        "website</a></dd>"
        "</dl>"
    )


# What `#cohorts` says each review label marks, as literals (CLAUDE.md section
# 4.38); the pending one is also inside the literal above.
_PENDING_KEY: Final = (
    "A record marked <strong>Pending review</strong> carries a question about it that this "
    "atlas leaves open for an expert: its note states what is unresolved, and the atlas "
    "takes no side."
)
_NOT_ACCESSIBLE_KEY: Final = (
    "A record marked <strong>Source not accessible</strong> rests on a source that could not "
    "be read: its note states what that leaves unchecked, and the atlas does not guess at it."
)


@pytest.mark.parametrize(
    ("marked", "key"),
    [
        pytest.param({}, "", id="none-marked"),
        pytest.param(
            {"ddd": ReviewStatus.NOT_ACCESSIBLE},
            f"<p>{_NOT_ACCESSIBLE_KEY}</p>",
            id="not-accessible-only",
        ),
        pytest.param(
            {
                "ukbb": ReviewStatus.NOT_ACCESSIBLE,
                "ddd": ReviewStatus.PENDING_REVIEW,
                "cnchd": ReviewStatus.NOT_ACCESSIBLE,
            },
            f"<p>{_PENDING_KEY} {_NOT_ACCESSIBLE_KEY}</p>",
            id="both",
        ),
    ],
)
def test_the_cohort_list_says_what_a_review_label_marks_only_for_labels_it_prints(
    tmp_path: Path, marked: dict[str, ReviewStatus], key: str
) -> None:
    """What "Pending review:" means, said where the label is printed (review, 2026-10-27, M8).

    The list printed the label after two descriptions, and every gene page's
    cohort fold links here, but no page said what it marks. A sentence
    between the intro and the list now does, one per status some record
    carries, in the vocabulary's order -- and none where no record carries
    one, since a key to a label the page does not print defines nothing a
    reader meets. `cohort_page` above pins the pending sentence alone, from a
    registry marking only `ukbb`. Here the same registry with nothing marked,
    with only the status the real registry has never carried, and with both,
    the pending record (`ddd`) between two not-accessible ones: `ukbb` comes
    first in the registry's insertion order and `cnchd` first by id, so a
    sentence order taken from either, rather than from the vocabulary, puts
    "Source not accessible" first and fails, and one sentence per record
    rather than per status prints it twice.
    """
    registry = {
        identifier: cohort.model_copy(
            update={
                "review_status": marked.get(identifier),
                "review_note": "What is unresolved." if identifier in marked else None,
            }
        )
        for identifier, cohort in _REGISTRY.items()
    }
    build_methods(Emitter(root=tmp_path), cohorts=registry)
    section = _section((tmp_path / METHODS_PAGE).read_text(encoding="utf-8"), "cohorts")
    intro_end = section.index("</p>") + len("</p>")
    assert section[intro_end : section.index('<dl class="cohort-list">')] == key


def test_the_cohort_list_says_what_its_rows_name_and_what_they_do_not(
    real_page: str, registry: dict[str, Cohort], mirrored: list[BurdenRow]
) -> None:
    """`#cohorts`'s intro, on the published page, held to the rows it describes.

    **It claimed a completeness the registry does not have** until 2026-10-26:
    "every sample collection a rare-variant study in this atlas draws on". The
    rows for PMID:34324492 named 4 collections, while that study's S2 Table
    lists 189 case sources holding 7,958 cases, of which those 4 hold 4,242
    (counted from the supplement's `s007.xlsx`, sheet "S2 Table": 189 source
    rows summing to its "N" row, 7,958). A reader told the list is complete
    concludes that two studies naming no collection in common share no people,
    which the next sentence must not license. It now says what the list is --
    every collection a row names -- and that a row may not name all its study
    drew on.

    **11 and 6,211 since owner decision R12 (2026-10-26)**, which registered
    the seven CHD case sources S2 identifies by PMID with 50 or more cases.
    Each description opens with S2's own count for its source, and the literal
    below is those counts as read from S2 on 2026-10-27 -- so a count mistyped
    in a description, or two counts swapped between sources, fails here,
    which the sum alone would not catch for a swap.

    Every figure in the sentence that the repository can check is checked:

    * "every sample collection a rare-variant row ... names" describes the list
      both ways: the registry is exactly the ids the mirror's rows cite;
    * the 11 is that study's distinct case cohorts in the mirror;
    * 6,211 is the sum of the counts those 11 descriptions open with, and
      7,958 is the total `taa_cases` states.

    189 is the supplement's alone; nothing committed carries it.
    """
    cited = {cohort for row in mirrored for cohort in (*row.case_cohorts, *row.control_cohorts)}
    assert cited == set(registry)

    named = sorted(
        {cohort for row in mirrored if row.study == "PMID:34324492" for cohort in row.case_cohorts}
    )
    opening = [re.match(r"([\d,]+) ", registry[cohort].description) for cohort in named]
    assert all(opening), f"a description for {named} does not open with its case count"
    counts = {
        cohort: int(match.group(1).replace(",", ""))
        for cohort, match in zip(named, opening, strict=True)
        if match
    }
    # S2 Table's "Number of Individuals" for each source, read from `s007.xlsx`.
    assert counts == {
        "decipher": 1221,
        "isca": 629,
        "pmid18535492_22349727_cases": 63,
        "pmid22318994_cases": 742,
        "pmid22912587_cases": 340,
        "pmid22929023_cases": 123,
        "pmid22939634_cases": 432,
        "pmid22969434_cases": 178,
        "pmid26655555_cases": 91,
        "signature": 1615,
        "taa_cases": 777,
    }
    held = sum(counts.values())
    assert "777 of its 7,958 cases" in registry["taa_cases"].description

    intro = re.search(r'<h2 id="cohorts">Cohorts</h2><p>(.*?)</p>', real_page, re.S)
    assert intro, "the published page has no cohort intro"
    assert (
        f"the rows for PMID:34324492 name {len(named)} of the 189 case sources its S2 Table "
        f"lists, and those {len(named)} hold {held:,} of its 7,958 cases."
    ) in intro.group(1)
    assert (len(named), held) == (11, 6211)


_CELL_TYPE_PARAGRAPHS = (
    "<p><strong>Cell types</strong> come from single-cell RNA sequencing of 8 hearts at 9 to "
    "15 weeks (Farah et al. 2024). For each heart, this atlas sums the counts of every cell "
    "of one type, and of that type in one region, and divides by all of that group's counts: "
    "counts per million (CPM). A page shows the median across hearts and the middle half of "
    "hearts. A heart counts only with 10 or more cells of that type there. Values below 1 CPM "
    "are below this atlas's detection floor, which borrows the Human Protein Atlas's "
    '"not detected" cut of 1. '
    "A dot's size is the share of that type's cells with at least one count.</p>"
    "<p>Regions pool all four ages. τ is computed across cell types as for organs above. "
    "Single-cell methods dissociate tissue, so how many cells of a type were captured "
    "reflects the method, not the heart. A high share of cells describes that cell type, "
    "not how common the type is.</p>"
)


def test_expression_states_what_every_gene_page_shortened_or_folded(page: str) -> None:
    """The expression reading notes, moved off every gene page (2026-10-26, D50).

    Eight paragraphs, pinned whole. The first is the bulk-dilution notice in
    full: a gene page keeps one line of it (`pages._BULK_TISSUE_NOTICE`,
    one of the ten caveats the spec keeps visible) and links here. The next
    two were the page's folded "How to read these figures" notes. The next
    two define tau and the developmental phases, which 118 and 119 gene pages
    print and nothing on the site defined until review (2026-10-26; each
    clause's source is recorded in `methods._EXPRESSION`). The sixth states in
    general the axes each chart's visible caption line states of its own
    figure (D57), and is where the charts' "How to read this chart" toggles
    link. The last two (`_CELL_TYPE_PARAGRAPHS`, at most 250 words) say how a
    cell-type figure is computed, and are where the cell-type block's fold
    links.

    Each was re-read for a page with no gene, dataset or figure on it, and
    reworded where its words assumed one: "below", "shown below", "on this
    page" and "elsewhere on this page" name nothing here. Read against the
    code and a build of 0fc451a (`methods._EXPRESSION` records each):

    * the bulk claim is scoped to the developmental *timeline*, the profile
      datasets whose stage series the gene page's section renders -- the
      cell-type paragraphs below it are single-cell -- and held to every
      profile dataset by
      `test_every_timeline_dataset_is_the_bulk_tissue_the_expression_notes_say`;
    * "the number of genes measured" is the number the grid was built from
      (`Placement.n_genes`: 58,735 in every cell of E-MTAB-6814's grid);
    * "a dataset's own detection floor" read as the source's, and the only
      floor in the corpus is this atlas's curated choice
      (`floor_source: "Curated, not mirrored"`). "its dataset's" read the same
      way, so since review (2026-10-26) the paragraph says whose a floor is
      and where `datasets.json` gives its value and basis, and the bulk
      paragraph says "the detection floor"; the dataset that sentence names is
      held to its record by
      `test_the_dataset_the_floor_sentence_names_declares_no_threshold_of_its_own`;
    * the reference distribution changes shape "over development", not "as
      the heart matures": every organ's grid moves -- liver's 90th percentile
      is 24 tpm at 4 wpc and 5 at adolescence, heart's 28 and 6;
    * the small multiples share one axis *within a dataset*
      (`_small_multiples` fits one per dataset block), and a trajectory follows
      one cardiac organ, its axis fitted to the gene's placed medians there;
    * tau's 0 and 1 are the formula's bounds, and no value between is glossed:
      D39(c) keeps every adjective about a tau off the site;
    * the phases' source is named in the form the trajectory's key renders
      `attributed_to`, and a stage matches a phase on `start <= wpc < end`.

    Pinned whole rather than paragraph by paragraph, for the reason
    `#rare-variants` is: a `findall` over `<p>` skips an appended `<p class=...>`
    or `<details>`.
    """
    assert _section(page, "expression") == (
        '<h2 id="expression">Developmental expression</h2>'
        '<p class="notice-inline"><strong>Every measurement on a developmental timeline in this '
        "atlas is from bulk tissue, never a single cell type.</strong> An embryonic heart "
        "sample mixes cardiomyocytes with endocardium, epicardium, cardiac neural crest, "
        "cushion mesenchyme and blood, and many congenital heart disease genes act in one rare "
        "lineage &mdash; the cardiac neural crest and outflow tract, the second heart field, "
        "the conduction system. A gene expressed intensely in one such lineage can still read "
        "low, or fall below the detection floor, once diluted across the sample: "
        "a low or absent bulk figure is not evidence that a gene is unimportant to heart "
        "development.</p>"
        "<p>A gene's median abundance in an organ is placed against that organ's own "
        "percentile grid for the same dataset and developmental stage, and the number of "
        "genes that grid was built from travels with every percentile a gene page shows. No "
        "percentile is published for a median below its dataset's detection floor, because "
        "the rank of an unreliable measurement is not itself a measurement. A detection floor "
        "is a threshold this atlas records per dataset, not necessarily one its source "
        "declares: where a source declares none, as the Expression Atlas declares none for "
        "E-MTAB-6814, the atlas chooses one, and every developmental dataset's record in "
        '<a href="datasets.json"><code>datasets.json</code></a> gives its floor and where '
        "that floor comes from.</p>"
        "<p>Percentiles in this atlas are <strong>not comparable across organs or across "
        "developmental stages</strong>: each organ transcribes a different fraction of the "
        "whole gene set, and an organ's reference distribution itself changes shape over "
        "development. A high percentile shown for one organ or stage says nothing about a "
        "lower one shown for another.</p>"
        "<p><strong>τ</strong> (tau; Yanai et al. 2005) measures how concentrated a gene's "
        "expression is across the organs a dataset sampled at one stage, computed on "
        "log2(median + 1): it is 0 when every organ's median is equal, and 1 when one organ's "
        "median is above zero and every other's is zero. It says nothing about <em>which</em> "
        "organ, so beside it a gene page names the organ with the highest median, or says the "
        "peak is tied. Every organ's median counts, even one below the detection floor; the "
        "floor decides only whether τ is published at all, and it is not when every organ "
        "sampled is below it. Where a stage has no τ, the page says why in its place, and each "
        "gene's JSON carries the formula and the medians beside every τ it publishes.</p>"
        "<p>A stage's <strong>developmental phases</strong> are the phases of heart "
        "morphogenesis under way at its age, from a vocabulary transcribed from Buijtendijk "
        "MFJ, Barnett P, van den Hoff MJB (2020), PMID:32048790. A phase runs from the start "
        "of one Carnegie stage the review names to the end of another. Each boundary is "
        "converted to weeks post conception from the day ranges the Human Developmental "
        "Stages ontology gives its Carnegie stage, not read off the review's prose, which "
        "counts weeks ordinally: its &ldquo;sixth week&rdquo; is five to six weeks elapsed. "
        "Where the review names no end, as for heart looping, the atlas derives one from the "
        "phase's Gene Ontology definition. A stage names every phase under way at its age, "
        "from the phase's start up to but not including its end, so one stage can name "
        "several; a stage no phase covers names none and says why, for example that it is "
        "post-natal or after the curated window. A chart following one cardiac organ shades "
        "each phase from the first stage "
        "that names it to the last, and its folded key names the phases and their "
        "source.</p>"
        "<p>Every expression chart lays its stages out left to right in curated developmental "
        "order, evenly spaced: the axis is <strong>order, not elapsed time</strong>. Within a "
        "dataset, a gene's small multiples &mdash; one panel per organ with a placed "
        "measurement &mdash; share <strong>one abundance axis</strong>, which is what makes a "
        "gene concentrated in one organ look different from a broadly expressed one; a "
        "per-organ axis would rescale every line to its own range and hide exactly that "
        "difference. A chart following one cardiac organ across a dataset's stages fits its "
        "vertical axis to the gene's own placed medians in that organ, so a curve's height "
        "compares nothing to another gene's.</p>" + _CELL_TYPE_PARAGRAPHS
    )
    added = re.sub(r"<[^>]+>", " ", _CELL_TYPE_PARAGRAPHS).split()
    assert len(added) <= 250, len(added)


def test_every_timeline_dataset_is_the_bulk_tissue_the_expression_notes_say(
    corpus: Corpus,
) -> None:
    """`#expression` opens "Every measurement on a developmental timeline in this
    atlas is from bulk tissue", and a gene page without cell-type data says
    "Every measurement below is from bulk tissue" above its timeline blocks.

    Both are unconditional by design -- a caveat beside only the bad news reads
    as an excuse -- so neither reads the data it describes. This is what does:
    a `profile`-design dataset of any other technology (single-cell,
    proteomics) would put a non-bulk measurement on a timeline and make both
    sentences false, and fails here rather than on a reader. Measured
    2026-10-26: one curated dataset, E-MTAB-6814, `bulk_rnaseq`.

    **Timeline datasets only, since the final review (2026-10-28).** It took
    every curated dataset from review (2026-10-26), while the sentence named
    every "developmental expression measurement", and a contrast dataset of
    fetal tissue is one. The cell-type paragraphs `#expression` now ends with
    made that sentence false, so it was scoped to the timeline, which only a
    `profile` dataset draws: a contrast dataset compares two conditions, and a
    `cell_type_profile` one pools 9 to 15 weeks. Checking every dataset would
    fail the moment the Farah record is curated, over a sentence that record
    does not contradict.
    """
    timeline = [d for d in corpus.datasets if d.design is DatasetDesign.PROFILE]
    assert timeline, "no profile dataset curated; the check would pass vacuously"
    other = sorted(
        f"{dataset.id} ({dataset.technology.value})"
        for dataset in timeline
        if dataset.technology is not Technology.BULK_RNASEQ
    )
    assert not other, (
        "re-read methods.html#expression's 'Every measurement on a developmental timeline in "
        "this atlas is from bulk tissue' and the gene page's 'Every measurement below "
        f"is from bulk tissue' (pages._BULK_TISSUE_NOTICE) against: {other}"
    )


def test_the_dataset_the_floor_sentence_names_declares_no_threshold_of_its_own(
    corpus: Corpus,
) -> None:
    """`#expression` says "the Expression Atlas declares none for E-MTAB-6814,
    the atlas chooses one" (2026-10-26), and that is a fact about one record.

    Unconditional, like the bulk claim, so it reads no data; this does. A
    curator who found the Atlas's own threshold and recorded it would make the
    sentence false, and the record's `floor_source` is where that would show.
    """
    record = next(dataset for dataset in corpus.datasets if dataset.id == "E-MTAB-6814")
    assert record.detection_floor is not None
    assert (record.floor_source or "").startswith(
        "Curated, not mirrored: the Expression Atlas declares no threshold for this experiment."
    ), record.floor_source


def test_the_phase_paragraph_names_the_source_and_the_derived_end_the_curation_records(
    corpus: Corpus,
) -> None:
    """`#expression` names the phases' source as the trajectory's key renders it
    and says the review names every end but heart looping's (2026-10-26).

    Both are facts about `curation/cardiac_phases.yaml`, which this page is
    built without, so this holds the sentence to the file: a new source, or a
    second phase whose end is derived or unstated, makes the paragraph false.
    """
    phases = corpus.cardiac_phases
    assert phases is not None and phases.phases, "no phase vocabulary curated"
    assert (phases.attributed_to, phases.citation) == (
        "Buijtendijk MFJ, Barnett P, van den Hoff MJB (2020)",
        "PMID:32048790",
    )
    assert {
        phase.id: phase.end_basis for phase in phases.phases if phase.end_basis != "stated"
    } == {"heart_looping": "derived"}


@pytest.mark.parametrize("built", ["page", "real_page"])
@pytest.mark.parametrize(
    "retired",
    [
        # Told 24 pages this atlas decides what counts as CHD (2026-08-06).
        "this atlas's CHD scope",
        # False since Costello and CFC entered scope (2026-10-25).
        "No gene here was disputed or refuted",
        # False for 5 genes with an out-of-scope `Disputed` (2026-08-06).
        "not published here at all",
        # There is no table on this page (D50).
        "the table below",
        # A one-submitter gate is not the gate (2026-08-06).
        "one or more Gene Curation",
        # A property of a disease, read as a sufficient condition for scope:
        # FLNA's periventricular nodular heterotopia has it and is not in
        # scope (2026-10-26). Written to catch the scope rule's form, the grade
        # key's and the card's alike.
        "in which Orphanet records",
        "in which <strong>Orphanet records",
        # The gene page's copy notice (`pages._CLINGEN_COPY_NOTICE`): "here" and
        # "the panel record above" name a table this page does not have (D50).
        "listed here only",
        # The refuted mechanism for an absent row (2026-08-05): false for two of
        # the three studies.
        "no variant of that class was seen in either group",
        # Present tense about a second study a one-study page does not have.
        "these cohorts overlap",
        # The tables round to three significant figures (2026-10-26).
        "exactly as published",
        # The card, the matrix and the forest fill apply a rule of the atlas's
        # own, which `#rare-variants` states (owner decision 2026-10-26, W2).
        "adds no assessment",
        # No study here treats a raw p below 0.05 as significant (review,
        # 2026-10-27): the threshold is the atlas's.
        "the studies themselves use",
        # `ddd`'s description published a plan a reader cannot check (owner
        # decision 2026-10-26, W4); `#cohorts` renders it on the real page.
        "plans to curate",
        # A shared collection may be adult controls (2026-10-26).
        "the same children",
        # True of one study, stated in general (2026-10-26).
        "uncorrected scale",
        # `all cases` is not a study's denominator (2026-10-26): PMID:34324492's
        # rows count 6,746 of its 7,958 cases, PMID:40127276's 3,887 trios or an
        # allele total against 11,555 probands.
        "whole case set",
        # The forest key's own deictics: no panel and no study is in sight here,
        # so its clauses are restated in general, never pasted (D57).
        "in this panel",
        "this study's",
        # The expression notes' deictics (2026-10-26): the gene page's bulk
        # notice, its reading notes and its percentile notice, and the band and
        # spark captions' -- no figure, organ or gene is in sight here (D50,
        # D57). "this dataset" is not among them: the inlined stylesheet says it.
        "Every measurement below",
        "shown below",
        "on this page",
        "this organ",
        "for this gene",
    ],
)
def test_no_sentence_the_move_retired_comes_back_here(
    request: pytest.FixtureRequest, built: str, retired: str
) -> None:
    """Whole page, both spellings of an apostrophe (section 4.35).

    A positive literal fails when a sentence is deleted; these fail when a
    deleted one comes back, and coming back -- an old constant pasted onto the
    page it was moved to -- is the failure nobody watches for. Checked on the
    unescaped text too, because the page spells one apostrophe `'` and another
    `&#x27;`, and a retired phrase must not pass by changing which.

    **On the page as published, too** (2026-10-26). It ran only on a page built
    with no cohorts, so it could not see `#cohorts` -- and the published page
    carried "the same children", in `pcgc`'s description, the whole time. No
    retired phrase is exempt from curated text: a curator's sentence on this
    page is held to what this page may say.
    """
    page = request.getfixturevalue(built)
    assert retired not in page
    assert retired not in html.unescape(page)


def test_the_key_follows_the_rule_its_dispute_sentence_points_back_to(page: str) -> None:
    """The rule renders above the key, because the key says "by either path above".

    The key's dispute sentence names no path; it points back at the rule, which
    names both. Rendered above the rule -- or anywhere the rule is not above it
    -- "above" points at nothing, the same defect "the table below" was.
    """
    section = _section(page, "scope")
    rule = section.index('<p class="scope-rule">')
    key = section.index('<div class="grade-key">')
    assert rule < key, "the grade key renders before the rule it says is above it"
