# tests/unit/test_build_pages.py
"""What a gene page must say, especially when the atlas has curated nothing."""

from __future__ import annotations

import html
import os
import re
import subprocess
import sys
from dataclasses import replace
from datetime import date
from enum import StrEnum
from pathlib import Path
from types import EllipsisType
from typing import Final

import pytest

from chd_atlas.build.burden import BurdenRow, CombinedRow
from chd_atlas.build.cell_types import CellTypeExpression
from chd_atlas.build.concordance import cohort_families, evidence_axes, gene_concordance
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.identity import GeneIdentity
from chd_atlas.build.pages import (
    _EM_DASH,
    _analysis_label,
    _chd_only,
    _count,
    _date,
    _forest_text,
    _glance_line,
    _inheritance,
    _phase_strip,
    _study_finding_summary,
    _tissue_percentile,
    build_gene_index_page,
    build_gene_pages,
)
from chd_atlas.build.paths import gene_page_path
from chd_atlas.build.profiles import (
    EMPTY_EXPRESSION_PROFILE,
    DatasetProfileEntry,
    ExpressionProfile,
    PhaseInfo,
    Placement,
    Specificity,
    StageProfileEntry,
    TissueProfileEntry,
)
from chd_atlas.build.render import Markup
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.models.assertion import Evidence, InTextLocator, LesionAssertion
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.dataset import Dataset, Stage
from chd_atlas.models.literature import Publication
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile, EndBasis
from chd_atlas.models.review import ReviewStatus
from chd_atlas.models.scope import ScopeEntry
from chd_atlas.orphanet_scope import Annotation, Warrant
from chd_atlas.vocab import (
    AtlasCuration,
    Classification,
    EvidenceClass,
    EvidenceStrength,
    LesionGroup,
    ValiditySource,
    ValidityState,
)
from tests.unit.test_build_cell_type_chart import dataset as cell_type_dataset
from tests.unit.test_build_cell_type_chart import payload as cell_type_payload
from tests.unit.test_build_cell_type_chart import specificity as cell_type_specificity

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"

# The rare-variant section's heading, as published (D49). One literal for every
# slice below, so a rename edits one line and a vanished heading raises rather
# than letting a slice fall through to the wrong section.
_RARE_VARIANT_H2: Final = '<h2 id="rare-variants">Rare-variant association</h2>'
_RARE_VARIANT_HEADING: Final = "Rare-variant association"


def _validity(
    disease_label: str = "Holt-Oram syndrome",
    report_url: str | None = "https://search.clinicalgenome.org/kb/gene-validity/x",
) -> GeneValidity:
    return GeneValidity(
        records=(
            ValidityRecord(
                source=ValiditySource.CLINGEN,
                classification=Classification.DEFINITIVE,
                classification_term="Definitive",
                disease="MONDO:0007732",
                disease_label=disease_label,
                moi="AD",
                sop="SOP11",
                gcep="Syndromic Disorders Gene Curation Expert Panel",
                report_url=report_url,
            ),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )


def _facts(
    gene: str,
    curation: AtlasCuration,
    groups: tuple[LesionGroup, ...] = (),
    functional_count: int = 0,
    no_association_reported_by: tuple[str, ...] = (),
) -> GeneFacts:
    """A published gene's facts. The third axis defaults off and is opt-in.

    `no_association_reported_by` carries both halves of the pair: the flag is
    derived from it here rather than passed separately, so a fixture cannot
    construct the impossible state -- an authority named while the flag reads
    false -- that `derive.gene_facts` refuses to produce.
    """
    return GeneFacts(
        gene=gene,
        headline_confidence=Classification.DEFINITIVE,
        validity_state=ValidityState.EXPERT_CURATED,
        has_conflicting_evidence=False,
        has_source_discordance=False,
        has_no_association_report=bool(no_association_reported_by),
        no_association_reported_by=no_association_reported_by,
        lesion_groups=groups,
        confidence_by_lesion_group={group: Classification.DEFINITIVE for group in groups},
        conflicting_lesion_groups=(),
        evidence_counts={EvidenceClass.GENETIC_CASE: 1} if groups else {},
        assertion_count=1 if curation is AtlasCuration.CURATED else 0,
        functional_count=functional_count,
        publications=("PMID:8988165",) if curation is AtlasCuration.CURATED else (),
        atlas_curation=curation,
    )


@pytest.fixture
def facts_uncurated() -> dict[str, GeneFacts]:
    return {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED)}


@pytest.fixture
def facts_curated() -> dict[str, GeneFacts]:
    return {TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,))}


@pytest.fixture
def facts_two() -> dict[str, GeneFacts]:
    return {
        TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,)),
        GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
    }


@pytest.fixture
def validity_two() -> dict[str, GeneValidity]:
    """Both genes of `facts_two`, each Definitive for a *different* in-scope disease.

    Different labels on purpose. The `definitive for` column and the rail row it
    mirrors are both per-gene lookups, and a fixture where every gene carried
    the same disease label could not tell a correct lookup from one wired to the
    first gene, to `next(iter(validity))`, or to a constant.
    """
    return {TBX5: _validity(), GATA4: _validity("structural congenital heart disease")}


@pytest.fixture
def assertions() -> dict[str, list[LesionAssertion]]:
    """One syndromic assertion, which is what `LesionAssertion` will accept.

    `syndromic_lists_extracardiac_features` refuses a `syndromic` record that
    names no extracardiac feature, so `extracardiac_features` is not optional
    here -- HP:0001156 is brachydactyly, the limb half of Holt-Oram.
    """
    return {
        TBX5: [
            LesionAssertion(
                id="CHDA:AST:0000001",
                gene=TBX5,
                phenotypes=["HP:0001631"],
                lesion_groups=[LesionGroup.SEPTAL],
                inheritance=["AD"],
                mechanism="haploinsufficiency",
                syndromic="syndromic",
                extracardiac_features=["HP:0001156"],
                curator="Enrique Audain",
                curated_on=date(2026, 7, 31),
                last_reviewed=date(2026, 7, 31),
                evidence=[
                    Evidence(
                        publication="PMID:8988165",
                        evidence_class=EvidenceClass.GENETIC_CASE,
                        locator=InTextLocator(element="Abstract"),
                        summary="A nonsense TBX5 mutation was found in affected members.",
                        strength=EvidenceStrength.STRONG,
                    )
                ],
            )
        ]
    }


@pytest.fixture
def publications() -> dict[str, Publication]:
    return {
        "PMID:8988165": Publication(
            id="PMID:8988165",
            title="Mutations in human TBX5 cause limb and cardiac malformation.",
            journal="Nature genetics",
            year=1997,
            authors=["Basson CT"],
            study_type="family_linkage",
            own_lab=False,
        )
    }


def _page(root: Path, name: str) -> str:
    return (root / "genes" / name).read_text(encoding="utf-8")


def _validity_table(page: str) -> str:
    """The mirrored-validity section alone: both authority groups.

    Ends at the next `<h2>` after the start. It also ended at the uncurated
    notice (`<p>The atlas has`) until D55 moved that notice from under this
    section into the card. The Orphanet notes are inside the section, above the
    first group, so they are not an end marker.

    Scoped rather than checked page-wide because the rail's "this gene as JSON"
    link and the shell's `<nav>` are both `<a href=`, so "no anchor anywhere"
    could never pass on any page. What the report column must not contain is one.
    Since D51 the section holds two tables, so it ends at the next heading of its
    own rank rather than at the first `</table>`.
    """
    start = page.index("<h2>Mirrored gene&ndash;disease validity</h2>")
    end = page.find("<h2>", start + 1)
    return page[start : end if end != -1 else len(page)]


def test_an_uncurated_gene_page_says_the_atlas_has_not_curated_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """118 of the 119 genes published today (measured 2026-10-26). The page must state it.

    An absent section is indistinguishable from "the atlas looked and found
    nothing", and a reader deciding what a gene means clinically must not have
    to infer which. The mirrored validity table is still present -- that is the
    content the page exists to show.

    **Re-pointed at the card (D55).** This asserted the main-column paragraph
    "The atlas has not yet curated a lesion assertion for this gene" until
    2026-10-26, when the card's curation sentence replaced it; the sentence is
    now asserted in the card's own element, and the old paragraph is asserted
    *absent*, so it cannot come back beside the card as a second copy.
    `test_the_card_says_in_plain_words_what_the_atlas_itself_curated` pins all
    four branches of the sentence; this one pins that a whole uncurated page
    carries it once, with the validity table it qualifies.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert re.findall(r'<p class="card-curation">(.*?)</p>', page) == [
        "Not yet curated by the atlas: every classification on this page is a ClinGen expert "
        "panel's or a named curator's."
    ]
    assert page.index('<p class="card-curation">') < page.index("</aside>")
    assert "The atlas has <strong>not yet curated</strong> a lesion assertion" not in page
    assert "#curated-evidence" not in page
    assert "Holt-Oram syndrome" in _validity_table(page)
    assert "GATA4" in page


def test_a_curated_gene_page_carries_its_evidence_quote_and_pmid(
    tmp_path: Path,
    facts_curated: dict[str, GeneFacts],
    assertions: dict[str, list[LesionAssertion]],
    publications: dict[str, Publication],
) -> None:
    """The evidence a curator wrote, on a page, plus the way back to the JSON.

    The bundle link is asserted here rather than in its own test because it is
    the same interpolation the escaping tests already cover; what is pinned is
    the shape -- one `../` out of `genes/` and back in, since the page and its
    bundle are siblings and `paths.py` names both from the site root.

    The PubMed href is pinned for a stronger reason: nothing observed the URL
    `_pubmed` constructs. Asserting `"PMID:8988165" in page` passes on the link
    *text*, which is the identifier verbatim, so the `href` beside it was
    unmeasured -- dropping `.removeprefix("PMID:")` publishes
    `https://pubmed.ncbi.nlm.nih.gov/PMID:8988165/`, a 404 on every evidence row
    of every curated gene page, and the whole suite stayed green (measured
    2026-08-04: 616 passed). With the href asserted below, that mutant fails
    here and nowhere else: 1 failed, 621 passed.

    Asserted as the full attribute rather than as the bare number, because
    "8988165" is a substring of "PMID:8988165" and a substring check on the
    digits alone cannot fail.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_curated,
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: _validity()},
        assertions=assertions,
        publications=publications,
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_11604.html")
    assert "PMID:8988165" in page
    assert 'href="https://pubmed.ncbi.nlm.nih.gov/8988165/"' in page
    assert "A nonsense TBX5 mutation was found in affected members." in page
    assert "Mutations in human TBX5 cause limb and cardiac malformation." in page
    # Scoped to the lesion-assertion axis rather than the bare phrase: this
    # gene carries no profiles data either (the call above passes no
    # `profiles=`), so its page also carries the unrelated "not yet curated a
    # developmental expression profile" notice -- a second, independent axis
    # legitimately reusing the template the lesion-assertion notice
    # established. A bare page-wide check cannot tell the two apart; CLAUDE.md
    # section 4.19 is this exact failure shape. The old main-column notice is
    # refused, and so are the card's two uncurated sentences (D55): a curated
    # gene's card says what was curated instead.
    assert "not yet curated</strong> a lesion assertion" not in page
    assert "Not yet curated by the atlas" not in page
    assert "No lesion assertion yet curated by the atlas" not in page
    assert "Not curated by the atlas" not in page
    assert "No lesion assertion curated by the atlas" not in page
    assert 'href="../genes/HGNC_11604.json"' in page
    # The section links its reading notes, in the pattern of every other
    # section's intro (owner decision 2026-10-26, E3); sliced to the section.
    section = page[page.index("<h2>Curated evidence</h2>") : page.index("<h2>Publications</h2>")]
    assert section.startswith(
        "<h2>Curated evidence</h2><p>Each lesion assertion below is this atlas's own "
        "curation, and changes no classification above. "
        '<a href="../methods.html#curated-evidence">How to read curated evidence</a>.</p>'
    )
    assert page.index("<h2>Mirrored gene&ndash;disease validity</h2>") < page.index(
        "<h2>Curated evidence</h2>"
    )


def test_a_symbol_carrying_markup_is_escaped(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Nothing curated or mirrored is reviewed the way page copy is.

    The symbol reaches the page twice -- the `<h1>` in the rail and the
    `<title>` through `document` -- and only the first is escaped by this
    module. Dropping that one `html.escape` fails the first assertion here
    while the second still passes on the title's own escaping, which is why
    both are asserted.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "<script>x</script>"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert "<script>x</script>" not in page
    assert "&lt;script&gt;" in page


@pytest.mark.parametrize(
    "report_url",
    ["Pseudoautosomal region, recessive", "javascript:alert(1)", "  https://evil.example"],
)
def test_a_non_http_report_url_never_becomes_a_link(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts], report_url: str
) -> None:
    """Escaping is not sanitising, and neither mirror schema constrains this column.

    The first case is real, not invented: `mirrors/gencc_submissions.tsv` carries
    one row -- SHOX, `MONDO:0009588`, Ambry Genetics -- whose `report_url` is the
    literal string `Pseudoautosomal region, recessive` (measured 2026-08-04, the
    only such row of the file). Langer mesomelic dysplasia is out of CHD scope
    today, so that row reaches no page and the defect is latent; the second case
    is the same column carrying a scheme that executes. `render.py` would escape
    either into a syntactically valid `href` and publish it.

    Both assertions matter. The first is the guard: no anchor at all for this
    row. The second is what stops a "fix" that merely escapes harder -- the raw
    string must not appear inside any attribute, quoted or not, because
    `href="javascript:alert(1)"` needs no unescaped character to run.

    The third case pins that leading whitespace does not smuggle a scheme past
    `str.startswith`; browsers strip it, `startswith` does not, so the URL never
    becomes a link there either. In every case the fixture's SOP renders as plain
    text in the report cell (`<td>SOP11</td>`): the SOP is still stated, and only
    the unusable URL is dropped.

    Measured 2026-08-04 with the guard removed (`Link(text="open", href=url) if
    url else _EM_DASH`): 3 failed, 619 passed -- one failure per case here, and
    no other test in the suite noticed.

    Since 2026-10-26 (D50) the section's intro links `methods.html#validity`,
    so "no anchor in the section" can no longer pass; the check reads the two
    authority groups, from the first `<h3>` on, which is where every row is.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity(report_url=report_url)},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    section = _validity_table(page)
    groups = section[section.index("<table") :]
    assert "<a href=" not in groups
    assert "<td>SOP11</td>" in groups
    assert report_url.strip() not in page


def test_panelapp_australias_homepage_link_opens_the_genes_own_page_and_nothing_else_moves(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """R8 (owner decision 2026-10-26): the one report URL the atlas builds.

    GenCC gives `https://panelapp-aus.org`, the bare homepage, as the report URL
    of every PanelApp Australia submission (5,478 of 5,478 in the mirror), so
    the row's link opened a page naming no gene. Exactly that URL from exactly
    that submitter links `/panels/entities/<symbol>` instead, the symbol
    percent-encoded (`@`, which HGNC's older group symbols carry, is the
    character that shows it). A PanelApp Australia URL that is not the bare
    homepage, and another submitter carrying that same homepage string, keep
    the URL GenCC gave. Rows are in the GenCC group, in mirror order.
    """
    homepage = "https://panelapp-aus.org"
    records = tuple(
        replace(_gencc(submitter, Classification.DEFINITIVE), report_url=url)
        for submitter, url in (
            ("PanelApp Australia", homepage),
            ("PanelApp Australia", "https://panelapp-aus.org/panels/250/"),
            ("G2P", homepage),
        )
    )
    build_gene_pages(
        facts_uncurated,
        Emitter(root=tmp_path),
        symbols={GATA4: "GATA4@"},
        validity={
            GATA4: GeneValidity(
                records=records,
                state=ValidityState.EXPERT_CURATED,
                has_source_discordance=False,
            )
        },
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    groups = _validity_table(_page(tmp_path, "HGNC_4173.html"))
    assert re.findall(r'<a href="(https://panelapp-aus[^"]*)">', groups) == [
        "https://panelapp-aus.org/panels/entities/GATA4%40",
        "https://panelapp-aus.org/panels/250/",
        "https://panelapp-aus.org",
    ]


def test_an_uncurated_notice_never_denies_functional_records_the_rail_counts(
    tmp_path: Path,
) -> None:
    """The card must not deny the atlas's own work in the sentence that counts it.

    `atlas_curation` is `curated` iff `assertion_count > 0` (`derive.gene_facts`),
    but functional records come from a population that needs no assertion at
    all, so a gene can carry curated functional records and still be
    `NOT_YET_CURATED`. The main-column notice used to answer that with "not yet
    curated evidence for this gene ... nothing on this page is the atlas's own
    assessment", one column away from a rail row reading `functional records |
    3`.

    **Re-pointed at the card (D55), 2026-10-26.** The rail no longer has a
    `functional records` row and the notice is now the card's curation
    sentence, so the count and the denial would sit in one sentence: "Not
    curated by the atlas ... The atlas has curated 3 functional records." What
    is pinned is the wording that avoids it -- the sentence names the *lesion
    assertion* as what is absent and counts the records -- and the denials in
    every form they could come back in, the new one included.

    `atlas_curation`'s meaning is deliberately unchanged -- it is tested,
    documented and published. Latent: the committed corpus has no
    `curation/functional/` directory, so every published gene counts zero
    functional records (measured 2026-08-04). Reachable only from this fixture.
    """
    facts = {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED, functional_count=3)}
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    curation = _slice_between(_rail_of(page), '<p class="card-curation">', ("</p>",))
    # The card counts them, and names the lesion assertion as what is absent.
    assert curation.startswith("No lesion assertion yet curated by the atlas: ")
    assert curation.endswith(
        "The atlas has curated 3 functional records, published in this gene's JSON."
    )
    # The sentences that would deny it, in the forms they could come back in:
    # the card's own unqualified openings, and the two the old notice used.
    assert "Not yet curated by the atlas" not in page
    assert "Not curated by the atlas" not in page
    assert "nothing on this page is the atlas's own assessment" not in page
    assert "not yet curated</strong> evidence for this gene" not in page


def test_a_mirrored_disease_label_carrying_markup_is_escaped(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The validity table is filled from a mirror, not from curated text."""
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity(disease_label="<img src=x onerror=alert(1)>")},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    assert "<img src=x" not in _page(tmp_path, "HGNC_4173.html")


# The `sources disagree` chip as published: its title is the only place a
# reader is told what the two words measure.
_SOURCES_DISAGREE_CHIP: Final = (
    '<span class="chip chip-warn" title="One source, ClinGen&#x27;s expert panels or the '
    "submissions to the Gene Curation Coalition, disputes or refutes this gene for a disease "
    "in scope while the other classifies it Limited or better for one, not necessarily the "
    'same disease; see the validity table.">sources disagree</span>'
)
_CONFLICTING_EVIDENCE_CHIP: Final = '<span class="chip chip-warn">conflicting evidence</span>'


@pytest.mark.parametrize(
    ("conflicting", "discordant", "shown", "hidden"),
    [
        (True, False, _CONFLICTING_EVIDENCE_CHIP, "sources disagree"),
        (False, True, _SOURCES_DISAGREE_CHIP, "conflicting evidence"),
    ],
)
def test_a_contested_gene_is_never_chipped_as_settled(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    conflicting: bool,
    discordant: bool,
    shown: str,
    hidden: str,
) -> None:
    """The green pill tracks the classification, not merely the fact of one.

        `chip-definitive` is the top of the stylesheet's grade ramp and the only
        class filled with the full success colour, so keying it on
        `headline_confidence is not None` would paint a refuted gene green while its
        own label read "refuted" -- the display failure `vocab.strongest` is
        documented against. Since 2026-08-06 the class comes from
        `render.GRADE_CHIP_KIND`, keyed on the grade's own value, and `refuted` maps
        to `warn`: the assertion below is that the rail's chip is the caution colour,
        not merely that it is not green.

        **The slice is the rail, not the page.** `render.grade_legend` rendered one
        chip of every class -- including `chip-definitive` -- inside a `<details>` on
        every page, so a page-wide `not in` could never fail again and passed only
        because the legend did not exist when it was written. CLAUDE.md section 4.19
        is this exact shape: a page-wide assertion is not a section assertion. The
        legend left the gene page for `methods.html#scope` on 2026-10-26 (D50); the
        slice stays, because the rail is still the element under test.

        Parametrised over the two warning flags in opposite states because they are
        separate axes: `has_conflicting_evidence` is any mixed evidence,
        `has_source_discordance` is specifically ClinGen and GenCC disagreeing with
        each other. Each case is the unique killer of dropping one of the two chips
        -- measured 2026-08-04 by replacing each `if` with `if False:` in `_rail`:
        dropping "sources disagree" failed only the second case, dropping
        "conflicting evidence" only the first, and nothing else in this file
        noticed either. The `hidden` assertion is what catches the two flags being
        read into the wrong chip.

    No published gene carries a contested headline -- the floor excludes those
        rungs and the ClinGen veto excludes the rest -- so the refuted chip is
        latent, reachable only from a fixture. The two flags are not: both are set
        on six published genes since the Orphanet scope path (BRAF, KRAS, MAP2K1,
        PTPN11, RAF1, SOS1; this read "neither flag is set on any of the 92
        genes" until 2026-10-26).

        **Each chip is asserted whole, title included (2026-10-26).** `sources
        disagree` had no title and `methods.html` never mentioned it, while on
        all six pages carrying it ClinGen and GenCC agree on every gene-disease
        pair they both classify: the flag is per gene, and fires on a dispute for
        one disease beside support for another. The title says so.

        **The gene index carries the same `conflicting evidence` chip, keyed on
        the same flag.** Its confidence column showed those six genes as a bare
        `definitive` until 2026-10-26, the settled reading `docs/data-api.md`
        forbids a client. Asserted on the gene's own row, beside the bare grade
        the confidence facet still matches on; the second case is the row that
        must not carry it.
    """
    facts = {
        GATA4: replace(
            facts_uncurated[GATA4],
            headline_confidence=Classification.REFUTED,
            has_conflicting_evidence=conflicting,
            has_source_discordance=discordant,
        )
    }
    # A flagged gene carries the contested record its flag was derived from:
    # the card names it, and refuses a flag with no record behind it.
    validity = _validity()
    if conflicting:
        refuted = replace(
            validity.records[0],
            classification=Classification.REFUTED,
            classification_term="Refuted",
            disease="MONDO:0009026",
            disease_label="Costello syndrome",
        )
        validity = replace(validity, records=(*validity.records, refuted))
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None, "the gene page no longer renders a rail"
    # The rendered chip, not the bare class name: `.chip-definitive` is also a
    # rule in the stylesheet every page inlines, so a substring check on the
    # class alone can never fail.
    assert 'class="chip chip-definitive"' not in rail.group()
    assert 'class="chip chip-warn">refuted</span>' in rail.group()
    assert shown in rail.group()
    assert hidden not in rail.group()

    build_gene_index_page(
        facts, emitter, symbols={GATA4: "GATA4"}, validity={GATA4: validity}, burden_counts={}
    )
    browse = _page(tmp_path, "index.html")
    row = re.search(r"<tr[^>]*>(?:(?!</tr>).)*HGNC_4173\.html(?:(?!</tr>).)*</tr>", browse, re.S)
    assert row is not None, "the gene index has no row for the gene"
    assert 'data-confidence="refuted"' in row.group()
    cells = re.findall(r"<td>(.*?)</td>", row.group(), re.S)
    assert cells[2] == (f"refuted {_CONFLICTING_EVIDENCE_CHIP}" if conflicting else "refuted")

    # **A contested gene no panel graded keeps its chip, beside the em dash**
    # (review, 2026-10-26). A GenCC dispute is no veto, so a gene admitted on
    # submitter agreement can carry the flag with a null headline; no published
    # gene does today, and suppressing the chip wherever the headline is null
    # survived the full suite. `data-confidence` stays empty: the facet has no
    # grade to match.
    ungraded = tmp_path / "ungraded"
    build_gene_index_page(
        {GATA4: replace(facts[GATA4], headline_confidence=None)},
        Emitter(root=ungraded),
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        burden_counts={},
    )
    row = re.search(
        r"<tr[^>]*>(?:(?!</tr>).)*HGNC_4173\.html(?:(?!</tr>).)*</tr>",
        _page(ungraded, "index.html"),
        re.S,
    )
    assert row is not None, "the gene index has no row for the ungraded gene"
    assert 'data-confidence=""' in row.group()
    cells = re.findall(r"<td>(.*?)</td>", row.group(), re.S)
    assert cells[2] == (f"— {_CONFLICTING_EVIDENCE_CHIP}" if conflicting else "—")


def test_every_gene_in_the_facts_gets_exactly_one_page(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Including a gene with no symbol in the registry -- `symbols` is empty here."""
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_two,
        emitter,
        symbols={},
        validity=validity_two,
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    assert sorted(path.name for path in (tmp_path / "genes").iterdir()) == [
        "HGNC_11604.html",
        "HGNC_4173.html",
    ]


def test_the_browse_table_is_complete_before_any_script_runs(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """D29. The rows are rendered by the build, not fetched.

    The mutant this kills is a page that ships an empty tbody and populates it
    from genes/index.json at runtime: it looks identical in a browser and serves
    nothing to curl, to a crawler, or to a reader with scripts disabled.

    Measured 2026-08-04 by rendering `<tbody>` empty and appending the rows from
    an inline `fetch('index.json')` instead: 3 failed, 609 passed. This test
    failed first, on `"TBX5" in page` -- the symbol never reaches the markup at
    all, so the later `fetch(` and `<script>`-ordering assertions are not what
    catches this particular mutant; they catch the narrower one that renders the
    rows *and* then re-fetches. The other two failures were the two browse tests
    below, which is expected: an empty table has no order and no links.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert "TBX5" in page and "GATA4" in page
    assert "fetch(" not in page
    assert page.index("HGNC_11604.html") < page.index("<script>")


def test_browse_rows_are_ordered_by_hgnc_id_against_a_literal(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Asserted against a literal, not by building twice and comparing.

    PYTHONHASHSEED is fixed for the life of an interpreter, so a same-process
    build-twice comparison cannot catch a dropped sort at any fixture size.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert page.index("HGNC_11604.html") < page.index("HGNC_4173.html")


def test_browse_facet_options_are_ordered_against_a_literal(
    tmp_path: Path, validity_two: dict[str, GeneValidity]
) -> None:
    """The `<option>` values are built from sets, so their order is seed-dependent.

    `build_gene_index_page` derives each facet's values with a set comprehension
    over `facts.values()`. Set iteration order for strings varies with
    `PYTHONHASHSEED`, and these options are part of the page's bytes and
    therefore of its checksum and its manifest entry -- so dropping the
    `sorted(values)` makes two builds of one commit differ. Measured 2026-08-04
    on the committed corpus, one process per seed: with the sort dropped,
    `genes/index.html` and `manifest.json` differ between `PYTHONHASHSEED=0` and
    `PYTHONHASHSEED=2`, and no other file does. The whole suite still passed --
    616 passed -- so nothing guarded this. With this test in place the mutant
    fails here and nowhere else: 1 failed, 621 passed.

    Asserted against a literal rather than by building twice: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so a same-process comparison cannot see
    this at any fixture size (CLAUDE.md §4.13).

    The fixture names **all nine** `LesionGroup` members on one gene, because
    this guard is probabilistic in the same way `test_cited_publications_come_back_sorted`
    is: *n* strings can land in an order that happens to equal sorted order, and
    the unsorted code then passes. Measured with the sort dropped, one process
    per seed, counting seeds on which the mutant survives:

    * the two-value `curation` facet the committed corpus would give
      (`curated`, `not_yet_curated`) -- **105/200**. Sizing this test to the real
      corpus would have produced a guard that waved the regression through on
      better than half of all runs, which §4.12 records as not being a guard.
    * the nine-member `lesion` facet below -- **0/200**, and **0/1000**.

    Nothing pins `PYTHONHASHSEED` in CI or in a `conftest.py`, so that difference
    is the whole value of the fixture. The `curation` assertion is kept beneath
    the `lesion` one anyway: it is the facet a reader actually uses to tell a
    curated gene from a mirrored one, and pinning its two values costs nothing.
    """
    facts = {
        TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=tuple(LesionGroup)),
        GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
    }
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    facets = dict(
        re.findall(r'<select name="([^"]+)"[^>]*>(.*?)</select>', _page(tmp_path, "index.html"))
    )
    assert re.findall(r'<option value="([^"]*)">', facets["lesion"]) == [
        "",
        "anomalous_venous_return",
        "avsd",
        "complex",
        "conotruncal",
        "heterotaxy",
        "lvoto",
        "other",
        "rvoto",
        "septal",
    ]
    assert re.findall(r'<option value="([^"]*)">', facets["curation"]) == [
        "",
        "curated",
        "not_yet_curated",
    ]


def test_the_browse_page_says_whose_classification_the_confidence_column_carries(
    tmp_path: Path,
) -> None:
    """23 rows reading `definitive` under a column headed `confidence`, unattributed.

    Every gene page says "Every classification below is an upstream panel's or
    submitter's" and the landing page's hero note says every gene is admitted
    on an upstream assessment, "never this atlas's own". This was the only page on the site where a
    panel's call could be read as the atlas's own, and the `atlas curation`
    column that answers the question was unexplained beside it.

    Both halves are asserted because they are separate claims: the attribution,
    and the pointer to the column that says whether the atlas has curated the
    gene at all. `_BROWSE_HEADERS` is checked in the same breath so the sentence
    cannot come to name a column the table stopped rendering.

    **The attribution itself was wrong, and is pinned as a literal since
    2026-10-26.** The paragraph said "Every confidence below is an upstream
    panel's or submitter's" and "A blank confidence means ...". The column is
    `headline_confidence`, the admitting ClinGen panel's grade and never a
    submitter's, and a gene no panel graded renders an em dash there, not a
    blank -- 21 of the 119 rows on a real build. So the fixture now carries one
    gene of each kind, and the em-dash cell is asserted beside the sentence
    that describes it: the two fail together.
    """
    graded, graded_validity = _limited(TBX5)
    admitted, admitted_validity = _submitter_admitted(GATA4)
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        {TBX5: graded, GATA4: admitted},
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity={TBX5: graded_validity, GATA4: admitted_validity},
        burden_counts={},
    )
    page = _page(tmp_path, "index.html")

    intro = _slice_between(page, "<p>Every confidence below", ("</p>",))
    assert intro == (
        " is a ClinGen expert panel's, republished with its provenance intact "
        "&mdash; the atlas authors no validity classification of its own. The "
        "<strong>atlas curation</strong> column is what says whether the atlas has "
        "curated a gene itself, and the <strong>graded for</strong> column names the "
        "disease each panel actually graded. A confidence shown as &mdash; means no "
        "ClinGen expert panel has graded the gene for any of the diseases in scope "
        "here; its page names the submitters whose agreement admitted it."
    )
    assert '<th scope="col">confidence</th>' in page
    row = re.search(r"<tr[^>]*>(?:(?!</tr>).)*HGNC_4173\.html(?:(?!</tr>).)*</tr>", page, re.S)
    assert row, "the submitter-admitted gene has no browse row"
    assert re.findall(r"<td>(.*?)</td>", row.group(0), re.S)[2] == "—"


def test_no_browse_row_states_a_bare_definitive_without_the_disease_it_is_for(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """`KMT2D | definitive` on a site called "CHD Atlas" is a wrong claim.

    The column is headed `graded for` since 2026-08-06, not `definitive for`: the
    widened gate publishes 43 `limited`, 9 `moderate` and 1 `strong` gene beside
    the 23 `definitive` ones, and a column head naming one rung above a cell
    naming another is a header that contradicts its own data.

    ClinGen's assertion is `KMT2D -- Definitive for Kabuki syndrome 1`, made by
    the SCID-CID GCEP. Measured 2026-08-04 against the committed mirrors, the
    browse page rendered 23 rows all reading `definitive`, with no disease column
    at all: 10 of the 23 are qualified by a panel that is not the Congenital
    Heart Disease GCEP, and 10 are definitive for a disease whose label names no
    cardiac feature (Kabuki 1 and 2, Mowat-Wilson, CHARGE, Alagille, Ellis-van
    Creveld twice, Holt-Oram, TARP, NR2F2-related multiple congenital anomalies).

    Asserted per row rather than as "both labels appear somewhere on the page",
    because the failure this guards is a lookup wired to the wrong gene: the
    fixture gives the two genes two different diseases, and a page naming
    Holt-Oram on both rows would satisfy a page-wide substring check while
    telling a reader GATA4 is definitive for Holt-Oram syndrome.

    The disease cell is required to sit in the same `<tr>` as the `definitive`
    cell, so moving the column to the far right -- or into a second table -- fails
    here even though every string is still on the page.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert '<th scope="col">graded for</th>' in page
    rows = re.findall(r"<tr(?: data-[^>]*)?><td>(.*?)</tr>", page)
    assert len(rows) == 2
    expected = {TBX5: "Holt-Oram syndrome", GATA4: "structural congenital heart disease"}
    for row in rows:
        gene = next(g for g in expected if g in row)
        assert f"<td>definitive</td><td>{expected[gene]}</td>" in row, row
        # The other gene's disease must not be on this row at all.
        other = next(expected[g] for g in expected if g != gene)
        assert other not in row


def test_a_gene_definitive_for_two_in_scope_diseases_names_both_in_a_fixed_order(
    tmp_path: Path,
) -> None:
    """Latent today, and sorted so it stays deterministic when it is not.

    Measured 2026-08-04 against the committed mirrors: all 23 published genes
    carry exactly one in-scope ClinGen `Definitive` record naming exactly one
    distinct disease label, so nothing in the real corpus exercises this. It is a
    fixture, not a hypothetical: `curation/chd_scope.yaml` is an editorial file a
    curator widens, and one added MONDO term can give a gene a second qualifying
    record without any code changing.

    Two properties, and they are separate. That **both** labels are named --
    dropping one would publish half the reason the gene is on the site. And that
    they come back in sorted order -- `_graded_diseases` de-duplicates through
    a `set`, whose iteration order for strings varies with `PYTHONHASHSEED`, so
    an unsorted return makes `genes/index.html` and its manifest checksum differ
    between two builds of one commit. Asserted against a literal in reverse
    alphabetical fixture order rather than by building twice: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so a same-process comparison cannot see
    a dropped sort at any fixture size (CLAUDE.md §4.13).

    Only a GenCC record calling the gene definitive is also present, and it is
    asserted absent: `_graded_diseases` names the diseases the *admitting panel*
    graded, and this gene was admitted on a ClinGen record, so naming a
    submitter's disease here would qualify the confidence with a disease that is
    not the one the gate turned on.

    That rationale read "GenCC admits no gene to this population (D21)" until
    2026-08-06, which the widening made false -- GenCC agreement now admits 16 of
    the 92. The assertion it justifies is unchanged and still correct; only the
    reason was stale. A GenCC-admitted gene has no headline grade at all, so
    `_graded_diseases` returns `()` for it and never reaches this branch.
    """
    two = GeneValidity(
        records=(
            _validity("Zellweger syndrome").records[0],
            _validity("Alagille syndrome").records[0],
            ValidityRecord(
                source=ValiditySource.GENCC,
                classification=Classification.DEFINITIVE,
                classification_term="Definitive",
                disease="MONDO:0000001",
                disease_label="a submitter's disease",
                moi="AD",
                report_url=None,
                submitter="Ambry Genetics",
            ),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )
    facts = {TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,))}
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts, emitter, symbols={TBX5: "TBX5"}, validity={TBX5: two}, burden_counts={}
    )
    build_gene_pages(
        facts,
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: two},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    browse = _page(tmp_path, "index.html")
    assert "<td>definitive</td><td>Alagille syndrome; Zellweger syndrome</td>" in browse
    assert "a submitter's disease" not in browse
    gene_page = _page(tmp_path, "HGNC_11604.html")
    assert "<dt>definitive for</dt><dd>Alagille syndrome; Zellweger syndrome</dd>" in gene_page


def test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Both page kinds say whose classification they show, and link the rule in full.

    No page on the site said what the published set is until 2026-08-04;
    `docs/data-api.md` did, and a reader of the HTML never sees that document.
    Until 2026-10-26 both page kinds carried the whole admission rule and a
    folded grade key. Both moved to `methods.html#scope` (D50), where
    `test_build_methods.py` pins them; each page kind now carries a one-line
    short form and a link to that section.

    What this asserts, on the browse page and on two gene pages: the short form
    as a literal, sliced to its `<p class="scope-rule">` element and present
    exactly once; the link inside it; that the full rule has not come back
    beside it (a page carrying both is the page D50 exists to slim); and the
    retired "this atlas's CHD scope" absent in both spellings of its
    apostrophe.

    **The short form must be true on a page whose panel graded CHD itself.**
    The spec's first wording -- "for the disease named, which is not the same
    as for CHD" -- contradicted the card on 53 of the 119 gene pages, whose
    panel graded congenital heart disease (ISL1: "definitive for congenital
    heart disease", one line above a sentence saying that is not CHD; measured
    2026-10-26). GATA4 here is graded for "structural congenital heart
    disease", so it is that page, independent of the real corpus: its card is
    asserted to name the disease, and its short form not to deny it.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )
    build_gene_pages(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    for name in ("index.html", "HGNC_11604.html", "HGNC_4173.html"):
        page = _page(tmp_path, name)
        # Both spellings: the build writes one apostrophe `'` and another
        # `&#x27;`, and a retired phrase must not pass by changing which.
        texts = (page, html.unescape(page))
        # **The published provenance claim, corrected 2026-08-06.** Every page
        # used to say a gene is definitive "for a disease in this atlas's CHD
        # scope", which told 24 published pages that this atlas decides what
        # counts as congenital heart disease. It does not, and its owner -- not
        # a clinician -- said so: scope is now a selection among terms external
        # authorities already use, and `curation/chd_scope.yaml` names the
        # authority for every term.
        #
        # The negative assertion is the load-bearing half. Without it a future
        # edit could restore the old sentence beside the new one and satisfy
        # every positive check here.
        assert all("this atlas's CHD scope" not in text for text in texts)
        # **Both warrants, or the sentence describes a gate the site does not
        # have**, and **sliced to the scope-rule element, because page-wide the
        # phrase checks could not fail** (CLAUDE.md section 4.19): both
        # recorded decisions moved with the full rule to `test_build_methods.py`
        # (2026-10-26, D50). The short form names neither warrant -- it says
        # whose classification it is and for what -- so what is pinned here is
        # the literal, and that the rule in full has not come back beside it.
        assert page.count('<p class="scope-rule">') == 1, f"{name}: one short form, once"
        rule = re.search(r'<p class="scope-rule">.*?</p>', page, re.S)
        assert rule, f"{name} carries no scope-rule element"
        sentence = rule.group(0)
        assert sentence == (
            '<p class="scope-rule">Published on an external authority&#x27;s classification '
            "for the disease it names, which for many genes here is a syndrome rather than "
            'CHD. <a href="../methods.html#scope">How a gene is admitted, and what each '
            "grade means</a>.</p>"
        ), f"{name}: the scope short form changed"
        assert 'href="../methods.html#scope"' in sentence
        assert all(
            "No disease is in scope on this atlas's own judgement" not in text for text in texts
        )

    # A page whose panel graded CHD itself: the card says so, and the short form
    # must not deny it.
    gata4 = _page(tmp_path, "HGNC_4173.html")
    assert "<dt>definitive for</dt><dd>structural congenital heart disease</dd>" in gata4
    short = _slice_between(gata4, '<p class="scope-rule">', ("</p>",))
    assert "not the same as for CHD" not in short


def test_a_gene_page_names_the_disease_beside_the_chip_that_says_definitive(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
) -> None:
    """The rail is the first screen, and the green chip on it said only `definitive`.

    `_rail` renders the symbol, the HGNC id and a `chip-definitive` pill, and on
    a viewport under 46rem `.layout` collapses to one column so those three fill
    the screen on their own. On a site titled "CHD Atlas" a green `definitive`
    with no disease reads as definitive for congenital heart disease; for CHD7 it
    means CHARGE syndrome, graded by the Hearing Loss GCEP.

    The order is asserted, not only the presence: the disease must come *after*
    the chip and before the rest of the rail's rows, so it is read as qualifying
    the chip rather than as one more fact three rows down.

    **Re-pointed at the card (D55), 2026-10-26.** "The rest of the rail's rows"
    was the `validity` row, which the card removed; it is now the card's first
    row, ahead of the gloss and the curation sentence that follow it.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity("CHARGE syndrome")},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    rail = _rail_of(page)
    assert _slice_between(rail, "<dl>", ("</dl>",)).startswith(
        "<dt>definitive for</dt><dd>CHARGE syndrome</dd>"
    )
    assert rail.index('chip chip-definitive">definitive') < rail.index("CHARGE syndrome")
    assert rail.index("CHARGE syndrome") < rail.index('<p class="card-gloss">')
    assert rail.index("CHARGE syndrome") < rail.index('<p class="card-curation">')


def test_every_browse_row_links_to_a_page_that_was_written(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_two,
        emitter,
        symbols={},
        validity={TBX5: _validity(), GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )
    build_gene_index_page(facts_two, emitter, symbols={}, validity=validity_two, burden_counts={})

    page = _page(tmp_path, "index.html")
    for name in ("HGNC_11604.html", "HGNC_4173.html"):
        assert name in page
        assert f"genes/{name}" in emitter.checksums


def test_a_facet_option_value_is_escaped_like_every_other_published_string(
    tmp_path: Path, facts_two: dict[str, GeneFacts]
) -> None:
    """The one interpolation on the site where escaping rested on provenance.

    Review on #15 measured this: every other value on every page reaches markup
    through a `render.py` primitive, and `<option value="{value}">` did not. It
    was safe only because all 21 members of `LesionGroup`, `Classification`,
    `ValidityState` and `AtlasCuration` happen to be `[a-z_]`. Driving the loop
    with `"><script>alert(1)</script>` emitted a live `<script>` element into
    the page while the same build's table cells stayed correctly escaped.

    The hostile value arrives as a `StrEnum` member rather than a bare `str`,
    because the facet loop reads `.value` off each field -- so a plain string
    never reaches the interpolation and a test using one would fail for the
    wrong reason, which is what happened when this test was first written.
    """

    class Hostile(StrEnum):
        X = '"><script>alert(1)</script>'

    facts = dict(facts_two)
    facts[GATA4] = replace(facts[GATA4], validity_state=Hostile.X)  # type: ignore[arg-type]
    emitter = Emitter(root=tmp_path)

    build_gene_index_page(facts, emitter, symbols={}, validity={}, burden_counts={})

    page = _page(tmp_path, "index.html")
    assert "<script>alert(1)</script>" not in page
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in page


def test_every_facet_names_a_data_attribute_the_filter_script_reads(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """A facet whose name matches no `data-*` attribute filters nothing, silently.

    `FILTER_SCRIPT` resolves each facet by string concatenation --
    `row.getAttribute('data-' + select.name)` -- so a `<select name="lesions">`
    over rows carrying `data-lesion` reads `null`, matches nothing, and hides
    every row the moment a reader touches that facet. Nothing in the type system
    couples the two spellings and no other test in this file reads either, so
    this is the coupling. Not in the plan; added because the mismatch is exactly
    the silent failure this project is written against.

    The empty `<option value="">` is asserted for the same reason: the script
    reads an empty select value as "no filter", so a facet lacking one has no
    way back to showing every gene. `name="q"` likewise -- the script does
    `form.querySelector('[name=q]').value`, which throws on a missing input and
    leaves every filter dead.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert 'name="q"' in page
    facets = re.findall(r'<select name="([^"]+)"[^>]*>(.*?)</select>', page)
    assert [name for name, _ in facets] == [
        "lesion",
        "confidence",
        "validity",
        "curation",
        "burden",
    ]
    for _, options in facets:
        assert options.startswith('<option value="">')
    for row in re.findall(r"<tr((?: data-[^>]*)?)>", page):
        if row:
            assert set(re.findall(r"data-([a-z]+)=", row)) == {"search", *(n for n, _ in facets)}


def test_every_browse_control_is_named_for_a_screen_reader(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Five form controls, five `aria-label`s. Raised by review on #14.

    Nothing else names any of them. The search box carries only a placeholder,
    which is announced inconsistently and disappears the moment a character is
    typed; each `<select>` carries only its own first option, and "any lesion"
    describes the control's *current value* rather than what it filters. A
    reader on a screen reader would meet five controls that announce nothing.

    Asserted against literals rather than against whatever the code emits,
    because a test that reads the label out of the page and compares it to
    itself would pass on an empty string. The labels differ deliberately from
    the neutral option text -- "Filter by lesion group" against "any lesion" --
    so an assertion that merely required *some* `aria-label` would not catch a
    label wired to the wrong facet.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert 'aria-label="Search by gene symbol or HGNC id"' in page
    assert re.findall(r'<select name="([^"]+)" aria-label="([^"]+)"', page) == [
        ("lesion", "Filter by lesion group"),
        ("confidence", "Filter by mirrored confidence"),
        ("validity", "Filter by validity state"),
        ("curation", "Filter by atlas curation"),
        ("burden", "Filter by whether rare-variant evidence exists"),
    ]
    # The facet's neutral option reads "any rare-variant result" (D49), while the
    # `burden` name and its yes/no values stay what the filter script reads.
    assert '<option value="">any rare-variant result</option>' in page


def _burden_row(**overrides: object) -> BurdenRow:
    """TAB2's syndromic loss-of-function row: the unbounded case, from real data."""
    payload: dict[str, object] = {
        "study": "PMID:42230622",
        "gene": GATA4,
        "cohort_stratum": "syndromic",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "maf_max": 0.001,
        "count_unit": "individuals",
        "n_case_carriers": 5,
        "n_cases": 1471,
        "comparator": "control_cohort",
        "n_control_carriers": 0,
        "n_controls": 45082,
        "expected_count": None,
        "effect": None,
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above",
        "ci_low": 28.1,
        "ci_high": None,
        "pvalue": 3.13e-08,
        "pvalue_test": "fisher_exact",
        "pvalue_adjusted": None,
        "pvalue_adjustment": None,
        "case_cohorts": ("cnchd", "ddd"),
        "control_cohorts": ("ukbb",),
        "method_note": None,
        "source": "audain2026_sd3",
    }
    payload.update(overrides)
    return BurdenRow(**payload)  # type: ignore[arg-type]


_COHORTS = {
    "cnchd": Cohort(id="cnchd", name="German Competence Network", description="German registry."),
    "ddd": Cohort(
        id="ddd",
        name="Deciphering Developmental Disorders (DDD)",
        description="Ascertained on developmental disorder, so enriched for syndromic CHD.",
    ),
    "ukbb": Cohort(
        id="ukbb",
        name="UK Biobank",
        description="Adults recruited at 40-69; survivorship bias against severe CHD.",
    ),
    "pcgc": Cohort(id="pcgc", name="PCGC", description="Pediatric Cardiac Genomics Consortium."),
}

_PUBLICATION = Publication(
    id="PMID:42230622",
    title="Assessing the contribution of rare variants to congenital heart disease",
    journal="NPJ genomic medicine",
    year=2026,
    authors=["Audain E", "Hitz MP"],
    study_type="case_control",  # type: ignore[arg-type]
    own_lab=True,
    tests_reported=138609,
)


_OTHER_PUBLICATIONS = {
    "PMID:34324492": Publication(
        id="PMID:34324492",
        title="Integrative analysis of genomic variants",
        journal="PLoS genetics",
        year=2021,
        authors=["Audain E", "Wilsdon A"],
        study_type="meta_analysis",  # type: ignore[arg-type]
        own_lab=True,
        tests_reported=11515,
    ),
}


_SIERANT_PUBLICATION = {
    "PMID:40127276": Publication(
        id="PMID:40127276",
        title="Landscape of rare variants in congenital heart disease",
        journal="Nature genetics",
        year=2025,
        authors=["Sierant MC"],
        study_type="case_control",  # type: ignore[arg-type]
        own_lab=False,
        tests_reported=4128,
    ),
}


def _burden_page(
    tmp_path: Path,
    facts: dict[str, GeneFacts],
    rows: list[BurdenRow],
    name: str = "HGNC_4173.html",
    publications: dict[str, Publication] | None = None,
    families: tuple[frozenset[str], ...] | None = None,
    axes: tuple[tuple[str, str], ...] | None = None,
    cohorts: dict[str, Cohort] | None = None,
) -> str:
    """Render one gene page.

    `families` and `axes` default to the derivation the runner performs, so a
    test that does not care about the matrix still renders a real one rather
    than an empty section -- and a test that does care can pin them. `cohorts`
    defaults to `_COHORTS`, and the families are derived over the same registry
    the page is handed, as the runner derives them (I3).
    """
    registry = _COHORTS if cohorts is None else cohorts
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity={GATA4: _validity(), TBX5: _validity()},
        assertions={},
        publications=publications or {_PUBLICATION.id: _PUBLICATION},
        burden={GATA4: rows},
        cohorts=registry,
        families=families if families is not None else cohort_families(rows, registry),
        axes=axes if axes is not None else evidence_axes(rows),
        cell_types={},
    )
    return _page(tmp_path, name)


def _burden_section_text(page: str) -> str:
    """The `Rare-variant association` section alone, up to the next `<h2>`.

    Sliced for the reason `_validity_table` and `_expression_section_text` are
    sliced (CLAUDE.md section 4.19). The developmental-expression section that
    follows this one draws a chart of its own on nearly every page, so a
    page-wide `"<svg" not in page` can be satisfied -- or defeated -- by a
    picture the burden section did not draw.
    """
    return _slice_between(page, _RARE_VARIANT_H2, ("<h2",))


def _sierant(**overrides: object) -> BurdenRow:
    """A PMID:40127276 de novo row: an enrichment ratio and a Benjamini-Hochberg q."""
    fields: dict[str, object] = {
        "study": "PMID:40127276",
        "comparator": "mutation_model",
        "effect": 41.9,
        "effect_measure": "enrichment_ratio",
        "effect_bound": None,
        "ci_low": None,
        "pvalue_adjustment": "benjamini_hochberg",
    }
    return _burden_row(**(fields | overrides))


@pytest.mark.parametrize(
    ("rows", "present", "absent"),
    [
        pytest.param(
            # SMAD2's shape: missense survives, LoF is nominal (p 0.022) but not.
            [
                _sierant(
                    consequence_class="missense_damaging", pvalue=5.8e-05, pvalue_adjusted=0.00257
                ),
                _sierant(consequence_class="lof", effect=8.0, pvalue=0.0221, pvalue_adjusted=0.145),
            ],
            ("Damaging Missense Significant", "q 0.00257", "Loss-of-function variants do not."),
            ("no enrichment", "ain-of-", "aploinsufficiency", "mechanism"),
            id="missense-survives-lof-nominal",
        ),
        pytest.param(
            [
                _sierant(
                    consequence_class="missense_damaging", pvalue=5.8e-05, pvalue_adjusted=0.00257
                )
            ],
            ("Damaging Missense Significant",),
            ("oss-of-function",),
            id="missense-survives-lof-untested",
        ),
        pytest.param(
            # PMID:34324492: deletions only, no effect size, a family-wise p.
            [
                _burden_row(
                    study="PMID:34324492",
                    effect=None,
                    effect_bound=None,
                    ci_low=None,
                    pvalue=9.999e-05,
                    pvalue_adjusted=0.0001,
                    pvalue_adjustment="familywise_permutation",
                )
            ],
            ("Loss-of-Function Significant", "corrected p 0.0001"),
            ("q 0.0001", "issense", "aploinsufficiency", _EM_DASH),
            id="cnv-familywise",
        ),
        pytest.param(
            # PMID:42230622 publishes no corrected p; this p clears even 0.05/138,609.
            [_burden_row()],
            ("publishes no corrected p",),
            ("does not survive", "does not reach"),
            id="no-correction-published",
        ),
        pytest.param(
            [_sierant(consequence_class="lof", pvalue=0.02, pvalue_adjusted=0.2)],
            ("does not survive this study's correction",),
            ("publishes no corrected p", "study-wide"),
            id="correction-published-not-survived",
        ),
    ],
)
def test_a_study_card_states_what_its_rows_show_and_no_mechanism(
    rows: list[BurdenRow], present: tuple[str, ...], absent: tuple[str, ...]
) -> None:
    """The finding card says which classes survive the study's own correction.

    It names no mechanism: no study here states one, and on 2026-11-02 31 of the
    93 published cards did (haploinsufficiency, gain-of-function). "Loss-of-function variants
    show no enrichment" was false on GDF1, SMAD2 and MYH6, whose LoF rows were
    nominal. A class the study never tested is not mentioned. The statistic is
    `_statistic`'s, so a family-wise p is never called q.
    """
    text = html.unescape(re.sub(r"<[^>]+>", "", _study_finding_summary(rows)))

    for phrase in present:
        assert phrase in text, (phrase, text)
    for phrase in absent:
        assert phrase not in text, (phrase, text)


def test_a_depleted_row_earns_no_enrichment_card() -> None:
    """A non-synonymous row with p < 0.05 below 1 is depletion, not "enrichment"."""
    depleted = _burden_row(effect=0.3, effect_bound=None, ci_low=0.1, ci_high=0.8, pvalue=0.01)

    assert _study_finding_summary([depleted]) == ""


def test_an_effect_size_is_never_rendered_without_the_measure_that_names_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The guard the single `effect` column was chosen against.

    One column carries odds ratios and de novo enrichments alike, which is what
    lets this schema absorb a fifth study without a migration -- and it is the
    one place two incomparable quantities could silently merge. An odds ratio of
    2.45 and a de novo enrichment of 2.45 are different claims, and a cell
    reading `2.45` under a header reading `effect` equates them.

    Both measures are rendered here from otherwise identical rows, so a
    `_effect` that dropped the label would produce two cells a reader could not
    tell apart -- and this test would see one string where it expects two.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1),
            _burden_row(
                cohort_stratum="all",
                comparator="mutation_model",
                n_control_carriers=None,
                n_controls=None,
                control_cohorts=(),
                expected_count=0.42,
                effect=2.45,
                effect_measure="enrichment_ratio",
                effect_bound=None,
                ci_low=1.2,
                ci_high=8.1,
                pvalue_test="poisson",
            ),
        ],
    )

    assert "OR 2.45 (95% CI 1.2–8.1)" in page
    assert "enrichment 2.45 (95% CI 1.2–8.1)" in page
    assert "<td>2.45</td>" not in page


def test_an_unbounded_odds_ratio_renders_its_lower_bound_rather_than_a_blank(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """34 rows of the committed mirror have no finite odds ratio, and they are
    the strongest results in the study.

    `Infinity` cannot be published -- `encode_json` sets `allow_nan=False`, and
    `JSON.parse` rejects it -- so the row carries a null effect and
    `effect_bound`. A page rendering that null as an em dash would throw away the
    finding: "at least 28.1" is what the study showed. The carrier counts beside
    it are what make it readable: 5 of 1,471 against 0 of 45,082.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "OR ∞ (95% CI 28.1–∞)" in page
    assert "<td>5 carriers / 1,471</td>" in page
    assert "<td>0 carriers / 45,082</td>" in page


def test_a_gene_with_no_burden_rows_gets_no_section_at_all(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Unlike the validity table, whose header names the columns and whose
    emptiness is itself an answer, an empty burden table would read as "this gene
    was studied and nothing was found" -- a claim no study made.
    """
    page = _burden_page(tmp_path, facts_uncurated, [])
    # The id, not the heading text: a heading renamed would make a text check
    # pass vacuously (it did, until 2026-10-26).
    assert 'id="rare-variants"' not in page
    assert _RARE_VARIANT_HEADING not in page


def test_two_studies_sharing_a_cohort_are_declared_not_independent(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Latent on the committed corpus, which carries one study.

    The atlas computes no pooled statistic because these cohorts overlap, and
    that decision is invisible to a reader comparing two tables by eye. The
    sentence is what makes it visible, and it names the collection rather than
    its id.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(), _burden_row(study="PMID:8988165", case_cohorts=("ddd", "pcgc"))],
    )

    assert "These two studies are not independent." in page
    assert "Deciphering Developmental Disorders (DDD)" in page


def test_two_studies_on_collections_declared_overlapping_are_named_and_counted_once(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Owner decision 2026-10-27 (I3): the declaration is the only link here.

    PMID:34324492 cites DECIPHER and PMID:42230622 cites DDD; no collection is
    shared. With the registry declaring the two overlapping, the page names
    the overlap the way it names a shared collection -- "may describe", for a
    likely overlap the registry records rather than one a shared id proves --
    and the matrix shows the two studies as one family row. The same rows over
    a registry without the declaration show neither, which is what every page
    carrying both studies showed until the declaration.
    """
    rows = [
        _burden_row(case_cohorts=("ddd",)),
        _audain2021(case_cohorts=("decipher",), control_cohorts=("gnomad",)),
    ]
    declared = {
        **_COHORTS,
        "ddd": _COHORTS["ddd"].model_copy(update={"overlaps_with": ["decipher"]}),
        "decipher": Cohort(
            id="decipher", name="DECIPHER", description="Submitted cases.", overlaps_with=["ddd"]
        ),
    }
    page = _burden_page(
        tmp_path, facts_uncurated, rows, publications=_ALL_STUDIES, cohorts=declared
    )

    notices = re.findall(r'<p class="notice-inline"><strong>These two studies.*?</p>', page)
    assert notices == [
        '<p class="notice-inline"><strong>These two studies are not independent.</strong> '
        "Audain E et al. 2021 draws on DECIPHER and Audain E et al. 2026 on Deciphering "
        "Developmental Disorders (DDD), collections this atlas records as overlapping, so "
        "their results may describe partly the same people and must not be combined.</p>"
    ]
    assert re.findall(r'<th scope="row">([^<]*)</th>', page) == [
        "Audain E et al. 2021 + Audain E et al. 2026"
    ]

    undeclared = _burden_page(
        tmp_path / "undeclared",
        facts_uncurated,
        rows,
        publications=_ALL_STUDIES,
        cohorts={
            **declared,
            "ddd": _COHORTS["ddd"],
            "decipher": declared["decipher"].model_copy(update={"overlaps_with": []}),
        },
    )
    assert "These two studies" not in undeclared
    assert re.findall(r'<th scope="row">([^<]*)</th>', undeclared) == [
        "Audain E et al. 2021",
        "Audain E et al. 2026",
    ]

    # **A study of the family with no row for this gene is marked** (review,
    # 2026-10-27): the row is still the family's, its cells come from the other
    # study alone, and an unmarked label named PMID:34324492 as reporting on a
    # gene it has no block for -- 46 of 118 matrix pages on a build of 5f7161b.
    # The families are the corpus's, as the runner derives them, not this
    # gene's: the gene carries PMID:42230622's row only. A one-study family
    # with no row is not marked: every cell of its row says `not tested`.
    lone = _burden_page(
        tmp_path / "lone",
        facts_uncurated,
        rows[:1],
        publications=_ALL_STUDIES,
        cohorts=declared,
        families=(*cohort_families(rows, declared), frozenset({"PMID:40127276"})),
    )
    assert "These two studies" not in lone
    assert re.findall(r'<th scope="row">([^<]*)</th>', lone) == [
        "Audain E et al. 2021 (no row for this gene) + Audain E et al. 2026",
        "Sierant MC et al. 2025",
    ]


# `test_the_section_says_what_a_missing_row_means_without_asserting_why` and
# `test_the_synonymous_negative_control_is_explained_only_where_one_exists` sat
# here until 2026-10-26, pinning the rare-variant reading notes a gene page
# folded under "How to read these numbers". The notes moved to
# `methods.html#rare-variants` (D50) and both tests' positives and recorded
# decisions moved with them, to
# `test_rare_variants_states_what_every_gene_page_folded_or_shortened` in
# `test_build_methods.py`. The refuted mechanism ("no variant of that class was
# seen in either group") is asserted absent on both pages: there by
# `test_no_sentence_the_move_retired_comes_back_here`, here by
# `test_the_general_reading_notes_are_on_methods_html_and_not_on_the_gene_page`.


def test_the_page_names_how_many_tests_the_study_ran(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The most serious finding of the 2026-08-05 review.

    Measured on the built site that day: of the 187 burden rows the 23 published
    gene pages render, **32 have p < 0.05 and 3 survive Bonferroni** over the
    study's own 138,609 comparisons. So 29 rows read as significant and are not,
    under a green `definitive` chip, on pages read by clinical geneticists.

    The atlas publishes no corrected p -- the supplement carries none, and
    computing one would be authoring a statistic (D12/D33). Naming the
    denominator is the honest alternative: it is a count the study made, and it
    is what lets a reader apply their own threshold. `0.05 / 138609` is stated
    as the arithmetic it is, not as a threshold the study endorsed.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    # Scoped to the method line. This asserted `"uncorrected" in page` until
    # 2026-08-05 and passed on the word from the *synonymous* paragraph, which is
    # a different sentence entirely -- the trap
    # `test_the_consequence_column_is_headed_for_the_column_it_renders` records.
    # When that paragraph became conditional the assertion failed, which is the
    # only reason the mis-scoping was found.
    method = page[page.index('<p class="method">') :]
    method = method[: method.index("</p>")]

    assert "138,609" in method
    assert "3.6e-07" in method
    assert "judge the raw one against the whole scan" in method
    # "138,609 such comparisons" until 2026-10-26: its antecedent, the design
    # clause, renders after it and folded under "How this study counted", so
    # "such" pointed at a sentence the reader had not met. Now self-contained.
    assert "This study reported 138,609 comparisons in all." in method
    assert "such comparisons" not in page


def test_the_cohort_caveats_a_curator_wrote_reach_the_reader(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """They reached no published byte until 2026-08-05.

    `curation/cohorts.yaml` records that UK Biobank's participants are adults
    while the CHD cases were largely enrolled in childhood -- a survivorship
    bias that inflates every odds ratio on every page -- and `models/cohort.py`
    says in as many words that `description` is where such caveats belong.
    `cohort_registry` returned `{id: name}` and dropped them.

    This is the project's characteristic failure in its purest form: the work
    was done, it was correct, it was validated, and no reader could reach it.

    **A record's open question reaches the same place** (owner's standing
    rule, 2026-10-27): after its description, escaped, on the one cohort of
    the three this row cites that carries one, and on no other. It is the
    `not_accessible` status here, so its label is pinned by literal too;
    `test_build_methods.py` pins `pending_review`'s.
    """
    registry = {
        **_COHORTS,
        "ukbb": _COHORTS["ukbb"].model_copy(
            update={
                "review_status": ReviewStatus.NOT_ACCESSIBLE,
                "review_note": "The <S2> sheet & its notes could not be read.",
            }
        ),
    }
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()], cohorts=registry)

    # In a `<details>`, so the caveats are one click from every table rather than
    # repeated at full length above 23 of them -- and still in the document for
    # a crawler, a `curl`, and a reader with no JavaScript. Sliced to the
    # cohort-notes fold itself (2026-10-26): `"<details" in page` passed on any
    # other fold on the page, the rare-variant reading notes among them, with
    # the caveats printed anywhere at all.
    notes = _slice_between(page, '<details class="cohort-notes">', ("</details>",))
    assert "survivorship bias against severe CHD" in notes
    # The title promises what the fold holds -- the collections, and the
    # caveats recorded for them -- and no more. It read "and what qualifies
    # every number above" until 2026-10-26, which several descriptions do not
    # do: some say only what a collection is.
    assert notes.startswith(
        "<summary>About these cohorts, and the caveats recorded for them</summary><ul>"
    )
    assert "qualifies every number" not in page
    assert "enriched for syndromic CHD" in notes
    assert (
        "<li><strong>UK Biobank</strong> &mdash; Adults recruited at 40-69; survivorship bias "
        "against severe CHD. <strong>Source not accessible:</strong> The &lt;S2&gt; sheet "
        "&amp; its notes could not be read.</li>"
    ) in notes
    assert notes.count("<strong>Source not accessible:</strong>") == 1
    # The fold stays (P4: narrowing where a cohort caveat reaches is what this
    # test exists to stop) and ends on a link to the whole registry on
    # methods.html, so this section links its anchor like every other
    # (2026-10-26, D50). Last in the fold, after the list it extends.
    assert notes.endswith(
        '</ul><p><a href="../methods.html#cohorts">Every cohort in this atlas</a></p>'
    )


def test_the_page_says_what_was_counted_not_only_how_many(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Four partition columns reached no page while populated on all 187 rows.

    A section headed "Rare variant burden" (now "Rare-variant association", D49)
    never said what *rare* meant, and
    `origin: any` -- a case-control count of variants regardless of inheritance
    -- was indistinguishable from a trio's de novo count. `vocab.VariantOrigin`
    warns about exactly that, and `tables.py` says the partition "exists to stop
    two incomparable rows from *looking* comparable".

    The origin label is asserted in full because the parenthetical is the whole
    point: "any inheritance" alone would still let a reader assume de novo.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "MAF below 0.001" in page
    assert "any inheritance (not a de novo test)" in page
    assert "SNVs and indels" in page
    assert "Fisher&#x27;s exact test" in page


def test_a_study_by_an_author_of_this_atlas_says_so(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`own_lab` reached `publications.json` and zero HTML files.

    Two of the three studies published today are by this repository's author, and
    a reader being asked to weigh their numbers should not have to
    cross-reference a JSON payload to learn that.

    This said "every burden table on the site today comes from one study whose
    first author is this repository's author" until 2026-08-05. PMID:40127276
    (Sierant MC et al., `own_lab: false`) made that false for 15 of the 23
    published gene pages, 8 of which now render three studies with the
    declaration on two of them. The sentence was the written record of when the
    line applies, so a stale version of it is worse than none -- it describes a
    corpus that no longer exists and invites the conclusion that the disclosure
    is unconditional.

    **Shortened to one sentence on 2026-10-26 (D50)**, and pinned whole, sliced
    to its own element; `methods.html#rare-variants` says the rest. The spec's
    approved short form ended "nothing added", and that is false where it sits:
    measured on a real build of 4233eca, 113 of the 190 own-lab study blocks
    also carry the atlas's own Bonferroni arithmetic (`0.05 / 138,609`) and
    draw a forest plot from the study's rows. It keeps the hedge the sentence
    it replaced had: no *assessment* added.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    disclosure = re.search(r'<p class="disclosure">.*?</p>', page, re.S)
    assert disclosure, "an own-lab study's block carries no disclosure"
    assert disclosure.group(0) == (
        '<p class="disclosure">Declaration: this study is by an author of this atlas; '
        "its results are republished as published and read by the same rule as every other "
        "study's.</p>"
    )
    for text in (page, html.unescape(page)):
        assert "nothing added" not in text
        assert "no assessment added" not in text

    outside = _PUBLICATION.model_copy(update={"own_lab": False})
    other = _burden_page(
        tmp_path, facts_uncurated, [_burden_row()], publications={outside.id: outside}
    )
    assert "by an author of this atlas" not in other


def test_a_de_novo_row_renders_the_expectation_it_was_compared_against(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_count(None, None)` rendered an em dash where the comparator belongs.

    On a `mutation_model` row the modelled expected count is the only thing the
    enrichment was computed against, and an em dash is indistinguishable from a
    control count nobody recorded. Latent until a trio study lands, which is
    precisely why it is pinned now.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                comparator="mutation_model",
                n_control_carriers=None,
                n_controls=None,
                control_cohorts=(),
                expected_count=0.42,
                effect=11.9,
                effect_measure="enrichment_ratio",
                effect_bound=None,
                ci_high=30.0,
                pvalue_test="poisson",
            )
        ],
    )

    assert "0.42 expected" in page
    assert "enrichment 11.9" in page


def test_a_method_note_is_rendered_as_the_footnote_the_schema_promises(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`tables.py` called `method_note` "rendered verbatim as a row footnote"
    from the day the column existed, and nothing read it.

    The column exists for the CNV case, where "carrier" means different things
    in different papers, so the first curator to need it would have got a green
    build, a correct bundle and no page text.

    A note of whitespace alone is no note: it rendered an empty `<li>` until
    2026-10-26, while the card's `caveats` label already read it as none.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(method_note="Carrier means any exonic overlap of the deletion."),
            _burden_row(cohort_stratum="all", method_note="   "),
        ],
    )

    assert re.findall(r'<ul class="footnotes">.*?</ul>', page) == [
        '<ul class="footnotes"><li>Carrier means any exonic overlap of the deletion.</li></ul>'
    ]


def test_the_provenance_line_unions_every_row_rather_than_reading_the_first(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`case_cohorts` is a per-row column and was reported per study block.

    A study whose strata drew on different collections published one row's
    cohorts as if they were all of them -- and `shared_cohorts` reads the same
    column with a *union* rule, so the "not independent" notice could name a
    cohort the provenance line did not. Both union now, so they cannot disagree.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(case_cohorts=("cnchd",)),
            _burden_row(cohort_stratum="all", case_cohorts=("ddd", "pcgc")),
        ],
    )

    # Scoped to the provenance paragraph. Asserted page-wide, this passed with
    # `_names` mutated to print bare ids, because `_cohort_notes` renders the
    # same names in the `<details>` block below the table.
    start = page.index('<p class="provenance">')
    line = page[start : page.index("</p>", start)]
    assert "German Competence Network" in line
    assert "Deciphering Developmental Disorders (DDD)" in line
    assert "PCGC" in line
    assert "cnchd" not in line


def test_a_study_heading_never_renders_a_nameless_author(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`Publication.authors` constrains the list, not its elements.

    A whitespace-only first author passes every validator and rendered
    "    et al. 2026" -- a study heading naming nobody. `authors[0]` was read
    only by `search.py` before this section existed, so the page is where the
    pre-existing model weakness first became visible. Same shape as the blank
    `symbol` defect `tables.cell_text` exists to catch.
    """
    nameless = _PUBLICATION.model_copy(update={"authors": ["   "]})
    page = _burden_page(
        tmp_path, facts_uncurated, [_burden_row()], publications={nameless.id: nameless}
    )

    assert "et al." not in page
    assert "PMID:42230622</a>" in page


def test_the_pooling_notice_is_the_matrix_caption_and_no_longer_conditional(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**Relocated, deliberately, and the condition went with it.**

    It used to render only where a gene carried two studies, because asserting
    "these cohorts overlap" above a single table sent a reader hunting for a
    second study that was not there.

    The evidence matrix changed what it qualifies. A summary that *counts
    datasets* is adjacent to pooling and invites exactly the arithmetic this
    sentence forbids -- so it is now the caption directly beneath the matrix,
    and it belongs there whenever the matrix does, one study or three. Folding
    it while introducing the thing it guards against would have been the worst
    available pairing.

    **Shortened to the spec's approved short form on 2026-10-26 (D50, P7)**,
    made explicit about what would double-count: the spec's "those sharing a
    collection would count" left the pooling implicit. Pinned whole, sliced to
    the section; the full sentence, with how a gene page names an overlap, is on
    `methods.html#rare-variants`, where "these cohorts overlap" is asserted
    absent too.

    **"people", not "children", since 2026-10-26.** `shared_cohorts` unions a
    study's case *and* control collections, and a shared collection may be a
    control cohort of adults -- UK Biobank's own note says "~500,000 UK adults
    recruited at ages 40-69" -- so the double-counted are not always children.
    The matrix caption and `docs/data-api.md` already said "people".
    """
    one = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert "no pooled statistic across studies" in one
    section = _burden_section_text(one)
    assert (
        "<p>The atlas computes <strong>no pooled statistic across studies</strong>: pooling "
        "two that share a sample collection would count the same people twice, and pooling "
        "two that draw on collections this atlas records as overlapping could.</p>" in section
    )
    assert "children" not in section
    # Immediately after the matrix, not four paragraphs above it.
    assert one.index("</table>") < one.index("no pooled statistic across studies")
    # **And it must be true of a one-study page.** It used to assert "these
    # cohorts overlap" in the present tense, which on this page sent a reader
    # hunting for a second study that is not there -- the defect the original
    # conditional existed to fix, reintroduced by relocating it. It is a
    # statement of policy now.
    assert "these cohorts overlap" not in one

    two = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(), _burden_row(study="PMID:8988165", case_cohorts=("ddd", "pcgc"))],
    )
    assert "no pooled statistic across studies" in two

    # **The caption also says the cell and the card lead differently** (owner
    # decision 2026-10-26, R4): the cell names the most significant row behind
    # its state, the card a study's surviving rows, a damaging union first.
    # Pinned as a literal, sliced to the caption -- from the matrix's end to
    # the first study block -- and rendered only where both can lead: with a
    # matrix, and where the card lists a surviving row (review, 2026-10-27).
    # On the 86 pages whose card lists none it was a sentence about a lead
    # that does not exist. `one`'s row publishes no corrected p, so nothing
    # survives; the same row with one does, and its card lists it.
    lead = (
        "<p>A cell names the most significant of the rows behind its state; the card lists a "
        "study's surviving rows with a damaging union before its parts, so the two can lead "
        "with different rows.</p>"
    )
    assert "A cell names the most significant" not in one
    survivor = _burden_row(pvalue_adjusted=1e-4, pvalue_adjustment="benjamini_hochberg")
    surviving = _burden_page(tmp_path / "surviving", facts_uncurated, [survivor])
    assert "surviving their own study's correction" in surviving, "the card lists no row"
    led = _burden_section_text(surviving)
    assert lead in led[led.index("</table></div>") : led.index("<h3")]
    # Two studies with no evidence axes draw no matrix, and keep the pooling
    # sentence without it -- and no lead note, though a row survives.
    bare = _burden_section_text(
        _burden_page(
            tmp_path / "bare",
            facts_uncurated,
            [survivor, _burden_row(study="PMID:8988165", case_cohorts=("ddd", "pcgc"))],
            axes=(),
        )
    )
    assert '<table class="matrix">' not in bare
    assert "no pooled statistic across studies" in bare
    assert "A cell names the most significant" not in bare


def test_the_consequence_column_is_headed_for_the_column_it_renders(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The header and the cell must name the same field.

    Raised by review on #17. `_BURDEN_HEADERS[1]` read "variant class" while the
    cell rendered `consequence_class`, and the method line above the table
    rendered the real `variant_class` -- so one phrase meant two things on one
    page and was wrong in the more prominent of them. The rename that introduced
    it was itself a fix for an ambiguous header.

    Asserted by pairing the header with a cell whose value could only have come
    from one field: `loss-of-function` is a `ConsequenceClass` member and
    `SNVs and indels` is a `VariantClass` label, so a header claiming the latter
    over a cell holding the former fails here.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    # Scoped to the burden section: `page.index("<thead>")` finds the mirrored
    # validity table, which is rendered first. Same trap as the provenance test.
    section = page[page.index(_RARE_VARIANT_HEADING) :]

    # Past the evidence matrix: it is now the first table in this section and
    # has its own `<thead>`. The comment above already records that
    # `page.index("<thead>")` finds the validity table; the matrix is the same
    # trap one layer in.
    section = section[section.index("</table>") :]
    header = section[section.index("<thead>") : section.index("</thead>")]
    assert '<th scope="col">consequence</th>' in header
    assert "variant class" not in header

    body = section[section.index("<tbody>") : section.index("</tbody>")]
    assert "loss-of-function" in body
    # The real `variant_class` is on the page exactly once, in the method line,
    # and never in this table.
    assert "SNVs and indels" not in body
    # Twice on the page, and neither is in this table: once as an evidence-matrix
    # column head, and once in the folded method line. It read `== 1` until the
    # matrix landed, when the column head made it 2 -- the number is asserted
    # rather than loosened so a third occurrence has to be a deliberate change.
    assert page.count("SNVs and indels") == 2


def test_the_browse_page_shows_how_much_burden_evidence_each_gene_has(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """`genes/index.json` published `burden_row_count`; the page showed nothing.

    Measured on the built site 2026-08-05: `genes/index.html` contained the word
    "burden" zero times -- no column, no filter, no count -- while every row of
    the payload behind it carried the number. A reader scanning 23 rows that all
    read `definitive` could not tell that four of them (ISL1, NR2F2, RBM10,
    SMAD2) carry no loss-of-function burden evidence at all.

    The em dash for a gene with none is deliberate and matches the `lesion
    groups` column beside it: `0` reads as a measured zero, while an em dash
    reads as "nothing here", which is what an absent row means.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={TBX5: 9},
    )

    page = _page(tmp_path, "index.html")

    assert '<th scope="col">rare variants across studies</th>' in page
    rows = {
        match.group(1): match.group(0)
        for match in re.finditer(r'<tr data-search="(hgnc:\d+)[^"]*".*?</tr>', page)
    }
    assert "9 tests" in rows["hgnc:11604"]
    assert "untested" in rows["hgnc:4173"]

    # The facet and the attribute must spell the same word; see
    # `test_every_facet_names_a_data_attribute_the_filter_script_reads`.
    assert 'data-burden="yes"' in rows["hgnc:11604"]
    assert 'data-burden="no"' in rows["hgnc:4173"]


def test_the_rail_states_rare_variant_evidence_and_no_bare_publication_count(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Two rail defects, both raised by review 2026-08-05, and what replaced the fixes.

    The rail read `publications | 0` on 22 of the 23 published genes while the
    burden section one column away cited a linked, PubMed-referenced study. The
    count was right -- `fact.publications` is what this atlas's *own* assertion
    evidence cites -- and the label was wrong, so the label became `curated
    publications`. And the rail carried no burden count at all, though
    `burden_row_count` exists and is published in every browse row, so it
    gained `rare-variant rows`.

    **Renamed and re-pointed at the card (D55), 2026-10-26** -- it was
    `test_the_rail_counts_burden_rows_and_names_the_publications_it_counts`.
    The card drops both rows: a zero count of curated publications says nothing
    a clinician can act on, and a row count is replaced by the rare-variant
    block, which says what the rows found. So the burden evidence is asserted
    in words, and every row label this test ever refused stays refused.
    """
    page = _burden_page(
        tmp_path, facts_uncurated, [_burden_row(), _burden_row(cohort_stratum="all")]
    )

    rail = _slice_between(page, '<aside class="rail">', ("</aside>",))
    card = _card_variants(rail)
    assert (
        "<p>Tested in 1 independent dataset (as far as the cohort registry records) (a study by "
        "an atlas author).</p>"
    ) in card
    assert '<p class="card-uncorrected"><a href="#study-42230622">' in card
    assert "<dt>rare-variant rows</dt>" not in rail
    assert "<dt>curated publications</dt>" not in rail
    assert "<dt>burden rows</dt>" not in rail
    assert "<dt>publications</dt>" not in page


def test_a_published_correction_is_rendered_and_names_its_method(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The column that answers the sharpest finding of the 2026-08-05 review.

    That review measured 187 uncorrected p-values on gene pages with nothing to
    judge them against. The atlas still computes no correction -- that would be
    authoring a statistic (D12/D33) -- but PMID:34324492 publishes one, and it
    **changes the conclusion**: CHD7 has a raw permutation p of 0.0068 and a
    family-wise corrected p of 0.991. Without this column the page shows a
    number that reads as significant for a gene the study found nothing for.

    The method is named in the cell because "0.991" alone does not say what it
    was corrected against, and a family-wise permutation correction and a
    Bonferroni factor are different claims.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                variant_class="cnv_deletion",
                effect=None,
                effect_measure=None,
                effect_bound=None,
                ci_low=None,
                ci_high=None,
                pvalue=0.0068,
                pvalue_test="permutation",
                pvalue_adjusted=0.991,
                pvalue_adjustment="familywise_permutation",
            )
        ],
    )

    assert '<th scope="col">corrected p</th>' in page
    assert "<td>0.991 (family-wise)</td>" in page
    assert "<td>0.0068</td>" in page


def test_the_method_line_tells_a_corrected_study_apart_from_an_uncorrected_one(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """A sentence true of one study is false of another.

    PMID:42230622 publishes no correction, so the reader needs the denominator
    and the advice to judge against the whole scan. PMID:34324492 corrects every
    row, so that same advice would point them at the column they should *not* be
    reading. The block says whichever is true of its own rows.

    **The corrected branch names no column since 2026-10-26 (P9).** It said
    "The corrected p column carries the study's own correction", and where the
    study has a forest plot that column is folded under "every figure from this
    study" while the plot names the same number `q` or `corrected p`. It now
    names the statistic, with the `q` it is shown as where it is a
    false-discovery rate (`_statistic`), and is pinned identical on a block with
    a forest and on one without, so it cannot come to depend on either layout.

    **And it states no count (2026-10-26).** It opened "This study reported N
    comparisons in all", and for PMID:40127276 N was `tests_reported: 248` --
    the study's panel *genes*, each tested in two arms over three consequence
    classes -- so the sentence was false on all 36 pages carrying that study.
    (That `tests_reported` is null since owner decision W1, 2026-10-26; the 248
    is its `correction_genes`.) A
    study that corrected every row needs no denominator from the reader; the
    count is shown only where the reader must apply it, beside the Bonferroni
    arithmetic, which is true of PMID:42230622 (138,609 is one Fisher test per
    row of its Supplementary Data 3). Pinned on a corrected study whose
    publication records no count at all, too: the sentence no longer reads
    one, so it must not depend on it. And pinned absent where no row carries a
    p-value at all -- a case series -- because there is then no corrected p to
    read: with the branch keyed on "no uncorrected row" alone it told such a
    block to read one (a mutant doing exactly that survived every other test,
    measured 2026-10-26).
    """
    read_corrected = (
        '<p class="method">Its <strong>corrected p</strong> (shown as <em>q</em> where it is a '
        "false-discovery rate) carries the study's own correction for multiple testing; "
        "read it rather than the raw p.</p>"
    )
    corrected = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(pvalue_adjusted=0.991, pvalue_adjustment="bonferroni")],
    )
    assert "<svg" in _burden_section_text(corrected), "this fixture must draw a forest"
    assert read_corrected in corrected
    assert "judge the raw one against the whole scan" not in corrected

    unplotted_row = _burden_row(
        study="PMID:34324492",
        variant_class="cnv_deletion",
        effect=None,
        effect_measure=None,
        effect_bound=None,
        ci_low=None,
        ci_high=None,
        pvalue_adjusted=0.991,
        pvalue_adjustment="familywise_permutation",
    )
    unplotted = _burden_page(
        tmp_path / "unplotted",
        facts_uncurated,
        [unplotted_row],
        publications=_OTHER_PUBLICATIONS,
    )
    assert "<svg" not in _burden_section_text(unplotted), "this fixture must draw no forest"
    assert read_corrected in unplotted

    uncounted = _OTHER_PUBLICATIONS["PMID:34324492"].model_copy(update={"tests_reported": None})
    no_count = _burden_page(
        tmp_path / "no_count",
        facts_uncurated,
        [unplotted_row],
        publications={uncounted.id: uncounted},
    )
    assert read_corrected in no_count
    for page in (corrected, unplotted, no_count):
        assert "comparisons in all" not in page
        assert "corrected p</strong> column" not in page

    series = _burden_page(
        tmp_path / "series",
        facts_uncurated,
        [
            _burden_row(
                comparator="none",
                n_control_carriers=None,
                n_controls=None,
                control_cohorts=(),
                effect=None,
                effect_measure=None,
                effect_bound=None,
                ci_low=None,
                ci_high=None,
                pvalue=None,
                pvalue_test=None,
            )
        ],
    )
    assert _RARE_VARIANT_H2 in series, "the case series must still render its block"
    assert "corrected p</strong> (shown as" not in series
    assert "comparisons in all" not in series

    plain = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert (
        '<p class="method">This study reported 138,609 comparisons in all. Where it published '
        "no corrected p, judge the raw one against the whole scan rather than against 0.05 "
        "&mdash; a Bonferroni threshold over that many is p &lt; 3.6e-07.</p>" in plain
    )
    assert "read it rather than the raw p" not in plain


def test_every_count_cell_names_what_it_counted(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_count`'s rule, and `_effect`'s rule applied one column over.

    Until 2026-08-05 the header read "cases (carriers / n)" and every cell was a
    bare "5 / 1,471", which was true only because every study curated so far
    counted people. PMID:40127276's case-control arm counts *alleles* -- measured
    against its own Dataset S4 in `CountUnit`, where its D-Mis `Obs` tracks the
    variant count in 245 of 248 genes and the distinct-proband count in only
    235 -- so under the old header its rows would have claimed 21,768 people
    sequenced where 11,555 were.

    The three units render in one table here because that is the real case: one
    study contributes a de novo arm counting mutations against trios and a
    case-control arm counting alleles, into the same table, and a per-column
    header cannot label either.

    Scoped to the burden `<tbody>`. A page-wide assertion passes with the unit
    dropped from the cells, because `_BURDEN_PREAMBLE` names all three words
    higher up the page -- the trap `test_the_consequence_column_is_headed_for_
    the_column_it_renders` records one section above.
    """
    rows = [
        _burden_row(consequence_class="lof", count_unit="individuals"),
        _burden_row(consequence_class="missense_damaging", count_unit="alleles"),
        _burden_row(
            consequence_class="missense_all",
            count_unit="de_novo_mutations",
            origin="de_novo",
            comparator="mutation_model",
            n_control_carriers=None,
            n_controls=None,
            control_cohorts=(),
            expected_count=0.166,
            effect=None,
            effect_measure=None,
            effect_bound=None,
            ci_low=None,
        ),
    ]
    section = _burden_page(tmp_path, facts_uncurated, rows)
    section = section[section.index(_RARE_VARIANT_HEADING) :]
    # Past the evidence matrix, whose cells also carry counts.
    section = section[section.index("</table>") :]
    body = section[section.index("<tbody>") : section.index("</tbody>")]

    assert "<td>5 carriers / 1,471</td>" in body
    assert "<td>5 alleles / 1,471</td>" in body
    # The denominator of a de novo row is trios, and the word travels with it:
    # 5 mutations in 1,471 *families* is not 5 in 1,471 alleles.
    assert "<td>5 de novo / 1,471 trios</td>" in body

    # The control column takes the same unit, so an allele row does not describe
    # its cases in alleles and its controls in people.
    assert "<td>0 alleles / 45,082</td>" in body
    assert "<td>0 carriers / 45,082</td>" in body

    # No cell anywhere in the table is a bare "n / n": that spelling is what the
    # column meant before the unit existed, and it is the mutant this kills.
    assert "<td>5 / 1,471</td>" not in body
    assert "<td>0 / 45,082</td>" not in body


def test_a_count_unit_nothing_has_taught_the_renderer_still_names_itself(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_COUNT_UNIT_LABEL` is hand-maintained, so it can fall behind `CountUnit`.

    A member added to the enum and forgotten here must degrade to its raw token
    rather than to silence: "5 genomes / 1,471" is ugly and true, while "5 /
    1,471" is the exact false claim this column was added to prevent. Same shape
    as `_study_label` falling back to the bare PMID rather than an em dash.

    Reached only behind a bypassed gate -- `validate_table` refuses a value
    outside `CountUnit`, and `build_site` refuses the corpus -- which is why it
    is pinned here rather than left to be discovered.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(count_unit="genomes")])
    section = page[page.index(_RARE_VARIANT_HEADING) :]
    # Past the evidence matrix, which is now the first table in this section and
    # has a `<tbody>` of its own. Slicing on the first one silently measured the
    # matrix instead -- the same trap this file records twice already.
    section = section[section.index("</table>") :]
    body = section[section.index("<tbody>") : section.index("</tbody>")]

    assert "<td>5 genomes / 1,471</td>" in body
    assert "<td>5 / 1,471</td>" not in body


def test_the_composite_row_is_named_as_a_union_of_the_two_below_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Three consequence rows must not read as three independent findings.

    PMID:40127276 reports `damaging (LOF + missense)` alongside its two
    components, and the composite is the analysis its 60 genes are defined by.
    Measured on CHD7: the damaging de novo row is 20 mutations and the two rows
    below it are 16 and 4 of *those same* 20, so a reader adding them up gets 40
    from 20 variants.

    The note is conditional, and both halves are asserted. A study reporting only
    components must not be told its rows decompose something that is not there --
    the defect `_POOLING_NOTICE` was made conditional for, where an
    unconditional sentence sent readers hunting for a second study that did not
    exist.
    """
    components = [
        _burden_row(consequence_class="lof"),
        _burden_row(consequence_class="missense_damaging"),
    ]
    composite = [_burden_row(consequence_class="damaging"), *components]

    with_composite = _burden_page(tmp_path, facts_uncurated, composite)
    assert "union</strong> of the loss-of-function" in with_composite
    # One sentence since 2026-10-26 (D50, P8), pinned whole and sliced to the
    # section. Why the composite is shown at all -- its p-value is not derivable
    # from its components' -- moved to `methods.html#rare-variants`.
    assert (
        '<p class="notice-inline">The <strong>damaging (LOF + missense)</strong> rows are '
        "the <strong>union</strong> of the loss-of-function and damaging-missense rows "
        "below them, not a third independent result &mdash; the same variants are counted "
        "in both.</p>" in _burden_section_text(with_composite)
    )
    assert "not derivable from the other two" not in with_composite

    without = _burden_page(tmp_path / "b", facts_uncurated, components)
    assert "union</strong> of the loss-of-function" not in without


def test_every_particular_caveat_survives_outside_the_fold(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**The fold rule needs a test, not good intentions.**

    A caveat may fold if it is *general* -- true on every gene page, teaching
    how to read the table. It must stay visible if it is *particular* -- it
    fired because of this gene's data and changes what a reader concludes about
    this gene.

    This project has twice measured that a caveat nobody reads is nearly a
    caveat that does not exist, so folding a particular one later has to turn
    this red. Asserted against the page with every `<details>` block removed,
    which is what "unfolded" means to the reader who never clicks.

    The multiple-testing warning is in this list deliberately. It was the
    sharpest finding of the 2026-08-05 review -- 32 of 187 published rows clear
    0.05 and 3 survive Bonferroni over the study's own 138,609 comparisons --
    so `_method_line` returns it separately from the design clause, which does
    fold. Reuniting them would bury it one release after it was added.

    The pooling caveat's wording became the spec's approved short form on
    2026-10-26 ("count the same children twice", D50, P7), and "children"
    became "people" the same day (a shared collection may be a control cohort
    of adults); the rule -- visible outside every fold -- is unchanged. The
    scan-size assertion is pinned to the method line's own words since the
    same day: the card's no-correction line also carries "138,609" ("judge its
    results against its scan of 138,609 tests"), so the bare number would pass
    with the method line folded away.

    **Two chart caveats joined the list on the owner's decision of 2026-10-26
    (E1)**, each particular to its figure and each folded under "How to read
    this chart" until then: a forest panel's warning that a `union`-tagged
    row is not fully split (10 panels on 10 pages, measured on a build of
    3b766ba), and the trajectory's note that its ticks are stages below the
    detection floor (32 of 112 trajectories). Each is required outside every
    fold, once, where it fires -- the `all cases` row below has its syndromic
    half drawn and no non-syndromic one, and the heart series has a stage
    below the floor -- and absent from the whole page where it does not, so a
    clause moved back into the toggle fails the first half and one rendered
    unconditionally fails the second.
    """
    rows = [
        _burden_row(consequence_class="lof"),
        _burden_row(
            study="PMID:34324492",
            consequence_class="lof",
            case_cohorts=("taa_cases",),
            method_note="Cases include 777 sporadic thoracic aortic aneurysm probands.",
        ),
    ]
    publications = {_PUBLICATION.id: _PUBLICATION, **_OTHER_PUBLICATIONS}
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [*rows, _burden_row(cohort_stratum="all", consequence_class="lof")],
        publications=publications,
    )
    unfolded = re.sub(r"<details.*?</details>", "", page, flags=re.S)
    split = (
        "The <strong>all cases · loss-of-function</strong> row is a union, and this study "
        "published no non-syndromic · loss-of-function row for this gene."
    )
    assert split in unfolded
    assert page.count("row is a union") == 1, "the split warning is on the page twice"
    # The same studies with every tagged row fully split -- here, no union at all.
    whole = _burden_page(tmp_path / "split", facts_uncurated, rows, publications=publications)
    assert "row is a union" not in whole

    # The trajectory's floor note, against a heart series with and without a
    # stage below the floor.
    def trajectory(*stages: tuple[str, float | None], into: str) -> str:
        return _expression_page(
            tmp_path / into,
            {GATA4: _heart_series(*stages)},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )

    ticks = "Ticks on the axis mark stages below this dataset's detection floor (1 tpm)"
    below = trajectory(("4wpc", 5.0), ("5wpc", None), ("6wpc", 40.0), into="below")
    assert ticks in re.sub(r"<details.*?</details>", "", below, flags=re.S)
    assert below.count("Ticks on the axis") == 1, "the floor note is on the page twice"
    above = trajectory(("4wpc", 5.0), ("5wpc", 12.0), ("6wpc", 40.0), into="above")
    assert "Ticks on the axis" not in above

    # The atlas's claim about itself.
    assert "computes none of them" in unfolded
    # The relocated pooling caption -- it qualifies the matrix directly above it.
    assert "count the same people twice" in unfolded
    # A row's own method note: this is where the TAA contamination, which is not
    # congenital heart disease, reaches a reader at all.
    assert "thoracic aortic aneurysm" in unfolded
    # The multiple-testing denominator, in the method line's own words.
    assert "138,609" in unfolded
    assert "138,609 comparisons in all" in unfolded


def test_the_general_reading_notes_are_on_methods_html_and_not_on_the_gene_page(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The other half of the rule, so "fold everything" fails as loudly as
    "fold nothing".

    These paragraphs are identical wherever they render and teach how to read
    the table rather than saying anything about this gene. Nothing is lost by
    moving them: the unit word is still in every count cell, and the matrix
    renders `not tested` as its own state rather than as a gap.

    **They left the gene page on 2026-10-26 (D50)**, and this test was named
    `test_the_general_reading_notes_are_folded_out_of_the_way` until then: they
    were folded under "How to read these numbers" on all 118 pages with a
    rare-variant section -- the synonymous paragraph on the 113 of them with a
    synonymous row -- and are now on `methods.html#rare-variants`, where
    `test_build_methods.py` pins them. What stays on the page is the section's
    one-line attribution, scoped to the study tables and linking there -- the
    evidence matrix beside it shows tallies the atlas derives ("n = N" on 86
    pages and "N rows" on 48, measured on a build of 4233eca), so "every count
    and statistic below" was false of the matrix -- and the pooling caption
    beneath the matrix (`test_the_pooling_notice_is_the_matrix_caption_and_no_longer_conditional`).

    The moved paragraphs are asserted absent from the **whole** page, in both
    apostrophe spellings, by a phrase unique to each, not by "negative
    control": the forest key may name the synonymous row's role, and a
    negative on that phrase would fail on a caption it was never about -- and
    since the forest key folded behind its own toggle (2026-10-26, D57) it
    does: its muted clause says "its negative control". The fixture is a
    synonymous row, so that clause renders, and it must not point at a fold
    that has left the page; the toggle it sits in closes on a link to the
    notes where they now are.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(consequence_class="synonymous")])
    section = _burden_section_text(page)

    assert section.startswith(
        "<p>Every count and statistic in the study tables below is an upstream study's, "
        "republished as published; the atlas computes none of them. "
        '<a href="../methods.html#rare-variants">How to read these results</a>.</p>'
    )
    assert '<details class="reading-notes"><summary>How to read these numbers' not in section
    figure = _forest_figures(section)[0]
    assert (
        "A <strong>muted</strong> row is one of this study's <strong>synonymous</strong> "
        "rows, its negative control, drawn apart from the results it calibrates.</p>"
        '<p><a href="../methods.html#rare-variants">How to read rare-variant results</a></p>'
        "</details>" in figure
    )
    for text in (page, html.unescape(page)):
        for moved in (
            "every cell names what it counted",
            "synonymous variants should show no enrichment",
            "Only a study that publishes one has it",
            "How to read these numbers",
            # The refuted mechanism for an absent row (2026-08-05).
            "no variant of that class was seen in either group",
            "every count and statistic below",
            "Every count and statistic below",
        ):
            assert moved not in text, f"{moved!r} is back on the gene page"


def test_the_matrix_names_an_untested_design_rather_than_leaving_a_blank(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The holes are the point, and a blank cell is indistinguishable from a bug.

    Columns come from the whole corpus, so a design nobody ran for this gene
    renders as absent instead of vanishing. Nobody has published a CNV de novo
    analysis for any gene, and that quadrant is empty on every page: the shape
    of the literature is itself information.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(consequence_class="lof")],
        families=(frozenset({"PMID:42230622"}), frozenset({"PMID:34324492"})),
        axes=(("snv_indel", "control_cohort"), ("cnv_deletion", "control_cohort")),
    )
    section = page[page.index(_RARE_VARIANT_HEADING) :]
    matrix = section[section.index('<table class="matrix">') : section.index("</table>")]

    # **The cell reads short and explains itself on hover.** The full sentence
    # set the column width for the whole table and pushed every cell wide enough
    # to crowd its neighbours, so it moved to the `title` -- where it is still
    # reachable, which is the half of this that must not regress.
    assert ">not tested</span>" in matrix
    assert 'title="this dataset did not test this gene"' in matrix

    # Every *data* cell says what it is in words and is never simply blank. The
    # one bare `<td>` is the corner of a table headed on both axes, which
    # carries no value by construction.
    assert matrix.count("<td></td>") == 1
    assert "<td></td>" not in matrix[matrix.index("<tbody>") :]


def test_the_gene_and_its_symbol_both_link_to_the_page_and_the_strip_does_not(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Where the affordance belongs, and where it does not.

    Both name cells link: `HGNC:11604` is the identifier a reader cites and
    `TBX5` is the name they scan for, and linking only the id put the affordance
    on the string nobody reaches for.

    **The strip is deliberately not a link.** It was one briefly, on the
    reasoning that it summarises the gene page's matrix -- but the row already
    links to that page twice, and a third link wrapped around a tally puts the
    affordance on a number rather than on a name. Each dot keeps its `title`,
    which is the detail the cell owes a reader.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={},
        concordance={
            gene: {"tested": 0, "enriched": 0, "corrected": 0, "families": []} for gene in facts_two
        },
    )
    page = _page(tmp_path, "index.html")
    href = f'href="../{gene_page_path(TBX5)}"'

    # Twice: once from the id, once from the symbol.
    assert page.count(href) == 2
    assert f"<a {href}>{TBX5}</a>" in page
    assert f"<a {href}>TBX5</a>" in page

    # And the strip carries no anchor at all.
    strip = page[page.index('<span class="strip">') :]
    strip = strip[: strip.index("</span></span>")]
    assert "<a " not in strip
    assert "strip-tally" in strip


def test_the_gene_index_prints_each_named_studys_summary_caveat_once(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Under the strip legend, in study-id order, for every study a strip names.

    The strips summarise each study's result for each gene -- "enriched, and
    survives that study's own correction" -- so the caveat the card prints
    beside a study's counts is owed here too. PMID:34324492 is named on both
    rows and noted once; the rows meet PMID:40127276 first, so the order is
    the sort's; PMID:42230622 carries a caveat and is named by no strip, and
    PMID:99999999 is named and has no record. Label and caveat are escaped.
    """

    def family(*studies: str, state: str = "corrected") -> dict[str, object]:
        return {"studies": list(studies), "state": state}

    publications = {
        "PMID:34324492": _OTHER_PUBLICATIONS["PMID:34324492"].model_copy(
            update={"summary_caveat": "Cases include <777> TAA probands & CHD."}
        ),
        "PMID:40127276": _SIERANT_PUBLICATION["PMID:40127276"].model_copy(
            update={"authors": ["D'Alessandro M"], "summary_caveat": "Counts alleles."}
        ),
        "PMID:42230622": _PUBLICATION.model_copy(update={"summary_caveat": "Named nowhere."}),
    }
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={},
        concordance={
            TBX5: {
                "tested": 2,
                "enriched": 2,
                "corrected": 1,
                "families": [
                    family("PMID:40127276"),
                    family("PMID:34324492", "PMID:99999999", state="nominal"),
                ],
            },
            GATA4: {
                "tested": 1,
                "enriched": 0,
                "corrected": 0,
                "families": [family("PMID:34324492", state="no_enrichment")],
            },
        },
        publications=publications,
    )
    page = _page(tmp_path, "index.html")
    notes = re.findall(r'<p class="strip-legend"><strong>.*?</p>', page[page.index("</table>") :])

    assert notes == [
        '<p class="strip-legend"><strong>Audain E et al. 2021</strong> (PMID:34324492): '
        "Cases include &lt;777&gt; TAA probands &amp; CHD.</p>",
        '<p class="strip-legend"><strong>D&#x27;Alessandro M et al. 2025</strong> '
        "(PMID:40127276): Counts alleles.</p>",
    ]
    assert "Named nowhere." not in page


def test_the_browse_strips_four_glyphs_are_all_named_in_its_legend(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """A key naming three of four states leaves the fourth as an unexplained mark.

    Measured by mutation 2026-08-05: replacing every label in `_STRIP_LEGEND`
    with a single letter survived this whole file. The gene page's matrix legend
    was guarded and the browse page's strip legend was not, which is the same
    asymmetry — and the same direction — as the missing matrix caption an
    adversarial review caught one release earlier.

    **Literals, and each caption asserted adjacent to its own glyph.** The
    original version of this test read the captions out of
    `EVIDENCE_STATE_LABELS` and checked membership, which passes under any
    permutation of that constant because the needles come from the thing being
    permuted. Measured 2026-08-06: swapping two entries survived all 787 tests
    and published a browse page whose key captions the untested glyph "tested,
    no enrichment detected" -- the `not_tested`-reads-as-`no_enrichment`
    collapse this project forbids, on a green build.

    **And the power caveat beside them** (review, 2026-10-27). The front page
    taught this key with the caveat until the doorway page moved both off it,
    pointing here; the caveat on this page was unguarded, and deleting it from
    `_STRIP_LEGEND` survived every test.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={},
        concordance={
            gene: {"tested": 0, "enriched": 0, "corrected": 0, "families": []} for gene in facts_two
        },
    )
    page = _page(tmp_path, "index.html")
    legend = page[page.index("</table>") :]

    for css, caption in (
        ("dot full", "enriched, and a result survives its own study&#x27;s correction"),
        ("dot half", "enriched nominally, or no correction published"),
        ("dot none", "tested, no enrichment detected"),
        ("dot untested", "not tested by that dataset"),
    ):
        pair = re.search(rf'<span class="{re.escape(css)}">[^<]*</span>\s*([^<]+)<', legend)
        assert pair, f"the strip legend no longer draws {css!r} with a caption"
        assert pair.group(1).strip() == caption, (
            f"{css!r} is captioned {pair.group(1).strip()!r}, expected {caption!r}"
        )
    assert (
        "<strong>No enrichment at these cohort sizes is not evidence against a gene</strong> "
        "&mdash; rare-variant tests routinely detect nothing for genes with overwhelming family "
        "and functional evidence."
    ) in legend


def test_the_strip_hover_names_each_study_of_a_family_with_its_own_state(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """A glyph is a family's; its hover says what each study in it says (review, 2026-10-27).

    Since I3 one family holds PMID:34324492 and PMID:42230622, and the hover
    read "<both ids>: <family state>". Measured on a build of 5f7161b: 10
    browse rows said both survive their own correction where only
    PMID:34324492's rows do -- PMID:42230622 publishes no corrected p -- and 46
    gave a state to a study with no row for the gene. Here TBX5's family is
    `corrected` on PMID:34324492's row alone and PMID:42230622's is `nominal`;
    GATA4's family is PMID:42230622's `no_enrichment`, and PMID:34324492 has no
    row for it, so it keeps the `not_tested` sentence. Rendered from the
    concordance the runner derives, so the payload's `members` and the hover
    are tested as one.
    """
    families = (
        frozenset({"PMID:34324492", "PMID:42230622"}),
        frozenset({"PMID:40127276"}),
    )
    rows = {
        TBX5: [_audain2021(gene=TBX5), _burden_row(gene=TBX5)],
        GATA4: [_burden_row(pvalue=0.6)],
    }
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={},
        concordance={gene: gene_concordance(rows[gene], families) for gene in facts_two},
        publications=_ALL_STUDIES,
    )
    page = _page(tmp_path, "index.html")
    strips = re.findall(r'<span class="strip">(.*?)<span class="strip-tally">', page)

    one_family = "One cohort family, counted once."
    assert [re.findall(r'title="([^"]*)"', strip) for strip in strips] == [
        [
            "PMID:34324492: enriched, and survives this study&#x27;s own correction. "
            "PMID:42230622: enriched nominally; not after correction, or no correction "
            f"published. {one_family}",
            "PMID:40127276: this dataset did not test this gene",
        ],
        [
            "PMID:34324492: this dataset did not test this gene. PMID:42230622: tested, no "
            f"enrichment detected. {one_family}",
            "PMID:40127276: this dataset did not test this gene",
        ],
    ]
    # The glyph and the tally stay the family's: TBX5's family is full.
    assert re.findall(r'<span class="dot (\w+)"', strips[0]) == ["full", "untested"]


# One matrix cell's state, `title`, visible effect and statistic, in order.
_MATRIX_CELL: Final = re.compile(
    r'class="cell (corrected|nominal)" title="([^"]*)">(.*?)<br><span class="sub">([^<]*)</span>'
)


@pytest.mark.parametrize(
    ("row", "shown"),
    [
        pytest.param(
            _burden_row(effect=6.53, effect_bound=None, ci_low=2.1, ci_high=20.3),
            "OR 6.53",
            id="bounded: the interval is dropped",
        ),
        pytest.param(
            _burden_row(ci_low=0.786, pvalue=0.0316),
            "OR ∞ (CI 0.786–∞)",
            id="TBX5: unbounded, the lower bound kept",
        ),
        pytest.param(
            _burden_row(
                cohort_stratum="nonsyndromic",
                consequence_class="missense_damaging",
                n_cases=2405,
                ci_low=3.5212,
                pvalue=0.0025639,
            ),
            "OR ∞ (CI 3.52–∞)",
            id="ZFPM1: the lower bound at three significant figures",
        ),
        pytest.param(
            _burden_row(ci_low=None),
            "OR ∞",
            id="unbounded with no lower bound to keep",
        ),
    ],
)
def test_the_matrix_cell_drops_the_interval_but_never_the_measure(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts], row: BurdenRow, shown: str
) -> None:
    """The one property `_effect_compact` trades width for, and the one it does not.

    The matrix cell is a summary, so it drops the confidence interval: that
    qualifies the estimate rather than naming it, it is in the row's own table a
    screen below, and carrying it ran a cell to 61 characters and set the column
    width for the whole grid.

    **The measure is not negotiable.** `_effect`'s rule is that no branch omits
    it, because an odds ratio of 3.1 and a de novo enrichment of 3.1 are
    different claims and a bare `3.1` equates them. A mutation matrix on
    2026-08-05 found that dropping it from the compact form survived the entire
    suite, which is the only reason this test exists.

    **An unbounded effect keeps its lower bound** (2026-10-26). It rendered a
    bare `OR ∞`, and on TBX5's and ZEB2's cells "OR ∞" over "p 0.0316" stood for
    an interval of 0.786–∞ that spans 1, visible only on hover -- the reading
    `_effect` already refuses, "infinitely enriched" with the bound that
    contradicts it out of sight. The bare form is left only where there is no
    bound to keep. Asserted on the cell's visible text, sliced from its `title`,
    which carries the interval wherever the row has one.

    **The bound is formatted, and TBX5's cannot show it**: 0.786 prints the same
    raw as at three significant figures, so printing `ci_low` unformatted
    survived the full suite (review, 2026-10-26). ZFPM1's non-syndromic
    missense row, 3.5212, is a real published bound whose two forms differ.
    """
    page = _burden_page(tmp_path, facts_uncurated, [row])
    section = _burden_section_text(page)
    matrix = _slice_between(section, '<table class="matrix">', ("</table>",))

    [(_, title, visible, _)] = _MATRIX_CELL.findall(matrix)
    assert visible == shown
    # The interval is on hover, and spelled in full only there and in the table.
    interval = row.ci_low is not None
    assert ("95% CI" in title) is interval
    assert "95% CI" not in re.sub(r'title="[^"]*"', "", matrix)
    assert ("95% CI" in section[section.index("</table>") :]) is interval


def test_the_browse_headers_and_cells_line_up(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Each column holds what its header says it holds.

    **`data_table` zips headers to cells by position and cannot detect a
    mismatch**, so inserting a header without moving its cell silently heads one
    column for a field two places away. That happened: "burden across studies"
    (now headed "rare variants across studies", D49) went in after
    `definitive for` while the strip stayed after `atlas
    curation`, and the whole suite passed -- 774 tests, none of which read a
    cell against its own header.

    This project already has `test_the_consequence_column_is_headed_for_the_
    column_it_renders` for exactly this failure on the burden table. The browse
    table had no equivalent, which is why the defect survived a build, a review
    of the rendered page, and a push.

    Asserted by pairing each header with a value only that column can produce,
    so a swap of any two columns fails rather than merely a shift of one.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={TBX5: 4},
        concordance={
            gene: {"tested": 2, "enriched": 1, "corrected": 1, "families": []} for gene in facts_two
        },
    )
    page = _page(tmp_path, "index.html")

    headers = [
        re.sub(r"<[^>]+>", "", cell)
        for cell in re.findall(r"<th[^>]*>.*?</th>", page[: page.index("</thead>")])
    ]
    row = next(r for r in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S) if f">{TBX5}<" in r)
    cells = [re.sub(r"<[^>]+>", "", cell) for cell in re.findall(r"<td>(.*?)</td>", row, re.S)]

    assert len(headers) == len(cells), "a header was added without a cell, or the reverse"
    by_header = dict(zip(headers, cells, strict=True))

    assert by_header["gene"] == TBX5
    assert by_header["symbol"] == "TBX5"
    assert by_header["confidence"] == "definitive"
    assert "tested" in by_header["rare variants across studies"]
    assert "tpm" in by_header["heart expression"].lower() or by_header["heart expression"] == "—"


def test_the_gene_page_matrix_carries_its_own_key_and_caption(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**The browse page had these and the gene page did not.**

    That is the wrong way round. A reader arriving from a search lands on the
    gene page, and KDM6A's matrix there is entirely hollow -- four `not tested`
    and `no enrichment` cells -- beside a green `definitive` chip. Without the
    caption, that reads as the data contradicting the classification. It does
    not: KDM6A causes Kabuki syndrome, and burden tests at these cohort sizes
    routinely detect nothing for genes with overwhelming family evidence.

    Found by an adversarial review of the promotion diff, on the page where the
    omission would have done the damage rather than on the one where it would
    not have.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(consequence_class="lof")])
    section = page[page.index(_RARE_VARIANT_HEADING) :]

    # The caption that stops an empty matrix reading as a verdict.
    assert "not evidence against a gene" in section
    assert "rare-variant tests routinely detect nothing" in section
    assert "burden tests routinely" not in section
    # And a key, so the four states are nameable rather than only coloured.
    # Each swatch beside its own caption, not merely both strings present: a
    # membership check passes with the two swapped, and a swapped key on this
    # page tells a reader an untested cell was tested and found nothing. No
    # assertion named `cell-key` existed at all until review 2026-08-06.
    for css, caption in (
        ("cell-key no-enrichment", "tested, no enrichment detected"),
        ("cell-key not-tested", "not tested by that dataset"),
    ):
        pair = re.search(rf'<span class="{re.escape(css)}">[^<]*</span>\s*([^<]+)<', section)
        assert pair, f"the matrix key no longer draws {css!r} with a caption"
        assert pair.group(1).strip() == caption, (
            f"{css!r} is captioned {pair.group(1).strip()!r}, expected {caption!r}"
        )
    # It sits with the matrix, not at the far end of the page.
    assert section.index("</table>") < section.index("not evidence against a gene")
    # Cohort families, not studies: two papers on one sample collection are one
    # row, never two replications. Literal, beside the matrix it qualifies.
    assert "Rows are <strong>independent cohort families</strong>, not studies" in section
    assert section.index("</table>") < section.index("independent cohort families")


def test_a_family_wise_corrected_p_is_not_labelled_q(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`q` denotes a false-discovery quantity, and one of the two methods is not.

    The mirror carries `benjamini_hochberg`, which is an FDR, and
    `familywise_permutation`, which is not. Labelling both `q` made the matrix
    contradict the same page's own table, which names the second "family-wise" a
    screen below -- so the page said two different things about one number.

    Measured 2026-08-05: 29 cells rendered `q`, some over a family-wise p.
    """
    # Different *studies*, so the two land in different cells: one cell shows one
    # row, and putting both in the same one would render only the smaller p.
    fdr = _burden_row(
        consequence_class="lof",
        pvalue=0.001,
        pvalue_adjusted=0.01,
        pvalue_adjustment="benjamini_hochberg",
    )
    familywise = _burden_row(
        study="PMID:34324492",
        consequence_class="lof",
        pvalue=0.001,
        pvalue_adjusted=0.02,
        pvalue_adjustment="familywise_permutation",
        case_cohorts=("decipher",),
        control_cohorts=("gain_controls",),
    )
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [fdr, familywise],
        publications={_PUBLICATION.id: _PUBLICATION, **_OTHER_PUBLICATIONS},
    )
    section = page[page.index(_RARE_VARIANT_HEADING) :]
    matrix = section[section.index('<table class="matrix">') : section.index("</table>")]

    assert "q 0.01" in matrix
    assert "q 0.02" not in matrix
    assert "corrected p 0.02" in matrix


# --- The effect-and-uncertainty forest ------------------------------------
#
# Every fixture below carries **both** values of whatever it is testing.
# CLAUDE.md section 4.36 is why that is spelled out rather than assumed: a
# fixture whose rows all share the value under test measures nothing, and it
# has cost this repository four separate defects -- the last of them in a page
# fixture exactly like these.


def _forest_figures(section: str) -> list[str]:
    """Each `<figure class="forest">` on the page, whole, in document order."""
    return re.findall(r'<figure class="forest".*?</figure>', section, flags=re.S)


def test_two_effect_measures_never_share_one_axis(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The merge `effect_measure` exists to prevent, and the order it is drawn in.

    779 corpus rows carry an odds ratio and 75 an enrichment ratio. An odds
    ratio of 3.1 and a de novo enrichment of 3.1 are different claims.

    **The panels' order is a published fact and was pinned by nothing.**
    `pages.py` builds the facet list as `sorted({row.effect_measure ...})`
    over a `set`, and removing that `sorted` survived the whole suite: the
    assertion below reads `sorted(set(figures))`, which sorts its own actual
    value and therefore cannot see order at all.
    `test_build_is_reproducible.py` cannot see it either -- it builds twice in
    one process, where `PYTHONHASHSEED` is fixed for the interpreter's life --
    so two builds of one commit on two machines would publish 25 gene pages
    with their forest panels in opposite orders while every checksum verified.

    **The real vocabulary is far too small to fixture safely, and this is
    the case CLAUDE.md section 4.42 says a synthetic set is for.** A probe
    constructing the exact `set` under 300 explicit `PYTHONHASHSEED` values
    measured the published pair `{enrichment_ratio, odds_ratio}` iterating
    pre-sorted **149/300 (49.7%)** -- a coin flip, so a literal-order
    assertion on the real pair is not a guard. All three members of
    `EffectMeasure` measured 73/300 (24.3%), and there is no fourth to widen
    into. Eight synthetic tokens measured **0/300**, matching
    `test_prf008_message_lists_every_offending_unit_in_sorted_order`'s fix for
    the identical constraint. `_effect` resolves an unrecognised measure to
    its own raw token by design, so the tokens render rather than being
    silently dropped, and `_burden_row` builds a `BurdenRow` directly -- which
    vocabulary a token belongs to is `tables.py`'s check, not this one's.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1),
                _burden_row(
                    cohort_stratum="all",
                    comparator="mutation_model",
                    n_control_carriers=None,
                    n_controls=None,
                    control_cohorts=(),
                    expected_count=0.42,
                    effect=2.45,
                    effect_measure="enrichment_ratio",
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                    pvalue_test="poisson",
                ),
            ],
        )
    )

    figures = re.findall(r'data-effect-measure="([^"]+)"', section)
    assert sorted(set(figures)) == ["enrichment_ratio", "odds_ratio"]
    assert len(figures) == len(set(figures)), "one panel carried two measures"
    # And the two panels are two pictures, not one picture drawn twice: the
    # rows are otherwise identical, so a facet keyed on `study` alone would
    # emit a single figure holding both numbers on one axis. The enrichment
    # ratio draws as the de novo structured table (`_denovo_table`) and the
    # odds ratio as the forest, so a single forest figure is expected here.
    assert len(_forest_figures(section)) == 1
    assert len(re.findall(r'<figure class="denovo"', section)) == 1

    # The order, on a set wide enough that a coincidentally-sorted iteration
    # is not what is being measured. Asserted as the published sequence, not
    # as `sorted(...)` of it.
    synthetic = ["golf", "alpha", "hotel", "charlie", "echo", "bravo", "foxtrot", "delta"]
    ordered = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    consequence_class=token,
                    effect=2.45,
                    effect_measure=token,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                )
                for token in synthetic
            ],
        )
    )
    assert re.findall(r'data-effect-measure="([^"]+)"', ordered) == sorted(synthetic)


def test_a_study_publishing_no_interval_still_gets_a_panel(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Otherwise the only picture on TBX5's page is of the null study.

    PMID:40127276 reports de novo LOF enriched 297x, q 5.6e-08, surviving its
    own correction, and publishes no intervals. It now draws as the de novo
    structured consequence table (`_denovo_table`) rather than a forest, so
    the surviving finding must appear there as a row -- not be demoted to
    text, which is what a blanket "draw only where intervals exist" rule would
    have done to the forest's point estimate.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="lof",
                    origin="de_novo",
                    comparator="mutation_model",
                    n_control_carriers=None,
                    n_controls=None,
                    control_cohorts=(),
                    expected_count=0.0128,
                    effect=297.419440226872,
                    effect_measure="enrichment_ratio",
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                    pvalue=1.35e-09,
                    pvalue_test="poisson",
                    pvalue_adjusted=5.58e-08,
                    pvalue_adjustment="benjamini_hochberg",
                ),
            ],
            publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        )
    )

    # The panel is drawn, not skipped: the surviving finding is a table row.
    figures = re.findall(r'<figure class="denovo".*?</figure>', section, flags=re.S)
    assert len(figures) == 1
    assert "denovo-sig" in figures[0]
    # The point estimate and its correction are in the cells, not demoted to
    # prose: the enrichment value as `×`, and the surviving `q`.
    assert "297" in figures[0]
    assert "q 5.58e-08" in figures[0]
    assert "no interval published" not in figures[0]


def test_an_unbounded_effect_draws_an_arrow_and_never_a_ceiling(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`effect_bound: unbounded_above` publishes a null effect and null
    ci_high. Those rows carry the strongest signals, and `ci_low` is the
    finding -- never a blank, an em dash, or an invented ceiling."""
    section = _burden_section_text(_burden_page(tmp_path, facts_uncurated, [_burden_row()]))

    # **And the bar reaches the axis edge rather than a number nobody
    # published.** Asserted against literals in the figure's own coordinate
    # system, not against the constants that produced them (CLAUDE.md section
    # 4.38): 470.0 is the right-hand edge of the plot area and 26.0 is the
    # first row's centre line. A renderer that invented a ceiling -- say
    # `ci_high = ci_low * 10` -- would end the bar short of the edge.
    #
    # The class is `chart-arrow-open`, not `chart-arrow`: this study published
    # no correction for this row, and an arrow carries the same fill rule as a
    # circle. `"chart-arrow" in section` was the first spelling of the check
    # above it and is now a substring of the hollow class, so it is asserted
    # here as a whole attribute instead.
    figure = _forest_figures(section)[0]
    assert '<polygon class="chart-arrow-open" points="470.0,26.0 463.0,22.5 463.0,29.5"/>' in figure
    assert 'class="chart-arrow"' not in figure
    # **Scoped to this row's own bar, not to the figure.** `'x2="470.0"' in
    # figure` was the first spelling of this and a mutant giving the row a
    # ceiling of `ci_low * 10` survived it: the horizontal axis line runs to
    # the same edge and satisfied the check on the bar's behalf (CLAUDE.md
    # section 4.19). Matched at `y1="26.0"`, the first row's centre.
    bar = re.search(
        r'<line class="chart-line" x1="[\d.]+" y1="26\.0" x2="([\d.]+)" y2="26\.0"/>', figure
    )
    assert bar is not None, "the unbounded row drew no interval bar"
    assert bar.group(1) == "470.0", "the bar stopped short of the axis edge"
    # The number the study did publish is still the whole finding.
    assert "OR ∞ (95% CI 28.1–∞)" in section


def test_a_union_row_is_marked_as_one(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`damaging` is LOF + damaging missense; `all` is syndromic +
    nonsyndromic. Stacked as sibling rows they would read as independent
    findings, which the section's prose already denies."""
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                # A union on the consequence axis and on the stratum axis at
                # once; a union on the stratum axis alone; and a row that is
                # neither, so "tag every row" fails here too.
                _burden_row(
                    cohort_stratum="all",
                    consequence_class="damaging",
                    effect=2.45,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                ),
                _burden_row(
                    cohort_stratum="all",
                    consequence_class="lof",
                    effect=3.1,
                    effect_bound=None,
                    ci_low=1.4,
                    ci_high=9.0,
                ),
                _burden_row(
                    cohort_stratum="syndromic",
                    consequence_class="lof",
                    effect=4.2,
                    effect_bound=None,
                    ci_low=1.9,
                    ci_high=11.0,
                ),
                # A second study whose every row is `all cases` and which
                # publishes neither stratum below it, and no composite.
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="lof",
                    effect=14.96,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                ),
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="missense_damaging",
                    effect=2.24,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                ),
            ],
        )
    )

    assert "chart-union" in section
    # Two panels, in study order: the second study's `all cases` rows are the
    # whole cohort rather than the sum of two rows on the page, because it
    # publishes no syndromic or non-syndromic row at all -- which is exactly
    # PMID:40127276's shape on all 25 genes it covers. Tagging them `union`
    # would point a reader at rows that are not there, so the tag is
    # conditional on the components being in the same panel, as
    # `_composite_note` is conditional on the study reporting both.
    figures = _forest_figures(section)
    assert [figure.count('class="chart-union"') for figure in figures] == [0, 2]
    # The two components are named, so the tag is a claim a reader can check.
    assert "loss-of-function and damaging-missense rows together" in figures[1]
    assert "syndromic and non-syndromic rows together" in figures[1]


def test_the_union_key_never_promises_a_component_row_the_panel_does_not_draw(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_union_kinds` tests the components panel-wide; a tag is worn per row.

    Measured 2026-08-20 on the built corpus: 288 rows carry the tag, **0 of
    them with neither component on the panel** -- so the tag itself is never
    false -- but on 10 panels a tagged row has only one of its two components
    drawn. TBX5's Audain odds-ratio panel is one: `all cases · loss-of-
    function` is tagged and the key says "`all cases` is the syndromic and
    non-syndromic rows together", while the panel carries `syndromic · loss-
    of-function` and no non-syndromic one. All 10 are a cell the study
    published no row for at all -- a missing (stratum, consequence) cell is
    not a null result -- so the row below the union is part of it, not all of
    it, and the key said otherwise.

    The tag stays: `all cases · loss-of-function` really does contain the
    syndromic row beneath it, and dropping the tag would publish the two as
    independent findings, which is the defect the tag exists to prevent. It
    is the key's reconciliation instruction that has to become conditional.

    Three panels, because the corpus can only show one of the two axes: all
    10 measured cases split on the **stratum** axis, so a fixture built from
    the corpus alone would leave the consequence branch unmeasured.

    **The visible line names the row and what is missing from it** (review,
    2026-10-27). It said "at least one of the rows it contains is not
    plotted; the table below carries every row this study published for this
    gene", pointing a reader at the table for a row that is in no table: on
    all 10 published panels the missing component is a cell the study never
    published for the gene. A fourth panel separates the two things an
    absence can be -- a component the study published and the panel cannot
    draw (no effect to place) is "not plotted here", never "published no" --
    which no published panel shows today, so the fixture is its only
    measurement. Each sentence is pinned as a literal on its panel's visible
    line.
    """
    placed: dict[str, object] = {"effect_bound": None, "ci_low": 1.2, "ci_high": 8.1}
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                # Panel one, Audain odds ratios: `all cases · lof` is tagged
                # and both of its components are drawn. Nothing to warn about.
                _burden_row(cohort_stratum="all", effect=2.45, **placed),
                _burden_row(cohort_stratum="syndromic", effect=3.1, **placed),
                _burden_row(cohort_stratum="nonsyndromic", effect=4.2, **placed),
                # Panel two, a synthetic odds-ratio study: `syndromic ·
                # damaging` is tagged on the consequence axis with only one of
                # its two component consequences drawn. No published panel is
                # in this state; every measured case splits on the other axis.
                # (Until 2026-11 this fixture was an `enrichment_ratio`
                # mutation-model forest; enrichment now renders as the de novo
                # structured table, which has no union reconciliation prose,
                # so the consequence branch is measured on a forest instead.)
                #
                # **The missing component is on the panel under a different
                # stratum**, and no row here is `all cases`, so the stratum
                # union is off and this measures the consequence branch alone.
                # Without it, dropping `_missing_components`' `other.cohort_stratum
                # == row.cohort_stratum` -- checking the union across strata
                # rather than within one -- survived all 1,050 tests: this
                # panel's rows all shared a stratum, so holding it fixed or
                # not made no difference to any of them (CLAUDE.md 4.36).
                _burden_row(
                    study="PMID:99999999",
                    cohort_stratum="syndromic",
                    consequence_class="damaging",
                    effect=5.0,
                    **placed,
                ),
                _burden_row(
                    study="PMID:99999999",
                    cohort_stratum="syndromic",
                    consequence_class="lof",
                    effect=6.0,
                    **placed,
                ),
                _burden_row(
                    study="PMID:99999999",
                    cohort_stratum="nonsyndromic",
                    consequence_class="missense_damaging",
                    effect=4.0,
                    **placed,
                ),
                # Panel three, Sierant odds ratios: TBX5's real shape -- a
                # stratum union with the non-syndromic side missing.
                _burden_row(study="PMID:40127276", cohort_stratum="all", effect=7.0, **placed),
                _burden_row(
                    study="PMID:40127276", cohort_stratum="syndromic", effect=8.0, **placed
                ),
                # Panel zero, PMID:34324492 odds ratios: the non-syndromic side
                # is published and has no effect to place, so it is in this
                # study's table and on no panel.
                _burden_row(study="PMID:34324492", cohort_stratum="all", effect=2.0, **placed),
                _burden_row(
                    study="PMID:34324492", cohort_stratum="syndromic", effect=3.0, **placed
                ),
                _burden_row(
                    study="PMID:34324492",
                    cohort_stratum="nonsyndromic",
                    effect=None,
                    effect_measure=None,
                    effect_bound=None,
                    ci_low=None,
                ),
            ],
            publications={
                _PUBLICATION.id: _PUBLICATION,
                **_SIERANT_PUBLICATION,
                **_OTHER_PUBLICATIONS,
            },
        )
    )

    # Studies in id order; each study now carries a single odds-ratio forest
    # (the enrichment panel that used to sit beside PMID:40127276 now renders
    # as the de novo table), and the synthetic consequence-union study sorts
    # last.
    unplotted, stratum_short, complete, consequence_short = _forest_figures(section)

    def visible(figure: str) -> str:
        return figure[figure.index("<figcaption><p>") : figure.index("</p>")]

    assert visible(unplotted).endswith(
        " The <strong>all cases · loss-of-function</strong> row is a union, and the "
        "non-syndromic · loss-of-function row it contains is not plotted here."
    )
    assert visible(consequence_short).endswith(
        " The <strong>syndromic · damaging (LOF + missense)</strong> row is a union, and this "
        "study published no syndromic · missense (damaging) row for this gene."
    )
    assert visible(stratum_short).endswith(
        " The <strong>all cases · loss-of-function</strong> row is a union, and this study "
        "published no non-syndromic · loss-of-function row for this gene."
    )
    assert "row is a union" not in complete
    # The tag itself is untouched on every one of them: the row really does
    # contain the row below it, and saying otherwise would publish two
    # overlapping rows as independent findings.
    assert all(
        figure.count('class="chart-union"') >= 1
        for figure in (unplotted, consequence_short, stratum_short, complete)
    )


def test_an_arrow_carries_the_same_correction_fill_as_a_circle(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """An arrow is a marker, so its fill has to mean what a circle's means.

    `_forest_estimate` bypassed `_mark_class` for both arrow branches and
    drew every one of them in `chart-arrow`, whose declaration was
    byte-identical to `chart-point`'s -- the fill that says "survived this
    study's own correction". Measured 2026-08-20 on the built corpus: **205
    of the 854 plotted marks were arrows and not one had survived any
    correction.** 30 were corrected and failed, with q as high as 1.0, and
    175 came from a study that published no correction at all; 13 of the
    corrected-and-failed sat on a ClinGen `definitive` gene.

    GATA4 was the worst of them and is the shape of the fixture below: three
    marks on one panel, all solid, no hollow mark anywhere to read them
    against, and the third was `p 1 · q 1`.

    The fixture carries all three fills in one panel, because a fixture whose
    rows share the value under test measures nothing (CLAUDE.md section 4.36)
    -- and the filled branch has **no row in the committed corpus at all**, so
    the surviving unbounded row here is the only thing that proves it is
    reachable rather than dead.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            # Unbounded above, no correction published: hollow.
            _burden_row(cohort_stratum="all", consequence_class="lof"),
            # Unbounded above and the study's own correction survived:
            # filled. No published row is in this state.
            _burden_row(
                cohort_stratum="syndromic",
                consequence_class="missense_damaging",
                pvalue_adjusted=0.0037,
                pvalue_adjustment="benjamini_hochberg",
            ),
            # An effect of exactly zero -- the other arrow branch, at the
            # opposite edge -- corrected and failed: hollow.
            _burden_row(
                cohort_stratum="nonsyndromic",
                consequence_class="lof",
                effect=0.0,
                effect_bound=None,
                ci_low=0.0,
                ci_high=3.2,
                pvalue_adjusted=1.0,
                pvalue_adjustment="benjamini_hochberg",
            ),
            # The study's own negative control, surviving and unbounded:
            # muted, never in the result colour.
            _burden_row(
                cohort_stratum="all",
                consequence_class="synonymous",
                pvalue_adjusted=0.0041,
                pvalue_adjustment="benjamini_hochberg",
            ),
        ],
    )

    figure = _forest_figures(_burden_section_text(page))[0]
    assert figure.count('class="chart-arrow-open"') == 2
    assert figure.count('class="chart-arrow"') == 1
    assert figure.count('class="chart-control-arrow"') == 1
    # And the hollow form is visibly hollow rather than a second name for the
    # same paint. Asserted on the published stylesheet -- the class is only
    # worth anything if it resolves to a different fill, and `chart-arrow`'s
    # declaration was byte-identical to `chart-point`'s for three releases.
    assert ".chart-arrow-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }" in page
    assert ".chart-control-arrow { fill: var(--muted); }" in page
    # A key that names only the hollow mark leaves the reader to infer the
    # complement, and for every arrow on the site that inference was wrong.
    assert "A <strong>filled</strong> marker is a row that <strong>survived</strong>" in figure


def test_each_fill_in_the_key_is_a_fill_the_reader_can_see(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_forest_caption` is conditional clause by clause; the fills were not.

    A legend entry for a glyph that is not on the panel sends a reader
    hunting for it, which is `_POOLING_NOTICE`'s recorded lesson (CLAUDE.md
    section 4.27). Measured 2026-08-20 on the built corpus: **11 of the 137
    panels printed the hollow clause with no hollow marker drawn**, across 10
    genes -- CHD7, GATA4, KMT2D, MYH7, NKX2-5, NODAL, RBFOX2, RNF40, WDR5,
    ZIC3, six of them ClinGen `definitive`. On 9 of the 11 the panel did draw
    arrows, so the key affirmatively told the reader that everything visible
    had failed its correction when the arrows carried no such claim.

    Both directions are asserted, from two panels that differ in exactly this
    -- an all-surviving panel and an all-failing one.
    """
    survived: dict[str, object] = {
        "pvalue_adjusted": 0.0037,
        "pvalue_adjustment": "benjamini_hochberg",
    }
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    cohort_stratum="all",
                    effect=2.45,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                    **survived,
                ),
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    effect=3.1,
                    effect_bound=None,
                    ci_low=1.4,
                    ci_high=9.0,
                ),
            ],
            publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        )
    )

    # PMID:40127276 sorts before PMID:42230622, so the failing panel is first.
    failing, surviving = _forest_figures(section)
    hollow = "A <strong>hollow</strong> marker is a row that <strong>did not survive</strong>"
    filled = "A <strong>filled</strong> marker is a row that <strong>survived</strong>"
    assert filled in surviving and hollow not in surviving
    assert hollow in failing and filled not in failing
    # The glyph each clause names really is the one drawn beside it.
    assert 'class="chart-point-open"' not in surviving
    assert 'class="chart-point"' not in failing


def test_marker_fill_encodes_survival_of_the_studys_own_correction(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Filled where the study's own correction was survived, hollow where not.

    Never a correction this atlas computed. A row whose study published no
    correction at all is hollow for the same honest reason -- nothing says it
    survived one.

    The fixture carries one surviving and one non-surviving row, because a
    fixture whose rows all share the value under test measures nothing.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                cohort_stratum="all",
                effect=2.45,
                effect_bound=None,
                ci_low=1.2,
                ci_high=8.1,
                pvalue_adjusted=0.0037,
                pvalue_adjustment="benjamini_hochberg",
            ),
            _burden_row(
                cohort_stratum="syndromic",
                effect=3.1,
                effect_bound=None,
                ci_low=1.4,
                ci_high=9.0,
                pvalue_adjusted=0.42,
                pvalue_adjustment="benjamini_hochberg",
            ),
            _burden_row(
                cohort_stratum="nonsyndromic",
                effect=4.2,
                effect_bound=None,
                ci_low=1.9,
                ci_high=11.0,
            ),
        ],
    )
    section = _burden_section_text(page)

    assert "chart-point-open" in section
    assert 'class="chart-point"' in section
    figure = _forest_figures(section)[0]
    # One survived; the other two did not, and one of those two is hollow for
    # the second honest reason -- its study published no correction at all.
    assert figure.count('class="chart-point"') == 1
    assert figure.count('class="chart-point-open"') == 2
    assert "did not survive" in figure
    # And the hollow form resolves to a visibly different paint, asserted on
    # the published stylesheet exactly as the arrow half asserts its own
    # (CLAUDE.md section 4.31: a guard on one of a pair is evidence about
    # one). `.chart-arrow-open` was pinned here and `.chart-point-open` was
    # not, in the release that cites this rule -- and pointing it at
    # `var(--link)` leaves all 604 hollow circles on the site indistinguishable
    # from the 250 filled ones.
    assert ".chart-point-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }" in page


def test_a_survivor_that_is_not_support_is_hollow_and_the_key_says_why(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The forest's fill reads the card's rule, and its key names no fill it does not draw.

    **A depleted row that cleared its correction is hollow (2026-10-26).** The
    fill read the corrected p alone, so this row would have been drawn filled
    while the card beside it said none survives a correction its study
    published: the card and the evidence matrix read
    `concordance.survives_own_correction`, which counts only support. Now every
    row but a synonymous one is filled by that rule. 0 published rows differed
    between the two rules that day, so this fixture is the only place the
    difference is drawn. The hollow clause, "did not survive", is false of
    this row, so the key says what it is instead -- and only that, because it
    is the panel's only hollow row.

    **And no filled clause where only the negative control survived.** A
    surviving synonymous row is drawn muted, never in the result colour, and a
    key promising a filled marker would send a reader hunting for one. Keying
    the clause on survival alone survived the release review's 634-test
    mutation run until this case (2026-10-26).
    """
    corrected: dict[str, object] = {"pvalue_adjustment": "benjamini_hochberg"}
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                cohort_stratum="all",
                effect=0.31,
                effect_bound=None,
                ci_low=0.12,
                ci_high=0.74,
                pvalue=0.001,
                pvalue_adjusted=0.012,
                **corrected,
            ),
            _burden_row(
                cohort_stratum="all",
                consequence_class="synonymous",
                effect=2.0,
                effect_bound=None,
                ci_low=1.3,
                ci_high=3.1,
                pvalue=0.0004,
                pvalue_adjusted=0.008,
                **corrected,
            ),
        ],
    )

    [figure] = _forest_figures(_burden_section_text(page))
    assert figure.count('class="chart-point-open"') == 1
    assert figure.count('class="chart-control-point"') == 1
    assert 'class="chart-point"' not in figure
    assert (
        "A row that <strong>survived</strong> its own study's correction without pointing to "
        "enrichment in cases is drawn <strong>hollow</strong>: it is not evidence for this "
        "gene." in figure
    )
    assert "<strong>did not survive</strong>" not in figure
    assert "A <strong>filled</strong> marker" not in figure
    assert "A <strong>muted</strong> row is one of this study's" in figure
    # The card on the same page reads the same rule.
    assert "none survives a correction its study published" in _card_variants(_rail_of(page))


def test_a_gene_with_no_plottable_row_gets_a_sentence_and_no_chart(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D42. Measured: KLF13, CFC1 and MYH11 carry only CNV rows with no effect
    measure at all; CRIPTO has no burden row. Four of 92."""
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    study="PMID:34324492",
                    variant_class="cnv_deletion",
                    cohort_stratum="all",
                    consequence_class="all_coding",
                    effect=None,
                    effect_measure=None,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                    pvalue=0.0068,
                    pvalue_test="binomial",
                    pvalue_adjusted=0.99,
                    pvalue_adjustment="familywise_permutation",
                ),
            ],
            publications={**_OTHER_PUBLICATIONS},
        )
    )

    assert "<svg" not in section
    assert "no effect measure" in section
    # And the figures it does have stay in plain sight rather than folding
    # behind a summary that promises a picture there is none of.
    assert "<table" in re.sub(r"<details.*?</details>", "", section, flags=re.S)


@pytest.mark.parametrize(
    ("measure", "label", "gloss"),
    [
        ("odds_ratio", "OR", "no association"),
        ("rate_ratio", "rate ratio", "no difference in rate"),
        ("hazard_ratio", "hazard_ratio", "no difference between the arms"),
    ],
)
def test_the_axis_key_glosses_the_null_line_in_its_own_panels_vocabulary(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    measure: str,
    label: str,
    gloss: str,
) -> None:
    """What 1 means depends on what is on the axis, and the key said one thing.

    `_forest_caption` made the measure label dynamic and then hardcoded the
    gloss as "no enrichment". Measured 2026-08-20 on the built site: 112 of
    the 137 panels, across 88 gene pages, read "Axis: **OR** ... the vertical
    line is **1**, no enrichment" -- and for an odds ratio 1 is *no
    association*. "No enrichment" is the de-novo-against-expectation
    vocabulary that belongs to the 25 `enrichment_ratio` panels. That is the
    conflation `effect_measure` exists to prevent (`_effect`'s docstring: an
    odds ratio of 3.1 and an enrichment of 3.1 are different claims), stated
    in the sentence whose whole job is telling the reader how to read the
    axis. The sibling `<title>` on the same figure already varied correctly.

    Each gloss is pinned against a **literal**, never against
    `_NULL_LINE_GLOSS`, per CLAUDE.md section 4.38: asserting the constant
    would pass for any rewording including a reintroduction of the single
    hardcoded gloss. `rate_ratio` is in `EffectMeasure` and carries 0 rows
    today; it is parametrised so the third vocabulary is written before a
    study needs it rather than after. The fourth case is an unrecognised
    token, which must still be glossed -- `_effect` has no branch that omits
    its label and this has none that omits its gloss -- with a phrase true of
    any ratio rather than one borrowed from a vocabulary the token may not
    belong to.

    Scoped to the `<figure class="forest">` and then to the axis sentence
    inside it (section 4.19): "no enrichment" also appears in this section's
    prose elsewhere, so a section-wide or even figure-wide assertion is
    answered by a sentence that is not the one under test.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    effect=2.45,
                    effect_measure=measure,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                )
            ],
        )
    )

    figure = _forest_figures(section)[0]
    # The axis line's own `<p>`, the caption's one visible line (D57). Until
    # 2026-10-26 the caption was one run of sentences and this matched up to
    # the first ". "; the axis sentence now ends at `</p>`, so that pattern
    # would run on into the folded key and check a slice it was not about.
    axis = re.search(r"<figcaption><p>(Axis: .*?)</p>", figure, re.S)
    assert axis is not None, "the panel prints no axis sentence"
    sentence = axis.group(1)

    assert f"<strong>{label}</strong>" in sentence
    assert f"the vertical line is <strong>1</strong>, {gloss}." in sentence
    for other in ("no association", "no enrichment", "no difference in rate"):
        if other != gloss:
            assert other not in sentence, f"the {measure} axis also claims 1 is {other!r}"


def test_the_forest_carries_no_summary_diamond(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D33: no pooled statistic, ever. These cohorts overlap -- DDD
    contributes cases to more than one cited paper -- so a pooled estimate
    would count the same children twice. The forest idiom invites a diamond
    and a reader expects one; its absence is deliberate and captioned.

    **Captioned in two places, and only one of them was guarded.** The
    figure's `<title>` is what a screen reader announces before any of the
    marks, and it makes the same claim in its own words. Measured 2026-08-20,
    rewriting it to announce "and a pooled summary" survived the whole suite:
    the assertion below is scoped to `<figure class="forest">`, which contains
    the caption as well, so the caption answered for the title. A reader using
    assistive technology would have been told this atlas publishes the one
    statistic D33 forbids -- and told it first.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [_burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1)],
        )
    )

    assert "chart-diamond" not in section
    assert "no pooled" in section
    # Scoped to the panel that invites the diamond. `_POOLING_NOTICE` already
    # carries "no pooled statistic across studies" as the matrix caption a
    # screen above, so the section-wide assertion above passes with or without
    # a word from this figure (CLAUDE.md section 4.19).
    figure = _forest_figures(section)[0]
    # The caption alone: the `<title>` below makes the same claim and answered
    # for the caption until 2026-10-26 (a deletion probe of the caption
    # survived every test).
    caption = re.sub(r"<title>.*?</title>", "", figure, flags=re.S)
    assert "A forest plot usually ends in a summary diamond" in caption
    assert "no pooled" in caption
    # "people", not "children" (2026-10-26): a shared collection may be an
    # adult control cohort, and `shared_cohorts` counts control collections too.
    assert (
        "two studies that share a sample collection describe partly the same people, and "
        "adding their rows would count them more than once." in caption
    )
    assert "children" not in caption
    # And the same scoping again, one level in: the caption and the title are
    # two sentences and a check over the figure is answered by either.
    announced = re.search(r"<title>(.*?)</title>", figure, re.S)
    assert announced is not None, "the figure announces nothing to a screen reader"
    assert "no pooled summary" in announced.group(1)
    assert "a pooled summary" not in announced.group(1).replace("no pooled summary", "")


def test_the_forest_key_folds_behind_one_visible_axis_line(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D57: a chart caption is one visible line plus a "How to read this chart" toggle.

    The key was one `<figcaption>` of up to twelve sentences under every panel
    -- 185 panels on 114 pages, measured on a build of 9f24190. The axis clause
    (the measure, its range, and what 1 means for that measure) stays in sight
    and every other clause folds but one. Most of them describe a mark in
    words that are the same wherever that mark is drawn; the not-fully-split
    clause does not -- it warns about this panel. It folded too, as the D57
    spec placed it (10 panels on 10 pages, measured 2026-10-26), and the owner
    moved it to the visible line the same day (E1). **Each clause
    keeps its condition**: the tests above that assert a clause renders only
    beside its mark slice the figure, which contains the toggle, so they hold
    unchanged; this one asserts where the clauses sit.

    One panel that draws every mark, so every conditional clause renders: a
    surviving composite (filled), an arrow to the right (unbounded, its
    correction failed: hollow), an arrow to the left (an effect of zero), a
    stratum union whose non-syndromic half is absent (the not-fully-split
    warning), and a synonymous row with no interval (hatched and muted).

    Both halves are pinned whole, as literals, so a clause that leaves the
    fold, drops out or changes fails here. That pins two clauses nothing
    pinned before (inventory finding 1) -- the left arrow and the muted row --
    and the muted clause changed on the fold: it ended "How to read these
    results says why", a link to the same section the toggle now closes on,
    so it names the row's role itself ("its negative control").

    **Two more changes after review (2026-10-26).** The muted clause was
    singular -- "The muted row is this study's synonymous row" -- while each of
    the 113 panels carrying it draws three, one per case group; it says "one
    of this study's synonymous rows". And the toggle's link read "How to read
    this atlas", the nav's text for the top of `methods.html`, while it led to
    one section; it names the section.
    """
    survived = {"pvalue_adjusted": 0.0037, "pvalue_adjustment": "benjamini_hochberg"}
    placed: dict[str, object] = {"effect_bound": None}
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                cohort_stratum="all",
                consequence_class="damaging",
                effect=2.45,
                ci_low=1.2,
                ci_high=8.1,
                **placed,
                **survived,
            ),
            _burden_row(
                cohort_stratum="all",
                consequence_class="lof",
                pvalue_adjusted=0.42,
                pvalue_adjustment="benjamini_hochberg",
            ),
            _burden_row(
                cohort_stratum="all",
                consequence_class="missense_damaging",
                n_case_carriers=0,
                effect=0.0,
                ci_low=0.0,
                ci_high=3.2,
                pvalue_adjusted=1.0,
                pvalue_adjustment="benjamini_hochberg",
                **placed,
            ),
            _burden_row(
                cohort_stratum="syndromic",
                consequence_class="lof",
                effect=3.1,
                ci_low=1.4,
                ci_high=9.0,
                **placed,
            ),
            _burden_row(
                cohort_stratum="all",
                consequence_class="synonymous",
                effect=1.1,
                ci_low=None,
                ci_high=None,
                pvalue_adjusted=0.9,
                pvalue_adjustment="benjamini_hochberg",
                **placed,
            ),
        ],
    )

    figure = _forest_figures(_burden_section_text(page))[0]
    caption = re.search(r"<figcaption>(.*?)</figcaption>", figure, re.S)
    assert caption is not None, "the forest carries no caption"
    visible, toggle, folded = caption.group(1).partition(
        '<details class="chart-notes"><summary>How to read this chart</summary>'
    )
    assert toggle, "the forest key is not folded behind a toggle"
    # The not-fully-split warning is on the visible line since the owner's
    # decision of 2026-10-26 (E1): it is particular to this panel. Since review
    # (2026-10-27) it names each unsplit row and each missing component, once,
    # in panel order -- here four of each, in its plural form.
    assert visible == (
        "<p>Axis: <strong>OR</strong> on a logarithmic scale, 0.562 to 50; the vertical "
        "line is <strong>1</strong>, no association. The <strong>all cases · damaging "
        "(LOF + missense)</strong>, <strong>all cases · loss-of-function</strong>, "
        "<strong>all cases · missense (damaging)</strong> and <strong>all cases · "
        "synonymous</strong> rows are unions, and this study published no syndromic · "
        "damaging (LOF + missense), non-syndromic · damaging (LOF + missense), "
        "non-syndromic · loss-of-function, syndromic · missense (damaging), non-syndromic · "
        "missense (damaging), syndromic · synonymous or non-syndromic · synonymous row for "
        "this gene.</p>"
    )
    assert folded == (
        "<p>A <strong>filled</strong> marker is a row that <strong>survived</strong> its own "
        "study's correction for multiple testing. "
        "A <strong>hollow</strong> marker is a row that <strong>did not survive</strong> its "
        "own study's correction for multiple testing, or whose study published no correction "
        "at all &mdash; the atlas computes none of its own. "
        "A forest plot usually ends in a summary diamond. There is <strong>no pooled</strong> "
        "estimate here and none is computed anywhere on this site: two studies that share a "
        "sample collection describe partly the same people, and adding their rows would "
        "count them more than once. "
        "A <strong>hatched</strong> bar means <strong>no interval published</strong> for "
        "that row, not a wide one. "
        "An arrow to the <strong>right</strong> is an effect unbounded above &mdash; no "
        "control carried the variant &mdash; so the study's lower bound is the whole "
        "finding and there is no point estimate to place. "
        "An arrow to the <strong>left</strong> is an effect of zero &mdash; no case "
        "carried the variant &mdash; which a logarithmic axis has no position for. "
        "A row tagged <strong>union</strong> contains at least one other row in this panel "
        "and is not an independent finding: <code>damaging</code> is the loss-of-function "
        "and damaging-missense rows together; <code>all cases</code> is the syndromic and "
        "non-syndromic rows together. "
        "A <strong>muted</strong> row is one of this study's <strong>synonymous</strong> "
        "rows, its negative control, drawn apart from the results it calibrates.</p>"
        '<p><a href="../methods.html#rare-variants">How to read rare-variant results</a></p>'
        "</details>"
    )


def test_a_forest_never_replaces_the_figures_it_summarises(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The evidence-loss mutant.

    A picture that took the table with it would remove every exact count,
    interval and p-value from the HTML -- this repository's characteristic
    defect. The table folds; it does not leave.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [_burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1)],
        )
    )

    assert "<svg" in section
    folded = re.search(r'<details class="study-figures">.*?</details>', section, re.S)
    assert folded is not None, "the study's table must be folded, not deleted"
    assert "<table" in folded.group(0)
    assert "5 carriers / 1,471" in folded.group(0)
    assert "OR 2.45 (95% CI 1.2–8.1)" in folded.group(0)
    # A section-wide `<details` check would not be enough: the how-this-study-
    # counted block is a `<details>` a few lines above, and the section's
    # reading notes were one until 2026-10-26 (CLAUDE.md section 4.19). Nor is
    # `"<table" not in unfolded`: the evidence matrix is a table too, and it is
    # outside the fold on purpose.
    # Asserted on a count only this study's own table renders -- the matrix
    # cell carries the effect and one statistic and no count at all.
    unfolded = re.sub(r"<details.*?</details>", "", section, flags=re.S)
    assert "5 carriers / 1,471" not in unfolded


def test_the_browse_page_refuses_a_gene_it_has_no_concordance_for(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """The guard `build_genes` has, applied to the layer that lacked it.

    `bundles._concordance_for` raises on a published gene the mapping does not
    cover, because defaulting it would publish `"tested": 0` -- byte-identical to
    a measured "no study reported this gene". The browse page did
    `(concordance or {}).get(gene, {})` and rendered exactly that zero, silently.
    `runner.py` says the two layers cannot disagree; the guard existed on one of
    them.

    `None` remains allowed and means "this build has no burden data at all",
    which is what most tests of this page pass. A *mapping* that omits a
    published gene is the error: it means the caller derived it over a different
    population.

    Found by an adversarial review of the promotion diff, and it survived the
    mutation matrix until this test existed.
    """
    emitter = Emitter(root=tmp_path)

    with pytest.raises(KeyError, match="no concordance derived for published gene"):
        build_gene_index_page(
            facts_two,
            emitter,
            symbols={GATA4: "GATA4", TBX5: "TBX5"},
            validity=validity_two,
            burden_counts={},
            concordance={GATA4: {"tested": 0, "enriched": 0, "corrected": 0, "families": []}},
        )


def _limited(gene: str) -> tuple[GeneFacts, GeneValidity]:
    """A gene whose only ClinGen record grades it `Limited`.

    **Not `Definitive`, and that is the whole point of the fixture.** Every other
    gene fixture in this file is graded `Definitive`, so "the headline's grade"
    and "`DEFINITIVE`" are the same value and a lookup keyed on the wrong one is
    invisible -- CLAUDE.md section 4.30, third recurrence. Measured 2026-08-06:
    three mutants survived the whole suite before this existed, including
    `_graded_diseases` keyed back on `DEFINITIVE` and the rail label pinned to
    the literal `definitive for`.

    `confidence_by_lesion_group` is emptied because `_facts` derives it from
    `groups`, which is empty here; leaving `_facts`' default would publish a map
    keyed on lesion groups the gene does not declare.
    """
    limited = ValidityRecord(
        source=ValiditySource.CLINGEN,
        classification=Classification.LIMITED,
        classification_term="Limited",
        disease="MONDO:0007732",
        disease_label="atrial septal defect",
        moi="AD",
        sop="SOP11",
        gcep="Congenital Heart Disease Gene Curation Expert Panel",
        report_url=None,
    )
    return (
        replace(
            _facts(gene, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=Classification.LIMITED,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=(limited,),
            state=ValidityState.EXPERT_CURATED,
            has_source_discordance=False,
        ),
    )


def _submitter_admitted(gene: str) -> tuple[GeneFacts, GeneValidity]:
    """A gene no panel graded, admitted by two GenCC submitters agreeing."""
    submissions = tuple(
        ValidityRecord(
            source=ValiditySource.GENCC,
            classification=Classification.LIMITED,
            classification_term="Limited",
            disease="MONDO:0007732",
            disease_label="atrial septal defect",
            moi="Autosomal dominant",
            report_url=None,
            submitter=submitter,
        )
        for submitter in ("Ambry Genetics", "G2P")
    )
    return (
        replace(
            _facts(gene, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=None,
            validity_state=ValidityState.SUBMITTER_CURATED,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=submissions,
            state=ValidityState.SUBMITTER_CURATED,
            has_source_discordance=False,
        ),
    )


def test_the_rail_names_the_grade_the_panel_gave_and_the_disease_it_gave_it_for(
    tmp_path: Path,
) -> None:
    """One claim in two halves, and before 2026-08-06 they could disagree.

    The rail's first row read `definitive for` unconditionally while
    `_graded_diseases` returned only `DEFINITIVE` records, so a `limited` gene
    rendered `definitive for | —`: a label naming a grade the panel did not give,
    over an em dash where the disease belonged. 69 of the 92 genes published
    today are in that position.

    Both halves now come from `fact.headline_confidence`, so this asserts the
    label, the disease and the chip class **as one string in one element**. A
    page-wide check on each separately would pass with the lookup keyed back on
    `DEFINITIVE`, because `render.grade_legend` put the word "definitive" on
    every page (CLAUDE.md section 4.19). It left the gene page for
    `methods.html#scope` on 2026-10-26 (D50); the slice stays, because nothing
    stops the next page-wide "definitive" from arriving.
    """
    fact, validity = _limited(TBX5)
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {TBX5: fact},
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_11604.html")
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None
    rail_html = rail.group()

    assert 'class="chip chip-limited">limited</span>' in rail_html
    assert "<dt>limited for</dt><dd>atrial septal defect</dd>" in rail_html
    assert "definitive for" not in rail_html
    assert _EM_DASH not in rail_html.split("<dl>")[1].split("</dd>")[0]


def test_a_gene_no_panel_graded_names_the_submitters_that_admitted_it(
    tmp_path: Path,
) -> None:
    """The other warrant, on the page. 16 of the 92 genes published are here.

    `headline_confidence` is `null`, so there is no grade to head the row with
    and no disease a panel named. The row becomes `admitted by` and carries the
    submitters, which is the only place on the page a reader can see *why* the
    gene is here at all -- and the chip must read `not classified` in the
    outline style rather than borrowing any rung of the grade ramp, because "no
    panel looked" is not "graded weakest".
    """
    fact, validity = _submitter_admitted(GATA4)
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {GATA4: fact},
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None
    rail_html = rail.group()

    assert 'class="chip chip-ungraded">not classified</span>' in rail_html
    assert "<dt>admitted by</dt><dd>Ambry Genetics; G2P</dd>" in rail_html
    for rung in ("definitive", "strong", "moderate", "limited"):
        assert f"chip-{rung}" not in rail_html


def test_both_page_kinds_leave_the_grade_key_to_methods_html(tmp_path: Path) -> None:
    """A colour ramp with no key is a ranking the reader has to guess at.

    Until 2026-08-06 the site needed no key: every published gene was
    `definitive` and the only chip was green. It now renders four rungs and an
    ungraded state, and `limited` -- a panel saying the case is not yet made --
    is the plurality at 43 of 92.

    Asserted on **both** page kinds, because a reader arriving from a search
    engine lands on the gene page and never sees the browse page. That is
    CLAUDE.md section 4.31 exactly: the browse page carried the evidence legend
    and the gene page did not, and the gene page is where the omission did the
    damage.

    **The label, the chip class and the gloss were asserted as one triple**,
    because checking only that seven phrases appeared somewhere on the page let
    `moderate` painted with `strong`'s chip survive all 796 tests (measured
    2026-08-06).

    **The key left both page kinds for `methods.html#scope` on 2026-10-26
    (D50, P3)**, and this test was named
    `test_both_page_kinds_gloss_every_grade_they_can_render` until then. The
    triple assertion moved with the key, to `test_build_methods.py`, and so did
    the reworded dispute sentence. Section 4.31 is still met, by a different
    route: both page kinds carry the scope
    short form, whose link names the key ("what each grade means"), and a gene
    page's card glosses its own headline grade in the key's words. What stays
    here is that neither page kind carries the key -- a copy left behind is a
    copy that drifts from the one on `methods.html` -- and the two retired
    denials, absent on both: a negative on published bytes is what fails when a
    deleted sentence comes back (section 4.35).
    """
    fact, validity = _limited(TBX5)
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        {TBX5: fact}, emitter, symbols={TBX5: "TBX5"}, validity={TBX5: validity}, burden_counts={}
    )
    build_gene_pages(
        {TBX5: fact},
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    for name in ("index.html", "HGNC_11604.html"):
        page = _page(tmp_path, name)
        assert 'class="grade-key"' not in page, f"{name} still carries the grade key"
        assert "What these classifications mean" not in page, f"{name} still folds a key"
        # Unqualified it was false for 5 published genes carrying an
        # out-of-scope ClinGen `Disputed`; and its in-scope half went false on
        # 2026-10-25, when six published RASopathy genes gained an in-scope
        # `Disputed` for Costello or CFC. Both spellings, as everywhere a
        # retired phrase is asserted absent.
        for text in (page, html.unescape(page)):
            assert "No gene here was disputed or refuted" not in text
            assert "not published here at all" not in text


# --- Developmental expression -----------------------------------------------
#
# Every fixture below builds `ExpressionProfile` TypedDicts directly, the same
# choice `_burden_row` makes for `BurdenRow`: this is a rendering test, and
# going through `build/profiles.py`'s real pipeline (a TSV on disk, a quantile
# grid, `gene_expression_profiles` itself) would test that module a second
# time rather than the English this one chooses for its output.


def _tissue_entry(
    tissue: str = "Heart",
    median: float = 100.0,
    unit: str = "rpkm",
    n_samples: int = 3,
    placement: Placement | None = None,
    not_placed_reason: str | None = None,
) -> TissueProfileEntry:
    return TissueProfileEntry(
        tissue=tissue,
        median_abundance=median,
        unit=unit,
        n_samples=n_samples,
        placement=placement,
        not_placed_reason=not_placed_reason,
    )


def _placement(
    median_percentile: int = 97,
    q25_percentile: int | None = 92,
    q75_percentile: int | None = 99,
    n_genes: int = 19842,
) -> Placement:
    return Placement(
        q25_percentile=q25_percentile,
        median_percentile=median_percentile,
        q75_percentile=q75_percentile,
        median_abundance=100.0,
        unit="rpkm",
        n_samples=3,
        n_genes=n_genes,
        method="lowest percentile of a tied breakpoint (bisect_left)",
    )


def _specificity(
    tau: float = 0.483,
    tissues: tuple[str, ...] = ("Brain", "Heart", "Liver"),
    highest_in: str | None = "Heart",
) -> Specificity:
    return Specificity(
        tau=tau,
        scale="log2(x+1)",
        method="tau (Yanai et al. 2005)",
        tissues=tissues,
        n_tissues=len(tissues),
        highest_in=highest_in,
        medians=dict.fromkeys(tissues, 10.0),
    )


def _phase_info(
    outcome: str | None = "matched",
    phase_ids: tuple[str, ...] = ("septation",),
    reason: str | None = None,
) -> PhaseInfo:
    return PhaseInfo(outcome=outcome, phase_ids=phase_ids, reason=reason)


def _stage_entry(
    stage: str | None = "7wpc",
    phase: PhaseInfo | None = None,
    specificity: Specificity | None = None,
    specificity_unavailable_reason: str | None = None,
    tissues: tuple[TissueProfileEntry, ...] = (),
) -> StageProfileEntry:
    return StageProfileEntry(
        stage=stage,
        phase=phase if phase is not None else _phase_info(),
        specificity=specificity,
        specificity_unavailable_reason=specificity_unavailable_reason,
        tissues=list(tissues),
    )


def _dataset_profile_entry(
    dataset: str = "E-MTAB-6814",
    quantile_shard: str | None = "omics/profile_quantiles/E-MTAB-6814.json",
    stages: tuple[StageProfileEntry, ...] = (),
) -> DatasetProfileEntry:
    return DatasetProfileEntry(dataset=dataset, quantile_shard=quantile_shard, stages=list(stages))


def _expression_profile(datasets: tuple[DatasetProfileEntry, ...] = ()) -> ExpressionProfile:
    return ExpressionProfile(datasets=list(datasets))


def _profile_dataset(
    accession: str = "E-MTAB-6814",
    *,
    cardiac_tissues: tuple[str, ...] = ("Heart",),
    detection_floor: float = 1.0,
) -> Dataset:
    """A minimal, validly-constructed `Dataset`, mirroring `test_build_profiles.

    py`'s own `_dataset()` factory field for field, so this fixture is proven
    to satisfy `Dataset`'s `a_profile_dataset_is_fully_declared` validator
    rather than guessed at independently.
    """
    return Dataset(
        id=accession,
        archive="arrayexpress",
        technology="bulk_rnaseq",
        design="profile",
        tissue="whole embryo",
        developmental_stage="embryonic",
        organism="NCBITaxon:9606",
        n_samples=3,
        licence="CC BY 4.0",
        contrasts=[],
        cardiac_tissues=cardiac_tissues,
        detection_floor=detection_floor,
        floor_source="source methods, section 4",
        quantile_estimator="linear",
        stages=(Stage(token="7wpc", wpc=7.0, order=1),),
    )


def _cardiac_phase(
    phase_id: str = "atrial_septum_morphogenesis",
    *,
    label: str = "Atrial septum morphogenesis",
    end_wpc: float | None = 6.29,
) -> CardiacPhase:
    """One curated phase, `end_basis` derived from whether an end is given.

    Derived rather than passed, so a fixture cannot construct the state
    `CardiacPhase.end_fields_match_end_basis` forbids and then be read as
    evidence about a state the real vocabulary can hold. `end_wpc=None` is the
    `NOT_STATED` case -- `heart_looping` in `curation/cardiac_phases.yaml`,
    the one phase whose end the source never states.
    """
    stated = end_wpc is not None
    return CardiacPhase(
        id=phase_id,
        go_id="GO:0060413",
        label=label,
        start_wpc=3.71,
        end_wpc=end_wpc,
        start_carnegie_stage="CS12",
        end_carnegie_stage="CS17" if stated else None,
        start_hsapdv_id="HsapDv:0000019",
        end_hsapdv_id="HsapDv:0000024" if stated else None,
        end_basis=EndBasis.STATED if stated else EndBasis.NOT_STATED,
    )


def _cardiac_phases(*phases: CardiacPhase) -> CardiacPhaseFile:
    """The curated phase vocabulary a page bands its trajectory with."""
    return CardiacPhaseFile(
        attributed_to="Buijtendijk MFJ, Barnett P, van den Hoff MJB (2020)",
        citation="PMID:32048790",
        phases=list(phases or (_cardiac_phase(),)),
    )


def _expression_page(
    tmp_path: Path,
    profiles: dict[str, ExpressionProfile],
    datasets: dict[str, Dataset] | None = None,
    facts: dict[str, GeneFacts] | None = None,
    name: str = "HGNC_4173.html",
    phases: CardiacPhaseFile | None = None,
    cell_types: dict[str, CellTypeExpression] | None = None,
) -> str:
    """One gene page, built only to exercise the expression section.

    Every other input is the smallest fixture `build_gene_pages` accepts: no
    burden, no curated assertions, one mirrored validity record. `facts`
    defaults to a single uncurated GATA4 so a caller testing one gene need not
    build a `GeneFacts` by hand.

    `phases` defaults to `None` -- the corpus with no curated phase
    vocabulary at all -- rather than to `_cardiac_phases()`, so a test that
    wants a banded trajectory has to ask for one. Every expression test
    written before the trajectory landed therefore still exercises exactly
    what it exercised then.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts or {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED)},
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity={GATA4: _validity(), TBX5: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        profiles=profiles,
        datasets=datasets or {},
        phases=phases,
        cell_types=cell_types or {},
    )
    return _page(tmp_path, name)


def _slice_between(text: str, start_marker: str, end_markers: tuple[str, ...]) -> str:
    """`text` after `start_marker`, up to the first of `end_markers` (or the end)."""
    start = text.index(start_marker) + len(start_marker)
    rest = text[start:]
    positions = [rest.index(marker) for marker in end_markers if marker in rest]
    return rest[: min(positions)] if positions else rest


def _expression_section_text(page: str) -> str:
    """The `Developmental expression` section alone, up to the next `<h2>`.

    Sliced for the reason `_validity_table` is sliced: a page-wide assertion
    cannot tell this section's own claim from an unrelated occurrence of the
    same words elsewhere on the page (CLAUDE.md section 4.19) -- the research
    notice, the scope rule and the burden section all sit on the same page.
    """
    return _slice_between(page, "<h2>Developmental expression</h2>", ("<h2",))


def test_a_gene_with_no_expression_data_says_the_atlas_has_not_curated_it(
    tmp_path: Path,
) -> None:
    """`{"datasets": []}` -- `EMPTY_EXPRESSION_PROFILE` itself, not a look-alike
    this test built by hand, so the page and `bundles.py`'s fallback are proven
    to agree on what "empty" means.

    **Chosen wording, and the rejected one.** "The atlas has not yet curated a
    profile for this gene" is the true claim: no dataset mentions this gene at
    all, so there is no specific dataset in view to say "does not cover it"
    about, and that alternative would invent a check this atlas never ran. The
    section still renders -- never an absent one, the reason the uncurated
    lesion-assertion notice existed (the card's curation sentence since D55):
    an absent section cannot be told apart from "the atlas looked and found
    nothing".
    """
    page = _expression_page(tmp_path, {GATA4: EMPTY_EXPRESSION_PROFILE})
    section = _expression_section_text(page)

    assert "not yet curated" in section
    assert "does not cover" not in section
    assert "<h3>" not in section, "an empty profile must name no dataset"
    assert "<table" not in section


# --- The stage tables (the fold under the bulk figure) ------------------------
#
# Two kinds of table sit in the fold: one per cardiac tissue (a row for every
# stage) and one for the other organs. Every assertion below is sliced to one
# table and, where it is about a stage, to one `<tr>` of it (CLAUDE.md 4.19):
# the stage labels repeat across tables, so a page-wide or fold-wide check
# could not tell which table carried a figure.


def _fold(page: str) -> str:
    """The first `stage-figures` fold of the expression section, summary included."""
    found = re.search(
        r'<details class="stage-figures">.*?</details>', _expression_section_text(page), re.S
    )
    assert found is not None, "the stage figures must be folded, not deleted"
    return found.group(0)


def _bulk_block(page: str) -> str:
    """The bulk dataset block: from its `<h3>` to the end of its stage-figures fold."""
    section = _expression_section_text(page)
    start = section.index("<h3>")
    return section[start : section.index("</details>", start) + len("</details>")]


def _stage_tables_by_label(page: str) -> dict[str, str]:
    """Each table of the fold keyed by its label (a tissue, "Other organs", "All stages")."""
    found = re.findall(
        r'<p class="method"><strong>([^<]*)</strong></p>'
        r'<div class="scroll">(<table.*?</table>)</div>',
        _fold(page),
        re.S,
    )
    return {label: table for label, table in found}


def _row_cells(table: str, stage: str) -> list[str]:
    """The cells of the one `<tr>` whose first cell is `stage`."""
    rows = re.findall(r"<tr><td>(.*?)</td>(.*?)</tr>", table, re.S)
    matching = [rest for label, rest in rows if label == stage]
    assert len(matching) == 1, f"expected one row for {stage!r}, found {len(matching)}"
    return [stage, *re.findall(r"<td>(.*?)</td>", matching[0], re.S)]


def _stage_page(tmp_path: Path, stages: tuple[StageProfileEntry, ...]) -> str:
    """A page whose one dataset declares Heart cardiac, so a heart table is drawn."""
    profile = _expression_profile((_dataset_profile_entry(stages=stages),))
    return _expression_page(tmp_path, {GATA4: profile}, {"E-MTAB-6814": _profile_dataset()})


def _headers(table: str) -> list[str]:
    return re.findall(r'<th scope="col">(.*?)</th>', table)


def test_every_figure_reaches_the_tables(tmp_path: Path) -> None:
    """Three stages by three organs, with a gap of every kind, and every figure found.

    Heart: placed at 4wpc, below the floor at 5wpc, a different gap at 6wpc.
    Liver: absent at 5wpc (a missing row, not a null). Brain: placed throughout.
    Each stage's row in its own table must carry its label, the median with its
    n, and the percentile (or the named reason), and the stage-level tau and
    highest-median cells must sit on the heart row. The Liver/Brain cells are
    asserted per `<tr>` of the other-organs table, so a figure printed in the
    wrong table or against the wrong stage fails.
    """
    gap = "no_quantile_grid"
    floor_reason = "below_detection_floor"
    stages = (
        _stage_entry(
            stage="4wpc",
            specificity=_specificity(tau=0.411, highest_in="Heart"),
            tissues=(
                _tissue_entry(
                    tissue="Heart", median=111.0, n_samples=3, placement=_placement(71, 60, 80)
                ),
                _tissue_entry(
                    tissue="Liver", median=22.0, n_samples=2, placement=_placement(33, 30, 36)
                ),
                _tissue_entry(
                    tissue="Brain", median=44.0, n_samples=4, placement=_placement(55, 50, 58)
                ),
            ),
        ),
        _stage_entry(
            stage="5wpc",
            specificity=_specificity(tau=0.222, highest_in="Brain"),
            tissues=(
                _tissue_entry(
                    tissue="Heart", median=0.3, placement=None, not_placed_reason=floor_reason
                ),
                _tissue_entry(
                    tissue="Brain", median=66.0, n_samples=5, placement=_placement(88, 80, 90)
                ),
            ),
        ),
        _stage_entry(
            stage="6wpc",
            specificity_unavailable_reason="undefined",
            tissues=(
                _tissue_entry(tissue="Heart", median=77.0, placement=None, not_placed_reason=gap),
                _tissue_entry(
                    tissue="Liver", median=0.4, placement=None, not_placed_reason=floor_reason
                ),
                _tissue_entry(tissue="Brain", median=99.0, placement=None, not_placed_reason=gap),
            ),
        ),
    )
    page = _expression_page(
        tmp_path,
        {GATA4: _expression_profile((_dataset_profile_entry(stages=stages),))},
        {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)},
    )
    tables = _stage_tables_by_label(page)
    assert list(tables) == ["Heart", "Other organs"]
    heart, other = tables["Heart"], tables["Other organs"]

    assert _headers(heart) == [
        "stage",
        "phases",
        "median rpkm (n)",
        "percentile (IQR)",
        "τ (log2(x+1))",
        "highest median",
    ]
    assert _headers(other) == ["stage", "Brain", "Liver"]

    four = _row_cells(heart, "4wpc")
    assert four[1:] == [
        "septation",
        "111 rpkm (n=3 samples)",
        "71 (IQR 60–80)",
        "0.411",
        "Heart",
    ]
    five = _row_cells(heart, "5wpc")
    assert five[1:] == [
        "septation",
        "0.3 rpkm (n=3 samples)",
        "below the 1 rpkm floor",
        "0.222",
        "Brain",
    ]
    six = _row_cells(heart, "6wpc")
    assert six[1:] == [
        "septation",
        "77 rpkm (n=3 samples)",
        "no complete percentile grid is published for this organ at this stage",
        "not available: it could not be computed for this gene at this stage",
        "–",
    ]

    assert _row_cells(other, "4wpc")[1:] == [
        "44 rpkm (n=4 samples); 55 (IQR 50–58)",
        "22 rpkm (n=2 samples); 33 (IQR 30–36)",
    ]
    assert _row_cells(other, "5wpc")[1:] == ["66 rpkm (n=5 samples); 88 (IQR 80–90)", "–"]
    assert _row_cells(other, "6wpc")[1:] == [
        "99 rpkm (n=3 samples); no complete percentile grid is published for this organ at "
        "this stage",
        "0.4 rpkm (n=3 samples); below the 1 rpkm floor",
    ]


def test_the_heart_table_keeps_a_row_for_an_unsampled_stage(tmp_path: Path) -> None:
    """Heart has no row at 5wpc. The row stays: `not sampled`, `–`, and still the
    stage's phases and tau, which belong to the stage and not to the tissue.
    """
    stages = (
        _stage_entry(
            stage="4wpc",
            specificity=_specificity(highest_in="Heart"),
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
        _stage_entry(
            stage="5wpc",
            phase=_phase_info(outcome="matched", phase_ids=("valve_formation",)),
            specificity=_specificity(tau=0.321, highest_in="Liver"),
            tissues=(_tissue_entry(tissue="Liver", placement=_placement()),),
        ),
    )
    tables = _stage_tables_by_label(_stage_page(tmp_path, stages))
    assert _row_cells(tables["Heart"], "5wpc") == [
        "5wpc",
        "valve formation",
        "not sampled",
        "–",
        "0.321",
        "Liver",
    ]


def test_a_dataset_with_no_cardiac_tissue_sampled_still_publishes_stage_facts(
    tmp_path: Path,
) -> None:
    """Phases, tau and the highest organ belong to a stage; with no heart table
    to carry them they would reach no byte of the page. One table of stages does.
    """
    stages = (
        _stage_entry(
            stage="4wpc",
            phase=_phase_info(outcome="matched", phase_ids=("a_b", "c_d")),
            specificity=_specificity(tau=0.5, highest_in="Liver"),
            tissues=(_tissue_entry(tissue="Liver", placement=_placement()),),
        ),
    )
    tables = _stage_tables_by_label(_stage_page(tmp_path, stages))
    assert list(tables) == ["All stages", "Other organs"]
    assert _headers(tables["All stages"]) == ["stage", "phases", "τ (log2(x+1))", "highest median"]
    assert _row_cells(tables["All stages"], "4wpc") == ["4wpc", "a b, c d", "0.5", "Liver"]


def test_the_percentile_caption_names_the_grid_once_when_uniform_and_per_cell_when_not(
    tmp_path: Path,
) -> None:
    """Uniform `n_genes` is said once in the caption and never in a cell; mixed
    values are said per cell and the caption names no population.
    """

    def fold_for(second_grid: int) -> str:
        stages = (
            _stage_entry(
                stage="4wpc",
                tissues=(
                    _tissue_entry(tissue="Heart", placement=_placement()),
                    _tissue_entry(tissue="Liver", placement=_placement()),
                ),
            ),
            _stage_entry(
                stage="5wpc",
                tissues=(
                    _tissue_entry(tissue="Heart", placement=_placement(n_genes=second_grid)),
                    _tissue_entry(tissue="Liver", placement=_placement(n_genes=second_grid)),
                ),
            ),
        )
        page = _expression_page(
            tmp_path / str(second_grid),
            {GATA4: _expression_profile((_dataset_profile_entry(stages=stages),))},
        )
        return _fold(page)

    (tmp_path / "19842").mkdir()
    (tmp_path / "100").mkdir()
    uniform = fold_for(19842)
    caption = uniform[: uniform.index("<p", uniform.index("<p") + 1)]
    assert (
        "A percentile places this gene among 19,842 genes in the same organ and stage, "
        "so it compares nothing across rows."
    ) in caption
    assert " of 19,842 genes" not in uniform

    mixed = fold_for(100)
    caption = mixed[: mixed.index("<p", mixed.index("<p") + 1)]
    assert (
        "A percentile compares this gene only with genes in the same organ and stage, "
        "never across rows."
    ) in caption
    assert "among" not in caption
    assert "97 (IQR 92–99) of 19,842 genes" in mixed
    assert "97 (IQR 92–99) of 100 genes" in mixed


def test_a_tied_peak_reads_tied_and_no_adjective_reaches_the_tables(tmp_path: Path) -> None:
    """`highest_in: null` reads `tied` on the heart row, beside a stage whose argmax is
    named; and the fold carries no adjective about tau (D39(c)).

    Both argmax cases are in one fixture so a gate that fired on "any tissue
    named" rather than on `highest_in is None` fails here. The retired wording
    ("heart-preferential") is asserted absent from the built bytes
    (CLAUDE.md 4.35), including the other organ's peak.
    """
    stages = (
        _stage_entry(
            stage="4wpc",
            specificity=_specificity(highest_in="Heart"),
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
        _stage_entry(
            stage="5wpc",
            specificity=_specificity(highest_in="Liver"),
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
        _stage_entry(
            stage="6wpc",
            specificity=_specificity(highest_in=None),
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
    )
    page = _stage_page(tmp_path, stages)
    heart = _stage_tables_by_label(page)["Heart"]
    assert _row_cells(heart, "4wpc")[-1] == "Heart"
    assert _row_cells(heart, "5wpc")[-1] == "Liver"
    assert _row_cells(heart, "6wpc")[-1] == "tied"
    fold = _fold(page).lower()
    for word in ("preferential", "specific", "broadly", "peaks in"):
        assert word not in fold


def test_mixed_tau_scales_refuse(tmp_path: Path) -> None:
    """The tau header names one scale, so two stages disagreeing must fail rather
    than publish a header false for some rows.
    """
    other = _specificity()
    other["scale"] = "linear"
    stages = (
        _stage_entry(stage="4wpc", specificity=_specificity()),
        _stage_entry(stage="5wpc", specificity=other),
    )
    with pytest.raises(ValueError, match="τ scales differ"):
        _stage_page(tmp_path, stages)


def test_the_old_per_stage_sentences_are_gone(tmp_path: Path) -> None:
    """Sliced to the bulk dataset block (its `<h3>` to the end of its fold), not the
    whole section: the cell-type block above it is not this layer's to constrain.
    """
    stages = (
        _stage_entry(
            stage="4wpc",
            specificity=_specificity(),
            tissues=(
                _tissue_entry(tissue="Heart", placement=_placement()),
                _tissue_entry(
                    tissue="Liver",
                    median=0.3,
                    placement=None,
                    not_placed_reason="below_detection_floor",
                ),
            ),
        ),
    )
    page = _stage_page(tmp_path, stages)
    block = _bulk_block(page)
    for retired in (
        "<h4>",
        "Developmental phase",
        "Expression peaks in",
        "below the detection floor in",
        "τ is not available",
        "τ =",
    ):
        assert retired not in block, retired


def test_a_stage_without_tau_names_its_reason_and_the_floor_with_its_unit(
    tmp_path: Path,
) -> None:
    """`specificity: null` beside a stage with a real tau, on one heart table.

    **And the floor's unit where the peak is below it (2026-10-26).** The old
    sentence read "(detection floor 1)" -- a number with no unit, directly above
    table cells quoting the same floor as "(detection floor 1 tpm)". The unit is
    the entry's own, so the fixture's `rpkm` (not the corpus's `tpm`) fails a
    hard-coded one. A dataset with no measurement to read a unit from renders the
    bare number rather than a dangling space.
    """
    below = "peak_below_detection_floor"
    unplaced = tuple(
        _tissue_entry(
            tissue=organ, median=0.3, placement=None, not_placed_reason="below_detection_floor"
        )
        for organ in ("Heart", "Liver")
    )
    stages = (
        _stage_entry(
            stage="one_organ",
            specificity_unavailable_reason="one_organ_sampled",
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
        _stage_entry(
            stage="two_organs",
            specificity=_specificity(tissues=("Heart", "Liver"), highest_in="Heart"),
            tissues=(_tissue_entry(tissue="Heart", placement=_placement()),),
        ),
        _stage_entry(stage="peak_below", specificity_unavailable_reason=below, tissues=unplaced),
    )
    page = _expression_page(
        tmp_path,
        {GATA4: _expression_profile((_dataset_profile_entry(stages=stages),))},
        {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)},
    )
    heart = _stage_tables_by_label(page)["Heart"]
    assert _row_cells(heart, "one_organ")[4:] == [
        "not available: only one organ was sampled at this stage, and τ needs at least two",
        "–",
    ]
    assert _row_cells(heart, "two_organs")[4:] == ["0.483", "Heart"]
    assert _row_cells(heart, "peak_below")[4:] == [
        "not available: every organ sampled is below the 1 rpkm floor",
        "–",
    ]

    no_unit = _stage_tables_by_label(
        _expression_page(
            tmp_path / "no_unit",
            {
                GATA4: _expression_profile(
                    (
                        _dataset_profile_entry(
                            stages=(
                                _stage_entry(stage="no_unit", specificity_unavailable_reason=below),
                            )
                        ),
                    )
                )
            },
            {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)},
        )
    )
    cells = _row_cells(no_unit["All stages"], "no_unit")
    assert cells[2] == "not available: every organ sampled is below the 1 floor"


def test_a_below_floor_tissue_names_the_floor_value_and_never_says_not_detected(
    tmp_path: Path,
) -> None:
    """A negative assertion on published bytes: grep the BUILT PAGE, never the source
    (`pages.py` splits copy across adjacent literals, CLAUDE.md 4.35). Both the
    heart cell and the other-organs cell, and a dataset with no declared floor.
    """
    stages = (
        _stage_entry(
            tissues=(
                _tissue_entry(
                    tissue="Heart",
                    median=0.3,
                    placement=None,
                    not_placed_reason="below_detection_floor",
                ),
                _tissue_entry(
                    tissue="Liver",
                    median=0.3,
                    placement=None,
                    not_placed_reason="below_detection_floor",
                ),
            )
        ),
    )
    profiles = {GATA4: _expression_profile((_dataset_profile_entry(stages=stages),))}
    page = _expression_page(tmp_path, profiles, {"E-MTAB-6814": _profile_dataset()})
    tables = _stage_tables_by_label(page)
    assert _row_cells(tables["Heart"], "7wpc")[3] == "below the 1 rpkm floor"
    assert (
        _row_cells(tables["Other organs"], "7wpc")[1]
        == "0.3 rpkm (n=3 samples); below the 1 rpkm floor"
    )
    assert "not detected" not in page.lower()

    # A floor the dataset never declared: the wording drops the number, not the state.
    unplaced = _tissue_entry(median=0.3, placement=None, not_placed_reason="below_detection_floor")
    assert _tissue_percentile(unplaced, None, with_grid=False) == "below the detection floor"


def test_n_samples_travels_with_every_figure_and_distinguishes_a_single_observation(
    tmp_path: Path,
) -> None:
    """The schema permits `n_samples` = 1, so "median" can be a single observation,
    and the page must say so -- and say something different for three. In both
    kinds of table.
    """
    stages = (
        _stage_entry(
            tissues=(
                _tissue_entry(tissue="Heart", n_samples=3, placement=_placement()),
                _tissue_entry(tissue="Liver", n_samples=1, placement=_placement()),
            )
        ),
    )
    tables = _stage_tables_by_label(_stage_page(tmp_path, stages))
    assert _row_cells(tables["Heart"], "7wpc")[2] == "100 rpkm (n=3 samples)"
    assert _row_cells(tables["Other organs"], "7wpc")[1].startswith("100 rpkm (n=1 sample);")


def test_every_stage_outcome_reads_as_a_state_never_a_bug_or_blank(tmp_path: Path) -> None:
    """All four `PhaseOutcome` members on one heart table, so a gene page is proven for
    a stage in each state -- including MATCHED with one phase, which the
    several-phases test deliberately does not also cover (a fixture where every
    MATCHED stage has the same number of phases could not tell "renders the one
    phase" from "renders `phase_ids[0]` and drops the rest"). The reason is
    rendered verbatim and escaped.
    """
    heart = (_tissue_entry(tissue="Heart", placement=_placement()),)
    stages = (
        _stage_entry(
            stage="s_matched",
            phase=_phase_info(outcome="matched", phase_ids=("outflow_tract_septation",)),
            tissues=heart,
        ),
        _stage_entry(
            stage="s_outside",
            phase=_phase_info(
                outcome="outside_window", phase_ids=(), reason="outside the curated window"
            ),
            tissues=heart,
        ),
        _stage_entry(
            stage="s_postnatal",
            phase=_phase_info(outcome="post_natal", phase_ids=(), reason="post-natal"),
            tissues=heart,
        ),
        _stage_entry(
            stage="s_undeclared",
            phase=_phase_info(
                outcome="undeclared", phase_ids=(), reason="stage not declared by this dataset"
            ),
            tissues=heart,
        ),
        _stage_entry(
            stage="s_empty", phase=_phase_info(outcome="matched", phase_ids=()), tissues=heart
        ),
        _stage_entry(
            stage="s_markup",
            phase=_phase_info(outcome="outside_window", phase_ids=(), reason="<b>&odd</b>"),
            tissues=heart,
        ),
        _stage_entry(
            stage="s_noreason",
            phase=_phase_info(outcome="outside_window", phase_ids=(), reason=None),
            tissues=heart,
        ),
    )
    table = _stage_tables_by_label(_stage_page(tmp_path, stages))["Heart"]
    expected = {
        "s_matched": "outflow tract septation",
        "s_outside": "outside the curated window",
        "s_postnatal": "post-natal",
        "s_undeclared": "stage not declared by this dataset",
        "s_empty": "unnamed phase",
        "s_markup": "&lt;b&gt;&amp;odd&lt;/b&gt;",
        "s_noreason": "not available",
    }
    for stage, text in expected.items():
        assert _row_cells(table, stage)[1] == text
    assert "<b>" not in table


def test_a_matched_stage_with_several_phases_names_every_one(tmp_path: Path) -> None:
    """A stage inside an overlap names ALL of its phases, not the first -- the
    unique killer of a renderer that reads only `phase_ids[0]`.

    Three phases, not two: with two, a truncating renderer showing one and one
    dropping the last of two could both pass a weaker assertion. Joined by ", ".
    """
    phase = _phase_info(
        outcome="matched",
        phase_ids=("atrial_septation", "outflow_tract_septation", "ventricular_septation"),
    )
    stages = (
        _stage_entry(phase=phase, tissues=(_tissue_entry(tissue="Heart", placement=_placement()),)),
    )
    table = _stage_tables_by_label(_stage_page(tmp_path, stages))["Heart"]
    assert _row_cells(table, "7wpc")[1] == (
        "atrial septation, outflow tract septation, ventricular septation"
    )


def test_the_bulk_dilution_caveat_is_unconditional_not_only_beside_bad_news(
    tmp_path: Path,
) -> None:
    """`_POOLING_NOTICE` is the recorded precedent: made conditional once, it
    read false in the position that motivated the fix. A caveat present only
    beside a below-floor tissue would read as an excuse for that gene's own
    numbers rather than as the general fact it is.

    **The short form since 2026-10-26 (D50, spec section 5).** The notice was
    95 words on every gene page; the gene page keeps the claim a reader acts
    on -- bulk tissue, and a low or absent figure is not evidence against the
    gene -- and links `methods.html#expression`, which carries the dilution
    argument in full. It is one of the ten caveats the spec keeps visible, so
    it is pinned whole, as a literal, outside every fold, and it is still
    unconditional.
    """
    all_detected = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(tissues=(_tissue_entry(placement=_placement()),)),)
            ),
        )
    )
    caveat = (
        '<p class="notice-inline"><strong>Every measurement below is from bulk '
        "tissue</strong>: a low or absent figure is not evidence a gene is unimportant to "
        'heart development. <a href="../methods.html#expression">How to read these '
        "figures</a>.</p>"
    )

    def unfolded(section: str) -> str:
        return re.sub(r"<details.*?</details>", "", section, flags=re.S)

    good_news_section = _expression_section_text(_expression_page(tmp_path, {GATA4: all_detected}))
    assert caveat in unfolded(good_news_section)

    below_floor = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(
                                median=0.3,
                                placement=None,
                                not_placed_reason="below_detection_floor",
                            ),
                        )
                    ),
                )
            ),
        )
    )
    bad_news_section = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: below_floor},
            {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)},
        )
    )
    assert caveat in unfolded(bad_news_section)

    # The combined banner replaces the bulk one on a gene with cell-type data;
    # it carries the same caveat as a clause, outside every fold.
    with_cells = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: below_floor},
            {"E-MTAB-6814": _profile_dataset(detection_floor=1.0), "hoc": cell_type_dataset()},
            cell_types={GATA4: CellTypeExpression(datasets=[cell_type_payload()])},
        )
    )
    assert (
        "a low or absent value is not evidence a gene is unimportant to heart development."
        in unfolded(with_cells)
    )


def test_the_expression_reading_notes_are_on_methods_html_and_not_on_the_gene_page(
    tmp_path: Path,
) -> None:
    """The section's reading notes left every gene page on 2026-10-26 (D50).

    They were folded under "How to read these figures" on all 119 gene pages:
    how a median is placed and why none is placed below the floor, and -- on
    the 118 pages laying out more than one placement, measured on a build of
    0fc451a -- rule 8, that percentiles are not comparable across organs or
    stages. Both teach how to read the figures rather than saying anything
    about this gene, and both are now on `methods.html#expression`, where
    `test_build_methods.py` pins them; the section's one-line bulk notice
    links there.

    This test was `test_percentiles_are_captioned_as_not_comparable_only_when_
    more_than_one_is_shown` until then. Its gate -- say rule 8 only where two
    placements sit side by side -- left with the sentence: stated once, in
    general, on a page with no figures, there is no position for it to be
    false in. The fixture is still the two-placement one the gate fired on,
    so a reinstated conditional fails here, and the denominator still travels
    with every percentile, which is the half of rule 8 that stays.

    Asserted absent from the **whole** page, in both apostrophe spellings, by
    a phrase unique to each moved paragraph (CLAUDE.md section 4.35), and the
    bulk notice's own full form by the phrases its short form dropped.
    """
    two_placements = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(tissue="Heart", placement=_placement()),
                            _tissue_entry(
                                tissue="Liver", placement=_placement(median_percentile=40)
                            ),
                        )
                    ),
                )
            ),
        )
    )
    page = _expression_page(tmp_path, {GATA4: two_placements})

    assert "97 (IQR 92–99)" in _fold(page)
    for text in (page, html.unescape(page)):
        for moved in (
            "<summary>How to read these figures</summary>",
            # The percentile notice; "on this page" was its deictic.
            "not comparable across organs",
            "elsewhere on this page",
            # The reading notes.
            "travels with every percentile",
            "the rank of an unreliable measurement",
            # The bulk notice in full: the short form keeps its claim and
            # drops the argument. Not "once diluted across the sample":
            # `_no_trajectory_sentence` says that of its own gene, rightly.
            "never a single cell type",
            "cardiac neural crest",
            "cushion mesenchyme",
        ):
            assert moved not in text, f"{moved!r} is back on the gene page"


def test_a_tissue_stage_and_dataset_token_carrying_markup_is_escaped(tmp_path: Path) -> None:
    """Tissue and stage tokens are mirrored third-party strings, and this is the
    one artifact kind where a `<` in an organ label is the difference between a
    string and a script tag.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                dataset="<i>E-MTAB-x</i>",
                quantile_shard=None,
                stages=(
                    _stage_entry(
                        stage="<b>7wpc</b>",
                        tissues=(
                            _tissue_entry(tissue="<script>xss</script>", placement=_placement()),
                        ),
                    ),
                ),
            ),
        )
    )
    page = _expression_page(tmp_path, {GATA4: profile})

    assert "<script>xss</script>" not in page
    assert "&lt;script&gt;xss&lt;/script&gt;" in page
    assert "<b>7wpc</b>" not in page
    assert "&lt;b&gt;7wpc&lt;/b&gt;" in page
    assert "<i>E-MTAB-x</i>" not in page
    assert "&lt;i&gt;E-MTAB-x&lt;/i&gt;" in page


def test_the_dataset_heading_links_to_its_own_percentile_grid(tmp_path: Path) -> None:
    """D39(b)'s auditability reaching a reader, not only a program: the same
    relative-path convention `_rail`'s "this gene as JSON" link uses.
    """
    profile = _expression_profile((_dataset_profile_entry(),))
    page = _expression_page(tmp_path, {GATA4: profile})

    assert '<a href="../omics/profile_quantiles/E-MTAB-6814.json">' in page


# --- The phase-banded trajectory (D43) --------------------------------------
#
# Every fixture below differs from every other in the one value the code under
# test reads: how many heart stages carry a percentile (three, one, none) and
# whether the phase a stage names has a curated end. CLAUDE.md section 4.36 is
# the reason that is spelled out here -- a page fixture whose genes all share
# the value under test measures nothing, and it has cost this repository four
# separate defects.

_HEART_DATASET = {"E-MTAB-6814": _profile_dataset(cardiac_tissues=("heart",))}


def _heart_series(
    *stages: tuple[str, float | None | EllipsisType | str],
    phase_ids: tuple[str, ...] = ("atrial_septum_morphogenesis",),
    unplaced_reason: str = "below_detection_floor",
) -> ExpressionProfile:
    """A heart series over `stages`, in the order given.

    Three states, because the page owes a reader three different sentences
    and `_tissue_medians`' docstring insists the last two "must not render the
    same way":

    * a `float` -- a median this dataset placed against its percentile grid;
    * `None` -- a stage sampled and **not placed**, for `unplaced_reason`,
      which `_trajectory` ticks on the axis;
    * `...` -- a stage this dataset has **no heart row for at all**, which
      gets no tick and nothing else, and is why the stage is still in
      `entry["stages"]` and still occupies an x position;
    * a `str` -- a stage not placed for *that* reason rather than
      `unplaced_reason`, so one figure can tick stages for two reasons.

    The order given is the published order -- `gene_expression_profiles` emits
    stages by `Stage.order`, so a fixture's own sequence is the trajectory's.

    `unplaced_reason` is a parameter rather than a constant because the reason
    changes what a page is entitled to say: below a floor is a low reading
    about the gene, and a missing percentile grid is a hole in the reference
    that says nothing about the gene at all.
    """

    def heart(median: float | None | str) -> TissueProfileEntry:
        if median is None or isinstance(median, str):
            reason = unplaced_reason if median is None else median
            return _tissue_entry(
                tissue="heart", median=0.3, unit="tpm", placement=None, not_placed_reason=reason
            )
        return _tissue_entry(
            tissue="heart",
            median=median,
            unit="tpm",
            placement=_placement(),
            not_placed_reason=None,
        )

    entries = [
        _stage_entry(
            stage=token,
            phase=_phase_info(phase_ids=phase_ids),
            tissues=() if median is ... else (heart(median),),
        )
        for token, median in stages
    ]
    return _expression_profile((_dataset_profile_entry(stages=tuple(entries)),))


def _dataset_lede(section: str) -> str:
    """One dataset block's opening: its chart, or the sentence instead of one.

    Sliced off the stage tables the way `_validity_table` is sliced, and for
    the same reason (CLAUDE.md section 4.19). `_BULK_TISSUE_NOTICE` sits
    above every dataset block, and until 2026-10-26 it contained "not evidence
    that", so a section-wide assertion on that phrase passed whether or not
    the sentence under test was rendered at all. Its short form (D50) says
    "not evidence a gene is unimportant", and the full notice, "not evidence
    that" included, is on `methods.html#expression`. The slice still starts
    below the notice, so a full notice pasted back onto the page cannot make
    a lede assertion pass on its words.

    It ends at the first `<details`. Every caller's fixture draws no chart, so
    that is the stage-figures fold; where a chart is drawn, it is the chart's
    own "How to read this chart" toggle (D57), and the lede is the chart and
    its one visible caption line.
    """
    lede = _slice_between(section, "</h3>", ("<details",))
    # The summary line (`_glance_line`) sits between the heading and the lede and
    # is pinned on its own; a lede assertion must not read its counts.
    if lede.startswith('<p class="glance">'):
        lede = lede[lede.index("</p>") + len("</p>") :]
    return lede


def _band_caption_text(section: str) -> str:
    """The figure's caption alone: its visible line and its toggle (D57).

    From the visible line's first words to the toggle's `</details>`, which
    comes before the stage-figures fold's. Scoped for the reason
    `_validity_table` is: `_percentile_cell` prints every gap clause the
    caption can, from the same map, in the stage table below, and the grid's
    legend prints them too, so a section-wide check on one is answered by
    another (CLAUDE.md section 4.19). The slice starts inside the first
    sentence, so it excludes the legend that precedes it.
    """
    return _slice_between(section, "per stage on a log scale", ("</details>",))


def test_every_chart_carries_an_atlas_specific_axis(tmp_path: Path) -> None:
    """D43: a chart earns its place by an axis the source's browser lacks.

    Replaces `test_the_section_never_renders_a_chart_only_text_and_tables`,
    which asserted the proxy "never an `<svg>`". The proxy forbade this
    trajectory -- which no external browser has -- while permitting the
    rendering that hid the alphabetical-ordering defect for three releases.

    A bare re-plot of the source's own curve, with no phase band, must fail
    here.

    **And the sentence that replaces the refused chart has to be true**, which
    is a second assertion and not the same one. There is no heart chart for
    three different reasons and this half exercises the third; a caller that
    routes all three into one sentence publishes, of a gene placed at every
    stage on the page below, that no measurement there is placed at all. That
    is what shipped until 2026-08-20: the refusal was guarded, the replacement
    was not.

    **The refusal draws no chart and still draws the grid** (2026-11-02): the
    figure outlives the heart chart, so "no `<svg>`" is no longer the
    assertion. What a refusal must not draw is the phase strip and the line.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<svg" in section
    assert 'class="chart-strip"' in section, (
        "a chart with no named phase is the re-plot D43 forbids"
    )
    assert 'class="chart-band"' in section
    assert 'class="chart-line"' in section

    # The same three measurements with no curated phase vocabulary behind them
    # are exactly the source's own curve, and are refused.
    bare = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=None)
        )
    )
    assert "chart-strip" not in bare
    assert "chart-line" not in bare
    assert "chart-band" not in bare
    assert 'class="heat-' in bare, "the grid is not the chart D43 refuses, and it still draws"
    # Every one of these three stages is placed, and the table below prints
    # its percentile. Neither of the other two refusals' sentences may be
    # published here.
    assert "no measurement there is placed" not in bare
    assert "detection floor at every" not in bare
    assert "curated cardiac phase" in bare, "the refusal has to name the reason it refused"
    # And it is never left dangling. The reason list is empty in this branch,
    # so a sentence built to end in one ends "&mdash; ." -- an em dash, a
    # space and a full stop, published on 55 pages by the mutant above.
    assert "&mdash; ." not in bare

    # The third fixture, and the one whose sentence was worst. A series with
    # a below-floor stage *between* two placed ones has exactly one recorded
    # gap reason, so a sentence keyed on the reason set alone reads it as
    # "every stage is below the floor" -- in bold, on a gene whose page
    # prints two percentiles immediately below it. "Nothing placed" and "one
    # gap, below the floor" are different facts and only the first earns that
    # sentence.
    mixed = _dataset_lede(
        _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _heart_series(("4wpc", 5.0), ("5wpc", None), ("6wpc", 275.0))},
                _HEART_DATASET,
                phases=None,
            )
        )
    )
    assert "chart-strip" not in mixed
    assert "chart-line" not in mixed
    assert "detection floor at every" not in mixed
    assert "no measurement there is placed" not in mixed


def _bulk_svg(section: str) -> str:
    """The bulk figure's one `<svg>`, out of the section it sits in.

    Sliced rather than asserted page-wide for the reason `_validity_table` is
    (CLAUDE.md section 4.19): the cell-type block draws `<svg>`s on the same
    page and carries its own text. The bulk figure is the 560-wide one, and
    there is exactly one per dataset.
    """
    found = re.findall(r'<svg class="chart" viewBox="0 0 560 \d+".*?</svg>', section, flags=re.S)
    assert len(found) == 1, f"expected one bulk figure, found {len(found)}"
    return found[0]


def _grid_cells(svg: str) -> list[str]:
    """Every grid cell's class, in drawing order: one entry per `<rect class="heat-...">`."""
    return re.findall(r'<rect class="(heat-[a-z0-9]+)"', svg)


def _organ_labels(svg: str) -> list[str]:
    """The grid's row labels, top to bottom.

    Sliced from the text just before the first grid cell: the phase strip's
    labels share `chart-tick-label` and, when flipped left of their bar,
    `text-anchor="end"`, and the stage axis's rotated labels share both too
    (but carry a `transform`). Only the grid's rows follow its first row label.
    """
    grid = svg[svg.rindex("<text", 0, svg.index('<rect class="heat-')) :]
    return re.findall(
        r'<text class="chart-tick-label(?: cardiac)?" x="[^"]*" y="[^"]*" '
        r'text-anchor="end">([^<]*)</text>',
        grid,
    )


def _strip_labels(svg: str) -> list[str]:
    """The phase strip's bar titles, top to bottom."""
    return re.findall(r'<rect class="chart-strip"[^>]*><title>([^<]*)</title>', svg)


def test_the_trajectory_says_what_its_vertical_axis_spans(tmp_path: Path) -> None:
    """A curve with no numbers on it is a shape, not a measurement.

    The trajectory contained **zero `<text>` elements** -- no ticks, no
    labels, no numbers -- and its vertical axis is fitted to this gene's own
    placed heart values, so every gene's curve filled the same 104 px whatever it
    spanned (104 viewBox units, unchanged by the 2026-11-02 redesign). Measured
    2026-08-20 on the built corpus: PKD1L1 runs 1 to 2
    tpm and TBX20 runs 1 to 526 tpm, and the two were drawn as the same
    full-height excursion; 8 of the 85 charted genes span under 5-fold and 39
    span 20-fold or more. The caption prints the range for that reason, and
    the grid below the chart (fixed bins, comparable across genes) is the
    figure that answers the question this axis cannot.

    Two genes with different spans, because a fixture whose genes share the
    value under test measures nothing (CLAUDE.md section 4.36): with one
    gene, a renderer that labelled the axis with any constant would pass.

    The numbers are `_fmt`'s, not `charts.coordinate`'s. `coordinate` is the
    *pixel* formatter -- it fixes to a tenth, so it would publish a median of
    275 tpm as `275.0` and one of 4,761 as `4761.0` beside `_fmt`'s `4,761`
    in the table below. The literals here are what `_fmt` produces and what
    `coordinate` does not.
    """
    wide = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )
    narrow = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 1.0), ("5wpc", 1.5), ("6wpc", 2.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )

    # Pinned as whole elements in the figure's own coordinate system: 60.0 is
    # the right edge of the margin `_PLOT_LEFT` reserves, 30.4 is the top of
    # the value scale and 134.4 its foot (one phase row above the plot puts its
    # top at 27.0), and `text-anchor="end"` is what keeps a four-digit label
    # off the plot.
    top = '<text class="chart-label" x="60.0" y="30.4" text-anchor="end">'
    foot = '<text class="chart-label" x="60.0" y="134.4" text-anchor="end">'
    assert f"{top}275</text>" in _bulk_svg(wide)
    assert f"{foot}5</text>" in _bulk_svg(wide)
    assert f"{top}2</text>" in _bulk_svg(narrow)
    assert f"{foot}1</text>" in _bulk_svg(narrow)

    # And the caption says the same span in words -- a reader who cannot see
    # the axis labels still gets it.
    assert "Median tpm in heart per stage on a log scale from 5 to 275, fitted to this gene" in wide
    assert "Median tpm in heart per stage on a log scale from 1 to 2, fitted to this gene" in narrow


def test_a_flat_heart_series_prints_the_axis_range_not_the_medians(tmp_path: Path) -> None:
    """Stated as the **axis's** range and not as the medians'.

    `_axis_bounds` widens a flat series by half a decade either side, because
    `LogScale` refuses a zero-span axis, and 3 of the 85 charted genes were
    flat when this was measured (2026-08-20): for TFAP2B, placed at 1 tpm and
    nowhere else, "from 1 to 1" would be wrong and "from 0.316 to 3.16"
    attributes to its medians a spread neither end of which anyone measured --
    which is why the sentence says the axis is *fitted*, and why `_FLAT_SERIES_SPAN`
    is the one place that widening lives.
    """
    flat = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 1.0), ("5wpc", 1.0), ("6wpc", 1.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )
    assert "per stage on a log scale from 0.316 to 3.16, fitted to this gene" in flat
    assert "from 1 to 1" not in flat


def test_the_heart_axis_is_fitted_to_the_heart_alone(tmp_path: Path) -> None:
    """A second organ's 5,000 must not stretch the heart's own axis.

    The axis used to span every other organ's dashed line, so a gene whose heart
    reads 10 to 100 beside a liver at 5,000 drew heart as a flat line along the
    floor. Heart and liver differ by two orders of magnitude here, so any axis
    that pools them cannot print "10 to 100" (CLAUDE.md section 4.36).
    """
    phase = _phase_info(phase_ids=("atrial_septum_morphogenesis",))
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=tuple(
                    _stage_entry(
                        stage=token,
                        phase=phase,
                        tissues=(
                            _tissue_entry(tissue="Heart", median=heart, placement=_placement()),
                            _tissue_entry(tissue="Liver", median=5000.0, placement=_placement()),
                        ),
                    )
                    for token, heart in (
                        ("4 week post conception", 10.0),
                        ("5 week post conception", 100.0),
                    )
                )
            ),
        )
    )
    section = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: profile},
            {"E-MTAB-6814": _profile_dataset(cardiac_tissues=("Heart",))},
            phases=_cardiac_phases(),
        )
    )
    assert "on a log scale from 10 to 100, fitted to this gene" in section
    assert "5,000 rpkm" not in _band_caption_text(section)
    assert "chart-context" not in section


def test_a_phase_the_source_never_ended_bands_nothing(tmp_path: Path) -> None:
    """`end_basis: not_stated` has no curated width, so it has no band.

    `heart_looping` is the live case: the review states when it starts and
    never when it stops, so `CardiacPhaseFile.phases_for` excludes it at every
    wpc. A band drawn to an invented edge would assert exactly the boundary
    that field exists to refuse -- and would do it in the one encoding a
    reader takes at a glance.

    The guard is here rather than trusted to `phases_for` upstream: a stage
    carrying an unended phase's id is a state the real pipeline never
    produces, which is why this fixture builds one by hand. CLAUDE.md section
    4.28 -- a guard added to one layer is not a guard.
    """
    vocabulary = _cardiac_phases(
        _cardiac_phase(),
        _cardiac_phase("heart_looping", label="Heart looping", end_wpc=None),
    )
    both = _heart_series(
        ("4wpc", 5.0),
        ("5wpc", 40.0),
        ("6wpc", 275.0),
        phase_ids=("atrial_septum_morphogenesis", "heart_looping"),
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: both}, _HEART_DATASET, phases=vocabulary)
    )
    assert section.count('class="chart-band"') == 1
    assert _strip_labels(_bulk_svg(section)) == ["Atrial septum morphogenesis"]
    assert "Heart looping" not in section

    # And with nothing else to band, there is no admissible chart at all.
    alone = _heart_series(
        ("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0), phase_ids=("heart_looping",)
    )
    unbanded = _expression_section_text(
        _expression_page(tmp_path, {GATA4: alone}, _HEART_DATASET, phases=vocabulary)
    )
    assert "chart-strip" not in unbanded
    assert "chart-line" not in unbanded
    assert "chart-band" not in unbanded
    assert "Heart looping" not in unbanded


def test_the_band_caption_names_only_the_breaks_this_figure_actually_draws(
    tmp_path: Path,
) -> None:
    """Two ways the heart chart's line can break, and the caption named one of
    them unconditionally -- including on the 55 pages that draw no tick at all.

    Measured 2026-08-20 on the built corpus: 85 pages carry a trajectory and
    every one of them said "the stage is ticked on the axis instead"; **55
    drew zero ticks**, and 55 of the 59 with a visibly broken line were among
    them. TBX5 is one -- two polylines, a visible gap, and the only sentence
    about gaps described a mark that is not on the page.

    The second cause was never named at all: this dataset has **no row** for
    an organ at some stages, identically for every gene (heart has none at
    `school age child` or `elderly`, forebrain none at `6 week post
    conception`, testis none at `neonate` or `school age child`), and that
    breaks the line with no tick under it. `_tissue_medians`' own docstring
    insists the two "must not render the same way"; in the picture they did.

    Four fixtures, one per corner: a fixture where every case shares the
    value under test measures nothing (CLAUDE.md section 4.36).

    **And a fifth, on the axis the four could not see.** A tick is drawn for
    every stage this dataset sampled and did not place, whatever the recorded
    reason; the clause explaining it was conditioned on `value is None` and
    so named the detection floor for all of them. A stage unplaced because no
    percentile grid was published is ticked and captioned as a low reading,
    three lines above a table cell reading "no complete percentile grid is
    published for this organ at this stage" -- the page contradicting itself,
    and in the direction that invents a measurement. `_heart_series` has
    carried `unplaced_reason` since the sentence tier was written and no test
    had pointed it at a caption.
    """

    def caption(
        *stages: tuple[str, float | None | EllipsisType | str],
        reason: str = "below_detection_floor",
    ) -> str:
        return _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _heart_series(*stages, unplaced_reason=reason)},
                _HEART_DATASET,
                phases=_cardiac_phases(),
            )
        )

    # The floor clause's words since E1 (2026-10-26), on the visible line.
    tick = "Ticks on the axis mark stages below this dataset's detection floor"
    no_row = "no row for this organ at that stage"

    # Ticked and unbroken: the floor clause, and nothing about a missing row.
    ticked = caption(("4wpc", None), ("5wpc", 5.0), ("6wpc", 40.0), ("7wpc", 275.0))
    assert tick in ticked
    assert no_row not in ticked

    # TBX5's shape: a line broken by a stage this dataset has no row for, and
    # not one tick anywhere. The floor clause described a mark that is absent.
    unsampled = caption(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", ...), ("7wpc", 275.0))
    assert tick not in unsampled
    assert no_row in unsampled
    # And it does not point back at a clause this figure does not carry.
    # "A break with no tick under it is different" beside no tick at all is
    # the same defect one sentence further on.
    assert "is different:" not in unsampled

    # Both causes at once -- the caption owes a reader both sentences.
    both = caption(
        ("4wpc", 5.0),
        ("5wpc", 40.0),
        ("6wpc", None),
        ("7wpc", 275.0),
        ("8wpc", ...),
        ("9wpc", 12.0),
    )
    assert tick in both
    assert no_row in both
    assert "is different:" in both

    # Neither: an unbroken, fully placed series says nothing about either.
    clean = caption(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    assert tick not in clean
    assert no_row not in clean

    # The same ticked shape, unplaced for a reason that is not the floor. The
    # mark is drawn either way, so the caption still owes the reader an
    # explanation of it -- it may not be the floor's.
    grid = caption(
        ("4wpc", None), ("5wpc", 5.0), ("6wpc", 40.0), ("7wpc", 275.0), reason="no_quantile_grid"
    )
    assert grid.count('class="chart-absent"') == 1, "the tick is drawn whatever the reason"
    # Scoped to the caption: `_percentile_cell` prints the same clause in the
    # stage table below, from the same map and deliberately so, and answers a
    # section-wide check on it whatever the caption said.
    caption_only = _band_caption_text(grid)
    assert tick not in caption_only
    assert "no complete percentile grid is published for this organ at this stage" in caption_only
    # And the mark is still accounted for rather than left on the axis with
    # nothing said about it, which is the defect the floor clause was written
    # for in the first place.
    assert "ticked" in caption_only

    # **Both reasons on one figure** (review, 2026-10-27). The visible line's
    # floor sentence is universal -- "Ticks on the axis mark stages below this
    # dataset's detection floor" -- and was false of the tick the missing grid
    # drew, while the clause saying a tick can mean something else sat in the
    # fold. With both present the visible sentence names both, and the fold
    # does not repeat the second. Latent on the committed corpus, where every
    # unplaced cell is below the floor.
    mixed = caption(("4wpc", None), ("5wpc", 5.0), ("6wpc", "no_quantile_grid"), ("7wpc", 275.0))
    assert mixed.count('class="chart-absent"') == 2, "one tick per unplaced stage"
    mixed_caption = _band_caption_text(mixed)
    shown, toggle, folded = mixed_caption.partition('<details class="chart-notes">')
    assert toggle, "the trajectory's key is not folded behind a toggle"
    assert shown.endswith(
        " Ticks on the axis, beneath the foot of the scale rather than at a value on it, mark "
        "stages below this dataset's detection floor (1 tpm) or where no complete percentile "
        "grid is published for this organ at this stage.</p>"
    )
    assert tick not in mixed_caption
    assert "percentile grid" not in folded

    # And the tick is placed off the bottom of the scale rather than at a
    # value on it. Measured: KIF20A's axis starts at 42 tpm against a 1 tpm
    # floor, so its 6 ticks sat where the axis reads 38.6 tpm; the words are
    # what stop a reader interpolating them.
    assert "beneath the foot of the scale rather than at a value on it" in ticked


def test_the_band_caption_folds_behind_one_visible_axis_line(tmp_path: Path) -> None:
    """D57, for the figure: its axes stay in sight and its key folds.

    The caption was one `<p class="method">` of four or five sentences under
    the trajectory on 112 pages (measured on a build of 0fc451a). The two axis
    sentences stay visible -- stage order rather than elapsed time, and a
    vertical axis fitted to this gene, with its range -- because they say what
    every point on the figure is. The floor, tick and no-row clauses and the
    shaded phases fold under "How to read this chart", which ends on a link to
    `methods.html#expression`, where the axes are stated once in general. The
    link reads "How to read expression figures", not the nav's "How to read
    this atlas", for the reason the forest's toggle names its section.
    **Each clause keeps its condition**: the tests above that assert a clause
    renders only beside its mark read the whole section, or the caption
    through `_band_caption_text`, which takes the toggle, so they hold
    unchanged; this one asserts where the clauses sit.

    The fixture is the one that earns every clause the floor allows at once: a
    below-floor stage between two placed ones (the floor clause) and a stage
    with no heart row (the no-row clause).

    **The floor clause is on the visible line since the owner's decision of
    2026-10-26 (E1)**, in a short form, because it is particular to this gene;
    `test_every_particular_caveat_survives_outside_the_fold` requires it there
    where it fires and nowhere where it does not. The no-row clause stays
    folded, and its first spelling reads "different" rather than "different
    again": the clause it contrasts with is now the line above the toggle, not
    the sentence before it in the same paragraph.

    Pinned whole, both halves, as literals.
    """
    section = _expression_section_text(
        _expression_page(
            tmp_path,
            {
                GATA4: _heart_series(
                    ("4wpc", 5.0),
                    ("5wpc", 40.0),
                    ("6wpc", None),
                    ("7wpc", 275.0),
                    ("8wpc", ...),
                    ("9wpc", 12.0),
                )
            },
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )
    start = section.index('<p class="method">Median ')
    caption = section[start : section.index("</details>", start) + len("</details>")]
    visible, toggle, folded = caption.partition(
        '<details class="chart-notes"><summary>How to read this chart</summary>'
    )

    assert toggle, "the figure's key is not folded behind a toggle"
    assert section.count("How to read this chart") == 1, "one fold per figure, not one per part"
    assert visible == (
        '<p class="method">Median tpm in heart per stage on a log scale from 5 to 275, fitted '
        "to this gene, so heights compare nothing across genes; stages are in curated order, "
        "evenly spaced, not elapsed time. Ticks on the axis mark stages below this dataset's "
        "detection floor (1 tpm), beneath the foot of the scale rather than at a value on "
        "it.</p>"
    )
    assert folded == (
        '<p class="method">A break in the line with no tick under it is different: this '
        "dataset has <strong>no row for this organ at that stage</strong>, so there is no "
        "measurement to place, low or otherwise. Phases above the chart: Buijtendijk MFJ, "
        "Barnett P, van den Hoff MJB (2020). The grid's shades are fixed bins, so they "
        "compare across genes.</p>"
        '<p class="method"><a href="../methods.html#expression">How to read expression '
        "figures</a></p></details>"
    )


def test_the_line_never_spans_a_stage_the_atlas_did_not_place(tmp_path: Path) -> None:
    """A below-floor stage breaks the line; it is not drawn straight over, is
    ticked rather than plotted at zero, and its tick never lands on a real
    measurement.

    Measured 2026-08-20 on the committed corpus: 20 of the 85 charted genes
    carry a below-floor stage strictly between two plotted points, and TBX1
    has 11 of them against 5 placed stages. A single line through the placed
    points runs smooth and high across exactly the stages where this dataset
    measured below its own floor -- the picture asserting a continuity the
    figures deny.

    The fixture has two runs and a lone leading point, so "one line through
    everything", "splits but drops a run" and "splits correctly" are three
    distinguishable outcomes rather than two.

    **Absorbs `test_a_stage_below_the_floor_is_an_axis_tick_and_never_a_zero`,
    whose fixture this one is a strict superset of.** Measured 2026-08-20:
    all three mutants that test could kill are killed here too, so it was
    paying for a build without being the unique killer of anything. Its
    rationale is the record and moves with it -- *below the floor*, *not
    sampled* and *zero* are three different claims: plotting a below-floor
    stage at zero puts it on the axis as the lowest measurement this atlas
    made, and leaving a gap makes it indistinguishable from a stage nobody
    sampled. The counts are asserted rather than the class name alone,
    because a renderer that ticked every stage, or that ticked the
    below-floor stage *and* plotted it, satisfies a bare `in` check.

    **And `_FLOOR_GAP` is what keeps the tick off the lowest real point.**
    That constant's whole stated purpose is to stop the atlas's weakest
    statement about a figure -- "we do not vouch for this one" -- landing on
    the same pixel row as its lowest real one. Nothing guarded it: the lowest
    placed point always sits at exactly the plot's bottom less `_FLOOR_GAP`,
    so ticking a below-floor stage there instead of on the axis puts the two
    marks on one row on every charted page -- cy 116.0 on HGNC:10249 before
    the phase strip moved the plot down (it is 131.0 against a tick at 139.0
    with one phase row above it) --
    which is the conflation `_percentile_cell` refuses in words.
    """
    profile = _heart_series(
        ("4wpc", 5.0),
        ("5wpc", None),
        ("6wpc", 40.0),
        ("7wpc", 275.0),
        ("8wpc", None),
        ("9wpc", 60.0),
        ("10wpc", 90.0),
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    drawn = re.findall(r'<polyline class="chart-line" points="([^"]+)"', section)
    assert len(drawn) == 2, "one line per run of adjacent placed stages"
    assert [len(run.split(" ")) for run in drawn] == [2, 2]
    # The stage that broke the line is still a point on the page, and the two
    # stages that broke it are still ticked.
    assert section.count('class="chart-point"') == 5
    assert section.count('class="chart-absent"') == 2

    # No tick shares a row with a measurement. Read off the published `cy`
    # rather than recomputed, so the assertion cannot agree with the renderer
    # by sharing its arithmetic.
    figure = _bulk_svg(section)
    rows = {
        css: {cy for cy in re.findall(rf'<circle class="{css}" cx="[^"]*" cy="([^"]*)"', figure)}
        for css in ("chart-point", "chart-absent")
    }
    assert rows["chart-absent"] == {"139.0"}, "a tick belongs on the axis"
    assert "131.0" in rows["chart-point"], "the lowest placed point sits one _FLOOR_GAP up"
    assert not rows["chart-absent"] & rows["chart-point"], (
        "a below-floor tick shares a row with a real measurement, which is the "
        "conflation _FLOOR_GAP exists to prevent"
    )


def test_a_gene_placed_at_one_stage_gets_markers_and_no_line(tmp_path: Path) -> None:
    """Tier 3. TFAP2B is placed at 1 of 19 stages in the committed corpus.

    A polyline through one point draws a trend the data does not contain.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", None), ("6wpc", None))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<polyline" not in section
    assert "chart-point" in section

    # Two adjacent placed stages is the case the constant actually gates, and
    # the one above cannot reach it: a run of a single point is dropped
    # whatever the threshold says, so `_STAGES_FOR_A_TRAJECTORY = 1` survives
    # a one-stage fixture. Two measurements draw a segment that claims
    # something about the interval between them.
    pair = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", None))
    boundary = _expression_section_text(
        _expression_page(tmp_path, {GATA4: pair}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<polyline" not in boundary
    assert boundary.count('class="chart-point"') == 2


@pytest.mark.parametrize(
    ("reason", "present", "absent"),
    [
        pytest.param(
            "below_detection_floor",
            ("detection floor at every", "not evidence that", "heart", "3 stages sampled"),
            (),
            id="below this dataset's own floor",
        ),
        pytest.param(
            "no_quantile_grid",
            ("no complete percentile grid is published for this organ at this stage",),
            ("detection floor at every",),
            id="a gap in the reference",
        ),
    ],
)
def test_nothing_placed_says_which_of_the_two_facts_it_is(
    tmp_path: Path, reason: str, present: tuple[str, ...], absent: tuple[str, ...]
) -> None:
    """D42, tier 4, and the reason it is not cosmetic.

    Measured 2026-08-20: seven published genes are below the floor in heart at
    every one of 19 stages -- SEMA3E, ZIC3, USP44, DAW1, FGF8, GDF1, NODAL.
    ZIC3 and NODAL are ClinGen definitive. An empty chart beside a definitive
    chip asserts the thing the bulk-dilution caveat denies.

    **"Nothing placed" is not one fact, and only one of them is about the
    gene.** `below_detection_floor` says this dataset measured the gene under
    its own floor, which is what earns the dilution argument beside it. Every
    other gap -- no percentile grid for this organ and stage, no floor
    declared, the dataset unregistered -- is a hole in the *reference*, and
    reporting it as a low reading asserts a measurement this dataset never
    made. The two cases were two tests with identical bodies differing only
    in this parameter, which the second's own docstring said; they are one
    parametrised test.

    Asserted against literals, not against the constant that produced them: a
    test that imports the string it asserts on compares the module to itself
    and passes for any rewording.

    Sliced to the dataset block's own lede before asserting, because
    `_BULK_TISSUE_NOTICE` carried "not evidence that" over every gene page
    on this site, charted or not, until its short form (2026-10-26); the slice
    stays, so the full notice pasted back cannot answer for this sentence.
    """
    profile = _heart_series(("4wpc", None), ("5wpc", None), ("6wpc", None), unplaced_reason=reason)
    section = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
        )
    )

    assert "chart-line" not in section
    assert "chart-strip" not in section
    assert 'class="heat-' in section, "the grid still draws where there is no heart chart"
    for phrase in present:
        assert phrase in section
    for phrase in absent:
        assert phrase not in section


def test_a_chart_never_replaces_the_figures_it_summarises(tmp_path: Path) -> None:
    """The evidence-loss mutant.

    A chart that removed the table would take every exact figure out of the
    HTML -- this repository's characteristic defect, curated work reaching no
    page. The table moves into `<details>`; it does not leave.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<svg" in section
    folded = re.search(r'<details class="stage-figures">.*?</details>', section, re.S)
    assert folded is not None, "the stage figures must be folded, not deleted"
    # Not merely present, and not merely somewhere on the page: every figure
    # the 21 blocks carried is inside the fold the chart summarises.
    assert "<table" in folded.group(0)
    assert "275 tpm (n=3 samples)" in folded.group(0)
    assert "97 (IQR 92–99)" in folded.group(0)
    assert "among 19,842 genes" in folded.group(0)
    # And the section-wide `<details` assertion this replaces is not enough on
    # its own: the trajectory's "How to read this chart" toggle renders one
    # a few lines above (D57; until 2026-10-26 it was the section's folded
    # reading notes), so a renderer that dropped the wrapper entirely would
    # still satisfy it (CLAUDE.md section 4.19).
    unfolded = re.sub(r"<details.*?</details>", "", section, flags=re.S)
    assert "<table" not in unfolded


# --- The organ grid and the stage axis (2026-11-02) ----------------------------
#
# The fixtures below use capitalised organs and rpkm, which the real corpus does
# not (tissues are lower-case, tpm), so a renderer that hard-coded either would
# fail here. `_HEART_PHASE` is the one phase id `_cardiac_phases()` names.

_HEART_PHASE = _phase_info(phase_ids=("atrial_septum_morphogenesis",))
_CAPITAL_HEART = {"E-MTAB-6814": _profile_dataset(cardiac_tissues=("Heart",))}


def _grid_page(
    tmp_path: Path,
    stages: tuple[tuple[str, tuple[TissueProfileEntry, ...]], ...],
    datasets: dict[str, Dataset] | None = None,
) -> str:
    """One dataset's expression section, for stages given as `(token, cells)`."""
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=tuple(
                    _stage_entry(stage=token, phase=_HEART_PHASE, tissues=cells)
                    for token, cells in stages
                )
            ),
        )
    )
    return _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: profile},
            _CAPITAL_HEART if datasets is None else datasets,
            phases=_cardiac_phases(),
        )
    )


def _placed(tissue: str, median: float) -> TissueProfileEntry:
    return _tissue_entry(tissue=tissue, median=median, placement=_placement())


def _unplaced(tissue: str, reason: str, median: float = 0.4) -> TissueProfileEntry:
    return _tissue_entry(tissue=tissue, median=median, placement=None, not_placed_reason=reason)


def test_every_stage_is_labelled_and_an_unsampled_one_is_greyed(tmp_path: Path) -> None:
    """The stage axis names every stage the dataset publishes, not only heart's.

    Heart has no row at the middle stage (the real dataset's shape at
    `school age child` and `elderly`), so the old axis, which labelled only the
    cardiac tissue's own stages, left that column unlabelled in a figure whose
    grid has a cell there. The prenatal run carries one shared
    "weeks post conception", centred, and the old axis title is gone.
    """
    section = _grid_page(
        tmp_path,
        (
            ("4 week post conception", (_placed("Heart", 10.0), _placed("Liver", 5.0))),
            ("5 week post conception", (_placed("Liver", 5.0),)),
            ("infant", (_placed("Heart", 20.0), _placed("Liver", 5.0))),
        ),
    )
    svg = _bulk_svg(section)

    assert ">4<" in svg
    assert ">infant<" in svg
    assert re.search(r'<text class="chart-tick-label unsampled"[^>]*>5</text>', svg)
    assert re.search(r'<text class="chart-tick-label"[^>]*>4</text>', svg)
    assert re.search(r'<text class="chart-tick-label"[^>]*>infant</text>', svg)
    assert svg.count("weeks post conception") == 1
    assert svg.count("curated order, not elapsed time)") == 0
    assert "developmental stage (curated" not in svg


@pytest.mark.parametrize(
    ("median", "placed", "reason", "expected"),
    [
        pytest.param(9.99, True, None, "heat-1", id="just under 10"),
        pytest.param(10.0, True, None, "heat-2", id="10 is the next bin"),
        pytest.param(100.0, True, None, "heat-3", id="100 is the next bin"),
        pytest.param(1000.0, True, None, "heat-4", id="1,000 is the top bin"),
        pytest.param(0.4, False, "below_detection_floor", "heat-floor", id="below the floor"),
        pytest.param(0.4, False, "no_quantile_grid", "heat-gap", id="any other reason"),
    ],
)
def test_each_grid_cell_takes_the_state_its_measurement_has(
    tmp_path: Path, median: float, placed: bool, reason: str | None, expected: str
) -> None:
    """Bins are `< 10`, `10-100`, `100-1,000`, `>= 1,000`, closed below and open above.

    One cell per page, so the class is asserted as the only grid cell there is
    and a renderer that classed every cell the same way cannot hide in a
    fixture whose cells all agree (CLAUDE.md section 4.36). An unplaced cell is
    never binned, whatever its median: 0.4 would be `heat-1` on the ramp.
    """
    cell = _placed("Heart", median) if placed else _unplaced("Heart", reason or "", median)
    svg = _bulk_svg(_grid_page(tmp_path, (("4 week post conception", (cell,)),)))
    assert _grid_cells(svg) == [expected]


def test_the_grid_draws_one_cell_per_row_and_none_where_there_is_no_row(tmp_path: Path) -> None:
    """3 organs x 3 stages with one missing row draws 8 cells, cardiac organ first.

    A missing row draws nothing: the dataset never looked, and a drawn empty cell
    would read as a measurement of nothing. Brain has no row at the middle stage.
    """
    svg = _bulk_svg(
        _grid_page(
            tmp_path,
            tuple(
                (
                    token,
                    tuple(
                        _placed(organ, value)
                        for organ, value in (("Liver", 5.0), ("Heart", 50.0), ("Brain", 500.0))
                        if not (organ == "Brain" and token == "5 week post conception")
                    ),
                )
                for token in (
                    "4 week post conception",
                    "5 week post conception",
                    "6 week post conception",
                )
            ),
        )
    )
    assert len(_grid_cells(svg)) == 8
    assert svg.count('<rect class="heat-') == 8
    assert _organ_labels(svg) == ["Heart", "Brain", "Liver"]
    assert svg.count('class="chart-tick-label cardiac"') == 1
    # Three stages would give 240-unit columns (`_half_step(3)` x 2); they are
    # capped at 40 px less a 2 px gutter.
    assert set(
        re.findall(r'<rect class="heat-[0-9]" x="[^"]*" y="[^"]*" width="([^"]*)"', svg)
    ) == {"38.0"}

    # A cell's title is escaped, whatever the organ is called.
    hostile = 'A<&"B'
    hostile_svg = _bulk_svg(
        _grid_page(
            tmp_path,
            (("4 week post conception", (_placed("Heart", 50.0), _placed(hostile, 100.0))),),
        )
    )
    assert "<title>A&lt;&amp;&quot;B \u00b7 4w \u00b7 100 rpkm</title>" in hostile_svg
    assert hostile not in hostile_svg


def test_the_grid_legend_names_only_the_states_drawn(tmp_path: Path) -> None:
    """A key for a mark that is not drawn sends the reader hunting (section 4.27).

    Sliced to the legend element: the same phrases appear in the stage tables
    and in cell titles, so a section-wide check would be answered by them.
    """

    def legend(section: str) -> str:
        return re.search(r'<p class="strip-legend heat-key">.*?</p>', section, re.S).group(0)  # type: ignore[union-attr]

    plain = legend(
        _grid_page(
            tmp_path,
            (
                ("4 week post conception", (_placed("Heart", 50.0), _placed("Liver", 5.0))),
                ("5 week post conception", (_placed("Heart", 60.0), _placed("Liver", 6.0))),
            ),
        )
    )
    assert "below the" not in plain
    assert "blank: not sampled" not in plain
    assert "no complete percentile grid" not in plain
    assert "1\u201310" in plain and "1,000+ rpkm" in plain

    full = legend(
        _grid_page(
            tmp_path,
            (
                (
                    "4 week post conception",
                    (_placed("Heart", 50.0), _unplaced("Liver", "below_detection_floor")),
                ),
                (
                    "5 week post conception",
                    (_placed("Heart", 60.0), _unplaced("Liver", "no_quantile_grid")),
                ),
                ("6 week post conception", (_placed("Heart", 70.0),)),
            ),
        )
    )
    assert "below the 1 rpkm floor" in full
    assert "no complete percentile grid is published for this organ at this stage" in full
    assert "blank: not sampled" in full

    # A below-floor cell with no *other* gap: the floor entry is present and
    # the gap entry (and its swatch) is not. Mutating `if grid.gap_reasons:`
    # to `if True:` adds "no reason was recorded for this gap" here.
    floor_only = legend(
        _grid_page(
            tmp_path,
            (
                (
                    "4 week post conception",
                    (_placed("Heart", 50.0), _unplaced("Liver", "below_detection_floor")),
                ),
            ),
        )
    )
    assert "below the 1 rpkm floor" in floor_only
    assert "no reason was recorded for this gap" not in floor_only
    assert "heat-gap" not in floor_only


def test_the_grey_stage_key_appears_only_when_a_stage_label_is_greyed(tmp_path: Path) -> None:
    """A greyed stage label is keyed, and only when one is drawn (section 4.27).

    The label is the only mark on the figure with no key, so the entry is
    present exactly when the heart has no row at some stage. Sliced to the
    legend element, as the test above does, because the phrase is not in the
    stage tables.
    """

    def legend(section: str) -> str:
        return re.search(r'<p class="strip-legend heat-key">.*?</p>', section, re.S).group(0)  # type: ignore[union-attr]

    greyed = legend(
        _grid_page(
            tmp_path,
            (
                ("4 week post conception", (_placed("Heart", 50.0), _placed("Liver", 5.0))),
                ("5 week post conception", (_placed("Liver", 5.0),)),
                ("6 week post conception", (_placed("Heart", 70.0),)),
            ),
        )
    )
    assert "grey stage: no Heart sample" in greyed
    assert "blank: not sampled" in greyed
    assert greyed.index("blank: not sampled") < greyed.index("grey stage: no Heart sample")

    complete = legend(
        _grid_page(
            tmp_path,
            (
                ("4 week post conception", (_placed("Heart", 50.0), _placed("Liver", 5.0))),
                ("5 week post conception", (_placed("Heart", 60.0), _placed("Liver", 6.0))),
            ),
        )
    )
    assert "grey stage" not in complete


def test_an_unregistered_dataset_draws_the_grid_alone(tmp_path: Path) -> None:
    """No `Dataset` means no cardiac tissue, floor or heart chart; the grid still draws."""
    section = _grid_page(
        tmp_path,
        (("4 week post conception", (_placed("Heart", 50.0), _placed("Liver", 5.0))),),
        datasets={},
    )
    assert _grid_cells(_bulk_svg(section)) == ["heat-2", "heat-1"]
    assert "chart-line" not in section
    assert "chart-strip" not in section
    assert 'class="glance"' not in section
    assert (
        '<p class="strip-legend heat-key"><span class="legend-item">'
        '<span class="heat-swatch heat-1"></span> under 10</span>'
    ) in section


def test_the_bulk_block_has_one_fold_and_no_dashed_lines_or_panels(tmp_path: Path) -> None:
    """The dashed other-organ lines, the sparkline panels and their second fold are gone.

    A negative assertion on published bytes, on the *page*: it is the one that
    fails when a deleted sentence comes back (CLAUDE.md section 4.35), and the
    stylesheet is part of the page.
    """
    emitter_page = _expression_page(
        tmp_path,
        {
            GATA4: _organ_panel(
                ("heart", (5.0, 40.0, 275.0)),
                ("liver", (1.0, 2.0, 4.0)),
                phase=_HEART_PHASE,
            )
        },
        _HEART_DATASET,
        phases=_cardiac_phases(),
    )
    section = _expression_section_text(emitter_page)

    assert section.count("How to read this chart") == 1
    for gone in ("spark", "chart-context", "Dashed grey lines", "One panel per organ"):
        assert gone not in emitter_page
    fold = section[section.index('<details class="chart-notes">') :]
    assert "The grid's shades are fixed bins" in fold[: fold.index("</details>")]


# --- The summary line (spec B1) ----------------------------------------------


def _glance_entry(
    *rows: tuple[str, float | None, str | None, str | None],
) -> DatasetProfileEntry:
    """One heart row per stage: `(label, median, not_placed_reason, highest_in)`.

    A `None` median is a stage sampled and not placed, for the reason given. A
    Liver row sits beside every heart row at 1 tpm, so `highest_in` is a fact a
    reader could check against a second organ and not a value that only heart
    reports.
    """
    return _dataset_profile_entry(
        stages=tuple(
            _stage_entry(
                stage=label,
                specificity=_specificity(tissues=("Heart", "Liver"), highest_in=highest_in),
                tissues=(
                    _tissue_entry(
                        tissue="Heart",
                        median=median if median is not None else 0.3,
                        unit="rpkm",
                        placement=_placement() if median is not None else None,
                        not_placed_reason=not_placed_reason,
                    ),
                    _tissue_entry(tissue="Liver", median=1.0, unit="rpkm"),
                ),
            )
            for label, median, not_placed_reason, highest_in in rows
        )
    )


def test_the_summary_line_reads_the_bundle_s_own_counts() -> None:
    """Every field of the sentence, from one fixture where each count differs.

    Placed 3 of 4 sampled stages, the fourth below the floor; peak 50 at a
    prenatal stage; `highest_in` is Heart at two stages, Liver at one and null
    at one. A fixture where any two of those counts agreed would let a mutant
    that published one in place of the other pass, so none of them agree.
    """
    line = _glance_line(
        _glance_entry(
            ("4 week post conception", 10.0, None, "Liver"),
            ("5 week post conception", 20.0, None, "Heart"),
            ("6 week post conception", 50.0, None, "Heart"),
            ("7 week post conception", None, "below_detection_floor", None),
        ),
        "Heart",
        1.0,
    )

    assert line == (
        '<p class="glance"><strong>Heart:</strong> above the 1 rpkm floor at 3 of 4 stages '
        "sampled; peak 50 rpkm at 6 weeks post conception; the highest of the organs sampled "
        "at 2 stages.</p>"
    )


@pytest.mark.parametrize(
    ("rows", "floor", "fragment"),
    [
        pytest.param(
            (
                ("5 week post conception", 50.0, None, "Heart"),
                ("6 week post conception", 50.0, None, "Heart"),
            ),
            1.0,
            "peak 50 rpkm at 5 weeks post conception and 1 other stage;",
            id="tied-peak-names-the-earliest-stage-and-the-other",
        ),
        pytest.param(
            (("5 week post conception", 50.0, None, None),),
            1.0,
            "at no stage.",
            id="zero-wins-reads-no-stage",
        ),
        pytest.param(
            (("5 week post conception", 50.0, None, "Heart"),),
            1.0,
            "at 1 stage.",
            id="one-win-reads-singular",
        ),
        pytest.param(
            (
                ("5 week post conception", 50.0, None, "Heart"),
                ("6 week post conception", None, "no_quantile_grid", None),
            ),
            1.0,
            "placed at 1 of 2 stages sampled; peak 50 rpkm",
            id="a-non-floor-gap-forbids-the-floor-clause",
        ),
        pytest.param(
            (
                ("5 week post conception", 50.0, None, None),
                ("6 week post conception", 40.0, None, None),
                ("7 week post conception", None, "no_quantile_grid", None),
            ),
            1.0,
            "placed at 2 of 3 stages sampled; peak 50 rpkm",
            id="the-placed-clause-counts-placed-and-sampled-separately",
        ),
        pytest.param(
            (
                ("5 week post conception", 50.0, None, None),
                ("6 week post conception", 40.0, None, None),
                ("7 week post conception", None, "below_detection_floor", None),
            ),
            None,
            "placed at 2 of 3 stages sampled; peak 50 rpkm",
            id="an-unknown-floor-clause-counts-placed-and-sampled-separately",
        ),
        pytest.param(
            (
                ("5 week post conception", 50.0, None, None),
                ("6 week post conception", 50.0, None, None),
                ("7 week post conception", 50.0, None, None),
            ),
            1.0,
            "peak 50 rpkm at 5 weeks post conception and 2 other stages;",
            id="three-way-tie-pluralises-the-other-stages",
        ),
        pytest.param(
            (("infant", 50.0, None, None), ("adult", 20.0, None, None)),
            1.0,
            "peak 50 rpkm at infant; the highest",
            id="a-non-prenatal-peak-is-named-by-its-own-label",
        ),
        pytest.param(
            (("<&", 50.0, None, None), ("adult", 20.0, None, None)),
            1.0,
            "peak 50 rpkm at &lt;&amp;; the highest",
            id="a-stage-label-is-escaped",
        ),
        pytest.param(
            (
                ("5 week post conception", 50.0, None, "Heart"),
                ("6 week post conception", None, "below_detection_floor", None),
            ),
            None,
            "placed at 1 of 2 stages sampled; peak 50 rpkm",
            id="an-unknown-floor-forbids-the-floor-clause",
        ),
    ],
)
def test_the_summary_line_s_branches(
    rows: tuple[tuple[str, float | None, str | None, str | None], ...],
    floor: float | None,
    fragment: str,
) -> None:
    """Each wording the sentence can take, against a fixture that forces it.

    Each row differs from the others in one respect, so each assertion is the
    only thing a mutant could break: a tie that is counted as a win, a zero
    that reads "1 stage", a grid gap that is called a low reading, and a floor
    that is announced without being known.
    """
    assert fragment in _glance_line(_glance_entry(*rows), "Heart", floor)


_ORGAN_ORDER_SCRIPT: Final = """
from chd_atlas.build.pages import _organ_order
names = ["Spleen", "Kidney", "Brain", "Liver", "Lung", "Colon", "Testis", "Ovary"]
entry = {"stages": [{"tissues": [{"tissue": n} for n in ["Heart", *names]]}]}
print("|".join(_organ_order(entry, frozenset({"Heart"}))))
"""


@pytest.mark.parametrize("seed", [0, 1, 7, 12345, 99991])
def test_the_organ_order_sorts_the_non_cardiac_organs_under_every_hash_seed(seed: int) -> None:
    """The grid's row order is published, so the non-cardiac group must be sorted.

    Run in a subprocess per explicit `PYTHONHASHSEED`, because one interpreter
    has one seed (CLAUDE.md section 4.13/4.42). Measured 2026-11-02 against the
    mutant that iterates `list({t for t in sampled})` instead of `sampled`:
    for these eight organ names the sorted order came out under **0 of 300**
    seeds (0..299, one probe interpreter per seed), so each seed here kills it
    with probability near 1. Were the true sorted-by-luck rate as high as 1%
    (about the 95% upper bound for 0/300), five seeds would leave a survival
    probability of 0.01^5 = 1e-10.
    """
    result = subprocess.run(
        [sys.executable, "-c", _ORGAN_ORDER_SCRIPT],
        env={**os.environ, "PYTHONHASHSEED": str(seed)},
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == "Heart|Brain|Colon|Kidney|Liver|Lung|Ovary|Spleen|Testis"


def test_no_summary_line_without_a_placed_heart_stage(tmp_path: Path) -> None:
    """Nothing placed means no summary line, and the lede still says why.

    Pins the wiring, not the wording: a dataset whose heart is below the floor
    at every stage must not gain a summary, and must keep the sentence that
    explains why there is no chart at all (D42).
    """
    section = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", None), ("5wpc", None))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )

    assert 'class="glance"' not in section
    assert "so there is no trajectory to draw" in section


def test_the_summary_line_reaches_the_page_before_the_chart(tmp_path: Path) -> None:
    """The line is the first thing under the dataset heading, as spec B1 places it.

    The literal is the one the unit test pins, so a page that drifted from the
    function (a line emitted after the figure, or a second copy) fails here.
    """
    section = _expression_section_text(
        _expression_page(
            tmp_path,
            {
                GATA4: _organ_panel(
                    ("Heart", (5.0, 40.0, 275.0)),
                    stages=(
                        "3 week post conception",
                        "4 week post conception",
                        "6 week post conception",
                    ),
                    highest_in="Heart",
                )
            },
            _CAPITAL_HEART,
            phases=_cardiac_phases(),
        )
    )
    heading_end = section.index("</h3>") + len("</h3>")
    glance_start = section.index('<p class="glance">')
    glance_end = section.index("</p>", glance_start) + len("</p>")

    assert section.count('<p class="glance">') == 1
    assert glance_start == heading_end
    assert section[glance_start:glance_end] == (
        '<p class="glance"><strong>Heart:</strong> above the 1 tpm floor at 3 of 3 stages '
        "sampled; peak 275 tpm at 6 weeks post conception; the highest of the organs sampled "
        "at 3 stages.</p>"
    )


# --- Organ fixtures (the small multiples were deleted 2026-11-02) -------------
#
# The fixtures below all sample at least one cardiac and one non-cardiac organ,
# at genuinely different magnitudes. CLAUDE.md section 4.36 is why that is
# spelled out: a fixture whose organs all share the value under test measures
# nothing, and it has cost this repository four separate defects. A panel where
# every organ rendered `chart-cardiac` and one where none did would be
# indistinguishable on a heart-only fixture.


def _organ_panel(
    *organs: tuple[str, tuple[float | None | EllipsisType, ...]],
    stages: tuple[str, ...] = ("4wpc", "5wpc", "6wpc"),
    highest_in: str | None = "heart",
    phase: PhaseInfo | None = None,
    unplaced_reason: str = "below_detection_floor",
) -> ExpressionProfile:
    """Several organs measured across the same stages, each at its own magnitude.

    The same three states `_heart_series` carries, for the same reason. A
    `float` is a median this dataset placed against its percentile grid; a
    `None` is a stage this organ was sampled at and not placed, below the
    detection floor; `...` is a stage this dataset has **no row for that
    organ** at, which is a third fact and breaks a panel's line with nothing
    to show for it -- E-MTAB-6814's real shape, identically for every gene.

    Every stage carries a real `specificity`, unlike `_heart_series`: tau is
    what this panel exists to make visible, and `_small_multiples` draws
    nothing for a dataset that reports none.

    `phase` defaults to `_phase_info()`, whose id is not in `_cardiac_phases()`
    -- so a panel built here draws no trajectory unless a caller asks for a
    phase the vocabulary actually names. That default is what every test
    written before the trajectory landed relies on; pass one to get bands.

    `unplaced_reason` is a parameter for the reason it is one on
    `_heart_series`: a `None` here is a stage this dataset sampled and did
    not place, and *why* decides what the caption beside the panel is
    entitled to say about it. Below a floor is a low reading about the gene;
    a missing percentile grid is a hole in the reference and says nothing
    about the gene at all.
    """
    names = tuple(name for name, _ in organs)
    entries = tuple(
        _stage_entry(
            stage=token,
            phase=phase,
            specificity=_specificity(tissues=names, highest_in=highest_in),
            tissues=tuple(
                _tissue_entry(
                    tissue=name,
                    median=values[index] if isinstance(values[index], float) else 0.3,
                    unit="tpm",
                    placement=_placement() if isinstance(values[index], float) else None,
                    not_placed_reason=(
                        None if isinstance(values[index], float) else unplaced_reason
                    ),
                )
                for name, values in organs
                if values[index] is not ...
            ),
        )
        for index, token in enumerate(stages)
    )
    return _expression_profile((_dataset_profile_entry(stages=entries),))


def test_one_placed_stage_with_no_band_says_why_in_words_and_draws_no_chart(
    tmp_path: Path,
) -> None:
    """One organ sampled means no specificity to show, and no phase means no chart.

    D42: where there is nothing to plot the page says so; it never renders an
    empty frame. (This was `test_small_multiples_are_absent_when_tau_is_
    undefined` until the small multiples were deleted on 2026-11-02; its
    figure assertions went with them and the lede assertions below are the
    ones that survive.)

    **The sentence the dataset block opens with is asserted too**, because
    this fixture is also the smallest case of D42's other refusal and the two
    share a page. One heart stage, placed, and no curated phase vocabulary:
    no trajectory can be drawn, and the reason is the missing band, not a
    missing measurement. Until 2026-08-20 the page said "no measurement there
    is placed against this dataset's percentile grid" of the one stage whose
    percentile it printed in the table below, and said it of "1 stages".
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        stage="4wpc",
                        specificity=None,
                        specificity_unavailable_reason="one_organ_sampled",
                        tissues=(
                            _tissue_entry(tissue="heart", unit="tpm", placement=_placement()),
                        ),
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET))

    assert "chart-line" not in section
    assert "chart-strip" not in section
    assert "only one organ was sampled" in section

    lede = _dataset_lede(section)
    assert "no measurement there is placed" not in lede
    assert "detection floor at every" not in lede
    assert "1 stages sampled" not in lede, "a count with no singular form"
    assert "1 stage sampled" in lede


# --- The order every cardiac organ is named in (CLAUDE.md section 4.36) ------
#
# `_dataset_block` builds `cardiac` as a `frozenset`, and three places walk it:
# `_charted_tissue` (which tissue the heart chart draws), `_organ_order` (the
# grid's rows) and `_sampled_cardiac`'s sentence. All sort it, and until
# 2026-08-20 nothing could tell (the third place was a spark caption's "this
# atlas declares" clause, deleted with the small multiples on 2026-11-02).
# `Dataset.cardiac_tissues` has exactly one member on the committed corpus and
# every fixture above declares one, so a one-element set iterates in one order
# under every seed. Two measurements, taken on the then three `sorted(cardiac)`
# calls: replacing all of them with a bare `cardiac` leaves the build
# byte-identical under PYTHONHASHSEED 0 and 12345, and mutating them one at a
# time was killed by 0 of the 1,050 tests that existed at 38fb290, the commit
# before this section.
#
# **Seven organs, and the number is measured rather than chosen.** A `frozenset`
# of *two* names iterates in sorted order about half the time, so a two-organ
# fixture would catch a dropped sort on roughly half of CI's random seeds --
# which is not a guard (section 4.12). Measured 2026-08-20 by constructing
# `frozenset(names)` under 400 explicit `PYTHONHASHSEED` values and counting how
# often it iterates pre-sorted:
#
#     2 names (atrium, ventricle)                 197/400   49.2%
#     3 names (aorta, heart, ventricle)            53/400   13.2%
#     4 names (+ atrium)                           17/400    4.2%
#     5 names (+ myocardium, endocardium)           6/400    1.5%
#     6 names (+ septum)                            2/400    0.5%
#     7 names (+ epicardium)   <- these seven       0/400    0.0%
#
# Widening the fixture is the cheap half of section 4.42's remedy and the whole
# of it here, because the rate reaches zero at a size a real dataset could
# plausibly declare. Re-measure rather than cite this if the names change: the
# rate is a property of the exact strings, not of the count.
_CARDIAC_ORGANS: tuple[str, ...] = (
    "ventricle",
    "aorta",
    "septum",
    "myocardium",
    "endocardium",
    "atrium",
    "epicardium",
)

# Declared above in a deliberately unsorted order, so a renderer walking
# `Dataset.cardiac_tissues` itself rather than the sorted frozenset fails too.
_CARDIAC_ORGANS_SORTED: tuple[str, ...] = (
    "aorta",
    "atrium",
    "endocardium",
    "epicardium",
    "myocardium",
    "septum",
    "ventricle",
)

_MANY_CARDIAC_DATASET = {"E-MTAB-6814": _profile_dataset(cardiac_tissues=_CARDIAC_ORGANS)}


_BAND_LABELS: Final = {"alpha": "Zulu", "bravo": "Yankee", "charlie": "Xray", "delta": "Whisky"}

# Deliberately neither sorted nor reverse-sorted: a walk in this order gives a
# fourth sequence again, so "sorted by id" is the only rule that produces the
# literal asserted below.
_BAND_WALK: Final = ("delta", "alpha", "charlie", "bravo")


def _phase_walk(*phase_ids: str) -> ExpressionProfile:
    """One placed heart stage per phase, in the order given."""
    return _expression_profile(
        (
            _dataset_profile_entry(
                stages=tuple(
                    _stage_entry(
                        stage=f"{index + 4}wpc",
                        phase=_phase_info(phase_ids=(phase_id,)),
                        tissues=(
                            _tissue_entry(
                                tissue="heart",
                                median=5.0 * (index + 1),
                                unit="tpm",
                                placement=_placement(),
                            ),
                        ),
                    )
                    for index, phase_id in enumerate(phase_ids)
                )
            ),
        )
    )


def test_the_phase_strip_names_each_banded_phase_once_in_chronological_order(
    tmp_path: Path,
) -> None:
    """The strip is a timeline: first stage, then last stage, then id.

    `_banded_phases` sorts by id and the strip must not inherit that. The four
    phases here have ids, labels, walk order and reversed walk order that are
    four different sequences, so only the rule "first stage, last stage, id"
    produces the literals below, and the second walk shows the order follows
    the *stages* rather than a rule that happens to be stable under reversal.
    The later phase of the pair sharing a stage has the alphabetically earlier
    id, so the id is shown to break a tie and nothing more.

    Each label is asserted as strip *text* exactly once: the bar's `<title>`
    carries it too, but a second text node (the faint band used to repeat it)
    would name a phase twice.
    """
    vocabulary = _cardiac_phases(
        *(
            _cardiac_phase(phase_id=phase_id, label=label)
            for phase_id, label in _BAND_LABELS.items()
        )
    )

    def strip(*phase_ids: str) -> tuple[list[str], str]:
        section = _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _phase_walk(*phase_ids)},
                _HEART_DATASET,
                phases=vocabulary,
            )
        )
        svg = _bulk_svg(section)
        return _strip_labels(svg), svg

    labels, svg = strip(*_BAND_WALK)
    # delta, alpha, charlie, bravo = Whisky, Zulu, Xray, Yankee. Id order would
    # be Zulu, Yankee, Xray, Whisky; label order Whisky, Xray, Yankee, Zulu.
    assert labels == ["Whisky", "Zulu", "Xray", "Yankee"]
    for label in labels:
        assert svg.count(f">{label}</text>") == 1
    assert strip(*reversed(_BAND_WALK))[0] == ["Yankee", "Xray", "Zulu", "Whisky"]

    # Two phases over one stage range: the id breaks the tie. alpha ("Zulu")
    # sorts before bravo ("Yankee") by id, and by nothing else.
    tied = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        stage="4wpc",
                        phase=_phase_info(phase_ids=("bravo", "alpha")),
                        tissues=(
                            _tissue_entry(tissue="heart", unit="tpm", placement=_placement()),
                        ),
                    ),
                )
            ),
        )
    )
    tied_section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: tied}, _HEART_DATASET, phases=vocabulary)
    )
    assert _strip_labels(_bulk_svg(tied_section)) == ["Zulu", "Yankee"]


def test_a_strip_label_that_would_leave_the_figure_goes_left_of_its_bar() -> None:
    """The label sits right of its bar, unless that would run off the figure.

    Called directly: no real phase reaches the right edge today (the vocabulary
    ends in the prenatal window, and E-MTAB-6814's post-natal stages are the
    last ten of 21), so no built page exercises it. A guard on a branch no page
    takes is still a guard (CLAUDE.md section 9), and an untested one is not.
    """
    phase = _cardiac_phase(label="A long phase label that cannot fit")
    early, _ = _phase_strip([(phase, 0, 1)], 21, 4.0)
    late, _ = _phase_strip([(phase, 19, 20)], 21, 4.0)
    assert 'text-anchor="start"' in early
    assert 'text-anchor="end"' in late


def test_the_cardiac_organs_are_charted_and_listed_in_one_order_whatever_the_frozenset_does(
    tmp_path: Path,
) -> None:
    """Seven cardiac organs: which one is charted, which are named, which rows lead the grid.

    Against literals rather than `sorted(_CARDIAC_ORGANS)`, for the reason
    `test_build_landing.py` pins published wording against literals (section
    4.38): a test that computes its own expectation from the same expression
    the code uses compares the module to itself and passes for any ordering
    both agree on.

    Three orders here come from the same `frozenset`: the tissue the chart
    draws (the first, sorted, with a placed median), the `Only ... is charted`
    clause's list of the rest, and the cardiac rows at the top of the grid.
    """
    profile = _organ_panel(
        *((organ, (5.0, 40.0, 275.0)) for organ in _CARDIAC_ORGANS),
        ("liver", (1.0, 2.0, 4.0)),
        highest_in="ventricle",
        phase=_phase_info(phase_ids=("atrial_septum_morphogenesis",)),
    )
    section = _expression_section_text(
        _expression_page(
            tmp_path, {GATA4: profile}, _MANY_CARDIAC_DATASET, phases=_cardiac_phases()
        )
    )
    svg = _bulk_svg(section)

    assert re.findall(r"<title>Median abundance in ([^,]+),", svg) == ["aorta"]
    assert (
        "Only aorta is charted; atrium, endocardium, epicardium, myocardium, septum, "
        "ventricle appear in the grid and the tables."
    ) in _band_caption_text(section)
    assert _organ_labels(svg) == [
        "aorta",
        "atrium",
        "endocardium",
        "epicardium",
        "myocardium",
        "septum",
        "ventricle",
        "liver",
    ]
    # The cardiac rows are the bold ones, and only they.
    assert svg.count('class="chart-tick-label cardiac"') == 7


def test_both_sentences_that_list_every_cardiac_organ_list_them_in_one_order(
    tmp_path: Path,
) -> None:
    """The two orders a page states in words rather than in geometry.

    `_no_trajectory_sentence` renders exactly when no cardiac organ has a
    heart chart. Every cardiac organ here is sampled and below the floor; liver
    is placed, so the grid still has a placed cell to draw (the second list,
    the grid's own rows, comes from `_sampled_tissues`' sort over a `set`).

    Both are asserted as one literal apiece rather than as a set of names,
    because the defect is the *order* and a membership check cannot see it.
    """
    profile = _organ_panel(
        *((organ, (None, None, None)) for organ in _CARDIAC_ORGANS),
        ("liver", (0.5, 1.2, 2.0)),
        highest_in="liver",
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _MANY_CARDIAC_DATASET)
    )

    assert (
        "aorta (3 stages sampled), atrium (3 stages sampled), "
        "endocardium (3 stages sampled), epicardium (3 stages sampled), "
        "myocardium (3 stages sampled), septum (3 stages sampled), "
        "ventricle (3 stages sampled)"
    ) in section
    # The second list, from `_sampled_tissues`' own sort over a `set` of the
    # organs the payload names -- a different frozen collection with the same
    # failure mode, and the same eight-name fixture measures 0/400 pre-sorted.
    # Cardiac organs lead the grid, sorted; liver follows.
    assert _organ_labels(_bulk_svg(section)) == [
        "aorta",
        "atrium",
        "endocardium",
        "epicardium",
        "myocardium",
        "septum",
        "ventricle",
        "liver",
    ]


def test_a_chart_label_refuses_an_anchor_outside_svgs_own_vocabulary() -> None:
    """The one attribute *value* the chart code builds from a parameter.

    `render.data_table` refuses an attribute *name* rather than escaping it,
    because `html.escape` rewrites the quotes and leaves the space and the `=`
    alone -- so a value reaching an attribute position can open a second
    attribute from a payload carrying no angle bracket. This is that rule one
    level down.

    Every call site passes a literal, so this is a guard on a bypassed gate
    (CLAUDE.md section 9): reaching it means a caller derived an attribute
    value from data, and that must fail rather than publish.
    """
    hostile = '"><script>alert(1)</script>'
    with pytest.raises(ValueError, match="text-anchor"):
        _forest_text(0.0, 0.0, Markup("x"), anchor=hostile)

    # and the vocabulary it does accept still works
    for anchor in ("start", "middle", "end"):
        assert f'text-anchor="{anchor}"' in _forest_text(0.0, 0.0, Markup("x"), anchor=anchor)


def _combined(**overrides: object) -> CombinedRow:
    """KMT2D's loss-of-function JL-FDR row from DS6, at the precision the mirror holds."""
    payload: dict[str, object] = {
        "study": "PMID:40127276",
        "gene": GATA4,
        "cohort_stratum": "all",
        "lesion_group": None,
        "consequence_class": "lof",
        "combines": ("mutation_model", "control_cohort"),
        "method": "joint_local_fdr",
        "pvalue": 1.6582359932982601e-38,
        "local_fdr": 2.6973631548392e-34,
        "method_note": "The study's own meta-analysis of its de novo and case-control arms.",
        "source": "sierant2025_ds6",
    }
    payload.update(overrides)
    return CombinedRow(**payload)  # type: ignore[arg-type]


def _study_block(section: str, pmid_number: str) -> str:
    """One study's block within the burden section: from its `<h3>` to the next.

    Split on the tag's opening, not on a bare `<h3>`: since the card links to
    each study (D55) the heading carries an `id`.
    """
    blocks = re.split(r"<h3[ >]", section)[1:]
    matches = [block for block in blocks if pmid_number in block.split("</h3>", 1)[0]]
    assert len(matches) == 1, f"expected one block for {pmid_number}"
    return matches[0]


def test_a_study_s_own_combination_is_stated_in_its_block_and_labelled_for_what_it_is(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The line, its label, and its caveat, pinned against literals (§4.38).

    Three things a reader would get wrong without the wording:

    * the number is the *authors'* combination within one cohort, not a pooling
      this atlas did -- D33 forbids the latter and the page must not look like it;
    * the second statistic is a **joint local FDR**, not a corrected p: it is
      smaller than its own p for about half the panel, so a reader who took it
      for one would conclude impossible things;
    * it combines two named arms, which the rows below it show.

    Scoped to the study's block (§4.19): the other study on the page carries no
    combination and must render none.
    """
    sierant = [
        _burden_row(
            study="PMID:40127276",
            cohort_stratum="all",
            comparator="mutation_model",
            origin="de_novo",
            count_unit="de_novo_mutations",
            n_case_carriers=21,
            n_cases=3887,
            n_control_carriers=None,
            n_controls=None,
            control_cohorts=(),
            expected_count=0.17,
            effect=126.8,
            effect_measure="enrichment_ratio",
            effect_bound=None,
            ci_low=None,
            pvalue=6.65e-37,
            pvalue_test="poisson",
            case_cohorts=("pcgc",),
        ),
    ]
    rows = [*sierant, _burden_row()]
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        burden={GATA4: rows},
        cohorts=_COHORTS,
        families=cohort_families(rows, _COHORTS),
        axes=evidence_axes(rows),
        combined={GATA4: [_combined()]},
        cell_types={},
    )
    section = _burden_section_text(_page(tmp_path, "HGNC_4173.html"))
    block = _study_block(section, "40127276")

    assert (
        "<strong>Combined across this study&#x27;s arms</strong> (de novo, trios + "
        "case-control, by the joint local FDR method): loss-of-function p 1.66e-38, "
        "joint local FDR 2.7e-34."
    ) in block
    assert "the authors&#x27; own meta-analysis of the rows below, within one cohort" in block
    assert "not a corrected p-value" in block
    # Never a bare "FDR" beside the number: every occurrence is "joint local FDR".
    assert block.count("FDR 2.7e-34") == block.count("joint local FDR 2.7e-34") == 1
    # The row's method note reaches the page with the arms' footnotes.
    assert "own meta-analysis of its de novo and case-control arms." in block
    assert "Combined across" not in _study_block(section, "42230622")


_PULMONIC = Annotation("HP:0001642", "Abnormal pulmonary valve morphology", "Frequent (79-30%)")


def _noonan_validity() -> GeneValidity:
    """One ClinGen-path disease and two Orphanet-scoped ones, one of them a subtype.

    Three because each kills a different mutant (CLAUDE.md section 4.36): a
    loop over every record's disease fails on Holt-Oram, which has no warrant; a
    note that ignored `via` reads the same for both Noonan terms; and an
    unsorted loop swaps them. A note for a warrant the gene holds no record for
    is guarded by the test's own extra warrant, not by these records: until
    2026-10-26 every warrant passed was one of them, and a loop over every
    warrant survived all 1,297 tests (built, it put all 36 notes on all 119
    pages).
    """
    base = _validity().records[0]
    return GeneValidity(
        records=(
            base,
            replace(
                base,
                source=ValiditySource.GENCC,
                disease="MONDO:0018997",
                disease_label="Noonan syndrome",
                gcep=None,
                sop=None,
                submitter="G2P",
            ),
            replace(base, disease="MONDO:0012547", disease_label="Noonan syndrome 4"),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )


def test_an_orphanet_scoped_syndrome_says_why_it_counts_as_chd_above_both_groups(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The note PTPN11's page needs, pinned against literals and sliced (§4.19, §4.38).

    A Noonan grade is a grade for the syndrome. Without this note a reader meets
    it in a CHD atlas with nothing saying why the syndrome counts, and takes it
    for a grade for heart disease. The scope table sits above both groups of
    classifications, because a warrant concerns the disease, which either group
    may grade. The scope rule at the top of the page named the Orphanet path
    too, and was checked here on the built element; it moved to
    `methods.html#scope` on 2026-10-26 (D50), and its Orphanet clause and
    `scope.json` parenthetical are pinned there, by `test_build_methods.py`.
    The page keeps a short form that names neither.

    **The whole-table lead-in names no one disease (2026-10-26).** The single
    sentences over the old paragraphs read "The classifications below are for the
    syndrome, not for heart disease", and above a table that was not only that
    syndrome's -- NOTCH1's lists an admitting record for a term a ClinGen expert
    panel scoped, and an Orphanet record for familial bicuspid aortic valve,
    beneath its two Adams-Oliver notes -- a single undirected sentence was false
    on that page. This fixture carries the same shape: a Holt-Oram record, scoped
    by an authority, below two Noonan rows. The lead-in says *each* disease's
    classifications are for *that* disease, which is true however many diseases
    the section below names; it names no parent, so it cannot be misread as one
    row's.

    **It ends on "that disease", not "the syndrome"** (review, 2026-10-26): in a
    subtype's row the "why" cell names the parent -- "Subtype of Noonan syndrome
    ..." -- so "the syndrome as a whole" would have the parent as its nearest
    antecedent. Asserted absent.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _noonan_validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        scope_warrants={
            "MONDO:0018997": Warrant("ORPHA:648", "Noonan syndrome", "direct", None, (_PULMONIC,)),
            "MONDO:0012547": Warrant(
                "ORPHA:648", "Noonan syndrome", "subtype", "OMIM:610733", (_PULMONIC,)
            ),
            # In scope, and not a disease this gene has a record for.
            _NSML: Warrant(
                "ORPHA:500", "Noonan syndrome with multiple lentigines", "direct", None, ()
            ),
        },
        scope_labels={
            "MONDO:0018997": "Noonan syndrome",
            "MONDO:0012547": "Noonan syndrome 4",
            _NSML: "Noonan syndrome with multiple lentigines",
        },
        cell_types={},
    )
    page = _page(tmp_path, "HGNC_4173.html")
    section = _slice_between(page, "<h2>Mirrored gene&ndash;disease validity</h2>", ("<h2",))

    direct_row = (
        "<td>Noonan syndrome</td>"
        "<td>Orphanet records a structural heart malformation in at least 30% of this "
        "disease's cases</td>"
        "<td>Abnormal pulmonary valve morphology at <em>Frequent (79-30%)</em> "
        "(Orphanet, ORPHA:648)</td>"
    )
    subtype_row = (
        "<td>Noonan syndrome 4</td>"
        "<td>Subtype of Noonan syndrome; Orphanet records a structural heart "
        "malformation in at least 30% of its cases (Orphanet maps OMIM:610733 under "
        "ORPHA:648)</td>"
        "<td>Abnormal pulmonary valve morphology at <em>Frequent (79-30%)</em></td>"
    )
    assert "Why these diseases are in scope" in section
    assert "Each disease&#x27;s classifications in the tables below are for that " in section
    assert direct_row in section and subtype_row in section
    # The lead-in names no one disease and no parent, because the validity table
    # below is not only these syndromes': Holt-Oram's record sits under both rows
    # here, as NOTCH1-related AOS spectrum disorder's does under NOTCH1's.
    assert "The classifications below are for the syndrome" not in section
    # And it ends on "that disease": after a subtype's "Subtype of Noonan
    # syndrome", "the syndrome as a whole" reads as the parent's.
    assert "for the syndrome as a whole" not in section
    # Sorted by disease id: MONDO:0012547 before MONDO:0018997.
    assert section.index(subtype_row) < section.index(direct_row)
    # The scope table sits above both classification groups, not under either.
    scope_table = section[section.index("<table") : section.index('id="validity-table"')]
    assert "Why these diseases are in scope" in section[: section.index("<table")]
    # Exactly the two Noonan rows -- a disease an expert panel curated as CHD
    # (Holt-Oram) and a scope term the gene has no record for (NSML) get no row.
    assert scope_table.count("<tr><td>") == 2
    assert "Holt-Oram" not in scope_table
    assert "multiple lentigines" not in page


# `test_the_grade_key_no_longer_denies_an_in_scope_dispute_it_now_displays` sat
# here until 2026-10-26, pinning the grade key's dispute sentence on a gene page.
# The key moved to `methods.html#scope` (D50) and the test moved with it, its
# docstring and its literal reworded for the new position:
# `test_scope_states_where_a_dispute_is_listed_without_pointing_below` in
# `test_build_methods.py`. Its gene-page negative ("No gene here was disputed or
# refuted" absent) is asserted on both page kinds by
# `test_both_page_kinds_leave_the_grade_key_to_methods_html`.


@pytest.mark.parametrize(
    ("moi", "shown"),
    [
        ("AD", "Autosomal dominant"),
        ("AR", "Autosomal recessive"),
        ("XL", "X-linked"),
        ("SD", "Semidominant"),
        ("MT", "Mitochondrial"),
        ("UD", "Unknown"),
        ("Autosomal dominant", "Autosomal dominant"),
        ("ZZ", "ZZ"),
    ],
)
def test_inheritance_reads_in_one_vocabulary_and_an_unknown_code_as_itself(
    moi: str, shown: str
) -> None:
    """D51: ClinGen's `AD` and GenCC's `Autosomal dominant` are one word on the page.

    Every ClinGen code in the committed mirror is a case here (measured
    2026-10-25: exactly these six). The last case is the one a mapping table
    gets wrong silently: a code nothing taught it must render as itself, never
    blank, because blank reads as "not stated".
    """
    assert _inheritance(moi) == shown


@pytest.mark.parametrize(
    ("timestamp", "shown"),
    [
        ("2024-11-07T17:00:00.000Z", "2024-11-07"),
        ("2018-05-31 20:04:32", "2018-05-31"),
        ("2025-01-17T00:00:00.000000Z", "2025-01-17"),
        ("2019-01-01T00:00:00", "2019-01-01"),
        (None, "—"),
        ("", "—"),
        ("before 2010", "before 2010"),
    ],
)
def test_a_date_shows_its_calendar_day_and_anything_else_verbatim(
    timestamp: str | None, shown: str
) -> None:
    """The timestamp shapes the two mirrors publish, and what is not one.

    The first four are shapes measured in the committed mirrors on 2026-10-25:
    ClinGen's `...T...000Z`, GenCC's space-separated form, GenCC's microsecond
    `Z` form, and ClinGen's zone-less `T` form (65 rows). GenCC's millisecond
    `Z` variant (10 rows) is not a separate case: it shares the first one's
    shape. The rest are a missing value and one that is not a date.
    """
    assert _date(timestamp) == shown


def _panel_and_copy(panel_date: str | None, **copy_overrides: object) -> GeneValidity:
    """A ClinGen record and ClinGen's own GenCC copy of it, the copy adjustable.

    The copy's timestamp is a different *shape* from the panel's on purpose, as
    in the mirrors: equality has to be on the calendar day, not the string.
    """
    panel = replace(_validity().records[0], classification_date=panel_date)
    copy = ValidityRecord(
        source=ValiditySource.GENCC,
        classification=panel.classification,
        classification_term=panel.classification_term,
        disease=panel.disease,
        disease_label=panel.disease_label,
        moi="Autosomal dominant",
        report_url=None,
        submitter="ClinGen",
        submitted_on="2024-11-07 17:00:00",
        sgc_id="SGC-000002",
    )
    return GeneValidity(
        records=(panel, replace(copy, **copy_overrides)),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )


@pytest.mark.parametrize(
    ("panel_date", "copy_overrides", "listed"),
    [
        ("2024-11-07T17:00:00.000Z", {}, False),
        (
            "2024-11-07T17:00:00.000Z",
            {"classification": Classification.STRONG, "classification_term": "Strong"},
            True,
        ),
        ("2024-11-07T17:00:00.000Z", {"submitted_on": "2019-02-01 00:00:00"}, True),
        ("2024-11-07T17:00:00.000Z", {"moi": "Autosomal recessive"}, True),
        ("2024-11-07T17:00:00.000Z", {"disease": "MONDO:0007733"}, True),
        ("2024-11-07T17:00:00.000Z", {"submitter": "Ambry Genetics"}, True),
        (None, {"submitted_on": None}, True),
    ],
    ids=[
        "exact copy",
        "grade differs",
        "date differs",
        "inheritance differs",
        "disease differs",
        "another submitter",
        "both undated",
    ],
)
def test_clingens_own_gencc_copy_is_listed_only_where_it_differs_from_the_panel(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    panel_date: str | None,
    copy_overrides: dict[str, object],
    listed: bool,
) -> None:
    """D51. ClinGen submits every classification to GenCC under its own name.

    Measured 2026-10-25: all 121 such copies among the 119 published genes match
    their direct record exactly, so every one was a duplicate row. Each case
    here changes one thing the match is on and must list the copy again -- a
    copy that disagrees is a different statement, and hiding it would make the
    page agree with itself by deleting the disagreement. "another submitter" is
    corroboration, never a copy. "both undated" is what stops two missing dates
    counting as agreement. The notice about ClinGen's copies must appear whenever
    the copy is ClinGen-submitted -- folded away or listed -- so a visible ClinGen
    row is never unexplained, and be absent only for "another submitter", which is
    no copy at all.

    **The notice's last sentence says "every copy" since 2026-10-26**, pinned
    whole. It said "This gene's JSON keeps every record", beneath an intro that
    now says the table may not be every classification the gene has -- so it
    read as though the JSON held the out-of-scope records the intro means. It
    holds none; it does hold every copy (`methods.html#validity` says the same).
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _panel_and_copy(panel_date, **copy_overrides)},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    others = _validity_table(page)
    assert others.count('data-authority="gencc"') == (1 if listed else 0)
    notice = (
        "<p>ClinGen also submits its own classifications to the Gene Curation Coalition; "
        "such a copy is listed here only where it differs from the panel record above, "
        "carries no date, or has no panel record to match. "
        "This gene's JSON keeps every copy.</p>"
    )
    assert (notice in others) is (copy_overrides.get("submitter", "ClinGen") == "ClinGen")
    assert "keeps every record" not in page


def test_a_copy_spelling_the_panels_rung_differently_is_still_folded(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D51 compares the rung, not the word -- and the corpus is why.

    ClinGen's panel records `Disputed` where its GenCC copy says `Disputed
    Evidence`: 9 real copies, all on the six published genes that carry
    `has_conflicting_evidence` (BRAF, KRAS, MAP2K1, PTPN11, RAF1, SOS1).
    Comparing `classification_term` would list every one of them as a second,
    apparently different ClinGen verdict -- measured on a mutant build
    2026-10-25 -- while every other fold case passed, because each reuses the
    panel's own word.
    """
    pair = _panel_and_copy(
        "2024-11-07T17:00:00.000Z",
        classification=Classification.DISPUTED,
        classification_term="Disputed Evidence",
    )
    panel = replace(
        pair.records[0], classification=Classification.DISPUTED, classification_term="Disputed"
    )
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: replace(pair, records=(panel, pair.records[1]))},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    others = _validity_table(_page(tmp_path, "HGNC_4173.html"))
    assert others.count('data-authority="gencc"') == 0


def test_the_validity_section_groups_by_authority_and_reads_one_vocabulary(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D51: panel first, every row dated, one word for one inheritance, SOP kept.

    The fixture's two records disagree on every surface the old table exposed
    raw -- ClinGen's `AD` against GenCC's phrase, an ISO timestamp against a
    space-separated one -- so each raw surface the old table exposed is pinned
    by at least one assertion.

    **The intro states the table's own truth condition since 2026-10-26 (D50).**
    The grade key said "a gene page lists only records for diseases in this
    atlas's remit", and the key left the page for `methods.html#scope`; so the
    intro now says it where the table is, and links `methods.html#validity`
    for how the records are listed. Pinned whole, sliced from the heading to
    the first group: nothing else sits there on a gene with no Orphanet note.
    It says the table **may not** be every classification the gene has, not
    that it is not: measured 2026-10-26 against the committed mirrors, 43 of
    the 119 published genes have no ClinGen or GenCC record for any disease
    outside scope, so "is not" was false on 43 pages.
    """
    panel = replace(_validity().records[0], classification_date="2024-11-07T17:00:00.000Z")
    submission = ValidityRecord(
        source=ValiditySource.GENCC,
        classification=Classification.STRONG,
        classification_term="Strong",
        disease=panel.disease,
        disease_label=panel.disease_label,
        moi="Autosomal dominant",
        report_url=None,
        submitter="Ambry Genetics",
        submitted_on="2018-05-31 20:04:32",
        sgc_id="SGC-000003",
    )
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={
            GATA4: GeneValidity(
                records=(panel, submission),
                state=ValidityState.EXPERT_CURATED,
                has_source_discordance=False,
            )
        },
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    section = _validity_table(_page(tmp_path, "HGNC_4173.html"))
    headers = re.findall(r'<th scope="col">([^<]*)</th>', section)
    assert headers == ["curated by", "disease", "inheritance", "classification", "date", "report"]
    assert section.count("<table") == 1
    # One table, the expert panel's row first and the only one marked gold standard.
    assert section.index("Syndromic Disorders Gene Curation Expert Panel") < section.index(
        "Ambry Genetics"
    )
    assert section.count("Gold standard</span>") == 1
    assert section.index("Gold standard") < section.index("Ambry Genetics")
    assert section.count('data-authority="clingen"') == 1
    assert section.count('data-authority="gencc"') == 1
    assert "<td>2024-11-07</td>" in section and "T17:00" not in section
    assert "<td>2018-05-31</td>" in section and "20:04:32" not in section
    intro = (
        "<p>Every classification below is an upstream panel's or submitter's, "
        "republished with its provenance intact. Only records for diseases in this "
        "atlas&#x27;s remit are listed, so this may not be every classification this gene "
        'has (<a href="../methods.html#validity">how these records are listed</a>).</p>'
    )
    head = section[: section.index("<table")]
    assert head == "<h2>Mirrored gene&ndash;disease validity</h2>" + intro
    # "below" is true only above both groups.
    assert section.index(intro) < section.index("<table")
    assert section.count("<td>Autosomal dominant</td>") == 2
    assert "<td>AD</td>" not in section
    assert ">SOP11 report</a>" in section
    assert "<td>clingen</td>" not in section and "<td>gencc</td>" not in section


@pytest.mark.parametrize(
    ("source", "absent_group_sentence"),
    [
        (ValiditySource.GENCC, "No ClinGen expert panel record for this gene"),
        (ValiditySource.CLINGEN, "No Gene Curation Coalition submitter other than ClinGen"),
    ],
)
def test_an_empty_authority_group_is_a_sentence_not_an_empty_table(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    source: ValiditySource,
    absent_group_sentence: str,
) -> None:
    """21 published genes have no panel record (measured 2026-10-25), and an empty
    "Expert panel" table reads as a rendering fault rather than as that fact."""
    record = replace(
        _validity().records[0],
        source=source,
        submitter="G2P" if source is ValiditySource.GENCC else None,
        submitted_on="2021-01-01 00:00:00" if source is ValiditySource.GENCC else None,
    )
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={
            GATA4: GeneValidity(
                records=(record,),
                state=ValidityState.EXPERT_CURATED,
                has_source_discordance=False,
            )
        },
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    section = _validity_table(_page(tmp_path, "HGNC_4173.html"))
    # One table holding the group that exists; the absent group is a sentence.
    assert section.count("<table") == 1
    assert absent_group_sentence in section
    assert section.count("Gold standard</span>") == (1 if source is ValiditySource.CLINGEN else 0)


# GATA4's row of `mirrors/genes.tsv`, measured 2026-10-25.
_GATA4_IDENTITY: Final = GeneIdentity(
    name="GATA binding protein 4",
    locus="8p23.1",
    ensembl_gene="ENSG00000136574",
    ncbi_gene=2626,
    uniprot="P43694",
    mane_select="NM_001308093.3",
)


def test_a_gene_page_names_its_gene_and_links_gnomad_and_clinvar_first(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D52, sliced to the rail and then to its link row (section 4.19)."""
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        identities={GATA4: _GATA4_IDENTITY},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    rail = _slice_between(_page(tmp_path, "HGNC_4173.html"), '<aside class="rail">', ("</aside>",))
    assert "<p>GATA binding protein 4</p>" in rail
    assert "<p>8p23.1 · MANE Select NM_001308093.3</p>" in rail
    # Under the symbol and above the HGNC id, where `_rail`'s docstring puts it.
    assert rail.index("</h1>") < rail.index("<p>GATA binding protein 4</p>") < rail.index("<code>")
    links = re.findall(
        r'<a href="([^"]+)">([^<]+)</a>', _slice_between(rail, '<p class="xrefs">', ("</p>",))
    )
    assert links[:2] == [
        ("https://gnomad.broadinstitute.org/gene/ENSG00000136574?dataset=gnomad_r4", "gnomAD"),
        ("https://www.ncbi.nlm.nih.gov/clinvar/?term=2626%5Bgeneid%5D", "ClinVar"),
    ]
    assert len(links) == 5


def test_the_rails_identity_text_is_escaped(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`name`, `locus` and `mane_select` are free text in the mirror schema.

    None carries a pattern in `tables.GENES`, so each reaches the rail as
    whatever HGNC published. Until the locus and MANE cases were added, a mutant
    dropping the escape on that line passed all 1218 tests (measured by review,
    2026-10-25).
    """
    build_gene_pages(
        facts_uncurated,
        Emitter(root=tmp_path),
        symbols={GATA4: "GATA4"},
        identities={
            GATA4: replace(
                _GATA4_IDENTITY, name="<b>x</b>", locus="<i>8p</i>", mane_select="<s>NM_1</s>"
            )
        },
        validity={GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )

    rail = _slice_between(_page(tmp_path, "HGNC_4173.html"), '<aside class="rail">', ("</aside>",))
    assert "<p>&lt;b&gt;x&lt;/b&gt;</p>" in rail
    assert "<p>&lt;i&gt;8p&lt;/i&gt; · MANE Select &lt;s&gt;NM_1&lt;/s&gt;</p>" in rail
    assert "<b>" not in rail and "<i>" not in rail and "<s>" not in rail


@pytest.mark.parametrize(
    "identities",
    [
        {GATA4: replace(_GATA4_IDENTITY, ensembl_gene=None)},
        {GATA4: replace(_GATA4_IDENTITY, ncbi_gene=None)},
        {},
    ],
    ids=["no Ensembl id", "no NCBI Gene id", "absent from the registry"],
)
def test_a_gene_page_that_cannot_carry_its_required_links_is_refused(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts], identities: dict[str, GeneIdentity]
) -> None:
    """D52: refused, not published without them. "absent from the registry" is the
    case a lookup returning `None` and skipping the row would get wrong."""
    with pytest.raises(ValueError, match="D52"):
        build_gene_pages(
            facts_uncurated,
            Emitter(root=tmp_path),
            symbols={GATA4: "GATA4"},
            identities=identities,
            validity={GATA4: _validity()},
            assertions={},
            publications={},
            burden={},
            cohorts={},
            cell_types={},
        )


def test_each_study_names_the_analyses_it_ran_under_its_heading(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D49: a de novo enrichment is not a burden test, and a reader is told which.

    One study running both arms (PMID:40127276 does) names both; a CNV study
    says so. Derived from `comparator` x `variant_class`, never typed per study.
    """
    rows = [
        _burden_row(),  # PMID:42230622, control_cohort, snv_indel
        _burden_row(study="PMID:34324492", variant_class="cnv_deletion", source="audain2021"),
        _burden_row(
            study="PMID:40127276",
            comparator="mutation_model",
            count_unit="de_novo_mutations",
            n_control_carriers=None,
            n_controls=None,
            expected_count=0.4,
            effect=12.5,
            effect_measure="enrichment_ratio",
            effect_bound=None,
            ci_low=None,
            pvalue_test="poisson",
            source="sierant2025",
        ),
        _burden_row(study="PMID:40127276", source="sierant2025"),
    ]
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            rows,
            publications={
                _PUBLICATION.id: _PUBLICATION,
                **_OTHER_PUBLICATIONS,
                **_SIERANT_PUBLICATION,
            },
        )
    )
    assert _subtitles_by_study(section) == {
        "42230622": "case-control burden test",
        "34324492": "CNV burden (case-control)",
        "40127276": "case-control burden test · de novo enrichment (trios vs mutation model)",
    }


def _subtitles_by_study(section: str) -> dict[str, str]:
    """{pmid: subtitle}, matched only where the subtitle sits directly under that
    study's heading, so a subtitle under the wrong study fails too."""
    return dict(
        re.findall(
            r'<h3 id="study-\d+"><a href="https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/?">[^<]*'
            r"</a></h3>"
            r'<p class="study-analyses">([^<]+)</p>',
            section,
        )
    )


def test_a_study_running_many_analyses_lists_them_in_sorted_order(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The subtitle joins a `set`, so its order is only deterministic if sorted.

    Eight synthetic comparators fall through `_analysis_label` to the raw
    `"<token> (snv_indel)"` form, giving eight distinct labels. Probed over
    `PYTHONHASHSEED` 0..299 in subprocesses: the set's own iteration order
    equalled the sorted join in 0 of 300 runs (a two-label set manages about
    half), so removing `sorted(` fails under essentially every seed.
    """
    tokens = ["zulu", "alpha", "mike", "echo", "kilo", "bravo", "tango", "delta"]
    rows = [_burden_row(comparator=t, cohort_stratum=t) for t in tokens]
    section = _burden_section_text(_burden_page(tmp_path, facts_uncurated, rows))
    assert _subtitles_by_study(section) == {
        "42230622": " · ".join(
            f"{t} (snv_indel)"
            for t in ["alpha", "bravo", "delta", "echo", "kilo", "mike", "tango", "zulu"]
        )
    }


@pytest.mark.parametrize(
    ("comparator", "variant_class", "label"),
    [
        ("mutation_model", "snv_indel", "de novo enrichment (trios vs mutation model)"),
        ("control_cohort", "snv_indel", "case-control burden test"),
        ("control_cohort", "cnv_deletion", "CNV burden (case-control)"),
        ("control_cohort", "sv_other", "case-control burden test (sv_other)"),
        ("none", "snv_indel", "case series (no comparison group)"),
        ("a_future_comparator", "snv_indel", "a_future_comparator (snv_indel)"),
    ],
)
def test_the_analysis_label_names_what_nothing_has_taught_it_rather_than_nothing(
    comparator: str, variant_class: str, label: str
) -> None:
    """Both raw fallbacks render the unknown token verbatim, as `_count` does."""
    assert _analysis_label(comparator, variant_class) == label


# --- The clinician card (D55, D56) -------------------------------------------
#
# The rail, rebuilt as the summary a clinician reads first. Every assertion
# below is sliced to the rail (CLAUDE.md section 4.19) and pinned against a
# literal rather than the constant that produced it (section 4.38).

_ALL_STUDIES: Final = {_PUBLICATION.id: _PUBLICATION, **_OTHER_PUBLICATIONS, **_SIERANT_PUBLICATION}

_CARD_CAVEAT: Final = (
    "<p><strong>No enrichment at these cohort sizes is not evidence against a gene</strong> "
    "&mdash; rare-variant tests routinely detect nothing for genes with overwhelming family and "
    "functional evidence.</p>"
)


def _rail_of(page: str) -> str:
    return _slice_between(page, '<aside class="rail">', ("</aside>",))


def _card_variants(rail: str) -> str:
    return _slice_between(rail, '<div class="card-variants">', ("</div>",))


def _sierant(**overrides: object) -> BurdenRow:
    """PMID:40127276's case-control arm, counting alleles: NKX2-5's damaging row."""
    payload: dict[str, object] = {
        "study": "PMID:40127276",
        "cohort_stratum": "all",
        "consequence_class": "damaging",
        "count_unit": "alleles",
        "n_case_carriers": 24,
        "n_cases": 20151,
        "n_control_carriers": 123,
        "n_controls": 237963,
        "effect": 2.31,
        "effect_bound": None,
        "ci_low": None,
        "pvalue": 4.58e-4,
        "pvalue_adjusted": 9.44e-3,
        "pvalue_adjustment": "benjamini_hochberg",
        "case_cohorts": ("pcgc",),
        "control_cohorts": (),
        "source": "sierant2025",
    }
    payload.update(overrides)
    return _burden_row(**payload)


def _sierant_de_novo(**overrides: object) -> BurdenRow:
    """PMID:40127276's de novo arm: KMT2D's damaging row, against a mutation model."""
    payload: dict[str, object] = {
        "comparator": "mutation_model",
        "count_unit": "de_novo_mutations",
        "n_case_carriers": 27,
        "n_cases": 3887,
        "n_control_carriers": None,
        "n_controls": None,
        "expected_count": 0.6071494,
        "effect": 44.5,
        "effect_measure": "enrichment_ratio",
        "pvalue": 7.21e-35,
        "pvalue_test": "poisson",
        "pvalue_adjusted": 1.79e-32,
    }
    payload.update(overrides)
    return _sierant(**payload)


def _audain2021(**overrides: object) -> BurdenRow:
    """PMID:34324492's CNV deletion test: no effect measure, a family-wise corrected p."""
    payload: dict[str, object] = {
        "study": "PMID:34324492",
        "variant_class": "cnv_deletion",
        "cohort_stratum": "all",
        "n_case_carriers": 10,
        "n_cases": 6746,
        "n_control_carriers": 0,
        "n_controls": 14024,
        "effect": None,
        "effect_measure": None,
        "effect_bound": None,
        "ci_low": None,
        "pvalue": 9.999e-05,
        "pvalue_test": "permutation",
        "pvalue_adjusted": 0.0121,
        "pvalue_adjustment": "familywise_permutation",
        "source": "audain2021_s6",
    }
    payload.update(overrides)
    return _burden_row(**payload)


def test_the_card_names_each_study_and_test_by_its_first_surviving_row(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D56: one line per (study, analysis) with a surviving row, never a ranking.

    The line shows the first surviving row in the table's own reading order
    (`burden._sort_key`: stratum, then the `damaging` union before its parts),
    never the most significant one -- the atlas authors no ranking -- then how
    many more in that group survive. Counts in their own unit through `_count`,
    and a de novo row against its expectation through `_controls`.

    The fixture is built so each rule has something to get wrong. Rows arrive
    in reverse study order and with the de novo arm before the case-control
    arm, so the groups must be sorted, not iterated as met. Within each group
    the later-ordered row arrives first (`lof` before `damaging`,
    `nonsyndromic` before `syndromic`), so the first row must come from the
    sort, and the CNV group's first survivor is a stratum other than `all`,
    which the line must name. Two rows that do not survive sit in surviving
    groups (corrected 0.536; corrected 0.2) and one study publishes no
    correction at all, so a row that is merely nominal must never reach a line.

    **Each line ends with a link to its study's block**, `#study-<PMID digits>`
    on that study's `<h3>`, labelled `caveats` where any of the group's
    *surviving* rows carries a `method_note` and `details` otherwise -- the
    note is where a caveat about those counts reaches a reader (PMID:34324492's
    cases include 777 thoracic aortic aneurysm probands). The CNV group's note
    is on its second survivor, not its first, so the label reads the whole
    group; the case-control group's note is on its one row that does not
    survive, so the label reads only survivors; and its first survivor's note
    is whitespace, which is no note. Every study block carries exactly one
    anchor, including the study with no line.

    **The rest of a group's survivors are named, not counted** ("also
    surviving: loss-of-function; missense (damaging)"), in the same reading
    order: "2 more survive" beside the `damaging` union counted its own two
    components as further findings. Each label carries its stratum and lesion
    group where they are not the whole cohort, so a per-lesion row is never
    shown unlabelled. A line from a study by an atlas author says so, as the
    study's block does; each link's accessible name names its study (WCAG
    2.4.9). And PMID:42230622, which publishes no corrected p for this gene,
    gets its own line saying so, with the scan size its block judges it by --
    never a threshold the atlas computed (D12/D33).
    """
    rows = [
        _burden_row(),  # PMID:42230622: p 3.13e-08, no correction published.
        _sierant_de_novo(
            consequence_class="missense_damaging",
            n_case_carriers=6,
            expected_count=0.44,
            pvalue_adjusted=5.84e-4,
        ),
        _sierant_de_novo(
            consequence_class="lof",
            n_case_carriers=21,
            expected_count=0.166,
            pvalue_adjusted=1.65e-34,
        ),
        _sierant_de_novo(),
        _sierant(
            consequence_class="missense_damaging",
            pvalue=0.0925,
            pvalue_adjusted=0.536,
            method_note="A note on a row that does not survive.",
        ),
        _sierant(
            consequence_class="lof",
            n_case_carriers=11,
            n_control_carriers=25,
            n_controls=222440,
            pvalue_adjusted=3.01e-3,
        ),
        _sierant(method_note="   "),
        _audain2021(pvalue_adjusted=0.2),
        _audain2021(
            cohort_stratum="nonsyndromic",
            lesion_group="conotruncal",
            n_case_carriers=6,
            n_cases=2405,
            method_note="Cases include 777 thoracic aortic aneurysm probands.",
        ),
        _audain2021(cohort_stratum="syndromic", n_case_carriers=4, n_cases=1471),
    ]
    page = _burden_page(tmp_path, facts_uncurated, rows, publications=_ALL_STUDIES)
    rail = _rail_of(page)
    card = _card_variants(rail)

    assert card.startswith(
        "<p><strong>Rare variants</strong> surviving their own study's correction</p><ul>"
    )
    assert re.findall(r"<li>(.*?)</li>", card) == [
        "Audain E et al. 2021 · CNV burden (case-control) — loss-of-function, syndromic: "
        "4 carriers / 1,471 vs 0 carriers / 14,024 · also surviving: loss-of-function, "
        "non-syndromic, lesion group conotruncal · study by an atlas author · "
        '<a href="#study-34324492" aria-label="Caveats for Audain E et al. 2021">caveats</a>',
        "Sierant MC et al. 2025 · case-control burden test — damaging (LOF + missense): "
        "24 alleles / 20,151 vs 123 alleles / 237,963 · also surviving: loss-of-function · "
        '<a href="#study-40127276" aria-label="Details for Sierant MC et al. 2025">details</a>',
        "Sierant MC et al. 2025 · de novo enrichment (trios vs mutation model) — damaging "
        "(LOF + missense): 27 de novo / 3,887 trios vs 0.607 expected · also surviving: "
        "loss-of-function; missense (damaging) · "
        '<a href="#study-40127276" aria-label="Details for Sierant MC et al. 2025">details</a>',
    ]
    assert "more survive" not in card
    after = card[card.index("</ul>") + len("</ul>") :]
    assert after == (
        '<p class="card-uncorrected"><a href="#study-42230622">Audain E et al. 2026</a> '
        "publishes no corrected p for this gene; judge its results against its scan of "
        f"138,609 tests (see its block) · study by an atlas author</p>{_CARD_CAVEAT}"
    )

    # One anchor per study block, on that study's own heading, and nowhere else.
    section = _burden_section_text(page)
    anchors = re.findall(
        r'<h3 id="([^"]+)"><a href="https://pubmed\.ncbi\.nlm\.nih\.gov/(\d+)/">', section
    )
    assert anchors == [
        ("study-34324492", "34324492"),
        ("study-40127276", "40127276"),
        ("study-42230622", "42230622"),
    ]
    for anchor, _ in anchors:
        assert page.count(f'id="{anchor}"') == 1, anchor


def test_the_card_states_a_studys_summary_caveat_once_beside_its_results(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """A curated `Publication.summary_caveat` follows the card's list, once per study.

    The card shows a study's counts and the caveat that says what they count
    sat one click away, in the study's block: PMID:34324492's cases include 777
    thoracic aortic aneurysm probands, and a carrier there is anyone with a
    rare deletion overlapping the gene. A count that reaches the card must not
    shed what qualifies it on the way.

    **Once per study with a line, in study-id order**: PMID:40127276 has two
    lines (two analyses) and one note, and its rows arrive first. **Never for a
    study with no line**: PMID:42230622 carries a caveat and only a
    no-correction line, which shows no counts to qualify, and a study whose
    rows all fail their correction shows no counts either. The caveat and
    the study's label are curator-written text, so both are escaped.

    **Then what S8 Table says of the CHD cases among the line's carriers**
    (owner decision R2, 2026-10-26), from the line's row's
    `n_chd_case_carriers`: "every carrier", and so no aneurysm proband, where
    it equals the row's case carriers (10 of 10 here), the two counts where it
    is fewer (7 of 10), and nothing where the row carries none -- the Sierant
    rows here, and 109 of PMID:34324492's 121 mirrored rows. The unequal case
    is on no published page (ANKRD11's 12 of 13 is the only one, and ANKRD11
    is unpublished), so this fixture is the only place "every" is held to
    equality. **From the line's own row**, the first of its group, never
    another the group holds; **once per study**, however many of its lines
    carry the same sentence; **never for a row with no case carrier**, of
    whom "every" says nothing; and **never without the caveat**, whose
    aneurysm probands it speaks of (all four review, 2026-10-27).
    """
    publications = {
        "PMID:42230622": _PUBLICATION.model_copy(
            update={"summary_caveat": "A caveat for a study with no line."}
        ),
        "PMID:34324492": _OTHER_PUBLICATIONS["PMID:34324492"].model_copy(
            update={"summary_caveat": "Cases include <777> TAA probands & CHD."}
        ),
        "PMID:40127276": _SIERANT_PUBLICATION["PMID:40127276"].model_copy(
            update={"authors": ["D'Alessandro M"], "summary_caveat": "Counts alleles."}
        ),
    }
    rows = [_sierant_de_novo(), _sierant(), _audain2021(n_chd_case_carriers=10), _burden_row()]
    page = _burden_page(tmp_path, facts_uncurated, rows, publications=publications)
    card = _card_variants(_rail_of(page))

    assert card.count("<li>") == 3
    after = card[card.index("</ul>") + len("</ul>") :]
    assert after == (
        '<p class="card-note"><a href="#study-34324492">Audain E et al. 2021</a>: '
        "Cases include &lt;777&gt; TAA probands &amp; CHD. For this gene, the study&#x27;s "
        "S8 Table counts every carrier among its CHD cases, so none is an aneurysm "
        "proband.</p>"
        '<p class="card-note"><a href="#study-40127276">D&#x27;Alessandro M et al. 2025</a>: '
        "Counts alleles.</p>"
        '<p class="card-uncorrected"><a href="#study-42230622">Audain E et al. 2026</a> '
        "publishes no corrected p for this gene; judge its results against its scan of "
        f"138,609 tests (see its block) · study by an atlas author</p>{_CARD_CAVEAT}"
    )

    failing = _burden_page(
        tmp_path / "failing",
        facts_uncurated,
        [_sierant(pvalue_adjusted=0.5)],
        publications=publications,
    )
    assert 'class="card-note"' not in _card_variants(_rail_of(failing))

    # Two lines of one study. The CNV line's own row is `damaging` (7 of 10),
    # which sorts before the `lof` row its group also holds (10 of 10); the
    # SNV line's row says what the first does, and the sentence appears once.
    fewer = _burden_page(
        tmp_path / "fewer",
        facts_uncurated,
        [
            _audain2021(n_chd_case_carriers=10),
            _audain2021(consequence_class="damaging", n_chd_case_carriers=7),
            _audain2021(variant_class="snv_indel", n_chd_case_carriers=7),
        ],
        publications=publications,
    )
    card = _card_variants(_rail_of(fewer))
    assert card.count("<li>") == 2
    assert re.findall(r'<p class="card-note">.*?</p>', card) == [
        '<p class="card-note"><a href="#study-34324492">Audain E et al. 2021</a>: '
        "Cases include &lt;777&gt; TAA probands &amp; CHD. For this gene, the study&#x27;s "
        "S8 Table counts 7 of the 10 carriers among its CHD cases.</p>"
    ]
    assert _chd_only([_audain2021(n_case_carriers=0, n_chd_case_carriers=0)]) == []
    uncaveated = _burden_page(
        tmp_path / "uncaveated",
        facts_uncurated,
        [_audain2021(n_chd_case_carriers=10)],
        publications={
            **publications,
            "PMID:34324492": publications["PMID:34324492"].model_copy(
                update={"summary_caveat": None}
            ),
        },
    )
    assert 'class="card-note"' not in _card_variants(_rail_of(uncaveated))
    assert "S8 Table" not in uncaveated

    # The table is named from `CHD_CASE_SPLITS` by the row's study. A count on a
    # row of a study it does not list is BUR024's to refuse, and past a
    # bypassed gate the card raises rather than name a table that study lacks.
    with pytest.raises(ValueError, match="a row of PMID:40127276 carries n_chd_case_carriers"):
        _burden_page(
            tmp_path / "other",
            facts_uncurated,
            [_sierant(n_chd_case_carriers=12)],
            publications=publications,
        )


def test_the_card_says_how_wide_a_survival_is_where_its_study_records_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The panel a q was corrected over, and the study's own genome-wide line (R3).

    Owner decision 2026-10-26: "surviving its own study's correction" read as
    genome-wide, and PMID:40127276's q is a Benjamini-Hochberg over its 248
    panel genes. Where a publication records `correction_genes`, the line's
    row carries its q and that scope; where it records `genome_wide_threshold`,
    every surviving row on the line strictly below it is marked, the named ones
    too. **The boundary is strict**, as the paper's "P-value < 2.6 x 10^-6"
    reads: the line's row here sits exactly on the threshold and is not
    marked, while the named row just below it is. Real rows (TBX5's card) are
    pinned in `test_built_site_is_consumable.py`; this is the fixture that can
    put a row on the line, which no real row is.

    A study recording neither field renders exactly as before, which the
    fixture literals in `test_the_card_names_each_study_and_test_by_its_first_
    surviving_row` hold: its PMID:40127276 publication records neither, and its
    de novo rows sit far below 2.6e-6.
    """
    publications = {
        "PMID:40127276": _SIERANT_PUBLICATION["PMID:40127276"].model_copy(
            update={"correction_genes": 248, "genome_wide_threshold": 2.6e-6}
        ),
    }
    rows = [
        _sierant_de_novo(pvalue=2.6e-6, pvalue_adjusted=1.2e-4),
        _sierant_de_novo(
            consequence_class="lof", n_case_carriers=21, pvalue=2.59e-6, pvalue_adjusted=1.1e-4
        ),
    ]
    page = _burden_page(tmp_path, facts_uncurated, rows, publications=publications)

    assert re.findall(r"<li>(.*?) · <a ", _card_variants(_rail_of(page))) == [
        "Sierant MC et al. 2025 · de novo enrichment (trios vs mutation model) — damaging "
        "(LOF + missense): 27 de novo / 3,887 trios vs 0.607 expected · q 0.00012 within its "
        "248-gene panel · also surviving: loss-of-function (genome-wide, p 2.59e-06)"
    ]


def test_the_card_and_the_matrix_read_one_rule_for_surviving_a_correction(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The same rows through both renderers (CLAUDE.md section 4.28).

    PMID:40127276's row clears its own Benjamini-Hochberg correction and the
    matrix fills its cell `corrected`; PMID:42230622's is nominally significant
    and its study publishes no correction, so the matrix draws it `nominal`. The
    card must list the first and not the second -- a card reading `_supports`
    alone would list both, and would then disagree with the matrix on the same
    page.
    """
    rows = [_burden_row(cohort_stratum="all", effect=3.1, effect_bound=None), _sierant()]
    page = _burden_page(tmp_path, facts_uncurated, rows, publications=_ALL_STUDIES)
    matrix = _slice_between(_burden_section_text(page), '<table class="matrix">', ("</table>",))
    by_family = dict(re.findall(r'<th scope="row">([^<]+)</th>.*?class="cell (\w[\w-]*)"', matrix))
    assert by_family == {"Audain E et al. 2026": "nominal", "Sierant MC et al. 2025": "corrected"}

    lines = re.findall(r"<li>(.*?)</li>", _card_variants(_rail_of(page)))
    assert [line.split(" · ")[0] for line in lines] == ["Sierant MC et al. 2025"]


@pytest.mark.parametrize(
    ("rows", "cell"),
    [
        pytest.param(
            [
                _sierant_de_novo(effect=15.8, pvalue=0.0613, pvalue_adjusted=0.200),
                _sierant_de_novo(
                    consequence_class="missense_damaging",
                    effect=20.8,
                    pvalue=0.0468,
                    pvalue_adjusted=0.270,
                ),
                _sierant_de_novo(
                    consequence_class="lof", effect=0.0, pvalue=1.0, pvalue_adjusted=1.0
                ),
            ],
            ("nominal", "enrichment 20.8; 0.27 (Benjamini-Hochberg)", "enrichment 20.8", "q 0.27"),
            id="GJA1: nominal by its missense row, not by the lower-q damaging union",
        ),
        pytest.param(
            [
                _sierant(),
                _sierant(
                    consequence_class="lof",
                    effect=4.2,
                    pvalue=1.2e-5,
                    pvalue_adjusted=None,
                    pvalue_adjustment=None,
                ),
            ],
            ("corrected", "OR 2.31; 0.00944 (Benjamini-Hochberg)", "OR 2.31", "q 0.00944"),
            id="corrected by its survivor, not by a lower raw p no correction cleared",
        ),
        pytest.param(
            [
                _sierant_de_novo(effect=40.45, pvalue=3.68e-06, pvalue_adjusted=9.1264e-05),
                _sierant_de_novo(
                    consequence_class="missense_damaging",
                    effect=54.1,
                    pvalue=1.17e-06,
                    pvalue_adjusted=0.00014508,
                ),
            ],
            (
                "corrected",
                "enrichment 40.5; 9.13e-05 (Benjamini-Hochberg)",
                "enrichment 40.5",
                "q 9.13e-05",
            ),
            id="RAF1: two survivors, named by the lower q, not the lower raw p",
        ),
    ],
)
def test_the_matrix_cell_names_a_row_its_state_is_true_of(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    rows: list[BurdenRow],
    cell: tuple[str, str, str, str],
) -> None:
    """The cell's glyph and the row it names come from one rule (2026-10-26).

    The cell named the most significant row in it -- corrected p where
    published, raw p otherwise -- whatever that row's role. GJA1's de novo cell
    in PMID:40127276 is `nominal` because its missense row is (p 0.0468, q
    0.270), and it named the `damaging` union (p 0.0613, q 0.200) under the
    "enriched nominally" glyph: a row that is not nominally significant at all.
    1 of 109 enriched cells on a real build. The fixture is GJA1's three rows.

    The second case is the same defect under the `corrected` glyph: a
    supporting row whose study published no correction sorts on its raw p, and
    1.2e-5 sorts before the survivor's q of 0.00944. It names the survivor.
    Both read `concordance.deciding_rows`, which `family_state` reads its tier
    from, so the glyph and the named row cannot come apart again.

    **Among rows that all decide, the corrected p ranks them, not the raw one.**
    The third case is RAF1's de novo cell: its `damaging` row has the lower q
    (9.13e-05 against 1.45e-04) and its missense row the lower raw p (1.17e-06
    against 3.68e-06), so a sort on raw p alone names the other row. Measured
    on a real build (review, 2026-10-26): that mutant changes the named row in
    8 of 109 enriched cells, and survived the full suite until this case.
    """
    page = _burden_page(tmp_path, facts_uncurated, rows, publications=_ALL_STUDIES)
    matrix = _slice_between(_burden_section_text(page), '<table class="matrix">', ("</table>",))

    assert _MATRIX_CELL.findall(matrix) == [cell]


_NO_TESTS_REPORTED: Final = {
    **_ALL_STUDIES,
    "PMID:34324492": _OTHER_PUBLICATIONS["PMID:34324492"].model_copy(
        update={"tests_reported": None}
    ),
}


@pytest.mark.parametrize(
    ("rows", "publications", "closing", "uncorrected"),
    [
        pytest.param(
            [
                _burden_row(cohort_stratum="all", effect=3.1, effect_bound=None),
                _audain2021(case_cohorts=("ddd",), pvalue_adjusted=0.2),
                _sierant(pvalue_adjusted=0.2),
                _sierant(consequence_class="synonymous", pvalue=1e-6, pvalue_adjusted=0.001),
            ],
            _ALL_STUDIES,
            "Tested in 2 independent datasets (as far as the cohort registry records) (2 of the "
            "3 studies by an atlas author); none survives a correction its study published.",
            [
                '<a href="#study-42230622">Audain E et al. 2026</a> publishes no corrected p '
                "for this gene; judge its results against its scan of 138,609 tests (see its "
                "block) · study by an "
                "atlas author"
            ],
            id="three studies in two families",
        ),
        pytest.param(
            [
                _sierant(pvalue=0.0934, pvalue_adjusted=0.256),
                _burden_row(
                    cohort_stratum="all",
                    consequence_class="missense_damaging",
                    effect=10.2,
                    effect_bound=None,
                    pvalue=1.6e-08,
                ),
                _burden_row(
                    consequence_class="missense_damaging",
                    effect=21.2,
                    effect_bound=None,
                    pvalue=2.47e-10,
                ),
            ],
            _ALL_STUDIES,
            "Tested in 2 independent datasets (as far as the cohort registry records) (1 of the "
            "2 studies by an atlas author); none survives a correction its study published.",
            [
                '<a href="#study-42230622">Audain E et al. 2026</a> publishes no corrected p '
                "for this gene; judge its results against its scan of 138,609 tests (see its "
                "block) · study by an "
                "atlas author"
            ],
            id="SMAD4: p 2.5e-10 in a study that publishes no correction",
        ),
        pytest.param(
            [_burden_row(cohort_stratum="all", effect=3.1, effect_bound=None)],
            _ALL_STUDIES,
            "Tested in 1 independent dataset (as far as the cohort registry records) (a study by "
            "an atlas author).",
            [
                '<a href="#study-42230622">Audain E et al. 2026</a> publishes no corrected p '
                "for this gene; judge its results against its scan of 138,609 tests (see its "
                "block) · study by an "
                "atlas author"
            ],
            id="one study",
        ),
        pytest.param(
            [_audain2021(case_cohorts=("ddd",), pvalue_adjusted=None)],
            _NO_TESTS_REPORTED,
            "Tested in 1 independent dataset (as far as the cohort registry records) (a study by "
            "an atlas author).",
            [
                '<a href="#study-34324492">Audain E et al. 2021</a> publishes no corrected p '
                "for this gene; see its block · study by an atlas author"
            ],
            id="no scan size recorded: no number invented",
        ),
        pytest.param(
            [_audain2021(case_cohorts=("ddd",), pvalue_adjusted=0.4)],
            _ALL_STUDIES,
            "Tested in 1 independent dataset (as far as the cohort registry records) (a study by "
            "an atlas author); none survives a correction its study published.",
            [],
            id="MYH11: one study by an atlas author, its correction failed",
        ),
        pytest.param(
            [
                _audain2021(case_cohorts=("ddd",), pvalue_adjusted=0.4),
                _burden_row(cohort_stratum="all", effect=3.1, effect_bound=None),
            ],
            _ALL_STUDIES,
            "Tested in 1 independent dataset (as far as the cohort registry records) (each study "
            "by an atlas author); none survives a correction its study published.",
            [
                '<a href="#study-42230622">Audain E et al. 2026</a> publishes no corrected p '
                "for this gene; judge its results against its scan of 138,609 tests (see its "
                "block) · study by an atlas author"
            ],
            id="two studies by an atlas author, one dataset: counted in studies",
        ),
        pytest.param(
            [_sierant(pvalue_adjusted=0.2)],
            _ALL_STUDIES,
            "Tested in 1 independent dataset (as far as the cohort registry records); none "
            "survives a correction its study published.",
            [],
            id="no study by an atlas author: no declaration",
        ),
    ],
)
def test_a_card_with_no_surviving_result_counts_the_datasets_that_tested_the_gene(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    rows: list[BurdenRow],
    publications: dict[str, Publication],
    closing: str,
    uncorrected: list[str],
) -> None:
    """Counted in cohort families, never in studies, and only those that tested it.

    PMID:42230622 and PMID:34324492 share a sample collection here (`ddd`), so
    three studies are two independent datasets; a third family that tested
    nothing for this gene is passed in and must not be counted. No row
    survives: one is nominal in a study publishing no correction, one fails its
    correction (0.2), and one is the synonymous negative control clearing a
    correction (0.001), which is not support. The power caveat stays: it is the
    sentence that stops "none survives" reading as evidence against the gene.

    **"None survives" may only speak of corrections a study published.**
    PMID:42230622 publishes no corrected p on any row, and the card read "none
    survives its own study's correction" for SMAD4 -- whose block tells the
    reader to judge raw p against the study's scan of 138,609 tests, which two
    SMAD4 rows (p 1.6e-08, 2.47e-10) would clear. Measured 2026-10-26: 24
    published genes are tested by that study alone. So the closing sentence
    says "a correction its study published", and each study that published
    none for this gene gets a line saying so, naming its scan size where
    `Publication.tests_reported` records it and inventing none where it does
    not. Whether a study published a correction is decided per gene, from its
    rows: any `pvalue_adjusted`. Where no study that tested the gene published
    one, the closing sentence says only how many tested it: "none survives a
    correction its study published" would be vacuously true there, and read as
    a correction applied and failed.

    **The closing sentence declares a study by an atlas author** (owner
    decision 2026-10-26, R7): it names no study, so MYH11's card, tested by
    PMID:34324492 alone, carried no declaration at all. Three forms, each
    pinned: one study, all of several, some of several -- counted in studies,
    not datasets, which the two-studies-one-family case is here to tell apart.
    A study by no atlas author adds nothing.

    **The count carries the front page's qualifier** (owner answer
    2026-10-27): "(as far as the cohort registry records)", straight after
    the count and before any declaration, in every form here -- singular and
    plural, with each declaration and with none. Pinned as literals, never as
    the constant that renders it.
    """
    families = (
        frozenset({"PMID:34324492", "PMID:42230622"}),
        frozenset({"PMID:40127276"}),
        frozenset({"PMID:99999999"}),
    )
    page = _burden_page(
        tmp_path, facts_uncurated, rows, publications=publications, families=families
    )
    card = _card_variants(_rail_of(page))

    assert card == (
        f"<p><strong>Rare variants</strong></p><p>{closing}</p>"
        + "".join(f'<p class="card-uncorrected">{line}</p>' for line in uncorrected)
        + _CARD_CAVEAT
    )
    assert "<li>" not in card
    assert "none survives its own study" not in page


def test_a_card_with_no_rare_variant_row_says_so_and_claims_no_test(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """CRIPTO's case: no study in this atlas publishes a row for the gene.

    Says what is known -- no study here reports a result -- and not "no study
    tested it", which a study whose supplement simply omits the gene cannot be
    said not to have done. No power caveat: there is no result for it to
    qualify, and no rare-variant section in the main column either.
    """
    page = _burden_page(tmp_path, facts_uncurated, [])
    card = _card_variants(_rail_of(page))

    assert card == (
        "<p><strong>Rare variants</strong></p>"
        "<p>No rare-variant study in this atlas reports a result for this gene.</p>"
    )
    assert "not evidence against a gene" not in _rail_of(page)
    assert _RARE_VARIANT_HEADING not in page


def _clingen(
    classification: Classification,
    disease: str,
    label: str,
    moi: str = "AD",
    gcep: str = "Congenital Heart Disease Gene Curation Expert Panel",
) -> ValidityRecord:
    return ValidityRecord(
        source=ValiditySource.CLINGEN,
        classification=classification,
        classification_term=classification.value.capitalize(),
        disease=disease,
        disease_label=label,
        moi=moi,
        report_url=None,
        sop="SOP11",
        gcep=gcep,
    )


def _gencc(
    submitter: str,
    classification: Classification | None,
    moi: str = "Autosomal dominant",
    disease: str = "MONDO:0007173",
    label: str = "atrial septal defect 7",
) -> ValidityRecord:
    return ValidityRecord(
        source=ValiditySource.GENCC,
        classification=classification,
        classification_term=classification.value if classification else "Supportive",
        disease=disease,
        disease_label=label,
        moi=moi,
        report_url=None,
        submitter=submitter,
        submitted_on="2020-08-20 13:02:04",
    )


def _card_page(
    tmp_path: Path,
    fact: GeneFacts,
    validity: GeneValidity,
    warrants: dict[str, Warrant] | None = None,
    scope_terms: dict[str, ScopeEntry] | None = None,
) -> str:
    """One gene page, rendered with no burden rows; returns its rail.

    `validity` also holds EVC, graded by the panel that scopes Ellis-van Creveld
    syndrome and never rendered: the gene a card names when that panel did not
    grade the page's own gene (`_card_scope`).
    """
    build_gene_pages(
        {GATA4: fact},
        Emitter(root=tmp_path),
        symbols={GATA4: "GATA4", _EVC_GENE: "EVC"},
        validity={GATA4: validity, _EVC_GENE: _EVC_VALIDITY},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        scope_warrants=warrants,
        scope_terms=scope_terms,
        cell_types={},
    )
    return _rail_of(_page(tmp_path, "HGNC_4173.html"))


def _graded(grade: Classification, *records: ValidityRecord) -> tuple[GeneFacts, GeneValidity]:
    """A ClinGen-admitted gene headlined `grade`, carrying `records`."""
    return (
        replace(
            _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=grade,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=records, state=ValidityState.EXPERT_CURATED, has_source_discordance=False
        ),
    )


def _agreed(*records: ValidityRecord) -> tuple[GeneFacts, GeneValidity]:
    """A gene no panel graded, admitted on GenCC submitters' agreement."""
    return (
        replace(
            _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=None,
            validity_state=ValidityState.SUBMITTER_CURATED,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=records, state=ValidityState.SUBMITTER_CURATED, has_source_discordance=False
        ),
    )


_NKX2_5_DISEASE: Final = "NKX2.5-related congenital, conduction and myopathic heart disease"


@pytest.mark.parametrize(
    ("gene", "rows"),
    [
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, "MONDO:0800441", _NKX2_5_DISEASE),
                # Below the headline: neither its disease nor its mode is the chip's.
                _clingen(Classification.LIMITED, "MONDO:0007173", "atrial septal defect 7", "XL"),
                _gencc("ClinGen", Classification.DEFINITIVE),
                _gencc("Ambry Genetics", Classification.STRONG),
                _gencc("G2P", Classification.DEFINITIVE),
                _gencc("Genomics England PanelApp", Classification.STRONG),
                _gencc("Labcorp Genetics (formerly Invitae)", Classification.STRONG),
                _gencc("Laboratory for Molecular Medicine", Classification.MODERATE),
                _gencc("Orphanet", None),
                _gencc("PanelApp Australia", Classification.STRONG, "Autosomal recessive"),
                _gencc("Illumina", Classification.NO_KNOWN_ASSOCIATION, "Unknown"),
            ),
            (
                ("definitive for", _NKX2_5_DISEASE),
                ("also listed by", "7 other curators (GenCC)"),
                ("inheritance", "Autosomal dominant"),
            ),
            id="NKX2-5: ClinGen and a dissent are not other curators",
        ),
        pytest.param(
            _graded(
                Classification.LIMITED,
                *(
                    _clingen(Classification.LIMITED, "MONDO:0005453", "congenital heart disease", m)
                    for m in ("XL", "AR", "SD", "UD", "AD", "MT")
                ),
            ),
            (
                ("limited for", "congenital heart disease"),
                (
                    "inheritance",
                    "Autosomal dominant / Autosomal recessive / Mitochondrial / Semidominant "
                    "/ Unknown / X-linked",
                ),
            ),
            id="ROBO1: one grade under several modes names every one",
        ),
        pytest.param(
            _agreed(
                _gencc("G2P", Classification.MODERATE, "Autosomal recessive"),
                _gencc("Ambry Genetics", Classification.LIMITED),
                _gencc("Orphanet", None, "X-linked"),
                _gencc("Illumina", Classification.NO_KNOWN_ASSOCIATION, "Unknown"),
            ),
            (
                ("admitted by", "Ambry Genetics; G2P"),
                ("also listed by", "1 other curator (GenCC)"),
                ("inheritance", "Autosomal dominant / Autosomal recessive"),
            ),
            id="GenCC-admitted: the modes the admitting records carry",
        ),
        pytest.param(
            _agreed(
                _gencc("Ambry Genetics", Classification.LIMITED),
                _gencc("Labcorp Genetics (formerly Invitae)", Classification.STRONG),
            ),
            (
                ("admitted by", "Ambry Genetics; Labcorp Genetics (formerly Invitae)"),
                ("inheritance", "Autosomal dominant"),
            ),
            id="TBX1: nobody else lists it, so no row says so",
        ),
    ],
)
def test_the_cards_rows_say_who_else_lists_the_gene_and_under_what_inheritance(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], rows: tuple[tuple[str, str], ...]
) -> None:
    """D55: the warrant row, then who else lists the gene, then its inheritance.

    **"Also listed by" counts `asserted_by.institutions` -- the list the bundle
    publishes -- minus ClinGen and minus the submitters that admitted the
    gene.** ClinGen is excluded under both names it appears as (its panel's own
    record, its GenCC copy), and a dissent is not a listing: Illumina's `No
    Known Disease Relationship` counts for nobody. `Supportive` does count --
    Orphanet asserts without grading.

    **Inheritance is the mode of the records that carry the warrant**: the
    ClinGen records at the headline grade, or the admitting submitters' records
    at or above the floor. Not every record's: the NKX2-5 fixture holds an
    `X-linked` record below its headline and an `Autosomal recessive` GenCC one,
    and neither is the mode of the classification the chip shows.

    **Every mode at the headline, not the first.** ROBO1, PRKD1 and LZTR1 each
    carry two ClinGen records at their headline grade for one disease, one
    `AD` and one `AR` (measured 2026-10-26), so the admitting record alone
    would state half of what the panel graded. All six of ClinGen's codes here,
    sorted by the words a reader sees. Six, not five, because the modes are
    de-duplicated through a `set`: probed as the six labels under
    `PYTHONHASHSEED` 0..499, one process each, the set iterated in sorted order
    on 0 of 500 seeds; the same five without `Unknown` did so on 2 of 300, a
    dropped sort surviving about one run in 150 (CLAUDE.md section 4.42).
    """
    fact, validity = gene
    rail = _card_page(tmp_path, fact, validity)
    listed = _slice_between(rail, "<dl>", ("</dl>",))
    assert listed == "".join(f"<dt>{dt}</dt><dd>{dd}</dd>" for dt, dd in rows)


@pytest.mark.parametrize(
    ("gene", "gloss"),
    [
        pytest.param(
            _graded(
                Classification.LIMITED,
                _clingen(Classification.LIMITED, "MONDO:0005453", "congenital heart disease"),
            ),
            '<p class="card-gloss"><em>Limited</em>: few probands, or evidence that is '
            "suggestive but not compelling. <strong>Limited is not a weak yes</strong> &mdash; "
            "it is a panel saying the case is not yet made.</p>",
            id="limited",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, "MONDO:0800441", _NKX2_5_DISEASE),
            ),
            '<p class="card-gloss"><em>Definitive</em>: repeatedly demonstrated in research and '
            "clinical settings and upheld over time, with no convincing evidence to the "
            "contrary.</p>",
            id="definitive",
        ),
        pytest.param(
            _agreed(
                _gencc("Ambry Genetics", Classification.LIMITED),
                _gencc("G2P", Classification.DEFINITIVE),
            ),
            None,
            id="ungraded: no grade, so no gloss",
        ),
    ],
)
def test_the_card_glosses_the_headline_grade_in_the_grade_keys_own_words(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], gloss: str | None
) -> None:
    """One line under the chip saying what the panel's rung means (D55).

    `limited` is the plurality (45 of 119, measured 2026-10-26) and the rung a
    reader is likeliest to misread as a weak yes, so it is pinned in full. The
    words are the grade key's, so the card and the key cannot come to gloss one
    rung two ways; a gene no panel graded has no rung to gloss, and its card
    says why it is published in the curation line instead.
    """
    fact, validity = gene
    rail = _card_page(tmp_path, fact, validity)
    if gloss is None:
        assert "card-gloss" not in rail
    else:
        assert gloss in rail
        assert rail.count("card-gloss") == 1


_RASOPATHY: Final = "RASopathy Gene Curation Expert Panel"


def _rasopathy(classification: Classification, disease: str, label: str) -> ValidityRecord:
    """A RASopathy-panel record, dated so ClinGen's own GenCC copy of it can fold."""
    return replace(
        _clingen(classification, disease, label),
        gcep=_RASOPATHY,
        classification_date="2024-11-07T17:00:00.000Z",
    )


@pytest.mark.parametrize(
    ("gene", "line"),
    [
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _rasopathy(Classification.DEFINITIVE, "MONDO:0018997", "Noonan syndrome"),
                _rasopathy(
                    Classification.DISPUTED, "MONDO:0015280", "cardiofaciocutaneous syndrome"
                ),
                _rasopathy(Classification.DISPUTED, "MONDO:0009026", "Costello syndrome"),
                # ClinGen's own GenCC copy of the Costello record: the table folds
                # it into the panel's row, so the card must not name it twice.
                replace(
                    _gencc(
                        "ClinGen",
                        Classification.DISPUTED,
                        disease="MONDO:0009026",
                        label="Costello syndrome",
                    ),
                    classification_term="Disputed Evidence",
                    submitted_on="2024-11-07 00:00:00",
                ),
                replace(
                    _gencc(
                        "Illumina",
                        Classification.DISPUTED,
                        disease="MONDO:0018997",
                        label="Noonan syndrome",
                    ),
                    classification_term="Disputed Evidence",
                ),
            ),
            "Disputed for Costello syndrome and cardiofaciocutaneous syndrome by the RASopathy "
            "Gene Curation Expert Panel; Disputed Evidence for Noonan syndrome by Illumina; see "
            "the validity table.",
            id="PTPN11: Definitive for Noonan, Disputed for Costello and CFC",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _rasopathy(Classification.DEFINITIVE, "MONDO:0018997", "Noonan syndrome"),
            ),
            None,
            id="graded, nothing disputed",
        ),
        pytest.param(
            _agreed(
                _gencc("Ambry Genetics", Classification.LIMITED),
                _gencc("G2P", Classification.DEFINITIVE),
                replace(
                    _gencc("Illumina", Classification.DISPUTED), classification_term="Disputed"
                ),
                replace(
                    _gencc(
                        "Ambry Genetics",
                        Classification.REFUTED,
                        disease="MONDO:0004933",
                        label="hypoplastic left heart syndrome",
                    ),
                    classification_term="Refuted Evidence",
                ),
            ),
            "Refuted Evidence for hypoplastic left heart syndrome by Ambry Genetics; Disputed for "
            "atrial septal defect 7 by Illumina; see the validity table.",
            id="no grade to gloss: the line takes the gloss's place",
        ),
    ],
)
def test_the_card_names_who_disputes_the_gene_and_for_what_beside_the_grade(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], line: str | None
) -> None:
    """`has_conflicting_evidence` on the card, in words, under the gloss (D55).

    The flag means exactly this (`vocab.has_conflicting_evidence`, over every
    in-scope record `derive.gene_facts` mirrors): some record is `Disputed` or
    `Refuted` *and* some record supports the gene. A line saying only that
    "an authority" disputes the gene told a reader neither who nor for what;
    PTPN11 is Definitive for Noonan syndrome and Disputed by the same panel for
    Costello and cardiofaciocutaneous syndrome, beside a Definitive gloss
    reading "no convincing evidence to the contrary". So the line names them:
    each authority's own term, verbatim, the diseases it applies it to, and
    the authority -- the panel for a ClinGen record, the submitter for a GenCC
    one. Panels first, then submitters, each by name. It names exactly the
    contested records the validity table lists: ClinGen's own GenCC copy of a
    panel record folds there (`_repeats_the_panel`) and is not named twice.
    The fixture lists cardiofaciocutaneous before Costello, so "Costello
    syndrome and cardiofaciocutaneous syndrome" can only come from the sort;
    the diseases are de-duplicated in record order, not through a `set`, so a
    dropped sort fails on every seed rather than on half of them.

    Renamed 2026-10-26 from
    `test_the_card_says_when_an_authority_disputes_the_gene_beside_the_grade`,
    when the line began naming its authorities.
    """
    fact, validity = gene
    rail = _card_page(tmp_path, replace(fact, has_conflicting_evidence=line is not None), validity)
    if line is None:
        assert "card-dispute" not in rail
        assert "Disputed" not in rail
        return
    expected = f'<p class="card-dispute">{line}</p>'
    assert rail.count('<p class="card-dispute">') == 1
    # Directly after the gloss, or where the gloss would sit.
    after = rail[rail.index("</dl>") + len("</dl>") :]
    if 'class="card-gloss"' in after:
        after = after[after.index("</p>") + len("</p>") :]
    assert after.startswith(expected)


def test_a_dispute_flag_with_no_contested_record_to_name_is_refused(tmp_path: Path) -> None:
    """The flag and the records come from one derivation; disagreeing, they are a bug.

    A card naming no authority would be the vague line this replaced, and one
    naming an authority the table does not list would be false, so the page is
    refused rather than published.
    """
    fact, validity = _graded(
        Classification.DEFINITIVE,
        _rasopathy(Classification.DEFINITIVE, "MONDO:0018997", "Noonan syndrome"),
    )
    with pytest.raises(ValueError, match="no contested record"):
        _card_page(tmp_path, replace(fact, has_conflicting_evidence=True), validity)


def _no_association(submitter: str, disease: str, label: str, on: str | None) -> ValidityRecord:
    """A GenCC submitter's `No Known Disease Relationship`, dated `on`."""
    return replace(
        _gencc(submitter, Classification.NO_KNOWN_ASSOCIATION, disease=disease, label=label),
        classification_term="No Known Disease Relationship",
        submitted_on=on,
    )


@pytest.mark.parametrize(
    ("gene", "line"),
    [
        pytest.param(
            _agreed(
                _gencc("G2P", Classification.DEFINITIVE),
                _gencc("Labcorp Genetics (formerly Invitae)", Classification.STRONG),
                _no_association(
                    "Illumina",
                    "MONDO:0013463",
                    "congenital heart defects, multiple types, 6",
                    "2021-12-10 00:00:00",
                ),
            ),
            "Illumina reports no known disease relationship for congenital heart defects, "
            "multiple types, 6 (2021); see the validity table.",
            id="GDF1: one submitter, one disease",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, "MONDO:0007732", "Holt-Oram syndrome"),
                # The panel's own record, and ClinGen's GenCC copy of it, which
                # folds in the table and so must not be named twice here.
                replace(
                    _clingen(
                        Classification.NO_KNOWN_ASSOCIATION,
                        "MONDO:0005453",
                        "congenital heart disease",
                        "UD",
                    ),
                    classification_term="No Known Disease Relationship",
                    classification_date="2024-10-01T16:00:00.000Z",
                ),
                replace(
                    _gencc(
                        "ClinGen",
                        Classification.NO_KNOWN_ASSOCIATION,
                        "UD",
                        disease="MONDO:0005453",
                        label="congenital heart disease",
                    ),
                    classification_term="No Known Disease Relationship",
                    submitted_on="2024-10-01 00:00:00",
                ),
                # Two submitters in reverse alphabetical order, the first with
                # two diseases in reverse order and one record with no date.
                _no_association("Illumina", "MONDO:0013463", "tetralogy of Fallot", None),
                _no_association(
                    "Illumina", "MONDO:0007173", "atrial septal defect 7", "2019-03-01 00:00:00"
                ),
                _no_association(
                    "Ambry Genetics", "MONDO:0007173", "atrial septal defect 7", "2020-08-20"
                ),
            ),
            "The Congenital Heart Disease Gene Curation Expert Panel reports no known disease "
            "relationship for congenital heart disease (2024); Ambry Genetics reports no known "
            "disease relationship for atrial septal defect 7 (2020); Illumina reports no known "
            "disease relationship for atrial septal defect 7 (2019) and tetralogy of Fallot; see "
            "the validity table.",
            id="a panel and two submitters, sorted, the panel's GenCC copy folded",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, "MONDO:0007732", "Holt-Oram syndrome"),
            ),
            None,
            id="no report",
        ),
    ],
)
def test_the_card_names_who_reports_no_known_disease_relationship_and_for_what(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], line: str | None
) -> None:
    """`has_no_association_report` on the card, in words (owner decision 2026-10-26, R6).

    The bundle has published the flag and `no_association_reported_by` since
    schema 2.10, and the validity table lists the record, but the card said
    nothing: GDF1's names G2P among the submitters that admit it, with no sign
    that Illumina reports no known disease relationship for one of the same
    diseases. The
    line names each authority -- the panel for a ClinGen record, the submitter
    for a GenCC one -- its diseases and the year of its record, panels first,
    then submitters by name, diseases sorted, from the records rather than
    from the bundle's name list. The second fixture lists its submitters and
    one submitter's diseases in reverse order and carries no `set`, so a
    dropped sort fails on every hash seed; its ClinGen copy of the panel's
    record folds as it does in the table. Directly after the gloss.
    """
    fact, validity = gene
    reporters = tuple(
        sorted(
            {
                record.submitter or record.gcep or ""
                for record in validity.records
                if record.classification is Classification.NO_KNOWN_ASSOCIATION
            }
        )
    )
    rail = _card_page(
        tmp_path,
        replace(
            fact,
            has_no_association_report=bool(reporters),
            no_association_reported_by=reporters,
        ),
        validity,
    )
    if line is None:
        assert "card-no-association" not in rail
        assert "no known disease relationship" not in rail
        return
    assert rail.count('<p class="card-no-association">') == 1
    after = rail[rail.index("</dl>") + len("</dl>") :]
    if 'class="card-gloss"' in after:
        after = after[after.index("</p>") + len("</p>") :]
    assert after.startswith(f'<p class="card-no-association">{line}</p>')


def test_a_no_association_flag_with_no_record_to_name_is_refused(tmp_path: Path) -> None:
    """The flag and the records come from one derivation; disagreeing, they are a bug."""
    fact, validity = _graded(
        Classification.DEFINITIVE,
        _clingen(Classification.DEFINITIVE, "MONDO:0007732", "Holt-Oram syndrome"),
    )
    flagged = replace(
        fact, has_no_association_report=True, no_association_reported_by=("Illumina",)
    )
    with pytest.raises(ValueError, match="no such record"):
        _card_page(tmp_path, flagged, validity)


_NOONAN: Final = "MONDO:0018997"
_NSML: Final = "MONDO:0007893"
_HOLT_ORAM: Final = "MONDO:0007732"
_AOS4: Final = "MONDO:0014124"
_ORPHANET_WARRANTS: Final = {
    _NOONAN: Warrant("ORPHA:648", "Noonan syndrome", "direct", None, (_PULMONIC,)),
    _NSML: Warrant("ORPHA:500", "Noonan syndrome with multiple lentigines", "direct", None, ()),
    _AOS4: Warrant("ORPHA:974", "Adams-Oliver syndrome", "subtype", "OMIM:615297", ()),
}
_IN_SCOPE_BY_ORPHANET: Final = (
    "in scope because Orphanet records a structural heart malformation in at least 30% of cases"
)
_HETEROTAXY: Final = "MONDO:1060197"
_TARP: Final = "MONDO:0010711"
_CHD: Final = "MONDO:0005453"
_SIBLING: Final = "MONDO:0009341"
_ELLIS_VAN_CREVELD: Final = "MONDO:0009162"
_SYNDROMIC_PANEL: Final = "Syndromic Disorders Gene Curation Expert Panel"
_CILIOPATHY_PANEL: Final = "Kidney Cystic and Ciliopathy Disorders Gene Curation Expert Panel"
_SKELETAL_PANEL: Final = "Skeletal Disorders Gene Curation Expert Panel"


def _term(
    mondo: str, label: str, admitted_by: str, attributed_to: str, *, syndrome: bool
) -> ScopeEntry:
    return ScopeEntry.model_validate(
        {
            "id": mondo,
            "label": label,
            "reason": "A fixture term.",
            "admitted_by": admitted_by,
            "attributed_to": attributed_to,
            "syndrome": syndrome,
            "admitted_on": "2026-10-26",
        }
    )


# R5's registry: two marked syndromes two ClinGen panels scoped, one the CHD
# panel's; the CHD-itself heterotaxy term, unmarked; the root CHD term; and a
# marked syndrome a GenCC submitter scoped, which names no ClinGen panel.
_SCOPE_TERMS: Final = {
    _HOLT_ORAM: _term(
        _HOLT_ORAM,
        "Holt-Oram syndrome",
        "clingen_expert_panel",
        "Syndromic Disorders Gene Curation Expert Panel",
        syndrome=True,
    ),
    _TARP: _term(
        _TARP,
        "TARP syndrome",
        "clingen_chd_panel",
        "Congenital Heart Disease Gene Curation Expert Panel",
        syndrome=True,
    ),
    _HETEROTAXY: _term(
        _HETEROTAXY,
        "congenital heart disease with heterotaxy syndrome",
        "clingen_chd_panel",
        "Congenital Heart Disease Gene Curation Expert Panel",
        syndrome=False,
    ),
    _CHD: _term(
        _CHD,
        "congenital heart disease",
        "clingen_chd_panel",
        "Congenital Heart Disease Gene Curation Expert Panel",
        syndrome=False,
    ),
    _SIBLING: _term(
        _SIBLING, "Mowat-Wilson syndrome", "gencc_submitter", "Ambry Genetics", syndrome=True
    ),
    # Scoped on the ciliopathy panel's curation, which graded EVC; EVC2's grade
    # for it is the skeletal panel's (the real mirror, 2026-10-27).
    _ELLIS_VAN_CREVELD: _term(
        _ELLIS_VAN_CREVELD,
        "Ellis-van Creveld syndrome",
        "clingen_expert_panel",
        _CILIOPATHY_PANEL,
        syndrome=True,
    ),
}
_ON_SYNDROMIC_CURATION: Final = (
    "in scope on its curation by ClinGen&#x27;s Syndromic Disorders Gene Curation Expert Panel"
)
# EVC: the gene the ciliopathy panel graded for Ellis-van Creveld syndrome.
# `_card_page` puts it in the validity mapping every card is rendered against.
_EVC_GENE: Final = "HGNC:3497"
_EVC_VALIDITY: Final = GeneValidity(
    records=(
        _clingen(
            Classification.DEFINITIVE,
            _ELLIS_VAN_CREVELD,
            "Ellis-van Creveld syndrome",
            "AR",
            gcep=_CILIOPATHY_PANEL,
        ),
    ),
    state=ValidityState.EXPERT_CURATED,
    has_source_discordance=False,
)


@pytest.mark.parametrize(
    ("gene", "line"),
    [
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _NOONAN, "Noonan syndrome"),
            ),
            f"Graded for Noonan syndrome as a whole syndrome; {_IN_SCOPE_BY_ORPHANET}.",
            id="Noonan",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE, _NSML, "Noonan syndrome with multiple lentigines"
                ),
                _clingen(Classification.DEFINITIVE, _NOONAN, "Noonan syndrome"),
            ),
            "Graded for Noonan syndrome and Noonan syndrome with multiple lentigines as "
            f"whole syndromes; {_IN_SCOPE_BY_ORPHANET}.",
            id="PTPN11: two Orphanet-scoped syndromes at the headline, both named",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _AOS4, "Adams-Oliver syndrome 4"),
            ),
            f"Graded for Adams-Oliver syndrome 4 as a whole syndrome; {_IN_SCOPE_BY_ORPHANET} of "
            "Adams-Oliver syndrome.",
            id="EOGT: a subtype is in scope on the syndrome it belongs to",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _NOONAN, "Noonan syndrome"),
                _clingen(Classification.DEFINITIVE, _AOS4, "Adams-Oliver syndrome 4"),
            ),
            "Graded for Adams-Oliver syndrome 4 and Noonan syndrome as whole syndromes; "
            f"{_IN_SCOPE_BY_ORPHANET} of Adams-Oliver syndrome and Noonan syndrome.",
            id="one direct and one subtype: the subtype's syndrome is still named",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _HOLT_ORAM, "Holt-Oram syndrome"),
                _clingen(Classification.DEFINITIVE, _NOONAN, "Noonan syndrome"),
            ),
            None,
            id="also graded for a disease scoped as CHD: the line would be false",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _CHD, "congenital heart disease"),
                _clingen(Classification.LIMITED, _NOONAN, "Noonan syndrome"),
            ),
            None,
            id="Orphanet-scoped only below the headline",
        ),
        pytest.param(
            _agreed(
                _gencc("Ambry Genetics", Classification.LIMITED, disease=_NOONAN, label="Noonan"),
                _gencc("G2P", Classification.DEFINITIVE, disease=_NOONAN, label="Noonan"),
            ),
            None,
            id="no panel graded it, so nothing was graded for anything",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE,
                    _HOLT_ORAM,
                    "Holt-Oram syndrome",
                    gcep=_SYNDROMIC_PANEL,
                ),
            ),
            f"Graded for Holt-Oram syndrome as a whole syndrome; {_ON_SYNDROMIC_CURATION}.",
            id="TBX5: a marked syndrome a ClinGen panel scoped",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE,
                    _ELLIS_VAN_CREVELD,
                    "Ellis-van Creveld syndrome",
                    "AR",
                    gcep=_SKELETAL_PANEL,
                ),
            ),
            "Graded for Ellis-van Creveld syndrome as a whole syndrome; in scope on "
            "ClinGen&#x27;s Kidney Cystic and Ciliopathy Disorders Gene Curation Expert "
            "Panel&#x27;s curation of it for EVC.",
            id="EVC2: the panel scoping the term graded another gene for it",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE,
                    _HOLT_ORAM,
                    "Holt-Oram syndrome",
                    gcep=_SYNDROMIC_PANEL,
                ),
                _clingen(Classification.DEFINITIVE, _TARP, "TARP syndrome"),
            ),
            "Graded for Holt-Oram syndrome and TARP syndrome as whole syndromes; in scope on "
            "ClinGen&#x27;s Syndromic Disorders Gene Curation Expert Panel&#x27;s curation of "
            "Holt-Oram syndrome and ClinGen&#x27;s Congenital Heart Disease Gene Curation "
            "Expert Panel&#x27;s curation of TARP syndrome.",
            id="two marked syndromes, each named with its own panel",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE,
                    _HOLT_ORAM,
                    "Holt-Oram syndrome",
                    gcep=_SYNDROMIC_PANEL,
                ),
                _clingen(
                    Classification.DEFINITIVE,
                    _ELLIS_VAN_CREVELD,
                    "Ellis-van Creveld syndrome",
                    "AR",
                    gcep=_SKELETAL_PANEL,
                ),
            ),
            "Graded for Ellis-van Creveld syndrome and Holt-Oram syndrome as whole syndromes; "
            "in scope on ClinGen&#x27;s Kidney Cystic and Ciliopathy Disorders Gene Curation "
            "Expert Panel&#x27;s curation of Ellis-van Creveld syndrome for EVC and "
            "ClinGen&#x27;s Syndromic Disorders Gene Curation Expert Panel&#x27;s curation of "
            "Holt-Oram syndrome.",
            id="two marked syndromes, one scoped on another gene's grade",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(
                    Classification.DEFINITIVE,
                    _HETEROTAXY,
                    "congenital heart disease with heterotaxy syndrome",
                ),
            ),
            None,
            id="ZIC3: CHD itself, unmarked, though its name says syndrome",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _HOLT_ORAM, "Holt-Oram syndrome"),
                _clingen(Classification.DEFINITIVE, _CHD, "congenital heart disease"),
            ),
            None,
            id="a marked syndrome beside an unmarked CHD term at the headline",
        ),
        pytest.param(
            _graded(
                Classification.DEFINITIVE,
                _clingen(Classification.DEFINITIVE, _SIBLING, "Mowat-Wilson syndrome"),
            ),
            None,
            id="a marked syndrome no ClinGen panel scoped names no panel",
        ),
    ],
)
def test_the_card_says_when_the_grade_is_for_a_syndrome_rather_than_for_heart_disease(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], line: str | None
) -> None:
    """PTPN11 is Definitive for Noonan syndrome, in a CHD atlas (D55).

    Since 2026-10-25 a syndrome can be in scope because Orphanet records a
    structural heart malformation in at least 30% of its cases. The validity
    section's notes say so at length; the card says it in one line beside the
    chip, because the chip is what a reader reads first -- and says *why* it
    is in scope, which "Graded for Noonan syndrome, not for heart disease"
    (the wording until 2026-10-26) left a reader to guess. Several diseases
    join with "and", never "; ", which the warrant row above already uses
    between diseases.

    **A subtype's warrant is its syndrome's.** Orphanet records the heart
    malformation for Adams-Oliver syndrome, and Adams-Oliver syndrome 4 (EOGT's
    headline disease, the one subtype among the 21 Orphanet lines measured
    2026-10-26) is in scope as its subtype; "in at least 30% of cases" alone
    would read as a frequency recorded for the subtype itself. **And still
    named beside a direct one**: a headline graded for Noonan syndrome and for
    Adams-Oliver syndrome 4 names both syndromes. Attributing only when every
    warrant is a subtype's would drop "of Adams-Oliver syndrome" there, and
    that mutant survived the release review's 634-test mutation run until this
    case was added (2026-10-26; latent, none of the 21 lines on a real build
    that day).

    **Rendered only when every disease the headline grade is for is
    Orphanet-scoped**, and naming all of them. A gene graded at its headline
    for a syndrome *and* for a disease an authority curated as CHD would read
    "not for heart disease" while it is -- latent today (0 of 98 graded genes,
    measured 2026-10-26), and pinned so it stays a sentence nobody can publish.
    PTPN11 is graded Definitive for two Orphanet-scoped syndromes, and a line
    naming one would state half of what the warrant row above it names.

    **A syndrome a ClinGen panel scoped gets the line too** (owner decision
    2026-10-26, R5), read off the curated `syndrome` mark: Holt-Oram is marked,
    and "congenital heart disease with heterotaxy syndrome" (ZIC3, NODAL) is not
    -- it is CHD itself, whatever its name says. Every headline disease must be
    marked and ClinGen-scoped, as every one must be Orphanet-scoped above: a
    marked syndrome beside an unmarked CHD term gets nothing, and a marked
    term a GenCC submitter scoped has no panel to name. It says the syndrome is
    in scope *on* the panel's curation, not *because* the panel curates it: the
    Syndromic Disorders panel curates 176 disease terms, 3 of them in scope.
    """
    fact, validity = gene
    rail = _card_page(
        tmp_path, fact, validity, warrants=_ORPHANET_WARRANTS, scope_terms=_SCOPE_TERMS
    )
    if line is None:
        assert "card-scope" not in rail
        assert "in scope because Orphanet" not in rail
        assert "in scope on" not in rail
        assert "not for heart disease" not in rail
    else:
        assert f'<p class="card-scope">{line}</p>' in rail


def test_a_syndrome_scoped_on_a_panel_that_graded_no_gene_for_it_is_refused(
    tmp_path: Path,
) -> None:
    """Names no gene rather than an empty one, and never "for another gene".

    A term attributed to a panel the mirror holds no record of for it is what
    SCP005 refuses, so the build cannot reach this; if it does, the gate was
    bypassed, and a card naming nobody, or a panel that graded nothing, is a
    provenance claim nothing supports.
    """
    fact, validity = _graded(
        Classification.DEFINITIVE,
        _clingen(Classification.DEFINITIVE, _TARP, "TARP syndrome", gcep=_SKELETAL_PANEL),
    )
    with pytest.raises(ValueError, match="graded no gene for it"):
        _card_page(tmp_path, fact, validity, scope_terms=_SCOPE_TERMS)


_NOT_CURATED_BY_PANEL: Final = (
    "Not yet curated by the atlas: every classification on this page is a ClinGen expert "
    "panel's or a named curator's."
)


@pytest.mark.parametrize(
    ("gene", "sentence"),
    [
        pytest.param(
            (
                _facts(GATA4, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,)),
                _validity(),
            ),
            "Curated by the atlas: 1 lesion assertion for the septal lesion group.",
            id="curated",
        ),
        pytest.param(
            (
                _facts(
                    GATA4,
                    AtlasCuration.CURATED,
                    groups=(LesionGroup.CONOTRUNCAL, LesionGroup.SEPTAL),
                    functional_count=2,
                ),
                _validity(),
            ),
            "Curated by the atlas: 1 lesion assertion for the conotruncal, septal lesion "
            "groups; 2 functional records.",
            id="curated, with functional records",
        ),
        pytest.param(
            (
                _facts(
                    GATA4, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,), functional_count=1
                ),
                _validity(),
            ),
            "Curated by the atlas: 1 lesion assertion for the septal lesion group; 1 functional "
            "record.",
            id="curated, with one functional record",
        ),
        pytest.param(
            (_facts(GATA4, AtlasCuration.NOT_YET_CURATED), _validity()),
            _NOT_CURATED_BY_PANEL,
            id="not curated, a panel graded it",
        ),
        pytest.param(
            (
                replace(
                    _agreed()[0],
                    validity_state=ValidityState.EXPERT_CURATED,
                ),
                GeneValidity(
                    records=(
                        _clingen(
                            Classification.NO_KNOWN_ASSOCIATION,
                            "MONDO:0007173",
                            "atrial septal defect 7",
                        ),
                        _gencc("Ambry Genetics", Classification.LIMITED),
                        _gencc("G2P", Classification.DEFINITIVE),
                    ),
                    state=ValidityState.EXPERT_CURATED,
                    has_source_discordance=False,
                ),
            ),
            _NOT_CURATED_BY_PANEL,
            id="not curated, admitted on agreement over a panel's no known association",
        ),
        pytest.param(
            _agreed(
                _gencc("Ambry Genetics", Classification.LIMITED),
                _gencc("G2P", Classification.DEFINITIVE),
            ),
            "Not yet curated by the atlas, and no ClinGen expert panel has graded this gene for "
            "any of the diseases in scope here: it is published because the Gene Curation "
            "Coalition submitters it is admitted by independently assert it.",
            id="not curated, admitted on agreement",
        ),
        pytest.param(
            (_facts(GATA4, AtlasCuration.NOT_YET_CURATED, functional_count=1), _validity()),
            "No lesion assertion yet curated by the atlas: every classification on this page is "
            "a ClinGen expert panel's or a named curator's. The atlas has curated 1 functional "
            "record, published in this gene's JSON.",
            id="not curated, with a functional record",
        ),
    ],
)
def test_the_card_says_in_plain_words_what_the_atlas_itself_curated(
    tmp_path: Path, gene: tuple[GeneFacts, GeneValidity], sentence: str
) -> None:
    """D55: one sentence, on the card, in place of a token and a main-column notice.

    Four branches, each pinned against its literal (CLAUDE.md section 4.38):

    * curated -- what the atlas curated, with the counts only where non-zero;
    * not curated, panel-graded -- every classification shown is someone
      else's;
    * not curated, admitted on agreement -- no panel graded it, and the
      submitters named in the card's own `admitted by` row are why it is here.
      **Scoped to diseases an external authority treats as CHD**, because
      unscoped it is false: ClinGen grades ACTC1 Definitive for hypertrophic
      cardiomyopathy and ELN Definitive for cutis laxa, both admitted here on
      agreement (measured 2026-10-26 against the ClinGen mirror). Its scope
      clause is the grade key's `ungraded` gloss's and the gene index's,
      word for word, so the three cannot describe the gate differently.
      **Scoped to the list, not to a property of a disease**: it read "for a
      disease an external authority treats as congenital heart disease, or
      in which Orphanet records a structural heart malformation in at least
      30% of cases" until 2026-10-26, false on FLNA's card -- ClinGen grades
      FLNA Definitive for periventricular nodular heterotopia, which would
      qualify through Orphanet and which no curator has nominated. That
      clause is asserted absent from the whole page, so it cannot return on
      another line of the card or in the validity section either;
    * not curated with functional records (latent: 0 genes today) -- the
      sentence names the lesion assertion as what is absent, because "not
      curated by the atlas" beside "the atlas has curated 1 functional record"
      denies the atlas's own work in the sentence that counts it.

    "Not *yet* curated", because the state is `not_yet_curated`: the atlas has
    not decided against curating these genes.

    The discriminator between the two uncurated sentences is `validity_state`,
    not the headline: a gene admitted on agreement that ClinGen recorded `No
    Known Disease Relationship` for in scope has a panel record, and the first
    sentence stays true of it where the second would not. That gene is a case
    here, because it is the only one that tells the two apart: every other
    fixture pairs a panel record with a headline, and keying the sentence on
    the headline survived all 1,297 tests until it was added (2026-10-26). Its
    page carries "no ClinGen expert panel has graded" nowhere.

    Exactly one curation paragraph per card, and none in the main column: the
    sentence moved here, and a page stating it twice is a page whose two copies
    can come to disagree.

    **The card ends on two links, after this sentence (P1, 2026-10-26):** "How
    to read this page", to the top of `methods.html`, and the gene's own JSON.
    The nav links the same file as "How to read this atlas"; the card's text
    differs on purpose, because from the card the reader's question is about
    the page in front of them. Every section a gene page renders has a section
    there, except the atlas's own curated evidence and its publications, which
    render on 1 of 119 pages (measured 2026-10-26). The link is to the top, not
    an anchor: the card abbreviates the scope, validity and rare-variant
    sections at once.
    """
    fact, validity = gene
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {GATA4: fact},
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        cell_types={},
    )
    page = _page(tmp_path, "HGNC_4173.html")
    rail = _rail_of(page)

    assert re.findall(r'<p class="card-curation">(.*?)</p>', rail) == [sentence]
    assert rail.endswith(
        f'<p class="card-curation">{sentence}</p>'
        '<p><a href="../methods.html">How to read this page</a> &middot; '
        '<a href="../genes/HGNC_4173.json">this gene as JSON</a></p>'
    )
    # Only where the card's own sentence says it, never on another line.
    denial = "no ClinGen expert panel has graded"
    assert (denial in page) is (denial in sentence)
    assert "in which Orphanet records" not in page
    main = page[page.index("</aside>") :]
    assert "curated by the atlas" not in main.lower()
    assert "not yet curated</strong> a lesion assertion" not in page


def test_every_value_the_card_interpolates_is_escaped(tmp_path: Path) -> None:
    """Disease labels and author names are mirrored or curated text, not page copy.

    Each card helper escapes what it interpolates, and each escape had no test
    until review measured that dropping it survived the suite: the Orphanet
    line's disease label, the dispute line's disease label, and the study
    label on a rare-variant line, its link's accessible name and the
    no-correction line. The no-association line's submitter and disease label
    joined them after the next review (2026-10-27): its escape survived all
    1,383 tests.

    **And the fixture carries all three lines a disputed, Orphanet-scoped gene
    can show, so their order is pinned here too**: the dispute line, then the
    no-association line, then the Orphanet line (`_rail`'s order, and
    `docs/data-api.md`'s). The R6 test pins the no-association line directly
    after the gloss, which only holds where no dispute line sits between.
    """
    marked = "<i>Noonan</i>"
    fact, validity = _graded(
        Classification.DEFINITIVE,
        _rasopathy(Classification.DEFINITIVE, _NOONAN, marked),
        _rasopathy(Classification.DISPUTED, "MONDO:0009026", "<b>Costello</b>"),
        _no_association("<q>Lab</q>", "MONDO:0007173", "<kbd>ASD7</kbd>", "2021-12-10 00:00:00"),
    )
    publications = {
        **_ALL_STUDIES,
        "PMID:40127276": _SIERANT_PUBLICATION["PMID:40127276"].model_copy(
            update={"authors": ["<s>Sierant</s>"]}
        ),
        "PMID:42230622": _PUBLICATION.model_copy(update={"authors": ["<u>Audain</u>"]}),
    }
    rows = [_sierant(), _burden_row(cohort_stratum="all", effect=3.1, effect_bound=None)]
    flagged = replace(
        fact,
        has_conflicting_evidence=True,
        has_no_association_report=True,
        no_association_reported_by=("<q>Lab</q>",),
    )
    build_gene_pages(
        {GATA4: flagged},
        Emitter(root=tmp_path),
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        assertions={},
        publications=publications,
        burden={GATA4: rows},
        cohorts=_COHORTS,
        families=cohort_families(rows, _COHORTS),
        axes=evidence_axes(rows),
        scope_warrants=_ORPHANET_WARRANTS,
        cell_types={},
    )
    rail = _rail_of(_page(tmp_path, "HGNC_4173.html"))

    scope = _slice_between(rail, '<p class="card-scope">', ("</p>",))
    assert scope.startswith("Graded for &lt;i&gt;Noonan&lt;/i&gt; as a whole syndrome;")
    dispute = _slice_between(rail, '<p class="card-dispute">', ("</p>",))
    assert dispute.startswith("Disputed for &lt;b&gt;Costello&lt;/b&gt; by the RASopathy")
    line = _slice_between(_card_variants(rail), "<li>", ("</li>",))
    assert line.startswith("&lt;s&gt;Sierant&lt;/s&gt; et al. 2025 · ")
    assert 'aria-label="Details for &lt;s&gt;Sierant&lt;/s&gt; et al. 2025"' in line
    assert '">&lt;u&gt;Audain&lt;/u&gt; et al. 2026</a> publishes no corrected p' in rail
    silence = _slice_between(rail, '<p class="card-no-association">', ("</p>",))
    assert silence == (
        "&lt;q&gt;Lab&lt;/q&gt; reports no known disease relationship for "
        "&lt;kbd&gt;ASD7&lt;/kbd&gt; (2021); see the validity table."
    )
    for raw in ("<i>", "<b>", "<s>", "<u>", "<q>", "<kbd>"):
        assert raw not in rail, raw
    lines = [
        rail.index(f'<p class="card-{kind}">') for kind in ("dispute", "no-association", "scope")
    ]
    assert lines == sorted(lines), "the dispute, no-association and Orphanet lines moved"


def test_two_studies_that_would_share_one_anchor_are_refused(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`PMID:123` and `123` both slug to `study-123`; one card link would land on
    the wrong study's caveats. Distinct PMIDs cannot collide, so this is a guard
    on a bypassed gate, and it raises rather than publishing."""
    rows = [_burden_row(study="PMID:123"), _burden_row(study="123")]
    with pytest.raises(ValueError, match="share a page anchor"):
        _burden_page(tmp_path, facts_uncurated, rows)


def test_a_study_no_cohort_family_holds_is_refused_rather_than_miscounted(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """ "Tested in N independent datasets" counts corpus families; a study absent
    from all of them means the families were derived over another corpus."""
    with pytest.raises(ValueError, match="no cohort family holds PMID:42230622"):
        _burden_page(tmp_path, facts_uncurated, [_burden_row()], families=(frozenset({"PMID:1"}),))


@pytest.mark.parametrize(
    ("carriers", "total", "unit", "shown"),
    [
        (1, 14040, "individuals", "1 carrier / 14,040"),
        (1, 20151, "alleles", "1 allele / 20,151"),
        (1, 3887, "de_novo_mutations", "1 de novo / 3,887 trios"),
        (0, 45082, "individuals", "0 carriers / 45,082"),
        (2, 237963, "alleles", "2 alleles / 237,963"),
    ],
)
def test_a_count_of_one_names_one_carrier(carriers: int, total: int, unit: str, shown: str) -> None:
    """ "1 carriers / 14,040" rendered in tables and on the card until 2026-10-26."""
    assert _count(carriers, total, unit) == shown


_COMBINED_BANNER = (
    '<p class="notice-inline">Cell-type values pool 8 hearts at 9 to 15 weeks and the '
    "timeline below is from bulk tissue; a low or absent value is not evidence a gene is "
    "unimportant to heart development.</p>"
)


def _bulk_profile() -> ExpressionProfile:
    return _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(tissues=(_tissue_entry(placement=_placement()),)),)
            ),
        )
    )


def test_a_gene_with_cell_type_data_carries_the_combined_banner_and_not_the_bulk_only_one(
    tmp_path: Path,
) -> None:
    page = _expression_page(
        tmp_path,
        {GATA4: _bulk_profile()},
        {"E-MTAB-6814": _profile_dataset(), "hoc": cell_type_dataset()},
        cell_types={GATA4: CellTypeExpression(datasets=[cell_type_payload()])},
    )
    section = _expression_section_text(page)

    assert _COMBINED_BANNER in section
    assert "Every measurement below is from bulk tissue" not in section
    # The cell-type block opens the section, before the first bulk dataset block.
    assert section.index("<h3>Cell types (Farah 2024)</h3>") < section.index("E-MTAB-6814")
    assert section.index(_COMBINED_BANNER) < section.index("<h3>Cell types")
    # The figures are named by the page's symbol, not by its id or the dataset's.
    assert 'aria-label="GATA4 expression by cell type"' in section


def test_a_gene_without_cell_type_data_keeps_the_bulk_banner_and_no_block(
    tmp_path: Path,
) -> None:
    # TBX5 has no entry in `cell_types`; GATA4 beside it does.
    page = _expression_page(
        tmp_path,
        {GATA4: _bulk_profile(), TBX5: _bulk_profile()},
        {"E-MTAB-6814": _profile_dataset(), "hoc": cell_type_dataset()},
        facts={
            GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
            TBX5: _facts(TBX5, AtlasCuration.NOT_YET_CURATED),
        },
        name="HGNC_11604.html",
        cell_types={GATA4: CellTypeExpression(datasets=[cell_type_payload()])},
    )
    section = _expression_section_text(page)

    assert "Every measurement below is from bulk tissue" in section
    assert "Cell-type values pool" not in section
    assert "Cell types (Farah" not in section
    assert '<svg class="ct-svg"' not in section


def test_a_gene_with_cell_type_data_and_no_bulk_profile_claims_no_timeline(
    tmp_path: Path,
) -> None:
    page = _expression_page(
        tmp_path,
        {GATA4: EMPTY_EXPRESSION_PROFILE},
        {"hoc": cell_type_dataset()},
        cell_types={GATA4: CellTypeExpression(datasets=[cell_type_payload()])},
    )
    section = _expression_section_text(page)

    assert (
        '<p class="notice-inline">Cell-type values pool 8 hearts at 9 to 15 weeks; a low or '
        "absent value is not evidence a gene is unimportant to heart development.</p>"
    ) in section
    assert "timeline" not in section
    assert "not yet curated" not in section
    assert "<h3>Cell types (Farah 2024)</h3>" in section


def test_the_cell_type_block_adds_at_most_41_visible_words_to_a_gene_page(
    tmp_path: Path,
) -> None:
    """Measured with the budget test's own counting (`_prose`, `_words`) on a fixture
    page (with its tau line) with and without the block. The spec's ceiling is 40;
    the measured growth was 40 (block 35 + banner 5, re-measured 2026-10-28 after
    the tau line moved below the chart) and is 41 since the owner's one-sentence
    banner of 2026-11-02 (the banner gained "and", which joins the two clauses), so the ceiling
    is raised by that word and any further visible word the block adds fails here.
    Words inside an SVG or a closed fold do not count, which is why the region grid's
    key costs none."""
    from tests.test_built_site_is_consumable import _prose, _words

    datasets = {"E-MTAB-6814": _profile_dataset(), "hoc": cell_type_dataset()}
    profiles = {GATA4: _bulk_profile()}
    plain = _expression_page(tmp_path / "a", profiles, datasets)
    with_cells = _expression_page(
        tmp_path / "b",
        profiles,
        datasets,
        cell_types={
            GATA4: CellTypeExpression(datasets=[cell_type_payload(spec=cell_type_specificity())])
        },
    )

    def words(page: str) -> int:
        return sum(_words(text) for text in _prose(page))

    assert words(with_cells) - words(plain) <= 41
