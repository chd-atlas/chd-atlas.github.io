# tests/unit/test_build_landing.py
"""`index.html`, the site's front door: the research-use statement and its counts."""

from __future__ import annotations

import html
import json
import re
import shutil
from datetime import date
from pathlib import Path

import pytest

from chd_atlas.build.burden import BurdenCensus, BurdenRow, burden_census
from chd_atlas.build.cell_types import CellTypeCensus
from chd_atlas.build.concordance import survives_own_correction
from chd_atlas.build.emit import Emitter
from chd_atlas.build.landing import build_landing
from chd_atlas.build.paths import LANDING
from chd_atlas.build.profiles import ProfileCensus
from chd_atlas.build.render import RESEARCH_USE_NOTICE, STYLESHEET
from chd_atlas.build.runner import build_site
from chd_atlas.build.validity import GeneValidity, ValidityRecord, uncurated
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator
from chd_atlas.models.literature import Publication
from chd_atlas.validate.sources import Source, SourceRegistry
from chd_atlas.vocab import Classification, ValiditySource, ValidityState


# What a corpus with no burden evidence at all publishes. Every test in this
# file that is about something else -- escaping, pluralisation, the shell --
# renders against it.
def _orphadata(
    source_id: str, file_version: str | None = None, file_date: str | None = None
) -> Source:
    """An Orphadata product's registry entry; only the header fields are read here."""
    return Source.model_validate(
        {
            "id": source_id,
            "name": "Orphadata Science",
            "version": "sha256:" + "0" * 64,
            "retrieved_on": "2026-10-25",
            "url": "https://www.orphadata.com/",
            "licence": "CC-BY-4.0",
            "redistribution": "permitted_with_attribution",
            "file_version": file_version,
            "file_date": file_date,
        }
    )


# Both products of one release, as `mirrors/sources.yaml` records them today:
# one version, one day, two times of day.
_SOURCES = SourceRegistry(
    sources=[
        _orphadata("orphanet_product4", "1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:57:18"),
        _orphadata("orphanet_product1", "1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:53:50"),
    ]
)

_NO_BURDEN = BurdenCensus(rows=0, genes=0, families=0, studies=(), surviving_genes=0)

# The developmental-expression counterpart to `_NO_BURDEN`, for the same tests.
_NO_PROFILES = ProfileCensus(genes=0, datasets=0)
_NO_CELL_TYPES = CellTypeCensus(genes=0, datasets=0)

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"
# Published, and carrying 9 rows in the committed mirror. The fixture below
# removes them to build a corpus where the published population is wider than
# the burden evidence.
TBX20 = "HGNC:11598"

# tests/unit/test_build_landing.py -> tests/unit -> tests -> repo root, the same
# climb `tests/unit/test_build_runner.py` uses to find the committed corpus.
REPO = Path(__file__).parent.parent.parent


def _burden_row(**overrides: object) -> BurdenRow:
    """A minimal published burden row. Only `study` is ever varied here.

    The census counts rows and reads `study` for nothing at all -- families are
    passed in -- so the statistical columns are inert in this file. They are
    filled with a real TAB2 row's values rather than zeros so that a future test
    reaching for one of them gets a plausible number instead of a fiction.
    """
    payload: dict[str, object] = {
        "study": "PMID:1",
        "gene": "HGNC:17075",
        "cohort_stratum": "all",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "maf_max": 0.001,
        "count_unit": "individuals",
        "n_case_carriers": 6,
        "n_cases": 3876,
        "comparator": "control_cohort",
        "n_control_carriers": 0,
        "n_controls": 45082,
        "expected_count": None,
        "effect": None,
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above",
        "ci_low": 13.7,
        "ci_high": None,
        "pvalue": 2.45e-07,
        "pvalue_test": "fisher_exact",
        "pvalue_adjusted": None,
        "pvalue_adjustment": None,
        "case_cohorts": ("cnchd",),
        "control_cohorts": ("ukbb",),
        "method_note": None,
        "source": "audain2026_sd3",
    }
    payload.update(overrides)
    return BurdenRow(**payload)  # type: ignore[arg-type]


def _evidence(**overrides: object) -> Evidence:
    payload: dict[str, object] = {
        "publication": "PMID:8988165",
        "evidence_class": "genetic_case",
        "locator": SupplementaryLocator(file="T.xlsx"),
        "summary": "s",
        "strength": "strong",
    }
    payload.update(overrides)
    return Evidence.model_validate(payload)


def _assertion(**overrides: object) -> LesionAssertion:
    payload: dict[str, object] = {
        "id": "CHDA:AST:0000001",
        "gene": TBX5,
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "both",
        "evidence": [_evidence()],
        "curator": "c",
        "curated_on": date(2026, 7, 1),
        "last_reviewed": date(2026, 7, 1),
    }
    payload.update(overrides)
    return LesionAssertion.model_validate(payload)


def _build(
    corpus: Corpus,
    symbols: dict[str, str],
    validity: dict[str, GeneValidity],
    tmp_path: Path,
    published: set[str] | None = None,
    census: BurdenCensus = _NO_BURDEN,
    profile_census: ProfileCensus = _NO_PROFILES,
    cell_type_census: CellTypeCensus = _NO_CELL_TYPES,
    sources: SourceRegistry = _SOURCES,
) -> str:
    """Render the page.

    `published` defaults to the genes the fixture asserts, which is what the
    page said before D21 and is *not* what it says now -- the default only
    keeps the tests that are about something else (escaping, pluralisation, the
    research-use notice) from having to name a population they do not care
    about. Every test that is about the published-gene figure passes it
    explicitly, and one of them passes a set that disagrees with the assertions
    on purpose.

    `census` defaults to an empty one for the same reason, and the tests about
    the burden figures build theirs from a real `burden_census` call over a
    fixture rather than handing this function three literals -- a page test that
    invented its own census could not see the derivation go wrong. `profile_census`
    defaults the same way, for the developmental-expression pair -- this
    default lives only here, never in `build_landing` itself, which requires
    the argument for the reason its own docstring gives.
    """
    emitter = Emitter(root=tmp_path)
    build_landing(
        corpus,
        symbols=symbols,
        validity=validity,
        published=(
            {assertion.gene for assertion in corpus.assertions} if published is None else published
        ),
        census=census,
        profile_census=profile_census,
        cell_type_census=cell_type_census,
        emitter=emitter,
        sources=sources,
    )
    return (tmp_path / LANDING).read_text(encoding="utf-8")


def _hero(text: str) -> str:
    """The hero band, whitespace-collapsed: from its `<div>` to the first `<h2>`."""
    return _prose(text[text.index('<div class="hero">') : text.index("<h2>")])


def _figure(value: int | str, label: str) -> str:
    """One hero figure as `_hero` renders it: the value immediately beside its label.

    Asserted with its own label, never as a bare number: three values in one
    band would otherwise let two of them swap unnoticed.
    """
    return f'"figure-value">{value}</span> <span class="figure-label">{label}<'


def _prose(text: str) -> str:
    """Whitespace-collapsed text, for asserting a sentence rather than a source line.

    The template wraps long sentences across several source lines the way every
    other line in this file does, and a browser collapses that whitespace when
    it renders — but `"in text"` does not, so a phrase split across a wrap would
    fail a substring check that has nothing to do with what a reader sees. Used
    for prose; the escaping and no-external-request tests below check the raw
    text instead, because a literal `<script` or an unescaped `<` is exactly
    what must not be normalised away.
    """
    return re.sub(r"\s+", " ", text)


def _section(text: str, heading: str) -> str:
    """The slice of the page from one `<h2>` heading up to (not including) the next.

    Lets a test pin *where* a row lives, not only whether the text appears
    somewhere on the page — the defect this file is guarding against was a true
    number in the wrong section, which a bare substring check cannot catch.
    """
    start = text.index(heading)
    nxt = text.find("<h2>", start + len(heading))
    return text[start : nxt if nxt != -1 else len(text)]


def _checkout(root: Path) -> Path:
    """The committed corpus, copied where a build can run over it.

    `curation` and `mirrors` are copied (cheap, and a caller may trim them);
    `ontologies` is symlinked, the same shortcut `test_build_runner.py` and
    `test_docs_match_the_build.py` use for the same 11 MB of pinned OBO releases.
    """
    root.mkdir(parents=True)
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    return root


@pytest.fixture(scope="module")
def build_with_an_untested_gene(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A build whose published population is wider than its burden evidence.

    **The committed corpus cannot tell two of the manifest's counts apart.** All
    23 published genes carry burden rows, so `len(published)` and "genes with at
    least one burden row" are both 23 and the mutant publishing the second under
    the name `genes` survives every test that reads a real build. Measured
    2026-08-05; the same degeneracy CLAUDE.md section 4.15b records for the
    one-gene corpus, in a new place.

    So this fixture drops one published gene's rows from the mirror and nothing
    else. TBX20 is chosen because it is published, carries rows in the committed
    mirror, and no other gene's row references it — the burden schema keys on
    (study, gene), so removing a gene's rows leaves every remaining row valid and
    every cohort, study and gene reference intact. `validate_repository` still
    passes, which is the point: this is a legitimate corpus, not a broken one,
    and it is what the atlas looks like the day a published gene has no burden
    study behind it.

    A second module-scoped build costs about a second and a half of validation.
    It buys the only distinction between two published figures that will diverge
    the first time D21 admits a gene no burden study covered.
    """
    base = tmp_path_factory.mktemp("landing-untested-gene")
    root = _checkout(base / "repo")
    mirror = root / "mirrors" / "burden.tsv"
    lines = mirror.read_text(encoding="utf-8").splitlines(keepends=True)
    header, rows = lines[0], lines[1:]
    gene = header.rstrip("\n").split("\t").index("gene")
    kept = [row for row in rows if row.split("\t")[gene] != TBX20]
    assert len(kept) < len(rows), "TBX20 has no rows in the mirror; pick another gene"
    mirror.write_text(header + "".join(kept), encoding="utf-8")

    out = base / "dist"
    build_site(root, out)
    return out


def test_the_manifest_gene_count_is_the_population_not_the_genes_carrying_burden(
    build_with_an_untested_gene: Path,
) -> None:
    """`counts.genes` sizes a consumer's fetch loop, so it must count bundles.

    One gene's burden rows are removed from the mirror and nothing else changes.
    The site still publishes every bundle — the publication gate is an external
    authority's classification, not the presence of burden evidence — so `genes`
    must stay at the bundle count while the genes carrying evidence falls. A
    consumer sizing a fetch loop on the smaller number silently skips a gene that
    has a bundle, an assertion and a page.

    Asserted as an inequality as well as an equality: `genes == len(bundles)`
    alone would pass on the corpus this atlas shipped on 2026-08-05, where the
    two numbers coincided, which is exactly how the mutant survived before this
    fixture existed.

    **The fixture's health check is TBX20's own bundle, not an arithmetic
    identity.** It was `carrying == len(bundles) - 1` until 2026-08-06, which
    silently assumed every other published gene carries burden rows. The widened
    gate admits genes no burden study covered — measured on the committed
    mirrors, one published gene has no rows before this fixture removes any — so
    that identity started failing for a reason that had nothing to do with what
    it was guarding. Naming the gene the fixture edits is the check that stays
    true however many others are untested, and it fails loudly if TBX20 ever
    stops being published.
    """
    manifest = json.loads((build_with_an_untested_gene / "manifest.json").read_text())
    counts = manifest["counts"]
    bundles = sorted((build_with_an_untested_gene / "genes").glob("HGNC_*.json"))
    carrying = sum(1 for path in bundles if json.loads(path.read_text(encoding="utf-8"))["burden"])

    emptied = build_with_an_untested_gene / "genes" / f"{TBX20.replace(':', '_')}.json"
    assert emptied.is_file(), "the fixture's gene is no longer published; pick another"
    assert not json.loads(emptied.read_text(encoding="utf-8"))["burden"], (
        "the fixture no longer removes that gene's rows"
    )

    assert counts["genes"] == len(bundles)
    assert counts["genes"] != carrying


@pytest.fixture(scope="module")
def real_build(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A real build of the committed corpus, for checking the page against `genes/index.json`.

    `curation` and `mirrors` are copied (cheap, and nothing here mutates them);
    `ontologies` is symlinked, the same shortcut `test_build_runner.py` and
    `test_docs_match_the_build.py` use for the same 11 MB of pinned OBO releases.
    Module-scoped because `build_site` costs real time (validation plus every
    builder) and every test in this file that needs a real build asks the same
    question of the same output.
    """
    base = tmp_path_factory.mktemp("landing-real-build")
    root = _checkout(base / "repo")
    out = base / "dist"
    build_site(root, out)
    return out


def test_the_page_is_published_through_write_text_and_reaches_the_checksums(
    tmp_path: Path,
) -> None:
    """`build_landing` writes exactly one file, at the path `LANDING` names.

    Through `emitter.checksums` rather than only the file on disk, because that
    mapping is what `write_manifest` publishes as `files` — a page written by a
    method that bypassed `_write` would still be on disk and invisible to the
    manifest, which is the exact failure `Emitter.write_text` exists to close.
    """
    emitter = Emitter(root=tmp_path)
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    build_landing(
        corpus,
        symbols={TBX5: "TBX5"},
        validity={TBX5: uncurated()},
        published={TBX5},
        census=_NO_BURDEN,
        profile_census=_NO_PROFILES,
        cell_type_census=_NO_CELL_TYPES,
        emitter=emitter,
        sources=_SOURCES,
    )

    assert list(emitter.checksums) == ["index.html"]
    assert (tmp_path / "index.html").is_file()


def _panel_graded(classification: Classification) -> GeneValidity:
    """A gene one ClinGen panel graded in scope, at `classification`."""
    return GeneValidity(
        records=(
            ValidityRecord(
                source=ValiditySource.CLINGEN,
                classification=classification,
                classification_term=classification.value.title(),
                disease="MONDO:0007732",
                disease_label="Holt-Oram syndrome",
                moi="AD",
                report_url=None,
                gcep="Congenital Heart Disease",
            ),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )


def _submitter_only() -> GeneValidity:
    """A gene two GenCC submitters assert and no ClinGen panel has graded."""
    return GeneValidity(
        records=tuple(
            ValidityRecord(
                source=ValiditySource.GENCC,
                classification=Classification.LIMITED,
                classification_term="Limited",
                disease="MONDO:0007732",
                disease_label="Holt-Oram syndrome",
                moi="Autosomal dominant",
                report_url=None,
                submitter=submitter,
            )
            for submitter in ("G2P", "PanelApp Australia")
        ),
        state=ValidityState.SUBMITTER_CURATED,
        has_source_discordance=False,
    )


def test_the_hero_note_splits_the_published_genes_by_which_warrant_admitted_them(
    tmp_path: Path,
) -> None:
    """The front page's one sentence about provenance, and it was false for 16 genes.

    It read "Every gene here is published on an upstream expert panel's
    classification" until 2026-08-06. That was true of all 23 genes published
    under the old rule and false the moment the gate widened: 16 of the 92 are
    admitted on GenCC submitter agreement and no ClinGen panel has graded them
    at all. The front door is where a reader forms the belief the rest of the
    site has to correct.

    **Measured 2026-08-06: reverting that sentence survived all 799 tests.** No
    test read the hero note, so the fix for it was unguarded -- a fix whose test
    passes with the fix removed is not a fix (CLAUDE.md section 8.3).

    The three figures are deliberately all different -- 3 panel-graded, 1
    submitter-admitted, 4 published -- because a fixture where two of them
    coincide cannot tell the sentence apart from one wired to the wrong
    population, and every gene fixture sharing one value is how three page
    mutants survived this suite once already (section 4.36). The panel-graded
    genes are graded at three different rungs for the same reason: the count is
    of genes a panel graded, not of genes it graded `Definitive`.
    """
    validity = {
        "HGNC:1": _panel_graded(Classification.DEFINITIVE),
        "HGNC:2": _panel_graded(Classification.MODERATE),
        "HGNC:3": _panel_graded(Classification.LIMITED),
        "HGNC:4": _submitter_only(),
    }

    text = _build(
        Corpus(root=Path("."), assertions=()),
        {},
        validity,
        tmp_path,
        published={"HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"},
    )

    note = re.search(r'<p class="hero-note">.*?</p>', text, re.S)
    assert note, "the front page carries no hero note"
    sentence = " ".join(note.group(0).split())

    assert "3 on a ClinGen expert panel's classification" in sentence
    assert "1 on two or more Gene Curation Coalition submitters" in sentence
    # The claim the widening made false must not come back (section 4.35: the
    # negative assertion is the half nobody is watching for).
    assert "Every gene here is published on an upstream expert panel" not in text
    # D49: the figures count every rare-variant test, and a de novo enrichment is
    # not a burden test, so the front page never calls them "burden statistics".
    # Per study, not "per independent cohort family" (review, 2026-10-27): rows
    # are each one study's, and since I3 one family holds two studies.
    assert (
        "Rare-variant statistics are reported per study, and none is pooled across studies."
        in sentence
    )
    assert "reported per independent cohort family" not in text
    assert "burden statistic" not in text.lower()
    # D49's definition of a rare-variant statistic, and of an independent
    # dataset, moved to `methods.html#rare-variants` with the doorway page
    # (2026-10-27); the note is how a reader of the figures above it gets there.
    # "these figures", not "them": the link must cover the dataset figure's
    # qualifier as well as the statistics the sentence before it names.
    assert sentence.endswith(
        'none is pooled across studies. <a href="methods.html#rare-variants">How to read '
        "these figures</a>.</p>"
    )


def test_the_published_count_is_derived_and_the_mirrored_count_is_not_stated(
    tmp_path: Path,
) -> None:
    """Different inputs must publish different numbers, and the mirror's count none.

    A landing page whose numbers were hardcoded strings would pass a test built
    against one fixture and only that one; varying the fixture and checking the
    rendered counts move with it is what tells a literal from a derivation. Two
    assertions on two genes, three published genes and five mirrored ones --
    none of which match the committed corpus's one assertion or its 119
    published genes -- so a hardcoded figure fails here.

    `published` is deliberately not the asserted genes. Since D21 those are
    separate populations, and the page used to count `{assertion.gene for ...}`
    under this very label, which would render 2 here instead of 3.

    **The mirrored-validity row is gone** (doorway page, 2026-10-27): every gene
    with mirrored validity in CHD scope, under a 32-word label explaining which
    of them are browsable. The rule it explained is `methods.html#scope`'s.
    Five mirrored genes against three published is what makes "5" mean that row.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002", gene=GATA4)),
    )
    validity = {gene: uncurated() for gene in ("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5")}

    text = _build(
        corpus,
        {TBX5: "TBX5", GATA4: "GATA4"},
        validity,
        tmp_path,
        published={TBX5, GATA4, "HGNC:1"},
    )

    assert "2 curated gene-disease assertions so far" in _prose(text)
    assert _figure(3, "genes published") in _hero(text)
    assert "mirrored validity" not in text
    assert ">5<" not in text


def test_the_published_gene_count_agrees_with_a_real_build_of_genes_index_json(
    real_build: Path,
) -> None:
    """The page's "genes published" figure must equal `len(genes/index.json["genes"])`.

    This is the regression the false "154 genes" claim was: `mirrored_gene_count
    = len(validity)` counted every gene either mirror curates, not the genes the
    build actually publishes. `genes/index.json` is keyed on `published`
    (`derive.gene_facts`). Checked here against the real built file rather than
    against a literal, so it fails the same way the original bug would if it
    recurred, and so that widening the gate again moves both numbers together or
    fails.
    """
    text = (real_build / "index.html").read_text(encoding="utf-8")
    genes_index = json.loads((real_build / "genes" / "index.json").read_text(encoding="utf-8"))
    published = len(genes_index["genes"])

    assert _figure(f"{published:,}", "genes published") in _hero(text)


def test_the_landing_page_distinguishes_published_from_curated(tmp_path: Path) -> None:
    """The promotion review forced this distinction into prose once already.

    92 genes are published on an upstream authority's classification -- 76 on a
    ClinGen expert panel's and 16 on two or more agreeing GenCC submitters --
    and 1 carries the atlas's own evidence. A single "genes" count would read
    as coverage the site does not have.

    The fixture separates all three numbers that could be wired to this row --
    two assertions, on one gene, against three published genes -- because the
    two plausible miswirings both render a number that is *true of something*:
    `len(corpus.assertions)` gives 2 and `len(published)` gives 3, and either
    would pass a fixture where they happened to coincide. Only 1 is the genes
    the atlas has curated.

    The sentence naming those genes is checked here too: it and the row are one
    derivation in `_render`, and a fixture with more assertions than genes is
    where a second derivation would show up as "2 ... so far: TBX5" beside a
    row saying 2.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002")),
    )

    text = _build(
        corpus,
        {TBX5: "TBX5"},
        {TBX5: uncurated()},
        tmp_path,
        published={TBX5, GATA4, "HGNC:1"},
    )
    published_section = _section(text, "<h2>What's published</h2>")

    assert _figure(3, "genes published") in _hero(text)
    assert re.search(r"<dt>Genes the atlas has curated</dt>\s*<dd>1</dd>", published_section)
    assert "2 curated gene-disease assertions so far: TBX5 (HGNC:11604)." in _prose(
        published_section
    )


def test_the_burden_census_reaches_the_page_once_per_figure(
    tmp_path: Path,
) -> None:
    """The hero band and the cards state one census, and no figure twice.

    This is the defect the whole change exists to fix, and it was an omission
    rather than a wrong number: measured 2026-08-05 against the live site, the
    front page's census published `assertions: 1, functional: 0, datasets: 0`
    for a build carrying 290 burden statistics over 23 genes from 3 independent
    datasets, and the word "burden" did not appear on it once. Every figure was
    true; the picture they composed was false.

    **All four figures the page states are deliberately distinct**, because a
    fixture where two of them coincide cannot see them swapped. The first version
    of this test gave the census 2 families and 2 genes-with-evidence, and the
    mutant putting `census.genes` in the hero where `census.families` belongs
    survived the whole file. Measured, then fixed: 6 statistics, 3 independent
    datasets, 2 genes with a reported result, 4 genes published -- and, since
    owner decision 2026-10-26 (R9), 1 gene with a result surviving its study's
    correction (`HGNC:2`'s one row; none of `HGNC:1`'s five), so the new figure
    differs from every other (CLAUDE.md section 4.30). Its label is not the
    owner's "Genes tested in rare-variant studies": the figure counts genes a
    study reports a result for, and the gene it leaves out on a real build
    (CRIPTO) is one a study may have tested and not reported.

    The third family reports on no published gene, which is what makes it
    different from the gene count. That is not contrived -- it is a study whose
    genes all sit outside the publication gate, and 127 of the burden mirror's
    150 genes are in exactly that position today.

    **The dataset count says what it can see** (owner decision 2026-10-26, R12):
    "(as far as the cohort registry records)", under the hero figure and on no
    other figure. The families are joined only through collections the registry
    names or declares overlapping, so case sources no row names are invisible to
    them. Pinned as a literal, never as the constant that renders it.

    **Each figure once** (doorway page, 2026-10-27): the three hero figures had
    a card each below them as well, so a reader met them twice. The cards now
    carry only what the hero does not.

    The census is built by `burden_census` from a real fixture rather than
    constructed here, so this test cannot pass against a derivation the module
    would get wrong; `test_build_burden.py` owns that derivation's own guard.
    """
    burden = {
        "HGNC:1": [_burden_row(study="PMID:1") for _ in range(5)],
        "HGNC:2": [
            _burden_row(
                study="PMID:2", pvalue_adjusted=0.01, pvalue_adjustment="benjamini_hochberg"
            )
        ],
    }
    published = {"HGNC:1", "HGNC:2", "HGNC:3", TBX5}
    families = (frozenset({"PMID:1"}), frozenset({"PMID:2"}), frozenset({"PMID:3"}))
    census = burden_census(burden, published, families, survives=survives_own_correction)
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path, published, census)
    hero = _hero(text)
    listed = _section(text, "<h2>What's published</h2>")

    for value, label in (
        (4, "genes published"),
        (3, "independent datasets"),
        (6, "rare-variant statistics"),
    ):
        assert _figure(value, label) in hero, f"the hero no longer states {value} beside {label!r}"
    qualifier = "(as far as the cohort registry records)"
    assert (
        f'{_figure(3, "independent datasets")}/span> <span class="figure-note">{qualifier}</span>'
    ) in hero
    assert text.count(qualifier) == 1
    assert re.search(r"<dt>Genes with a reported rare-variant result</dt>\s*<dd>2</dd>", listed)
    assert re.search(
        r"<dt>Genes with a rare-variant result surviving its study's correction</dt>"
        r"\s*<dd>1</dd>",
        listed,
    )
    assert "rare-variant evidence" not in listed
    assert "Genes tested" not in listed
    for repeated in ("Genes published", "Rare-variant statistics", "Independent datasets"):
        assert f"<dt>{repeated}</dt>" not in listed, f"{repeated!r} is in the hero already"


def _study(pmid: str, *, own_lab: bool) -> Publication:
    """A curated rare-variant study; only its id and `own_lab` are read here."""
    return Publication(
        id=pmid,
        title=f"Study {pmid}",
        journal="J",
        year=2026,
        authors=["A"],
        study_type="case_control",  # type: ignore[arg-type]
        own_lab=own_lab,
    )


@pytest.mark.parametrize(
    ("tested", "by_an_author", "sentence"),
    [
        pytest.param(
            ("PMID:1", "PMID:2", "PMID:3"),
            {"PMID:1"},
            "Of the 3 rare-variant studies whose results appear here, 1 is by an author of this "
            "atlas, and a gene page says so wherever it shows its results.",
            id="one of several",
        ),
        pytest.param(
            ("PMID:1", "PMID:2", "PMID:3"),
            {"PMID:1", "PMID:2"},
            "Of the 3 rare-variant studies whose results appear here, 2 are by an author of this "
            "atlas, and a gene page says so wherever it shows their results.",
            id="several",
        ),
        pytest.param(
            ("PMID:1", "PMID:2", "PMID:3"),
            set(),
            "Of the 3 rare-variant studies whose results appear here, none is by an author of "
            "this atlas.",
            id="none",
        ),
        pytest.param(
            ("PMID:1",),
            {"PMID:1"},
            "Of the 1 rare-variant study whose results appear here, 1 is by an author of this "
            "atlas, and a gene page says so wherever it shows its results.",
            id="a single study",
        ),
        pytest.param((), set(), None, id="no study"),
    ],
)
def test_the_front_page_declares_the_rare_variant_studies_by_an_atlas_author(
    tmp_path: Path, tested: tuple[str, ...], by_an_author: set[str], sentence: str | None
) -> None:
    """One sentence, both counts computed (owner decision 2026-10-26, R7).

    The cards and the study blocks declare a study by an author of this atlas
    and the front page said nothing. The population is the studies whose rows
    reach a published bundle (`BurdenCensus.studies`), never every curated
    publication: two decoys are by an atlas author and are not among them --
    `PMID:9` tests only a gene outside the gate, and `PMID:8` tests nothing --
    so counting `corpus.publications` instead reads 3 of 5, 4 of 5 and 2 of 5
    in the first three cases, where the sentence must read 1, 2 and none of 3
    (CLAUDE.md section 4.30: two figures equal on the real corpus would let
    the wrong one pass). Pinned in each of its forms, as literals.

    **And at both ends of the count** (review, 2026-10-27): a single study
    takes the singular noun, and a build publishing no study's rows carries no
    sentence at all -- there is nothing to declare, and "Of the 0 rare-variant
    studies" would say there is. Both survived every test until then. "by an
    author", in every form, as each block's declaration says it.

    **Shortened with the doorway page, then partly restored** (2026-10-27).
    "published here" read as though this atlas had published the studies,
    which is the one misreading a conflict-of-interest line cannot carry, and
    "its block" named a thing the front page never defines. "whose results
    appear here" and "wherever it shows its results" say what is true of each.
    """
    burden = {
        "HGNC:1": [_burden_row(study=pmid) for pmid in tested[:2]],
        "HGNC:2": [_burden_row(study=pmid) for pmid in tested[2:]],
        # Outside the gate: mirrored, never published.
        "HGNC:9": [_burden_row(study="PMID:9")],
    }
    published = {"HGNC:1", "HGNC:2", TBX5}
    families = tuple(frozenset({f"PMID:{n}"}) for n in (1, 2, 3, 9))
    census = burden_census(burden, published, families, survives=survives_own_correction)
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(),),
        publications=tuple(
            _study(pmid, own_lab=pmid in by_an_author or pmid in {"PMID:8", "PMID:9"})
            for pmid in ("PMID:1", "PMID:2", "PMID:3", "PMID:8", "PMID:9")
        ),
    )

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path, published, census)
    listed = _prose(_section(text, "<h2>What's published</h2>"))

    if sentence is None:
        assert "rare-variant stud" not in listed
    else:
        assert f"<p>{sentence}</p>" in listed


@pytest.mark.parametrize(
    ("panels", "note"),
    [
        pytest.param(
            {"PMID:2": 248},
            "For PMID:2, the correction a result must survive ran over its 248-gene panel, not "
            "the genome.",
            id="one published study corrected over a panel",
        ),
        pytest.param(
            {"PMID:1": 1200, "PMID:3": 248},
            "For PMID:1 and PMID:3, the correction a result must survive ran over each study's "
            "gene panel (1,200 genes for PMID:1, 248 genes for PMID:3), not the genome.",
            id="two, each with its own panel",
        ),
        pytest.param({"PMID:9": 248}, None, id="the only such study is outside the gate"),
        pytest.param({}, None, id="no study corrected over a panel"),
    ],
)
def test_the_front_page_says_whose_correction_ran_over_a_panel_not_the_genome(
    tmp_path: Path, panels: dict[str, int], note: str | None
) -> None:
    """The surviving count's qualifier (review, 2026-10-27; owner decision R3).

    "Genes with a rare-variant result surviving its study's correction" is the
    owner's label and reads as genome-wide, and 22 of the 32 genes it counted on
    a build of 5f7161b survive only PMID:40127276's Benjamini-Hochberg over its
    248 panel genes. So the page names each published study whose correction
    ran over a panel, from its curated `correction_genes`. Only the studies
    whose rows reach a published bundle: `PMID:9` tests only a gene outside the
    gate, so a note naming it would qualify a count it is not in. The note's
    link to `methods.html#rare-variants` went with the doorway page
    (2026-10-27): the hero note links there already.
    """
    burden = {
        "HGNC:1": [_burden_row(study="PMID:1"), _burden_row(study="PMID:2")],
        "HGNC:2": [_burden_row(study="PMID:3")],
        "HGNC:9": [_burden_row(study="PMID:9")],
    }
    published = {"HGNC:1", "HGNC:2", TBX5}
    families = tuple(frozenset({f"PMID:{n}"}) for n in (1, 2, 3, 9))
    census = burden_census(burden, published, families, survives=survives_own_correction)
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(),),
        publications=tuple(
            _study(pmid, own_lab=False).model_copy(update={"correction_genes": panels.get(pmid)})
            for pmid in ("PMID:1", "PMID:2", "PMID:3", "PMID:9")
        ),
    )

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path, published, census)
    listed = _prose(_section(text, "<h2>What's published</h2>"))

    if note is None:
        assert "not the genome" not in text
    else:
        assert f"<p>{note}</p>" in listed


def test_the_landing_card_names_the_layer_the_build_publishes(tmp_path: Path) -> None:
    """Adding a whole evidence layer without touching the front page is how
    `index.html` once advertised an empty atlas for three releases with the
    word "burden" appearing on it zero times (CLAUDE.md section 4.29). Pinned
    against a literal, never against a constant a rewrite could carry along
    with it.

    **Neither profile count is distinguishable on the committed corpus.** This
    said "every profile count is 0 -- there is no `profiles` mirror to derive
    one from" until 2026-08-20; `mirrors/profiles/E-MTAB-6814.tsv` landed on
    2026-08-19 and the real build now publishes `profile_genes: 92` and
    `profile_datasets: 1`. Those are also `genes: 92` and `datasets: 1`, so a
    card wired to the wrong count still reads correctly there. `genes=3,
    datasets=2` here is what tells a real card from one reading a neighbour or
    a literal. `test_every_figure_the_page_states_is_the_manifests_own` is what
    proves this same wiring holds on the real build; this is what proves it is
    not two coincidences wired to each other.

    The datasets card went with the doorway page (2026-10-27): one dataset
    today, and `manifest.json`'s `profile_datasets` states it.

    **The cell-type card is held the same way** (final review, 2026-10-28).
    Both cell-type counts are 0 on the committed corpus, so a card reading
    `datasets` instead of `genes` survived every test; `genes=7, datasets=3`
    separates them, and from the profile card's 3 and 2.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(
        corpus,
        {TBX5: "TBX5"},
        {TBX5: uncurated()},
        tmp_path,
        profile_census=ProfileCensus(genes=3, datasets=2),
        cell_type_census=CellTypeCensus(genes=7, datasets=3),
    )
    listed = _section(text, "<h2>What's published</h2>")

    assert re.search(r"<dt>Genes with developmental expression data</dt>\s*<dd>3</dd>", listed)
    assert re.search(r"<dt>Genes with cell-type expression data</dt>\s*<dd>7</dd>", listed)


def test_the_front_page_is_a_doorway_and_its_reading_notes_live_elsewhere(
    tmp_path: Path,
) -> None:
    """What the doorway page moved off the front page stays off it (2026-10-27).

    The page carried 768 words and 17 cards, and a first-time reader had to get
    through them to reach the genes. Each block it dropped is stated in one
    other place, which is where it is now pinned:

    * the glyph key and its power caveat, on the browse page and every gene page
      beside the glyphs themselves (`test_build_pages.py`);
    * the definitions of a rare-variant statistic and an independent dataset,
      on `methods.html#rare-variants` (`test_build_methods.py`);
    * "What this is", which restated the tagline and the hero note;
    * thirteen of its 17 cards -- three repeating the hero -- in `manifest.json`'s
      `counts`.

    Asserted absent here because a pasted old block is how they would come back
    (CLAUDE.md section 4.35: the negative assertion is the half nobody watches).
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert 'class="legend-item"' not in text
    assert "not evidence against a gene" not in text
    assert "<strong>rare-variant statistic</strong>" not in text
    assert "<strong>independent dataset</strong>" not in text
    assert "<h2>What this is</h2>" not in text
    assert re.findall(r"<h2>([^<]+)</h2>", text) == ["What's published", "Sources", "Data"]


def test_the_browse_page_is_the_pages_call_to_action(tmp_path: Path) -> None:
    """The browse page is the product; it used to be a bullet under two sections.

    Pinned inside the hero rather than anywhere on the page, because the `<nav>`
    links to the same file and would satisfy a bare substring check while the
    call to action disappeared.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert '<a class="cta" href="genes/index.html">' in _hero(text)


def test_every_figure_the_page_states_is_the_manifests_own(
    real_build: Path,
) -> None:
    """Three artifacts, one set of numbers, checked against a real build.

    The page states seven keys of `counts` and the manifest all of them
    (`docs/data-api.md`; the page stated every key until the doorway page,
    2026-10-27): six as figures -- `cell_type_genes` joined them with schema
    2.20 -- and `assertions` in the sentence naming the curated genes. They are
    derived once in `build_site` (the burden figures from a single
    `burden_census` call) and handed to both precisely so they cannot disagree
    -- and a test that only read fixtures could not see them start to, because
    the drift would be in the wiring rather than in either derivation.

    `burden_rows` is checked against the gene bundles too, which is the claim
    that actually matters: the figure must count rows a consumer can *fetch*, not
    rows in `mirrors/burden.tsv`. Measured 2026-08-20: the mirror holds 1,475
    rows over 150 genes, 59 of those genes publish no page, and 915 rows reach
    a bundle -- so a census counting the mirror would advertise 1.6 times the
    evidence the site serves. Summed over the published bundles rather than
    over the mirror for that reason. (This read "127 of those genes" and "five
    times" until 2026-08-20: both were measured on the 23-gene corpus and both
    moved when the gate widened to 92 on 2026-08-06.)

    **Every figure the page states, not a sample.** This checked three keys
    until review 2026-08-06 pointed out that the invariant was claimed for the
    whole census and guarded for a third of it -- so a cross-wire between two
    cards publishes a false front-page figure with the manifest still correct
    and nothing failing. Every card and hero figure on the page must be one of
    the labels below, each with a value derived without `landing.py`, so a
    figure added without one fails here rather than going unguarded. The two
    cards that are not manifest keys are held to what the build published or
    read instead (review, 2026-10-27: they were listed here and never checked):
    the genes `genes/index.json` gives a rare-variant row, and the distinct
    genes the curation's assertions name, published or not, as the card counts.
    """
    text = (real_build / "index.html").read_text(encoding="utf-8")
    manifest = json.loads((real_build / "manifest.json").read_text(encoding="utf-8"))
    counts = manifest["counts"]
    bundles = sorted((real_build / "genes").glob("HGNC_*.json"))
    fetchable = sum(len(json.loads(path.read_text(encoding="utf-8"))["burden"]) for path in bundles)

    assert counts["burden_rows"] == fetchable
    assert counts["genes"] == len(bundles)
    index = json.loads((real_build / "genes" / "index.json").read_text(encoding="utf-8"))
    curation, issues = load_curation(real_build.parent / "repo")
    assert issues == [], issues

    # The label the page prints for each figure, and the value it must state.
    hero = {
        "genes published": "genes",
        "independent datasets": "cohort_families",
        "rare-variant statistics": "burden_rows",
    }
    cards = {
        "Genes with a rare-variant result surviving its study's correction": counts[
            "burden_surviving_genes"
        ],
        "Genes with developmental expression data": counts["profile_genes"],
        "Genes with cell-type expression data": counts["cell_type_genes"],
        "Genes with a reported rare-variant result": sum(
            1 for row in index["genes"] if row["burden_row_count"]
        ),
        "Genes the atlas has curated": len({assertion.gene for assertion in curation.assertions}),
    }
    figures = re.findall(
        r'"figure-value">([\d,]+)</span>\s*<span class="figure-label">([^<]+)<', text
    )
    stated = re.findall(r"<dt>([^<]+)</dt>\s*<dd>([\d,]+)</dd>", text)
    assert {label for _, label in figures} == set(hero)
    assert {html.unescape(label) for label, _ in stated} == set(cards)
    for value, label in figures:
        assert int(value.replace(",", "")) == counts[hero[label]], label
    for label, value in stated:
        assert int(value.replace(",", "")) == cards[html.unescape(label)], label
    assert re.search(
        rf"<p>{counts['assertions']} curated gene-disease assertions? so far: ", text
    ), "the assertion sentence does not state the manifest's `assertions`"


def test_the_page_shares_the_shell_and_the_stylesheet_every_other_page_uses(
    tmp_path: Path,
) -> None:
    """The landing page carries no `<style>` of its own.

    It had one until the gene pages arrived, and three page kinds each carrying
    their own copy is three stylesheets that drift. `render.STYLESHEET` is now
    the only one, so this asserts the rendered page contains it verbatim rather
    than merely that some `<style>` block exists -- a landing page that grew a
    second block back would still satisfy the weaker check.

    The `<nav>` is the other half of the shell, and `root=""` is what makes its
    links right from the site root: a landing page rendered with `"../"` would
    point every visitor's "Genes" link above the deploy root.

    Its third link, to `methods.html` (D50, 2026-10-26), is what puts the
    site's reading notes one click from every page, the gene pages included.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert f"<style>{STYLESHEET}</style>" in text
    assert text.count("<style>") == 1
    assert '<a href="index.html">CHD Atlas</a>' in text
    assert '<a href="genes/index.html">Genes</a>' in text
    assert '<a href="methods.html">How to read this atlas</a></nav>' in text


def test_the_data_line_links_each_payload_a_consumer_starts_from(tmp_path: Path) -> None:
    """One line where "Browse the data" was a five-item list (doorway page, 2026-10-27).

    Its first item, `genes/index.html`, is the hero's call to action and the
    `<nav>`'s second link; `sources.json` is linked from "Sources". What is left
    is for a program, and each link is pinned inside the section, because the
    footer and "Sources" link nearby files and would satisfy a page-wide check.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)
    data = _section(text, "<h2>Data</h2>")
    data = data[: data.index("<footer>")]

    for href in ("genes/index.json", "manifest.json", "https://github.com/chd-atlas"):
        assert f'href="{href}"' in data, href


def test_the_front_page_names_orphanet_and_attributes_it_as_its_licence_requires(
    tmp_path: Path,
) -> None:
    """Orphanet appeared on the front page zero times until 2026-10-26.

    Its phenotype frequencies are what bring Noonan, Costello, Adams-Oliver and
    the other nominated syndromes into scope, and for some published genes
    those syndromes are the only route in; `sources.json` records both Orphadata
    products as CC-BY-4.0, attribution required, and Orphadata asks to be cited
    in a stated form. Pinned against a literal in "Sources", whitespace-collapsed
    (the template wraps it), so the citation cannot be shortened: the citation
    is what the licence asks for, not the name.

    **Pinned whole, since the doorway page (2026-10-27)**, when the section went
    from 193 words to 123 in one paragraph: "The atlas authors no validity
    classification of its own" is the hero note's claim already, the mirrored-
    validity row is gone, and the other sources each take one clause.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))
    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    sources_section = _prose(_section(text, "<h2>Sources</h2>"))
    assert sources_section == (
        "<h2>Sources</h2> <p>Gene-disease validity is mirrored from <strong>ClinGen</strong> "
        "and <strong>GenCC</strong>, and phenotype terms from the <strong>Human Phenotype "
        "Ontology</strong>. Mondo's disease mappings are mirrored under CC BY 4.0. "
        "<strong>Orphanet</strong>'s cardiovascular phenotype annotations, which decide "
        "whether a nominated syndrome is in scope, and its OMIM mappings are mirrored under "
        "CC BY 4.0 — <em>Orphadata Science: Free access data from Orphanet. &copy; INSERM "
        '1999. Available on <a href="http://sciences.orphadata.com/">http://sciences.orphadata.com/'
        "</a>. Data version 1.3.42 / 4.1.8 [2025-03-03], 2026-06-23.</em> "
        # Every other attribution-required source (owner decision 2026-10-26,
        # R10); `test_built_site_is_consumable.py` holds the set it names to
        # `sources.json`.
        "Developmental expression is the <strong>EMBL-EBI Expression Atlas</strong>'s "
        "(E-MTAB-6814), under EMBL-EBI's Terms of Use. Rare-variant results come from the "
        "supplements of Audain et al. 2026 and Sierant et al. 2025, both CC BY 4.0, and of "
        'Audain et al. 2021, CC0. <a href="sources.json">sources.json</a> identifies each '
        "mirrored file and its licence terms.</p> "
    )
    # ClinGen's, GenCC's and HGNC's entries are dates ("FILE CREATED: 2026-08-03",
    # "retrieved 2026-08-03", "2026-10-16"), so "by version or checksum" was
    # false of three sources (review, 2026-10-26).
    assert "by version or checksum" not in sources_section


@pytest.mark.parametrize(
    ("sources", "stated"),
    [
        pytest.param(
            SourceRegistry(
                sources=[
                    _orphadata(
                        "orphanet_product4", "1.3.43 / 4.1.8 [2026-01-01]", "2026-12-01 01:00:00"
                    ),
                    _orphadata(
                        "orphanet_product1", "1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:53:50"
                    ),
                ]
            ),
            "Data version 1.3.43 / 4.1.8 [2026-01-01], 2026-12-01 (product 4); "
            "1.3.42 / 4.1.8 [2025-03-03], 2026-06-23 (product 1).</em>",
            id="two releases: each product named",
        ),
        pytest.param(
            SourceRegistry(
                sources=[_orphadata("orphanet_product4"), _orphadata("orphanet_product1")]
            ),
            "http://sciences.orphadata.com/</a>.</em>",
            id="no header recorded: nothing to cite",
        ),
        pytest.param(
            SourceRegistry(
                sources=[
                    _orphadata(
                        "orphanet_product4", "1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:57:18"
                    ),
                    _orphadata("orphanet_product1"),
                ]
            ),
            "Data version 1.3.42 / 4.1.8 [2025-03-03], 2026-06-23 (product 4).</em>",
            id="one product's header: named, never cited as both",
        ),
        pytest.param(
            SourceRegistry(
                sources=[
                    _orphadata("orphanet_product4", "1.3 <b> & 4", "2026<12>01 01:00:00"),
                    _orphadata("orphanet_product1", "1.3 <b> & 4", "2026<12>01 02:00:00"),
                ]
            ),
            "Data version 1.3 &lt;b&gt; &amp; 4, 2026&lt;12&gt;01.</em>",
            id="markup in a header is escaped, the date as the version",
        ),
        pytest.param(
            SourceRegistry(
                sources=[
                    _orphadata("orphanet_product4", "1.3 <b>", "2026<12>01 01:00:00"),
                    _orphadata("orphanet_product1", "1.2 & 4", "2025<06>23 02:00:00"),
                ]
            ),
            "Data version 1.3 &lt;b&gt;, 2026&lt;12&gt;01 (product 4); 1.2 &amp; 4, "
            "2025&lt;06&gt;23 (product 1).</em>",
            id="and escaped per product where the releases differ",
        ),
    ],
)
def test_the_orphadata_citation_states_the_release_each_product_was_mirrored_from(
    tmp_path: Path, sources: SourceRegistry, stated: str
) -> None:
    """Orphadata's citation asks for its "Data version" (owner decision 2026-10-26, R11).

    The page printed none and pointed to `sources.json`, which identified both
    files by sha256 alone. Each product's own `<JDBOR>` header version and date
    are now in `mirrors/sources.yaml`, written there by the converter, and the
    citation prints them: once where both products are one release on one day
    (the literal in the test above), each by name where they differ, and
    nothing where no header was recorded -- never a version one product
    carries attributed to both, which is also why a release only one product
    records is printed against that product's name (review, 2026-10-27: the
    `len(stated) == len(_ORPHADATA)` guard that keeps it so survived every
    test). Both header fields are a curator's transcription interpolated into
    HTML, so both are escaped; the date reached the page unescaped until then.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))
    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path, sources=sources)

    cited = _prose(_section(text, "<h2>Sources</h2>"))
    assert stated in cited
    if "Data version" not in stated:
        assert "Data version" not in cited


def test_a_single_assertion_is_singular_not_plural(tmp_path: Path) -> None:
    """The wording matters as much as the number: "1 assertions" reads as unreviewed."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert "1 curated gene-disease assertion so far" in text
    assert "1 curated gene-disease assertions" not in text


def test_no_assertions_yet_reads_as_none_rather_than_as_an_empty_list(tmp_path: Path) -> None:
    """A corpus with nothing curated still has to publish a page, and an honest one."""
    corpus = Corpus(root=Path("."), assertions=())

    text = _build(corpus, {}, {}, tmp_path)

    assert "0 curated gene-disease assertions so far" in text
    assert "none yet" in text
    assert _figure(0, "genes published") in _hero(text)


def test_the_research_use_statement_reaches_the_front_page_exactly_once(
    tmp_path: Path,
) -> None:
    """The development-status notice: exactly two sentences, once, and no more.

    The paragraph this notice used to carry beneath them -- "Today it holds {N}
    curated gene-disease assertion..." -- editorialised with a claim false for
    100% of described genes ("almost every gene this site describes carries an
    upstream expert panel's or submitter's classification and nothing this
    atlas has independently curated"): the committed corpus describes one gene,
    and it is the curated one. That paragraph is gone, and this page has no
    editorial prose left between the notice and the hero.

    The notice itself now belongs to `render.document`, which puts it on all 25
    pages instead of this one (see `test_build_render.py`). What is still this
    file's business is that the landing page did not keep a second copy: the
    section that used to hold it was removed here rather than left to render
    alongside the shared one, and two copies of the site's only research-use
    statement is how one of them comes to be edited and the other not.

    Counted on `RESEARCH_USE_NOTICE` rather than on the substring "not a
    clinical decision-support tool", because the emphasis markup is part of what
    the owner specified and a page carrying the sentence with the `<strong>`
    dropped would satisfy the looser check.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert text.count(RESEARCH_USE_NOTICE) == 1
    assert "<h2>Development status &amp; research use</h2>" not in text
    assert "Today it holds" not in text
    assert "independently curated" not in text


def test_a_gene_symbol_is_escaped_before_it_reaches_the_page(tmp_path: Path) -> None:
    """A mirror cell is curated text, not markup — the same discipline `bundles.py`
    and `search.py` already apply to a symbol before it reaches a payload.

    `mirrors/genes.tsv` is reviewed for scientific accuracy, not for HTML safety,
    so nothing upstream of this page guarantees a symbol contains no `<`. Without
    `html.escape` this fixture's symbol would open a real `<script>` tag in the
    page every visitor's browser parses.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "<script>alert(1)</script>"}, {TBX5: uncurated()}, tmp_path)

    assert "<script>alert(1)</script>" not in text
    assert "&lt;script&gt;" in text


def test_a_gene_absent_from_the_registry_falls_back_to_its_hgnc_id(tmp_path: Path) -> None:
    """The same registry fallback `build_genes` and `search.py` already apply."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {}, {TBX5: uncurated()}, tmp_path)

    assert TBX5 in text


def test_the_page_makes_no_external_request(tmp_path: Path) -> None:
    """Self-contained: no script, no stylesheet or font pulled from a CDN.

    A `<a href="https://…">` is fine — that is a link a reader clicks, not a
    request the page issues on load — so this checks for the tags and
    attributes that *do* fetch automatically rather than for the substring
    "http" anywhere in the document.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert "<script" not in text
    assert "<link" not in text
    assert "@import" not in text
