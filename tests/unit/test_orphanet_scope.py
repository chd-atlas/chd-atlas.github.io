# tests/unit/test_orphanet_scope.py
"""The Orphanet frequency rule, on a synthetic ontology small enough to read.

Each test below exists to kill one way the rule could quietly admit the wrong
syndrome. Measured on the real data before the rule was written (spec
2026-10-24): applied without the structural restriction, the same floor admits
Fabry disease and Hurler syndrome through their cardiomyopathy and valve
findings, so the restriction is the part most worth pinning.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from chd_atlas.orphanet_scope import (
    Annotation,
    OrphanetTables,
    Warrant,
    load_orphanet_tables,
    pinned_hpo,
    pinned_structural_phenotypes,
    qualifying,
    qualifying_targets,
    reach,
    structural_phenotypes,
    warrant,
)
from chd_atlas.tables import MONDO_XREFS, ORPHANET_OMIM, ORPHANET_PHENOTYPES
from chd_atlas.validate.sources import load_sources

# One term per branch the rule has to tell apart: the root, a septal child with
# an alt id, the myocardium branch it excludes, and two vascular-branch terms it
# adds back. The Typedef is there because a stanza that is not a Term must reset
# the reader, or its lines would attach to the last term.
_OBO = """format-version: 1.2

[Term]
id: HP:0000118
name: Phenotypic abnormality

[Term]
id: HP:0001627
name: Abnormal heart morphology
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001671
name: Abnormal cardiac septum morphology
is_a: HP:0001627 ! Abnormal heart morphology

[Term]
id: HP:0001631
name: Atrial septal defect
alt_id: HP:0001630
is_a: HP:0001671 ! Abnormal cardiac septum morphology

[Term]
id: HP:0001637
name: Abnormal myocardium morphology
is_a: HP:0001627 ! Abnormal heart morphology

[Term]
id: HP:0001638
name: Cardiomyopathy
is_a: HP:0001637 ! Abnormal myocardium morphology

[Term]
id: HP:0030962
name: Abnormal morphology of the great vessels
is_a: HP:0000118 ! Phenotypic abnormality

[Term]
id: HP:0001643
name: Patent ductus arteriosus
is_a: HP:0030962 ! Abnormal morphology of the great vessels

[Term]
id: HP:0001680
name: Coarctation of aorta
is_a: HP:0030962 ! Abnormal morphology of the great vessels

[Typedef]
id: part_of
is_a: HP:0001627
"""

ASD = Annotation("HP:0001631", "Atrial septal defect", "Frequent (79-30%)")
HCM = Annotation("HP:0001638", "Cardiomyopathy", "Very frequent (99-80%)")
PDA = Annotation("HP:0001643", "Patent ductus arteriosus", "Frequent (79-30%)")


@pytest.fixture
def structural(tmp_path: Path) -> frozenset[str]:
    obo = tmp_path / "hp.obo"
    obo.write_text(_OBO, encoding="utf-8")
    return structural_phenotypes(obo)


def _tables(
    annotations: dict[str, tuple[Annotation, ...]] | None = None,
    omim_links: dict[tuple[str, str], tuple[str, str]] | None = None,
    xrefs: dict[str, frozenset[str]] | None = None,
) -> OrphanetTables:
    return OrphanetTables(
        annotations=annotations or {},
        names={"ORPHA:648": "Noonan syndrome"},
        omim_links=omim_links or {},
        xrefs=xrefs or {},
    )


def test_the_structural_set_excludes_cardiomyopathy_and_adds_pda(
    structural: frozenset[str],
) -> None:
    """The definition the owner approved, on one term per branch."""
    assert {"HP:0001627", "HP:0001671", "HP:0001631", "HP:0001643", "HP:0001680"} <= structural
    # An annotation may carry a retired id; it must still count.
    assert "HP:0001630" in structural
    assert not {"HP:0001637", "HP:0001638", "HP:0030962", "HP:0000118"} & structural
    # The Typedef's `is_a` must not have attached `part_of` to the heart.
    assert "part_of" not in structural


@pytest.mark.parametrize(
    ("frequency", "admitted"),
    [
        pytest.param("Obligate (100%)", True, id="obligate"),
        pytest.param("Very frequent (99-80%)", True, id="very-frequent"),
        pytest.param("Frequent (79-30%)", True, id="frequent"),
        pytest.param("Occasional (29-5%)", False, id="occasional"),
        pytest.param(None, False, id="no-band"),
    ],
)
def test_only_frequent_or_stronger_qualifies(
    structural: frozenset[str], frequency: str | None, admitted: bool
) -> None:
    """The floor is Frequent (30%), HPO's own word for "frequent"."""
    tables = _tables({"ORPHA:648": (Annotation("HP:0001631", "ASD", frequency),)})
    assert bool(qualifying("ORPHA:648", tables, structural)) is admitted


def test_a_disorder_qualifying_only_by_cardiomyopathy_does_not_qualify(
    structural: frozenset[str],
) -> None:
    """The Fabry/Hurler case: a cardiac finding that is not a malformation."""
    assert qualifying("ORPHA:648", _tables({"ORPHA:648": (HCM,)}), structural) == ()


def test_patent_ductus_arteriosus_qualifies_though_hpo_files_it_outside_the_heart(
    structural: frozenset[str],
) -> None:
    assert qualifying("ORPHA:648", _tables({"ORPHA:648": (PDA,)}), structural) == (PDA,)


def test_qualifying_annotations_are_ordered_strongest_first(structural: frozenset[str]) -> None:
    """Rank first, then HPO id: the note on the page reads the strongest first."""
    obligate = Annotation("HP:0001680", "Coarctation of aorta", "Obligate (100%)")
    tables = _tables({"ORPHA:648": (PDA, ASD, HCM, obligate)})
    assert qualifying("ORPHA:648", tables, structural) == (obligate, ASD, PDA)


def test_a_term_reaches_its_own_orphanet_code_directly() -> None:
    tables = _tables(xrefs={"MONDO:0018997": frozenset({"ORPHA:648", "OMIM:163950"})})
    assert reach("MONDO:0018997", "ORPHA:648", tables) == ("direct", None)


@pytest.mark.parametrize(
    ("relation", "validation", "reached"),
    [
        pytest.param("BTNT", "Validated", True, id="broader-validated"),
        pytest.param("E", "Validated", True, id="exact-validated"),
        pytest.param("NTBT", "Validated", False, id="narrower"),
        pytest.param("ND", "Validated", False, id="undecided"),
        pytest.param("BTNT", "Not yet validated", False, id="unvalidated"),
    ],
)
def test_a_subtype_reaches_its_parent_only_through_a_validated_exact_or_broader_mapping(
    relation: str, validation: str, reached: bool
) -> None:
    """Noonan syndrome 4's path to ORPHA:648, under each relation product 1 uses."""
    tables = _tables(
        omim_links={("ORPHA:648", "OMIM:610733"): (relation, validation)},
        xrefs={"MONDO:0012547": frozenset({"OMIM:610733"})},
    )
    route = reach("MONDO:0012547", "ORPHA:648", tables)
    assert route == (("subtype", "OMIM:610733") if reached else None)


def test_a_warrant_names_which_half_of_the_rule_failed(structural: frozenset[str]) -> None:
    """SCP006 prints this reason, so a curator knows what to fix."""
    reachable = _tables(
        annotations={"ORPHA:648": (HCM,)},
        xrefs={"MONDO:0018997": frozenset({"ORPHA:648"})},
    )
    unreachable = _tables(annotations={"ORPHA:648": (ASD,)})

    assert "does not reach ORPHA:648" in str(
        warrant("MONDO:0018997", "ORPHA:648", unreachable, structural)
    )
    assert "no structural heart phenotype" in str(
        warrant("MONDO:0018997", "ORPHA:648", reachable, structural)
    )
    good = _tables(
        annotations={"ORPHA:648": (ASD,)}, xrefs={"MONDO:0018997": frozenset({"ORPHA:648"})}
    )
    assert warrant("MONDO:0018997", "ORPHA:648", good, structural) == Warrant(
        "ORPHA:648", "Noonan syndrome", "direct", None, (ASD,)
    )


def test_qualifying_targets_follows_both_routes(structural: frozenset[str]) -> None:
    """SCP007's search: a direct code and a validated subtype parent, both found."""
    tables = OrphanetTables(
        annotations={"ORPHA:648": (ASD,), "ORPHA:500": (PDA,), "ORPHA:9": (HCM,)},
        names={},
        omim_links={("ORPHA:500", "OMIM:1"): ("BTNT", "Validated")},
        xrefs={"MONDO:0000001": frozenset({"ORPHA:648", "ORPHA:9", "OMIM:1"})},
    )
    assert qualifying_targets("MONDO:0000001", tables, structural) == ["ORPHA:500", "ORPHA:648"]


def _write(root: Path, name: str, columns: tuple[str, ...], rows: list[tuple[str, ...]]) -> None:
    mirrors = root / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(columns), *("\t".join(row) for row in rows)]
    (mirrors / name).write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_the_three_mirrors_load_into_the_shape_the_rule_reads(tmp_path: Path) -> None:
    """And a missing one loads as None, never as empty tables (§1: data wrong vs never loaded)."""
    _write(
        tmp_path,
        "orphanet_phenotypes.tsv",
        ORPHANET_PHENOTYPES.column_names,
        [
            (
                "ORPHA:648",
                "Noonan syndrome",
                "HP:0001631",
                "Atrial septal defect",
                "Frequent (79-30%)",
                "s",
            )
        ],
    )
    _write(
        tmp_path,
        "orphanet_omim.tsv",
        ORPHANET_OMIM.column_names,
        [("ORPHA:648", "OMIM:610733", "BTNT", "Validated", "s")],
    )
    assert load_orphanet_tables(tmp_path) is None  # mondo_xrefs.tsv still absent

    _write(
        tmp_path,
        "mondo_xrefs.tsv",
        MONDO_XREFS.column_names,
        [("MONDO:0012547", "OMIM:610733", "s")],
    )
    tables = load_orphanet_tables(tmp_path)
    assert tables is not None
    assert tables.annotations == {"ORPHA:648": (ASD,)}
    assert tables.names == {"ORPHA:648": "Noonan syndrome"}
    assert tables.omim_links == {("ORPHA:648", "OMIM:610733"): ("BTNT", "Validated")}
    assert tables.xrefs == {"MONDO:0012547": frozenset({"OMIM:610733"})}


@pytest.mark.parametrize(
    ("pinned", "content"),
    [
        pytest.param({}, None, id="unpinned"),
        pytest.param({"HP": "ontologies/hp.obo"}, None, id="pinned-but-absent"),
        pytest.param({"HP": "ontologies/hp.obo"}, b"id: HP:\xff\n", id="not-utf-8"),
    ],
)
def test_an_hpo_release_that_cannot_be_read_is_none_never_an_empty_set(
    tmp_path: Path, pinned: dict[str, str], content: bytes | None
) -> None:
    """None means "the release never loaded"; an empty set means "nothing is structural".

    The second is a claim about every syndrome. Measured 2026-10-26 by handing
    `validate_orphanet_warrants` the committed scope and Orphanet mirrors: with
    None it reports one SCP006 naming the pinned HPO release and one SCP008
    skip; with `frozenset()` it reports 36 SCP006 refusals, one per
    Orphanet-scoped term, each saying the syndrome "carries no structural heart
    phenotype" -- 36 false findings about Orphanet's data and not one about
    the file that did not load. A `frozenset()` in place of the first `None`
    survived the release review's 634-test mutation set.
    """
    if content is not None:
        (tmp_path / "ontologies").mkdir()
        (tmp_path / "ontologies" / "hp.obo").write_bytes(content)
    assert pinned_structural_phenotypes(tmp_path, pinned) is None
    assert pinned_hpo(tmp_path, pinned) is None


def test_the_pinned_release_gives_the_structural_set_and_each_terms_name(
    tmp_path: Path, structural: frozenset[str]
) -> None:
    """`pinned_hpo`, which the build reads once for the scope warrants and the variant table.

    On the synthetic release: one label per `[Term]`, keyed by its id -- never by
    an alt id (HP:0001630), and never the Typedef's -- and the same structural set
    `structural_phenotypes` reads. On the pinned release, the two labels the
    variant table's tests print.
    """
    (tmp_path / "ontologies").mkdir()
    (tmp_path / "ontologies" / "hp.obo").write_text(_OBO, encoding="utf-8")
    synthetic = pinned_hpo(tmp_path, {"HP": "ontologies/hp.obo"})
    assert synthetic is not None
    assert synthetic.structural == structural
    assert synthetic.labels == {
        "HP:0000118": "Phenotypic abnormality",
        "HP:0001627": "Abnormal heart morphology",
        "HP:0001671": "Abnormal cardiac septum morphology",
        "HP:0001631": "Atrial septal defect",
        "HP:0001637": "Abnormal myocardium morphology",
        "HP:0001638": "Cardiomyopathy",
        "HP:0030962": "Abnormal morphology of the great vessels",
        "HP:0001643": "Patent ductus arteriosus",
        "HP:0001680": "Coarctation of aorta",
    }
    repo = Path(__file__).resolve().parents[2]
    registry, _ = load_sources(repo)
    pinned = pinned_hpo(repo, registry.ontology_files())
    assert pinned is not None, "the pinned HPO release did not load"
    assert pinned.labels["HP:0001631"] == "Atrial septal defect"
    assert pinned.labels["HP:0011722"] == "Mixed total anomalous pulmonary venous connection"


def test_the_structural_set_the_build_applies_is_the_one_methods_html_states() -> None:
    """291 terms over the pinned HPO release, alt ids included (measured 2026-10-26).

    An HPO upgrade that moves this changes which syndromes the rule admits
    without any line of this repository changing, so it must fail here and be
    re-measured on purpose rather than drift in. Read through
    `pinned_structural_phenotypes` and the source registry, as `build_site` and
    SCP006 read it, since 2026-10-26; it opened `hp-2026-06-23.obo` by name
    until then, so re-pinning a new release in `mirrors/sources.yaml` would have
    changed the rule while this went on measuring the old file.

    **`methods.html#scope` defines a structural heart malformation in words**
    (`methods.SCOPE_RULE`, pinned as a literal in `test_build_methods.py`), and
    every term below is a clause of that sentence, measured against this
    release: eight lesions are counted though HPO files them outside Abnormal
    heart morphology, each with every subtype filed under it; cardiomyopathy,
    tumours and the pericardium are excluded; the findings counted "still
    include some that are not malformations, such as a calcified or prolapsing
    valve"; and peripheral pulmonary artery stenosis and an aortopulmonary
    window are not counted. So widening or narrowing the rule where the
    sentence names it fails here, and the sentence is reworded on purpose rather
    than left describing a rule the build no longer applies.

    **Widened by owner decision 2026-10-26 (R1)** from 264 terms: six of the
    eight lesions were named "not counted" (three of them) or not named at all
    until then. Measured that day over the committed mirrors, the widening makes
    11 more Orphanet disorders qualify (Noonan syndrome-like disorder with loose
    anagen hair, ORPHA:2701, among them) and takes SCP007's digest from 239
    unnominated terms to 248; it admits no term, because nomination stays a
    curator's. Each lesion's subtypes are listed so that dropping any one of the
    eight from `ADDED_BRANCHES` fails here naming it.
    """
    repo = Path(__file__).resolve().parents[2]
    registry, _ = load_sources(repo)
    structural = pinned_structural_phenotypes(repo, registry.ontology_files())
    assert structural is not None, "the pinned HPO release did not load"
    assert len(structural) == 291
    assert "HP:0001631" in structural  # atrial septal defect
    # Counted, though HPO files them outside Abnormal heart morphology, each with
    # every subtype (and alt id) this release files under it.
    added = {
        "coarctation": {"HP:0001680", "HP:0004974", "HP:0005151", "HP:0011646", "HP:0011647"},
        "patent ductus": {"HP:0001643", "HP:0011648", "HP:0011649", "HP:0011650"},
        "pulmonic stenosis": {"HP:0001642", "HP:0034348", "HP:0034349", "HP:0034350"},
        "aortic valve stenosis": {"HP:0001650", "HP:0005140"},
        "supravalvular aortic stenosis": {"HP:0004381"},
        "truncus": {"HP:0001660", "HP:0004384", "HP:0011608", "HP:0011609", "HP:0011610"},
        "interrupted arch": {"HP:0011611", "HP:0006680", "HP:0011612", "HP:0011613", "HP:0011614"},
        "anomalous pulmonary venous return": {
            "HP:0010772",
            "HP:0005160",  # total
            "HP:0010773",  # partial
            "HP:0011626",  # scimitar anomaly, a partial return
            "HP:0011719",
        },
    }
    for lesion, terms in added.items():
        assert terms <= structural, f"{lesion} is no longer counted: {sorted(terms - structural)}"
    # Excluded: cardiomyopathy, a neoplasm of the heart, the pericardium.
    assert not {"HP:0001638", "HP:0100544", "HP:0001697"} & structural
    # Kept, and not malformations: aortic valve calcification, mitral valve prolapse.
    assert {"HP:0004380", "HP:0001634"} <= structural
    # Not counted: peripheral pulmonary artery stenosis, an aortopulmonary window.
    assert not {"HP:0004969", "HP:0011604"} & structural
