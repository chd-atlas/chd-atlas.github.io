# tests/unit/test_build_bundles.py
import json
from datetime import date
from pathlib import Path
from typing import Any, Final

import pytest

from chd_atlas.build.bundles import build_genes
from chd_atlas.build.burden import BurdenRow
from chd_atlas.build.concordance import cohort_families, gene_concordance
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.identity import GeneIdentity
from chd_atlas.build.omics import ModalitySummary
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.vocab import Classification, ValiditySource, ValidityState

# Two HGNC ids whose lexical order is the reverse of their symbols': "HGNC:11604"
# sorts before "HGNC:4173", while "GATA4" sorts before "TBX5". Exactly one test
# below relies on that disagreement to tell one ordering rule from the other —
# `test_the_index_is_ordered_by_hgnc_id_rather_than_by_symbol`. Every other use
# of GATA4 here is just a second gene.
TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"
SYMBOLS = {TBX5: "TBX5", GATA4: "GATA4"}


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


def _functional(**overrides: object) -> FunctionalEvidence:
    payload: dict[str, object] = {
        "id": "CHDA:FUN:0000001",
        "gene": TBX5,
        "organism": "NCBITaxon:10090",
        "perturbation": "knockout",
        "zygosity": "heterozygous",
        "cardiac_phenotype": ["MP:0000284"],
        "phenocopies_human": "partial",
        "rescue_outcome": "not_attempted",
        "publication": "PMID:8988165",
    }
    payload.update(overrides)
    return FunctionalEvidence.model_validate(payload)


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {"root": Path("."), "assertions": (_assertion(),)}
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def _summary(**overrides: object) -> ModalitySummary:
    payload: dict[str, Any] = {
        "count": 3,
        "shards": ["omics/expression/GSE1000.json"],
        "top": [{"gene": TBX5, "fdr": 0.01}],
    }
    payload.update(overrides)
    return ModalitySummary(**payload)


def _variant(**overrides: object) -> dict[str, Any]:
    payload: dict[str, Any] = {"vrs_id": "ga4gh:VA.x", "gene": TBX5, "chrom": "12", "pos": 100}
    payload.update(overrides)
    return payload


def _validity_record(**overrides: object) -> ValidityRecord:
    payload: dict[str, object] = {
        "source": ValiditySource.CLINGEN,
        "classification_term": "Definitive",
        "classification": Classification.DEFINITIVE,
        "disease": "MONDO:0007732",
        "disease_label": "Holt-Oram syndrome",
        "moi": "AD",
        "report_url": None,
    }
    payload.update(overrides)
    return ValidityRecord(**payload)  # type: ignore[arg-type]


def _gene_validity(**overrides: object) -> GeneValidity:
    payload: dict[str, object] = {
        "state": ValidityState.EXPERT_CURATED,
        "records": (_validity_record(),),
        "has_source_discordance": False,
    }
    payload.update(overrides)
    return GeneValidity(**payload)  # type: ignore[arg-type]


_EMPTY_CONCORDANCE: Final[dict[str, Json]] = {
    "tested": 0,
    "enriched": 0,
    "corrected": 0,
    "families": [],
}


def _no_burden() -> dict[str, Json]:
    """The concordance mapping for a corpus where no study reported anything.

    `build_genes` takes it as a required argument rather than deriving it, so
    the runner's single derivation reaches the bundles, the gene pages and the
    browse page alike. A default would let a caller forget it and publish zeros
    indistinguishable from a measured "no dataset tested this" -- this project's
    characteristic failure -- so the tests supply it explicitly.

    A real `dict` over the genes these fixtures publish, **not a
    `defaultdict`**. `_concordance_for` refuses a gene the mapping does not
    contain, and a `defaultdict` reports `__contains__` as `False` for exactly
    the keys it would happily invent -- so it would trip the guard while looking
    like it satisfied it. That is the guard working, and it is why this returns
    the population rather than a fallback.
    """
    return {gene: dict(_EMPTY_CONCORDANCE) for gene in SYMBOLS}


def _placement(**overrides: object) -> dict[str, Any]:
    """A `build.profiles.Placement`, hand-built the way `_summary` stands in
    for a `ModalitySummary` -- `build_genes` treats `profiles` as an
    already-assembled input, the same dependency-injection precedent `omics`
    already follows in this file, so the assembly itself
    (`build.profiles.gene_expression_profiles`) is exercised in
    `test_build_profiles.py`, not here.
    """
    payload: dict[str, Any] = {
        "q25_percentile": 40,
        "median_percentile": 50,
        "q75_percentile": 60,
        "median_abundance": 12.5,
        "unit": "rpkm",
        "n_samples": 3,
        "n_genes": 19842,
        "method": "lowest percentile of a tied breakpoint (bisect_left)",
    }
    payload.update(overrides)
    return payload


def _specificity_payload(**overrides: object) -> dict[str, Any]:
    """A `build.profiles.Specificity`, hand-built -- see `_placement`."""
    payload: dict[str, Any] = {
        "tau": 0.48,
        "scale": "log2(x+1)",
        "method": "tau (Yanai et al. 2005)",
        "tissues": ["Brain", "Heart", "Liver"],
        "n_tissues": 3,
        "highest_in": "Heart",
        "medians": {"Brain": 10.0, "Heart": 100.0, "Liver": 10.0},
    }
    payload.update(overrides)
    return payload


def _profile_payload(**overrides: object) -> dict[str, Any]:
    """A `build.profiles.ExpressionProfile`, hand-built -- see `_placement`."""
    payload: dict[str, Any] = {
        "datasets": [
            {
                "dataset": "E-MTAB-6814",
                "quantile_shard": "omics/profile_quantiles/E-MTAB-6814.json",
                "stages": [
                    {
                        "stage": "7wpc",
                        "phase": {
                            "outcome": "matched",
                            "phase_ids": ["septation"],
                            "reason": None,
                        },
                        "specificity": _specificity_payload(),
                        "specificity_unavailable_reason": None,
                        "tissues": [
                            {
                                "tissue": "Heart",
                                "median_abundance": 100.0,
                                "unit": "rpkm",
                                "n_samples": 3,
                                "placement": _placement(median_abundance=100.0),
                                "not_placed_reason": None,
                            }
                        ],
                    }
                ],
            }
        ]
    }
    payload.update(overrides)
    return payload


def _read(root: Path, relative: str) -> dict[str, Any]:
    """One emitted artifact, read back from the file that was written."""
    payload: dict[str, Any] = json.loads((root / relative).read_text())
    return payload


def _entries(root: Path) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = _read(root, "genes/index.json")["genes"]
    return entries


def _by_gene(root: Path) -> dict[str, dict[str, Any]]:
    return {str(entry["gene"]): entry for entry in _entries(root)}


def test_an_index_row_is_exactly_what_the_browser_filters_on(tmp_path: Path) -> None:
    """The index is downloaded by every visitor before they pick a gene.

    Its size is the cost of the browse page, so what belongs in a row is what
    ranks or filters it plus the path to fetch the rest. An assertion, its
    evidence, a variant list or a publication list added here would be
    downloaded by everyone to be read by nobody — the bundle already carries all
    four.

    Asserted as one whole-row equality rather than as a key set and a list of
    field values: equality fails on a key added, a key removed and any wrong
    value alike, so it is the two checks in one and cannot drift apart from
    them. `bundle` is pinned as a literal string for the reason the key exists —
    the client must never build a bundle path itself, so this exact value is the
    contract, not just its shape.
    """
    emitter = Emitter(root=tmp_path)
    validity = {TBX5: _gene_validity()}

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _entries(tmp_path) == [
        {
            "gene": TBX5,
            "symbol": "TBX5",
            "headline_confidence": "definitive",
            "validity_state": "expert_curated",
            "atlas_curation": "curated",
            "has_conflicting_evidence": False,
            "has_source_discordance": False,
            "has_no_association_report": False,
            "no_association_reported_by": [],
            "lesion_groups": ["septal"],
            "confidence_by_lesion_group": {"septal": "definitive"},
            "conflicting_lesion_groups": [],
            "evidence_counts": {"genetic_case": 1},
            "assertion_count": 1,
            "functional_count": 0,
            "variant_count": 0,
            "burden_row_count": 0,
            # Every key present even where no study reported the gene, so a
            # consumer reads one shape rather than guarding for a missing one.
            # `families` is empty rather than absent for the same reason
            # `burden` is an empty array: "no study reported this" must not be
            # indistinguishable from "the build dropped it".
            "independent_datasets": {
                "tested": 0,
                "enriched": 0,
                "corrected": 0,
                "families": [],
            },
            "bundle": "genes/HGNC_11604.json",
        }
    ]
    # `is False`, not the `== False` the equality above settles for: `0 == False`
    # in Python, so a flag written as an int would publish `0` and slip through
    # that comparison while a consumer testing `=== false` reads it as unset.
    assert _entries(tmp_path)[0]["has_conflicting_evidence"] is False
    assert _entries(tmp_path)[0]["has_source_discordance"] is False


def test_a_published_gene_the_atlas_has_not_curated_gets_a_row_and_a_bundle(
    tmp_path: Path,
) -> None:
    """The other 91 of the 92 genes published today, end to end.

    D21 publishes a gene on an upstream authority's call -- a ClinGen record at
    or above `Limited`, or two agreeing GenCC submitters -- which is not the
    atlas's own work; `atlas_curation` is the one key a consumer filtering a
    browse list reads to tell the two apart. Both genes here are published and
    only TBX5 is asserted, so a `_headline` that hardcoded either value — or
    that dropped the key and let a reader's `.get(..., "curated")` fill it in —
    fails on one of the two rows rather than passing on a fixture that only ever
    shows one answer.

    The uncurated gene's whole payload is checked, not only its flag: it has to
    carry a real headline (`definitive`, from the mirror, which is the entire
    reason it is published at all) and *empty* evidence containers, because a
    page that rendered a missing key as a broken section is the failure this
    replaces. The bundle is asserted to exist because the index row advertises
    it; a row without a file behind it is the dead link nothing downstream
    catches.
    """
    emitter = Emitter(root=tmp_path)
    validity = {TBX5: _gene_validity(), GATA4: _gene_validity()}

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    rows = _by_gene(tmp_path)
    assert rows[TBX5]["atlas_curation"] == "curated"
    assert rows[GATA4]["atlas_curation"] == "not_yet_curated"
    assert rows[GATA4]["assertion_count"] == 0
    assert rows[GATA4]["headline_confidence"] == "definitive"
    assert rows[GATA4]["lesion_groups"] == []
    assert rows[GATA4]["confidence_by_lesion_group"] == {}
    assert rows[GATA4]["evidence_counts"] == {}

    bundle = _read(tmp_path, str(rows[GATA4]["bundle"]))
    assert bundle["atlas_curation"] == "not_yet_curated"
    assert bundle["assertions"] == []
    assert bundle["functional"] == []
    assert bundle["publications"] == []
    assert bundle["validity"]["records"], "the panel call it is published on must be on the page"


def test_a_bundle_carries_the_whole_gene_page_and_nothing_more(tmp_path: Path) -> None:
    """A gene detail page is one fetch, so every section it renders is here.

    `admitted_by` and `asserted_by` joined the set in schema 2.8. They are the
    only fields here that answer "why is this gene on the site at all", which
    stopped being a question with one answer when the gate widened past ClinGen
    `Definitive` -- 16 of the 92 genes published are admitted by GenCC submitter
    agreement and carry no ClinGen record at all.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert set(_read(tmp_path, "genes/HGNC_11604.json")) == {
        "gene",
        "symbol",
        "identity",
        "admitted_by",
        "asserted_by",
        "headline_confidence",
        "validity_state",
        "atlas_curation",
        "has_conflicting_evidence",
        "has_source_discordance",
        "has_no_association_report",
        "no_association_reported_by",
        "independent_datasets",
        "validity",
        "lesion_groups",
        "publications",
        "assertions",
        "functional",
        "variants",
        "omics",
        "burden",
        "burden_combined",
        "expression_profile",
        "cell_type_expression",
    }
    # Always present, empty when the gene has none -- an absent key would make
    # "no study combined its arms" indistinguishable from "the build dropped it".
    assert _read(tmp_path, "genes/HGNC_11604.json")["burden_combined"] == []


def test_a_gene_with_no_mirrored_validity_publishes_the_uncurated_shape(
    tmp_path: Path,
) -> None:
    """`validity={}` is the same fallback `gene_facts` documents for this gene.

    The bundle's `validity` block must read exactly what `build.validity.
    uncurated()` returns -- an explicit `"uncurated"`/empty-records shape --
    rather than an absent key or some other guess at what "nothing curated
    this" should look like.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _read(tmp_path, "genes/HGNC_11604.json")["validity"] == {
        "state": "uncurated",
        "has_source_discordance": False,
        "records": [],
    }


def test_every_validity_record_carries_the_same_key_set_regardless_of_source(
    tmp_path: Path,
) -> None:
    """A ClinGen record and a GenCC record publish identical keys.

    `sop`, `classification_date` and `gcep` are ClinGen-only and `submitter`,
    `submitted_on` and `sgc_id` are GenCC-only on `ValidityRecord` (see its
    docstring); the field a record's own source does not carry must still
    appear, as `null`, rather than being missing from that record's object.
    Publishing a mix of shapes is a trap for a consumer that reads
    `record["sop"]` off every entry of the array — it would have to check
    `source` first, which the published contract does not require.
    """
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN,
                    sop="SOP11",
                    gcep="G",
                    classification_date="2025-03-25T16:00:00.000Z",
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.MODERATE,
                    submitter="Ambry Genetics",
                    submitted_on="2018-05-31 20:04:32",
                    sgc_id="SGC-102815",
                ),
            )
        )
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    records = _read(tmp_path, "genes/HGNC_11604.json")["validity"]["records"]
    clingen = next(r for r in records if r["source"] == "clingen")
    gencc = next(r for r in records if r["source"] == "gencc")
    # A literal, not `set(clingen) == set(gencc)` alone: two records that both
    # dropped `submitted_on` would still agree with each other (section 4.38).
    assert (
        set(clingen)
        == set(gencc)
        == {
            "source",
            "classification_term",
            "classification",
            "disease",
            "disease_label",
            "moi",
            "sop",
            "classification_date",
            "gcep",
            "report_url",
            "submitter",
            "submitted_on",
            "sgc_id",
        }
    )
    assert gencc["submitted_on"] == "2018-05-31 20:04:32"
    assert gencc["sgc_id"] == "SGC-102815"
    assert clingen["submitted_on"] is None
    assert clingen["sgc_id"] is None
    assert clingen["sop"] == "SOP11"
    assert clingen["submitter"] is None
    assert clingen["classification_date"] == "2025-03-25T16:00:00.000Z"
    assert gencc["classification_date"] is None
    assert gencc["gcep"] is None
    assert gencc["submitter"] == "Ambry Genetics"
    assert gencc["sop"] is None


def test_the_bundle_does_not_resort_the_records_validity_py_already_ordered(
    tmp_path: Path,
) -> None:
    """`records` is published in whatever order `GeneValidity.records` holds.

    Handed GenCC-then-ClinGen here — the reverse of the `(source, disease, moi,
    submitter, sgc_id)` order `build.validity.gene_validity` itself sorts by —
    so a `bundles.py` that quietly re-sorted, or rebuilt the tuple through a
    `dict` or a `set` on the way to JSON, would publish ClinGen first and this
    would fail. Not re-sorting is the only way it passes.
    """
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(source=ValiditySource.GENCC, submitter="Zeta"),
                _validity_record(source=ValiditySource.CLINGEN),
            )
        )
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    records = _read(tmp_path, "genes/HGNC_11604.json")["validity"]["records"]
    assert [r["source"] for r in records] == ["gencc", "clingen"]


def test_every_bundle_path_the_index_advertises_was_written(tmp_path: Path) -> None:
    """Checked over the whole index rather than on one hardcoded name.

    `bundle` is the URL a detail page is fetched from, so an index entry whose
    path was never written is a dead link that nothing else in the build
    notices: `emitter.checksums` records what was written, never what a payload
    promised. The two come from one `gene_bundle_path` call per gene, and this
    is what holds them together if that stops being true.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    advertised = [str(entry["bundle"]) for entry in _entries(tmp_path)]
    assert len(advertised) == 2
    for relative in advertised:
        assert relative in emitter.checksums, f"{relative} is in the index but was never written"
        assert (tmp_path / relative).is_file(), f"{relative} is in the manifest but not on disk"


def test_a_contested_gene_is_flagged_in_both_the_index_and_the_bundle(
    tmp_path: Path,
) -> None:
    """Spec 5.2: a contested gene must never be displayed as settled.

    `strongest()` ranks on one linear scale where definitive outranks refuted,
    so a gene carrying both resolves to definitive and the refutation is
    invisible in `headline_confidence` alone. The flag is the other half of that
    pair, and this module is its first consumer.

    Both genes are asserted here, contested and uncontested, with the contest
    carried by their *mirrored* records rather than by anything the curated
    assertion says. A writer that hardcoded either value, or that dropped the
    key and let a reader's `.get(..., False)` fill it in, fails on one of the
    assertions below rather than passing on a fixture that only ever shows one
    answer. The `headline_confidence` assertions are what make the point: TBX5
    still reads `definitive`, so the flag is carrying information that the
    confidence cannot — and GATA4 reads `moderate`, so a headline hardcoded to
    the value a contested gene happens to resolve to fails here too.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
        )
    )
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN, classification=Classification.DEFINITIVE
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.REFUTED,
                    submitter="s",
                ),
            )
        ),
        GATA4: _gene_validity(records=(_validity_record(classification=Classification.MODERATE),)),
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    entries = _by_gene(tmp_path)
    assert entries[TBX5]["headline_confidence"] == "definitive"
    assert entries[TBX5]["has_conflicting_evidence"] is True
    assert entries[GATA4]["headline_confidence"] == "moderate"
    assert entries[GATA4]["has_conflicting_evidence"] is False

    contested = _read(tmp_path, "genes/HGNC_11604.json")
    settled = _read(tmp_path, "genes/HGNC_4173.json")
    assert contested["headline_confidence"] == "definitive"
    assert contested["has_conflicting_evidence"] is True
    assert settled["headline_confidence"] == "moderate"
    assert settled["has_conflicting_evidence"] is False


def test_confidence_is_broken_down_by_lesion_group(tmp_path: Path) -> None:
    """Every declared group carries the breakdown, keyed the same way the headline is.

    ClinGen and GenCC classify a gene against a disease, never against a lesion,
    so the mirrored classifications carry no per-group signal of their own.
    TBX5 here declares two lesion groups across two assertions; both publish the
    *same* collapsed value — `strongest()` of the gene's mirrored records — which
    is `definitive` even though one of the two mirrored records is `refuted`.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", lesion_groups=["septal"]),
            _assertion(id="CHDA:AST:0000002", lesion_groups=["conotruncal"]),
        )
    )
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN, classification=Classification.DEFINITIVE
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.REFUTED,
                    submitter="s",
                ),
            )
        )
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    entry = _entries(tmp_path)[0]
    assert entry["confidence_by_lesion_group"] == {
        "conotruncal": "definitive",
        "septal": "definitive",
    }
    # Sorted by value, and the same set of groups the breakdown keys on.
    assert entry["lesion_groups"] == ["conotruncal", "septal"]
    assert _read(tmp_path, "genes/HGNC_11604.json")["lesion_groups"] == ["conotruncal", "septal"]


def test_evidence_counts_are_carried_per_evidence_class(tmp_path: Path) -> None:
    """Two classes with different counts, so a total or a length cannot pass."""
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(evidence_class="genetic_case"),
                    _evidence(evidence_class="genetic_case", publication="PMID:9"),
                    _evidence(
                        evidence_class="functional_model",
                        functional_evidence="CHDA:FUN:0000001",
                    ),
                ]
            ),
        ),
        functional=(_functional(),),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _entries(tmp_path)[0]["evidence_counts"] == {"functional_model": 1, "genetic_case": 2}


def test_the_browse_counts_match_the_bundle_they_link_to(tmp_path: Path) -> None:
    """A count smaller than the page it links to reads as missing data.

    The three counts are the browse layer's promise about what opening a gene
    will show, and each is derived separately from the list the bundle carries —
    `assertion_count` and `functional_count` in `derive.py`, `variant_count`
    here. The two genes differ on all three, so a count wired to the wrong
    collection, or to a constant, fails rather than coinciding.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=TBX5),
            _assertion(id="CHDA:AST:0000003", gene=GATA4),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=TBX5),
            _functional(id="CHDA:FUN:0000003", gene=TBX5),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={
            TBX5: [_variant(vrs_id="ga4gh:VA.1"), _variant(vrs_id="ga4gh:VA.2")],
            GATA4: [_variant(vrs_id="ga4gh:VA.3", gene=GATA4)],
        },
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    counted = {
        str(entry["gene"]): (
            entry["assertion_count"],
            entry["functional_count"],
            entry["variant_count"],
        )
        for entry in _entries(tmp_path)
    }
    assert counted == {TBX5: (2, 3, 2), GATA4: (1, 0, 1)}
    for entry in _entries(tmp_path):
        bundle = _read(tmp_path, str(entry["bundle"]))
        assert entry["assertion_count"] == len(bundle["assertions"])
        assert entry["functional_count"] == len(bundle["functional"])
        assert entry["variant_count"] == len(bundle["variants"])


def test_a_gene_absent_from_the_registry_keeps_its_hgnc_id_as_its_label(
    tmp_path: Path,
) -> None:
    """`mirrors/genes.tsv` is a mirror, so a gene can be asserted before it lands.

    The fallback is the id rather than the empty string because the browse row
    and the page heading are rendered from this value: blank reads as a bug in
    the site, while `HGNC:4173` reads as a gene whose symbol is not in yet and
    is still enough to search for. TBX5 is in the registry in the same call, so
    a writer that always falls back fails here too.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols={TBX5: "TBX5"},
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    entries = _by_gene(tmp_path)
    assert entries[TBX5]["symbol"] == "TBX5"
    assert entries[GATA4]["symbol"] == GATA4
    assert _read(tmp_path, "genes/HGNC_4173.json")["symbol"] == GATA4


def test_a_bundle_carries_its_assertions_in_full(tmp_path: Path) -> None:
    """Including the evidence, which is the record a curator is judged on."""
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assertion = bundle["assertions"][0]
    assert assertion["id"] == "CHDA:AST:0000001"
    assert assertion["evidence"][0]["summary"] == "s"
    assert assertion["evidence"][0]["locator"] == {
        "kind": "supplementary",
        "file": "T.xlsx",
        "sheet": None,
        "row": None,
    }
    # Dumped in JSON mode, so a date is a string rather than a `date` that
    # `json.dumps` refuses part-way through a build.
    assert assertion["curated_on"] == "2026-07-01"
    assert "classification" not in assertion, (
        "the curated record no longer carries a validity call (D12) -- that "
        "field is now published only from the mirror, as headline_confidence"
    )


def test_a_bundle_carries_every_functional_record_about_the_gene(tmp_path: Path) -> None:
    """Every record about the gene, not only those an assertion cites.

    `referential.py` requires a *cited* functional record to be about the
    asserting gene, but never requires a curated record to be cited at all, so
    the two populations differ. `derive.functional_count` counts this same
    population, and the browse count would otherwise exceed the page.
    """
    corpus = _corpus(
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=TBX5),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [record["id"] for record in bundle["functional"]] == [
        "CHDA:FUN:0000001",
        "CHDA:FUN:0000002",
    ]
    assert bundle["functional"][0]["organism"] == "NCBITaxon:10090"


def test_a_bundle_carries_the_omics_summaries_verbatim(tmp_path: Path) -> None:
    """The whole `ModalitySummary`, not a re-listing of the fields known today.

    The summary here carries a fourth key that `ModalitySummary` does not
    declare, standing in for the field someone adds to it next. With only the
    three declared fields, rebuilding the summary key by key produces byte-identical
    output and this assertion cannot tell the two implementations apart; the
    undeclared key is what makes whole-value equality a real guard rather than a
    stated intention. `_summaries` passes the value through, so it survives; a
    rebuild drops it silently, and the symptom would be a gene page missing data
    the shard it links to still holds.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={TBX5: {"expression": _summary(median_log2fc=1.5)}},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _read(tmp_path, "genes/HGNC_11604.json")["omics"] == {
        "expression": {
            "count": 3,
            "shards": ["omics/expression/GSE1000.json"],
            "top": [{"gene": TBX5, "fdr": 0.01}],
            "median_log2fc": 1.5,
        }
    }


def test_a_bundle_embeds_its_variants_rather_than_linking_them(tmp_path: Path) -> None:
    """Affordable because curation bounds the count, not because the data does.

    The order is the one `build_variants` produced — karyotype order across
    shards, then each shard's own validated sort order — and nothing here
    reorders it.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={TBX5: [_variant(vrs_id="ga4gh:VA.2"), _variant(vrs_id="ga4gh:VA.1")]},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [row["vrs_id"] for row in bundle["variants"]] == ["ga4gh:VA.2", "ga4gh:VA.1"]
    assert bundle["variants"][0]["pos"] == 100


def test_a_gene_with_no_omics_or_variants_gets_empty_containers(tmp_path: Path) -> None:
    """A consumer should read `bundle["omics"]`, never guard for its absence."""
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert bundle["omics"] == {}
    assert bundle["variants"] == []


def test_a_gene_with_no_profile_data_still_carries_the_key(tmp_path: Path) -> None:
    """The rule `bundles.py` already states for `omics` and `burden`.

    An absent key would make "no study reported this" indistinguishable from
    "the build dropped it".
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert "expression_profile" in bundle
    assert bundle["expression_profile"] == {"datasets": []}


def test_every_bundle_carries_cell_type_expression_even_when_empty(tmp_path: Path) -> None:
    """The rule `omics` and `expression_profile` keep: always present.

    A gene absent from the mapping publishes exactly `{"datasets": []}`; one
    covered gene is carried verbatim. The mapping itself has no default since
    review (2026-10-28), so a build handed none fails rather than publishing
    the empty shape on every bundle.
    """
    emitter = Emitter(root=tmp_path)
    covered: Any = {"datasets": [{"dataset": "hoc", "unit": "cpm"}]}

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={TBX5: covered},
    )

    assert _read(tmp_path, "genes/HGNC_11604.json")["cell_type_expression"] == covered
    assert _read(tmp_path, "genes/HGNC_4173.json")["cell_type_expression"] == {"datasets": []}


def test_the_bundle_carries_the_medians_tau_was_computed_from(tmp_path: Path) -> None:
    """D39(b): tau's inputs must be reachable from the payload publishing tau.

    The omics slice is truncated (`omics.select_top`), so the per-organ
    values cannot be recovered from `top` in general.
    """
    emitter = Emitter(root=tmp_path)
    profile = _profile_payload()

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={TBX5: profile},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    specificity = bundle["expression_profile"]["datasets"][0]["stages"][0]["specificity"]
    assert set(specificity["medians"]) == set(specificity["tissues"])
    assert specificity["n_tissues"] == len(specificity["medians"])
    # The whole payload is carried verbatim, not rebuilt field by field --
    # the same guarantee `test_a_bundle_carries_the_omics_summaries_verbatim`
    # pins for `omics`, and for the same reason: a rebuild drops a field
    # silently, while a passthrough fails loudly the moment the shapes stop
    # matching.
    assert bundle["expression_profile"] == profile


def test_a_merely_registered_genes_profile_facts_never_reach_a_bundle(tmp_path: Path) -> None:
    """`profiles` may carry more genes than `published` -- the same
    registry-vs-published asymmetry `mirrors/genes.tsv` already has (154
    registered, 92 published). `gene_expression_profiles` takes no
    `published` argument at all, so filtering happens exactly once, here, by
    the same iteration over `published` that already restricts `omics` and
    `variants`.
    """
    emitter = Emitter(root=tmp_path)
    profile = _profile_payload()

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},  # GATA4 is in `profiles` below but NOT published
        burden={},
        concordance=_no_burden(),
        profiles={TBX5: profile, GATA4: _profile_payload()},
        cell_types={},
    )

    assert not (tmp_path / "genes" / "HGNC_4173.json").exists()
    assert _read(tmp_path, "genes/HGNC_11604.json")["expression_profile"] == profile


def test_a_bundle_holds_only_its_own_genes_evidence(tmp_path: Path) -> None:
    """Every collection in the bundle is per-gene, and each is filtered separately.

    Both genes carry all four kinds of evidence here, so a collection joined on
    the wrong key — or not joined at all — shows up as another gene's records in
    this bundle rather than as an empty one.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(
                id="CHDA:AST:0000002", gene=GATA4, evidence=[_evidence(publication="PMID:9")]
            ),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=GATA4),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={TBX5: {"expression": _summary()}, GATA4: {"profiles": _summary(count=9)}},
        variants={TBX5: [_variant()], GATA4: [_variant(gene=GATA4, vrs_id="ga4gh:VA.y")]},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_4173.json")
    assert bundle["gene"] == GATA4
    assert [item["id"] for item in bundle["assertions"]] == ["CHDA:AST:0000002"]
    assert [item["id"] for item in bundle["functional"]] == ["CHDA:FUN:0000002"]
    assert [row["vrs_id"] for row in bundle["variants"]] == ["ga4gh:VA.y"]
    assert list(bundle["omics"]) == ["profiles"]
    assert bundle["publications"] == ["PMID:9"]


def test_the_index_is_ordered_by_hgnc_id_rather_than_by_symbol(tmp_path: Path) -> None:
    """The index is a JSON array, so its order is published and must be fixed.

    `encode_json` sorts dict keys and never array elements, so whatever order
    this list reaches the emitter in is the order every consumer downloads. The
    two genes are given in the reverse of both candidate orders, and their
    symbols rank them the other way round from their ids, so a build that
    emitted encounter order or sorted on the symbol fails here.

    This pins the published contract end to end. It cannot attribute the sort:
    `gene_facts` already returns its genes sorted by id, so `bundles.py`
    sorting them again is defence against that changing rather than what
    produces the order today.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert [entry["gene"] for entry in _entries(tmp_path)] == [TBX5, GATA4]
    assert [entry["symbol"] for entry in _entries(tmp_path)] == ["TBX5", "GATA4"]


def test_bundle_assertions_and_functional_records_are_ordered_by_id(tmp_path: Path) -> None:
    """Both are JSON arrays whose order survives `encode_json` untouched.

    Given as ids 3, 1, 2 — deliberately neither id order nor its reverse. A
    build that kept corpus order fails on every run rather than on the fraction
    a hash-ordered container would, and so does one that reversed corpus order
    instead of sorting, which a 3, 2, 1 fixture would hand the right answer to.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000003"),
            _assertion(id="CHDA:AST:0000001"),
            _assertion(id="CHDA:AST:0000002"),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000003"),
            _functional(id="CHDA:FUN:0000001"),
            _functional(id="CHDA:FUN:0000002"),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [item["id"] for item in bundle["assertions"]] == [
        "CHDA:AST:0000001",
        "CHDA:AST:0000002",
        "CHDA:AST:0000003",
    ]
    assert [item["id"] for item in bundle["functional"]] == [
        "CHDA:FUN:0000001",
        "CHDA:FUN:0000002",
        "CHDA:FUN:0000003",
    ]


def test_a_bundle_lists_the_publications_its_evidence_cites(tmp_path: Path) -> None:
    """Sorted lexically, matching `publications.json` rather than contradicting it.

    Six of them because `derive.py` accumulates them through a `set`, whose
    iteration order for strings varies with `PYTHONHASHSEED`: at three PMIDs an
    unsorted tuple still comes out sorted on roughly one seed in four. The
    ordering rule is `derive.py`'s; what this pins is that the bundle publishes
    it rather than re-deriving one of its own.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(publication="PMID:9"),
                    _evidence(publication="PMID:8988165"),
                    _evidence(publication="PMID:11729"),
                    _evidence(publication="PMID:3"),
                    _evidence(publication="PMID:40404"),
                    _evidence(publication="PMID:777"),
                ]
            ),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _read(tmp_path, "genes/HGNC_11604.json")["publications"] == [
        "PMID:11729",
        "PMID:3",
        "PMID:40404",
        "PMID:777",
        "PMID:8988165",
        "PMID:9",
    ]


def test_a_gene_with_evidence_but_outside_the_published_set_is_not_published(
    tmp_path: Path,
) -> None:
    """No expert panel calls it definitive, so it gets no bundle and no index row.

    `gene_facts` omits such a gene deliberately — D21 makes `published` the
    population and nothing else admits a gene to it — and this module publishes
    what that returns. Nothing here raises over it: the corpus below is legal
    under every validator in the project, and which genes exist is `gene_facts`'
    decision, not this module's.

    What the rule costs differs by evidence kind, and the asymmetry is the
    reason this test carries all three. GATA4's omics rows and variants are
    still published, in the omics shards and in `variants/<chrom>.json.gz`; they
    are merely unreachable from the gene index, so getting to them means already
    knowing the HGNC id. Its functional record is not published anywhere at all
    — no other build module reads `corpus.functional`, `derive.py` reads it only
    to count, and `build_genes` writes it only inside the loop over published
    genes. So the second half of this test is not about a missing link but about
    curated evidence that reaches no file: the record's id appears in none of
    the bytes the build wrote.

    Widening the gate to 23 genes did not close this. It moved which genes fall
    outside it — the cost now lands on a gene no panel has graded definitive
    rather than on one nobody has asserted — and D37 is the decision that keeps
    such candidate genes out of the definitive set rather than beside it.

    Pinned here so that either half changing is a decision rather than a side
    effect. The obvious one-line fix — keying the loop on the union of
    `published` and the functionally-evidenced genes — is what the index
    assertion below refuses; before this test existed, that mutation passed the
    whole suite.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(functional=(_functional(id="CHDA:FUN:0000009", gene=GATA4),)),
        emitter,
        symbols=SYMBOLS,
        omics={GATA4: {"expression": _summary()}},
        variants={GATA4: [_variant(gene=GATA4)]},
        validity={},
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert [entry["gene"] for entry in _entries(tmp_path)] == [TBX5]
    assert set(emitter.checksums) == {"genes/index.json", "genes/HGNC_11604.json"}
    written = "".join((tmp_path / relative).read_text() for relative in emitter.checksums)
    assert "CHDA:FUN:0000009" not in written
    assert GATA4 not in written


def test_a_bundle_that_cannot_be_written_leaves_no_index_at_all(tmp_path: Path) -> None:
    """The index is written last, so it never advertises a file nothing wrote.

    `Emitter` refuses to write one path twice — two genes normalising to one
    filename is a curation bug, not a build that should quietly publish one page
    where two were expected — and the refusal is what a pre-filled `checksums`
    reproduces here without needing that exotic pair of identifiers.

    Emitting the index first would leave `genes/index.json` on disk pointing at
    a bundle the build then failed to produce. The build fails either way; what
    this pins is which artifacts exist when it does. `Emitter` writes each file
    as it reaches it, so everything written before the failure stays where it
    landed.
    """
    emitter = Emitter(root=tmp_path, checksums={"genes/HGNC_11604.json": "sha256:0"})

    with pytest.raises(ValueError, match="written twice"):
        build_genes(
            _corpus(),
            emitter,
            symbols=SYMBOLS,
            omics={},
            variants={},
            validity={},
            published={TBX5},
            burden={},
            concordance=_no_burden(),
            profiles={},
            cell_types={},
        )

    assert not (tmp_path / "genes" / "index.json").exists()
    assert "genes/index.json" not in emitter.checksums


def test_the_index_is_emitted_for_an_empty_corpus(tmp_path: Path) -> None:
    """Emitted with nothing to list, rather than left absent.

    A consumer that has to tell "no genes curated yet" from "wrong URL" by
    reading a 404 will get it wrong, and the browse page fetches this before it
    can know either.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        Corpus(root=Path(".")),
        emitter,
        symbols={},
        omics={},
        variants={},
        validity={},
        published=set(),
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    assert _read(tmp_path, "genes/index.json") == {"genes": []}
    assert set(emitter.checksums) == {"genes/index.json"}


def test_a_contested_gene_names_every_declared_lesion_group_as_conflicting(
    tmp_path: Path,
) -> None:
    """The gene-level pairing, one level down — issue #4.

    `confidence_by_lesion_group` collapses each group with `strongest()` of the
    gene's *mirrored* records, so a gene carrying both a definitive and a
    refuted record publishes `definitive` for every group and the refutation
    vanishes exactly as it would from the headline. `has_conflicting_evidence`
    is computed gene-wide and `conflicting_lesion_groups` is its per-group
    counterpart — but the mirrors classify per gene-disease, never per lesion,
    so there is no finer signal left to divide the groups with: a contested
    gene names *every* group it declares, never a subset.

    TBX5 here declares two lesion groups across two assertions, septal and
    conotruncal, while its mirrored records are definitive and refuted. A build
    that named only one of the two groups, or that let the gene-level flag
    stand in without naming any group, fails here.

    The collapsed values are asserted alongside, because the list is only
    meaningful as a qualifier on them: both groups must still read
    `definitive`, or the list is describing something the map no longer says.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", lesion_groups=["septal"]),
            _assertion(id="CHDA:AST:0000002", lesion_groups=["conotruncal"]),
        )
    )
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN, classification=Classification.DEFINITIVE
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.REFUTED,
                    submitter="s",
                ),
            )
        )
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    entry = _entries(tmp_path)[0]
    assert entry["has_conflicting_evidence"] is True, "the gene is contested"
    assert entry["conflicting_lesion_groups"] == ["conotruncal", "septal"], (
        "every declared group, not a subset of them"
    )
    assert entry["confidence_by_lesion_group"] == {
        "conotruncal": "definitive",
        "septal": "definitive",
    }


def test_the_conflicting_groups_are_a_subset_of_the_groups_that_carry_confidence(
    tmp_path: Path,
) -> None:
    """A list naming a group the map does not carry would describe nothing.

    Both are built from the same `ordered_groups`, so this cannot diverge today;
    it is pinned because the list is published beside the map and a reader joins
    the two by key. An uncontested gene must publish an empty list rather than
    omit the field, so a consumer reads one shape.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5, lesion_groups=["septal"]),
            _assertion(id="CHDA:AST:0000002", gene=GATA4, lesion_groups=["septal"]),
        )
    )
    validity = {
        TBX5: _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN, classification=Classification.REFUTED
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.MODERATE,
                    submitter="s",
                ),
            )
        ),
        GATA4: _gene_validity(
            records=(_validity_record(classification=Classification.DEFINITIVE),)
        ),
    }
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity=validity,
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    for entry in _entries(tmp_path):
        named = set(entry["conflicting_lesion_groups"])
        assert named <= set(entry["confidence_by_lesion_group"])
    rows = _by_gene(tmp_path)
    assert rows[TBX5]["conflicting_lesion_groups"] == ["septal"]
    assert rows[GATA4]["conflicting_lesion_groups"] == []


def test_the_bundle_carries_a_genes_burden_rows_and_the_browse_row_counts_them(
    tmp_path: Path,
) -> None:
    """Both halves, because neither is guarded by the other.

    Measured 2026-08-05: with `burden={}` at every call site in this file, a
    mutant emitting `"burden": []` unconditionally and a mutant emitting
    `"burden_row_count": 0` unconditionally both survived the whole suite. An
    empty fixture cannot tell a builder that publishes from one that drops.

    The count is asserted against the same list the bundle embeds. A browse row
    promising nine rows of burden evidence on a page carrying none is the browse
    layer lying about the page -- the reason `variant_count` is derived the same
    way one line above it.
    """
    emitter = Emitter(root=tmp_path)
    rows = [
        BurdenRow(
            study="PMID:42230622",
            gene=TBX5,
            cohort_stratum=stratum,
            lesion_group=None,
            variant_class="snv_indel",
            consequence_class="lof",
            origin="any",
            maf_max=0.001,
            count_unit="individuals",
            n_case_carriers=5,
            n_cases=1471,
            comparator="control_cohort",
            n_control_carriers=0,
            n_controls=45082,
            expected_count=None,
            effect=None,
            effect_measure="odds_ratio",
            effect_bound="unbounded_above",
            ci_low=28.1,
            ci_high=None,
            pvalue=3.13e-08,
            pvalue_test="fisher_exact",
            pvalue_adjusted=None,
            pvalue_adjustment=None,
            case_cohorts=("cnchd", "ddd"),
            control_cohorts=("ukbb",),
            method_note=None,
            source="audain2026_sd3",
        )
        for stratum in ("all", "syndromic")
    ]

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={TBX5: rows},
        # Derived the way the runner derives it, so this asserts the real
        # shape rather than a hand-written one.
        concordance={
            **_no_burden(),
            TBX5: gene_concordance(rows, cohort_families(rows, {})),
        },
        profiles={},
        cell_types={},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [row["cohort_stratum"] for row in bundle["burden"]] == ["all", "syndromic"]
    assert bundle["burden"][0]["case_cohorts"] == ["cnchd", "ddd"]
    # The unbounded row publishes no number and keeps its lower bound, which is
    # the whole finding. `is None`, not `== None`: a consumer branches on it.
    assert bundle["burden"][0]["effect"] is None
    assert bundle["burden"][0]["effect_bound"] == "unbounded_above"
    assert bundle["burden"][0]["ci_low"] == 28.1

    rows_by_gene = _by_gene(tmp_path)
    assert rows_by_gene[TBX5]["burden_row_count"] == 2
    # A published gene no study reported carries the empty shape, never a
    # missing key: "no study reported this gene" must not read as a dropped join.
    assert rows_by_gene[GATA4]["burden_row_count"] == 0
    assert _read(tmp_path, "genes/HGNC_4173.json")["burden"] == []


def test_a_published_gene_with_no_concordance_derived_fails_loudly(tmp_path: Path) -> None:
    """`build_genes` refuses rather than defaulting, and the refusal names the gap.

    `concordance` is a required argument precisely so a caller cannot forget it.
    But a caller can still build it over the wrong population -- the published
    set changed, the mapping did not -- and the failure a silent default would
    produce is the worst kind this project has: a bundle publishing
    `"tested": 0`, which is indistinguishable from the true statement that no
    study reported that gene.

    So the guard `raise`s, and it says which gene and why. `raise` rather than
    `assert`: `-O` strips `assert`.

    A mutation matrix on 2026-08-05 found this guard unguarded -- disabling it
    survived the whole suite -- which is the only reason this test exists.
    """
    emitter = Emitter(root=tmp_path)

    with pytest.raises(KeyError, match="no concordance derived for published gene"):
        build_genes(
            _corpus(),
            emitter,
            symbols=SYMBOLS,
            omics={},
            variants={},
            validity={},
            published={TBX5},
            burden={},
            concordance={},
            profiles={},
            cell_types={},
        )


def test_a_bundle_carries_the_identity_its_registry_row_holds(tmp_path: Path) -> None:
    """D52. TBX5's real registry row, and GATA4 with none: six nulls, never a missing key."""
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        identities={
            TBX5: GeneIdentity(
                name="T-box transcription factor 5",
                locus="12q24.21",
                ensembl_gene="ENSG00000089225",
                ncbi_gene=6910,
                uniprot="Q99593",
                mane_select="NM_181486.4",
            )
        },
        omics={},
        variants={},
        validity={},
        published={TBX5, GATA4},
        burden={},
        concordance=_no_burden(),
        profiles={},
        cell_types={},
    )

    tbx5 = _read(tmp_path, "genes/HGNC_11604.json")["identity"]
    assert tbx5["ensembl_gene"] == "ENSG00000089225"
    assert tbx5["ncbi_gene"] == 6910 and type(tbx5["ncbi_gene"]) is int
    assert _read(tmp_path, "genes/HGNC_4173.json")["identity"] == dict.fromkeys(tbx5)
