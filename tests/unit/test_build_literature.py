# tests/unit/test_build_literature.py
import json
from collections.abc import Collection
from datetime import date
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.literature import build_literature
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import FeaturedManuscript, PhenotypeTerm, Publication


def _publication(pmid: str = "PMID:8988165", **overrides: object) -> Publication:
    payload: dict[str, object] = {
        "id": pmid,
        "title": "Mutations in human TBX5 cause Holt-Oram syndrome",
        "journal": "Nature Genetics",
        "year": 1997,
        "authors": ["Li QY"],
        "study_type": "family_linkage",
    }
    payload.update(overrides)
    return Publication.model_validate(payload)


def _phenotype(
    term: str = "HP:0001631", label: str = "Atrial septal defect", **overrides: object
) -> PhenotypeTerm:
    payload: dict[str, object] = {"id": term, "label": label, "lesion_group": "septal"}
    payload.update(overrides)
    return PhenotypeTerm.model_validate(payload)


def _dataset(accession: str = "GSE1000", **overrides: object) -> Dataset:
    payload: dict[str, object] = {
        "id": accession,
        "archive": "geo",
        "technology": "bulk_rnaseq",
        "design": "contrast",
        "tissue": "right ventricular myocardium",
        "developmental_stage": "infant",
        "organism": "NCBITaxon:9606",
        "n_samples": 6,
        "licence": "CC-BY-4.0",
        "contrasts": [
            {
                "id": "tof_vs_control",
                "description": "Tetralogy of Fallot versus structurally normal heart",
                "case_group": "tof",
                "control_group": "control",
                "statistical_method": "Wald test",
                "software": "DESeq2 1.44.0",
            }
        ],
    }
    payload.update(overrides)
    return Dataset.model_validate(payload)


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {
        "root": Path("."),
        "publications": (_publication(),),
        "featured": (
            FeaturedManuscript(publication="PMID:8988165", order=1, blurb="b", topic="genomics"),
        ),
        "phenotypes": (_phenotype(),),
    }
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def _payload(root: Path, name: str, key: str) -> list[dict[str, object]]:
    """One emitted payload's array, read back from the file that was written."""
    written = json.loads((root / name).read_text())
    records: list[dict[str, object]] = written[key]
    return records


def _cited(entry: dict[str, object]) -> str:
    """The PMID of the publication one featured entry carries inline."""
    publication = entry["publication"]
    assert isinstance(publication, dict)
    return str(publication["id"])


def test_emits_all_four_payloads(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    build_literature(_corpus(), emitter, burden_studies=())

    assert set(emitter.checksums) == {
        "publications.json",
        "featured.json",
        "phenotypes.json",
        "datasets.json",
    }


def test_featured_manuscripts_are_emitted_in_display_order(tmp_path: Path) -> None:
    """The front page renders this list as given; order is the curator's decision.

    Nothing else in these two entries produces this result: their PMIDs, blurbs,
    topics and corpus order each rank the pair the other way round, so a build
    that sorted on any of those — or that emitted the corpus order untouched —
    fails here rather than passing on a fixture that was already sorted.
    """
    corpus = _corpus(
        publications=(_publication("PMID:1"), _publication("PMID:2")),
        featured=(
            FeaturedManuscript(
                publication="PMID:1", order=2, blurb="a broader review", topic="genomics"
            ),
            FeaturedManuscript(
                publication="PMID:2", order=1, blurb="zebrafish knockdown", topic="review"
            ),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    featured = _payload(tmp_path, "featured.json", "featured")
    assert [entry["order"] for entry in featured] == [1, 2]


def test_the_inline_publication_is_the_one_the_entry_cites(tmp_path: Path) -> None:
    """The landing page is one fetch: resolving a PMID against publications.json
    client-side would make first paint depend on a second request.

    Asserted against two curated publications rather than one, because with one
    "the publication it cites" and "any publication" are the same record. A build
    that inlined the first — or the wrong element of a list it had sorted
    separately — puts a manuscript on the front page under another's blurb.
    """
    corpus = _corpus(
        publications=(
            _publication("PMID:1", title="Not this one"),
            _publication("PMID:2", title="The featured one"),
        ),
        featured=(FeaturedManuscript(publication="PMID:2", order=1, blurb="b", topic="genomics"),),
    )
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    featured = _payload(tmp_path, "featured.json", "featured")
    assert featured[0]["publication"] == {
        "id": "PMID:2",
        "title": "The featured one",
        "journal": "Nature Genetics",
        "year": 1997,
        "authors": ["Li QY"],
        "study_type": "family_linkage",
        "doi": None,
        "pmcid": None,
        "own_lab": False,
        "cohort_size": None,
        "ancestry": [],
        "tests_reported": None,
        "correction_genes": None,
        "genome_wide_threshold": None,
        "summary_caveat": None,
    }


def test_featured_entries_sharing_a_display_order_are_ranked_by_pmid(tmp_path: Path) -> None:
    """`FeaturedFile.display_order_is_unique` forbids this pair, and featured.yaml is
    the only file the list is loaded from, so a tie is reachable only by building a
    `Corpus` directly — as this does.

    Pinned anyway because the alternative is a tie broken by the order the two
    records happen to arrive in, which is not content: nothing about the records
    themselves would decide which manuscript leads the page.
    """
    corpus = _corpus(
        publications=(_publication("PMID:1"), _publication("PMID:2")),
        featured=(
            FeaturedManuscript(publication="PMID:2", order=1, blurb="b", topic="genomics"),
            FeaturedManuscript(publication="PMID:1", order=1, blurb="b", topic="genomics"),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    featured = _payload(tmp_path, "featured.json", "featured")
    assert [_cited(entry) for entry in featured] == ["PMID:1", "PMID:2"]


def test_a_featured_entry_keeps_the_curators_own_words(tmp_path: Path) -> None:
    """The blurb says why a manuscript leads the front page and exists nowhere else
    in the corpus, so a payload that inlined the publication and dropped everything
    around it would lose the only writing the entry contributes."""
    corpus = _corpus(
        featured=(
            FeaturedManuscript(
                publication="PMID:8988165",
                order=1,
                blurb="The founding demonstration of TBX5 haploinsufficiency.",
                topic="review",
            ),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    entry = _payload(tmp_path, "featured.json", "featured")[0]
    assert entry["blurb"] == "The founding demonstration of TBX5 haploinsufficiency."
    assert entry["topic"] == "review"
    assert entry["order"] == 1


def test_a_featured_entry_citing_an_unknown_publication_is_refused(tmp_path: Path) -> None:
    """`validate_references` reports this as REF002 and `build_site` refuses to
    build on any error, so reaching here means validation was bypassed.

    Publishing anyway is the costly branch: the entry's `publication` would be a
    bare PMID string where every other entry carries an object, so a page reading
    `entry.publication.title` gets `undefined` with nothing raised anywhere — in
    the one payload the landing page renders without a second fetch.
    """
    corpus = _corpus(
        publications=(),
        featured=(
            FeaturedManuscript(publication="PMID:404", order=1, blurb="b", topic="genomics"),
        ),
    )
    emitter = Emitter(root=tmp_path)

    with pytest.raises(ValueError, match="PMID:404"):
        build_literature(corpus, emitter, burden_studies=())

    # Refused before anything is written, rather than after half the payloads are.
    assert emitter.checksums == {}
    assert not (tmp_path / "publications.json").exists()


@pytest.mark.parametrize(
    ("corpus", "name", "key", "expected"),
    [
        pytest.param(
            _corpus(
                publications=(
                    _publication("PMID:9", title="third"),
                    _publication("PMID:10", title="second"),
                    _publication("PMID:100", title="first"),
                ),
                featured=(),
            ),
            "publications.json",
            "publications",
            ["PMID:10", "PMID:100", "PMID:9"],
            id="publications",
        ),
        pytest.param(
            _corpus(datasets=(_dataset("GSE9"), _dataset("GSE10"), _dataset("GSE100"))),
            "datasets.json",
            "datasets",
            ["GSE10", "GSE100", "GSE9"],
            id="datasets",
        ),
    ],
)
def test_identifier_arrays_are_ordered_lexically(
    tmp_path: Path, corpus: Corpus, name: str, key: str, expected: list[str]
) -> None:
    """A published array's order is the order every consumer downloads.

    Lexical on the identifier, which puts PMID:10 and PMID:100 ahead of PMID:9.
    That is deliberate: a PMID is issued when a record is indexed, so PMID order
    is neither chronological nor topical and nothing displays this array in it —
    a reader is shown publications by year, or in the curator's featured order.
    What the sort has to do is agree between two builds of one commit, and agree
    with `derive.py`, which orders a gene's publication list the same way. Making
    only this array numeric would leave the two disagreeing about the order of
    the same PMIDs. For an accession there is not even a numeric order to choose:
    `GSE`, `PXD`, `E-` and `EGAS` would each need a parsing rule, and that rule
    is itself a thing that can be wrong.

    Neither the corpus order nor the titles produce these results.
    """
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    published = _payload(tmp_path, name, key)
    assert [record["id"] for record in published] == expected


@pytest.mark.parametrize(
    ("corpus", "name", "key", "field", "expected"),
    [
        pytest.param(
            _corpus(
                publications=(
                    _publication("PMID:1", title="first transcription"),
                    _publication("PMID:1", title="second transcription"),
                ),
                featured=(),
            ),
            "publications.json",
            "publications",
            "title",
            ["first transcription", "second transcription"],
            id="publications",
        ),
        pytest.param(
            _corpus(datasets=(_dataset("GSE1000", n_samples=6), _dataset("GSE1000", n_samples=99))),
            "datasets.json",
            "datasets",
            "n_samples",
            [6, 99],
            id="datasets",
        ),
    ],
)
def test_records_colliding_on_one_identifier_are_all_published(
    tmp_path: Path,
    corpus: Corpus,
    name: str,
    key: str,
    field: str,
    expected: list[object],
) -> None:
    """Keying an emitted array on the identifier drops the second of two records
    sharing one, and drops it silently: the build stays green, the manifest counts
    what it wrote, and nothing anywhere says a curated record was discarded. A
    record published twice is at least visible to whoever reads the file.

    The two cases are not equally reachable, which is the fact worth recording.
    `PublicationFile.ids_are_unique` rejects the whole file, so a duplicate PMID
    reaches `load_curation` as SCHEMA001 and leaves `corpus.publications` empty —
    this module never sees that pair. A dataset is one file per accession, so two
    files both declaring GSE1000 load with no issue at all and both records arrive
    here; only `REF011` reports them, and that is a referential check run beside
    the loader, not a guard inside it.
    """
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    assert [record[field] for record in _payload(tmp_path, name, key)] == expected


def test_phenotypes_are_ordered_by_identifier_and_keep_their_synonyms(tmp_path: Path) -> None:
    """Synonyms are what a search box matches "VSD" on, so dropping them is silent.

    The labels rank these three differently from their identifiers, and the
    corpus order differs from both.
    """
    corpus = _corpus(
        phenotypes=(
            _phenotype("MONDO:0005453", "Atrioventricular septal defect", lesion_group="avsd"),
            _phenotype("HP:0001631", "Atrial septal defect", synonyms=["ASD"]),
            _phenotype("HP:0001629", "Ventricular septal defect", synonyms=["VSD"]),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    published = _payload(tmp_path, "phenotypes.json", "phenotypes")
    assert [term["id"] for term in published] == [
        "HP:0001629",
        "HP:0001631",
        "MONDO:0005453",
    ]
    assert published[0]["synonyms"] == ["VSD"]


def test_a_dataset_is_published_with_the_contrasts_that_interpret_it(tmp_path: Path) -> None:
    """A fold change is uninterpretable without knowing what was compared.

    The omics shards carry the numbers and name a contrast id; this payload is
    the only place that id is defined, so a dataset published without its
    contrasts makes every shard row unreadable.
    """
    emitter = Emitter(root=tmp_path)

    build_literature(_corpus(datasets=(_dataset(),)), emitter, burden_studies=())

    published = _payload(tmp_path, "datasets.json", "datasets")
    contrasts = published[0]["contrasts"]
    assert isinstance(contrasts, list)
    assert contrasts[0]["id"] == "tof_vs_control"
    assert contrasts[0]["case_group"] == "tof"
    assert published[0]["n_samples"] == 6


def test_vocabularies_are_published_as_their_values(tmp_path: Path) -> None:
    """What the front end reads off these fields, pinned as a wire contract.

    Each is a `StrEnum`, which `json.dumps` writes as its value whether the model
    was dumped in json or in python mode — so this does not pin that choice, and
    nothing here can. What it pins is that the value reaches the payload at all,
    as a bare string rather than renamed or nested. `lesion_group` and `archive`
    are asserted nowhere else in this file, and they are what a browse facet and
    an archive filter switch on.
    """
    corpus = _corpus(datasets=(_dataset(),))
    emitter = Emitter(root=tmp_path)

    build_literature(corpus, emitter, burden_studies=())

    assert _payload(tmp_path, "publications.json", "publications")[0]["study_type"] == (
        "family_linkage"
    )
    assert _payload(tmp_path, "featured.json", "featured")[0]["topic"] == "genomics"
    assert _payload(tmp_path, "phenotypes.json", "phenotypes")[0]["lesion_group"] == "septal"
    assert _payload(tmp_path, "datasets.json", "datasets")[0]["archive"] == "geo"


def test_an_empty_corpus_still_emits_every_payload(tmp_path: Path) -> None:
    """A consumer fetching phenotypes.json must not 404 before curation starts."""
    emitter = Emitter(root=tmp_path)

    build_literature(Corpus(root=Path(".")), emitter, burden_studies=())

    assert json.loads((tmp_path / "phenotypes.json").read_text()) == {"phenotypes": []}
    assert json.loads((tmp_path / "publications.json").read_text()) == {"publications": []}
    assert json.loads((tmp_path / "featured.json").read_text()) == {"featured": []}
    assert json.loads((tmp_path / "datasets.json").read_text()) == {"datasets": []}


# --- withholding a publication only an unapproved literature finding cites ---


def _finding(pmid: str, approved: bool, number: int = 1) -> dict[str, object]:
    finding: dict[str, object] = {
        "id": f"CHDA:LIT:{number:07d}",
        "publication": pmid,
        "design": "cohort_screen",
        "methods": ["exome"],
        "count_unit": "probands",
        "screened": 7,
        "carriers": {"de_novo": 2, "inherited": 1, "origin_unknown": 3},
        "variant_classes": ["snv_indel"],
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "syndromic": "both",
        "shares_probands_with": "unknown",
        "locators": [{"kind": "in_text", "element": "Table 2"}],
        "summary": "A summary.",
        "extracted_by": "agent:extract/test",
    }
    if approved:
        finding["approved_by"] = "A. Curator"
        finding["approved_on"] = "2026-11-05"
    return finding


def _literature(*findings: dict[str, object]) -> tuple[tuple[Path, GeneLiteratureFile], ...]:
    file = GeneLiteratureFile.model_validate(
        {
            "gene": "HGNC:11604",
            "searches": [
                {
                    "id": 1,
                    "source": "pubmed",
                    "query": "q",
                    "run_on": "2026-11-02",
                    "hits": 1,
                    "new_hits": 1,
                }
            ],
            "triage": [],
            "findings": list(findings),
        }
    )
    return ((Path("TBX5.yaml"), file),)


def _published(
    tmp_path: Path, corpus: Corpus, burden_studies: Collection[str] = ()
) -> list[object]:
    build_literature(corpus, Emitter(root=tmp_path), burden_studies=burden_studies)
    return [p["id"] for p in _payload(tmp_path, "publications.json", "publications")]


def test_a_publication_only_an_unapproved_finding_cites_is_withheld(tmp_path: Path) -> None:
    """Owner decision 3: the agent's study_type guess reaches no byte before approval."""
    corpus = _corpus(
        publications=(_publication("PMID:1"), _publication("PMID:2"), _publication("PMID:3")),
        featured=(),
        literature=_literature(
            _finding("PMID:1", approved=False, number=1),
            _finding("PMID:2", approved=True, number=2),
        ),
    )

    assert _published(tmp_path, corpus) == ["PMID:2", "PMID:3"]


def test_an_approved_finding_among_unapproved_ones_publishes_the_publication(
    tmp_path: Path,
) -> None:
    corpus = _corpus(
        publications=(_publication("PMID:1"),),
        featured=(),
        literature=_literature(
            _finding("PMID:1", approved=False, number=1),
            _finding("PMID:1", approved=True, number=2),
        ),
    )

    assert _published(tmp_path, corpus) == ["PMID:1"]


@pytest.mark.parametrize("citer", ["featured", "burden", "dataset", "assertion", "functional"])
def test_another_citer_keeps_a_publication_an_unapproved_finding_also_cites(
    tmp_path: Path, citer: str
) -> None:
    """Every citer `validate_references` and the burden validators resolve a PMID for."""
    extra: dict[str, object] = {}
    burden: set[str] = set()
    if citer == "featured":
        extra["featured"] = (
            FeaturedManuscript(publication="PMID:1", order=1, blurb="b", topic="genomics"),
        )
    elif citer == "burden":
        burden = {"PMID:1"}
    elif citer == "dataset":
        extra["datasets"] = (_dataset(publication="PMID:1"),)
    elif citer == "assertion":
        extra["assertions"] = (
            LesionAssertion.model_validate(
                {
                    "id": "CHDA:AST:0000001",
                    "gene": "HGNC:11604",
                    "phenotypes": ["HP:0001631"],
                    "lesion_groups": ["septal"],
                    "inheritance": ["AD"],
                    "mechanism": "haploinsufficiency",
                    "syndromic": "both",
                    "evidence": [
                        Evidence.model_validate(
                            {
                                "publication": "PMID:1",
                                "evidence_class": "genetic_case",
                                "locator": SupplementaryLocator(file="T.xlsx"),
                                "summary": "s",
                                "strength": "strong",
                            }
                        )
                    ],
                    "curator": "c",
                    "curated_on": date(2026, 7, 1),
                    "last_reviewed": date(2026, 7, 1),
                }
            ),
        )
    else:
        extra["functional"] = (
            FunctionalEvidence.model_validate(
                {
                    "id": "CHDA:FUN:0000009",
                    "gene": "HGNC:11604",
                    "organism": "NCBITaxon:10090",
                    "perturbation": "knockout",
                    "zygosity": "heterozygous",
                    "cardiac_phenotype": ["MP:0010402"],
                    "phenocopies_human": "partial",
                    "rescue_outcome": "not_attempted",
                    "publication": "PMID:1",
                }
            ),
        )
    corpus = _corpus(
        publications=(_publication("PMID:1"),),
        featured=extra.pop("featured", ()),
        literature=_literature(_finding("PMID:1", approved=False)),
        **extra,
    )

    assert _published(tmp_path, corpus, burden) == ["PMID:1"]


def test_a_publication_no_literature_finding_cites_publishes_as_before(tmp_path: Path) -> None:
    corpus = _corpus(publications=(_publication("PMID:1"),), featured=(), literature=())

    assert _published(tmp_path, corpus) == ["PMID:1"]
