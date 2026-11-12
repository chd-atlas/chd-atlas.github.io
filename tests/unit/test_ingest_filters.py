"""The two tables `chd-atlas ingest` keeps only part of (spec F23).

An Orphanet product 4 package carries every phenotype annotation and a Mondo
SSSOM package every exact match; ingest keeps the annotations under HP:0001626
in the HPO release `mirrors/sources.yaml` pins, and the Mondo terms either
validity mirror names, so the committed mirrors keep the rows the converter
(`scripts/convert_orphanet.py`) kept. The ingest tests run on `small_repo`
(`genes.tsv` and `sources.yaml` copied, validation stubbed), so the filter is
the only thing between a package row and the mirror.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Callable
from pathlib import Path

import pytest

from chd_atlas.ingest.apply import IngestResult, ingest_package
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.orphanet_scope import (
    pinned_hpo_release,
    pinned_structural_phenotypes,
    structural_phenotypes,
)
from chd_atlas.tables import TABLE_SCHEMAS
from chd_atlas.validate.sources import load_sources
from tests.unit.test_ingest_apply import REPO_ROOT, small_repo, tree, write_tables

COMMITTED_HPO = "ontologies/hp-2026-06-23.obo"
ORPHANET = "orphanet_product4"
MONDO = "mondo_sssom"


def _obo(*terms: tuple[str, str | None, str | None]) -> str:
    """An OBO of ``(id, parent, alt_id)`` terms, each a `[Term]` stanza."""
    stanzas = ["format-version: 1.2\n"]
    for term, parent, alt in terms:
        lines = ["[Term]", f"id: {term}", f"name: {term}"]
        lines += [f"alt_id: {alt}"] if alt else []
        lines += [f"is_a: {parent} ! parent"] if parent else []
        stanzas.append("\n".join(lines) + "\n")
    return "\n".join(stanzas)


# HP:0001626 -> HP:0001627 -> HP:0001631 (alt id HP:0009999), and HP:0000256
# under HP:0000118, outside the cardiovascular branch.
FIRST = _obo(
    ("HP:0000118", None, None),
    ("HP:0001626", "HP:0000118", None),
    ("HP:0001627", "HP:0001626", None),
    ("HP:0001631", "HP:0001627", "HP:0009999"),
    ("HP:0000256", "HP:0000118", None),
)
# The same release but for HP:0000256, which this one files under HP:0001626.
SECOND = _obo(
    ("HP:0000118", None, None),
    ("HP:0001626", "HP:0000118", None),
    ("HP:0001627", "HP:0001626", None),
    ("HP:0001631", "HP:0001627", "HP:0009999"),
    ("HP:0000256", "HP:0001626", None),
)


def orphanet_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, obo: str | bytes) -> Path:
    """`small_repo` with ``obo`` at the path `sources.yaml` pins for HP."""
    repo = small_repo(tmp_path, monkeypatch)
    path = repo / COMMITTED_HPO
    path.parent.mkdir(parents=True)
    path.write_bytes(obo.encode("utf-8") if isinstance(obo, str) else obo)
    return repo


def phenotype_rows(*hpo: str) -> list[dict[str, str]]:
    return [
        {
            "orpha": "ORPHA:648",
            "disorder_name": "Noonan syndrome",
            "hpo": term,
            "hpo_label": f"label of {term}",
            "frequency": "Frequent (79-30%)",
            "source": ORPHANET,
        }
        for term in hpo
    ]


PACKAGE_HPO = ("HP:0001631", "HP:0009999", "HP:0000256", "HP:7777777")


def ingest_phenotypes(repo: Path, tmp_path: Path) -> IngestResult:
    tables = [("orphanet_phenotypes", None, phenotype_rows(*PACKAGE_HPO))]
    package = write_tables(tmp_path / "pkg", tables, package=ORPHANET)
    return ingest_package(repo, ORPHANET, package)


def mirrored_hpo(repo: Path) -> list[str]:
    header, *lines = (repo / "mirrors" / "orphanet_phenotypes.tsv").read_text().splitlines()
    assert header.split("\t") == list(TABLE_SCHEMAS["orphanet_phenotypes"].column_names)
    return [line.split("\t")[2] for line in lines]


def test_an_orphanet_package_keeps_only_the_cardiovascular_branch_of_the_pinned_hpo(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """An alt id of a term in the branch is kept; an id the release lacks is named, not counted.

    HP:7777777 is in no release, so whether it is cardiovascular is unknown:
    it is listed by name and its row is not written, but it is not one of the
    rows "outside the filter", which would read as a fact about the term.
    """
    repo = orphanet_repo(tmp_path, monkeypatch, FIRST)

    result = ingest_phenotypes(repo, tmp_path)

    assert mirrored_hpo(repo) == ["HP:0001631", "HP:0009999"]
    assert (result.rows_in, result.rows_written) == (4, 2)
    assert (result.rows_outside_filter, result.rows_unknown_term) == (1, 1)
    assert result.unknown_terms == ("HP:7777777",)
    assert result.filter_reasons == (f"hpo not under HP:0001626 in {COMMITTED_HPO}",)
    assert (result.rows_outside_registry, result.rows_unresolved) == (0, 0)


def test_the_orphanet_filter_reads_the_hpo_release_sources_yaml_pins(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One HPO release for the mirror and for the rule that reads it.

    SCP006 and the build decide what is structural over the release
    `mirrors/sources.yaml` pins (`pinned_structural_phenotypes`), so the
    filter must read that one too: filtering through another would leave an
    annotation the pinned release places under the heart unmirrored, and the
    syndrome it qualifies failing the rule with nothing saying why. The first
    release stays at the path the committed file pins, so a filter reading a
    fixed path instead of the pin keeps the first release's answer and fails.
    (This replaces `test_convert_orphanet.py`'s coupling test, which compared
    the converter's literal path with the pin.)
    """
    repo = orphanet_repo(tmp_path, monkeypatch, FIRST)
    (repo / "ontologies" / "hp-second.obo").write_text(SECOND, encoding="utf-8")
    sources = repo / "mirrors" / "sources.yaml"
    text = sources.read_text(encoding="utf-8")
    pin = f"ontology_file: {COMMITTED_HPO}\n"
    assert text.count(pin) == 1
    sources.write_text(text.replace(pin, "ontology_file: ontologies/hp-second.obo\n"))

    result = ingest_phenotypes(repo, tmp_path)

    assert mirrored_hpo(repo) == ["HP:0000256", "HP:0001631", "HP:0009999"]
    assert (result.rows_outside_filter, result.rows_unknown_term) == (0, 1)
    assert result.filter_reasons == ("hpo not under HP:0001626 in ontologies/hp-second.obo",)


def validity_row(relative: str, disease: str) -> dict[str, str]:
    """The committed mirror's first row with its ``disease`` replaced."""
    with (REPO_ROOT / relative).open(encoding="utf-8") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    return {**row, "disease": disease}


def write_mirror(repo: Path, relative: str, rows: list[dict[str, str]]) -> None:
    with (REPO_ROOT / relative).open(encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    lines = ["\t".join(header), *("\t".join(row[c] for c in header) for row in rows)]
    (repo / relative).write_text("\n".join(lines) + "\n", encoding="utf-8")


CLINGEN = "mirrors/clingen_gene_validity.tsv"
GENCC = "mirrors/gencc_submissions.tsv"
# A scope naming only the term no validity mirror names, so a filter keyed on
# the scope rather than the mirrors keeps exactly the wrong row.
SCOPE = """diseases:
  - id: MONDO:0000003
    label: a term only the scope names
    reason: a fixture
    admitted_by: gencc_submitter
    attributed_to: Ambry Genetics
    admitted_on: 2026-11-11
"""


def mondo_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """ClinGen names MONDO:0000001, GenCC MONDO:0000002, the scope MONDO:0000003."""
    repo = small_repo(tmp_path, monkeypatch)
    write_mirror(repo, CLINGEN, [validity_row(CLINGEN, "MONDO:0000001")])
    write_mirror(repo, GENCC, [validity_row(GENCC, "MONDO:0000002")])
    (repo / "curation").mkdir()
    (repo / "curation" / "chd_scope.yaml").write_text(SCOPE, encoding="utf-8")
    return repo


def ingest_xrefs(repo: Path, tmp_path: Path) -> IngestResult:
    rows = [
        {"mondo": "MONDO:0000001", "target": "ORPHA:1", "source": MONDO},
        {"mondo": "MONDO:0000002", "target": "OMIM:100002", "source": MONDO},
        {"mondo": "MONDO:0000003", "target": "ORPHA:3", "source": MONDO},
    ]
    package = write_tables(tmp_path / "pkg", [("mondo_xrefs", None, rows)], package=MONDO)
    return ingest_package(repo, MONDO, package)


def test_a_mondo_package_keeps_the_terms_a_validity_mirror_names(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Either mirror, every row of it whatever its classification, as it stands when ingest runs.

    `scripts/convert_orphanet.py::_mirrored_terms` kept the same set (the
    `disease` of every row of both mirrors). Measured 2026-11-12 over the
    committed mirrors: all 10,011 `mondo_xrefs.tsv` rows name a disease of
    one of them, over 7,418 distinct MONDO terms.
    """
    repo = mondo_repo(tmp_path, monkeypatch)

    result = ingest_xrefs(repo, tmp_path)

    lines = (repo / "mirrors" / "mondo_xrefs.tsv").read_text(encoding="utf-8").splitlines()
    assert lines == [
        "mondo\ttarget\tsource",
        f"MONDO:0000001\tORPHA:1\t{MONDO}",
        f"MONDO:0000002\tOMIM:100002\t{MONDO}",
    ]
    assert (result.rows_in, result.rows_written, result.rows_outside_filter) == (3, 2, 1)
    assert (result.rows_unknown_term, result.unknown_terms) == (0, ())
    assert result.filter_reasons == (f"mondo named by neither {CLINGEN} nor {GENCC}",)


Setup = Callable[[Path, pytest.MonkeyPatch], tuple[Path, str]]


def _undecodable_obo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, str]:
    return orphanet_repo(tmp_path, monkeypatch, b"[Term]\nid: HP:0001626\xff\n"), ORPHANET


def _obo_without_the_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, str]:
    return orphanet_repo(tmp_path, monkeypatch, _obo(("HP:0000118", None, None))), ORPHANET


def _unreadable_clingen(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, str]:
    repo = mondo_repo(tmp_path, monkeypatch)
    (repo / CLINGEN).write_text("gene\tdisease\nHGNC:1\tMONDO:0000001\n", encoding="utf-8")
    return repo, MONDO


def _no_validity_mirror(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, str]:
    repo = mondo_repo(tmp_path, monkeypatch)
    (repo / CLINGEN).unlink()
    (repo / GENCC).unlink()
    return repo, MONDO


@pytest.mark.parametrize(
    ("setup", "message"),
    [
        (
            _undecodable_obo,
            f"the orphanet_phenotypes filter cannot be built: the pinned HPO release "
            f"{COMMITTED_HPO} cannot be read",
        ),
        (
            _obo_without_the_root,
            f"the orphanet_phenotypes filter cannot be built: the pinned HPO release "
            f"{COMMITTED_HPO} has no term HP:0001626",
        ),
        (
            _unreadable_clingen,
            f"the mondo_xrefs filter cannot be built: {CLINGEN} fails its schema",
        ),
        (
            _no_validity_mirror,
            f"the mondo_xrefs filter cannot be built: neither {CLINGEN} nor {GENCC} names a "
            "disease, so the filter would drop every row",
        ),
    ],
    ids=["undecodable-obo", "obo-without-hp-0001626", "unreadable-clingen", "no-validity-mirror"],
)
def test_a_filter_that_cannot_be_built_is_ing010_before_any_write(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, setup: Setup, message: str
) -> None:
    """Refused rather than filtered to nothing, which would read as a fact about the package.

    An absent validity mirror counts as one with no rows, as it does for the
    registry (`chd_atlas.registry`); only both absent leaves nothing to keep.
    """
    repo, source = setup(tmp_path, monkeypatch)
    if source == ORPHANET:
        rows = phenotype_rows(*PACKAGE_HPO)
        tables = [("orphanet_phenotypes", None, rows)]
    else:
        rows = [{"mondo": "MONDO:0000001", "target": "ORPHA:1", "source": MONDO}]
        tables = [("mondo_xrefs", None, rows)]
    package = write_tables(tmp_path / "pkg", tables, package=source)
    before = tree(repo)

    with pytest.raises(IngestRefused) as excinfo:
        ingest_package(repo, source, package)

    assert excinfo.value.code == "ING010"
    assert excinfo.value.message.startswith(message)
    assert tree(repo) == before


def test_every_committed_orphanet_phenotype_is_in_the_pinned_cardiovascular_branch() -> None:
    """The filter keeps every row the converter kept: 7,902 rows over 730 terms.

    Read through `pinned_hpo_release` over the release `sources.yaml` pins, as
    the ingest reads it. Measured 2026-11-12 over the committed mirror and
    `hp-2026-06-23.obo`: 561 of its 730 terms (5,514 rows) are outside both
    `structural_phenotypes` and HP:0001627's closure, so a closure started at
    HP:0001627 would drop them; HP:0000822 (Hypertension, 201 rows) is the one
    named below. 5 terms (9 rows) are alt ids (HP:0007850, HP:0008014,
    HP:0031014, HP:0031297, HP:0031349), so a closure without alt ids would
    drop those.
    """
    registry, _ = load_sources(REPO_ROOT)
    pinned = registry.ontology_files()
    release = pinned_hpo_release(REPO_ROOT, pinned)
    structural = pinned_structural_phenotypes(REPO_ROOT, pinned)
    assert release is not None and structural is not None
    with (REPO_ROOT / "mirrors" / "orphanet_phenotypes.tsv").open(encoding="utf-8") as handle:
        terms = [row["hpo"] for row in csv.DictReader(handle, delimiter="\t")]

    assert len(terms) == 7902
    assert sorted(set(terms) - release.cardiovascular) == []
    assert "HP:0000822" in terms
    assert "HP:0000822" in release.cardiovascular
    assert "HP:0000822" not in structural
    assert "HP:0031014" in terms
    assert "HP:0031014" in release.cardiovascular


def test_structural_phenotypes_is_unchanged_by_the_shared_reader() -> None:
    """The structural set the rule applies did not move when its OBO reader was shared.

    Measured 2026-11-12 at 4b1f8f6, before the reader was factored out, over
    the pinned `hp-2026-06-23.obo`: 291 terms, and the sha256 below of their
    sorted ids joined by newlines.
    """
    structural = structural_phenotypes(REPO_ROOT / COMMITTED_HPO)

    assert len(structural) == 291
    digest = hashlib.sha256("\n".join(sorted(structural)).encode("utf-8")).hexdigest()
    assert digest == "5dd4bf2b8bab2344dac1729018ed07ac854bffb5ed9858aa46a7f8db9335a7de"
