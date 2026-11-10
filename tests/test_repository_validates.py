import shutil
from pathlib import Path

from chd_atlas.validate.runner import validate_repository

REPO_ROOT = Path(__file__).parent.parent


def test_the_committed_repository_validates_cleanly() -> None:
    """0 errors and the baseline 3 warnings (two SCP003, one SCP007), measured 2026-11-10."""
    report = validate_repository(REPO_ROOT)
    assert report.ok, report.render()
    assert report.warning_count == 3, report.render()


def test_the_seed_corpus_is_not_empty() -> None:
    from chd_atlas.corpus import load_curation

    corpus, _ = load_curation(REPO_ROOT)
    assert corpus.assertions
    assert corpus.publications
    assert corpus.featured


def test_the_exported_schemas_are_committed() -> None:
    """test_committed_schemas_match_the_models returns early when schemas/ is
    absent, so without this it passes vacuously and never checks for drift."""
    from chd_atlas.schema_export import EXPORTED_MODELS

    committed = {path.name for path in (REPO_ROOT / "schemas").glob("*.schema.json")}
    assert committed == {f"{name}.schema.json" for name in EXPORTED_MODELS}

    from chd_atlas.tables import TABLE_SCHEMAS

    tables = {path.name for path in (REPO_ROOT / "schemas" / "tables").glob("*.json")}
    assert tables == {f"{name}.v1.json" for name in TABLE_SCHEMAS}


def test_a_misnamed_assertion_in_the_real_corpus_fails_validation(tmp_path: Path) -> None:
    """The case that motivated the stray-entry sweep: renaming an assertion to
    .yml made it invisible to every check while validate still exited 0."""
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, tmp_path / item)
    assert validate_repository(tmp_path).ok, "the copied corpus should start clean"

    assertions = tmp_path / "curation" / "assertions"
    original = next(assertions.glob("*.yaml"))
    original.rename(original.with_suffix(".yml"))

    report = validate_repository(tmp_path)

    assert report.ok is False
    assert "CUR001" in [issue.code for issue in report.issues]


def test_missing_validity_mirrors_fails_validation(tmp_path: Path) -> None:
    """Reproduces the gap found reviewing Task 10's wiring: with both validity
    mirrors gone and everything else in the corpus valid, `_mirrored_validity`
    returns None and only SCP000 -- a WARNING -- was reported. `ok` ignores
    warnings, so `build_site` (which gates on `report.ok`) would have
    published a site where scope, and so which genes exist at all, could not
    be checked, while `chd-atlas validate` still exited 0. TBL012 makes the
    missing mirrors an ERROR, matching how TBL008 covers a missing gene
    registry the same way.

    `is False`, not `== False`: `{"f": 0} == {"f": False}` is `True` in
    Python, and this project pins boolean contracts with `is`.
    """
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, tmp_path / item)
    assert validate_repository(tmp_path).ok, "the copied corpus should start clean"

    (tmp_path / "mirrors" / "clingen_gene_validity.tsv").unlink()
    (tmp_path / "mirrors" / "gencc_submissions.tsv").unlink()

    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert report.ok is False
    assert "TBL012" in codes
    assert "SCP000" in codes


# Every `design="profile"` dataset's whole stage sequence, in curated order.
#
# Keyed by accession and asserted exhaustively, because the version this
# replaced was hardcoded to `id == "E-MTAB-6814"` and the slice `ordered[13:]`
# -- so a second profile dataset inherited no pin at all, and one committed
# with an inverted block validated 0 errors / 4 warnings, exactly like the
# clean baseline. A dataset added without an entry here now fails rather than
# passing silently.
#
# The prenatal half is pinned too, not only the post-natal block PRF012 cannot
# reach: it costs nothing, and a slice boundary is one more number to get
# wrong.
_EXPECTED_STAGE_ORDER: dict[str, list[str]] = {
    "E-MTAB-6814": [
        "4 week post conception",
        "5 week post conception",
        "6 week post conception",
        "7 week post conception",
        "8 week post conception",
        "9 week post conception",
        "10 week post conception",
        "11 week post conception",
        "12 week post conception",
        "13 week post conception",
        "16 week post conception",
        "18 week post conception",
        "19 week post conception",
        "neonate",
        "infant",
        "toddler",
        "school age child",
        "adolescent",
        "young adult",
        "middle adult",
        "elderly",
    ],
}


def test_every_profile_dataset_pins_its_stage_sequence_against_a_literal() -> None:
    """Post-natal tokens carry no `wpc`, so `order` is their only chronology.

    PRF014 now refuses a post-natal stage ordered before a prenatal one, but
    nothing can check the post-natal block's *internal* sequence: the model
    accepts any ascending integers and no validator can contradict a null
    `wpc`. Re-alphabetising that block would silently restore the exact defect
    `Stage.order` removes -- `elderly` published second of eight -- with every
    gate green.

    Asserted against a literal spelled out in `_EXPECTED_STAGE_ORDER`, never
    against `sorted()` of the tokens themselves: sorting them is what the
    defect did.

    **Iterates every profile dataset, and refuses an unpinned one.** The
    version this replaced named one accession and sliced at 13, so it guarded
    exactly one dataset's tail; a second dataset with an inverted block
    measured 0 errors and 4 warnings, byte-identical codes to the clean
    baseline, and this test would not have run against it at all.
    """
    from chd_atlas.corpus import load_curation
    from chd_atlas.models.dataset import DatasetDesign

    corpus, _ = load_curation(REPO_ROOT)
    profiles = {
        dataset.id: dataset
        for dataset in corpus.datasets
        if dataset.design is DatasetDesign.PROFILE
    }

    assert profiles, "no profile dataset in the corpus; this test measures nothing"
    assert set(profiles) == set(_EXPECTED_STAGE_ORDER), (
        "a profile dataset has no pinned stage sequence (or a pin names no dataset): "
        f"{sorted(set(profiles) ^ set(_EXPECTED_STAGE_ORDER))}"
    )

    for accession, dataset in sorted(profiles.items()):
        ordered = [stage.token for stage in sorted(dataset.stages, key=lambda item: item.order)]
        assert ordered == _EXPECTED_STAGE_ORDER[accession], accession
