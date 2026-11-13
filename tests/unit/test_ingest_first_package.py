"""A new source's first package ingests with SRC007 and SRC008 live.

`tests/unit/test_ingest_apply.py` switches both codes off for its whole module
(its fixtures keep one pin), `test_cli_ingest.py` stubs ingest, so before this
file nothing ran `ingest_package` through the real validator on a repository that
KEEPS the committed pins and manifests. That is the curator's workflow for a new
source, and the one SRC007/SRC008 sit in: the new source's pin and manifest must
arrive in the same transaction as its rows, or the ingest refuses itself.
"""

import json
import shutil
from pathlib import Path

from chd_atlas.ingest.apply import ingest_package
from chd_atlas.validate.runner import validate_repository
from tests.unit.package_fixtures import burden_row, write_package

REPO_ROOT = Path(__file__).parent.parent.parent
SOURCE = "toy_burden"
SLICE_COMMIT = "0123456789abcdef0123456789abcdef01234567"


def test_a_new_sources_first_package_is_accepted_beside_the_committed_pins(
    tmp_path: Path,
) -> None:
    repo = tmp_path / "repo"
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, repo / item)
    committed_pins = sorted(p.name for p in (repo / "mirrors" / "packages").glob("*.json"))
    assert len(committed_pins) >= 10, "the committed manifests must be kept, not removed"
    package = write_package(
        tmp_path,
        [burden_row(source=SOURCE, maf_max="0.0001")],
        ids="hgnc_id",
        coverage={"registry_slice": {"atlas_commit": SLICE_COMMIT, "genes": 185}},
    )

    result = ingest_package(repo, SOURCE, package)

    assert result.rows_written == 1
    assert sorted(p.name for p in (repo / "mirrors" / "packages").glob("*.json")) == sorted(
        [*committed_pins, f"{SOURCE}.json"]
    )
    report = validate_repository(repo)
    assert [
        i.code for i in report.issues if i.code in {"SRC005", "SRC006", "SRC007", "SRC008"}
    ] == []
    assert report.ok
    manifest = json.loads((repo / "mirrors" / "packages" / f"{SOURCE}.json").read_text())
    assert manifest["package"] == SOURCE
