"""`report-changes` over two real builds: the corpus, and a copy without TBX5 burden rows.

The synthetic suite proves each rule on a fixture; this proves the report
reads what `build_site` actually writes. The pair is built once (two builds,
about 25 s: measured 12.5 s per build on 2026-11-10).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import pytest
from typer.testing import CliRunner

from chd_atlas.build.runner import build_site
from chd_atlas.changes import SiteChanges, render_markdown, site_changes
from chd_atlas.cli import app

REPO_ROOT = Path(__file__).parent.parent
TBX5 = "HGNC:11604"


def copy_repo(target: Path) -> Path:
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, target / item)
    return target


@pytest.fixture(scope="module")
def pair(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path, Path, Path]:
    work = tmp_path_factory.mktemp("changes")
    head_root = copy_repo(work / "head")
    # burden_combined rows name their arms in burden.tsv (BUR020), so TBX5's go too.
    for name in ("burden.tsv", "burden_combined.tsv"):
        table = head_root / "mirrors" / name
        lines = table.read_text(encoding="utf-8").splitlines(keepends=True)
        kept = [line for line in lines if TBX5 not in line.split("\t")[:3]]
        assert len(kept) < len(lines), (
            f"{name}: TBX5 must carry rows for this fixture to mean anything"
        )
        table.write_text("".join(kept), encoding="utf-8")
    base_site, head_site = work / "base-site", work / "head-site"
    build_site(REPO_ROOT, base_site)
    build_site(head_root, head_site)
    return REPO_ROOT, base_site, head_root, head_site


@pytest.fixture(scope="module")
def changes(pair: tuple[Path, Path, Path, Path]) -> SiteChanges:
    base_root, base_site, head_root, head_site = pair
    return site_changes(base_root, base_site, head_root, head_site)


def test_the_report_reads_a_real_build(changes: SiteChanges) -> None:
    [tbx5] = [c for c in changes.mirror_rows if c.table == "burden" and c.gene == TBX5]
    assert tbx5.before > 0 and tbx5.after == 0
    [burden_rows] = [c for c in changes.counts if c[0] == "burden_rows"]
    assert isinstance(burden_rows[1], int) and isinstance(burden_rows[2], int)
    assert burden_rows[2] < burden_rows[1]
    # TBX5 stays published: it loses evidence, it does not leave.
    assert TBX5 not in [gene for gene, _ in changes.left]
    assert TBX5 in [gene.gene for gene in changes.genes]

    text = render_markdown(changes)
    mirror = text.split("### Mirror rows per gene")[1].split("###")[0]
    assert f"- burden {TBX5}: {tbx5.before} → **0**" in mirror
    assert "genes/HGNC_11604.json" in text.split("files changed")[1]


def test_the_command_writes_the_same_report(
    pair: tuple[Path, Path, Path, Path], changes: SiteChanges, tmp_path: Path
) -> None:
    base_root, base_site, head_root, head_site = pair
    out = tmp_path / "report.md"
    result = CliRunner().invoke(
        app,
        [
            "report-changes",
            "--base-root", str(base_root),
            "--base-site", str(base_site),
            "--head-root", str(head_root),
            "--head-site", str(head_site),
            "--out", str(out),
        ],
    )  # fmt: skip
    assert result.exit_code == 0, result.output
    assert out.read_text(encoding="utf-8") == render_markdown(changes)


def test_a_site_without_a_manifest_is_a_bad_argument(
    pair: tuple[Path, Path, Path, Path], tmp_path: Path
) -> None:
    base_root, base_site, head_root, _ = pair
    result = CliRunner().invoke(
        app,
        [
            "report-changes",
            "--base-root", str(base_root),
            "--base-site", str(base_site),
            "--head-root", str(head_root),
            "--head-site", str(tmp_path),
        ],
    )  # fmt: skip
    assert result.exit_code == 2
    assert "has no manifest.json" in result.output
