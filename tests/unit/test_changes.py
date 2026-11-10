"""`site_changes` and `render_markdown` over synthetic sites and checkouts.

Each fixture separates the quantities a mutant would confuse (CLAUDE.md §4.30):
the gene that leaves (G_LEFT) is not the gene whose grade changes (G_GRADE), and
the manifest count that changes (`burden_rows`) is not `genes`. Strings are
asserted against literals, and a rendered assertion slices to the section it is
about (§4.19), because the same gene id appears in several sections.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
from typer.testing import CliRunner

from chd_atlas.changes import REPORT_MARKER, MirrorRowChange, render_markdown, site_changes
from chd_atlas.cli import app
from tests.unit.package_fixtures import manifest_dict

G_STAY = "HGNC:1"
G_GRADE = "HGNC:2"
G_LEFT = "HGNC:3"
G_NEW = "HGNC:4"
SHA_A = "a" * 64
SHA_B = "b" * 64


def row(**cells: Any) -> dict[str, Any]:
    base: dict[str, Any] = {
        "study": "PMID:1",
        "cohort_stratum": "all",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "inherited",
        "maf_max": None,
        "pvalue": 0.5,
    }
    base.update(cells)
    return base


def entry(gene: str, symbol: str, **cells: Any) -> dict[str, Any]:
    found: dict[str, Any] = {
        "gene": gene,
        "symbol": symbol,
        "bundle": f"genes/{gene.replace(':', '_')}.json",
        "headline_confidence": "definitive",
        "has_conflicting_evidence": False,
        "has_no_association_report": False,
    }
    found.update(cells)
    return found


def write_site(
    site: Path,
    genes: list[tuple[dict[str, Any], dict[str, Any]]],
    *,
    counts: dict[str, int] | None = None,
    files: dict[str, str] | None = None,
    schema_version: str = "2.22",
    source_commit: str | None = None,
) -> Path:
    """`genes` is (index entry, bundle) pairs; the bundle is written verbatim."""
    (site / "genes").mkdir(parents=True)
    for index_entry, bundle in genes:
        (site / index_entry["bundle"]).write_text(json.dumps(bundle), encoding="utf-8")
    (site / "genes" / "index.json").write_text(
        json.dumps({"genes": [index_entry for index_entry, _ in genes]}), encoding="utf-8"
    )
    manifest = {
        "counts": counts if counts is not None else {"genes": len(genes)},
        "files": files if files is not None else {},
        "schema_version": schema_version,
        "source_commit": source_commit,
        "status": "ok",
    }
    (site / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return site


def write_root(
    root: Path,
    burden_genes: list[str] = (),
    pin: tuple[str, str] | None = None,  # type: ignore[assignment]
    redistribution: str = "permitted",
) -> Path:
    (root / "mirrors").mkdir(parents=True)
    # Only the gene column matters to the report, which reads nothing else.
    header = "study\tgene\tsource"
    lines = [header, *(f"PMID:1\t{gene}\ttoy" for gene in burden_genes)]
    (root / "mirrors" / "burden.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    package = (
        f"\n    package:\n      version: {pin[0]}\n      manifest_sha256: {pin[1]}" if pin else ""
    )
    (root / "mirrors" / "sources.yaml").write_text(
        "sources:\n"
        "  - id: toy_burden\n    name: Toy\n    version: v1\n    retrieved_on: 2026-11-01\n"
        "    url: https://example.org\n    licence: CC-BY-4.0\n"
        f"    redistribution: {redistribution}{package}\n",
        encoding="utf-8",
    )
    return root


def write_genes(path: Path, genes: list[str]) -> None:
    """A mirror holding only its gene column, one line per row: all the report reads."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(["gene", *genes]) + "\n", encoding="utf-8")


def pair(tmp_path: Path, base_genes: list[Any], head_genes: list[Any], **head: Any) -> Any:
    base_site = write_site(tmp_path / "bs", base_genes)
    head_site = write_site(tmp_path / "hs", head_genes, **head)
    return site_changes(
        write_root(tmp_path / "br"), base_site, write_root(tmp_path / "hr"), head_site
    )


def test_entered_and_left_are_not_swapped(tmp_path: Path) -> None:
    changes = pair(
        tmp_path,
        [(entry(G_STAY, "STAY"), {}), (entry(G_LEFT, "LEFT"), {})],
        [(entry(G_STAY, "STAY"), {}), (entry(G_NEW, "NEW"), {})],
    )
    assert changes.entered == ((G_NEW, "NEW"),)
    assert changes.left == ((G_LEFT, "LEFT"),)
    assert changes.genes == ()
    text = render_markdown(changes)
    entered = text.split("### Genes entering the published set")[1].split("###")[0]
    left = text.split("### Genes leaving the published set")[1].split("###")[0]
    assert "- HGNC:4 NEW" in entered and "HGNC:3" not in entered
    assert "- HGNC:3 LEFT" in left and "HGNC:4" not in left


def test_a_grade_and_a_flag_change_are_reported_per_gene(tmp_path: Path) -> None:
    changes = pair(
        tmp_path,
        [(entry(G_GRADE, "GRADE"), {}), (entry(G_STAY, "STAY"), {}), (entry(G_LEFT, "LEFT"), {})],
        [
            (
                entry(
                    G_GRADE,
                    "GRADE",
                    headline_confidence="limited",
                    has_conflicting_evidence=True,
                    has_no_association_report=True,
                ),
                {},
            ),
            (entry(G_STAY, "STAY"), {}),
        ],
    )
    assert changes.left == ((G_LEFT, "LEFT"),)
    assert [g.gene for g in changes.genes] == [G_GRADE]
    assert changes.genes[0].lines == (
        "headline_confidence: definitive → limited",
        "has_conflicting_evidence: False → True",
        "has_no_association_report: False → True",
    )


def test_burden_rows_are_added_removed_and_changed_by_study_stratum_and_consequence(
    tmp_path: Path,
) -> None:
    kept = row(origin="de_novo")
    base = {
        "burden": [
            kept,
            row(consequence_class="lof", origin="inherited", pvalue=0.5),
            row(consequence_class="missense", origin="inherited"),
            row(consequence_class="missense", origin="de_novo"),
        ]
    }
    head = {
        "burden": [
            kept,
            row(consequence_class="lof", origin="inherited", pvalue=0.01),  # same key, new value
            row(consequence_class="missense", origin="inherited"),
            row(study="PMID:2", consequence_class="lof"),  # added
        ]
    }
    changes = pair(tmp_path, [(entry(G_GRADE, "GRADE"), base)], [(entry(G_GRADE, "GRADE"), head)])
    assert changes.genes[0].lines == (
        "burden PMID:1 / all / lof: 1 changed",
        "burden PMID:1 / all / missense: 1 removed",
        "burden PMID:2 / all / lof: 1 added",
    )


def test_a_gene_whose_rows_go_to_zero_is_named_in_bold(tmp_path: Path) -> None:
    base_site = write_site(tmp_path / "bs", [])
    head_site = write_site(tmp_path / "hs", [])
    changes = site_changes(
        write_root(tmp_path / "br", [G_GRADE, G_GRADE, G_STAY, G_STAY, G_STAY]),
        base_site,
        write_root(tmp_path / "hr", [G_STAY]),
        head_site,
    )
    # A gene going to zero comes first, although HGNC:1 sorts before HGNC:2.
    assert [(c.table, c.gene, c.before, c.after) for c in changes.mirror_rows] == [
        ("burden", G_GRADE, 2, 0),
        ("burden", G_STAY, 3, 1),
    ]
    section = render_markdown(changes).split("### Mirror rows per gene")[1].split("###")[0]
    assert "- burden HGNC:2: 2 → **0**" in section
    assert "- burden HGNC:1: 3 → 1\n" in section


def test_a_gene_going_to_zero_is_never_folded_into_a_capped_table(tmp_path: Path) -> None:
    """Sixty GenCC genes lose a row and a sixty-first loses its only one; TBX5's
    nine `profiles` rows go to none.

    By (table, gene) `profiles` sorts after `gencc_submissions`, so a single cap
    over the section would fold TBX5 into the remainder (the reviewer's
    measurement), and HGNC:2999 sorts after every other GenCC gene, so a cap per
    table alone would fold it. Zeros come first and are never capped; the rest
    is capped per table, so `profiles`' own non-zero change survives too.
    """
    many = [f"HGNC:{2000 + i}" for i in range(60)]
    base_root, head_root = write_root(tmp_path / "br"), write_root(tmp_path / "hr")
    write_genes(base_root / "mirrors" / "gencc_submissions.tsv", [*many, *many, "HGNC:2999"])
    write_genes(head_root / "mirrors" / "gencc_submissions.tsv", many)
    shard = Path("mirrors") / "profiles" / "E-MTAB-6814.tsv"
    write_genes(base_root / shard, ["HGNC:11604"] * 9 + ["HGNC:7"] * 2)
    write_genes(head_root / shard, ["HGNC:7"])

    changes = site_changes(
        base_root, write_site(tmp_path / "bs", []), head_root, write_site(tmp_path / "hs", [])
    )

    assert changes.mirror_rows[:2] == (
        MirrorRowChange("gencc_submissions", "HGNC:2999", 1, 0),
        MirrorRowChange("profiles", "HGNC:11604", 9, 0),
    )
    section = render_markdown(changes).split("### Mirror rows per gene")[1].split("###")[0]
    lines = [line for line in section.splitlines() if line]
    assert lines[:2] == [
        "- gencc_submissions HGNC:2999: 1 → **0**",
        "- profiles HGNC:11604: 9 → **0**",
    ]
    assert sum(line.startswith("- gencc_submissions ") for line in lines[2:]) == 50
    assert "… and 10 more in gencc_submissions" in lines
    assert "- profiles HGNC:7: 2 → 1" in lines
    assert not any(line.startswith("… and") and "profiles" in line for line in lines)


def test_a_base_bundle_without_literature_is_read_as_empty(tmp_path: Path) -> None:
    head = {
        "burden": [row()],
        "literature": {"findings": [{}, {}], "pending_approval": 1},
        "expression_profile": {"datasets": [{"dataset": "E-MTAB-6814"}]},
        "cell_type_expression": {"datasets": [{"dataset": "farah2024_hoc"}]},
    }
    # The base is a bundle built before any of these keys existed.
    changes = pair(tmp_path, [(entry(G_GRADE, "GRADE"), {})], [(entry(G_GRADE, "GRADE"), head)])
    assert changes.genes[0].lines == (
        "burden PMID:1 / all / lof: 1 added",
        "expression_profile datasets added: E-MTAB-6814",
        "cell_type_expression datasets added: farah2024_hoc",
        "literature findings: 0 → 2",
        "literature pending approval: 0 → 1",
    )


def test_identical_sites_report_no_published_file_changes(tmp_path: Path) -> None:
    files = {"genes/index.json": f"sha256:{SHA_A}"}
    genes = [(entry(G_STAY, "STAY"), {})]
    base_site = write_site(tmp_path / "bs", genes, files=files, source_commit="1" * 40)
    # source_commit differs between any two commits and says nothing about data.
    head_site = write_site(tmp_path / "hs", genes, files=files, source_commit="2" * 40)
    changes = site_changes(
        write_root(tmp_path / "br"), base_site, write_root(tmp_path / "hr"), head_site
    )
    assert changes.files_changed == ()
    assert changes.files_total == 2
    assert changes.counts == ()
    text = render_markdown(changes)
    assert "No published file changes." in text
    # Both commits head the report, so a stale comment can be told from a current one.
    assert f"\nBase built from {'1' * 40}; head built from {'2' * 40}.\n" in text
    assert text.index("Base built from") < text.index("No published file changes.")


def test_a_count_and_a_file_change_are_listed_and_genes_is_not(tmp_path: Path) -> None:
    genes = [(entry(G_STAY, "STAY"), {})]
    base_site = write_site(
        tmp_path / "bs",
        genes,
        counts={"genes": 1, "burden_rows": 10, "old_key": 5},
        files={
            "a.json": f"sha256:{SHA_A}",
            "b.json": f"sha256:{SHA_A}",
            "gone.json": f"sha256:{SHA_A}",
        },
    )
    head_site = write_site(
        tmp_path / "hs",
        genes,
        counts={"genes": 1, "burden_rows": 8, "new_key": 2},
        files={
            "a.json": f"sha256:{SHA_A}",
            "b.json": f"sha256:{SHA_B}",
            "new.json": f"sha256:{SHA_B}",
        },
    )
    changes = site_changes(
        write_root(tmp_path / "br"), base_site, write_root(tmp_path / "hr"), head_site
    )
    assert changes.counts == (
        ("burden_rows", 10, 8),
        ("new_key", None, 2),
        ("old_key", 5, None),
    )
    assert changes.files_changed == ("b.json", "gone.json", "manifest.json", "new.json")
    assert changes.files_total == 5  # a, b, gone, new + manifest.json
    text = render_markdown(changes)
    assert "- burden_rows: 10 → 8" in text.split("### Manifest counts")[1].split("###")[0]
    assert "4 of 5 files changed" in text


def test_a_changed_pin_reports_its_dropped_rows(tmp_path: Path) -> None:
    head_root = write_root(tmp_path / "hr", pin=("2026-11-10.1", SHA_B))
    (head_root / "mirrors" / "packages").mkdir()
    (head_root / "mirrors" / "packages" / "toy_burden.json").write_text(
        json.dumps(
            manifest_dict(
                dropped=[{"reason": "symbol not in HGNC", "count": 2, "examples": ["FOO1"]}]
            )
        ),
        encoding="utf-8",
    )
    changes = site_changes(
        write_root(tmp_path / "br", pin=("2026-10-01.1", SHA_A)),
        write_site(tmp_path / "bs", []),
        head_root,
        write_site(tmp_path / "hs", []),
    )
    [package] = changes.packages
    assert package.source_id == "toy_burden"
    assert package.before == "2026-10-01.1 (sha256:aaaaaaaaaaaa)"
    assert package.after == "2026-11-10.1 (sha256:bbbbbbbbbbbb)"
    assert [(d.reason, d.count) for d in package.dropped] == [("symbol not in HGNC", 2)]
    assert package.input_rows == 10
    section = render_markdown(changes).split("### Packages")[1]
    assert "dropped 2: symbol not in HGNC (e.g. FOO1)" in section
    assert "10 input rows" in section


def test_a_repin_shows_the_licence_and_a_redistribution_change_on_its_own_line(
    tmp_path: Path,
) -> None:
    head_root = write_root(
        tmp_path / "hr", pin=("v2", SHA_B), redistribution="permitted_with_attribution"
    )
    (head_root / "mirrors" / "packages").mkdir()
    (head_root / "mirrors" / "packages" / "toy_burden.json").write_text(
        json.dumps(manifest_dict(licence_evidence="Supplement footer, CC BY 4.0")),
        encoding="utf-8",
    )
    changes = site_changes(
        write_root(tmp_path / "br", pin=("v1", SHA_A)),
        write_site(tmp_path / "bs", []),
        head_root,
        write_site(tmp_path / "hs", []),
    )
    [package] = changes.packages
    assert (package.licence_before, package.licence_after) == ("CC-BY-4.0", "CC-BY-4.0")
    assert (package.redistribution_before, package.redistribution_after) == (
        "permitted",
        "permitted_with_attribution",
    )
    assert package.licence_evidence == "Supplement footer, CC BY 4.0"
    lines = render_markdown(changes).split("### Packages")[1].split("###")[0].splitlines()
    assert "  - licence: CC-BY-4.0" in lines
    assert "  - redistribution: permitted → permitted_with_attribution" in lines
    assert "  - licence evidence (head manifest): Supplement footer, CC BY 4.0" in lines


def test_an_unreadable_head_manifest_is_named_not_raised(tmp_path: Path) -> None:
    changes = site_changes(
        write_root(tmp_path / "br"),
        write_site(tmp_path / "bs", []),
        write_root(tmp_path / "hr", pin=("v2", SHA_B)),  # no mirrors/packages/ file
        write_site(tmp_path / "hs", []),
    )
    [package] = changes.packages
    assert package.before is None and package.dropped == () and package.input_rows is None
    assert "head manifest unreadable" in render_markdown(changes).split("### Packages")[1]


def _n_of(section: str, n: int, tmp_path: Path) -> tuple[Any, str]:
    """A pair of sites with ``n`` items in ``section``, and that section's heading."""
    genes = [entry(f"HGNC:{1000 + i}", f"G{i:03d}") for i in range(n)]
    if section == "entering":
        return pair(tmp_path, [], [(g, {}) for g in genes]), "Genes entering the published set"
    if section == "leaving":
        return pair(tmp_path, [(g, {}) for g in genes], []), "Genes leaving the published set"
    if section == "genes changed":
        limited = [(dict(g, headline_confidence="limited"), {}) for g in genes]
        return pair(tmp_path, [(g, {}) for g in genes], limited), "Genes changed"
    base_site = write_site(tmp_path / "bs", [], counts={f"k{i:03d}": 1 for i in range(n)})
    head_site = write_site(tmp_path / "hs", [], counts={f"k{i:03d}": 2 for i in range(n)})
    changes = site_changes(
        write_root(tmp_path / "br"), base_site, write_root(tmp_path / "hr"), head_site
    )
    return changes, "Manifest counts"


@pytest.mark.parametrize("section", ["entering", "leaving", "genes changed", "counts"])
@pytest.mark.parametrize("n", [50, 51, 60])
def test_lists_are_capped_at_fifty_with_the_remainder_counted(
    tmp_path: Path, n: int, section: str
) -> None:
    changes, heading = _n_of(section, n, tmp_path)
    text = render_markdown(changes).split(f"### {heading}\n")[1].split("###")[0]
    # Top-level items only: a changed gene's own lines are indented under it.
    assert sum(line.startswith("- ") for line in text.splitlines()) == min(n, 50)
    assert ("… and " in text) == (n > 50)
    if n > 50:
        assert f"… and {n - 50} more" in text


def test_the_report_opens_with_its_marker(tmp_path: Path) -> None:
    text = render_markdown(pair(tmp_path, [], []))
    assert text.startswith("<!-- chd-atlas:site-change-report -->\n")
    assert REPORT_MARKER == "<!-- chd-atlas:site-change-report -->"


@pytest.mark.parametrize(
    ("unreadable", "label"),
    [
        ("study\tsource\nPMID:1\ttoy\n", "lacking the gene column"),
        ("", "empty"),
    ],
)
def test_an_unreadable_mirror_is_named_not_reported_as_zero(
    tmp_path: Path, unreadable: str, label: str
) -> None:
    base_root = write_root(tmp_path / "br", [G_GRADE, G_GRADE, G_STAY])
    (base_root / "mirrors" / "burden.tsv").write_text(unreadable, encoding="utf-8")
    head_root = write_root(tmp_path / "hr", [G_STAY])
    changes = site_changes(
        base_root, write_site(tmp_path / "bs", []), head_root, write_site(tmp_path / "hs", [])
    )
    assert changes.mirror_rows == ()
    assert changes.unreadable_mirrors == (("burden", "base"),)
    section = render_markdown(changes).split("### Mirror rows per gene")[1].split("###")[0]
    assert "- burden: unreadable on base, rows not compared" in section
    assert "→ **0**" not in section, label


def test_an_unwritable_out_exits_2_with_its_message(tmp_path: Path) -> None:
    """As `build` treats an unwritable --out: exit 2, `error: …`, no traceback."""
    out = tmp_path / "missing" / "report.md"
    paths = {
        "--base-root": write_root(tmp_path / "br"),
        "--base-site": write_site(tmp_path / "bs", []),
        "--head-root": write_root(tmp_path / "hr"),
        "--head-site": write_site(tmp_path / "hs", []),
        "--out": out,
    }
    args = ["report-changes", *(item for flag, path in paths.items() for item in (flag, str(path)))]

    result = CliRunner().invoke(app, args)

    assert result.exit_code == 2
    assert result.stdout.startswith(f"error: could not write the report to {out}: ")
    assert not out.parent.exists()
