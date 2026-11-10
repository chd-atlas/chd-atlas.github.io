"""What a change does to the published site: the report an owner reads instead of a row diff.

It compares two *built* sites and the two checkouts they were built from, so CI
can build the base with the base's own code. Everything is read as the files
say, with no knowledge of which build wrote them: a bundle from older code may
lack a key (``literature`` before schema 2.21), and a missing key is read as
empty rather than raised on.

Every list in `SiteChanges` is sorted by the caller-visible key (gene id, count
key, path, source id), so two runs over the same pair render the same bytes.
`mirror_rows` is the one list with a different order: every (table, gene) whose
rows went to zero comes first, then the rest, each part sorted by (table, gene),
because a gene losing all its rows is the change the report must never fold
away (spec §4.4: "a gene going from 9 rows to 0 is named").
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import polars as pl

from chd_atlas.models.package import DroppedRows, ManifestError, load_manifest
from chd_atlas.tables import BURDEN, PACKAGES_DIR, TABLE_SCHEMAS, mirror_paths
from chd_atlas.validate.sources import PackagePin, load_sources

REPORT_MARKER: Final = "<!-- chd-atlas:site-change-report -->"

# A rendered list stops here and counts the remainder: genes entering, leaving
# and changed (capped by gene, so a gene's own lines are never cut), manifest
# counts, packages and changed files. Mirror rows are capped per table, and a
# gene whose rows went to zero is never capped. A changed gene's lines, a
# package's licence and dropped lines, and the unreadable-mirror lines are not
# capped.
LIST_CAP: Final = 50

_GRADE_FIELDS: Final = (
    "headline_confidence",
    "has_conflicting_evidence",
    "has_no_association_report",
)
# The full burden row key is `BURDEN.sort_key`; a row is the same row when all
# of it but the gene matches, and the report groups by these three of them.
_ROW_KEY: Final = tuple(column for column in BURDEN.sort_key if column != "gene")
_GROUP_KEY: Final = ("study", "cohort_stratum", "consequence_class")


@dataclass(frozen=True)
class MirrorRowChange:
    table: str
    gene: str
    before: int
    after: int


@dataclass(frozen=True)
class PackageChange:
    source_id: str
    before: str | None
    after: str | None
    dropped: tuple[DroppedRows, ...]
    input_rows: int | None
    # Set when the head pin's committed manifest could not be read; `dropped`,
    # `input_rows` and `licence_evidence` are then unknown, and the report says so.
    manifest_error: str | None = None
    # From the base and head `sources.yaml` entries; None where a side has no
    # entry for the source.
    licence_before: str | None = None
    licence_after: str | None = None
    redistribution_before: str | None = None
    redistribution_after: str | None = None
    # The head pin's committed manifest's `licence_evidence`.
    licence_evidence: str | None = None


@dataclass(frozen=True)
class GeneChange:
    gene: str
    symbol: str
    lines: tuple[str, ...]


@dataclass(frozen=True)
class SiteChanges:
    schema_version: tuple[str, str]
    entered: tuple[tuple[str, str], ...]
    left: tuple[tuple[str, str], ...]
    genes: tuple[GeneChange, ...]
    counts: tuple[tuple[str, object, object], ...]
    files_changed: tuple[str, ...]
    files_total: int
    mirror_rows: tuple[MirrorRowChange, ...]
    packages: tuple[PackageChange, ...]
    # (table, side) pairs whose mirror could not be read; those tables are
    # left out of `mirror_rows` rather than counted as zero.
    unreadable_mirrors: tuple[tuple[str, str], ...] = ()
    # Each manifest's `source_commit` (base, head); None where it is null.
    source_commits: tuple[str | None, str | None] = (None, None)


def _read_json(path: Path) -> dict[str, Any]:
    loaded = json.loads(path.read_text(encoding="utf-8"))
    return loaded if isinstance(loaded, dict) else {}


def _as_dict(value: object) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _as_list(value: object) -> list[Any]:
    return value if isinstance(value, list) else []


def _index(site: Path) -> dict[str, dict[str, Any]]:
    entries = _as_list(_read_json(site / "genes" / "index.json").get("genes"))
    return {str(entry["gene"]): entry for entry in entries if isinstance(entry, dict)}


def _bundle(site: Path, entry: dict[str, Any]) -> dict[str, Any]:
    return _read_json(site / str(entry["bundle"]))


def _show(value: object) -> str:
    return "none" if value is None else str(value)


def _burden_lines(base: dict[str, Any], head: dict[str, Any]) -> list[str]:
    def keyed(bundle: dict[str, Any]) -> dict[tuple[Any, ...], dict[str, Any]]:
        rows = [row for row in _as_list(bundle.get("burden")) if isinstance(row, dict)]
        return {tuple(row.get(column) for column in _ROW_KEY): row for row in rows}

    before, after = keyed(base), keyed(head)

    def group(key: tuple[Any, ...]) -> tuple[str, str, str]:
        found = dict(zip(_ROW_KEY, key, strict=True))
        return tuple(_show(found[column]) for column in _GROUP_KEY)  # type: ignore[return-value]

    added = Counter(group(key) for key in after.keys() - before.keys())
    removed = Counter(group(key) for key in before.keys() - after.keys())
    changed = Counter(
        group(key) for key in before.keys() & after.keys() if before[key] != after[key]
    )
    lines = []
    for label in sorted(added.keys() | removed.keys() | changed.keys()):
        parts = [
            f"{count} {word}"
            for word, count in (
                ("added", added[label]),
                ("removed", removed[label]),
                ("changed", changed[label]),
            )
            if count
        ]
        lines.append(f"burden {' / '.join(label)}: {', '.join(parts)}")
    return lines


def _datasets(bundle: dict[str, Any], key: str) -> set[str]:
    datasets = _as_list(_as_dict(bundle.get(key)).get("datasets"))
    return {
        str(item["dataset"]) for item in datasets if isinstance(item, dict) and "dataset" in item
    }


def _gene_lines(
    base_entry: dict[str, Any],
    head_entry: dict[str, Any],
    base: dict[str, Any],
    head: dict[str, Any],
) -> list[str]:
    lines = [
        f"{field}: {_show(base_entry.get(field))} → {_show(head_entry.get(field))}"
        for field in _GRADE_FIELDS
        if base_entry.get(field) != head_entry.get(field)
    ]
    lines += _burden_lines(base, head)
    for key in ("expression_profile", "cell_type_expression"):
        was, now = _datasets(base, key), _datasets(head, key)
        for word, names in (("added", now - was), ("removed", was - now)):
            if names:
                lines.append(f"{key} datasets {word}: {', '.join(sorted(names))}")
    was_lit, now_lit = _as_dict(base.get("literature")), _as_dict(head.get("literature"))
    found = (len(_as_list(was_lit.get("findings"))), len(_as_list(now_lit.get("findings"))))
    if found[0] != found[1]:
        lines.append(f"literature findings: {found[0]} → {found[1]}")
    pending = (was_lit.get("pending_approval", 0), now_lit.get("pending_approval", 0))
    if pending[0] != pending[1]:
        lines.append(f"literature pending approval: {pending[0]} → {pending[1]}")
    return lines


def _manifest_files(manifest: dict[str, Any]) -> dict[str, str]:
    return {str(k): str(v) for k, v in _as_dict(manifest.get("files")).items()}


def _without_commit(manifest: dict[str, Any]) -> dict[str, Any]:
    return {k: v for k, v in manifest.items() if k != "source_commit"}


def _counts(base: dict[str, Any], head: dict[str, Any]) -> tuple[tuple[str, object, object], ...]:
    was, now = _as_dict(base.get("counts")), _as_dict(head.get("counts"))
    return tuple(
        (key, was.get(key), now.get(key))
        for key in sorted(was.keys() | now.keys())
        if was.get(key) != now.get(key)
    )


# What polars raises for a TSV that cannot give the gene column, measured on
# polars 1.43.1 by reading small files with the call below: the header lacks
# the column, including when a quote opening a header cell swallows the rest of
# the file into one name (ColumnNotFoundError); the file is empty
# (NoDataError); a gene cell opens a quote that never closes (ComputeError).
# Two shapes raise nothing. A header-only file with the gene column reads as
# zero rows, which is what it holds. A quote opening a cell of another column
# swallows the rest of the file into that cell, so that row and every row
# after it read as one row with a null gene. This reader does not detect that;
# no committed mirror contains a quote character (measured 2026-11-10).
_UNREADABLE: Final = (
    pl.exceptions.ColumnNotFoundError,
    pl.exceptions.NoDataError,
    pl.exceptions.ComputeError,
)


def _mirror_counts(
    root: Path, relative: dict[Path, str]
) -> tuple[Counter[tuple[str, str]], set[str]]:
    """Rows per (table, gene), the gene column read as text and nothing else.

    Also the tables with a file that could not be read. Those are not counted
    as zero: a file that gave no gene column says nothing about its genes.
    """
    found: Counter[tuple[str, str]] = Counter()
    unreadable: set[str] = set()
    for rel, table in relative.items():
        path = root / rel
        gene_column = TABLE_SCHEMAS[table].gene_column
        if gene_column is None or not path.is_file():
            continue
        # Text, one column: a column added or retyped between base and head
        # must not fail the report.
        try:
            frame = pl.read_csv(path, separator="\t", infer_schema_length=0, columns=[gene_column])
        except _UNREADABLE:
            unreadable.add(table)
            continue
        for gene in frame[gene_column].drop_nulls():
            found[(table, gene)] += 1
    return found, unreadable


def _mirror_changes(
    base_root: Path, head_root: Path
) -> tuple[tuple[MirrorRowChange, ...], tuple[tuple[str, str], ...]]:
    relative: dict[Path, str] = {}
    for root in (head_root, base_root):
        for path, table in mirror_paths(root):
            relative[path.relative_to(root)] = table
    was, was_bad = _mirror_counts(base_root, relative)
    now, now_bad = _mirror_counts(head_root, relative)
    skipped = was_bad | now_bad
    changes = tuple(
        sorted(
            (
                MirrorRowChange(table, gene, was[(table, gene)], now[(table, gene)])
                for table, gene in was.keys() | now.keys()
                if was[(table, gene)] != now[(table, gene)] and table not in skipped
            ),
            key=lambda change: (change.after != 0, change.table, change.gene),
        )
    )
    unreadable = tuple(sorted([(t, "base") for t in was_bad] + [(t, "head") for t in now_bad]))
    return changes, unreadable


def _pin_text(pin: PackagePin | None) -> str | None:
    return None if pin is None else f"{pin.version} (sha256:{pin.manifest_sha256[:12]})"


def _package_changes(base_root: Path, head_root: Path) -> tuple[PackageChange, ...]:
    was = load_sources(base_root)[0].by_id
    now = load_sources(head_root)[0].by_id
    changes = []
    for source_id in sorted(was.keys() | now.keys()):
        old, new = was.get(source_id), now.get(source_id)
        before = None if old is None else old.package
        after = None if new is None else new.package
        if before == after:
            continue
        dropped: tuple[DroppedRows, ...] = ()
        input_rows: int | None = None
        evidence: str | None = None
        error: str | None = None
        if after is not None:
            try:
                manifest = load_manifest(head_root / "mirrors" / PACKAGES_DIR / f"{source_id}.json")
            except ManifestError as exc:
                error = str(exc)
            else:
                dropped = tuple(manifest.dropped)
                input_rows = manifest.input_counts.rows
                evidence = manifest.licence_evidence
        changes.append(
            PackageChange(
                source_id,
                _pin_text(before),
                _pin_text(after),
                dropped,
                input_rows,
                error,
                licence_before=None if old is None else old.licence,
                licence_after=None if new is None else new.licence,
                redistribution_before=None if old is None else old.redistribution.value,
                redistribution_after=None if new is None else new.redistribution.value,
                licence_evidence=evidence,
            )
        )
    return tuple(changes)


def site_changes(base_root: Path, base_site: Path, head_root: Path, head_site: Path) -> SiteChanges:
    """Compare two built sites, and the two checkouts they came from."""
    base_manifest, head_manifest = (
        _read_json(base_site / "manifest.json"),
        _read_json(head_site / "manifest.json"),
    )
    base_files, head_files = _manifest_files(base_manifest), _manifest_files(head_manifest)
    changed = {
        path
        for path in base_files.keys() | head_files.keys()
        if base_files.get(path) != head_files.get(path)
    }
    # manifest.json cannot carry its own checksum, and `source_commit` differs
    # between any two commits by construction, so it is listed only when
    # something else in it differs.
    if _without_commit(base_manifest) != _without_commit(head_manifest):
        changed.add("manifest.json")

    base_index, head_index = _index(base_site), _index(head_site)
    genes = []
    for gene in sorted(base_index.keys() & head_index.keys()):
        lines = _gene_lines(
            base_index[gene],
            head_index[gene],
            _bundle(base_site, base_index[gene]),
            _bundle(head_site, head_index[gene]),
        )
        if lines:
            genes.append(GeneChange(gene, str(head_index[gene].get("symbol", "")), tuple(lines)))

    def named(index: dict[str, dict[str, Any]], ids: set[str]) -> tuple[tuple[str, str], ...]:
        return tuple((gene, str(index[gene].get("symbol", ""))) for gene in sorted(ids))

    mirror_rows, unreadable_mirrors = _mirror_changes(base_root, head_root)
    return SiteChanges(
        schema_version=(
            str(base_manifest.get("schema_version")),
            str(head_manifest.get("schema_version")),
        ),
        entered=named(head_index, head_index.keys() - base_index.keys()),
        left=named(base_index, base_index.keys() - head_index.keys()),
        genes=tuple(genes),
        counts=_counts(base_manifest, head_manifest),
        files_changed=tuple(sorted(changed)),
        files_total=len(base_files.keys() | head_files.keys()) + 1,
        mirror_rows=mirror_rows,
        packages=_package_changes(base_root, head_root),
        unreadable_mirrors=unreadable_mirrors,
        source_commits=(
            _commit_or_none(base_manifest.get("source_commit")),
            _commit_or_none(head_manifest.get("source_commit")),
        ),
    )


def _commit_or_none(value: object) -> str | None:
    return value if isinstance(value, str) else None


def _capped(items: list[str]) -> list[str]:
    shown = items[:LIST_CAP]
    if len(items) > LIST_CAP:
        shown.append(f"… and {len(items) - LIST_CAP} more")
    return shown


def _section(title: str, items: list[str], *, cap: bool = True) -> list[str]:
    if not items:
        return []
    return ["", f"### {title}", "", *(_capped(items) if cap else items)]


def _mirror_lines(changes: SiteChanges) -> list[str]:
    """Unreadable tables, then every gene going to zero, then the rest capped per table."""
    lines = [
        f"- {table}: unreadable on {side}, rows not compared"
        for table, side in changes.unreadable_mirrors
    ]
    lines += [
        f"- {c.table} {c.gene}: {c.before} → **0**" for c in changes.mirror_rows if c.after == 0
    ]
    rest: dict[str, list[str]] = {}
    for c in changes.mirror_rows:
        if c.after != 0:
            rest.setdefault(c.table, []).append(f"- {c.table} {c.gene}: {c.before} → {c.after}")
    for table in sorted(rest):
        lines += rest[table][:LIST_CAP]
        if len(rest[table]) > LIST_CAP:
            lines.append(f"… and {len(rest[table]) - LIST_CAP} more in {table}")
    return lines


def _commit_text(value: str | None) -> str:
    return "an unrecorded commit" if value is None else value


def render_markdown(changes: SiteChanges) -> str:
    """One Markdown comment. Only non-empty sections are printed."""
    base_commit, head_commit = changes.source_commits
    lines = [
        REPORT_MARKER,
        "",
        "## Site-change report",
        "",
        f"Base built from {_commit_text(base_commit)}; "
        f"head built from {_commit_text(head_commit)}.",
    ]
    before, after = changes.schema_version
    if before != after:
        lines += ["", f"Schema version: {before} → {after}"]
    lines += _section(
        "Genes entering the published set", [f"- {g} {s}" for g, s in changes.entered]
    )
    lines += _section("Genes leaving the published set", [f"- {g} {s}" for g, s in changes.left])
    # Capped by gene, not by line, so a gene's lines are never cut in half.
    lines += _section_by_gene(changes)
    lines += _section(
        "Manifest counts",
        [f"- {key}: {_show(was)} → {_show(now)}" for key, was, now in changes.counts],
    )
    # Capped per table inside `_mirror_lines`, never as a whole: see LIST_CAP.
    lines += _section("Mirror rows per gene", _mirror_lines(changes), cap=False)
    lines += _section("Packages", [_package_text(p) for p in changes.packages])
    if changes.files_changed:
        lines += _section(
            f"{len(changes.files_changed)} of {changes.files_total} files changed",
            [f"- {path}" for path in changes.files_changed],
        )
    else:
        lines += ["", "No published file changes."]
    return "\n".join(lines) + "\n"


def _section_by_gene(changes: SiteChanges) -> list[str]:
    if not changes.genes:
        return []
    out = ["", "### Genes changed", ""]
    for gene in changes.genes[:LIST_CAP]:
        out.append(f"- {gene.gene} {gene.symbol}")
        out += [f"  - {line}" for line in gene.lines]
    if len(changes.genes) > LIST_CAP:
        out.append(f"… and {len(changes.genes) - LIST_CAP} more")
    return out


def _package_text(change: PackageChange) -> str:
    text = f"- {change.source_id}: {_show(change.before)} → {_show(change.after)}"
    if change.manifest_error is not None:
        text += f"; head manifest unreadable ({change.manifest_error})"
    elif change.input_rows is not None:
        text += f"; {change.input_rows} input rows"
    # Each on its own line, so a redistribution change cannot hide in a long one.
    for label, was, now in (
        ("licence", change.licence_before, change.licence_after),
        ("redistribution", change.redistribution_before, change.redistribution_after),
    ):
        shown = _show(now) if was == now else f"{_show(was)} → {_show(now)}"
        text += f"\n  - {label}: {shown}"
    if change.licence_evidence is not None:
        text += f"\n  - licence evidence (head manifest): {change.licence_evidence}"
    for dropped in change.dropped:
        examples = f" (e.g. {', '.join(dropped.examples)})" if dropped.examples else ""
        text += f"\n  - dropped {dropped.count}: {dropped.reason}{examples}"
    return text
