"""What a change does to the published site: the report an owner reads instead of a row diff.

It compares two *built* sites and the two checkouts they were built from, so CI
can build the base with the base's own code. Everything is read as the files
say, with no knowledge of which build wrote them: a bundle from older code may
lack a key (``literature`` before schema 2.21), and a missing key is read as
empty rather than raised on.

Every list in `SiteChanges` is sorted by the caller-visible key (gene id, count
key, path, source id), so two runs over the same pair render the same bytes.
"""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

import polars as pl

from chd_atlas.models.package import DroppedRows, ManifestError, load_manifest
from chd_atlas.tables import BURDEN, TABLE_SCHEMAS, mirror_paths
from chd_atlas.validate.sources import PackagePin, load_sources

REPORT_MARKER: Final = "<!-- chd-atlas:site-change-report -->"

# Every rendered list stops here and counts the remainder.
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
    # Set when the head pin's committed manifest could not be read; `dropped`
    # and `input_rows` are then unknown, and the report says so.
    manifest_error: str | None = None


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


# What polars raises for a TSV that cannot give the gene column: the column is
# absent (ColumnNotFoundError, measured also for a header-only file and for an
# unterminated quote), the file is empty (NoDataError), or it will not parse
# (ComputeError).
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
        MirrorRowChange(table, gene, was[(table, gene)], now[(table, gene)])
        for table, gene in sorted(was.keys() | now.keys())
        if was[(table, gene)] != now[(table, gene)] and table not in skipped
    )
    unreadable = tuple(sorted([(t, "base") for t in was_bad] + [(t, "head") for t in now_bad]))
    return changes, unreadable


def _pin_text(pin: PackagePin | None) -> str | None:
    return None if pin is None else f"{pin.version} (sha256:{pin.manifest_sha256[:12]})"


def _package_changes(base_root: Path, head_root: Path) -> tuple[PackageChange, ...]:
    was = {source.id: source.package for source in load_sources(base_root)[0].sources}
    now = {source.id: source.package for source in load_sources(head_root)[0].sources}
    changes = []
    for source_id in sorted(was.keys() | now.keys()):
        before, after = was.get(source_id), now.get(source_id)
        if before == after:
            continue
        dropped: tuple[DroppedRows, ...] = ()
        input_rows: int | None = None
        error: str | None = None
        if after is not None:
            try:
                manifest = load_manifest(head_root / "mirrors" / "packages" / f"{source_id}.json")
            except ManifestError as exc:
                error = str(exc)
            else:
                dropped = tuple(manifest.dropped)
                input_rows = manifest.input_counts.rows
        changes.append(
            PackageChange(
                source_id, _pin_text(before), _pin_text(after), dropped, input_rows, error
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
    )


def _capped(items: list[str]) -> list[str]:
    shown = items[:LIST_CAP]
    if len(items) > LIST_CAP:
        shown.append(f"… and {len(items) - LIST_CAP} more")
    return shown


def _section(title: str, items: list[str]) -> list[str]:
    return ["", f"### {title}", "", *_capped(items)] if items else []


def render_markdown(changes: SiteChanges) -> str:
    """One Markdown comment. Only non-empty sections are printed."""
    lines = [REPORT_MARKER, "", "## Site-change report"]
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
    lines += _section(
        "Mirror rows per gene",
        [
            f"- {table}: unreadable on {side}, rows not compared"
            for table, side in changes.unreadable_mirrors
        ]
        + [
            f"- {c.table} {c.gene}: {c.before} → {'**0**' if c.after == 0 else c.after}"
            for c in changes.mirror_rows
        ],
    )
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
        return f"{text}; head manifest unreadable ({change.manifest_error})"
    if change.input_rows is not None:
        text += f"; {change.input_rows} input rows"
    for dropped in change.dropped:
        examples = f" (e.g. {', '.join(dropped.examples)})" if dropped.examples else ""
        text += f"\n  - dropped {dropped.count}: {dropped.reason}{examples}"
    return text
