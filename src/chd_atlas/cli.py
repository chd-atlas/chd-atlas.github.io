from __future__ import annotations

import hashlib
import json
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import date
from pathlib import Path
from typing import Any, Final

import typer

from chd_atlas.build.questions import (
    QUESTION_LABELS,
    SETTLED_LABELS,
    SETTLED_ORDER,
    open_questions,
)
from chd_atlas.build.runner import BuildRefused, build_site
from chd_atlas.changes import render_markdown, site_changes
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.fs import write_bytes_atomically
from chd_atlas.gene_literature_diff import diff_findings
from chd_atlas.gene_literature_io import (
    LiteratureFileError,
    literature_path,
    parse_findings_file,
    parse_literature,
    with_findings_updated,
    with_synthesis_approval,
)
from chd_atlas.ingest.apply import (
    REGISTRY_INPUTS,
    IngestAllResult,
    IngestResult,
    ingest_all,
    ingest_package,
    package_pins,
    pinned_package_dir,
)
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.issues import ValidationIssue
from chd_atlas.literature_query import (
    TargetError,
    literature_frame,
    literature_target,
    resolve_symbol,
)
from chd_atlas.models.assertion import InTextLocator, Locator
from chd_atlas.models.gene_literature import (
    Finding,
    GeneLiteratureFile,
    Variant,
    finding_digest,
    shown_digest,
    synthesis_digest,
)
from chd_atlas.models.literature import Publication
from chd_atlas.models.review import ReviewStatus
from chd_atlas.models.review_decisions import Decision, ReviewDecisions, ReviewerRole
from chd_atlas.schema_export import TABLE_SCHEMA_DIR, export_schemas, export_table_schemas
from chd_atlas.validate.gene_synthesis import validate_gene_synthesis
from chd_atlas.validate.runner import validate_repository
from chd_atlas.vocab import SettledBy

app = typer.Typer(help="CHD Atlas curation tooling.", no_args_is_help=True)
schemas_app = typer.Typer(help="JSON Schema utilities.", no_args_is_help=True)
literature_app = typer.Typer(help="Literature review helpers.", no_args_is_help=True)
app.add_typer(schemas_app, name="schemas")
app.add_typer(literature_app, name="literature")


@app.command()
def validate(
    root: Path = typer.Option(Path("."), help="Repository root to validate."),
) -> None:
    """Validate the curation corpus and mirror tables."""
    # Every validator reports its own target as missing, so a --root that does
    # not exist produces output byte-identical to a real but empty repository:
    # a confident list of content errors about a repository that was never
    # read. Exit 2 rather than 1 so CI can distinguish "you pointed me at the
    # wrong place" from "this repository has errors".
    if not root.is_dir():
        typer.echo(f"error: --root {root} is not a directory")
        raise typer.Exit(code=2)

    report = validate_repository(root)
    typer.echo(report.render())
    raise typer.Exit(code=0 if report.ok else 1)


@app.command()
def build(
    root: Path = typer.Option(Path("."), help="Repository root to build from."),
    out: Path = typer.Option(Path("dist"), help="Directory to write the site into."),
) -> None:
    """Build the published data API into a directory."""
    # Checked before the build, for the same reason `validate` checks it: a
    # --root that does not exist otherwise reaches the validator, which reports
    # every target as missing and refuses — so a typo produces a confident list
    # of content errors about a repository that was never read. Exit 2 rather
    # than 1 so CI can tell "you pointed me at the wrong place" from "this
    # repository has errors"; the first is a pipeline bug, the second a
    # curation one, and they are fixed by different people.
    if not root.is_dir():
        typer.echo(f"error: --root {root} is not a directory")
        raise typer.Exit(code=2)

    # `--out` may legitimately not exist, since the build creates it — but if
    # something is already there it has to be a directory. `is_symlink()` is
    # tested first because `exists()` follows the link, so a dangling symlink
    # answers False to both `exists()` and `is_dir()` and would slip past a
    # check written the obvious way.
    #
    # Checked here rather than left to fail during the build for the same
    # reason `--root` is: an argument that cannot work should not cost a second
    # and a half of validation before saying so.
    if (out.is_symlink() or out.exists()) and not out.is_dir():
        typer.echo(f"error: --out {out} exists and is not a directory")
        raise typer.Exit(code=2)

    try:
        written = build_site(root, out)
    except BuildRefused as exc:
        # The rendered report, not just the fact of failure: a curator who ran
        # `build` needs to know what to fix without being sent to a second
        # command to find out. `build_site` writes nothing when it refuses, so
        # there is no partial `dist/` here for a deploy step to pick up.
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc
    except OSError as exc:
        # A destination that cannot be written is a bad argument, not a defect,
        # and the same exit code `schemas export` uses for the same reason. The
        # guard above catches the shapes visible before the build; this catches
        # the rest — a read-only parent, a full disk, a path that becomes
        # unwritable while the build runs. Reads are not in scope: `read_table`
        # and `_read_yaml` turn their own failures into validation issues, which
        # the gate above already acted on.
        typer.echo(f"error: could not write the site to {out}: {exc}")
        raise typer.Exit(code=2) from exc
    except ValueError as exc:
        # Every refusal inside `build/` is a `ValueError`: the emitter's path,
        # duplicate, case and seal guards, the variant shard-name and foreign-row
        # guards, the unresolvable featured PMID, the unsplit alias cell.
        #
        # These were documented as firing only on a bypassed gate. That is not
        # true yet and this is where it showed: `validate/` checks a variant
        # shard's `chrom` column but not its filename, so `mirrors/variants/
        # chr12.tsv` validates clean at 0 errors and then raises here. Measured.
        # Without this clause typer printed a rich traceback and exited 1, so a
        # curator got a stack trace under the code that means "your data is
        # wrong" — right about the data, useless about which file.
        #
        # Exit 1 is correct: it *is* a repository problem, the same class as a
        # refusal. What was wrong was the presentation. Unlike `BuildRefused`,
        # this fires after writing has begun, so `out` may hold a partial site —
        # said plainly, because a deploy step must not upload it.
        typer.echo(f"error: the build refused this repository: {exc}")
        typer.echo(
            f"note: {out} may hold a partial site and must not be published; "
            f"`chd-atlas validate` does not yet report this class of problem"
        )
        raise typer.Exit(code=1) from exc

    typer.echo(f"wrote {len(written)} file(s) to {out}")


@app.command()
def ingest(
    source_id: str | None = typer.Argument(
        None, help="The source's id in mirrors/sources.yaml (none with --all)."
    ),
    package: Path | None = typer.Option(
        None,
        "--package",
        help="The package directory: manifest.json and the tables it lists. "
        "Default: the pinned version's cache, local/packages/SOURCE_ID/VERSION.",
    ),
    every: bool = typer.Option(
        False,
        "--all",
        help="Re-ingest every pinned package from its cache, remove the rows of genes that "
        "left the registry, and validate once at the end.",
    ),
    repin: list[str] | None = typer.Option(
        None,
        "--repin",
        metavar="SOURCE=DIR",
        help="With --all: ingest SOURCE from the package directory DIR instead of its pinned "
        "cache, and pin DIR's version, in the same transaction. Repeatable. The route for a "
        "ClinGen or GenCC package, whatever it does to the gene registry; ingested alone "
        "it is refused when it brings a gene into the registry.",
    ),
    root: Path = typer.Option(Path("."), help="Repository root to ingest into."),
) -> None:
    """Merge a verified source package into mirrors/ and pin it in sources.yaml.

    With --all, also removes from every registry-filtered mirror the rows of genes that
    left the registry (the hgnc step), printing each file, gene and row count.
    """
    # Exit 2 for an argument that cannot work or a filesystem error, 1 for a
    # package or repository the ingest refuses (every file restored): the same
    # split `validate` and `build` make.
    one_of = "give a SOURCE_ID, or --all with no SOURCE_ID or --package"
    if not root.is_dir():
        raise _fail(f"--root {root} is not a directory", 2)
    if repin and not every:
        raise _fail("--repin needs --all (a single ingest takes its package with --package)", 2)
    if every:
        if source_id is not None or package is not None:
            raise _fail(one_of, 2)
        repins = _repins(repin or [])
        with _ingest_errors():
            every_result = ingest_all(root, repins)
        for line in _ingest_all_lines(root, every_result):
            typer.echo(line)
        return
    if source_id is None:
        raise _fail(one_of, 2)
    if package is not None and not package.is_dir():
        raise _fail(f"--package {package} is not a directory", 2)
    with _ingest_errors():
        if package is None:
            pin = package_pins(root).get(source_id)
            if pin is None:
                raise _fail(
                    f"source '{source_id}' has no package pin in mirrors/sources.yaml; "
                    "give the package directory with --package",
                    2,
                )
            package = pinned_package_dir(root, source_id, pin)
        result = ingest_package(root, source_id, package)
    for line in _ingest_lines(root, result):
        typer.echo(line)


def _repins(values: Sequence[str]) -> dict[str, Path]:
    """Each ``--repin SOURCE=DIR`` by source; exit 2 for one that cannot work."""
    found: dict[str, Path] = {}
    for value in values:
        source_id, equals, directory = value.partition("=")
        if not (source_id and equals and directory):
            raise _fail(f"--repin {value} is not SOURCE=DIR", 2)
        if source_id in found:
            raise _fail(f"--repin names {source_id} more than once", 2)
        path = Path(directory)
        if not path.is_dir():
            raise _fail(f"--repin {value}: {directory} is not a directory", 2)
        found[source_id] = path
    return found


@contextmanager
def _ingest_errors() -> Iterator[None]:
    """Exit 1 printing a refusal (every file restored), 2 on a filesystem error."""
    try:
        yield
    except IngestRefused as exc:
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc
    except OSError as exc:
        # A write that failed, its files restored, or a restore that failed,
        # whose message names every path it left changed (`apply._restore`).
        raise _fail(str(exc), 2) from exc


@app.command("report-changes")
def report_changes(
    base_root: Path = typer.Option(..., "--base-root", help="The base checkout."),
    base_site: Path = typer.Option(..., "--base-site", help="The site built from the base."),
    head_root: Path = typer.Option(..., "--head-root", help="The head checkout."),
    head_site: Path = typer.Option(..., "--head-site", help="The site built from the head."),
    out: Path | None = typer.Option(None, "--out", help="Write the report here, not to stdout."),
) -> None:
    """Report what a change does to the published site, as one Markdown comment."""
    for flag, path in (("--base-root", base_root), ("--head-root", head_root)):
        if not path.is_dir():
            raise _fail(f"{flag} {path} is not a directory", 2)
    for flag, path in (("--base-site", base_site), ("--head-site", head_site)):
        if not (path / "manifest.json").is_file():
            raise _fail(f"{flag} {path} has no manifest.json", 2)
    report = render_markdown(site_changes(base_root, base_site, head_root, head_site))
    if out is None:
        typer.echo(report, nl=False)
        return
    try:
        write_bytes_atomically(out, report.encode("utf-8"))
    except OSError as exc:
        # As `build` reports an unwritable --out.
        raise _fail(f"could not write the report to {out}: {exc}", 2) from exc


def _ingest_lines(root: Path, result: IngestResult) -> list[str]:
    """What an ingest wrote and, by name, everything it did not; then its validation."""
    lines = _summary_lines(root, result)
    if result.source_id in REGISTRY_INPUTS:
        lines.append(ALL_NOTE)
    return lines + _validation_lines(result.warnings)


def _ingest_all_lines(root: Path, every: IngestAllResult) -> list[str]:
    """Each source's summary, every pruned row by file and gene, then the one validation."""
    lines = [line for result in every.results for line in _summary_lines(root, result)]
    if not every.results:
        lines.append(
            "no source in mirrors/sources.yaml is pinned to a package; nothing was ingested"
        )
    if every.pruned:
        genes = sorted({item.gene for item in every.pruned})
        total = sum(item.rows for item in every.pruned)
        lines.append(
            f"removed {total} row(s) naming {len(genes)} gene(s) that left the registry: "
            + ", ".join(genes)
        )
        by_path: dict[str, list[str]] = {}
        for item in every.pruned:
            by_path.setdefault(item.path, []).append(f"{item.gene} ({item.rows})")
        lines += [f"removed from {path}: {', '.join(found)}" for path, found in by_path.items()]
    return lines + _validation_lines(every.warnings)


# How many unknown ids the summary names. An every-gene package can leave
# thousands unresolved (E-MTAB-6814: 16,868 of 58,735 Ensembl ids, measured by
# the phase 3 plan), which as one line buries everything else; the full list
# goes to a file outside the package cache instead.
UNRESOLVED_SHOWN: Final = 25

# Printed after a successful single ClinGen or GenCC ingest, which writes only
# its own mirror (`ingest.apply` module docstring).
ALL_NOTE: Final = (
    "note: a single ingest writes only this mirror; `chd-atlas ingest --all` (with "
    "`--repin <id>=<dir>` for a new version) re-filters the registry, the "
    "registry-filtered mirrors and mondo_xrefs"
)


def _write_unresolved(root: Path, result: IngestResult) -> str:
    """Write every unresolved id to local/ingest/; return its path under ``root``.

    Raises OSError; `_summary_lines` turns that into a line.

    Outside ``local/packages/`` so the cache's bytes (ING017) are untouched.
    """
    relative = f"local/ingest/{result.source_id}-{result.version}.unresolved.tsv"
    rows = ["id\tresolution\tcandidates"] + [
        f"{item.symbol}\t{item.resolution.value}\t{'|'.join(item.candidates)}"
        for item in sorted(result.unresolved, key=lambda item: item.symbol)
    ]
    (root / relative).parent.mkdir(parents=True, exist_ok=True)
    write_bytes_atomically(root / relative, ("\n".join(rows) + "\n").encode("utf-8"))
    return relative


def _summary_lines(root: Path, result: IngestResult) -> list[str]:
    """What an ingest wrote and, by name, everything it did not."""
    lines = [
        f"ingested {result.source_id} {result.version}: {result.rows_in} row(s) in the "
        f"package, {result.rows_written} written",
        f"{result.rows_outside_registry} row(s) name a gene outside the registry and were "
        f"not written",
        f"{result.rows_unresolved} row(s) name a gene id that did not resolve and were not written",
    ]
    if result.unresolved:
        path = f"local/ingest/{result.source_id}-{result.version}.unresolved.tsv"
        # The ingest has committed by now, so a file that cannot be written is
        # a line of its own: raising here would hide the summary of a change
        # already made (and, for --all, every source's).
        write_failure: str | None = None
        try:
            _write_unresolved(root, result)
        except OSError as exc:
            write_failure = f"could not write the complete list of unresolved ids to {path}: {exc}"
        ordered = sorted(result.unresolved, key=lambda item: item.symbol)
        # An ambiguous id is a collision a curator must read, so all are named;
        # only the plain unknowns are capped.
        ambiguous = [item for item in ordered if item.candidates]
        unknown = [item for item in ordered if not item.candidates]
        shown = ambiguous + unknown[:UNRESOLVED_SHOWN]
        names = "; ".join(
            f"{item.symbol} ({item.resolution.value}"
            + (f": {', '.join(item.candidates)}" if item.candidates else "")
            + ")"
            for item in shown
        )
        more = len(unknown) - UNRESOLVED_SHOWN
        if more > 0:
            where = f" in {path}" if write_failure is None else ""
            names += f"; … and {more} more (all {len(ordered)}{where})"
        lines.append(
            f"{len(result.unresolved)} gene id(s) did not resolve and their rows were not "
            f"written: {names}"
        )
        if write_failure is not None:
            lines.append(write_failure)
    # Only for a package a table filter applied to; an unknown term is its own
    # line, because it is not a row the filter judged (`ingest.filters`).
    if result.filter_reasons:
        lines.append(
            f"{result.rows_outside_filter} row(s) are outside the filter "
            f"({'; '.join(result.filter_reasons)})"
        )
    if result.unknown_terms:
        lines.append(
            f"{result.rows_unknown_term} row(s) name a term the pinned release does not know: "
            + ", ".join(result.unknown_terms)
        )
    absent = result.registry_genes_absent
    if absent:
        genes = ", ".join(absent)
        # Facts only: why a gene has no row (the source has none, its id did
        # not resolve, the registry changed after a slice was cut) is not
        # something the ingest measures.
        if result.coverage == "every_gene":
            lines.append(
                f"{len(absent)} registry gene(s) have no resolved row in the package: {genes}"
            )
        else:
            lines.append(
                f"{len(absent)} registry gene(s) have no resolved row in this registry slice "
                f"(cut at {result.coverage.atlas_commit} over {result.coverage.genes} genes; "
                f"the registry now lists {result.registry_genes}): {genes}"
            )
    for entry in result.dropped:
        examples = f" (e.g. {', '.join(entry.examples)})" if entry.examples else ""
        lines.append(f"dropped before packaging: {entry.reason}: {entry.count}{examples}")
    return lines


def _validation_lines(warnings: Sequence[ValidationIssue]) -> list[str]:
    # An error refuses the ingest, so the count is 0 whenever this prints.
    return [
        f"validate: 0 error(s), {len(warnings)} warning(s)",
        *(issue.format() for issue in warnings),
    ]


@schemas_app.command("export")
def schemas_export(
    target: Path = typer.Option(Path("schemas"), help="Directory to write schemas into."),
) -> None:
    """Regenerate committed JSON Schema files from the Pydantic models."""
    # A target that is an existing file, or a directory that cannot be written,
    # is a bad argument rather than a defect: report the path and exit 2 rather
    # than showing the user a traceback. FileExistsError and PermissionError
    # are both OSError.
    try:
        written = export_schemas(target)
        tables = export_table_schemas(target / TABLE_SCHEMA_DIR)
    except OSError as exc:
        typer.echo(f"error: could not write schemas to {target}: {exc}")
        raise typer.Exit(code=2) from exc

    typer.echo(f"wrote {len(written)} schema(s) to {target}")
    typer.echo(f"wrote {len(tables)} table schema(s) to {target / TABLE_SCHEMA_DIR}")


def _fail(message: str, code: int) -> typer.Exit:
    typer.echo(f"error: {message}")
    return typer.Exit(code=code)


def _require_root(root: Path) -> None:
    if not root.is_dir():
        raise _fail(f"--root {root} is not a directory", 2)


def _resolve_gene(root: Path, symbol: str) -> tuple[str, str, Sequence[str]]:
    """`(HGNC id, approved symbol, previous symbols)` for a typed symbol, or exit 2."""
    try:
        return resolve_symbol(root, symbol)
    except TargetError as exc:
        raise _fail(str(exc), 2) from None


@literature_app.command("target")
def literature_target_command(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Print the atlas's answer about a gene as one JSON line."""
    _require_root(root)
    try:
        target = literature_target(root, symbol)
    except TargetError as exc:
        raise _fail(str(exc), 2) from None
    typer.echo(json.dumps(target.as_json(), sort_keys=True, ensure_ascii=False))


@literature_app.command("frame")
def literature_frame_command(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Print a gene's pinned MANE Select transcript and protein as one JSON line."""
    _require_root(root)
    try:
        frame = literature_frame(root, symbol)
    except TargetError as exc:
        raise _fail(str(exc), 2) from None
    typer.echo(json.dumps(frame.as_json(), sort_keys=True, ensure_ascii=False))


@literature_app.command("digests")
def literature_digests_command(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Print each finding's digests and the file's sha256 as one JSON line. Read-only.

    `finding_digest` is what an approval binds to; `shown_digest` also covers
    the review marks and the paper's study type and own-lab mark, so a review
    page can tell that a card changed after it was shown. `approved` is true for
    any finding carrying an approval, void or not; `approval_void` is true only
    when the approval no longer matches the finding (LIT015). Exit 2 for an
    unknown symbol or a gene with no literature file, 1 for an unparsable
    literature file or a `curation/publications.yaml` that does not load.
    """
    _require_root(root)
    approved_symbol, path = _require_literature_file(root, symbol)
    raw = path.read_bytes()
    try:
        record = parse_literature(raw)
    except LiteratureFileError as exc:
        raise _fail(f"{path}: {exc}", 1) from None
    papers = _publications_by_id(root)
    findings = [
        {
            "id": f.id,
            "finding_digest": finding_digest(f),
            "shown_digest": shown_digest(f, papers.get(f.publication)),
            "approved": f.approved_by is not None,
            "approval_void": _approval_void(f),
            "verified": f.verified_by is not None,
        }
        for f in sorted(record.findings, key=lambda f: f.id)
    ]
    out = {
        "format": 1,
        "gene": str(record.gene),
        "symbol": approved_symbol,
        "literature_sha256": "sha256:" + hashlib.sha256(raw).hexdigest(),
        "findings": findings,
    }
    typer.echo(json.dumps(out, sort_keys=True, ensure_ascii=False))


def _publications_by_id(root: Path) -> dict[str, Publication]:
    """`curation/publications.yaml` by paper id; exit 1 naming the file if it did not load.

    `load_curation` reports an unreadable file as an issue and returns no
    publications. Continuing would hash every finding with `publication: null`,
    so an edit to a paper's study type would go unnoticed.
    """
    corpus, issues = load_curation(root)
    broken = [i for i in issues if "publications.yaml" in i.location]
    if broken:
        raise _fail(f"{broken[0].location}: {broken[0].message}", 1)
    return {p.id: p for p in corpus.publications}


def _require_literature_file(root: Path, symbol: str) -> tuple[str, Path]:
    """`(approved symbol, path)` of an existing literature file, or exit 2 naming the path."""
    _, approved, _ = _resolve_gene(root, symbol)
    path = literature_path(root, approved)
    if not path.is_file():
        raise _fail(f"{path} does not exist", 2)
    return approved, path


def _ids(values: Sequence[str]) -> str:
    return ", ".join(values) or "none"


def _is_agent(name: str) -> bool:
    return name.strip().lower().startswith("agent:")


def _locator_text(locator: Locator) -> str:
    if isinstance(locator, InTextLocator):
        return locator.element
    parts = [locator.file, locator.sheet, f"row {locator.row}" if locator.row else None]
    return " / ".join(p for p in parts if p)


def _counts_text(finding: Finding) -> str:
    groups: list[str] = []
    if finding.screened is not None:
        unit = finding.count_unit.value if finding.count_unit else "(unit not stated)"
        groups.append(f"screened {finding.screened} {unit}")
    if finding.carriers is not None:
        c = finding.carriers
        groups.append(
            f"carriers {c.total} (de novo {c.de_novo}, inherited {c.inherited},"
            f" origin unknown {c.origin_unknown})"
        )
    tail = [
        text
        for text in (
            f"cases {finding.cases}" if finding.cases is not None else "",
            f"controls {finding.controls}" if finding.controls is not None else "",
            f"p {finding.p_value!r}" if finding.p_value is not None else "",
            f"OR {finding.odds_ratio!r}" if finding.odds_ratio is not None else "",
        )
        if text
    ]
    if tail:
        groups.append(", ".join(tail))
    return "; ".join(groups) or "no counts"


def _today() -> date:
    """Today's date; a test replaces it, so the 'not in the future' rule is measurable."""
    return date.today()


def _variant_line(variant: Variant) -> str:
    """One variant as the person approving reads it, every digest-covered field but one.

    Shown: the paper's own text, placement, the checked names, zygosity and
    links, carriers, phenotypes, `kind`, the non-null `reported` names (labelled
    transcript, c, p, rsid) and the locator (its kind, then its element or file).
    Omitted on purpose: `normalisation`, the id of the normaliser run, which is
    bookkeeping a person approving the finding has nothing to read in.
    """
    names = f"{variant.hgvs_c or '-'} {variant.hgvs_p or '-'} {variant.grch38 or '-'}"
    zygosity = variant.zygosity.value
    if variant.with_:
        phase = f" ({variant.phase.value})" if variant.phase is not None else ""
        zygosity += f" with {', '.join(str(i) for i in variant.with_)}{phase}"
    c = variant.carriers
    carriers = f"{c.de_novo}/{c.inherited}/{c.origin_unknown}" if c is not None else "not stated"
    r = variant.reported
    reported = (
        ", ".join(
            f"{label} {value}"
            for label, value in (
                ("transcript", r.transcript),
                ("c", r.c),
                ("p", r.p),
                ("rsid", r.rsid),
            )
            if value is not None
        )
        or "none"
    )
    line = (
        f"  variant {variant.id}: {variant.as_reported} | {variant.placement.value} | {names}"
        f" | {zygosity} | carriers {carriers} | {', '.join(variant.phenotypes) or 'no phenotypes'}"
        f" | kind {variant.kind.value} | reported {reported}"
        f" | locator {variant.locator.kind} {_locator_text(variant.locator)}"
    )
    return line


def _mark_text(
    status: ReviewStatus | None, note: str | None, settled_by: SettledBy | None
) -> str | None:
    """`review: pending_review (settled by the authors): <note>`; a `not_accessible` mark
    names no settler (the full text settles it by definition); None when unmarked."""
    if status is None:
        return None
    if settled_by is None:
        return f"review: {status.value}: {note}"
    return f"review: {status.value} (settled by {SETTLED_LABELS[settled_by.value]}): {note}"


def _variant_lines(variant: Variant) -> list[str]:
    """The variant's line, then (indented) its own mark and one line per unclear field."""
    lines = [_variant_line(variant)]
    mark = _mark_text(variant.review_status, variant.review_note, variant.review_settled_by)
    if mark is not None:
        lines.append(f"    {mark}")
    for field in ("zygosity", "phase"):
        entry = variant.unclear.get(field)
        if entry is not None:
            lines.append(f"    unclear {field} ({entry.settled_by.value}): {entry.reason}")
    return lines


def _listing(finding: Finding, paper: Publication | None) -> list[str]:
    """Every field the approval digest covers, except the bookkeeping ones.

    Omitted on purpose: `id` and the link (the heading line), `approved_*`
    (a listed finding is unapproved, or its approval is void and `APPROVAL
    VOID` stands in for it), `verified_by` (the heading line) and each
    variant's `normalisation` id (see `_variant_line`). Floats go
    through `repr`, so 1.23456789e-9 is not shown as 1.23457e-09 (final review
    I2). The paper's `study_type` and `own_lab` come from publications.yaml,
    where an agent's reading of the design may be the only copy (R32).
    """
    link = f"https://pubmed.ncbi.nlm.nih.gov/{finding.publication.removeprefix('PMID:')}/"
    if paper is None:
        paper_line = "paper: not in curation/publications.yaml"
    else:
        own_lab = "yes" if paper.own_lab else "no"
        paper_line = f"paper: study_type {paper.study_type.value}; own_lab {own_lab}"
    shared = finding.shares_probands_with
    void = "  APPROVAL VOID (LIT015)" if _approval_void(finding) else ""
    lines = [
        f"{finding.id}  {link}  {finding.design.value}  {finding.verified_by or 'UNVERIFIED'}"
        f"{void}",
        f"  {paper_line}",
        f"  methods: {', '.join(m.value for m in finding.methods)}",
        f"  comparator: {finding.comparator.value if finding.comparator else 'none'}",
        f"  burden study: {finding.burden_study or 'none'}",
        f"  variant classes: {_ids([v.value for v in finding.variant_classes])}",
        *(line for v in finding.variants for line in _variant_lines(v)),
        f"  {_counts_text(finding)}",
        f"  phenotypes: {', '.join(finding.phenotypes)}",
        f"  lesion groups: {', '.join(g.value for g in finding.lesion_groups)}",
        f"  syndromic: {finding.syndromic.value}",
        f"  cohorts: {_ids(finding.cohorts)}",
        f"  locators: {'; '.join(_locator_text(loc) for loc in finding.locators)}",
        f"  shares probands with: {shared if isinstance(shared, str) else _ids(shared)}",
    ]
    mark = _mark_text(finding.review_status, finding.review_note, finding.review_settled_by)
    if mark is not None:
        lines.append(f"  {mark}")
    lines.append(f"  extracted by: {finding.extracted_by}")
    lines.append(f"  {finding.summary}")
    return lines


def _first_search_run_on(record: GeneLiteratureFile, finding: Finding) -> date | None:
    """The run date of the search that first found the finding's paper (the LIT008 rule)."""
    first = next((e.search for e in record.triage if e.pmid == finding.publication), None)
    return next((s.run_on for s in record.searches if s.id == first), None)


def _own_assertions(curation: Any, record: GeneLiteratureFile) -> list[Any]:
    """The loaded curation's lesion assertions for the file's own gene.

    The one filter behind the synthesis listing and `--approve-synthesis`: the
    synthesis digest covers the gene's own assertions, so both must pass the
    same list or an approval is stale the moment it is written.
    """
    return [a for a in curation.assertions if str(a.gene) == str(record.gene)]


def _approval_void(finding: Finding) -> bool:
    """The LIT015 rule: an approved finding whose content no longer matches its approval."""
    return finding.approved_by is not None and finding.approved_digest != finding_digest(finding)


def _approval_state(finding: Finding) -> str:
    """` [unapproved]`, ` [approval void: LIT015]`, or nothing for a sound approval."""
    if finding.approved_by is None:
        return " [unapproved]"
    if _approval_void(finding):
        return " [approval void: LIT015]"
    return ""


def _cited_state(finding: Finding) -> str:
    """The mark after a cited finding: ` [unapproved]`, ` [approval void: LIT015]`, or
    ` [<review mark>]`, in that precedence, else nothing."""
    state = _approval_state(finding)
    if state:
        return state
    if finding.review_status is not None:
        return f" [{finding.review_status.value}]"
    return ""


# Hex characters of the synthesis digest the listing prints; `--digest` must give
# at least this many, each matching.
_LISTED_HEX: Final = 12


def _listed_digest(digest: str) -> str:
    """What the listing prints for `--digest`: the digest's first `_LISTED_HEX` hex characters."""
    return digest.removeprefix("sha256:")[:_LISTED_HEX]


def _synthesis_listing(
    symbol: str, record: GeneLiteratureFile, assertions: Sequence[Any]
) -> list[str]:
    """The synthesis's approval state, its sentences and what each cites.

    An unapproved or stale synthesis is followed by the command that approves
    it, carrying the digest of what is listed (`--digest`). A cited assertion
    is shown with its lesion groups, phenotypes and evidence, so the person
    approving reads what the sentence rests on.
    """
    synthesis = record.synthesis
    if synthesis is None:
        return []
    digest = synthesis_digest(synthesis, record.findings, assertions)
    approve = (
        f"  approve with: chd-atlas literature review {symbol} --approve-synthesis "
        f"--by <name> --digest {_listed_digest(digest)}"
    )
    if synthesis.approved_by is None:
        lines = ["synthesis: unapproved", approve]
    else:
        who = f"approved by {synthesis.approved_by} on {synthesis.approved_on}"
        if synthesis.approved_digest == digest:
            lines = [f"synthesis: {who}"]
        else:
            lines = [
                f"synthesis: STALE ({who}; the synthesis or the records it covers changed since)",
                approve,
            ]
    findings: dict[str, Finding] = {f.id: f for f in record.findings}
    own: dict[str, Any] = {a.id: a for a in assertions}
    for number, sentence in enumerate(synthesis.sentences, 1):
        lines.append(f"  {number}. {sentence.text}")
        for cite in sentence.cites:
            if cite in findings:
                cited = findings[cite]
                lines.append(f"     {cite}: {cited.summary}{_cited_state(cited)}")
            elif cite in own:
                assertion = own[cite]
                lines.append(f"     {cite}: lesion assertion for {assertion.gene}")
                groups = ", ".join(g.value for g in assertion.lesion_groups)
                lines.append(f"       lesion groups: {groups}")
                lines.append(f"       phenotypes: {', '.join(assertion.phenotypes)}")
                for item in assertion.evidence:
                    lines.append(f"       evidence: {item.publication} — {item.summary}")
            else:
                lines.append(f"     {cite}: (not found)")
    return lines


def _write_unless_changed(path: Path, read: bytes, content: bytes) -> None:
    """Write `content` unless the file's bytes differ from `read`, those this command read.

    Re-read immediately before the write: another writer (a person approving, an
    agent redrafting) whose bytes landed since would otherwise be overwritten
    with both commands reporting success (K-1). The window left is the one
    between this read and the rename.
    """
    if path.read_bytes() != read:
        raise _fail(f"{path} changed since it was read; nothing was written", 2)
    write_bytes_atomically(path, content)


def _approval_updates(
    record: GeneLiteratureFile,
    ids: Sequence[str],
    by: str,
    when: date,
    allow_unverified: bool,
    path: Path,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    """The approval fields for `ids` and one message per finding; exit 2 if any is refused.

    Refused: an unknown id, a finding already soundly approved (approving again
    would erase who approved it), a date before the search that found the paper
    (LIT008), and an unverified finding unless `allow_unverified`. Nothing is
    written here; the caller makes the one write.
    """
    by_id: dict[str, Finding] = {f.id: f for f in record.findings}
    ordered = sorted(set(ids))
    for finding_id in ordered:
        target = by_id.get(finding_id)
        if target is None:
            raise _fail(f"no such finding {finding_id} in {path}", 2)
        if target.approved_by is not None and not _approval_void(target):
            raise _fail(
                f"{finding_id} is already approved by {target.approved_by}; "
                f"approving again would erase who approved it",
                2,
            )
        found_on = _first_search_run_on(record, target)
        if found_on is not None and when < found_on:
            raise _fail(
                f"--on {when.isoformat()} is before search {found_on.isoformat()} "
                f"found {target.publication}; an approval cannot predate it (LIT008)",
                2,
            )
        if target.verified_by is None and not allow_unverified:
            raise _fail(
                f"{finding_id} is unverified (no verified_by); "
                f"run `chd-curator diff` or pass --allow-unverified",
                2,
            )
    updates = {
        i: {
            "approved_by": by.strip(),
            "approved_on": when,
            "approved_digest": finding_digest(by_id[i]),
        }
        for i in ordered
    }
    messages = [
        f"re-approved {i} ({by.strip()}, {when.isoformat()}; "
        f"was {by_id[i].approved_by} on {by_id[i].approved_on})"
        if by_id[i].approved_by is not None
        else f"approved {i} ({by.strip()}, {when.isoformat()})"
        for i in ordered
    ]
    return updates, messages


@literature_app.command("questions")
def literature_questions(
    symbol: str | None = typer.Argument(
        None, help="Gene symbol (approved, alias or previous); omit for every literature file."
    ),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """List every open question in a gene's findings, approved or not. Read-only.

    Per gene: `<SYMBOL>: N open questions`, then a heading per settler in the
    order the atlas ranks them, `settled by <label> (K)`; N and every K count
    questions, not records (a `with` group's shared unclear phase is one question
    over several variants). A line names the finding (and `v<id>`), its PMID and the
    question; a finding that is unapproved or whose approval is void is marked
    exactly as `literature review` marks it. Exit 2 for an unknown symbol or a gene
    with no literature file.
    """
    _require_root(root)
    if symbol is not None:
        approved, path = _require_literature_file(root, symbol)
        files = [(approved, path)]
    else:
        directory = literature_path(root, "X").parent
        files = (
            [(p.stem, p) for p in sorted(directory.glob("*.yaml"))] if directory.is_dir() else []
        )
    for name, path in files:
        try:
            record = parse_literature(path.read_bytes())
        except LiteratureFileError as exc:
            raise _fail(f"{path}: {exc}", 1) from None
        findings = sorted(record.findings, key=lambda f: f.id)
        by_id: dict[str, Finding] = {f.id: f for f in findings}
        questions = open_questions([f.model_dump(mode="json", by_alias=True) for f in findings])
        count = len(questions)
        typer.echo(f"{name}: {count} open question{'' if count == 1 else 's'}")
        for settler in SETTLED_ORDER:
            group = [q for q in questions if q["settled_by"] == settler]
            if not group:
                continue
            typer.echo(f"settled by {SETTLED_LABELS[settler]} ({len(group)})")
            for question in group:
                for rec in question["records"]:
                    where = (
                        rec["finding"]
                        if rec["variant"] is None
                        else (f"{rec['finding']} v{rec['variant']}")
                    )
                    typer.echo(
                        f"  {where}  {rec['publication']}  {QUESTION_LABELS[question['kind']]}"
                        f"{_approval_state(by_id[rec['finding']])}"
                    )
                if question["kind"] != "not_accessible":
                    typer.echo(f"    {question['reason']}")


@literature_app.command("review")
def literature_review(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    approve: list[str] = typer.Option(
        [],
        "--approve",
        help="A finding id to approve, or to approve again when its approval is void "
        "(LIT015); repeat the option for several.",
    ),
    by: str = typer.Option("", "--by", help="The approving person's name."),
    on: str = typer.Option("", "--on", help="Approval date, YYYY-MM-DD (default today)."),
    allow_unverified: bool = typer.Option(
        False, "--allow-unverified", help="Approve a finding no verifier has compared."
    ),
    against: Path | None = typer.Option(
        None,
        "--against",
        help="Listing only: the blind second pass; prints each finding's field differences.",
    ),
    pending_review: str = typer.Option(
        "", "--pending-review", help="A finding id to mark pending_review."
    ),
    approve_synthesis: bool = typer.Option(
        False,
        "--approve-synthesis",
        help="Approve the gene's synthesis (a stale approval is re-approved).",
    ),
    digest: str = typer.Option(
        "",
        "--digest",
        help="With --approve-synthesis: the digest the listing printed (12 or more hex digits).",
    ),
    note: str = typer.Option("", "--note", help="The measured disagreement, for --pending-review."),
    settled_by: SettledBy | None = typer.Option(
        None,
        "--settled-by",
        help="What would settle the question, for --pending-review: "
        "full_text, authors, expert or atlas_model.",
    ),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """List the findings awaiting approval (and the synthesis, if any), approve some, or mark one.

    `--approve-synthesis` approves the gene's synthesis, or re-approves a stale
    one, only with `--digest` matching the synthesis now on disk: the listing
    prints the command with it, so the person approves the text they read. It is
    refused while any approved finding of the file fails LIT015. A review mark is
    published with the finding; it does not withdraw an approval, and a marked
    finding can still be approved. A finding carries one mark: a second is
    refused rather than overwriting the first note. A finding whose approval is
    void (LIT015) is listed, marked `APPROVAL VOID`, and may be approved again;
    the output names the approval it replaced. Every write is refused if the
    file changed since this command read it.
    """
    _require_root(root)
    if sum(bool(a) for a in (approve, pending_review, approve_synthesis)) > 1:
        raise _fail(
            "--approve, --pending-review and --approve-synthesis are separate actions; "
            "use only one of them",
            2,
        )
    if approve_synthesis and allow_unverified:
        raise _fail("--allow-unverified does not apply to --approve-synthesis", 2)
    if not approve and not approve_synthesis and (by or on or allow_unverified):
        raise _fail(
            "--by and --on only apply with --approve or --approve-synthesis, "
            "and --allow-unverified only with --approve",
            2,
        )
    if not pending_review and note:
        raise _fail("--note only applies with --pending-review", 2)
    if not pending_review and settled_by is not None:
        raise _fail("--settled-by only applies with --pending-review", 2)
    if not approve_synthesis and digest:
        raise _fail("--digest only applies with --approve-synthesis", 2)
    if against is not None and (approve or pending_review or approve_synthesis):
        raise _fail("--against only applies to the listing, not to an action", 2)
    if approve or approve_synthesis:
        flag = "--approve" if approve else "--approve-synthesis"
        if not by.strip():
            raise _fail(f"--by (the approving person's name) is required with {flag}", 2)
        if _is_agent(by):
            raise _fail(f"--by {by!r} names an agent; a person must approve", 2)
        try:
            when = date.fromisoformat(on) if on else _today()
        except ValueError:
            raise _fail(f"--on {on!r} is not a YYYY-MM-DD date", 2) from None
        if when > _today():
            raise _fail(f"--on {when.isoformat()} is in the future", 2)
    if approve_synthesis and not digest.strip():
        raise _fail(
            "--digest is required with --approve-synthesis: run "
            f"`chd-atlas literature review {symbol}` and copy the command it prints",
            2,
        )
    if pending_review and not note.strip():
        raise _fail("--note (the measured disagreement) is required with --pending-review", 2)
    if pending_review and settled_by is None:
        raise _fail(
            "--settled-by (what would settle the question: "
            + ", ".join(v.value for v in SettledBy)
            + ") is required with --pending-review",
            2,
        )

    approved_symbol, path = _require_literature_file(root, symbol)
    raw = path.read_bytes()
    try:
        record = parse_literature(raw)
    except LiteratureFileError as exc:
        raise _fail(f"{path}: {exc}", 1) from None
    by_id: dict[str, Finding] = {f.id: f for f in record.findings}

    if not approve and not pending_review and not approve_synthesis:
        waiting = sorted(
            (f for f in record.findings if f.approved_by is None or _approval_void(f)),
            key=lambda f: f.id,
        )
        if not waiting:
            typer.echo(f"{approved_symbol}: no findings awaiting approval")
        curation = load_curation(root)[0]
        papers = {p.id: p for p in curation.publications}
        differences: dict[str, list[tuple[str, Any, Any]]] = {}
        unpaired: set[str] = set()
        if against is not None:
            if not against.is_file():
                raise _fail(f"--against {against} is not a file", 2)
            try:
                second = parse_findings_file(against.read_bytes())
            except LiteratureFileError as exc:
                raise _fail(str(exc), 1) from None
            diff = diff_findings(record.findings, second)
            differences = {k: list(v) for k, v in diff.fields.items()}
            unpaired = set(diff.only_first)
        for f in waiting:
            for line in _listing(f, papers.get(f.publication)):
                typer.echo(line)
            if against is not None:
                if f.id in unpaired:
                    typer.echo("  second pass: no partner finding")
                for name, mine, theirs in differences.get(f.id, []):
                    typer.echo(
                        f"  differs: {name}: first {json.dumps(mine, sort_keys=True)}"
                        f" | second {json.dumps(theirs, sort_keys=True)}"
                    )
            typer.echo("")
        if record.synthesis is not None:
            own = _own_assertions(curation, record)
            for line in _synthesis_listing(approved_symbol, record, own):
                typer.echo(line)
        return

    if approve_synthesis:
        synthesis = record.synthesis
        if synthesis is None:
            raise _fail(f"{path} has no synthesis to approve", 2)
        own = _own_assertions(load_curation(root)[0], record)
        current = synthesis_digest(synthesis, record.findings, own)
        if synthesis.approved_by is not None and synthesis.approved_digest == current:
            raise _fail(
                f"the synthesis is already approved by {synthesis.approved_by}; "
                f"approving again would erase who approved it",
                2,
            )
        given = digest.strip()
        if len(given) < _LISTED_HEX or not current.removeprefix("sha256:").startswith(given):
            raise _fail(
                f"the synthesis is not the one listed (its digest is {_listed_digest(current)}); "
                "list it again before approving",
                2,
            )
        void = sorted(f.id for f in record.findings if _approval_void(f))
        if void:
            raise _fail(
                f"cannot approve the synthesis: {', '.join(void)} fail(s) LIT015 "
                "(approved content changed since approval); re-approve or correct them first",
                2,
            )
        errors = sorted(
            {
                i.code
                for i in validate_gene_synthesis(
                    Corpus(root=root, assertions=tuple(own), literature=((path, record),))
                )
                if i.code in {"SYN001", "SYN002", "SYN003"}
            }
        )
        if errors:
            raise _fail(
                f"the synthesis fails validation ({', '.join(errors)}); "
                f"fix it before approving: `chd-atlas validate`",
                2,
            )
        approved_dates = [f.approved_on for f in record.findings if f.approved_on is not None]
        if approved_dates and when < max(approved_dates):
            latest = max(approved_dates)
            raise _fail(
                f"--on {when.isoformat()} is before {latest.isoformat()}, when a finding "
                f"this synthesis covers was approved; an approval cannot predate what it covers",
                2,
            )
        try:
            content = with_synthesis_approval(
                raw,
                {"approved_by": by.strip(), "approved_on": when, "approved_digest": current},
            )
        except LiteratureFileError as exc:
            raise _fail(str(exc), 1) from None
        _write_unless_changed(path, raw, content)
        typer.echo(f"approved synthesis ({by.strip()}, {when.isoformat()})")
        return

    if pending_review:
        target = by_id.get(pending_review)
        if target is None:
            raise _fail(f"no such finding {pending_review} in {path}", 2)
        if target.review_status is not None:
            raise _fail(
                f"{pending_review} already carries a review mark ({target.review_status.value}); "
                f"changing it is a hand edit",
                2,
            )
        updates: dict[str, dict[str, Any]] = {
            pending_review: {
                "review_status": ReviewStatus.PENDING_REVIEW.value,
                "review_note": note.strip(),
                "review_settled_by": settled_by.value if settled_by else None,
            }
        }
        messages = [f"marked {pending_review} pending_review"]
    else:
        updates, messages = _approval_updates(
            record, approve, by, when, allow_unverified=allow_unverified, path=path
        )
    try:
        content = with_findings_updated(raw, updates)
    except LiteratureFileError as exc:
        raise _fail(str(exc), 1) from None
    _write_unless_changed(path, raw, content)
    for message in messages:
        typer.echo(message)


@literature_app.command("apply-review")
def literature_apply_review(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    decisions: Path = typer.Option(..., "--decisions", help="The review page's decisions file."),
    by: str = typer.Option("", "--by", help="The approving person's name; the reviewer's."),
    on: str = typer.Option("", "--on", help="Approval date, YYYY-MM-DD (default today)."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Approve the findings a review page's owner approved, refusing any changed since it was built.

    Reads a decisions file (`chd-atlas literature digests` prints the digests a
    page is built from). Every refusal below exits 2 and writes nothing: an
    unreadable or invalid file, a file for another gene, a session never
    finished (`session.finished_at` empty), a reviewer who is a second opinion
    (never applied as approvals), `--by` differing from the reviewer's name or
    naming an agent, a finding id absent from the literature file (on any
    decision), `--on` in the future or before the date the page was built (local
    time: an approval cannot predate the page that showed it), and, found only
    for a finding about to be approved, after the stale and already-approved
    sort, `--on` before the search that found its paper (LIT008). A
    `curation/publications.yaml` that does not load exits 1 before any approval is written.
    Then per `approve`: already approved and sound, printed and skipped; its current
    `shown_digest` differs from the file's, printed as `stale` and skipped (exit 1 once the
    others are written); else approved, as `review --approve --allow-unverified` would, so a
    finding no verifier has compared is approved. All approvals go in one write. Written
    nowhere, only printed: the `next round:` lines for a `fix`, a `later`, and a triage
    review that is not `approve` (with its flips), `note on <id>` lines for a note left on
    an `approve` decision (whatever became of it), and `papers_seen`.
    """
    _require_root(root)
    try:
        review = ReviewDecisions.model_validate_json(decisions.read_bytes())
    except (OSError, ValueError) as exc:
        raise _fail(f"{decisions}: cannot be applied: {exc}", 2) from None
    approved_symbol, path = _require_literature_file(root, symbol)
    raw = path.read_bytes()
    try:
        record = parse_literature(raw)
    except LiteratureFileError as exc:
        raise _fail(f"{path}: {exc}", 1) from None
    if str(review.gene) != str(record.gene):
        raise _fail(
            f"{decisions} is for {review.gene} ({review.symbol}), not {approved_symbol} "
            f"({record.gene}); nothing was applied",
            2,
        )
    if review.session.finished_at is None:
        raise _fail(
            "the reviewer has not pressed I'm done (session.finished_at is empty); "
            "nothing was applied",
            2,
        )
    if review.reviewer.role is not ReviewerRole.OWNER:
        raise _fail("a second opinion is never applied as approvals; nothing was applied", 2)
    if not by.strip() or by.strip() != review.reviewer.name:
        raise _fail(
            f"--by {by!r} is not the reviewer's name ({review.reviewer.name!r}); "
            "nothing was applied",
            2,
        )
    if _is_agent(by):
        raise _fail(f"--by {by!r} names an agent; a person must approve", 2)
    by_id: dict[str, Finding] = {f.id: f for f in record.findings}
    unknown = sorted(str(d.id) for d in review.findings if str(d.id) not in by_id)
    if unknown:
        raise _fail(f"no such finding {', '.join(unknown)} in {path}; nothing was applied", 2)
    try:
        when = date.fromisoformat(on) if on else _today()
    except ValueError:
        raise _fail(f"--on {on!r} is not a YYYY-MM-DD date", 2) from None
    if when > _today():
        raise _fail(f"--on {when.isoformat()} is in the future", 2)
    built = review.page.built_on
    built_on = (built.astimezone() if built.tzinfo is not None else built).date()
    if when < built_on:
        raise _fail(
            f"--on {when.isoformat()} is before the page was built on {built_on.isoformat()}; "
            "an approval cannot predate the page that showed it",
            2,
        )

    papers = _publications_by_id(root)
    to_approve: list[str] = []
    lines: list[str] = []
    stale = 0
    for item in review.findings:
        if item.decision is not Decision.APPROVE:
            continue
        finding = by_id[str(item.id)]
        if finding.approved_by is not None and not _approval_void(finding):
            lines.append(f"already approved {finding.id} by {finding.approved_by}")
        elif shown_digest(finding, papers.get(finding.publication)) != item.shown_digest:
            lines.append(f"stale {finding.id}: changed since the page was built; review it again")
            stale += 1
        else:
            to_approve.append(finding.id)
    updates, messages = _approval_updates(
        record, to_approve, review.reviewer.name, when, allow_unverified=True, path=path
    )
    if updates:
        try:
            content = with_findings_updated(raw, updates)
        except LiteratureFileError as exc:
            raise _fail(str(exc), 1) from None
        _write_unless_changed(path, raw, content)
    lines.extend(messages)
    for item in review.findings:
        if item.decision is not Decision.APPROVE:
            lines.append(f"next round: {item.decision.value} {item.id}: {item.note.strip()}")
        elif item.note.strip():
            lines.append(f"note on {item.id}: {item.note.strip()}")
    if review.triage is not None:
        if review.triage.decision is not Decision.APPROVE:
            lines.append(
                f"next round: triage {review.triage.decision.value}: {review.triage.note.strip()}"
            )
        for flip in review.triage.flips:
            lines.append(f"next round: triage {flip.pmid} -> {flip.to}: {flip.note}")
    if review.papers_seen:
        lines.append("papers seen: " + ", ".join(str(p) for p in review.papers_seen))
    for line in lines:
        typer.echo(line)
    if stale:
        raise typer.Exit(1)
