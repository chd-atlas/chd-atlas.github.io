from __future__ import annotations

import json
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import date
from pathlib import Path
from typing import Any

import typer

from chd_atlas.build.runner import BuildRefused, build_site
from chd_atlas.changes import render_markdown, site_changes
from chd_atlas.corpus import load_curation
from chd_atlas.fs import write_bytes_atomically
from chd_atlas.gene_literature_diff import diff_findings
from chd_atlas.gene_literature_io import (
    LiteratureFileError,
    literature_path,
    parse_findings_file,
    parse_literature,
    with_findings_updated,
)
from chd_atlas.ingest.apply import (
    IngestAllResult,
    IngestResult,
    ingest_all,
    ingest_package,
    package_pins,
    pinned_package_dir,
)
from chd_atlas.ingest.package import IngestRefused
from chd_atlas.issues import ValidationIssue
from chd_atlas.literature_query import TargetError, literature_target, resolve_symbol
from chd_atlas.models.assertion import InTextLocator, Locator
from chd_atlas.models.gene_literature import Finding, GeneLiteratureFile, finding_digest
from chd_atlas.models.literature import Publication
from chd_atlas.models.review import ReviewStatus
from chd_atlas.schema_export import TABLE_SCHEMA_DIR, export_schemas, export_table_schemas
from chd_atlas.validate.runner import validate_repository

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
    if every:
        if source_id is not None or package is not None:
            raise _fail(one_of, 2)
        with _ingest_errors():
            every_result = ingest_all(root)
        for line in _ingest_all_lines(every_result):
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
    for line in _ingest_lines(result):
        typer.echo(line)


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


def _ingest_lines(result: IngestResult) -> list[str]:
    """What an ingest wrote and, by name, everything it did not; then its validation."""
    return _summary_lines(result) + _validation_lines(result.warnings)


def _ingest_all_lines(every: IngestAllResult) -> list[str]:
    """Each source's summary, every pruned row by file and gene, then the one validation."""
    lines = [line for result in every.results for line in _summary_lines(result)]
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


def _summary_lines(result: IngestResult) -> list[str]:
    """What an ingest wrote and, by name, everything it did not."""
    lines = [
        f"ingested {result.source_id} {result.version}: {result.rows_in} row(s) in the "
        f"package, {result.rows_written} written",
        f"{result.rows_outside_registry} row(s) name a gene outside the registry and were "
        f"not written",
        f"{result.rows_unresolved} row(s) name a gene id that did not resolve and were not written",
    ]
    if result.unresolved:
        names = "; ".join(
            f"{item.symbol} ({item.resolution.value}"
            + (f": {', '.join(item.candidates)}" if item.candidates else "")
            + ")"
            for item in result.unresolved
        )
        lines.append(
            f"{len(result.unresolved)} gene id(s) did not resolve and their rows were not "
            f"written: {names}"
        )
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


def _listing(finding: Finding, paper: Publication | None) -> list[str]:
    """Every field of a finding a person approves, except the bookkeeping ones.

    Omitted on purpose: `id` and the link (the heading line), `approved_*`
    (a listed finding has none) and `verified_by` (the heading line). Floats go
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
    lines = [
        f"{finding.id}  {link}  {finding.design.value}  {finding.verified_by or 'UNVERIFIED'}",
        f"  {paper_line}",
        f"  methods: {', '.join(m.value for m in finding.methods)}",
        f"  comparator: {finding.comparator.value if finding.comparator else 'none'}",
        f"  burden study: {finding.burden_study or 'none'}",
        f"  variant classes: {_ids([v.value for v in finding.variant_classes])}",
        f"  {_counts_text(finding)}",
        f"  phenotypes: {', '.join(finding.phenotypes)}",
        f"  lesion groups: {', '.join(g.value for g in finding.lesion_groups)}",
        f"  syndromic: {finding.syndromic.value}",
        f"  cohorts: {_ids(finding.cohorts)}",
        f"  locators: {'; '.join(_locator_text(loc) for loc in finding.locators)}",
        f"  shares probands with: {shared if isinstance(shared, str) else _ids(shared)}",
    ]
    if finding.review_status is not None:
        lines.append(f"  review: {finding.review_status.value}: {finding.review_note}")
    lines.append(f"  extracted by: {finding.extracted_by}")
    lines.append(f"  {finding.summary}")
    return lines


def _first_search_run_on(record: GeneLiteratureFile, finding: Finding) -> date | None:
    """The run date of the search that first found the finding's paper (the LIT008 rule)."""
    first = next((e.search for e in record.triage if e.pmid == finding.publication), None)
    return next((s.run_on for s in record.searches if s.id == first), None)


@literature_app.command("review")
def literature_review(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    approve: list[str] = typer.Option(
        [], "--approve", help="A finding id to approve; repeat the option for several."
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
    note: str = typer.Option("", "--note", help="The measured disagreement, for --pending-review."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """List the findings awaiting approval, approve some, or mark one pending review.

    A review mark is published with the finding; it does not withdraw an
    approval, and a marked finding can still be approved. A finding carries one
    mark: a second is refused rather than overwriting the first note.
    """
    _require_root(root)
    if approve and pending_review:
        raise _fail("--approve and --pending-review are two actions; use one, not both", 2)
    if not approve and (by or on or allow_unverified):
        raise _fail("--by, --on and --allow-unverified only apply with --approve", 2)
    if not pending_review and note:
        raise _fail("--note only applies with --pending-review", 2)
    if against is not None and (approve or pending_review):
        raise _fail("--against only applies to the listing, not to the two actions", 2)
    if approve:
        if not by.strip():
            raise _fail("--by (the approving person's name) is required with --approve", 2)
        if _is_agent(by):
            raise _fail(f"--by {by!r} names an agent; a person must approve", 2)
        try:
            when = date.fromisoformat(on) if on else _today()
        except ValueError:
            raise _fail(f"--on {on!r} is not a YYYY-MM-DD date", 2) from None
        if when > _today():
            raise _fail(f"--on {when.isoformat()} is in the future", 2)
    if pending_review and not note.strip():
        raise _fail("--note (the measured disagreement) is required with --pending-review", 2)

    approved_symbol, path = _require_literature_file(root, symbol)
    raw = path.read_bytes()
    try:
        record = parse_literature(raw)
    except LiteratureFileError as exc:
        raise _fail(f"{path}: {exc}", 1) from None
    by_id: dict[str, Finding] = {f.id: f for f in record.findings}

    if not approve and not pending_review:
        waiting = sorted((f for f in record.findings if f.approved_by is None), key=lambda f: f.id)
        if not waiting:
            typer.echo(f"{approved_symbol}: no findings awaiting approval")
        papers = {p.id: p for p in load_curation(root)[0].publications}
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
            }
        }
        messages = [f"marked {pending_review} pending_review"]
    else:
        ids = sorted(set(approve))
        for finding_id in ids:
            target = by_id.get(finding_id)
            if target is None:
                raise _fail(f"no such finding {finding_id} in {path}", 2)
            if target.approved_by is not None:
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
            for i in ids
        }
        messages = [f"approved {i} ({by.strip()}, {when.isoformat()})" for i in ids]
    try:
        content = with_findings_updated(raw, updates)
    except LiteratureFileError as exc:
        raise _fail(str(exc), 1) from None
    write_bytes_atomically(path, content)
    for message in messages:
        typer.echo(message)
