from __future__ import annotations

import os
import re
import urllib.error
from collections.abc import Sequence
from datetime import date
from pathlib import Path

import typer

from chd_atlas import pubmed
from chd_atlas.build.runner import BuildRefused, build_site
from chd_atlas.build.validity import gene_validity
from chd_atlas.corpus import load_curation
from chd_atlas.fs import write_bytes_atomically
from chd_atlas.gene_literature_io import (
    LiteratureFileError,
    literature_path,
    parse_literature,
    with_publication,
    with_search,
)
from chd_atlas.pubmed import Fetch, PubMedError, urllib_fetch
from chd_atlas.schema_export import export_schemas
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.runner import _gate_published_genes, _symbol_registry, validate_repository
from chd_atlas.vocab import StudyType, TriageDecision

app = typer.Typer(help="CHD Atlas curation tooling.", no_args_is_help=True)
schemas_app = typer.Typer(help="JSON Schema utilities.", no_args_is_help=True)
literature_app = typer.Typer(help="Literature search and curation helpers.", no_args_is_help=True)
app.add_typer(schemas_app, name="schemas")
app.add_typer(literature_app, name="literature")

# The one network hook. Commands read it at call time, so a test replaces it
# with `monkeypatch.setattr(cli, "_fetch", fake)` and nothing reaches PubMed.
_fetch: Fetch = urllib_fetch


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
    except OSError as exc:
        typer.echo(f"error: could not write schemas to {target}: {exc}")
        raise typer.Exit(code=2) from exc

    typer.echo(f"wrote {len(written)} schema(s) to {target}")


def _fail(message: str, code: int) -> typer.Exit:
    typer.echo(f"error: {message}")
    return typer.Exit(code=code)


def _api_key() -> str | None:
    """The NCBI key from the environment, or None when unset or blank."""
    return (os.environ.get("NCBI_API_KEY") or "").strip() or None


def _redact(message: str, key: str | None) -> str:
    return message.replace(key, "***") if key else message


def _require_root(root: Path) -> None:
    if not root.is_dir():
        raise _fail(f"--root {root} is not a directory", 2)


def _resolve_gene(root: Path, symbol: str) -> tuple[str, str, Sequence[str]]:
    """`(HGNC id, approved symbol, previous symbols)` for a typed symbol, or exit 2.

    Any unambiguous resolution is accepted -- approved, alias or previous
    (ruling R19) -- because curators type old names; the caller writes the
    approved one.
    """
    registry = _symbol_registry(root)
    if registry is None:
        raise _fail("mirrors/genes.tsv is missing or unreadable", 2)
    resolved = registry.resolve(symbol)
    if resolved.gene is None:
        detail = (
            f"{symbol!r} is claimed by {', '.join(resolved.candidates)}"
            if resolved.candidates
            else f"{symbol!r} is not a known gene symbol"
        )
        raise _fail(detail, 2)
    approved = registry.approved_symbol(resolved.gene)
    if approved is None:
        raise _fail(f"{resolved.gene} has no approved symbol in mirrors/genes.tsv", 2)
    return resolved.gene, approved, registry.previous_symbols(resolved.gene)


@literature_app.command("search")
def literature_search(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    run_on: str = typer.Option(
        "", "--run-on", help="Date of the search, YYYY-MM-DD (default today)."
    ),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Run the recorded PubMed search for a published gene into its literature file."""
    _require_root(root)
    try:
        when = date.fromisoformat(run_on) if run_on else date.today()
    except ValueError:
        raise _fail(f"--run-on {run_on!r} is not a YYYY-MM-DD date", 2) from None

    gene, approved, previous = _resolve_gene(root, symbol)
    corpus, _ = load_curation(root)
    if gene not in _gate_published_genes(root, corpus):
        raise _fail(
            f"{approved} ({gene}) is outside the publication gate; "
            f"its literature would reach no page",
            2,
        )

    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    if clingen is None or gencc is None:  # the gate above would have been empty
        raise _fail("the validity mirrors are unreadable", 2)
    scope = {str(entry.id) for entry in corpus.chd_scope}
    records = gene_validity(clingen, gencc, in_scope=scope)[gene].records
    diseases = sorted({record.disease_label for record in records})
    query = pubmed.build_query(approved, previous, diseases)

    path = literature_path(root, approved)
    existing = path.read_bytes() if path.is_file() else None
    key = _api_key()
    # Every network call finishes, and the new file is serialised and checked,
    # before anything is written: a failure leaves the file as it was.
    try:
        hits, pmids = pubmed.esearch_all(query, _fetch, api_key=key)
        content, search_id, new_hits = with_search(
            existing, gene=gene, query=query, run_on=when, hits=hits, pmids=pmids
        )
    except (PubMedError, urllib.error.URLError, LiteratureFileError) as exc:
        raise _fail(_redact(str(exc), key), 1) from None
    path.parent.mkdir(parents=True, exist_ok=True)
    write_bytes_atomically(path, content)
    typer.echo(f"{approved}: search {search_id} found {hits} hit(s), {new_hits} new; wrote {path}")


@literature_app.command("abstracts")
def literature_abstracts(
    symbol: str = typer.Argument(..., help="Gene symbol (approved, alias or previous)."),
    decision: str = typer.Option("pending", help="pending, included or excluded."),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Print PMID, year, title and abstract for each triage entry with a decision."""
    _require_root(root)
    try:
        wanted = TriageDecision(decision)
    except ValueError:
        raise _fail(
            f"--decision {decision!r} is not one of pending, included, excluded", 2
        ) from None
    _, approved, _ = _resolve_gene(root, symbol)
    path = literature_path(root, approved)
    if not path.is_file():
        raise _fail(f"{path} does not exist; run `chd-atlas literature search {approved}` first", 1)
    key = _api_key()
    try:
        record = parse_literature(path.read_bytes())
        pmids = [e.pmid.removeprefix("PMID:") for e in record.triage if e.decision is wanted]
        titles = pubmed.efetch_abstracts(pmids, _fetch, api_key=key)
        meta = pubmed.esummary(pmids, _fetch, api_key=key)
    except (PubMedError, urllib.error.URLError, LiteratureFileError) as exc:
        raise _fail(_redact(str(exc), key), 1) from None
    for pmid in pmids:
        found = titles.get(pmid)
        paper = meta.get(pmid)
        year = str(paper.year) if paper else "year unknown"
        if found is None:
            typer.echo(f"PMID:{pmid} · {year} · (PubMed returned no record for this PMID)")
        else:
            title, abstract = found
            typer.echo(f"PMID:{pmid} · {year} · {title}")
            typer.echo(abstract or "(no abstract in PubMed)")
        typer.echo("")


@literature_app.command("add-publication")
def literature_add_publication(
    pmid: str = typer.Argument(..., help="PubMed id, with or without the PMID: prefix."),
    study_type: str = typer.Option(
        ..., "--study-type", help="The curator's reading of the design."
    ),
    root: Path = typer.Option(Path("."), help="Repository root."),
) -> None:
    """Append a Publication record, built from PubMed, to curation/publications.yaml."""
    _require_root(root)
    digits = pmid.removeprefix("PMID:")
    if not re.fullmatch(r"[1-9]\d*", digits):
        raise _fail(f"{pmid!r} is not a PubMed id", 2)
    try:
        kind = StudyType(study_type)
    except ValueError:
        allowed = ", ".join(t.value for t in StudyType)
        raise _fail(f"--study-type {study_type!r} is not one of: {allowed}", 2) from None
    path = root / "curation" / "publications.yaml"
    if not path.is_file():
        raise _fail(f"{path} does not exist", 1)
    key = _api_key()
    try:
        found = pubmed.esummary([digits], _fetch, api_key=key).get(digits)
        if found is None:
            raise PubMedError(f"PubMed returned no record for PMID:{digits}")
        content = with_publication(path.read_bytes(), pubmed.publication_from_meta(found, kind))
    except (PubMedError, urllib.error.URLError, LiteratureFileError, ValueError) as exc:
        raise _fail(_redact(str(exc), key), 1) from None
    write_bytes_atomically(path, content)
    typer.echo(f"added PMID:{digits} to {path}")
