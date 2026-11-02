# src/chd_atlas/gene_literature_io.py
"""Reading and writing the curator-edited literature files without losing a byte
of what the curator wrote.

Two files, two strategies, chosen by what each needs to preserve:

- ``curation/literature/<SYMBOL>.yaml`` is rewritten through ruamel's
  round-trip loader, so a curator's comments and the existing entries survive a
  re-run of ``search``. New entries are appended; nothing is reordered.
- ``curation/publications.yaml`` is *appended to as text*: its comments sit
  between records, folded titles are wrapped by hand, and a round trip would
  re-flow them. Every existing byte stays; one record is added at the end
  (ruling R7 -- the file is in curation order, not sorted, and no validator
  sorts it).

Both return **bytes** and write nothing: the caller serialises, validates and
only then writes (`fs.write_bytes_atomically`), so a failure anywhere leaves the
file as it was.
"""

from __future__ import annotations

from datetime import date
from io import BytesIO
from pathlib import Path
from typing import Any

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import Publication, PublicationFile


class LiteratureFileError(Exception):
    """A literature or publications file cannot be safely extended."""


def literature_path(root: Path, symbol: str) -> Path:
    return root / "curation" / "literature" / f"{symbol}.yaml"


def _round_trip_yaml() -> YAML:
    yaml = YAML()
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.preserve_quotes = True
    # Wide enough that a long recorded query is never folded, and an existing
    # line is never re-wrapped on the way through.
    yaml.width = 4096
    return yaml


def dump_bytes(data: Any) -> bytes:
    buffer = BytesIO()
    _round_trip_yaml().dump(data, buffer)
    return buffer.getvalue()


def load_round_trip(raw: bytes) -> CommentedMap:
    try:
        data = _round_trip_yaml().load(raw)
    except Exception as exc:  # ruamel raises several YAMLError subclasses
        raise LiteratureFileError(f"not valid YAML: {exc}") from exc
    if not isinstance(data, CommentedMap):
        raise LiteratureFileError("the top level is not a mapping")
    return data


def parse_literature(raw: bytes) -> GeneLiteratureFile:
    """The file as the model sees it, or `LiteratureFileError` naming why not."""
    try:
        data = YAML(typ="safe").load(raw)
        return GeneLiteratureFile.model_validate(data)
    except Exception as exc:  # YAMLError or pydantic.ValidationError
        raise LiteratureFileError(str(exc)) from exc


def pmid_key(pmid: str) -> int:
    """Numeric order of a `PMID:123` string (string order puts 9 after 10)."""
    return int(pmid.removeprefix("PMID:"))


def with_search(
    existing: bytes | None,
    *,
    gene: str,
    query: str,
    run_on: date,
    hits: int,
    pmids: list[str],
) -> tuple[bytes, int, int]:
    """The file's new content after recording one search: `(bytes, search_id, new_hits)`.

    `pmids` are the bare PMIDs esearch returned. Only those not already in the
    file's triage become entries, `pending` and attributed to the new search,
    appended in numeric PMID order (ruling R8). The search id is the largest
    existing id plus one. The result is checked through `GeneLiteratureFile`
    before it is returned.
    """
    if existing is None:
        data: CommentedMap = CommentedMap()
        data["gene"] = gene
        data["searches"] = CommentedSeq()
        data["triage"] = CommentedSeq()
    else:
        data = load_round_trip(existing)
        if data.get("gene") != gene:
            raise LiteratureFileError(f"the file's gene is {data.get('gene')!r}, not {gene}")
    searches = data.setdefault("searches", CommentedSeq())
    triage = data.setdefault("triage", CommentedSeq())
    for seq in (searches, triage):
        seq.fa.set_block_style()

    search_id = max((int(s["id"]) for s in searches), default=0) + 1
    known = {str(t["pmid"]) for t in triage}
    fresh = sorted((f"PMID:{p}" for p in set(pmids)), key=pmid_key)
    fresh = [p for p in fresh if p not in known]

    searches.append(
        {
            "id": search_id,
            "source": "pubmed",
            "query": query,
            "run_on": run_on,
            "hits": hits,
            "new_hits": len(fresh),
        }
    )
    for pmid in fresh:
        triage.append({"pmid": pmid, "search": search_id, "decision": "pending"})

    content = dump_bytes(data)
    parse_literature(content)
    return content, search_id, len(fresh)


def with_publication(existing: bytes, publication: Publication) -> bytes:
    """`publications.yaml` with one record appended; every existing byte kept."""
    try:
        before = PublicationFile.model_validate(YAML(typ="safe").load(existing))
    except Exception as exc:
        raise LiteratureFileError(f"publications.yaml does not validate: {exc}") from exc
    if any(p.id == publication.id for p in before.publications):
        raise LiteratureFileError(f"{publication.id} is already in publications.yaml")

    record: dict[str, Any] = {
        "id": publication.id,
        "title": publication.title,
        "journal": publication.journal,
        "year": publication.year,
        "authors": list(publication.authors),
        "study_type": publication.study_type.value,
    }
    for field in ("doi", "pmcid"):
        value = getattr(publication, field)
        if value is not None:
            record[field] = value
    text = dump_bytes({"publications": [record]}).decode("utf-8")
    snippet = text.split("\n", 1)[1]  # drop the "publications:" line, keep the item's indent

    separator = b"" if existing.endswith(b"\n") else b"\n"
    content = existing + separator + snippet.encode("utf-8")

    try:
        after = PublicationFile.model_validate(YAML(typ="safe").load(content))
    except Exception as exc:
        raise LiteratureFileError(f"the appended record does not parse: {exc}") from exc
    if (
        len(after.publications) != len(before.publications) + 1
        or after.publications[-1] != publication
    ):
        raise LiteratureFileError("the appended record did not land as the last publication")
    return content
