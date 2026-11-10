"""Shared fixtures for the package-manifest tests (feeds phase 1)."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

from chd_atlas.tables import TABLE_SCHEMAS


def manifest_dict(**overrides: Any) -> dict[str, Any]:
    """A complete, valid manifest as a fresh dict; top-level keys can be overridden."""
    manifest: dict[str, Any] = {
        "package": "toy_burden",
        "version": "2026-11-10.1",
        "ids": "symbol",
        "coverage": "every_gene",
        "tables": [
            {
                "name": "burden",
                "format_version": 1,
                "file": "burden.tsv",
                "sha256": "b" * 64,
                "rows": 3,
                "genes": 2,
                "shard": None,
            }
        ],
        "inputs": [
            {
                "name": "supplement_s1",
                "url": "https://example.org/s1.xlsx",
                "sha256": "c" * 64,
                "retrieved_on": "2026-11-01",
                "upstream_release": None,
            }
        ],
        "input_counts": {"rows": 10, "genes": 4},
        "dropped": [{"reason": "symbol not in HGNC", "count": 2, "examples": ["FOO1"]}],
        "built_by": {
            "recipe_commit": "a" * 40,
            "hvantk_release": None,
            "image_digest": None,
            "host": "workstation",
            "job": None,
        },
        "source": {
            "name": "Toy burden study",
            "version": "2026-11-10.1",
            "retrieved_on": "2026-11-01",
            "url": "https://example.org/toy",
            "licence": "CC-BY-4.0",
            "redistribution": "permitted",
            "file_version": None,
            "file_date": None,
        },
        "licence_evidence": "Supplement licence statement, https://example.org/licence",
    }
    manifest.update(overrides)
    return manifest


_METHOD_NOTE = "Cases are the Pediatric Cardiac Genomics Consortium CHD GENES cohort."


def burden_row(**cells: str) -> dict[str, str]:
    """One burden row as ``{column: text}``; unnamed columns are left blank.

    Values are copied from the first PMID:40127276 rows of ``mirrors/burden.tsv``
    (a de novo and a case-control arm), with the note shortened.
    """
    row = {
        "study": "PMID:40127276",
        "gene": "HGNC:10023",
        "cohort_stratum": "all",
        "variant_class": "snv_indel",
        "consequence_class": "damaging",
        "origin": "de_novo",
        "count_unit": "de_novo_mutations",
        "n_case_carriers": "2",
        "n_cases": "3887",
        "comparator": "mutation_model",
        "expected_count": "0.0357604",
        "effect": "55.9277860426617",
        "effect_measure": "enrichment_ratio",
        "pvalue": "0.000624",
        "pvalue_test": "poisson",
        "pvalue_adjusted": "0.00583259259259259",
        "pvalue_adjustment": "benjamini_hochberg",
        "case_cohorts": "pcgc",
        "method_note": _METHOD_NOTE,
        "source": "toy_burden",
    }
    row.update(cells)
    return row


def default_burden_rows() -> list[dict[str, str]]:
    """Three rows over two genes (the manifest fixture's ``rows: 3, genes: 2``)."""
    return [
        burden_row(),
        burden_row(consequence_class="lof", n_case_carriers="0", effect="0.0", pvalue="1.0"),
        burden_row(
            gene="HGNC:11604",
            origin="transmitted_or_unphased",
            maf_max="1e-05",
            count_unit="alleles",
            n_case_carriers="9",
            n_cases="22558",
            comparator="control_cohort",
            n_control_carriers="41",
            n_controls="270661",
            effect_measure="odds_ratio",
            effect="2.63442751348609",
            pvalue_test="fisher_exact",
            control_cohorts="gnomad_controls",
        ),
    ]


def tsv_bytes(rows: Sequence[Mapping[str, str]], columns: Sequence[str]) -> bytes:
    lines = ["\t".join(columns)]
    lines += ["\t".join(row.get(column, "") for column in columns) for row in rows]
    return ("\n".join(lines) + "\n").encode("utf-8")


def write_package(
    tmp_path: Path,
    rows: Sequence[Mapping[str, str]] | None = None,
    *,
    gz: bool = False,
    columns: Sequence[str] | None = None,
    **manifest_overrides: Any,
) -> Path:
    """Write a burden package (TSV plus ``manifest.json``) under ``tmp_path / "pkg"``.

    The manifest's table entry carries the sha256 of the bytes on disk, the data
    row count and the distinct non-blank ``gene`` values, so the package is valid
    as written; a test then corrupts exactly one thing with ``edit_table``,
    ``replace_table_bytes`` or by appending to the file. ``columns`` overrides
    the header (default: the burden schema's, in schema order).
    """
    chosen = list(rows) if rows is not None else default_burden_rows()
    header = list(columns) if columns is not None else list(TABLE_SCHEMAS["burden"].column_names)
    directory = tmp_path / "pkg"
    directory.mkdir()
    name = "burden.tsv.gz" if gz else "burden.tsv"
    raw = tsv_bytes(chosen, header)
    on_disk = gzip.compress(raw, mtime=0) if gz else raw
    (directory / name).write_bytes(on_disk)
    manifest = manifest_dict(**manifest_overrides)
    manifest["tables"] = manifest_overrides.get("tables") or [
        {
            **manifest["tables"][0],
            "file": name,
            "sha256": hashlib.sha256(on_disk).hexdigest(),
            "rows": len(chosen),
            "genes": len({row["gene"].strip() for row in chosen if row.get("gene", "").strip()}),
        }
    ]
    (directory / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return directory


def edit_table(directory: Path, **fields: Any) -> None:
    """Rewrite fields of the first table entry in ``manifest.json``."""
    path = directory / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["tables"][0].update(fields)
    path.write_text(json.dumps(manifest), encoding="utf-8")


def replace_table_bytes(directory: Path, data: bytes) -> None:
    """Replace the table file's bytes and re-hash, so only the content is wrong."""
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    (directory / manifest["tables"][0]["file"]).write_bytes(data)
    edit_table(directory, sha256=hashlib.sha256(data).hexdigest())


def write_hgnc_root(
    root: Path,
    rows: Sequence[Mapping[str, str]],
    *,
    version: str = "2026-10-16",
    gz: bool = False,
    pinned: bool = True,
    cached: bool = True,
    edit_manifest: Callable[[dict[str, Any]], None] | None = None,
    raw_manifest: str | None = None,
    pin_sha256: str | None = None,
) -> Path:
    """A temp repository root holding an ``hgnc`` source, its pin and a cached HGNC package.

    Writes ``local/packages/hgnc/<version>/`` holding ``manifest.json`` and a
    ``genes`` table of ``rows`` (when ``cached``), then ``mirrors/sources.yaml``
    (one entry, with a ``package`` pin when ``pinned``). The pin is the real
    sha256 of the manifest bytes written, so a test damages the cache *with* a
    matching pin; ``pin_sha256`` overrides it to model a stale or swapped cache.
    ``edit_manifest`` mutates the manifest dict before writing and
    ``raw_manifest`` replaces its text. Never touches the real ``local/``: ``root``
    is the caller's ``tmp_path``.
    """
    manifest_sha = "0" * 64
    if cached:
        directory = root / "local" / "packages" / "hgnc" / version
        directory.mkdir(parents=True)
        raw = tsv_bytes(rows, TABLE_SCHEMAS["genes"].column_names)
        name = "genes.tsv.gz" if gz else "genes.tsv"
        on_disk = gzip.compress(raw, mtime=0) if gz else raw
        (directory / name).write_bytes(on_disk)
        manifest = manifest_dict(
            package="hgnc",
            version=version,
            ids="hgnc_id",
            tables=[
                {
                    "name": "genes",
                    "format_version": 1,
                    "file": name,
                    "sha256": hashlib.sha256(on_disk).hexdigest(),
                    "rows": len(rows),
                    "genes": len(rows),
                    "shard": None,
                }
            ],
        )
        if edit_manifest is not None:
            edit_manifest(manifest)
        text = raw_manifest if raw_manifest is not None else json.dumps(manifest)
        (directory / "manifest.json").write_bytes(text.encode("utf-8"))
        manifest_sha = hashlib.sha256(text.encode("utf-8")).hexdigest()
    mirrors = root / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    pin = (
        f'    package: {{version: "{version}", manifest_sha256: "{pin_sha256 or manifest_sha}"}}\n'
        if pinned
        else ""
    )
    (mirrors / "sources.yaml").write_text(
        "sources:\n"
        "  - id: hgnc\n"
        "    name: HUGO Gene Nomenclature Committee\n"
        f'    version: "{version}"\n'
        "    retrieved_on: 2026-10-25\n"
        "    url: https://www.genenames.org/\n"
        "    licence: CC0-1.0\n"
        "    redistribution: permitted\n" + pin,
        encoding="utf-8",
    )
    return root
