"""Shared fixtures for the package-manifest tests (feeds phase 1)."""

from __future__ import annotations

from typing import Any


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
