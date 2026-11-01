#!/usr/bin/env python3
"""Flatten `hvantk expression summarize` outputs into the pseudobulk TSV.

Runs on the cluster under the hvantk venv (it needs anndata). Reads one or both
per-group AnnData files hvantk writes -- layers `sum` and `count_nonzero`
(groups x genes), `obs['n_cells']`, the group-by columns in `obs`, gene symbols
as the var index -- and writes one row per (symbol, heart, cell type, region)
for the approved symbols in `--registry`. Aliases and previous symbols are not
read: the source carries no Ensembl id, and an alias can be an unrelated
gene's approved symbol (`HGF` is SOS1's alias), so keeping one would carry the
wrong gene's column into the pseudobulk.

`--by-region` is the run grouped by Age, Rep, major_cell_class and Region; its
rows take each group's own Region. `--pooled` is the run without Region; its
rows are tagged `all`. `library_size` is the group's `sum` row summed over every
gene in the file, before the registry filter.

The obs column names are those hvantk's documentation gives; they have not been
checked against a real output, so a missing column stops the run naming it.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys
from collections.abc import Mapping, Sequence
from collections.abc import Set as AbstractSet
from pathlib import Path
from typing import Any

GROUP_COLUMNS = ("Age", "Rep", "major_cell_class")
REGION_COLUMN = "Region"
OUTPUT_COLUMNS = (
    "symbol",
    "heart",
    "cell_type",
    "region",
    "counts",
    "count_nonzero",
    "n_cells",
    "library_size",
)
CHUNK = 128


def _whole(value: object, what: str) -> int:
    number = float(str(value))
    if number != int(number):
        raise SystemExit(f"ERROR: {what} {value!r} is not a whole number")
    return int(number)


def extract_rows(
    groups: Sequence[Mapping[str, object]],
    genes: Sequence[str],
    sums: Sequence[Sequence[float]],
    nonzero: Sequence[Sequence[float]],
    wanted: AbstractSet[str],
    region: str | None,
) -> list[dict[str, object]]:
    """Rows for the wanted symbols of each group.

    `region` is the token for every row, or None to take each group's own
    `Region`. `library_size` sums the group's whole `sums` row, so `genes` must
    be the complete gene axis, not the wanted subset.
    """
    columns = [(j, symbol) for j, symbol in enumerate(genes) if symbol in wanted]
    rows: list[dict[str, object]] = []
    for group, sum_row, nonzero_row in zip(groups, sums, nonzero, strict=True):
        age = _whole(group["Age"], "Age")
        rep = _whole(group["Rep"], "Rep")
        library = float(sum(sum_row))
        tag = str(group[REGION_COLUMN]) if region is None else region
        for j, symbol in columns:
            rows.append(
                {
                    "symbol": symbol,
                    "heart": f"{age}pcw-{rep}",
                    "cell_type": str(group["major_cell_class"]),
                    "region": tag,
                    "counts": float(sum_row[j]),
                    "count_nonzero": _whole(nonzero_row[j], "count_nonzero"),
                    "n_cells": _whole(group["n_cells"], "n_cells"),
                    "library_size": library,
                }
            )
    return rows


def read_wanted(registry: Path) -> set[str]:
    """Every approved symbol in the registry TSV, and nothing from its alias or
    previous-symbol columns."""
    wanted: set[str] = set()
    with registry.open(newline="", encoding="utf-8") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            wanted.add((record.get("symbol") or "").strip())
    wanted.discard("")
    return wanted


def _dense(block: Any) -> list[list[float]]:
    to_dense = getattr(block, "toarray", None)
    dense = to_dense() if callable(to_dense) else block
    return list(map(list, dense.tolist()))


def _rows_of(path: Path, wanted: set[str], region: str | None) -> list[dict[str, object]]:
    import anndata  # only the hvantk venv has it

    data = anndata.read_h5ad(path)
    needed = [*GROUP_COLUMNS, "n_cells", *([REGION_COLUMN] if region is None else [])]
    missing = [c for c in needed if c not in data.obs.columns]
    absent_layers = [layer for layer in ("sum", "count_nonzero") if layer not in data.layers]
    if missing or absent_layers:
        raise SystemExit(
            f"ERROR: {path}: obs columns missing {missing}, layers missing {absent_layers}; "
            f"obs has {list(data.obs.columns)}, layers {list(data.layers)}"
        )
    genes = [str(g) for g in data.var_names]
    records = data.obs[needed].to_dict("records")
    rows: list[dict[str, object]] = []
    for start in range(0, data.n_obs, CHUNK):
        stop = min(start + CHUNK, data.n_obs)
        rows.extend(
            extract_rows(
                records[start:stop],
                genes,
                _dense(data.layers["sum"][start:stop]),
                _dense(data.layers["count_nonzero"][start:stop]),
                wanted,
                region,
            )
        )
    print(f"{path}: {data.n_obs} groups x {len(genes)} genes -> {len(rows)} rows", file=sys.stderr)
    return rows


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--by-region", type=Path, required=True)
    parser.add_argument("--pooled", type=Path, required=True)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    wanted = read_wanted(args.registry)
    rows = _rows_of(args.by_region, wanted, None) + _rows_of(args.pooled, wanted, "all")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.out.with_name(args.out.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, args.out)
    print(f"wrote {len(rows)} rows to {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
