# src/chd_atlas/build/cell_types.py
"""A gene's cell-type expression, assembled from the `cell_type_profiles` mirror.

One entry per `cell_type_profile` dataset that has rows for the gene. Each
entry lists the dataset's declared cell types that have data, in the order the
dataset declares them, with the pooled (`region == "all"`) values and one
entry per declared region that has a row, also in declared order. Both orders
come from the dataset record, never from a set or from file order, so the
output does not move with `PYTHONHASHSEED` (CLAUDE.md section 4.42).

**Absence is not a low value.** A cell type with no pooled row for a gene was
not measured in any heart with at least 10 cells of it, so it is omitted from
`cell_types`; a region with no row is omitted from `regions` the same way. The
only entries that say "low" are those with `below_floor` true. A page that
wants every declared cell type lists the dataset's own declaration and writes
"no data" where this payload has no entry.

Tau is computed over the pooled medians of the cell types present, through
`profiles.specificity`, the same function the bulk layer uses, so the two
cannot disagree on the arithmetic. Two things differ: its `method` string says
"organs", which is false for cell types, and its `tissues` are sorted, where
`specificity.cell_types` lists the present cell types in declared order, the
order every other list in this payload and on the page follows.
"""

from __future__ import annotations

from collections.abc import Collection, Mapping, Sequence
from pathlib import Path
from typing import Any, Final, TypedDict

from chd_atlas.build.profiles import TAU_METHOD, specificity
from chd_atlas.models.dataset import Dataset
from chd_atlas.tables import CELL_TYPE_PROFILES, POOLED_REGION, mirror_paths, read_table
from chd_atlas.vocab import DatasetDesign

# `TAU_METHOD` with "organs"/"organ's" replaced; pinned in full by
# `test_the_tau_method_names_cell_types_not_organs`.
CELL_TYPE_TAU_METHOD: Final = TAU_METHOD.replace("organs", "cell types").replace(
    "organ's", "cell type's"
)


class CellTypeValues(TypedDict):
    """The numbers for one cell type, pooled or in one region."""

    median: float
    q25: float | None
    q75: float | None
    n_donors: int
    n_cells: int
    fraction_expressing: float
    below_floor: bool


class RegionEntry(CellTypeValues):
    """`CellTypeValues` for one region."""

    region: str


class CellTypeEntry(TypedDict):
    """One cell type: its pooled values and its per-region values."""

    cell_type: str
    label: str
    pooled: CellTypeValues
    regions: list[RegionEntry]


class CellTypeSpecificity(TypedDict):
    """Tau over the pooled medians, with what is needed to re-derive it."""

    tau: float
    scale: str
    method: str
    cell_types: tuple[str, ...]
    n_cell_types: int
    highest_in: str | None
    medians: dict[str, float]


class CellTypeDataset(TypedDict):
    """One dataset's cell-type expression for one gene."""

    dataset: str
    unit: str
    detection_floor: float
    n_donors: int
    specificity: CellTypeSpecificity | None
    cell_types: list[CellTypeEntry]


class CellTypeExpression(TypedDict):
    """One gene's cell-type expression; `datasets` is empty if no mirror covers it."""

    datasets: list[CellTypeDataset]


EMPTY_CELL_TYPE_EXPRESSION: Final[CellTypeExpression] = CellTypeExpression(datasets=[])


class CellTypeCensus(TypedDict):
    """What the cell-type layer this site publishes amounts to."""

    genes: int
    datasets: int


def _values(record: Mapping[str, Any], floor: float) -> CellTypeValues:
    median = float(record["median_abundance"])
    return CellTypeValues(
        median=median,
        q25=float(record["q25"]) if record["q25"] is not None else None,
        q75=float(record["q75"]) if record["q75"] is not None else None,
        n_donors=int(record["n_donors"]),
        n_cells=int(record["n_cells"]),
        fraction_expressing=float(record["fraction_expressing"]),
        below_floor=median < floor,
    )


def _dataset_entry(
    dataset: Dataset, unit: str, rows: Mapping[tuple[str, str], Mapping[str, Any]]
) -> CellTypeDataset:
    """One gene's entry for one dataset; `rows` is keyed by (cell_type, region)."""
    floor = dataset.detection_floor
    if floor is None:
        # The model refuses a cell_type_profile dataset without one, so only a
        # bypassed gate reaches here; `raise`, because `-O` strips an `assert`.
        raise ValueError(f"dataset {str(dataset.id)!r} has no detection_floor")
    entries: list[CellTypeEntry] = []
    medians: dict[str, float] = {}
    for declared in dataset.cell_types:
        pooled = rows.get((declared.token, POOLED_REGION))
        if pooled is None:
            continue
        values = _values(pooled, floor)
        medians[declared.token] = values["median"]
        regions = [
            RegionEntry(region=region, **_values(rows[declared.token, region], floor))
            for region in dataset.regions
            if (declared.token, region) in rows
        ]
        entries.append(
            CellTypeEntry(
                cell_type=declared.token, label=declared.label, pooled=values, regions=regions
            )
        )
    tau = specificity(medians, floor)
    return CellTypeDataset(
        dataset=str(dataset.id),
        unit=unit,
        detection_floor=floor,
        n_donors=dataset.n_samples,
        specificity=(
            None
            if tau is None
            else CellTypeSpecificity(
                tau=tau["tau"],
                scale=tau["scale"],
                method=CELL_TYPE_TAU_METHOD,
                cell_types=tuple(c.token for c in dataset.cell_types if c.token in medians),
                n_cell_types=tau["n_tissues"],
                highest_in=tau["highest_in"],
                medians=tau["medians"],
            )
        ),
        cell_types=entries,
    )


def gene_cell_type_expression(
    root: Path, datasets: Sequence[Dataset]
) -> dict[str, CellTypeExpression]:
    """Every gene the `cell_type_profiles` mirror mentions, "published or not".

    A gene's `datasets` are ordered by dataset id. Rows for a dataset that is
    not a registered `cell_type_profile` dataset, and rows naming a cell type
    or region the dataset does not declare, are skipped: `validate_table` and
    `validate_cell_type_profiles` report them against the same files. A key
    repeated across rows raises: TBL007 and CTP011 refuse it, and keeping
    either row would publish one value and drop the other without a word.
    """
    by_id = {
        str(dataset.id): dataset
        for dataset in datasets
        if dataset.design == DatasetDesign.CELL_TYPE_PROFILE
    }
    # gene -> dataset -> (cell_type, region) -> row. `units` collects every
    # distinct unit over all of a dataset's rows, across genes. TBL004 holds
    # every row to `cpm`, so a mixed dataset reaches here only past a bypassed
    # validator; the check below refuses it rather than publish one unit for
    # values measured in another (CLAUDE.md section 9: a bypassed gate fails).
    units: dict[str, set[str]] = {}
    rows: dict[str, dict[str, dict[tuple[str, str], Mapping[str, Any]]]] = {}
    for path, schema_name in mirror_paths(root):
        if schema_name != CELL_TYPE_PROFILES.name:
            continue
        frame, _ = read_table(path, CELL_TYPE_PROFILES)
        if frame is None:
            continue
        for record in frame.iter_rows(named=True):
            dataset = by_id.get(str(record["dataset"]))
            if dataset is None or record["gene"] is None:
                continue
            units.setdefault(str(dataset.id), set()).add(str(record["unit"]))
            tokens = {c.token for c in dataset.cell_types}
            region = str(record["region"])
            if record["cell_type"] not in tokens or (
                region != POOLED_REGION and region not in dataset.regions
            ):
                continue
            cells = rows.setdefault(str(record["gene"]), {}).setdefault(str(dataset.id), {})
            key = (str(record["cell_type"]), region)
            if key in cells:
                raise ValueError(
                    f"dataset {str(dataset.id)!r} repeats {record['gene']}/{key[0]}/{key[1]} "
                    f"in cell_type_profiles (again in {path})"
                )
            cells[key] = record

    for dataset_id, found in sorted(units.items()):
        if len(found) > 1:
            raise ValueError(
                f"dataset {dataset_id!r} mixes units in cell_type_profiles: {sorted(found)}"
            )

    result: dict[str, CellTypeExpression] = {}
    for gene, per_dataset in rows.items():
        entries = []
        for dataset_id in sorted(per_dataset):
            cells = per_dataset[dataset_id]
            (unit,) = units[dataset_id]
            entries.append(_dataset_entry(by_id[dataset_id], unit, cells))
        result[gene] = CellTypeExpression(datasets=entries)
    return result


def cell_type_census(
    expression: Mapping[str, CellTypeExpression], published: Collection[str]
) -> CellTypeCensus:
    """What the cell-type layer *this site publishes* amounts to.

    Both figures are restricted to `published`. `genes` counts published genes
    whose `datasets` list is non-empty, never `len(expression)`:
    `gene_cell_type_expression` returns every gene the mirror mentions, whether
    published or not, so counting the mapping would advertise cell-type
    evidence behind genes this site publishes no bundle for. `datasets` counts
    the distinct dataset ids named by those same genes' entries, never the
    number of registered cell-type datasets.
    """
    genes = 0
    datasets: set[str] = set()
    for gene in published:
        entry = expression.get(gene)
        if entry is None or not entry["datasets"]:
            continue
        genes += 1
        datasets.update(item["dataset"] for item in entry["datasets"])
    return CellTypeCensus(genes=genes, datasets=len(datasets))
