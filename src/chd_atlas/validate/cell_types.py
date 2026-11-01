"""Cross-checks between ``mirrors/cell_type_profiles/`` and the dataset records.

`validate_table` checks each shard against its own column declarations and
cannot see the dataset a row belongs to. A cell-type row is only meaningful
against that dataset's declarations: the cell types it names, the regions it
sampled, and how many donors it had. Eleven ways the mirror can disagree with
them, or with the rules the cell-type block prints of every value:

- **CTP001** the row's dataset has no `cell_type_profile` record.
- **CTP002** the cell type is not one the dataset declares.
- **CTP003** the region is neither one the dataset declares nor `all`.
- **CTP004** a region row has no `all` row for the same
  (dataset, gene, cell_type); `all` is the pooled value every page falls back to.
- **CTP005** `q25 <= median_abundance <= q75` does not hold, or exactly one of
  q25 and q75 is null (an interval with one end is not an interval).
- **CTP006** `n_donors` exceeds the dataset's `n_samples`.
- **CTP007** the gene is not in the gene registry. This function skips it only
  when `known_genes` is None, and flags every gene against an empty set. The
  runner passes `known_genes or None`, so there a missing, unreadable or
  header-only registry skips CTP007 and TBL008 or REF001 reports why.
- **CTP008** a `cell_type_profile` dataset declares a cell type, or a region,
  that no row of that dataset carries: a declaration with no data, as PRF004 is
  for the bulk layer. One issue per missing token, or one for the dataset when
  no row carries it at all. A missing `mirrors/cell_type_profiles/` is exactly
  that case, not a reason to skip.
- **CTP009** quartiles that contradict the hearts behind them: q25 and q75
  both present with fewer than `MIN_HEARTS_FOR_QUARTILES` (3) hearts, or both
  null with 3 or more. Exactly one null is CTP005's, so the two never fire on
  one row.
- **CTP010** `n_cells` below `MIN_CELLS_PER_HEART` (10) times `n_donors`. Every
  heart counted had at least 10 cells there, so the pool has at least ten per
  heart. Necessary, not sufficient: a pool can clear it with one heart under
  10 cells (hearts of 25 and 5 cells pool 30, over 20), and only the per-heart
  groups the converter reads can show that (`convert_farah2024.py` refuses
  them).
- **CTP011** one (dataset, gene, cell_type, region) in more than one shard. The
  build would keep one row and drop the other in silence. TBL007 already
  refuses a repeat within one shard on the same key, so this covers repeats
  across shards only.

CTP002, CTP003 and CTP006 need the dataset's record (its cell types, regions
and `n_samples`), so they are skipped for a row whose dataset is unrecorded
(CTP001 already names that). A shard that cannot be read, or lacks a column a
rule needs, is skipped: TBL000/TBL001 already report it. CTP008 is then skipped
for every dataset, since the rows it looks for may be in that shard. The other
rules report once per offending value per shard, not once per row; CTP008 once
per missing token, and CTP011 once per key at each later shard carrying it.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Collection, Sequence
from pathlib import Path

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.dataset import Dataset
from chd_atlas.tables import (
    CELL_TYPE_PROFILES,
    MIN_CELLS_PER_HEART,
    MIN_HEARTS_FOR_QUARTILES,
    POOLED_REGION,
    SHARDED_TABLES,
    mirror_paths,
    read_table,
)
from chd_atlas.vocab import DatasetDesign

_COLUMNS = [
    "dataset",
    "gene",
    "cell_type",
    "region",
    "median_abundance",
    "q25",
    "q75",
    "n_donors",
    "n_cells",
]

_Key = tuple[str, str, str, str]


def validate_cell_type_profiles(
    root: Path,
    datasets: Sequence[Dataset],
    known_genes: Collection[str] | None,
) -> list[ValidationIssue]:
    """CTP001-CTP011 over every `cell_type_profiles` shard under `root`."""
    records: dict[str, Dataset] = {
        str(d.id): d for d in datasets if d.design is DatasetDesign.CELL_TYPE_PROFILE
    }
    issues: list[ValidationIssue] = []
    complete = True
    shards: list[tuple[str, list[dict[str, object]]]] = []
    for path, name in mirror_paths(root):
        if name != CELL_TYPE_PROFILES.name:
            continue
        frame, _ = read_table(path, CELL_TYPE_PROFILES)
        if frame is None or not set(_COLUMNS) <= set(frame.columns):
            complete = False
            continue
        rows = list(frame.select(_COLUMNS).iter_rows(named=True))
        shards.append((str(path), rows))
        issues.extend(_shard_issues(str(path), rows, records, known_genes))
    issues.extend(_repeated_across_shards(shards))
    if complete:
        issues.extend(_declared_without_rows(root, shards, records))
    return issues


def _shard_issues(
    location: str,
    rows: list[dict[str, object]],
    records: dict[str, Dataset],
    known_genes: Collection[str] | None,
) -> list[ValidationIssue]:
    pooled = {
        (r["dataset"], r["gene"], r["cell_type"]) for r in rows if r["region"] == POOLED_REGION
    }
    found: dict[tuple[str, str], str] = {}

    def report(code: str, key: str, message: str) -> None:
        found.setdefault((code, key), message)

    for row in rows:
        dataset, gene = str(row["dataset"]), str(row["gene"])
        cell_type, region = str(row["cell_type"]), str(row["region"])
        where = f"{dataset}/{gene}/{cell_type}/{region}"
        record = records.get(dataset)
        if record is None:
            report(
                "CTP001",
                dataset,
                f"dataset '{dataset}' has no cell_type_profile record in datasets/",
            )
        else:
            if cell_type not in {c.token for c in record.cell_types}:
                report(
                    "CTP002",
                    f"{dataset}/{cell_type}",
                    f"cell type '{cell_type}' is not declared by dataset '{dataset}'",
                )
            if region != POOLED_REGION and region not in record.regions:
                report(
                    "CTP003",
                    f"{dataset}/{region}",
                    f"region '{region}' is not declared by dataset '{dataset}' "
                    f"(declared: {sorted(record.regions)}, plus '{POOLED_REGION}')",
                )
            n_donors = row["n_donors"]
            if isinstance(n_donors, int) and n_donors > record.n_samples:
                report(
                    "CTP006",
                    where,
                    f"{gene}/{cell_type}/{region}: n_donors {n_donors} exceeds "
                    f"dataset '{dataset}' n_samples {record.n_samples}",
                )
        if region != POOLED_REGION and (dataset, gene, cell_type) not in pooled:
            report(
                "CTP004",
                f"{dataset}/{gene}/{cell_type}",
                f"{gene}/{cell_type} in dataset '{dataset}' has region rows but no "
                f"'{POOLED_REGION}' row",
            )
        broken = _quartile_problem(row)
        if broken:
            report(
                "CTP005",
                where,
                f"{gene}/{cell_type}/{region} in dataset '{dataset}': {broken}",
            )
        unsupported = _quartiles_against_hearts(row)
        if unsupported:
            report(
                "CTP009",
                where,
                f"{gene}/{cell_type}/{region} in dataset '{dataset}': {unsupported}",
            )
        n_donors, n_cells = row["n_donors"], row["n_cells"]
        if (
            isinstance(n_donors, int)
            and isinstance(n_cells, int)
            and n_cells < MIN_CELLS_PER_HEART * n_donors
        ):
            report(
                "CTP010",
                where,
                f"{gene}/{cell_type}/{region} in dataset '{dataset}': n_cells {n_cells} is "
                f"below {MIN_CELLS_PER_HEART} per heart over n_donors {n_donors}",
            )
        if known_genes is not None and gene not in known_genes:
            report("CTP007", gene, f"gene '{gene}' is not in the gene registry")
    return [
        ValidationIssue(code, Severity.ERROR, location, message)
        for (code, _), message in sorted(found.items())
    ]


def _quartile_problem(row: dict[str, object]) -> str | None:
    """What is wrong with the row's q25/median/q75, or None if it is fine."""
    q25, median, q75 = row["q25"], row["median_abundance"], row["q75"]
    if (q25 is None) != (q75 is None):
        return "exactly one of q25 and q75 is null"
    if q25 is None or q75 is None or median is None:
        return None
    if not (q25 <= median <= q75):  # type: ignore[operator]
        return f"q25 {q25} <= median {median} <= q75 {q75} does not hold"
    return None


def _quartiles_against_hearts(row: dict[str, object]) -> str | None:
    """CTP009's finding for the row, or None.

    Only a consistent pair is judged here: exactly one null quartile is
    CTP005's, and judging it again would report one fault twice.
    """
    q25, q75, n_donors = row["q25"], row["q75"], row["n_donors"]
    if (q25 is None) != (q75 is None) or not isinstance(n_donors, int):
        return None
    enough = n_donors >= MIN_HEARTS_FOR_QUARTILES
    if q25 is not None and not enough:
        return (
            f"quartiles from {n_donors} hearts; they are published from "
            f"{MIN_HEARTS_FOR_QUARTILES} or more"
        )
    if q25 is None and enough:
        return (
            f"no quartiles from {n_donors} hearts; they are published from "
            f"{MIN_HEARTS_FOR_QUARTILES} or more"
        )
    return None


def _repeated_across_shards(
    shards: list[tuple[str, list[dict[str, object]]]],
) -> list[ValidationIssue]:
    """CTP011: a key some earlier shard already carries, reported at the later one.

    Each shard's keys are deduplicated first, so a repeat within one shard is
    left to TBL007.
    """
    first: dict[_Key, str] = {}
    issues: list[ValidationIssue] = []
    for location, rows in shards:
        keys = sorted(
            {
                (str(r["dataset"]), str(r["gene"]), str(r["cell_type"]), str(r["region"]))
                for r in rows
            }
        )
        for key in keys:
            earlier = first.setdefault(key, location)
            if earlier != location:
                dataset, gene, cell_type, region = key
                issues.append(
                    ValidationIssue(
                        "CTP011",
                        Severity.ERROR,
                        location,
                        f"{gene}/{cell_type}/{region} in dataset '{dataset}' is also in "
                        f"{earlier}; one key, one row",
                    )
                )
    return issues


def _declared_without_rows(
    root: Path,
    shards: list[tuple[str, list[dict[str, object]]]],
    records: dict[str, Dataset],
) -> list[ValidationIssue]:
    """CTP008: each declared cell type and region no row of its dataset carries.

    Sorted by dataset id, then in the dataset's declared order, so the output
    does not follow `PYTHONHASHSEED`. The location is the first shard carrying
    the dataset, or the shard its id would name by convention when none does --
    the file a curator would create, as PRF004 does.
    """
    cell_types: dict[str, set[str]] = defaultdict(set)
    regions: dict[str, set[str]] = defaultdict(set)
    located: dict[str, str] = {}
    for location, rows in shards:
        for row in rows:
            dataset = str(row["dataset"])
            cell_types[dataset].add(str(row["cell_type"]))
            regions[dataset].add(str(row["region"]))
            located.setdefault(dataset, location)
    directory = root / "mirrors" / SHARDED_TABLES[CELL_TYPE_PROFILES.name]
    issues: list[ValidationIssue] = []
    for dataset_id in sorted(records):
        record = records[dataset_id]
        if dataset_id not in located:
            messages = [
                f"dataset '{dataset_id}' is a cell_type_profile dataset and no "
                "cell_type_profiles row carries it"
            ]
        else:
            messages = [
                f"dataset '{dataset_id}' declares cell type '{c.token}', carried by none "
                "of its rows"
                for c in record.cell_types
                if c.token not in cell_types[dataset_id]
            ] + [
                f"dataset '{dataset_id}' declares region '{region}', carried by none of its rows"
                for region in record.regions
                if region not in regions[dataset_id]
            ]
        location = located.get(dataset_id, str(directory / f"{dataset_id}.tsv"))
        issues.extend(
            ValidationIssue("CTP008", Severity.ERROR, location, message) for message in messages
        )
    return issues
