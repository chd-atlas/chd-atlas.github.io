# src/chd_atlas/validate/profiles.py
"""Cross-table rules for ``mirrors/profiles/`` and ``mirrors/profile_quantiles/``.

`validate_table` checks each shard against its own column declarations and
cannot see that the two tables are two views of one measurement: a gene's
`median_abundance` only means something once it is placed inside the
distribution `profile_quantiles` describes for that same
`(dataset, tissue, stage)`, and that placement is only valid when both tables
agree on what they are measuring. Four ways that agreement can silently break:

- **The unit.** `unit` is `allowed` per cell (TBL005), not per dataset, so
  nothing stops one dataset reporting `rpkm` for one organ and `tpm` for
  another (**PRF008**), or `profiles` and `profile_quantiles` naming the same
  dataset in two different units (**PRF001**). Either way, a median gets read
  against a ruler marked in the wrong quantity -- the 40th percentile of a
  TPM distribution is not the 40th percentile of an RPKM one.
- **The cell.** A gene's median is placed against the breakpoints for its own
  `(dataset, tissue, stage)`. A triple with `profiles` rows and no
  `profile_quantiles` rows (**PRF002**) has nothing to be placed against.
- **The ruler.** `profile_quantiles` publishes 101 breakpoints per triple,
  meant to be read as a monotone lookup table. A breakpoint that decreases
  (**PRF003**) means the 40th percentile is a smaller number than the 30th,
  and a gene whose value straddles the inversion gets a percentile that
  depends on which side of a coin flip the lookup lands on.

**`profiles.stage` is nullable; `profile_quantiles.stage` is not** -- a
quantile grid always names the stage it summarises, so a null-stage `profiles`
row can never have a quantile partner *by construction*. PRF002 skips those
rows rather than reporting a design fact as a curation error once per tissue
on every dataset that has one. The row itself is not this function's to drop;
whatever later joins the two tables is expected to left-join so the row still
publishes, with a null percentile. Both halves matter for the same reason a
skipped check and a dropped row are different failures: skip the check alone
and a real gap still reads the same as this permitted one; keep the row alone
and there is nothing here to tell the two apart, which is how a curator learns
to stop reading PRF002 at all.

**A `profile_quantiles` mirror that is missing or entirely unreadable is not
reported here.** Comparing `profiles` against a table that never loaded would
flag every `(dataset, tissue, stage)` in `profiles` as a gap and every
dataset's unit as disagreeing -- the same cascade `REF000`/`SRC000`/`ONT000`
exist to prevent for a curated registry that failed to load. When no
`profile_quantiles` shard could be read at all, PRF001 and PRF002 are silently
skipped here rather than cascading; naming that skip so a report cannot go
quietly green over it is a separate check, added together with the runner
wiring that calls this module.

Reference checks against other registries (a gene id, a declared cardiac
tissue or stage token, a curated phase vocabulary), and the two checks a
curated dataset's own stage chronology needs (**PRF011**, **PRF012**,
**PRF014**), live in
`validate_profile_references` below, kept separate for the same reason
`validate_burden_references` is separate from `validate_burden`: one
unreadable registry must not be able to take these table-internal checks
down with it.

Every check here reports rather than raises, per the project's split:
validators collect issues, and `build_site` is what refuses.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile
from chd_atlas.tables import PROFILE_QUANTILES, PROFILES, TableSchema, mirror_paths, read_table

# A private alias purely for readability below: every check in this module is
# keyed on the same triple.
_Triple = tuple[str, str, str]

# One merged run of cardiac-phase coverage: (start_wpc, end_wpc, the sorted
# ids of every phase that merged into it). See `_coverage_spans`.
_CoverageSpan = tuple[float, float, tuple[str, ...]]

_PROFILES_COLUMNS: tuple[str, ...] = ("dataset", "tissue", "stage", "unit")
_QUANTILES_COLUMNS: tuple[str, ...] = ("dataset", "tissue", "stage", "percentile", "value", "unit")
# Columns `validate_profile_references` needs from each table -- distinct from
# the two constants above (`validate_profiles`'s own needs): it reads `gene`,
# which internal consistency never touches, and does not need `unit`/
# `percentile`/`value` at all.
_PROFILES_REFERENCE_COLUMNS: tuple[str, ...] = ("dataset", "gene", "tissue", "stage")
_QUANTILES_REFERENCE_COLUMNS: tuple[str, ...] = ("dataset", "tissue", "stage")


def validate_profiles(root: Path) -> list[ValidationIssue]:
    """Check `profiles` and `profile_quantiles` for internal consistency.

    Reads every shard under both directories via `mirror_paths` and
    `read_table`. A shard that cannot be read, or that is missing a column
    this function needs, is skipped rather than aborting the rest of the
    corpus -- `validate_table` already reports it as TBL000/TBL001 against
    the same path, so the build refuses either way and a second report here
    would only say it twice, per `validate_burden`'s precedent.

    Distinct `(dataset, tissue, stage)` triples and distinct
    `(dataset, unit)` facts are checked once, not once per row: a
    15,000-row shard with one bad triple yields one issue, matching
    `validate_burden_references`.
    """
    # dataset -> distinct profiles units, and a shard path to blame for each
    # dataset (the first, in `mirror_paths`' sorted order, that named it).
    profile_units: dict[str, set[str]] = defaultdict(set)
    profile_unit_path: dict[str, Path] = {}
    # (dataset, tissue, stage) seen in a profiles row whose stage is not null,
    # and a shard path to blame for each triple.
    profile_triples: set[_Triple] = set()
    profile_triple_path: dict[_Triple, Path] = {}

    for path, row in _shard_rows(root, PROFILES, "profiles", _PROFILES_COLUMNS):
        dataset, tissue, stage, unit = row["dataset"], row["tissue"], row["stage"], row["unit"]
        if dataset is None:
            # A null on a non-nullable column is TBL003's to report; nothing
            # useful can be attributed to a row with no dataset at all.
            continue
        if unit is not None:
            profile_units[dataset].add(unit)
            profile_unit_path.setdefault(dataset, path)
        # The stage-nullability asymmetry with `profile_quantiles` (see the
        # module docstring) means a null-stage row can never have a quantile
        # partner by construction, so it must never enter `profile_triples`.
        # Its unit still counts towards PRF001/PRF008 above: which stages a
        # dataset covers and which units it reports are independent facts.
        if tissue is not None and stage is not None:
            triple = (dataset, tissue, stage)
            profile_triples.add(triple)
            profile_triple_path.setdefault(triple, path)

    # dataset -> distinct profile_quantiles units.
    quantile_units: dict[str, set[str]] = defaultdict(set)
    # (dataset, tissue, stage) seen in any profile_quantiles row.
    quantile_triples: set[_Triple] = set()
    # (dataset, tissue, stage) -> [(percentile, value), ...], and a shard path
    # to blame for each triple.
    breakpoints: dict[_Triple, list[tuple[int, float]]] = defaultdict(list)
    breakpoint_path: dict[_Triple, Path] = {}
    # True the moment any profile_quantiles shard yields at least one usable
    # row. Stays False for "no shard exists", "every shard was unreadable or
    # missing a needed column", and the edge case of a shard that read cleanly
    # but has zero data rows -- all three mean there is nothing here to
    # compare `profiles` against, which PRF001/PRF002 must not mistake for
    # "compared, and everything is missing".
    quantile_available = False

    for path, row in _shard_rows(root, PROFILE_QUANTILES, "profile_quantiles", _QUANTILES_COLUMNS):
        quantile_available = True
        dataset, tissue, stage = row["dataset"], row["tissue"], row["stage"]
        percentile, value, unit = row["percentile"], row["value"], row["unit"]
        if dataset is None or tissue is None or stage is None:
            continue
        triple = (dataset, tissue, stage)
        quantile_triples.add(triple)
        if unit is not None:
            quantile_units[dataset].add(unit)
        if percentile is not None and value is not None:
            breakpoints[triple].append((percentile, value))
            breakpoint_path.setdefault(triple, path)

    issues: list[ValidationIssue] = []
    issues.extend(_prf008_issues(profile_units, profile_unit_path))
    issues.extend(
        _prf001_issues(profile_units, profile_unit_path, quantile_units, quantile_available)
    )
    issues.extend(
        _prf002_issues(profile_triples, profile_triple_path, quantile_triples, quantile_available)
    )
    issues.extend(_prf003_issues(breakpoints, breakpoint_path))
    return issues


def validate_profile_references(
    root: Path,
    datasets: tuple[Dataset, ...],
    known_genes: set[str] | None,
    published_genes: set[str],
    phases: CardiacPhaseFile | None,
) -> list[ValidationIssue]:
    """Check profiles rows against the registries that give them meaning.

    Separate from `validate_profiles` for the reason BUR009-011 are separate
    from BUR001-008: one unreadable registry must not silently take the
    consistency checks with it. This function does its own pass over
    `profiles`/`profile_quantiles` rather than sharing state with
    `validate_profiles`, matching how `validate_burden_references` re-reads
    `mirrors/burden.tsv` independently of `validate_burden`.

    `known_genes` is `None` when the gene registry did not load, and every
    gene check (PRF007) is then skipped rather than reporting each row as
    dangling -- the cascade REF000/SRC000/ONT000 exist to prevent. It is
    deliberately distinct from an empty-but-loaded registry (a header-only
    `mirrors/genes.tsv`), which correctly flags every gene as unknown; the
    skip is an `is not None` check for exactly that reason, never a truthy
    one, since an empty `set` is falsy.

    `published_genes` is a different set from `known_genes` and the two must
    not be conflated: `known_genes` is the whole gene registry (every gene
    this atlas has an id for), `published_genes` is the population that
    clears the publication gate (`build.validity.published_genes()`'s
    return). PRF007 checks the profiles column against the registry, because
    an unregistered id is a shape problem regardless of what publishes;
    PRF009 checks a gap against the published population, because a gene
    absent from a cell only matters to a reader once it has a page to be
    absent *from*.

    Eight codes:

    - **PRF004** (ERROR) -- a dataset's declared `cardiac_tissues` token
      matched by no `profiles` row for that dataset.
    - **PRF005** (ERROR) -- a `profiles.stage` token a dataset's own record
      does not declare in `stages`.
    - **PRF006** (WARNING) -- an interior gap no declared cardiac phase
      covers; see `_prf006_issues` for why only interior gaps, and only a
      warning.
    - **PRF007** (ERROR) -- a `profiles.gene` id in no gene registry.
    - **PRF009** (WARNING) -- a published gene missing a `profiles` row in a
      `(dataset, tissue, stage)` cell that `profile_quantiles` shows was
      actually assayed; see `_prf009_issues` for the grouping choice.
    - **PRF011** (ERROR) -- two stages of one dataset claiming one `order`.
    - **PRF012** (ERROR) -- a dataset's `order` contradicting its own `wpc`
      about which of two stages is earlier.
    - **PRF014** (ERROR) -- a stage with no `wpc` ordered before a stage that
      has one; the half of the chronology PRF012 cannot reach, since it
      compares only across stages carrying both.

    PRF011, PRF012 and PRF014 read no mirror at all -- they compare a curated
    dataset record against itself -- and are here rather than in
    `validate_profiles` because this is the function that receives
    `datasets`. Deliberately *one*
    call site, added to the one the runner already reaches on both the
    mirror-readable and mirrors-missing branches (see the runner's own
    comment, and section 4.34): a second call site would be a second thing to
    forget on one of them.

    PRF004/005/007 are errors because each names a claim the corpus makes
    that its own mirrors do not back up (a declaration with no data, data
    with no declaration, an id naming nothing). PRF011, PRF012 and PRF014 are
    errors for a different reason -- no mirror is involved at all: a curated
    record contradicts itself, so there is no reading of it that is true. PRF009 is
    a warning because a gene missing from one source matrix is ordinary --
    the burden layer's nine registered genes absent from the Audain
    supplement are the precedent -- and the point is to name the gap, not to
    block the build over it.
    """
    # One pass over `profiles`: which tissues, stages and genes each dataset's
    # rows actually carry, and which shard (or its first-seen file) to blame.
    tissues_by_dataset: dict[str, set[str]] = defaultdict(set)
    stages_by_dataset: dict[str, set[str]] = defaultdict(set)
    shard_path_by_dataset: dict[str, Path] = {}
    profile_genes: set[str] = set()
    gene_path: dict[str, Path] = {}
    # (dataset, tissue, stage) -> genes with a profiles row there. Restricted
    # to rows with a non-null stage, for the same reason `validate_profiles`'s
    # `profile_triples` is: `profile_quantiles.stage` is never null, so a
    # null-stage row could never belong to an assayed cell in the first place.
    genes_by_triple: dict[_Triple, set[str]] = defaultdict(set)

    for path, row in _shard_rows(root, PROFILES, "profiles", _PROFILES_REFERENCE_COLUMNS):
        dataset, gene, tissue, stage = row["dataset"], row["gene"], row["tissue"], row["stage"]
        if dataset is None:
            # A null on a non-nullable column is TBL003's to report.
            continue
        shard_path_by_dataset.setdefault(dataset, path)
        if tissue is not None:
            tissues_by_dataset[dataset].add(tissue)
        if stage is not None:
            stages_by_dataset[dataset].add(stage)
        if gene is not None:
            profile_genes.add(gene)
            gene_path.setdefault(gene, path)
        if tissue is not None and stage is not None and gene is not None:
            genes_by_triple[(dataset, tissue, stage)].add(gene)

    # One pass over `profile_quantiles`: which (dataset, tissue, stage) cells
    # were actually assayed, and a shard path to blame for each.
    quantile_triples: set[_Triple] = set()
    quantile_triple_path: dict[_Triple, Path] = {}
    for path, row in _shard_rows(
        root, PROFILE_QUANTILES, "profile_quantiles", _QUANTILES_REFERENCE_COLUMNS
    ):
        dataset, tissue, stage = row["dataset"], row["tissue"], row["stage"]
        if dataset is None or tissue is None or stage is None:
            continue
        triple = (dataset, tissue, stage)
        quantile_triples.add(triple)
        quantile_triple_path.setdefault(triple, path)

    issues: list[ValidationIssue] = []
    issues.extend(_prf004_issues(root, datasets, tissues_by_dataset, shard_path_by_dataset))
    issues.extend(_prf005_issues(root, datasets, stages_by_dataset, shard_path_by_dataset))
    issues.extend(_prf011_issues(root, datasets))
    issues.extend(_prf012_issues(root, datasets))
    issues.extend(_prf014_issues(root, datasets))
    issues.extend(_prf006_issues(root, phases))
    issues.extend(_prf007_issues(profile_genes, gene_path, known_genes))
    issues.extend(
        _prf009_issues(quantile_triples, quantile_triple_path, genes_by_triple, published_genes)
    )
    return issues


def _shard_rows(
    root: Path, schema: TableSchema, schema_name: str, columns: tuple[str, ...]
) -> Iterator[tuple[Path, dict[str, Any]]]:
    """Every row of every shard for `schema_name`, restricted to `columns`.

    Yields `(path, row)` so a caller can blame the file a row came from. A
    shard that could not be read, or that is missing one of `columns`, is
    skipped rather than raising -- see `validate_profiles`'s docstring for why
    that is safe here.
    """
    for path, name in mirror_paths(root):
        if name != schema_name:
            continue
        frame, _ = read_table(path, schema)
        if frame is None or not set(columns) <= set(frame.columns):
            continue
        yield from ((path, row) for row in frame.select(list(columns)).iter_rows(named=True))


def _prf008_issues(
    units_by_dataset: dict[str, set[str]], path_by_dataset: dict[str, Path]
) -> list[ValidationIssue]:
    """PRF008 -- one dataset, more than one `profiles` unit.

    Sorted by dataset, and each message's unit list is sorted too: both are
    built from a `set`, whose iteration order is not insertion order and can
    follow `PYTHONHASHSEED`, so an unsorted report would let two validate runs
    over one identical repository disagree on the order issues arrive in, or
    on the order units are named within one message.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(units_by_dataset):
        units = units_by_dataset[dataset]
        if len(units) > 1:
            issues.append(
                ValidationIssue(
                    "PRF008",
                    Severity.ERROR,
                    str(path_by_dataset[dataset]),
                    f"dataset '{dataset}' reports profiles.median_abundance in more "
                    f"than one unit: {sorted(units)}",
                )
            )
    return issues


def _prf001_issues(
    profile_units: dict[str, set[str]],
    profile_unit_path: dict[str, Path],
    quantile_units: dict[str, set[str]],
    quantile_available: bool,
) -> list[ValidationIssue]:
    """PRF001 -- one dataset's `profiles` and `profile_quantiles` units must match.

    Skipped entirely when no `profile_quantiles` shard could be read at all:
    see the module docstring for why comparing against an unloaded table would
    cascade rather than inform. Independent of PRF008 above -- a dataset can
    trip both at once, and each names a different fact worth fixing on its
    own: PRF008 that `profiles` alone is inconsistent, this that whatever it
    reports does not fully match `profile_quantiles` either.
    """
    if not quantile_available:
        return []
    issues: list[ValidationIssue] = []
    for dataset in sorted(set(profile_units) & set(quantile_units)):
        p_units, q_units = profile_units[dataset], quantile_units[dataset]
        if p_units != q_units:
            issues.append(
                ValidationIssue(
                    "PRF001",
                    Severity.ERROR,
                    str(profile_unit_path[dataset]),
                    f"dataset '{dataset}' reports profiles in unit(s) {sorted(p_units)} "
                    f"but profile_quantiles in unit(s) {sorted(q_units)}",
                )
            )
    return issues


def _prf002_issues(
    profile_triples: set[_Triple],
    profile_triple_path: dict[_Triple, Path],
    quantile_triples: set[_Triple],
    quantile_available: bool,
) -> list[ValidationIssue]:
    """PRF002 -- a `(dataset, tissue, stage)` with `profiles` rows and no
    `profile_quantiles` rows.

    Skipped entirely when no `profile_quantiles` shard could be read at all,
    for the same reason as PRF001: every triple in `profiles` would otherwise
    report "missing", naming a hundred symptoms of one unread file rather
    than the file itself.

    `sorted(... - ...)` rather than iterating the set difference directly:
    the difference is a fresh `set`, and its iteration order can follow
    `PYTHONHASHSEED` rather than any property of the data.
    """
    if not quantile_available:
        return []
    issues: list[ValidationIssue] = []
    for triple in sorted(profile_triples - quantile_triples):
        dataset, tissue, stage = triple
        issues.append(
            ValidationIssue(
                "PRF002",
                Severity.ERROR,
                str(profile_triple_path[triple]),
                f"dataset '{dataset}' has profiles rows for tissue '{tissue}' stage "
                f"'{stage}', but profile_quantiles has none",
            )
        )
    return issues


def _prf003_issues(
    breakpoints: dict[_Triple, list[tuple[int, float]]], breakpoint_path: dict[_Triple, Path]
) -> list[ValidationIssue]:
    """PRF003 -- breakpoints must be non-decreasing along the percentile axis.

    Non-decreasing, not strictly increasing: roughly half a bulk library sits
    at zero, so long runs of tied breakpoints at the floor are the ordinary
    case, not a defect.

    One issue per offending triple, naming the first inversion found, not one
    per decreasing pair -- the same distinct-entity aggregation as PRF002,
    chosen so a triple with several inversions does not drown a report in
    near-duplicate lines.

    Sorted by `percentile` before comparison rather than trusting row order:
    nothing guarantees a shard arrives with its rows in percentile order by
    the time this function sees it (`validate_sort_order` is a separate
    check, and this function is exercised standalone in every test in this
    module), so trusting file order could read an unrelated transcription
    order as a real inversion, or miss one hidden by it.
    """
    issues: list[ValidationIssue] = []
    for triple in sorted(breakpoints):
        dataset, tissue, stage = triple
        ordered = sorted(breakpoints[triple], key=lambda item: item[0])
        # strict=False: `ordered[1:]` is one element shorter than `ordered` by
        # construction -- this is the standard pairwise-zip idiom, not a
        # length mismatch to guard against.
        for (prev_percentile, prev_value), (percentile, value) in zip(
            ordered, ordered[1:], strict=False
        ):
            if value < prev_value:
                issues.append(
                    ValidationIssue(
                        "PRF003",
                        Severity.ERROR,
                        str(breakpoint_path[triple]),
                        f"dataset '{dataset}' tissue '{tissue}' stage '{stage}': "
                        f"percentile {percentile} has value {value}, below percentile "
                        f"{prev_percentile}'s value {prev_value}",
                    )
                )
                break
    return issues


def _prf004_issues(
    root: Path,
    datasets: tuple[Dataset, ...],
    tissues_by_dataset: dict[str, set[str]],
    shard_path_by_dataset: dict[str, Path],
) -> list[ValidationIssue]:
    """PRF004 -- a declared `cardiac_tissues` token matched by no profiles row.

    Checked from the curated side, one dataset at a time: for each dataset
    this atlas actually has a record for, its own declared tokens are
    compared against the tissues its own profiles rows use. A `profiles` row
    naming a dataset with no curated record at all is out of scope for this
    check (and for PRF005 below) -- there is no declaration here to compare
    it against, and inventing one would be guessing rather than checking.

    A dataset with zero profiles rows at all still fires this, once, naming
    every declared token: a curated record with nothing behind it yet is a
    claim its own mirror has not backed up, not a tolerated mid-curation
    state. The location then falls back to the shard path the dataset's
    accession would use by convention (`mirrors/profiles/<id>.tsv`), even
    though the file may not exist -- still the file a curator would create.

    Sorted by dataset id explicitly, rather than trusting `datasets`' own
    tuple order (which follows directory-listing order in `corpus.py`, not
    id) -- and each message's missing-token list is sorted too, since it is
    built from a `set` difference whose unsorted iteration order can follow
    PYTHONHASHSEED. Pinned together by
    `test_prf004_and_prf005_are_reported_in_sorted_order_by_dataset_and_within_message`.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(datasets, key=lambda item: item.id):
        observed = tissues_by_dataset.get(dataset.id, set())
        missing = sorted(set(dataset.cardiac_tissues) - observed)
        if missing:
            path = shard_path_by_dataset.get(
                dataset.id, root / "mirrors" / "profiles" / f"{dataset.id}.tsv"
            )
            issues.append(
                ValidationIssue(
                    "PRF004",
                    Severity.ERROR,
                    str(path),
                    f"dataset '{dataset.id}' declares cardiac_tissues {missing}, "
                    f"matched by no profiles row for this dataset",
                )
            )
    return issues


def _prf005_issues(
    root: Path,
    datasets: tuple[Dataset, ...],
    stages_by_dataset: dict[str, set[str]],
    shard_path_by_dataset: dict[str, Path],
) -> list[ValidationIssue]:
    """PRF005 -- a `profiles.stage` token a dataset's own record does not declare.

    The reverse direction from PRF004, over the same curated datasets: PRF004
    asks whether a declaration is backed up by data, this asks whether the
    data is backed up by a declaration. A dataset with zero profiles rows
    contributes nothing here (`stages_by_dataset.get(..., set())` is empty),
    correctly -- there is no observed token to call undeclared.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(datasets, key=lambda item: item.id):
        declared = {stage.token for stage in dataset.stages}
        extra = sorted(stages_by_dataset.get(dataset.id, set()) - declared)
        if extra:
            path = shard_path_by_dataset.get(
                dataset.id, root / "mirrors" / "profiles" / f"{dataset.id}.tsv"
            )
            issues.append(
                ValidationIssue(
                    "PRF005",
                    Severity.ERROR,
                    str(path),
                    f"dataset '{dataset.id}' has profiles rows naming stage(s) {extra} "
                    f"that its own dataset record does not declare in 'stages'",
                )
            )
    return issues


def _dataset_location(root: Path, dataset: Dataset) -> str:
    """The curated file a stage-ordering issue belongs to.

    By convention rather than by lookup: `corpus.py` loads every `*.yaml`
    under `curation/datasets/` and does not record which file a record came
    from, so this names the file a curator would open. Same convention
    `_prf004_issues` uses for its shard-path fallback.
    """
    return str(root / "curation" / "datasets" / f"{dataset.id}.yaml")


def _prf011_issues(root: Path, datasets: tuple[Dataset, ...]) -> list[ValidationIssue]:
    """PRF011 -- two stages of one dataset claiming the same `order`.

    An error, not a warning: with a duplicate position the sort's tie-break
    (the token) decides the sequence, which is exactly the alphabetical
    ordering `Stage.order` exists to remove -- reappearing inside one
    dataset, silently, under a green build.

    **Numbered PRF011, not PRF010.** `validate/runner.py` already issues
    PRF010 for an unrelated failure (profiles rows present with no
    `profile_quantiles` data), where it is PRF000's partner error. A skip
    warning is only safe because it arrives with *the* error that caused it,
    so one code must not name two failures.

    Scanned in `(order, token)` order rather than in the curated list's own
    order: which of a colliding pair gets named must not depend on which
    record a curator happened to type first, and `Stage`'s own docstring says
    the list's order is deliberately not load-bearing. Datasets are sorted by
    id for the same reason `_prf004_issues` sorts them -- `corpus.datasets`
    follows directory-listing order, not id.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(datasets, key=lambda item: item.id):
        seen: dict[int, str] = {}
        for stage in sorted(dataset.stages, key=lambda item: (item.order, item.token)):
            first = seen.get(stage.order)
            if first is None:
                seen[stage.order] = stage.token
                continue
            issues.append(
                ValidationIssue(
                    "PRF011",
                    Severity.ERROR,
                    _dataset_location(root, dataset),
                    f"dataset '{dataset.id}' stages {first!r} and {stage.token!r} both "
                    f"declare order {stage.order}; a position must name one stage",
                )
            )
    return issues


def _prf012_issues(root: Path, datasets: tuple[Dataset, ...]) -> list[ValidationIssue]:
    """PRF012 -- `order` disagrees with `wpc` about which stage is earlier.

    This is the failure `Stage.order` introduces by existing: every record is
    individually valid -- an integer and a positive float -- so pydantic
    cannot see it, and a stage transcribed into the wrong slot would publish
    a developmental trajectory with two of its points swapped.

    Checked only across stages carrying **both** values. A post-natal stage's
    null `wpc` is a fact about the stage, not a missing value, and comparing
    against it would report every profile dataset this atlas will curate.

    **It cannot check the post-natal block's internal order at all** -- those
    eight stages have no `wpc`, which is precisely why they needed `order` --
    so `tests/test_repository_validates.py` pins that sequence against a
    literal instead. Do not read a clean PRF012 as "the chronology is
    checked".

    What it also could not check, until PRF014, is where that block sits
    relative to the prenatal one: see `_prf014_issues`, which enforces the
    one implication `Stage.wpc`'s own semantics already carry.

    The comparison is `>=`, not `>`: two tokens at one `wpc` are two names for
    one point in time, so `order` would be claiming a sequence its own
    evidence does not support, and `wpc` is the key `assign_phase` joins on,
    so both would resolve to identical phases while the axis drew them apart.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(datasets, key=lambda item: item.id):
        dated = sorted(
            (stage for stage in dataset.stages if stage.wpc is not None),
            key=lambda item: item.order,
        )
        # Adjacent pairs, so the two sequences differ in length by design --
        # `strict=False` is the intent, not a silenced B905.
        for earlier, later in zip(dated, dated[1:], strict=False):
            # The two `is not None` tests are already true of everything in
            # `dated`; they are here so mypy can narrow `float | None` to
            # `float`, not because a null can reach this line.
            if earlier.wpc is not None and later.wpc is not None and earlier.wpc >= later.wpc:
                issues.append(
                    ValidationIssue(
                        "PRF012",
                        Severity.ERROR,
                        _dataset_location(root, dataset),
                        f"dataset '{dataset.id}' orders stage {earlier.token!r} "
                        f"(order {earlier.order}, {earlier.wpc} wpc) before "
                        f"{later.token!r} (order {later.order}, {later.wpc} wpc)",
                    )
                )
    return issues


def _prf014_issues(root: Path, datasets: tuple[Dataset, ...]) -> list[ValidationIssue]:
    """PRF014 -- a stage with no `wpc` ordered before a stage that has one.

    The half of the chronology PRF012 cannot reach. PRF012 compares `order`
    against `wpc` only across stages carrying **both**, and every post-natal
    stage has `wpc: null` by definition -- so nothing constrained where the
    post-natal block sorted relative to the prenatal one, and the two checks
    together left a hole exactly the shape of the defect `Stage.order` was
    added to fix.

    Measured 2026-08-20 on the committed corpus: renumbering E-MTAB-6814's
    eight post-natal tokens to `order` 1-8 and its thirteen prenatal tokens
    to 9-21 produced **0 errors, 4 warnings** -- codes byte-identical to the
    clean baseline -- and TBX5's chart published "Median abundance in whole
    heart, neonate to 19 week post conception", septation bands on the right
    half of the axis. The only thing standing against it was one literal in
    `tests/test_repository_validates.py` hardcoded to `id == "E-MTAB-6814"`
    and the slice `ordered[13:]`; a second profile dataset inherited no pin
    at all, and one with an inverted block also validated 0/4.

    This asserts nothing new about the data -- it enforces what the model
    already says. `Stage.wpc`'s docstring: a null is post-natal, "which is
    not a gap: a post-natal stage is outside every cardiac morphogenetic
    phase by definition". Post-natal is after prenatal, so a null-`wpc` stage
    ordering before a dated one contradicts the record's own semantics. An
    ERROR for PRF011/PRF012's reason: no mirror is involved, a curated record
    contradicts itself, and there is no reading of it that is true.

    One issue per offending post-natal stage rather than one per dataset, and
    each names the **latest** dated stage it claims to precede: a curator
    fixing this needs to know which token moved. Sorted by `order` within a
    dataset and by dataset id across them, both comparison-based over tuples
    and lists, so no PYTHONHASHSEED dependence is in play (section 4.42).

    A dataset with no dated stage at all cannot violate this, and neither can
    one with no post-natal stage -- both fall out of `max`/comparison rather
    than needing a guard, and `test_prf014_is_silent_when_every_post_natal_
    stage_orders_last` pins the ordinary shape.
    """
    issues: list[ValidationIssue] = []
    for dataset in sorted(datasets, key=lambda item: item.id):
        dated = [stage for stage in dataset.stages if stage.wpc is not None]
        if not dated:
            continue
        latest = max(dated, key=lambda item: item.order)
        undated = sorted(
            (stage for stage in dataset.stages if stage.wpc is None),
            key=lambda item: item.order,
        )
        for stage in undated:
            if stage.order < latest.order:
                issues.append(
                    ValidationIssue(
                        "PRF014",
                        Severity.ERROR,
                        _dataset_location(root, dataset),
                        f"dataset '{dataset.id}' orders stage {stage.token!r} "
                        f"(order {stage.order}, no wpc) before {latest.token!r} "
                        f"(order {latest.order}, {latest.wpc} wpc); a stage with no "
                        "wpc is post-natal by construction and must order after "
                        "every stage that has one",
                    )
                )
    return issues


def _coverage_spans(phases: Sequence[CardiacPhase]) -> list[_CoverageSpan]:
    """Merge overlapping or touching phase intervals into maximal coverage spans.

    **Why a merge, not a pairwise-adjacent comparison.** Before phases could
    overlap, "sort by start, compare each phase only to the very next one"
    was already a correct interior-gap test, because no phase could ever
    reach past its immediate successor. Once phases may overlap
    (`models/phases.py`'s module docstring), that stops being true: a long
    phase can span past several shorter ones nested or overlapping inside it,
    and comparing only *adjacent-in-sort-order* phases would then report a
    false gap between two short phases that a third, wider phase already
    covers. Worked example: `A=[1,10)`, `B=[2,3)`, `C=[8,9)` sorts as
    `A, B, C`; comparing `B` to `C` directly sees `C.start (8) > B.end (3)`
    and reports a gap that does not exist, because `A` already covers all of
    `[1,10)`, `B` and `C` included. Merging into spans first is what keeps
    the check correct once overlap is legal, not merely unchanged for the
    non-overlapping case it used to be the only case.

    A phase merges into the running span when its own start falls at or
    before the span's current end (`<=`, matching the half-open-interval
    "touching is not a gap" rule the old pairwise check already used); the
    span's end then becomes the *later* of the two ends, never simply
    replaced, so a short phase nested inside a longer one cannot shrink the
    span back down. Each span's `phase_ids` names every phase merged into it,
    sorted by `(start_wpc, id)` via the traversal order -- `_prf006_issues`
    re-sorts them again before rendering, so this function does not need to
    guarantee that ordering on its own return.

    **A phase with no stated end (`end_wpc is None`, `EndBasis.NOT_STATED`)
    contributes no coverage at all** -- skipped before it can be merged, the
    same exclusion `CardiacPhaseFile.phases_for` applies and for the same
    reason: there is no source-given point at which such a phase's coverage
    would end, so it must not be treated as closing a gap it was never shown
    to close. This is precisely why the real vocabulary now has a genuine
    interior gap where an earlier, invented end for `heart_looping` used to
    paper over one -- see `curation/cardiac_phases.yaml`'s own header comment.
    """
    ordered = sorted(phases, key=lambda phase: (phase.start_wpc, phase.id))
    spans: list[_CoverageSpan] = []
    for phase in ordered:
        if phase.end_wpc is None:
            continue
        if spans and phase.start_wpc <= spans[-1][1]:
            start, end, ids = spans[-1]
            spans[-1] = (start, max(end, phase.end_wpc), (*ids, phase.id))
        else:
            spans.append((phase.start_wpc, phase.end_wpc, (phase.id,)))
    return spans


def _prf006_issues(root: Path, phases: CardiacPhaseFile | None) -> list[ValidationIssue]:
    """PRF006 -- an interior gap no declared cardiac phase covers.

    WARNING, not ERROR, and interior gaps only -- the spec's single-line code
    table hid that this splits in two.

    **Redefined for overlap.** Phases may now overlap by design
    (`models/phases.py`'s module docstring), so "gap" can no longer mean "the
    space between one phase and the very next in sorted order" -- that
    pairwise reading misfires once a wider phase can span past a nested,
    narrower one; see `_coverage_spans` for a worked counter-example. A gap is
    instead a wpc region strictly between the earliest start and the latest
    end that **no phase, individually, reaches** -- computed by merging every
    phase into maximal coverage spans first (`_coverage_spans`) and reporting
    only between two *spans*, never between two raw phases.

    Never the region before the first phase or after the last: a curated
    window legitimately stops before the post-natal stages (a null-wpc stage
    is post-natal *by construction*, per `Stage.wpc`), so firing there would
    trigger on every correct vocabulary -- and a check that fires on every
    correct input is a check a curator learns to ignore. `phases=None` (no
    vocabulary curated yet) skips this entirely, the same treatment
    `validate_profiles` gives an unread `profile_quantiles` mirror; an empty
    `phases.phases` list needs no separate guard, because `_coverage_spans`
    of an empty sequence is `[]` and the pairwise loop below is then a no-op.

    Each message names every phase id bordering the gap on each side, not
    just one, because a merged span can carry several overlapping phases --
    naming only one would let a curator "fix" the wrong phase. Pinned by
    `test_prf006_reports_one_issue_per_interior_gap_in_wpc_order` (the
    non-overlapping case, phases declared out of wpc order so a dropped sort
    compares the wrong pairs rather than merely reordering the report) and by
    `test_prf006_does_not_report_a_gap_a_wider_overlapping_phase_already_
    covers` (the overlap case `_coverage_spans`'s own docstring works through).
    """
    if phases is None:
        return []
    issues: list[ValidationIssue] = []
    location = str(root / "curation" / "cardiac_phases.yaml")
    spans = _coverage_spans(phases.phases)
    # strict=False: `spans[1:]` is one element shorter than `spans` by
    # construction, the standard pairwise-zip idiom -- not a length mismatch
    # to guard against. Consecutive spans only, which is what makes this
    # interior: the region before `spans[0]` and after `spans[-1]` never
    # appears as either half of a pair.
    for earlier, later in zip(spans, spans[1:], strict=False):
        _, earlier_end, earlier_ids = earlier
        later_start, _, later_ids = later
        if later_start > earlier_end:
            issues.append(
                ValidationIssue(
                    "PRF006",
                    Severity.WARNING,
                    location,
                    f"a gap [{earlier_end}, {later_start}) falls between phase(s) "
                    f"{sorted(earlier_ids)} (covering up to {earlier_end}) and phase(s) "
                    f"{sorted(later_ids)} (from {later_start}); a stage in this range "
                    f"falls in no phase",
                )
            )
    return issues


def _prf007_issues(
    profile_genes: set[str], gene_path: dict[str, Path], known_genes: set[str] | None
) -> list[ValidationIssue]:
    """PRF007 -- a `profiles.gene` id in no gene registry.

    `known_genes=None` means the registry did not load, and this whole check
    is skipped rather than reporting each row as dangling -- the cascade
    REF000/SRC000/ONT000 exist to prevent. Deliberately `is not None` rather
    than a truthy check: an empty-but-loaded registry (a header-only
    `mirrors/genes.tsv`) must still flag every gene as unknown, and an empty
    `set` is falsy in Python.

    One issue per distinct gene id, matching `validate_burden_references`'s
    BUR011 -- a shard with one bad id yields one issue, not one per row that
    cites it.
    """
    if known_genes is None:
        return []
    issues: list[ValidationIssue] = []
    for gene in sorted(profile_genes - known_genes):
        issues.append(
            ValidationIssue(
                "PRF007",
                Severity.ERROR,
                str(gene_path[gene]),
                f"gene {gene} is not in mirrors/genes.tsv",
            )
        )
    return issues


def _prf009_issues(
    quantile_triples: set[_Triple],
    quantile_triple_path: dict[_Triple, Path],
    genes_by_triple: dict[_Triple, set[str]],
    published_genes: set[str],
) -> list[ValidationIssue]:
    """PRF009 -- a published gene missing a profiles row in a cell that was assayed.

    Inverts PRF002's question. `profile_quantiles` records which
    `(dataset, tissue, stage)` cells were actually assayed; for each one, a
    published gene with no `profiles` row there was *dropped*, not merely
    unsampled -- a distinction `tau` cannot make, since it only ever sees the
    rows that exist. WARNING, not ERROR: a gene absent from a source matrix is
    ordinary on its own (the burden layer's nine registered genes missing
    from the Audain supplement are the precedent), and the point is to name
    the gap, not to block the build over it.

    Checked against `published_genes`, never `known_genes`: the two are
    different sets and this function must not conflate them. A merely
    registered gene missing from one cell is unremarkable -- 67 of the
    186-gene registry are not published at all (2026-11-14) -- while a
    *published* one is a
    gap on a page a reader can already reach.

    **Grouped by assayed cell, not by gene or by (gene, cell) pair.** The
    realistic failure this guards against -- a converter that drops a
    scattered subset of rows within one tissue/stage slice, such as dropped
    NA rows or a bad range in a copy-paste -- leaves every dropped gene
    sharing the *same* cell, so one issue per cell names the whole defect
    once, listing every gene missing there. Grouping by gene instead would
    split that single root cause into as many issues as genes it happened to
    drop; grouping by (gene, cell) pair would size the report at the product
    of published genes and assayed cells (92 x roughly 150 today -- up to
    13,800), which is the spam this choice exists to avoid. Bounded instead
    by the number of assayed cells alone. Pinned by
    `test_prf009_groups_by_assayed_cell_not_by_gene_or_by_row`.

    No `quantile_available`-style guard is needed here, unlike PRF001/PRF002
    in `validate_profiles`. Those subtract *from* the profiles side, so
    treating an unread quantile mirror as merely empty would make that
    subtraction return everything in `profiles` -- the false cascade the flag
    exists to stop. This function instead iterates `quantile_triples`
    *directly*: when no `profile_quantiles` shard could be read, that set is
    empty by construction, and the loop below is then simply a no-op. Pinned
    by `test_prf009_reports_nothing_when_no_quantile_shard_can_be_read`
    rather than left as reasoning nobody checked.

    `quantile_triples` is itself a `set` of 3-tuples, so its unsorted
    iteration order can follow PYTHONHASHSEED exactly like PRF002's
    `profile_triples - quantile_triples` in `validate_profiles` -- wrapping a
    value in a tuple does not remove that risk, it only changes the specific
    coincidence rate (measured: a five-word tissue set that never coincided
    with sorted order as bare strings measured a real, nonzero rate once
    wrapped as `(dataset, tissue, stage)`). `sorted(quantile_triples)` fixes
    that regardless of table layout, comparing the full triple
    lexicographically. Each message's missing-gene list is sorted for the
    same reason PRF004's is. See
    `test_prf009_issues_and_their_message_lists_are_sorted` for the measured
    rates behind the fixture size chosen there.
    """
    issues: list[ValidationIssue] = []
    for triple in sorted(quantile_triples):
        missing = sorted(published_genes - genes_by_triple.get(triple, set()))
        if missing:
            dataset, tissue, stage = triple
            issues.append(
                ValidationIssue(
                    "PRF009",
                    Severity.WARNING,
                    str(quantile_triple_path[triple]),
                    f"dataset '{dataset}' tissue '{tissue}' stage '{stage}' was assayed "
                    f"(profile_quantiles has data for it), but {len(missing)} published "
                    f"gene(s) have no profiles row here: {missing}",
                )
            )
    return issues
