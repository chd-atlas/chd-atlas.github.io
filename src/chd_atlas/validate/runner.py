# src/chd_atlas/validate/runner.py
"""Runs every validator and aggregates the result into one deterministic report."""

from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from pathlib import Path

from chd_atlas.build.identity import read_gene_identities
from chd_atlas.build.validity import gene_validity, published_genes
from chd_atlas.corpus import Corpus, load_curation, unexpected_curation_entries
from chd_atlas.genes import GeneRegistry
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.dataset import Dataset
from chd_atlas.orphanet_scope import load_orphanet_tables, pinned_structural_phenotypes
from chd_atlas.tables import (
    CHD_CASE_SPLITS,
    PROFILE_QUANTILES,
    PROFILES,
    TABLE_SCHEMAS,
    TableSchema,
    mirror_paths,
    read_table,
    unexpected_mirror_entries,
    validate_table,
)
from chd_atlas.validate.burden import (
    validate_burden,
    validate_burden_combined,
    validate_burden_combined_references,
    validate_burden_references,
)
from chd_atlas.validate.cell_types import validate_cell_type_profiles
from chd_atlas.validate.genes import validate_identity_ids, validate_mirror_symbols
from chd_atlas.validate.ids import load_id_registry, validate_ids
from chd_atlas.validate.ontology import (
    OntologyRegistry,
    validate_labels,
    validate_terms,
)
from chd_atlas.validate.profiles import validate_profile_references, validate_profiles
from chd_atlas.validate.referential import (
    validate_mirror_references,
    validate_ptm_evidence_is_reachable,
    validate_references,
)
from chd_atlas.validate.scope import (
    orphanet_candidates,
    scope_candidates,
    validate_curation_is_in_scope,
    validate_orphanet_warrants,
    validate_scope_attribution,
    validate_scope_terms,
)
from chd_atlas.validate.sort_order import validate_sort_order
from chd_atlas.validate.sources import load_sources, validate_source_references


@dataclass
class ValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    def __post_init__(self) -> None:
        # Sorted by location first so everything wrong with one file reads together.
        # ValidationIssue's own ordering is code-first, which groups by error type and
        # scatters a single file's problems across the report — fine as a stable
        # default, wrong for a human fixing one file at a time.
        self.issues = sorted(
            set(self.issues),
            key=lambda issue: (issue.location, issue.code, issue.severity, issue.message),
        )

    @property
    def error_count(self) -> int:
        return sum(1 for issue in self.issues if issue.severity is Severity.ERROR)

    @property
    def warning_count(self) -> int:
        return sum(1 for issue in self.issues if issue.severity is Severity.WARNING)

    @property
    def ok(self) -> bool:
        # Warnings are deliberately ignored, so a warnings-only report exits 0.
        # Two kinds of warning exist, and the reasoning differs for each.
        #
        # SKIP warnings (REF000, SRC000, ONT000) mean "a check did not run".
        # Ignoring one would let CI go green on a repository whose checks were
        # silently skipped, so each is emitted in the same branch as the error
        # that caused the skip — a failed corpus load, or a failed source
        # registry load — and a skip can therefore never be the only thing
        # wrong. A skip warning that could arrive alone must be raised to an
        # error. Pinned by
        # test_every_skip_warning_arrives_with_the_error_that_caused_it.
        #
        # GAP warnings (REF013) mean "a check ran and found curated data that
        # will not reach a reader". They arrive alone, deliberately: the record
        # is sound and the site is publishable, and what is missing is one
        # mirror cell. Blocking a whole deploy on that would be the wrong trade
        # for a resource curated incrementally. This is a real weakening — the
        # site can publish with a known gap — and it is accepted because the
        # alternative before REF013 existed was the same gap reported nowhere at
        # all. The safety they rest on is different from the skip warnings': the
        # check ran, so the report names the gene, the assertion and the file.
        # Pinned by test_a_gap_warning_is_reported_without_blocking_the_build.
        return self.error_count == 0

    def render(self) -> str:
        lines = [issue.format() for issue in self.issues]
        # Blank line so the summary reads as a summary rather than as one more
        # issue running straight on from the last.
        lines.append("")
        lines.append(f"{self.error_count} error(s), {self.warning_count} warning(s)")
        return "\n".join(lines)


def _used_sources(root: Path) -> dict[str, set[str]]:
    """Collect the distinct ``source`` values each mirror table relies on.

    The burden mirror relies on one more for each study in `CHD_CASE_SPLITS`
    whose rows carry `n_chd_case_carriers`: the file that column is transcribed
    from, which no row's `source` names (owner decision R2). Collected here so
    it is held to the registry like every other, and the table a gene page's
    card cites is one whose licence is recorded and permits redistribution.
    """
    used: dict[str, set[str]] = {}
    for path, schema_name in mirror_paths(root):
        schema = TABLE_SCHEMAS[schema_name]
        if not any(column.name == "source" for column in schema.columns):
            continue
        frame, _ = read_table(path, schema)
        if frame is None or "source" not in frame.columns:
            continue
        values = {str(value) for value in frame["source"].to_list() if value is not None}
        if {"study", "n_chd_case_carriers"} <= set(frame.columns):
            values |= {
                CHD_CASE_SPLITS[study].source
                for study, chd in zip(
                    frame["study"].to_list(), frame["n_chd_case_carriers"].to_list(), strict=True
                )
                if chd is not None and study in CHD_CASE_SPLITS
            }
        if values:
            used[str(path)] = values
    return used


def _known_genes(root: Path) -> set[str] | None:
    """The gene registry, or None when it could not be read.

    None is deliberately distinct from an empty set. An empty registry makes
    every assertion and every functional record report REF001, burying the one
    real problem — the registry itself — under one error per record, with
    nothing anywhere naming the cause.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return None
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns:
        return None
    return {value for value in frame["hgnc_id"].to_list() if value is not None}


def _symbol_registry(root: Path) -> GeneRegistry | None:
    """`mirrors/genes.tsv` as a symbol lookup, or None when it cannot be read.

    None rather than an empty registry, for the reason `_gene_registry` above
    gives: an empty one makes every symbol in the corpus report GEN001, burying
    the single real problem under one warning per row.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return None
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns or "symbol" not in frame.columns:
        return None
    return GeneRegistry(frame.iter_rows(named=True))


def _mirror_gene_labels(root: Path) -> list[tuple[str, str, str]]:
    """`(hgnc_id, symbol, where)` for every distinct pair the validity mirrors print.

    Deduplicated by the caller's `set`, not here: ClinGen carries one row per
    gene-disease pair, so a gene with four curated diseases would otherwise
    report the same naming disagreement four times.
    """
    labelled: list[tuple[str, str, str]] = []
    for name, schema in (
        ("clingen_gene_validity.tsv", "clingen_validity"),
        ("gencc_submissions.tsv", "gencc_submissions"),
    ):
        path = root / "mirrors" / name
        if not path.is_file():
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema])
        if frame is None or "gene" not in frame.columns or "gene_symbol" not in frame.columns:
            continue
        for gene, symbol in zip(
            frame["gene"].to_list(), frame["gene_symbol"].to_list(), strict=True
        ):
            if gene is not None and symbol is not None:
                labelled.append((str(gene), str(symbol), str(path)))
    return labelled


def _mirrored_validity(
    root: Path,
) -> (
    tuple[
        dict[str, str],
        dict[str, set[str]],
        dict[str, set[str]],
        dict[str, frozenset[str]],
        dict[str, frozenset[str]],
    ]
    | None
):
    """Disease labels and gene/disease cross-links, from every validity mirror.

    Returns three mappings the scope validators need, built once from both
    mirror tables so each validator stays a pure check over already-loaded
    data:

    - a MONDO term -> the label its authority (ClinGen or GenCC) published
    - a MONDO term -> every HGNC id curated under it
    - an HGNC id -> every MONDO term it is curated under
    - a MONDO term -> every ClinGen expert panel that curated a record naming it
    - a MONDO term -> every GenCC submitter that did

    The last two are kept APART rather than merged into one set of names, and
    that separation is the whole of SCP005's strength. Merged, the check can
    only ask "did some authority use this term", so a curator may take a term
    whose sole warrant is a commercial laboratory's GenCC submission, label it
    `admitted_by: clingen_chd_panel`, and pass -- publishing on 24 pages that
    ClinGen's chartered panel scoped a disease it never curated. Found by
    adversarial review of #30, on the guard added in the same commit.

    None is deliberately distinct from an empty result, for exactly the reason
    `_known_genes` returns None: checking every scope term against an empty
    mapping would report SCP001 once per term — 68 failures naming the
    symptom, none naming the cause — so a mirror that could not be read at all
    is reported once, as SCP000, and nothing else runs. Unlike `_known_genes`,
    which needs one file, this needs *neither* mirror to be readable before
    giving up — one working mirror still gives real, if partial, coverage.
    """
    labels: dict[str, str] = {}
    genes_by_disease: dict[str, set[str]] = {}
    diseases_by_gene: dict[str, set[str]] = {}
    clingen_panels: dict[str, set[str]] = {}
    gencc_submitters: dict[str, set[str]] = {}
    read_any = False
    # GenCC first, so ClinGen's label overwrites it below: ClinGen is the
    # primary source for gene-disease validity and GenCC harmonises rather
    # than authors, so preferring GenCC's label would attribute it to the
    # resource that did not write it.
    for schema_name, filename in (
        ("gencc_submissions", "gencc_submissions.tsv"),
        ("clingen_validity", "clingen_gene_validity.tsv"),
    ):
        path = root / "mirrors" / filename
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        # All three columns this loop selects, not only `disease`: an upstream
        # rename can drop `disease_label` or `gene` just as easily, and
        # `frame.select` raises `ColumnNotFoundError` on any one of them
        # missing — `validate_table` (called separately, over the same file)
        # already reports that as TBL001, so nothing is lost by skipping the
        # mirror here rather than crashing the whole validation run on it.
        if frame is None or not {"disease", "disease_label", "gene"}.issubset(frame.columns):
            continue
        read_any = True
        # The column naming the body that made each call: ClinGen records carry
        # the expert panel, GenCC rows the submitter. Selected separately from
        # the three above because either mirror may lack it without the other
        # checks being affected -- `validate_table` reports a missing column as
        # TBL001 over the same file, so skipping it here loses nothing.
        clingen = schema_name == "clingen_validity"
        who = "gcep" if clingen else "submitter"
        into = clingen_panels if clingen else gencc_submitters
        curators = frame[who].to_list() if who in frame.columns else [None] * frame.height
        rows = frame.select(["disease", "disease_label", "gene"]).iter_rows()
        for (disease, label, gene), curator in zip(rows, curators, strict=True):
            if disease is None:
                continue
            if label is not None:
                labels[str(disease)] = str(label)
            if gene is not None:
                genes_by_disease.setdefault(str(disease), set()).add(str(gene))
                diseases_by_gene.setdefault(str(gene), set()).add(str(disease))
            if curator is not None:
                into.setdefault(str(disease), set()).add(str(curator))
    if not read_any:
        return None
    return (
        labels,
        genes_by_disease,
        diseases_by_gene,
        {term: frozenset(names) for term, names in clingen_panels.items()},
        {term: frozenset(names) for term, names in gencc_submitters.items()},
    )


def _mirror_has_rows(root: Path, schema_name: str, schema: TableSchema) -> bool:
    """True when some shard registered under `schema_name` has a readable row.

    Shared by the PRF010/PRF000 check below for both `profiles` and
    `profile_quantiles`. Three inputs all read as "no rows" here, deliberately
    conflated the way `validate_profiles`'s own `quantile_available` flag
    conflates them (see that module's docstring): a directory that was never
    created (`mirror_paths` yields a shard only when `path.is_file()`, per
    `tables.py`), a shard that exists but fails to parse (`read_table` returns
    `None` — already reported once, as TBL000, by the unconditional
    `mirror_paths` loop in `validate_repository`), and a shard that parses
    cleanly but is header-only (`frame.height == 0`). All three leave nothing
    for a `profiles` row's percentile to be read against, which is the one
    fact PRF010 exists to name — the finer question of *why* is TBL000/TBL001's
    to answer, not this function's.
    """
    for path, name in mirror_paths(root):
        if name != schema_name:
            continue
        frame, _ = read_table(path, schema)
        if frame is not None and frame.height > 0:
            return True
    return False


def _phase_matchable_stages(root: Path, datasets: tuple[Dataset, ...]) -> set[tuple[str, str]]:
    """Every `(dataset, stage)` in `profiles` that a cardiac phase could match.

    "Could match" is narrower than "is a profiles row", and the difference is
    section 4.41's whole point. A stage with a null `wpc` is post-natal *by
    construction* (`Stage.wpc`) and falls outside every cardiac morphogenetic
    phase whether or not a vocabulary is curated; a stage token the dataset's
    own record never declares has no `wpc` to compare at all and is PRF005's
    to report. Neither is work a phase vocabulary would have done, so neither
    may trigger the skip below -- a lone PRF000 in an otherwise clean report
    is exactly the failure `test_every_skip_warning_arrives_with_the_error_
    that_caused_it` exists to catch, and `GEN000` broke this on its first run
    by asking only whether the registry was absent.

    Returns the pairs rather than a bool so PRF013's message can say how much
    is at stake; the count is what tells a curator whether they are looking at
    a whole layer or one stray row.
    """
    wpc_by_dataset: dict[str, set[str]] = {
        dataset.id: {stage.token for stage in dataset.stages if stage.wpc is not None}
        for dataset in datasets
    }
    matchable: set[tuple[str, str]] = set()
    for path, name in mirror_paths(root):
        if name != "profiles":
            continue
        frame, _ = read_table(path, PROFILES)
        if frame is None:
            continue
        for dataset, stage in zip(
            frame["dataset"].to_list(), frame["stage"].to_list(), strict=True
        ):
            if dataset is None or stage is None:
                continue
            if stage in wpc_by_dataset.get(dataset, frozenset()):
                matchable.add((dataset, stage))
    return matchable


def _gate_published_genes(root: Path, corpus: Corpus) -> set[str]:
    """The D21 population `PRF009` checks `profiles` rows against.

    `validate_profile_references` needs the *published* gene set, never the
    whole registry `_known_genes` returns — profiles.py's own docstring is
    explicit the two must not be conflated, since most of the registry is not
    published at all. The only correct source for that population is
    `build.validity.published_genes`, the same function `build_site` calls to
    decide what to publish. Reimplementing the gate's rules a second time here
    (which submitters count, the agreement floor, ClinGen's veto) would risk
    the two copies drifting apart the next time the gate widens — the same
    class of defect CLAUDE.md section 9 keeps `headline_confidence` and
    `has_conflicting_evidence` written together in one place to prevent.

    `validate/burden.py`'s BUR009-011 docstring declined to add the equivalent
    check for burden rows, because deriving the published set inside
    `validate/` "needs build.validity.published_genes and inverts the
    layering: the validator would depend on the builder it exists to gate."
    That general concern does not reproduce as an import cycle here:
    `chd_atlas.build.validity` imports only `chd_atlas.build.emit` and
    `chd_atlas.vocab`, and neither imports anything under `chd_atlas.validate`
    — `build.runner` already imports both `validate.runner` and
    `build.validity` today, in that order, with no cycle, and this module
    still imports cleanly on its own with the new import added. The choice
    actually being made here is between reusing the one correct
    implementation and duplicating a classification this involved; for a
    project this dependent on a single source of truth, duplication is the
    larger risk.

    **Never lets `gene_validity`/`published_genes` raise out of this
    function.** They are written for post-gate data — `build_site` calls them
    only after `validate_repository` has already reported 0 errors — and they
    are not defensive the way this module's own mirror readers are.
    Reproduced directly against synthetic frames: a column an upstream rename
    dropped (the same shape Copilot found in `_mirrored_validity`) raises
    `KeyError`, and a classification term the atlas vocabulary has not mapped
    yet raises `ValueError` by design (`gene_validity`'s own docstring: "a
    term either mirror adds in a later release must reach a human"). Both are
    the right behaviour for `build_site`, which must refuse rather than
    publish a confidence nobody vetted, and neither is acceptable here: this
    call computes one WARNING-level check's population, not a
    correctness-critical result, so degrading to "nothing published, this
    attempt" is the safe failure, and a raised traceback out of `chd-atlas
    validate` is not. `build_site`'s own call is untouched by this and still
    raises exactly as before.

    Returns an empty set, never `None`, when either validity mirror cannot be
    read or the computation above fails. `published_genes` has no optional
    form the way `known_genes` does — it must always be a `set[str]` — so an
    empty set is the fail-safe equivalent of the `None`-means-skip convention
    used elsewhere in this module: PRF009 only ever subtracts published genes
    *from* what a cell already has rows for, so an empty left-hand side can
    under-report a gap but can never invent a false one.

    **And the degradation cannot happen silently, which is what makes it
    acceptable without a skip code of its own.** Every way this returns empty
    is already an ERROR somewhere else in the same report, so the rule that a
    skip must arrive with the error that caused it holds structurally rather
    than by a paired warning. Measured 2026-08-15:

    - an unreadable mirror is TBL000, and `_mirrored_validity` pairs TBL012
      with it;
    - a dropped column is TBL001;
    - an unmapped classification term — the one case where the mirrors are
      readable *and* schema-valid — cannot reach `gene_validity` at all,
      because both `classification` columns carry `allowed=` sets that
      `tests/unit/test_vocab.py::test_every_mirrored_classification_maps_and_
      an_unknown_one_raises` pins equal to the vocab maps this function's
      callee dispatches on. Measured: 7/7 and 8/8 terms, no drift either way.

    **The catch is narrow for the same reason.** Only `KeyError` and
    `ValueError` are absorbed -- the two data failures above. A programmer
    error has no accompanying validation error to make it visible, so it must
    propagate rather than leave PRF009 silently dead in a green report.

    That last one is the load-bearing link and the only one that is not
    self-evident from this file. If that test is ever deleted, the two
    hand-written literals can drift, a term valid to the schema but unmapped
    by the vocab reaches `gene_validity`, this `except` swallows it, and
    PRF009 silently checks nothing in an otherwise green report.
    """
    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    if clingen is None or gencc is None:
        return set()
    scope_terms = {str(entry.id) for entry in corpus.chd_scope}
    try:
        return published_genes(gene_validity(clingen, gencc, in_scope=scope_terms))
    except (KeyError, ValueError):
        # Narrow deliberately. These two are the *data* failure modes named
        # above -- a column an upstream rename dropped raises `KeyError`, and a
        # classification term the atlas vocabulary has not mapped raises
        # `ValueError` by design -- and each already arrives with its own ERROR
        # from `validate_table`, so degrading here cannot hide anything.
        #
        # A blanket `except Exception` also swallowed programmer errors, and
        # those have no accompanying error to make them visible. Measured
        # 2026-08-19 on this branch: injecting a `TypeError` into
        # `gene_validity` -- the shape of an ordinary refactor regression --
        # produced `0 error(s), 4 warning(s)` with `ok=True` and PRF009
        # silently checking nothing. The docstring's guarantee covered data
        # failures and not code ones. Raised by review on PR #39.
        return set()


def _relative_to_root(issues: list[ValidationIssue], root: Path) -> list[ValidationIssue]:
    """Rewrite absolute locations as repo-relative.

    Locations are embedded verbatim in the rendered report, so absolute paths
    would make output differ between a curator's machine and CI, and would sort
    a hardcoded relative literal away from the same file's other issues.

    The separator comes from `os.sep` rather than a literal `/` because every
    location is built as `str(path)`, which renders `C:\\repo\\mirrors\\genes.tsv`
    on Windows. A hardcoded `/` would match no prefix there, leave every location
    absolute, and silently void the guarantee this function exists to provide.
    """
    prefix = f"{root}{os.sep}"
    return [replace(issue, location=issue.location.removeprefix(prefix)) for issue in issues]


def validate_repository(root: Path) -> ValidationReport:
    issues: list[ValidationIssue] = []

    corpus, corpus_issues = load_curation(root)
    issues.extend(corpus_issues)

    registry, source_issues = load_sources(root)
    issues.extend(source_issues)

    ontologies = OntologyRegistry.from_files(root, registry.ontology_files())
    issues.extend(ontologies.load_issues)

    for path, schema_name in mirror_paths(root):
        schema = TABLE_SCHEMAS[schema_name]
        issues.extend(validate_table(path, schema))
        issues.extend(validate_sort_order(path, schema))
    # Unconditional: the loop above iterates what the schemas expect, so an entry
    # nothing claims — a shard directory lost to a typo — is invisible to it.
    issues.extend(unexpected_mirror_entries(root))
    # Also unconditional, and deliberately outside the corpus branch below: the
    # burden table's cross-column rules read nothing but the table itself, so a
    # corpus that failed to load must not silently take these checks with it.
    issues.extend(validate_burden(root))
    # Unconditional for the same reason: it reads the two burden mirrors and
    # nothing else.
    issues.extend(validate_burden_combined(root))
    # Same reasoning again: `profiles`/`profile_quantiles` internal consistency
    # reads nothing but the two mirrors themselves.
    issues.extend(validate_profiles(root))
    # PRF010/PRF000 -- the skip machinery `validate_profiles` deliberately does
    # not provide itself (see its module docstring): comparing `profiles`
    # against a `profile_quantiles` mirror that is missing or entirely
    # unreadable would otherwise cascade the way an unread `genes.tsv` or
    # validity mirror does elsewhere in this function -- every
    # (dataset, tissue, stage) in `profiles` would report PRF002 and every
    # dataset's unit would report PRF001, naming a hundred symptoms of one gap
    # rather than the gap itself.
    #
    # `mirror_paths` yields a table only when `path.is_file()` (`tables.py`),
    # so an ABSENT `mirrors/profile_quantiles/` directory raises no TBL000 at
    # all -- and that absence is the *ordinary* state mid-curation (profiles
    # committed, quantiles not yet), not a rare failure. Without a partner
    # error here, PRF000 would be a lone warning; `ValidationReport.ok` ignores
    # warnings, so the build would publish every gene with its percentile band
    # silently missing. Precedents: the hand-written TBL008 below (missing
    # gene registry) and TBL012 (missing validity mirrors).
    #
    # Conditioned on `profiles_have_rows`, per section 4.41: emit a skip only
    # when there was work to skip. A repository with no `profiles` data at all
    # -- true of the committed corpus today -- must fire neither code, exactly
    # like GEN000's own guard.
    profiles_have_rows = _mirror_has_rows(root, "profiles", PROFILES)
    if profiles_have_rows and not _mirror_has_rows(root, "profile_quantiles", PROFILE_QUANTILES):
        issues.append(
            ValidationIssue(
                "PRF010",
                Severity.ERROR,
                str(root / "mirrors" / "profile_quantiles"),
                "profiles rows are present but mirrors/profile_quantiles/ holds "
                "no shard with any row; every percentile would be silently absent",
            )
        )
        issues.append(
            ValidationIssue(
                "PRF000",
                Severity.WARNING,
                str(root / "mirrors" / "profile_quantiles"),
                "skipped percentile checks: profiles rows exist but no profile_quantiles data does",
            )
        )
    # Same reasoning for the interpretive tree, which holds the curator
    # judgement the atlas exists to record.
    issues.extend(unexpected_curation_entries(root))

    # Referential checks over a corpus that failed to load would compare against a
    # knowingly incomplete set. The file wrappers validate their whole record list
    # as one model, so a single malformed publication empties `corpus.publications`
    # entirely — and every assertion citing a perfectly valid PMID would then be
    # reported as a dangling reference. Hundreds of false positives would bury the
    # one real error. Report the skip instead.
    if corpus_issues:
        issues.append(
            ValidationIssue(
                "REF000",
                Severity.WARNING,
                str(root / "curation"),
                f"skipped referential checks: {len(corpus_issues)} record(s) failed to load",
            )
        )
    else:
        known_genes = _known_genes(root)
        if known_genes is None:
            issues.append(
                ValidationIssue(
                    "TBL008",
                    Severity.ERROR,
                    str(root / "mirrors" / "genes.tsv"),
                    "gene registry is missing or unreadable, so gene references cannot be checked",
                )
            )
        # Without the phenotype vocabulary, REF007/REF009/REF010 silently no-op:
        # `lesion_group_of` is empty so every lookup misses. An error, not a
        # warning, because `ok` ignores warnings and a check stopped running.
        if corpus.assertions and not corpus.phenotypes:
            issues.append(
                ValidationIssue(
                    "CUR002",
                    Severity.ERROR,
                    str(root / "curation" / "phenotypes.yaml"),
                    "phenotype vocabulary is missing or empty, so lesion-group "
                    "consistency checks cannot run",
                )
            )
        # Otherwise every citation reports REF002 and nothing names the cause.
        if (corpus.assertions or corpus.functional) and not corpus.publications:
            issues.append(
                ValidationIssue(
                    "CUR003",
                    Severity.ERROR,
                    str(root / "curation" / "publications.yaml"),
                    "publication registry is missing or empty, so every citation "
                    "will be reported as dangling",
                )
            )
        issues.extend(validate_references(corpus, known_genes=known_genes))
        issues.extend(validate_mirror_references(root, corpus))
        # Deliberately not folded into the `validate_burden` call above: these
        # need registries a failed corpus load empties, and the comparator rules
        # must keep running when that happens.
        #
        # **Each registry is passed as None when it is empty, not as an empty
        # set**, and this is the load-bearing part. The first version passed the
        # sets unconditionally and broke the guarantee in §1 of CLAUDE.md that a
        # missing registry must not report as hundreds of dangling references.
        # Measured 2026-08-05 against the committed corpus:
        #
        #   deleting curation/cohorts.yaml   -> 4 BUR009 errors, none naming the cause
        #   header-only mirrors/genes.tsv    -> 146 errors (1 REF001 + 145 BUR011),
        #                                       against 1 error before this table existed
        #
        # `_known_genes` returns None only for a missing or unreadable file; a
        # header-only mirror reads fine and yields an empty set, which is why the
        # emptiness test lives here rather than there.
        burden_path = root / "mirrors" / "burden.tsv"
        cohort_ids = {str(cohort.id) for cohort in corpus.cohorts}
        study_ids = {str(publication.id) for publication in corpus.publications}
        if burden_path.is_file() and not cohort_ids:
            issues.append(
                ValidationIssue(
                    "CUR004",
                    Severity.ERROR,
                    str(root / "curation" / "cohorts.yaml"),
                    "cohort registry is missing or empty, so every cohort a burden "
                    "row cites will be reported as dangling",
                )
            )
        # CUR003 above fires only when an assertion or a functional record needs
        # the registry. A burden row cites a publication too, so the same cause
        # needs the same report when burden is the only thing citing one.
        if burden_path.is_file() and not study_ids and not (corpus.assertions or corpus.functional):
            issues.append(
                ValidationIssue(
                    "CUR003",
                    Severity.ERROR,
                    str(root / "curation" / "publications.yaml"),
                    "publication registry is missing or empty, so every citation "
                    "will be reported as dangling",
                )
            )
        issues.extend(
            validate_burden_references(
                root,
                known_cohorts=cohort_ids or None,
                known_genes=known_genes or None,
                known_studies=study_ids or None,
            )
        )
        issues.extend(
            validate_burden_combined_references(
                root,
                known_genes=known_genes or None,
                known_studies=study_ids or None,
            )
        )
        # Inside the same branch as the reference checks above, and for the same reason: it
        # reads the gene registry, so on a corpus that failed to load it would
        # report a missing accession for every gene at once.
        issues.extend(validate_ptm_evidence_is_reachable(root, corpus))

        # Same branch, same reasoning again: `corpus.datasets` and
        # `corpus.cardiac_phases` would both be empty on a corpus that failed
        # to load, making every declared cardiac tissue and stage look
        # unbacked by data and vice versa.
        #
        # Called once here rather than duplicated inside the `mirrored`
        # if/else immediately below -- unlike the scope validators there,
        # `datasets`, `known_genes` and `phases` do not come from either
        # validity mirror, so the same call answers both branches. Deliberate,
        # not an oversight (CLAUDE.md section 4.34: a validator called on two
        # branches is tested on the branch you were thinking about) --
        # `test_prf_references_run_on_both_scope_branches` pins that this
        # single call site is in fact reached on both the mirror-readable and
        # mirrors-missing (TBL012) branches.
        # The gate's population, once, for PRF009 and GEN004 below.
        published = _gate_published_genes(root, corpus)
        issues.extend(
            validate_profile_references(
                root,
                datasets=corpus.datasets,
                # `or None` for the reason the burden call above uses it: an
                # unreadable registry must skip the gene check, not report every
                # row as dangling. Passing the empty set bare made PRF007 fire
                # 154 times -- once per mirrored gene -- which is precisely the
                # cascade REF000/SRC000/ONT000 exist to prevent. Measured when
                # the first profiles mirror landed; the unit test passed `None`
                # explicitly and so never exercised this call site.
                known_genes=known_genes or None,
                published_genes=published,
                phases=corpus.cardiac_phases,
            )
        )
        # CTP001-011: cell-type rows against their dataset record and the rules
        # the cell-type block prints of them. The same `or None` as above, for the same
        # reason: an unreadable or header-only registry skips CTP007 rather
        # than flagging every gene; TBL008 or REF001 reports why. The other
        # checks still run, which
        # `test_a_header_only_gene_registry_skips_only_the_cell_type_gene_check`
        # pins: with this call deleted the whole suite once passed (review,
        # 2026-10-28).
        issues.extend(
            validate_cell_type_profiles(
                root, datasets=corpus.datasets, known_genes=known_genes or None
            )
        )
        # GEN004 (owner decision 2026-10-26, R13): the identity ids every
        # published gene page links by, read by the rule the build refuses a
        # page on, through the reader the build uses (`read_gene_identities`).
        # In this branch because the published population is the gate's,
        # computed from the curated scope a failed corpus load empties.
        #
        # It passed `_gene_identities(root) or None` until review (2026-10-27),
        # so that a header-only registry -- which reads, as empty, and which
        # TBL008 does not report -- would not report all 119 published genes'
        # ids missing. That also switched the check off: with no curated
        # assertion to raise REF001, such a registry validated with 0 errors
        # and `build_site` raised from `external_links` after writing 133
        # files (measured on a copy of 5f7161b). An empty registry is now one
        # GEN004 naming the file, from `validate_identity_ids`.
        issues.extend(
            validate_identity_ids(
                read_gene_identities(root), published, str(root / "mirrors" / "genes.tsv")
            )
        )

        # PRF013/PRF000 -- the same pairing PRF010/PRF000 gives the quantile
        # mirror, for the reference vocabulary rather than the reference
        # table, and here for the same reason: distinguishing "the data is
        # wrong" from "the reference data never loaded" is this project's
        # core rule, and the phase vocabulary was the one input that had no
        # such guard.
        #
        # `_prf006_issues` returns `[]` the moment `phases is None`, so an
        # absent `curation/cardiac_phases.yaml` skipped every phase check in
        # silence. Measured on the committed corpus: deleting the file took
        # `validate` from 4 warnings to **3** -- the count went *down* -- and
        # the build stayed green while 460 stage entries flipped from
        # `matched` to `outside_window` and all 85 charted pages printed
        # "none of this dataset's developmental stages falls inside a curated
        # cardiac phase this atlas can band". That sentence's subject is the
        # dataset's stages, and its cause was the atlas's own vocabulary.
        #
        # An empty `phases` list is treated identically and is not the same
        # input as a malformed file: `phases: []` alone fails
        # `CardiacPhaseFile`'s required `attributed_to`/`citation` and is
        # already SCHEMA001, but a file carrying both and declaring no phases
        # parses cleanly and is exactly as silent as an absent one. A guard
        # written as `is None` would catch only the first.
        #
        # Conditioned on `_phase_matchable_stages`, per section 4.41: a corpus
        # whose profiles rows are all post-natal (or all undeclared) has no
        # phase assignment to lose, and firing there would put a lone PRF000
        # in a report `ValidationReport.ok` calls clean.
        phases = corpus.cardiac_phases
        if phases is None or not phases.phases:
            matchable = _phase_matchable_stages(root, corpus.datasets)
            if matchable:
                issues.append(
                    ValidationIssue(
                        "PRF013",
                        Severity.ERROR,
                        str(root / "curation" / "cardiac_phases.yaml"),
                        f"{len(matchable)} profiles (dataset, stage) pair(s) carry a "
                        "wpc that a cardiac phase would be matched against, but no "
                        "phase vocabulary is curated; every one would publish as "
                        "'outside_window' and every chart would lose its bands",
                    )
                )
                issues.append(
                    ValidationIssue(
                        "PRF000",
                        Severity.WARNING,
                        str(root / "curation" / "cardiac_phases.yaml"),
                        "skipped cardiac phase checks: profiles rows would be "
                        "phase-matched but no phase vocabulary is curated",
                    )
                )

        # Same branch, same reasoning again: a corpus that failed to load empties
        # `corpus.chd_scope`, and every scope term would then look absent.
        scope_location = str(root / "curation" / "chd_scope.yaml")
        # SCP006/SCP008. One call site, before the validity-mirror if/else and
        # reached on both of its branches, because an Orphanet warrant is checked
        # against the Orphanet mirrors and the pinned HPO release, neither of which
        # is a validity mirror (CLAUDE.md section 4.34: a call that appears in an
        # `if` and an `else` is two call sites).
        orphanet_tables = load_orphanet_tables(root)
        structural = pinned_structural_phenotypes(root, registry.ontology_files())
        issues.extend(
            validate_orphanet_warrants(
                corpus.chd_scope, orphanet_tables, structural, scope_location
            )
        )
        mirrored = _mirrored_validity(root)
        if mirrored is None:
            # SCP000 alone is a WARNING, and `ValidationReport.ok` ignores
            # warnings — so without an accompanying ERROR, a repository missing
            # both validity mirrors would validate "clean" while every scope
            # check silently could not run, and `build_site` would publish it
            # anyway. Same shape as `TBL008` for a missing gene registry: name
            # the cause as an error, located at the mirrors directory because
            # two files are implicated and neither alone explains it.
            issues.append(
                ValidationIssue(
                    "TBL012",
                    Severity.ERROR,
                    str(root / "mirrors"),
                    "the validity mirrors are missing or unreadable, so no gene "
                    "can be placed in scope",
                )
            )
            issues.extend(validate_scope_terms(corpus.chd_scope, None, scope_location))
            issues.extend(validate_scope_attribution(corpus.chd_scope, None, None, scope_location))
        else:
            labels, genes_by_disease, diseases_by_gene, panels, submitters = mirrored
            issues.extend(validate_scope_terms(corpus.chd_scope, labels, scope_location))
            issues.extend(
                validate_scope_attribution(corpus.chd_scope, panels, submitters, scope_location)
            )
            scope_terms = {str(entry.id) for entry in corpus.chd_scope}
            in_scope_genes = {
                gene for term in scope_terms for gene in genes_by_disease.get(term, set())
            }
            issues.extend(
                scope_candidates(
                    scope_terms, labels, diseases_by_gene, in_scope_genes, scope_location
                )
            )
            # SCP007 needs the mirrored labels, so it lives on this branch only,
            # beside SCP003, the other "terms a curator should look at" digest.
            issues.extend(
                orphanet_candidates(
                    scope_terms, labels, orphanet_tables, structural, scope_location
                )
            )
            # `Corpus` does not record which file each assertion was loaded from
            # (see `corpus.py::Corpus`), so every curated gene is attributed to
            # the assertions directory rather than a specific file within it —
            # coarser than a curator fixing SCP004 would want, but not a
            # precision this module has and is pretending not to.
            curated_genes = {
                str(assertion.gene): str(root / "curation" / "assertions")
                for assertion in corpus.assertions
            }
            issues.extend(validate_curation_is_in_scope(curated_genes, in_scope_genes))

    # Issue #33: every cross-source join keys on HGNC id, and the registry that
    # makes that possible has to stay able to resolve the symbols the mirrors
    # actually print. Independent of scope -- a stale registry entry is a
    # naming fact about a gene, not about which diseases it is curated for --
    # so this runs on its own rather than inside either scope branch, which is
    # also what keeps it from being a call site that exists on one branch and
    # not the other (CLAUDE.md section 4.34).
    issues.extend(
        validate_mirror_symbols(
            _mirror_gene_labels(root),
            _symbol_registry(root),
            str(root / "mirrors" / "genes.tsv"),
        )
    )

    # Same reasoning as the referential skip above: on a failed registry load
    # `registry` is empty, so every source every mirror table uses would report
    # SRC003 and bury the one real error.
    if source_issues:
        issues.append(
            ValidationIssue(
                "SRC000",
                Severity.WARNING,
                str(root / "mirrors" / "sources.yaml"),
                "skipped source-reference checks: the source registry did not load",
            )
        )
    else:
        issues.extend(validate_source_references(_used_sources(root), registry))

    # The ontology pins live in sources.yaml, so a registry that failed to load
    # pins nothing and every term reports ONT003 "no pinned ontology" — one per
    # term, burying the single SRC001 that caused them. Same cascade as above.
    # ...and likewise when a pinned release itself failed to load: ONT003 reads
    # "no pinned ontology for prefix 'HP'", which is false when the prefix IS
    # pinned and the file is simply unreadable. A curator would go add a pin
    # that already exists. The ONT004 naming the real cause is already reported.
    if source_issues or ontologies.load_issues:
        cause = (
            "the source registry did not load, so no releases are pinned"
            if source_issues
            else "a pinned ontology release failed to load"
        )
        issues.append(
            ValidationIssue(
                "ONT000",
                Severity.WARNING,
                str(root / "mirrors" / "sources.yaml"),
                f"skipped ontology checks: {cause}",
            )
        )
    else:
        for assertion in corpus.assertions:
            issues.extend(
                validate_terms(
                    [*assertion.phenotypes, *assertion.extracardiac_features],
                    ontologies,
                    f"assertion {assertion.id}",
                )
            )
        for record in corpus.functional:
            issues.extend(
                validate_terms(
                    [record.organism, *record.cardiac_phenotype],
                    ontologies,
                    f"functional evidence {record.id}",
                )
            )
        phenotypes_file = str(root / "curation" / "phenotypes.yaml")
        for term in corpus.phenotypes:
            issues.extend(validate_terms([term.id], ontologies, phenotypes_file))
        # The file declares its labels transcribed from the pinned release, and
        # those labels are what the atlas displays. Pinning makes the claim
        # checkable, so check it.
        issues.extend(
            validate_labels(
                [(term.id, term.label) for term in corpus.phenotypes],
                ontologies,
                phenotypes_file,
            )
        )

    id_registry, id_issues = load_id_registry(root)
    issues.extend(id_issues)
    # None means the counter file is missing or malformed, which ID001/ID004
    # already named. Checking against an empty registry would additionally make
    # every atlas ID exceed its ceiling, burying that one error under one per record.
    if id_registry is not None:
        issues.extend(
            validate_ids(
                [assertion.id for assertion in corpus.assertions]
                + [record.id for record in corpus.functional],
                id_registry,
            )
        )

    # Relativised before the report is built, so `__post_init__` sorts the
    # locations the reader will actually see.
    return ValidationReport(issues=_relative_to_root(issues, root))
