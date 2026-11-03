# src/chd_atlas/validate/referential.py
"""Cross-entity reference integrity across the curation corpus and mirror tables."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Final

from chd_atlas.corpus import Corpus
from chd_atlas.duplicates import duplicates
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table
from chd_atlas.vocab import EvidenceClass, LesionGroup

# Mirror tables whose rows point back at a curated dataset (and usually a contrast).
_DATASET_LINKED_TABLES: Final[frozenset[str]] = frozenset(
    {"expression", "profiles", "proteomics", "phospho"}
)


def lesion_group_faults(
    phenotypes: Iterable[str],
    lesion_groups: Iterable[LesionGroup],
    lesion_group_of: Mapping[str, LesionGroup | None],
) -> list[tuple[str, str]]:
    """REF007 and REF009 for one record's phenotypes against its lesion groups.

    Shared by lesion assertions and literature findings (final review I3): both
    carry `phenotypes` and `lesion_groups`, and the mapping lives in
    curation/phenotypes.yaml, which no single-record model can see. REF007 is a
    mapped phenotype whose group is not declared; REF009 is a declared group no
    phenotype justifies, checked only when every phenotype is mapped. An unmapped
    phenotype is exempt from both (the assertion path reports it as REF012).
    Returns (code, message) pairs, REF007s first.
    """
    declared_groups = list(lesion_groups)
    terms = list(phenotypes)
    faults: list[tuple[str, str]] = []
    for phenotype in terms:
        expected = lesion_group_of.get(phenotype)
        if expected is not None and expected not in declared_groups:
            declared = sorted(group.value for group in declared_groups)
            faults.append(
                (
                    "REF007",
                    f"phenotype {phenotype} belongs to lesion group "
                    f"'{expected.value}', which is not among {declared}",
                )
            )
    if all(phenotype in lesion_group_of for phenotype in terms):
        justified = {lesion_group_of[p] for p in terms}
        unjustified = sorted(group.value for group in declared_groups if group not in justified)
        if unjustified:
            faults.append(
                (
                    "REF009",
                    f"lesion groups {unjustified} are not justified by any declared phenotype",
                )
            )
    return faults


def validate_references(corpus: Corpus, known_genes: set[str] | None) -> list[ValidationIssue]:
    """Check that every reference in the corpus resolves to something real.

    ``known_genes`` comes from ``mirrors/genes.tsv``; it is passed in rather
    than read here so this function stays a pure check over loaded data.

    ``None`` means that registry could not be read, so the gene checks are
    skipped. An empty set is indistinguishable from a missing file, and a
    missing file would otherwise report one REF001 per assertion and per
    functional record — hundreds of failures naming the cause nowhere. The
    caller reports the unreadable registry once instead.
    """
    issues: list[ValidationIssue] = []
    known_publications = {publication.id for publication in corpus.publications}
    functional_by_id = {record.id: record for record in corpus.functional}
    known_datasets = {dataset.id for dataset in corpus.datasets}
    lesion_group_of: dict[str, LesionGroup | None] = {
        term.id: term.lesion_group for term in corpus.phenotypes
    }

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    for assertion in corpus.assertions:
        location = f"assertion {assertion.id}"

        if known_genes is not None and assertion.gene not in known_genes:
            error("REF001", location, f"gene {assertion.gene} is not in mirrors/genes.tsv")

        # An assertion pairing an atrial septal defect with lesion_groups [conotruncal]
        # is internally inconsistent. The mapping lives in curation/phenotypes.yaml, so
        # the single-record model cannot see it — it has to be checked here.
        for code, message in lesion_group_faults(
            assertion.phenotypes, assertion.lesion_groups, lesion_group_of
        ):
            if code == "REF007":
                error(code, location, message)

        # REF009 below exempts an assertion with any unmapped phenotype, which is
        # sound but silent: one term missing from the register disables all
        # lesion-group checking for that record with nothing saying so. Naming the
        # term restores the report.
        #
        # Guarded on a non-empty register for the same reason `known_genes is
        # None` guards REF001: when phenotypes.yaml is absent or empty every term
        # is unmapped, and one issue per phenotype per assertion would bury the
        # single CUR002 that names the cause.
        if lesion_group_of:
            for phenotype in assertion.phenotypes:
                if phenotype not in lesion_group_of:
                    error(
                        "REF012",
                        location,
                        f"phenotype {phenotype} is not in curation/phenotypes.yaml, "
                        f"so its lesion group cannot be checked",
                    )

        # The converse of REF007: a declared lesion group that no phenotype accounts for.
        # It files the assertion under a browse facet nothing in the record supports.
        # Sound only when every phenotype is mapped — an unmapped one could legitimately
        # justify the extra group, which is why REF007 exempts unmapped phenotypes too.
        for code, message in lesion_group_faults(
            assertion.phenotypes, assertion.lesion_groups, lesion_group_of
        ):
            if code == "REF009":
                error(code, location, message)

        # curation/phenotypes.yaml also registers extracardiac terms, purely so
        # their labels get checked (`PhenotypeTerm.lesion_group` is `None` for
        # those), so `lesion_group_of.get(feature)` is `None` for a genuinely
        # extracardiac term and only non-`None` for a term that is *also* a
        # cardiac lesion -- one that cannot be what makes the assertion syndromic.
        for feature in assertion.extracardiac_features:
            feature_group = lesion_group_of.get(feature)
            if feature_group is not None:
                error(
                    "REF010",
                    location,
                    f"extracardiac feature {feature} is a cardiac lesion in "
                    f"group '{feature_group.value}'",
                )

        for index, evidence in enumerate(assertion.evidence):
            evidence_location = f"{location} evidence[{index}]"

            if evidence.publication not in known_publications:
                error(
                    "REF002",
                    evidence_location,
                    f"publication {evidence.publication} is not in curation/publications.yaml",
                )
            if evidence.functional_evidence is not None:
                cited = functional_by_id.get(evidence.functional_evidence)
                if cited is None:
                    error(
                        "REF003",
                        evidence_location,
                        f"functional evidence {evidence.functional_evidence} does not exist",
                    )
                # Functional evidence feeds ClinGen-style scoring, so a record about
                # another gene can carry an assertion to 'definitive' on model data
                # that says nothing about the gene being asserted.
                elif cited.gene != assertion.gene:
                    error(
                        "REF008",
                        evidence_location,
                        f"functional evidence {evidence.functional_evidence} is about "
                        f"gene {cited.gene}, not {assertion.gene}",
                    )
            if evidence.dataset is not None and evidence.dataset not in known_datasets:
                error(
                    "REF004",
                    evidence_location,
                    f"dataset {evidence.dataset} is not in curation/datasets/",
                )

    for record in corpus.functional:
        location = f"functional evidence {record.id}"
        if known_genes is not None and record.gene not in known_genes:
            error("REF001", location, f"gene {record.gene} is not in mirrors/genes.tsv")
        if record.publication not in known_publications:
            error(
                "REF002",
                location,
                f"publication {record.publication} is not in curation/publications.yaml",
            )

    for entry in corpus.featured:
        if entry.publication not in known_publications:
            error(
                "REF002",
                f"featured entry {entry.order}",
                f"publication {entry.publication} is not in curation/publications.yaml",
            )

    # Each dataset is its own file under curation/datasets/, so no single model
    # can see a collision between two of them. `validate_ids` does not cover it
    # either: it checks only atlas-minted CHDA identifiers, and an accession is a
    # third-party one. Two files declaring PXD012345 with different sample counts
    # both load, and every mirror row citing that accession then resolves to
    # whichever the filesystem happened to yield first.
    for repeated in sorted(duplicates(dataset.id for dataset in corpus.datasets)):
        error(
            "REF011",
            f"dataset {repeated}",
            "accession is declared by more than one file under curation/datasets/",
        )

    for dataset in corpus.datasets:
        if dataset.publication is not None and dataset.publication not in known_publications:
            error(
                "REF002",
                f"dataset {dataset.id}",
                f"publication {dataset.publication} is not in curation/publications.yaml",
            )

    return issues


def validate_mirror_references(root: Path, corpus: Corpus) -> list[ValidationIssue]:
    """Check that mirror rows point at datasets and contrasts that actually exist.

    Distinct (dataset, contrast) pairs are checked once rather than per row, so a
    100,000-row table with one bad contrast yields one issue, not 100,000.
    """
    known_datasets = {dataset.id for dataset in corpus.datasets}
    known_contrasts = {
        (dataset.id, contrast.id) for dataset in corpus.datasets for contrast in dataset.contrasts
    }
    issues: list[ValidationIssue] = []

    for path, schema_name in mirror_paths(root):
        if schema_name not in _DATASET_LINKED_TABLES:
            continue

        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None or "dataset" not in frame.columns:
            continue

        has_contrast = "contrast" in frame.columns
        selected = ["dataset", "contrast"] if has_contrast else ["dataset"]
        # A null dataset cell is a schema violation `validate_table` already reports as
        # TBL003 — both columns are non-nullable in every linked schema. Carrying it
        # here would add "dataset None is not in curation/datasets/", naming a dataset
        # that was never written down. A null contrast is skipped below, by the same
        # guard that skips tables having no contrast column at all.
        pairs = {
            (row[0], row[1] if has_contrast else None)
            for row in frame.select(selected).rows()
            if row[0] is not None
        }

        for dataset, contrast in sorted(pairs, key=lambda pair: (str(pair[0]), str(pair[1]))):
            if dataset not in known_datasets:
                issues.append(
                    ValidationIssue(
                        "REF005",
                        Severity.ERROR,
                        str(path),
                        f"dataset {dataset} is not in curation/datasets/",
                    )
                )
            elif contrast is not None and (dataset, contrast) not in known_contrasts:
                issues.append(
                    ValidationIssue(
                        "REF006",
                        Severity.ERROR,
                        str(path),
                        f"contrast '{contrast}' is not defined by dataset {dataset}",
                    )
                )

    return issues


def validate_ptm_evidence_is_reachable(root: Path, corpus: Corpus) -> list[ValidationIssue]:
    """Warn when a gene asserts PTM evidence but has no accession to join it by.

    `mirrors/phospho/` carries no gene column at all, so the only route from a
    phospho row to a gene is `mirrors/genes.tsv:uniprot`, and that column is
    nullable. A gene whose cell is blank has its PTM sites published in a shard
    and summarised in no bundle: `"omics": {}` on the gene page while the browse
    row beside it advertises `evidence_counts: {"ptm": 1}` from the assertion.
    Reproduced end to end on a build reporting 0 errors and 0 warnings, which is
    what makes it this project's characteristic failure rather than a nuisance.

    Keyed on the curated side, not the mirror. Scanning the omics tables for
    unmatched accessions would report every protein the atlas has no gene for,
    which for a real proteomics release is most of them — the cascade REF000 and
    SRC000 exist to prevent. An assertion that *claims* PTM evidence is a small,
    curated set where every hit is actionable: the curator wrote the claim, so
    the missing cell is theirs to fill.

    Exact rather than heuristic. REF001 already guarantees the asserting gene is
    in the registry, so a null `uniprot` is the only remaining hole.

    A WARNING, not an error. The evidence is real and the assertion is sound; what
    is missing is one mirror cell, and refusing to build the site over it would be
    the wrong trade for a resource curated incrementally. `proteomics` is
    deliberately not covered — its rows carry their own nullable `gene` column, so
    an accession is one of two routes there rather than the only one.
    """
    asserted_ptm = sorted(
        {
            (str(assertion.gene), str(assertion.id))
            for assertion in corpus.assertions
            for evidence in assertion.evidence
            if evidence.evidence_class is EvidenceClass.PTM
        }
    )
    if not asserted_ptm:
        return []

    path = root / "mirrors" / "genes.tsv"
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns or "uniprot" not in frame.columns:
        # A registry that did not load is TBL008 and a missing column is TBL002.
        # Reporting "no accession" for every gene on top of either is the cascade
        # REF000 exists to avoid.
        return []

    accessions = {
        str(row["hgnc_id"]): row["uniprot"]
        for row in frame.select(["hgnc_id", "uniprot"]).to_dicts()
        if row["hgnc_id"] is not None
    }
    return [
        ValidationIssue(
            "REF013",
            Severity.WARNING,
            str(path),
            f"{gene} asserts PTM evidence in {assertion_id} but has no uniprot accession, "
            f"so its phospho sites can reach no gene page",
        )
        for gene, assertion_id in asserted_ptm
        if not accessions.get(gene)
    ]
