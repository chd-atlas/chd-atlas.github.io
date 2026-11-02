# src/chd_atlas/corpus.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from pydantic import BaseModel, ValidationError
from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.fs import list_dir
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.assertion import AssertionFile, LesionAssertion
from chd_atlas.models.cohort import Cohort, CohortFile
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalEvidence, FunctionalFile
from chd_atlas.models.gene_literature import GeneLiteratureFile
from chd_atlas.models.literature import (
    FeaturedFile,
    FeaturedManuscript,
    PhenotypeFile,
    PhenotypeTerm,
    Publication,
    PublicationFile,
)
from chd_atlas.models.phases import CardiacPhaseFile
from chd_atlas.models.scope import ChdScopeFile, ScopeEntry


@dataclass(frozen=True)
class Corpus:
    """Every interpretive record loaded from ``curation/``."""

    root: Path
    assertions: tuple[LesionAssertion, ...] = ()
    functional: tuple[FunctionalEvidence, ...] = ()
    publications: tuple[Publication, ...] = ()
    featured: tuple[FeaturedManuscript, ...] = ()
    phenotypes: tuple[PhenotypeTerm, ...] = ()
    datasets: tuple[Dataset, ...] = ()
    chd_scope: tuple[ScopeEntry, ...] = ()
    cohorts: tuple[Cohort, ...] = ()
    cardiac_phases: CardiacPhaseFile | None = None
    # One entry per `curation/literature/*.yaml`, filename-sorted. The path is
    # kept because the file name is part of the record (LIT009) and names the
    # location of an error.
    literature: tuple[tuple[Path, GeneLiteratureFile], ...] = ()


@dataclass
class _Accumulator:
    issues: list[ValidationIssue] = field(default_factory=list)

    def error(self, code: str, location: Path | str, message: str) -> None:
        self.issues.append(
            ValidationIssue(
                code=code,
                severity=Severity.ERROR,
                location=str(location),
                message=message,
            )
        )


# Distinguishes "could not read the file" from "the file parsed to null".
# Returning None for both would let an empty YAML file silently load as no
# records and no issues.
_UNREADABLE: Final = object()


def _read_yaml(path: Path, acc: _Accumulator) -> Any:
    yaml = YAML(typ="safe")
    try:
        return yaml.load(path.read_text(encoding="utf-8"))
    # ruamel builds a `datetime.date` directly and raises a bare ValueError on
    # an out-of-range literal such as `curated_on: 2026-13-45`, which YAMLError
    # does not cover. UnicodeDecodeError is already a ValueError subclass; it is
    # kept listed to document that non-UTF-8 bytes are handled here too.
    # Pydantic's ValidationError is also a ValueError but is raised in `_parse`,
    # so it is not swallowed by this guard.
    #
    # OSError covers everything `read_text` can raise before parsing begins: an
    # unreadable mode, a dangling symlink, a directory where a file was
    # expected. One unreadable file must be reported as one issue rather than
    # aborting the whole run with a traceback, exactly as `read_table` does for
    # the TSV side.
    except (YAMLError, UnicodeDecodeError, ValueError, OSError) as exc:
        acc.error("YAML001", path, f"could not read YAML: {exc}")
        return _UNREADABLE


def _parse[ModelT: BaseModel](model: type[ModelT], path: Path, acc: _Accumulator) -> ModelT | None:
    raw = _read_yaml(path, acc)
    if raw is _UNREADABLE:
        return None
    try:
        return model.model_validate(raw)
    except ValidationError as exc:
        for error in exc.errors():
            location = ".".join(str(part) for part in error["loc"])
            acc.error("SCHEMA001", f"{path}:{location}", error["msg"])
        return None


def _record_files(directory: Path) -> list[Path]:
    """Every ``*.yaml`` file in one record directory, in a stable order.

    The `is_file` guard matters because ``glob`` matches on name alone: a
    *directory* named ``TBX5.yaml`` would otherwise be handed to `read_text`.
    Skipping it here is safe only because `unexpected_curation_entries` reports
    it — silently skipping is what this whole family of checks exists to stop.
    """
    return sorted(path for path in directory.glob("*.yaml") if path.is_file())


def unexpected_curation_entries(root: Path) -> list[ValidationIssue]:
    """Anything under ``curation/`` that no loader claims.

    The loaders glob for what they expect, so anything misnamed is simply never
    seen: an assertion saved as ``TBX5.yml``, a directory misspelled
    ``curation/assertion/``, a directory where a record file belongs. Each makes
    curator judgement vanish from validation while the gate still passes, which
    is the worst failure this tool can have. Naming the stray entry catches it.

    An *absent* expected file or directory stays legitimate — there is no
    ``functional/`` before any functional evidence has been curated.
    """
    curation = root / "curation"
    if not curation.is_dir():
        return []

    expected_files = {
        ".id_registry.yaml",
        "cardiac_phases.yaml",
        "chd_scope.yaml",
        "cohorts.yaml",
        "featured.yaml",
        "phenotypes.yaml",
        "publications.yaml",
    }
    expected_dirs = {"assertions", "datasets", "functional", "literature"}
    issues: list[ValidationIssue] = []

    def error(entry: Path, message: str) -> None:
        issues.append(ValidationIssue("CUR001", Severity.ERROR, str(entry), message))

    entries, listing_issues = list_dir(curation, "CUR001")
    issues.extend(listing_issues)

    for entry in entries:
        if entry.is_dir():
            if entry.name not in expected_dirs:
                error(
                    entry,
                    f"unexpected directory under curation/; expected one of "
                    f"{sorted(expected_dirs)}",
                )
                continue
            records, record_listing_issues = list_dir(entry, "CUR001")
            issues.extend(record_listing_issues)
            for record in records:
                if not record.is_file():
                    error(record, f"'{record.name}' should be a YAML file, not a directory")
                elif record.suffix != ".yaml":
                    error(record, f"record files must end .yaml; '{record.name}' does not")
        # Checked before the unexpected-file case: the name sets are disjoint, so
        # testing expected_files first would report a file named `assertions` as
        # merely unnamed rather than as the wrong kind of entry.
        elif entry.name in expected_dirs:
            error(entry, f"'{entry.name}' should be a directory of records, not a file")
        elif entry.name not in expected_files:
            error(
                entry,
                f"unexpected file under curation/; expected one of {sorted(expected_files)}",
            )
    return issues


def load_curation(root: Path) -> tuple[Corpus, list[ValidationIssue]]:
    """Load every YAML record under ``root/curation``.

    Errors are accumulated rather than raised so that a single run reports
    every problem in the corpus, not just the first one encountered.
    """
    acc = _Accumulator()
    curation = root / "curation"
    if not curation.is_dir():
        acc.error("CORPUS001", curation, "curation directory not found")
        return Corpus(root=root), acc.issues

    assertions: list[LesionAssertion] = []
    for path in _record_files(curation / "assertions"):
        parsed = _parse(AssertionFile, path, acc)
        if parsed is not None:
            assertions.extend(parsed.assertions)

    functional: list[FunctionalEvidence] = []
    for path in _record_files(curation / "functional"):
        parsed_functional = _parse(FunctionalFile, path, acc)
        if parsed_functional is not None:
            functional.extend(parsed_functional.functional_evidence)

    datasets: list[Dataset] = []
    for path in _record_files(curation / "datasets"):
        parsed_dataset = _parse(Dataset, path, acc)
        if parsed_dataset is not None:
            datasets.append(parsed_dataset)

    literature: list[tuple[Path, GeneLiteratureFile]] = []
    for path in _record_files(curation / "literature"):
        parsed_literature = _parse(GeneLiteratureFile, path, acc)
        if parsed_literature is not None:
            literature.append((path, parsed_literature))

    publications: tuple[Publication, ...] = ()
    publications_path = curation / "publications.yaml"
    if publications_path.is_file():
        parsed_publications = _parse(PublicationFile, publications_path, acc)
        if parsed_publications is not None:
            publications = tuple(parsed_publications.publications)

    featured: tuple[FeaturedManuscript, ...] = ()
    featured_path = curation / "featured.yaml"
    if featured_path.is_file():
        parsed_featured = _parse(FeaturedFile, featured_path, acc)
        if parsed_featured is not None:
            featured = tuple(parsed_featured.featured)

    phenotypes: tuple[PhenotypeTerm, ...] = ()
    phenotypes_path = curation / "phenotypes.yaml"
    if phenotypes_path.is_file():
        parsed_phenotypes = _parse(PhenotypeFile, phenotypes_path, acc)
        if parsed_phenotypes is not None:
            phenotypes = tuple(parsed_phenotypes.phenotypes)

    chd_scope: tuple[ScopeEntry, ...] = ()
    chd_scope_path = curation / "chd_scope.yaml"
    if chd_scope_path.is_file():
        parsed_chd_scope = _parse(ChdScopeFile, chd_scope_path, acc)
        if parsed_chd_scope is not None:
            chd_scope = tuple(parsed_chd_scope.diseases)

    cohorts: tuple[Cohort, ...] = ()
    cohorts_path = curation / "cohorts.yaml"
    if cohorts_path.is_file():
        parsed_cohorts = _parse(CohortFile, cohorts_path, acc)
        if parsed_cohorts is not None:
            cohorts = tuple(parsed_cohorts.cohorts)

    cardiac_phases: CardiacPhaseFile | None = None
    cardiac_phases_path = curation / "cardiac_phases.yaml"
    if cardiac_phases_path.is_file():
        cardiac_phases = _parse(CardiacPhaseFile, cardiac_phases_path, acc)

    corpus = Corpus(
        root=root,
        assertions=tuple(assertions),
        functional=tuple(functional),
        publications=publications,
        featured=featured,
        phenotypes=phenotypes,
        datasets=tuple(datasets),
        chd_scope=chd_scope,
        cohorts=cohorts,
        cardiac_phases=cardiac_phases,
        literature=tuple(literature),
    )
    return corpus, acc.issues
