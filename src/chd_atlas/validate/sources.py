"""The source registry: which external resource, which release, which licence.

Licence provenance is a correctness requirement rather than documentation,
because the atlas redistributes third-party data.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping
from datetime import date
from pathlib import Path, PurePosixPath

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    ValidationError,
    field_validator,
    model_validator,
)
from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.duplicates import duplicates
from chd_atlas.fs import list_dir
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.package import SEGMENT, ManifestError, load_manifest
from chd_atlas.tables import FLAT_TABLES, PACKAGES_DIR, mirror_paths
from chd_atlas.vocab import Redistribution


class PackagePin(BaseModel):
    """The data package a source's mirrored rows were ingested from (spec F11).

    `version` is the package's own version, so it takes the manifest's
    ``version`` pattern (never ``.`` or ``..``: it names the package's cache
    directory); `manifest_sha256` is the sha256 of its `manifest.json`, which
    is what makes the pin checkable.
    """

    model_config = ConfigDict(extra="forbid")

    version: str = Field(pattern=SEGMENT)
    manifest_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")


class Source(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str = Field(min_length=1)
    name: str = Field(min_length=1)
    version: str = Field(min_length=1)
    retrieved_on: date
    url: str = Field(min_length=1)
    licence: str = Field(min_length=1)
    redistribution: Redistribution
    ontology_prefix: str | None = None
    ontology_file: str | None = None
    # What the file itself says it is, beside a `version` that is a checksum
    # (owner decision 2026-10-26, R11). Orphadata's XML opens on
    # `<JDBOR date="..." version="...">`, and its citation asks for that data
    # version; a sha256 proves which file was mirrored but tells a reader
    # nothing they can cite. Transcribed verbatim from the header by the
    # chd-atlas/feeds `orphanet_product4` / `orphanet_product1` recipes
    # (`feeds.orphadata.jdbor_release`) into each package's manifest, and
    # recorded here beside the sha256 so the three describe one file.
    # Both or neither.
    file_version: str | None = None
    file_date: str | None = None
    # `None` when the source's rows were not ingested from a package (spec F11):
    # every source as of schema 2.22, and for good a source no mirror table
    # holds rows from, such as `hpo`.
    package: PackagePin | None = None

    @field_validator("licence")
    @classmethod
    def licence_is_not_blank(cls, value: str) -> str:
        """`min_length=1` accepts "   ", which records no provenance at all.

        The licence vocabulary itself is deliberately free text: the atlas
        mirrors resources whose terms have no SPDX identifier.
        """
        if not value.strip():
            raise ValueError("licence must not be blank")
        return value

    @model_validator(mode="after")
    def ontology_prefix_and_file_are_paired(self) -> Source:
        """Either both are set or neither is.

        `ontology_files()` is keyed on prefix and valued on file, so a source
        declaring one without the other is dropped from the mapping and the
        pin it was meant to express is never applied.
        """
        if self.ontology_prefix is not None and self.ontology_file is None:
            raise ValueError(
                f"source '{self.id}' declares ontology_prefix "
                f"'{self.ontology_prefix}' but no ontology_file"
            )
        if self.ontology_file is not None and self.ontology_prefix is None:
            raise ValueError(
                f"source '{self.id}' declares ontology_file "
                f"'{self.ontology_file}' but no ontology_prefix"
            )
        return self

    @model_validator(mode="after")
    def file_version_and_date_are_paired_and_not_blank(self) -> Source:
        """A file's own version means little without the date it states, and the reverse."""
        for label, value in (("file_version", self.file_version), ("file_date", self.file_date)):
            if value is not None and not value.strip():
                raise ValueError(f"source '{self.id}' {label} must not be blank")
        if (self.file_version is None) != (self.file_date is None):
            raise ValueError(
                f"source '{self.id}' declares one of file_version and file_date; "
                f"a file's header states both"
            )
        return self

    @model_validator(mode="after")
    def ontology_file_stays_inside_the_repository(self) -> Source:
        """The pin is resolved as ``root / ontology_file``.

        `pathlib` discards the left operand entirely when the right one is
        absolute, and ``..`` walks out of the repository, so either would
        resolve against the curator's machine rather than the checkout.
        """
        if self.ontology_file is None:
            return self
        # Repo-relative paths are declared posix-style regardless of platform.
        path = PurePosixPath(self.ontology_file)
        if path.is_absolute() or ".." in path.parts:
            raise ValueError(
                f"source '{self.id}' ontology_file '{self.ontology_file}' must be "
                f"a relative path inside the repository"
            )
        return self


class SourceRegistry(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sources: list[Source] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> SourceRegistry:
        found = duplicates(source.id for source in self.sources)
        if found:
            raise ValueError(f"duplicate source ids: {found}")
        return self

    @model_validator(mode="after")
    def ontology_prefixes_are_unique(self) -> SourceRegistry:
        """Two sources claiming one prefix would resolve by YAML line order.

        `ontology_files()` is keyed on prefix, so leaving a superseded release
        in place would silently change which release the whole corpus is
        validated against — the opposite of pinning it.
        """
        found = duplicates(
            source.ontology_prefix for source in self.sources if source.ontology_prefix is not None
        )
        if found:
            raise ValueError(f"duplicate ontology prefixes: {found}")
        return self

    @property
    def by_id(self) -> dict[str, Source]:
        return {source.id: source for source in self.sources}

    def ontology_files(self) -> dict[str, str]:
        """Map ontology prefix (``HP``) to its pinned OBO path."""
        return {
            source.ontology_prefix: source.ontology_file
            for source in self.sources
            if source.ontology_prefix is not None and source.ontology_file is not None
        }


def _empty_registry() -> SourceRegistry:
    """A fresh registry per failed load.

    A module-level singleton would be shared by every failed load in the
    process, so one accidental append to its mutable `sources` list would
    contaminate all the others.
    """
    return SourceRegistry.model_construct(sources=[])


def load_sources(root: Path) -> tuple[SourceRegistry, list[ValidationIssue]]:
    path = root / "mirrors" / "sources.yaml"
    if not path.is_file():
        return _empty_registry(), [
            ValidationIssue("SRC001", Severity.ERROR, str(path), "source registry not found")
        ]

    yaml = YAML(typ="safe")
    try:
        raw = yaml.load(path.read_text(encoding="utf-8"))
    # ruamel builds a `datetime.date` directly and raises a bare ValueError on
    # an out-of-range literal such as `2026-13-45`, which YAMLError does not
    # cover. UnicodeDecodeError is already a ValueError subclass; it is kept
    # listed to document that non-UTF-8 bytes are handled here too. Pydantic's
    # ValidationError is also a ValueError but is raised in the try below, so
    # it is not swallowed by this guard.
    # OSError covers what read_text can raise before parsing: an unreadable
    # mode, a dangling symlink, a directory where a file was expected. One
    # bad file must be one issue, not a traceback aborting the whole run.
    except (YAMLError, UnicodeDecodeError, ValueError, OSError) as exc:
        return _empty_registry(), [
            ValidationIssue("YAML001", Severity.ERROR, str(path), f"could not read YAML: {exc}")
        ]

    try:
        return SourceRegistry.model_validate(raw), []
    except ValidationError as exc:
        issues = [
            ValidationIssue(
                "SRC002",
                Severity.ERROR,
                str(path),
                f"{'.'.join(str(part) for part in error['loc'])}: {error['msg']}",
            )
            for error in exc.errors()
        ]
        return _empty_registry(), issues


def validate_source_references(
    used: dict[str, set[str]], registry: SourceRegistry
) -> list[ValidationIssue]:
    """Check every ``source`` value used by a mirror table may be redistributed.

    A record whose source is unregistered has no recorded licence, which is a
    redistribution risk rather than a cosmetic omission. A source whose licence
    is recorded as forbidding redistribution is the same risk with the answer
    already known, so mirroring it is an error rather than a warning.
    """
    by_id = registry.by_id
    issues: list[ValidationIssue] = []
    for location, sources in sorted(used.items()):
        for source in sorted(sources):
            entry = by_id.get(source)
            if entry is None:
                issues.append(
                    ValidationIssue(
                        "SRC003",
                        Severity.ERROR,
                        location,
                        f"source '{source}' has no entry in mirrors/sources.yaml (licence unknown)",
                    )
                )
            elif entry.redistribution is Redistribution.PROHIBITED:
                issues.append(
                    ValidationIssue(
                        "SRC004",
                        Severity.ERROR,
                        location,
                        f"source '{source}' is recorded as "
                        f"redistribution: {Redistribution.PROHIBITED.value} "
                        f"and must not be mirrored",
                    )
                )
    return issues


def validate_package_pins(root: Path, registry: SourceRegistry) -> list[ValidationIssue]:
    """Check each source's `package` pin against the manifest committed for it.

    The pin is only checkable because the manifest travels with the repository
    (`mirrors/packages/<source_id>.json`, byte-for-byte the package's own), so
    this reads that file: its sha256 is taken over the bytes on disk, never a
    re-serialisation, and its `package` and `version` must be the pin's. A file
    that exists but is not a valid manifest is SRC005 carrying the
    `ManifestError` text, not a raise: validators report, the build refuses.
    Any `.json` whose stem is no pinned source is SRC006; other entries under
    the directory are TBL009's (`unexpected_mirror_entries`).
    """
    directory = root / "mirrors" / PACKAGES_DIR
    issues: list[ValidationIssue] = []
    pinned: set[str] = set()

    for source in registry.sources:
        pin = source.package
        if pin is None:
            continue
        pinned.add(source.id)
        path = directory / f"{source.id}.json"

        def refuse(message: str, *, path: Path = path) -> None:
            issues.append(ValidationIssue("SRC005", Severity.ERROR, str(path), message))

        try:
            data = path.read_bytes()
        except OSError as exc:
            refuse(f"source '{source.id}' pins a package but its manifest cannot be read: {exc}")
            continue
        digest = hashlib.sha256(data).hexdigest()
        if digest != pin.manifest_sha256:
            refuse(
                f"manifest sha256 is {digest}, but the pin for '{source.id}' says "
                f"{pin.manifest_sha256}"
            )
            continue
        try:
            manifest = load_manifest(path)
        except ManifestError as exc:
            refuse(str(exc))
            continue
        if manifest.package != source.id:
            refuse(f"manifest is for package '{manifest.package}', not '{source.id}'")
        if manifest.version != pin.version:
            refuse(
                f"manifest version is '{manifest.version}', but the pin for "
                f"'{source.id}' says '{pin.version}'"
            )

    # An unreadable directory is reported once, as TBL009, by
    # `unexpected_mirror_entries`; repeating it here would double-count it.
    entries = list_dir(directory, "TBL009")[0] if directory.is_dir() else []
    for entry in entries:
        if entry.is_file() and entry.suffix == ".json" and entry.stem not in pinned:
            issues.append(
                ValidationIssue(
                    "SRC006",
                    Severity.ERROR,
                    str(entry),
                    f"'{entry.name}' is a package manifest but no source in "
                    f"mirrors/sources.yaml pins package '{entry.stem}'",
                )
            )
    return issues


# A source no table's rows are, but one column of another source's rows is
# (plan Decision 8): `audain2021_s8` supplies `n_chd_case_carriers` on the
# `audain2021_s6` package's burden rows. It keeps a hand-written registry entry
# and no package; it is accepted because its `version` sha256 is an input of
# the named package's committed manifest.
COLUMN_ONLY_SOURCES: dict[str, str] = {"audain2021_s8": "audain2021_s6"}


def validate_pins_cover_mirrors(
    root: Path, registry: SourceRegistry, used: Mapping[str, set[str]]
) -> list[ValidationIssue]:
    """Every mirror file is listed by a pinned manifest, and every source a mirror
    relies on is pinned and that source's own manifest lists the table (SRC007,
    SRC008).

    This does NOT verify a row against a package's bytes: packages are not
    committed, so a row filed under the correct pinned source cannot be checked
    here. What it catches is structural.

    SRC007 has two forms. (1) A mirror file (`mirror_paths`) that no pinned
    source's committed manifest lists by `(table, shard)`. (2) A mirror that
    relies on a pinned source (`used`, as `runner._used_sources` returns it)
    whose own manifest does not list that file's `(table, shard)`: a row
    relabelled to some other pinned source is caught by form (2) and not by
    form (1), because another manifest lists the file. A `COLUMN_ONLY_SOURCES`
    entry is held to its owner's manifest. SRC008: a source a mirror relies on
    that has no `package` pin and is not a `COLUMN_ONLY_SOURCES` entry whose
    `version` sha256 AND `url` match one input of its package's manifest. A
    source id the registry does not hold is SRC003's, not reported here.

    Runs only when some source pins a package, so pin-less repositories (most
    unit fixtures) stay valid; `tests/test_committed_mirrors.py` holds the real
    corpus without that gate. A pinned manifest that cannot be read is SRC005's
    (`validate_package_pins`): listing nothing, it would report every file it
    would have listed, so SRC007 is skipped when any pinned manifest is
    unreadable, and a source whose manifest is unreadable is not reported
    either.
    """
    pinned = [source for source in registry.sources if source.package is not None]
    if not pinned:
        return []

    inputs: dict[str, set[tuple[str, str]]] = {}
    listed_by: dict[str, set[tuple[str, str | None]]] = {}
    unreadable = False
    for source in pinned:
        try:
            manifest = load_manifest(root / "mirrors" / PACKAGES_DIR / f"{source.id}.json")
        except ManifestError:
            unreadable = True
            continue
        inputs[source.id] = {(entry.sha256, entry.url) for entry in manifest.inputs}
        listed_by[source.id] = {(table.name, table.shard) for table in manifest.tables}
    listed = set().union(*listed_by.values()) if listed_by else set()

    issues: list[ValidationIssue] = []
    keys: dict[str, tuple[str, str | None]] = {}
    for path, schema_name in mirror_paths(root):
        shard = None if schema_name in FLAT_TABLES else path.stem
        keys[str(path)] = (schema_name, shard)
        if not unreadable and (schema_name, shard) not in listed:
            issues.append(
                ValidationIssue(
                    "SRC007",
                    Severity.ERROR,
                    str(path),
                    f"no pinned package manifest lists table '{schema_name}'"
                    + (f" shard '{shard}'" if shard else "")
                    + "; mirror rows must be ingested from a pinned package",
                )
            )

    for location in sorted(used):
        key = keys.get(location)
        for source_id in sorted(used[location]):
            entry = registry.by_id.get(source_id)
            if entry is None:
                continue
            owner = COLUMN_ONLY_SOURCES.get(source_id)
            if entry.package is not None or owner is not None:
                holder = source_id if entry.package is not None else owner
                if (
                    key is not None
                    and key in listed
                    and holder in listed_by
                    and key not in listed_by[holder]
                ):
                    table, shard = key
                    others = sorted(i for i, tables in listed_by.items() if key in tables)
                    issues.append(
                        ValidationIssue(
                            "SRC007",
                            Severity.ERROR,
                            location,
                            f"this mirror relies on source '{source_id}', but the manifest "
                            f"of package '{holder}' does not list table '{table}'"
                            + (f" shard '{shard}'" if shard else "")
                            + f"; listed by: {', '.join(others)}",
                        )
                    )
            if entry.package is not None:
                continue
            if owner is not None:
                if owner not in inputs and any(s.id == owner for s in pinned):
                    continue  # the owner's manifest is unreadable: SRC005's
                digest = entry.version.removeprefix("sha256:")
                if entry.version.startswith("sha256:") and (digest, entry.url) in inputs.get(
                    owner, set()
                ):
                    continue
            if owner is not None:
                listed_for_url = sorted(
                    sha for sha, url in inputs.get(owner, set()) if url == entry.url
                )
                detail = (
                    f" and whose version is not an input of package '{owner}' with "
                    "the entry's url; "
                    + (
                        f"that manifest lists sha256:{', sha256:'.join(listed_for_url)} for it; "
                        if listed_for_url
                        else "that manifest lists no input at that url; "
                    )
                    + "edit the entry's `version` and `retrieved_on` by hand to match"
                )
            else:
                detail = ""
            issues.append(
                ValidationIssue(
                    "SRC008",
                    Severity.ERROR,
                    location,
                    f"this mirror relies on source '{source_id}', which pins no package" + detail,
                )
            )
    return issues
