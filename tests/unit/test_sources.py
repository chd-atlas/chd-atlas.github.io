import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.validate.sources import (
    load_sources,
    validate_package_pins,
    validate_pins_cover_mirrors,
)
from tests.unit.package_fixtures import manifest_dict

SOURCES_YAML = """\
sources:
  - id: clinvar
    name: ClinVar
    version: "2026-07"
    retrieved_on: 2026-07-15
    url: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
    licence: public-domain
    redistribution: permitted
  - id: hpo
    name: Human Phenotype Ontology
    version: "2026-06-01"
    retrieved_on: 2026-07-15
    url: https://hpo.jax.org/
    licence: HPO-1.0
    redistribution: permitted
    ontology_prefix: HP
    ontology_file: ontologies/hp-2026-06-01.obo
"""


def _write_sources(root: Path, text: str = SOURCES_YAML) -> None:
    (root / "mirrors").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "sources.yaml").write_text(text)


def test_loads_the_source_registry(tmp_path: Path) -> None:
    _write_sources(tmp_path)

    registry, issues = load_sources(tmp_path)

    assert issues == []
    assert registry.by_id["clinvar"].licence == "public-domain"
    assert registry.ontology_files() == {"HP": "ontologies/hp-2026-06-01.obo"}


def test_missing_licence_is_an_error(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML.replace("    licence: public-domain\n", ""))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]
    assert "licence" in issues[0].message


def test_missing_registry_file_is_an_error(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir()

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC001"]


def test_duplicate_source_ids_are_rejected(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML + SOURCES_YAML.split("sources:\n")[1])

    _, issues = load_sources(tmp_path)

    assert any(i.code == "SRC002" and "duplicate source" in i.message for i in issues)


def test_non_utf8_bytes_are_reported_rather_than_raised(tmp_path: Path) -> None:
    """A Latin-1 byte must not abort the whole validation run."""
    (tmp_path / "mirrors").mkdir(parents=True, exist_ok=True)
    (tmp_path / "mirrors" / "sources.yaml").write_bytes(
        SOURCES_YAML.replace("ClinVar", "Clin\xe9Var").encode("latin-1")
    )

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["YAML001"]


def test_reports_every_duplicate_source_id(tmp_path: Path) -> None:
    _write_sources(
        tmp_path,
        SOURCES_YAML + SOURCES_YAML.split("sources:\n")[1] + SOURCES_YAML.split("sources:\n")[1],
    )

    _, issues = load_sources(tmp_path)

    message = " ".join(i.message for i in issues)
    assert "clinvar" in message
    assert "hpo" in message


def test_unparseable_yaml_is_reported(tmp_path: Path) -> None:
    _write_sources(tmp_path, "sources: [ unclosed\n")

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["YAML001"]


def test_ontology_prefix_requires_a_pinned_file(tmp_path: Path) -> None:
    _write_sources(
        tmp_path,
        SOURCES_YAML.replace("    ontology_file: ontologies/hp-2026-06-01.obo\n", ""),
    )

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]
    assert "ontology_file" in issues[0].message


def test_out_of_range_date_is_reported_rather_than_raised(tmp_path: Path) -> None:
    """Unquoted dates are the house style; ruamel raises a bare ValueError on month 13."""
    _write_sources(tmp_path, SOURCES_YAML.replace("2026-07-15", "2026-13-45"))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["YAML001"]


def test_schema_violation_still_reports_src002_not_yaml001(tmp_path: Path) -> None:
    """ValidationError is a ValueError; it must not be swallowed by the YAML guard."""
    _write_sources(tmp_path, SOURCES_YAML.replace("    licence: public-domain\n", ""))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]


def test_duplicate_ontology_prefix_is_rejected(tmp_path: Path) -> None:
    extra = (
        "  - id: hpo_old\n"
        "    name: Human Phenotype Ontology\n"
        '    version: "2020-01-01"\n'
        "    retrieved_on: 2026-07-15\n"
        "    url: https://hpo.jax.org/\n"
        "    licence: HPO-1.0\n"
        "    redistribution: permitted\n"
        "    ontology_prefix: HP\n"
        "    ontology_file: ontologies/hp-2020-01-01.obo\n"
    )
    _write_sources(tmp_path, SOURCES_YAML + extra)

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]
    assert "ontology prefix" in issues[0].message


def test_ontology_file_without_a_prefix_is_rejected(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML.replace("    ontology_prefix: HP\n", ""))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]


def test_absolute_ontology_file_is_rejected(tmp_path: Path) -> None:
    _write_sources(
        tmp_path,
        SOURCES_YAML.replace("ontologies/hp-2026-06-01.obo", "/etc/passwd"),
    )

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]


def test_whitespace_licence_is_rejected(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML.replace("licence: public-domain", 'licence: "   "'))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]


def test_source_referenced_by_a_record_must_have_a_registry_entry(tmp_path: Path) -> None:
    from chd_atlas.validate.sources import validate_source_references

    _write_sources(tmp_path)
    registry, _ = load_sources(tmp_path)

    issues = validate_source_references({"mirrors/variants/12.tsv": {"gnomad"}}, registry)

    assert [i.code for i in issues] == ["SRC003"]
    assert "licence unknown" in issues[0].message


def test_registered_sources_pass_the_reference_check(tmp_path: Path) -> None:
    from chd_atlas.validate.sources import validate_source_references

    _write_sources(tmp_path)
    registry, _ = load_sources(tmp_path)

    assert validate_source_references({"mirrors/variants/12.tsv": {"clinvar"}}, registry) == []


def test_prohibited_source_used_by_a_mirror_is_reported(tmp_path: Path) -> None:
    from chd_atlas.validate.sources import validate_source_references

    _write_sources(
        tmp_path,
        SOURCES_YAML.replace("redistribution: permitted", "redistribution: prohibited", 1),
    )
    registry, _ = load_sources(tmp_path)

    issues = validate_source_references({"mirrors/variants/12.tsv": {"clinvar"}}, registry)

    assert [i.code for i in issues] == ["SRC004"]
    assert "clinvar" in issues[0].message


@pytest.mark.parametrize(
    ("fields", "accepted"),
    [
        pytest.param('    file_version: "1.3.42"\n    file_date: "2026-06-23"\n', True, id="both"),
        pytest.param('    file_version: "1.3.42"\n', False, id="version without date"),
        pytest.param('    file_date: "2026-06-23"\n', False, id="date without version"),
        pytest.param('    file_version: "  "\n    file_date: "2026-06-23"\n', False, id="blank"),
    ],
)
def test_a_files_own_version_and_date_are_declared_together(
    tmp_path: Path, fields: str, accepted: bool
) -> None:
    """R11: a file's header version cites as little without its date as the reverse.

    `file_version` and `file_date` are what an Orphadata file states in its
    `<JDBOR>` header, written by the converter as a pair; one alone, or one
    that is blank, is a malformed entry reported as SRC002 rather than
    published as half a citation.
    """
    _write_sources(
        tmp_path,
        SOURCES_YAML.replace(
            "    licence: public-domain\n", fields + "    licence: public-domain\n"
        ),
    )

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ([] if accepted else ["SRC002"])


@pytest.mark.parametrize(
    ("digest", "accepted"),
    [
        pytest.param("a" * 64, True, id="64 hex"),
        pytest.param("a" * 63, False, id="63 hex"),
        pytest.param("A" * 64, False, id="uppercase"),
    ],
)
def test_a_malformed_pin_is_refused(tmp_path: Path, digest: str, accepted: bool) -> None:
    """A package pin names a manifest by its sha256; anything else is SRC002."""
    pin = f'    package:\n      version: "2026.1"\n      manifest_sha256: "{digest}"\n'
    _write_sources(
        tmp_path,
        SOURCES_YAML.replace("    licence: public-domain\n", pin + "    licence: public-domain\n"),
    )

    registry, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ([] if accepted else ["SRC002"])
    if accepted:
        assert registry.by_id["clinvar"].package is not None
        assert registry.by_id["clinvar"].package.manifest_sha256 == digest


# --- package pins agree with their committed manifests (feeds phase 1, task 4) ---


def _pinned_repo(
    root: Path, *, manifest: dict[str, object] | None = None, raw: bytes | None = None
) -> Path:
    """A repo whose `clinvar` source pins package `toy_burden`; returns the manifest path.

    The pin is computed from the bytes written, so a test changes exactly one thing.
    """
    data = raw if raw is not None else json.dumps(manifest or manifest_dict()).encode()
    pin = (
        "    package:\n"
        f'      version: "{manifest_dict()["version"]}"\n'
        f'      manifest_sha256: "{hashlib.sha256(data).hexdigest()}"\n'
    )
    text = SOURCES_YAML.replace("  - id: clinvar\n", "  - id: toy_burden\n", 1).replace(
        "    redistribution: permitted\n", "    redistribution: permitted\n" + pin, 1
    )
    _write_sources(root, text)
    path = root / "mirrors" / "packages" / "toy_burden.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(data)
    return path


def _pin_issues(root: Path) -> list[tuple[str, str]]:
    registry, load_issues = load_sources(root)
    assert load_issues == []
    return [(i.code, i.location) for i in validate_package_pins(root, registry)]


@pytest.mark.parametrize(
    ("damage", "fragment"),
    [
        ({"missing": True}, "manifest cannot be read"),
        ({"manifest": manifest_dict(package="other_package")}, "'other_package'"),
        ({"manifest": manifest_dict(version="2026-11-11.1")}, "'2026-11-11.1'"),
        ({"raw": b"{not json"}, "cannot read manifest"),
    ],
    ids=["missing", "another-package", "another-version", "not-a-manifest-is-not-a-raise"],
)
def test_a_pin_its_committed_manifest_does_not_honour_is_src005(
    tmp_path: Path, damage: dict[str, Any], fragment: str
) -> None:
    path = _pinned_repo(tmp_path, **{k: v for k, v in damage.items() if k != "missing"})
    if damage.get("missing"):
        path.unlink()

    registry, _ = load_sources(tmp_path)
    issues = validate_package_pins(tmp_path, registry)

    assert [(i.code, i.location) for i in issues] == [("SRC005", str(path))]
    assert fragment in issues[0].message


def test_a_pin_whose_manifest_bytes_changed_is_src005(tmp_path: Path) -> None:
    """The same manifest re-indented: it still parses and still agrees with the pin's
    package and version, so only the sha256 over the bytes on disk can refuse it."""
    path = _pinned_repo(tmp_path)
    assert _pin_issues(tmp_path) == []
    pinned = hashlib.sha256(path.read_bytes()).hexdigest()
    path.write_text(json.dumps(manifest_dict(), indent=2))
    on_disk = hashlib.sha256(path.read_bytes()).hexdigest()
    assert on_disk != pinned

    registry, _ = load_sources(tmp_path)
    issues = validate_package_pins(tmp_path, registry)

    assert [(i.code, i.location) for i in issues] == [("SRC005", str(path))]
    assert issues[0].message == (
        f"manifest sha256 is {on_disk}, but the pin for 'toy_burden' says {pinned}"
    )


def test_a_manifest_no_pin_claims_is_src006(tmp_path: Path) -> None:
    _write_sources(tmp_path)
    assert _pin_issues(tmp_path) == [], "no pins and no packages directory"
    stray = tmp_path / "mirrors" / "packages" / "orphan_pkg.json"
    stray.parent.mkdir(parents=True)
    stray.write_text("{}")
    (stray.parent / "notes.txt").write_text("ignored here; TBL009 owns it")

    assert _pin_issues(tmp_path) == [("SRC006", str(stray))]


# --- every mirror row comes from a pinned package (feeds phase 3, task 18) ---


def _cover(root: Path, used: dict[str, set[str]]) -> list[tuple[str, str]]:
    registry, load_issues = load_sources(root)
    assert load_issues == []
    return [(i.code, i.location) for i in validate_pins_cover_mirrors(root, registry, used)]


def _touch(root: Path, relative: str) -> Path:
    path = root / "mirrors" / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x\n")
    return path


def test_a_mirror_file_no_manifest_lists_is_src007(tmp_path: Path) -> None:
    """`toy_burden`'s manifest lists `burden` (flat); the genes file and a
    `profiles` shard are listed by nobody. A listed table in a different shard
    does not cover another shard of the same table."""
    _pinned_repo(tmp_path)
    burden = _touch(tmp_path, "burden.tsv")
    genes = _touch(tmp_path, "genes.tsv")
    shard = _touch(tmp_path, "profiles/E-MTAB-6814.tsv")

    issues = _cover(tmp_path, {str(burden): {"toy_burden"}})

    assert sorted(issues) == sorted([("SRC007", str(genes)), ("SRC007", str(shard))])


def test_a_listed_shard_is_covered_by_table_and_shard_together(tmp_path: Path) -> None:
    manifest = manifest_dict()
    manifest["tables"][0].update(name="profiles", file="profiles.tsv", shard="E-MTAB-6814")
    _pinned_repo(tmp_path, manifest=manifest)
    listed = _touch(tmp_path, "profiles/E-MTAB-6814.tsv")
    other = _touch(tmp_path, "profiles/hoc.tsv")

    assert _cover(tmp_path, {}) == [("SRC007", str(other))]
    assert listed.exists()


def test_a_source_a_mirror_relies_on_without_a_pin_is_src008(tmp_path: Path) -> None:
    _pinned_repo(tmp_path)
    burden = _touch(tmp_path, "burden.tsv")

    issues = _cover(tmp_path, {str(burden): {"toy_burden", "hpo", "not_in_registry"}})

    # `hpo` is registered and unpinned; an unregistered id is SRC003's, not ours.
    assert issues == [("SRC008", str(burden))]


_S8_SHA = "d" * 64
_S6_SOURCES = """\
sources:
  - id: audain2021_s6
    name: S6
    version: "sha256:{s6}"
    retrieved_on: 2026-08-05
    url: https://x.org/s6
    licence: CC0-1.0
    redistribution: permitted
    package:
      version: "2026-11-10.1"
      manifest_sha256: "{pin}"
  - id: audain2021_s8
    name: S8
    version: "sha256:{s8}"
    retrieved_on: 2026-08-05
    url: https://x.org/s8
    licence: CC0-1.0
    redistribution: permitted
"""


def _s8_repo(root: Path, *, s8_version_sha: str, in_package: str = "audain2021_s6") -> Path:
    """`audain2021_s6` pins a manifest whose second input is `_S8_SHA`; the S8
    source's `version` carries `s8_version_sha`. A second package, `other`,
    also lists `_S8_SHA` when `in_package` is not s6 (the exemption must not
    accept another package's inputs)."""
    manifest = manifest_dict(package="audain2021_s6", version="2026-11-10.1")
    manifest["inputs"].append(dict(manifest["inputs"][0], name="s013.xlsx", sha256=_S8_SHA))
    if in_package != "audain2021_s6":
        manifest["inputs"][1]["sha256"] = "e" * 64
    data = json.dumps(manifest).encode()
    _write_sources(
        root,
        _S6_SOURCES.format(s6="c" * 64, s8=s8_version_sha, pin=hashlib.sha256(data).hexdigest()),
    )
    packages = root / "mirrors" / "packages"
    packages.mkdir(parents=True, exist_ok=True)
    (packages / "audain2021_s6.json").write_bytes(data)
    return _touch(root, "burden.tsv")


def test_a_column_source_is_accepted_only_as_an_input_of_its_rows_package(
    tmp_path: Path,
) -> None:
    burden = _s8_repo(tmp_path, s8_version_sha=_S8_SHA)
    used = {str(burden): {"audain2021_s6", "audain2021_s8"}}
    assert _cover(tmp_path, used) == []

    # One hex digit off: no longer an input of the package.
    one_off = tmp_path / "off"
    burden = _s8_repo(one_off, s8_version_sha="d" * 63 + "e")
    assert _cover(one_off, {str(burden): {"audain2021_s6", "audain2021_s8"}}) == [
        ("SRC008", str(burden))
    ]

    # The right hex, but the package does not list it (the exemption names a
    # package's own inputs, not "any sha anywhere").
    elsewhere = tmp_path / "elsewhere"
    burden = _s8_repo(elsewhere, s8_version_sha=_S8_SHA, in_package="other")
    assert _cover(elsewhere, {str(burden): {"audain2021_s6", "audain2021_s8"}}) == [
        ("SRC008", str(burden))
    ]


def test_a_repository_with_no_pin_is_not_checked(tmp_path: Path) -> None:
    """Decision 9's gate: many unit fixtures build pin-less repositories."""
    _write_sources(tmp_path)
    burden = _touch(tmp_path, "burden.tsv")

    assert _cover(tmp_path, {str(burden): {"clinvar"}}) == []


def test_an_unreadable_manifest_adds_nothing_beyond_its_src005(tmp_path: Path) -> None:
    """SRC005 already names the broken manifest; every file it would have listed
    must not also report SRC007."""
    path = _pinned_repo(tmp_path, raw=b"{not json")
    burden = _touch(tmp_path, "burden.tsv")
    assert path.exists()

    assert _cover(tmp_path, {str(burden): {"toy_burden"}}) == []


def test_a_column_source_whose_package_manifest_is_unreadable_is_not_reported(
    tmp_path: Path,
) -> None:
    burden = _s8_repo(tmp_path, s8_version_sha=_S8_SHA)
    (tmp_path / "mirrors" / "packages" / "audain2021_s6.json").write_bytes(b"{not json")

    assert _cover(tmp_path, {str(burden): {"audain2021_s6", "audain2021_s8"}}) == []
