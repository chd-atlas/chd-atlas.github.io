# tests/unit/test_dataset_model.py
import pytest
from pydantic import ValidationError

from chd_atlas.models.dataset import _ARCHIVE_TECHNOLOGIES, Contrast, Dataset, Stage
from chd_atlas.vocab import Archive, Technology


def _contrast(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "tof_vs_control",
        "description": "Tetralogy of Fallot right ventricle versus non-CHD control",
        "case_group": "TOF RV",
        "control_group": "Non-CHD RV",
        "statistical_method": "limma-voom",
        "software": "limma 3.58.1",
    }
    base.update(overrides)
    return base


def _dataset(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "PXD012345",
        "archive": "pride",
        "technology": "phosphoproteomics",
        "design": "contrast",
        "tissue": "right ventricle",
        "developmental_stage": "infant",
        "organism": "NCBITaxon:9606",
        "n_samples": 24,
        "licence": "CC-BY-4.0",
        "contrasts": [_contrast()],
    }
    base.update(overrides)
    return base


def _a_contrast(**overrides: object) -> Contrast:
    base: dict[str, object] = {
        "id": "tof_vs_control",
        "description": "Tetralogy of Fallot right ventricle versus non-CHD control",
        "case_group": "TOF RV",
        "control_group": "Non-CHD RV",
        "statistical_method": "limma-voom",
        "software": "limma 3.58.1",
    }
    base.update(overrides)
    return Contrast(**base)


def _profile_dataset(**overrides: object) -> Dataset:
    base: dict[str, object] = {
        "id": "E-MTAB-6814",
        "archive": "arrayexpress",
        "technology": "bulk_rnaseq",
        "design": "profile",
        "tissue": "Heart",
        "developmental_stage": "embryonic",
        "organism": "NCBITaxon:9606",
        "n_samples": 297,
        "licence": "CC BY 4.0",
        "contrasts": [],
        "cardiac_tissues": ["Heart"],
        "detection_floor": 1.0,
        "floor_source": "source methods, section 4",
        "quantile_estimator": "linear",
        "stages": [
            {"token": "7wpc", "wpc": 7.0, "order": 1},
            {"token": "senior", "wpc": None, "order": 2},
        ],
    }
    base.update(overrides)
    return Dataset(**base)


def test_parses_a_well_formed_dataset() -> None:
    dataset = Dataset.model_validate(_dataset())
    assert dataset.contrasts[0].id == "tof_vs_control"
    assert dataset.publication is None


def test_contrast_id_must_be_snake_case() -> None:
    with pytest.raises(ValidationError):
        Contrast.model_validate(_contrast(id="TOF vs Control"))


def test_dataset_requires_at_least_one_contrast() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(contrasts=[]))


def test_dataset_rejects_duplicate_contrast_ids() -> None:
    with pytest.raises(ValidationError, match="duplicate contrast ids:.*tof_vs_control"):
        Dataset.model_validate(_dataset(contrasts=[_contrast(), _contrast()]))


def test_dataset_requires_a_licence() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(licence=""))


def test_dataset_rejects_an_accession_that_does_not_match_its_archive() -> None:
    with pytest.raises(ValidationError, match="accession GSE123456 does not match archive 'pride'"):
        Dataset.model_validate(_dataset(id="GSE123456"))


def test_archive_technology_map_is_total() -> None:
    """A missing key would raise KeyError instead of a clean ValidationError."""
    assert set(_ARCHIVE_TECHNOLOGIES) == set(Archive)
    covered: set[Technology] = set()
    for technologies in _ARCHIVE_TECHNOLOGIES.values():
        covered |= technologies
    assert covered == set(Technology)


def test_rejects_mass_spec_technology_in_a_sequencing_archive() -> None:
    with pytest.raises(ValidationError, match="phosphoproteomics"):
        Dataset.model_validate(_dataset(id="GSE123456", archive="geo"))


def test_rejects_sequencing_technology_in_pride() -> None:
    with pytest.raises(ValidationError, match="bulk_rnaseq"):
        Dataset.model_validate(_dataset(technology="bulk_rnaseq"))


def test_accepts_sequencing_technology_in_geo() -> None:
    dataset = Dataset.model_validate(
        _dataset(id="GSE123456", archive="geo", technology="bulk_rnaseq")
    )
    assert dataset.technology == "bulk_rnaseq"


def test_dataset_with_contrasts_needs_at_least_two_samples() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(n_samples=1))


def test_rejects_a_self_comparing_contrast() -> None:
    with pytest.raises(ValidationError, match="case_group"):
        Contrast.model_validate(_contrast(case_group="RV", control_group="rv "))


def test_reports_every_duplicate_contrast_id() -> None:
    with pytest.raises(ValidationError) as exc:
        Dataset.model_validate(
            _dataset(
                contrasts=[
                    _contrast(id="a"),
                    _contrast(id="a"),
                    _contrast(id="b"),
                    _contrast(id="b"),
                ]
            )
        )
    message = str(exc.value)
    assert "a" in message
    assert "b" in message


def test_rejects_a_whitespace_only_licence() -> None:
    """`min_length=1` accepts "   ", which records no provenance at all.

    `Source.licence` already rejects it. A dataset redistributes third-party
    measurements exactly as a source does, so the two must agree on what counts
    as a recorded licence.
    """
    with pytest.raises(ValidationError, match="licence must not be blank"):
        Dataset.model_validate(_dataset(licence="   "))


def test_a_profile_dataset_declares_no_contrast_and_a_contrast_dataset_must() -> None:
    """`design` is explicit so an empty list cannot mean two different things.

    Relaxing `contrasts` to "may be empty" would make a curator who simply
    forgot the contrasts indistinguishable from an abundance atlas that
    correctly has none.
    """
    assert _profile_dataset().contrasts == []

    with pytest.raises(ValidationError, match="design 'profile' declares no contrasts"):
        _profile_dataset(contrasts=[_a_contrast()])

    with pytest.raises(ValidationError, match="design 'contrast' requires at least one contrast"):
        _profile_dataset(design="contrast", contrasts=[])


def test_a_profile_dataset_must_name_a_cardiac_tissue_it_actually_has() -> None:
    """D40: the cardiac token is declared, never matched on the string "heart".

    A set, not a scalar: a dataset with two cardiac tokens would otherwise file
    the second as a comparison organ, and tau would measure specificity against
    the heart's own other chamber.
    """
    assert _profile_dataset().cardiac_tissues == ("Heart",)
    assert _profile_dataset(cardiac_tissues=["Heart", "Atrium"]).cardiac_tissues == (
        "Heart",
        "Atrium",
    )

    with pytest.raises(ValidationError, match="cardiac_tissues must not be empty"):
        _profile_dataset(cardiac_tissues=[])

    with pytest.raises(ValidationError, match="duplicate cardiac_tissues"):
        _profile_dataset(cardiac_tissues=["Heart", "Heart"])


@pytest.mark.parametrize("blank_tissues", [[""], ["  "]])
def test_a_blank_cardiac_tissue_is_rejected(blank_tissues: list[str]) -> None:
    """Same reasoning as `licence_is_not_blank`: a blank entry names no tissue,
    and `cardiac_tissues_are_unique` only catches two identical blanks, not one
    stray blank from a YAML slip.
    """
    with pytest.raises(ValidationError, match="cardiac_tissues must not contain a blank"):
        _profile_dataset(cardiac_tissues=blank_tissues)


@pytest.mark.parametrize(
    ("override", "match"),
    [
        pytest.param({"stages": []}, "stages must not be empty", id="stages"),
        pytest.param(
            {"detection_floor": None}, "detection_floor and floor_source", id="detection_floor"
        ),
        pytest.param({"floor_source": ""}, "detection_floor and floor_source", id="floor_source"),
        pytest.param({"quantile_estimator": ""}, "quantile_estimator", id="quantile_estimator"),
    ],
)
def test_a_profile_dataset_must_declare_everything_the_derivation_needs(
    override: dict[str, object], match: str
) -> None:
    """D41 and the fields beside it: each is the input to a published figure
    (see `a_profile_dataset_is_fully_declared`), so a profile dataset missing
    any one of them is rejected rather than silently computed against a
    default nobody curated. Each case changes exactly one field away from the
    otherwise-valid `_profile_dataset()` baseline, so it is unambiguous which
    branch caught it.
    """
    with pytest.raises(ValidationError, match=match):
        _profile_dataset(**override)


def test_stage_tokens_are_unique_and_wpc_may_be_null_postnatally() -> None:
    stages = {stage.token: stage.wpc for stage in _profile_dataset().stages}
    assert stages["7wpc"] == 7.0
    assert stages["senior"] is None

    with pytest.raises(ValidationError, match="duplicate stage tokens"):
        _profile_dataset(
            stages=[
                {"token": "7wpc", "wpc": 7.0, "order": 1},
                {"token": "7wpc", "wpc": 9.0, "order": 2},
            ]
        )


@pytest.mark.parametrize(
    "fields",
    [
        pytest.param({"wpc": 0, "order": 1}, id="wpc of zero"),
        pytest.param({"wpc": -1, "order": 1}, id="negative wpc"),
        pytest.param({"wpc": None}, id="no order at all"),
        pytest.param({"wpc": None, "order": 0}, id="order of zero"),
        pytest.param({"wpc": None, "order": -1}, id="negative order"),
    ],
)
def test_a_stage_refuses_a_position_it_cannot_place(fields: dict[str, object]) -> None:
    """Both bounds and the requiredness, on one baseline.

    The two bounds are spelled differently and this docstring said they were
    one: `wpc` is `gt=0` (a float, so any positive value) and `order` is
    `ge=1` (an int, so 1-based). For an integer the two are the same set, and
    that coincidence is exactly why the wrong name survived -- a reader
    checking the model finds `ge=1` and a claim about `gt=0`.

    Nullability (a postnatal stage carrying no `wpc`) is covered above and is
    a fact about the stage, not a missing value.

    **`order` is required, not defaulted**, and that is the case with a
    reason worth keeping: a stage with no declared position cannot be placed
    on a time axis, and a default would silently put every undeclared stage
    in the same slot -- reintroducing the alphabetical scramble this field
    exists to remove. It was its own one-line test beside two more, all three
    differing only in which field they moved off a valid baseline.
    """
    with pytest.raises(ValidationError):
        Stage(token="neonate", **fields)  # type: ignore[arg-type]


def _cell_type_dataset(**overrides: object) -> Dataset:
    base: dict[str, object] = {
        "id": "hoc",
        "archive": "ucsc_cell_browser",
        "technology": "scrnaseq",
        "design": "cell_type_profile",
        "tissue": "Heart",
        "developmental_stage": "adult",
        "organism": "NCBITaxon:9606",
        "n_samples": 10,
        "licence": "CC BY 4.0",
        "cell_types": [
            {"token": "Cardiomyocyte", "label": "Cardiomyocyte", "compartment": "myocardial"},
            {"token": "Fibroblast", "label": "Fibroblast", "compartment": "stromal"},
        ],
        "regions": ["LV", "RV"],
        "detection_floor": 0.1,
        "floor_source": "source methods, section 2",
        "quantile_estimator": "linear",
    }
    base.update(overrides)
    return Dataset(**base)


def test_a_ucsc_dataset_id_is_its_lowercase_browser_name() -> None:
    """The UCSC Cell Browser names a dataset by a lowercase slug, not an accession."""
    assert _cell_type_dataset().id == "hoc"
    assert _cell_type_dataset(id="heart-atlas").id == "heart-atlas"

    with pytest.raises(ValidationError, match="does not match archive 'ucsc_cell_browser'"):
        _cell_type_dataset(id="GSE1")
    with pytest.raises(ValidationError, match="does not match archive 'geo'"):
        _cell_type_dataset(archive="geo")


@pytest.mark.parametrize(
    ("override", "match"),
    [
        pytest.param({"cell_types": []}, "cell_types must not be empty", id="cell_types"),
        pytest.param({"regions": []}, "regions must not be empty", id="regions"),
        pytest.param(
            {"detection_floor": None}, "detection_floor and floor_source", id="detection_floor"
        ),
        pytest.param({"floor_source": ""}, "detection_floor and floor_source", id="floor_source"),
        pytest.param({"quantile_estimator": ""}, "quantile_estimator", id="quantile_estimator"),
        pytest.param(
            {
                "cell_types": [
                    {"token": "Fibroblast", "label": "Fibroblast", "compartment": "stromal"},
                    {"token": "Fibroblast", "label": "Fibro", "compartment": "stromal"},
                ]
            },
            "duplicate cell_type tokens",
            id="duplicate_tokens",
        ),
        pytest.param(
            {"contrasts": [_contrast()]},
            "design 'cell_type_profile' declares no contrasts",
            id="contrasts",
        ),
        # Measured before the check (review, 2026-10-28): the model accepted
        # ['LV', 'LV', 'all', 'Apex', ' '], and "Apex" would draw a column of
        # dashes captioned "no heart had 10 or more cells there".
        pytest.param({"regions": ["LV", "LV"]}, "duplicate regions", id="duplicate_regions"),
        pytest.param({"regions": ["LV", "all"]}, "'all'", id="pooled_region"),
        pytest.param({"regions": ["LV", "Apex"]}, "'Apex'", id="unknown_region"),
        pytest.param({"regions": ["LV", " "]}, "' '", id="blank_region"),
        *(
            pytest.param(
                {
                    "cell_types": [
                        {"token": "Fibroblast", "label": "Fibroblast", "compartment": "stromal"}
                        | {field: "   "}
                    ]
                },
                f"{field} must not be blank",
                id=f"blank_{field}",
            )
            for field in ("token", "label", "compartment")
        ),
    ],
)
def test_a_cell_type_profile_dataset_is_fully_declared(
    override: dict[str, object], match: str
) -> None:
    """Each case moves exactly one field off the valid baseline, so the branch
    that caught it is unambiguous (same shape as the profile-dataset test).
    """
    assert _cell_type_dataset().regions == ("LV", "RV")
    with pytest.raises(ValidationError, match=match):
        _cell_type_dataset(**override)


def test_a_cell_type_profile_dataset_needs_no_stages() -> None:
    """Stages and cardiac_tissues belong to the bulk profile design; a cell-type
    dataset is indexed by cell type and region instead -- Farah et al. 2024's
    fetal hearts pool 9 to 15 weeks rather than follow a stage series.
    """
    dataset = _cell_type_dataset()
    assert dataset.stages == ()
    assert dataset.cardiac_tissues == ()


@pytest.mark.parametrize("design", ["contrast", "profile"])
@pytest.mark.parametrize(
    "field",
    [
        pytest.param(
            {"cell_types": [{"token": "T", "label": "T", "compartment": "c"}]}, id="cell_types"
        ),
        pytest.param({"regions": ["LV"]}, id="regions"),
    ],
)
def test_cell_types_and_regions_are_refused_on_any_other_design(
    design: str, field: dict[str, object]
) -> None:
    """The `[]`-unless-`cell_type_profile` claim in the model comment and in
    docs/data-api.md is enforced here, not only stated.
    """
    base = _profile_dataset().model_dump() if design == "profile" else _dataset()
    base.update(field)
    with pytest.raises(ValidationError, match=f"design '{design}' declares no cell_types"):
        Dataset.model_validate(base)
