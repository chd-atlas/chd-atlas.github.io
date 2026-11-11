# tests/unit/test_identifiers.py
import pytest
from pydantic import TypeAdapter, ValidationError

from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    ContrastId,
    DiseaseId,
    Doi,
    FunctionalId,
    HgncId,
    ModelPhenotypeId,
    ModificationId,
    PhenotypeId,
    Pmcid,
    Pmid,
    SequenceOntologyId,
    TaxonId,
    UniprotAccession,
)


@pytest.mark.parametrize(
    ("type_", "value"),
    [
        (HgncId, "HGNC:11604"),
        (Pmid, "PMID:8988165"),
        (PhenotypeId, "HP:0001631"),
        (PhenotypeId, "MONDO:0007259"),
        (ModelPhenotypeId, "MP:0010402"),
        (ModelPhenotypeId, "ZP:0000706"),
        (ModelPhenotypeId, "HP:0001631"),
        (AssertionId, "CHDA:AST:0000001"),
        (FunctionalId, "CHDA:FUN:0000042"),
        (UniprotAccession, "Q99593"),
        (UniprotAccession, "A0A0B4J2D5"),
        (AccessionId, "PXD012345"),
        (AccessionId, "GSE123456"),
        (AccessionId, "E-MTAB-1234"),
        (TaxonId, "NCBITaxon:9606"),
        (ContrastId, "tof_vs_control"),
        (ContrastId, "a1"),
        (SequenceOntologyId, "SO:0001587"),
        (ModificationId, "MOD:00046"),
        (AccessionId, "EGAS00001000123"),
        (AccessionId, "PXD1234567"),
        (Doi, "10.1038/ng0197-30"),
        (Pmcid, "PMC1234567"),
    ],
)
def test_accepts_well_formed_identifiers(type_: type, value: str) -> None:
    assert TypeAdapter(type_).validate_python(value) == value


@pytest.mark.parametrize(
    ("type_", "value"),
    [
        (HgncId, "11604"),
        (HgncId, "hgnc:11604"),
        (Pmid, "8988165"),
        (PhenotypeId, "HP:1631"),
        (PhenotypeId, "OMIM:142900"),
        (AssertionId, "CHDA:AST:1"),
        (AssertionId, "CHDA:FUN:0000001"),
        (UniprotAccession, "NOTANACC"),
        (AccessionId, "PXD12345"),
        (TaxonId, "9606"),
        (TaxonId, "NCBITaxon:"),
        (ContrastId, "TOF_vs_control"),
        (ContrastId, "_foo"),
        (ContrastId, "foo_"),
        (SequenceOntologyId, "SO:12345"),
        (ModificationId, "MOD:123456"),
        (FunctionalId, "CHDA:AST:0000001"),
        (AccessionId, "EGAS0001"),
        (Doi, "not-a-doi"),
        (Doi, "10.1038"),
        (Pmcid, "1234567"),
    ],
)
def test_rejects_malformed_identifiers(type_: type, value: str) -> None:
    with pytest.raises(ValidationError):
        TypeAdapter(type_).validate_python(value)


def test_disease_id_accepts_mondo_and_rejects_hp() -> None:
    """DiseaseId is narrower than PhenotypeId, which admits both prefixes.

    A disease entity and a phenotypic feature are different claims: MONDO:0007732
    is Holt-Oram syndrome, HP:0001631 is an atrial septal defect. The mirror keys
    on the first and must not silently accept the second.
    """
    adapter = TypeAdapter(DiseaseId)
    assert adapter.validate_python("MONDO:0007732") == "MONDO:0007732"
    for rejected in ("HP:0001631", "MONDO:7732", "MONDO:00077321", "mondo:0007732"):
        with pytest.raises(ValidationError):
            adapter.validate_python(rejected)


def test_the_uniprot_list_pattern_takes_a_pipe_separated_list_and_nothing_looser() -> None:
    import re

    from chd_atlas.identifiers import UNIPROT_LIST_PATTERN

    for good in ("P98175", "P0DW28|P98175", "P98175-2|P0DW28"):
        assert re.match(UNIPROT_LIST_PATTERN, good), good
    for bad in ("P0DW28|", "|P98175", "P0DW28||P98175", "P0DW28|nope", ""):
        assert not re.match(UNIPROT_LIST_PATTERN, bad), bad
