# tests/unit/test_genes.py
"""Symbol -> HGNC id: resolved, or reported. Never guessed, never dropped. Issue #33."""

from __future__ import annotations

import polars as pl
import pytest

from chd_atlas.build.identity import GeneIdentity
from chd_atlas.genes import GeneRegistry, Resolution
from chd_atlas.issues import Severity
from chd_atlas.validate.genes import (
    validate_identity_ids,
    validate_mirror_symbols,
    validate_one_uniprot_per_gene,
    validate_symbol_resolution,
)

# Three genes, and no two of them carry the same combination of columns, so a
# lookup cannot pass by accident: `POPDC1` has both aliases and a retired name,
# `TBX5` an alias and no retired name, `GATA4` neither. A fixture whose rows all
# populate the same columns cannot distinguish "read the alias index" from "read
# the previous index" (CLAUDE.md 4.36).
#
# The values are the real registry's, and `POPDC1` is the whole issue in one
# row: ClinGen labels HGNC:1152 `BVES`, which is **not** among its aliases
# (`HBVES` is a different string) and **is** its previous symbol.
_ROWS = (
    {
        "hgnc_id": "HGNC:1152",
        "symbol": "POPDC1",
        "aliases": "POP1|HBVES",
        "prev_symbols": "BVES",
    },
    {"hgnc_id": "HGNC:11604", "symbol": "TBX5", "aliases": "HOS", "prev_symbols": ""},
    {"hgnc_id": "HGNC:4173", "symbol": "GATA4", "aliases": "", "prev_symbols": ""},
)


def test_a_symbol_resolves_by_approved_name_then_alias_then_retired_name() -> None:
    """All three tiers, each reporting which one answered.

    The `Resolution` is not decoration. A caller that only got the id back could
    not tell a current symbol from a retired one, and those want different
    treatment when a source disagrees with the registry about a gene's name --
    a source still writing `BVES` is using a name HGNC withdrew, which is worth
    knowing even though the row resolves.

    Each tier is asserted on the gene that *only* that tier can answer for, so
    a mutant collapsing two indexes into one fails here rather than passing on
    a row both would have resolved.
    """
    registry = GeneRegistry(_ROWS)
    assert len(registry) == 3

    approved = registry.resolve("TBX5")
    assert approved.resolution is Resolution.APPROVED
    assert approved.gene == "HGNC:11604"
    assert approved.ok is True

    alias = registry.resolve("HOS")
    assert alias.resolution is Resolution.ALIAS
    assert alias.gene == "HGNC:11604"

    retired = registry.resolve("BVES")
    assert retired.resolution is Resolution.PREVIOUS
    assert retired.gene == "HGNC:1152"
    assert retired.ok is True
    assert retired.candidates == ("HGNC:1152",)

    # Surrounding whitespace is a formatting artefact of a source file, not a
    # different gene; a bare dict lookup would miss it and drop the row.
    assert registry.resolve("  TBX5  ").gene == "HGNC:11604"


def test_an_unknown_symbol_is_reported_and_never_resolved_to_a_guess() -> None:
    """The load-bearing rule: a symbol that maps to nothing must not vanish.

    **`BVES` was this test's subject and is no longer**, which is worth
    recording rather than quietly rewriting. From 2026-08-06 to 2026-08-24 it
    was the live case: ClinGen labels HGNC:1152 `BVES`, the mirror called it
    `POPDC1` with aliases `POP1|HBVES`, and `HBVES` is a *different string*, so
    a symbol-keyed join on `BVES` resolved to nothing. It resolves now -- via
    the previous tier, asserted above -- because the mirror gained the
    `prev_symbols` column that issue #33 asked for and the resolver gained a
    tier to read it with.

    The rule it guarded is unchanged and still needs a subject, so the subject
    is a string no tier claims. A gene whose symbol resolves to nothing and is
    then dropped leaves a green build with the gene simply absent, which is this
    project's characteristic failure.
    """
    registry = GeneRegistry(_ROWS)

    unknown = registry.resolve("NOSUCHSYMBOL")
    assert unknown.resolution is Resolution.UNKNOWN
    assert unknown.gene is None
    assert unknown.ok is False
    assert unknown.candidates == ()

    assert registry.resolve("").resolution is Resolution.UNKNOWN
    assert registry.resolve("   ").resolution is Resolution.UNKNOWN


def test_a_symbol_two_genes_claim_is_refused_rather_than_picked() -> None:
    """Ambiguity names every candidate and resolves to none of them.

    The retired Audain 2026 converter declined to match aliases at all, and its
    stated reason was this exact hazard: "an alias that is another gene's current symbol
    would silently attribute one gene's burden to another". Measured 2026-08-14
    the collision does not occur in the committed registry -- 154 approved
    symbols, 309 aliases, zero overlap -- so only a fixture can reach it, and it
    is guarded structurally rather than on the bet that the next upstream
    refresh keeps that true.

    **Approved beats alias, and that is not the same case.** HGNC retires a
    symbol by making it a *previous* symbol of the gene that kept the name, so a
    string that is one gene's approved symbol and another's alias names the
    first gene today -- resolving to the alias would attribute the row to a gene
    the source did not mean. That must resolve, not report.
    """
    shared = GeneRegistry(
        (
            {"hgnc_id": "HGNC:1", "symbol": "AAA", "aliases": "SHARED"},
            {"hgnc_id": "HGNC:2", "symbol": "BBB", "aliases": "SHARED"},
        )
    )
    ambiguous = shared.resolve("SHARED")
    assert ambiguous.resolution is Resolution.AMBIGUOUS
    assert ambiguous.gene is None
    assert ambiguous.candidates == ("HGNC:1", "HGNC:2")

    # One gene's approved symbol is another's alias: approved wins outright.
    contested = GeneRegistry(
        (
            {"hgnc_id": "HGNC:1", "symbol": "NAME", "aliases": ""},
            {"hgnc_id": "HGNC:2", "symbol": "OTHER", "aliases": "NAME"},
        )
    )
    resolved = contested.resolve("NAME")
    assert resolved.resolution is Resolution.APPROVED
    assert resolved.gene == "HGNC:1"


def test_a_name_the_alias_and_previous_tiers_give_to_different_genes_is_refused() -> None:
    """Owner decision F27 (2026-11-13): the alias and previous tiers answer
    together, so a name they give to more than one gene is refused, exactly as
    a name two genes retired is.

    **`ODD` is the case that decided it.** It is a current alias of GJA1
    (HGNC:4274) and the withdrawn name of OSR1 (HGNC:8111), and the only string
    the committed registry's alias and previous tiers both claim (re-measured
    2026-11-13 over its 185 genes). Until F27 the live alias won. Over the
    whole HGNC release the pins use (2026-10-06) 356 names sit in both tiers,
    and neither order is right for them: the old name means the previous gene
    for `MLL2` (KMT2D) and `HRAS1` (HRAS), the alias gene or nobody for `CBP`
    (CREBBP) and `SIP1` (ZEB2). Refusing names the collision in the ingest
    report, where an order would land the rows on a gene silently.

    **`CBP` is the case a reordering reaches.** Two genes claim it as an alias
    (CREBBP, PAG1) and a third retired it (OPN1LW). Every claimant is a
    candidate: stopping at the alias tier hides OPN1LW, and reading the previous
    tier first answers OPN1LW, a colour-vision gene, for a name that usually
    means CREBBP.

    **`HIGM1` is the case that must still resolve.** It is CD40LG's alias *and*
    its withdrawn name, so one gene claims it: it resolves, labelled an alias,
    so `previous_symbols` leaves it out.

    **`MADH7` is unchanged.** SMAD6 (HGNC:6772) and SMAD7 (HGNC:6773) both
    retired it, so it is ambiguous *within* the previous tier.
    """
    odd = GeneRegistry(
        (
            {"hgnc_id": "HGNC:4274", "symbol": "GJA1", "aliases": "ODD", "prev_symbols": ""},
            {"hgnc_id": "HGNC:8111", "symbol": "OSR1", "aliases": "", "prev_symbols": "ODD"},
        )
    ).resolve("ODD")
    assert odd.resolution is Resolution.AMBIGUOUS
    assert odd.gene is None
    assert odd.candidates == ("HGNC:4274", "HGNC:8111")

    cbp = GeneRegistry(
        (
            {"hgnc_id": "HGNC:2348", "symbol": "CREBBP", "aliases": "CBP", "prev_symbols": ""},
            {"hgnc_id": "HGNC:30043", "symbol": "PAG1", "aliases": "CBP", "prev_symbols": ""},
            {"hgnc_id": "HGNC:9936", "symbol": "OPN1LW", "aliases": "", "prev_symbols": "CBP"},
        )
    ).resolve("CBP")
    assert cbp.resolution is Resolution.AMBIGUOUS
    assert cbp.gene is None
    assert cbp.candidates == ("HGNC:2348", "HGNC:30043", "HGNC:9936")

    one_gene = GeneRegistry(
        (
            {
                "hgnc_id": "HGNC:11935",
                "symbol": "CD40LG",
                "aliases": "CD40L|HIGM1",
                "prev_symbols": "HIGM1|TNFSF5",
            },
        )
    )
    higm1 = one_gene.resolve("HIGM1")
    assert higm1.resolution is Resolution.ALIAS
    assert higm1.gene == "HGNC:11935"
    assert higm1.candidates == ("HGNC:11935",)
    assert one_gene.previous_symbols("HGNC:11935") == ("TNFSF5",)

    madh7 = GeneRegistry(
        (
            {"hgnc_id": "HGNC:6772", "symbol": "SMAD6", "aliases": "", "prev_symbols": "MADH7"},
            {"hgnc_id": "HGNC:6773", "symbol": "SMAD7", "aliases": "", "prev_symbols": "MADH7"},
        )
    ).resolve("MADH7")
    assert madh7.resolution is Resolution.AMBIGUOUS
    assert madh7.gene is None
    assert madh7.candidates == ("HGNC:6772", "HGNC:6773")


def test_a_row_with_no_id_is_skipped_without_taking_the_registry_with_it() -> None:
    """A blank `hgnc_id` has nothing to key on; every other row still loads."""
    registry = GeneRegistry(({"hgnc_id": "  ", "symbol": "GHOST", "aliases": ""}, *_ROWS))
    assert len(registry) == 3
    assert registry.resolve("GHOST").resolution is Resolution.UNKNOWN
    assert registry.resolve("TBX5").gene == "HGNC:11604"


def test_gen003_reports_a_mirror_symbol_the_registry_cannot_resolve() -> None:
    """The check that fires on the committed data, and what it must not report.

    Four rows, one finding.

    - `GATA4`/`GATA4-OLD` is a **genuinely stale registry**: a name no tier
      claims. That is the failure mode GEN003 is for now that a retired symbol
      resolves on its own.
    - `POPDC1`/`BVES` **must stay silent**, and it is the case that changed.
      Until 2026-08-24 this was the finding -- the only one on the committed
      corpus, and the whole standing warning on the validate baseline. `BVES`
      now resolves through the previous tier, so a mirror using it is
      consistent, not stale. Kept as a row here precisely so a mutant deleting
      the previous tier re-reports it and fails.
    - `TBX5` agrees outright and must stay silent.
    - HGNC:9999 is a gene the registry does not carry at all, which is **not** a
      naming disagreement -- `mirrors/genes.tsv` is deliberately narrower than
      the validity mirrors, and reporting those would bury the real finding
      under thousands.
    """
    issues = validate_mirror_symbols(
        (
            ("HGNC:4173", "GATA4-OLD", "mirrors/clingen_gene_validity.tsv"),
            ("HGNC:1152", "BVES", "mirrors/clingen_gene_validity.tsv"),
            ("HGNC:11604", "TBX5", "mirrors/clingen_gene_validity.tsv"),
            ("HGNC:9999", "NOTINREGISTRY", "mirrors/gencc_submissions.tsv"),
        ),
        GeneRegistry(_ROWS),
        "mirrors/genes.tsv",
    )

    assert [issue.code for issue in issues] == ["GEN003"]
    assert issues[0].severity is Severity.WARNING
    assert "HGNC:4173" in issues[0].message
    assert "'GATA4-OLD'" in issues[0].message
    assert (
        "pin a newer HGNC package from chd-atlas/feeds "
        "(`chd-atlas ingest hgnc --package <its directory>`)"
    ) in issues[0].message
    assert "re-ingesting the pinned version changes nothing" in issues[0].message


def test_losing_the_prev_symbols_column_brings_gen003_back_rather_than_going_quiet() -> None:
    """GEN003 is the *second* detector for a vanished column, and independent.

    **The first is `TBL001`, and this docstring claimed otherwise for one
    revision.** It said a mirror without `prev_symbols` "loads with zero
    issues", which was measured against `read_table` -- the wrong entry point.
    `validate_table` is what the repository actually runs, and it reports
    `TBL001 missing required column 'prev_symbols'` (measured 2026-08-24). A
    nullable column may hold nulls; it may not be absent. So the column cannot
    disappear quietly, and the test that would have shipped asserting it can
    was wrong about its own subject.

    What survives the correction is worth keeping, because the two guards fail
    for different reasons and neither implies the other. `TBL001` reads the
    *file* and fires on a header. GEN003 reads the *resolutions* and fires on a
    consequence -- a symbol our own validity mirrors print that the registry
    can no longer answer for. A `GeneRegistry` built from rows that never went
    through `validate_table` (a bypassed gate, which CLAUDE.md section 9 keeps
    guards for deliberately) reaches the second and not the first.

    The demonstration is the nicest available one: the warning that cleared
    when this column landed is the warning that comes back if it leaves.
    Asserted in both directions so a mutant reporting GEN003 unconditionally
    fails the first half.
    """
    stripped = GeneRegistry(
        tuple({key: value for key, value in row.items() if key != "prev_symbols"} for row in _ROWS)
    )
    assert stripped.resolve("BVES").resolution is Resolution.UNKNOWN

    labelled = (("HGNC:1152", "BVES", "mirrors/clingen_gene_validity.tsv"),)
    assert validate_mirror_symbols(labelled, GeneRegistry(_ROWS), "mirrors/genes.tsv") == []

    reported = validate_mirror_symbols(labelled, stripped, "mirrors/genes.tsv")
    assert [issue.code for issue in reported] == ["GEN003"]
    assert "'BVES'" in reported[0].message


def test_a_mirror_symbol_naming_a_different_gene_is_reported_as_such() -> None:
    """The serious case: the symbol resolves, but to somebody else.

    Unreachable from a correctly generated mirror -- both columns come from one
    upstream row -- so this exists to stop a hand-edited mirror introducing it
    undetected, and the message must name the gene it actually resolved to
    rather than saying only that something is wrong.
    """
    issues = validate_mirror_symbols(
        (("HGNC:1152", "TBX5", "mirrors/clingen_gene_validity.tsv"),),
        GeneRegistry(_ROWS),
        "mirrors/genes.tsv",
    )
    assert [issue.code for issue in issues] == ["GEN003"]
    assert "resolves to HGNC:11604" in issues[0].message


def test_an_unreadable_registry_reports_one_skip_and_not_one_per_symbol() -> None:
    """GEN000, the REF000/SRC000/ONT000/SCP000 precedent.

    A missing reference table must say "checks skipped" once. Reporting every
    symbol as unresolvable would bury the real error -- the registry itself --
    under one warning per row, with nothing naming the cause. Both entry points
    take the same branch, because a caller that got a clean report from one and
    a cascade from the other would trust the wrong one.
    """
    skipped = validate_mirror_symbols(
        (("HGNC:1152", "BVES", "x"), ("HGNC:11604", "TBX5", "y")), None, "mirrors/genes.tsv"
    )
    assert [issue.code for issue in skipped] == ["GEN000"]
    assert skipped[0].severity is Severity.WARNING
    assert "skipped" in skipped[0].message

    resolution = validate_symbol_resolution(
        (("BVES", "x"), ("HOS", "y")), None, "mirrors/genes.tsv"
    )
    assert [issue.code for issue in resolution] == ["GEN000"]


def test_gen001_and_gen002_report_an_unmapped_and_an_ambiguous_symbol() -> None:
    """The entry point a symbol-keyed source will use. No caller in the build yet.

    **Stated plainly because it matters:** nothing in `build/` calls
    `validate_symbol_resolution` today, since every source currently mirrored --
    ClinGen and GenCC -- is already keyed on HGNC id. It exists so the machinery
    is in place before the first symbol-only source lands (CHDgene is the likely
    one, #31), and it is exercised here rather than left as a guard nobody has
    watched fail.

    A resolvable symbol produces no issue at all: a validator that reported
    every row would be noise, and the report has to be readable to be read.
    """
    issues = validate_symbol_resolution(
        (
            ("NOSUCHSYMBOL", "chdgene.csv:row 4"),
            ("SHARED", "chdgene.csv:row 9"),
            ("TBX5", "ok"),
            # Resolvable only through the previous tier, and therefore silent.
            # A source still writing a withdrawn name is not an unmapped row.
            ("BVES", "chdgene.csv:row 12"),
        ),
        GeneRegistry(
            (
                *_ROWS,
                {"hgnc_id": "HGNC:1", "symbol": "AAA", "aliases": "SHARED"},
                {"hgnc_id": "HGNC:2", "symbol": "BBB", "aliases": "SHARED"},
            )
        ),
        "mirrors/genes.tsv",
    )

    assert sorted(issue.code for issue in issues) == ["GEN001", "GEN002"]
    unmapped = next(issue for issue in issues if issue.code == "GEN001")
    ambiguous = next(issue for issue in issues if issue.code == "GEN002")
    assert unmapped.location == "chdgene.csv:row 4"
    assert "would be dropped" in unmapped.message
    assert "HGNC:1, HGNC:2" in ambiguous.message
    assert "refusing to guess" in ambiguous.message


@pytest.mark.parametrize(
    ("identity", "problem"),
    [
        pytest.param(
            GeneIdentity(ensembl_gene="ENSG0089225", ncbi_gene=6910),
            "'ENSG0089225' is not an Ensembl gene id (ENSG + 11 digits)",
            id="malformed Ensembl",
        ),
        pytest.param(
            GeneIdentity(ensembl_gene="ENSG00000089225", ncbi_gene=0),
            "0 is not an NCBI Gene id (a positive integer)",
            id="NCBI Gene id below 1",
        ),
        pytest.param(
            GeneIdentity(ensembl_gene="ENSG00000089225"), "no NCBI Gene id", id="missing NCBI"
        ),
        pytest.param(None, "no Ensembl gene id or NCBI Gene id", id="no registry row"),
        pytest.param(GeneIdentity(ensembl_gene="ENSG00000089225", ncbi_gene=6910), None, id="ok"),
    ],
)
def test_a_published_genes_identity_ids_are_held_to_the_rule_the_build_refuses_on(
    identity: GeneIdentity | None, problem: str | None
) -> None:
    """GEN004 (owner decision 2026-10-26, R13): validate and build read one rule.

    `identity.identity_id_problems` is what `external_links` refuses a page
    on, so a malformed or missing id on a published gene is one ERROR here
    rather than a build that refuses after writing has begun. The same row on
    a gene the gate does not publish reports nothing -- the build never links
    it -- and an unreadable registry reports nothing, because TBL008 already
    names that error. An empty one is reported once.
    """
    identities = {"HGNC:11604": identity} if identity is not None else {}
    unpublished = {"HGNC:2": GeneIdentity(ensembl_gene="ENSG1", ncbi_gene=-3)}

    issues = validate_identity_ids(
        {**identities, **unpublished}, {"HGNC:11604"}, "mirrors/genes.tsv"
    )

    if problem is None:
        assert issues == []
    else:
        assert [(issue.code, issue.severity) for issue in issues] == [("GEN004", Severity.ERROR)]
        assert issues[0].message.startswith("HGNC:11604 is published")
        assert f"these ids: {problem}; the build refuses" in issues[0].message
    assert validate_identity_ids(None, {"HGNC:11604"}, "mirrors/genes.tsv") == []
    # A registry that reads and holds no row is one GEN004 naming it -- never
    # nothing, which is what `or None` at the call site made it until review
    # (2026-10-27), and never one per published gene.
    [empty] = validate_identity_ids({}, {"HGNC:11604", "HGNC:2"}, "mirrors/genes.tsv")
    assert (empty.code, empty.severity) == ("GEN004", Severity.ERROR)
    assert empty.message.startswith("the gene registry holds no gene rows, and 2 published")
    assert validate_identity_ids({}, set(), "mirrors/genes.tsv") == []


def test_previous_symbols_drop_the_names_the_resolver_refuses() -> None:
    """`MADH7` (retired by two genes) and `ODD` (OSR1's retired name and GJA1's
    alias, refused under F27) are not names a search may attach to one gene; a
    clean retired name is, sorted."""
    registry = GeneRegistry(
        [
            {
                "hgnc_id": "HGNC:6772",
                "symbol": "SMAD6",
                "aliases": "",
                "prev_symbols": "MADH6|MADH7",
            },
            {"hgnc_id": "HGNC:6773", "symbol": "SMAD7", "aliases": "", "prev_symbols": "MADH7"},
            {"hgnc_id": "HGNC:4274", "symbol": "GJA1", "aliases": "ODD", "prev_symbols": ""},
            {
                "hgnc_id": "HGNC:8111",
                "symbol": "OSR1",
                "aliases": "",
                "prev_symbols": "ZZ9|ODD|AA1",
            },
        ]
    )
    assert registry.previous_symbols("HGNC:6772") == ("MADH6",)
    assert registry.previous_symbols("HGNC:6773") == ()
    assert registry.previous_symbols("HGNC:8111") == ("AA1", "ZZ9")
    assert registry.previous_symbols("HGNC:4274") == ()
    assert registry.previous_symbols("HGNC:0") == ()


def test_approved_symbol_of_an_id() -> None:
    registry = GeneRegistry(_ROWS)
    assert registry.approved_symbol("HGNC:11604") == "TBX5"
    assert registry.approved_symbol("HGNC:0") is None


def test_a_registry_gene_naming_two_proteins_is_refused() -> None:
    """GEN005: RBM10 as HGNC lists it is an error; one accession or none is silent.

    The fixture mixes a single accession, a null and a blank beside the list, so
    "any non-null cell" and "contains a pipe" are distinguishable.
    """
    frame = pl.DataFrame(
        {
            "hgnc_id": ["HGNC:9896", "HGNC:11604", "HGNC:4173", "HGNC:1152"],
            "uniprot": ["P0DW28|P98175", "Q99593", None, ""],
        }
    )
    issues = validate_one_uniprot_per_gene(frame)
    assert [(i.code, i.severity) for i in issues] == [("GEN005", Severity.ERROR)]
    message = issues[0].message
    assert "HGNC:9896" in message
    assert "P0DW28" in message and "P98175" in message
    assert "pin the gene's own protein in the HGNC recipe" in message
    assert validate_one_uniprot_per_gene(frame.head(0)) == []
