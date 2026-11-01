# src/chd_atlas/vocab.py
from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from typing import Final


class Classification(StrEnum):
    """ClinGen-style gene-disease validity classification."""

    DEFINITIVE = "definitive"
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    DISPUTED = "disputed"
    REFUTED = "refuted"
    NO_KNOWN_ASSOCIATION = "no_known_association"


CLASSIFICATION_RANK: Final[dict[Classification, int]] = {
    Classification.DEFINITIVE: 6,
    Classification.STRONG: 5,
    Classification.MODERATE: 4,
    Classification.LIMITED: 3,
    Classification.DISPUTED: 2,
    Classification.REFUTED: 1,
    Classification.NO_KNOWN_ASSOCIATION: 0,
}


def strongest(classifications: Iterable[Classification]) -> Classification:
    """Return the highest-ranked classification, used to derive headline confidence.

    This does not flag conflicting evidence: DISPUTED and REFUTED still sit on the
    same linear rank as the rest, so a REFUTED classification will not surface if a
    higher-ranked classification is also present for the same gene. Callers that
    need to know whether a gene's evidence is contested must also call
    `has_conflicting_evidence()`.
    """
    items = list(classifications)
    if not items:
        raise ValueError("strongest() requires at least one classification")
    return max(items, key=lambda c: CLASSIFICATION_RANK[c])


CONTESTED: Final[frozenset[Classification]] = frozenset(
    {Classification.DISPUTED, Classification.REFUTED}
)


def has_conflicting_evidence(classifications: Iterable[Classification]) -> bool:
    """True when a gene carries both a supportive and a contesting classification.

    ClinGen treats disputed and refuted as a separate axis rather than weaker
    rungs of the definitive-to-limited ladder, so they cannot be compared on the
    single rank `strongest()` uses. Callers displaying headline confidence must
    call this alongside `strongest()`, so a contested gene is never presented as
    settled.
    """
    items = set(classifications)
    contesting = items & CONTESTED
    supportive = items - CONTESTED - {Classification.NO_KNOWN_ASSOCIATION}
    return bool(contesting and supportive)


def reports_no_association(classifications: Iterable[Classification]) -> bool:
    """True when some authority reported no known association *and* another supports one.

    **The third axis, and issue #13's resolution.** `NO_KNOWN_ASSOCIATION` takes
    neither side of `has_conflicting_evidence` -- it is subtracted from the
    supportive set and is not in `CONTESTED` -- so a gene one submitter calls
    `Definitive` and another calls `No Known Disease Relationship` published
    `has_conflicting_evidence: false` and said nothing else at all. That
    exclusion is right and stays: "a panel looked and found no reported
    evidence" is not "a panel disputes this", and folding the two together
    would give one laboratory's null result the force of a chartered panel's
    refutation. What was wrong was that the disagreement then reached no
    published byte.

    So this is a *separate* flag rather than a widening of the other. ClinGen
    treats the assertion the same way -- a distinct verdict, not a rung of the
    definitive-to-limited ladder -- which is the argument the issue itself makes
    for a third axis over either of the two that exist.

    **Symmetric with `has_conflicting_evidence`: it needs both sides.** A gene
    whose *only* in-scope record is `no_known_association` is not in
    disagreement with anything, and flagging it would tell a reader two
    authorities differ where one spoke. Measured 2026-08-14 over the committed
    mirrors: 9 in-scope genes are in exactly that position and none of them is
    published, so the distinction is currently invisible in the output and is
    pinned by fixture instead.

    GDF1 is the live case and the only one: G2P `Definitive`, Labcorp `Strong`,
    Illumina `No Known Disease Relationship`.
    """
    items = set(classifications)
    supportive = items - CONTESTED - {Classification.NO_KNOWN_ASSOCIATION}
    return Classification.NO_KNOWN_ASSOCIATION in items and bool(supportive)


# The single place an authority's vocabulary becomes the atlas's. Mirrors store
# what was published; this maps it. There is deliberately no `.get(term,
# default)` anywhere downstream -- a term absent from these dicts is an error,
# because a classification silently coerced to an adjacent rung is a confidence
# nobody asserted.
CLINGEN_CLASSIFICATIONS: Final[dict[str, Classification]] = {
    "Definitive": Classification.DEFINITIVE,
    "Strong": Classification.STRONG,
    "Moderate": Classification.MODERATE,
    "Limited": Classification.LIMITED,
    "Disputed": Classification.DISPUTED,
    "Refuted": Classification.REFUTED,
    "No Known Disease Relationship": Classification.NO_KNOWN_ASSOCIATION,
}

# `Supportive` maps to None rather than to a member. It is GenCC's mapping
# exception for submitters that assert an association without grading evidence,
# so it is not a weaker rung of the ladder and must not be placed on one.
# Callers publish the verbatim term and omit the mapped value, which is honest:
# the submitter said "associated", not "how strongly".
GENCC_CLASSIFICATIONS: Final[dict[str, Classification | None]] = {
    "Definitive": Classification.DEFINITIVE,
    "Strong": Classification.STRONG,
    "Moderate": Classification.MODERATE,
    "Limited": Classification.LIMITED,
    "Disputed Evidence": Classification.DISPUTED,
    "Refuted Evidence": Classification.REFUTED,
    "No Known Disease Relationship": Classification.NO_KNOWN_ASSOCIATION,
    "Supportive": None,
}


class ValidityState(StrEnum):
    """How well curated a gene's validity is, published as a fact.

    `UNCURATED` is not an absence. A consumer cannot tell "no authority has
    assessed this gene" from "the build dropped it" by looking at a missing
    field, and that distinction is exactly the silent evidence loss this project
    treats as its characteristic failure.
    """

    EXPERT_CURATED = "expert_curated"
    SUBMITTER_CURATED = "submitter_curated"
    UNCURATED = "uncurated"


class AtlasCuration(StrEnum):
    """Whether the atlas itself has curated evidence for a gene it publishes.

    D21 publishes a gene on an external authority's assessment, which is not the
    atlas's own work. Nearly every published gene carries no
    `LesionAssertion` at all (118 of 119, measured 2026-10-26), so a consumer
    that read `headline_confidence` alone would take a mirrored ClinGen call for
    curated atlas content -- and since the 2026-08-06 widening, some of those
    genes have no ClinGen record either (16 that day, 21 on 2026-10-26), so
    `headline_confidence` is `null` and the same consumer sees a gene with no
    stated confidence from anyone. `admitted_by` in the bundle is what
    says who did admit it.

    Published as a field rather than left inferable from `assertion_count == 0`
    for the reason `validity_state` is: a consumer filtering a browse list needs
    one key to read, not an arithmetic rule to reimplement, and the two must not
    be able to drift apart.
    """

    CURATED = "curated"
    NOT_YET_CURATED = "not_yet_curated"


class ValiditySource(StrEnum):
    CLINGEN = "clingen"
    GENCC = "gencc"


class LesionGroup(StrEnum):
    CONOTRUNCAL = "conotruncal"
    SEPTAL = "septal"
    LVOTO = "lvoto"
    RVOTO = "rvoto"
    AVSD = "avsd"
    HETEROTAXY = "heterotaxy"
    ANOMALOUS_VENOUS_RETURN = "anomalous_venous_return"
    COMPLEX = "complex"
    OTHER = "other"


class Inheritance(StrEnum):
    AD = "AD"
    AR = "AR"
    XL = "XL"
    DE_NOVO = "de_novo"
    SOMATIC = "somatic"
    OLIGOGENIC = "oligogenic"
    UNKNOWN = "unknown"


class Mechanism(StrEnum):
    LOSS_OF_FUNCTION = "loss_of_function"
    GAIN_OF_FUNCTION = "gain_of_function"
    DOMINANT_NEGATIVE = "dominant_negative"
    HAPLOINSUFFICIENCY = "haploinsufficiency"
    UNKNOWN = "unknown"


class SyndromicStatus(StrEnum):
    ISOLATED = "isolated"
    SYNDROMIC = "syndromic"
    BOTH = "both"


class EvidenceClass(StrEnum):
    GENETIC_CASE = "genetic_case"
    GENETIC_SEGREGATION = "genetic_segregation"
    GENETIC_CASE_CONTROL = "genetic_case_control"
    FUNCTIONAL_MODEL = "functional_model"
    EXPRESSION = "expression"
    PROTEOMIC = "proteomic"
    PTM = "ptm"
    REGULATORY = "regulatory"


class EvidenceStrength(StrEnum):
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    CONTRADICTORY = "contradictory"


class SourceTier(StrEnum):
    OWN_CURATION = "own_curation"
    THIRD_PARTY_LITERATURE = "third_party_literature"
    DATABASE_MIRROR = "database_mirror"


class Redistribution(StrEnum):
    """What the source's licence permits the atlas to republish.

    Read by `validate_source_references`, so a source marked `prohibited` is
    reported rather than silently mirrored.
    """

    PERMITTED = "permitted"
    PERMITTED_WITH_ATTRIBUTION = "permitted_with_attribution"
    NON_COMMERCIAL_ONLY = "non_commercial_only"
    PROHIBITED = "prohibited"


class StudyType(StrEnum):
    TRIO_WES = "trio_wes"
    TRIO_WGS = "trio_wgs"
    FAMILY_LINKAGE = "family_linkage"
    CASE_CONTROL = "case_control"
    FUNCTIONAL = "functional"
    PROTEOMICS = "proteomics"
    TRANSCRIPTOMICS = "transcriptomics"
    REVIEW = "review"
    META_ANALYSIS = "meta_analysis"


class BurdenComparator(StrEnum):
    """What a burden row's statistic compared the case count against.

    This is the discriminator the whole burden schema turns on. Every published
    burden analysis asks the same question -- is this gene hit more often than
    expected? -- and differs only in where "expected" comes from:

    - `CONTROL_COHORT`: an observed rate in sequenced controls (case-control).
    - `MUTATION_MODEL`: a per-gene expectation from a trinucleotide mutation
      rate model, the standard for de novo enrichment in trios, where there is
      no control group at all.
    - `NONE`: a case series. A numerator and a denominator, and nothing to
      compare them with -- so no effect size and no p-value may be published,
      which `validate_burden` enforces rather than trusts.

    A row's design is a property of the *row*, not of the study: one paper
    routinely contributes de novo rows and case-control rows, so this cannot be
    hoisted onto the publication record.
    """

    CONTROL_COHORT = "control_cohort"
    MUTATION_MODEL = "mutation_model"
    NONE = "none"


class EffectMeasure(StrEnum):
    """Which quantity a burden row's `effect` column holds.

    Mandatory whenever an effect is present, because an odds ratio of 3.1 and a
    de novo enrichment of 3.1 are different claims and a column headed "effect"
    holding a bare number silently equates them. `EFFECT_MEASURES` below pins
    which measures each comparator can produce; the renderer additionally must
    never print the number without this label.
    """

    ODDS_RATIO = "odds_ratio"
    ENRICHMENT_RATIO = "enrichment_ratio"
    RATE_RATIO = "rate_ratio"


# Which effect measures each comparator can legitimately yield. A mutation
# model cannot produce an odds ratio -- it has no control odds -- and a control
# cohort cannot produce an "enrichment over expectation" because its expectation
# is observed, not modelled. Read by `validate_burden`; the empty set for `NONE`
# is what makes a case series unable to carry an effect at all.
EFFECT_MEASURES: Final[dict[BurdenComparator, frozenset[EffectMeasure]]] = {
    BurdenComparator.CONTROL_COHORT: frozenset(
        {EffectMeasure.ODDS_RATIO, EffectMeasure.RATE_RATIO}
    ),
    BurdenComparator.MUTATION_MODEL: frozenset({EffectMeasure.ENRICHMENT_RATIO}),
    BurdenComparator.NONE: frozenset(),
}


class EffectBound(StrEnum):
    """Why a row reports no finite effect size though the test ran.

    Fisher's exact test returns an unbounded odds ratio when no control carries
    the variant class, and the published upper confidence bound is likewise
    infinite. Measured 2026-08-04 against Supplementary Data 3 of PMID 42230622:
    927 of 138,609 rows across its three strata report `Infinity`, always in
    `fet.odds_ratio` and `fet.ci_95_upper` together and nowhere else.

    Those rows carry the *strongest* signal in the study -- TAB2, 5 syndromic
    carriers of 1,471 against 0 of 45,082 controls -- so dropping them would
    discard the result the reader most wants. (This said "6 syndromic carriers"
    until 2026-08-05; 6 is the *all-cases* count, of 3,876. Verified against the
    mirror, the published bundle and the paper's own Table 1.)

    They cannot be published as a number either:
    `encode_json` sets `allow_nan=False`, and `Infinity` is accepted by Python's
    `json.loads` but rejected by `JSON.parse`, so it would break every page
    while the manifest checksum still verified.

    The row therefore carries a null `effect`, a null `ci_high`, and this flag,
    and a page renders "unbounded" rather than a blank or an invented ceiling.

    There is deliberately no `unbounded_below` member. An odds ratio of zero --
    13,430 rows in the same sheet, where no case carries -- is finite, publishes
    as `0.0`, and comes with a finite upper bound, so it needs no flag.
    """

    UNBOUNDED_ABOVE = "unbounded_above"


class CohortStratum(StrEnum):
    """Which slice of a study's cases a burden row counted.

    Separate from `SyndromicStatus`, which records what a curator asserts about
    a gene. This records how a *published analysis* partitioned its cohort, and
    its `ALL` is not that enum's `BOTH`: `ALL` means the undivided case set,
    while `BOTH` means a gene causes isolated and syndromic disease alike.
    """

    ALL = "all"
    SYNDROMIC = "syndromic"
    NONSYNDROMIC = "nonsyndromic"


class VariantClass(StrEnum):
    """What kind of DNA change a burden row counted.

    Kept separate from `ConsequenceClass` rather than merged into one
    "variant category": deletion-and-loss-of-function and SNV-and-loss-of-function
    are both real cells, and merging them would make "all loss-of-function
    evidence for this gene, whatever the variant type" unaskable.
    """

    SNV_INDEL = "snv_indel"
    CNV_DELETION = "cnv_deletion"
    CNV_DUPLICATION = "cnv_duplication"
    SV_OTHER = "sv_other"


class ConsequenceClass(StrEnum):
    """What the counted variants do to the protein.

    `SYNONYMOUS` is a first-class member, not an artifact to be filtered out on
    the way in. Synonymous variants should show no case-control enrichment, so a
    synonymous row is the study's own negative control -- and several published
    genes carry a nominally significant one. Publishing it beside the
    loss-of-function row for the same gene, study and stratum is what lets a
    reader judge how well that comparison was calibrated, which no other CHD
    resource shows. `build/bundles.py` groups on that key for exactly this
    reason; nothing anywhere may drop a row for being synonymous.
    """

    LOF = "lof"
    MISSENSE_DAMAGING = "missense_damaging"
    MISSENSE_ALL = "missense_all"
    SYNONYMOUS = "synonymous"
    ALL_CODING = "all_coding"
    # Loss-of-function and damaging missense counted together, which is the
    # composite most CHD burden papers report as their *primary* analysis --
    # PMID:40127276's 60 genes are defined by it. Distinct from `ALL_CODING`,
    # which would additionally include tolerated missense and synonymous
    # variants and so is a different denominator, and not derivable from the two
    # component rows: the composite carries its own p-value, which is not a
    # function of theirs.
    #
    # A page showing all three must say that this row is the union of the other
    # two rather than a third independent result. `build/burden.py` orders it
    # first for that reason -- it is the headline, and its components read as a
    # breakdown of it rather than as separate findings.
    DAMAGING = "damaging"


class VariantOrigin(StrEnum):
    """Whether the counted variants were required to be de novo.

    Without it, a trio study's de novo loss-of-function count and a
    case-control study's rare-inherited loss-of-function count differ in no
    published column and read as the same measurement.
    """

    DE_NOVO = "de_novo"
    INHERITED = "inherited"
    ANY = "any"
    # PMID:40127276's "transmitted/unphased variants" (TUVs): its case-control
    # arm removes de novo mutations, then pools variants whose transmission is
    # *known* (in its 3,887 trios) with variants whose transmission is *unknown*
    # (in its 7,668 singletons). Neither `INHERITED` nor `ANY` is true of that
    # set -- `ANY` because de novo variants were deliberately excluded, and
    # `INHERITED` because most of the set was never phased at all. Measured in
    # its Dataset S4: of 36,054 variant rows, 27,429 are `Unphased`, 8,221
    # `Trans` and 404 `DNM`, so the unphased majority is the whole difficulty.
    TRANSMITTED_OR_UNPHASED = "transmitted_or_unphased"


class CountUnit(StrEnum):
    """What a burden row's numerator and denominator actually count.

    **A count without its unit is the same defect as an effect without its
    measure** -- see `EffectMeasure`. "12 / 21,768" and "6 / 6,746" are both
    plausible burden cells, and read side by side they claim to be the same kind
    of measurement. They are not, and no other column in this schema
    distinguishes them.

    Each member fixes *both* halves of the pair, which is why one column
    suffices where the numerator and denominator are counted in different things:

    - `INDIVIDUALS`: people carrying at least one qualifying variant, out of
      people sequenced. Someone with two qualifying variants counts once.
    - `ALLELES`: qualifying alleles observed, out of alleles called. Someone with
      two qualifying variants counts twice, and the denominator is roughly twice
      the sample size and varies gene by gene with coverage.
    - `DE_NOVO_MUTATIONS`: de novo mutations observed, out of *trios*. The
      denominator is families, not alleles, because a trio is what it takes to
      call one.

    The distinction was found, not assumed. Measured 2026-08-05 against
    PMID:40127276's Dataset S6 and Dataset S4: its D-Mis case-control `Obs`
    equals the number of qualifying variant rows in 245 of 248 genes but the
    number of distinct probands in only 235 -- CACNA1A 80 variants across 79
    probands, TSC1 36 across 34, LRP1 135 across 132 -- and `Obs` follows the
    variant every time. Publishing that column under this atlas's previous
    header, "cases (carriers / n)", would have asserted 21,768 people sequenced
    where 11,555 were, and called alleles carriers.

    Nothing here is inferable from `comparator`. A case-control study may count
    people (PMID:42230622) or alleles (PMID:40127276), and both are ordinary
    practice, so `validate_burden` asserts only the one implication that cannot
    fail: counting de novo mutations requires a de novo variant set.
    """

    INDIVIDUALS = "individuals"
    ALLELES = "alleles"
    DE_NOVO_MUTATIONS = "de_novo_mutations"


class StatisticalTest(StrEnum):
    """Which test produced a burden row's p-value.

    Deliberately not constrained against `BurdenComparator` the way
    `EFFECT_MEASURES` constrains the measure. A binomial test appears against
    both a control rate and a modelled rate in the published literature, so a
    compatibility table here would encode a rule this project has not measured.
    The one rule that *is* asserted -- a `NONE` comparator carries no test at
    all -- lives in `validate_burden`.
    """

    FISHER_EXACT = "fisher_exact"
    POISSON = "poisson"
    BINOMIAL = "binomial"
    # PLINK's `--mperm` CNV association test, which draws its null by permuting
    # case/control labels rather than from a closed-form distribution. Added
    # 2026-08-05 with PMID:34324492, whose whole CNV analysis uses it.
    PERMUTATION = "permutation"
    # DeNovoWEST and MuPIT both compare an observed de novo count against a
    # per-gene mutation-rate model. Named for the family rather than the tool:
    # PMID:34324492 reports `p_dnv` as the minimum of the two, so no single tool
    # name would be true of the column.
    MUTATION_RATE = "mutation_rate"


class PvalueAdjustment(StrEnum):
    """How a published corrected p-value was corrected.

    **The atlas computes no correction and never will** -- that would be
    authoring a statistic (D12/D33). This records one a study published.

    It exists because the 2026-08-05 review found the sharpest gap in this
    layer: gene pages showed 187 uncorrected p-values, of which 32 clear 0.05
    and 3 survive Bonferroni over the study's own 138,609 tests. That study
    published no corrected column, so the page names the denominator instead.
    PMID:34324492 *does* publish one -- PLINK's EMP2 beside EMP1, and a
    Bonferroni-adjusted p beside the raw de novo p -- and dropping it would
    discard the one number that answers the question the review raised.

    `FAMILYWISE_PERMUTATION` is PLINK's EMP2: the max(T) family-wise error rate
    over the same permutations that produced EMP1, not a Bonferroni factor.
    """

    BONFERRONI = "bonferroni"
    FAMILYWISE_PERMUTATION = "familywise_permutation"
    BENJAMINI_HOCHBERG = "benjamini_hochberg"


class CombinationMethod(StrEnum):
    """How a study combined two of its own burden arms into one statistic.

    Not a correction, and not a pooling across studies. PMID:40127276 merges
    each gene's de novo Poisson p and case-control Fisher p by the Joint-Local
    False Discovery Rate method (Jiang & Yu 2017, PMID:28011772), within one
    cohort, and publishes a meta-analysis p and a local FDR. The atlas computes
    no combination of its own (D12/D33); this records one a study published.
    """

    JOINT_LOCAL_FDR = "joint_local_fdr"


class FeaturedTopic(StrEnum):
    GENOMICS = "genomics"
    VARIANTS = "variants"
    EXPRESSION = "expression"
    PROTEOMICS = "proteomics"
    REVIEW = "review"


class Perturbation(StrEnum):
    KNOCKOUT = "knockout"
    KNOCKIN = "knockin"
    CONDITIONAL_KO = "conditional_ko"
    MORPHOLINO = "morpholino"
    CRISPR = "crispr"
    OVEREXPRESSION = "overexpression"
    PATIENT_VARIANT_KNOCKIN = "patient_variant_knockin"


class Zygosity(StrEnum):
    HOMOZYGOUS = "homozygous"
    HETEROZYGOUS = "heterozygous"
    COMPOUND_HET = "compound_het"
    HEMIZYGOUS = "hemizygous"
    NOT_APPLICABLE = "n/a"


class PhenocopyAssessment(StrEnum):
    YES = "yes"
    PARTIAL = "partial"
    NO = "no"
    NOT_ASSESSED = "not_assessed"


class RescueOutcome(StrEnum):
    NOT_ATTEMPTED = "not_attempted"
    NO_RESCUE = "no_rescue"
    PARTIAL_RESCUE = "partial_rescue"
    FULL_RESCUE = "full_rescue"


class Archive(StrEnum):
    PRIDE = "pride"
    GEO = "geo"
    ARRAYEXPRESS = "arrayexpress"
    EGA = "ega"
    UCSC_CELL_BROWSER = "ucsc_cell_browser"


class Technology(StrEnum):
    BULK_RNASEQ = "bulk_rnaseq"
    SCRNASEQ = "scrnaseq"
    MICROARRAY = "microarray"
    TMT_MS = "tmt_ms"
    LFQ_MS = "lfq_ms"
    PHOSPHOPROTEOMICS = "phosphoproteomics"


class DatasetDesign(StrEnum):
    """Whether a dataset's rows are contrasts or abundances.

    Declared rather than inferred from `len(contrasts)`, so a curator who
    forgot the contrasts is distinguishable from an abundance atlas that
    correctly has none.
    """

    CONTRAST = "contrast"
    PROFILE = "profile"
    # One abundance per cell type and region (Farah et al. 2024), not per
    # tissue and stage: it declares cell_types and regions where PROFILE
    # declares cardiac_tissues and stages.
    CELL_TYPE_PROFILE = "cell_type_profile"


MODEL_ORGANISMS: Final[dict[str, str]] = {
    "NCBITaxon:9606": "Homo sapiens",
    "NCBITaxon:10090": "Mus musculus",
    "NCBITaxon:10116": "Rattus norvegicus",
    "NCBITaxon:7955": "Danio rerio",
    "NCBITaxon:8355": "Xenopus laevis",
    "NCBITaxon:8364": "Xenopus tropicalis",
    "NCBITaxon:9031": "Gallus gallus",
    "NCBITaxon:9823": "Sus scrofa",
}
