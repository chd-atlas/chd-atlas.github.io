# src/chd_atlas/build/burden.py
"""Per-gene rare-variant burden, loaded from `mirrors/burden.tsv`.

Loading and grouping only. The bundle shape is `burden_payload`'s and the page
markup is `pages._burden_section`'s, the same split `omics.py` has with
`bundles.py`: this module knows the mirror, and nothing here knows HTML.

**Nothing filters a row out.** Not the synonymous rows, which are the study's own
negative control and the reason this layer is worth publishing at all; not the
rows with no case carrier, which are evidence *against* enrichment; not the rows
whose odds ratio came back unbounded, which carry the strongest signals in the
data. A gene's rows reach the page exactly as the study computed them, and the
page's job is to label them, not to choose between them.

The one thing that is *not* published is anything this atlas computed. There is
no pooled p-value across studies and no per-gene summary verdict, because the
CHD literature reuses cohorts -- DDD contributed to more than one CHD exome
study, though only one cited today names it (`curation/cohorts.yaml`) -- so a pooled
statistic would count the same children twice. `shared_cohorts` is the
alternative: it names the overlap instead of pretending it away, and it
derives nothing beyond a set intersection over what each study already declared.
"""

from __future__ import annotations

from collections.abc import Callable, Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Final, NamedTuple

from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import COHORTS
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.review import check_review
from chd_atlas.tables import BURDEN, BURDEN_COMBINED, read_table

# The order strata are presented in, everywhere. Not alphabetical, which would
# put `nonsyndromic` before `syndromic` and both before the whole-cohort figure a
# reader looks at first; not the mirror's sort order either, which is
# alphabetical for the same reason. This is the reading order.
STRATUM_ORDER: Final[tuple[str, ...]] = ("all", "syndromic", "nonsyndromic")

# Within a stratum: the result, then the other result, then the control that
# calibrates them. `synonymous` last is a deliberate reading order rather than an
# ordering of importance -- it is the row that tells you whether to believe the
# two above it.
#
# `damaging` is first because it is the *composite* of `lof` and
# `missense_damaging`, and the primary analysis of the study that reports it
# (PMID:40127276 defines its 60 genes by it). Putting the union above its two
# components makes them read as a breakdown of the headline; putting it after
# them, or alphabetically -- where it would land first anyway, by coincidence --
# would make three rows read as three independent findings.
CONSEQUENCE_ORDER: Final[tuple[str, ...]] = (
    "damaging",
    "lof",
    "missense_damaging",
    "missense_all",
    "all_coding",
    "synonymous",
)


@dataclass(frozen=True)
class BurdenRow:
    """One published burden statistic, as `mirrors/burden.tsv` records it.

    A flat mirror of the schema rather than a reshaped view of it. The columns
    that constrain each other -- `comparator` against the control counts,
    `effect` against `effect_measure` and `effect_bound` -- are already checked
    by `validate/burden.py`, and reshaping them here into, say, a union type
    would move that check into the build, where a violation raises mid-write
    instead of being reported by the gate.

    `n_chd_case_carriers` sits last, out of the schema's order, only because a
    field with a default must follow those without; `load_burden` always
    passes it.
    """

    study: str
    gene: str
    cohort_stratum: str
    lesion_group: str | None
    variant_class: str
    consequence_class: str
    origin: str
    maf_max: float | None
    count_unit: str
    n_case_carriers: int
    n_cases: int
    comparator: str
    n_control_carriers: int | None
    n_controls: int | None
    expected_count: float | None
    effect: float | None
    effect_measure: str | None
    effect_bound: str | None
    ci_low: float | None
    ci_high: float | None
    pvalue: float | None
    pvalue_test: str | None
    pvalue_adjusted: float | None
    pvalue_adjustment: str | None
    case_cohorts: tuple[str, ...]
    control_cohorts: tuple[str, ...]
    method_note: str | None
    source: str
    n_chd_case_carriers: int | None = None


def _cohorts(value: object) -> tuple[str, ...]:
    """The `;`-joined cohort column as a tuple, in the order it was written.

    Not sorted: the package that wrote the mirror lists them in the order the study's own
    cohort description names them (the recipe's documented rule, not re-measured),
    and a reader meeting "cnchd, ddd, nottingham"
    should meet the study's ordering rather than the alphabet's. The mirror cell
    is a single string, so this order is stable across builds by construction.
    """
    if value is None:
        return ()
    return tuple(part for part in str(value).split(";") if part)


def _stratum_rank(stratum: str) -> int:
    """Index in `STRATUM_ORDER`; an unknown stratum sorts last rather than raising."""
    return STRATUM_ORDER.index(stratum) if stratum in STRATUM_ORDER else len(STRATUM_ORDER)


def _consequence_rank(consequence: str) -> int:
    """Index in `CONSEQUENCE_ORDER`; an unknown class sorts last rather than raising."""
    return (
        CONSEQUENCE_ORDER.index(consequence)
        if consequence in CONSEQUENCE_ORDER
        else len(CONSEQUENCE_ORDER)
    )


def _sort_key(row: BurdenRow) -> tuple[object, ...]:
    """Reading order within one gene: study, then stratum, then consequence.

    `STRATUM_ORDER` and `CONSEQUENCE_ORDER` are looked up by index rather than
    compared as strings, so the sequence a reader sees is the one declared above
    and not the alphabet's. A value outside either tuple sorts last rather than
    raising: it is a vocabulary member `validate_table` accepted and this module
    has not been taught to order, and dropping it or dying on it would be worse
    than putting it at the end.
    """
    return (
        row.study,
        _stratum_rank(row.cohort_stratum),
        row.cohort_stratum,
        _consequence_rank(row.consequence_class),
        row.consequence_class,
        row.origin,
    )


def load_burden(root: Path) -> dict[str, list[BurdenRow]]:
    """Every burden row, grouped by HGNC id and ordered for reading.

    Returns an empty mapping when the mirror is absent -- there is none before
    any is curated -- and when it is unreadable or missing a column, which
    `validate_table` reports as TBL000/TBL001 and `build_site` refuses on, so
    neither is reachable behind the gate.

    The tie-break in `_sort_key` runs to `origin` and stops there, which is not
    the mirror's full key: `lesion_group`, `variant_class` and `maf_max` are
    excluded from it.

    **Both futures this paragraph used to anticipate have arrived.** It said
    "today every row shares one value of each" and described what would happen
    "when a second frequency threshold or a CNV row arrives"; PMID:34324492
    brought `cnv_deletion` and PMID:40127276 brought a second `maf_max` plus a
    null one, in the same release. The conclusion still holds -- reading order is
    total -- but for a different reason than the one stated: `study` and `origin`
    separate the rows that now share a gene and stratum, not the accident of a
    single-valued corpus.

    Where two rows do compare equal, `sorted` is stable, so they hold the
    mirror's order -- which `validate_sort_order` pins. Order stays deterministic
    either way; it is only the *reading* order that would need extending.
    """
    path = root / "mirrors" / "burden.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, BURDEN)
    if frame is None or any(column not in frame.columns for column in BURDEN.column_names):
        return {}

    grouped: dict[str, list[BurdenRow]] = {}
    for record in frame.iter_rows(named=True):
        row = BurdenRow(
            study=str(record["study"]),
            gene=str(record["gene"]),
            cohort_stratum=str(record["cohort_stratum"]),
            lesion_group=record["lesion_group"],
            variant_class=str(record["variant_class"]),
            consequence_class=str(record["consequence_class"]),
            origin=str(record["origin"]),
            maf_max=record["maf_max"],
            count_unit=str(record["count_unit"]),
            n_case_carriers=int(record["n_case_carriers"]),
            n_cases=int(record["n_cases"]),
            comparator=str(record["comparator"]),
            n_control_carriers=record["n_control_carriers"],
            n_controls=record["n_controls"],
            expected_count=record["expected_count"],
            effect=record["effect"],
            effect_measure=record["effect_measure"],
            effect_bound=record["effect_bound"],
            ci_low=record["ci_low"],
            ci_high=record["ci_high"],
            pvalue=record["pvalue"],
            pvalue_test=record["pvalue_test"],
            pvalue_adjusted=record["pvalue_adjusted"],
            pvalue_adjustment=record["pvalue_adjustment"],
            case_cohorts=_cohorts(record["case_cohorts"]),
            control_cohorts=_cohorts(record["control_cohorts"]),
            method_note=record["method_note"],
            source=str(record["source"]),
            n_chd_case_carriers=record["n_chd_case_carriers"],
        )
        grouped.setdefault(row.gene, []).append(row)

    return {gene: sorted(rows, key=_sort_key) for gene, rows in grouped.items()}


class Overlap(NamedTuple):
    """Why two studies are not independent: what they share, or are declared to.

    `shared` is every collection both studies cite. `declared` is every pair
    `(left's collection, right's collection)` that `curation/cohorts.yaml`
    records as overlapping (`Cohort.overlaps_with`), the left study's first --
    two ids for partly the same people, which no shared id can show.
    """

    shared: tuple[str, ...]
    declared: tuple[tuple[str, str], ...]


def shared_cohorts(
    rows: Iterable[BurdenRow], cohorts: Mapping[str, Cohort]
) -> dict[tuple[str, str], Overlap]:
    """Which sample collections each pair of studies has in common, or overlapping.

    The reason `curation/cohorts.yaml` exists. The atlas publishes per-study rows
    and no pooled statistic, because pooling studies that share samples inflates
    the result -- but a reader looking at two rows will compare them anyway, and
    nothing on the page would tell them the two are not independent. This names
    the overlap so the page can.

    Keyed on the ordered pair with the smaller study id first, so one pair yields
    one entry however the rows are encountered, and the mapping iterates in a
    stable order. Pairs sharing nothing are omitted: the interesting fact is the
    overlap, and an entry per independent pair would bury it.

    **No longer latent, and live by declaration rather than by a shared id**
    (owner decision 2026-10-27, I3). The three studies cite disjoint ids, so
    this returned an empty mapping, and `concordance.cohort_families`, which
    walks these edges, made them three families. But PMID:34324492 cites
    `decipher` and PMID:42230622 `ddd`, and DDD links its participants into
    DECIPHER: an author of the first answered that its case set does or likely
    does include DDD participants. `cohorts` is the registry, and a pair of
    collections it declares overlapping (`Cohort.overlaps_with`) is an edge as
    a shared id is, so the two studies are one family. Required, never
    defaulted, so no consumer can read families without the declarations
    (CLAUDE.md section 4.28): the families, every card's dataset count, the
    matrix, the browse strip and the page's notice all come through here.
    """
    by_study: dict[str, set[str]] = {}
    for row in rows:
        by_study.setdefault(row.study, set()).update(row.case_cohorts, row.control_cohorts)
    declared_with = {
        identifier: set(cohort.overlaps_with) for identifier, cohort in cohorts.items()
    }

    overlaps: dict[tuple[str, str], Overlap] = {}
    studies = sorted(by_study)
    for index, left in enumerate(studies):
        for right in studies[index + 1 :]:
            common = by_study[left] & by_study[right]
            declared = sorted(
                (mine, theirs)
                for mine in by_study[left]
                for theirs in by_study[right]
                if theirs in declared_with.get(mine, ())
            )
            if common or declared:
                overlaps[(left, right)] = Overlap(tuple(sorted(common)), tuple(declared))
    return overlaps


def burden_payload(rows: Iterable[BurdenRow]) -> list[Json]:
    """One gene's burden rows, JSON-ready, in the order they will be rendered.

    Every key is present on every row, `null` where the comparator does not
    populate it, for the reason `bundles._validity_record` publishes ClinGen-only
    and GenCC-only fields as `null` rather than omitting them: an array whose
    objects differ in shape is a trap for a consumer that reads a field off
    `rows[0]` and expects the same key on `rows[1]`.

    Cohorts publish as arrays rather than as the mirror's `;`-joined string. The
    separator is a storage detail of a TSV that keeps one value per cell; a
    consumer should not have to know it, and a client splitting on `;` itself is
    a client that will one day meet a cohort id containing one.
    """
    return [
        {
            "study": row.study,
            "cohort_stratum": row.cohort_stratum,
            "lesion_group": row.lesion_group,
            "variant_class": row.variant_class,
            "consequence_class": row.consequence_class,
            "origin": row.origin,
            "maf_max": row.maf_max,
            # Qualifies all four count keys below. A consumer reading
            # `n_case_carriers` off two studies and comparing them is doing
            # arithmetic on alleles and people unless this key is there to stop
            # it -- the same defect the gene page carried until 2026-08-05.
            "count_unit": row.count_unit,
            "n_case_carriers": row.n_case_carriers,
            "n_cases": row.n_cases,
            # How many of the case carriers are CHD cases, where the study
            # counted that apart from a case set holding another condition;
            # `null` elsewhere (owner decision R2, 2026-10-26).
            "n_chd_case_carriers": row.n_chd_case_carriers,
            "comparator": row.comparator,
            "n_control_carriers": row.n_control_carriers,
            "n_controls": row.n_controls,
            "expected_count": row.expected_count,
            "effect": row.effect,
            "effect_measure": row.effect_measure,
            "effect_bound": row.effect_bound,
            "ci_low": row.ci_low,
            "ci_high": row.ci_high,
            "pvalue": row.pvalue,
            "pvalue_test": row.pvalue_test,
            "pvalue_adjusted": row.pvalue_adjusted,
            "pvalue_adjustment": row.pvalue_adjustment,
            "case_cohorts": list(row.case_cohorts),
            "control_cohorts": list(row.control_cohorts),
            "method_note": row.method_note,
            "source": row.source,
        }
        for row in rows
    ]


@dataclass(frozen=True)
class CombinedRow:
    """One study's own combination of two of its burden arms (`burden_combined.tsv`).

    `local_fdr` is a local false discovery rate and never a corrected p; see
    `BURDEN_COMBINED` and `pages._combined_line`, which is why it is not called
    `pvalue_adjusted` anywhere.
    """

    study: str
    gene: str
    cohort_stratum: str
    lesion_group: str | None
    consequence_class: str
    combines: tuple[str, ...]
    method: str
    pvalue: float
    local_fdr: float
    method_note: str | None
    source: str


def _combined_sort_key(row: CombinedRow) -> tuple[object, ...]:
    """The arm rows' reading order: study, stratum, consequence."""
    return (
        row.study,
        _stratum_rank(row.cohort_stratum),
        row.cohort_stratum,
        row.lesion_group or "",
        _consequence_rank(row.consequence_class),
        row.consequence_class,
    )


def load_combined(root: Path) -> dict[str, list[CombinedRow]]:
    """Every combined row, grouped by HGNC id and ordered for reading.

    Empty when the mirror is absent, unreadable or missing a column. The last two
    are TBL000/TBL001 errors the gate refuses on, so behind the gate only
    "absent" is reachable.
    """
    path = root / "mirrors" / "burden_combined.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, BURDEN_COMBINED)
    if frame is None or any(c not in frame.columns for c in BURDEN_COMBINED.column_names):
        return {}

    grouped: dict[str, list[CombinedRow]] = {}
    for record in frame.iter_rows(named=True):
        row = CombinedRow(
            study=str(record["study"]),
            gene=str(record["gene"]),
            cohort_stratum=str(record["cohort_stratum"]),
            lesion_group=record["lesion_group"],
            consequence_class=str(record["consequence_class"]),
            # The `;` convention the cohort columns use, split by the same
            # function and kept in the order the study's arms were written in.
            combines=_cohorts(record["combines"]),
            method=str(record["method"]),
            pvalue=float(record["pvalue"]),
            local_fdr=float(record["local_fdr"]),
            method_note=record["method_note"],
            source=str(record["source"]),
        )
        grouped.setdefault(row.gene, []).append(row)
    return {gene: sorted(rows, key=_combined_sort_key) for gene, rows in grouped.items()}


def combined_payload(rows: Iterable[CombinedRow]) -> list[Json]:
    """One gene's combined rows, JSON-ready, every key on every object.

    Every key on every object for the reason `burden_payload` gives; `combines`
    publishes as an array for the reason its cohort columns do.
    """
    return [
        {
            "study": row.study,
            "cohort_stratum": row.cohort_stratum,
            "lesion_group": row.lesion_group,
            "consequence_class": row.consequence_class,
            "combines": list(row.combines),
            "method": row.method,
            "pvalue": row.pvalue,
            "local_fdr": row.local_fdr,
            "method_note": row.method_note,
            "source": row.source,
        }
        for row in rows
    ]


def combined_census(
    combined: Mapping[str, Sequence[CombinedRow]], published: Collection[str]
) -> int:
    """Combined results reaching a published bundle, never the mirror's row count."""
    return sum(len(combined.get(gene, ())) for gene in published)


@dataclass(frozen=True)
class BurdenCensus:
    """The burden figures, derived once and published in two places.

    `index.html`, in its hero band and its cards, and `manifest.json`'s
    `counts`. They come from one object for the same reason `published` is
    threaded down from `build_site` rather than recomputed: two derivations are
    two things that drift, and a front page whose figure contradicts the payload
    a program reads instead of the page is worse than one carrying neither.

    That is also why this lives in `burden.py` rather than in `landing.py`, where
    it was first written. A census with two consumers in different layers belongs
    to neither of them.
    """

    rows: int
    genes: int
    families: int
    # The studies whose rows reach a published bundle, sorted: the population
    # the front page's own-lab sentence counts (owner decision 2026-10-26,
    # R7). Not published in `manifest.json`, which names its `counts` keys one
    # by one.
    studies: tuple[str, ...]
    # Published genes with at least one row that survives its own study's
    # correction: the genes whose card lists a rare-variant result (owner
    # decision 2026-10-26, R9). Published as `counts.burden_surviving_genes`.
    surviving_genes: int


def burden_census(
    burden: Mapping[str, Sequence[BurdenRow]],
    published: Collection[str],
    families: tuple[frozenset[str], ...],
    *,
    survives: Callable[[BurdenRow], bool],
) -> BurdenCensus:
    """What the burden evidence *this site publishes* amounts to.

    **`published` is the restriction, and dropping it would overstate the atlas
    by a factor of five.** `mirrors/burden.tsv` is deliberately wider than the
    publication gate — measured 2026-08-05, 1,475 rows over 150 genes, of which
    127 genes publish no page (CLAUDE.md section 10 records that asymmetry as
    intended, so re-mirroring is not needed the day D21 widens). A front page
    counting the mirror would advertise 1,475 statistics a reader cannot reach
    from any link on it, and `manifest.json` would promise a consumer five times
    the rows it can fetch. It counts the rows that reach a published bundle --
    290 that day.

    **`families` is the dataset count, never `len({row.study for row in rows})`.**
    A cohort family is a connected component over shared sample collections, so
    two papers drawing on one collection describe the same people and count once
    (`concordance.cohort_families`). The two were equal over the committed
    corpus — three studies citing disjoint collections — which is precisely why
    they were separated here and separated again in the tests. They parted on
    owner decision 2026-10-27 (I3): DECIPHER and DDD are declared overlapping
    (`Cohort.overlaps_with`), so PMID:34324492 and PMID:42230622 are one family
    and the census reads 2 families for 3 studies. Taken as a parameter rather than
    recomputed, so this page and the strip on every browse row count the same
    families.

    Genes are counted as those with at least one *published* row rather than as
    `len(published)`. They were equal at 23 of 23 when this was written, and
    stating the published-gene count under a burden label would have kept
    reading as 23 on the day a gene was published with no burden evidence at
    all -- which has since happened: they differ by every published gene no
    rare-variant study covered.

    **`families` is the one figure here NOT restricted to `published`, and that
    asymmetry is deliberate.** `rows` and `genes` count what a reader can reach;
    `families` counts the cohort families in the corpus, including any whose
    studies tested no published gene. Restricting it would be worse, not better:
    a family exists because two studies share a sample collection, and dropping
    a member because none of its genes cleared the publication gate can split
    one family into what then presents as two independent datasets -- inventing
    independence D33 exists to deny. So the front page can state a family count
    one higher than the number of families any published gene was tested by.
    Reviewed 2026-08-06 and kept; the alternative was measured and inverts a
    real overlap.

    **`studies` is restricted to `published`, like `rows`** (owner decision
    2026-10-26, R7): the front page says how many of the rare-variant studies
    whose results the site publishes are by an author of this atlas, and a
    study whose rows all sit outside the gate has no block on any page to say
    so. Sorted, because it is derived from a set.

    **`surviving_genes` counts the published genes whose card lists a
    result** (owner decision 2026-10-26, R9): a gene with at least one row
    `survives` passes. `survives` is `concordance.survives_own_correction`,
    the predicate the card's list reads (`pages._card_rare_variants` groups
    exactly the rows it passes, with no exclusion of its own), handed in by
    the caller rather than imported, because `concordance` imports this
    module; required, so no caller can count by a rule of its own. Measured on
    a real build the day it was written: 32 of the 118 genes with a published
    row, and every one of the 32 cards lists a line.
    """
    rows = 0
    genes = 0
    studies: set[str] = set()
    for gene in published:
        count = len(burden.get(gene, ()))
        rows += count
        if count:
            genes += 1
        studies.update(row.study for row in burden.get(gene, ()))
    surviving_genes = sum(
        1 for gene in published if any(survives(row) for row in burden.get(gene, ()))
    )
    return BurdenCensus(
        rows=rows,
        genes=genes,
        families=len(families),
        studies=tuple(sorted(studies)),
        surviving_genes=surviving_genes,
    )


def cohort_registry(cohorts: Iterable[Cohort]) -> dict[str, Cohort]:
    """Cohort id -> the curated record, from `corpus.cohorts`.

    Returns the whole record rather than just the name. The first version
    returned `{id: name}` and dropped `description` -- which is where
    `models/cohort.py` says the caveats qualifying every number drawn from a
    collection belong, and where a curator had written that UK Biobank's
    controls are adults recruited at 40-69 while the CHD cases were largely
    enrolled in childhood. Measured 2026-08-05: no cohort description reached a
    single published byte. The curator wrote the caveat, the model documented
    where it goes, and no reader could reach it.

    A cohort absent from the registry keeps its id as its label at the render
    site -- the same fallback `build_genes` applies to a gene missing from
    `mirrors/genes.tsv`. `validate_burden_references` reports the absence as
    BUR009 and the gate refuses on it, so that fallback is unreachable behind
    `build_site`.
    """
    return {str(cohort.id): cohort for cohort in cohorts}


def build_cohorts(cohorts: Iterable[Cohort], emitter: Emitter) -> None:
    """Publish `cohorts.json`: the table a bare cohort id resolves against.

    **Every burden row names its collections by bare id and nothing mapped those
    ids to anything.** Measured 2026-08-06 over the published corpus: 13 distinct
    ids appear across 915 rows -- `taa_cases`, `gnomad_controls`, `ukbb` and ten
    more -- and `curation/cohorts.yaml` reached zero published bytes, so every
    programmatic consumer got the numbers without the caveats that qualify them.
    `Cohort.url` reached no byte at all. The descriptions rendered only inside a
    `<details>` on the gene page, which is a fact about the HTML and not about
    the API.

    The caveats are the point, not decoration. `taa_cases` is 777 sporadic
    thoracic aortic aneurysm probands, a condition **outside this atlas's
    scope**, in a case set the study reports as CHD (whether any of them also
    has a congenital heart lesion is recorded nowhere this atlas has read);
    `ukbb`'s controls are adults recruited at 40-69 against cases a substantial
    proportion of whom were enrolled in childhood; `ddd` ascertains on developmental disorder, so
    its contribution is enriched for syndromic CHD. A consumer computing anything
    from a row without those sentences is computing something else.

    **The whole registry is published, not only the ids the published rows
    cite.** The two sets were identical when measured -- 20 and 20 on
    2026-10-27 -- and identical figures are one figure to every test (CLAUDE.md
    section 4.30), so this says which one it is: this is a *resolution table*,
    and a table that omitted an id
    some row cites would 404 exactly where a reader needed it. The direction is
    also the safe one -- `validate_burden_references` (BUR009) already refuses a
    row citing a cohort this file does not carry, so publishing the registry
    entire means every id in every row resolves by construction, while
    publishing the cited subset would make that guarantee depend on the gate.

    Sorted by id: `encode_json`'s `sort_keys` orders dict keys and never list
    elements, so an unsorted array here checksums differently between two builds
    of one commit. `url` is emitted as `null` rather than omitted where a
    collection has no public page (11 of the 20 on 2026-10-27), because an
    object whose shape varies is a trap for a consumer reading a field off one
    record and expecting it on the next -- the rule `admission_provenance` and
    `burden_payload` already follow. `review_status` and `review_note` follow
    it too: both `null` on a record no open question is recorded for.

    **Both or neither, checked again here** (review, 2026-10-27, N2): the
    model refuses one without the other where the file loads, but a record
    built without validation (`model_construct`, `model_copy`) reached this
    emitter and published `"pending_review"` beside a `null` note -- a
    question with no statement of what it is. `models.review.check_review`
    raises before anything is written, a guard on a bypassed gate like the
    pages' `render.review_status_note`.
    """
    ordered = sorted(cohorts, key=lambda item: str(item.id))
    for cohort in ordered:
        check_review(str(cohort.id), cohort.review_status, cohort.review_note)
    emitter.write_json(
        COHORTS,
        {
            "cohorts": [
                {
                    "id": str(cohort.id),
                    "name": cohort.name,
                    "description": cohort.description,
                    "url": cohort.url,
                    # The registry's declared overlaps (owner decision
                    # 2026-10-27, I3), sorted, `[]` where none: the edges
                    # `cohort_families` merges on besides a shared id.
                    "overlaps_with": sorted(str(other) for other in cohort.overlaps_with),
                    # A question about the record this atlas leaves open
                    # (owner's standing rule, 2026-10-27; `models/review.py`),
                    # `null` and `null` where none is recorded -- never
                    # omitted, like `url`.
                    "review_status": (
                        None if cohort.review_status is None else cohort.review_status.value
                    ),
                    "review_note": cohort.review_note,
                }
                for cohort in ordered
            ]
        },
    )
