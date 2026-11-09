# src/chd_atlas/build/pages.py
"""One HTML page per published gene, plus the browse page over all of them.

Everything a gene page shows has been published as JSON since the validity
backbone landed. `genes/HGNC_11604.json` carries seven mirrored validity
records -- four of them with a report URL, and all four render as links --
alongside a curated assertion whose evidence quote, class, strength and
publication this module renders, under a summary of its lesion groups,
inheritance, mechanism and syndromic status (measured 2026-08-04 against a real
build of the committed corpus). All of it sat at a URL, with a verifying
checksum, and with nothing rendering it. This module is what makes that content
reachable by a person rather than only by a program.

It does not make *all* of it reachable, and this paragraph claimed it did until
the claim was measured. On the same real build, `HP:0001631` appears once in
`genes/HGNC_11604.json` and **zero** times in `genes/HGNC_11604.html`; the
evidence locator's `Abstract` likewise, 1 and 0. The assertion's phenotype terms
and the in-text locator on every evidence item reach no page. Rendering them is
a product decision rather than a defect to be patched, so what is recorded here
is only that the JSON, not the page, is the complete artifact.

Written through `Emitter.write_text`, never `write_json`. `Json`'s union
includes `str`, so `write_json(path, page)` type-checks and then runs the whole
document through `json.dumps`, publishing it quoted and escaped as a JSON string
literal. mypy cannot catch that; the mistake is which encoder to call.

**The uncurated sentence is load-bearing.** 118 of the 119 genes published
carry no `LesionAssertion` (measured 2026-10-26 on a real build). Their pages
state that the atlas has not curated them rather than leaving the evidence
section silently absent, because an absent section is indistinguishable from
"the atlas looked and found nothing" and a reader deciding what a gene means
clinically must not have to infer which. Since D55 the sentence is on the card,
in the rail beside the chip, rather than a paragraph in the main column; its
wording must not deny the functional records the atlas may have curated for the
same gene, and `_card_curation` is where that wording and the reason for it are
recorded.

**Neither page kind states a bare grade.** `_SCOPE_SHORT` says in one line whose
classification it is and for what, and links `methods.html#scope`, where the
rule that admits a gene to this site at all and `render.grade_key`'s gloss of
each rung in the grading authority's own words are stated in full (D50,
2026-10-26); and `_graded_diseases` names the disease each
ClinGen panel actually graded -- in the browse table's `graded for` column and in
the gene page's rail, one line under the chip. Both read the same `validity`
mapping the bundles were written from, so the two pages cannot name different
diseases for one gene.

Every value interpolated into a page goes through `render.py`, which escapes
text, links and attributes -- except the bare interpolations in `_rail` (three,
plus its identity block) and in the card helpers it calls, which escape
themselves and say so. Nothing curated or
mirrored is reviewed the way page copy is, and this is the one artifact kind
where a `<` in a disease label is the difference between a string and a script
tag. A mirrored disease label now reaches two more places for that reason: the
rail row and the browse cell, both of which go through `definition_list` and
`data_table` respectively.

The research-use notice is `render.document`'s, not this module's, so it is on
every gene page and on the browse page without either builder mentioning it.

Nothing here reads a clock and every value comes from the caller, so two builds
of one commit render byte-identical pages. `build_gene_pages` iterates its
mapping in sorted order for the reason every other loop in this build does.
"""

from __future__ import annotations

import html
import math
import re
from collections.abc import Collection, Iterable, Mapping, Sequence
from typing import Final, NamedTuple
from urllib.parse import quote

from chd_atlas.build.burden import BurdenRow, CombinedRow, Overlap, _sort_key, shared_cohorts
from chd_atlas.build.cell_type_chart import cell_type_block
from chd_atlas.build.cell_types import EMPTY_CELL_TYPE_EXPRESSION, CellTypeExpression
from chd_atlas.build.charts import (
    LogScale,
    Shape,
    coordinate,
    glyph,
    marker,
    polyline,
    svg_figure,
)
from chd_atlas.build.concordance import (
    FamilyState,
    deciding_rows,
    family_state,
    is_upward,
    survives_own_correction,
)
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.gene_literature import (
    EMPTY_GENE_LITERATURE,
    GeneLiterature,
    SearchRecord,
    pmid_number,
)
from chd_atlas.build.identity import EMPTY_IDENTITY, GeneIdentity, external_links
from chd_atlas.build.paths import GENE_INDEX_PAGE, gene_bundle_path, gene_page_path
from chd_atlas.build.profiles import (
    EMPTY_EXPRESSION_PROFILE,
    DatasetProfileEntry,
    ExpressionProfile,
    PhaseInfo,
    PhaseOutcome,
    ProfileGap,
    StageProfileEntry,
    TissueProfileEntry,
)
from chd_atlas.build.render import (
    EVIDENCE_POWER_CAVEAT,
    FILTER_SCRIPT,
    GRADE_CHIP_KIND,
    GRADE_GLOSS,
    REGISTRY_QUALIFIER,
    Cell,
    Link,
    Markup,
    Row,
    chip,
    data_table,
    definition_list,
    document,
    evidence_legend,
    review_status_note,
)
from chd_atlas.build.validity import (
    GeneValidity,
    ValidityRecord,
    agreeing_records,
    agreeing_submitters,
    asserting_institutions,
    uncurated,
)
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.literature import Publication
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile
from chd_atlas.models.review import ReviewStatus
from chd_atlas.models.scope import ScopeAuthority, ScopeEntry
from chd_atlas.orphanet_scope import Warrant
from chd_atlas.tables import CHD_CASE_SPLITS
from chd_atlas.vocab import (
    CONTESTED,
    AtlasCuration,
    Classification,
    FindingDesign,
    ValiditySource,
    ValidityState,
)

# What the published set *is*, in one line, on both page kinds that show it: the
# short form of the scope rule, which lives in full on `methods.html#scope`
# (`methods.SCOPE_RULE`, D50, spec section 5's table) together with the grade
# key. Both page kinds carried the whole rule and the key until 2026-10-26; the
# history of why the rule exists at all moved with it.
#
# What the short form keeps is the half a reader needs on every page, because it
# is the half whose absence did the damage: whose classification this is (an
# external authority's, not the atlas's) and what it is for (the disease it
# names). A green `definitive` under a title reading "CHD Atlas" read as
# definitive for CHD, which measured 2026-08-04 was false for most of the first
# 23 genes.
#
# **It is hedged, because the plain form was false on 53 pages.** The spec's
# short form read "for the disease named, which is not the same as for CHD".
# On 53 of the 119 gene pages the panel graded congenital heart disease itself
# -- ISL1's card reads "definitive for congenital heart disease" one line above
# it (measured 2026-10-26) -- so the sentence denied the card. "Which for many
# genes here is a syndrome rather than CHD" is true on every page: it says what
# the disease often is, not what it is on this one, and the card beside it says
# that.
#
# Both page kinds sit in `genes/`, hence `../`. One constant for both, for the
# reason the full rule was one: two copies of one claim are two things that
# drift.
_SCOPE_SHORT: Final = (
    '<p class="scope-rule">Published on an external authority&#x27;s classification for '
    "the disease it names, which for many genes here is a syndrome rather than CHD. "
    '<a href="../methods.html#scope">How a gene is admitted, and what each grade '
    "means</a>.</p>"
)

# **Two notices, because the sentence was true of 23 genes and is now true of
# 76.** Until 2026-08-06 every published gene carried a ClinGen grade, so "the
# classification above is an expert panel's" was unconditionally true. The
# widened gate admits 16 genes no panel has graded at all -- measured on the
# committed corpus, all 16 shipped that sentence, telling a reader a commercial
# laboratory's submission was a chartered panel's classification, one screen
# below a legend saying no panel had graded it. A sentence's truth conditions
# travel with the population it describes.
#
# `validity_state` is the discriminator rather than `headline_confidence is
# None`, because the two are equal on the committed corpus (16 and 16) and equal
# figures are one figure to every test.
#
# Reworded 2026-10-25 because the validity section became two groups (D51):
# "the classification above is an expert panel's" then sat under the
# "Other curators (GenCC)" table on 59 pages. A sentence's truth conditions
# travel with its position.
#
# **Moved into the card on 2026-10-26 (D55)**, as one sentence each. They were
# two main-column paragraphs, `_NOT_CURATED` ("The atlas has not yet curated a
# lesion assertion for this gene. Each classification above is a ClinGen
# expert panel's or a named curator's ...") and `_NOT_CURATED_UNGRADED`, under
# the validity section; the card states the same claims beside the chip, and
# the main column no longer repeats them. `_card_curation` composes the
# sentence. Two claims had to survive the move, and each is narrower than the
# obvious wording:
#
# * The submitter-admitted sentence keeps its scope. "No ClinGen expert panel
#   has graded this gene" is false unscoped for ACTC1 (Definitive for
#   hypertrophic cardiomyopathy), ELN (Definitive for cutis laxa), FLNA, SMAD6
#   and MNS1 (Disputed for primary ciliary dyskinesia) -- five of the 21
#   admitted on agreement carry a ClinGen record for a disease out of scope,
#   measured 2026-10-26 against `mirrors/clingen_gene_validity.tsv`. Those
#   records reach no page, so the page could not even show the reader the
#   record that contradicted it. The scope clause -- "for any of the diseases
#   in scope here" -- is the grade key's `ungraded` gloss's (`render.GRADE_GLOSS`,
#   on `methods.html#scope`) and the gene index's em-dash sentence's, word for
#   word, so the three scope the claim alike. Only that clause is shared: the
#   rest of each sentence is its own (this one says the submitters
#   "independently assert it", the gloss that they "classify it Limited or
#   better").
# * **"In scope here" names the list, not a property of a disease.** The clause
#   read "for a disease an external authority treats as congenital heart
#   disease, or in which Orphanet records a structural heart malformation in at
#   least 30% of cases" until 2026-10-26, which reads as a sufficient condition
#   and was false on FLNA's card: ClinGen grades FLNA Definitive for
#   periventricular nodular heterotopia (MONDO:0020341), which SCP007 reports
#   would qualify through Orphanet (ORPHA:98892) and which no curator has
#   nominated. A disease is in scope only once it is in
#   `curation/chd_scope.yaml`; `methods.SCOPE_RULE` says how a term gets there.
# * The opening names the lesion assertion as what is absent whenever the
#   atlas has curated functional records for the gene; see `_card_curation`.
# * **"the submitters it is admitted by", not "named above"** (review,
#   2026-10-27). Above this sentence the card named only the `admitted by`
#   row's submitters until R6 put a no-association line between them: on
#   GDF1 it names Illumina, whose `No Known Disease Relationship` asserts
#   nothing, so "the submitters named above independently assert it" took
#   in the one submitter that dissents. The clause now names the admitting
#   submitters by their role, wherever the card puts them.
_CURATION_BY_PANEL: Final = (
    "every classification on this page is a ClinGen expert panel's or a named curator's."
)

_CURATION_BY_AGREEMENT: Final = (
    "and no ClinGen expert panel has graded this gene for any of the diseases in scope "
    "here: it is published because the Gene Curation Coalition submitters it is admitted by "
    "independently assert it."
)

# The `sources disagree` chip's tooltip: what `has_source_discordance` measures
# (`validity._has_source_discordance`), which its two words do not say. The flag
# is per gene by the owner's decision (2026-10-25), and on every page carrying
# it today ClinGen and GenCC agree on each gene-disease pair they both classify
# (BRAF, KRAS, MAP2K1, PTPN11, RAF1 and SOS1, measured 2026-10-26): it fires
# because one source disputes the gene for Costello or cardiofaciocutaneous
# syndrome while the other supports it for another disease (Noonan syndrome
# among them on all six). "Sources disagree" alone reads as two authorities
# contradicting each other about one disease, so the tooltip says "not
# necessarily the same disease". "The other" because there are exactly two
# sources, and the second is GenCC's submissions as mirrored -- ClinGen's own
# copies included, which `methods.html#validity` states.
_SOURCE_DISCORDANCE_TITLE: Final = (
    "One source, ClinGen's expert panels or the submissions to the Gene Curation Coalition, "
    "disputes or refutes this gene for a disease in scope while the other classifies it "
    "Limited or better for one, not necessarily the same disease; see the validity table."
)

_EM_DASH: Final = "—"

_VALIDITY_HEADERS: Final = (
    "curated by",
    "disease",
    "inheritance",
    "classification",
    "date",
    "report",
)

# `_orphanet_notes`'s scope table: one row per disease a gene carries that is in
# scope because Orphanet records a structural heart malformation at or above the
# frequency floor. Display-only wording; see `orphanet_scope.warrant`.
_ORPHANET_SCOPE_HEADERS: Final = (
    "disease",
    "why it is in scope",
    "structural heart malformations Orphanet records",
)

# ClinGen's mode-of-inheritance codes in GenCC's words, so both groups of the
# validity section read one vocabulary (D51). Measured 2026-10-25 against the
# committed mirrors: ClinGen uses exactly these six codes (AR 1,833 rows, AD
# 1,471, XL 207, SD 59, MT 54, UD 29), and GenCC spells the first five exactly
# so. `UD` is ClinGen's "undetermined", which GenCC records as `Unknown`. The
# bundle's `moi` stays verbatim; this is display only.
_INHERITANCE_LABEL: Final[dict[str, str]] = {
    "AD": "Autosomal dominant",
    "AR": "Autosomal recessive",
    "XL": "X-linked",
    "SD": "Semidominant",
    "MT": "Mitochondrial",
    "UD": "Unknown",
}

# The calendar date every timestamp either mirror publishes begins with.
# Measured 2026-10-25, by replacing digits with 9 and counting. ClinGen's
# `classification_date`: `2024-11-07T17:00:00.000Z` (3,588 rows) and a zone-less
# `2017-11-21T00:00:00` (65 rows). GenCC's `submitted_on`: space-separated
# `2018-05-31 20:04:32` (20,543 rows), `2025-01-17T00:00:00.000000Z` (9,857) and
# the millisecond `...00.000Z` (10).
_CALENDAR_DATE: Final = re.compile(r"\d{4}-\d{2}-\d{2}")

_GOLD_STANDARD_BADGE: Final = (
    ' <span class="gold-standard" title="Written by a ClinGen expert panel, the primary '
    'authority for gene&ndash;disease validity">Gold standard</span>'
)

_NO_PANEL_RECORD: Final = (
    "<p>No ClinGen expert panel record for this gene and an in-scope disease is mirrored here.</p>"
)

_NO_OTHER_CURATOR: Final = (
    "<p>No Gene Curation Coalition submitter other than ClinGen is mirrored here for this "
    "gene and an in-scope disease.</p>"
)

# Rendered when `_repeats_the_panel` folded at least one copy away, or when any
# listed record is ClinGen's own submission -- so a ClinGen row under "Other
# curators" is always explained, and the notice never describes an omission that
# did not happen.
#
# "keeps every copy", not "every record", since 2026-10-26: beneath an intro
# saying the table may not be every classification the gene has, "every record"
# read as though the JSON held the out-of-scope records. It holds none, and
# every copy (591 bundle records, 0 out of scope, all 121 copies).
_CLINGEN_COPY_NOTICE: Final = (
    "<p>ClinGen also submits its own classifications to the Gene Curation Coalition; "
    "such a copy is listed here only where it differs from the panel record above, "
    "carries no date, or has no panel record to match. "
    "This gene's JSON keeps every copy.</p>"
)

_EVIDENCE_HEADERS: Final = ("class", "strength", "summary", "publication")

# Three of these read as something else on the same page, which review measured
# on 2026-08-05:
#
# * `cohort` held `all cases / syndromic / non-syndromic` -- strata of the case
#   set -- six lines below a provenance paragraph using "cohort" for the sample
#   *collections* (`cnchd`, `ddd`, `ukbb`), which is also what
#   `curation/cohorts.yaml` is named for. `case group` is the stratum.
# * `variants` held `loss-of-function / synonymous`; in a row whose other four
#   cells are numbers, that header invites reading it as a count. It became
#   "variant class" on 2026-08-05 and that was worse: this column renders
#   `consequence_class`, while the method line above the table renders the real
#   `variant_class` ("SNVs and indels"). One phrase meant two things on one page
#   and was wrong in the more prominent place. Raised by review on #17.
#   `test_the_variant_column_is_headed_for_the_column_it_renders` pins it.
# * `cases` and `controls` hold "carriers / total", and nothing said so. The
#   denominator must travel with the numerator -- six carriers is a different
#   claim in 3,876 people than in 45,082 -- and the header is where a reader
#   looks for what a cell means. `methods.html#rare-variants` states the rule
#   as well (`_BURDEN_PREAMBLE` until 2026-10-26).
_BURDEN_HEADERS: Final = (
    "case group",
    "consequence",
    # Not "cases (carriers / n)". The unit is a property of the *row*, not of
    # the column -- one study's table can hold a de novo arm counting mutations
    # and a case-control arm counting alleles -- so it is rendered per cell by
    # `_count` and the header names only what the column is.
    "cases",
    "controls",
    "effect",
    "p",
    "corrected p",
)

# The word for each count unit, as (numerator, denominator suffix). The
# denominator of a de novo row is *trios*, not alleles or people, because a trio
# is what it takes to call one -- so this pairs both halves rather than labelling
# the numerator and leaving the reader to assume the other matches.
_COUNT_UNIT_LABEL: Final[dict[str, tuple[str, str]]] = {
    "individuals": ("carriers", ""),
    "alleles": ("alleles", ""),
    "de_novo_mutations": ("de novo", " trios"),
}

_SINGULAR_COUNT_NOUN: Final[dict[str, str]] = {"carriers": "carrier", "alleles": "allele"}

# How a published correction is named in the cell, so a reader is never left to
# guess what "0.99" was corrected against.
_ADJUSTMENT_LABEL: Final[dict[str, str]] = {
    "bonferroni": "Bonferroni",
    "familywise_permutation": "family-wise",
    "benjamini_hochberg": "Benjamini-Hochberg",
}

_STRATUM_LABEL: Final[dict[str, str]] = {
    "all": "all cases",
    "syndromic": "syndromic",
    "nonsyndromic": "non-syndromic",
}

# The four states, as a reader meets them. `not_tested` is deliberately not a
# synonym for `no result`: KDM6A is ClinGen definitive -- it causes Kabuki
# syndrome -- and shows nothing in either dataset that tested it, so a page that
# merged the two would report "0 of 3" beside a green definitive chip and read
# as the data contradicting the classification.
# The design each comparator names, for the matrix column heads. Reuses the
# discriminator the whole burden schema turns on rather than inventing a
# parallel vocabulary for the page.
# What the glyphs mean, and the sentence that stops the column reading as a
# verdict. **The caption is not optional.** KDM6A is ClinGen definitive -- it
# causes Kabuki syndrome -- and shows nothing in either dataset that tested it,
# so its row reads "0 of 2 tested" beside a green `definitive` chip. Without the
# caption a clinical geneticist reads that as the data contradicting the
# classification; with it, they read it as the power statement it is.
# The matrix's own key and caption. **The browse page carried these and the gene
# page did not**, which is the wrong way round: a reader arriving from a search
# lands on the gene page, and KDM6A's matrix there is entirely hollow beside a
# green `definitive` chip with nothing to explain it. Found by review before the
# deploy, on the page where the omission would have done the damage.
_MATRIX_LEGEND: Final = (
    evidence_legend(swatches=True)
    + '<p class="strip-legend">Rows are <strong>independent cohort families</strong>, not '
    "studies: two papers that share a sample collection, or draw on collections this atlas "
    "records as overlapping, can describe the same people and appear once. "
    + EVIDENCE_POWER_CAVEAT
    + "</p>"
)

_STRIP_LEGEND: Final = (
    evidence_legend(swatches=False)
    + '<p class="strip-legend">One glyph per <strong>independent cohort family</strong>, not '
    "per study: two papers that share a sample collection, or draw on collections this atlas "
    "records as overlapping, can describe the same people and count once. "
    + EVIDENCE_POWER_CAVEAT
    + " The tally reads <em>enriched</em> of <em>tested</em>: the first "
    "number counts datasets that found something, the second counts datasets that looked. "
    "Neither is a verdict on the gene.</p>"
)

_COMPARATOR_LABEL: Final[dict[str, str]] = {
    "control_cohort": "case-control",
    "mutation_model": "de novo, trios",
    "none": "case series",
}


# What kind of analysis each arm of a study ran, named as the literature names it
# (D49): a case-control burden test, a CNV burden test, or a de novo enrichment
# against a mutation model -- which is not a burden test at all. A combination
# nothing has taught this table names itself rather than nothing, the rule
# `_count` applies to an unknown unit.
def _analysis_label(comparator: str, variant_class: str) -> str:
    if comparator == "mutation_model":
        return "de novo enrichment (trios vs mutation model)"
    if comparator == "control_cohort":
        if variant_class.startswith("cnv_"):
            return "CNV burden (case-control)"
        if variant_class == "snv_indel":
            return "case-control burden test"
        return f"case-control burden test ({variant_class})"
    if comparator == "none":
        return "case series (no comparison group)"
    return f"{comparator} ({variant_class})"


# How a combination method is named on the page, what its second statistic is
# called, and the caveat that must sit beside it. An unknown method renders its
# raw token and no caveat rather than nothing -- the rule `_count` applies to an
# unknown unit -- because a missing name is worse than an unpolished one.
_COMBINATION_LABEL: Final[dict[str, tuple[str, str, str]]] = {
    "joint_local_fdr": (
        "the joint local FDR method",
        "joint local FDR",
        "A joint local FDR is <strong>not a corrected p-value</strong>: it can be "
        "smaller than the p beside it, so it is not comparable with the "
        "Benjamini-Hochberg column.",
    ),
}

_STATE_GLYPH: Final[dict[str, str]] = {
    "corrected": "full",
    "nominal": "half",
    "no_enrichment": "none",
    "not_tested": "untested",
}

# What a matrix cell *reads*, as against what it says on hover. The full
# sentence -- "this dataset did not test this gene" -- set the column width for
# the whole table and pushed the cells wide enough to crowd each other. The
# short form keeps the grid tight; `_STATE_TITLE` still carries the sentence on
# the cell's `title`, so nothing is lost to a reader who wants it.
_STATE_CELL: Final[dict[str, str]] = {
    "not_tested": "not tested",
    "no_enrichment": "tested",
}

_STATE_TITLE: Final[dict[str, str]] = {
    "corrected": "enriched, and survives this study's own correction",
    "nominal": "enriched nominally; not after correction, or no correction published",
    "no_enrichment": "tested, no enrichment detected",
    "not_tested": "this dataset did not test this gene",
}

_CONSEQUENCE_LABEL: Final[dict[str, str]] = {
    # Names its two components rather than saying "damaging", so a reader who
    # meets this row first is told what it is the union of before they reach the
    # two rows below it. `_composite_note` says the rest.
    "damaging": "damaging (LOF + missense)",
    "lof": "loss-of-function",
    "missense_damaging": "missense (damaging)",
    "missense_all": "missense (all)",
    "all_coding": "all coding",
    "synonymous": "synonymous",
}

# Spelled out beside every number, because the abbreviation is the difference
# between two claims. See `_effect`.
_MEASURE_LABEL: Final[dict[str, str]] = {
    "odds_ratio": "OR",
    "enrichment_ratio": "enrichment",
    "rate_ratio": "rate ratio",
}

# What the value 1 *means*, which is not the same sentence for every measure --
# and `_forest_caption` printed one of them on all of them.
#
# Measured 2026-08-20 on the built site: 112 of the 137 forest panels, across
# 88 gene pages, read "Axis: **OR** ... the vertical line is **1**, no
# enrichment". For an odds ratio 1 is *no association*; "no enrichment" is the
# de-novo-against-expectation vocabulary belonging to the 25
# `enrichment_ratio` panels. That is exactly the merge `effect_measure` exists
# to prevent (see `_effect`), committed in the one sentence whose job is
# telling a reader how to read the axis.
#
# Keyed on the same token as `_MEASURE_LABEL` and consulted in the same place,
# so the label and its gloss cannot drift apart again: a fourth measure that
# lands here with a label and no gloss gets `_NULL_LINE_DEFAULT`, which is
# true of any ratio, rather than one vocabulary's words applied to another's
# number. `rate_ratio` carries 0 rows today and is written anyway, for the
# reason `_effect` has no branch that omits its measure.
_NULL_LINE_GLOSS: Final[dict[str, str]] = {
    "odds_ratio": "no association",
    "enrichment_ratio": "no enrichment",
    "rate_ratio": "no difference in rate",
}

# Deliberately not a blank and not a refusal. `_effect` renders an unknown
# measure's raw token rather than nothing, on the rule that a missing label is
# worse than an unfamiliar one; the same rule applies to the gloss, and every
# member of `EffectMeasure` is a ratio, for which 1 is the arms agreeing.
_NULL_LINE_DEFAULT: Final = "no difference between the arms"

# The partition columns that reached no page until 2026-08-05. Without them a
# section headed "Rare variant burden" (now "Rare-variant association", D49)
# never says what *rare* means and a
# case-control count of rare inherited variants reads exactly like a trio's de
# novo count.
#
# This said "the four partition columns... each is populated on all 187 published
# rows" and named `lesion_group` among them. **`lesion_group` is null on all 290
# published rows and always has been** -- no curated study stratifies by lesion
# group -- and `maf_max` is null on 45 of them, the de novo rows that applied no
# frequency filter. So the claim was false for two of its four columns.
# `_LESION_LABEL` stays empty rather than being deleted: the column is live in
# the schema, and the first study that uses it needs that map filled.
_VARIANT_CLASS_LABEL: Final[dict[str, str]] = {
    "snv_indel": "SNVs and indels",
    "cnv_deletion": "CNV deletions",
    "cnv_duplication": "CNV duplications",
    "sv_other": "other structural variants",
}

_ORIGIN_LABEL: Final[dict[str, str]] = {
    "de_novo": "de novo only",
    "inherited": "inherited only",
    "any": "any inheritance (not a de novo test)",
    # Says what was excluded and what is unknown, because that is the whole
    # content of the value: de novo variants were removed, and most of what
    # remains was never phased. "inherited" would claim transmission the study
    # did not observe for the singleton majority.
    # No internal semicolon: `_method_line` joins its parts with "; " and the
    # origin labels within a part with ", ", so a label containing one produced
    # "de novo excluded; transmission otherwise known or unknown, de novo only"
    # -- which reads as three clauses at two nesting levels and parses as none.
    "transmitted_or_unphased": "de novo excluded (transmission otherwise known or unknown)",
}

_LESION_LABEL: Final[dict[str, str]] = {}

_TEST_LABEL: Final[dict[str, str]] = {
    "fisher_exact": "Fisher's exact test",
    "poisson": "Poisson test",
    "binomial": "binomial test",
}

# Attribution, in the voice `_validity_section` and `landing.py` already use for
# republished content. The `<h2>` alone said "Rare variant burden" (now
# "Rare-variant association", D49) while its
# sibling says "Mirrored gene-disease validity", and the uncurated notice that
# sat above it (on the card since D55) scopes itself to *classifications* -- so
# nothing between the heading and the table said whose numbers these are.
#
# One line since 2026-10-26 (D50), linking `methods.html#rare-variants`, which
# states it in full. **Scoped to the study tables, not to "below".** The
# evidence matrix directly under this line shows tallies the atlas derives --
# "n = N" in a cell on 86 pages, "N rows" on 48 (measured on a build of
# 4233eca) -- so "every count and statistic below ... the atlas computes none
# of them" was false of the first thing below it, and had been since the matrix
# arrived. "exactly" went too: the tables round to three significant figures.
_MIRRORED_NOTICE: Final = (
    "<p>Every count and statistic in the study tables below is an upstream study's, "
    "republished as published; the atlas computes none of them. "
    '<a href="../methods.html#rare-variants">How to read these results</a>.</p>'
)

# Rendered only when a gene carries rows from more than one study. Rendered
# unconditionally it asserted "these cohorts overlap" in the present tense on a
# page showing a single table, sending a reader to look for a second study that
# is not there.
# Now the matrix's caption, so it is rendered on every page with a matrix rather
# than only where a gene carries two studies. That change made the old wording
# false: it asserted "these cohorts overlap" in the present tense, which on a
# page showing one study sent a reader hunting for a second -- the very defect
# the original conditional was added to fix. It is a statement of policy now,
# true whether the page shows one study or three, and `shared_cohorts` still
# names an actual overlap where one exists.
#
# The spec's approved short form since 2026-10-26 (D50, P7), made explicit about
# what would double-count: the spec's "those sharing a collection would count"
# left the pooling implicit. The full sentence, with how a gene page names an
# overlap, is on `methods.html#rare-variants`. "people", not the spec's
# "children" (2026-10-26): `shared_cohorts` unions case and control collections,
# and a shared collection may be a control cohort of adults -- UK Biobank's
# registry note says "~500,000 UK adults recruited at ages 40-69".
#
# The declared-overlap clause since review (2026-10-27): the reason it gave was
# a shared collection alone, and on the 72 pages carrying both PMID:34324492
# and PMID:42230622 the reason is a declared overlap (`Cohort.overlaps_with`,
# owner decision I3), for which "may describe partly the same people" is all
# the registry knows -- so "would" for a shared collection, "could" for the
# declared kind, as `methods.html#rare-variants` words the two.
_POOLING_NOTICE: Final = (
    "<p>The atlas computes <strong>no pooled statistic across studies</strong>: pooling "
    "two that share a sample collection would count the same people twice, and pooling "
    "two that draw on collections this atlas records as overlapping could.</p>"
)

# The matrix and the card each lead with one row, by two rules the owner kept
# (owner decision 2026-10-26, R4), so the caption says they differ rather than
# leaving a reader to find the two disagreeing. A cell names the most
# significant of the rows that give it its state (`_matrix_cell`: among
# `concordance.deciding_rows`, by corrected p where published and raw p
# otherwise); the card lists a study's surviving rows in its table's reading
# order (`_card_rare_variants`, `burden._sort_key`), in which the `damaging`
# union sorts before its components. Measured on a build of 3b766ba: of the 45
# card lines on 32 pages, 15 lead with a different row from their matrix cell,
# on 13 pages. The owner's draft said "its strongest supporting row": the cell
# names a surviving row in a `corrected` cell, ranked by p, not by effect --
# "most significant", the word `docs/data-api.md` uses for it.
#
# **Rendered only where the card lists a surviving row** (review,
# 2026-10-27): with the matrix, and only if some row of this gene survives
# its own study's correction (`concordance.survives_own_correction`, the
# predicate the card's list reads). It rendered with every matrix, 118 pages,
# and on the 86 whose card lists no row there is no lead for a cell to differ
# from; on 55 of those no cell names a row at all (measured on a build of
# ad03128). On the 32 where it renders, card and matrix can lead differently.
# And it claims the card's order only as far as it holds: the card groups a
# study's surviving rows into one line per analysis, ordered by analysis
# label, so "it keeps each study table's order" was true within a line and
# not across a study's lines; that a damaging union comes before its parts is
# true on every line.
_MATRIX_LEAD_NOTE: Final = (
    "<p>A cell names the most significant of the rows behind its state; the card lists a "
    "study's surviving rows with a damaging union before its parts, so the two can lead "
    "with different rows.</p>"
)

# `_BURDEN_PREAMBLE` and `_SYNONYMOUS_NOTICE` sat here until 2026-10-26: the
# rare-variant reading notes every gene page folded under "How to read these
# numbers". They moved, with their histories, to `methods._RARE_VARIANTS` (D50).


def _heart_expression_summary(
    profile: ExpressionProfile | None,
) -> tuple[float, int | None, str | None]:
    """Returns (max_tpm, best_percentile, best_stage) for embryonic/fetal heart."""
    if not profile or not profile.get("datasets"):
        return (0.0, None, None)
    max_tpm = 0.0
    best_pct = None
    best_stage = None
    for ds in profile["datasets"]:
        for st in ds.get("stages", []):
            for t in st.get("tissues", []):
                if t.get("tissue") == "heart":
                    val = float(t.get("median_abundance", 0.0))
                    if val > max_tpm:
                        max_tpm = val
                        placement = t.get("placement")
                        best_pct = placement.get("median_percentile") if placement else None
                        best_stage = st.get("stage")
    return (max_tpm, best_pct, best_stage)


def _heart_expression_cell(profile: ExpressionProfile | None) -> tuple[Cell, str, str]:
    max_tpm, best_pct, best_stage = _heart_expression_summary(profile)
    if max_tpm <= 0.0:
        return (_EM_DASH, "none", "0")
    level = (
        "high"
        if max_tpm >= 50.0
        else "moderate"
        if max_tpm >= 5.0
        else "low"
        if max_tpm >= 1.0
        else "absent"
    )
    if max_tpm < 1.0:
        content = (
            '<span class="expression-cell" '
            f'title="Peak at {html.escape(best_stage or "")}: {max_tpm:.1f} TPM">'
            '&lt;1 TPM <span class="expression-pct">(low)</span></span>'
        )
    else:
        top_val = max(1, 100 - best_pct) if best_pct is not None else None
        pct_str = f"top {top_val}%" if top_val is not None else ""
        pct_span = f' <span class="expression-pct">({pct_str})</span>' if pct_str else ""
        content = (
            '<span class="expression-cell" '
            f'title="Peak at {html.escape(best_stage or "")}: {max_tpm:.1f} TPM">'
            f"{int(round(max_tpm))} TPM{pct_span}</span>"
        )
    return (Markup(content), level, str(max_tpm))


def _pubmed(pmid: str) -> str:
    """The canonical article URL for a `PMID:` identifier.

    A literal rather than a mirrored field: `Publication` records a DOI and a
    PMCID but no PubMed URL, and this is the one link on the page the atlas
    constructs rather than republishes.
    """
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid.removeprefix('PMID:')}/"


def _graded_diseases(gene_validity: GeneValidity, grade: Classification | None) -> tuple[str, ...]:
    """The in-scope diseases a ClinGen panel graded this gene at `grade`.

    **Keyed on the gene's own headline, not on `DEFINITIVE`.** It was the latter
    until 2026-08-06, when the gate widened past ClinGen `Definitive` and 69 of
    the 92 published genes stopped having one: a function that only ever looks
    for definitive records returns nothing for them, and the caller renders an em
    dash beside a `limited` chip -- the disease the panel actually graded silently
    absent from the one row that exists to name it. That is this project's
    characteristic failure in page copy rather than in a join.

    `grade` is `GeneFacts.headline_confidence`, which `derive.gene_facts` derives
    from ClinGen records alone, so what comes back is always the disease behind
    the chip the reader is looking at. `None` -- a gene admitted on GenCC
    agreement, which no panel graded -- returns empty, and the caller names the
    submitters instead.

    A GenCC submitter's own grade is deliberately never here. It admits a gene
    only by agreement with a second submitter (D21) and is never the headline, so
    naming its disease would qualify a chip with a record that did not produce
    it.

    Returns labels rather than MONDO ids: this is the string a reader reads
    beside the grade, and the id is already on the gene page in
    `_validity_section`, which renders both.

    **Sorted, and de-duplicated through a `set` first.** A gene may hold more
    than one record at its headline grade -- two panels grading it for two
    in-scope diseases, or one disease reached through two rows -- and a `set` of
    strings iterates in an order that varies with `PYTHONHASHSEED`, which would
    make the browse page's bytes and therefore its checksum differ between two
    builds of one commit. Measured 2026-08-06 against the committed mirrors: all
    76 panel-graded genes have exactly one distinct label at their headline
    grade, so the multi-label case is latent and its ordering is guarded by
    fixture alone
    (`test_a_gene_definitive_for_two_in_scope_diseases_names_both_in_a_fixed_order`).
    That name was wrong here until 2026-08-06 -- it cited a
    `..._graded_for_two_...` test that has never existed, which is the one thing
    a cross-reference must not do: it sends the next reader looking for a guard
    they will not find, and deleting the real one then costs nothing.
    """
    if grade is None:
        return ()
    return tuple(
        sorted(
            {
                record.disease_label
                for record in gene_validity.records
                if record.source is ValiditySource.CLINGEN and record.classification is grade
            }
        )
    )


def _card_curation(fact: GeneFacts, literature: GeneLiterature = EMPTY_GENE_LITERATURE) -> str:
    """The card's atlas-curation sentence, and what it must not deny. D55.

    Four branches, in the order the code takes them:

    * **curated** -- "Curated by the atlas: 1 lesion assertion (septal)." with
      the functional-record count appended only when non-zero (TBX5 and NR2F2
      are the curated genes today; one assertion each). With approved literature
      findings the sentence adds "and 24 findings from 21 papers" (owner
      decision 2); a gene with only those reads "Curated by the atlas: 5
      findings from 4 papers." Papers are the distinct PMIDs among the
      findings. `literature` is the bundle's, which carries approved findings
      only, so a gene whose findings all await approval is not counted here;
    * **not curated, a ClinGen panel holds a record** -- every classification
      on the page is a panel's or a named curator's (97 genes);
    * **not curated, no panel record** -- no ClinGen expert panel has graded
      the gene for any of the diseases in scope here, and it is published
      because the submitters named in the card's `admitted by` row assert it
      (21 genes). Scoped, because unscoped it is false; and scoped to the
      list, because scoped to a property of a disease it was false too; see
      `_CURATION_BY_AGREEMENT`;
    * **functional records but no lesion assertion** (latent: the committed
      corpus has no `curation/functional/` directory) -- appends "The atlas
      has curated N functional record(s)" and opens on "No lesion assertion
      yet curated by the atlas" instead of "Not yet curated by the atlas".

    "Not *yet*" in both openings, because the state is `not_yet_curated`; it
    read "Not curated by the atlas" for one commit.

    **That last opening is the lesson this function inherits.** `atlas_curation`
    is derived from `LesionAssertion`s alone -- `curated` iff `assertion_count >
    0`, `derive.gene_facts` -- while functional records come from a population
    `referential.py` never requires an assertion to cite, so a gene can carry
    curated functional records and still resolve to `NOT_YET_CURATED`. The
    notice this replaced said "not yet curated evidence for this gene ...
    nothing on this page is the atlas's own assessment" until 2026-08-04, beside
    a rail counting `functional records | 3`: denying curated work on the page
    that renders it, this project's characteristic failure in prose instead of
    a dropped join. "Not yet curated by the atlas" followed by "the atlas has
    curated 1 functional record" would be the same defect inside one sentence.
    `atlas_curation`'s meaning is deliberately untouched; what was wrong was
    the sentence, not the field. Literature is the same case: a gene with
    approved findings and no assertion is `not_yet_curated` in the JSON and
    "Curated by the atlas" on the card, because the card must not deny work
    the page below it renders.

    **`validity_state` is the discriminator between the two uncurated
    sentences, not `headline_confidence is None`**, for the reason the notices
    it replaced used it: the two agree on every gene published today (98 and
    21), and equal figures are one figure to every test. They can differ: a
    gene admitted on agreement that ClinGen recorded `No Known Disease
    Relationship` for in scope has a panel record and no headline, and the
    panel sentence stays true of it where "no ClinGen expert panel has graded
    this gene" would not.

    It was the main-column paragraph `_not_curated`, rendered after the
    validity section, until D55 moved it here; the page now states it once.
    """
    findings = literature["findings"]
    if fact.atlas_curation is AtlasCuration.CURATED or findings:
        parts = []
        if fact.atlas_curation is AtlasCuration.CURATED:
            noun = "assertion" if fact.assertion_count == 1 else "assertions"
            groups = [group.value for group in fact.lesion_groups]
            scope = f" ({html.escape(', '.join(groups))})" if groups else ""
            parts.append(f"{fact.assertion_count} lesion {noun}{scope}")
        if findings:
            papers = len({str(finding["publication"]) for finding in findings})
            parts.append(
                f"{len(findings)} {'finding' if len(findings) == 1 else 'findings'} "
                f"from {papers} {'paper' if papers == 1 else 'papers'}"
            )
        functional = (
            f"; {fact.functional_count} functional record{'s' if fact.functional_count > 1 else ''}"
            if fact.functional_count
            else ""
        )
        sentence = f"Curated by the atlas: {' and '.join(parts)}{functional}."
        return f'<p class="card-curation">{sentence}</p>'

    # "Not yet": the state is `not_yet_curated`, and the atlas has not decided
    # against curating these genes.
    opening = (
        "Not yet curated by the atlas"
        if fact.functional_count == 0
        else "No lesion assertion yet curated by the atlas"
    )
    if fact.validity_state is ValidityState.EXPERT_CURATED:
        sentence = f"{opening}: {_CURATION_BY_PANEL}"
    else:
        sentence = f"{opening}, {_CURATION_BY_AGREEMENT}"
    if fact.functional_count:
        noun = "record" if fact.functional_count == 1 else "records"
        sentence += (
            f" The atlas has curated {fact.functional_count} functional {noun}, "
            "published in this gene's JSON."
        )
    # Nothing here is mirrored or curated text: counts, constants, and the
    # lesion groups' own enum values, escaped where they are interpolated.
    return f'<p class="card-curation">{sentence}</p>'


def _inheritance(moi: str) -> str:
    """A record's mode of inheritance as the page shows it.

    A code nothing has taught this table renders as itself -- never blank and
    never coerced to a neighbouring phrase -- for the reason `_count` renders an
    unknown count unit verbatim: a missing word reads as "not stated", which is
    a different claim from "stated in a code this page does not know".
    """
    return _INHERITANCE_LABEL.get(moi, moi)


def _date(timestamp: str | None) -> str:
    """The `YYYY-MM-DD` a mirror timestamp starts with; anything else verbatim.

    Truncated, not parsed: every shape either mirror publishes begins with the
    calendar date (`_CALENDAR_DATE`), and the time of day is upstream
    bookkeeping no reader acts on. A value that does not start with a date is
    shown as the authority wrote it rather than dropped.
    """
    if not timestamp:
        return _EM_DASH
    match = _CALENDAR_DATE.match(timestamp)
    return match.group(0) if match else timestamp


def _record_date(record: ValidityRecord) -> str | None:
    """The record's own date: ClinGen's `classification_date`, GenCC's `submitted_on`."""
    if record.source is ValiditySource.CLINGEN:
        return record.classification_date
    return record.submitted_on


def _repeats_the_panel(record: ValidityRecord, panel: Sequence[ValidityRecord]) -> bool:
    """True iff `record` is ClinGen's own GenCC copy of a panel record already shown. D51.

    ClinGen submits its classifications to GenCC under its own name, so the
    mirror holds each ClinGen verdict twice. Measured 2026-10-25: all 121 such
    copies among the 119 published genes match their direct record on disease,
    mapped grade, inheritance and calendar date -- the 9 whose
    `classification_term` differs say `Disputed Evidence` where the panel says
    `Disputed`, the same rung -- so on that corpus every copy is a duplicate.

    **Hidden only on an exact match, never on the submitter's name alone.** A
    copy that disagrees on any of the four is a different statement, and hiding
    it would make the page agree with itself by deleting the disagreement. A
    copy with no date matches nothing: two missing dates are not evidence that
    two records agree.

    Only ClinGen's own copies are candidates. Another submitter agreeing with a
    panel is corroboration and is always listed, the principle
    `test_a_clingen_record_marks_expert_curated_and_gencc_records_still_publish`
    records for the bundle -- which keeps every record, these copies included.
    Folding them is the page's choice alone.
    """
    if record.source is not ValiditySource.GENCC or record.submitter != "ClinGen":
        return False
    copied = _record_date(record)
    if not copied:
        return False
    return any(
        original.disease == record.disease
        and original.classification is record.classification
        and _inheritance(original.moi) == _inheritance(record.moi)
        and _date(_record_date(original)) == _date(copied)
        for original in panel
    )


def _report_link(url: str | None, sop: str | None = None) -> Cell:
    """A link only for an `http`/`https` report URL; an em dash for anything else.

    Neither mirror schema constrains this column, and the failure is not
    hypothetical: `mirrors/gencc_submissions.tsv` carries one row -- SHOX,
    `MONDO:0009588` Langer mesomelic dysplasia, submitted by Ambry Genetics --
    whose `report_url` is the literal string `Pseudoautosomal region, recessive`
    (measured 2026-08-04, the only such row of the file). That disease is out of
    CHD scope today, so the row reaches no page; the same column would carry a
    `javascript:` URL into an `href` just as readily, and `render.py` escapes an
    attribute value without sanitising the scheme it names.

    Fixed here rather than in `tables.py` on purpose. The mirror is a faithful
    copy of what upstream published and has to stay one -- rejecting the row
    would drop a real GenCC submission over a field no reader needs -- so the
    scheme is decided at the render boundary, where a URL becomes a link.

    A guard on a bypassed gate, in the same idiom as `encode_json`'s
    `allow_nan=False` and `variants.py`'s shard checks: nothing upstream
    promises this, so the page refuses rather than publishes.

    `sop` labels a ClinGen report with the framework version it was curated
    under (`SOP9 report`), which is where the page states the SOP since D51
    folded its column away. Without a usable URL the SOP still renders, as
    text: a ClinGen classification shown without its SOP is an unqualified
    claim (`docs/data-api.md`).
    """
    if url and url.startswith(("http://", "https://")):
        return Link(text=f"{sop} report" if sop else "open", href=url)
    return sop or _EM_DASH


def _and_join(items: Sequence[str]) -> str:
    """ "a", "a and b", "a, b and c" -- for names in a sentence, in the order given."""
    if len(items) <= 2:
        return " and ".join(items)
    return f"{', '.join(items[:-1])} and {items[-1]}"


def _row_label(row: BurdenRow) -> str:
    """A burden row's result, named as the table names it: consequence, then partitions.

    The stratum where it is not the whole cohort, and the lesion group where
    there is one -- null on every published row today, and rendered anyway so
    the first per-lesion row is never shown as if it were the whole cohort's.
    """
    label = _CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class)
    if row.cohort_stratum != "all":
        label += f", {_STRATUM_LABEL.get(row.cohort_stratum, row.cohort_stratum)}"
    if row.lesion_group is not None:
        label += f", lesion group {_LESION_LABEL.get(row.lesion_group, row.lesion_group)}"
    return label


# Appended to a card line from a study by an author of this atlas, which the
# study's block declares (`_disclosure`): a result that reaches the card must
# not shed the declaration on the way.
_OWN_LAB_MARK: Final = " · study by an atlas author"


def _correction_scope(row: BurdenRow, publication: Publication | None) -> str:
    """What a card line's row survived, where its study says (owner decision 2026-10-26, R3).

    "Surviving its own study's correction" read as genome-wide on every line,
    and for PMID:40127276 it is not: that study's q is a Benjamini-Hochberg
    correction over its 248 panel genes (`Publication.correction_genes`,
    measured against Dataset S6), and the same study draws a genome-wide line
    of its own at p < 2.6e-6 (`Publication.genome_wide_threshold`, from its
    Fig. 1 legend; `_genome_wide`). Measured on a build of e0393f4: 67
    non-synonymous rows on 25 published genes survive at q < 0.05, and 25 of
    them, on 13 genes, are below the study's genome-wide line -- two strengths
    of finding the card printed identically. TBX5's case-control damaging row
    (p 1.34e-4, q 0.00369) is the first kind and its de novo damaging row (p
    1.07e-8) the second.

    So a line's row gains its corrected p and the panel it was corrected over.
    Only the line's own row: the rest of the group are named, not stated, and
    carry no corrected p. The statistic is named by `_statistic`, so it is `q`
    only where the correction is a false-discovery rate, as on the matrix and
    the forest panel. Renders only where the field is curated, so every other
    study today renders exactly as before.
    """
    if publication is None or publication.correction_genes is None:
        return ""
    if row.pvalue_adjusted is None:
        return ""
    return f" · {_statistic(row)} within its {publication.correction_genes:,}-gene panel"


def _chd_only(rows: Sequence[BurdenRow]) -> list[str]:
    """What a study's own table says of the CHD cases among a card's carriers (R2).

    PMID:34324492's case set combines CHD with a series of 777 sporadic
    thoracic aortic aneurysm probands, and the study's summary caveat on the
    card says so. For the genes its S8 Table counts, `n_chd_case_carriers`
    says how many of a row's case carriers are CHD cases. Where that is all of
    them: "For this gene, the study's S8 Table counts every carrier among its
    CHD cases, so none is an aneurysm proband." The second clause is what
    stops the first reading as a remark on how S8 counts, true of any gene
    (review, 2026-10-27). It is the consequence, not a second figure read off
    S8: the table's bands partition its carriers (measured 2026-10-27, "only
    CHD" plus "only TAA" is "CHD+TAA" for all 528 genes, a blank read as
    none), and of the 10 published genes the line shows, S8's TAA-only `case`
    is 0 for TBX1, EFNB2, SOX7 and ETS1 and blank for the other six -- so the
    sentence never says S8 "counts none" there. It does not retire the
    caveat: the line's denominator is the case set the study analysed,
    aneurysm probands among it (S8 has a TAA-only carrier for 144 of its 528
    genes).

    Where it is fewer, the sentence gives both counts, never "every", and no
    aneurysm count: the mirror carries none, and a difference of two counts is
    not a figure S8 prints. Measured on a build of this tree (2026-10-27): all
    10 published genes with a card line from that study carry the column and
    are equal (ELN, TBX1, GATA4, TAB2, EFNB2, SOX7, ROBO4, ETS1, NKX2-5,
    NR2F2); the one row where it is smaller, ANKRD11's 12 of 13, is on no
    published page, so that form is reachable from a fixture only.

    `rows` are the rows the card's lines show, one per line. The table, and
    the case the second clause rules out, are the ones `CHD_CASE_SPLITS`
    records for the row's study: BUR024 refuses the column on a row of any
    study it does not list, so a row that reaches here without one means the
    gate was bypassed, and that raises rather than names a table the study
    never had. A row with no case carrier says nothing,
    since "every carrier" of none is vacuous. Deduplicated in order, so two
    lines of one study cannot repeat a sentence.
    """
    sentences: list[str] = []
    for row in rows:
        chd = row.n_chd_case_carriers
        if chd is None:
            continue
        split = CHD_CASE_SPLITS.get(row.study)
        if split is None:
            raise ValueError(
                f"a row of {row.study} carries n_chd_case_carriers, which only a study in "
                f"CHD_CASE_SPLITS may (BUR024); the validation gate was bypassed"
            )
        if row.n_case_carriers == 0:
            continue
        sentence = (
            f"For this gene, the study's {split.table} counts every carrier among its CHD "
            f"cases, so none is {split.other_case}."
            if chd == row.n_case_carriers
            else f"For this gene, the study's {split.table} counts {chd:,} of the "
            f"{row.n_case_carriers:,} carriers among its CHD cases."
        )
        if sentence not in sentences:
            sentences.append(sentence)
    return sentences


def _genome_wide(row: BurdenRow, publication: Publication | None) -> str:
    """`p <raw p>` where a row clears its study's own genome-wide line, else "" (R3).

    Strictly below, as the paper's "P-value < 2.6 x 10^-6" reads; formatted as
    the study tables format a p. Marked on **every** surviving row of a card
    line, the named ones too, not only the line's own row: measured on a real
    build, marking the line's row alone left 15 of the 25 genome-wide rows
    unmarked and 6 of the 13 genes with one -- FLT4's case-control
    loss-of-function row (p 1.75e-15) sat named under a damaging row that
    survives within the panel only. Independent of `correction_genes`: a study
    may state a genome-wide line without correcting over a panel.
    """
    if publication is None or publication.genome_wide_threshold is None:
        return ""
    if row.pvalue is None or not row.pvalue < publication.genome_wide_threshold:
        return ""
    return f"p {row.pvalue:.3g}"


def _card_rare_variants(
    rows: Sequence[BurdenRow],
    families: tuple[frozenset[str], ...],
    publications: Mapping[str, Publication],
) -> str:
    """The card's rare-variant block: what survives its own study's correction (D56).

    **One line per (study, analysis) with a surviving row**, never one per row:
    KMT2D has five surviving rows over two analyses of one study, and five
    lines in a 15rem column would bury the rest of the card. The line shows the
    group's first surviving row in the table's own reading order
    (`burden._sort_key`, which puts the `damaging` union -- PMID:40127276's
    primary analysis -- before its components), then **names** the group's
    other survivors in the same order: "also surviving: loss-of-function;
    missense (damaging)". It counted them ("2 more survive") until 2026-10-26,
    which beside the `damaging` union counted the union's own components as
    further findings; the union's label already says what it is the union of.
    **Not the "most significant" row**: a line keeps its study table's reading
    order. The evidence matrix's cell is where a row is chosen by p
    (`_matrix_cell`, among the rows behind its state), and its caption says
    the two can lead with different rows (`_MATRIX_LEAD_NOTE`). This read
    "choosing one by p-value would be a ranking, and the atlas authors none
    (D12)" until review (2026-10-27), beside a matrix that chooses by p.
    Groups in study-id order, then analysis label, the order the section's
    own subtitles use.

    "Survives" is `concordance.survives_own_correction`, the predicate the
    evidence matrix's `corrected` tier reads, so the card and the matrix on the
    same page cannot disagree about a row (CLAUDE.md section 4.28).

    **A study that published no correction for this gene gets its own line**,
    whether or not another study's results survive: "<study> publishes no
    corrected p; judge its results against its scan of 138,609 tests in its
    block". Decided per gene, from its rows -- any `pvalue_adjusted` means the
    study published one. PMID:42230622 publishes none on any row and tests 113
    published genes, 24 of them alone (measured 2026-10-26); its block tells a
    reader to judge raw p against its scan ("a Bonferroni threshold over that
    many is p < 3.6e-07"), which SMAD4's p of 2.47e-10 clears. The card read
    "none survives its own study's correction" of it -- a correction that study
    never made. So
    the closing sentence now speaks only of "a correction its study published",
    and the line names the scan size from `Publication.tests_reported` where
    recorded and none where not. No threshold is applied and no statistic
    computed (D12/D33): the reader is pointed at the count the study made.

    Counts in their own unit through `_count`, and a de novo row against its
    expectation through `_controls`, so the card says "24 alleles / 20,151"
    where the table does (D46). A line from a study by an atlas author carries
    `_OWN_LAB_MARK`, as its block carries the declaration.

    **Each result line ends with a link to its study's block** (`_study_anchor`,
    the `id` on that study's `<h3>`): `caveats` where any of the group's
    surviving rows carries a `method_note`, `details` otherwise, with an
    accessible name naming the study (WCAG 2.4.9: the visible word alone does
    not say whose caveats). The counts on a line can need the note to be read
    right -- PMID:34324492's case set includes 777 thoracic aortic aneurysm
    probands, PMID:40127276's case-control arm counts a proband with two
    qualifying variants twice -- and the full note is footnoted in the block;
    only a curated summary of it follows the lines (below). Measured
    2026-10-26: every surviving row of both studies carries one, so all 45
    card lines read `caveats` and `details` is reachable from a
    fixture only. A no-correction line links its study's name instead.

    **Beneath the lines, a study's curated `Publication.summary_caveat`**,
    once per study with a line and in the lines' study order: the one caveat
    the curator says changes what that study's counts mean, so it travels with
    them rather than staying one click away. Measured 2026-10-26: only
    PMID:34324492 carries one (the 777 thoracic aortic aneurysm probands in its
    case set, and a carrier being anyone with a deletion overlapping the
    gene), and it reaches every card that shows its counts.

    With no surviving row it counts the **cohort families** with any row for
    the gene -- independent datasets, never studies, because two papers on one
    sample collection describe the same people -- and, where some study
    published a correction for the gene, says none survives one; where none
    did, it says only how many tested it, and the no-correction lines say the
    rest.

    **That line named no study, so it shed the declaration** (owner decision
    2026-10-26, R7): MYH11, HRAS, KLF13 and CFC1 are tested by PMID:34324492
    alone, nothing survives, and their cards said nothing of its being by an
    atlas author. Where any study that tested the gene is one, the line says
    so: "(a study by an atlas author)" where one study tested it, "(each
    study by an atlas author)" where several did and all are, and "(k of the
    n studies by an atlas author)" where only some are -- counted in studies,
    because the line's own count is of datasets and one family may hold both
    kinds. Measured on a real build: all 86 such lines carry it (24 + 4 of
    the first form, 48 of the second, 7 + 3 of the third). Where a
    no-correction line follows for the same study, it carries the mark too
    (24 cards); each line keeps its own declaration. With no row at all it says no
    study here reports a result, and not that no study *tested* the gene: a
    supplement that omits a gene (CRIPTO, in PMID:42230622) does not say the
    study never ran it. The power caveat follows whenever the gene has rows,
    so "none survives" never stands alone beside a definitive chip.

    **The count carries the front page's qualifier** (owner answer
    2026-10-27), `render.REGISTRY_QUALIFIER`, "(as far as the cohort registry
    records)": the families are joined only through collections the registry
    names or declares overlapping, so this count, like the front page's, sees
    no overlap the registry does not record. It follows the count directly
    and precedes the declaration: placed after it, it would read as
    qualifying who wrote the studies, which the registry does not record at
    all. Measured on a real build the day it was added: all 86 lines carry
    it, in these five forms: "Tested in 1 independent dataset (as far
    as the cohort registry records) (a study by an atlas author)." on 24, the
    same ending "; none survives a correction its study published." on 4,
    "(each study by an atlas author)" on 48, and two datasets with "(1 of the
    2 studies ...)" on 7 and "(2 of the 3 studies ...)" on 3.

    `families` is the corpus-wide grouping `build_site` derives once; a study
    of this gene's rows that no family holds means it was derived over another
    corpus, and the count would be wrong, so that raises rather than renders.
    """
    heading = "<p><strong>Rare variants</strong></p>"
    if not rows:
        return (
            f'<div class="card-variants">{heading}'
            "<p>No rare-variant study in this atlas reports a result for this gene.</p></div>"
        )

    def by_an_author(study: str) -> bool:
        publication = publications.get(study)
        return publication is not None and publication.own_lab

    def own_lab(study: str) -> str:
        return _OWN_LAB_MARK if by_an_author(study) else ""

    groups: dict[tuple[str, str], list[BurdenRow]] = {}
    for row in rows:
        if survives_own_correction(row):
            key = (row.study, _analysis_label(row.comparator, row.variant_class))
            groups.setdefault(key, []).append(row)

    studies = {row.study for row in rows}
    published_correction = {row.study for row in rows if row.pvalue_adjusted is not None}
    uncorrected: list[str] = []
    for study in sorted(studies - published_correction):
        publication = publications.get(study)
        tests = publication.tests_reported if publication is not None else None
        # "for this gene": the decision above is made from this gene's rows, so
        # the sentence claims no more than that, even for a study that corrects
        # some genes and not others.
        advice = (
            f"judge its results against its scan of {tests:,} tests (see its block)"
            if tests
            else "see its block"
        )
        name = html.escape(_study_label(study, publications))
        uncorrected.append(
            f'<p class="card-uncorrected"><a href="#{html.escape(_study_anchor(study))}">'
            f"{name}</a> publishes no corrected p for this gene; {advice}"
            f"{own_lab(study)}</p>"
        )

    if groups:
        items: list[str] = []
        for study, analysis in sorted(groups):
            surviving = sorted(groups[(study, analysis)], key=_sort_key)
            first, others = surviving[0], surviving[1:]
            publication = publications.get(study)
            named = []
            for row in others:
                genome_wide = _genome_wide(row, publication)
                named.append(
                    f"{_row_label(row)} (genome-wide, {genome_wide})"
                    if genome_wide
                    else _row_label(row)
                )
            also = f" · also surviving: {'; '.join(named)}" if named else ""
            genome_wide = _genome_wide(first, publication)
            name = _study_label(study, publications)
            line = (
                f"{name} · {analysis} — {_row_label(first)}: "
                f"{_count(first.n_case_carriers, first.n_cases, first.count_unit)} vs "
                f"{_controls(first)}{_correction_scope(first, publication)}"
                f"{f' · genome-wide ({genome_wide})' if genome_wide else ''}"
                f"{also}{own_lab(study)}"
            )
            # `caveats` where a surviving row's `method_note` qualifies these
            # counts -- it is footnoted in the study's block, outside every
            # fold -- and `details` otherwise. Whitespace is no note, as
            # `_footnotes` reads it.
            noted = any((row.method_note or "").strip() for row in surviving)
            label = "caveats" if noted else "details"
            purpose = html.escape(f"{label.capitalize()} for {name}")
            items.append(
                f"<li>{html.escape(line)} · "
                f'<a href="#{html.escape(_study_anchor(study))}" aria-label="{purpose}">'
                f"{label}</a></li>"
            )
        # A study's curated summary caveat, once per study with a line, in the
        # lines' own study order: it qualifies the counts those lines show, so
        # a study with no line -- a no-correction line shows no counts -- has
        # nothing here for it to qualify. After it, what the study's own table
        # says of the CHD cases among the carriers those lines show
        # (`_chd_only`, owner decision R2): whether the caveat reaches those
        # carriers, never whether it reaches the denominator. Only after it: a
        # study with no caveat gets neither, because that sentence speaks of
        # aneurysm probands only the caveat introduces (review, 2026-10-27).
        notes: list[str] = []
        for study in dict.fromkeys(study for study, _ in sorted(groups)):
            publication = publications.get(study)
            if publication is None or publication.summary_caveat is None:
                continue
            shown = [
                sorted(members, key=_sort_key)[0]
                for (owner, _), members in sorted(groups.items())
                if owner == study
            ]
            text = " ".join([publication.summary_caveat, *_chd_only(shown)])
            notes.append(
                f'<p class="card-note"><a href="#{html.escape(_study_anchor(study))}">'
                f"{html.escape(_study_label(study, publications))}</a>: "
                f"{html.escape(text)}</p>"
            )
        body = (
            "<p><strong>Rare variants</strong> surviving their own study's correction</p>"
            f"<ul>{''.join(items)}</ul>{''.join(notes)}"
        )
    else:
        unplaced = sorted(studies - {study for family in families for study in family})
        if unplaced:
            raise ValueError(
                f"no cohort family holds {', '.join(unplaced)}; the families were derived over "
                f"a different corpus, and counting them would misstate how many independent "
                f"datasets tested this gene"
            )
        tested = sum(1 for family in families if family & studies)
        noun = "dataset" if tested == 1 else "datasets"
        # "None survives" only where some study published a correction to
        # survive. With none published (24 genes, all tested by PMID:42230622
        # alone) it would be vacuously true and read as a correction applied
        # and failed; the no-correction line below says what is true instead.
        verdict = (
            "; none survives a correction its study published." if published_correction else "."
        )
        # The declaration this line names no study to carry (owner decision
        # 2026-10-26, R7). Counted in studies, never in datasets: a cohort
        # family can hold a study by an atlas author beside one that is not.
        own = sum(1 for study in studies if by_an_author(study))
        declared = ""
        if len(studies) == 1 and own:
            declared = " (a study by an atlas author)"
        elif own == len(studies) > 1:
            declared = " (each study by an atlas author)"
        elif own:
            declared = f" ({own} of the {len(studies)} studies by an atlas author)"
        # What the count can see (owner answer 2026-10-27): the front page's
        # qualifier, from one constant, straight after the count it qualifies
        # and before the declaration, which is about studies and which the
        # registry does not record.
        body = (
            f"{heading}<p>Tested in {tested} independent {noun} {REGISTRY_QUALIFIER}"
            f"{declared}{verdict}</p>"
        )
    return (
        f'<div class="card-variants">{body}{"".join(uncorrected)}'
        f"<p>{EVIDENCE_POWER_CAVEAT}</p></div>"
    )


def _warrant_records(fact: GeneFacts, gene_validity: GeneValidity) -> list[ValidityRecord]:
    """The records whose classification the card's warrant row states.

    For a graded gene, every ClinGen record at the headline grade -- the
    population `_graded_diseases` names in that row. For a gene admitted on
    agreement, the admitting submitters' qualifying records
    (`validity.agreeing_records`). Never every record: a GenCC submitter's
    `Autosomal recessive` row beside a panel's `AD` Definitive is not the mode
    of the classification the chip shows.
    """
    headline = fact.headline_confidence
    if headline is None:
        return list(agreeing_records(gene_validity))
    return [
        record
        for record in gene_validity.records
        if record.source is ValiditySource.CLINGEN and record.classification is headline
    ]


def _card_rows(fact: GeneFacts, gene_validity: GeneValidity) -> list[tuple[str, str]]:
    """The card's definition list: the warrant, who else lists the gene, inheritance. D55.

    **The warrant row is unchanged** -- `<grade> for <diseases>`, or `admitted
    by <submitters>` for a gene no panel graded -- and stays first; `_rail`'s
    docstring records why.

    **"Also listed by N other curators (GenCC)"** counts
    `validity.asserting_institutions` -- the `asserted_by.institutions` the
    bundle publishes -- minus ClinGen and minus the submitters that admitted
    the gene, so a ClinGen panel is never counted as "another curator" of its
    own grade, under either name it appears as. Omitted at zero (TBX1: both of
    its asserting submitters are the two that admitted it), because "0 other
    curators" reads as a measured absence of support. A count, never a score:
    eight laboratories listing a gene is not eight times the evidence
    (`admission_provenance`).

    **Inheritance is every distinct mode of `_warrant_records`**, sorted by the
    words a reader sees and joined " / ". Every mode, not the admitting
    record's alone: ROBO1, PRKD1 and LZTR1 each carry two ClinGen records at
    their headline grade for one disease, one `AD` and one `AR` (measured
    2026-10-26), so one mode would state half of what the panel graded. For a
    gene admitted on agreement, GATA5 and GDF1 carry three modes across their
    admitting records, `Unknown` among them -- the submitter's own word,
    rendered as it is in the validity table. Omitted when no record states a
    mode, rather than rendered blank.
    """
    headline = fact.headline_confidence
    # The submitters that admitted the gene: `admitted_by.submitters`, which is
    # empty where a ClinGen record admitted it. Not `agreeing_submitters` on a
    # graded gene -- NKX2-5's GenCC submitters agree at or above the floor too,
    # and subtracting them would leave "1 other curator" where seven list it.
    submitters = agreeing_submitters(gene_validity) if headline is None else ()
    rows = [
        (
            f"{headline.value} for",
            "; ".join(_graded_diseases(gene_validity, headline)) or _EM_DASH,
        )
        if headline is not None
        else ("admitted by", "; ".join(submitters) or _EM_DASH)
    ]
    others = set(asserting_institutions(gene_validity)) - {"ClinGen"} - set(submitters)
    if others:
        noun = "curator" if len(others) == 1 else "curators"
        rows.append(("also listed by", f"{len(others)} other {noun} (GenCC)"))
    modes = sorted(
        {_inheritance(record.moi) for record in _warrant_records(fact, gene_validity) if record.moi}
    )
    if modes:
        rows.append(("inheritance", " / ".join(modes)))
    return rows


def _card_gloss(fact: GeneFacts) -> str:
    """One line saying what the headline grade means, in the grade key's words. D55.

    `render.GRADE_GLOSS` verbatim -- trusted constant markup, and the same
    text the grade key renders -- so the card and the key cannot gloss one rung
    two ways. `limited` is the plurality (45 of 119 genes, measured 2026-10-26)
    and the rung most easily read as a weak yes; its gloss says it is not.

    Nothing for a gene no panel graded: there is no rung to gloss, and the
    card's curation line says why it is published instead. Nothing for a grade
    the key does not gloss either (`disputed`, `refuted`): no published gene is
    headlined by one, and a guard on a bypassed gate renders no claim rather
    than raising mid-build.
    """
    headline = fact.headline_confidence
    if headline is None:
        return ""
    gloss = GRADE_GLOSS.get(headline.value)
    if gloss is None:
        return ""
    name = html.escape(headline.value.replace("_", " ").capitalize())
    return f'<p class="card-gloss"><em>{name}</em>: {gloss}.</p>'


def _card_dispute(fact: GeneFacts, gene_validity: GeneValidity) -> str:
    """Who disputes or refutes the gene, and for what -- under the gloss. D55.

    Rendered exactly where `has_conflicting_evidence` is true. The flag means
    (`vocab.has_conflicting_evidence` over every in-scope record
    `derive.gene_facts` mirrors) that some record is `Disputed` or `Refuted`
    and some record supports the gene; it does not say which disease, and on
    all 6 flagged genes today (BRAF, KRAS, MAP2K1, PTPN11, RAF1, SOS1,
    measured 2026-10-26) the dispute is for a different disease from the
    headline's -- beside a gloss reading "with no convincing evidence to the
    contrary". A line saying only that "an authority" disputes the gene (its
    wording for one commit) told a reader neither who nor for what.

    So the line names them from the contested records themselves (`CONTESTED`):
    grouped by authority -- the expert panel for a ClinGen record, the
    submitter for a GenCC one -- and by the authority's own term, verbatim, with
    the diseases it applies it to: "Disputed for Costello syndrome and
    cardiofaciocutaneous syndrome by the RASopathy Gene Curation Expert Panel;
    see the validity table." Panels first, then submitters, each by name;
    diseases sorted. It names exactly the contested records the validity table
    lists: ClinGen's own GenCC copy of a panel record folds there
    (`_repeats_the_panel`) and is not named a second time, under ClinGen's
    name, here.

    A flag with no contested record to name means the facts and the records
    were derived apart, and the page is refused rather than published with a
    line naming nobody. The gloss's place when the gene has no gloss: a gene
    admitted on agreement can carry the flag through a submitter's dispute.
    """
    if not fact.has_conflicting_evidence:
        return ""
    panel = [record for record in gene_validity.records if record.source is ValiditySource.CLINGEN]
    contested = [
        record
        for record in gene_validity.records
        if record.classification in CONTESTED and not _repeats_the_panel(record, panel)
    ]
    if not contested:
        raise ValueError(
            f"{fact.gene} is flagged has_conflicting_evidence with no contested record to name; "
            "its facts and its validity records were derived apart"
        )
    # Diseases kept in `records` order (already sorted, `validity._sort_key`) and
    # de-duplicated through a dict, not a `set`: the sort below is then for the
    # reader, and dropping it changes the text deterministically rather than on
    # one hash seed in two.
    groups: dict[tuple[int, str, str], dict[str, None]] = {}
    for record in contested:
        if record.source is ValiditySource.CLINGEN:
            key = (
                0,
                f"the {record.gcep}" if record.gcep else "ClinGen",
                record.classification_term,
            )
        else:
            key = (
                1,
                record.submitter or "a Gene Curation Coalition submitter",
                record.classification_term,
            )
        groups.setdefault(key, {})[record.disease_label] = None
    clauses = [
        f"{term} for {_and_join(sorted(diseases))} by {authority}"
        for (_, authority, term), diseases in sorted(groups.items())
    ]
    sentence = "; ".join(clauses) + "; see the validity table."
    return f'<p class="card-dispute">{html.escape(sentence)}</p>'


def _card_no_association(fact: GeneFacts, gene_validity: GeneValidity) -> str:
    """Who reports no known disease relationship, for what, and when (owner decision R6).

    Rendered exactly where `has_no_association_report` is true: some authority
    recorded `No Known Disease Relationship` in scope while another supports
    the gene (`vocab.reports_no_association`). The bundle published the flag
    and `no_association_reported_by` since schema 2.10, and the validity table
    listed the record, but the card -- the summary a reader takes in first --
    said nothing. Measured on a build of 3b766ba: one published gene, GDF1,
    admitted on GenCC agreement, so its card's chip reads `not classified`
    and its `admitted by` row names four submitters, G2P (`Definitive`) among
    them; the line beneath reads "Illumina reports no known disease
    relationship for congenital heart defects, multiple types, 6 (2021); see
    the validity table."

    Built from the records, never from a literal or from the bundle's name
    list: grouped by authority -- the expert panel for a ClinGen record, the
    submitter for a GenCC one, as `_card_dispute` names them -- with each
    authority's diseases and the year of its own record (`_record_date`, the
    date the validity table shows). Panels first, then submitters, each by
    name; diseases sorted. ClinGen's own GenCC copy of a panel record folds as
    it does in the table (`_repeats_the_panel`), so it is not named twice.
    "no known disease relationship" is the one term either mirror maps to
    `NO_KNOWN_ASSOCIATION` (`vocab`), said in words rather than in the
    authority's capitals. A record with no calendar date names no year.

    A flag with no record to name means the facts and the records were
    derived apart, and the page is refused, as `_card_dispute` refuses.
    """
    if not fact.has_no_association_report:
        return ""
    panel = [record for record in gene_validity.records if record.source is ValiditySource.CLINGEN]
    reports = [
        record
        for record in gene_validity.records
        if record.classification is Classification.NO_KNOWN_ASSOCIATION
        and not _repeats_the_panel(record, panel)
    ]
    if not reports:
        raise ValueError(
            f"{fact.gene} is flagged has_no_association_report with no such record to name; "
            "its facts and its validity records were derived apart"
        )
    # De-duplicated through a dict, not a `set`, so the sorts below are for the
    # reader and dropping one changes the text on every hash seed.
    groups: dict[tuple[int, str], dict[str, None]] = {}
    for record in reports:
        if record.source is ValiditySource.CLINGEN:
            key = (0, f"the {record.gcep}" if record.gcep else "ClinGen")
        else:
            key = (1, record.submitter or "a Gene Curation Coalition submitter")
        dated = _date(_record_date(record))
        year = f" ({dated[:4]})" if _CALENDAR_DATE.match(dated) else ""
        groups.setdefault(key, {})[f"{record.disease_label}{year}"] = None
    clauses = [
        f"{authority} reports no known disease relationship for {_and_join(sorted(diseases))}"
        for (_, authority), diseases in sorted(groups.items())
    ]
    sentence = "; ".join(clauses) + "; see the validity table."
    sentence = sentence[0].upper() + sentence[1:]
    return f'<p class="card-no-association">{html.escape(sentence)}</p>'


# The two scope classes whose authority is a ClinGen expert panel, for which a
# marked syndrome's card line names the panel (`_card_scope`, R5).
_CLINGEN_SCOPES: Final = frozenset(
    {ScopeAuthority.CLINGEN_CHD_PANEL, ScopeAuthority.CLINGEN_EXPERT_PANEL}
)


def _panel_graders(
    validity: Mapping[str, GeneValidity], symbols: Mapping[str, str]
) -> dict[tuple[str, str], tuple[str, ...]]:
    """`(disease, panel)` -> the genes that ClinGen panel graded for that disease.

    Over every gene the validity mirrors curate in scope, published or not, by
    registry symbol (the HGNC id where the registry names none), sorted. For
    `_card_scope`, which must name them when the panel a syndrome is scoped on
    is not the one that graded the page's own gene.
    """
    graders: dict[tuple[str, str], set[str]] = {}
    for gene, gene_validity in validity.items():
        for record in gene_validity.records:
            if record.source is ValiditySource.CLINGEN and record.gcep is not None:
                graders.setdefault((record.disease, record.gcep), set()).add(
                    symbols.get(gene, gene)
                )
    return {key: tuple(sorted(names)) for key, names in graders.items()}


def _card_scope(
    fact: GeneFacts,
    gene_validity: GeneValidity,
    warrants: Mapping[str, Warrant],
    scope_terms: Mapping[str, ScopeEntry] | None = None,
    panel_graders: Mapping[tuple[str, str], tuple[str, ...]] | None = None,
) -> str:
    """Why a syndrome's grade is on a CHD atlas, beside the chip. D55.

    "Graded for Noonan syndrome as a whole syndrome; in scope because Orphanet
    records a structural heart malformation in at least 30% of cases." A
    syndrome can be in scope on Orphanet's phenotype frequencies
    (`orphanet_scope`), and a grade for it is a grade for the syndrome. The
    validity section's notes say so at length; this says it beside the chip,
    with the reason -- which "Graded for Noonan syndrome, not for heart
    disease", its wording until 2026-10-26, left out.

    **Only when every disease the headline grade is for is Orphanet-scoped**,
    naming all of them, joined with "and" rather than the warrant row's "; ".
    A gene graded at its headline for a syndrome and for a disease an
    authority curated as CHD would read as in scope through Orphanet while it
    is not only that -- 0 of the 98 graded genes today (measured 2026-10-26),
    so the condition is latent and fixture-pinned. And all of them, not the
    admitting record's alone: PTPN11 is Definitive for Noonan syndrome *and*
    Noonan syndrome with multiple lentigines. Nothing for a gene no panel
    graded: nothing was graded.

    **A subtype's warrant is its syndrome's** (`Warrant.via == "subtype"`):
    Orphanet records the malformation for Adams-Oliver syndrome, and EOGT's
    Adams-Oliver syndrome 4 is in scope as its subtype -- the one such headline
    among the 21 Orphanet lines measured 2026-10-26. The clause then names the
    Orphanet syndromes the frequency was recorded for, so it is never read as
    a frequency recorded for the subtype itself.

    **A syndrome a ClinGen expert panel scoped gets the line too** (owner
    decision 2026-10-26, R5): "Graded for Holt-Oram syndrome as a whole
    syndrome; in scope on its curation by ClinGen's Syndromic Disorders Gene
    Curation Expert Panel." Ten genes are graded at their headline for such a
    syndrome and their cards said nothing. Which terms are syndromes is the
    curated `ScopeEntry.syndrome` mark, because nothing else says it: ZIC3 and
    NODAL are graded for "congenital heart disease with heterotaxy syndrome",
    which is CHD itself. Rendered only when **every** headline disease is a
    marked term an authority of a ClinGen class scoped, for the reason the
    Orphanet line needs all of them.

    **The reason clause is not the owner's "because <panel> curates it"**,
    which reads as sufficient and is not: measured 2026-10-26 on the ClinGen
    mirror, the Syndromic Disorders panel curates 176 disease terms and 3 are
    in scope, Hearing Loss 63 and 1. What `admitted_by` records, and SCP005
    checks against the mirror, is that the term is listed on that panel's
    curation of it -- so the line says it is in scope *on* that curation. The
    authority is the term's `attributed_to`, verbatim, after "ClinGen's".

    **That panel need not be the one that graded this gene** (review,
    2026-10-27). `attributed_to` names the panel whose curation scopes the
    *term*. EVC2's grade for Ellis-van Creveld syndrome is the Skeletal
    Disorders panel's, and the term is scoped on the Kidney Cystic and
    Ciliopathy Disorders panel's curation of it for EVC; "in scope on its
    curation by" that panel read as the panel having graded EVC2. Where no
    ClinGen record of this gene for the term is that panel's, the clause names
    the genes that panel did grade for it (`panel_graders`, from the validity
    mirror the build already holds): "in scope on ClinGen's <panel>'s curation
    of it for EVC." Measured on a build of 5f7161b: 1 of the 10 lines (EVC2);
    the other nine panels graded the page's own gene. None found is a scope
    term SCP005 refuses, so it raises rather than inventing a gene.

    **Each disease is paired with its own panel.** With two diseases the
    clause sorted the panels apart from the diseases and joined the two lists,
    so Holt-Oram syndrome and TARP syndrome read as curated by the Congenital
    Heart Disease and the Syndromic Disorders panels -- the reverse of the
    terms' attributions. Several diseases now get one "ClinGen's <panel>'s
    curation of <disease>" each, in the order the line names the diseases.
    Latent: every one of the ten lines names one disease.
    """
    headline = fact.headline_confidence
    if headline is None:
        return ""
    warranting = _warrant_records(fact, gene_validity)
    graded = {record.disease for record in warranting}
    if not graded:
        return ""
    diseases = _graded_diseases(gene_validity, headline)
    noun = "a whole syndrome" if len(diseases) == 1 else "whole syndromes"
    if graded <= warrants.keys():
        found = [warrants[disease] for disease in sorted(graded)]
        recorded_for = (
            ""
            if all(warrant.via == "direct" for warrant in found)
            else f" of {_and_join(sorted(dict.fromkeys(warrant.name for warrant in found)))}"
        )
        reason = (
            "in scope because Orphanet records a structural heart malformation in at least 30% "
            f"of cases{recorded_for}"
        )
    else:
        terms = [(scope_terms or {}).get(disease) for disease in sorted(graded)]
        marked = [
            term
            for term in terms
            if term is not None and term.syndrome and term.admitted_by in _CLINGEN_SCOPES
        ]
        if len(marked) != len(terms):
            return ""
        labels = {record.disease: record.disease_label for record in warranting}
        clauses: list[tuple[str, str, tuple[str, ...]]] = []
        for term in sorted(marked, key=lambda term: (labels[term.id], term.id)):
            own = {
                record.gcep
                for record in gene_validity.records
                if record.source is ValiditySource.CLINGEN and record.disease == term.id
            }
            others: tuple[str, ...] = ()
            if term.attributed_to not in own:
                others = (panel_graders or {}).get((term.id, term.attributed_to), ())
                if not others:
                    raise ValueError(
                        f"scope term {term.id} is attributed to {term.attributed_to!r}, which "
                        f"graded no gene for it in the validity mirror; SCP005 refuses this"
                    )
            clauses.append((f"ClinGen's {term.attributed_to}", labels[term.id], others))
        if len(clauses) == 1:
            [(panel, _, others)] = clauses
            reason = (
                f"in scope on {panel}'s curation of it for {_and_join(others)}"
                if others
                else f"in scope on its curation by {panel}"
            )
        else:
            reason = "in scope on " + _and_join(
                [
                    f"{panel}'s curation of {label}"
                    + (f" for {_and_join(others)}" if others else "")
                    for panel, label, others in clauses
                ]
            )
    sentence = f"Graded for {_and_join(diseases)} as {noun}; {reason}."
    return f'<p class="card-scope">{html.escape(sentence)}</p>'


def _rail(
    gene: str,
    symbol: str,
    fact: GeneFacts,
    gene_validity: GeneValidity,
    rows: Sequence[BurdenRow],
    *,
    families: tuple[frozenset[str], ...] = (),
    publications: Mapping[str, Publication] | None = None,
    warrants: Mapping[str, Warrant] | None = None,
    identity: GeneIdentity | None = None,
    scope_terms: Mapping[str, ScopeEntry] | None = None,
    panel_graders: Mapping[tuple[str, str], tuple[str, ...]] | None = None,
    literature: GeneLiterature = EMPTY_GENE_LITERATURE,
) -> str:
    """The clinician card: the summary column, for every gene whether curated here or not.

    **D55: the rail is the card.** Top to bottom: the symbol and the gene's
    identity; the HGNC id and the chips; a definition list (`_card_rows`: the
    warrant, who else lists the gene, its inheritance); the headline grade's
    gloss (`_card_gloss`); a disputed-gene line (`_card_dispute`); a
    no-association line (`_card_no_association`, R6); the Orphanet line
    (`_card_scope`); the rare-variant results that survive their own
    study's correction (`_card_rare_variants`, D56), each linking to its
    study's block; the atlas-curation sentence (`_card_curation`); and a last
    line linking the top of `methods.html` as "How to read this page" (P1,
    2026-10-26) beside the gene's JSON. The nav's link to the same file reads
    "How to read this atlas"; from the card the question is about this page.
    Every row from data the page already has -- nothing computed.

    It replaced a definition list of raw tokens and zero counts -- `validity:
    expert_curated`, `atlas curation: not_yet_curated`, `assertions: 0`,
    `functional records: 0`, `curated publications: 0` -- plus a burden-row
    count, which say nothing a clinician can act on and read a measured zero
    where they mean "not curated here". Every one of them is still in the
    gene's JSON. The main column's "not yet curated" notice moved here too, as
    the curation sentence, so the page states it once.

    **The chip class is looked up in `GRADE_CHIP_KIND`, keyed on the grade's own
    value.** Until 2026-08-06 it was `"definitive" if headline is DEFINITIVE else
    ""`, which was right while the gate admitted nothing else and became wrong
    the moment it widened: 53 of the 92 genes published carry `strong`,
    `moderate` or `limited`, and a bare `.chip` renders all three in the same
    grey as `not classified`. A reader cannot tell a panel's weakest supportive
    rung from no panel having looked, which is the distinction the widened gate
    turns on. The definitions are the other half -- the colours are ordinal and
    mean nothing without ClinGen's own words beside them. Since 2026-10-26 (D50)
    the card's gloss line gives the headline rung's, and `render.grade_key` on
    `methods.html#scope`, linked from the page, gives every rung's.

    **The first row's label carries the grade, so the label and the disease
    cannot come apart.** It read `definitive for` unconditionally until the same
    date, and `_graded_diseases` was returning the *definitive* records, so on a
    `limited` gene the row was headed `definitive for` beside a `limited` chip.
    Now the label is the headline's own value -- `limited for: atrial septal
    defect` -- and for a gene no panel graded it becomes `admitted by`, naming
    the submitters whose agreement admitted it. Both come from `fact`, so the
    chip and the row below it are one claim.

    It stays first in the list for the reason the browse column sits beside
    `confidence`: the chip one line above says a bare grade, and on this site's
    title and tagline that reads as "definitive for congenital heart disease".
    For KMT2D it means Kabuki syndrome 1. The row is what stops the chip standing
    alone on the first screen.

    **The gene's identity sits under its symbol (D52)**: approved name, locus,
    MANE Select transcript, and a link row with gnomAD and ClinVar first.
    `external_links` raises for a gene missing the ids those two are keyed on,
    so a page is refused rather than published without them. `identity=None`
    renders none of it; only a caller holding no registry at all passes that,
    and `build_site` always passes one.
    """
    headline = fact.headline_confidence
    chips = [
        chip(
            headline.value if headline is not None else "not classified",
            kind=GRADE_CHIP_KIND[headline.value] if headline is not None else "ungraded",
        )
    ]
    if fact.has_conflicting_evidence:
        chips.append(chip("conflicting evidence", kind="warn"))
    if fact.has_source_discordance:
        chips.append(chip("sources disagree", kind="warn", title=_SOURCE_DISCORDANCE_TITLE))

    identity_block = ""
    if identity is not None:
        links = external_links(gene, identity)
        details = [identity.locus] if identity.locus else []
        if identity.mane_select:
            details.append(f"MANE Select {identity.mane_select}")
        identity_block = (
            (f"<p>{html.escape(identity.name)}</p>" if identity.name else "")
            + (f"<p>{html.escape(' · '.join(details))}</p>" if details else "")
            + '<p class="xrefs">'
            + " · ".join(
                f'<a href="{html.escape(link.href)}">{html.escape(link.label)}</a>'
                for link in links
            )
            + "</p>"
        )
    bundle = gene_bundle_path(HgncId(gene))
    # `html.escape` here rather than a `render.py` primitive: these three, and
    # the identity block above, are bare interpolations into markup this
    # function assembles. The card's rows go through `definition_list`, which
    # escapes its arguments, and each card helper escapes what it interpolates;
    # `GRADE_GLOSS` and `EVIDENCE_POWER_CAVEAT` are trusted constant markup.
    return (
        f'<aside class="rail"><h1>{html.escape(symbol)}</h1>{identity_block}'
        f"<p><code>{html.escape(gene)}</code></p><p>{''.join(chips)}</p>"
        f"{definition_list(_card_rows(fact, gene_validity))}"
        f"{_card_gloss(fact)}"
        f"{_card_dispute(fact, gene_validity)}"
        f"{_card_no_association(fact, gene_validity)}"
        f"{_card_scope(fact, gene_validity, warrants or {}, scope_terms, panel_graders)}"
        f"{_card_rare_variants(rows, families, publications or {})}"
        f"{_card_curation(fact, literature)}"
        '<p><a href="../methods.html">How to read this page</a> &middot; '
        f'<a href="../{html.escape(bundle)}">this gene as JSON</a></p></aside>'
    )


# GenCC's `report_url` for every PanelApp Australia submission: the site's bare
# homepage, on all 5,478 of them in the mirror (measured 2026-10-26), so the
# validity table's "open" link landed on a page that does not name the gene.
_PANELAPP_AUSTRALIA: Final = "PanelApp Australia"
_PANELAPP_AUSTRALIA_HOME: Final = "https://panelapp-aus.org"


def _report_url(record: ValidityRecord, symbol: str | None) -> str | None:
    """Where a validity row's report link goes: GenCC's URL, with one exception (R8).

    Owner decision 2026-10-26. A PanelApp Australia record whose `report_url`
    is exactly the bare homepage links that site's page for the gene instead,
    `https://panelapp-aus.org/panels/entities/<HGNC symbol>`: verified
    2026-10-26 to answer 200 for TBX5, whose PanelApp Australia API record
    lists it on 11 panels. **The one report URL this atlas builds rather than
    republishes** (the identity links are built too, from ids, by
    `identity.external_links`), because the one GenCC supplies names no gene and so cannot
    lead to the record the row reports. Every other submitter's URL, and a
    PanelApp Australia URL that is anything but the homepage, is GenCC's as
    given; the mirror and the bundle keep GenCC's value. The symbol is the gene
    registry's, percent-encoded; a gene the registry does not name keeps the
    homepage rather than a URL built on its HGNC id.
    """
    if (
        record.source is ValiditySource.GENCC
        and record.submitter == _PANELAPP_AUSTRALIA
        and record.report_url == _PANELAPP_AUSTRALIA_HOME
        and symbol
    ):
        return f"{_PANELAPP_AUSTRALIA_HOME}/panels/entities/{quote(symbol, safe='')}"
    return record.report_url


def _validity_row(record: ValidityRecord, symbol: str | None = None) -> Row:
    """One mirrored classification as a row, in the section's one vocabulary.

    `classification_term` is what the authority published, not the rung it maps
    onto: GenCC's `Supportive` maps to no rung at all, and a table showing the
    mapped value would render that submitter's row as blank or, worse, as a
    grade nobody gave. `symbol` is the gene's registry symbol, for
    `_report_url`.
    """
    curator = record.gcep if record.source is ValiditySource.CLINGEN else record.submitter
    term = record.classification_term
    kind = record.classification.value if record.classification else "ungraded"
    chip_html = chip(term, kind=kind)
    gold = record.source is ValiditySource.CLINGEN
    # The one row-level distinction on the table: a chartered expert panel wrote
    # this row, not a submitting laboratory. Marked on the row, not with a
    # separate table, so a reader compares every classification of a disease in
    # one place and still sees which one is the gold standard.
    curated_by: Cell = (
        Markup(f"{html.escape(curator or _EM_DASH)}{_GOLD_STANDARD_BADGE}")
        if gold
        else curator or _EM_DASH
    )
    return Row(
        attributes=(("authority", "clingen" if gold else "gencc"),),
        cells=(
            curated_by,
            f"{record.disease_label} ({record.disease})",
            _inheritance(record.moi),
            Markup(chip_html),
            _date(_record_date(record)),
            _report_link(_report_url(record, symbol), record.sop),
        ),
    )


def _validity_section(
    gene_validity: GeneValidity, notes: str = "", symbol: str | None = None
) -> str:
    """The mirrored classifications: the expert panel's first, then everyone else's. D51.

    `notes` (the Orphanet scope notes) is placed after the intro and before the
    first group, because a scope warrant concerns the disease and the disease
    may be graded in either group, or both.

    Two groups under one heading rather than a `source` column, because the
    column made a reader decode `clingen` and `gencc` to learn the one thing
    that matters about a row -- whether a chartered panel or a submitting
    laboratory wrote it. Both groups share one header and one vocabulary
    (`_inheritance`, `_date`), so they read across.

    An empty group is a sentence, not a header over no rows: 21 published genes
    have no panel record at all (measured 2026-10-25), and an empty "Expert
    panel" table reads as a rendering fault rather than as that fact.

    `_repeats_the_panel` decides which of ClinGen's own GenCC copies are
    listed. Measured on a real build 2026-10-25: 591 records, 121 copies folded,
    470 rows rendered.

    **The intro states the table's own truth condition (2026-10-26, D50).** The
    grade key's last sentence -- a gene page lists only records for diseases in
    this atlas's remit -- left the page with the key for `methods.html#scope`,
    and it is a fact about *this* table, so a short form of it now sits where
    the table is, linking `methods.html#validity` for how the records are
    listed. Both groups are filtered, not only ClinGen's: `gene_validity` drops
    an out-of-scope row from either mirror before a gene is keyed. It says the
    table **may not** be every classification the gene has, not that it is
    not: measured 2026-10-26 against the committed mirrors, 43 of the 119
    published genes have no ClinGen or GenCC record for any disease outside
    scope, so "is not" would be false on 43 pages.
    """
    panel = [r for r in gene_validity.records if r.source is ValiditySource.CLINGEN]
    submitted = [r for r in gene_validity.records if r.source is ValiditySource.GENCC]
    listed = [r for r in submitted if not _repeats_the_panel(r, panel)]
    panel_rows = [_validity_row(record, symbol) for record in panel]
    listed_rows = [_validity_row(record, symbol) for record in listed]
    folded = len(listed) < len(submitted)
    clingen_listed = any(record.submitter == "ClinGen" for record in listed)
    rows = panel_rows + listed_rows
    return (
        "<h2>Mirrored gene&ndash;disease validity</h2>"
        "<p>Every classification below is an upstream panel's or submitter's, "
        "republished with its provenance intact. Only records for diseases in this "
        "atlas&#x27;s remit are listed, so this may not be every classification this gene "
        'has (<a href="../methods.html#validity">how these records are listed</a>).</p>'
        + notes
        + ("" if panel_rows else _NO_PANEL_RECORD)
        + ("" if listed_rows else _NO_OTHER_CURATOR)
        + (data_table(_VALIDITY_HEADERS, rows, table_id="validity-table") if rows else "")
        + (_CLINGEN_COPY_NOTICE if folded or clingen_listed else "")
    )


def _orphanet_notes(
    gene_validity: GeneValidity, warrants: Mapping[str, Warrant], labels: Mapping[str, str]
) -> str:
    """Why each Orphanet-scoped disease in this gene's records counts as CHD here.

    Since 2026-10-25 (issue #8 part 1) a syndrome can be in scope because Orphanet
    records structural heart malformations in at least 30% of its cases --
    `orphanet_scope.FREQUENCY_FLOOR`, `Frequent (79-30%)` -- rather than because
    an authority curated it as CHD. Without this note PTPN11 would sit in a CHD
    atlas with a Noonan grade and nothing saying why Noonan counts, and a reader
    would take the grade for a grade for heart disease.

    One row per distinct Orphanet-scoped disease among the gene's records,
    sorted by disease id, in a table above both groups of classifications: it is
    particular to this gene. The table's third column names the malformations
    Orphanet records for the disease and the band it records them at, both
    verbatim; the second says how the disease is reached -- its own cases, or a
    validated subtype of a parent syndrome, with the OMIM-to-Orphanet mapping
    that certificates it.

    **The lead-in does not name any one disease for the whole table's sake.**
    It read "The classifications below are for the syndrome, not for heart
    disease" until 2026-10-26, above a table that is not only that syndrome's:
    NOTCH1's lists an admitting record for a term a ClinGen expert panel scoped
    (NOTCH1-related AOS spectrum disorder) and an Orphanet record for familial
    bicuspid aortic valve beneath its two Adams-Oliver notes, so a single
    undirected sentence was false on that page. The whole-table lead-in says
    *each* disease's classifications are for *that* disease, which is true no
    matter how many diseases the section below names; the caveat names no
    parent, so "the syndrome" cannot read as one row's parent.

    **The lead-in ends on "that disease", not "the syndrome"** (review,
    2026-10-26): a subtype's row names its parent syndrome in the "why" cell --
    "subtype of Noonan syndrome ..." -- so "the syndrome as a whole" would read
    as the parent's.
    """
    diseases = sorted({record.disease for record in gene_validity.records} & set(warrants))
    if not diseases:
        return ""
    rows: list[Row] = []
    for disease in diseases:
        found = warrants[disease]
        label = html.escape(labels.get(disease, disease))
        malformation_cells = "; ".join(
            f"{html.escape(a.label)} at <em>{html.escape(a.frequency or '')}</em>"
            for a in found.phenotypes
        ) + (f" (Orphanet, {html.escape(found.orpha)})" if found.via == "direct" else "")
        if found.via == "subtype":
            why = (
                "Subtype of "
                + html.escape(found.name)
                + "; Orphanet records a structural heart malformation in at least 30% of its "
                "cases (Orphanet maps "
                + html.escape(found.omim or "")
                + " under "
                + html.escape(found.orpha)
                + ")"
            )
        else:
            why = (
                "Orphanet records a structural heart malformation in at least 30% of this "
                "disease's cases"
            )
        rows.append(Row(cells=(label, Markup(why), Markup(malformation_cells))))
    return (
        "<h3>Why these diseases are in scope</h3>"
        '<p class="notice-inline">These diseases count as congenital heart disease in '
        "this atlas because Orphanet records a structural heart malformation in at least "
        "30% of their cases. Each disease&#x27;s classifications in the tables below are "
        "for that disease as a whole, not for heart disease alone.</p>"
        + data_table(_ORPHANET_SCOPE_HEADERS, rows)
    )


def _evidence_section(
    records: Sequence[LesionAssertion], publications: Mapping[str, Publication]
) -> str:
    """Each assertion and its evidence, each item naming the paper it cites.

    The paper is named in the row, "Basson CT et al. 1997" linked to PubMed with
    its title on hover (`_paper_link`), as the rare-variant and literature
    sections name theirs. A "Publications" table listing the cited papers by id,
    title and year followed this section until the owner removed it
    (2026-11-09): it repeated what the rows linked, and the Literature section
    takes its place on the page.
    """
    blocks: list[str] = []
    for assertion in records:
        summary = definition_list(
            [
                ("lesion groups", ", ".join(group.value for group in assertion.lesion_groups)),
                ("inheritance", ", ".join(assertion.inheritance)),
                ("mechanism", assertion.mechanism),
                ("syndromic", assertion.syndromic),
            ]
        )
        items = [
            Row(
                cells=(
                    item.evidence_class.value,
                    item.strength.value,
                    item.summary,
                    _paper_link(item.publication, publications),
                )
            )
            for item in assertion.evidence
        ]
        blocks.append(
            f"<h3>{html.escape(assertion.id)}</h3>{summary}{data_table(_EVIDENCE_HEADERS, items)}"
        )
    # The section's link to `methods.html#curated-evidence` (owner decision
    # 2026-10-26, E3), in the pattern of the other sections' one-line intros.
    # "above" is true in this position: the validity table renders first.
    return (
        "<h2>Curated evidence</h2>"
        "<p>Each lesion assertion below is this atlas's own curation, and changes no "
        'classification above. <a href="../methods.html#curated-evidence">How to read '
        "curated evidence</a>.</p>" + "".join(blocks)
    )


def _paper_link(pmid: str, publications: Mapping[str, Publication]) -> Markup:
    """A cited paper as "Basson CT et al. 1997", linked to PubMed, its title on hover.

    A paper missing from `publications` keeps its link and is labelled by its
    PMID (`_study_label`'s fallback) with no hover title, rather than dropping
    the citation: the row is the curator's, and a citation that vanishes because
    a lookup missed is the silent evidence loss this project guards against.
    `referential.py` refuses a corpus citing an unknown publication, so this is
    a guard on a bypassed gate.
    """
    publication = publications.get(pmid)
    hover = f' title="{html.escape(publication.title)}"' if publication is not None else ""
    return Markup(
        f'<a href="{html.escape(_pubmed(pmid))}"{hover}>'
        f"{html.escape(_study_label(pmid, publications))}</a>"
    )


def _study_label(pmid: str, publications: Mapping[str, Publication]) -> str:
    """ "Audain et al. 2026", or the bare PMID when the record is missing.

    The fallback is a guard on a bypassed gate -- BUR010 reports a burden row
    citing an unregistered study and `build_site` refuses on it -- and it renders
    the PMID rather than an em dash for the reason `build_genes` falls back to an
    HGNC id: an identifier is something a reader can still look up.
    """
    publication = publications.get(pmid)
    if publication is None:
        return pmid
    # `Publication.authors` constrains the *list* to `min_length=1`, not its
    # elements, so a whitespace-only first author passes every validator and
    # renders "    et al. 2026" -- a study heading naming nobody. Same defect
    # `tables.cell_text` exists to catch on the mirror side, reaching a page for
    # the first time here: `authors[0]` was previously read only by `search.py`.
    author = publication.authors[0].strip()
    if not author:
        return pmid
    return f"{author} et al. {publication.year}"


def _study_anchor(study: str) -> str:
    """`study-34324492`: the `id` of a study's block, and the card's link to it (D55).

    The PMID's digits, so the anchor is stable across builds and readable in a
    URL. One function for both ends, so the heading and the link that points at
    it cannot be spelled two ways.
    """
    return f"study-{study.removeprefix('PMID:')}"


def _count(carriers: int | None, total: int | None, unit: str) -> str:
    """ "6 carriers / 3,876" -- a numerator, its denominator, and their unit.

    Never the numerator alone. Six carriers is a different claim in 3,876 cases
    than in 45,082, and the two columns of this table hold exactly that contrast.
    The separator is a literal `,` via `:,`, which is locale-independent, so two
    builds on two machines render the same bytes.

    **And never without the unit.** This is `_effect`'s rule applied to the
    count columns, for the same reason and after the same near miss: until
    2026-08-05 the header read "cases (carriers / n)" and every cell was assumed
    to be people. PMID:40127276's case-control arm counts *alleles* -- measured,
    not assumed, in `CountUnit` -- so that header would have claimed 21,768
    people sequenced where 11,555 were, and called an allele count a carrier
    count, on a page whose checksum verified. An unknown unit therefore renders
    its raw token rather than nothing: there is no branch here that omits it.
    """
    if carriers is None or total is None:
        return _EM_DASH
    numerator, denominator = _COUNT_UNIT_LABEL.get(unit, (unit, ""))
    # "1 carrier", not "1 carriers": published in tables and on the card until
    # 2026-10-26. Only the two plural nouns change; "de novo" has no plural,
    # and an unknown unit's raw token is never rewritten.
    if carriers == 1:
        numerator = _SINGULAR_COUNT_NOUN.get(numerator, numerator)
    return f"{carriers:,} {numerator} / {total:,}{denominator}"


def _effect(row: BurdenRow) -> str:
    """The effect size, its measure, and its interval -- never a bare number.

    **This is the guard the single `effect` column was chosen against.** One
    column holding both odds ratios and de novo enrichments is what lets this
    schema absorb a fifth study without a migration, and it is also the one place
    two incomparable quantities could silently merge: an odds ratio of 3.1 and an
    enrichment of 3.1 are different claims, and a cell reading "3.1" under a
    header reading "effect" equates them. `_MEASURE_LABEL` is therefore consulted
    on every row and there is no branch that omits it.

    `∞` rather than a number, and rather than a blank. Fisher's exact test
    returns an infinite odds ratio where no control carries, which
    `encode_json`'s `allow_nan=False` refuses to publish -- and which says
    nothing about how strong the result is: 20 of the 37 such mirrored rows
    have a lower bound below 1 (measured 2026-10-26). The lower bound survives
    and is the whole finding: TAB2's syndromic row says the true odds ratio is
    at least 28.1, which a blank cell would have thrown away.
    """
    if row.effect_measure is None:
        return _EM_DASH
    measure = _MEASURE_LABEL.get(row.effect_measure, row.effect_measure)
    if row.effect_bound == "unbounded_above":
        # `∞ (95% CI 0.298–∞)`, not `unbounded (95% CI ≥0.298)`.
        #
        # Measured 2026-08-05 on the built site: 19 cells render an unbounded
        # odds ratio and **11 of them have a lower bound below 1** -- GATA6 and
        # TBX5 at ≥0.298 with p = 0.079, among others. The old wording led with
        # the word "unbounded", which reads as "infinitely enriched", and put
        # the bound that contradicts it in a parenthesis at lower salience, in a
        # `≥` notation that is not how a Fisher interval is written. Eleven times
        # out of nineteen the parenthetical was the entire message.
        #
        # The conventional form puts the interval in the reader's usual shape,
        # so an interval spanning 1 is visible as one. `∞` is safe here and not
        # in the payload: this is HTML text, while `allow_nan=False` governs the
        # JSON, which carries `null` plus `effect_bound`.
        low = _fmt(row.ci_low) if row.ci_low is not None else None
        return f"{measure} ∞" + (f" (95% CI {low}–∞)" if low else "")
    if row.effect is None:
        return _EM_DASH
    value = f"{measure} {_fmt(row.effect)}"
    if row.ci_low is None or row.ci_high is None:
        return value
    return f"{value} (95% CI {_fmt(row.ci_low)}–{_fmt(row.ci_high)})"


def _fmt(value: float) -> str:
    """A statistic at three significant figures, without scientific notation.

    `f"{1810.0:.3g}"` is `1.81e+03`. Measured on the built site, that spelling
    reached four confidence-interval upper bounds -- `1.81e+03`, `2.35e+03`,
    `2.18e+03`, `1.36e+03` -- sitting mid-column beside plain decimals like
    `1.3`. The upper bound is the number that says "do not act on this", and in
    the least legible glyph in the table a skimming reader can take `1.81e+03`
    for `1.81`, turning an uninformative interval into a narrow-looking one.

    Values at or above 1,000 render as separated integers instead: `1,810`.
    Three significant figures are already far more precision than an interval
    that wide carries. `:,` is locale-independent, like `_count`'s.

    Not used for p-values: there an exponent is correct and expected, and
    `2.45e-07` is how a reader wants to see it.
    """
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    return f"{value:.3g}"


def _corrected(row: BurdenRow) -> str:
    """A published corrected p-value, named for the correction that produced it.

    **This column is the answer to the sharpest finding of the 2026-08-05
    review**, which measured 187 uncorrected p-values on gene pages with nothing
    to judge them against. The atlas still computes no correction -- that would
    be authoring a statistic (D12/D33) -- but where a study publishes one, it
    reaches the reader.

    It changes conclusions, not just presentation. CHD7 in PMID:34324492 has a
    raw permutation p of 0.0068 and a family-wise corrected p of 0.99: without
    this column the page shows a number that reads as significant for a gene the
    study found nothing for.

    An em dash where a study published no correction, never a blank: the method
    line above the table then names how many comparisons the study ran, which is
    what a reader has instead.
    """
    if row.pvalue_adjusted is None:
        return _EM_DASH
    label = _ADJUSTMENT_LABEL.get(row.pvalue_adjustment or "", row.pvalue_adjustment or "")
    return f"{row.pvalue_adjusted:.3g} ({label})" if label else f"{row.pvalue_adjusted:.3g}"


def _finding_detail(row: BurdenRow) -> str:
    """A card's " (effect, statistic)": `_effect` and `_statistic`, so the card
    names the measure and the correction as the matrix cell and the forest do.
    `_statistic` calls a family-wise p "corrected p", never q."""
    parts = [part for part in (_effect(row), _statistic(row)) if part != _EM_DASH]
    return f" ({html.escape(', '.join(parts))})" if parts else ""


def _study_finding_summary(study_rows: Sequence[BurdenRow]) -> str:
    """The study's finding at this gene: which consequence classes survive its
    own correction, as its rows state it, and nothing more.

    **No mechanism.** No study here states one, and inferring haploinsufficiency
    or gain of function from which class survives is a conclusion no source
    reached. "Loss-of-function variants show no enrichment" was false on GDF1,
    SMAD2 and MYH6, whose LoF rows were nominal (p 0.009 to 0.023); the card now
    says only that they do not survive, and says it only of a class the study
    tested -- PMID:34324492 tests deletions alone, and a missing class is not a
    null result. Measured 2026-11-02, before this: 31 of 93 cards named a
    mechanism, 10 called a family-wise p "q", and 35 said a PMID:42230622 row
    "does not reach" a correction that study never published.
    """
    lof_surv = [r for r in study_rows if r.consequence_class == "lof" and _survived_correction(r)]
    ms_surv = [
        r
        for r in study_rows
        if r.consequence_class == "missense_damaging" and _survived_correction(r)
    ]
    union_surv = [
        r for r in study_rows if r.consequence_class == "damaging" and _survived_correction(r)
    ]
    tested = {r.consequence_class for r in study_rows}

    if lof_surv and not ms_surv:
        return (
            '<div class="finding-card finding-lof">'
            '<div class="finding-badge">★ Loss-of-Function Significant</div>'
            '<div class="finding-summary">Enrichment of <strong>loss-of-function variants</strong> '
            f"survives this study's correction{_finding_detail(lof_surv[0])}."
            f"{' Damaging missense variants do not.' if 'missense_damaging' in tested else ''}"
            "</div></div>"
        )
    if ms_surv and not lof_surv:
        return (
            '<div class="finding-card finding-missense">'
            '<div class="finding-badge">★ Damaging Missense Significant</div>'
            '<div class="finding-summary">Enrichment of '
            "<strong>damaging missense variants</strong> "
            f"survives this study's correction{_finding_detail(ms_surv[0])}."
            f"{' Loss-of-function variants do not.' if 'lof' in tested else ''}"
            "</div></div>"
        )
    if lof_surv and ms_surv:
        return (
            '<div class="finding-card finding-both">'
            '<div class="finding-badge">★ Both LoF & Missense Significant</div>'
            '<div class="finding-summary">Both <strong>loss-of-function</strong> and '
            "<strong>damaging missense</strong> variants survive this study's correction.</div>"
            "</div>"
        )
    if union_surv:
        return (
            '<div class="finding-card finding-union">'
            '<div class="finding-badge">★ Damaging Composite Significant</div>'
            '<div class="finding-summary">Enrichment of combined '
            "<strong>damaging variants</strong> survives this study's correction"
            f"{_finding_detail(union_surv[0])}.</div>"
            "</div>"
        )

    # "Enrichment" only of a row that points to it (`is_upward`), so a
    # depleted p below 0.05 earns no card.
    nominals = [
        r
        for r in study_rows
        if r.pvalue is not None
        and r.pvalue < 0.05
        and r.consequence_class != "synonymous"
        and is_upward(r)
    ]
    if not nominals:
        return ""
    if all(r.pvalue_adjusted is None for r in nominals):
        return (
            '<div class="finding-card finding-nominal">'
            '<div class="finding-badge">Uncorrected p only</div>'
            '<div class="finding-summary">Nominal enrichment (<em>p</em> &lt; 0.05); this '
            "study publishes no corrected p for this gene.</div>"
            "</div>"
        )
    return (
        '<div class="finding-card finding-nominal">'
        '<div class="finding-badge">Nominal Enrichment (Uncorrected)</div>'
        '<div class="finding-summary">Nominal enrichment (<em>p</em> &lt; 0.05) that does not '
        "survive this study's correction.</div>"
        "</div>"
    )


# The de novo (mutation-model) structured consequence table's columns. A
# forest plot is the wrong idiom for these rows: Poisson de novo enrichment
# carries a point estimate and no symmetric confidence interval, and a null
# result -- three consequence rows all at `p 1 · q 1` -- draws as three
# identical-looking flat bars that each appear to carry a message. This table
# puts the numbers a geneticist actually reads in a row each.
_DENOVO_HEADERS: Final = (
    "consequence",
    "observed de novo",
    "expected",
    "\u00d7 enrichment",
    "corrected p",
    "reading",
)


def _denovo_enrichment(row: BurdenRow) -> str:
    """A de novo row's fold enrichment, as a bare `\u00d7` value: the column head
    names the measure, so a lone number here left the reader no cheaper way to
    read it -- `_effect`'s no-bare-number rule governs a *mixed* column and has
    no force in one whose head is the measure."""
    if row.effect is None:
        return _EM_DASH
    return f"{_fmt(row.effect)}\u00d7"


def _denovo_reading(row: BurdenRow) -> Markup:
    """The one-line reading of a de novo row, as a tinted chip."""
    if _survived_correction(row) and row.consequence_class != "synonymous":
        return Markup('<span class="denovo-chip denovo-sig">Significant</span>')
    if row.pvalue is not None and row.pvalue < 0.05:
        return Markup('<span class="denovo-chip denovo-nominal">Nominal (uncorrected)</span>')
    return Markup('<span class="denovo-chip denovo-none">No enrichment</span>')


def _denovo_table(rows: Sequence[BurdenRow]) -> str:
    """A study's de novo (mutation-model) enrichment rows, as a structured table.

    This is the de novo panel after the 2026-11 review: the sparse forest that
    used to draw these rows put a point estimate and no whisker on a
    logarithmic axis, and a fully null panel rendered three identical flat bars
    reading `p 1 · q 1`. The consequence rows are the whole reading for a
    biologist -- which consequence class is enriched, how many de novos were
    observed against the modelled expectation, by how much, and whether it
    clears the study's own correction -- so they now render as a table, the
    idiom the owner approved for them, and the forest survives only for
    case-control odds ratios and rate ratios, which have real confidence
    intervals to draw. `data-effect-measure` is kept on the figure so the two
    panels are never conflated and the measure-pair guard still sees both.

    `expected` is the mutation-model expectation the observed count was tested
    against; where a study publishes none (`None`), an em dash -- never a
    recomputed value (D12/D33).
    """
    ordered = sorted(rows, key=_sort_key)
    table = data_table(
        _DENOVO_HEADERS,
        [
            Row(
                cells=(
                    _CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class),
                    _count(row.n_case_carriers, row.n_cases, row.count_unit),
                    f"{_fmt(row.expected_count)}" if row.expected_count is not None else _EM_DASH,
                    _denovo_enrichment(row),
                    _statistic(row),
                    _denovo_reading(row),
                )
            )
            for row in ordered
        ],
    )
    return (
        f'<figure class="denovo" data-effect-measure="enrichment_ratio">'
        f'<div class="table-scroll">{table}</div>'
        "<figcaption><p>De novo enrichment against a modelled expectation (Poisson). "
        "The corrected <em>q</em> carries the study's own multiple-testing correction; "
        "where the study published none, the raw <em>p</em> is shown and "
        '<a href="../methods.html#rare-variants">the scan it came from</a> is the '
        "reading a null result is judged against.</p></figcaption></figure>"
    )


def _burden_section(
    rows: Sequence[BurdenRow],
    publications: Mapping[str, Publication],
    cohorts: Mapping[str, Cohort],
    families: tuple[frozenset[str], ...] = (),
    axes: tuple[tuple[str, str], ...] = (),
    combined: Sequence[CombinedRow] = (),
) -> str:
    """Published rare-variant burden, one table per study.

    Three things this section has to say that the numbers alone do not, each of
    which a reader would otherwise get wrong:

    * **The rows are per study and are never pooled.** The CHD literature reuses
      cohorts, so a pooled p-value counts the same people twice. Declining to
      compute one does not stop a reader doing it by eye, so where two studies
      share a collection the section names it (`shared_cohorts`).
    * **A consequence class with no row is one the study did not report** -- not
      "was tested and found nothing", and *why* it is absent is the study's own
      rule rather than a fact about the data. This bullet claimed the stronger
      thing until 2026-08-05: "had no carrier in either group... measured over
      the committed mirror: zero of its 1,192 rows have no case carrier *and* no
      control carrier". That measurement was over PMID:42230622 alone, and both
      halves have since failed. The mirror now holds 1,475 rows and **one of
      them does have zero carriers on both sides**; and PMID:40127276 observed
      14,364 synonymous variants (its Dataset S4) while publishing no synonymous
      row at all, because its gene-level table covers only damaging classes. The
      matrix is genuinely sparse, so a reader meets a gap often enough for the
      distinction to matter -- which is why the sentence must be weak enough to
      be true of every study on the page.
    * **The synonymous row is the study's own negative control**, and only a
      study that publishes one has it. Synonymous variants should show no
      enrichment; where one does, that gene's comparison is poorly calibrated.
      It is sorted last within each stratum so it reads as what it is -- the row
      that says whether to believe the two above it. 69 of the 290 published
      rows are synonymous and every one is from PMID:42230622, so
      `_SYNONYMOUS_NOTICE` was conditional: a page with no such row must not
      promise a check it cannot offer.

    The last two bullets' sentences left the page on 2026-10-26 (D50): they
    were the folded reading notes, identical wherever they rendered -- the
    absent-row sentence on all 118 pages with a rare-variant section, the
    synonymous one on the 113 of them with a synonymous row -- and are now on
    `methods.html#rare-variants` (`methods._RARE_VARIANTS`), stated once with
    their conditions in their own words. The first stays, as the pooling
    caption and the `shared_cohorts` notice.

    Returns `""` for a gene with no burden rows rather than an empty section:
    unlike the validity section, whose empty group is a sentence stating that no
    record from that authority is mirrored (D51) -- a fact about the mirror, and
    an answer in itself -- an empty burden table would say "this gene was
    studied and nothing was found", which is a claim no study made.
    """
    if not rows:
        return ""

    overlaps = shared_cohorts(rows, cohorts)
    studies = sorted({row.study for row in rows})
    # One block per study, and each block's heading is the card's link target
    # (D55): two studies sharing an anchor would send one study's caveats link
    # to the other's block. Distinct PMIDs cannot collide; this refuses the
    # bypassed gate where an id is not one.
    if len({_study_anchor(study) for study in studies}) != len(studies):
        raise ValueError(f"two studies share a page anchor: {studies}")
    blocks: list[str] = []
    for study in studies:
        study_rows = [row for row in rows if row.study == study]
        study_combined = [row for row in combined if row.study == study]
        table = data_table(
            _BURDEN_HEADERS,
            [
                Row(
                    cells=(
                        _STRATUM_LABEL.get(row.cohort_stratum, row.cohort_stratum),
                        _CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class),
                        _count(row.n_case_carriers, row.n_cases, row.count_unit),
                        _controls(row),
                        _effect(row),
                        f"{row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH,
                        _corrected(row),
                    )
                )
                for row in study_rows
            ],
        )
        design, warning = _method_line(study_rows, publications.get(study))
        # `design` and `_provenance` fold together: both are identical for a
        # given study on all 91 gene pages carrying one. `warning`, the disclosure, the
        # composite note and the footnotes do not -- each changes what a reader
        # concludes from the table right below it.
        counted = (
            f'<details class="reading-notes"><summary>How this study counted'
            f"</summary>{design}{_provenance(study_rows, cohorts)}</details>"
            if design or _provenance(study_rows, cohorts)
            else ""
        )
        # One panel per effect measure, never one per study: an odds ratio of
        # 3.1 and a de novo enrichment of 3.1 are different claims, and 75 of
        # the 915 published rows are the second kind. `sorted` over a set, for
        # the reason every other loop here is sorted. De novo enrichment
        # (mutation_model) rows draw as the structured consequence table rather
        # than a forest, which suits point estimates without confidence
        # intervals and keeps the all-null panel (`p 1 · q 1` three times)
        # legible instead of three identical flat bars.
        placed = [row for row in study_rows if _plottable(row)]
        measures = sorted({row.effect_measure for row in placed if row.effect_measure is not None})
        # Every (stratum, consequence) cell this study published for this gene,
        # plotted or not: what a panel's unsplit-union line reads to say a
        # missing component was never published rather than not drawn.
        cells = frozenset((row.cohort_stratum, row.consequence_class) for row in study_rows)
        measure_rows = {
            measure: [row for row in placed if row.effect_measure == measure]
            for measure in measures
        }
        panels = "".join(
            _denovo_table(measure_rows[measure])
            if measure == "enrichment_ratio"
            else _forest(
                measure_rows[measure],
                measure,
                _study_label(study, publications),
                cells,
            )
            for measure in measures
        )
        # **Folded only where a panel above it shows the same rows.** THE FOLD
        # RULE says a caveat may fold when it is general and must stay visible
        # when it is particular; a table of this gene's own counts is as
        # particular as anything on the page, and what earns the fold is not
        # that it repeats but that the picture directly above it is drawn from
        # exactly these rows. Where no row is plottable there is no picture, so
        # the table is the only record and stays in plain sight -- every
        # PMID:34324492 block is in that state (77 pages, measured 2026-10-26;
        # for KLF13, CFC1, HRAS and MYH11 it is the only study), and folding
        # away the only figures it has would be evidence loss dressed as
        # tidiness.
        figures = (
            f'<details class="study-figures"><summary>every figure from this study'
            f"</summary>{table}</details>"
            if panels
            else table
        )
        analyses = html.escape(
            " · ".join(sorted({_analysis_label(r.comparator, r.variant_class) for r in study_rows}))
        )
        blocks.append(
            f'<div class="study-card">'
            f'<h3 id="{html.escape(_study_anchor(study))}">'
            f'<a href="{html.escape(_pubmed(study))}">'
            f"{html.escape(_study_label(study, publications))}</a></h3>"
            f'<p class="study-analyses">{analyses}</p>'
            f"{_study_finding_summary(study_rows)}"
            f"{_disclosure(study, publications)}"
            f"{warning}"
            f"{counted}"
            f"{_composite_note(study_rows)}"
            f"{_combined_line(study_combined)}"
            f"{panels or _NO_EFFECT_MEASURE}"
            f"{figures}"
            f"{_footnotes(study_rows, study_combined)}"
            f"</div>"
        )

    shared = "".join(
        _overlap_notice(left, right, overlap, publications, cohorts)
        for (left, right), overlap in overlaps.items()
    )

    # The pooling sentence is conditional. Rendered unconditionally it asserted
    # "these cohorts overlap" in the present tense on a page showing one table,
    # sending a reader to look for a second study that is not there and
    # undermining the caveats that *are* live.
    pooling = _POOLING_NOTICE if len({row.study for row in rows}) > 1 else ""

    matrix = _evidence_matrix(rows, families, axes, publications)
    # The card lists a row exactly when one survives (`_card_rare_variants`).
    lead_differs = any(survives_own_correction(row) for row in rows)

    # THE FOLD RULE. A caveat may fold if it is *general* -- true on every gene
    # page, teaching how to read the table. It must stay visible if it is
    # *particular* -- it fired because of this gene's data and changes what a
    # reader concludes about this gene.
    #
    # The general reading notes -- the units glossary, what an absent cell means,
    # and the synonymous negative control -- were folded here under "How to read
    # these numbers": the first two on all 118 pages carrying a rare-variant
    # section, the third on the 113 of them with a synonymous row.
    # They moved to `methods.html#rare-variants` on 2026-10-26 (D50), which
    # `_MIRRORED_NOTICE` links; none is lost -- the unit word is still in every
    # cell and the matrix renders `not tested` as its own state.
    #
    # What stays visible is unchanged, outside every fold: `_MIRRORED_NOTICE`
    # (one line, and it is the atlas's claim about itself), the pooling caption,
    # the `shared_cohorts` non-independence notice, the own-lab disclosure, the
    # method line's multiple-testing clause, and every row's `method_note`
    # footnote -- which is where PMID:34324492's thoracic-aortic-aneurysm
    # contamination reaches a reader.
    # `test_every_particular_caveat_survives_outside_the_fold` is what stops one
    # of those being tidied inside later.

    return (
        '<h2 id="rare-variants">Rare-variant association</h2>'
        + _MIRRORED_NOTICE
        + matrix
        # RELOCATED, not folded. A summary that counts datasets is adjacent to
        # pooling and invites exactly the arithmetic this sentence forbids, so it
        # belongs directly beneath the matrix rather than four paragraphs above
        # it. Rendered whenever there is a matrix, because the matrix is what it
        # qualifies -- unlike its previous home, where it was conditional on
        # there being more than one study to pool. The lead note beside it is
        # conditional on the card listing a surviving row; see its constant.
        + (_POOLING_NOTICE + (_MATRIX_LEAD_NOTE if lead_differs else "") if matrix else pooling)
        + shared
        + "".join(blocks)
        + _cohort_notes(rows, cohorts)
    )


def _controls(row: BurdenRow) -> str:
    """The control column, which a de novo row fills with an expectation.

    `_count(None, None)` renders an em dash, and on a `mutation_model` row that
    dash is where the *entire comparator* belongs: the modelled expected count
    is the only thing the enrichment was computed against, and an em dash is
    indistinguishable from a control count nobody recorded.
    """
    if row.comparator == "mutation_model" and row.expected_count is not None:
        return f"{_fmt(row.expected_count)} expected"
    return _count(row.n_control_carriers, row.n_controls, row.count_unit)


def _disclosure(study: str, publications: Mapping[str, Publication]) -> str:
    """Say so when the study is by an author of this atlas.

    `Publication.own_lab` exists to record exactly this and reached
    `publications.json` and zero HTML files. Two of the three studies published
    today -- `PMID:42230622` and `PMID:34324492` -- are by this repository's
    author, and a reader weighing those numbers should be told without having to
    cross-reference a JSON payload.

    **The conditional is now load-bearing rather than latent.** This docstring
    said "every burden table on the site today is from PMID:42230622" until
    2026-08-05, when `PMID:40127276` (`own_lab: false`) began rendering on 15 of
    the 23 gene pages published then. Re-measured 2026-08-06 against the widened
    population: 91 of the 92 pages carry a burden table, drawn from three studies
    -- PMID:42230622 on 87, PMID:34324492 on 61, PMID:40127276 on 25 -- and 13
    carry all three, with the own-lab declaration on **two** of the three
    (PMID:42230622 and PMID:34324492; PMID:40127276 is not this atlas's). That is
    the arrangement that makes the line informative: an unconditional disclosure
    says nothing, and a stale claim that it is unconditional invites someone to
    make it one. This sentence read "one of the three" until 2026-08-06, which
    contradicted this docstring's own opening paragraph two screens up.

    **One sentence since 2026-10-26 (D50)**; `methods.html#rare-variants` says
    the rest. The spec's approved short form ended "nothing added", and that is
    false where it sits: measured on a build of 4233eca, 113 of the 190 own-lab
    study blocks also carry the atlas's own Bonferroni arithmetic
    (`_method_line`'s `0.05 / tests_reported`) and a forest plot the atlas drew
    from the study's rows.

    **"no assessment added" went too (review, 2026-10-27).** The same day,
    `methods.html#rare-variants` began stating the one rule the atlas applies
    to every study's rows (owner decision W2): which support the gene and
    which survive their study's correction, the rule the card, the matrix and
    the forest fill read. A block saying the atlas added no assessment
    contradicted that page on all 190 own-lab blocks. The sentence says what
    is true instead: the rule is applied, and it is the same rule
    (`concordance.survives_own_correction` reads no `own_lab`).
    """
    publication = publications.get(study)
    if publication is None or not publication.own_lab:
        return ""
    return (
        '<p class="disclosure">Declaration: this study is by an author of this atlas; '
        "its results are republished as published and read by the same rule as every other "
        "study's.</p>"
    )


def _labelled(values: Iterable[str | None], labels: Mapping[str, str]) -> list[str]:
    """Distinct values of one column, mapped to display text, sorted.

    Sorted because these come from a `set`, whose iteration order varies with
    `PYTHONHASHSEED` -- and this string is part of a page's bytes and therefore
    its checksum. Sorted on the label rather than the raw value so the rendered
    line reads in the order it is shown.
    """
    return sorted({labels.get(value, value) for value in values if value is not None})


def _method_line(rows: Sequence[BurdenRow], publication: Publication | None) -> tuple[str, str]:
    """What was counted, and against how many other tests.

    **Three of the eight partition columns reached no page.** `variant_class`,
    `origin` and `maf_max` were populated and appeared nowhere, so a section
    headed "Rare variant burden" (now "Rare-variant association", D49) never said
    what *rare* meant, and a
    case-control count of rare-inherited variants was indistinguishable from a
    trio's de novo count. `tables.py` states the case against exactly this: the
    partition "exists to stop two incomparable rows from *looking* comparable".

    This said "four" and named `lesion_group` as the fourth until 2026-08-05.
    **`lesion_group` is null on all 290 published rows and always has been** --
    no curated study stratifies by lesion group yet -- so a quarter of the
    stated rationale was for a column that had nothing to render. `_LESION_LABEL`
    is deliberately left empty rather than removed: the column is live in the
    schema and a study that uses it should render, at which point that map needs
    filling. Measured, not assumed; the claim it replaces was not.

    Rendered per study rather than per row, from the distinct values across the
    block. **A field that is multi-valued within a study renders every value**,
    which says what the block contains without claiming which row is which --
    the JSON is row-precise and the preamble says so. That is no longer
    hypothetical: PMID:40127276 contributes two origins and two tests to one
    block, and its `maf_max` is populated on the case-control rows and null on
    the de novo rows, which is why the frequency clause below distinguishes
    "filtered at this threshold" from "not filtered at all" rather than dropping
    the null.

    The multiple-testing sentence is the other half, and the more serious one.
    Measured 2026-08-05 over the built site: 32 of PMID:42230622's 187 published
    rows have p < 0.05 and **3** survive Bonferroni over that study's own
    138,609 tests, so 29 rows read as significant and are not. That study
    publishes no corrected p, and computing one would be authoring a statistic
    (D12/D33) -- but naming the denominator is a fact the study supplies, and it
    is what lets a reader apply their own threshold. The other two studies do
    publish corrections, so they get the other branch below.
    """
    parts: list[str] = []
    for values, labels in (
        ([row.variant_class for row in rows], _VARIANT_CLASS_LABEL),
        ([row.origin for row in rows], _ORIGIN_LABEL),
        ([row.lesion_group for row in rows], _LESION_LABEL),
        ([row.pvalue_test for row in rows], _TEST_LABEL),
    ):
        rendered = _labelled(values, labels)
        if rendered:
            parts.append(", ".join(rendered))
    # A null `maf_max` means "this analysis applied no frequency filter", which
    # is a fact about the rows and not an absence of one. Dropping it and then
    # stating the surviving threshold over the whole block attributed
    # PMID:40127276's MAF < 1e-05 to its 45 de novo rows, which had no filter at
    # all -- a de novo mutation is defined by absence from both parents, not by
    # a population frequency. Measured 2026-08-05: that sentence was wrong on 15
    # of the 23 published gene pages.
    mafs = sorted({row.maf_max for row in rows if row.maf_max is not None})
    unfiltered = any(row.maf_max is None for row in rows)
    if mafs and unfiltered:
        parts.insert(
            1,
            "MAF below "
            + ", ".join(f"{value:g}" for value in mafs)
            + " for the rows that applied a frequency filter, and no filter on the rest",
        )
    elif mafs:
        parts.insert(1, "MAF below " + ", ".join(f"{value:g}" for value in mafs))
    elif unfiltered:
        parts.insert(1, "no frequency filter")

    sentence = html.escape("; ".join(parts)) + "." if parts else ""

    # Three cases, because a sentence true of one study is false of another.
    # PMID:42230622 publishes no correction at all, so the reader needs the
    # denominator; PMID:34324492 corrects every row, so telling them to judge
    # the raw p against the scan would be pointing at a column they should not
    # be reading. A study with any uncorrected row gets the denominator, worded
    # for the rows that lack a correction.
    #
    # **The count is shown only where the reader must apply it** (2026-10-26).
    # Both branches opened "This study reported N comparisons in all", and
    # `tests_reported` did not always count comparisons: PMID:40127276's was
    # 248, its panel *genes*, each tested in two arms over three consequence
    # classes, so the sentence was false on all 36 pages carrying that study.
    # (Its `tests_reported` is null since owner decision 2026-10-26, W1, and
    # the 248 is `correction_genes`.) Where a study corrected every row the reader
    # applies no threshold of their own, so no count is needed, and none is
    # stated; where it did not, the count is the Bonferroni denominator the
    # sentence after it uses -- true of PMID:42230622, whose 138,609 is one
    # Fisher test per row of the three sheets of its Supplementary Data 3
    # (`curation/publications.yaml`). The corrected sentence therefore no
    # longer reads `tests_reported` at all, and renders whether or not the
    # publication records one -- but only where a corrected p exists to read.
    correction = ""
    uncorrected = [row for row in rows if row.pvalue is not None and row.pvalue_adjusted is None]
    if uncorrected and publication is not None and publication.tests_reported:
        total = publication.tests_reported
        threshold = 0.05 / total
        # "{total:,} such comparisons" until 2026-10-26: "such" pointed back at
        # the design clause, which renders after this one and folded under "How
        # this study counted", so its antecedent was a sentence the reader had
        # not met.
        correction = (
            f" This study reported {total:,} comparisons in all. Where it published no "
            f"corrected p, judge the raw one against the whole scan rather than against "
            f"0.05 &mdash; a Bonferroni threshold over that many is p &lt; {threshold:.2g}."
        )
    elif uncorrected:
        correction = " p-values are <strong>uncorrected</strong> as published."
    elif any(row.pvalue_adjusted is not None for row in rows):
        # Names the statistic, not "the corrected p column" (P9, 2026-10-26):
        # where the study has a forest plot that column is folded under "every
        # figure from this study", and the plot names the same number `q` where
        # it is a false-discovery rate (`_statistic`). True with a forest and
        # without one.
        correction = (
            " Its <strong>corrected p</strong> (shown as <em>q</em> where it is a "
            "false-discovery rate) carries the study's own correction for multiple testing; "
            "read it rather than the raw p."
        )

    if not sentence and not correction:
        return ("", "")

    # **The two halves are returned separately because they fold differently.**
    #
    # The design clause ("SNVs and indels; MAF below 0.001; Fisher's exact
    # test") is *general*: identical for a given study on all 23 gene pages, and
    # it teaches how to read the block rather than saying anything about this
    # gene. It folds, with the provenance line, under "How this study counted".
    #
    # The multiple-testing clause does not. It was the sharpest finding of the
    # 2026-08-05 review -- 32 of PMID:42230622's 187 published rows clear 0.05
    # and 3 survive Bonferroni over its own 138,609 comparisons -- so it changes
    # what a reader concludes from every p-value in the block below it. Folding
    # it would reintroduce, one release later, exactly the defect that review
    # found.
    #
    # Returned as a pair rather than as one string the caller splits: this
    # module builds HTML and must never parse it back, and a split on "</p>"
    # would break silently the first time either half gained a paragraph.
    design = f'<p class="method">{sentence}</p>' if sentence else ""
    warning = f'<p class="method">{correction.strip()}</p>' if correction else ""
    return (design, warning)


def _provenance(rows: Sequence[BurdenRow], cohorts: Mapping[str, Cohort]) -> str:
    """The collections a study's rows drew on, unioned across all of them.

    Read off `study_rows[0]` until 2026-08-05, which stated one row's cohorts
    over the whole block. `case_cohorts` and `control_cohorts` are per-row
    columns, so a study whose strata draw on different collections published a
    provenance line the data contradicted -- and `shared_cohorts` reads the same
    columns with a *union* rule, so the "not independent" notice could name a
    cohort the provenance line did not. Both now union, so they cannot disagree.
    """
    cases = sorted({identifier for row in rows for identifier in row.case_cohorts})
    controls = sorted({identifier for row in rows for identifier in row.control_cohorts})
    return (
        f'<p class="provenance">Cases: {html.escape(_names(cases, cohorts))}. '
        f"Controls: {html.escape(_names(controls, cohorts) or 'none (see the method above)')}.</p>"
    )


def _concordance_for(
    concordance: Mapping[str, Mapping[str, Json]] | None, gene: str
) -> Mapping[str, Json]:
    """One gene's concordance, or a refusal naming the gap.

    **`bundles._concordance_for` raised on a missing gene and this page defaulted
    it**, so a mapping built over the wrong population would have made the bundle
    fail loudly and the browse row render "0 of 0 tested" -- a tally byte-identical
    to a measured "no study reported this gene". The guard existed on one of the
    two layers `runner.py` says cannot disagree.

    `None` is still allowed and means "this build has no burden data at all",
    which is what every unit test of this page passes. A *mapping* that omits a
    gene is the error: it means the caller derived it over a different population.
    """
    if concordance is None:
        return {"tested": 0, "enriched": 0, "corrected": 0, "families": []}
    if gene not in concordance:
        raise KeyError(
            f"no concordance derived for published gene {gene}; the mapping was built over "
            f"a different population, and defaulting it would render 'no dataset tested "
            f"this' as if it had been measured"
        )
    return concordance[gene]


def _strip_notes(
    genes: Iterable[str],
    concordance: Mapping[str, Mapping[str, Json]] | None,
    publications: Mapping[str, Publication],
) -> str:
    """Each curated `Publication.summary_caveat` of a study the strips name, once.

    The strips summarise every study's result for every gene -- PMID:34324492
    sits in a family on all of them, and a hover reading "PMID:34324492:
    enriched, and survives this study's own correction" is a stronger summary
    than any count -- so the caveat the
    gene page's card prints beside that study's counts is printed here too,
    under the legend, in study-id order. The page has no study blocks to link
    to, so each note names its study by label and PMID.
    """
    # Keys of a dict, not a set: met in row order, so the `sorted` below is what
    # orders them, and a test can tell (CLAUDE.md section 4.42).
    studies: dict[str, None] = {}
    for gene in genes:
        families = _concordance_for(concordance, gene).get("families") or ()
        if not isinstance(families, Sequence):
            continue
        for entry in families:
            named = entry.get("studies") if isinstance(entry, Mapping) else None
            if isinstance(named, Sequence) and not isinstance(named, str):
                studies.update(dict.fromkeys(str(study) for study in named))
    notes: list[str] = []
    for study in sorted(studies):
        publication = publications.get(study)
        if publication is None or publication.summary_caveat is None:
            continue
        notes.append(
            f'<p class="strip-legend"><strong>'
            f"{html.escape(_study_label(study, publications))}</strong> "
            f"({html.escape(study)}): {html.escape(publication.summary_caveat)}</p>"
        )
    return "".join(notes)


def _dot_strip(concordance: Mapping[str, Json]) -> str:
    """One glyph per cohort family, plus a tally naming both statistics.

    **Fill encodes the correction**, so a single glyph carries two facts a
    reader would otherwise need two columns for: that a dataset showed
    enrichment, and whether that survived the study's own correction. TAB2 is
    the case it was designed against -- every family pointing the same way,
    and fewer surviving correction: it read "3 of 3 tested · 1 corrected"
    over three families until owner decision 2026-10-27 (I3) made two of its
    three studies one family, and reads "2 of 2 tested · 1 corrected" since.

    The tally counts *tested* families and says both numbers, because either
    alone misleads. `2 of 2` hides that one does not survive correction;
    `1 corrected` hides that PMID:42230622 publishes no correction on any row
    and so can never earn one.

    A family that did not test the gene gets its own glyph and is excluded from
    the denominator. That is not presentation: KDM6A is ClinGen definitive and
    shows nothing in either dataset that tested it, and a denominator that also
    counted a family which never looked would add one more dataset to those
    finding nothing, beside a green chip. KLF13, which one of the two families
    tested, reads "1 of 1 tested" (measured on a build of 5f7161b).

    **The hover names each study with its own state** (`members`, review
    2026-10-27). The glyph and the tally are the family's; the `title` is
    where a reader asks which study said what, and it built "<studies>:
    <family state>" until then, which attributed a family's state to every
    study in it -- "enriched, and survives this study's own correction" of a
    study that publishes no corrected p, and a state at all of a study with no
    row for the gene. A family of two or more says it is one, counted once.
    """
    families = concordance.get("families") or ()
    if not isinstance(families, Sequence):
        return _EM_DASH
    glyphs: list[str] = []
    for entry in families:
        if not isinstance(entry, Mapping):
            continue
        state = str(entry.get("state", ""))
        members = entry.get("members") or ()
        said = [
            f"{member.get('study', '')}: "
            f"{_STATE_TITLE.get(str(member.get('state', '')), str(member.get('state', '')))}"
            for member in (members if isinstance(members, Sequence) else ())
            if isinstance(member, Mapping)
        ]
        # ". " between studies, not "; ": the `nominal` state's own sentence
        # carries a semicolon.
        title = ". ".join(said) + (". One cohort family, counted once." if len(said) > 1 else "")
        glyph = _STATE_GLYPH.get(state, "untested")
        body = "&ndash;" if glyph == "untested" else ""
        glyphs.append(f'<span class="dot {glyph}" title="{html.escape(title)}">{body}</span>')
    tested = concordance.get("tested", 0)
    corrected = concordance.get("corrected", 0)
    enriched = concordance.get("enriched", 0)
    tally = f"{enriched} of {tested} tested &middot; {corrected} corrected"
    # Not a link. The row already links to the gene page twice -- from the id and
    # from the symbol -- and a third link wrapped around a tally puts the
    # affordance on a number rather than on a name. Each dot still carries a
    # `title` naming its studies and state, which is the detail this cell owes a
    # reader; the numbers behind it are one click away by the name they scanned
    # for.
    return f'<span class="strip">{"".join(glyphs)}<span class="strip-tally">{tally}</span></span>'


def _evidence_matrix(
    rows: Sequence[BurdenRow],
    families: tuple[frozenset[str], ...],
    axes: tuple[tuple[str, str], ...],
    publications: Mapping[str, Publication],
) -> str:
    """Cohort family against evidence design, with the holes left visible.

    **The empty cells are the point.** Columns come from the whole corpus rather
    than from this gene, so a dataset that did not test this gene renders as an
    absence rather than vanishing -- the same confusion between "not tested" and
    "tested and found nothing" that `FamilyState` has four members for.

    **The guarantee is narrower than an earlier version of this docstring
    claimed.** It said the CNV de novo quadrant "is empty on every page". It is
    on *no* page: `evidence_axes` derives columns from the rows that exist, so a
    design *nobody* ran produces no column at all. What the matrix shows is a
    hole where one family did not run a design another family did; it is silent
    about a design the whole corpus lacks. Measured 2026-08-05: three columns,
    none of them CNV de novo.

    The cell shows the most significant of the rows that give it its state
    (`concordance.deciding_rows`) -- the surviving rows of a `corrected` cell,
    the supporting rows of a `nominal` one -- ordered by corrected p where
    published and by raw p otherwise, tie-broken on
    `(consequence_class, cohort_stratum)` so the choice is deterministic across
    builds. It chose among every row in the cell until 2026-10-26, and GJA1's
    de novo cell in PMID:40127276 named its `damaging` union (p 0.0613, q
    0.200) under the `nominal` glyph its missense row (p 0.0468) had earned:
    measured that day on a real build, 1 of the 109 enriched cells named a row
    its state is not true of. Where a cell holds more than one row the count is
    stated, so a reader knows the cell summarises and the table below is the
    record.

    **In a family of two or more studies, a study with no row for this gene is
    marked in the row's label** ("(no row for this gene)"), because the
    family's cells cannot say which of its studies they come from. Until
    review (2026-10-27) the label joined every member's name, so 46 of the 118
    pages with a matrix (measured on a build of 5f7161b) named a study as part
    of the row whose cells it contributes nothing to, and which has no block
    on the page. A one-study family is not marked: its cells already say
    `not tested` in every column.
    """
    if not rows or not families or not axes:
        return ""

    by_study: dict[str, list[BurdenRow]] = {}
    for row in rows:
        by_study.setdefault(row.study, []).append(row)

    header = "".join(
        f'<th scope="col">{html.escape(_VARIANT_CLASS_LABEL.get(variant, variant))}'
        f'<br><span class="sub">{html.escape(_COMPARATOR_LABEL.get(comparator, comparator))}'
        f"</span></th>"
        for variant, comparator in axes
    )

    body: list[str] = []
    for family in families:
        members = [row for study in sorted(family) for row in by_study.get(study, ())]
        label = " + ".join(
            _study_label(study, publications)
            + (" (no row for this gene)" if len(family) > 1 and study not in by_study else "")
            for study in sorted(family)
        )
        cells: list[str] = []
        for variant, comparator in axes:
            in_cell = [
                row
                for row in members
                if row.variant_class == variant and row.comparator == comparator
            ]
            cells.append(_matrix_cell(in_cell))
        body.append(f'<tr><th scope="row">{html.escape(label)}</th>{"".join(cells)}</tr>')

    return (
        '<div class="matrix-card">'
        '<div class="scroll"><table class="matrix">'
        f"<thead><tr><td></td>{header}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div>" + _MATRIX_LEGEND + "</div>"
    )


def _matrix_cell(rows: Sequence[BurdenRow]) -> str:
    """One family's result for one evidence design.

    An empty cell says `not tested` in words rather than rendering blank: a
    blank cell is indistinguishable from a rendering bug, and this one carries a
    claim -- that this dataset did not run this design on this gene.
    """
    if not rows:
        return (
            f'<td><span class="cell not-tested" title="{_STATE_TITLE["not_tested"]}">'
            f"{_STATE_CELL['not_tested']}</span></td>"
        )

    state = family_state(rows)
    if state is FamilyState.NO_ENRICHMENT:
        detail = f"n = {len(rows)}" if len(rows) > 1 else "&nbsp;"
        return (
            f'<td><span class="cell no-enrichment" title="{_STATE_TITLE["no_enrichment"]}">'
            f'{_STATE_CELL["no_enrichment"]}<br><span class="sub">{detail}</span></span></td>'
        )

    # Most significant first among the rows the state is true of -- never a
    # row that is not support, or not a survivor under a `corrected` glyph --
    # by corrected p where published, raw p otherwise, then a stable tie-break
    # so two builds pick the same row.
    best = sorted(
        deciding_rows(rows),
        key=lambda row: (
            row.pvalue_adjusted if row.pvalue_adjusted is not None else (row.pvalue or 1.0),
            row.consequence_class,
            row.cohort_stratum,
        ),
    )[0]

    # **The cell is a summary; the table below is the record.** So it carries the
    # effect and one statistic and stops there. The full spelling -- the
    # confidence interval, the name of the correction -- goes on the `title` and
    # is in the table in full a screen further down.
    #
    # This is a width decision with a correctness edge. The interval and the
    # correction name ran a cell to 61 characters, which set the column width for
    # the whole grid and crowded the cells against each other. What may *not* be
    # dropped is the measure: `_effect` has no branch that omits it, because an
    # odds ratio of 3.1 and a de novo enrichment of 3.1 are different claims, and
    # `_effect_compact` keeps that property.
    kind = "corrected" if state is FamilyState.CORRECTED else "nominal"
    # `_statistic` is shared with the forest panel above, so the two pictures on
    # one page cannot name the same number two different ways; the rule it
    # keeps -- `q` only for a real false-discovery rate -- is recorded there.
    statistic = _statistic(best)
    more = f'<span class="sub"> &middot; {len(rows)} rows</span>' if len(rows) > 1 else ""
    full = f"{_effect(best)}; {_corrected(best) if best.pvalue_adjusted is not None else ''}".strip(
        "; "
    )
    return (
        f'<td><span class="cell {kind}" title="{html.escape(full)}">'
        f"{html.escape(_effect_compact(best))}<br>"
        f'<span class="sub">{html.escape(statistic)}</span>{more}</span></td>'
    )


def _effect_compact(row: BurdenRow) -> str:
    """The effect and its measure, without the interval. For the matrix cell only.

    **The measure survives; the interval is dropped, except for an unbounded
    effect.** `_effect`'s rule -- that no branch omits the measure, because an
    odds ratio of 3.1 and a de novo enrichment of 3.1 are different claims --
    is the one property this must not trade for width. The interval is a
    different kind of thing: it qualifies the estimate rather than naming it,
    and it is in the row's own table a screen below, plus on this cell's
    `title`.

    **An unbounded effect keeps its lower bound: `OR ∞ (CI 0.786–∞)`**, never a
    bare `OR ∞` and never a `≥` (see `_effect`). Fisher returns an infinite
    odds ratio wherever no control carries, which says nothing about how large
    the effect is: the lower bound is the finding, and alone `∞` reads as
    "infinitely enriched". This rendered the bare form until 2026-10-26, and
    on TBX5's and ZEB2's cells "OR ∞" over "p 0.0316" stood for an interval of
    0.786–∞ that spans 1, visible only on hover. Measured that day: 37 of the
    1,805 mirrored rows are unbounded above and 20 of them have a lower bound
    below 1 (29 and 14 of the 1,210 published), and 9 matrix cells name one.
    The bare form is the fallback for an unbounded row with no lower bound, of
    which there are none.
    """
    if row.effect_measure is None:
        return _EM_DASH
    measure = _MEASURE_LABEL.get(row.effect_measure, row.effect_measure)
    if row.effect_bound == "unbounded_above":
        low = f" (CI {_fmt(row.ci_low)}–∞)" if row.ci_low is not None else ""
        return f"{measure} ∞{low}"
    if row.effect is None:
        return _EM_DASH
    return f"{measure} {_fmt(row.effect)}"


# --- Effect and uncertainty, drawn ------------------------------------------
#
# The evidence matrix a screen above answers *which datasets tested this gene
# and which found enrichment*. It cannot answer **how large** and **how
# certain**, because a swatch has no room for a number and no room for an
# interval. That is what this panel is for, and it is the half of the burden
# layer a reader was previously left to assemble from a nine-row table by eye.
#
# **Every (study, effect measure) pair gets a panel, including the two studies
# that publish no interval at all.** Measured 2026-08-20 over the 915 published
# rows: PMID:42230622 publishes an interval on 100% of its 704 rows and
# PMID:40127276 on 0% of its 150 -- and PMID:40127276 is the study with the
# findings that survive correction. On TBX5 it reports de novo loss-of-function
# enriched 297x at q 5.6e-08, while none of PMID:42230622's eight TBX5 rows
# carries a corrected p at all and seven of their intervals cross 1. The eighth
# does not, and this comment said all eight did until it was measured: the
# non-syndromic synonymous row sits entirely below 1, at OR 0.251 (95% CI
# 0.030-0.924, p 0.0375). It is that study's own negative control, so it is the
# one interval a sentence about this panel least ought to misdescribe.
#
# Drawing only where an interval exists would put a picture of the uncorrected
# result on the page and leave the surviving finding as text -- curated
# evidence visually demoted, which is this repository's characteristic failure
# in a new medium.

# One panel's geometry, in the SVG's own coordinates. The label column is wide
# because the longest label this vocabulary produces -- "non-syndromic ·
# damaging (LOF + missense) union" -- is 47 characters, and a label running
# under the plot area would overlap the very interval it names.
_FOREST_WIDTH: Final = 600
_FOREST_ROW: Final = 20.0
_FOREST_TOP: Final = 26.0
_FOREST_LABEL_X: Final = 2.0
_FOREST_PLOT_LEFT: Final = 244.0
_FOREST_PLOT_RIGHT: Final = 470.0
_FOREST_STAT_X: Final = 478.0
_FOREST_AXIS_GAP: Final = 10.0
_FOREST_FOOT: Final = 16.0

# Text is dropped by an explicit offset rather than by `dominant-baseline`,
# whose initial value differs between renderers: a chart whose labels sit half
# a row too high in one browser and correctly in another is a chart nobody can
# check against the table below it.
_FOREST_TEXT_DROP: Final = 3.4
_FOREST_ARROW: Final = 7.0
_FOREST_ARROW_HALF: Final = 3.5

# A quarter decade of clearance beyond the extreme values, so that no measured
# number lands on an edge of the axis. The edges are **reserved**: an arrow
# there means "past the end of this axis", and a real value sharing that pixel
# would read as one. Same reason `_FLOOR_GAP` exists on the trajectory chart --
# the atlas's weakest statement about a figure and its most extreme real one
# must not land on the same pixel column.
_EFFECT_AXIS_PAD: Final = 10.0**0.25

# What "survived the correction" means, and the only threshold this module
# applies to any number. It is applied to a p-value the **study** corrected and
# published; the atlas corrects nothing (D12/D33), which is why a row whose
# study published no correction at all is drawn hollow rather than filled --
# nothing says it survived one. A non-synonymous row's fill reads
# `concordance.survives_own_correction` instead, against the same 0.05
# (`concordance.ALPHA`); see `_survived_correction`.
_CORRECTED_ALPHA: Final = 0.05

# The two unions this schema contains, each as (the union's own value, the
# components it is the union of). `_composite_note` states the first in words
# and `_CONSEQUENCE_LABEL["damaging"]` names it in the table; the second is
# stated by the counts themselves (measured on TBX5: 3,876 = 1,471 + 2,405).
# The plot has to mark both, for the reason the note exists: bars stacked as
# siblings read as independent findings, and here two of them are inside a
# third.
_CONSEQUENCE_UNION: Final[tuple[str, tuple[str, ...]]] = (
    "damaging",
    ("lof", "missense_damaging"),
)
_STRATUM_UNION: Final[tuple[str, tuple[str, ...]]] = (
    "all",
    ("syndromic", "nonsyndromic"),
)

# What a page says instead of drawing an empty frame (D42). An empty axis under
# a heading naming a study asserts that the study measured something and found
# it to be nothing, which is not what happened.
#
# It renders in PMID:34324492's block on **every** page carrying that study,
# because its CNV rows publish no effect measure of any kind -- not only on the
# genes whose *only* rows are that study's. This comment said "3 published
# genes -- KLF13, CFC1, MYH11" (measured 2026-08-20) until 2026-10-26, which
# counted the second set and read as the first. Measured that day on a build of
# 4233eca: **77 pages**, one block each, all PMID:34324492; 4 of them (KLF13,
# CFC1, HRAS, MYH11) carry no other study. CRIPTO carries no burden row at all
# and reaches this section never.
_NO_EFFECT_MEASURE: Final = (
    '<p class="method">This study publishes <strong>no effect measure</strong> for this '
    "gene &mdash; no odds ratio, enrichment or rate ratio &mdash; so there is nothing to "
    "place on an effect axis. Its counts and p-values are in the table below.</p>"
)


def _plottable(row: BurdenRow) -> bool:
    """Whether this row can be placed on an effect axis at all.

    An `unbounded_above` row **is** plottable and renders as an arrow: some of
    those rows carry the strongest signals in the mirror (TAB2: 5 syndromic
    carriers, 0 of 45,082 controls), and `ci_low` is the finding. Excluding
    them was an error in this plan's own first measurements.

    A row with an effect measure, no effect and no unbounded flag is refused
    rather than drawn as a bar with nothing on it. It does not occur in the
    mirror -- `validate/burden.py` pairs `effect` with `effect_measure` -- and
    a panel row with no estimate and no arrow would read as an estimate of
    zero, which is a different and much stronger claim.
    """
    return row.effect_measure is not None and (
        row.effect is not None or row.effect_bound == "unbounded_above"
    )


def _union_kinds(rows: Sequence[BurdenRow]) -> tuple[bool, bool]:
    """Whether each union is *shown as one* in this panel: (consequence, stratum).

    Conditional on the components being in the same panel, exactly as
    `_composite_note` is conditional on the study reporting both. Measured
    2026-08-20: every PMID:40127276 row is `cohort_stratum: all` and that study
    publishes no syndromic or non-syndromic row at all, so its `all` rows are
    the whole cohort rather than the sum of two rows a reader can see. Tagging
    them `union` would point at rows that are not on the page.

    **This answers for the panel, and a tag is worn by a row.** A panel can
    turn a union on for one consequence and then hang the tag on a row whose
    own components are not all drawn -- 10 panels do. `_missing_components`
    is what asks the per-row question, and `_unsplit_line` is where the answer
    reaches the panel's visible line.
    """
    composite, parts = _CONSEQUENCE_UNION
    consequence = any(row.consequence_class == composite for row in rows) and any(
        row.consequence_class in parts for row in rows
    )
    whole, strata = _STRATUM_UNION
    stratum = any(row.cohort_stratum == whole for row in rows) and any(
        row.cohort_stratum in strata for row in rows
    )
    return consequence, stratum


def _is_union(row: BurdenRow, kinds: tuple[bool, bool]) -> bool:
    """Whether this row contains other rows in its own panel."""
    consequence, stratum = kinds
    return (consequence and row.consequence_class == _CONSEQUENCE_UNION[0]) or (
        stratum and row.cohort_stratum == _STRATUM_UNION[0]
    )


def _missing_components(
    row: BurdenRow, rows: Sequence[BurdenRow], kinds: tuple[bool, bool]
) -> list[tuple[str, str]]:
    """The `(stratum, consequence)` cells this union contains and `rows` does not draw.

    Empty when **every** row this union contains is itself drawn on this
    panel, and for a row that is not a union on it. It was `_fully_split`, a
    `bool`, until review (2026-10-27) asked the panel's visible line to name
    what is missing (`_unsplit_line`). In reading order: the consequence
    union's parts, then the stratum union's, each in its `_..._UNION` tuple's
    order.

    `_union_kinds` asks its question of the panel; a tag is worn by a row, and
    the two are not the same question. A panel carrying `all cases · missense`
    beside `syndromic · missense` and `non-syndromic · missense` turns the
    stratum union on for the whole panel -- and then `all cases · loss-of-
    function`, whose non-syndromic half the study never published, wears the
    same tag with only one component under it.

    Measured 2026-08-20 on the built corpus: 288 rows carry the tag, **0 of
    them with neither component present**, so the tag is never false; but on
    10 panels a tagged row has one of its two components missing. All 10 are
    a (stratum, consequence) cell absent from the mirror entirely -- and a
    missing cell is not a null result, so nothing here may read as one.

    Both axes are checked with the *other* axis held fixed, because that is
    what the union is: `damaging` for `all cases` is the loss-of-function and
    damaging-missense rows **for `all cases`**, not for any stratum.
    """
    consequence, stratum = kinds
    missing: list[tuple[str, str]] = []
    if consequence and row.consequence_class == _CONSEQUENCE_UNION[0]:
        same_stratum = [other for other in rows if other.cohort_stratum == row.cohort_stratum]
        drawn = {other.consequence_class for other in same_stratum}
        parts = _CONSEQUENCE_UNION[1]
        missing += [(row.cohort_stratum, part) for part in parts if part not in drawn]
    if stratum and row.cohort_stratum == _STRATUM_UNION[0]:
        same_class = [other for other in rows if other.consequence_class == row.consequence_class]
        drawn = {other.cohort_stratum for other in same_class}
        strata = _STRATUM_UNION[1]
        missing += [(part, row.consequence_class) for part in strata if part not in drawn]
    return missing


def _forest_label(stratum: str, consequence: str) -> str:
    """A forest row's name, as its label and the panel's visible line both print it."""
    stratum_label = _STRATUM_LABEL.get(stratum, stratum)
    return f"{stratum_label} · {_CONSEQUENCE_LABEL.get(consequence, consequence)}"


def _unsplit_line(
    rows: Sequence[BurdenRow], kinds: tuple[bool, bool], published: Collection[tuple[str, str]]
) -> str:
    """One sentence naming each `union`-tagged row whose components are not all drawn.

    On the panel's visible line, where it fires (owner decision 2026-10-26,
    E1). It read "A row tagged union here is not fully split: at least one of
    the rows it contains is not plotted; the table below carries every row
    this study published for this gene" until review (2026-10-27), which
    pointed a reader at the table for a row that is in no table: on all 10
    panels it fires on (measured on a build of ad03128) the missing component
    is a cell the study never published for the gene, and TBX5's panel tags
    three rows of which only `all cases · loss-of-function` is unsplit. So it
    names the row, and what is missing from it, and says which of two things
    the absence is: `published` holds the cells of every row this study
    published for this gene, plotted or not, so a component the study
    published and this panel does not draw -- no effect to place, or another
    effect measure -- is said to be "not plotted here", never unpublished.
    Latent today: on none of the 10 panels is a missing component published.

    Several such rows share one sentence, "The A and B rows are unions, and
    this study published no C or D row for this gene", each component named
    once, in panel order: the 10 published panels each have one, and the
    forest key's fixture, with four, reads as a list rather than as four
    repetitions of one sentence.
    """
    unions: list[str] = []
    absent: list[str] = []
    undrawn: list[str] = []
    for row in rows:
        if not _is_union(row, kinds):
            continue
        missing = _missing_components(row, rows, kinds)
        if not missing:
            continue
        # Escaped: an unknown token renders as itself (`_STRATUM_LABEL.get`).
        label = html.escape(_forest_label(row.cohort_stratum, row.consequence_class))
        unions.append(f"<strong>{label}</strong>")
        for cell in missing:
            name = html.escape(_forest_label(*cell))
            named = undrawn if cell in published else absent
            if name not in named:
                named.append(name)
    if not unions:
        return ""
    several = len(unions) > 1
    clauses: list[str] = []
    if absent:
        names = absent[0] if len(absent) == 1 else f"{', '.join(absent[:-1])} or {absent[-1]}"
        clauses.append(f"this study published no {names} row for this gene")
    if undrawn:
        plural = len(undrawn) > 1
        clauses.append(
            f"the {_and_join(undrawn)} row{'s' if plural else ''} "
            f"{'they contain' if several else 'it contains'} "
            f"{'are' if plural else 'is'} not plotted here"
        )
    subject = (
        f"The {_and_join(unions)} rows are unions" if several else f"The {unions[0]} row is a union"
    )
    return f" {subject}, and {', and '.join(clauses)}."


def _survived_correction(row: BurdenRow) -> bool:
    """Whether this row's forest marker is drawn as having cleared **its own
    study's** correction for multiple testing.

    Never a correction this atlas computed -- there is none (D12/D33). A row
    whose study published no correction is `False` for the same honest
    reason a row that failed one is: nothing says it cleared anything.

    **For every row but a synonymous one this is
    `concordance.survives_own_correction`**, the rule the card and the
    evidence matrix on the same page read, so the three cannot disagree about
    a row (CLAUDE.md section 4.28). It read the corrected p alone until
    2026-10-26, so a depleted row with q below 0.05 would have been drawn
    filled while the card said none survives a correction its study
    published. 0 published rows differed between the two rules that day.

    **A synonymous row keeps the corrected p alone**, because the shared rule
    never passes one -- a negative control is not support -- and a surviving
    one is a warning about the comparison that `_mark_class` draws muted
    rather than hollow. One predicate, because it decides the fill of two
    different glyphs and the key's clauses about them.
    """
    if row.consequence_class == "synonymous":
        return row.pvalue_adjusted is not None and row.pvalue_adjusted < _CORRECTED_ALPHA
    return survives_own_correction(row)


def _survived_against_the_gene(row: BurdenRow) -> bool:
    """A non-synonymous row that cleared its study's correction without supporting the gene.

    Drawn hollow by `_survived_correction`, so the key's hollow clause --
    "did not survive" -- is false of it, and `_forest_caption` adds a clause
    wherever one is drawn. On a validated corpus that means a row whose
    statistic does not point to enrichment in cases: BUR018 refuses a
    corrected p smaller than its raw one, so a corrected p below 0.05 implies
    a raw p below it too. Latent: 0 published rows (measured 2026-10-26).
    """
    return (
        row.consequence_class != "synonymous"
        and row.pvalue_adjusted is not None
        and row.pvalue_adjusted < _CORRECTED_ALPHA
        and not _survived_correction(row)
    )


def _mark_class(row: BurdenRow, glyph: str) -> str:
    """Filled where the study's own correction was survived, hollow where not.

    "Survived" is `_survived_correction`: for every row but a synonymous one,
    survived *as support for the gene* -- the card's and the matrix's rule --
    so a depleted row that cleared its correction is hollow (2026-10-26).

    **`glyph` is the shape, and the fill rule is the same for every shape.**
    A panel draws two of them -- `point` for a placed estimate, `arrow` for a
    row the axis has no position for -- and each has the same three states,
    so the class names are composed here rather than written out per shape.
    `_forest_estimate` bypassed this function entirely for its two arrow
    branches and drew every arrow in `chart-arrow`, whose declaration was
    byte-identical to `chart-point`'s. Measured 2026-08-20 on the built
    corpus: **205 of the 854 plotted marks were arrows, not one of them had
    survived any correction, and every one was drawn in the fill that says it
    had.** 30 were corrected and failed, with q as high as 1.0, and 175 came
    from a study publishing no correction at all; **13 of the 30 sat on a
    ClinGen `definitive` gene**, GATA4 among them, whose Sierant enrichment
    panel drew three solid marks of which the third was `p 1 · q 1`. Two
    families of class name, written twice, is how that happened; one function
    is the fix.

    Measured 2026-08-20 over the 915 published rows: 55 carry an adjusted p
    below 0.05 and 45 of those are plottable, all from PMID:40127276.

    A synonymous row that *does* survive is drawn muted rather than in the
    result colour, because it is the study's own negative control and a
    surviving one is a warning about the comparison rather than a finding
    (`methods.html#rare-variants` says the rest). No published row was in that
    state when measured on 2026-10-26 -- every mirrored synonymous row then came
    from PMID:42230622, which publishes no correction -- so the branch is
    exercised by its test rather than by the corpus. Hollowness is what may not
    be traded away here: it carries the correction, which the muted colour does
    not.
    """
    if not _survived_correction(row):
        return f"chart-{glyph}-open"
    return f"chart-control-{glyph}" if row.consequence_class == "synonymous" else f"chart-{glyph}"


def _statistic(row: BurdenRow) -> str:
    """The one statistic that speaks for a row: the published correction where
    there is one, and the raw p otherwise.

    `q` only where the correction really is a false-discovery rate. The mirror
    carries two methods -- `benjamini_hochberg`, which is an FDR, and
    `familywise_permutation`, which is not -- and labelling the second `q`
    contradicted the same page's own table, which names it "family-wise" a
    screen below. Measured 2026-08-05: 29 matrix cells rendered `q`, some of
    them over a family-wise p.

    One function for the matrix cell and the forest panel, so the two pictures
    on one page cannot name the same number two different ways.
    """
    if row.pvalue_adjusted is None:
        return f"p {row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH
    if row.pvalue_adjustment == "benjamini_hochberg":
        return f"q {row.pvalue_adjusted:.3g}"
    return f"corrected p {row.pvalue_adjusted:.3g}"


def _forest_statistic(row: BurdenRow) -> str:
    """Both statistics for a plotted row: the raw p, and the correction if any.

    The matrix cell shows one because it has one line; a panel row has a whole
    column, and the pair is what the fill of its marker is derived from. A
    reader who sees a hollow marker beside `p 0.045 · q 0.29` can check the
    encoding against the numbers, which is the only way a legend is worth
    anything.
    """
    raw = f"p {row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH
    if row.pvalue_adjusted is None:
        return raw
    return f"{raw} · {_statistic(row)}"


def _effect_bounds(rows: Sequence[BurdenRow]) -> tuple[float, float]:
    """The axis one panel's rows are placed on, padded clear of both edges.

    **1.0 is always among the values**, so the null line is always on the
    axis: a reader must always be able to see which side of "no enrichment" a
    bar falls on, and a panel whose rows are all above 1 would otherwise draw
    no null at all and read as if every row were enriched relative to the
    weakest of them.

    Only positive values are collected. A zero effect and a zero lower bound
    have no position on a logarithmic axis -- `LogScale` refuses them rather
    than clamping, because a clamped zero reads as "the smallest number here"
    -- and `_forest_estimate` draws them as an arrow off the left edge
    instead. Measured 2026-08-20: 182 of the 915 published rows have an effect
    of exactly zero, every one of them with no case carrier, and 152 of those
    also publish a zero lower bound.

    `_axis_bounds` is what makes a panel whose values are all equal drawable at
    all; the padding is applied after it, so the flat case is widened and then
    padded rather than either alone.
    """
    values = [1.0]
    for row in rows:
        values.extend(
            value
            for value in (row.effect, row.ci_low, row.ci_high)
            if value is not None and value > 0
        )
    low, high = _axis_bounds(values)
    return low / _EFFECT_AXIS_PAD, high * _EFFECT_AXIS_PAD


def _forest_line(x1: float, x2: float, y: float, css_class: str) -> str:
    """One horizontal rule at `y`, from `x1` to `x2`."""
    return (
        f'<line class="{html.escape(css_class)}" x1="{coordinate(x1)}" y1="{coordinate(y)}" '
        f'x2="{coordinate(x2)}" y2="{coordinate(y)}"/>'
    )


def _forest_arrow(x: float, y: float, *, pointing: float, css_class: str) -> str:
    """A triangle at an axis edge: this row continues past the end of the axis.

    `pointing` is +1 for the right-hand edge and -1 for the left. The apex sits
    exactly on the edge, which is why `_EFFECT_AXIS_PAD` keeps every measured
    value away from it.

    `css_class` is required rather than defaulted, so that adding a third
    arrow branch cannot quietly inherit the fill of the first two -- which is
    exactly how every arrow on the site came to be drawn as a surviving
    result. `_mark_class` is what callers pass.
    """
    base = x - pointing * _FOREST_ARROW
    return (
        f'<polygon class="{html.escape(css_class)}" points="{coordinate(x)},{coordinate(y)} '
        f"{coordinate(base)},{coordinate(y - _FOREST_ARROW_HALF)} "
        f'{coordinate(base)},{coordinate(y + _FOREST_ARROW_HALF)}"/>'
    )


# SVG's own closed vocabulary for `text-anchor`. Named rather than
# interpolated freely because this is the one *attribute value* the chart code
# builds from a parameter, and `render.data_table` already records why that is
# the dangerous shape: `html.escape` rewrites the quote characters and leaves
# the space and the `=` alone, so a value reaching an attribute outside its
# quotes opens a second attribute. That primitive *refuses* an attribute name
# rather than escaping it, for exactly this reason; this is the same rule one
# level down.
_TEXT_ANCHORS: Final = frozenset({"start", "middle", "end"})


def _forest_text(x: float, y: float, content: Markup, *, anchor: str) -> str:
    """One label.

    `content` is `Markup` rather than `str`, so the obligation is visible at
    the call site instead of living in this sentence: a row label carries a
    `<tspan>` for its union tag, and escaping here would publish the tag as
    literal text. `render.py` expresses the identical contract the identical
    way -- taking a bare `str` here was an asymmetry a review named, not a
    considered difference.

    `anchor` is checked against SVG's own closed set. Every call site passes a
    literal, so this is a guard on a bypassed gate in `data_table`'s idiom:
    reaching it means a caller derived an attribute value from data, and that
    must fail rather than publish. `raise`, never `assert` -- `-O` strips
    `assert`.
    """
    if anchor not in _TEXT_ANCHORS:
        raise ValueError(
            f"text-anchor {anchor!r} is not one of {sorted(_TEXT_ANCHORS)}; escaping "
            f"cannot make an attribute value safe, because a space in one opens a "
            f"second attribute"
        )
    return (
        f'<text class="chart-label" x="{coordinate(x)}" '
        f'y="{coordinate(y + _FOREST_TEXT_DROP)}" text-anchor="{anchor}">{content.html}</text>'
    )


def _forest_interval(row: BurdenRow, y: float, scale: LogScale) -> str:
    """The row's published interval, or the statement that there is none.

    A row with neither bound gets a hatched bar spanning the whole axis, not a
    short one and not nothing. Measured 2026-08-20 over the 915 published rows:
    211 publish neither bound, and every one of the 55 rows carrying an
    adjusted p below 0.05 is among them. A short bar would invent a precision
    the study never claimed, and no bar at all would make the strongest rows on
    the site the faintest marks on it.

    Only 150 of those 211 reach a panel at all -- the other 61 are
    PMID:34324492's, which publish no effect measure either, so `_plottable`
    refuses them before this is called. The two figures are not the same
    quantity and must not be written as one: 10 of the 55 surviving rows are in
    the 61, so "the plottable subset holds every surviving finding" is false
    while "the rows with no interval do" is true.

    A bound that is absent or zero resolves to the axis edge, where
    `_forest_estimate` has already drawn the arrow that says the row runs past
    it. Measured 2026-08-20 on the published corpus, the two coincide exactly:
    every row with a zero lower bound has an effect of zero, and every row with
    no upper bound is flagged `unbounded_above`.
    """
    if row.ci_low is None and row.ci_high is None:
        return _forest_line(_FOREST_PLOT_LEFT, _FOREST_PLOT_RIGHT, y, "chart-nointerval")
    low = scale.x(row.ci_low) if row.ci_low is not None and row.ci_low > 0 else _FOREST_PLOT_LEFT
    high = (
        scale.x(row.ci_high) if row.ci_high is not None and row.ci_high > 0 else _FOREST_PLOT_RIGHT
    )
    control = row.consequence_class == "synonymous"
    return _forest_line(low, high, y, "chart-control" if control else "chart-line")


def _forest_estimate(row: BurdenRow, y: float, scale: LogScale) -> str:
    """The point estimate, or an arrow where the axis has no position for it.

    **Every branch takes its fill from `_mark_class`**, arrows included. The
    arrow branches carried a hardcoded `chart-arrow` until 2026-08-20 and it
    resolved to the same paint as a survived-its-correction circle; see
    `_mark_class` for what that published.

    Two rows have no point to draw and both are drawn as arrows rather than as
    blanks:

    * `effect_bound: unbounded_above` -- Fisher returns an infinite odds ratio
      where no control carries, and `allow_nan=False` refuses to publish it.
      All 23 such published rows have zero control carriers. The lower bound is
      the finding, and it is the left end of the bar this arrow terminates.
    * an effect of exactly zero -- no case carried, on all 182 published rows
      in that state. Zero has no position on a logarithmic axis, and clamping
      it onto the axis floor would read as "the smallest effect measured here"
      rather than as "none".
    """
    if row.effect is None:
        return _forest_arrow(
            _FOREST_PLOT_RIGHT, y, pointing=1.0, css_class=_mark_class(row, "arrow")
        )
    if row.effect <= 0:
        return _forest_arrow(
            _FOREST_PLOT_LEFT, y, pointing=-1.0, css_class=_mark_class(row, "arrow")
        )
    return marker(scale.x(row.effect), y, css_class=_mark_class(row, "point"))


def _forest_row(row: BurdenRow, index: int, scale: LogScale, *, union: bool) -> str:
    """One burden row: its name, its interval, its estimate and its statistics."""
    y = _FOREST_TOP + index * _FOREST_ROW
    label = _forest_label(row.cohort_stratum, row.consequence_class)
    # A `<tspan>` rather than a second `<text>` at a computed x: the label's
    # width depends on the reader's font, which this build cannot measure, and
    # a tag positioned by arithmetic would land on top of the longest labels.
    tag = '<tspan class="chart-union"> union</tspan>' if union else ""
    return (
        _forest_text(_FOREST_LABEL_X, y, Markup(html.escape(label) + tag), anchor="start")
        + _forest_interval(row, y, scale)
        + _forest_estimate(row, y, scale)
        + _forest_text(
            _FOREST_STAT_X, y, Markup(html.escape(_forest_statistic(row))), anchor="start"
        )
    )


def _forest_caption(
    rows: Sequence[BurdenRow],
    measure: str,
    bounds: tuple[float, float],
    kinds: tuple[bool, bool],
    published: Collection[tuple[str, str]],
) -> str:
    """The key. Every clause names a mark that is actually on this panel.

    **One visible line, and the rest behind a toggle (D57, 2026-10-26).** The
    key was one run of up to twelve sentences under every panel -- 185 panels
    on 114 pages, measured on a build of 9f24190. The axis clause stays in
    sight: the measure, its range, and what 1 means for that measure. Every
    other clause folds under "How to read this chart", which closes on a link
    to `methods.html#rare-variants`, where the key is stated once in general.
    **Not every clause is the same wherever it renders.** The mark clauses
    are: each says what a mark means in words that do not change with the
    panel. The not-fully-split clause is not -- it is a warning about this
    panel, that a tagged row's components are not all drawn on it. It folded
    with the rest as the approved D57 spec placed it, on 10 panels on 10 pages
    (TBX5, TBX20, ZEB2, CHD7, ZFPM1, GATA6, FLNA, FLT4, KMT2D, NODAL, measured
    2026-10-26), and **the owner moved it to the visible line the same day
    (E1)**, in a short form that names the `union` tag itself, since the
    clause explaining the tag stays folded. The owner's draft said "some rows
    it contains are not plotted; the table below has every row": on every one
    of those 10 panels exactly one component is missing, and it is a cell the
    study never published, so it is in no table either -- the clause says "at
    least one" and "every row this study published for this gene". **Neither
    move changed a clause's condition**: each still renders only where it is
    true. The muted clause alone was reworded, to name the
    row's role ("its negative control") instead of linking the section the
    toggle already links -- and again the same day, to "one of this study's
    synonymous rows", because each of the 113 panels that draws one draws
    three, one per case group. The toggle's link reads "How to read
    rare-variant results": it read "How to read this atlas", the site nav's
    text for the top of `methods.html`, while it led to one section of it.

    Conditional clause by clause, for the reason `_composite_note` is
    conditional (and `_SYNONYMOUS_NOTICE` was) and `_POOLING_NOTICE` had to be
    reworded when it moved (CLAUDE.md section 4.27): a sentence's truth
    conditions travel with its position, and a legend entry for a glyph the
    reader cannot see sends them hunting for it.

    **The hollow clause was the one exception, in the function whose docstring
    says there are none.** Measured 2026-08-20 on the built corpus: 11 of the
    137 panels printed it with no hollow marker drawn, across 10 genes --
    CHD7, GATA4, KMT2D, MYH7, NKX2-5, NODAL, RBFOX2, RNF40, WDR5, ZIC3, six
    of them ClinGen `definitive` -- and on 9 of those the panel did draw
    arrows, so the key told the reader that everything visible had failed a
    correction. Every other clause measured 0 captioned-but-not-drawn.

    **And there was no clause for a filled mark at all**, so its meaning had
    to be inferred as the complement of the hollow one. That inference was
    wrong for all 205 arrows on the site until `_mark_class` began deciding
    an arrow's fill; it is stated now rather than left to be worked out.
    The filled clause excludes a surviving *synonymous* row, which is drawn
    muted rather than in the result colour and has its own clause below.

    **The axis clause was wrong in a different way: it named a mark that is
    on the panel and then said the wrong thing about it.** The measure label
    was already dynamic; the gloss beside it was the literal "no enrichment"
    for every measure, so 112 of the 137 panels told the reader that 1 on an
    odds-ratio axis means no enrichment. Both halves now resolve from the same
    token -- see `_NULL_LINE_GLOSS` for the measurement and the third
    vocabulary. The sibling `<title>` on the same figure had always varied
    correctly, which is why reading the panel did not show it.

    **The survived-but-hollow clause (2026-10-26)** renders only where a
    non-synonymous row cleared its study's correction without supporting the
    gene (`_survived_against_the_gene`), which `_mark_class` draws hollow --
    no published panel today. The hollow clause's "did not survive" is false
    of that row, so the hollow clause renders only for the other hollow rows
    and the new clause stands on its own where that row is the panel's only
    one.

    **The union clause is conditional twice**: once on a union being shown as
    one at all, and once more on whether every tagged row is fully split.
    Measured 2026-08-20: 10 panels tag a row with only one of its two
    components drawn, TBX5's Audain odds-ratio panel among them, and the key
    told the reader to reconcile it against a row that is not there.

    The pooling clause is **not** conditional, and it is worded as policy
    rather than as an observation about this page. It is here at all because
    the forest idiom ends in a summary diamond and a reader who knows the
    idiom will look for one; its absence has to be an answer rather than a
    gap. It deliberately does not reuse `_POOLING_NOTICE`'s own sentence:
    that constant is the matrix's caption, and a test asserts it survives
    outside the fold, which a second copy of its wording here would satisfy on
    its behalf.
    """
    low, high = bounds
    label = _MEASURE_LABEL.get(measure, measure)
    # Both halves resolve from the same token. The gloss was hardcoded beside
    # this already-dynamic label, which is how 112 odds-ratio panels came to
    # call 1 "no enrichment" -- see `_NULL_LINE_GLOSS`.
    gloss = _NULL_LINE_GLOSS.get(measure, _NULL_LINE_DEFAULT)
    parts = [
        f"Axis: <strong>{html.escape(label)}</strong> on a logarithmic scale, "
        f"{_fmt(low)} to {_fmt(high)}; the vertical line is <strong>1</strong>, "
        f"{html.escape(gloss)}. ",
    ]
    if any(_survived_correction(row) and row.consequence_class != "synonymous" for row in rows):
        parts.append(
            "A <strong>filled</strong> marker is a row that <strong>survived</strong> its own "
            "study's correction for multiple testing. "
        )
    if any(not _survived_correction(row) and not _survived_against_the_gene(row) for row in rows):
        parts.append(
            "A <strong>hollow</strong> marker is a row that <strong>did not survive</strong> its "
            "own study's correction for multiple testing, or whose study published no correction "
            "at all &mdash; the atlas computes none of its own. "
        )
    # The one hollow row the clause above is false of; see
    # `_survived_against_the_gene`. Self-contained, because on a panel where it
    # is the only hollow row the clause above does not render.
    if any(_survived_against_the_gene(row) for row in rows):
        parts.append(
            "A row that <strong>survived</strong> its own study's correction without pointing "
            "to enrichment in cases is drawn <strong>hollow</strong>: it is not evidence for "
            "this gene. "
        )
    parts.append(
        "A forest plot usually ends in a summary diamond. There is <strong>no pooled</strong> "
        "estimate here and none is computed anywhere on this site: two studies that share a "
        "sample collection describe partly the same people, and adding their rows would "
        "count them more than once. "
    )
    if any(row.ci_low is None and row.ci_high is None for row in rows):
        parts.append(
            "A <strong>hatched</strong> bar means <strong>no interval published</strong> for "
            "that row, not a wide one. "
        )
    if any(row.effect is None for row in rows):
        parts.append(
            "An arrow to the <strong>right</strong> is an effect unbounded above &mdash; no "
            "control carried the variant &mdash; so the study's lower bound is the whole "
            "finding and there is no point estimate to place. "
        )
    if any(row.effect is not None and row.effect <= 0 for row in rows):
        parts.append(
            "An arrow to the <strong>left</strong> is an effect of zero &mdash; no case "
            "carried the variant &mdash; which a logarithmic axis has no position for. "
        )
    consequence, stratum = kinds
    reasons = []
    if consequence:
        reasons.append(
            "<code>damaging</code> is the loss-of-function and damaging-missense rows together"
        )
    if stratum:
        reasons.append("<code>all cases</code> is the syndromic and non-syndromic rows together")
    if reasons:
        parts.append(
            "A row tagged <strong>union</strong> contains at least one other row in this panel "
            "and is not an independent finding: " + "; ".join(reasons) + ". "
        )
    if any(row.consequence_class == "synonymous" for row in rows):
        # What this clause owes the reader is what the *muting* means, which
        # nothing else on the panel says, and the row's role, in two words.
        # Why a synonymous row calibrates the others is
        # `methods.html#rare-variants`'s, which the toggle this clause sits in
        # closes on a link to. It pointed at "<em>How to read these
        # numbers</em>, above" until 2026-10-26, a fold that left the page that
        # day (D50), and then for one revision at that section by a link of its
        # own -- a second link to the same place, one sentence before the
        # toggle's.
        parts.append(
            "A <strong>muted</strong> row is one of this study's <strong>synonymous</strong> "
            "rows, its negative control, drawn apart from the results it calibrates."
        )
    # The axis clause is first and unconditional, and the pooling clause is
    # unconditional too, so `rest` is never empty.
    # Particular to this panel, so on the visible line beside the axis and only
    # where it fires (owner decision 2026-10-26, E1); the union clause it
    # qualifies stays folded, and this names the row and the tag itself so it
    # reads without it. See `_unsplit_line`.
    split = _unsplit_line(rows, kinds, published)
    axis, rest = parts[0], "".join(parts[1:]).strip()
    return (
        f"<figcaption><p>{axis.strip()}{split}</p>"
        '<details class="chart-notes"><summary>How to read this chart</summary>'
        f"<p>{rest}</p>"
        '<p><a href="../methods.html#rare-variants">How to read rare-variant results</a></p>'
        "</details></figcaption>"
    )


def _forest(
    rows: Sequence[BurdenRow],
    measure: str,
    study_label: str,
    published: Collection[tuple[str, str]],
) -> str:
    """One study's rows for one effect measure, as effect and interval.

    **One panel per (study, effect measure), never one per study.** 779 of the
    915 published rows carry an odds ratio and 75 an enrichment ratio, and one
    study publishes both: an odds ratio of 3.1 and a de novo enrichment of 3.1
    are different claims, so placing them on one axis would equate them at the
    only scale a reader actually reads. `_effect` has no branch that omits the
    measure's name for the same reason; this is that rule in geometry.

    Rows are in `_sort_key` order, so the panel reads top to bottom in the
    order the table below it reads -- the composite above its components, the
    synonymous negative control last.

    Raises on an empty panel rather than emitting an axis with nothing on it.
    The caller filters on `_plottable` and never calls this with nothing, so
    this is a guard on a bypassed gate: a bypassed gate must fail rather than
    publish (`raise`, never `assert` -- `-O` strips `assert`).
    """
    if not rows:
        raise ValueError("a forest panel needs at least one row to place")
    ordered = sorted(rows, key=_sort_key)
    bounds = _effect_bounds(ordered)
    scale = LogScale(
        low=bounds[0],
        high=bounds[1],
        left=_FOREST_PLOT_LEFT,
        width=_FOREST_PLOT_RIGHT - _FOREST_PLOT_LEFT,
    )
    kinds = _union_kinds(ordered)
    axis_y = _FOREST_TOP + (len(ordered) - 1) * _FOREST_ROW + _FOREST_AXIS_GAP
    top = _FOREST_TOP - _FOREST_ROW / 2
    null_x = scale.x(1.0)
    body = (
        _forest_line(_FOREST_PLOT_LEFT, _FOREST_PLOT_RIGHT, axis_y, "chart-axis")
        + f'<line class="chart-axis" x1="{coordinate(null_x)}" y1="{coordinate(top)}" '
        f'x2="{coordinate(null_x)}" y2="{coordinate(axis_y)}"/>'
        + _forest_text(null_x, top - 6.0, Markup("1"), anchor="middle")
        + _forest_text(_FOREST_PLOT_LEFT, axis_y + 8.0, Markup(_fmt(bounds[0])), anchor="start")
        + _forest_text(_FOREST_PLOT_RIGHT, axis_y + 8.0, Markup(_fmt(bounds[1])), anchor="end")
        + "".join(
            _forest_row(row, index, scale, union=_is_union(row, kinds))
            for index, row in enumerate(ordered)
        )
    )
    noun = "statistic" if len(ordered) == 1 else "statistics"
    title = (
        f"{study_label}: {len(ordered)} published "
        f"{_MEASURE_LABEL.get(measure, measure)} {noun} for this gene on a logarithmic axis, "
        f"with 95% confidence intervals where this study published them, and no pooled summary."
    )
    figure = svg_figure(
        width=_FOREST_WIDTH,
        height=round(axis_y + _FOREST_FOOT),
        title=title,
        body=body,
    )
    return (
        f'<figure class="forest" data-effect-measure="{html.escape(measure)}">'
        f"{figure}{_forest_caption(ordered, measure, bounds, kinds, published)}</figure>"
    )


def _composite_note(rows: Sequence[BurdenRow]) -> str:
    """Say that the `damaging` rows are the union of the two below them.

    Without it a reader meets three consequence rows and reads three findings,
    then adds up carrier counts that already include each other. Measured on
    CHD7 in PMID:40127276: the damaging de novo row is 20 mutations, and the
    loss-of-function and damaging-missense rows below it are 16 and 4 of *those
    same* 20.

    The composite is not droppable in favour of its parts -- it carries its own
    p-value, which is not a function of theirs, and it is the analysis the study
    defines its results by -- so the relationship has to be stated rather than
    designed away.

    Conditional on the study actually reporting both a composite and a
    component. Rendered unconditionally it would assert a decomposition that
    does not exist for the two studies whose tables carry no `damaging` row at
    all, which is the defect `_POOLING_NOTICE` was made conditional for.

    **One sentence since 2026-10-26 (D50, P8).** Its second -- the composite is
    shown because the study defines its findings by it and its p-value is not
    derivable from the other two -- is general, true of every block that
    carries the note, and moved to `methods.html#rare-variants`. The forest's
    union tag stays in its key, and the row's `method_note` footnote, which is
    mirror data, is unchanged.
    """
    consequences = {row.consequence_class for row in rows}
    if "damaging" not in consequences or not consequences & {"lof", "missense_damaging"}:
        return ""
    return (
        '<p class="notice-inline">The <strong>damaging (LOF + missense)</strong> rows are '
        "the <strong>union</strong> of the loss-of-function and damaging-missense rows "
        "below them, not a third independent result &mdash; the same variants are counted "
        "in both.</p>"
    )


def _combined_line(rows: Sequence[CombinedRow]) -> str:
    """A study's own combination of its arms, stated above the arms' figures.

    PMID:40127276's primary result is neither arm but their combination, and
    until schema 2.14 it reached no page, because it could not go where a
    corrected p goes (BUR018 refuses it there, correctly). Rendered outside every
    fold because it is particular to this gene.

    One paragraph per distinct (arms, method) pair, in sorted order; every study
    today has at most one. Conditional, like `_composite_note`: a study with no
    combination renders nothing, because an absent combination is not a claim.
    """
    paragraphs = []
    for combines, method in sorted({(row.combines, row.method) for row in rows}):
        group = [row for row in rows if (row.combines, row.method) == (combines, method)]
        arms = " + ".join(_COMPARATOR_LABEL.get(c, c) for c in combines)
        name, statistic, caveat = _COMBINATION_LABEL.get(method, (method, method, ""))
        results = " &middot; ".join(
            (
                f"{html.escape(_STRATUM_LABEL.get(row.cohort_stratum, row.cohort_stratum))}, "
                if row.cohort_stratum != "all"
                else ""
            )
            + html.escape(_CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class))
            + f" p {row.pvalue:.3g}, {html.escape(statistic)} {row.local_fdr:.3g}"
            for row in group
        )
        paragraphs.append(
            '<p class="notice-inline"><strong>'
            + html.escape("Combined across this study's arms")
            + f"</strong> ({html.escape(arms)}, by {html.escape(name)}): {results}. "
            + html.escape(
                "This is the authors' own meta-analysis of the rows below, within one cohort"
            )
            + " &mdash; not a statistic this atlas computed, and not a pooling across "
            + f"studies. {caveat}</p>"
        )
    return "".join(paragraphs)


def _footnotes(rows: Sequence[BurdenRow], combined: Sequence[CombinedRow] = ()) -> str:
    """`method_note`, rendered -- which `tables.py` said it was and it was not.

    The column's comment reads "rendered verbatim as a row footnote", and no
    code path read it. It exists for the CNV case, where "carrier" means
    different things in different papers (any overlap, exonic, whole gene), so
    the first curator to need it would have got a green build, a correct bundle
    and no page text.

    A combined row's note joins the arms', so the method it names is footnoted
    under the figures it combines.

    A note that is empty once stripped is no note. One of whitespace alone
    rendered an empty `<li>` until 2026-10-26 -- 0 such rows in the mirror, so
    latent -- while the card's `caveats` label already read it as none.
    """
    raw = [row.method_note for row in rows] + [row.method_note for row in combined]
    notes = sorted({note.strip() for note in raw if note and note.strip()})
    if not notes:
        return ""
    items = "".join(f"<li>{html.escape(note)}</li>" for note in notes)
    return f'<ul class="footnotes">{items}</ul>'


def _cohort_notes(rows: Sequence[BurdenRow], cohorts: Mapping[str, Cohort]) -> str:
    """What each collection is, and the caveats its curator recorded.

    **Titled "About these cohorts, and the caveats recorded for them"**, which
    read "..., and what qualifies every number above" until 2026-10-26 on all
    118 pages carrying it. Several descriptions qualify nothing -- `cnchd`'s,
    `nottingham`'s and `isca`'s say what the collection is and where a study
    names it -- so the title claimed of every collection what only some carry;
    `methods.html#cohorts` made the same correction that day.

    **These reached no published byte.** `curation/cohorts.yaml` records, for
    UK Biobank, that its participants are adults recruited at 40-69 while a
    substantial proportion of the CHD cases were enrolled in childhood -- a
    survivorship bias that inflates every odds ratio on every page carrying a
    UK Biobank comparison -- and that
    the combined cohort is ~92% European. `models/cohort.py` says in as many
    words that `description` is where such caveats belong. `cohort_labels`
    returned `{id: name}` and dropped them, so the curator wrote the caveat, the
    model documented it, and no reader could reach it.

    In a `<details>` element: the caveats are long, they repeat on every page
    carrying a burden table, and
    a reader who has read them once should not have to scroll past them again.
    `<details>` is native HTML, needs no script, and its contents are in the
    document for a crawler, a `curl` and a reader with JavaScript disabled.

    **The fold stays, and ends on a link to `methods.html#cohorts`** (P4,
    2026-10-26), which lists every cohort in the registry. Only this gene's
    cohorts are here; the full list is a second place a caveat reaches, never
    a replacement for this one (CLAUDE.md section 9).

    **A record's open question follows its description here too** (owner's
    standing rule, 2026-10-27): its status's label and the curator's note,
    escaped, by `render.review_status_note`, the function the methods page's
    list uses. It reaches the same pages the description does, and a
    collection with no question recorded gets nothing.
    """
    cited = sorted(
        {identifier for row in rows for identifier in (*row.case_cohorts, *row.control_cohorts)}
    )
    items = [
        f"<li><strong>{html.escape(cohort.name)}</strong> &mdash; "
        f"{html.escape(cohort.description)}"
        f"{review_status_note(cohort.id, cohort.review_status, cohort.review_note)}</li>"
        for identifier in cited
        if (cohort := cohorts.get(identifier)) is not None
    ]
    if not items:
        return ""
    return (
        '<details class="cohort-notes"><summary>About these cohorts, and the caveats '
        "recorded for them</summary><ul>"
        + "".join(items)
        + '</ul><p><a href="../methods.html#cohorts">Every cohort in this atlas</a></p></details>'
    )


def _overlap_notice(
    left: str,
    right: str,
    overlap: Overlap,
    publications: Mapping[str, Publication],
    cohorts: Mapping[str, Cohort],
) -> str:
    """Why two studies on one page are not independent, the way the page names a shared one.

    A shared collection reads as it always has: both "draw on" it, so their
    results "describe partly the same people". A declared overlap (owner
    decision 2026-10-27, I3; `Cohort.overlaps_with`) is two collections, one
    per study, so it says which each draws on, that this atlas records them as
    overlapping, and "may describe", because the registry records a likely
    overlap, not a counted one: PMID:34324492's author answered that its case
    set does or likely does include DDD participants. Measured on a real build
    the day it was written: 72 pages carry both PMID:34324492 and PMID:42230622,
    and each now carries this notice once.
    """
    first = html.escape(_study_label(left, publications))
    second = html.escape(_study_label(right, publications))
    clauses = []
    if overlap.shared:
        clauses.append(
            f"{first} and {second} both draw on {html.escape(_names(overlap.shared, cohorts))}"
        )
    if overlap.declared:
        mine = sorted({collection for collection, _ in overlap.declared})
        theirs = sorted({collection for _, collection in overlap.declared})
        clauses.append(
            f"{first} draws on {html.escape(_names(mine, cohorts))} and {second} on "
            f"{html.escape(_names(theirs, cohorts))}, collections this atlas records as "
            "overlapping"
        )
    verb = "describe" if overlap.shared else "may describe"
    return (
        '<p class="notice-inline"><strong>These two studies are not independent.</strong> '
        f"{', and '.join(clauses)}, so their results {verb} partly the same people and "
        "must not be combined.</p>"
    )


def _names(ids: Sequence[str], cohorts: Mapping[str, Cohort]) -> str:
    """Cohort ids rendered as the collections they name, comma-separated.

    An id with no record renders as the id, which a reader can still look up;
    BUR009 reports the absence and the gate refuses, so this is a guard on a
    bypassed gate.
    """
    record = {identifier: cohorts.get(identifier) for identifier in ids}
    return ", ".join(
        value.name if value is not None else identifier for identifier, value in record.items()
    )


# --- Developmental expression -----------------------------------------------
#
# `ExpressionProfile` (`build/profiles.py`) already carries every figure this
# section renders -- a percentile band, tau, a phase placement -- with every
# qualifier D39 requires attached to it in the JSON. Nothing below computes
# anything; it only chooses the English for a payload whose own module
# docstring says exactly that is a later task's job.
#
# D43 governs this section: **never a panel redundant with the source's own
# browser**. That is D32's actual test -- "if a panel would look the same on
# the Heart Cell Atlas, it does not belong here" -- and it replaces the proxy
# this comment carried until 2026-08-20, which read "sentences and small
# per-organ tables, never a plot".
#
# The proxy was wrong in both directions. It forbade the phase-banded
# trajectory, which no external browser has, while permitting the rendering
# that hid a live defect for the whole life of schema 2.11: stages published
# in alphabetical order, 4 wpc eighth, because 21 individually-correct blocks
# cannot show a trajectory and therefore cannot show a scrambled one. (This
# read "for three releases" until 2026-08-21; the profiles mirror landed on
# 2026-08-19 and nothing published a stage array before it, so 2.11 is the
# only version that ever carried the defect. See `manifest.py`'s 2.12 note.)
#
# A bare heart curve IS what the source's own browser shows and remains
# forbidden. The same curve banded by this atlas's own cardiac-phase
# vocabulary is not. `test_every_chart_carries_an_atlas_specific_axis` is what
# watches for the difference.

# What each `ProfileGap` means, in English. Split into three maps rather than
# one, because `ProfileGap`'s own docstring says only the first two members
# mean exactly the same fact wherever they appear -- `dataset_not_registered`
# and `detection_floor_undeclared` are dataset-level gaps with one meaning in
# both fields that publish them, while `no_quantile_grid` is placement-only and
# `one_organ_sampled`/`undefined` are tau-only. Sharing only the two that
# genuinely share a meaning is `EVIDENCE_STATE_LABELS`'s discipline applied
# here: one vocabulary where the fact is one fact, never a false shared
# spelling for two different ones.
_DATASET_GAP_CLAUSE: Final[dict[str, str]] = {
    ProfileGap.DATASET_NOT_REGISTERED.value: (
        "this dataset is not registered in the atlas's own records"
    ),
    ProfileGap.FLOOR_UNDECLARED.value: "this dataset has not declared a detection floor",
}

_PLACEMENT_GAP_CLAUSE: Final[dict[str, str]] = {
    **_DATASET_GAP_CLAUSE,
    ProfileGap.NO_QUANTILE_GRID.value: (
        "no complete percentile grid is published for this organ at this stage"
    ),
    ProfileGap.NOT_ON_A_LOG_AXIS.value: (
        "this atlas placed the measurement and a logarithmic axis has no "
        "position for it: the median is zero or below"
    ),
}

_SPECIFICITY_GAP_CLAUSE: Final[dict[str, str]] = {
    **_DATASET_GAP_CLAUSE,
    ProfileGap.ONE_ORGAN_SAMPLED.value: (
        "only one organ was sampled at this stage, and τ needs at least two"
    ),
    ProfileGap.UNDEFINED.value: "it could not be computed for this gene at this stage",
}

# The shape `bundles._expression_profile` publishes for a gene no curated bulk
# RNA-seq dataset mentions at all (`EMPTY_EXPRESSION_PROFILE`, imported rather
# than reconstructed here, so the two cannot drift on what "empty" means).
#
# **Chosen over "this dataset does not cover it", and that is a different
# claim.** `datasets == []` means no `profiles.tsv` row anywhere names this
# gene: there is no specific dataset in view to say did not cover it, and
# saying so would invent a check this atlas never ran. This is the idiom of the
# lesion-assertion notice (`_not_curated` until D55 moved it into the card as
# `_card_curation`), not `_burden_section`'s: an absent section here would read as
# "the atlas looked and found nothing", indistinguishable from what this
# sentence actually means, so the section renders regardless.
_NOT_CURATED_EXPRESSION: Final = (
    "<p>The atlas has <strong>not yet curated</strong> a developmental expression "
    "profile for this gene: no curated bulk RNA-seq dataset mentions it today. That is "
    "a different claim from a specific dataset having looked at this gene and found "
    "nothing.</p>"
)

# Unconditional -- rendered beside a page of nothing but detected, well-placed
# figures exactly as beside a gene with every organ below the floor. The
# precedent for both the wording and the placement decision is
# `_POOLING_NOTICE`: made conditional once, it read false in the position that
# motivated the fix, because a caveat that appears only next to the bad news
# reads as an excuse for it. This one is a general fact about bulk RNA-seq and
# CHD gene biology, true whatever this gene's own numbers say, so it carries no
# clause that could be false in either position.
#
# This is the burden layer's "a missing cell is not a null result" discipline,
# needed more here: that layer leaves a gap as an absence, while this one turns
# a low number into a sentence, which is the more persuasive -- and more
# dangerous -- thing to get wrong.
#
# The spec's short form since 2026-10-26 (D50, spec section 5), one of the ten
# caveats it keeps visible: the claim a reader acts on -- bulk tissue, and a low
# or absent figure is not evidence the gene is unimportant to heart development
# -- and a link. The dilution argument, which lineages a heart sample averages
# over, is on `methods.html#expression` in full (`methods._EXPRESSION`). The
# short form drops "never a single cell type", which "bulk tissue" already
# says, and "bulk" before "figure", which the bold clause before it already
# scopes. "below" stays true in this position: the notice opens the section,
# above every dataset block.
#
# "whole, bulk tissue" until the owner's decision of 2026-10-26 (E2): the
# dataset's own curation record sets the Expression Atlas's organ tokens
# (forebrain, hindbrain) against the paper's (brain, cerebellum), so "whole"
# need not hold of every organ and stage, and the dilution argument rests on
# bulk -- mixed cell types -- alone.
#
# Named `_BULK_DILUTION_NOTICE` until review (2026-10-26): the one line it keeps
# says nothing of dilution, so the name described text that had left it.
_BULK_TISSUE_NOTICE: Final = (
    '<p class="notice-inline"><strong>Every measurement below is from bulk '
    "tissue</strong>: a low or absent figure is not evidence a gene is unimportant to "
    'heart development. <a href="../methods.html#expression">How to read these '
    "figures</a>.</p>"
)

# The banner of a gene that has cell-type data, in place of `_BULK_TISSUE_NOTICE`
# (never beside it: that notice says every measurement below is bulk, which is
# false once a cell-type block opens the section). It keeps the same caveat,
# unconditional, in one sentence that also says what the two layers are.
_CELL_TYPE_NOTICE: Final = (
    '<p class="notice-inline">Cell-type values pool 8 hearts at 9 to 15 weeks and the '
    "timeline below is from bulk tissue; a low or absent value is not evidence a gene is "
    "unimportant to heart development.</p>"
)

# The same, for a gene with cell-type data and no bulk profile: there is no
# timeline below, so the banner must not claim one.
_CELL_TYPE_ONLY_NOTICE: Final = (
    '<p class="notice-inline">Cell-type values pool 8 hearts at 9 to 15 weeks; a low or '
    "absent value is not evidence a gene is unimportant to heart development.</p>"
)

# `_EXPRESSION_READING_NOTES` and `_PERCENTILE_COMPARABILITY_NOTICE` sat here
# until 2026-10-26: the reading notes every gene page folded under "How to read
# these figures". They moved, with their histories, to `methods._EXPRESSION`
# (D50).


def _abundance(tissue: TissueProfileEntry) -> str:
    """ "100 rpkm (n=3 samples)" -- never the number alone.

    `n_samples` travels with every figure on this page, the `_count`/
    `count_unit` discipline applied here: the schema permits `n_samples` = 1,
    so "median" can be one observation, and a reader must be told rather than
    left to assume a study-sized sample sits behind every row.
    """
    noun = "sample" if tissue["n_samples"] == 1 else "samples"
    unit = tissue["unit"]
    return f"{_fmt(tissue['median_abundance'])} {unit} (n={tissue['n_samples']} {noun})"


def _phases_cell(phase: PhaseInfo) -> str:
    """Where a stage falls in the curated phase vocabulary, or why not -- never blank.

    - **Every id, never the first.** Human cardiac morphogenesis runs several
      processes at once (at 6 wpc an AVSD gene needs atrial septation and a
      TGA gene outflow-tract septation, both true together), so a MATCHED stage
      names every phase id, joined by ", " with "_" read as a space and no other
      rewrite: an id is this atlas's own slug, not mirrored third-party text.
    - **An empty MATCHED stage reads "unnamed phase"**, a named state rather
      than a blank that reads as a rendering bug.
    - **Any other outcome renders `reason` verbatim** (`profiles.py` wrote it as
      a short reader-facing clause), or "not available"; a second vocabulary
      here could only repeat or drift from it. `data_table` escapes the cell.
    """
    if phase["outcome"] == PhaseOutcome.MATCHED.value:
        ids = phase["phase_ids"]
        if not ids:
            return "unnamed phase"
        return ", ".join(phase_id.replace("_", " ") or "unnamed phase" for phase_id in ids)
    return phase["reason"] or "not available"


def _entry_unit(entry: DatasetProfileEntry) -> str:
    """The first unit this entry's measurements carry, or "".

    One dataset quotes every organ in one unit because validation makes it:
    PRF008 is an ERROR for a dataset whose `median_abundance` comes in more than
    one unit, and the build refuses a corpus with any error. It names the unit
    in a header and in the floor clauses, where no single cell supplies it.
    """
    return next(
        (cell["unit"] for stage in entry["stages"] for cell in stage["tissues"] if cell["unit"]),
        "",
    )


def _floor_words(floor: float | None, unit: str) -> str:
    """ "1 tpm" for the floor value with its unit, "" where the floor is unknown."""
    if floor is None:
        return ""
    return " ".join(part for part in (_fmt(floor), unit) if part)


def _tissue_percentile(tissue: TissueProfileEntry, floor: float | None, *, with_grid: bool) -> str:
    """The percentile, or the named reason there is none -- never "not detected".

    **`below_detection_floor` has its own required wording (rule 5).** "Not
    detected" is a negative assertion a clinical reader takes as evidence
    against a gene, and this atlas measured no such thing: a value below a bulk
    assay's floor is one it does not vouch for. The floor's value and unit sit
    beside it, because a reader cannot judge how far below without them. Every
    other gap resolves through `_PLACEMENT_GAP_CLAUSE`, the one place that says
    why a placement is missing.
    """
    placed = tissue["placement"]
    if placed is not None:
        text = f"{placed['median_percentile']}"
        if placed["q25_percentile"] is not None and placed["q75_percentile"] is not None:
            text += f" (IQR {placed['q25_percentile']}–{placed['q75_percentile']})"
        if with_grid:
            text += f" of {placed['n_genes']:,} genes"
        return text
    reason = tissue["not_placed_reason"]
    if reason == ProfileGap.BELOW_DETECTION_FLOOR.value:
        words = _floor_words(floor, tissue["unit"])
        return f"below the {words} floor" if words else "below the detection floor"
    return _PLACEMENT_GAP_CLAUSE.get(
        reason or "", "no percentile is available for this measurement"
    )


def _tau_cell(stage: StageProfileEntry, floor: float | None, unit: str) -> str:
    """tau, or "not available: " and the reason it is absent -- never blank.

    **No adjective and no band on the number** (D39(c) reserves that judgement
    from the atlas): the scale is in the column header and the organs are in
    the highest-median column, and the reader weighs one against the other.
    `peak_below_detection_floor` is kept out of `_SPECIFICITY_GAP_CLAUSE`
    because it is the one reason with a number to show -- the floor, with its
    unit, as rule 5 owes a reader told a figure is missing because of one.
    """
    spec = stage["specificity"]
    if spec is not None:
        return _fmt(spec["tau"])
    reason = stage["specificity_unavailable_reason"]
    if reason == ProfileGap.PEAK_BELOW_DETECTION_FLOOR.value:
        words = _floor_words(floor, unit)
        qualifier = f" {words}" if words else " detection"
        return f"not available: every organ sampled is below the{qualifier} floor"
    clause = _SPECIFICITY_GAP_CLAUSE.get(reason or "", "no reason was recorded for this gap")
    return f"not available: {clause}"


def _highest_cell(stage: StageProfileEntry) -> str:
    """The organ with the highest median, "tied", or "-" where there is no tau.

    tau measures concentration, not location, so an organ-naming cell keyed on
    tau alone can state the opposite of the truth for a gene concentrated
    elsewhere: the organ is `highest_in`, the argmax. A tie (`highest_in` null)
    reads "tied" rather than naming an arbitrary winner.
    """
    spec = stage["specificity"]
    if spec is None:
        return "–"
    return spec["highest_in"] or "tied"


def _other_organ_cell(
    cell: TissueProfileEntry | None, floor: float | None, *, with_grid: bool
) -> str:
    """One other-organ cell: the median with its n, then the percentile or the reason."""
    if cell is None:
        return "–"
    figure = _abundance(cell)
    return f"{figure}; {_tissue_percentile(cell, floor, with_grid=with_grid)}"


def _stage_tables(entry: DatasetProfileEntry, cardiac: frozenset[str], floor: float | None) -> str:
    """Every figure of one dataset: a table per cardiac tissue, one for the other organs.

    This replaced 21 near-identical per-stage blocks per dataset (59% of the
    section's text was verbatim repetition, measured on TBX5 by the method in
    `_dataset_block`'s docstring) without taking a figure out of the
    HTML: a chart that replaced them would be the characteristic defect, curated
    work reaching no page, in a redesign. Decisions the blocks made, kept here:

    - **The caption restates rule 8.** Rows laid out across stages invite
      comparing percentiles that come from different grids, so it says a
      percentile places this gene among N genes in the same organ and stage.
      Where `n_genes` differs between cells it says so and every percentile cell
      carries its own " of N genes".
    - **The heart table has a row for every stage**, sampled or not, so the
      stage-level facts (phases, tau, highest median) survive for a stage the
      tissue lacks; its median cell reads "not sampled", its percentile "–".
    - **Phases:** every id, never the first (`_phases_cell`). **tau:** no
      adjective; the scale is in the header and a tie reads "tied"
      (`_highest_cell`). **Every gap is a named state**, never blank, and
      rule 5's below-floor wording never reads "not detected".
    - **Stage-level facts never depend on a tissue.** With no cardiac tissue
      sampled, one table of stages carries phases, tau and the highest median.
    - **One tau scale per table.** Its header names the scale, so two stages
      disagreeing raise: a bypassed gate must fail rather than publish a header
      false for some rows.
    """
    stages = entry["stages"]
    unit = _entry_unit(entry)
    scales = {stage["specificity"]["scale"] for stage in stages if stage["specificity"]}
    if len(scales) > 1:
        raise ValueError(f"τ scales differ across stages: {sorted(scales)}")
    tau_header = f"τ ({next(iter(scales))})" if scales else "τ"
    grids = {
        cell["placement"]["n_genes"]
        for stage in stages
        for cell in stage["tissues"]
        if cell["placement"] is not None
    }
    uniform = len(grids) == 1
    with_grid = not uniform
    if uniform:
        caption = (
            f"A percentile places this gene among {next(iter(grids)):,} genes in the same "
            "organ and stage, so it compares nothing across rows or columns."
        )
    else:
        caption = (
            "A percentile compares this gene only with genes in the same organ and stage, "
            "never across rows or columns."
        )
    organs = _organ_order(entry, cardiac)
    heart = [organ for organ in organs if organ in cardiac]
    others = [organ for organ in organs if organ not in cardiac]
    median_header = f"median {unit} (n)" if unit else "median (n)"
    parts = [f'<p class="method">{html.escape(caption)}</p>']
    for organ in heart:
        rows = []
        for stage in stages:
            cell = next((c for c in stage["tissues"] if c["tissue"] == organ), None)
            rows.append(
                Row(
                    cells=(
                        _stage_label(stage),
                        _phases_cell(stage["phase"]),
                        _abundance(cell) if cell is not None else "not sampled",
                        _tissue_percentile(cell, floor, with_grid=with_grid)
                        if cell is not None
                        else "–",
                        _tau_cell(stage, floor, unit),
                        _highest_cell(stage),
                    )
                )
            )
        headers = (
            "stage",
            "phases",
            median_header,
            "percentile (IQR)",
            tau_header,
            "highest median",
        )
        parts.append(
            f'<p class="method"><strong>{html.escape(organ)}</strong></p>'
            f'<div class="scroll">{data_table(headers, rows)}</div>'
        )
    if not heart:
        # No cardiac tissue was sampled (or the dataset declares none), so no
        # heart table carries the stage-level facts -- the phase, tau and the
        # highest organ belong to the stage, not to a tissue, and would reach no
        # byte of the page. One table of stages holds them instead.
        rows = [
            Row(
                cells=(
                    _stage_label(stage),
                    _phases_cell(stage["phase"]),
                    _tau_cell(stage, floor, unit),
                    _highest_cell(stage),
                )
            )
            for stage in stages
        ]
        stage_headers = ("stage", "phases", tau_header, "highest median")
        parts.append(
            '<p class="method"><strong>All stages</strong></p>'
            f'<div class="scroll">{data_table(stage_headers, rows)}</div>'
        )
    if others:
        rows = [
            Row(
                cells=(
                    _stage_label(stage),
                    *(
                        _other_organ_cell(
                            next((c for c in stage["tissues"] if c["tissue"] == organ), None),
                            floor,
                            with_grid=with_grid,
                        )
                        for organ in others
                    ),
                )
            )
            for stage in stages
        ]
        parts.append(
            '<p class="method"><strong>Other organs</strong></p>'
            f'<div class="scroll">{data_table(("stage", *others), rows)}</div>'
        )
    return "".join(parts)


# --- The trajectory chart ----------------------------------------------------

# Below this many placed stages a trajectory is not a trajectory. Two points
# make a line segment, and a line segment is a claim about the interval
# between them that two measurements do not support. Measured over the
# committed corpus 2026-08-20: 82 of 92 published genes clear this in heart,
# 3 fall to markers, and 7 have nothing placed at all.
_STAGES_FOR_A_TRAJECTORY: Final = 3

_CHART_WIDTH: Final = 560
# Wide enough for the grid's organ labels, which sit in the margin left of the
# first column (`_grid`): "forebrain" and "hindbrain", the longest organ names
# in E-MTAB-6814, are 9 characters. The value axis's own labels and title share
# this margin.
_PLOT_LEFT: Final = 64.0
_PLOT_RIGHT: Final = 544.0
# The heart chart's height. Its top and bottom are not constants: the phase
# strip above it has one row per banded phase, so where the plot starts depends
# on how many phases this dataset's stages fall inside (`_bulk_figure`).
_PLOT_HEIGHT: Final = 112.0

# The lowest placed measurement stops this far above the axis, and the gap is
# load-bearing rather than padding. `_bulk_figure` ticks a below-floor stage ON
# the axis; without the gap, the atlas's weakest statement about a figure ("we
# do not vouch for this one") and its lowest real one ("this is the smallest
# number we measured") would land on the same pixel row, which is the exact
# conflation `_tissue_percentile` refuses to make in words.
_FLOOR_GAP: Final = 8.0

# The value axis's own two labels, in the margin `_PLOT_LEFT` already
# reserves. Without them the 104 px between the plot's top and the foot of the
# scale carried no number at all, and every gene's curve filled all of it:
# measured 2026-08-20, PKD1L1 spans 1 to 2 tpm and TBX20 spans 1 to 526, and
# the two were drawn as the same full-height excursion. Right-anchored, so a
# four-digit median (MYH6 reaches 4,761 tpm in heart) grows leftward into the
# margin rather than across the plot. The drop centres the glyphs on the
# point they label, the way `_FOREST_TEXT_DROP` centres a forest row's.
_AXIS_LABEL_X: Final = _PLOT_LEFT - 4.0
_AXIS_LABEL_DROP: Final = 3.4

# A series whose placed values are all equal has no range, and `LogScale`
# refuses a zero-span axis rather than silently placing every point on one
# pixel. One decade centred on the value is the substitute: half a decade
# either side puts a flat series in the middle of the plot, where it reads as
# "one value, no range", instead of pinned to the floor or the ceiling, where
# it would read as low or high against a scale nothing else occupies.
#
# **The condition is equal values, not few of them**, and the two are easy to
# write as one because they overlap. Measured 2026-08-20 on the committed
# corpus, three genes reach this branch on their heart trajectory: FOXH1 (one
# placed stage, 2 tpm), TFAP2B (one placed stage, 1 tpm) and MMP21 -- which is
# placed at *three* stages, all of them 1.0 tpm. A count of genes with few
# placed stages would be a different set: 9 of the 92 are placed at fewer than
# two heart stages, and 7 of those 9 place nothing at all, so `_charted_tissue`
# draws no chart and there is no axis to scale. No gene reaches this branch
# through any other path: the heart chart is the only figure here with a value
# axis, and it is fitted to the cardiac tissue's own placed medians alone.
_FLAT_SERIES_SPAN: Final = 10.0**0.5


def _even_x(index: int, count: int, left: float, right: float) -> float:
    """The `index`-th of `count` evenly spaced positions between `left` and `right`.

    Kept as its own function from when two pictures on this page plotted the
    same stage sequence at different widths (the trajectory and a row of organ
    small multiples, deleted 2026-11-02). Today one figure holds the phase
    strip, the heart chart and the organ grid, and `_stage_x` is the one place
    that fixes where a stage sits, so a column is one stage from strip to grid
    by construction.
    """
    if count < 2:
        return (left + right) / 2
    return left + (right - left) * index / (count - 1)


def _stage_x(index: int, count: int) -> float:
    """Where the `index`-th of `count` stages sits on the horizontal axis.

    The axis is the curated stage **order**, never `wpc`, and the two are
    different claims. Eight of this dataset's 21 stages are post-natal and
    carry no `wpc` at all, so a time-proportional axis could not place them;
    and even among the prenatal ones the sampling is uneven (13, then 16, 18,
    19 wpc). Evenly spaced positions say only "these came in this order",
    which is the whole of what `Stage.order` supports. `_band_caption` says so
    in words, because a reader cannot see it from the picture.
    """
    return _even_x(index, count, _PLOT_LEFT, _PLOT_RIGHT)


def _half_step(count: int) -> float:
    """Half the gap between two adjacent stages, in pixels.

    A phase band, and a bar in the phase strip, runs half a step outside the
    first and last stage it covers, so that a phase covering exactly one stage
    draws a visible band rather than a zero-width rectangle nobody sees. One
    published phase does exactly that: `outflow_tract_septum_morphogenesis` covers 6 wpc alone.
    """
    if count < 2:
        return (_PLOT_RIGHT - _PLOT_LEFT) / 2
    return (_PLOT_RIGHT - _PLOT_LEFT) / (count - 1) / 2


def _axis_bounds(values: Sequence[float]) -> tuple[float, float]:
    """A log axis wide enough to place every value, and never zero-wide."""
    low, high = min(values), max(values)
    if high > low:
        return low, high
    return low / _FLAT_SERIES_SPAN, high * _FLAT_SERIES_SPAN


def _axis_label(value: float, y: float) -> str:
    """One end of the value axis, named.

    `_fmt`, never `charts.coordinate`. `coordinate` is the *pixel* formatter
    and fixes to a tenth of a pixel, so it would publish a median of 275 tpm
    as `275.0` and one of 4,761 as `4761.0`, beside `_fmt`'s `4,761` in the
    table below the chart -- the same figure spelled two ways on one page.
    """
    return (
        f'<text class="chart-label" x="{coordinate(_AXIS_LABEL_X)}" '
        f'y="{coordinate(y + _AXIS_LABEL_DROP)}" text-anchor="end">'
        f"{html.escape(_fmt(value))}</text>"
    )


# The figure's axes. The value axis gets a line, a gridline at every whole
# decade the data spans (beside the two end labels, which stay), and a title
# naming the quantity and its unit. The stage axis labels **every** stage the
# dataset publishes, because one x position per stage is shared by the phase
# strip above, the heart chart and the organ grid below, and the order of the
# stages is the only thing the horizontal position says.
_X_TICK: Final = 4.0
_X_LABEL_DROP: Final = 12.0
_Y_TITLE_X: Final = 9.0
_GRID_CLEARANCE: Final = 9.0
_PRENATAL_STAGE: Final = re.compile(r"^(\d+) week post conception$")
_PRENATAL_TITLE: Final = "weeks post conception"
# A tick label's rotation, and the width one character of it takes. The width
# is an estimate used only to size the margin under the axis and to decide
# which side of a strip bar its label goes on. It was set by eye against a
# headless-Chrome screenshot of the TBX5 page (the longest label, "school age
# child", clears the first grid row with room to spare), never derived from
# font metrics: a different font would move it.
_TICK_ANGLE: Final = 40.0
_CHAR_WIDTH: Final = 4.6


def _short_stage(label: str) -> str:
    """A stage's short label: "4w" for "4 week post conception", else the name itself.

    Used where there is no room for the axis's own convention (a week number
    under a shared "weeks post conception" title): a hover title, a grid cell's
    `<title>`.
    """
    prenatal = _PRENATAL_STAGE.fullmatch(label)
    return f"{prenatal.group(1)}w" if prenatal else label


def _decades(low: float, high: float) -> list[float]:
    """Every power of ten inside `[low, high]`, ascending."""
    first, last = math.ceil(math.log10(low)), math.floor(math.log10(high))
    return [10.0**exponent for exponent in range(first, last + 1)]


def _value_axis(
    scale: LogScale, low: float, high: float, unit: str, top: float, bottom: float
) -> str:
    """The value axis: its line, a gridline per whole decade, and its title.

    `top` and `bottom` are the plot's own, because the phase strip above it
    makes where the plot starts a function of how many phases there are.
    """
    base = bottom - _FLOOR_GAP
    ends = (base - scale.x(low), base - scale.x(high))
    out = (
        f'<line class="chart-axis" x1="{coordinate(_PLOT_LEFT)}" y1="{coordinate(top)}" '
        f'x2="{coordinate(_PLOT_LEFT)}" y2="{coordinate(bottom)}"/>'
    )
    for decade in _decades(low, high):
        y = base - scale.x(decade)
        if any(abs(y - end) < _GRID_CLEARANCE for end in ends):
            continue
        out += (
            f'<line class="chart-grid" x1="{coordinate(_PLOT_LEFT)}" y1="{coordinate(y)}" '
            f'x2="{coordinate(_PLOT_RIGHT)}" y2="{coordinate(y)}"/>' + _axis_label(decade, y)
        )
    centre = (top + bottom) / 2
    title = f"median abundance ({unit}, log scale)" if unit else "median abundance (log scale)"
    return out + (
        f'<text class="chart-axis-title" x="{coordinate(_Y_TITLE_X)}" y="{coordinate(centre)}" '
        f'text-anchor="middle" transform="rotate(-90 {coordinate(_Y_TITLE_X)} '
        f'{coordinate(centre)})">{html.escape(title)}</text>'
    )


def _axis_depth(labels: Sequence[str]) -> float:
    """How far below the axis line its labels reach, in pixels.

    A prenatal stage prints its week number on one row, and one shared
    "weeks post conception" sits on the row below; any other stage prints its
    full name rotated, and the longest name decides how far down those reach.
    """
    prenatal = any(_PRENATAL_STAGE.fullmatch(label) for label in labels)
    others = [len(label) for label in labels if not _PRENATAL_STAGE.fullmatch(label)]
    depth = _X_LABEL_DROP + 3.0
    if prenatal:
        depth = _X_LABEL_DROP + 15.0
    if others:
        slanted = max(others) * _CHAR_WIDTH * math.sin(math.radians(_TICK_ANGLE))
        depth = max(depth, _X_LABEL_DROP + slanted + 4.0)
    return depth


def _stage_axis(labels: Sequence[str], unsampled: frozenset[int], baseline: float) -> str:
    """The stage axis: a tick and a label under **every** stage the dataset publishes.

    Until 2026-11-02 only the cardiac tissue's own sampled stages were labelled
    (as "4w", rotated) and the axis carried a title, "developmental stage
    (curated order, not elapsed time)". Two stages heart has no row at
    (`school age child`, `elderly`) therefore had no tick at all, in a figure
    whose columns are now shared with a grid that does have a cell there.

    A prenatal stage ("N week post conception") shows its number, unrotated,
    with exactly one "weeks post conception" centred under the prenatal run;
    any other stage shows its full name rotated and anchored at its end. A
    stage in `unsampled` -- one where the charted cardiac tissue has no row --
    takes the lighter `unsampled` class: it is labelled, because the grid
    below has a column there, and greyed, because the line above has nothing.
    """
    count = len(labels)
    out = ""
    prenatal_x: list[float] = []
    for index, label in enumerate(labels):
        x = _stage_x(index, count)
        y = baseline + _X_LABEL_DROP
        css = "chart-tick-label unsampled" if index in unsampled else "chart-tick-label"
        out += (
            f'<line class="chart-axis" x1="{coordinate(x)}" y1="{coordinate(baseline)}" '
            f'x2="{coordinate(x)}" y2="{coordinate(baseline + _X_TICK)}"/>'
        )
        prenatal = _PRENATAL_STAGE.fullmatch(label)
        if prenatal:
            prenatal_x.append(x)
            out += (
                f'<text class="{css}" x="{coordinate(x)}" y="{coordinate(y)}" '
                f'text-anchor="middle">{prenatal.group(1)}</text>'
            )
        else:
            out += (
                f'<text class="{css}" x="{coordinate(x)}" y="{coordinate(y)}" '
                f'text-anchor="end" transform="rotate(-{coordinate(_TICK_ANGLE)} '
                f'{coordinate(x)} {coordinate(y)})">{html.escape(label)}</text>'
            )
    if prenatal_x:
        centre = (min(prenatal_x) + max(prenatal_x)) / 2
        out += (
            f'<text class="chart-axis-title" x="{coordinate(centre)}" '
            f'y="{coordinate(baseline + _X_LABEL_DROP + 12.0)}" text-anchor="middle">'
            f"{_PRENATAL_TITLE}</text>"
        )
    return out


def _stage_label(stage: StageProfileEntry) -> str:
    """A stage's own token, or the state a null one is in.

    The one spelling of a null stage, taken from one place so a chart's axis
    and the tables' stage column cannot spell it two different ways.
    """
    return stage["stage"] if stage["stage"] is not None else "no stage recorded"


class _Median(NamedTuple):
    """One organ's figure at one stage, and -- if it is not plottable -- why not.

    A `NamedTuple` rather than a bare tuple because the fourth field is the
    one two captions were reading off the wrong axis. `value is None` says a
    stage was sampled and not placed; it does not say the gene reads low
    there, and two captions both said so (`_band_caption` and the since-deleted
    spark caption). The reason
    now travels *with* the point, so a caption cannot look it up on a
    different pass over the data than the one that decided the point was
    unplaceable -- the discipline `_adjacent_runs` keeps by taking its pixel
    points beside its stage indices.

    `gap` is `None` exactly when `value` is not, and is the raw `ProfileGap`
    token the payload recorded rather than a sentence: what the wording
    should be is `_PLACEMENT_GAP_CLAUSE`'s job, in one place, and this module
    already has three maps saying so.

    `position` rather than `index`, which is `tuple.index` and cannot be a
    field name on a `NamedTuple` at all -- mypy refuses it outright, which is
    the only reason this is not a silently shadowed method.
    """

    position: int
    label: str
    value: float | None
    gap: ProfileGap | None


def _tissue_medians(entry: DatasetProfileEntry, tissue: str) -> list[_Median]:
    """One organ's series: `(stage index, stage label, plottable median, gap)`.

    `None` is a stage this dataset **sampled and did not place** -- below the
    detection floor, or with no complete percentile grid behind it. Those
    stages stay in the list rather than being filtered out, because a gap in
    the sequence and a stage nobody sampled must not render the same way:
    `_bulk_figure` ticks the first on the axis and the second is simply not
    here at all.

    A stage where this organ has no row is not in the list, for that reason.
    Measured 2026-08-20 on the committed corpus: heart has a row at 19 of the
    21 stages this dataset publishes, missing at `school age child` and
    `elderly`, for all 92 published genes.

    **What is plotted is exactly what the atlas placed.** A below-floor
    figure has a median too, and plotting it would put a number this atlas
    declines to vouch for on the same line as ones it does -- the discipline
    `_tissue_percentile` already keeps in words. A non-positive median folds
    into the same `None` for a second, narrower reason: it cannot be located
    on a log axis at all. Neither loses a figure; every one of them is in the
    table inside the `<details>`.
    """
    series: list[_Median] = []
    for index, stage in enumerate(entry["stages"]):
        for measured in stage["tissues"]:
            if measured["tissue"] != tissue:
                continue
            plottable = measured["median_abundance"]
            placed = measured["placement"] is not None and plottable > 0
            series.append(
                _Median(
                    position=index,
                    label=_stage_label(stage),
                    value=plottable if placed else None,
                    # `or` rather than a bare read: a cell the payload placed
                    # whose median is `<= 0` carries no `not_placed_reason` by
                    # `TissueProfileEntry`'s own contract, and would otherwise
                    # reach `_Median` with neither a value nor a gap -- the one
                    # pairing its docstring rules out.
                    gap=(
                        None
                        if placed
                        else (measured["not_placed_reason"] or ProfileGap.NOT_ON_A_LOG_AXIS)
                    ),
                )
            )
    return series


def _tissue_unit(entry: DatasetProfileEntry, tissue: str) -> str:
    """The unit this organ's medians are quoted in, or `""` if it has none.

    Read off the measurements rather than off the `Dataset`, for the same
    reason `_abundance` reads it there: the unit is a property of the mirrored
    row, and a page must not quote a floor in one unit beside a median in
    another.
    """
    for stage in entry["stages"]:
        for measured in stage["tissues"]:
            if measured["tissue"] == tissue:
                return measured["unit"]
    return ""


def _glance_line(entry: DatasetProfileEntry, tissue: str, floor: float | None) -> str:
    """The one-sentence summary a reader meets before any chart, for one cardiac tissue.

    Every number is read off the bundle's own fields and nothing is re-derived
    (spec B1), so the sentence says no more than the payload publishes. Returns
    `""` when no stage of `tissue` is placed: there is then no peak to state,
    and `_dataset_block` says why in the lede instead.

    **The floor clause is conditional, and the condition is read off the gaps.**
    "Above the floor" is true only when every stage this tissue sampled and did
    not place was below the floor. A stage left unplaced for a reference reason
    (no percentile grid) is a hole in the atlas, not a low reading, so one such
    stage forces the neutral "placed at" wording.

    **The highest-count carries no denominator, deliberately.** It counts the
    stages whose published `highest_in` is this tissue and does not say "of N".
    tau is absent at some sampled stages (measured 2026-11-02 on the committed
    corpus: 28 of 119 gene bundles carry a dataset stage with no tau), so the
    stages where the tissue could have been highest are not the stages
    sampled; and a tie is not a win (`highest_in` is `None` there: 120 of
    the 2,499 dataset stages), so the stages where it was not highest are not
    the complement either. "Of 19" would be false for both reasons.

    **The peak is named at its earliest stage.** Stages arrive in developmental
    order, so the first placed stage carrying the peak value is the earliest;
    any later equal value is counted as a tie, not as a further peak.
    """
    series = _tissue_medians(entry, tissue)
    placed = [measured for measured in series if measured.value is not None]
    if not placed:
        return ""
    unit = _tissue_unit(entry, tissue)
    below_floor_only = all(
        measured.gap == ProfileGap.BELOW_DETECTION_FLOOR
        for measured in series
        if measured.value is None
    )
    if floor is not None and below_floor_only:
        words = html.escape(_floor_words(floor, unit), quote=False)
        clause = f"at or above the {words} floor at {len(placed)} of {len(series)} stages sampled"
    else:
        clause = f"placed at {len(placed)} of {len(series)} stages sampled"

    values = [measured.value for measured in placed if measured.value is not None]
    peak_value = max(values)
    at_peak = [measured for measured in placed if measured.value == peak_value]
    earliest = at_peak[0]
    other_peaks = len(at_peak) - 1
    if other_peaks == 0:
        tie_clause = ""
    elif other_peaks == 1:
        tie_clause = " and 1 other stage"
    else:
        tie_clause = f" and {other_peaks} other stages"
    prenatal = _PRENATAL_STAGE.fullmatch(earliest.label)
    stage = f"{prenatal.group(1)} weeks post conception" if prenatal else earliest.label

    wins = sum(
        1
        for stage_entry in entry["stages"]
        if stage_entry["specificity"] is not None
        and stage_entry["specificity"]["highest_in"] == tissue
    )
    if wins == 0:
        highest = "no stage"
    elif wins == 1:
        highest = "1 stage"
    else:
        highest = f"{wins} stages"

    name = html.escape(tissue[:1].upper() + tissue[1:], quote=False)
    return (
        f'<p class="glance"><strong>{name}:</strong> {clause}; '
        f"peak {_fmt(peak_value)} {html.escape(unit, quote=False)} at "
        f"{html.escape(stage, quote=False)}{tie_clause}; "
        f"the highest of the organs sampled at {highest}.</p>"
    )


def _banded_phases(
    entry: DatasetProfileEntry, phases: CardiacPhaseFile | None
) -> list[tuple[CardiacPhase, int, int]]:
    """Each curated phase this dataset's stages fall inside, and the span it covers.

    **A phase whose `end_basis` is `NOT_STATED` is skipped, unconditionally.**
    It has no curated width, so a band drawn for it would have to run to an
    edge nobody stated -- asserting exactly the boundary that field exists to
    refuse, and asserting it in the one encoding a reader takes at a glance
    rather than reads. `heart_looping` is the live case: the source states
    when it begins and never when it ends.

    That duplicates a rule `CardiacPhaseFile.phases_for` already keeps, and
    the duplication is deliberate. `phases_for` is why no stage the real
    pipeline builds ever names an unended phase; this function is why that
    stays true of the *picture* even if some later caller assembles
    `phase_ids` another way. A guard added to one layer is not a guard.

    Sorted by phase id -- and **not as a determinism guard**, which is what
    this sentence claimed until it was measured. `covered` is a dict, filled
    by walking a list of stages and, inside each, `phase_ids` in the
    `(start_wpc, id)` order `phases_for` fixes; no `set` is iterated anywhere,
    so its own key order is already reproducible. Measured 2026-08-20 over six
    explicit `PYTHONHASHSEED` values (0, 1, 7, 42, 12345, 99991): one distinct
    build with the sort and one without.

    What the sort does buy is that the band order is a function of the covered
    ids alone rather than of the stage walk that found them, so it does not
    move when a dataset's stage list does. The two orders are genuinely
    different: dropping the sort changed 85 gene pages and the manifest, and
    reversing this dataset's 21 stages reorders the unsorted one while leaving
    the sorted one alone.
    """
    if phases is None:
        return []
    stated = {phase.id: phase for phase in phases.phases if phase.end_wpc is not None}
    covered: dict[str, list[int]] = {}
    for index, stage in enumerate(entry["stages"]):
        for phase_id in stage["phase"]["phase_ids"]:
            if phase_id in stated:
                covered.setdefault(phase_id, []).append(index)
    return [
        (stated[phase_id], min(covered[phase_id]), max(covered[phase_id]))
        for phase_id in sorted(covered)
    ]


def _phase_bands(
    entry: DatasetProfileEntry,
    phases: CardiacPhaseFile | None,
    count: int,
    top: float,
    bottom: float,
) -> str:
    """The faint shaded stage ranges behind the heart chart.

    **This is what makes the chart admissible under D43.** The curve itself is
    what the source's own browser already draws; the phases are a claim no
    external browser makes, transcribed from PMID:32048790 into
    `curation/cardiac_phases.yaml` and carried onto every gene page here.
    `_bulk_figure` draws no heart chart when this returns `""`, so the rule is
    enforced by construction rather than by a comment asking a later author to
    keep it.

    The bands carry no `<title>` since the phase strip above the plot was added
    (2026-11-02): each phase is now named once, in the strip's own row, and a
    second `<title>` per band would name it twice. They overlap, deliberately
    -- three cardiac processes run concurrently at 6 wpc -- and share one fill,
    merging into a single shaded region marking the morphogenetic window.
    """
    edge = _half_step(count)
    rects = []
    for _, first, last in _banded_phases(entry, phases):
        left = max(_stage_x(first, count) - edge, 0.0)
        right = min(_stage_x(last, count) + edge, float(_CHART_WIDTH))
        rects.append(
            f'<rect class="chart-band" x="{coordinate(left)}" y="{coordinate(top)}" '
            f'width="{coordinate(right - left)}" height="{coordinate(bottom - top)}"/>'
        )
    return "".join(rects)


_STRIP_TOP: Final = 4.0
_STRIP_ROW: Final = 13.0
_STRIP_BAR: Final = 8.0
_STRIP_TO_PLOT: Final = 10.0


def _phase_strip(
    rows: Sequence[tuple[CardiacPhase, int, int]], count: int, top: float
) -> tuple[str, float]:
    """One named bar per banded phase, in chronological order, and the strip's height.

    **Chronological, not `_banded_phases`' id order.** That function sorts by
    phase id so the order does not move with the stage walk, and an id order
    is alphabetical, which says nothing about when a phase happens. A reader
    follows a strip as a timeline, so the rows are ordered by first stage, then
    last stage, then id -- the id only breaking a tie, so the order is still a
    function of the data alone. `_banded_phases` is left alone: other callers
    rely on its order.

    Each bar spans the same half-step-padded range the faint band behind the
    plot does, so the two cannot disagree; its label sits right of the bar, or
    left of it when the label would run off the figure. A bar carries a
    `<title>` with the same label, for a reader who hovers.
    """
    edge = _half_step(count)
    out = ""
    ordered = sorted(rows, key=lambda row: (row[1], row[2], row[0].id))
    for position, (phase, first, last) in enumerate(ordered):
        y = top + position * _STRIP_ROW
        left = max(_stage_x(first, count) - edge, 0.0)
        right = min(_stage_x(last, count) + edge, float(_CHART_WIDTH))
        label = html.escape(phase.label)
        if right + 4.0 + len(phase.label) * _CHAR_WIDTH <= _CHART_WIDTH:
            text_x, anchor = right + 4.0, "start"
        else:
            text_x, anchor = left - 4.0, "end"
        out += (
            f'<rect class="chart-strip" x="{coordinate(left)}" y="{coordinate(y)}" '
            f'width="{coordinate(right - left)}" height="{coordinate(_STRIP_BAR)}" rx="3">'
            f"<title>{label}</title></rect>"
            f'<text class="chart-tick-label" x="{coordinate(text_x)}" '
            f'y="{coordinate(y + _STRIP_BAR / 2 + _AXIS_LABEL_DROP)}" '
            f'text-anchor="{anchor}">{label}</text>'
        )
    return out, len(ordered) * _STRIP_ROW


def _other_gap_clause(reasons: frozenset[str]) -> str:
    """Every gap reason *except* the detection floor's, spelled out.

    The floor is excluded because it is the one reason with its own required
    wording -- `_tissue_percentile`'s rule 5 and the short form `_band_caption`
    carries on its visible line (E1). Everything else
    resolves through `_PLACEMENT_GAP_CLAUSE`, so a caption and the table cell
    three lines below it cannot spell one gap two ways.

    Sorted because `reasons` is a `frozenset` and an unsorted join over one
    would reorder a published sentence between two builds of one commit --
    **and, measured 2026-08-20, that cannot happen today.** `_tissue_entry`
    assigns `dataset_not_registered` and `detection_floor_undeclared` from
    `_dataset_gap`, which reads the `Dataset` alone, so either one applies to
    every cell of that dataset or to none; only `no_quantile_grid` is
    per-cell, and it is reached only where `_dataset_gap` returned nothing.
    Two distinct non-floor reasons therefore cannot co-occur inside one
    `DatasetProfileEntry`, this set holds at most one member, and removing
    the sort changes no published byte. The sort stays as the defensive kind
    section 9 keeps deliberately: this module renders whatever the payload
    says and cannot see an invariant that lives in `profiles.py`. It is not
    guarded by a test, and that is the reason -- `_banded_phases`' own
    docstring records the same correction, made after its sort was claimed as
    a determinism guard and measured not to be one.

    Returns `""` when the floor is the only reason present, which is every
    published page today -- measured 2026-08-20, all 2,817 unplaced cells on
    the committed corpus are below-floor.
    """
    others = sorted(reason for reason in reasons if reason != ProfileGap.BELOW_DETECTION_FLOOR)
    if not others:
        return ""
    return "; ".join(
        _PLACEMENT_GAP_CLAUSE.get(reason, "no reason was recorded for this gap")
        for reason in others
    )


def _chart_notes(notes: str) -> str:
    """An expression chart's folded key, under "How to read this chart" (D57).

    One implementation for every expression chart's toggle, so the summary and
    the link cannot be worded two ways. (Until 2026-11-02 a page carried two
    of them, the trajectory's and the small multiples'; there is one fold per
    figure now.) The link names its section, as the forest's toggle names `#rare-variants`: the
    nav's "How to read this atlas" and the gene card's "How to read this page"
    are the two texts that open the top of `methods.html`.
    Both paragraphs are `method`, because unlike the forest's these sit outside
    a `<figcaption>`, and a bare `<p>` would render at body size under the
    caption line it continues.
    """
    return (
        '<details class="chart-notes"><summary>How to read this chart</summary>'
        f'<p class="method">{notes}</p>'
        '<p class="method"><a href="../methods.html#expression">How to read expression '
        "figures</a></p></details>"
    )


_GRID_BINS_NOTE: Final = "The grid's shades are fixed bins, so they compare across genes."


def _band_caption(
    entry: DatasetProfileEntry,
    phases: CardiacPhaseFile | None,
    tissue: str,
    floor: float | None,
    low: float,
    high: float,
    series: Sequence[_Median],
    *,
    lined: bool,
    uncharted: Sequence[str] = (),
) -> str:
    """What the figure cannot say about itself, in words beside it.

    The horizontal axis is stage *order*, not elapsed time, so the gap between
    13 and 16 wpc is drawn the same width as the gap between 4 and 5. Saying
    so is the difference between an ordinal axis and a false linear one.

    **And the vertical axis is fitted to this gene alone**, which is the
    second trap. Every gene's curve spans the same 104 units whatever its range
    (the 112-unit plot less its 8-unit floor gap):
    measured 2026-08-20, PKD1L1 runs 1 to 2 tpm in heart and TBX20 runs 1 to
    526, and the two were drawn as the same full-height excursion under a
    caption that named neither number. The range is printed so a reader can
    tell which axis they are looking at -- and the grid below the chart is
    the opposite, fixed bins, which the fold says in its last sentence.

    **It fits the cardiac tissue alone since 2026-11-02.** The axis also
    spanned every other organ's dashed line for as long as those were drawn,
    and the caption said "across every organ drawn"; the dashed lines are gone
    (the grid carries the other organs), so the range is the tissue's own
    placed medians, and the clause is gone with them.

    **Both break clauses are conditional on the break being in this figure.**
    The tick clause named a mark most pages do not draw: measured 2026-08-20,
    all 85 trajectory pages said "the stage is ticked on the axis instead"
    while **55 drew zero ticks**. A line also breaks where this dataset has no
    row for that organ at that stage, which for E-MTAB-6814 is the same set of
    stages for every gene (heart has none at `school age child` or `elderly`).
    `_tissue_medians`' docstring insists that state and a below-floor one "must
    not render the same way", and in the picture they did.

    **It is stated as the axis's range, never as the medians'**, and the two
    differ on exactly the genes where the wording matters most. `_axis_bounds`
    widens a flat series by half a decade either side, because `LogScale`
    refuses a zero-span axis; 3 of the 85 charted genes were flat then (TFAP2B
    and MMP21 at 1 tpm, FOXH1 at 2), and "median abundance in heart, 0.316 to
    3.16 tpm" would attribute to their medians a spread neither end of which
    anyone measured.

    **One visible sentence and the rest behind a toggle (D57, 2026-10-26).**
    The floor clause came back out of the fold (owner decision 2026-10-26,
    E1), in a short form on the visible line and only where it fires. A tick is
    a `marker` centred on the axis line, so it says "on the axis", and it
    keeps "beneath the foot of the scale rather than at a value on it", which
    stops a reader reading a tick as a measurement. The fold carries the tick
    clause, the no-row clause, the phases' attribution and the grid's bins; it
    closes on a link to `methods.html#expression`. **The fold changed no
    clause's condition.**

    `uncharted` names any further cardiac tissue with a placed median: the
    chart draws one tissue, and a second one's line would be silently absent
    from it. The corpus has one cardiac tissue per dataset (E-MTAB-6814:
    heart), so this is latent and fixture-tested; the grid and the tables
    carry the others in full, and the fold says so.
    """
    attribution = html.escape(phases.attributed_to) if phases is not None else ""
    unit = html.escape(_tissue_unit(entry, tissue))
    # Every stage this dataset sampled and did not place is ticked, whatever
    # the recorded reason -- so the clause explaining the mark is keyed on the
    # reasons that actually occur, never on `value is None`. Read off the
    # whole series rather than off `_line_breaks`, which sees only the gaps
    # that interrupt a *line*: a leading or trailing unplaced stage draws a
    # tick and breaks nothing, and it still needs explaining.
    ticked = frozenset(point.gap or "" for point in series if point.value is None)
    # Visible, on the line above the toggle, and only where it fires (E1).
    #
    # **Where a tick also stands for another reason, the visible line names it
    # too** (review, 2026-10-27). "Ticks on the axis mark stages below this
    # dataset's detection floor" is universal, and beside a tick drawn because
    # no percentile grid was published it is false of that tick. So with both
    # reasons present the sentence is a disjunction true of every tick.
    # Latent: all 2,830 unplaced cells in the published bundles are below the
    # floor (measured 2026-10-27), so no published page takes this form.
    other_clause = _other_gap_clause(ticked)
    floor_clause = ""
    tick_clause = ""
    below_floor = ProfileGap.BELOW_DETECTION_FLOOR.value in ticked
    if floor is not None and below_floor and other_clause:
        floor_clause = (
            " Ticks on the axis, beneath the foot of the scale rather than at a value on it, "
            f"mark stages below this dataset's detection floor ({_floor_words(floor, unit)}) or "
            f"where {other_clause}."
        )
    elif floor is not None and below_floor:
        floor_clause = (
            f" Ticks on the axis mark stages below this dataset's detection floor "
            f"({_floor_words(floor, unit)}), beneath the foot of the scale rather than at a "
            "value on it."
        )
    elif other_clause:
        tick_clause = (
            f" A ticked stage is one this atlas did not place: {other_clause}. The mark "
            "says nothing about how high the gene reads there."
        )
    _, no_row = _line_breaks(series)
    # Two spellings, because the sentence contrasts itself with the floor
    # clause and that clause is not always there. "A break with no tick under
    # it is different" on a figure drawing no tick at all sends a reader
    # hunting for one (CLAUDE.md section 4.27).
    gap_clause = ""
    if lined and no_row:
        gap_clause = (
            " A break in the line with no tick under it is different: this dataset has "
            "<strong>no row for this organ at that stage</strong>, so there is no measurement "
            "to place, low or otherwise."
            if floor_clause
            else " The line breaks where this dataset has <strong>no row for this organ at "
            "that stage</strong>: there is no measurement to place, low or otherwise."
        )
    uncharted_clause = ""
    if uncharted:
        names = ", ".join(html.escape(name) for name in uncharted)
        uncharted_clause = (
            f" Only {html.escape(tissue)} is charted; {names} appear in the grid and the tables."
        )
    visible = (
        f'<p class="method">Median {unit or "abundance"} in {html.escape(tissue)} per stage on a '
        f"log scale from {_fmt(low)} to {_fmt(high)}, fitted to this gene, so heights "
        "compare nothing across genes; stages are in curated order, evenly spaced, "
        f"not elapsed time.{floor_clause}</p>"
    )
    notes = (
        f"{tick_clause}{gap_clause}{uncharted_clause} Phases above the chart: {attribution}. "
        f"{_GRID_BINS_NOTE}"
    ).strip()
    return visible + _chart_notes(notes)


def _adjacent_runs(
    placed: Sequence[tuple[int, float]], points: Sequence[tuple[float, float]]
) -> list[list[tuple[float, float]]]:
    """`points` split wherever the stages behind them are not consecutive.

    A break is any stage the atlas did not place sitting between two it did --
    below this dataset's detection floor, or with no percentile grid behind
    it -- and also a stage where this organ has no row at all. The two are
    different facts and `_bulk_figure` renders them differently (the first gets
    an axis tick, the second nothing), but they are the same fact *for the
    line*: there is no measurement here to draw through.

    Takes the pixel points beside the stage indices rather than recomputing
    them, so the run boundaries and the markers cannot come from two different
    passes over the series.
    """
    runs: list[list[tuple[float, float]]] = []
    previous: int | None = None
    for (index, _), point in zip(placed, points, strict=True):
        if previous is None or index != previous + 1:
            runs.append([])
        runs[-1].append(point)
        previous = index
    return runs


def _line_breaks(series: Sequence[_Median]) -> tuple[frozenset[str], bool]:
    """Which causes actually interrupt this organ's drawn line:
    `(the recorded reasons for a stage sampled and not placed, whether a stage
    has no row at all)`.

    **The first half is a set of reasons and was a bare `bool`.** "A stage
    sampled and not placed" is not one fact -- below a detection floor is a
    low reading about the gene, a missing percentile grid is a hole in the
    reference -- and a caller handed only `True` can say nothing but the
    first. A caption (since deleted) said exactly that, of every break,
    whatever the payload recorded.

    `_adjacent_runs` splits the line on both and cannot tell them apart,
    because for the *line* they are the same fact. For the *caption* they are
    not: the first is a measurement this atlas declines to vouch for and gets
    a tick on the axis, the second is the dataset never having looked, and
    `_tissue_medians`' docstring is explicit that they "must not render the
    same way". They did -- to a reader, a break with a tick and a break
    without one differed only in a mark the caption promised on every page
    and drew on 30 of 85.

    Only a gap *between two placed stages* counts, for either cause. A stage
    the dataset skipped before the first placed one or after the last breaks
    no line, so nothing on the figure needs explaining.
    """
    recorded = {point.position: point.gap for point in series}
    placed = [point.position for point in series if point.value is not None]
    unplaced: set[str] = set()
    unsampled = False
    for first, second in zip(placed, placed[1:], strict=False):
        for gap in range(first + 1, second):
            if gap in recorded:
                unplaced.add(recorded[gap] or "")
            else:
                unsampled = True
    return frozenset(unplaced), unsampled


def _charted_tissue(
    entry: DatasetProfileEntry, cardiac: frozenset[str], phases: CardiacPhaseFile | None
) -> str | None:
    """The cardiac tissue the heart chart draws, or `None` where there is no chart.

    `None` is the caller's cue to say something in words instead
    (`_unbanded_sentence`, `_no_trajectory_sentence`), in each of the cases
    where a chart would assert more than the data does: nothing placed for any
    cardiac tissue, and no curated phase band to draw, which is the bare re-plot
    D43 forbids. **They are different facts and `_dataset_block` owes a reader
    a different sentence for each**; it asks `_anything_placed` to tell them
    apart.

    The first cardiac tissue, sorted, with a placed median. Sorted because
    `cardiac` is a `frozenset`: an unsorted pick would draw a different organ
    between two builds of one commit.
    """
    if not _banded_phases(entry, phases):
        return None
    for tissue in sorted(cardiac):
        if any(point.value is not None for point in _tissue_medians(entry, tissue)):
            return tissue
    return None


def _organ_order(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> list[str]:
    """Every organ this dataset sampled for this gene: cardiac first, then the rest.

    Each group sorted, because both come from a `set` and the grid's row order
    is published. Cardiac first because the rows below the heart chart should
    read as "the organs this atlas is about, then the context" -- and because
    the declaration is this atlas's own (`Dataset.cardiac_tissues`), not the
    source's.
    """
    sampled = _sampled_tissues(entry)
    return [t for t in sampled if t in cardiac] + [t for t in sampled if t not in cardiac]


# The grid's fixed bins, in the dataset's own unit: a cell's shade means the
# same abundance on every gene's page, which is what the heart chart's
# per-gene axis cannot offer. Edges are 10, 100 and 1,000.
_HEAT_EDGES: Final = (10.0, 100.0, 1000.0)
_GRID_ROW: Final = 14.0
_GRID_CELL: Final = 12.0
_GRID_GAP_ABOVE: Final = 6.0
# A column is as wide as its stage's slot, capped: a dataset with a handful of
# stages would otherwise draw cells wider than the grid is tall.
_CELL_HALF_MAX: Final = 20.0


def _grid_state(cell: TissueProfileEntry) -> str:
    """The class one grid cell takes: a shade for a placed median, else why not.

    **A placed cell is binned on its median; an unplaced one is never binned.**
    A below-floor figure has a median too (0.4 tpm), and shading it `heat-1`
    would put a number this atlas declines to vouch for on the same ramp as
    ones it does -- the discipline `_tissue_medians` keeps for the line and
    `_tissue_percentile` keeps in words. `below_detection_floor` is the one
    reason with its own state (rule 5); any other reason is a hole in the
    reference, not a low reading, and gets the other.
    """
    if cell["placement"] is not None:
        value = cell["median_abundance"]
        for number, edge in enumerate(_HEAT_EDGES, start=1):
            if value < edge:
                return f"heat-{number}"
        return "heat-4"
    if cell["not_placed_reason"] == ProfileGap.BELOW_DETECTION_FLOOR.value:
        return "heat-floor"
    return "heat-gap"


def _cell_title(
    organ: str, stage: str, cell: TissueProfileEntry, state: str, floor: float | None
) -> str:
    """What hovering one grid cell says: organ, stage, value, and its state if unplaced."""
    unit = cell["unit"]
    text = f"{organ} \u00b7 {_short_stage(stage)} \u00b7 " + " ".join(
        part for part in (_fmt(cell["median_abundance"]), unit) if part
    )
    if state == "heat-floor":
        words = _floor_words(floor, unit)
        text += f" \u00b7 below the {words} floor" if words else " \u00b7 below the detection floor"
    elif state == "heat-gap":
        reason = cell["not_placed_reason"] or ""
        text += " \u00b7 " + _PLACEMENT_GAP_CLAUSE.get(
            reason, "no reason was recorded for this gap"
        )
    return html.escape(text)


class _Grid(NamedTuple):
    """The grid's markup, and which of its optional states it drew.

    The legend names only the states that occur (CLAUDE.md section 4.27: a key
    for a mark that is not drawn sends the reader hunting), so the grid says
    what it drew rather than the legend re-reading the payload on a second
    pass that could disagree.
    """

    svg: str
    floor: bool
    gap_reasons: frozenset[str]
    blank: bool
    shaded: bool


def _grid(
    entry: DatasetProfileEntry,
    organs: Sequence[str],
    cardiac: frozenset[str],
    floor: float | None,
    top: float,
) -> _Grid:
    """One row per organ, one column per stage, one cell per row the payload has.

    **A stage where an organ has no row draws nothing.** That is not a cell
    with a blank state: the dataset never looked, and a drawn empty cell would
    read as a measurement of nothing (`_tissue_medians` keeps the same line
    for the heart series). Cardiac organ labels are bold.
    """
    count = len(entry["stages"])
    half = min(_half_step(count), _CELL_HALF_MAX)
    out = ""
    drawn_floor = False
    drawn_shade = False
    gap_reasons: set[str] = set()
    cells = 0
    for row, organ in enumerate(organs):
        y = top + row * _GRID_ROW
        css = "chart-tick-label cardiac" if organ in cardiac else "chart-tick-label"
        out += (
            f'<text class="{css}" x="{coordinate(_PLOT_LEFT - half - 4.0)}" '
            f'y="{coordinate(y + _GRID_CELL / 2 + _AXIS_LABEL_DROP)}" text-anchor="end">'
            f"{html.escape(organ)}</text>"
        )
        for index, stage in enumerate(entry["stages"]):
            for cell in stage["tissues"]:
                if cell["tissue"] != organ:
                    continue
                cells += 1
                state = _grid_state(cell)
                drawn_floor = drawn_floor or state == "heat-floor"
                drawn_shade = drawn_shade or state[5:].isdigit()
                if state == "heat-gap":
                    gap_reasons.add(cell["not_placed_reason"] or "")
                x = _stage_x(index, count)
                title = _cell_title(organ, _stage_label(stage), cell, state, floor)
                out += (
                    f'<rect class="{state}" x="{coordinate(x - half + 1.0)}" y="{coordinate(y)}" '
                    f'width="{coordinate(2 * half - 2.0)}" height="{coordinate(_GRID_CELL)}" '
                    f'rx="2"><title>{title}</title></rect>'
                )
                if state == "heat-floor":
                    out += (
                        f'<circle class="cell-dot" cx="{coordinate(x)}" '
                        f'cy="{coordinate(y + _GRID_CELL / 2)}" r="1.6"/>'
                    )
    return _Grid(out, drawn_floor, frozenset(gap_reasons), cells < len(organs) * count, drawn_shade)


def _heat_legend(grid: _Grid, floor: float | None, unit: str, greyed: str | None) -> str:
    """The grid's key, in HTML: a swatch per state that occurs, and no others.

    The four shade labels are the bins' edges in the dataset's unit; the first
    starts at the floor ("1-10") because a placed median cannot be below it,
    and says "under 10" where the floor is unknown or not below the bin's edge.

    `greyed` is the charted cardiac tissue when the figure greys at least one
    of its stage labels (a stage that tissue has no row for). The grey label
    has no swatch of its own, so its entry names the tissue in text alone.
    """
    unit_text = f" {html.escape(unit)}" if unit else ""
    first = f"{_fmt(floor)}\u201310" if floor is not None and floor < _HEAT_EDGES[0] else "under 10"
    # The four shades are keyed only where the grid draws at least one placed
    # cell: a key for a shade that is not drawn sends the reader hunting
    # (GDF1 reads below the floor everywhere and draws no shaded cell).
    entries = (
        [
            ("heat-1", first),
            ("heat-2", "10\u2013100"),
            ("heat-3", "100\u20131,000"),
            ("heat-4", f"1,000+{unit_text}"),
        ]
        if grid.shaded
        else []
    )
    if grid.floor:
        words = html.escape(_floor_words(floor, unit), quote=False)
        entries.append(
            ("heat-floor", f"below the {words} floor" if words else "below the detection floor")
        )
    if grid.gap_reasons:
        entries.append(
            (
                "heat-gap",
                _other_gap_clause(grid.gap_reasons) or "no reason was recorded for this gap",
            )
        )
    if grid.blank:
        entries.append(("heat-blank", "blank: not sampled"))
    items = "".join(
        f'<span class="legend-item"><span class="heat-swatch {state}"></span> {text}</span> '
        for state, text in entries
    )
    if greyed is not None:
        items += f'<span class="legend-item">grey stage: no {html.escape(greyed)} sample</span> '
    return f'<p class="strip-legend heat-key">{items.rstrip()}</p>'


def _bulk_figure(
    entry: DatasetProfileEntry, dataset: Dataset | None, phases: CardiacPhaseFile | None
) -> str:
    """One dataset's figure: a named phase strip, the heart over time, every organ by stage.

    One `svg_figure`, one `<title>`. Every part reads the same `_stage_x`, so
    a column is one stage from the strip's bars through the heart chart's
    points and the axis labels to the grid's cells, by construction rather than
    by three pieces of arithmetic agreeing.

    Top to bottom: the **phase strip** (one bar per banded phase, named,
    chronological), the **heart chart** (the cardiac tissue's placed medians on
    a log axis fitted to those medians alone; a line only through adjacent
    placed stages, a tick on the axis for a stage sampled and not placed), the
    **stage axis** (every stage labelled), and the **organ grid** (cardiac
    organs first, shades in fixed bins). Then the grid's HTML legend and the
    caption with its one fold.

    **The figure survives a missing chart.** Where nothing is placed in the
    cardiac tissue, or no phase band exists (`_charted_tissue` is `None`), the
    strip and the plot are not drawn, the stage axis moves to the top and the
    grid is drawn in full; `_dataset_block` keeps its lede sentence for that
    case. Returns `""` only where no organ has any row at all: there is then
    no stage to label and nothing to say.

    **What was deleted, and why it is not missed.** The dashed other-organ
    lines (the heart chart's axis had to span them, so it compared nothing for
    a gene whose heart reads low) and the five sparkline panels (a second
    figure re-asking the question the grid now answers with a column shared
    with the heart chart). Their conditions were not carried over: a mark not
    drawn gets no clause.
    """
    cardiac = frozenset(dataset.cardiac_tissues) if dataset is not None else frozenset()
    organs = _organ_order(entry, cardiac)
    if not organs:
        return ""
    stages = entry["stages"]
    count = len(stages)
    labels = [_stage_label(stage) for stage in stages]
    floor = dataset.detection_floor if dataset is not None else None
    tissue = _charted_tissue(entry, cardiac, phases)

    body = ""
    caption = ""
    unsampled: frozenset[int] = frozenset()
    summary = f"{entry['dataset']}, {labels[0]} to {labels[-1]}"
    if tissue is None:
        baseline = _STRIP_TOP
        # Completed once the grid is built: whether any cell is shaded is the
        # grid's to say, and a title promising shades on a grid with none
        # (GDF1) would describe a mark the figure does not draw.
        title = f"Median abundance by organ and developmental stage in {summary}: every organ "
    else:
        series = _tissue_medians(entry, tissue)
        placed = [(point.position, point.value) for point in series if point.value is not None]
        unsampled = frozenset(range(count)) - {point.position for point in series}
        rows = _banded_phases(entry, phases)
        strip, strip_height = _phase_strip(rows, count, _STRIP_TOP)
        top = _STRIP_TOP + strip_height + _STRIP_TO_PLOT
        bottom = top + _PLOT_HEIGHT
        baseline = bottom
        low, high = _axis_bounds([value for _, value in placed])
        scale = LogScale(low=low, high=high, left=0.0, width=bottom - _FLOOR_GAP - top)
        points = [
            (_stage_x(index, count), bottom - _FLOOR_GAP - scale.x(value))
            for index, value in placed
        ]
        unit = _tissue_unit(entry, tissue)
        body += (
            _phase_bands(entry, phases, count, top, bottom)
            + strip
            + f'<line class="chart-axis" x1="{coordinate(_PLOT_LEFT)}" '
            f'y1="{coordinate(bottom)}" x2="{coordinate(_PLOT_RIGHT)}" '
            f'y2="{coordinate(bottom)}"/>'
            + _axis_label(high, bottom - _FLOOR_GAP - scale.x(high))
            + _axis_label(low, bottom - _FLOOR_GAP - scale.x(low))
            + _value_axis(scale, low, high, unit, top, bottom)
            + "".join(
                marker(_stage_x(point.position, count), bottom, css_class="chart-absent")
                for point in series
                if point.value is None
            )
        )
        enough = len(points) >= _STAGES_FOR_A_TRAJECTORY
        lines = (
            "".join(
                polyline(run, css_class="chart-line")
                for run in _adjacent_runs(placed, points)
                if len(run) >= 2
            )
            if enough
            else ""
        )
        body += (
            lines
            + "".join(marker(x, y, css_class="chart-point") for x, y in points)
            + _hover_dots(tissue, placed, points, labels, unit)
        )
        floor_words = (
            f", at or above a detection floor of {_floor_words(floor, unit)}"
            if floor is not None
            else ""
        )
        title = (
            f"Median abundance in {tissue}, {series[0].label} to {series[-1].label}: "
            f"{len(placed)} of {len(series)} sampled developmental stages placed against this "
            f"dataset's percentile grid{floor_words}, with this atlas's own curated cardiac "
            "phases above the chart and every organ sampled in a grid below it."
        )
        uncharted = [
            other
            for other in sorted(cardiac)
            if other != tissue and any(p.value is not None for p in _tissue_medians(entry, other))
        ]
        # `lined` rather than a second look at the data: the caption's clause
        # about a break in the line may only be published if there *is* a line,
        # and only this function knows whether one was drawn.
        caption = _band_caption(
            entry, phases, tissue, floor, low, high, series, lined=bool(lines), uncharted=uncharted
        )

    body += _stage_axis(labels, unsampled, baseline)
    if tissue is None:
        body += (
            f'<line class="chart-axis" x1="{coordinate(_PLOT_LEFT)}" y1="{coordinate(baseline)}" '
            f'x2="{coordinate(_PLOT_RIGHT)}" y2="{coordinate(baseline)}"/>'
        )
    grid_top = baseline + _axis_depth(labels) + _GRID_GAP_ABOVE
    grid = _grid(entry, organs, cardiac, floor, grid_top)
    if tissue is None:
        title += (
            "sampled, one cell per stage sampled, shaded in fixed bins."
            if grid.shaded
            else "sampled, one cell per stage sampled; no median is placed, so no cell is shaded."
        )
    height = int(math.ceil(grid_top + len(organs) * _GRID_ROW + 4.0))
    figure = svg_figure(width=_CHART_WIDTH, height=height, title=title, body=body + grid.svg)
    if tissue is None and grid.shaded:
        # The fold's only sentence is about the shades, so it needs one drawn.
        caption = _chart_notes(_GRID_BINS_NOTE)
    legend = _heat_legend(
        grid,
        floor,
        _dominant_unit(entry, set(organs)),
        greyed=tissue if unsampled else None,
    )
    return f'<div class="chart-scroll">{figure}</div>{legend}{caption}'


def _hover_dots(
    name: str,
    placed: Sequence[tuple[int, float]],
    points: Sequence[tuple[float, float]],
    labels: Sequence[str],
    unit: str,
) -> str:
    """An invisible, easy-to-hit circle on every placed point, titled with what it is.

    The marks and lines are a pixel or two wide, which is too small to hover;
    these are a few pixels wider and carry the organ, the stage and the value,
    so a reader can name any point without leaving the chart. They draw
    nothing (`fill: none`, `pointer-events: all`), so the picture is unchanged.
    """
    suffix = f" {html.escape(unit)}" if unit else ""
    return "".join(
        f'<circle class="chart-hit" cx="{coordinate(x)}" cy="{coordinate(y)}" r="5">'
        f"<title>{html.escape(name)} \u00b7 {html.escape(_short_stage(labels[index]))} \u00b7 "
        f"{html.escape(_fmt(value))}{suffix}</title></circle>"
        for (index, value), (x, y) in zip(placed, points, strict=True)
    )


def _sampled_cardiac(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> list[tuple[str, int]]:
    """Each cardiac organ this dataset sampled for this gene, and how often.

    Sorted, because `cardiac` is a `frozenset` and an unsorted iteration would
    reorder a published sentence between two builds of one commit.
    """
    return [
        (tissue, len(_tissue_medians(entry, tissue)))
        for tissue in sorted(cardiac)
        if _tissue_medians(entry, tissue)
    ]


def _sampled_organs_clause(sampled: Sequence[tuple[str, int]]) -> str:
    """ "heart (19 stages sampled)" -- the organ clause both refusals open with.

    One implementation, because the two sentences name the same organs for
    two different reasons and a reader who sees one never sees the other; two
    copies of this could spell the same fact two ways with nothing failing.

    **`stage` has a singular form here and did not before.** The count is a
    real 1 on a live corpus shape -- a dataset sampling one cardiac stage for
    a gene -- and "1 stages sampled" is the `_abundance`/`_count` discipline
    lapsing in the one place nothing was reading.
    """
    return ", ".join(
        f"{html.escape(tissue)} ({count} stage{'' if count == 1 else 's'} sampled)"
        for tissue, count in sampled
    )


def _unbanded_sentence(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> str:
    """D43's refusal, in words: placed measurements but no phase to band them with.

    The third of the three reasons there is no heart chart, and the only one that is not
    about the gene at all. `_phase_bands` returns nothing when no curated
    cardiac phase with a stated end covers any of this dataset's stages, and
    a curve drawn without one is the re-plot of the source's own browser that
    D43 forbids. That is a fact about this atlas's vocabulary, so this
    sentence claims nothing about the measurements -- every one of which is
    placed, and printed in the tables below.

    Reached only where something *is* placed. Where nothing is,
    `_no_trajectory_sentence` has the better sentence whatever the bands say:
    "no trajectory" is then a fact about the gene, and the missing band is a
    second reason a reader does not need.

    Returns `""` when this dataset sampled no cardiac organ for this gene,
    the same degenerate case `_no_trajectory_sentence` returns `""` for and
    for the same reason: there is no organ to make a statement about.
    """
    sampled = _sampled_cardiac(entry, cardiac)
    if not sampled:
        return ""
    return (
        f'<p class="notice-inline">No trajectory is drawn for '
        f"{_sampled_organs_clause(sampled)}: none of this dataset's developmental stages "
        "falls inside a curated cardiac phase this atlas can band, and an unbanded curve "
        "is the source's own plot rather than anything this atlas adds. Every per-stage "
        "figure is in the tables below.</p>"
    )


def _anything_placed(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> bool:
    """Does any cardiac organ carry a median this dataset placed?

    The one question that tells the second reason there is no heart chart from the third,
    and it is asked of every cardiac organ at once because the answer decides
    one sentence for the whole dataset block.
    """
    return any(
        point.value is not None
        for tissue in sorted(cardiac)
        for point in _tissue_medians(entry, tissue)
    )


def _no_trajectory_sentence(
    entry: DatasetProfileEntry, cardiac: frozenset[str], floor: float | None
) -> str:
    """D42: what a page says instead of drawing an empty chart.

    Measured 2026-08-20 against a real build: seven published genes are below
    the floor in every cardiac tissue at every sampled stage -- SEMA3E, ZIC3,
    USP44, DAW1, FGF8, GDF1, NODAL. Their headline grades are `definitive`
    (ZIC3, NODAL), `moderate` (DAW1, FGF8), `limited` (USP44) and none at all
    (SEMA3E, GDF1). They are laterality and outflow-tract genes acting in rare
    lineages, which is exactly the dilution `methods.html#expression` states
    in full. `_BULK_TISSUE_NOTICE`, the line above every dataset block
    (`_CELL_TYPE_NOTICE` where the gene has cell-type data, which says no more
    of it), has said nothing of dilution since 2026-10-26 (bulk tissue, and a low figure
    is not evidence the gene is unimportant, then the link), so on these
    pages this sentence is the only statement of the mechanism: measured on a
    build of 802d694, "diluted" appears on 7 gene pages, these seven, once
    each, all in this sentence. An empty chart beside a definitive chip
    asserts the thing that caveat denies, so this tier states the fact and
    the mechanism rather than showing a frame with nothing in it.

    **Two branches, because "nothing placed" is not one fact.** All seven of
    those genes are below the floor, and that is what earns the dilution
    argument; a cell left unplaced because no percentile grid was published,
    or because the dataset declares no floor, is a gap in the *reference*
    rather than a low reading, and saying "below the detection floor" of it
    would be false. The second branch names the recorded reason instead and
    claims nothing about the gene.

    Returns `""` when this dataset sampled no cardiac organ for this gene at
    all: there is then no organ to make either statement about, and every
    figure it did publish is in the tables below regardless.

    **Reached only where nothing at all is placed, and the first branch is
    false the moment that stops holding.** "At every stage sampled" is read
    off the set of recorded gap *reasons*, which says what the gaps were and
    not how many stages had one; a series placed at 4 and 6 wpc and below the
    floor at 5 has exactly one reason in that set. It published the bold
    claim of a gene whose page printed two percentiles three lines below it
    (2026-08-20). `_dataset_block` is what holds the precondition, by asking
    `_anything_placed` before choosing this sentence over
    `_unbanded_sentence`; the guard is stated here because this is where a
    later author would have to break it.
    """
    sampled = _sampled_cardiac(entry, cardiac)
    if not sampled:
        return ""
    organs = _sampled_organs_clause(sampled)
    named = {tissue for tissue, _ in sampled}
    reasons = {
        measured["not_placed_reason"]
        for stage in entry["stages"]
        for measured in stage["tissues"]
        if measured["tissue"] in named and measured["placement"] is None
    }
    if reasons != {ProfileGap.BELOW_DETECTION_FLOOR.value}:
        clauses = "; ".join(
            _PLACEMENT_GAP_CLAUSE.get(reason or "", "no reason was recorded for this gap")
            for reason in sorted(reason or "" for reason in reasons)
        )
        return (
            f'<p class="notice-inline">No trajectory is drawn for {organs}: no measurement '
            f"there is placed against this dataset's percentile grid &mdash; {clauses}. Every "
            "per-stage figure is in the tables below.</p>"
        )
    # The floor's own value is a separate sentence rather than a parenthesis
    # inside the first, because the organ clause already ends in one -- "in
    # heart (19 stages sampled) (detection floor 1 tpm)" is what nesting
    # them produced, and a reader stops reading at the second bracket.
    qualifier = (
        f" This dataset's detection floor is {_fmt(floor)} "
        f"{html.escape(_dominant_unit(entry, named))}."
        if floor is not None
        else ""
    )
    return (
        f'<p class="notice-inline">This gene reads <strong>below this dataset\'s detection '
        f"floor at every stage sampled</strong> in {organs}, so there is no trajectory to "
        f"draw.{qualifier} That is "
        "<strong>not evidence that</strong> the gene is unimportant to heart development, and it "
        "is not a measured zero: every figure here is a median over bulk tissue, and a gene "
        "transcribed intensely in one rare cardiac lineage still reads below a bulk floor once "
        "diluted across the sample. Every per-stage figure is in the tables below.</p>"
    )


def _dominant_unit(entry: DatasetProfileEntry, tissues: set[str]) -> str:
    """The unit the named organs' medians are quoted in.

    A floor quoted without one is the `_count`/`count_unit` defect in a new
    place; the organs here all come from one dataset, so one unit answers for
    all of them, and the first is read rather than a set assembled.
    """
    for tissue in sorted(tissues):
        unit = _tissue_unit(entry, tissue)
        if unit:
            return unit
    return ""


def _sampled_tissues(entry: DatasetProfileEntry) -> list[str]:
    """Every organ this dataset sampled for this gene, sorted.

    Sorted because it is derived from a `set`: `sort_keys` in the JSON encoder
    has nothing to say about the order the grid's rows are laid out in, and two
    builds of one commit must emit the same bytes.
    """
    return sorted({cell["tissue"] for stage in entry["stages"] for cell in stage["tissues"]})


def _dataset_block(
    entry: DatasetProfileEntry, dataset: Dataset | None, phases: CardiacPhaseFile | None
) -> str:
    """One dataset's whole contribution: its figure, then every figure in tables.

    The figure (`_bulk_figure`) answers one question at a glance -- when the
    cardiac tissue reads, how high, and whether it is where the gene reads
    highest -- and it is drawn for every dataset that sampled any organ. Where
    the cardiac tissue has no heart chart to draw it returns the grid alone,
    and this keeps the lede sentence (`_unbanded_sentence`,
    `_no_trajectory_sentence`) in front of it: the chart's absence is a fact a
    reader is owed in words (D42).

    The link is D39(b)'s other half reaching a reader rather than only a
    program: `build_profile_quantiles` publishes the 101-point grid a
    percentile above was read against, and this is the one place a person
    reading the page -- not a script reading the bundle -- can reach it.
    Omitted when `quantile_shard` is `None` (a corpus mid-curation, or the
    gate bypassed) rather than linking a file `Emitter` never wrote.

    `dataset` is looked up by the caller and may be `None` --
    `ProfileGap.DATASET_NOT_REGISTERED` is a real, named state, not a caller
    error, so this degrades to no cardiac tissues and no floor rather than
    raising: the gap is already stated on every cell it affects, and this
    function's job is to render what the payload says, not to police a
    registration mismatch a validator already reports elsewhere.

    **The stage tables move into a `<details>`; they do not leave.** This
    section rendered 21 near-identical blocks per dataset, one per stage, and
    59% of its text was verbatim repetition -- "below the detection floor in
    whole liver at this stage (detection floor 1 tpm)" appeared 20 times on
    the TBX5 page.

    **A bare percentage is not reproducible, so here is the method that
    produced that one.** Build `dev`, take the TBX5 page from
    `<h2>Developmental expression</h2>` to the end, replace every tag with a
    newline (which puts each table cell on its own line), strip and drop
    blanks, and count the words in lines that duplicate an earlier line:
    1,854 of 3,136 words over 513 lines, 146 of them distinct -- 59.1%. Other
    ways of asking the same page give very different answers, none of them
    wrong: 12% of its characters sit inside a repeated *sentence*, and 89% of
    its tokens sit inside a repeated 4-gram. The figure is about *lines*, and
    the lines are mostly table cells.

    A chart that *replaced* them would take every exact figure out of the
    HTML, which is this repository's characteristic defect (curated work
    reaching no page) wearing a redesign. The summary is what a reader sees
    first; the record is one click away and still complete.
    """
    heading = f"<h3>{html.escape(entry['dataset'])}"
    shard = entry["quantile_shard"]
    if shard:
        heading += f' <a href="../{html.escape(shard)}">percentile grid as JSON</a>'
    heading += "</h3>"
    cardiac = frozenset(dataset.cardiac_tissues) if dataset is not None else frozenset()
    floor = dataset.detection_floor if dataset is not None else None
    # **Three refusals, three sentences.** There is no heart chart for three
    # different reasons -- nothing sampled, nothing placed, and no curated
    # phase band -- and a caller that routes them all into one sentence
    # publishes a claim that is false for two of them. Measured 2026-08-20
    # with `runner.py` passing `phases=None`: 55 pages said "no measurement
    # there is placed against this dataset's percentile grid &mdash; ." with
    # the reason list empty, and 30 said, in bold, that the gene reads below
    # the detection floor at every stage sampled -- one of them with 18 of its
    # 19 heart stages placed and printed in the tables below. Nothing was
    # published in either state, because `runner.py` passes the real
    # vocabulary; the sentences were wrong all the same, and a dataset whose
    # stages fall in no *ended* cardiac phase reaches them.
    #
    # Placement first, because it is the fact about the gene. A missing band
    # is a fact about this atlas's own vocabulary, and where nothing is
    # placed a reader is better served by the dilution caveat than by it.
    if _charted_tissue(entry, cardiac, phases) is not None:
        lede = ""
    elif _anything_placed(entry, cardiac):
        lede = _unbanded_sentence(entry, cardiac)
    else:
        lede = _no_trajectory_sentence(entry, cardiac, floor)
    glance = "".join(_glance_line(entry, tissue, floor) for tissue in sorted(cardiac))
    figure = _bulk_figure(entry, dataset, phases)
    stages = _stage_tables(entry, cardiac, floor)
    total = len(entry["stages"])
    noun = "stage" if total == 1 else "stages"
    folded = (
        f'<details class="stage-figures"><summary>every figure, all {total} {noun}'
        f"</summary>{stages}</details>"
    )
    return heading + glance + lede + figure + folded


def _expression_section(
    profile: ExpressionProfile,
    datasets: Mapping[str, Dataset],
    phases: CardiacPhaseFile | None,
    cell_types: CellTypeExpression = EMPTY_CELL_TYPE_EXPRESSION,
    *,
    gene: str,
) -> str:
    """Developmental expression: percentile, tau and phase, per organ and stage.

    **Renders a section even for a gene with no data, unlike `_burden_
    section`.** A burden table's absence would misstate a null result as
    evidence against a gene it was never tested for, so that section returns
    `""`. There is no equivalent risk here: `datasets == []` states a fact this
    atlas already knows -- no curated dataset mentions this gene -- and
    omitting the section would make that indistinguishable from the atlas
    never having built this feature at all. `_NOT_CURATED_EXPRESSION` follows
    the idiom of the lesion-assertion notice (now the card's curation sentence,
    `_card_curation`), not `_burden_section`'s.

    `datasets` resolves each entry's own accession to the curated `Dataset` it
    came from, for `cardiac_tissues` (the argmax gate) and `detection_floor`
    (rule 5's required value) -- both dataset-level facts this per-gene
    payload does not repeat. A gene absent from the caller's own `profiles`
    mapping is the caller's job to default to `EMPTY_EXPRESSION_PROFILE`,
    matching `bundles._expression_profile`'s own fallback, so the page and the
    bundle cannot disagree about which genes carry data.

    **The reading notes left on 2026-10-26 (D50).** They were folded under
    "How to read these figures" on all 119 gene pages -- how a median is
    placed and why none is placed below the floor, and, where the section laid
    out more than one placement (118 pages, measured on a build of 0fc451a),
    that percentiles are not comparable across organs or stages. Both teach
    how to read the figures and say nothing about this gene, so both are on
    `methods.html#expression`, where the conditional went too: stated once in
    general, the comparability sentence has no position to be false in. What
    stays is `_BULK_TISSUE_NOTICE`'s short form, visible and unconditional,
    which links there.

    **Cell types open the section.** A gene with a non-empty
    `cell_types["datasets"]` gets `_CELL_TYPE_NOTICE` instead of
    `_BULK_TISSUE_NOTICE` and a `cell_type_block` per dataset before the bulk
    blocks; a gene without keeps the bulk notice and no block. Each cell-type
    entry's dataset must be in `datasets` (`KeyError` otherwise): the block
    draws the dataset's own declared cell types and regions, and a missing
    record would leave nothing to draw them from. A gene with cell-type data
    but no bulk profile gets the block and `_CELL_TYPE_ONLY_NOTICE` (no
    timeline is claimed) in place of `_NOT_CURATED_EXPRESSION`, whose "not yet
    curated" would contradict it. `gene` is the symbol (or id) the block's
    figures are named by for a screen reader.
    """
    entries = profile["datasets"]
    cell_entries = cell_types["datasets"]
    if not entries and not cell_entries:
        return "<h2>Developmental expression</h2>" + _NOT_CURATED_EXPRESSION

    blocks = "".join(
        _dataset_block(entry, datasets.get(entry["dataset"]), phases) for entry in entries
    )
    if not cell_entries:
        return "<h2>Developmental expression</h2>" + _BULK_TISSUE_NOTICE + blocks
    cells = "".join(
        cell_type_block(entry, datasets[entry["dataset"]], gene) for entry in cell_entries
    )
    notice = _CELL_TYPE_NOTICE if entries else _CELL_TYPE_ONLY_NOTICE
    return "<h2>Developmental expression</h2>" + notice + cells + blocks


# P7 (plan decision): a gene with no literature file says so in three words,
# which keeps the page under the gene-page budget test's six-word boilerplate
# threshold on all 119 pages. Rendered unconditionally for such a gene, for the
# reason `_NOT_CURATED_EXPRESSION` is: omitting the section would make "the
# atlas has not searched for this gene" read the same as "the atlas has not
# built this feature". It is not a claim that the literature holds nothing.
_NOT_CURATED_LITERATURE: Final = "<p>Not yet curated.</p>"

# What the exclusion tally calls each reason (`vocab.ExclusionReason`), in
# that enum's order. A reason with no entry renders its raw token, so a reason
# added to the vocabulary is shown rather than dropped.
_EXCLUSION_LABELS: Final = (
    ("not_peer_reviewed", "not peer reviewed"),
    ("review_article", "review article"),
    ("not_human_genetics", "not human genetics"),
    ("no_finding_in_gene", "no finding in this gene"),
    ("not_chd", "not about congenital heart disease"),
    ("same_data_as", "same data as another paper"),
    ("full_text_not_accessible", "full text not accessible"),
    ("other", "other"),
)

_FINDING_HEADERS: Final = (
    "year",
    "paper",
    "design",
    "methods",
    "counts",
    "carriers by origin",
    "result",
    "lesions",
    "independence",
    "finding",
    "approved by",
)

# Tokens whose words are not the token with its underscores removed: an acronym
# keeps its capitals (`karyotype_or_fish` read "karyotype or fish").
_TOKEN_LABELS: Final = {"karyotype_or_fish": "karyotype or FISH"}


def _words(token: str) -> str:
    """A vocabulary token as words: `family_segregation` -> `family segregation`."""
    return _TOKEN_LABELS.get(token, token.replace("_", " "))


def _number(n: int, singular: str, plural: str) -> str:
    """`1 proband`, `2 probands`: the count with its noun in agreement."""
    return f"{n:,} {singular if n == 1 else plural}"


# What one counted person or family is called, by `count_unit`: the singular of
# the unit, and the carrier noun in that unit (R30: with `families`, each origin
# count is carrier families, so "carriers" would be a different quantity).
_UNIT_NOUNS: Final = {
    "probands": ("proband", "probands", "carrier", "carriers"),
    "families": ("family", "families", "carrier family", "carrier families"),
}


def _finding_counts(finding: Mapping[str, Json]) -> str:
    """A finding's counts, always with what was counted (the `count_unit` rule).

    Raises `ValueError` for a finding with neither shape of count: the model
    refuses one, so reaching it means a bypassed gate, and a blank cell would
    pass for a finding with no denominator.
    """
    screened, unit, carriers = finding["screened"], finding["count_unit"], finding["carriers"]
    if isinstance(screened, int) and isinstance(unit, str) and isinstance(carriers, dict):
        total = sum(v for v in carriers.values() if isinstance(v, int))
        # An unknown unit renders its raw token rather than nothing (the
        # `count_unit` rule), in the plural for both nouns.
        one, many, carrier_one, carrier_many = _UNIT_NOUNS.get(
            unit, (unit, unit, "carrier", "carriers")
        )
        return f"{_number(total, carrier_one, carrier_many)} among {_number(screened, one, many)}"
    cases, controls = finding["cases"], finding["controls"]
    if isinstance(cases, int) and isinstance(controls, int):
        return f"{_number(cases, 'case', 'cases')}, {_number(controls, 'control', 'controls')}"
    raise ValueError(f"{finding['id']}: a published finding has no counts")


def _finding_carriers(finding: Mapping[str, Json]) -> str:
    carriers = finding["carriers"]
    if not isinstance(carriers, dict):
        return "\u2014"
    return (
        f"{carriers['de_novo']} de novo, {carriers['inherited']} inherited, "
        f"{carriers['origin_unknown']} origin unknown"
    )


def _finding_result(finding: Mapping[str, Json]) -> str:
    """What the finding measured against, and the statistics exactly as stored.

    A rare-variant association names its comparator (a control cohort and a
    mutation model are different questions); `p` and `OR` print through `repr`,
    unrounded, so 1.23456789e-09 is not shown as 0. `\u2014` where the finding
    states none of them (a case report, a cohort screen).
    """
    parts = []
    comparator = finding["comparator"]
    if isinstance(comparator, str):
        parts.append(f"compared with {_words(comparator)}")
    p_value, odds_ratio = finding["p_value"], finding["odds_ratio"]
    if p_value is not None:
        parts.append(f"p {p_value!r}")
    if odds_ratio is not None:
        parts.append(f"OR {odds_ratio!r}")
    return ", ".join(parts) or "\u2014"


def _finding_independence(finding: Mapping[str, Json]) -> str:
    """Which collections the finding drew on and whether its probands are shared (R31, R37).

    `unknown` is a value the curator wrote and `[]` is an assertion (the paper
    states its cohort shares no probands with any other finding), so they are
    told apart rather than both rendering as blank. A named id may belong to a
    finding that is not yet approved and so not on any page.
    """
    parts = []
    cohorts = _as_list(finding["cohorts"])
    if cohorts:
        parts.append("cohorts: " + ", ".join(str(c) for c in cohorts))
    shared = finding["shares_probands_with"]
    if isinstance(shared, list):
        parts.append(
            "probands shared with: " + ", ".join(str(i) for i in shared)
            if shared
            else "probands: none shared"
        )
    else:
        parts.append(f"probands: {shared}")
    return "; ".join(parts)


def _literature_searches(searches: Sequence[SearchRecord]) -> str:
    """Every recorded search, folded: its date, its arithmetic and the exact query.

    The query is the published statement of what was searched; it is escaped and
    shown in `<code>` unaltered. Nothing here is visible until opened.
    """
    items = "".join(
        f"<li>{html.escape(str(search['run_on']))}: "
        f"{_number(int(str(search['hits'])), 'hit', 'hits')}, {search['new_hits']} new; "
        f"<code>{html.escape(str(search['query']))}</code></li>"
        for search in searches
    )
    return f"<details><summary>Searches run</summary><ul>{items}</ul></details>"


def _literature_exclusions(excluded: Mapping[str, int]) -> str:
    """The exclusion tally by reason, folded; `""` where nothing was excluded."""
    if not excluded:
        return ""
    known = [token for token, _ in _EXCLUSION_LABELS]
    labels = dict(_EXCLUSION_LABELS)
    order = [t for t in known if t in excluded] + sorted(t for t in excluded if t not in labels)
    items = "".join(
        f"<li>{html.escape(labels.get(token, token))}: {excluded[token]}</li>" for token in order
    )
    return f"<details><summary>Excluded papers, by reason</summary><ul>{items}</ul></details>"


# One mark per paper (owner, 2026-11-03). A paper with one design keeps that design's
# shape; a paper with several designs is a ring, so the ring is never a design. Every
# member of FindingDesign must have a shape: `_paper_shape` refuses a token it does not
# know rather than drawing a default.
_DESIGN_SHAPES: Final[dict[str, Shape]] = {
    "case_report": "circle",
    "family_segregation": "square",
    "cohort_screen": "triangle_up",
    "rare_variant_association": "triangle_down",
    "common_variant_association": "star",
}
_SEVERAL_DESIGNS: Final[Shape] = "ring"
_DESIGN_ORDER: Final[dict[str, int]] = {d.value: i for i, d in enumerate(FindingDesign)}


class _PaperGroup(NamedTuple):
    """One publication's approved findings, with the distinct designs they report."""

    pmid: str
    publication: Publication
    findings: tuple[Mapping[str, Json], ...]
    designs: tuple[str, ...]


def _paper_groups(
    findings: Sequence[Mapping[str, Json]],
    publications: Mapping[str, Publication],
) -> list[_PaperGroup]:
    """The findings grouped by paper, ordered by year then numeric PMID.

    Grouping keys on the PMID, never on adjacency, so a paper's findings need not be
    contiguous in the input. The PMID sort is numeric (`PMID:10` after `PMID:2`). A
    finding whose paper is unregistered raises: validation refuses it (LIT011), so
    reaching here means the gate was bypassed.
    """
    by_pmid: dict[str, list[Mapping[str, Json]]] = {}
    for finding in findings:
        pmid = str(finding["publication"])
        if pmid not in publications:
            raise ValueError(f"finding {finding['id']} cites {pmid}, which is not registered")
        by_pmid.setdefault(pmid, []).append(finding)
    groups = []
    for pmid, members in by_pmid.items():
        designs = sorted(
            {str(f["design"]) for f in members},
            key=lambda d: (_DESIGN_ORDER.get(d, len(_DESIGN_ORDER)), d),
        )
        groups.append(_PaperGroup(pmid, publications[pmid], tuple(members), tuple(designs)))
    groups.sort(key=lambda g: (g.publication.year, pmid_number(g.pmid)))
    return groups


def _paper_shape(group: _PaperGroup) -> Shape:
    """The paper's mark: its design's shape, or the ring when it reports several designs."""
    for design in group.designs:
        if design not in _DESIGN_SHAPES:
            raise ValueError(f"no shape for finding design {design!r} ({group.pmid})")
    if len(group.designs) > 1:
        return _SEVERAL_DESIGNS
    return _DESIGN_SHAPES[group.designs[0]]


# Timeline geometry, in pixels. Stacking is 11 px: a glyph is a 10 px box, so two marks of
# one year never touch (the star runs from y-5 to about y+4.05 and fits).
_TL_PAD: Final = 24.0
_TL_STACK: Final = 11.0
_TL_TOP: Final = 10.0
_TL_KEY_ROW: Final = 16.0
_SEVERAL_DESIGNS_LABEL: Final = "several designs"


def _timeline_year_labels(first: int, last: int) -> list[int]:
    """Both ends, and each multiple of 10 between them at least 4 years from both."""
    if first == last:
        return [first]
    decades = [
        y for y in range(first + 1, last) if y % 10 == 0 and y - first >= 4 and last - y >= 4
    ]
    return [first, *decades, last]


def _own_lab_mark(publication: Publication) -> str:
    """`_OWN_LAB_MARK` for a paper an atlas author wrote, else the empty string.

    The timeline hover and the fold summary both carry it, so a reader who meets
    a paper on either surface meets the declaration (owner decision 5).
    """
    return _OWN_LAB_MARK if publication.own_lab else ""


def _literature_timeline(
    groups: Sequence[_PaperGroup],
    publications: Mapping[str, Publication],
) -> str:
    """One SVG: each approved paper a mark on a year axis, shaped by design.

    A paper's mark is `_paper_shape`; papers of one year stack upward, in group order,
    so none hides another. The hover names the paper and gives one clause per finding
    -- the counts of different findings are never added (T6). Raises `ValueError` on no
    papers: the caller renders the no-findings line instead of an empty axis.
    """
    if not groups:
        raise ValueError("a timeline needs at least one paper; there are no papers to draw")
    years = [g.publication.year for g in groups]
    first, last = min(years), max(years)
    span = _CHART_WIDTH - 2 * _TL_PAD

    def x_of(year: int) -> float:
        return (
            _CHART_WIDTH / 2 if first == last else _TL_PAD + (year - first) / (last - first) * span
        )

    tallest = max(years.count(year) for year in set(years))
    baseline = _TL_TOP + 5 + _TL_STACK * (tallest - 1)
    axis_y = baseline + 10
    seen: dict[int, int] = {}
    marks = []
    shapes_used = set()
    for group in groups:
        year = group.publication.year
        k = seen.get(year, 0)
        seen[year] = k + 1
        shape = _paper_shape(group)
        shapes_used.add(shape)
        clauses = "; ".join(
            f"{_words(str(f['design']))}: {_finding_counts(f)}" for f in group.findings
        )
        label = _study_label(group.pmid, publications) + _own_lab_mark(group.publication)
        hover = " · ".join([label, group.pmid, group.publication.title, clauses])
        marks.append(glyph(x_of(year), baseline - _TL_STACK * k, shape=shape, title=hover))
    body = (
        f'<line class="chart-axis" x1="{coordinate(_TL_PAD)}" y1="{coordinate(axis_y)}" '
        f'x2="{coordinate(_CHART_WIDTH - _TL_PAD)}" y2="{coordinate(axis_y)}"/>'
        + "".join(
            f'<text class="chart-tick-label" x="{coordinate(x_of(y))}" '
            f'y="{coordinate(axis_y + 12)}" text-anchor="middle">{y}</text>'
            for y in _timeline_year_labels(first, last)
        )
        + f'<g class="lit-marks">{"".join(marks)}</g>'
    )
    # The key: only the shapes drawn here, design order then the ring, wrapped to rows.
    entries = [
        (_DESIGN_SHAPES[d.value], _words(d.value))
        for d in FindingDesign
        if _DESIGN_SHAPES[d.value] in shapes_used
    ]
    if _SEVERAL_DESIGNS in shapes_used:
        entries.append((_SEVERAL_DESIGNS, _SEVERAL_DESIGNS_LABEL))
    key = []
    x, row = _TL_PAD, 0
    for shape, label in entries:
        width = 14 + len(label) * _CHAR_WIDTH + 18
        if x > _TL_PAD and x + width > _CHART_WIDTH - _TL_PAD:
            x, row = _TL_PAD, row + 1
        y = axis_y + 30 + _TL_KEY_ROW * row
        key.append(glyph(x + 5, y, shape=shape, title=label))
        key.append(
            f'<text class="chart-label" x="{coordinate(x + 14)}" y="{coordinate(y + 3)}">'
            f"{html.escape(label)}</text>"
        )
        x += width
    height = int(math.ceil(axis_y + 30 + _TL_KEY_ROW * row + 10))
    n = len(groups)
    when = f"in {first}" if first == last else f"{first} to {last}"
    title = (
        f"{n} {'publication' if n == 1 else 'publications'} with approved findings, "
        f"{when}, by year and design"
    )
    return svg_figure(
        width=_CHART_WIDTH,
        height=height,
        title=title,
        body=body + f'<g class="lit-key">{"".join(key)}</g>',
    )


_PAPER_FINDING_HEADERS: Final = _FINDING_HEADERS[2:]


def _finding_row(finding: Mapping[str, Json]) -> Row:
    """A finding's cells from `design` to `approved by`; the caller adds neither year nor paper.

    The review mark (R24) sits beside the summary, so it qualifies the finding
    and never hides it.
    """
    status, note = finding["review_status"], finding["review_note"]
    mark = review_status_note(
        str(finding["id"]),
        ReviewStatus(status) if isinstance(status, str) else None,
        note if isinstance(note, str) else None,
    )
    return Row(
        (
            _words(str(finding["design"])),
            ", ".join(_words(str(m)) for m in _as_list(finding["methods"])),
            _finding_counts(finding),
            _finding_carriers(finding),
            _finding_result(finding),
            ", ".join(_words(str(g)) for g in _as_list(finding["lesion_groups"])),
            _finding_independence(finding),
            Markup(html.escape(str(finding["summary"])) + mark),
            str(finding["approved_by"]),
        )
    )


def _literature_publications(
    groups: Sequence[_PaperGroup],
    publications: Mapping[str, Publication],
) -> str:
    """The approved papers folded, each opening to its own findings.

    The PMID link is in the opened block, not in the `<summary>`, which holds no
    interactive content. A paper's findings are listed, never summed (T6).
    """
    n = len(groups)
    blocks = []
    for group in groups:
        k = len(group.findings)
        designs = " + ".join(_words(d) for d in group.designs)
        summary = (
            f"{_study_label(group.pmid, publications)} \u00b7 {group.publication.journal} "
            f"\u00b7 {designs} \u00b7 {k} {'finding' if k == 1 else 'findings'}"
            f"{_own_lab_mark(group.publication)}"
        )
        table = data_table(_PAPER_FINDING_HEADERS, [_finding_row(f) for f in group.findings])
        blocks.append(
            f'<details class="lit-paper"><summary>{html.escape(summary)}</summary>'
            f"<p>{html.escape(group.publication.title)} "
            f'<a href="{html.escape(_pubmed(group.pmid))}">{html.escape(group.pmid)}</a></p>'
            f'<div class="scroll">{table}</div></details>'
        )
    noun = "publication" if n == 1 else "publications"
    return f"<details><summary>{n} {noun}</summary>{''.join(blocks)}</details>"


def _literature_section(
    literature: GeneLiterature,
    publications: Mapping[str, Publication],
) -> str:
    """The gene's curated literature, in this order: the summary line, the timeline of
    approved papers and their folded list (or `No approved findings yet.` and no axis),
    the exclusion tally, every search.

    No search recorded: `_NOT_CURATED_LITERATURE` and nothing else. Otherwise a
    visible line counting what the searches returned and what was decided (one
    search: `searched D`; several: `N searches, latest D` and `in all`) -- the
    papers, included and excluded counts always print, so `0 included` or
    `0 excluded` can render; only the optional clauses (`, N not yet triaged`,
    `; N findings await approval`) are omitted when their count is 0 -- then
    the exclusion tally and every search folded with its exact query. The
    timeline and the list cover approved findings only, which is all the bundle
    carries; the `await approval` clause is what stops the unapproved ones being
    invisible. A finding's `review_status` and `review_note` are rendered beside
    its summary (ruling R24): the mark qualifies a published finding and never hides it.

    `publications` supplies each paper's year; `_paper_groups` raises on an
    unregistered paper, because validation refuses one (LIT011).
    """
    searches = literature["searches"]
    if not searches:
        return "<h2>Literature</h2>" + _NOT_CURATED_LITERATURE
    triage = literature["triage"]
    papers = sum(search["new_hits"] for search in searches)
    excluded_total = sum(triage["excluded"].values())
    latest = html.escape(max(search["run_on"] for search in searches))
    noun = "paper" if papers == 1 else "papers"
    if len(searches) == 1:
        opening = f"PubMed, searched {latest}: {papers} {noun}; "
    else:
        opening = f"PubMed, {len(searches)} searches, latest {latest}: {papers} {noun} in all; "
    line = f"{opening}{triage['included']} included, {excluded_total} excluded"
    if triage["pending"]:
        line += f", {triage['pending']} not yet triaged"
    waiting = literature["pending_approval"]
    if waiting:
        line += f"; {waiting} {'finding awaits' if waiting == 1 else 'findings await'} approval"
    findings = literature["findings"]
    if findings:
        groups = _paper_groups(findings, publications)
        body = _literature_timeline(groups, publications) + _literature_publications(
            groups, publications
        )
    else:
        body = "<p>No approved findings yet.</p>"
    return (
        f'<h2>Literature</h2><p class="method">{line}.</p>'
        + body
        + _literature_exclusions(triage["excluded"])
        + _literature_searches(searches)
    )


def _as_list(value: Json) -> list[Json]:
    return value if isinstance(value, list) else []


def build_gene_pages(
    facts: Mapping[str, GeneFacts],
    emitter: Emitter,
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    assertions: Mapping[str, list[LesionAssertion]],
    publications: Mapping[str, Publication],
    burden: Mapping[str, Sequence[BurdenRow]],
    cohorts: Mapping[str, Cohort],
    families: tuple[frozenset[str], ...] = (),
    axes: tuple[tuple[str, str], ...] = (),
    profiles: Mapping[str, ExpressionProfile] | None = None,
    datasets: Mapping[str, Dataset] | None = None,
    phases: CardiacPhaseFile | None = None,
    combined: Mapping[str, Sequence[CombinedRow]] | None = None,
    scope_warrants: Mapping[str, Warrant] | None = None,
    scope_labels: Mapping[str, str] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
    scope_terms: Mapping[str, ScopeEntry] | None = None,
    cell_types: Mapping[str, CellTypeExpression],
    literature: Mapping[str, GeneLiterature],
) -> None:
    """Emit one HTML page per gene in `facts`.

    `facts` is `build_genes`' return, so the pages render from the same
    derivation the bundles were written from rather than deriving it twice --
    a page and its bundle disagreeing about a gene's confidence is the failure
    this shared input removes.

    `symbols` is the same `{gene: labels.symbol}` projection `build_genes` and
    `build_landing` take; a gene absent from it reads as its own HGNC id, the
    same registry fallback those two apply.

    `identities` is `identity.read_gene_identities`' return (D52). Given a mapping,
    every gene must carry the ids its required links are keyed on: a gene
    absent from it gets `EMPTY_IDENTITY`, which `external_links` refuses. `None`
    renders no identity block at all -- the unit suites that do not exercise
    identity pass nothing -- and `build_site` always passes the registry.

    A gene absent from `validity` gets `uncurated()` -- the section under its own
    heading, each authority group stating it holds no record, never a missing
    section.

    `profiles` and `datasets` both default to `None` -- and to `{}` inside --
    rather than to `{}` in the signature, the same reason `concordance` on
    `build_gene_index_page` does: an empty *literal* dict default is one object
    shared across every call that omits the argument, which is harmless only as
    long as nobody ever mutates it, and this project does not rely on that
    holding forever. A gene absent from `profiles` gets `EMPTY_EXPRESSION_
    PROFILE`, matching `bundles._expression_profile`'s own fallback exactly, so
    the page and the bundle cannot disagree about which genes carry a
    developmental expression profile. `datasets` resolves an entry's own
    accession to the curated record `_expression_section` reads
    `cardiac_tissues` and `detection_floor` from; see that function's docstring
    for why a missing dataset degrades rather than raises.

    `phases` is the same `CardiacPhaseFile` the profiles were derived against,
    threaded down rather than re-read, so a band on a chart and the phase
    sentence under the stage it covers cannot name different vocabularies.
    `None` -- no curated phase file at all -- costs the heart chart, not the
    section: `_charted_tissue` refuses to chart an unbanded curve, because a curve
    with no band is the panel the source's own browser already publishes, and
    every figure it would have summarised is in the tables regardless.

    `cell_types` is `gene_cell_type_expression`'s return; a gene absent from it
    gets `EMPTY_CELL_TYPE_EXPRESSION` and keeps the bulk-only section. Required,
    with no default, for the reason `build_genes` gives: with one, deleting it
    from the runner's call passed the whole suite (review, 2026-10-28). A
    caller with no cell-type data passes `{}`.

    `literature` is `gene_literature`'s return; a gene absent from it gets
    `EMPTY_GENE_LITERATURE` and the "Not yet curated" line. Required, with no
    default, for the reason `cell_types` is: the committed corpus carries no
    literature file, so a real build cannot tell a dropped mapping from an
    empty one. A caller with no literature passes `{}`.

    Sorted, like every loop in this build that iterates a mapping: `sort_keys`
    orders dict keys in a JSON payload and has nothing to say about the order
    files are written in, and a differing write order would move nothing here
    today but is the habit this project keeps.
    """
    profile_by_gene = profiles or {}
    dataset_registry = datasets or {}
    combined_by_gene = combined or {}
    warrants = scope_warrants or {}
    labels = scope_labels or {}
    graders = _panel_graders(validity, symbols)
    for gene in sorted(facts):
        fact = facts[gene]
        symbol = symbols.get(gene, gene)
        records = assertions.get(gene, [])
        gene_validity = validity.get(gene, uncurated())
        identity = None if identities is None else identities.get(gene, EMPTY_IDENTITY)
        gene_literature_bundle = literature.get(gene, EMPTY_GENE_LITERATURE)
        body = (
            '<div class="layout">'
            + _rail(
                gene,
                symbol,
                fact,
                gene_validity,
                burden.get(gene, ()),
                families=families,
                publications=publications,
                warrants=warrants,
                identity=identity,
                scope_terms=scope_terms,
                panel_graders=graders,
                literature=gene_literature_bundle,
            )
            + "<div>"
            + _SCOPE_SHORT
            # The Orphanet notes sit above both groups, not under either: a
            # scope warrant concerns the disease, which may be graded in the
            # panel group, the GenCC group, or both (measured: for 14 of 27
            # pages' notes the disease has no GenCC row), and a note placed
            # under one group read as qualifying only that group's rows.
            + _validity_section(
                gene_validity,
                notes=_orphanet_notes(gene_validity, warrants, labels),
                symbol=symbols.get(gene),
            )
            # The atlas's own curated evidence, for a curated gene. An uncurated
            # gene's page says so on the card (`_card_curation`, D55), not here:
            # the main-column notice that sat in this slot until 2026-10-26 is
            # the sentence the card now carries, and a page stating it twice is
            # a page whose two copies can drift apart.
            + (
                _evidence_section(records, publications)
                if fact.atlas_curation is AtlasCuration.CURATED
                else ""
            )
            # The atlas's literature beside its lesion assertions (owner,
            # 2026-11-09), on every page, "Not yet curated." included: both are
            # the atlas's own curation, so they sit together after the mirrored
            # validity they are distinguished from. It was the page's last
            # section until then, appended after expression.
            + _literature_section(gene_literature_bundle, publications)
            # After the atlas's own section, not between it and the mirrored
            # validity table, which keeps the atlas's own curation beside the
            # classifications it is distinguished from. The old main-column
            # notice said "each classification above is a ClinGen expert
            # panel's or a named curator's", and a burden table above it would
            # have put rows of statistics inside what it referred to.
            + _burden_section(
                burden.get(gene, ()),
                publications,
                cohorts,
                families,
                axes,
                combined_by_gene.get(gene, ()),
            )
            # Last: burden and expression are both mirrored/derived evidence
            # layers rather than the atlas's own curation, and expression is
            # the newer of the two.
            + _expression_section(
                profile_by_gene.get(gene, EMPTY_EXPRESSION_PROFILE),
                dataset_registry,
                phases,
                cell_types.get(gene, EMPTY_CELL_TYPE_EXPRESSION),
                gene=symbol,
            )
            + "</div></div>"
        )
        emitter.write_text(
            gene_page_path(HgncId(gene)),
            document(title=f"{symbol} — CHD Atlas", root="../", body=body),
        )


_BROWSE_HEADERS: Final = (
    "gene",
    "symbol",
    "confidence",
    "graded for",
    "rare variants across studies",
    "heart expression",
)


def render_combined_table(
    facts: Mapping[str, GeneFacts],
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    burden_counts: Mapping[str, int],
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
    profiles: Mapping[str, ExpressionProfile] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
    prefix: str = "../",
) -> tuple[str, str]:
    """Render (filters_html, table_html) for the published gene population."""
    rows: list[Row] = []
    for gene in sorted(facts):
        fact = facts[gene]
        symbol = symbols.get(gene, gene)
        name = (
            identities[gene].name
            if identities and gene in identities and identities[gene].name
            else ""
        )

        conf_val = fact.headline_confidence.value if fact.headline_confidence else ""
        confidence_cell = (
            Markup(
                f"{html.escape(conf_val or _EM_DASH)} " + chip("conflicting evidence", kind="warn")
            )
            if fact.has_conflicting_evidence
            else conf_val or _EM_DASH
        )

        diseases = _graded_diseases(validity.get(gene, uncurated()), fact.headline_confidence)
        disease_text = "; ".join(diseases) or _EM_DASH

        groups = " ".join(group.value for group in fact.lesion_groups)

        n_tests = burden_counts.get(gene, 0)
        tests_label = f"{n_tests} test{'s' if n_tests != 1 else ''}" if n_tests else "untested"
        dot_html = _dot_strip(_concordance_for(concordance, gene))
        cohort_content = f'<div>{dot_html}</div><div class="cohort-tests">{tests_label}</div>'

        heart_cell, _, _ = _heart_expression_cell(profiles.get(gene) if profiles else None)

        rows.append(
            Row(
                cells=(
                    Link(text=gene, href=f"{prefix}{gene_page_path(HgncId(gene))}"),
                    Link(text=symbol, href=f"{prefix}{gene_page_path(HgncId(gene))}"),
                    confidence_cell,
                    disease_text,
                    Markup(cohort_content),
                    heart_cell,
                ),
                attributes=(
                    ("search", f"{gene} {symbol} {name} {disease_text} {groups}".lower()),
                    ("lesion", groups),
                    ("confidence", conf_val),
                    ("validity", fact.validity_state.value),
                    ("curation", fact.atlas_curation.value),
                    ("burden", "yes" if n_tests else "no"),
                ),
            )
        )

    selects = "".join(
        f'<select name="{name}" aria-label="Filter by {aria}">'
        f'<option value="">{label}</option>'
        + "".join(
            f'<option value="{html.escape(value)}">{html.escape(value)}</option>'
            for value in sorted(values)
        )
        + "</select>"
        for name, label, aria, values in (
            (
                "lesion",
                "any lesion",
                "lesion group",
                {g.value for f in facts.values() for g in f.lesion_groups},
            ),
            (
                "confidence",
                "any confidence",
                "mirrored confidence",
                {f.headline_confidence.value for f in facts.values() if f.headline_confidence},
            ),
            (
                "validity",
                "any source",
                "validity state",
                {f.validity_state.value for f in facts.values()},
            ),
            (
                "curation",
                "any curation",
                "atlas curation",
                {f.atlas_curation.value for f in facts.values()},
            ),
            (
                "burden",
                "any rare-variant result",
                "whether rare-variant evidence exists",
                {"yes" if burden_counts.get(gene) else "no" for gene in facts},
            ),
        )
    )

    filters = (
        '<form id="filters" class="filters">'
        '<input name="q" type="search" aria-label="Search by gene symbol or HGNC id" '
        'placeholder="symbol or HGNC id">'
        f"{selects}</form>"
    )
    table = data_table(_BROWSE_HEADERS, rows, table_id="gene-table")
    return filters, table


def render_browse_table(
    facts: Mapping[str, GeneFacts],
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    burden_counts: Mapping[str, int],
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
    profiles: Mapping[str, ExpressionProfile] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
    prefix: str = "../",
) -> tuple[str, str]:
    return render_combined_table(
        facts,
        symbols=symbols,
        validity=validity,
        burden_counts=burden_counts,
        concordance=concordance,
        profiles=profiles,
        identities=identities,
        prefix=prefix,
    )


def build_gene_index_page(
    facts: Mapping[str, GeneFacts],
    emitter: Emitter,
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    burden_counts: Mapping[str, int],
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
    publications: Mapping[str, Publication] | None = None,
    profiles: Mapping[str, ExpressionProfile] | None = None,
    identities: Mapping[str, GeneIdentity] | None = None,
) -> None:
    """Emit `genes/index.html`: every published gene, filterable in the browser.

    `burden_counts` is `{gene: len(rows)}`, the same number `genes/index.json`
    publishes as `burden_row_count`. The browse page did not mention burden at
    all -- no column, no filter, no count -- while the payload behind it did, so
    a reader scanning 23 identical `definitive` rows could not tell that four of
    them (ISL1, NR2F2, RBM10, SMAD2) carry no loss-of-function evidence at all.
    Raised by review 2026-08-05 and measured on the built site.

    Every row is rendered here, by the build. The inline script only *hides*
    rows, so the page is complete to `curl`, to a crawler, and to a reader with
    JavaScript disabled -- the alternative, an empty `<tbody>` populated from
    `genes/index.json` at runtime, looks identical in a browser and serves
    nothing to any of the three.

    Sorted by HGNC id, the order `genes/index.json` publishes in, so the browse
    page and the payload behind it agree on which gene comes first without
    either reading the other.

    **The attribution line is not decoration.** This page renders 23 rows
    reading `definitive` under a column headed `confidence` (measured
    2026-08-04), and it was the only page on the site where that could be read
    as the atlas's own call: `_validity_section` says "Every classification
    below is an upstream panel's or submitter's" on each gene page, and
    `landing.py`'s hero note says every gene is admitted on an upstream
    assessment, "never this atlas's own". The sentence below says the same thing in the same voice,
    and names `atlas curation` as the column that answers the question the
    confidence column does not -- whether the atlas has curated the gene at all.

    **Whose call it is was only half of it.** All 23 of those rows read
    `definitive` with no disease anywhere on the page, and the scope rule
    (`_SCOPE_SHORT` here since D50, in full on `methods.html#scope`) plus the
    `graded for` column are the other half: a row reading `KMT2D | definitive` on
    a site titled "CHD Atlas" claims KMT2D is definitive for congenital heart
    disease, and ClinGen's actual assertion is `KMT2D -- Definitive for Kabuki
    syndrome 1`, made by the SCID-CID GCEP. `validity` is taken for that column
    alone; `facts` carries no disease, and deriving one from a second read of the
    mirrors would let this page and the gene page it links to name different
    diseases for the same gene.

    Since the gate widened the column answers a second question too. An empty
    `confidence` -- rendered as an em dash, like every empty cell here -- now
    means something specific: no ClinGen panel graded the gene, which is why
    `_graded_diseases` returns nothing for it either, and that pair of dashes is
    the only one on this table that is not a gap in the data. The paragraph
    under the heading says so (it said "a blank confidence" until 2026-10-26),
    and `render.grade_key` glosses it on `methods.html#scope`, which both page
    kinds link from their short form (it was a folded key on both pages until
    2026-10-26, D50).

    **A contested gene is not shown as settled here either.** Until 2026-10-26
    the confidence cell was the bare grade, so BRAF, KRAS, MAP2K1, PTPN11, RAF1
    and SOS1 read as ordinary `definitive` rows while `genes/index.json` behind
    this page flagged each `has_conflicting_evidence` -- the reading
    `docs/data-api.md` tells a client never to commit. The cell now carries the
    gene card's own `conflicting evidence` chip, wherever the flag is set;
    `data-confidence` stays the bare grade, because the facet matches it.
    `test_every_contested_gene_is_chipped_on_the_gene_index_and_no_other_is`
    holds it against a real build.

    Two couplings to `render.FILTER_SCRIPT` that nothing in the type system
    holds, both pinned by
    `test_every_facet_names_a_data_attribute_the_filter_script_reads`:

    * each `<select>`'s `name` is resolved into an attribute by concatenation,
      `row.getAttribute('data-' + select.name)`, so the four names below must
      spell the four `attributes` keys exactly. A mismatch reads `null` and
      hides every row the moment that facet is used, while the page looks
      correct until someone touches it.
    * the script matches a facet with `.split(' ').indexOf(want)`, so a value
      carrying a space could never be selected. None can: measured 2026-08-04,
      no member of `LesionGroup`, `Classification`, `ValidityState` or
      `AtlasCuration` contains one. That is what lets `data-lesion` hold a whole
      space-joined list and still match a single group.

    Every `<option>` value comes from a vocabulary member -- an enum in
    `vocab.py`, not a curated or mirrored string -- so nothing reaches the
    `selects` markup below unescaped. Everything else on the page goes through
    `render.py`, which escapes cells and attributes itself.
    """
    filters, table = render_combined_table(
        facts,
        symbols=symbols,
        validity=validity,
        burden_counts=burden_counts,
        concordance=concordance,
        profiles=profiles,
        identities=identities,
        prefix="../",
    )
    body = (
        "<h1>Genes</h1>"
        + '<div style="display: none;"><p class="scope-rule">'
        + _SCOPE_SHORT.removeprefix('<p class="scope-rule">').removesuffix("</p>")
        + "</p><p>Every confidence below is a ClinGen expert panel's, republished "
        "with its provenance intact &mdash; the atlas authors no validity classification "
        "of its own. The <strong>atlas curation</strong> column is what says whether the "
        "atlas has curated a gene itself, and the <strong>graded for</strong> column "
        "names the disease each panel actually graded. A confidence shown as &mdash; means "
        "no ClinGen expert panel has graded the gene for any of the diseases in scope here; "
        "its page names the submitters whose agreement "
        "admitted it.</p></div>" + f'<div style="margin-top: 1rem;"><span id="shown" '
        f'style="display: none;">{len(facts)}</span>'
        + filters
        + f'<div class="table-scroll">{table}</div>'
        + f'<div style="display: none;">{_STRIP_LEGEND}'
        + f"{_strip_notes(facts, concordance, publications or {})}</div>"
        + """<div class="concordance-key">
      <span><strong>Replication across independent cohorts:</strong></span>
      <span><span class="dot full"></span> Multi-cohort replicated</span>
      <span><span class="dot half"></span> Nominal</span>
      <span><span class="dot none"></span> No enrichment</span>
      <span><span class="dot untested">&ndash;</span> Not tested</span>
    </div>"""
        + "</div>"
    )
    emitter.write_text(
        GENE_INDEX_PAGE,
        document(title="Genes — CHD Atlas", root="../", body=body, script=FILTER_SCRIPT),
    )
