# src/chd_atlas/build/methods.py
"""`methods.html`, "How to read this atlas": site-wide explanation, written once.

D50. Every gene page repeated the same explanation -- the rule that admits a
gene, the grade key and, in later sections, the reading notes for the validity,
rare-variant and expression evidence -- so a reader met the same paragraphs on
every page and the part particular to the gene was buried under them. This page
holds that explanation once; each gene page keeps a short form and a link to the
section here that states it in full.

**Every sentence below the introduction either was moved off a gene page, the
gene index or the front page, or was written here to state in full what a gene page now says in
one line or prints without defining -- and each was read in this position.** A
sentence's truth conditions travel with its position (CLAUDE.md section 4.27):
"the table below" was true on a gene page, false on the gene index and false
here, and "this page" here names a page with no gene on it.

* `#scope` was moved: the scope rule and the grade key. Those that were false
  here were reworded -- the scope rule's last clause, and the grade key's
  dispute and completeness sentences. The rest are unchanged because they are
  true anywhere: the key's opening sentences, and the rung glosses
  (`render.GRADE_GLOSS`), which stay verbatim on purpose because a gene page's
  card glosses its headline grade from the same constant, and one rung must not
  be glossed two ways.
* `#validity` was written new (2026-10-26). Its gene-page sentences did not
  leave: the ClinGen-copy notice stays on the page (P5), and the table's intro
  keeps a short form of its scope. The section states in full what those
  abbreviate, for a page with no table on it -- so "a gene page lists", never
  "listed here" or "the panel record above". Its third paragraph states in
  full what the card's two warning chips abbreviate to two words each.
* `#curated-evidence` was written new (2026-10-26, owner decision E3): what
  a curated gene's two extra sections hold -- the atlas's own lesion
  assertions and the publications their evidence cites -- which no page
  explained while every card's "How to read this page" linked here. Each
  clause was read off `models/assertion.py`, `pages._evidence_section`,
  `derive.gene_facts` and `validity.published_genes`; it names no count, so
  nothing in it expires when a second gene is curated.
* `#rare-variants` opens with the front page's two definitions, moved here with
  the doorway page (2026-10-27). After them it is both: three paragraphs moved
  off every gene page (the folded reading notes and the composite note's
  second sentence), and two are the full forms of notices a gene page keeps in
  one line. Two more state in
  general what each forest plot's folded key says of its own marks (D57), and
  what a block with no effect measure shows instead of a plot; one more, what
  the `all cases` case group is, which a forest key says in one clause
  (2026-10-26).
* `#cohorts` is curated data under one written paragraph: every record of
  `curation/cohorts.yaml`, the same ones `cohorts.json` publishes. A gene page
  keeps the notes of its own cohorts, folded, and links here (P4).
* `#expression` is both (2026-10-26): two paragraphs moved off the gene pages
  (the folded reading notes: the placement note from every page with data,
  the comparability notice from every page that laid out more than one
  placement), one is the full form of the bulk-tissue notice a gene page
  keeps in one line, and one states in general the chart axes each chart's
  visible caption line states of its own figure (D57). Two more define what
  the pages print and nothing defined (review, 2026-10-26): tau, and the
  developmental phases each stage names and each trajectory shades.

`tests/unit/test_build_methods.py` pins the scope rule and the key's opening
paragraph against literals, sliced to the section; the glosses by a phrase
unique to each, as one label-chip-gloss triple per rung; `#validity`,
`#rare-variants`, `#expression` and, on a fixture registry, `#cohorts` whole, as
literals; the cohort intro's figures against the mirror and the committed
registry; the bulk-tissue claim, and the dataset the floor sentence names,
against the curated datasets; and the phrases the move retired as absent from
the whole page, built with no cohorts and with the committed registry.
`tests/test_built_site_is_consumable.py` checks `#cohorts` against the
`cohorts.json` of a real build.

Written through `Emitter.write_text`, for the reason `pages.py` gives: `Json`'s
union includes `str`, so `write_json` would type-check and publish the page as a
quoted JSON string. Nothing here reads a clock or iterates a set, so two builds
of one commit write the same bytes.
"""

from __future__ import annotations

import html
from collections.abc import Mapping
from typing import Final, assert_never

from chd_atlas.build.concordance import ALPHA
from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import COHORTS, DATASETS, METHODS_PAGE
from chd_atlas.build.render import (
    document,
    grade_key,
    review_status_label,
    review_status_note,
)
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.review import ReviewStatus

# What the published set *is*, stated in full here and in one line on both page
# kinds that show it (`pages._SCOPE_SHORT`, which links this section).
#
# No page on this site said it. `docs/data-api.md` states the rule and a reader
# of the HTML never sees it, so a browse page of 23 rows all reading
# `definitive`, under a title reading "CHD Atlas" and a tagline about congenital
# heart disease, said "these 23 genes are definitive for CHD". Measured
# 2026-08-04 against the committed mirrors, that was false for most of them: 10
# of the 23 were qualified by a panel that is not the Congenital Heart Disease
# GCEP -- SCID-CID for KMT2D and KDM6A, Hearing Loss for CHD7, Kidney Cystic and
# Ciliopathy Disorders for NOTCH2 and EVC, Intellectual Disability and Autism
# for ZEB2, Skeletal Disorders for EVC2, Syndromic Disorders for TBX5 and
# TFAP2B, Prenatal for FLT4 -- and 10 were definitive for a disease whose label
# names no cardiac feature at all (Kabuki syndrome 1 and 2, Mowat-Wilson,
# CHARGE, Alagille, Ellis-van Creveld twice, Holt-Oram, TARP, and NR2F2-related
# multiple congenital anomalies).
#
# The 2026-08-06 widening makes the sentence carry more, not less. Measured on
# the widened corpus: 92 genes, of which 76 are admitted by a ClinGen record at
# or above `Limited` and 16 by two or more GenCC submitters with no ClinGen
# record at all. A reader meeting `limited` needs to know both that a chartered
# panel wrote it and that this atlas did not.
#
# Moved from `pages._SCOPE_RULE` (2026-10-26, D50), where it was one constant
# for both page kinds. Its last clause said the disease "is shown beside the
# gene", true on a gene page and on the gene index and unanchored here, so it
# now names both places -- and says "for a gene a ClinGen panel graded", since
# the 21 genes admitted on submitter agreement have no panel-named disease to
# show (their card names the submitters instead).
#
# Its first part was rewritten on 2026-10-26, after review measured it false
# read as a rule (each measured that day on the committed mirrors and the
# pinned HPO release):
#
# * It named a *property of a disease* where the code reads a *list*. A disease
#   is in scope only once it is in `curation/chd_scope.yaml`, on either path:
#   SCP007 reports 239 unnominated terms that would qualify on Orphanet's
#   frequencies (FLNA's periventricular nodular heterotopia among them), and
#   Orphanet's own GenCC submissions name Ebstein anomaly (MONDO:0009144),
#   which is not a scope term. So the sentence says what a term needs to be
#   listed, and that a qualifying disease is not in scope until it is.
# * "independently assert it" was weaker than the gate, which counts a GenCC
#   submitter only at `Limited` or better (`validity.agreeing_records`).
# * "Structural heart malformation" is now defined, from
#   `orphanet_scope.structural_phenotypes` -- and two clauses are narrower
#   than the text this was drafted from. The excluded branches are named as
#   *excluded*, not as all the branches that are not malformations: the kept
#   ones still hold aortic valve calcification, mitral valve prolapse and left
#   ventricular hypertrophy, and the sentence says so. And a subtype qualifies
#   on its syndrome's frequency (`Warrant.via == "subtype"`), which the
#   parenthesis states. `test_the_structural_set_the_build_applies_is_the_one_
#   methods_html_states` (`test_orphanet_scope.py`) ties every term named here
#   to the set the build computes, so changing the rule where this names it
#   fails there until this is reworded.
#
# Widened by owner decision 2026-10-26 (R1): pulmonic stenosis, aortic valve
# stenosis, supravalvular aortic stenosis, truncus arteriosus, an interrupted
# aortic arch and anomalous pulmonary venous return joined coarctation and the
# ductus (`orphanet_scope.ADDED_BRANCHES`), so the clause that named three of
# them "not counted" now lists all eight. Two clauses are narrower than the
# owner's preview, each measured that day on the pinned release: the set still
# leaves lesions out -- peripheral pulmonary artery stenosis (HP:0004969) and an
# aortopulmonary window (HP:0011604) are filed elsewhere and not added -- so the
# sentence names two rather than implying every lesion now counts; and aortic
# valve stenosis is not always a malformation, which the sentence says beside
# the calcified valve.
SCOPE_RULE: Final = (
    '<p class="scope-rule">A gene is published here when an external authority already '
    "treats it as a congenital heart disease gene: either a <strong>ClinGen expert "
    "panel</strong> classifies it <strong>Limited or better</strong> for a disease in scope "
    "here, or <strong>two or more Gene Curation Coalition submitters</strong> classify it "
    "Limited or better for one and no ClinGen panel disputes it. A disease is in scope when "
    "it is one of the terms this atlas lists, and a term is listed only if an external "
    "authority treats it as congenital heart disease, or if a curator nominates a syndrome "
    "and <strong>Orphanet records a structural heart malformation in at least 30% of its "
    "cases</strong> (for a subtype, of the syndrome it belongs to); a disease that would "
    "qualify is not in scope until it is listed. Here a structural heart malformation is a "
    "finding under the Human Phenotype Ontology's <em>Abnormal heart morphology</em> "
    "(HP:0001627), less branches of it excluded as not malformations (cardiomyopathy, "
    "tumours, the pericardium and the like), plus eight lesions the ontology files "
    "elsewhere, with any subtypes: coarctation of the aorta, a patent ductus "
    "arteriosus, pulmonic stenosis, aortic valve stenosis, supravalvular aortic stenosis, "
    "truncus arteriosus, an interrupted aortic arch and anomalous pulmonary venous return. "
    "The findings counted still include some that are not malformations, such as a "
    "calcified or prolapsing valve, and aortic valve stenosis can be acquired; other "
    "lesions the ontology files elsewhere, among them peripheral pulmonary artery "
    "stenosis and an aortopulmonary window, are not counted. "
    "<strong>No disease is in scope on this atlas's own judgement</strong>, and neither "
    "is any gene: every scope term names the authority it came from "
    "(<code>scope.json</code>, from <code>curation/chd_scope.yaml</code>), and every gene "
    "page names the authority that admitted the gene. A classification is not the same "
    "as definitive for congenital heart disease: for many of these genes the disease "
    "graded is a syndrome of which a cardiac lesion is one feature, and for a gene a "
    "ClinGen panel graded, the disease the panel named is shown beside the gene, on "
    "its page and in the gene index.</p>"
)

# How a gene page's validity table is listed (2026-10-26, D50). The page's own
# intro states the short form -- only in-scope records are listed -- and links
# here. The copy sentence describes `pages._repeats_the_panel` and is copied
# from `docs/data-api.md`'s account of it, read against the function that day:
# a copy is hidden only on a match on disease, mapped grade, inheritance and
# calendar date, and a copy with no date matches nothing.
#
# Two clauses are narrower than the text this was drafted from, each measured
# 2026-10-26 on the committed mirrors and a real build:
#
# * "Each gene's JSON keeps every record" read, beside the next paragraph's
#   "that record is not shown", as though the JSON held the out-of-scope
#   records. It holds none -- 591 bundle records, 0 for a disease outside
#   scope -- and all 121 copies, so it says "every copy".
# * "a gene page is not every classification the gene has" was false for 43 of
#   the 119 published genes, which have no ClinGen or GenCC record for any
#   disease outside scope. "may not be", as the gene page says.
#
# Reworded after review the same day. It opened "Every classification in this
# atlas is an upstream panel's or submitter's, republished with its provenance
# intact", which the grade key's first sentence already says two paragraphs up;
# it is said once, there. "the ClinGen expert panel's own records" became
# "ClinGen expert panels'", and "A panel may have classified" became "A panel or
# submitter": the scope filter drops GenCC rows too, and most out-of-scope
# records are GenCC's (76 published genes have one, 27 an out-of-scope ClinGen
# record).
#
# The third paragraph glosses the card's two flags (2026-10-26), which no page
# explained: `methods.html` mentioned "sources disagree" zero times, and on all
# six pages carrying it ClinGen and GenCC agree on every gene-disease pair they
# both classify. Each clause is read off the code that sets the flag:
# `vocab.has_conflicting_evidence` over every in-scope record with a mapped
# classification (`derive.gene_facts`), and `validity._has_source_discordance`,
# whose sources are `ValiditySource.CLINGEN` and `ValiditySource.GENCC` -- so a
# ClinGen copy submitted to GenCC is GenCC's there. `NO_KNOWN_ASSOCIATION` is
# subtracted from the supporting side and is not in `CONTESTED`; `Supportive`
# maps to no classification and is read by neither. The flag is per gene by the
# owner's decision (2026-10-25), so the paragraph says what that permits rather
# than narrowing it. "Narrower" is exact: a gene with `has_source_discordance`
# has a contesting and a supporting record, which is `has_conflicting_evidence`.
_VALIDITY: Final = (
    '<h2 id="validity">Mirrored validity records</h2>'
    "<p>A gene page lists its classifications in two groups: ClinGen expert panels' own "
    "records first, then other curators' submissions to the Gene Curation Coalition "
    "(GenCC). ClinGen also submits its own classifications to GenCC; a gene page lists "
    "such a copy only where its mapped grade, date, inheritance or disease differs from "
    "the panel's own record, or where it carries no date or has no panel record to "
    "match. Each gene's JSON keeps every copy, listed or not.</p>"
    "<p>Only records for diseases in this atlas's remit are listed. A panel or submitter "
    "may have classified the same gene for a disease outside it, and that record is not "
    "shown, so a gene page may not be every classification the gene has.</p>"
    "<p>A gene page's card may carry two flags beside its grade. <em>conflicting "
    "evidence</em> means that among the gene's records for diseases in scope, at least one "
    "is <em>Disputed</em> or <em>Refuted</em> and at least one classifies the gene Limited "
    "or better, whoever made each and whatever disease each is for; the gene index marks "
    "the same genes beside their confidence. <em>sources disagree</em> is narrower: a "
    "record that disputes or refutes the gene and one that classifies it Limited or better "
    "come from different sources, one a ClinGen expert panel's own record and the other a "
    "submission to GenCC "
    "(ClinGen's own copies count as GenCC's here). It too is a flag on the gene, not on a "
    "disease: the two records may be for different diseases, and the two sources may agree "
    "on every disease they both classify. Neither flag counts a <em>No Known Disease "
    "Relationship</em> record on either side, nor GenCC's <em>Supportive</em>, which "
    "grades no evidence.</p>"
)

# What a curated gene's page adds (owner decision 2026-10-26, E3). A gene page
# renders "Curated evidence" and "Publications" only for a gene with a lesion
# assertion (`pages._evidence_section`, gated on `AtlasCuration.CURATED`), and
# this page explained every other section but those two. Each clause read off
# the code that day:
#
# * An assertion is `LesionAssertion`: a curator's record that the gene belongs
#   to one or more lesion groups (`lesion_groups`, at least one), with
#   `inheritance`, `mechanism` and `syndromic` (isolated, syndromic or both) --
#   the four rows the page's definition list renders.
# * Each `Evidence` item cites exactly one publication (`publication`, a
#   required PMID) and carries `evidence_class`, `strength` (strong, moderate,
#   supporting or contradictory) and the curator's own `summary` -- the
#   table's four columns. The owner's draft said "quoted evidence": `summary`
#   is the curator's summary, not a quotation, and the record locates the
#   passage instead (`locator`), which the page does not render.
# * It changes nothing mirrored, and that is D37: `validity.published_genes`
#   reads the validity records alone, and `derive.gene_facts` headlines the
#   strongest ClinGen record whatever the assertions say; an assertion feeds
#   only `lesion_groups` and the curation fields.
# * "Publications" lists the distinct publications the evidence cites, sorted,
#   with their id, title and year (`_CITATION_HEADERS`), and renders exactly
#   where "Curated evidence" does. Named as a curated gene's page's own: it
#   read "The Publications section that follows" until review (2026-10-27),
#   and on this page what follows is the rare-variant section.
_CURATED_EVIDENCE: Final = (
    '<h2 id="curated-evidence">Curated evidence</h2>'
    "<p>Where this atlas has curated a gene itself, its page carries a <strong>Curated "
    "evidence</strong> section. Each entry is a <strong>lesion assertion</strong>: a "
    "curator's record that the gene belongs to one or more congenital heart disease lesion "
    "groups, with its inheritance, its mechanism and whether it is isolated, syndromic or "
    "both. Each item of evidence behind an assertion cites one publication, and gives the "
    "class of the evidence, its strength as the curator grades it (strong, moderate, "
    "supporting or contradictory) and the curator's summary of what the publication "
    "shows.</p>"
    "<p>An assertion is this atlas's own work, and it <strong>changes nothing "
    "mirrored</strong>: a gene is published, and its grade is headlined, on the upstream "
    "validity records alone, so no assertion admits a gene or moves its grade. A curated "
    "gene's page also has a <strong>Publications</strong> section, which lists every "
    "publication its evidence cites, with its title and year.</p>"
)

# How to read a gene page's rare-variant section (2026-10-26, D50). The page
# keeps a one-line short form of the two paragraphs after the front page's
# definitions (`pages._MIRRORED_NOTICE` and `pages._POOLING_NOTICE`) and links
# here; the three after those were the page's folded reading notes and its
# composite note's second sentence, and left it.
#
# Narrower than the text this was drafted from, each measured 2026-10-26 on a
# real build of 4233eca:
#
# * "Every count and statistic in a gene page's rare-variant section" -> "in a
#   gene page's study tables". The evidence matrix in that section shows tallies
#   the atlas derives -- "n = N" in a cell on 86 pages, "N rows" on 48 -- and
#   each PMID:42230622 block carries the atlas's own Bonferroni arithmetic
#   (`0.05 / 138,609`); neither is an upstream study's number.
# * "republished exactly as published" -> "as published". The tables round to
#   three significant figures: 376 of the 1,210 published raw p-values render
#   rounded (PMID:40127276's 0.00583259259259259 reads 0.00583). Each gene's
#   JSON carries them in full.
# * The units sentence said `cases` and `controls` "each show observed / total".
#   On a row tested against a mutation model the `controls` cell is the model's
#   expected count (`pages._controls`, which keys on the comparator): 108 cells
#   on 36 pages read "0.0358 expected" and the like. It says so.
#
# Reworded after review the same day: "children" -> "people", because
# `shared_cohorts` counts shared control collections too and a control cohort
# may be adults (UK Biobank's note: "~500,000 UK adults recruited at ages
# 40-69"); the overlap clause says what the page can know -- collections shared
# *as this atlas has recorded them* -- since `shared_cohorts` sees only shared
# registry ids and declared overlaps, and two studies may share individuals no
# entry names (review raised DDD and DECIPHER, which the registry has declared
# overlapping since owner decision 2026-10-27, I3); the synonymous paragraph
# lost "the same uncorrected scale" (true of the one study that publishes a
# synonymous row, stated in general) and "a scan of that size" (no antecedent
# here).
# `concordance.ALPHA` as the rule's sentence prints it: "0.05".
_ALPHA: Final = f"{ALPHA:g}"

_RARE_VARIANTS: Final = (
    '<h2 id="rare-variants">Rare-variant association</h2>'
    # The front page's two definitions, moved here with the doorway page
    # (2026-10-27): its hero note links this section rather than stating them.
    # D49: a statistic names all three analyses, and "copy-number variants"
    # rather than deletions, so a future non-deletion CNV analysis does not
    # make it false.
    "<p>A <strong>rare-variant statistic</strong> is one study's rare-variant test of one "
    "gene under one design: a case-control burden test of rare SNVs and indels or of "
    "copy-number variants, or a de novo enrichment test, which are different analyses and "
    "are named as such on each gene's page. An <strong>independent dataset</strong> is a "
    "cohort family rather than a paper: two studies drawing on the same sample collection, "
    "or on collections this atlas records as overlapping, can describe the same people and "
    "count once. How many datasets agree about a gene is shown on every row of the browse "
    "page, and it is a count of evidence &mdash; never a verdict on the gene.</p>"
    # The one rule the atlas applies (owner decision 2026-10-26, W2). It said
    # "adds no assessment of them" while the card, the matrix and the forest
    # fill each decide what survives by `concordance.survives_own_correction`.
    # Read off the code that day: `_supports` (not synonymous, a raw p below
    # `ALPHA` -- strictly -- and `is_upward`, which reads the effect, an
    # unbounded bound, or the published counts), then a study-published
    # `pvalue_adjusted` below `ALPHA`, strictly; a row with no corrected p does
    # not survive; `family_state` reads the matrix tier off the same rows;
    # `pages._survived_correction` fills a forest marker by it for every row but
    # a synonymous one, which is drawn muted (the forest paragraph below says
    # so). The threshold is rendered from `ALPHA`, so this page cannot state a
    # number the rule does not use. It read "at the 0.05 the studies themselves
    # use" until review (2026-10-27): no study here treats a raw p below 0.05 as
    # significant -- PMID:42230622's block tells a reader to judge raw p against
    # its whole scan -- so the number is stated as the atlas's, unattributed.
    #
    # The sentence after "survives none" (owner decision 2026-10-26, R3): a
    # survival is only as wide as the correction, and PMID:40127276's is a
    # Benjamini-Hochberg over its 248 panel genes (`Publication.
    # correction_genes`, reproduced from Dataset S6 to 7e-15), which the card
    # now names beside the q (`pages._correction_scope`), marking every
    # surviving row below the study's own genome-wide line
    # (`pages._genome_wide`). Stated in general; true today of the one study
    # whose correction ran over a panel, the only one curated with either field.
    "<p>Every count and statistic in a gene page's study tables is an upstream study's, "
    "republished as published. The atlas computes none of them. It applies one rule of its "
    f"own, at {_ALPHA}: a row that is not synonymous "
    "<strong>supports</strong> the gene when its result points to enrichment in cases and "
    f"its p is below {_ALPHA}, and <strong>survives</strong> its study's correction when the "
    f"corrected p its study published is below {_ALPHA} as well; a row whose study published "
    "no corrected p survives none. A study's correction covers only the genes it was run "
    "over: where those were a gene panel rather than the genome, the card says so beside "
    "the corrected p, and marks any surviving row whose raw p is also below the threshold "
    "the study itself states for genome-wide significance. A gene page's card lists the "
    "rows that survive, the "
    "evidence matrix marks a dataset whose rows survive apart from one whose rows only "
    "support, and a forest plot draws a surviving row filled. Where a study is by an "
    "author of this atlas, its block says so.</p>"
    # was pages._POOLING_NOTICE in full. The overlap clauses since owner
    # decision 2026-10-27 (I3): DECIPHER and DDD are declared overlapping
    # (`Cohort.overlaps_with`), and the page's notice names such a pair.
    # "could count", not "would" (review, 2026-10-27): the people are counted
    # twice only if the overlap the clause before calls a "may" is real.
    "<p>The atlas computes <strong>no pooled statistic across studies</strong>. Where two "
    "studies draw on the same sample collection they describe partly the same people, and "
    "where they draw on collections this atlas records as overlapping they may, so "
    "combining their results could count those people twice; a gene page names any "
    "sample collection its studies share, and any overlap between theirs, as this atlas "
    "has recorded them.</p>"
    # was pages._BURDEN_PREAMBLE; "In each study's table," anchors what it describes
    "<p>In each study's table, <code>cases</code> and <code>controls</code> each show "
    "<em>observed / total</em> (for a row tested against a mutation model, "
    "<code>controls</code> is instead the count that model expected), and <strong>every "
    "cell names what it counted</strong> &mdash; carriers, alleles, or de novo mutations. "
    "The three are not interchangeable: a study counting alleles counts a person with two "
    "qualifying variants twice, and its denominator is roughly twice the number of people "
    "sequenced. A consequence class or case group with no row is one the study "
    "<strong>did not report</strong>, which is not a null result and does not always "
    "mean no such variant was seen &mdash; studies differ in which classes their "
    "gene-level tables cover.</p>"
    # was pages._SYNONYMOUS_NOTICE; "above it" -> "beside it"
    "<p>The <strong>synonymous</strong> row is a study's own negative control &mdash; "
    "synonymous variants should show no enrichment. Read it on the same scale as the rows "
    "beside it: a <em>strongly</em> enriched synonymous row is a warning about that "
    "gene's comparison, while a nominally significant one is what a scan over thousands "
    "of genes produces by chance. <strong>Only a study that publishes one has "
    "it</strong>; a table with no synonymous row offers no such check.</p>"
    # was the second sentence of pages._composite_note
    "<p>A <strong>damaging (LOF + missense)</strong> row is the <strong>union</strong> of "
    "the study's loss-of-function and damaging-missense rows, not a third independent "
    "result &mdash; the same variants are counted in both. It is shown because a study "
    "that defines its findings by the composite publishes a p-value for it that is not "
    "derivable from the other two.</p>"
    # The case-group union, which nothing here stated until 2026-10-26: the
    # paragraph above explained `damaging` while a table's `all cases` row sat
    # beside the rows it contains. Measured that day on the committed mirror:
    # all three studies publish an `all` stratum and only PMID:42230622 also
    # splits its cases, into syndromic (1,471) and non-syndromic (2,405); on all
    # 470 cells where it publishes all three, the `all` row's cases and carriers
    # are the sum of the other two's and its controls are theirs. Stated in
    # general, so it names no study; `test_an_all_cases_row_is_the_union_of_
    # the_case_groups_its_study_splits` holds it to the mirror.
    #
    # The last clause said `all cases` was a study's "whole case set" until
    # review on 2026-10-26: it is not a denominator. PMID:34324492's rows count
    # 6,746 cases of the 7,958 its own note gives the study, and
    # PMID:40127276's count 3,887 trios or an allele total against 11,555
    # probands. What holds of every study is the narrower claim: not a subgroup.
    "<p>Where a study splits its cases into <strong>syndromic</strong> and "
    "<strong>non-syndromic</strong> groups, its <strong>all cases</strong> group is "
    "likewise the union of the two, and an <strong>all cases</strong> row is not a result "
    "independent of theirs. Where a study does not split its cases, <strong>all "
    "cases</strong> is its only case group: every case the row counts, not a subgroup.</p>"
    # The forest key in general; each panel's own key, conditional clause by
    # clause, is folded under "How to read this chart" (`pages._forest_caption`,
    # D57). Each clause measured against a real build of 9f24190 before it was
    # stated in general: all 29 right arrows are odds ratios with 0 control
    # carriers and a published lower bound; all 257 left arrows have 0 case
    # carriers; every plotted study corrects every row or none. "contains other
    # rows in the same panel" was false for 10 union-tagged rows, on 10 panels,
    # that contain exactly one other drawn row, so it says "at least one", as
    # each panel's key does.
    #
    # Two clauses reworded after review (2026-10-26): "these have none" had no
    # antecedent nearer than "a forest plot", so it names this atlas's; and the
    # muted row was singular, while every panel that draws one draws three
    # (measured on a real build: 113 panels, one synonymous row per case group
    # of PMID:42230622), so it is "one of a study's synonymous rows", true for
    # one or many.
    #
    # One sentence added 2026-10-26, when a marker's fill began reading
    # `concordance.survives_own_correction` for every row but a synonymous one
    # (`pages._survived_correction`): a depleted row that cleared its
    # correction is drawn hollow, and "did not" survive would be false of it.
    # 0 published rows are in that state (measured that day). "Non-synonymous"
    # since review the same day: a synonymous row that survived its correction is
    # drawn filled and muted (`pages._mark_class`), never hollow -- latent, as no
    # mirrored synonymous row carries a corrected p.
    "<p>A forest plot shows each row's effect and interval on a logarithmic axis. A "
    "<strong>filled</strong> marker is a row that survived its own study's correction for "
    "multiple testing; a <strong>hollow</strong> marker is a row that did not, or whose "
    "study published no correction at all &mdash; the atlas computes none of its own. A "
    "non-synonymous row that survived its correction without pointing to enrichment in "
    "cases is drawn hollow too: it is not evidence for the gene. A forest plot usually "
    "ends in a summary diamond; this atlas's have none, because no pooled estimate is "
    "computed anywhere on this site. A <strong>hatched</strong> bar "
    "means no interval was published for that row, not a wide one. An arrow to the "
    "<strong>right</strong> is an effect unbounded above &mdash; no control carried the "
    "variant &mdash; so the study's lower bound is the whole finding; an arrow to the "
    "<strong>left</strong> is an effect of zero &mdash; no case carried the variant "
    "&mdash; which a logarithmic axis has no position for. A row tagged "
    "<strong>union</strong> contains at least one other row in the same panel and is not "
    "an independent finding. A <strong>muted</strong> row is one of a study's synonymous "
    "rows, drawn apart from the results it calibrates.</p>"
    # The general form of `pages._NO_EFFECT_MEASURE`, which a block shows in
    # place of a plot (77 pages, every one in PMID:34324492's block).
    "<p>A study that publishes <strong>no effect measure</strong> for a gene &mdash; no "
    "odds ratio, enrichment or rate ratio &mdash; has nothing to place on an effect axis, "
    "so its block shows its counts and p-values in a table instead of a plot.</p>"
)

# Every collection the rare-variant rows cite, with its curated description
# (P4, 2026-10-26). A gene page folds the notes of only the collections its own
# rows cite (`pages._cohort_notes`, whose fold ends on a link here) and keeps
# them: CLAUDE.md section 9 forbids narrowing where a cohort caveat reaches, so
# this is a second place, not a move. It is the whole registry -- the records
# `cohorts.json` publishes -- and every registered collection is cited by a
# mirror row, so "every sample collection a rare-variant row names" describes
# the list both ways; `test_the_cohort_list_says_what_its_rows_name_and_what_
# they_do_not` holds that against the mirror.
#
# It read "every sample collection a rare-variant study in this atlas draws on,
# with what qualifies every number drawn from it" until review (2026-10-26)
# measured both halves false. PMID:34324492's rows named 4 collections
# (`decipher`, `isca`, `signature`, `taa_cases`), and the study's S2 Table lists
# 189 case sources holding 7,958 cases, of which those 4 hold 4,242 -- counted
# from the supplement's `s007.xlsx` (sheet "S2 Table": 189 source rows whose
# individuals sum to 7,958, matching the sheet's own "N" row; Signature 1,615,
# Decipher_database 1,221, PMID 21092924 777, ISCA 629). A reader told the list
# is complete concludes that two studies naming no collection in common share
# no people, which is the inference the next sentence must not license; the
# gene page's pooling sentence was already hedged ("as this atlas has recorded
# them"). And some descriptions qualify nothing (`isca`'s is a case count), so
# the list carries "the caveats its curator recorded", whatever they are. The
# sources the rows do not name are not named here either: an unregistered
# source has no record to resolve against.
#
# 11 and 6,211 since owner decision R12 (2026-10-26) registered the seven
# CHD case sources S2 identifies by PMID with 50 or more cases, 1,969 between
# them (S2's own counts, each recorded in its `curation/cohorts.yaml`
# description). The 178 S2 sources still unnamed hold 1,747 cases, so "two
# that name neither may still share people this registry does not record"
# still covers them. The 11 is pinned against the mirror, 6,211 against the 11
# descriptions' own counts and 7,958 against `taa_cases`'s; 189 is the
# supplement's alone and nothing here can check it.
# The declared-overlap clause since review (2026-10-27): "two that name none in
# common may still share people this registry does not record" stopped being
# the whole truth when the registry recorded one such pair (`ddd` and
# `decipher`, owner decision I3).
_COHORTS_INTRO: Final = (
    '<h2 id="cohorts">Cohorts</h2>'
    "<p>Every sample collection a rare-variant row in this atlas names, with the caveats "
    "its curator recorded. A row may not name every collection its study drew on: the "
    "rows for PMID:34324492 name 11 of the 189 case sources its S2 Table lists, and those "
    "11 hold 6,211 of its 7,958 cases. Two studies that name the same collection describe "
    "partly the same people; two that name collections this registry records as "
    "overlapping may; and two that name neither may still share people this registry "
    "does not record. The same records are published as "
    f'<a href="{COHORTS}"><code>{COHORTS}</code></a>.</p>'
)


def _review_gloss(status: ReviewStatus) -> str:
    """What a record carrying `status` is, said after its label (`_review_key`)."""
    match status:
        case ReviewStatus.PENDING_REVIEW:
            return (
                "carries a question about it that this atlas leaves open for an expert: its "
                "note states what is unresolved, and the atlas takes no side."
            )
        case ReviewStatus.NOT_ACCESSIBLE:
            return (
                "rests on a source that could not be read: its note states what that leaves "
                "unchecked, and the atlas does not guess at it."
            )
        case _:
            assert_never(status)


def _review_key(cohorts: Mapping[str, Cohort]) -> str:
    """What each review label in the list marks, for the labels it prints; `""` for none.

    Review (2026-10-27, M8): the list printed "Pending review:" after two
    descriptions, and every gene page's cohort fold links here, but no page
    said what the label means or that the atlas takes no side. One
    sentence per status the registry carries, in the vocabulary's order,
    and nothing when no record carries one: a key for a label the page does
    not print defines nothing a reader meets. On the committed registry, 2 of
    20 records carry `pending_review` and none `not_accessible` (measured
    2026-10-27), so the page prints the first sentence alone.
    """
    present = {cohort.review_status for cohort in cohorts.values()}
    sentences = [
        f"A record marked <strong>{review_status_label(status)}</strong> {_review_gloss(status)}"
        for status in ReviewStatus
        if status in present
    ]
    return f"<p>{' '.join(sentences)}</p>" if sentences else ""


def _cohorts(cohorts: Mapping[str, Cohort]) -> str:
    """The registry as a definition list: name and id, then the description.

    Sorted by id here, because `burden.cohort_registry` returns the YAML's order
    and an array derived from a mapping is the caller's to sort (CLAUDE.md
    section 9) -- the same order `cohorts.json` publishes. Every curated field
    is escaped: names, descriptions and urls are free text. A url is the
    collection's public page, linked inside its own description; a cohort with
    none (11 of 20, measured 2026-10-27) gets no link rather than an empty
    one. The link reads "website", not "source" (2026-10-26): it leads to the
    collection's home page, and "source" beside a description reads as the
    citation for it.

    A record carrying an open question (owner's standing rule, 2026-10-27)
    says so right after its description, before any link: its status's label
    ("Pending review:" or "Source not accessible:") and the curator's note,
    escaped (`render.review_status_note`, which a gene page's cohort notes
    share). A record with none gets nothing. Before the list, where any record
    carries one, a sentence says what each printed label marks (`_review_key`).

    `cohort-list` is layout only. The base `dl` is a two-column grid whose label
    column grows to its widest label, which suits the grade key's chips and
    would give a cohort's full name most of the line.
    """
    items = "".join(
        f"<dt>{html.escape(cohort.name)} (<code>{html.escape(cohort.id)}</code>)</dt>"
        f"<dd>{html.escape(cohort.description)}"
        + review_status_note(cohort.id, cohort.review_status, cohort.review_note)
        + (f' <a href="{html.escape(cohort.url)}">website</a>' if cohort.url else "")
        + "</dd>"
        for cohort in sorted(cohorts.values(), key=lambda cohort: cohort.id)
    )
    return f'{_COHORTS_INTRO}{_review_key(cohorts)}<dl class="cohort-list">{items}</dl>'


# How to read a gene page's developmental expression section (2026-10-26, D50,
# D57). The page keeps a one-line short form of the first paragraph
# (`pages._BULK_TISSUE_NOTICE`, one of the ten caveats the spec keeps
# visible) and links here; the next two were its folded "How to read these
# figures" notes, and left it; the next two define tau and the developmental
# phases, which the pages print and nothing defined (review, the same day, each
# recorded where it is written); the sixth states in general the axes each
# chart's visible caption line states of its own figure, and is where every
# chart's "How to read this chart" toggle links (`pages._chart_notes`); the
# last two, added 2026-10-28, say how a cell-type figure is computed and what
# it does not say.
#
# Each was re-read for a page with no gene, dataset or figure on it, and read
# against the code and a build of 0fc451a that day. Narrower or different from
# the text this was drafted from:
#
# * The bulk claim names the developmental *timeline*: the profile datasets
#   whose stage series the gene page's section renders, all `bulk_rnaseq`
#   today (E-MTAB-6814 alone). It named every "developmental expression
#   measurement" until the final review (2026-10-28); the cell-type
#   paragraphs at the end of this section made that false, being single-cell
#   values that pool 9 to 15 weeks rather than follow a timeline.
#   Unconditional by design, it reads no data, so
#   `test_every_timeline_dataset_is_the_bulk_tissue_the_expression_notes_say`
#   holds it to every `profile`-design dataset instead; it held every curated
#   dataset from review (2026-10-26) until the sentence was scoped, since a
#   contrast dataset of fetal tissue was a developmental measurement and is
#   on no timeline. Its "the detection floor" was "a dataset's" for one
#   revision and is "the" again since review (2026-10-26): "a dataset's" read
#   as the source's, and the second paragraph now says whose a floor is.
# * "Each organ's median abundance" had a gene in view; it is "a gene's median
#   abundance in an organ". "the number of genes measured" is the number the
#   grid was built from (`Placement.n_genes`, 58,735 in every cell of
#   E-MTAB-6814's grid), and it travels with every percentile "a gene page
#   shows" -- "shown below" named nothing here. Measured: all 11,331
#   percentiles on the 119 pages are stage-table cells reading "N of 58,735
#   genes".
# * "a dataset's own detection floor" read as the source's, and the one floor
#   in the corpus is this atlas's curated choice (E-MTAB-6814's `floor_source`:
#   "Curated, not mirrored: the Expression Atlas declares no threshold for this
#   experiment"). Dropping "own" did not fix that -- "its dataset's" still
#   credits the floor to the dataset, which a reader takes for the source's,
#   and no HTML page said otherwise (review, 2026-10-26). So the paragraph now
#   says it: a floor is a threshold this atlas records per dataset, not
#   necessarily its source's; E-MTAB-6814 is named as a dataset whose source
#   declares none; and `datasets.json`, which publishes `detection_floor` and
#   `floor_source` for every dataset, is linked for the value and its basis.
#   "every developmental dataset's record" holds by construction: a
#   developmental dataset here is a profile one, as in the bulk claim above,
#   and a profile dataset declares both fields or does not load
#   (`Dataset.a_profile_dataset_is_fully_declared`). No number: a floor is per
#   dataset, and the record carries it. No percentile is published below it: `profiles.placement`
#   returns `None` for a median under the floor, the bundle publishes that
#   `null`, and `percentile_annotations` skips it for the omics preview.
# * The comparability sentence was conditional on a page laying out more than
#   one placement (118 of 119 did); stated once in general it has no position
#   to be false in. "on this page" became "in this atlas" and "elsewhere on
#   this page" "for another". "as the heart matures" became "over development",
#   because every organ's grid moves, not the heart's alone: liver's 90th
#   percentile is 24 tpm at 4 wpc and 5 at adolescence, heart's 28 and 6.
# * The axes sentence was checked against `_even_x` (both charts space a
#   dataset's stages evenly, in `Stage.order`), `_small_multiples` (one
#   `LogScale` per dataset block, over every organ drawn) and `_trajectory`
#   (one chart per declared cardiac organ, its axis `_axis_bounds` over that
#   organ's placed medians). "share one abundance axis" is "within a dataset":
#   a second dataset's grid would draw its own.
_EXPRESSION: Final = (
    '<h2 id="expression">Developmental expression</h2>'
    # was pages._BULK_TISSUE_NOTICE in full, under its old name
    # `_BULK_DILUTION_NOTICE`; "below" -> "in this atlas". "whole" left it on
    # the owner's decision of 2026-10-26 (E2): E-MTAB-6814's curation record
    # sets the Expression Atlas's organ tokens (forebrain, hindbrain) against
    # the paper's (brain, cerebellum), so a sample need not be a whole organ;
    # the dilution argument rests on bulk, mixed cell types, alone.
    '<p class="notice-inline"><strong>Every measurement on a developmental timeline in this '
    "atlas is from bulk tissue, never a single cell type.</strong> An embryonic heart "
    "sample mixes cardiomyocytes with endocardium, epicardium, cardiac neural crest, "
    "cushion mesenchyme and blood, and many congenital heart disease genes act in one rare "
    "lineage &mdash; the cardiac neural crest and outflow tract, the second heart field, "
    "the conduction system. A gene expressed intensely in one such lineage can still read "
    "low, or fall below the detection floor, once diluted across the sample: "
    "a low or absent bulk figure is not evidence that a gene is unimportant to heart "
    "development.</p>"
    # was pages._EXPRESSION_READING_NOTES
    "<p>A gene's median abundance in an organ is placed against that organ's own "
    "percentile grid for the same dataset and developmental stage, and the number of "
    "genes that grid was built from travels with every percentile a gene page shows. No "
    "percentile is published for a median below its dataset's detection floor, because "
    "the rank of an unreliable measurement is not itself a measurement. A detection floor "
    "is a threshold this atlas records per dataset, not necessarily one its source "
    "declares: where a source declares none, as the Expression Atlas declares none for "
    "E-MTAB-6814, the atlas chooses one, and every developmental dataset's record in "
    f'<a href="{DATASETS}"><code>{DATASETS}</code></a> gives its floor and where that '
    "floor comes from.</p>"
    # was pages._PERCENTILE_COMPARABILITY_NOTICE
    "<p>Percentiles in this atlas are <strong>not comparable across organs or across "
    "developmental stages</strong>: each organ transcribes a different fraction of the "
    "whole gene set, and an organ's reference distribution itself changes shape over "
    "development. A high percentile shown for one organ or stage says nothing about a "
    "lower one shown for another.</p>"
    # tau, which 118 gene pages print once per stage and nothing on the site
    # defined (review, 2026-10-26; counted on a build of 802d694). Each clause
    # read off the code that day: `profiles.specificity` and `TAU_METHOD`
    # (log2(median + 1), raw medians even below the floor, the peak gated
    # `< floor`), `pages._specificity_sentence` (the argmax organ, or "the
    # peak is tied") and `pages._specificity_gap_sentence` (a reason for every
    # gap). 0 and 1 are the formula's own bounds, not a banding of the values
    # between: D39(c) keeps every adjective about a tau off the site. A
    # negative median is clamped to zero (`TAU_METHOD` says so); none of the
    # 22,015 rows of `mirrors/profiles/E-MTAB-6814.tsv` is negative, so the
    # "1" clause names zero alone, a sufficient condition and a true one.
    "<p><strong>τ</strong> (tau; Yanai et al. 2005) measures how concentrated a gene's "
    "expression is across the organs a dataset sampled at one stage, computed on "
    "log2(median + 1): it is 0 when every organ's median is equal, and 1 when one organ's "
    "median is above zero and every other's is zero. It says nothing about <em>which</em> "
    "organ, so beside it a gene page names the organ with the highest median, or says the "
    "peak is tied. Every organ's median counts, even one below the detection floor; the "
    "floor decides only whether τ is published at all, and it is not when every organ "
    "sampled is below it. Where a stage has no τ, the page says why in its place, and each "
    "gene's JSON carries the formula and the medians beside every τ it publishes.</p>"
    # The phases every stage block names and the trajectory shades (review,
    # 2026-10-26: "Developmental phase" on 119 pages, "Shaded:" on 112, on a
    # build of 802d694, and this section said nothing of either). Read against
    # `curation/cardiac_phases.yaml` (its header, and each phase's Carnegie
    # stages and `end_basis`), `CardiacPhaseFile.phases_for` (`start <= wpc <
    # end`), `profiles.assign_phase` (the reasons for no match) and
    # `pages._banded_phases` (a band from the first stage naming a phase to the
    # last, read off the same `phase_ids`). The source is named in the form
    # the chart's key renders `attributed_to`. Every stated boundary is the
    # start or end of a Carnegie stage; heart looping's end is the one
    # `derived`, from GO:0001947's definition read against septation's stated
    # start, and it is named because the vocabulary holds it, though no
    # published stage falls in that phase.
    "<p>A stage's <strong>developmental phases</strong> are the phases of heart "
    "morphogenesis under way at its age, from a vocabulary transcribed from Buijtendijk "
    "MFJ, Barnett P, van den Hoff MJB (2020), PMID:32048790. A phase runs from the start "
    "of one Carnegie stage the review names to the end of another. Each boundary is "
    "converted to weeks post conception from the day ranges the Human Developmental "
    "Stages ontology gives its Carnegie stage, not read off the review's prose, which "
    "counts weeks ordinally: its &ldquo;sixth week&rdquo; is five to six weeks elapsed. "
    "Where the review names no end, as for heart looping, the atlas derives one from the "
    "phase's Gene Ontology definition. A stage names every phase under way at its age, "
    "from the phase's start up to but not including its end, so one stage can name "
    "several; a stage no phase covers names none and says why, for example that it is "
    "post-natal or after the curated window. A chart following one cardiac organ shades "
    "each phase from the first stage "
    "that names it to the last, and its folded key names the phases and their "
    "source.</p>"
    # The chart axes in general; each chart's own visible line states them of
    # its figure (`pages._band_caption`, `pages._spark_caption`).
    "<p>Every expression chart lays its stages out left to right in curated developmental "
    "order, evenly spaced: the axis is <strong>order, not elapsed time</strong>. Within a "
    "dataset, a gene's small multiples &mdash; one panel per organ with a placed "
    "measurement &mdash; share <strong>one abundance axis</strong>, which is what makes a "
    "gene concentrated in one organ look different from a broadly expressed one; a "
    "per-organ axis would rescale every line to its own range and hide exactly that "
    "difference. A chart following one cardiac organ across a dataset's stages fits its "
    "vertical axis to the gene's own placed medians in that organ, so a curve's height "
    "compares nothing to another gene's.</p>"
    # The single-cell layer (Farah et al. 2024): how a cell-type figure is
    # computed and what it does not say.
    "<p><strong>Cell types</strong> come from single-cell RNA sequencing of 8 hearts at 9 to "
    "15 weeks (Farah et al. 2024). For each heart, this atlas sums the counts of every cell "
    "of one type, and of that type in one region, and divides by all of that group's counts: "
    "counts per million (CPM). A page shows the median across hearts and the middle half of "
    "hearts. A heart counts only with 10 or more cells of that type there. Values below 1 CPM "
    "are below this atlas's detection floor, which borrows the Human Protein Atlas's "
    '"not detected" cut of 1. '
    "A dot's size is the share of that type's cells with at least one count.</p>"
    "<p>Regions pool all four ages. τ is computed across cell types as for organs above. "
    "Single-cell methods dissociate tissue, so how many cells of a type were captured "
    "reflects the method, not the heart. A high share of cells describes that cell type, "
    "not how common the type is.</p>"
)


# The page's own account of what it holds (CLAUDE.md section 4.29); the cohort
# list joined it on 2026-10-26. It names only sections the page renders, so the
# expression clause joined it only with `#expression` (2026-10-26) -- it was
# removed once for promising a section the page did not yet hold. "the
# rare-variant studies draw on" became "the rare-variant results name" the same
# day, for the reason `_COHORTS_INTRO` records: the list is what the rows name,
# not all a study drew on. The curated-evidence clause joined it with
# `#curated-evidence` (E3), in gene-page order.
_INTRO: Final = (
    "<h1>How to read this atlas</h1>"
    "<p>What every gene page assumes, written once: which genes are here and what each "
    "classification means, how the validity records are listed, what the atlas's own "
    "curated evidence is, how to read the rare-variant results, every sample collection "
    "the rare-variant results name, and how to read the developmental expression "
    "figures.</p>"
)


def build_methods(emitter: Emitter, *, cohorts: Mapping[str, Cohort]) -> None:
    """Emit `methods.html`, the site's reading notes (D50).

    `cohorts` is the registry `runner.build_site` already hands the gene pages
    (`burden.cohort_registry`); `#cohorts` lists all of it.
    """
    sections = (
        '<h2 id="scope">Which genes are here, and what a classification means</h2>'
        + SCOPE_RULE
        + grade_key()
        + _VALIDITY
        + _CURATED_EVIDENCE
        + _RARE_VARIANTS
        + _cohorts(cohorts)
        + _EXPRESSION
    )
    emitter.write_text(
        METHODS_PAGE,
        document(title="How to read this atlas — CHD Atlas", root="", body=_INTRO + sections),
    )
