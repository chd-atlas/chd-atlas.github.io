# src/chd_atlas/build/manifest.py
"""The manifest: what was built, from which commit, and a checksum per file.

This is what a consumer verifies a download against and what a cache keys on,
so `files` maps each published URL to the sha256 of the exact bytes served at
it — the compressed bytes for a `.gz`, not the JSON inside. A checksum that
disagrees with the file is worse than no manifest at all: a consumer verifying
a download rejects a file that is correct.

Nothing here varies between two builds of the same commit. A build timestamp
was considered and rejected: it would make the output irreproducible, so no
test could assert byte-identical rebuilds and the checksums would stop being a
pure function of the inputs. A consumer wanting a publication date reads the
commit's.

The manifest cannot list itself — its own checksum would have to be computed
from bytes containing that checksum. The exclusion is not a filter but an
ordering: `write_json` encodes the payload before it records the new entry, so
`manifest.json` enters `checksums` a moment after the bytes naming every other
file are already fixed.
"""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.corpus import Corpus

# What a consumer can rely on across two builds, in the usual major.minor sense
# — spelled out here because the field was published for three releases without
# anyone saying what it meant.
#
# MINOR rises when a field is added and nothing is removed or repurposed: a
# consumer written against the previous version keeps working, and one that
# wants the new field can test for the version rather than for the key.
# MAJOR rises when a field changes shape or leaves, which breaks such a reader.
#
# 1.1 added `genes` to every omics shard row (issue #3) and
# `conflicting_lesion_groups` to every gene index row (issue #4). Both are
# additive, so 1.0 readers are unaffected.
#
# 2.0 removed `classification` and `source_tier` from the curated assertion and
# moved gene-disease validity to mirrored, attributed records (design decision
# D12: the atlas no longer authors a validity call of its own). A 1.x reader
# looking for a classification on the assertion finds none, which is a
# breaking change and so MAJOR rather than MINOR, even though the release's
# other change — the gene bundle's new `validity` object, carrying every
# mirrored record and its provenance — is purely additive on its own.
#
# 2.1 added `status`, so a programmatic consumer can read the atlas's own
# readiness without scraping `index.html`'s prose for it. Nothing existing
# changed shape or left, so this is MINOR despite following a MAJOR release —
# each release's letter is decided against what it changes, never against the
# rung before it.
#
# 2.2 added `atlas_curation` to every gene index row and every gene bundle
# (design decision D21/D22). The same release widened the published population
# from the genes carrying a curated assertion to the genes a ClinGen expert
# panel calls definitive — 1 row to 23 — which is not a schema change at all:
# more rows of an unchanged shape. The new field is what stops that widening
# being a silent one, since 22 of the 23 rows now describe a gene this atlas has
# not curated. Additive, so MINOR: a 2.1 reader keeps working.
#
# That widening is also where the rule above had to be extended, because the
# rule speaks only about fields and a release can change which *rows* appear
# without touching one. A population change inside an unchanged shape stays
# MINOR: a consumer's parsing is unaffected, since every row still carries every
# field 2.1 published and none of them has changed shape. A consumer's
# *display* is very much affected — 22 of the 23 rows describe a gene a 2.1
# reader had never seen and would render as this atlas's own curated content —
# so the warning belongs in `docs/data-api.md`, where a reader meets the file,
# rather than in a version letter that cannot carry it. Read the letter for
# "will my parser still work"; read the document for "does my page still say
# something true".
#
# 2.2 also added `genes/index.html` and `genes/<slug>.html`, taking the site
# from the single page 2.1 shipped — `index.html`, this build's first HTML
# artifact — to 25 (measured 2026-08-04 against a build of the committed
# corpus). Additive in a stronger sense than a field is: no consumer of the data
# API is affected by a page appearing beside a payload, and every page is
# checksummed into `files` like any other published byte.
#
# 2.3 adds `burden` to every gene bundle and `burden_row_count` to every
# `genes/index.json` row: the per-study rare-variant burden statistics from
# `mirrors/burden.tsv` -- 1,192 rows over 145 genes *at that version*, from one
# study. The mirror has grown twice since (1,475 rows, 150 genes, three studies);
# these figures describe what 2.3 shipped and are deliberately not updated,
# because a changelog entry that tracks the present tells a reader nothing about
# the version it names. Current counts live in `docs/data-api.md`, where
# `tests/test_docs_match_the_build.py` measures them against a real build.
#
# Additive, so MINOR — a 2.2 reader keeps working, and both keys are always
# present (an empty array and a zero for a gene no study reported), so a
# consumer never has to guard for a missing one.
#
# The display warning that accompanied 2.2's population change applies here in a
# sharper form, and `docs/data-api.md` carries it: `effect` must never be
# rendered without `effect_measure` beside it, because an odds ratio of 3.1 and
# a de novo enrichment of 3.1 are different claims sharing one column; and
# `effect` is `null` on 34 rows whose odds ratio the study published as
# infinite, where `effect_bound` is `unbounded_above` and `ci_low` carries the
# whole finding. A consumer treating that null as "not tested" would drop the
# strongest results in the data.
#
# 2.4 adds `pvalue_adjusted` and `pvalue_adjustment` to every burden object,
# and `mirrors/burden.tsv` gains its second study (PMID:34324492, CNV deletion
# burden by permutation). Additive, so MINOR: both keys are always present,
# `null` where the study published no correction.
#
# The display warning matters more than the parsing one, again. A raw p and a
# corrected p can point opposite ways -- CHD7 in PMID:34324492 is 0.0068 raw and
# 0.991 family-wise corrected -- so a consumer rendering `pvalue` alone will
# show a result as significant that the study reported as null. Where
# `pvalue_adjusted` is present it is the number to read.
# 2.5 adds `count_unit` to every burden object. Additive, so MINOR: the key is
# always present, and every row published before it carried `individuals`.
#
# It is the third display warning in a row, and the most load-bearing of them,
# because unlike the two above it changes what the *existing* keys mean.
# `n_case_carriers`, `n_cases`, `n_control_carriers` and `n_controls` are not
# comparable across rows with different `count_unit`: `individuals` counts
# people, `alleles` counts alleles -- so a person with two qualifying variants
# counts twice and the denominator is roughly twice the sample size -- and
# `de_novo_mutations` counts mutations against a denominator of *trios*. A
# consumer that divides `n_case_carriers` by `n_cases` across studies without
# reading this key is comparing three different quantities. That is not
# hypothetical: PMID:42230622 and PMID:40127276 both publish case-control rows
# for the same genes, in different units.
# 2.6 adds `independent_datasets` to every gene bundle and every
# `genes/index.json` row. Additive, so MINOR: the key is always present, and
# `families` always carries one entry per cohort family in the corpus -- a gene
# no study reported gets `not_tested` entries, never a shorter list.
#
# **It is a count of datasets, not a validity call, and the display obligation
# is that it must never be rendered as one.** `headline_confidence` is a
# mirrored ClinGen classification; this object sits beside it and says how many
# independent cohort families showed enrichment. A consumer that renders
# "0 of 2" as a verdict beside a green `definitive` chip will tell a clinician
# the data contradict the classification. They do not -- KDM6A causes Kabuki
# syndrome and shows nothing in either dataset that tested it, because burden
# tests at these cohort sizes routinely miss genes with overwhelming family
# evidence.
#
# Two consequences for any consumer: read `tested` as the denominator, never
# `len(families)`; and render `not_tested` distinguishably from
# `no_enrichment`, because collapsing them is what turns "nobody looked" into
# "somebody looked and found nothing".
#
# 2.10 adds `has_no_association_report` and `no_association_reported_by` to
# every gene bundle and browse row. Additive, and a *third axis* rather than a
# widening of `has_conflicting_evidence`, which keeps meaning exactly
# `CONTESTED`. Issue #13: a gene one submitter calls `Definitive` and another
# calls `No Known Disease Relationship` published `has_conflicting_evidence:
# false` and said nothing else, so the disagreement reached no byte. GDF1 is the
# live case. Folding the two together instead would give one laboratory's null
# result the weight of a chartered panel's refutation.
#
# 2.9 adds `cohorts.json` and a `cohorts` count. Additive, so MINOR: no existing
# key changes meaning and no payload loses one.
#
# It closes the gap 2.3 opened. Every burden row has named its sample collections
# by bare id since `burden` was published -- `["taa_cases"]`, `["ukbb"]` -- and
# nothing resolved those ids to anything: measured 2026-08-06, `curation/
# cohorts.yaml` reached zero published bytes and `Cohort.url` reached none at
# all. So every programmatic consumer had the numbers and none of the caveats
# that qualify them, including that `taa_cases` is 777 probands who do not have
# congenital heart disease. The descriptions rendered only inside a `<details>`
# on the gene page, which is a property of the HTML and not of this API.
#
# 2.7 adds `genes`, `burden_rows` and `cohort_families` to `counts`. Additive,
# so MINOR: every 2.6 key is present and unchanged.
#
# It is a correctness fix, not a convenience. Every `counts` key through 2.6
# counted a `corpus` collection, so the object described *the curation* while
# reading as a census of the build -- measured 2026-08-05 it published
# `assertions: 1, datasets: 0, functional: 0` for a site publishing 23 genes and
# 290 burden statistics, and named neither. A consumer summarising the atlas from
# this object reported an empty one. `index.html` had the same omission and is
# fixed in the same release, because `docs/data-api.md` then described the page
# as stating "every key of `counts`" -- one census in two formats, which must not
# be able to disagree. Since the doorway page (2026-10-27) the page states six of
# them, and this object all of them.
#
# So `counts` is now a census of what the build *published*, not of what the
# curation holds. The distinction matters for reading a diff against 2.6:
# `genes` is D21's published population (23), not `len(corpus.assertions)`'s
# genes (1) and not the 154 genes either validity mirror covers.
#
# `burden_rows` counts rows reaching a published bundle -- 290 -- not the 1,475
# in `mirrors/burden.tsv`. The mirror is deliberately wider than the publication
# gate, so 127 of its 150 genes publish no page; a count of the mirror here would
# advertise evidence no consumer of this API can fetch. It is exactly the sum of
# every `genes/index.json` row's `burden_row_count`, and a consumer may check it
# that way.
#
# `cohort_families`, and deliberately **not** `independent_datasets` even though
# that is the key every gene bundle carries. `counts` already has `datasets`,
# which means *omics* datasets and is 0; two adjacent keys reading `datasets: 0`
# and `independent_datasets: 3` is a misreading waiting to happen. This is the
# length of any bundle's `independent_datasets.families` array, which is the same
# for every gene by construction.
#
# 2.8 adds `admitted_by` and `asserted_by` to every gene bundle, and widens the
# published population from 23 genes to 92. Additive, so MINOR -- and the letter
# is decided the same way 2.2's was, which is the precedent that matters here
# because this release repeats its shape exactly.
#
# **MINOR was checked against the rule, not assumed.** The one candidate for
# MAJOR is `headline_confidence`, which is `null` on 16 of the 92 genes now
# published and was `"definitive"` on all 23 before. That is not a shape change:
# the field has been `string | null` since 2.0, `docs/data-api.md` has documented
# the null since then, and a parser written against 2.7 keeps working. What
# changed is which rows appear -- MINOR by the paragraph above, which was written
# for exactly this case.
#
# **The display obligation is the real cost of this release, and it is larger
# than 2.2's.** A 2.7 consumer could render `headline_confidence` verbatim and be
# right every time, because every published gene was `definitive`. Measured
# 2026-08-06 on the widened corpus that becomes 23 `definitive`, 1 `strong`,
# 9 `moderate`, 43 `limited` and 16 `null`. Three consequences, none of which a
# version letter can carry:
#
# * a consumer that treats presence in `genes/index.json` as "an expert panel
#   called this a CHD gene" is now wrong for 16 genes, which no panel graded at
#   all. `validity_state` is what distinguishes them **in that file**:
#   `"expert_curated"` on 76, `"submitter_curated"` on 16, and measured
#   2026-08-06 the `"submitter_curated"` set is exactly the `null`-headline set
#   and exactly the `gencc_agreement` set. The richer `admitted_by.authority`
#   naming the warrant is in the gene bundle, not in the index -- this comment
#   pointed a browse-row consumer at it until 2026-08-06, and no row has ever
#   carried it.
# * a consumer that renders the chip without the grade's meaning turns
#   `limited` -- a panel saying the case is not yet made -- into a weak yes.
#   The grade key is what the site itself does about this (`render.grade_key`,
#   on `methods.html#scope` since 2.18; before it a folded `grade_legend` on
#   every gene page and the gene index).
# * `asserted_by.count` is a count of institutions, never a score. Eight
#   authorities asserting a gene frequently means it sits on more commercial
#   panels, not that it is eight times better supported. D12 applies: a rank
#   derived from it would be a validity call the atlas authored.
#
# 2.11 adds `profile_datasets` and `profile_genes` to `counts`, and a
# "Developmental expression" pair of cards to `index.html`'s published list --
# the bulk developmental RNA expression layer built across Tasks 1-14, landing
# here with nothing published mentioning it yet. Additive, so MINOR: no
# existing key changes meaning and no payload loses one.
#
# Both are restricted to `published`, for the same reason `genes` and
# `burden_rows` already are (2.7's own paragraph above): `mirrors/profiles/*.tsv`
# can name a gene no external authority has admitted, the same 154-vs-92
# asymmetry `mirrors/genes.tsv` carries for burden, so a count reading the
# mirror instead of the population would advertise developmental-expression
# evidence no consumer of this API can fetch. `profile_datasets` counts the
# distinct dataset accessions named by a *published* gene's own
# `expression_profile.datasets` -- never `len(corpus.datasets)`, which already
# means something else and is unrestricted: every registered omics dataset,
# `profile`-design and `contrast`-design alike. Two keys both named "datasets"
# and both meaning something different is exactly the confusion
# `cohort_families` was kept apart from `independent_datasets` to avoid (2.7's
# paragraph again), so this one counts only what a published gene's own bundle
# actually carries.
#
# **Both were 0 when this paragraph was written and are not now.** It said "no
# `profiles` mirror has ever been committed" -- true on 2026-08-15 and false
# from 04a1077 four days later, which committed `mirrors/profiles/E-MTAB-6814.tsv`.
# Measured 2026-08-20 on the committed corpus: `profile_datasets` is 1 and
# `profile_genes` is 92.
#
# The reason the sentence gave still holds, in a sharper form. A real build
# cannot tell these two figures from the ones beside them: `profile_genes` is
# 92 and so is `genes`, `profile_datasets` is 1 and so is `datasets`. Two
# figures that are equal today are one figure to every test (CLAUDE.md section
# 4.30), so `profiles.profile_census`'s own tests build a fixture that
# separates them rather than trusting a real build to.
#
# 2.12 adds `order` to every `stages` entry in `datasets.json`, and changes the
# order of `expression_profile.datasets[].stages[]` in all 92 gene bundles from
# alphabetical to chronological. MINOR -- but the second half is a correction to
# already-published data, so the letter was checked against the rule above
# rather than assumed.
#
# `order` is the curated chronological position of one of a dataset's own stage
# tokens: 1-based, ascending, unique within a dataset (PRF011, and PRF012
# refuses an `order` that contradicts its own `wpc`). Required and always
# present, so a 2.11 parser is unaffected and a 2.12 reader never guards for a
# missing key. It exists because `wpc` cannot do this job: `wpc` is null for
# every post-natal stage -- 8 of E-MTAB-6814's 21 -- so the tokens that most
# need ordering carry no number at all, and what the build sorted on instead was
# the token string.
#
# **The reordering is the release, and it is a defect being fixed rather than a
# feature being added.** Through 2.11 a gene bundle's stage array was sorted
# alphabetically: `4 week post conception` published eighth, after `19 week post
# conception`, and `elderly` second of the eight post-natal stages.
# Deterministic, reproducible, byte-identical between builds, and wrong -- a
# developmental time series published in dictionary order.
#
# **Through 2.11 means under 2.11 and no other version.** This paragraph read
# "it survived three releases" until 2026-08-21, in the same breath as the
# sentence above that names the single version. Measured: the profiles mirror
# was first committed on 2026-08-19 (`git log --diff-filter=A -- mirrors/
# profile*`), at which commit `SCHEMA_VERSION` was already 2.11; before it,
# every gene bundle's `expression_profile` was `{"datasets": []}` and no stage
# array was published at all. The order was wrong for one schema version,
# promoted to `main` on the day the mirror landed and corrected the next.
#
# It survived that version because nothing asserted a chronology and no
# published artifact drew one -- which is the part worth keeping, and does not
# need an inflated count to be alarming.
#
# MINOR by the rule, for the reason 2.2's population change was MINOR: no field
# is added, removed or reshaped, every entry still carries every key 2.11
# published, and a parser is unaffected. The paragraph 2.2 added for *rows*
# applies to *positions* -- what changed is which entry sits where, not what an
# entry is.
#
# **The display obligation is heavier than the parsing one, again, and it is
# unusual: a correct consumer's output changes.** Anything that plotted `stages`
# in array order -- the obvious thing to do with a developmental series -- drew
# an axis in the wrong sequence and will now draw a different picture. That is
# the correction arriving, not a regression to investigate. A consumer that
# worked around it by sorting the tokens itself must *remove* the workaround,
# because sorting the tokens is precisely the bug. No version letter can carry
# that, so `docs/data-api.md` states it where a reader meets the array.
#
# `datasets.json`'s own `stages` array **is** sorted by the build, on `order`,
# since 2026-08-21. It was not: the file is serialised generically from the
# curated record (`literature._dump`), so its sequence was the curation file's
# declaration order and agreed with `order` only because the curator typed it
# that way. Measured -- moving the `elderly` block to the top of the YAML while
# leaving `order: 21` alone published `elderly` first with `validate` at the
# exact documented baseline. PRF011 and PRF014 constrain the values of `order`,
# never the sequence a file lists them in, so the guarantee lives in the
# emitter where it holds by construction.
#
# **That array moved in this release too, and measurement is the only reason
# this sentence is here.** A diff of two real builds (2.11 against 2.12) shows
# `datasets.json` changing in *two* ways, not the one that was expected: every
# stage gains `order`, and the eight post-natal entries come back in a
# different sequence -- `neonate` first rather than `adolescent` -- because the
# curator renumbered and re-sequenced that block in the same commit. The token
# *set* is unchanged and every `wpc` is unchanged; the sequence is not. It was
# first recorded as "gained exactly one key, with identical stage tokens",
# which is true of the set and false of the order (CLAUDE.md section 4.24).
# ---------------------------------------------------------------------------
# 2.13 -- a gene's retired HGNC symbols become search terms.
#
# `mirrors/genes.tsv` gained `prev_symbols` with issue #33, for the join: the
# resolver in `chd_atlas/genes.py` needs it so a source still writing a
# withdrawn name lands on the right gene. That change published nothing. This
# one does -- `search/index.json.gz` now carries those symbols in a gene
# record's `terms`.
#
# MINOR, by the same rule as 2.2's population change and 2.12's reordering: no
# field is added, removed or reshaped, every record still carries every key
# 2.12 published, and `terms` was already a variable-length array of strings a
# consumer must not depend on the length of.
#
# Measured over the committed corpus, one real build against the previous one:
# 35 of the 92 published gene records gain a term, 46 terms in total, and
# **zero** of them already named anything else in the index -- so no query that
# returned one gene before now returns two by accident. `TBX5` gains `HOS`,
# `VEGFA` gains `VEGF`, `ZIC3` gains `HTX1`.
#
# **One term is deliberately published on two genes: `MADH7`, on SMAD6
# (HGNC:6772) and SMAD7 (HGNC:6773), which both retired it.** `genes.py`
# *refuses* that symbol -- `Resolution.AMBIGUOUS`, no id -- and the two are not
# in conflict, because they answer different questions. The resolver is asked
# "which single gene does this string denote?", where guessing attributes one
# gene's evidence to another. Search is asked "what might this reader be looking
# for?", where showing both is the only honest answer and the reader picks.
# Filtering ambiguous retired symbols out of the index to make the two agree
# would make `MADH7` unfindable, which serves nobody.
#
# 2.14 adds `burden_combined` to every gene bundle and `burden_combined` to
# `counts`: one study's own combination of two of its burden arms, from
# `mirrors/burden_combined.tsv`. Additive, so MINOR: the key is always present,
# `[]` where the gene has none. Measured on the committed corpus: 75 objects over
# 25 published genes, all PMID:40127276's JL-FDR meta-analysis, against 90
# mirrored rows for 30 registered genes. Its `local_fdr` is a local false
# discovery rate -- smaller than its own `pvalue` for about half that study's
# panel -- and a consumer must never present it as a corrected p.
#
# 2.15 publishes `scope.json` and adds `scope_terms` to `counts`: the CHD scope
# registry, every term with the authority it is attributed to and, for a term
# admitted on Orphanet's phenotype frequencies (issue #8 part 1), the warrant
# re-derived from the mirrors. Additive, so MINOR: no field of 2.14 is removed
# or reshaped. It lands with a population change -- 92 -> 119 genes, the
# Orphanet-scoped syndromes (Noonan, Costello, Adams-Oliver, Myhre,
# Rubinstein-Taybi and their kin) plus NOTCH1 on its own ClinGen term -- which
# carries the display obligation: a Noonan grade is a grade for the syndrome,
# and a consumer showing it must say why the syndrome counts, which is what
# `scope.json`'s `orphanet` object is for.
#
# 2.16 adds `submitted_on` and `sgc_id` to every record of a bundle's
# `validity.records`, `null` on a ClinGen record. Additive, so MINOR. Until it,
# a GenCC classification published no date at all -- 470 of 591 records on
# 2026-10-25 -- and no identifier: `(source, disease, moi, submitter)` is not
# unique, and PKD1L1 publishes two Ambry submissions that share it.
#
# 2.17 adds `identity` to every gene bundle: the gene's approved name, locus,
# Ensembl and NCBI Gene ids, UniProt accession and MANE Select transcript from
# `mirrors/genes.tsv`, every key always present. Additive, so MINOR. All six
# were mirrored for every published gene; before 2.17 `locus`, `ensembl_gene`,
# `ncbi_gene`, `uniprot` and `mane_select` reached no published byte, and
# `name` only search's `terms`.
#
# 2.18 adds `summary_caveat` to every record of `publications.json` (and so
# to `featured.json`'s embedded publications), `null` where the curator
# recorded none. Additive, so MINOR. It is the one sentence a gene page's
# card prints beneath a study's results, and the gene index under its strip
# legend; on the committed corpus only PMID:34324492 carries one, naming the
# 777 thoracic aortic aneurysm probands in its case set.
#
# It also publishes `methods.html`, the site's reading notes: an HTML page,
# additive in the sense 2.2 records, so no bump of its own.
#
# 2.19 (owner decisions of 2026-10-26, one bump for the batch). MINOR, every
# change additive or a population change in an existing field:
#
# * The structural set Orphanet scoping reads gains six lesions HPO files
#   outside Abnormal heart morphology (`orphanet_scope.ADDED_BRANCHES`). No
#   field changes; ten `scope.json` warrants gain pulmonic stenosis among their
#   `phenotypes`, and no term or gene enters or leaves scope.
# * Every record of `publications.json` (and `featured.json`'s embedded ones)
#   gains `correction_genes` and `genome_wide_threshold`, `null` where not
#   curated: the gene count a study's own correction ran over, and the p it
#   states for genome-wide significance. PMID:40127276's `tests_reported` goes
#   from 248 to `null`, because 248 counted its panel genes, not comparisons.
# * Every `scope.json` term gains `syndrome`, a curated boolean, `true` on the
#   nine ClinGen-scoped syndromes published genes are graded for.
# * `counts` gains `burden_surviving_genes`: the published genes with a row
#   surviving its own study's correction, the genes whose card lists one.
# * Every `sources.json` entry gains `file_version` and `file_date`, what the
#   file states about itself, `null` but for the two Orphadata products.
# * Every `cohorts.json` record gains `overlaps_with`, the declared overlaps
#   (`ddd` and `decipher`), `[]` elsewhere. A population change too: the two
#   studies citing them become one cohort family, so `counts.cohort_families`
#   and every bundle's `independent_datasets.families` go from 3 to 2.
# * Every `independent_datasets.families` entry (bundle and `genes/index.json`
#   row) gains `members`: one `{study, state}` per study in the family, sorted
#   by study, each study's own state by the family rule over its rows alone. A
#   family of two studies made the family's state no longer every member's.
# * `cohorts.json` gains seven records (owner decision R12, 2026-10-26): the
#   CHD case series PMID:34324492's S2 Table identifies by PMID with 50 or more
#   cases, which every one of that study's rows now names in `case_cohorts`.
#   A population change: `counts.cohorts` goes from 13 to 20, and
#   `counts.cohort_families` stays 2, since no other study names one of them.
# * Every object of a bundle's `burden` array gains `n_chd_case_carriers`
#   (owner decision R2, 2026-10-26): how many of `n_case_carriers` are CHD
#   cases, where a study whose case set also holds another condition counted
#   them apart; `null` elsewhere. Only PMID:34324492's S8 Table gives one, on
#   the 12 mirrored rows (10 published) for the genes it lists.
# * Every `cohorts.json` record and every `scope.json` term gains
#   `review_status` and `review_note` (owner's answer, 2026-10-27), both `null`
#   where none is recorded: `pending_review`, a question about the record that
#   this atlas leaves open for an expert (on the two cohorts marked, the
#   sources they rest on conflict; on the two scope terms, a clinical
#   judgement), or `not_accessible`, a source the record rests on could not be
#   read, with a note stating what is unresolved (`models/review.py`).
#
# 2.20 (single-cell cell-type expression, Farah et al. 2024). MINOR, additive:
#
# * Every gene bundle gains `cell_type_expression`, always present, and
#   `{"datasets": []}` for a gene no cell-type dataset covers. A profile of
#   cell types, never a contrast: a low or absent value is not evidence
#   against a gene.
# * `counts` gains `cell_type_genes` and `cell_type_datasets`, the published
#   genes with a non-empty `cell_type_expression` and the distinct datasets
#   those bundles name (`cell_types.cell_type_census`). Both are 0 until
#   `mirrors/cell_type_profiles/` is committed.
# * Every `datasets.json` record gains `cell_types` (each `token`, `label`,
#   `compartment`) and `regions`, both `[]` on every design but
#   `cell_type_profile`.
# * `design` gains the value `cell_type_profile`, and `archive` the value
#   `ucsc_cell_browser` (`vocab.DatasetDesign`, `vocab.Archive`).
# * A dataset id may be a lowercase slug (`hoc`, `heart-atlas`), the UCSC Cell
#   Browser's name for a dataset, as well as an archive accession: in
#   `datasets.json`'s `id` and in the `dataset` an assertion's evidence names
#   (`identifiers.AccessionId`). No curated record uses one yet.
SCHEMA_VERSION: Final = "2.20"

# What `status` publishes today. A literal rather than something derived from
# the corpus, unlike every field in `counts`: there is no measurement of "is
# this atlas ready for clinical use" to compute, only a decision to record, and
# recording it as a constant is what stops a future release forgetting to flip
# it — the two-word status and the paragraph on `index.html` explaining what it
# concretely means are edited in the same commit, or neither is trustworthy.
STATUS: Final = "in-development"


def source_commit(root: Path) -> str | None:
    """The commit the build was made from, or None when there is not one.

    None rather than a raise: building from an unpacked tarball is legitimate
    and should produce a complete site, just one that cannot state its
    provenance. The same answer covers a `root` that does not exist and a
    checkout with no commits yet — `git -C` failing to change directory records
    exactly as much provenance as a missing `.git` does, so neither earns a
    branch of its own here.

    The `returncode` check is load-bearing rather than defensive, which is not
    obvious: in a repository with no commits `git rev-parse HEAD` exits 128 but
    prints the literal string `HEAD` on *stdout* (measured, git 2.50). Without
    the check this would publish `"source_commit": "HEAD"` — a provenance claim
    resolving to nothing, which is worse than the null a consumer can test for.

    `check=False` is not a second decision: `check=True` raises
    `CalledProcessError`, a `SubprocessError`, which the handler below already
    catches, so the two spellings behave identically. `timeout` is a real one —
    without it a git that blocks on a dead network mount or a credential prompt
    hangs the build with nothing to read.

    The toplevel comparison is what makes "not a checkout" mean it. `git -C X`
    searches *upward* from X, so a `--root` that is merely unpacked inside some
    unrelated repository resolves that repository's HEAD, and this published a
    real, verifiable sha naming a commit containing none of the built data —
    strictly worse than the `"HEAD"` the returncode check exists to prevent,
    because a false provenance claim that resolves is one nobody thinks to
    check. Two environment variables do the same thing more quietly: `GIT_DIR`
    and `GIT_WORK_TREE` override `-C` outright, so they are cleared rather than
    trusted. Reproduced before fixing, and pinned by
    `test_a_root_inside_an_unrelated_checkout_has_no_provenance`.
    """
    # `-C` still points git at `root`; `--show-toplevel` then reports which
    # repository it actually found, which is the answer that has to be checked.
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    try:
        located = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
            env=environment,
        )
    # OSError is git absent from PATH, which a build must survive for the same
    # reason a tarball build must.
    except (OSError, subprocess.SubprocessError):
        return None
    if located.returncode != 0:
        return None

    # Line-based, not whitespace-based. `--show-toplevel` prints a path, and a
    # path may contain spaces — `~/My Projects/atlas` is ordinary on macOS.
    # Splitting on whitespace made that three tokens instead of two and returned
    # `None`, so a perfectly valid checkout published no provenance and did it
    # by failing into the same answer that means "not a checkout at all".
    # Reported by review on #7 and reproduced before fixing.
    #
    # The commit is the last line and never contains whitespace; everything
    # before it is the path, rejoined, so even a newline in a directory name
    # compares correctly rather than silently answering `None`.
    lines = located.stdout.splitlines()
    if len(lines) < 2:
        return None
    toplevel, commit = "\n".join(lines[:-1]), lines[-1]
    try:
        if Path(toplevel).resolve() != root.resolve():
            return None
    except OSError:
        return None
    # `or None`, because "" is falsy but still reads as a recorded commit to a
    # consumer testing the key against null.
    return commit.strip() or None


def write_manifest(
    corpus: Corpus,
    emitter: Emitter,
    commit: str | None,
    counts: Mapping[str, int],
) -> None:
    """Write `manifest.json` and seal the emitter. The build ends here.

    `emitter.checksums` is complete only once every other builder has run, so
    this must be the last thing a build emits. Sealing is what makes that a rule
    rather than an intention — see `Emitter.seal`.

    `counts` carries the figures this module cannot derive from a `Corpus`: the
    published gene population, the burden rows reaching a bundle, the number
    of cohort families, and the developmental-expression layer's own
    `profile_datasets`/`profile_genes` (`profiles.profile_census`) and the
    single-cell layer's `cell_type_datasets`/`cell_type_genes`
    (`cell_types.cell_type_census`). They are
    passed in rather than recomputed for the reason `runner.py` threads
    `published` to three builders — the front page, the browse payload and this
    manifest must state one census, and the only way they cannot disagree is to
    count the same objects. Merged over the corpus counts rather than nested
    under a key of their own, because a consumer reading `counts` wants the
    census, not a lesson about which half of it came from where.
    """
    corpus_counts = {
        "assertions": len(corpus.assertions),
        # The registry `cohorts.json` publishes entire, not the subset the
        # published rows cite. Those were 20 and 20 on 2026-10-27; the key is
        # named for the file so the two cannot be confused if they ever part.
        "cohorts": len(corpus.cohorts),
        "datasets": len(corpus.datasets),
        "featured": len(corpus.featured),
        "functional": len(corpus.functional),
        "phenotypes": len(corpus.phenotypes),
        "publications": len(corpus.publications),
        # The whole registry `scope.json` publishes: every disease term in
        # scope, on whichever authority's warrant. 105 since 2026-10-25.
        "scope_terms": len(corpus.chd_scope),
    }
    # A build count silently shadowing a corpus count would republish one figure
    # under another's name, and `sort_keys` would leave the output looking
    # perfectly ordinary. `datasets` is the live hazard: it means *omics*
    # datasets here and a caller reaching for the cohort-family count is one
    # plausible keystroke from overwriting it. `raise`, never `assert` -- `-O`
    # strips `assert`, and this is a guard on a mistake no test would otherwise
    # see.
    collisions = sorted(set(counts) & set(corpus_counts))
    if collisions:
        raise ValueError(
            f"build counts {collisions} would shadow the corpus counts of the same name; "
            f"rename them rather than overwriting a figure derived from the corpus"
        )

    emitter.write_json(
        "manifest.json",
        {
            "schema_version": SCHEMA_VERSION,
            "source_commit": commit,
            "status": STATUS,
            "counts": {**corpus_counts, **counts},
            # A snapshot, not an ordering. Sorting here was specified and then
            # dropped: `encode_json` passes `sort_keys=True`, which orders keys
            # at every level, so an unsorted mapping encodes to the same bytes.
            # Measured rather than reasoned — restoring the sort changes no
            # published byte and no test can see it — because a sort that buys
            # nothing has already been documented twice in this build as the
            # thing that keeps it reproducible.
            #
            # The copy is equally undetectable and is kept as hygiene, which is
            # the honest description: the payload is encoded before `write_json`
            # records anything, so aliasing `checksums` here would not put
            # `manifest.json` into its own `files` today. It would if that
            # ordering ever changed, and a mapping handed out is one a caller
            # can mutate.
            "files": dict(emitter.checksums),
        },
    )
    emitter.seal()
