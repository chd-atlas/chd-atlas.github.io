# House rules for every curator skill

Each rule names the code that enforces it where one does. `SCHEMA001` is a per-record
model rule in `src/chd_atlas/models/gene_literature.py`. A rule with no code is a
convention only a person can check.

## Principles

- A finding is evidence, never a gene-disease validity classification. Validity stays
  mirrored. Do not touch `curation/assertions/`, `mirrors/` or any file your skill's Steps
  do not name.
- Record what the authors claim. "Carrier" means the authors' own qualifying variant.
  Never re-classify a variant.
- Every number has a locator and a unit. A finding with no locator is refused
  (`SCHEMA001`); `count_unit` is required for the rare-variant designs (`SCHEMA001`).
  That each number has its own locator is a convention.
- "Unknown" is a value. Where a source disagrees with its paper, or cannot be read, set
  `review_status` and `review_note` (`SCHEMA001` refuses one without the other). State
  the measured disagreement only. Never guess which side is right, never drop the record.
- A summary is one or two factual sentences and is never stronger than the paper.
- Validators report; the build refuses. Run `.venv/bin/chd-atlas validate` before every
  hand-off.

## Sources

- PMIDs come only from `chd-atlas literature search` and `abstracts`, never from memory.
  A finding citing a PMID not included in the file's triage is refused (`LIT004`).
- Preprints and reviews are excluded by the query. One that slips through is excluded by
  hand with `not_peer_reviewed` or `review_article`.
- Never write the NCBI key into any file. It is optional and read from `NCBI_API_KEY`.

## Files and ids

- The file is `curation/literature/<SYMBOL>.yaml`, named by the gene's approved symbol
  (`LIT009`). `search` accepts an old name but writes the approved one.
- A finding id is `CHDA:LIT:<7-digit ordinal>`, taken from the `LIT` counter in
  `curation/.id_registry.yaml`, which you bump in the same edit. An id beyond the counter
  is refused (`ID002`); an id repeated within one file is refused (`SCHEMA001`), and one
  used in two files is refused (`ID003`).
- Write every value out on the record it belongs to. A literature file with a YAML anchor
  (`&a`), alias (`*a`) or merge key (`<<`) is refused whole (`YAML002`), because `<<: *a`
  would copy one finding's approval onto a finding nobody approved.
- Every cited PMID needs a record in `curation/publications.yaml`, added with
  `add-publication` (`LIT011`). A cohort must be in `curation/cohorts.yaml` (`LIT012`).
- A literature file for a gene outside the publication gate reaches no page (`LIT013`,
  a warning).

## Triage

- Every hit is triaged once (`LIT005`), to a search the file records (`LIT010`). Each
  search's `new_hits` must equal the triage entries attributed to it (`LIT001`). Never
  change a `pmid` or a `search` on a triage entry.
- A hit still `pending` is a warning (`LIT002`). A triage hand-off leaves none; the
  search hand-off ends with every new hit pending, by design.
- Exclude only for a reason the abstract shows. When unsure, include.
- An excluded paper needs a `reason` (`SCHEMA001`). `same_data_as` needs `same_as`, which
  must name a PMID included in the same file (`SCHEMA001`, `LIT006`). `other` and
  `full_text_not_accessible` need a `note` (`SCHEMA001`).
- A paper that cannot be read is included, and its findings record what the abstract
  supports, marked `not_accessible`. Where the abstract supports no finding, the curator
  excludes the paper as `full_text_not_accessible`, with a `note` naming what the abstract
  lacks (no abstract, no count, no carrier's heart lesion). It is a gap in what was read,
  not a judgement on the paper.
- DNA methylation or RNA expression measured in human samples is not human genetic
  evidence: exclude it as `not_human_genetics`. A paper that also reports a variant in
  the gene is included, and only its variant result is extracted.
- An erratum or corrigendum is excluded as `other` with `note: erratum to PMID:<the
  corrected paper>`. The extractor reads it with the corrected paper and writes a number
  it corrects as corrected, with a locator naming the erratum.
- A negative result: a CHD cohort screened for the gene with no carrier is included (its
  finding is a `cohort_screen` with every carrier count 0); a single patient or family in
  whom the gene was tested and ruled out, the cause lying elsewhere, is excluded as
  `no_finding_in_gene`.
- Cell, protein or in-silico work on variants that reports no carrier of its own is
  `not_human_genetics`, even when the variants came from patients in earlier papers.
- A trait that is not a congenital heart lesion (atrial fibrillation alone, adult heart
  structure, myocardial infarction) is `not_chd`; a review whose only new analysis is of
  such a trait is `review_article`.
- An included paper needs at least one finding (`LIT003`). Never write a finding only to
  clear `LIT003`: say why the paper has none, and the curator re-triages it.

## Findings

- One finding per distinct result: a paper reporting two designs or two separate cohorts
  yields one finding each (`LIT004` refuses a finding for a paper that is not included).
- The design fixes the counts. Carriers summed across origins never exceed `screened`; a
  `case_report` screens exactly one; a `rare_variant_association` has a comparator other
  than `none`; a `common_variant_association` has `cases` and `controls` and no carrier
  counts; the genotyping methods appear only with `common_variant_association`. All
  `SCHEMA001`. An infinite or NaN `p_value` or `odds_ratio` (`.inf`, `.nan`) is refused
  (`SCHEMA001`): write the bound the paper gives, or leave the field out.
- `carriers` and `screened` are counted in `count_unit`. With `probands` a carrier is a
  person. With `families`, `screened` is families and each origin count is carrier
  families, so how many relatives carry the variant goes in the `summary` with a locator,
  not in a count.
- Origin. A carrier family is `inherited` when the paper shows the variant passing from a
  carrier parent to a child, `de_novo` when it arose in the family's affected child
  (absent in both parents), and `origin_unknown` only when no transmission is shown. A
  proband is `de_novo` if absent in both parents, `inherited` if found in a parent, and
  otherwise `origin_unknown`.
- A paper that names no detection method has `methods: [not_stated]`, alone (`SCHEMA001`).
  Never write the nearest method.
- `methods` lists every method that found the variant in a counted carrier, proband or
  relative (capillary sequencing of relatives counts). Never a step that only confirms a
  variant already found, nor one that found nothing.
- A pooled analysis of cohorts the paper also reports one by one is not a finding of its
  own: one finding per cohort.
- A number the paper does not print (a total summed from its stages) is written only with
  the summary saying how it was derived.
- A finding needs a congenital heart lesion. Carriers with only rhythm or conduction
  disease, a cardiomyopathy (such as LVNC) or no heart defect get no finding of their
  own; name them in the summary of the paper's CHD finding where it matters. A paper left
  with no finding is re-triaged `not_chd` by the curator.
- `phenotypes` are the carriers' own congenital heart lesions, each as the most specific
  HPO term the paper names (the paper says secundum: the secundum term). The cohort's
  lesions only for a screen with no carrier, or where the paper never gives the carriers'
  lesions; then the summary says the lesions are the cohort's. No extracardiac features.
- `lesion_groups` agree with the `phenotypes` that `curation/phenotypes.yaml` maps to a
  group (`REF007`, `REF009`); a phenotype it does not map is not checked. `other` is only
  for a structural lesion with no group (patent ductus arteriosus, Ebstein anomaly, cor
  triatriatum) and for CHD whose type the paper does not give (HP:0001627).
- `burden_study` links a study `mirrors/burden.tsv` already holds for this gene; with no
  such row the finding is refused (`LIT014`).
- Independence: list `cohorts` from the registry. Set `shares_probands_with` to the
  finding ids that share probands. `[]` asserts that no proband is shared with any other
  finding: write it only where the paper states its cohort is new. Otherwise write
  `unknown`. Leaving the field out is refused (`SCHEMA001`); not guessing is a convention
  nothing enforces. A named id must exist and must not be the finding's own (`LIT007`);
  it may name a finding not yet approved. The curator checks independence by
  reading this field in the `chd-atlas literature review` listing; `diff` does not
  compare it. The blind verifier is the one exception to naming ids: it always writes
  `unknown`, because its ids are local to its own file.

## Attribution and approval

- `extracted_by` and `verified_by` are `agent:<skill>/<model>`. The verifier is a
  different model from the extractor, in a fresh context.
- `approved_by` is a human's name. A value starting `agent:` is refused by the model
  (`SCHEMA001`) and by `chd-atlas literature review`. `approved_by`, `approved_on` and
  `approved_digest` are set together (`SCHEMA001`).
- `approved_digest` is a hash of the finding as it was approved; `literature review
  --approve` writes it. Editing any field it covers (the id, counts, summary, locators,
  anything but the approval, `verified_by` and the review mark) leaves the approval
  stale: the validator reports `LIT015` and the build refuses until the three approval
  fields are deleted and a person approves again. Never write or copy a digest by hand.
- An agent never runs `literature review --approve` and never writes `approved_by`. That
  is a convention: the command cannot tell a person's name from an agent's invented one.
- Approval is not dated before the search that found the paper (`LIT008`).
- A review mark (`pending_review` or `not_accessible`, with its `review_note`) is
  published with the finding. It does not withdraw an approval.
