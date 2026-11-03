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
  must name a PMID included in the same file (`SCHEMA001`, `LIT006`). `other` needs a
  `note` (`SCHEMA001`).
- A paper that cannot be read is not excluded. Include it, and record what the abstract
  supports, marked `not_accessible`.
- An included paper needs at least one finding (`LIT003`).

## Findings

- One finding per distinct result: a paper reporting two designs or two separate cohorts
  yields one finding each (`LIT004` refuses a finding for a paper that is not included).
- The design fixes the counts. Carriers summed across origins never exceed `screened`; a
  `case_report` screens exactly one; a `rare_variant_association` has a comparator other
  than `none`; a `common_variant_association` has `cases` and `controls` and no carrier
  counts; the genotyping methods appear only with `common_variant_association`. All
  `SCHEMA001`. An infinite or NaN `p_value` or `odds_ratio` (`.inf`, `.nan`) is refused
  (`SCHEMA001`): write the bound the paper gives, or leave the field out.
- `carriers` and `screened` are counted in `count_unit`. With `probands` or `individuals`
  a carrier is a person. With `families`, `screened` is families and each origin count is
  carrier families, so how many relatives carry the variant goes in the `summary` with a
  locator, not in a count.
- `lesion_groups` agree with the `phenotypes` that `curation/phenotypes.yaml` maps to a
  group (`REF007`, `REF009`); a phenotype it does not map is not checked.
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
