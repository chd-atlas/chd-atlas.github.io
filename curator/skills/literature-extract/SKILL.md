---
name: literature-extract
description: Use when a gene's included papers need findings written, one per distinct result, each number with a locator, and a Publication record added.
---

# literature-extract

Read `curator/curator-rules.md` first. Model tier: a mid-tier model (for example Sonnet).

## Inputs

- A gene symbol with no `pending` triage entries (the `literature-triage` skill is done).
- The model name you run as, for `extracted_by`.

## Steps

1. List the papers to extract:

   ```bash
   .venv/bin/chd-atlas literature abstracts <SYMBOL> --decision included
   ```

   Read each paper in full where PubMed Central or Europe PMC has open-access text. Extract
   only the paper itself, not results it cites from other papers.
2. For each paper write one finding per distinct result under `findings:` in
   `curation/literature/<SYMBOL>.yaml`. Two designs, or two separate cohorts, are two
   findings. Choose the design (`case_report`, `family_segregation`, `cohort_screen`,
   `rare_variant_association`, `common_variant_association`) by what the study observed.
3. Take the id: read the `LIT` counter under `prefixes` in `curation/.id_registry.yaml`,
   use the next ordinal as `CHDA:LIT:<7 digits>`, and set the counter to the last ordinal
   you used, in the same edit.
4. Fill the counts by design. Rare-variant designs need `count_unit` (`probands` or
   `families`), `screened`, `carriers` (`de_novo`, `inherited`, `origin_unknown`) and
   `variant_classes`. A `rare_variant_association` also needs a `comparator`
   (`control_cohort` or `mutation_model`); a case series is a `cohort_screen`. A
   `common_variant_association` needs `cases` and `controls`, takes `p_value` and
   `odds_ratio` when reported, and has no carrier counts. If the study's rows are already
   in `mirrors/burden.tsv`, set `burden_study` to its PMID.
5. Add a `locators` entry for every number (`kind: in_text` with `element`, or
   `kind: supplementary` with `file`, `sheet`, `row`). Set `cohorts` from
   `curation/cohorts.yaml`, and `shares_probands_with` to the ids of findings that share
   probands with this one, or `unknown`. The curator checks independence by reading this
   field in the `review` listing.
6. A paper you cannot read stays in. Record only what the abstract supports, use the
   locator element `Abstract`, and set `review_status: not_accessible` with a
   `review_note` saying what could not be read.
7. Where the paper disagrees with a source the finding rests on, set
   `review_status: pending_review` and state only the measured disagreement in
   `review_note`.
8. Set `extracted_by: agent:literature-extract/<model>`. Leave `verified_by`,
   `approved_by` and `approved_on` unset.
9. Add each paper's record the first time it gets a finding, using the paper's first
   finding for the study type:

   | Design | StudyType |
   |---|---|
   | `case_report` | `case_report` |
   | `family_segregation` | `family_linkage` |
   | `cohort_screen` | `cohort_sequencing` |
   | `rare_variant_association` | `case_control` |
   | `common_variant_association` with `genotyping_array` | `gwas` |
   | `common_variant_association` otherwise | `candidate_variant` |

   ```bash
   .venv/bin/chd-atlas literature add-publication <PMID> --study-type <StudyType>
   ```

   It exits 1 if the PMID is already in `curation/publications.yaml`; then there is
   nothing to add.
10. Run `.venv/bin/chd-atlas validate` and fix every error it names.

An illustrative finding (every value is made up):

```yaml
- id: CHDA:LIT:0000001
  publication: PMID:10000001
  design: case_report
  methods: [targeted_sequencing]
  count_unit: probands
  screened: 1
  carriers: {de_novo: 1, inherited: 0, origin_unknown: 0}
  variant_classes: [snv_indel]
  phenotypes: [HP:0001631]
  lesion_groups: [septal]
  syndromic: isolated
  shares_probands_with: unknown
  locators:
    - kind: in_text
      element: Table 1
  summary: One proband with an atrial septal defect carried a de novo missense variant.
  extracted_by: agent:literature-extract/<model>
```

## Done when

- Every included paper has a finding and a Publication record: validate reports no
  `LIT003` or `LIT011`.
- `.venv/bin/chd-atlas validate` reports no error. Hand off to the `literature-verify`
  skill, then the curator runs `.venv/bin/chd-atlas literature review <SYMBOL>`.

## Never

- Never state more than the paper does, or a number without its locator and unit.
- Never re-classify a variant or call a gene definitive.
- Never fill `verified_by` or `approved_by`, or run `literature review --approve`.
- Never invent an id outside the counter, a PMID, a cohort or an HPO term.
