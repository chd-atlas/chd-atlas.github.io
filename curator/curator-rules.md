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

- PMIDs come only from `chd-curator search` and `chd-curator abstracts`, never from memory.
  A finding citing a PMID not included in the file's triage is refused (`LIT004`).
- Preprints and reviews are excluded by the query. One that slips through is excluded by
  hand with `not_peer_reviewed` or `review_article`.
- Never write the NCBI key into any file. It is optional and read from `NCBI_API_KEY`.

## Files and ids

- The file is `curation/literature/<SYMBOL>.yaml`, named by the gene's approved symbol
  (`LIT009`). `chd-curator search` accepts an old name but writes the approved one.
- A finding id is `CHDA:LIT:<7-digit ordinal>`, taken from the `LIT` counter in
  `curation/.id_registry.yaml`, which you bump in the same edit. An id beyond the counter
  is refused (`ID002`); an id repeated within one file is refused (`SCHEMA001`), and one
  used in two files is refused (`ID003`).
- Write every value out on the record it belongs to. A literature file with a YAML anchor
  (`&a`), alias (`*a`) or merge key (`<<`) is refused whole (`YAML002`), because `<<: *a`
  would copy one finding's approval onto a finding nobody approved.
- Every cited PMID needs a record in `curation/publications.yaml`, added with
  `chd-curator add-publication` (`LIT011`). A cohort must be in `curation/cohorts.yaml` (`LIT012`).
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

The include rule, the reason table and the unreadable-paper rule below moved here from curator's `literature-triage` skill (its steps 3 and 4), so every triage role cites one copy.

Include when the paper reports human genetic evidence on this gene in CHD. When unsure,
include. A negative result is included when it is a CHD cohort screened for the gene with no
carrier (its finding is a `cohort_screen` with every carrier count 0). Exclude only for a reason the abstract shows:

| reason | use when the abstract shows |
|---|---|
| `not_peer_reviewed` | a preprint, thesis, a conference or meeting abstract, or other unreviewed item |
| `review_article` | a review or commentary with no new data; a review whose only new analysis is of a trait that is not a congenital heart lesion |
| `not_human_genetics` | animal, cell or protein work with no human genetic finding; or DNA methylation or RNA expression measured in human samples with no variant in the gene reported (a paper that also reports a variant in the gene is included); cell, protein or in-silico work on variants that reports no carrier of its own, even when the variants came from patients in earlier papers |
| `no_finding_in_gene` | the gene named only as background, or nothing reported about a variant in it (a gene among a study's reported results is a finding: include it); or a single patient or family in whom the gene was tested and ruled out, the cause lying elsewhere (a CHD cohort screened with no carrier is included, not this) |
| `not_chd` | no congenital heart phenotype; a trait that is not a congenital heart lesion (atrial fibrillation alone, adult heart structure, myocardial infarction) |
| `same_data_as` | the same patients and the same result reported again (add `same_as`); a paper with new results on the same patients is included, and its finding names the other in `shares_probands_with` |
| `other` | none of the above (add `note` saying why); an erratum or corrigendum is excluded here with `note: erratum to PMID:<the corrected paper>` |

A paper whose abstract is missing or whose full text you cannot reach is not excluded
for that. Include it; the extract step records it as `not_accessible`. The one
exception is an erratum, which PubMed usually prints with no abstract and the corrected
paper's title. Exclude it only when you know it is one (its title says so, or the
curator told you); a hit with no abstract whose title repeats another hit's is
included and named in your hand-off as unsure.

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
  reading this field in the `chd-atlas literature review` listing; `chd-curator diff` does not
  compare it. The blind verifier is the one exception to naming ids: it always writes
  `unknown`, because its ids are local to its own file.

### Variants

A finding's `variants` list the variants its paper names, each placed on the gene's MANE
Select frame, which `chd-atlas literature frame <SYMBOL>` prints. `LIT016`, `LIT017`,
`LIT018`, `LIT019`, `LIT020` and `LIT021` check an entry's `hgvs_c`, `hgvs_p` and `grch38`
against that frame. `chd-curator normalise` writes `hgvs_c`, `grch38` and `normalisation`
of a `mane` variant, so those checks matter once it has run; a `protein` variant's
`hgvs_p` is the agent's, and is checked at the first `validate`.
`LIT022` checks an entry's shape, `LIT023` its carriers, `LIT024` its phenotypes and
`LIT025` that `hgvs_p` agrees with the paper's own `p.` name. The model's own pairings
(a `pending_review` without `review_settled_by`, an `unclear` value without its entry
and the reverse) are `SCHEMA001`.

- One entry per distinct variant per genotype context. The context is the partners it is
  carried with and the `phase`: one variant in two contexts is two entries, and one
  variant in one context is one entry (`LIT022` refuses the same name twice in one
  context, and twice in one `with` group). `kind` is among the finding's
  `variant_classes` (`LIT022`).
- `as_reported` is verbatim: the text at the locator. Where the paper prints `c.` and `p.`
  in separate columns, take both, in the order printed, joined by one space. `reported`
  holds each name the paper gives (`c`, `p`, `rsid`, and `transcript` only when the paper
  states one) as written. Never choose between the paper's own names, even when they
  disagree: the checks above report it.
- Placement, by the first rule that fits:
  - `structural` for a CNV or chromosomal change, even written in `c.` notation, with
    `kind` `cnv` or `chromosomal`;
  - `mane` for a `c.` name within the exons ±50 bp, with or without an rsID, with `kind`
    `snv_indel`. A coding-DNA name printed without its `c.` prefix ("456delC", "998G>T")
    is still a `c.` name: write it `mane` with `reported.c` exactly as printed; `normalise`
    decides, and an agent never marks it `unplaced` itself;
  - `outside_region` for a variant the paper places outside the exons ±50 bp (an intron
    beyond 50, upstream or the promoter, downstream, an enhancer), whatever names it: a
    `c.`, an rsID, a `g.` position or a description, with `kind` `snv_indel`. Position
    beats naming, and a position is never looked up in dbSNP, a genome browser or the
    agent's own knowledge;
  - `protein` for a `p.` name alone, with `hgvs_p` on the MANE protein in three-letter
    code in parentheses;
  - `unplaced` for a variant named by an rsID or a `g.` name alone that the paper does not
    place (inside the exons, or of unstated position), or a `protein` variant whose residue validate then refuses (or whose `p.` names
    no residue), which moves to `unplaced`.
- What each placement holds and leaves null (`LIT022`):
  - `mane` holds `hgvs_c`, `grch38` and `normalisation`, all written by `normalise`; an
    agent writes none of them, nor `hgvs_p` of a `mane` variant;
  - `protein` holds `hgvs_p` and leaves `hgvs_c`, `grch38` and `normalisation` null;
  - `structural` and `outside_region` leave `hgvs_c`, `hgvs_p` and `grch38` null;
  - `unplaced` holds `review_status` and leaves `hgvs_c` and `grch38` null; it may keep an
    `hgvs_p`, which is still checked, but a refused one is removed and the note states
    only the mismatch.
- `carriers` per variant, in the finding's `count_unit`. A variant needs
  them on every design but `common_variant_association`, and none the paper gives no
  carrier for is named; a variant with 0 is refused (`LIT022`). The variants' `de_novo` and `inherited` carriers may not exceed
  the finding's, nor may their total, a `with` group counted once (`LIT023`). A
  `common_variant_association` names each variant it tested and has no `carriers`
  (`LIT022`).
- `zygosity` as the paper states it. A variant carried with another variant of this gene
  lists their ids in `with`. Links are mutual, each variant is `heterozygous`, every
  entry of a group carries the same total, and `phase` is set exactly when `with` is
  (`trans`, `cis`, `not_stated` or `unclear`); a group of more than two is never
  `trans` (`LIT022`). A group whose phase is `unclear` carries one identical
  `unclear.phase` entry on every member.
- `phenotypes` only where the paper ties the lesion to these carriers, each among the
  finding's own (`LIT024`).
- Three ways to leave a value plain. `not_stated` only when the paper is silent about it.
  `unclear`, with its `reason` and `settled_by`, when the paper says something one value
  cannot hold (heterozygous in 2 carriers, homozygous in 1) or two incompatible things:
  only `zygosity` and `phase` take it, and the model refuses the value without its entry
  and the entry without the value. A review mark (`pending_review`) is a question about
  the record as a whole: ask whether it would go away if this one field had a value the
  atlas can hold, and if so use `unclear`. Every `pending_review` names
  `review_settled_by`, and `not_accessible` takes none.
- A settler names what would settle the question: `full_text` (only the abstract was
  read), `authors` (the paper contradicts itself or the reference sequence), `expert` (a
  domain judgement the paper cannot supply) or `atlas_model` (the paper is clear and the
  atlas's record cannot hold it). An `unplaced` variant named by an rsID or a `g.` name alone that the paper does not
  place is `atlas_model`; a `protein` variant whose residue is
  refused is `authors`, or `atlas_model` when its `p.` names no residue;
  a `mane` variant `normalise` cannot place becomes `unplaced` with `expert`.

One entry in each placement state, as an agent writes it. A `mane` entry (here and in
the `with` example) is before `normalise`, so `validate` reports "a mane variant needs
hgvs_c", "a mane variant needs grch38" and "a mane variant needs normalisation"
(`LIT022`), and, where its `reported.p` names a residue, "the paper names p.R279X but the
variant has no hgvs_p" (`LIT025`). Every other entry validates as written, on TBX5's
frame. Values are illustrative.

```yaml
variants:
  - id: 1
    as_reported: c.835C>T (p.R279X)
    locator: {kind: in_text, element: Table 1}
    kind: snv_indel
    reported: {c: c.835C>T, p: p.R279X, rsid: rs10000001}
    placement: mane
    zygosity: heterozygous
    carriers: {de_novo: 1, inherited: 0, origin_unknown: 0}
    phenotypes: [HP:0001631]
  - id: 2
    as_reported: p.R279X
    locator: {kind: in_text, element: Table 2}
    kind: snv_indel
    reported: {p: p.R279X}
    placement: protein
    hgvs_p: NP_852259.1:p.(Arg279Ter)
    zygosity: heterozygous
    carriers: {de_novo: 0, inherited: 0, origin_unknown: 1}
  - id: 3
    as_reported: 12q24.21 deletion including the gene
    locator: {kind: in_text, element: Case 3}
    kind: cnv
    placement: structural
    zygosity: heterozygous
    carriers: {de_novo: 0, inherited: 0, origin_unknown: 1}
  - id: 4
    as_reported: c.1234+51A>G
    locator: {kind: in_text, element: Table 2}
    kind: snv_indel
    reported: {c: c.1234+51A>G}
    placement: outside_region
    carriers: {de_novo: 0, inherited: 1, origin_unknown: 0}
  - id: 5
    as_reported: rs10000002
    locator: {kind: in_text, element: Table 3}
    kind: snv_indel
    reported: {rsid: rs10000002}
    placement: unplaced
    carriers: {de_novo: 0, inherited: 0, origin_unknown: 1}
    review_status: pending_review
    review_note: The paper names the variant by rsID alone and does not give its position.
    review_settled_by: atlas_model
```

Two variants carried together, one on each allele (the same total on both):

```yaml
variants:
  - id: 1
    as_reported: c.100A>G
    locator: {kind: in_text, element: Table 1}
    kind: snv_indel
    reported: {c: c.100A>G}
    placement: mane
    zygosity: heterozygous
    with: [2]
    phase: trans
    carriers: {de_novo: 0, inherited: 1, origin_unknown: 0}
  - id: 2
    as_reported: c.200C>T
    locator: {kind: in_text, element: Table 1}
    kind: snv_indel
    reported: {c: c.200C>T}
    placement: mane
    zygosity: heterozygous
    with: [1]
    phase: trans
    carriers: {de_novo: 0, inherited: 1, origin_unknown: 0}
```

A variant whose carriers differ in zygosity is one entry, `unclear`, not two entries for
one name:

```yaml
variants:
  - id: 1
    as_reported: p.N174S
    locator: {kind: in_text, element: Table 2}
    kind: snv_indel
    reported: {p: p.N174S}
    placement: protein
    hgvs_p: NP_852259.1:p.(Asn174Ser)
    zygosity: unclear
    unclear:
      zygosity:
        reason: Table 2 lists two carriers as heterozygous and one as homozygous; the atlas records one zygosity per variant.
        settled_by: atlas_model
    carriers: {de_novo: 0, inherited: 0, origin_unknown: 3}
```

## Attribution and approval

- `extracted_by` and `verified_by` are `agent:<skill>/<model>`. The verifier is a
  different model from the extractor, in a fresh context.
- `approved_by` is a human's name. A value starting `agent:` is refused by the model
  (`SCHEMA001`) and by `chd-atlas literature review` and `literature apply-review`. `approved_by`, `approved_on` and
  `approved_digest` are set together (`SCHEMA001`).
- `approved_digest` is a hash of the finding as it was approved; `literature review
  --approve` and `literature apply-review` both write it. Editing any field it covers (the id, counts, summary, locators,
  anything but the approval, `verified_by` and the review mark) leaves the approval
  stale: the validator reports `LIT015` and the build refuses until the three approval
  fields are deleted and a person approves again. Never write or copy a digest by hand.
- An agent never runs `literature review --approve` and never writes `approved_by` by
  hand. `literature apply-review` writes it, so an agent runs it only on a decisions file
  the review page downloaded, or one curator converted from the page's store, after the
  reviewer pressed "I'm done", with `--by` the reviewer's name; it never writes or edits
  that file by hand. That is a convention: `apply-review` refuses a `--by` that is not the
  file's reviewer or that names an agent, an unfinished session and a second opinion, but
  cannot tell who produced the file, and `review --approve` cannot tell a person's name
  from an agent's invented one. `apply-review` passes `allow_unverified=True`, so it
  approves findings no verifier has compared; the reviewer's decision is the check.
- Approval is not dated before the search that found the paper (`LIT008`).
- A review mark (`pending_review` or `not_accessible`, with its `review_note`) is
  published with the finding. It does not withdraw an approval.

## Synthesis

An optional `synthesis` block in a gene's literature file: the "In short" paragraph of
its page.

- Reader: a clinical geneticist. Cover the lesion spectrum, inheritance, syndromic
  features, genotype–phenotype observations, and how much rests on abstracts alone. It
  synthesises this gene's approved findings and lesion assertions only. The card states
  the validity grade, the burden results and the curated mechanism above it; the
  synthesis does not restate them.
- Sources: atlas records only. Every sentence rests on a record published on the same
  page: an approved finding of this file or a lesion assertion of this gene.
- Every sentence cites at least one such record in `cites`. A cite that is not an
  approved finding of this file or a lesion assertion of this gene is `SYN001`.
- No numbers (`SYN002`); counts live in the card's literature line. It refuses a
  standalone digit run; a number word (one to twenty, the tens from thirty to ninety,
  dozen, hundred, thousand, million, half, quarter, majority, most, twice, thrice,
  percent, per cent, and dozens, hundreds, thousands, millions, halves, quarters,
  majorities), alone or glued to "fold" ("twofold"); a digit run read as an ordinal
  ("3rd"), a decade ("1990s", "mid-2000s") or a multiple ("3x"); a numeric symbol
  ("½", "³"); and `>` or `<` before a digit. Gene names, ontology ids and variants such
  as NKX2-5, 22q11.2, HP:0001631 and p.Gly80Arg pass. Terms that trip it, and what to
  write instead: "half-sibling" ("a sibling sharing only the mother", or "only the
  father"), "trisomy 21" ("Down syndrome"), "45,X" ("Turner syndrome"), "46,XX" ("an
  XX karyotype"), a numbered subtype ("Noonan syndrome 1": name the syndrome without
  its number), "one" as a pronoun ("one of the reports": rephrase), a frameshift with
  its terminal count ("p.Gly35Argfs*75": name the variant class), a transcript id
  ("NM_181486.4": omit it) and an en dash in a gene name ("NKX2–5": use a hyphen).
- One to five sentences, and at most 120 words in all (`SYN003` is the 120 words; the
  sentence count is the model's rule).
- `drafted_by` is `agent:<skill>/<model>`; the model enforces it (a name without the `agent:`
  prefix is a schema error in `validate`). An agent only drafts. `approved_by` is a
  person's name, and `approved_by`, `approved_on` and `approved_digest` are set together.
- The approval covers the text and its cites, `drafted_by`, every approved finding of
  the file (its content, its review mark and its approval: marking it `pending_review`
  or `not_accessible`, or re-approving it, voids the synthesis's approval), which
  findings are approved, and every lesion assertion of the gene. Changing any of them
  makes the validator report `SYN004`, and the synthesis is not published until a
  person approves it again.
- Withdrawing the approval of a finding a sentence cites, or deleting a cited
  assertion, is `SYN001`, an error: the build refuses until the finding is approved
  again or the sentence stops citing it.
- A person lists the gene with `chd-atlas literature review <SYMBOL>`, which shows the
  text, what each sentence cites and the command that approves it, with the digest of
  what was listed: `chd-atlas literature review <SYMBOL> --approve-synthesis --by <name>
  --digest <digest>`. It refuses a synthesis that changed since the listing, and any
  approved finding of the file that fails `LIT015`.

```yaml
synthesis:
  sentences:
    - text: Reported lesions include atrial and ventricular septal defects.
      cites: [CHDA:LIT:0000001, CHDA:AST:0000001]
  drafted_by: agent:gene-synthesis/model
```
