# Published data API

The atlas is served as static JSON. There is no server, no query language and no
pagination: a consumer fetches whole files and filters them client-side.

Build it with `uv run chd-atlas build --root . --out dist`. The command refuses
to write anything if `chd-atlas validate` would report an error, so a published
site is always one that passed its own gate. Every example below is copied from
a real build of the committed corpus rather than written by hand, except where
the example itself says otherwise: a long array or string is abridged with `…`,
and a shape the committed corpus holds no instance of — a contrast omics row, a
populated variant index, a contested gene curated for a lesion group — is
marked *illustrative* where it appears.

Two builds of one commit are byte-identical. Nothing in the output carries a
timestamp — a consumer wanting a publication date should read the commit named
in `manifest.json`.

## Reading this API

**Never construct a path.** Every payload that refers to another file carries
the path to fetch. `genes/index.json` gives each gene its `bundle`, and the
search index gives each record its `path`. The slug rule that turns `HGNC:11604`
into `HGNC_11604` is an implementation detail and is not part of this contract.

**`.json.gz` is not transparently decompressed.** GitHub Pages serves a `.json.gz`
as `application/gzip` and sets no `Content-Encoding`, so the browser hands you
raw compressed bytes. Decompress them yourself:

```js
const response = await fetch("search/index.json.gz");
const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
const { records } = await new Response(stream).json();
```

Plain `.json` files are gzipped in transit by Pages and need no such handling.

**A contested gene must never be displayed as settled.** See the note at the end
of this document; it is the one consumer obligation this API imposes, and
`atlas_curation` extends it to the genes this atlas has not yet curated — 117 of
the 119 published today.

---

## `manifest.json`

What the build produced, and a checksum for every file in it.

```json
{
  "counts": {
    "assertions": 2, "burden_combined": 108, "burden_rows": 1210,
    "burden_surviving_genes": 32, "cell_type_datasets": 0, "cell_type_genes": 0,
    "cohort_families": 2, "cohorts": 20, "datasets": 1, "featured": 1, "functional": 0,
    "genes": 119, "phenotypes": 8, "profile_datasets": 1, "profile_genes": 119,
    "publications": 10, "scope_terms": 105
  },
  "files": {
    "genes/index.json": "sha256:<64 hex>",
    "publications.json": "sha256:<64 hex>"
  },
  "schema_version": "2.20",
  "source_commit": "<40-hex commit sha, or null outside a git checkout>",
  "status": "in-development"
}
```

The structure above is copied from a real build; the checksum and commit values
are shown as placeholders on purpose. Both are derived from content that changes
with every commit, so pinning real ones here would guarantee this document is
wrong by the next one.

- `files` maps a relative URL to the sha256 of **the bytes actually served** at
  it. For a `.json.gz` that is the digest of the compressed file, which is what
  a consumer can verify against what it downloaded.
- `manifest.json` is absent from its own `files`. It cannot be present: the
  value would have to be computed from bytes containing it.
- `source_commit` is `null` when the build was made outside a git checkout — an
  unpacked tarball still produces a complete site, just one that cannot state
  its provenance.
- `counts` is a census of **what this build published**, which is not the same as
  a count of files and — since `2.7` — no longer only a count of curated records.
  Eight of its keys count a curated file, whole: `assertions`, `datasets`
  (omics), `featured`, `functional`, `phenotypes` and `publications` count
  curated records, and `cohorts` and `scope_terms` count the two registries
  [`cohorts.json`](#cohortsjson) and [`scope.json`](#scopejson) publish entire
  — every sample collection and every in-scope disease term, not the subset the
  published genes happen to cite. A gene the atlas has not curated contributes
  to none of them and still gets a bundle — 117 of the 119 genes published today
  are in exactly that position, published on an external authority's
  classification rather than on curation done here. Read `atlas_curation` on the
  browse row to tell the two apart.

  The other nine keys are computed by the build, and the first three below were
  added because the object read as a census of the atlas while describing only
  its curation: through `2.6` it published `assertions: 1, datasets: 0,
  functional: 0` for a site serving 23 genes and 290 burden statistics — the
  corpus as it stood that day — and named neither.

  - `genes` is the published population — the number of gene bundles, and the
    length of `genes/index.json`'s `genes` array.
  - `burden_rows` is the number of rare-variant statistics reaching a bundle, which is
    exactly the sum of every browse row's `burden_row_count`. It counts what a
    consumer can fetch, **not** what `mirrors/burden.tsv` holds: that file is
    deliberately wider than the publication gate, and today 63 of its 181 genes
    publish no page.
  - `cohort_families` is the number of independent cohort families, the length
    of any bundle's `independent_datasets.families` array — the same for every
    gene by construction. It is deliberately not named `independent_datasets`,
    because `datasets` above already means omics datasets and two adjacent keys
    reading `datasets: 0` and `independent_datasets: 3` invite a misreading.
  - `burden_surviving_genes` (since `2.19`) is the number of published genes
    with at least one rare-variant row that survives its own study's
    correction — the genes whose gene-page card lists a rare-variant result,
    by the card's own rule (`methods.html#rare-variants`). A survival is only
    as wide as its study's correction, which for PMID:40127276 is within its
    248-gene panel ([`publications.json`](#publicationsjson)'s
    `correction_genes`).
  - `burden_combined` is the number of a study's own combined results reaching a
    published bundle — the total of every bundle's
    [`burden_combined`](#the-bundles-burden_combined-array-one-studys-own-combination-of-its-arms)
    array, not the mirror's row count.
  - `profile_genes` is the number of published genes whose
    [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
    is non-empty.
  - `profile_datasets` is the number of profile datasets those bundles name —
    not `datasets`, which counts every curated omics dataset, published gene or
    not.
  - `cell_type_genes` (since `2.20`) is the number of published genes whose
    [`cell_type_expression`](#the-bundles-cell_type_expression-object-a-cell-type-profile-never-a-contrast)
    is non-empty. It is 0 until a cell-type mirror is committed.
  - `cell_type_datasets` (since `2.20`) is the number of cell-type datasets
    those bundles name, 0 likewise.

  Each of the nine counts only what the published genes' bundles carry, except
  `cohort_families`, which groups the whole burden mirror's cohorts: restricting
  it to the published genes could split one family into what then reads as two
  independent ones. `index.html` states these six keys too: `genes`,
  `burden_rows`, `cohort_families`, `burden_surviving_genes`, `profile_genes`
  and `cell_type_genes`. They are derived once in the build and handed to
  both, so the page and this object cannot disagree.
- `schema_version` is `major.minor`. **Minor** rises when a field is added and
  nothing is removed or repurposed, so a consumer written against an earlier
  version keeps working; **major** rises when a field changes shape or leaves.
  `1.1` added `genes` to omics shard rows and `conflicting_lesion_groups` to
  gene index rows. `2.0` removed `classification` and `source_tier` from the
  curated assertion — the atlas no longer authors a gene-disease validity call
  of its own — and added the gene bundle's `validity` object in their place.
  The removal is what makes it major: a 1.x reader looking for a
  classification on the assertion now finds none. `2.1` added `status`,
  purely additively. `2.2` added `atlas_curation` to every gene index row and
  every gene bundle, also additively — in the release that widened the
  published gene set from the genes this atlas has curated to the genes a
  ClinGen expert panel calls definitive. More rows of an unchanged shape is not
  a schema change; the new field is, and it is what tells the two kinds of row
  apart. `2.3` added `burden` to every gene bundle and `burden_row_count` to
  every gene index row, again additively; both keys are always present, so a
  2.2 reader keeps working and a 2.3 reader never guards for a missing one.
  `2.4` added `pvalue_adjusted` and `pvalue_adjustment` to every burden object,
  and `2.5` added `count_unit`. Both additive, both always present. Both carry a
  *display* obligation heavier than the parsing one, described with the
  [`burden` array](#the-bundles-burden-array-per-study-never-pooled) below: a raw
  and a corrected p can point opposite ways, and two rows counting different
  units are not comparable at all. `2.6` added `independent_datasets`, also
  additive and always present, and its display obligation is the heaviest of the
  three: it is a **count of datasets, not a validity call**, and rendering it as
  a verdict beside a mirrored ClinGen `definitive` tells a clinician something
  the data do not say. See
  [`independent_datasets`](#independent_datasets-a-count-of-datasets-never-a-verdict).
  `2.7` added `genes`, `burden_rows` and `cohort_families` to `counts`, described
  above. Additive: every `2.6` key is present and unchanged, so the only thing a
  `2.6` reader misses is the half of the census that was never there.
  `2.8` added [`admitted_by` and `asserted_by`](#admitted_by-and-asserted_by-why-this-gene-is-here)
  to every gene bundle, and **widened the population from 23 genes to 92**. Both
  keys are additive and always present, so a 2.7 parser keeps working — but the
  population change carries the heaviest display obligation on this list, and it
  is not additive. Until 2.8 every published gene was `definitive`, so a
  consumer could render `headline_confidence` verbatim and be right every time.
  That is now 23 `definitive`, 1 `strong`, 9 `moderate`, 43 `limited` and 16
  `null`. A consumer that renders the chip without the grade's meaning turns
  `limited` — a panel saying the case is **not yet made** — into a weak yes, and
  one that treats presence in `genes/index.json` as "a panel called this a CHD
  gene" is now wrong for 16 genes. This is a population change inside an
  unchanged shape, which is why it is MINOR by the rule above; the obligation it
  creates is real regardless of the version letter.
  `2.10` added `has_no_association_report` and `no_association_reported_by` to
  every gene bundle and browse row, described above. Additive, and a third axis
  rather than a widening of `has_conflicting_evidence`.
  `2.9` added [`cohorts.json`](#cohortsjson) and a `cohorts` count. Additive: no
  existing key changes meaning and no payload loses one. It closes a gap open
  since `2.3` — every burden row has named its sample collections by bare id
  since `burden` was published, and nothing resolved those ids to anything, so a
  consumer had the numbers and none of the caveats that qualify them.
  `2.11` added `profile_datasets` and `profile_genes` to `counts`, restricted
  to `published` for the same reason `genes` and `burden_rows` already are,
  plus a "Developmental expression" pair of cards on `index.html`. Additive:
  both keys are always present, and both were `0` when the version was minted,
  because no `profiles` mirror had been committed yet. `E-MTAB-6814` has
  landed since, and the example above carries the figures a build produces
  today. See
  [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  for the whole layer these two counts summarise.
  `2.12` added `order` to every entry of a dataset record's
  [`stages`](#datasetsjson) array, and **changed the order of the
  `expression_profile.datasets[].stages[]` array in every gene bundle from
  alphabetical to chronological**. The added key is additive and always
  present, so a 2.11 parser keeps working; the reordering is a correction to
  data that was already published, and it is the release. Through `2.11` a
  gene's developmental series came back in dictionary order — `4 week post
  conception` eighth, *after* `19 week post conception`, and `elderly` second
  of the eight post-natal stages. It was deterministic and reproducible and it
  was not a chronology. **If you plotted `stages` in array order you drew an
  axis in the wrong sequence, and you will now draw a different picture; that
  is the fix arriving, not a regression.** If you worked around it by sorting
  the stage tokens yourself, remove the workaround — sorting the tokens is
  exactly what the bug did. Like `2.2`'s and `2.8`'s population changes this is
  MINOR by the rule above, because no field is added, removed or reshaped and
  every entry still carries every key `2.11` published; the display obligation
  is real regardless of the letter.

  `2.13` added a gene's **withdrawn** HGNC symbols to its `terms` in
  [`search/index.json.gz`](#searchindexjsongz). Measured against `2.12`: 35 of
  the 92 indexed genes gain a term, 46 terms in total, and none of them already
  named anything else — so a query answered by one gene before is still answered
  by one gene, and `HOS`, `VEGF` and `HTX1` now find TBX5, VEGFA and ZIC3
  instead of nothing. MINOR by the rule above: no field is added, removed or
  reshaped, and `terms` was already a variable-length array no consumer may
  depend on the length of. **One term now names two records** — `MADH7`,
  withdrawn by both SMAD6 and SMAD7 — which `terms` already permitted through a
  shared alias and which a client must present rather than collapse.
  `2.14` added [`burden_combined`](#the-bundles-burden_combined-array-one-studys-own-combination-of-its-arms)
  to every gene bundle and to `counts`. It is additive and always present, `[]`
  where a gene has none. Its display obligation: `local_fdr` is never a
  corrected p.
  `2.15` added [`scope.json`](#scopejson) and `scope_terms` to `counts`, additively,
  and **widened the population from 92 genes to 119**: syndromes where structural
  CHD is frequent but not part of the diagnosis entered scope on Orphanet's
  phenotype frequencies (Noonan syndrome and Noonan syndrome with multiple
  lentigines, Costello, cardiofaciocutaneous, Adams-Oliver, Myhre and
  Rubinstein-Taybi syndromes, and their subtypes), and NOTCH1 on its own ClinGen
  term. Its display obligation is that a Noonan grade is a grade for the
  syndrome: say why the syndrome counts,
  which is what `scope.json`'s `orphanet` object is for. It also brought the first
  published genes that set `has_conflicting_evidence` — see
  [contested genes](#contested-genes-the-one-consumer-obligation).
  `2.16` added `submitted_on` and `sgc_id` to every record in a bundle's
  [`validity.records`](#the-bundles-validity-object-mirrored-attributed-never-authored-here),
  `null` on a ClinGen record. Additive. Before it, a GenCC classification
  carried no date at all.
  `2.17` added [`identity`](#genesslugjson) to every gene bundle — name, locus,
  Ensembl and NCBI Gene ids, UniProt accession and MANE Select transcript.
  Additive.
  `2.18` added `summary_caveat` to every record of
  [`publications.json`](#publicationsjson) and so to `featured.json`'s
  embedded publications, `null` where none is curated, and added
  [`methods.html`](#methodshtml), the site's reading notes. Additive.
  `2.19` widened the structural heart phenotypes [`scope.json`](#scopejson)'s
  Orphanet warrants are read against by six lesions; ten warrants gained
  pulmonic stenosis, and no term or gene entered or left scope. It added
  `correction_genes` and `genome_wide_threshold` to every record of
  [`publications.json`](#publicationsjson) (and `featured.json`'s embedded
  ones), `null` where not curated, and set PMID:40127276's `tests_reported` to
  `null`: its 248 counted panel genes, not comparisons. It declared DDD and
  DECIPHER overlapping (`overlaps_with` on every [`cohorts.json`](#cohortsjson)
  record, `[]` where none), which merged PMID:34324492 and PMID:42230622 into
  one cohort family: `counts.cohort_families` went from 3 to 2, and every
  `independent_datasets` array from three families to two, each family gaining
  `members`, every study's own state. It added `file_version`
  and `file_date` to every [`sources.json`](#sourcesjson) entry, and `syndrome` to
  every [`scope.json`](#scopejson) term, `false` unless curated, and
  `burden_surviving_genes` to `counts`. It registered seven more case
  collections, the CHD case series PMID:34324492's S2 Table identifies by PMID
  with 50 or more cases, and that study's rows name them in `case_cohorts`:
  `counts.cohorts` went from 13 to 20, and `counts.cohort_families` did not
  move, since no other study names any of them. And it added
  [`n_chd_case_carriers`](#the-bundles-burden-array-per-study-never-pooled) to
  every object of a bundle's `burden` array, `null` but where a study counted
  the CHD cases among its case carriers apart. It added `review_status` and
  `review_note` to every [`cohorts.json`](#cohortsjson) record and every
  [`scope.json`](#scopejson) term: `"pending_review"`, a question about the
  record that this atlas leaves open for an expert, or `"not_accessible"`, a
  source the record rests on could not be read, with a note stating what is
  unresolved; `null` and `null` where none is recorded. Additive.
  `2.20` added [`cell_type_expression`](#the-bundles-cell_type_expression-object-a-cell-type-profile-never-a-contrast)
  to every gene bundle, `{"datasets": []}` where none is covered, and
  `cell_type_genes` and `cell_type_datasets` to `counts`. It added
  `cell_types` and `regions` to every [`datasets.json`](#datasetsjson) record,
  `[]` on every design but the new `"cell_type_profile"`, and the archive
  `"ucsc_cell_browser"`, whose datasets are named by a lowercase slug (`hoc`)
  rather than an accession, so a dataset id may now be one. Additive.
- `status` is the atlas's own readiness, so a program can read it without
  scraping `index.html`'s prose. Today it is always `"in-development"` — one
  curated gene-disease assertion alongside mirrored ClinGen/GenCC validity for
  many more genes than that. This is a research resource, not a clinical
  decision-support tool; see `index.html` at the site root for the statement
  in full, and the note on [contested genes](#contested-genes-the-one-consumer-obligation)
  below for what the mirrored validity fields do and do not assert.

## The site root: `index.html`

The page a person opens directly rather than fetches as JSON, and the entry
point to the other 121. It carries the same development-status and
research-use statement `status` above is the machine-readable half of, as every
page does, and is a doorway rather than the census (since 2026-10-27): it says
what the atlas is, states the figures a first-time reader needs, and sends
them to `genes/index.html`. It states seven keys of [`counts`](#manifestjson):
`genes`, `burden_rows` and `cohort_families` in the hero,
`burden_surviving_genes`, `profile_genes` and `cell_type_genes` as cards, and `assertions` in the
sentence naming the curated genes. Each equals its key in `manifest.json`,
which states the rest. Two cards are not manifest keys: the genes with a
reported rare-variant result (any row at all; the label read "Genes with
rare-variant evidence" until `2.19`) and the genes this atlas has itself
curated. Nor is the hero's split of the published genes by which upstream
assessment admitted them. The independent-datasets figure
carries "(as far as the cohort registry records)": it counts cohort families,
which join two studies only through a collection both name or collections
[`cohorts.json`](#cohortsjson) records as overlapping, so two studies that
share people the registry does not record count as two. A gene page's card
qualifies its own dataset count the same way. What one rare-variant statistic
and one independent dataset are is [`methods.html`](#methodshtml)'s
`#rare-variants`, which the hero links.

Below the cards it says, for any published study whose correction ran over a
gene panel, that the correction was not genome-wide, and how many of the
rare-variant studies whose results the site publishes are by an author of this
atlas — both counts computed, from the studies with a row on a published gene
and their `own_lab` in [`publications.json`](#publicationsjson) — and that a
gene page says so wherever it shows such a study's results. Its *Sources* names every source
[`sources.json`](#sourcesjson) records as `permitted_with_attribution`, with
the licence where it is CC BY, and gives Orphadata's requested citation. It
links to `genes/index.html`, `genes/index.json`, `manifest.json`,
`sources.json`, [`methods.html`](#methodshtml) (from the navigation bar every
page shares, and from the hero) and the repository. Self-contained: no external
request, no build timestamp, byte-identical between two builds of one commit
like everything else here.

**Every page is checksummed in `manifest.json` exactly like a payload.** The
build behind this document publishes 122 of them — this one, `genes/index.html`,
[`methods.html`](#methodshtml) and one per published gene — and each has an entry
in `files` giving the sha256 of the bytes served at it. A page is published
output, so it is verifiable output.

## `methods.html`

"How to read this atlas": the explanation every gene page assumes, written once
at the site root rather than repeated on each page (design decision D50), and
linked from the navigation bar on every page and from the foot of every gene
page's card, as *How to read this page*. It has six sections, in this order,
and each carries an `id`, so a page can link the paragraph it shortens.

- `#scope` states in full the rule that admits a gene to this atlas and the rule
  that puts a disease in scope, and carries the grade key: what each rung means,
  glossed from ClinGen's Gene-Disease Validity SOP, and where a disputed or
  refuted record is listed. The gene index and every gene page carry a one-line
  short form of the rule that links here.
- `#validity` says how a gene page lists the mirrored validity records: in two
  groups, which of ClinGen's own GenCC copies are listed, and that a record for
  a disease outside this atlas's remit is not shown. Every gene page's validity
  table's intro carries a one-sentence short form of that last point and links
  here. It also says what the card's two flags mean, `conflicting evidence` and
  `sources disagree`: each is a flag on the gene, not on a disease, and the
  records behind it may be for different diseases.
- `#curated-evidence` says what a curated gene's page adds: a lesion
  assertion is a curator's record that the gene belongs to one or more
  lesion groups, with its inheritance, mechanism and syndromic status, and each
  item of evidence behind it cites one publication with its class, the
  curator's grading of its strength and the curator's summary; an assertion
  changes nothing mirrored — a gene is published and its grade headlined on
  the upstream validity records alone; and the *Publications* section lists
  every publication that evidence cites. A curated gene's evidence section
  opens with a line linking here, *How to read curated evidence*.
- `#rare-variants` says how to read a gene page's rare-variant results: what
  one rare-variant statistic and one independent dataset are (the front page's
  figures, whose hero links here), that
  the counts and statistics in its study tables are the studies' own, and the
  one rule the atlas applies to them — a non-synonymous row whose result points
  to enrichment in cases with p below 0.05 supports the gene, and survives its
  study's correction when the corrected p its study published is below 0.05
  too, which decides the card's list, the matrix's states and a forest
  marker's fill for every row but a synonymous one (a synonymous row, never
  support, is filled where its own corrected p is below 0.05, and drawn
  muted); that no
  statistic is pooled across studies, what each count unit means and what an
  absent row does not mean, the synonymous negative control, what a
  damaging composite row is, that a study's `all cases` group is the union of
  its syndromic and non-syndromic groups where it splits its cases (and its
  only case group, not a subgroup, where it does not), what each mark on a
  forest plot means, and what a study block with no effect measure shows
  instead of a plot. A gene page keeps one-line forms of the paragraphs on the
  studies' own numbers and on pooling, and links here, and each forest plot's folded key ends with a link here.
- `#cohorts` lists every sample collection a rare-variant row names, sorted by
  id: its name, id, description and public page where it has one — the same
  records as [`cohorts.json`](#cohortsjson), each followed, where it carries a
  `review_status`, by that open question ("Pending review:" or "Source not
  accessible:", then its `review_note`). Before the list, only where some
  record carries a status, one sentence per status present says what its
  label marks and that the atlas does not settle the question. It says too
  that a row may not name
  every collection its study drew on, so two studies naming no collection in
  common may still share people. A gene page keeps the notes of the cohorts its
  own rows cite, open questions included, folded at the end of its rare-variant
  section under *About these cohorts, and the caveats recorded for them*, and
  that fold ends with a link here.
- `#expression` says how to read a gene page's developmental expression
  section: that every measurement on a developmental timeline here is from
  bulk tissue, so a low or absent figure is not evidence a gene is unimportant
  to heart development; how a median is placed against its organ's percentile
  grid, and that none is placed below its dataset's detection floor, a
  threshold this atlas records per dataset and not necessarily one its source
  declares (a profile dataset's record in [`datasets.json`](#datasetsjson)
  gives the floor and its basis); that percentiles are not comparable across
  organs or developmental stages; what τ measures and when none is published;
  what a stage's developmental phases are, where their boundaries come from
  and how a stage is matched to them; what the charts' axes are; and, for the
  single-cell layer, how a cell-type value is computed (counts per million,
  the median across hearts, a heart counting only with 10 or more cells of a
  type), its detection floor, what a dot's size is, and what a cell-type
  figure does not say. A gene page keeps a one-line form of the first and
  links here, and each expression chart's folded key, and the cell-type
  block's reading fold, end with a link here.

Like the other pages it is for people: its markup is not a contract, and no
payload depends on it.

## `genes/index.json`

The browse payload. Downloaded by every visitor before they pick a gene, so it
carries what ranks or filters a row and the path to fetch the rest — and none of
the evidence itself.

**Which genes are listed.** One row per gene an external authority already
treats as a congenital heart disease gene — 119 genes today, on either of two
warrants:

- a **ClinGen** record classifying it `Limited` or better for a disease in scope
  here (98 genes), or
- **two or more Gene Curation Coalition submitters** independently asserting it
  at `Limited` or better for one, where no ClinGen panel disputes it (21 genes).

A disease is in scope when it is one of the terms this atlas lists, and a term
is listed only if an external authority treats it as congenital heart disease,
or if a curator nominates a syndrome and Orphanet records a structural heart
malformation in at least 30% of its cases (for a subtype, of the syndrome it
belongs to). A disease that would qualify is not in scope until it is listed.

**No disease is in scope on this atlas's own judgement, and neither is any
gene.** `curation/chd_scope.yaml` records, for every term, which authority uses
it, and a validator (`SCP005`) checks that attribution against the mirrors
rather than taking it on trust. Since `2.15` a term can also be in scope on
Orphanet's own phenotype-frequency annotation: a curator nominates a syndrome,
and `SCP006` admits it only if Orphanet records a structural heart phenotype for
it at `Frequent (79-30%)` or stronger. Every term and its warrant is published in
[`scope.json`](#scopejson). Every bundle carries `admitted_by`, naming the
one warrant that admitted that gene. That makes the whole chain a mirrored
decision rather than a curated one: 117 of the 119 genes carry no assertion
authored here, and `atlas_curation` is the field that says which.

**The two warrants are not symmetric, deliberately.** ClinGen admits alone
because it is the only source here with chartered expert panels and published
SOPs; GenCC harmonises submissions rather than adjudicating between them, so it
needs two. ClinGen's own GenCC submissions are excluded from that count — it is
GenCC's largest in-scope submitter (111 rows over 109 genes, measured
2026-08-06), and counting them would let one body vote twice. A ClinGen
`Disputed` or `Refuted` record vetoes the GenCC route entirely: a source trusted
to admit a gene alone must be trusted to refuse one.

**Read `atlas_curation` before presenting a row as curated content.** Every
`headline_confidence` is an upstream expert panel's call and never this atlas's,
and on 117 of the 119 rows there is no assessment by this atlas behind it at all.

**This population widened on 2026-08-06, from 23 genes to 92.** A consumer
written against the earlier data will have seen only `"definitive"` in
`headline_confidence` and only `"expert_curated"` in `validity_state`. Both
assumptions are now false — see the next section.

```json
{
  "genes": [
    {
      "assertion_count": 1,
      "atlas_curation": "curated",
      "bundle": "genes/HGNC_11604.json",
      "confidence_by_lesion_group": { "septal": "definitive" },
      "conflicting_lesion_groups": [],
      "evidence_counts": { "genetic_case": 1 },
      "functional_count": 0,
      "gene": "HGNC:11604",
      "has_conflicting_evidence": false,
      "has_no_association_report": false,
      "has_source_discordance": false,
      "headline_confidence": "definitive",
      "no_association_reported_by": [],
      "lesion_groups": ["septal"],
      "symbol": "TBX5",
      "validity_state": "expert_curated",
      "variant_count": 0,
      "burden_row_count": 14,
      "independent_datasets": {
        "tested": 2,
        "enriched": 2,
        "corrected": 1,
        "families": [
          { "studies": ["PMID:34324492", "PMID:42230622"], "state": "nominal",
            "members": [ { "study": "PMID:34324492", "state": "not_tested" },
                         { "study": "PMID:42230622", "state": "nominal" } ] },
          { "studies": ["PMID:40127276"], "state": "corrected",
            "members": [ { "study": "PMID:40127276", "state": "corrected" } ] }
        ]
      }
    }
  ]
}
```

- The array is ordered by HGNC id as a string, not by symbol and not
  numerically — so `HGNC:10023` precedes `HGNC:1008`, the way `publications.json`
  orders PMIDs. JSON arrays keep their order, so this is part of the contract.
- `bundle` is the path to fetch for the detail page. Do not build it yourself.
- `symbol` falls back to the HGNC id for a gene not yet in `mirrors/genes.tsv`,
  so it is always a non-empty string you can render and search.
- `headline_confidence` and the rest of this row come from the mirrored
  ClinGen/GenCC validity records for the gene, never from a curated assertion —
  this atlas mirrors gene-disease validity, it does not author it. The example
  above is TBX5's real row: its one in-scope ClinGen record — the **Syndromic
  Disorders** Gene Curation Expert Panel, for **Holt-Oram syndrome** — makes it
  `"expert_curated"` with a `"definitive"` headline. It is not a Congenital
  Heart Disease GCEP record, which this paragraph claimed until 2026-08-06 and
  which the mirror has never held for TBX5; `admitted_by.panel` and
  `admitted_by.disease_label` in the bundle now name both, so a reader can check
  this sentence rather than trust it. That a gene's headline can come from a
  panel whose remit is a syndrome, for a disease whose label names no cardiac
  feature, is the ordinary case here and not an exception: 10 of the 23
  originally published genes were graded by some panel other than the CHD GCEP. `headline_confidence` is `null` for a gene no authority has
  assessed. It is never `"no_known_association"` for that case: that
  classification is itself an assessed verdict ("a panel looked and found
  nothing"), and asserting it for a gene nobody has assessed would state a
  conclusion no authority reached.
- `atlas_curation` is `"curated"` when the atlas holds at least one curated
  assertion for the gene and `"not_yet_curated"` when it holds none — 117 of the
  119 rows published today. It is about *this* resource's own work and says
  nothing about the gene's validity: an uncurated row still carries whatever
  `headline_confidence` an expert panel assigned, which is exactly why the gene
  is listed. Read this rather than testing
  `assertion_count == 0`, so a browse filter does not have to reimplement the
  rule. It appears on the browse row and in the bundle, written in one place so
  the two cannot disagree.
- `validity_state` says how well curated the gene is: `"expert_curated"` (an
  in-scope ClinGen record exists), `"submitter_curated"` (only GenCC has
  assessed it) or `"uncurated"` (neither mirror has). `headline_confidence` is
  `null` **iff `validity_state` is not `"expert_curated"`** — the headline is
  the strongest grade among the gene's in-scope ClinGen records, which for a
  gene ClinGen admitted is the admitting record's grade, so a gene no ClinGen
  panel has graded in scope has no headline, however many GenCC submitters
  graded it and however highly. Measured 2026-08-06: the 16 `null` rows are
  exactly the 16 `"submitter_curated"` rows, as sets and not merely as counts.
  This rule read "`null` iff the gene is `"uncurated"` *or* every in-scope
  record maps to no rung" until 2026-08-06, which is false for every one of
  those 16 — GDF1 carries a G2P `Definitive` and a Labcorp `Strong`, and
  headlines `null` because neither is ClinGen's. Six genes in the
  committed mirrors (HGNC:1541, HGNC:24595, HGNC:4317, HGNC:6188, HGNC:9380,
  HGNC:9381) would resolve to `null` while `"submitter_curated"`: each carries
  exactly one in-scope record, an Orphanet `Supportive` submission, which
  `vocab.GENCC_CLASSIFICATIONS` maps to `None` because the submitter asserted an
  association without grading its evidence, not because nobody looked. Do not
  infer "no authority has assessed this gene" from `headline_confidence: null`
  alone — check `validity_state` for that.
- **`headline_confidence: null` and `"submitter_curated"` are now in the
  payload, and a consumer written before 2026-08-06 will not have seen them.**
  Measured against a real build (2026-10-25): `headline_confidence` is
  `"definitive"` on 42 rows, `"limited"` on 45, `"moderate"` on 10, `"strong"` on
  1 and `null` on 21; `validity_state` is `"expert_curated"` on 98 and
  `"submitter_curated"` on 21. (23/43/9/1/16 and 76/16 on 2026-08-06.) Every `null` headline is on a
  `"submitter_curated"` row — no panel graded the gene, and `admitted_by` names
  the submitters instead. Until 2026-08-05 all 23 published rows carried
  `"definitive"` and `"expert_curated"` and those were the only values either
  field took, so any rendering tested against that data is now incomplete.
  `"uncurated"` still does not occur, and the six `Supportive`-only genes above
  are still in the mirrors rather than the output.

- **`limited` is not a weak yes.** It is a ClinGen expert panel saying the
  gene-disease case is *not yet made* — few probands, or evidence that is
  suggestive but not compelling. 45 of the 119 rows carry it, so a consumer that
  renders the string with no gloss, or that colours it on a good-to-bad ramp
  without saying what the rungs mean, presents the plurality of this atlas as
  better supported than the grading authority does. The site's own pages carry
  ClinGen's definitions beside every chip for this reason.
- `has_source_discordance` is `true` when one mirrored source contests the gene
  while the *other* supports it. It is narrower than `has_conflicting_evidence`,
  which also fires when a single source is internally split across diseases or
  panels — see [Contested genes](#contested-genes-the-one-consumer-obligation).
- `has_no_association_report` is **a third axis, not a third rung**, added in
  schema `2.10`. It is `true` when some authority reported *no known
  association* in scope while another asserted one, and
  `no_association_reported_by` names who — sorted, and empty exactly when the
  flag is `false`, so the two can never say different things.

  It is deliberately **not** folded into `has_conflicting_evidence`, which means
  exactly `disputed`/`refuted`. "A panel looked and found no reported evidence"
  is not "a panel disputes this", and merging them would give one laboratory's
  null result the weight of a chartered panel's refutation. ClinGen treats the
  assertion the same way — a distinct verdict rather than a rung of the
  definitive-to-limited ladder.

  Before `2.10` that disagreement reached no published byte at all. **GDF1 is
  the live case and currently the only one:** G2P grades it `Definitive`,
  Labcorp `Strong`, and Illumina reports `No Known Disease Relationship`. It
  published `has_conflicting_evidence: false` with nothing beside it, so a
  consumer was told the evidence did not conflict while two bodies disagreed
  about whether an association exists. A gene whose *only* in-scope record is
  `no_known_association` does **not** set this flag — it is not in disagreement
  with anything, and there are 9 such genes in the mirrors, none published.
- `confidence_by_lesion_group` applies the gene's mirrored `headline_confidence`
  to every lesion group its curated assertions declare — ClinGen and GenCC
  classify a gene against a disease, never against a specific lesion, so there
  is no finer-grained signal to divide the groups with. Its keys are exactly the
  gene's `lesion_groups`, so it is empty whenever that array is empty — which is
  every gene this atlas has not curated yet — and empty as well when
  `headline_confidence` is `null`.
- **An empty `confidence_by_lesion_group` does not mean "not assessed."** It
  means this atlas has recorded no lesion group for the gene, which is a fact
  about this resource's curation queue and says nothing about the gene's
  validity. Measured against a real build (2026-10-26), 96 of the 119 published
  rows carry a non-null `headline_confidence` together with
  `"confidence_by_lesion_group": {}` — an expert panel graded each of them, and
  no curator here has written an assertion for them yet. TBX5 and NR2F2 are the
  only rows with a non-empty map today (`{"septal": "definitive"}` for TBX5;
  `{"septal","avsd","lvoto","heterotaxy"}` all `"definitive"` for NR2F2). A consumer that
  renders an empty map as "not assessed" therefore understates every one of
  those expert-panel calls as unassessed while still counting them among the
  genes it presents. Read `validity_state` for whether any authority has
  assessed the gene, and `atlas_curation` for
  whether this atlas has curated it; this field answers neither question.
- `conflicting_lesion_groups` names every group in `confidence_by_lesion_group`
  when the gene is contested, and none when it is not — the same "no per-group
  signal" reasoning applies, so this can never name a proper subset of the
  gene's declared groups. It is the per-group counterpart of
  `has_conflicting_evidence` — see [Contested genes](#contested-genes-the-one-consumer-obligation),
  which is the one obligation this API places on a consumer.
- **The two lesion-group fields appear here and nowhere else.** The gene bundle
  carries `headline_confidence`, `validity_state`, `has_conflicting_evidence` and
  `has_source_discordance`, but neither of them, so a detail page that needs
  group-level confidence must carry it over from the browse row it was opened
  from. It cannot be derived from the bundle: that would mean reimplementing
  the classification ranking and the contested test, neither of which is
  published.
- The three counts describe what the bundle contains, so a browse row never
  promises more than the page delivers.

## `genes/index.html`

The browse page: the same 119 rows `genes/index.json` publishes, rendered as a
table a person can read and filter. Each row carries the HGNC id and the
symbol — both linked to that gene's page — `headline_confidence`, the disease
it is graded for, the gene's `independent_datasets` as a strip (one glyph
per cohort family, whose hover names each of the family's studies with that
study's own state, from `members`),
`validity_state`, `atlas_curation`, `burden_row_count` and the gene's lesion
groups. Above the table sit a text box
matching id or symbol and five menus (lesion group, confidence, validity state,
atlas curation, rare-variant evidence), whose options are the values actually present
in the build rather than every value the vocabulary allows. Above them, one
line says each classification is an external authority's, for the disease it
names, which for many genes here is a syndrome rather than CHD, and links
[`methods.html#scope`](#methodshtml) for the admission rule and the grade key.
21 of the 119 genes have no headline at all — they are admitted on submitter
agreement and carry `null` — and the browse page renders an em dash (—) in
their confidence and graded-for cells, as it does for every empty cell. The
confidence menu offers no value for them: only *any confidence* shows them. A
grade is plain text, not a chip; a contested gene's confidence cell adds a
`conflicting evidence` chip after it — the same genes `has_conflicting_evidence`
marks — while the row's `data-confidence` attribute, which the menu filters on,
stays the bare grade.

**Every row is rendered by the build, and the inline script only hides rows.**
There is no empty `<tbody>` filled in by a fetch, so `curl`, a crawler and a
reader with JavaScript disabled all get the complete table of 119 genes,
unfiltered — which is what the page shows before anyone touches a control in
any case. Nothing is loaded from anywhere: the stylesheet and the script are
inline, and the page makes no external request.

**`genes/index.json` is the machine-readable contract; this page is not.** The
column set, the markup and the `data-` attributes the filter reads may change
without a `schema_version` bump. Write against the payload and leave the page
to people.

## `genes/<slug>.json`

One gene's whole detail page, in one fetch. TBX5's, below, with every
top-level key the bundle carries; the arrays and objects marked `…` are
abridged, and each is documented in full further down.

```json
{
  "gene": "HGNC:11604",
  "symbol": "TBX5",
  "identity": {
    "name": "T-box transcription factor 5",
    "locus": "12q24.21",
    "ensembl_gene": "ENSG00000089225",
    "ncbi_gene": 6910,
    "uniprot": "Q99593",
    "mane_select": "NM_181486.4"
  },
  "headline_confidence": "definitive",
  "validity_state": "expert_curated",
  "atlas_curation": "curated",
  "has_conflicting_evidence": false,
  "has_source_discordance": false,
  "has_no_association_report": false,
  "no_association_reported_by": [],
  "lesion_groups": ["septal"],
  "validity": { "state": "expert_curated", "has_source_discordance": false, "records": [ … ] },
  "admitted_by": { "authority": "clingen", "classification": "definitive",
                   "disease": "MONDO:0007732", "disease_label": "Holt-Oram syndrome",
                   "panel": "Syndromic Disorders Gene Curation Expert Panel",
                   "submitters": [] },
  "asserted_by": { "count": 6, "institutions": ["Ambry Genetics", "ClinGen", "G2P", … ] },
  "publications": ["PMID:8988165"],
  "assertions": [ { "id": "CHDA:AST:0000001", "lesion_groups": ["septal"], "evidence": [ … ], … } ],
  "functional": [],
  "variants": [],
  "omics": { "profiles": { "count": 119, "shards": ["omics/profiles/E-MTAB-6814.json"], "top": [ … ] } },
  "burden": [ { "study": "PMID:40127276", "cohort_stratum": "all", "consequence_class": "damaging",
                "comparator": "mutation_model", … }, … ],
  "burden_combined": [ { "study": "PMID:40127276", "consequence_class": "damaging", … }, … ],
  "independent_datasets": {
    "tested": 2, "enriched": 2, "corrected": 1,
    "families": [
      { "studies": ["PMID:34324492", "PMID:42230622"], "state": "nominal", "members": [ … ] },
      { "studies": ["PMID:40127276"], "state": "corrected", "members": [ … ] }
    ]
  },
  "expression_profile": { "datasets": [ { "dataset": "E-MTAB-6814", … } ] },
  "cell_type_expression": { "datasets": [] }
}
```

- `identity` is the gene's row of the HGNC registry this atlas mirrors
  (`mirrors/genes.tsv`), with one cell the atlas chooses: RBM10's `uniprot`,
  where HGNC lists two accessions (see the `mirrors/genes.tsv` notes below).
  Every key is present on every gene and `null` only where the registry holds no
  value; `ncbi_gene` is an integer. The gene page links gnomAD by
  `ensembl_gene` and ClinVar by `ncbi_gene`, and the build fails rather than
  publish a gene missing either, or carrying an `ensembl_gene` that is not
  `ENSG` and 11 digits or an `ncbi_gene` below 1, so a build that completes
  never publishes either as `null` or malformed. Since `2.19` `chd-atlas
  validate` reports the same rule as `GEN004`, so the gate refuses such a
  registry before the build begins. The bundle carries the ids, not the URLs.
- `assertions` carry their full `evidence` array, including each item's
  `locator`, `strength` and `summary` — the record a curator is judged on.
  Assertion fields: `curated_on`, `curator`, `evidence`, `extracardiac_features`,
  `gene`, `id`, `inheritance`, `last_reviewed`, `lesion_groups`, `mechanism`,
  `notes`, `phenotypes`, `syndromic`. **No `classification` or `source_tier`**:
  the atlas mirrors gene-disease validity rather than curating its own, so an
  assertion says only what the curator is the authority for -- which lesions a
  gene is claimed for, and on what evidence. `headline_confidence` above is
  where the mirrored classification lives.
- `functional` holds **every** functional record about the gene, not only those
  an assertion cites.
- `atlas_curation` reads the same here as on the browse row. On the 117 genes
  published today without curation here it is `"not_yet_curated"`, and
  `assertions`, `publications` and `functional` are then empty arrays: the page
  is the panel's classification plus whatever this atlas has recorded, which
  may be nothing at all. An empty `assertions` array is a curation gap, never a
  fetch that failed.
- `omics`, `variants` and `burden` are always present and may be empty. Read
  them without guarding for a missing key.
- **`variants` are embedded; omics rows are linked.** That asymmetry is a
  curation policy, not a property of the data — this atlas curates variants by
  hand, so the count per gene is bounded by effort. Omics tables are not, so a
  bundle carries per-modality summaries with `shards` to fetch. The omics section
  below gives their shape and how to select a gene's rows out of a fetched
  shard.
- `publications` lists the PMIDs the gene's assertion evidence cites, in lexical
  order. It does not include PMIDs cited only by its functional records.
- `assertions` and `functional` are ordered by id.

### The bundle's `validity` object: mirrored, attributed, never authored here

**This atlas publishes no gene-disease validity classification of its own.**
`headline_confidence`, `validity_state`, `has_conflicting_evidence` and
`has_source_discordance` at the top of the bundle come entirely from ClinGen
and GenCC, mirrored and attributed rather than asserted by a curator — see
[Contested genes](#contested-genes-the-one-consumer-obligation) for how to
read them safely. The bundle's `validity` object is where the mirrored records
behind those fields live:

- `state` repeats `validity_state`. `has_source_discordance` repeats the
  top-level field of the same name. Both are published twice on purpose: the
  flat fields are what the browse row (`genes/index.json`) and the bundle
  publish identically, and this object is the self-contained provenance
  record for a consumer that wants only "who curated this gene, and what did
  each of them say" without cross-referencing the fields beside it.
- `records` is one entry per in-scope mirrored classification, ClinGen's and
  GenCC's alike, in no particular order a consumer should rely on beyond what
  is published. Every record carries the same key set regardless of source:

  ```json
  {
    "source": "clingen",
    "classification_term": "Definitive",
    "classification": "definitive",
    "disease": "MONDO:0007732",
    "disease_label": "Holt-Oram syndrome",
    "moi": "AD",
    "sop": "SOP11",
    "classification_date": "2025-03-25T16:00:00.000Z",
    "gcep": "Syndromic Disorders Gene Curation Expert Panel",
    "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_24e6c85a-33cf-4248-be1f-6431c7c6b1e5-2025-03-25T160000.000Z",
    "submitter": null,
    "submitted_on": null,
    "sgc_id": null
  }
  ```

  - `classification_term` is the authority's own word, verbatim — "Definitive",
    "Disputed Evidence", "Supportive". `classification` maps that term onto
    this atlas's `Classification` scale (the same values `headline_confidence`
    publishes), or `null` where the term is not a rung on it at all — GenCC's
    `Supportive`, a submitter asserting an association without grading its
    strength. A consumer that wants to render exactly what the authority said
    reads `classification_term`; one that wants to filter or rank reads
    `classification` and must handle `null`.
  - `sop`, `classification_date` and `gcep` are populated only on a ClinGen
    record; `submitter`, `submitted_on` and `sgc_id` only on a GenCC record.
    The field a record's own source does not carry is published as `null`
    rather than omitted, so every object in `records` has the same shape and a
    consumer never has to check `source` before it can look a key up.
  - **A GenCC record's date is `submitted_on`, not `classification_date`.**
    `classification_date` is ClinGen's and is `null` on every GenCC record.
    `submitted_on` is the date the submitter attached to its submission
    (GenCC's `submitted_as_date`), and `sgc_id` is GenCC's identifier for the
    submission; both since `2.16`, both verbatim from the mirror. Every GenCC
    record's `submitted_on` is either `null` or begins with an ISO calendar
    date: the mirror's schema enforces the prefix and permits a missing date,
    and all 30,410 mirrored submissions carry one today. It is `null` on every
    ClinGen record. When it is not `null`, take its first ten characters: what
    follows varies by submitter, and occasionally within one
    (`2018-05-31 20:04:32`, `2025-01-17T00:00:00.000000Z`), and the
    space-separated shape states no time zone, so do not parse it as an
    instant. Key a GenCC record on `sgc_id`.
  - **`report_url` does not always lead back to the submission.** Of the 470
    GenCC records across the 119 gene bundles, 189 publish `report_url: null`
    and 76 more point at `https://panelapp-aus.org` — a homepage, not a
    submission — so 265 of the 470 offer a reader no published route back to
    the record behind the classification (measured 2026-10-25; 157 of 266 on
    2026-08-06, before the population widened). Attribute these to their
    `submitter` and present them as classifications this atlas mirrors; do not
    render them as though every one resolves to a citable report. The bundle
    publishes PanelApp Australia's homepage as GenCC gave it. Since `2.19` the
    site's gene pages link such a row to
    `https://panelapp-aus.org/panels/entities/<symbol>` instead, the gene's own
    page on that site (built from the HGNC symbol, percent-encoded) — the one
    report URL this atlas constructs rather than republishes; a client may do
    the same. It leads to the gene's PanelApp Australia panels, not to the one
    submission, which the URL GenCC gave does not identify.
  - `sop` is published because ClinGen's committed mirror spans SOP4 through
    SOP12 with no crosswalk between framework versions published anywhere. A
    classification attributed to ClinGen without its SOP version is an
    unqualified claim — the same rule under a different framework applies
    equally to a rung on the classification ladder itself.
  - `moi` and `disease`/`disease_label` are the authority's own mode of
    inheritance and disease term for that specific record — ClinGen and GenCC
    do not always agree on either, which is part of why two records for one
    gene can differ.
- `state` is one of three values, the same ones `validity_state` publishes:
  `"expert_curated"` (at least one in-scope ClinGen record exists),
  `"submitter_curated"` (only GenCC has assessed the gene) or `"uncurated"`
  (neither mirror has). An uncurated gene carries `records: []` and
  `has_source_discordance: false` — the explicit empty shape, not an absent
  `validity` key, so a consumer can tell "no authority has assessed this
  gene" from "the field is missing" without guessing. **Two of the three are
  reachable in the payload since 2026-08-06**: measured on a real build
  (2026-10-25), 98 bundles say `"expert_curated"` and 21 say `"submitter_curated"`. Before the
  gate widened only the first occurred, because a gene was published on an
  in-scope ClinGen `Definitive` call and that made it `"expert_curated"` by
  construction. `"uncurated"` remains unreachable — a published gene has at
  least one in-scope record by definition — so no bundle carries an empty
  `records` array, and a consumer that treats `records: []` as a failed fetch
  would still break on the first gene that had one.

`sources.json` (below) carries the licence terms this atlas mirrors ClinGen
and GenCC under, the same way it does for HPO.

### `admitted_by` and `asserted_by`: why this gene is here

Added in schema 2.8. The gate is a claim about external authorities, and these
two keys are what let a consumer **check** it rather than trust it.

`admitted_by` is the single warrant that cleared the gate. Exactly one warrant
admits a gene, and ClinGen is checked first, so a gene with both is published on
ClinGen's word:

```json
{ "authority": "clingen", "classification": "definitive",
  "disease": "MONDO:0007732", "disease_label": "Holt-Oram syndrome",
  "panel": "Syndromic Disorders Gene Curation Expert Panel", "submitters": [] }
```

```json
{ "authority": "gencc_agreement", "classification": null, "disease": null,
  "disease_label": null, "panel": null,
  "submitters": ["Ambry Genetics", "Labcorp Genetics (formerly Invitae)"] }
```

`authority` is `"clingen"` on 98 of the 119 genes and `"gencc_agreement"` on 21.
**Every key is present on every gene** — `submitters` is `[]` rather than absent
where ClinGen admitted the gene, and `classification`/`disease`/`disease_label`/
`panel` are `null` rather than absent where GenCC did — because an object whose
shape varies is a trap for a consumer reading a field off one gene and expecting
it on the next.

`admitted_by.classification` is the same value as `headline_confidence` on
every gene published today, but the two are computed separately and are not
one figure by construction. `admitted_by.classification` is the grade of the
record that cleared the gate; `headline_confidence` is the strongest grade
among the gene's in-scope ClinGen records, and for a gene ClinGen admitted
those are the same record. They would part on a gene admitted on submitter
agreement that also carries an in-scope ClinGen `No Known Disease Relationship`
record — `null` here, `"no_known_association"` in the headline — and no such
gene is published. A `"gencc_agreement"` gene otherwise has `null` for both,
however highly its submitters graded it.

`asserted_by` is every distinct institution **asserting** the gene in scope,
deduped by institution, with a count:

```json
{ "count": 5, "institutions": ["Ambry Genetics", "G2P",
  "Labcorp Genetics (formerly Invitae)", "Orphanet", "PanelApp Australia"] }
```

Two things it deliberately is not. It is **not a count of records** — ClinGen
submits to GenCC under its own name, so counting `gcep` and `submitter` values
naively overstates by exactly one per gene (134 against 111 over the 23 genes
published before the widening). And a **dissent is not an assertion**: a record
of `no_known_association`, `disputed` or `refuted` is excluded, because counting
the institutions that disagree as institutions that agree is the opposite of
what the field's name says. GDF1 is the case — it published `"count": 6`
including Illumina's `No Known Disease Relationship` until 2026-08-06, and GDF1
is the gene that motivates requiring two submitters in the first place. GenCC's
`Supportive` **is** counted: it asserts an association without grading its
evidence, which is an assertion, merely ungraded.

**`count` is not a score, and must never be rendered as one.** Eight authorities
asserting a gene is not evidence it is eight times better supported than a gene
with one; it frequently means it sits on more commercial test panels. The atlas
publishes no validity call of its own, and a rank derived from this count would
be exactly that.

Neither key appears in `genes/index.json` — the browse rows carry
`validity_state` instead, which separates the same two populations
(`"expert_curated"` 98, `"submitter_curated"` 21). Fetch the bundle for the
warrant itself.

### The bundle's `burden` array: per study, never pooled

Published rare-variant statistics for the gene, one object per
(study, cohort stratum, variant class, consequence class, frequency threshold).
**1210 rows reach the API**, across 118 of the 119 published genes, from
3 studies (PMID:34324492, PMID:40127276, PMID:42230622), contributing
77, 216, 917 rows respectively. `mirrors/burden.tsv` holds 1,805 rows for 181 genes; the other 595,
covering 63 genes, are for genes the site does not publish and reach no
bundle and no page. They are held
so that widening the publication gate later needs no re-mirroring — which is
exactly what happened on 2026-08-06, when the gate widened from 23 genes to 92
and 625 already-mirrored rows became fetchable with no change to the mirror.
The same reason `mirrors/genes.tsv` registers 185 genes and 119 publish.

**Terms.** A *burden test* compares how often cases and controls carry
qualifying rare variants in a gene (case-control). A *CNV burden* test does the
same for copy-number variants (deletions, today). A *de novo enrichment* test compares de novo
mutations in affected trios with a mutation-rate model's expectation; it is not
a burden test. The bundle keeps the historical key `burden` for every
rare-variant analysis, case series included; the gene page heads them "Rare-variant association" and names each study's
analysis.

**One published gene carries no burden rows at all** and its `burden` array is
`[]`. That is not a fetch failure: the publication gate is an authority's
classification, not the presence of rare-variant evidence, and a gene no rare-variant study
covered is exactly what a widened gate admits.
`burden_row_count` on the browse row is this array's length.

Every count in this section is asserted against a real build by
`tests/test_docs_match_the_build.py`, so it fails rather than rots when a study
lands. It rotted once — this paragraph described two studies and 200 rows for
one commit after the third study shipped.

```json
{
  "study": "PMID:42230622",
  "cohort_stratum": "syndromic",
  "lesion_group": null,
  "variant_class": "snv_indel",
  "consequence_class": "lof",
  "origin": "any",
  "maf_max": 0.001,
  "count_unit": "individuals",
  "n_case_carriers": 5,
  "n_cases": 1471,
  "n_chd_case_carriers": null,
  "comparator": "control_cohort",
  "n_control_carriers": 0,
  "n_controls": 45082,
  "expected_count": null,
  "effect": null,
  "effect_measure": "odds_ratio",
  "effect_bound": "unbounded_above",
  "ci_low": 28.1,
  "ci_high": null,
  "pvalue": 3.13e-08,
  "pvalue_test": "fisher_exact",
  "pvalue_adjusted": null,
  "pvalue_adjustment": null,
  "case_cohorts": ["cnchd", "ddd", "nottingham"],
  "control_cohorts": ["ukbb"],
  "method_note": null,
  "source": "audain2026_sd3"
}
```

Every key is present on every object, `null` where the row's comparator does
not populate it, so a consumer never has to guard for a missing one.

**`count_unit` qualifies all four count keys, and they are not comparable
without it.** The four studies this schema was designed against do not count the
same thing:

| `count_unit` | numerator | denominator |
| --- | --- | --- |
| `individuals` | people carrying at least one qualifying variant | people sequenced |
| `alleles` | qualifying alleles observed | alleles called — roughly twice the people, and varying per gene with coverage |
| `de_novo_mutations` | de novo mutations observed | **trios**, not alleles and not people |

A consumer dividing `n_case_carriers` by `n_cases` across studies without
reading this key is doing arithmetic on three different quantities. Someone
carrying two qualifying variants counts once under `individuals` and twice under
`alleles`.

**`comparator` is the field the whole array turns on.** The published burden
literature answers one question — is this gene hit more often than expected? —
and differs only in what "expected" meant:

| `comparator` | expectation from | populated | `effect_measure` |
| --- | --- | --- | --- |
| `control_cohort` | `n_control_carriers` / `n_controls` | both control fields | `odds_ratio`, `rate_ratio` |
| `mutation_model` | `expected_count` | `expected_count` | `enrichment_ratio` |
| `none` | nothing (a case series) | neither | none, and no `pvalue` |

**Both `control_cohort` (1102 rows) and `mutation_model` (108 rows) appear today;
`none` does not.** PMID:40127276 contributes the `mutation_model` rows — de novo
mutations in 3,887 trios against a mutability-based expectation — and they carry
null control fields, a non-null `expected_count`, and an `enrichment_ratio`
rather than an odds ratio. A consumer that reads only `n_control_carriers` will
see `null` on 108 of the 1,210 rows, several of which are the strongest results the
atlas publishes. The build refuses to publish a row whose statistic contradicts
its comparator.

**`consequence_class` is not a partition: `damaging` is the union of `lof` and
`missense_damaging`.** PMID:40127276 reports its primary analysis over
loss-of-function and damaging missense together, and publishes that composite
alongside its two components — 72 of the 1210 rows. Its `n_case_carriers` is
exactly the sum of the two component rows for the same gene, stratum and
comparator, verified on every one. **Summing `n_case_carriers` across
`consequence_class` therefore double-counts those variants.** Group by
`consequence_class` and pick, or exclude `damaging`; do not aggregate over it.
The composite is published rather than dropped because its p-value is not a
function of its components' and is the statistic that study defines its results
by.

Six obligations, each of which is a wrong claim if you get it wrong:

1. **Never render `effect` without `effect_measure` beside it.** One column
   holds odds ratios and de novo enrichments alike, and an odds ratio of 3.1 and
   an enrichment of 3.1 are different claims. A cell reading `3.1` under a header
   reading "effect" equates them.
2. **`effect: null` with `effect_bound: "unbounded_above"` is a result, not a
   missing one.** Fisher's exact test returns an infinite odds ratio where no
   control carries; `Infinity` is accepted by `JSON.parse` nowhere, so the
   number cannot be published and `ci_low` carries the finding instead — TAB2's
   syndromic loss-of-function row has a 95% interval that starts at 28.1. 29
   published rows today. An infinite point estimate says nothing about how
   strong the result is: 14 of those 29 rows have a `ci_low` below 1, an
   interval that includes no association at all. Render `ci_low` with it, as
   the gene pages do (`OR ∞ (95% CI 28.1–∞)`); a blank, "not tested" or a bare
   "∞" each say something the row does not.
3. **An absent (stratum, consequence) cell is not a null result — but why it is
   absent is the study's rule, not the atlas's.** Read it per study rather than
   as one law. For `PMID:42230622`, zero of its 1,438 mirror rows have no case
   carrier *and* no control carrier, so a 2×2 of all zeros supports no test and
   that study emitted none. The other two do not follow the same rule:
   `PMID:34324492` tests one consequence class by construction (CNV deletions),
   and `PMID:40127276` **observed 14,364 synonymous variants and still published
   no synonymous row**, because its gene-level table reports only damaging
   classes. So an absent cell means "this study did not report one here", and
   only for `PMID:42230622` does it additionally mean "no carrier on either
   side". `PMID:40127276` does publish a row with zeros on both sides: WDR5's
   loss-of-function case-control row.
4. **Do not filter out `consequence_class: "synonymous"`.** It is a study's own
   negative control — synonymous variants should show no enrichment — and where
   one is significant, that gene's comparison is poorly calibrated, which a
   reader can only see if you show it. **It is not available everywhere:** 339
   of the 1210 published rows are synonymous and all of them come from
   `PMID:42230622`. The other two studies publish no synonymous row at all, so
   their results have no negative control on this page and must be read without
   one. (This obligation once called synonymous "the most numerous class, 435
   rows against 345 loss-of-function", measured against a one-study corpus; then
   "92 loss-of-function and 69 synonymous", measured against the 23-gene
   population. Then "306 and 261" over the 92-gene population. All three went
   stale the same way. Today loss-of-function is 401 and synonymous 339, beside
   398 damaging missense and 72 damaging composite —
   re-derive these rather than citing them.)
5. **Read `pvalue_adjusted` where it is present.** A raw and a corrected
   p-value can point opposite ways: CHD7 in `PMID:34324492` is `0.0068` raw and
   `0.991` after the study's own family-wise permutation correction, so a
   consumer rendering `pvalue` alone shows as significant a result the study
   reported as null. `pvalue_adjustment` names the correction, because a
   family-wise permutation correction and a Bonferroni factor are different
   claims. Both keys are `null` where the study published none — and there the
   study's own comparison count, `tests_reported` in `publications.json`, is
   what a reader has instead. **The atlas computes no correction itself.**
6. **Do not pool across studies.** The CHD literature reuses cohorts, so a
   combined p-value counts the same people twice. `case_cohorts` and
   `control_cohorts` name the collections each row drew on precisely so overlap
   is visible.

   **These ids resolve against [`cohorts.json`](#cohortsjson), and you should
   resolve them.** `case_cohorts` and `control_cohorts` publish bare strings like
   `"taa_cases"`; that file maps each to a name, a URL where one exists, and the
   caveats that qualify every number drawn from it. `taa_cases` is 777 probands
   with **thoracic aortic aneurysm, a condition outside this atlas's scope** — a
   row citing it is not a plain CHD case set, and nothing in the row says so but
   its `method_note`; the description is the caveat in full.

   Until schema `2.9` these ids resolved to nothing at all: the descriptions
   lived in the repository and reached a reader only through the gene *page*,
   under "About these cohorts". A consumer reading the JSON got the numbers
   without any of the qualifications, which is exactly the failure the cohort
   columns exist to prevent.

`n_cases` and `n_controls` are the row's own denominators: the cases and
controls its study's statistic was computed on. They may differ from the figures
a paper's abstract reports — for `PMID:42230622` they are 3,876 cases (1,471 syndromic +
2,405 non-syndromic) against 45,082 controls, while the abstract gives 4,747 and
52,881. The paper does not reconcile the two. Use the row's own denominators:
for `PMID:42230622` they are the ones its statistic was computed from.

For `PMID:34324492` they are the 6,746 cases and 14,024 controls that the
article's methods say remained after filtering, the set its permutation test was
run on. The study's S1 Table, which the article cites for the cohort as
assembled, gives 7,107 cases and 14,040 controls in its deletion column; earlier
releases published those figures as the rows' denominators.

**`n_chd_case_carriers`** (since `2.19`) says how many of `n_case_carriers` are
congenital heart disease cases, where the study's case set also holds cases of
a condition outside this atlas's scope and the study counted the CHD cases'
carriers apart; it is `null` everywhere else, which is not a count of zero.
Same unit as `n_case_carriers`, and never above it. Only `PMID:34324492` gives
one (validation refuses one on a row of any other study, `BUR024`). Its cases
combine congenital heart disease with a series of 777 sporadic thoracic aortic
aneurysm probands (the `taa_cases` cohort), and the study does not report how
many of the 777 are in the case set its rows count. Its S8 Table counts the
CHD-only carriers for 521 of the 528 genes that reached significance and
leaves the count blank for the other 7: a row carries the count where S8 gives
one for its gene, and `null` otherwise. Every other number on those rows is
from S6 Table, the file the row's `source` names; S8 Table is recorded in
[`sources.json`](#sourcesjson) as `audain2021_s8`. Where the count equals
`n_case_carriers`, every carrier the row counts is a CHD case, and a gene
page's card says so in the paragraph that holds the study's caveat about its
aneurysm probands, after the caveat.

### The bundle's `burden_combined` array: one study's own combination of its arms

On every gene bundle, `[]` where the gene has none. Each object is a study's
**own** combination of two of its burden arms, for one stratum and consequence
class. It is never a pooling across studies, which this atlas does not do (see
the `burden` array above). Today it holds PMID:40127276's JL-FDR meta-analysis
of its de novo arm (Poisson, 3,887 trios) and its case-control arm (Fisher,
alleles against gnomAD), by the joint local false discovery rate method of
Jiang & Yu 2017 (PMID:28011772).

**108 combined results reach the API**, across 36 published genes.
`mirrors/burden_combined.tsv` holds 123 rows for 41 genes: one per arm pair for
every panel gene in `mirrors/genes.tsv`, wider than the publication gate for the
reason `burden.tsv` is. Every row's arms are `burden` rows (BUR020 refuses one
that is not), matched on `study`, `gene`, `cohort_stratum`, `lesion_group` and
`consequence_class`; `combines` names their comparators, in the order the page
names them.

```json
{
  "study": "PMID:40127276",
  "cohort_stratum": "all",
  "lesion_group": null,
  "consequence_class": "damaging",
  "combines": ["mutation_model", "control_cohort"],
  "method": "joint_local_fdr",
  "pvalue": 9.72491826673901e-12,
  "local_fdr": 3.09825857398814e-10,
  "method_note": "The study's own meta-analysis of its de novo and case-control arms, …",
  "source": "sierant2025_ds6"
}
```

- **`local_fdr` is not a corrected p.** It is a local false discovery rate.
  Measured over the study's 248 panel genes, it is smaller than its own
  `pvalue` for 118 damaging, 135 loss-of-function and 129 damaging-missense
  genes, which a step-up correction can never be. Never display it under a
  "corrected p" or "q" heading, and never compare it with a burden row's
  `pvalue_adjusted`.
- **`pvalue` is the study's, transcribed.** For many of the study's genes it
  does not reproduce as Fisher's method on the two printed arm p-values (within
  10% for only 116, 133 and 144 of the 248 genes, by consequence class), so no
  formula is named for it, and it is not always stronger than the stronger arm:
  NODAL's damaging value is 7.57e-10 against a case-control 9.86e-14. The
  study's genome-wide threshold is p < 2.6 × 10⁻⁶.
- `counts.burden_combined` in [`manifest.json`](#manifestjson) is this array's
  total over the published genes.

### `independent_datasets`: a count of datasets, never a verdict

On every gene bundle and every `genes/index.json` row. It answers one question —
**how many independent datasets have looked at this gene, and what did they
find** — and it answers it by counting, never by combining. No pooled statistic
is computed here; see the previous section for why.

```json
"independent_datasets": {
  "tested": 2,
  "enriched": 2,
  "corrected": 1,
  "families": [
    { "studies": ["PMID:34324492", "PMID:42230622"], "state": "nominal",
      "members": [ { "study": "PMID:34324492", "state": "not_tested" },
                   { "study": "PMID:42230622", "state": "nominal" } ] },
    { "studies": ["PMID:40127276"], "state": "corrected",
      "members": [ { "study": "PMID:40127276", "state": "corrected" } ] }
  ]
}
```

**A family is not a study.** It is a set of studies that share at least one
sample collection, or draw on collections [`cohorts.json`](#cohortsjson) records
as overlapping (`overlaps_with`), walked transitively — so two papers that can
describe the same people count **once**. The three curated studies cite
disjoint collection ids, but since `2.19` DECIPHER (cited by PMID:34324492) and
DDD (cited by PMID:42230622) are declared overlapping — DDD links its
participants into DECIPHER, and an author of PMID:34324492 confirms its case set
does or likely does include DDD participants — so there are two families, one of
them those two studies. A family's `state` reads all of its studies' rows; a
gene page names the overlap above the two studies' blocks.

**`members` says what each study in the family says on its own** (since
`2.19`): one entry per study, sorted by study id, on every family, a one-study
family included, each with that study's `state` by the same rule over its rows
alone. A family's `state` is the strongest of its members'. It is there because
a family's state is not every member's: in the example above, TBX5's, the
family is `nominal` on PMID:42230622's rows, and PMID:34324492 has no row for
the gene at all. Attribute a state to a study from `members`, never from the
family; the glyph and the tally are the family's.

| `state` | meaning |
| --- | --- |
| `corrected` | enriched, and a row survives its own study's published correction |
| `nominal` | enriched, but nominal only — or the study published no correction at all |
| `no_enrichment` | that dataset tested the gene and detected nothing |
| `not_tested` | that dataset did not test this gene |

A family counts as enriched only if a non-synonymous row has `pvalue < 0.05`
**and points upward**. Direction is read from what the study published — the
`unbounded_above` flag, then `effect` against 1, then the published rates where
a study reports no effect at all. A significantly *depleted* row is not
agreement, and a synonymous row is never support: it is the negative control.

**Three obligations, and the first one matters most.**

1. **This is not a validity call, and it must never be rendered as one.**
   `headline_confidence` beside it is a mirrored ClinGen classification. A
   consumer that renders "0 of 2" as a verdict next to a green `definitive` chip
   tells a clinician the data contradict the classification. They do not:
   **KDM6A causes Kabuki syndrome and shows nothing in either dataset that
   tested it**, because burden tests at these cohort sizes routinely detect
   nothing for genes with overwhelming family and functional evidence. If you
   display this, display alongside it that no enrichment here is not evidence
   against a gene.
2. **Read `tested` as the denominator, never `len(families)`.** A gene absent
   from a study's panel was not examined and found wanting; it was not examined.
3. **`not_tested` must render distinguishably from `no_enrichment`.** Collapsing
   them is what turns "nobody looked" into "somebody looked and found nothing".

`families` is always present, ordered deterministically, and **always carries one
entry per cohort family in the whole corpus** — two today. A gene no study
reported gets two `not_tested` entries, not an empty array; `[]` occurs only if
the corpus has no burden data at all. That is deliberate: every gene shows the
same slots, so a dataset that did not test this gene is visible as an absence
rather than as a shorter list.

**The browse page heads this column "rare variants across studies", not
`independent_datasets`.** The key names what is counted for a program; the
header names what it is for a reader. Both avoid a verdict word on purpose —
"replicated in" was the shorter candidate and was rejected, because for a gene
showing 0 of 2 it reads as "not replicated", which is a claim the data do not
make.

### The bundle's `expression_profile` object: a developmental transcriptome, never a contrast

On every gene bundle. Always present, and empty (`{"datasets": []}`) for a
gene no profile dataset's mirror covers — the same rule `omics` and `burden`
above already keep, so "no profile dataset says anything about this gene"
cannot be confused with "the build dropped it".

```json
{ "expression_profile": { "datasets": [] } }
```

**That is still the real shape for a gene no profile dataset's rows mention.**
`mirrors/profiles/E-MTAB-6814.tsv` is committed, and `manifest.json`'s
`counts.profile_genes`/`counts.profile_datasets` say how far that reaches:
119 genes carry a developmental expression profile, across 1 dataset. Most
published genes carry the populated shape today, not the empty one above.
The nested shape below is illustrative — restructured for
readability (a real gene's `stages` runs to a dozen or more entries, one per
token the dataset declares) rather than copied verbatim from one bundle — but
every value in its `phase` block is real: it is `assign_phase`'s actual,
reproducible answer for a 7-elapsed-week stage against the boundaries
`curation/cardiac_phases.yaml` curates today. The `stage` token, the `unit`
and `n_genes` are real too: stage tokens are the dataset's own long-form
strings, this dataset reports `tpm`, and its grids were built from a 58,735-gene
transcriptome. The `tau`, the abundances and the percentiles are the invented
part.

```json
{
  "expression_profile": {
    "datasets": [
      {
        "dataset": "E-MTAB-6814",
        "quantile_shard": "omics/profile_quantiles/E-MTAB-6814.json",
        "stages": [
          {
            "stage": "7 week post conception",
            "phase": { "outcome": "matched",
                       "phase_ids": ["ventricular_septum_morphogenesis",
                                     "heart_valve_morphogenesis"],
                       "reason": null },
            "specificity": { "tau": 0.62, "scale": "log2(x+1)", "method": "…",
                              "tissues": ["heart", "kidney", "liver"],
                              "n_tissues": 3, "highest_in": "heart",
                              "medians": { "heart": 42.0, "kidney": 8.0, "liver": 6.0 } },
            "specificity_unavailable_reason": null,
            "tissues": [
              { "tissue": "heart", "median_abundance": 42.0, "unit": "tpm",
                "n_samples": 3,
                "placement": { "q25_percentile": 44, "median_percentile": 50,
                                "q75_percentile": 56, "median_abundance": 42.0,
                                "unit": "tpm", "n_samples": 3, "n_genes": 58735,
                                "method": "lowest percentile of a tied breakpoint (bisect_left)" },
                "not_placed_reason": null }
            ]
          }
        ]
      }
    ]
  }
}
```

**Two mirror tables feed this, and each reaches a consumer a different way.**
`mirrors/profiles/<accession>.tsv` (`dataset`, `gene`, `tissue`, `stage`,
`median_abundance`, `unit`, `q25`, `q75`, `n_samples`) is one of the four
tables [`omics/<modality>/<accession>.json`](#omicsmodalityaccessionjson)
above shards like any other — `profiles` is one of that section's four
modalities — so its raw rows reach `omics/profiles/<accession>.json`
regardless of gene publication, and `median_abundance`/`unit`/`n_samples`
below are read straight from it. `mirrors/profile_quantiles/<accession>.tsv`
never goes through that mechanism at all: it has no gene column, so
`build_omics` skips it outright, and it is published only as
[`omics/profile_quantiles/<accession>.json`](#omicsprofile_quantilesaccessionjson)
(below), one full shard per dataset.

- `datasets` is one entry per profile dataset whose mirror mentions this
  gene. `quantile_shard` names that accession's
  `omics/profile_quantiles/<accession>.json` payload, or `null` if the build
  wrote no shard for it.
- `stages` is one entry per stage token the mirror's rows use for this
  (gene, dataset) pair — `stage` is `null` for a measurement with no
  developmental stage recorded at all.

  **The array is in chronological order, and was not before schema `2.12`.**
  It is sorted on the dataset record's own curated
  [`order`](#datasetsjson) — not on `wpc`, which is `null` for every
  post-natal stage, and not on the stage token. A token the dataset does not
  declare sorts after every one it does; the `null` stage sorts last. The
  entries themselves do not carry `order`; the array's sequence is what
  publishes it, so read the series as it arrives.

  Through schema `2.11` this array came back **alphabetically**: `4 week post
  conception` published eighth, after `19 week post conception`, and `elderly`
  second of the eight post-natal stages. That was deterministic and
  reproducible between builds — and it was not a chronology. If you plotted
  these entries in array order against an earlier release you drew a
  developmental trajectory with its axis shuffled, and the same code now draws
  a different, correct picture. If you sorted the stage tokens yourself to
  work around it, stop: sorting the tokens is what the defect did.
- `phase` places the stage in the curated cardiac morphogenetic window
  (`curation/cardiac_phases.yaml`):

  | `phase.outcome` | meaning |
  | --- | --- |
  | `matched` | the stage's own developmental age (`wpc`) falls inside one or more curated cardiac phases; `phase_ids` names all of them |
  | `outside_window` | a real age exists, but no curated phase covers it |
  | `post_natal` | the stage has no developmental age at all |
  | `undeclared` | the dataset's own record does not declare this stage token |
  | `null` | the row itself carries no stage token; `reason` says so |

  **`phase_ids` may name more than one phase, and often does.** Human cardiac
  morphogenesis runs several processes concurrently — at 6 elapsed weeks post
  conception, atrial septation, ventricular septation and outflow tract
  septation are all underway at once — so a stage legitimately matches every
  phase whose window contains its `wpc`, never just the nearest or the first
  declared. `curation/cardiac_phases.yaml` curates six phases transcribed from
  Buijtendijk et al. 2020 (PMID:32048790), each carrying the Gene Ontology
  term and the Carnegie stage(s)/HsapDv id(s) its boundaries were read from —
  not published in this bundle field, only in the curation file itself, so a
  consumer auditing a boundary starts there. `phase_ids` is empty exactly
  when `outcome` is not `matched`.

  **Each phase's end carries an `end_basis` in the curation file, and one of
  the six ends where this atlas concluded rather than where the source
  says.** `stated` is a boundary the source gives. `derived` is one the atlas
  concluded from the phase's own Gene Ontology definition read against a
  boundary the source does give, and is kept distinct so a reader can tell
  the two apart. `not_stated` is a phase whose end nothing gives: its end is
  left `null` and the phase is excluded from every stage's match at every
  `wpc`, rather than asserted to run on past a point no source supports.
  `heart_looping` is the one `derived` phase: the source states when looping
  starts (CS10) and never when it ends, but GO:0001947, the term the phase
  carries, defines looping as ending "preceding septation", and the source
  dates septation's start at CS12 — so looping's end is CS12, 3.71 wpc
  (`HsapDv:0000019`). It carried `not_stated` until 2026-08-24. None of the
  six is `not_stated` today.

  **Two of the six never appear in a published `phase_ids` today, and that
  is a fact about the data, not a rule.** `embryonic_heart_tube_morphogenesis`
  ends at 3.29 wpc and `heart_looping` at 3.71, and the earliest stage any
  published dataset samples is 4 wpc, so no published stage falls inside
  either. A dataset sampling earlier stages would be the first to show them;
  a consumer reading `curation/cardiac_phases.yaml` sees all six.
- `specificity` is Yanai's τ (tau) over every organ a dataset sampled at one
  stage, or `null` with `specificity_unavailable_reason` naming why:

  | `specificity_unavailable_reason` | meaning |
  | --- | --- |
  | `dataset_not_registered` | the mirror names a dataset accession with no curated record at all |
  | `detection_floor_undeclared` | the dataset record exists but declares no detection floor |
  | `one_organ_sampled` | fewer than two organs were sampled at this stage — τ's denominator is `n − 1` |
  | `peak_below_detection_floor` | every sampled organ's median is below the dataset's own floor |
  | `undefined` | a residual this atlas cannot itself further diagnose |

  `tissue`-level entries carry a **related but not identical** vocabulary
  through `not_placed_reason` — the two share `dataset_not_registered` and
  `detection_floor_undeclared`, but a single tissue can never be short of
  organs or undefined the way a whole stage's τ can, and gains two reasons of
  its own instead: `no_quantile_grid` (no complete 101-point breakpoint grid
  exists for this cell) and `below_detection_floor` (a grid and a floor both
  exist; the median itself falls below it). Four possible values, never the
  stage-level five.
- **τ is published bare — no adjective, no band.** A value of 0.71 is not
  glossed "intermediate" or "specific": choosing a threshold would be a
  classification this atlas authors about how tissue-specific a gene is,
  which its governing rule forbids for every derived figure in this layer, the
  same way D12 forbids the atlas authoring a gene-disease validity call. If
  you display τ, display the number.
- **τ is blind to *where* a gene peaks — read `highest_in`, never the number
  alone.** τ measures concentration, not location. Measured: a gene at heart
  20, liver 200 and five other organs at 5 (rpkm; seven organs sampled) scores
  τ = 0.963 computed on the raw values and 0.623 on the `log2(x+1)` scale this
  atlas actually publishes — and in both, the peak organ is **liver**, not
  heart. A consumer who reads a moderately high τ as "heart-preferential" for
  a gene concentrated in liver states the opposite of the truth. `highest_in`
  is `null` when two or more organs tie exactly at the peak, where naming one
  of them would be arbitrary.
- `scale` (`"log2(x+1)"`) and `method` travel beside every τ so a consumer can
  re-derive it without reading this atlas's source rather than trusting a
  paraphrase of it. `method` reads exactly: tau (Yanai et al. 2005): the sum
  over organs of (1 - x_i/x_max) divided by (n - 1), where n is the number of
  organs and x = log2(median+1); a negative median is clamped to 0 before the
  transform, and every organ's raw median is used even below the dataset's
  detection floor. `medians` is τ's own input, published exactly
  as measured — including any organ a gene bundle's own `omics` preview may
  have already dropped from `top` — so τ's inputs stay reachable from the one
  payload that publishes τ (design decision D39(b)).
- `placement` is one gene's percentile band in one (dataset, tissue, stage)
  cell: `q25_percentile`/`median_percentile`/`q75_percentile`, the outer two
  `null` below `n_samples = 3`, where a quartile of two points is not a
  quartile. `method`
  (`"lowest percentile of a tied breakpoint (bisect_left)"`) names the
  tie-break — a value tied with a run of identical breakpoints reports the
  *lowest* percentile in that run, so an unexpressed gene never reads above
  the median it shares with every other unexpressed one. `n_genes` is the
  size of the *source* transcriptome that cell's grid was built from — tens
  of thousands of genes, never the 119 this atlas publishes — because "top 4%"
  names no population without its denominator.

**Percentiles are not comparable across organs, and — more importantly — not
comparable across developmental stages.** Each organ transcribes a different
fraction of the whole gene universe, so "top 10% in liver" and "top 10% in
heart" describe two different reference distributions, not one ranking read
two ways. And a reference distribution's own shape changes over
development, in every organ, so a gene with flat, unchanging absolute
abundance across development can still show a *moving* percentile from one
stage to the next, purely because everything else in the transcriptome is
moving around it. The cross-stage caveat matters more here, because this
layer's headline question — is this gene's cardiac expression distinctive,
and does that change over development — is asked *within* one organ, across
its own stages: reading a percentile trend without this caveat can report
development where none occurred, or miss it where it did.

**Why the quantile grid is published at all, and exactly how far that goes.**
A percentile's input is the whole source transcriptome — tens of thousands of
genes — and this atlas publishes 119 of them; design decision D32 forbids
re-hosting the matrix a percentile was read against. Without a published grid,
`median_percentile` would be a number nobody outside this atlas could check.
[`omics/profile_quantiles/<accession>.json`](#omicsprofile_quantilesaccessionjson)
(below) publishes exactly the 101 breakpoints instead — which is a
**one-level guarantee, not an unbounded one**: a gene's percentile is
re-derivable from its own `median_abundance` and the published breakpoints;
the breakpoints themselves are not re-derivable, because they were read off
the whole matrix D32 keeps unpublished. State it as the trade it is —
auditable one level down, not provably correct all the way back to the source
data.

### The bundle's `cell_type_expression` object: a cell-type profile, never a contrast

On every gene bundle. Always present, and `{"datasets": []}` for a gene no
cell-type dataset covers; no cell-type mirror is committed yet, so that is
every bundle today and `counts.cell_type_genes` and `counts.cell_type_datasets`
are both 0.

A low or absent value is not evidence against a gene: it describes the cells
sampled, not the gene's role. Per-region values pool 9 to 15 weeks of
development, so they are not a stage series. A cell type or region that is
absent means no heart had at least 10 cells there. `below_floor` marks a
median under the dataset's `detection_floor`. `specificity` is tau over the
pooled medians, with its inputs; its `cell_types` are those with a pooled value,
in the dataset's declared order. It is `null` when fewer than two cell types
have a pooled value, or when every pooled median is below the floor.

```json
{"datasets": [{"dataset": "hoc", "unit": "cpm", "detection_floor": 1.0, "n_donors": 8,
  "specificity": {"tau": 0.591, "scale": "log2(x+1)", "method": "<cell_types.CELL_TYPE_TAU_METHOD>",
    "cell_types": ["aCM", "vCM", "ncCM", "Fibro", "SMC", "Epicardial", "Endocardial",
                   "BEC", "LEC", "P-RBC", "WBC", "Neuronal"],
    "n_cell_types": 12, "highest_in": "aCM",
    "medians": {"aCM": 220.67, "vCM": 41.48, "ncCM": 159.96, "Fibro": 11.16, "SMC": 4.26,
                "Epicardial": 20.75, "Endocardial": 8.3, "BEC": 2.7, "LEC": 0.41,
                "P-RBC": 0.2, "WBC": 7.38, "Neuronal": 6.84}},
  "cell_types": [{"cell_type": "aCM", "label": "Atrial cardiomyocytes",
    "pooled": {"median": 220.67, "q25": 194.18, "q75": 241.99, "n_donors": 8,
               "n_cells": 20268, "fraction_expressing": 0.837, "below_floor": false},
    "regions": [{"region": "LA", "median": 228.06, "q25": 198.51, "q75": 241.66,
                 "n_donors": 8, "n_cells": 7437, "fraction_expressing": 0.823,
                 "below_floor": false}]}]}]}
```

## `genes/<slug>.html`

One page per published gene, rendering that gene's bundle for a reader.

**The summary column is a card**, written for a clinician reading it first.
Top to bottom:

1. The gene's symbol, its approved name, locus and MANE Select transcript, then
   a row of links out: gnomAD and ClinVar first and ClinGen and Ensembl after
   them, all four on every page, and UniProt where the registry has an
   accession.
2. The HGNC id and the headline classification as a chip — with a `conflicting
   evidence` chip beside it when `has_conflicting_evidence` is set and a
   `sources disagree` chip when `has_source_discordance` is. Both are flags on
   the gene, not on a disease: the disputing record and the supporting one may
   be for different diseases, as they are on every gene that sets them today,
   and [`methods.html#validity`](#methodshtml) says so.
3. Up to three rows. The disease that grade is for (or, for a gene admitted on
   submitter agreement, the submitters who admitted it). **Also listed by** —
   how many other institutions in `asserted_by.institutions` list the gene,
   ClinGen and the admitting submitters excluded; omitted when there are none.
   **Inheritance** — every distinct mode of the records that carry the
   warrant: the ClinGen records at the headline grade, or the admitting
   submitters' records at `Limited` or better, in the validity table's
   vocabulary. Every mode, not one: a panel can grade one gene at one rung under
   two modes.
4. For a graded gene, one line glossing its headline grade in the grade key's
   own words — on a `limited` gene, that *Limited is not a weak yes*. Where
   `has_conflicting_evidence` is true, a line under it (or in its place, for a
   gene no panel graded) naming who disputes or refutes the gene and for what,
   from the contested records the validity table lists: each authority's own
   term, the diseases it applies it to, and the authority — the expert panel
   for a ClinGen record, the submitter for a GenCC one — e.g. "Disputed for
   Costello syndrome and cardiofaciocutaneous syndrome by the RASopathy Gene
   Curation Expert Panel; see the validity table." The flag itself does not
   say which disease; on every gene that sets it today the dispute is for a
   different disease from the one the headline grade is for. Where
   `has_no_association_report` is true, a line after it names each authority
   that reports no known disease relationship, the diseases, and the year of
   its record, from the records the validity table lists — GDF1's reads
   "Illumina reports no known disease relationship for congenital heart
   defects, multiple types, 6 (2021); see the validity table."
5. Where every disease the headline grade is for is in scope through Orphanet's
   phenotype frequencies, the reason: "Graded for Noonan syndrome as a whole
   syndrome; in scope because Orphanet records a structural heart malformation
   in at least 30% of cases." ("… as whole syndromes" where the grade is for
   more than one.) For a subtype, whose warrant is its syndrome's, the clause
   names the syndrome the frequency was recorded for ("… of cases of
   Adams-Oliver syndrome"). Where instead every such disease is a term
   `scope.json` marks `syndrome` and a ClinGen expert panel scoped, the line
   names the panel: "Graded for Holt-Oram syndrome as a whole syndrome; in
   scope on its curation by ClinGen's Syndromic Disorders Gene Curation Expert
   Panel." (since `2.19`; ten cards). The panel is the term's `attributed_to`,
   the one whose curation scopes the term; where it is not a panel that graded
   this gene for the disease, the line names the genes it did grade for it —
   EVC2's: "… in scope on ClinGen's Kidney Cystic and Ciliopathy Disorders Gene
   Curation Expert Panel's curation of it for EVC." Several diseases each name
   their own panel.
6. **Rare variants**: one line per study and analysis with a result that
   survives its own study's correction — the rule the evidence matrix's
   `corrected` state uses, so the two cannot disagree — naming the first such
   row in the table's own order (the `damaging` union first) with its counts in
   their own unit, then the group's other surviving results by name ("also
   surviving: loss-of-function; missense (damaging)") — named, not counted,
   because the union's components are not further findings. A line from a
   study by an author of this atlas says so ("study by an atlas author"), as
   the study's block does. Each line ends with a link to that study's block on
   the same page (`#study-<PMID digits>`, the `id` of the study's heading):
   `caveats` where any of the group's surviving rows carries a `method_note`,
   `details` otherwise, with an accessible name naming the study. Not the
   most significant row: a line keeps its study table's reading order, while
   a matrix cell names the most significant of the rows behind its state (see
   the rare-variant section below). Where the study's publication records
   `correction_genes`, the line's row also carries its corrected p and what it
   was corrected over ("q 0.00369 within its 248-gene panel"), because a q over
   a panel is not genome-wide; where it records `genome_wide_threshold`, every
   surviving row on the line whose raw p is below it is marked, the line's own
   row as "genome-wide (p 1.07e-08)" and a named one as "loss-of-function
   (genome-wide, p 1.35e-09)" (TBX5's de novo line). Today only PMID:40127276
   records either. Beneath the lines, once
   per study with a line, that study's curated `summary_caveat` from
   `publications.json`, where it has one; the full caveat stays in the
   study's block. After it, and never without it, where the line's row
   carries `n_chd_case_carriers`, what that count says of the carriers the
   line shows: "For this gene, the study's S8 Table counts every carrier among
   its CHD cases, so none is an aneurysm proband." where it equals the row's
   case carriers, and both counts where it is fewer. Only PMID:34324492's rows
   carry it, so only its caveat is ever followed by one. It speaks of the
   carriers alone: the line's denominator is the case set the study analysed,
   aneurysm probands among it.

   **A study that published no corrected p for the gene gets a line of its
   own**, whether or not another study's results survive: "<study> publishes
   no corrected p for this gene; judge its results against its scan of 138,609
   tests (see its block)", with the scan size from the publication's
   `tests_reported`; where none is recorded the line ends "see its block".
   Under the one rule `methods.html#rare-variants` states, such a study's rows
   survive nothing, so none of them is listed on the card; the evidence matrix
   still marks one as supporting where it points to enrichment in cases and
   its p is below 0.05. Whether a study published a correction is read per
   gene, from its rows.

   With no surviving row, how many independent datasets (cohort families, not
   studies) tested the gene, qualified as the front page qualifies its count —
   "(as far as the cohort registry records)", straight after the count — then
   declaring any study among them by an author of
   this atlas — "(a study by an atlas author)" where one study tested it,
   "(each study by an atlas author)", or "(1 of the 2 studies by an atlas
   author)", counted in studies — and, where some study published a
   correction for it, that none survives a correction its study published: "Tested in 2
   independent datasets (as far as the cohort registry records) (1 of the 2
   studies by an atlas author); none survives a correction its study
   published." With no row at
   all, that no rare-variant study here reports a result. The power caveat — no enrichment at these cohort sizes is not
   evidence against a gene — follows whenever the gene has rows.
7. One sentence saying what the atlas itself has curated (below).
8. *How to read this page*, a link to the top of [`methods.html`](#methodshtml),
   and a link to the gene's own JSON.

The page no longer prints `validity_state`, `atlas_curation`, the lesion groups
or the assertion, functional-record, curated-publication and burden-row counts
as rows; every one of them is in the bundle, which remains the record.

Beside the card, the main column opens with the same one-line short form of the
scope rule as the gene index, linking [`methods.html#scope`](#methodshtml),
which states the rule in full and carries the grade key; before this release,
every gene page carried both. Below it, the mirrored validity table in two groups —
**Expert panel (ClinGen)** first, then **Other curators (GenCC)** — each row
giving who curated it, the disease, the mode of inheritance in one vocabulary
(ClinGen's `AD` reads as GenCC's "Autosomal dominant"), the authority's own
`classification_term` **verbatim** rather than the rung this atlas maps it onto,
the date as
`YYYY-MM-DD`, and a link to the upstream report where one is published, labelled
with its SOP for a ClinGen record. ClinGen also submits its classifications to
GenCC; the page lists such a copy only where its mapped grade, date, inheritance or
disease differs from the panel's own record, or where it carries no date or has
no panel record to match. The bundle keeps every record. The table's intro says
that only records for in-scope diseases are listed, so the table may not be every
classification the gene has — "may not", because many published genes have no
ClinGen or GenCC record for any other disease — and links
[`methods.html#validity`](#methodshtml). Neither the page nor the bundle carries
a record for a disease outside scope. A curated gene's page
follows it with the curated evidence — a line saying each assertion is this
atlas's own curation and changes no classification above, linking
[`methods.html#curated-evidence`](#methodshtml); each assertion, its evidence
items with their class, strength and summary; and the publications that
evidence cites.

The rare-variant section opens with one line attributing every count and
statistic in its study tables to the upstream study — "Every count and statistic
in the study tables below is an upstream study's, republished as published; the
atlas computes none of them." — and a link to
[`methods.html#rare-variants`](#methodshtml), *How to read these results*. The
line is scoped to the study
tables because the evidence matrix beneath it shows tallies the atlas derives
(how many rows a cell summarises). A matrix row is a cohort family and names its
studies; in a family of two or more, a study with no row for the gene is
marked "(no row for this gene)", since the family's cells cannot say which
study they come from. A `corrected` or `nominal` matrix cell also
names one row: the most significant of the rows that give the cell that state —
the surviving rows under `corrected`, the supporting rows under `nominal` — with
its effect and one p-value. An unbounded odds ratio keeps its lower bound there,
e.g. `OR ∞ (CI 0.786–∞)`, because a bare `∞` reads as a strong result whatever
the interval. The reading notes this section folded before this release — the
count units, what an absent row means, the synonymous
negative control — are on that page now, as is why a damaging composite row is
shown. What stays visible beside the matrix, the plots and the tables, outside
every fold:

- under the matrix, its legend — rows are independent cohort families, and no
  enrichment at these cohort sizes is not evidence against a gene — and that
  no statistic is pooled across studies: "The atlas computes no pooled
  statistic across studies: pooling two that share a sample collection would
  count the same people twice, and pooling two that draw on collections this
  atlas records as overlapping could." — and, where the card lists a surviving row,
  that a cell and the card can lead with different rows: "A cell names the
  most significant of the rows behind its state; the card lists a study's
  surviving rows with a damaging union before its parts, so the two can lead
  with different rows.";
- for two of the gene's studies that are not independent, a notice naming why,
  as this atlas has recorded it: a sample collection both draw on ("… so their
  results describe partly the same people and must not be combined"), or
  collections [`cohorts.json`](#cohortsjson) records as overlapping
  (`overlaps_with`), one per study ("… collections this atlas records as
  overlapping, so their results may describe partly the same people and must
  not be combined"). On today's corpus no two studies share a collection, and
  every such notice is the second kind: PMID:34324492's DECIPHER beside
  PMID:42230622's DDD, on every page carrying both;
- on a block from a study by an author of this atlas, a one-line declaration:
  "Declaration: this study is by an author of this atlas; its results are
  republished as published and read by the same rule as every other
  study's.";
- the multiple-testing line. For a study with any row it published no
  corrected p for, the study's number of tests (`tests_reported` in
  `publications.json`, where one is recorded) and the Bonferroni threshold
  over it. `tests_reported` is named on a gene page only for such a study,
  here and in the card's line for it. For a study that
  corrected every row, that its corrected p carries its own correction for
  multiple testing, and no count;
- where a block has a damaging composite, that it is the union of the two rows
  below it;
- a study's own combination of its arms, where it publishes one
  ([`burden_combined`](#the-bundles-burden_combined-array-one-studys-own-combination-of-its-arms));
- where a study publishes no effect measure for the gene, a sentence saying so
  in place of a plot, with its counts and p-values in the table;
- every row's `method_note` footnote.

Each forest plot — one per study and effect measure — shows one line of its
caption: the axis, naming the effect measure, its range on a logarithmic scale
and what the line at 1 means for that measure. Where a `union`-tagged row
contains a row the panel does not plot, the same line names that row and the
one missing from it — "The **all cases · loss-of-function** row is a union,
and this study published no non-syndromic · loss-of-function row for this
gene." — or, for a row the study did publish and the panel does not draw,
that it is not plotted here. The rest of its key is folded under *How to read this chart*:
what each mark means — filled and hollow markers, hatched bars, arrows, muted
and `union`-tagged rows — each only where the panel draws one, and why there
is no summary diamond, ending with a link, *How to read rare-variant results*,
to [`methods.html#rare-variants`](#methodshtml), which states the key in
general.

The developmental expression section, last on the page, opens with one visible
line: "Every measurement below is from bulk tissue: a low or absent
figure is not evidence a gene is unimportant to heart development." It links
[`methods.html#expression`](#methodshtml) as *How to read these figures*; that
section states the argument in full with the reading notes this section folded
before this release — how a median is placed against its organ's percentile
grid, why none is placed below a dataset's detection floor, and why percentiles
are not comparable across organs or developmental stages — and defines τ and
the developmental phases each stage names. Each chart — a phase-banded
trajectory per cardiac organ where one can be drawn, and the organ small
multiples — shows one line of its caption: its axes, in stage order rather than
elapsed time, and the range of its abundance axis. Under a trajectory that ticks
a stage below the dataset's detection floor, the same line says what the ticks
mark, with the floor's value, and names any other reason a tick on it was
drawn for. The rest of its key is folded under *How to read
this chart*, ending with a link, *How to read expression figures*, to
`methods.html#expression`. Under a trajectory the fold always names the shaded
phases and their source; why the line breaks where this dataset has no row for
the organ, and what a tick means on a figure with no tick for the floor, it
adds only where the figure shows one. Under the small multiples it always says why
the panels share one axis, which organ is drawn heavier as this atlas's
declared cardiac tissue (or, where this dataset placed nothing there, that no
panel is drawn heavier), that only medians placed against the percentile grid
are drawn, and that every figure is in the table below; why a line breaks, and
which sampled organs have no panel, it adds only where the panels show one.

For the 117 genes published today with no curation here, the card's curation
sentence says exactly that:

> Not yet curated by the atlas: every classification on this page is a ClinGen
> expert panel's or a named curator's.

That wording renders on 96 of those 117 pages, the genes a ClinGen panel holds
an in-scope record for. The 21 admitted on submitter agreement have no panel
classification to describe, and say so instead:

> Not yet curated by the atlas, and no ClinGen expert panel has graded this
> gene for any of the diseases in scope here: it is published because the Gene
> Curation Coalition submitters it is admitted by independently assert it.

The scope in that second sentence is part of the claim. Unscoped it is false for
five of the 21: ClinGen grades ACTC1 Definitive for hypertrophic cardiomyopathy
and ELN Definitive for cutis laxa, for instance, and those diseases are outside
this atlas's remit, so their records are not on the page to contradict it.

`validity_state` selects between the two sentences, not `headline_confidence is
null`: the two agree on every gene published today, and two figures that are
equal are one figure to every test.

The sentence is there instead of the curation simply going unmentioned. A
missing evidence section is indistinguishable from "the atlas looked and found
nothing", and a reader deciding what a gene means clinically must not have to
infer which.

Read it precisely, because it is narrower than it first appears. **A lesion
assertion** is what decides it, not evidence in general: `atlas_curation` is
derived from curated `LesionAssertion` records alone, so a gene can carry
functional-evidence records the atlas curated and still report
`not_yet_curated`. Such a card opens "No lesion assertion yet curated by the
atlas" instead, and adds "The atlas has curated N functional record(s),
published in this gene's JSON." — naming those records as the atlas's own work
rather than denying them in the sentence that counts them. No published gene is
in that state today.

A gene page carries no script at all — every value on it was rendered at build
time — and, like the browse page, makes no external request. The bundle remains
the machine-readable contract: the page's markup is not one.

**No payload carries the path to a gene's page.** A browse row's `bundle` is
the JSON; the route to the page is the link on `genes/index.html`, which is why
the slug rule stays an implementation detail (see
[Reading this API](#reading-this-api)).

## `omics/<modality>/<accession>.json`

A gene bundle's `omics` maps a modality — `expression`, `profiles`, `proteomics`
or `phospho` — to a summary of that gene's rows. The committed corpus mirrors
`profiles` alone today, so that is the only key a bundle carries; TBX5's, with
`top` abridged to its first row:

```json
{
  "omics": {
    "profiles": {
      "count": 119,
      "shards": ["omics/profiles/E-MTAB-6814.json"],
      "top": [
        { "dataset": "E-MTAB-6814", "gene": "HGNC:11604", "genes": ["HGNC:11604"],
          "tissue": "heart", "stage": "16 week post conception",
          "median_abundance": 184.0, "q25": 122.0, "q75": 279.0, "unit": "tpm",
          "n_samples": 3, "percentile": 100 },
        …
      ]
    }
  }
}
```

- `count` is every row about the gene, across every shard listed.
- `top` holds the same row objects the shard does, `genes` included.
- `shards` are the files holding them. Each shard is
  `{"table": "<modality>", "rows": [ … ]}` — the mirror rows, each with one field
  added by the build, `genes` (see below). A `profiles` row gains a second,
  `percentile`: that row's median placed on its cell's quantile grid — the same
  figure as `placement.median_percentile` in the gene's
  [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  — or `null` where the build could not place it.
- **`top` is capped at 25 rows.** For `expression`, `proteomics` and `phospho`
  it is a preview, ranked by significance (ascending FDR, ties broken by the
  table's own sort order), not a page of results. `count` is frequently
  larger, and the cap is not carried in the payload, so do not infer
  completeness from `len(top)`.
- **`profiles` ranks `top` differently, because that table has no
  significance column at all** — it is an abundance table, not a contrast.
  The slice is *stratified*: the cardiac series leads, ranked by the derived
  `placement.median_percentile`
  (see [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  above), but slots are **reserved for the non-cardiac tissues**, up to half of `top`.

  With `top` capped at 25 that reservation is at most 12, so every other
  tissue present is represented only while there are 12 or fewer of them —
  measured, and true of this dataset, which has six. Past that the cardiac
  series can consume its budget and the remaining tissues share what is left,
  so do not read the reservation as a guarantee that every tissue appears.

  **The reservation exists because of τ.** A gene's `tau` figure (also
  documented there) is computed over every organ a dataset sampled at one
  developmental stage, and design decision D39(b) requires τ's own inputs to
  stay reachable from the same payload that publishes τ. Without the
  reservation, a bulk developmental atlas sampling many stages per organ
  fills every one of the 25 rows with the cardiac series alone — measured:
  ranking `profiles` the way the other three tables are ranked put 0 of 14
  heart rows in a 25-row slice at 14 stages per organ (every row ties at "no
  FDR", so the tie-break falls back to an alphabetical tissue order that
  heart loses), and ranking cardiac-first with no reservation at all
  over-corrected, leaving 0 of 6 comparison organs in a 25-row slice at 23
  stages per organ. **A consumer who does not know this will read the
  reserved rows as a ranking bug**: they are exactly the organs a gene's own
  `tau` was computed over, shown so that figure is auditable from the same
  payload that publishes it.

**To get the rows a bundle counted, filter the shard on `genes`.**

Every shard row carries `genes`, a list of the HGNC ids that row is evidence
about. *Illustrative* — no proteomics or phospho mirror is committed, so this
row is constructed to show the case the field exists for, a row keyed on a
protein rather than a gene:

```json
{ "dataset": "PXD012345", "protein": "Q99593", "position": 100, "genes": ["HGNC:11604"] }
```

```js
const shard = await (await fetch(summary.shards[0])).json();
const mine = shard.rows.filter(row => row.genes.includes("HGNC:11604"));
// mine.length === summary.count, when the gene has one shard
```

That equality is the point of the field, and it holds by construction: `count` is
derived from the same attribution the rows publish, computed once. It is a
**list** because one protein accession can belong to several genes — a histone
cluster, for instance — and a single-valued field would silently drop all but one.

`genes` is present on every modality, including `expression` and `profiles` whose
rows already carry their own `gene` column, so a consumer filters one way
everywhere. Where the two exist side by side they agree; `genes` is the one
`count` is built from.

`mirrors/genes.tsv` holds 185 rows — exactly the genes ClinGen or GenCC curates
within the scope `curation/chd_scope.yaml` declares, of which 119 clear the
publication gate — across ten columns. The six identity columns (`name`,
`locus`, `ensembl_gene`, `ncbi_gene`, `uniprot`, `mane_select`) are populated on
all 185 rows (re-measured 2026-10-25). **All ten columns reach a published byte
since `2.17`**; before it, `ensembl_gene`, `ncbi_gene`, `locus`, `uniprot` and
`mane_select` reached none. `hgnc_id` and `symbol` are published on every gene
bundle and index row; `name`, `aliases` and — since `2.13` — `prev_symbols` are
published as `terms` in `search/index.json.gz`; and `name`, `locus`,
`ensembl_gene`, `ncbi_gene`, `uniprot` and `mane_select` are the bundle's
[`identity`](#genesslugjson) object, for each published gene.

**`prev_symbols` was a lookup input only until `2.13`, and this paragraph said
so.** It carries the HGNC symbols a gene has had withdrawn, and it arrived for
the join: the resolver in `chd_atlas/genes.py` ranks an approved symbol above a
current alias above a retired name, so a source still writing an old name lands
on the right gene. Search did not read it, which meant **35 of the 92 published
genes held a retired symbol nothing on this site could find** — `TBX5`←`HOS`,
`VEGFA`←`VEGF`, `ZIC3`←`HTX1`. `2.13` publishes them: 46 terms across those 35
records, and **zero of them already named anything else in the index**, so no
query that returned one gene before returns two by accident.

**One does return two, deliberately.** `MADH7` is published on both SMAD6
(HGNC:6772) and SMAD7 (HGNC:6773), which each retired it — and the resolver
*refuses* that same symbol rather than picking a gene. The two are not in
conflict. Resolution asks which single gene a string denotes, where a guess
attributes one gene's evidence to another; search asks what a reader might be
looking for, where showing both is the only honest answer. A consumer building
its own index from `terms` should expect a term to name more than one record,
as it already could for a shared alias.

**Since `2.17` each bundle carries four cross-references for its gene, and no
others.** Its [`identity`](#genesslugjson) object holds the Ensembl gene id, the
NCBI Gene id, the UniProt accession and the MANE Select transcript, plus the
gene's cytogenetic band, mirrored from HGNC rather than mapped by this atlas —
with one exception the atlas resolves itself. HGNC lists two UniProt accessions
for RBM10 (HGNC:9896), `P0DW28|P98175`, in alphabetical order, and the first is
MINAS-60, a 130-residue microprotein read from an alternative frame inside
RBM10. The converter pins RBM10's own protein, `P98175`, and refuses to run if
any other gene gains a second accession, until a curator has chosen one. Until
`2.17` the API carried none of these. No other identifier for the gene — an
OMIM number, a transcript other than MANE Select — is published, and a consumer
that needs one must resolve `hgnc_id` against HGNC itself. The atlas identifies
a gene by its HGNC id and leaves every further mapping to the authority that
maintains it.

## `omics/profile_quantiles/<accession>.json`

The percentile grid a `placement`
(see [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
above) was read against — the other half of the trade design decision D39(b)
makes. This is **not** one of the four modalities in the section above: it has
no gene column at all — there is no gene to attribute one breakpoint row to —
so `build_omics` skips it outright (that function shards only the tables its
own `_GENE_COLUMN` map names). A separate builder
(`profiles.build_profile_quantiles`) emits it instead, to this separate,
sibling path under `omics/`.

```json
{
  "table": "profile_quantiles",
  "rows": [
    { "dataset": "E-MTAB-6814", "tissue": "heart",
      "stage": "7 week post conception", "percentile": 50,
      "value": 0.0, "unit": "tpm", "n_genes": 58735 }
  ]
}
```

Exactly the `mirrors/profile_quantiles/<accession>.tsv` rows for that dataset,
republished verbatim — `dataset`, `tissue`, `stage`, `percentile`, `value`,
`unit`, `n_genes` — sorted by `(dataset, tissue, stage, percentile)`. A
complete grid is 101 rows for one (dataset, tissue, stage) cell, `percentile`
running 0 through 100; `n_genes` is the size of the source transcriptome that
grid was built from, the same figure a `placement` in that cell publishes
under its own `n_genes`.

A gene bundle reaches this file through
`expression_profile.datasets[].quantile_shard`, never by constructing the path
itself — the "never construct a path" rule this document opens with. Without
this file, `median_percentile` would be an unauditable number: the whole point
of publishing it is so a consumer can look up the median a gene reported and
confirm which percentile it lands in, the same lookup this atlas's own build
performed once.

## `variants/index.json` and `variants/<chrom>.json.gz`

The committed corpus curates no variants, so the real index is empty:

```json
{ "shards": [] }
```

*Illustrative* — the same file once variants are curated on chromosomes 1 and X:

```json
{ "shards": ["variants/1.json.gz", "variants/X.json.gz"] }
```

- The index exists so a consumer can enumerate chromosomes without probing for
  404s, and is emitted even when empty — as it is in the committed corpus today.
- `shards` is in karyotype order (1…22, X, Y, MT), not lexical, so it can drive
  a chromosome picker directly.
- Each shard is `{"chrom": "12", "rows": [ … ]}` and holds only rows on the
  chromosome it is named for.
- Gzipped: see the decompression note above.

## `publications.json`

```json
{
  "publications": [
    {
      "id": "PMID:8988165",
      "title": "Mutations in human TBX5 [corrected] cause limb and cardiac malformation in Holt-Oram syndrome.",
      "journal": "Nature genetics",
      "year": 1997,
      "authors": ["Basson CT", "Bachinsky DR", "…"],
      "study_type": "family_linkage",
      "doi": "10.1038/ng0197-30",
      "pmcid": null,
      "own_lab": false,
      "cohort_size": null,
      "ancestry": [],
      "tests_reported": null,
      "correction_genes": null,
      "genome_wide_threshold": null,
      "summary_caveat": null
    }
  ]
}
```

`summary_caveat` is one curated sentence (at most 240 characters) saying what
a reader must know to read the study's counts — that PMID:34324492's cases were
drawn from 7,958 individuals, 777 of them thoracic aortic aneurysm probands,
outside this atlas's scope. A gene page's card prints it beneath that study's
results, and the gene index under its strip legend. It is `null` where no curator has
recorded one, which does not mean the study's results need no qualification:
its rows' `method_note` and the [`cohorts.json`](#cohortsjson) descriptions
remain the full caveat. A consumer showing the study's counts or results away
from those should show this sentence too.

`tests_reported` is the number of comparisons a study ran, where it states
one; it is `null` for PMID:40127276 since `2.19`, which recorded 248 until then:
that is the number of its panel *genes*, each tested in two arms over several
consequence classes, not a count of comparisons. That study corrects every row,
so no page needs a scan size for it.

`correction_genes` (since `2.19`) is the number of genes a study's own
multiple-testing correction ran over, recorded for a study whose correction was
over a prespecified gene panel rather than the genome: 248 for PMID:40127276,
whose per-arm Benjamini-Hochberg q reproduces over exactly its 248 panel genes.
`genome_wide_threshold` (since `2.19`) is the p-value the study itself states
for genome-wide significance, transcribed: `2.6e-06` for PMID:40127276, from its
Fig. 1 legend. Both are `null` where not recorded. A row of such a study with
`pvalue_adjusted` below 0.05 survives a correction **within its panel**; whether
its raw `pvalue` also clears the study's genome-wide line is a separate, stronger
claim, and a consumer showing "surviving" should say which.

Ordered by PMID lexically, not numerically — so `PMID:10` precedes `PMID:9`. A
PMID is issued at indexing time, so numeric order ranks by nothing a reader
asked for; pages should rank by year or by the curated featured list.

## `featured.json`

The landing page's manuscripts, in curator-chosen `order`.

```json
{
  "featured": [
    {
      "order": 1,
      "topic": "…",
      "blurb": "One of two back-to-back 1997 reports identifying TBX5 mutations …",
      "publication": { "id": "PMID:8988165", "title": "…", "journal": "Nature genetics", "…": "…" }
    }
  ]
}
```

`publication` is the **resolved object**, not a PMID string, so the landing page
renders without a second fetch. It is never a bare string.

## `phenotypes.json`

```json
{
  "phenotypes": [
    { "id": "HP:0001629", "label": "Ventricular septal defect",
      "lesion_group": "septal", "synonyms": ["VSD"] }
  ]
}
```

`lesion_group` is the facet the gene index's `lesion_groups` and
`confidence_by_lesion_group` key on, which is what lets a phenotype filter drive
a gene list. It is `null` for a term that is not itself a cardiac lesion — an
extracardiac feature (e.g. a limb malformation cited by a syndromic assertion)
registered here only so its label is checked against the pinned HPO release,
the same guarantee every cardiac term gets.

## `datasets.json`

```json
{
  "datasets": [
    {
      "id": "E-MTAB-6814", "archive": "arrayexpress", "design": "profile",
      "technology": "bulk_rnaseq", "organism": "NCBITaxon:9606",
      "tissue": "forebrain, heart, hindbrain, kidney, liver, ovary, testis",
      "developmental_stage": "4 weeks post conception to elderly adult (13 prenatal, 8 post-natal stages)",
      "n_samples": 287, "licence": "EMBL-EBI Terms of Use",
      "publication": "PMID:31243369", "contrasts": [],
      "cardiac_tissues": ["heart"],
      "detection_floor": 1.0, "floor_source": "…", "quantile_estimator": "linear",
      "stages": [
        { "token": "4 week post conception", "wpc": 4.0, "order": 1 },
        { "token": "19 week post conception", "wpc": 19.0, "order": 13 },
        { "token": "neonate", "wpc": null, "order": 14 },
        { "token": "elderly", "wpc": null, "order": 21 }
      ]
    }
  ]
}
```

`stages` and `floor_source` are abridged above — the real record declares all
21 stage tokens and spells the floor's provenance out in full. Every other
value is copied from a build of the committed corpus.

One record per omics dataset, serialised generically from the curated model
(`build/literature.py`'s `_dump`, the same function `publications.json`,
`featured.json` and `phenotypes.json` use) — every field it declares reaches
this file, with nothing filtered out: `id`, `archive`, `technology`, `tissue`,
`developmental_stage`, `organism`, `n_samples`, `licence`, `contrasts`,
`design`, `cardiac_tissues`, `stages`, `cell_types`, `regions`,
`detection_floor`, `floor_source`, `quantile_estimator` and `publication`.
This is what an omics row's `dataset` column resolves against, the way
`publications.json` resolves a PMID.

**The eight from `design` to `quantile_estimator` are what distinguish a
developmental expression dataset, bulk or single-cell, from every other kind
mirrored here, and the two from each other.** `design` is one of
`"contrast"` (a differential-expression comparison — `expression`,
`proteomics`, `phospho`), `"profile"` (a raw abundance series across
development — `profiles`, see
[`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
above) or `"cell_type_profile"` (single-cell abundance by cell type and
region, described below). `cardiac_tissues` and `stages` are `[]` and
`detection_floor`/`floor_source`/`quantile_estimator` are `null` on a
`"contrast"` record — never omitted — so every record has one shape whichever
design it declares. `cardiac_tissues` names which of a profile dataset's own
tissue tokens this atlas reads as "the heart" for that dataset.

`cell_types` (each `token`, `label`, `compartment`) and `regions` are `[]` on
every record except a `"cell_type_profile"`-design one, the single-cell design
that is indexed by cell type and region instead of by stage; the curation model
rejects them on any other design. On that design each region appears once and
is one of `LA`, `RA`, `LV`, `RV` and `IVS`; the pooled `all` is never declared.
No such record is curated yet.

`stages` is one entry per developmental stage the dataset declares. Each
carries the stage's own `token` — the literal string the mirror's rows use,
which is what an `expression_profile` stage joins against — and two positions
on the developmental axis:

- `wpc` is the stage's age in weeks post conception, and is `null` for every
  post-natal stage — 8 of E-MTAB-6814's 21. That null is not a gap in the
  curation; a stage after birth has no such age.
- `order` (added in schema `2.12`) is the **curated chronological position**:
  1-based, ascending, and unique within a dataset. It covers the post-natal
  stages, which is the whole reason it exists — `wpc` is `null` for all of
  them, so those tokens carry no other chronology, and the only remaining way
  to sequence them is the stage *string*, which puts `elderly` before
  `infant`. It is what orders the developmental series published on every gene
  bundle: see
  [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  above.

**This array is sorted by the build, on `order`** — the same sequence the
gene bundle's [`stages`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
array arrives in.

Through schema `2.11` it was not, and the paragraph here said so. Every record
in `datasets.json` is serialised generically from the curated file, so `stages`
arrived in whatever sequence a curator typed and agreed with `order` only
because they typed it in sequence. No validator could have caught a re-typed
block: `PRF011` (two stages claiming one `order`), `PRF012` (`order`
contradicting a stage's own `wpc`) and `PRF014` (a post-natal stage ordered
before a prenatal one) all constrain the *values* of `order`, never the
sequence a file lists them in. Measured — moving the `elderly` block to the top
of the curation file while leaving `order: 21` alone published `elderly` first
with `validate` at its exact documented baseline. The guarantee therefore lives
in the emitter, where it holds by construction rather than by care.

**No published byte changed when it landed**, which is why it carries no
schema letter: the curated file was already in order, so this array is
identical before and after and no consumer can tell the two builds apart. What
changed is that it can no longer stop being in order without the build putting
it back.

**This array's own sequence did move in `2.12`.** The same commit that numbered
the stages re-sequenced the eight post-natal ones in the curated file, so they
arrive `neonate` first rather than `adolescent` first. Every token and every
`wpc` is unchanged — only the order is different, and only for the post-natal
block.

The committed corpus holds one dataset today, the `"profile"`-design
`E-MTAB-6814` shown above. No `"contrast"`-design dataset has been curated
yet, so every `contrasts` array a consumer meets here is empty.

## `cohorts.json`

The sample collections every burden row names, and the caveats that qualify
them. Added in schema `2.9`.

```json
{
  "cohorts": [
    {
      "id": "taa_cases",
      "name": "Sporadic thoracic aortic aneurysm case series (PMID 21092924)",
      "description": "777 individuals with sporadic thoracic aortic aneurysm, included in the CNV case set of PMID:34324492 alongside the CHD cases; the study's abstract describes its probands as \"probands with CHD including cases with sporadic thoracic aortic aneurysm\". Thoracic aortic aneurysm is outside this atlas's scope …",
      "url": null,
      "overlaps_with": [],
      "review_status": null,
      "review_note": null
    }
  ]
}
```

**This file is the reason `case_cohorts` and `control_cohorts` are worth
reading.** A burden row names its collections by bare id; without this table
those ids resolve to nothing, and a consumer computing anything from the row is
working with numbers stripped of the sentences that say what they count.

- `id` is exactly the string a burden row's `case_cohorts` / `control_cohorts`
  array carries. Every id in every published row appears here — the build
  publishes the whole curated registry rather than the subset the published rows
  happen to cite, and a validator (`BUR009`) refuses a row naming a collection
  this file does not carry, so resolution cannot fail.
- `name` is the collection's full name, always a non-empty string.
- `description` carries the caveats, and it is the field this file exists for.
  Examples from the current registry: `ukbb`'s controls are adults recruited at
  40–69 while a substantial proportion of the cases were enrolled in childhood;
  `ddd` ascertains on developmental disorder rather than heart disease, so its
  contribution is
  enriched for **syndromic** CHD; `gnomad_controls` is not screened for
  congenital heart disease. None of these is expressible as a column.
- `url` is the collection's public page, or **`null`** where it has none (11
  of the 20 today). Always present, never omitted.
- `overlaps_with` (since `2.19`) lists, sorted, the collections this one is
  declared to share people with although no row names both — `ddd` and
  `decipher` list each other — and is `[]` elsewhere. Declared on both entries,
  checked when the registry loads. `cohort_families` merges studies across these
  as across a shared id, so it is part of what makes two studies one dataset.
- `review_status` and `review_note` (since `2.19`) record a question about
  the record that this atlas leaves open rather than settles. The rule, the
  atlas owner's: such a record says so rather than choosing a side,
  speculating, or being dropped. `review_status` is `"pending_review"`, a
  question about the record left open for an expert, or `"not_accessible"`, a
  source the record rests on could not be read; `review_note` states exactly
  what is unresolved and nothing more. The vocabulary is
  [`scope.json`](#scopejson)'s too, where the question can be a clinical
  judgement rather than a conflict between sources. Both are set or both are
  `null`, never omitted, and `null` means no question has been recorded, not
  that the record was checked and found clean. Two records carry one today,
  both `"pending_review"`: `pmid22939634_cases` and `pmid22969434_cases`, on
  each of which PMID:34324492's S2 Table and the cited paper's abstract
  conflict. The note renders after the description on
  [`methods.html#cohorts`](#methodshtml) and in a gene page's cohort notes.

The array is sorted by `id`. There is no `role` field: a collection is cases in
one study and could be controls in another, so which it was is a property of the
burden row that cites it — the column it appears in — not of the collection.

**A row may not name every collection its study drew on.** PMID:34324492's S2
Table lists case sources that no row names;
[`methods.html#cohorts`](#methodshtml) states how many of its sources the rows
do name, and how many of its cases those hold. So two studies that name no
collection in common, and none this file records as overlapping, may still
share people, and `cohort_families` cannot see it. The seven
`CHD case series of PMID:…` records, added in `2.19`, are the CHD case sources
that table identifies by PMID with 50 or more cases (the one other source it
identifies by PMID with that many is the study's aneurysm series,
`taa_cases`); each description says what the table says of its series and
nothing more.

`counts.cohorts` in `manifest.json` is the length of this array. It is **not**
`counts.cohort_families`, which counts how many independent collections the
rare-variant evidence groups into (2 today, 3 before `2.19` declared DDD and
DECIPHER overlapping) and is what
[`independent_datasets`](#independent_datasets-a-count-of-datasets-never-a-verdict)
is derived from.

## `scope.json`

The registry of disease terms this atlas treats as congenital heart disease,
and for each one the external authority it is attributed to. A bundle names
every validity record's disease by bare MONDO id; this is the table those ids
resolve against, in the same way [`cohorts.json`](#cohortsjson) resolves a
burden row's cohort ids.

**105 terms**, sorted by `id`. 69 are attributed to an authority that curates
them (a ClinGen expert panel, or a named GenCC submitter); **36** are admitted
on Orphanet's phenotype frequencies (`admitted_by:
"orphanet_phenotype_frequency"`), the path added in `2.15` for syndromes where
structural CHD is frequent but not part of the diagnosis.

```json
{
  "id": "MONDO:0012547",
  "label": "Noonan syndrome 4",
  "admitted_by": "orphanet_phenotype_frequency",
  "attributed_to": "Orphanet",
  "orphanet": {
    "code": "ORPHA:648",
    "name": "Noonan syndrome",
    "via": "subtype",
    "omim": "OMIM:610733",
    "phenotypes": [
      {
        "hpo": "HP:0001641",
        "label": "Abnormal pulmonary valve morphology",
        "frequency": "Frequent (79-30%)"
      }
    ]
  },
  "syndrome": false,
  "review_status": null,
  "review_note": null
}
```

- `syndrome` (since `2.19`) is a curated mark, always a boolean: `true` where
  the term names a syndrome — a disease with a heart malformation among its
  features — rather than heart disease itself, so a grade for it is a grade
  for the syndrome as a whole. Set on the nine ClinGen-scoped syndrome terms
  that published genes are graded for at their headline (Holt-Oram, CHARGE,
  Alagille, Ellis-van Creveld, Kabuki 1 and 2, Mowat-Wilson, TARP, and
  NR2F2-related multiple congenital anomalies/dysmorphic syndrome). `false`
  says no curator has marked the term, which is not a finding that it is not
  a syndrome: the Orphanet-path terms are syndromes by their warrant and are
  not marked, and "congenital heart disease with heterotaxy syndrome" is CHD
  itself whatever its name says.
- `review_status` and `review_note` (since `2.19`) record an open question
  about the term, in the vocabulary [`cohorts.json`](#cohortsjson) uses: both
  set or both `null`, never omitted. Two terms carry `"pending_review"` today,
  MONDO:1010098 and MONDO:1060150, marking a question left open for an expert
  rather than a conflict between sources: whether each names a syndrome,
  which decides whether a gene page calls a grade for it a grade for a
  syndrome as a whole, is left open for clinical review, and `syndrome` stays
  `false` meanwhile. Published here only: no page renders a term's question.
- `admitted_by` is a string token naming the kind of authority that puts the
  term in scope, one of four:

  | `admitted_by` | the term is in scope because | `attributed_to` names |
  | --- | --- | --- |
  | `clingen_chd_panel` | ClinGen's Congenital Heart Disease Gene Curation Expert Panel uses it | that panel |
  | `clingen_expert_panel` | another ClinGen expert panel uses it | that panel |
  | `gencc_submitter` | a body submitting to GenCC uses it | the submitter, never GenCC itself |
  | `orphanet_phenotype_frequency` | a curator nominated the syndrome and Orphanet's phenotype frequencies warrant it | `Orphanet` |

  Read the token, not the name: Orphanet also submits to GenCC, so
  `attributed_to: "Orphanet"` appears on `gencc_submitter` terms too.
- **This `admitted_by` is not a gene bundle's `admitted_by`, though the two share
  a name.** The bundle's is an object — `authority` is `"clingen"` or
  `"gencc_agreement"` — saying which warrant admitted a *gene*; this one is a
  string saying which kind of authority put a *disease* in scope. They answer
  different questions and do not join: PTPN11's bundle says `"clingen"`, for a
  disease this file says is in scope on `"orphanet_phenotype_frequency"`.
- `orphanet` is `null` on every term an authority curated, never absent.
- For an Orphanet-path term it is the warrant. It names the Orphanet disorder,
  how the term reaches it (`via: "direct"`, a MONDO exact match, or
  `"subtype"`, with `omim` naming the OMIM entry Orphanet maps under its code
  by a *validated* exact-or-broader mapping; `omim` is `null` on every
  `"direct"` term), and every structural heart phenotype Orphanet records for
  it at `Frequent (79-30%)` or stronger.
  Frequencies are Orphanet's band names verbatim.
- **The warrant is re-derived from the mirrors on every build**, by the same
  rule the validator `SCP006` refuses a nomination on, so it cannot drift from
  what Orphanet records. "Structural" is HPO's `Abnormal heart morphology`
  branch minus its cardiomyopathy, cardiomegaly, tumour, pericardium,
  endocardium and imaging-finding branches, plus eight lesions HPO files
  elsewhere, with any subtypes: coarctation of the aorta, patent ductus
  arteriosus, pulmonic stenosis, aortic valve stenosis, supravalvular aortic
  stenosis, truncus arteriosus, interrupted aortic arch and anomalous pulmonary
  venous return, over the pinned HPO release. The last six were added in `2.19`,
  which added pulmonic stenosis to the published phenotypes of ten terms (the
  cardiofaciocutaneous, Costello and Noonan-with-multiple-lentigines families)
  and admitted no term by itself: a curator still nominates.
- **A grade for one of these syndromes is a grade for the syndrome, not for
  heart disease.** A client showing PTPN11's `Definitive` for Noonan syndrome
  should say why Noonan syndrome counts here, which is what this object is
  for. The site's gene pages do so at the top of the validity section, above both
  groups of classifications.

## `sources.json`

What the atlas mirrors, and on whose terms.

```json
{
  "sources": [
    {
      "id": "hpo",
      "name": "Human Phenotype Ontology",
      "version": "hp/releases/2026-06-23",
      "retrieved_on": "2026-07-31",
      "url": "https://hpo.jax.org/",
      "licence": "https://hpo.jax.org/app/license",
      "redistribution": "permitted_with_attribution",
      "ontology_prefix": "HP",
      "ontology_file": "ontologies/hp-2026-06-23.obo",
      "file_version": null,
      "file_date": null
    }
  ]
}
```

**Read this before redistributing anything from this site.** Phenotype labels
and synonyms in `phenotypes.json` and in the search index are transcribed from
the pinned HPO release, whose terms are `permitted_with_attribution` — so a
consumer republishing them carries the same obligation, and this file is where
the attribution to satisfy it comes from.

The repository's `LICENSE` (Apache-2.0) covers the **code**. It does not govern
mirrored third-party content, whose terms are the ones recorded here.

`version` is the upstream release identifier and `retrieved_on` the date it was
taken, so a claim can be traced to the exact release it rests on.

`file_version` and `file_date` (since `2.19`) are what a mirrored file states
about itself, verbatim, where `version` is a checksum that identifies the file
but cites as nothing. Both Orphadata products carry them, from their `<JDBOR>`
header — `"1.3.42 / 4.1.8 [2025-03-03]"`, dated `2026-06-23` — written beside
the sha256 by the converter on every run; `null` elsewhere. Orphadata's
requested citation ends "Data version [date]", and the site's front page prints
both there.

The registry also carries `clingen` and `gencc` — the two sources every
mirrored record in a gene bundle's `validity.records` is attributed to.
Both are recorded `"redistribution": "permitted"` (CC0-1.0), unlike HPO's
`"permitted_with_attribution"`: neither authority's licence requires an
attribution notice to redistribute their content, though GenCC's terms
*request* one. This file is where that distinction is recorded, rather than
assumed from the presence of a `validity` object.

**`scope.json` republishes Orphanet content, and that content requires
attribution.** Each Orphanet-path term's `orphanet` object carries Orphanet's
disorder code and name, its OMIM reference and its phenotype frequency bands,
from the entries `orphanet_product1` and `orphanet_product4` record here as
`CC-BY-4.0`, `"permitted_with_attribution"`; the MONDO-to-Orphanet match behind
`via` is read from `mondo_sssom`, recorded the same way, and the phenotype
labels are HPO term names as Orphanet records them. A consumer republishing
that object carries the attribution obligation, as for HPO. Other entries carry their own terms — every source's
`redistribution` is in this file, and it, not this paragraph, is the list.

## `search/index.json.gz`

A flat array of records over genes, publications and phenotypes. Deliberately
not an inverted index: at this corpus size a client filters the whole array in a
fraction of a frame.

```json
{
  "records": [
    { "kind": "gene", "id": "HGNC:11604", "label": "TBX5",
      "path": "genes/HGNC_11604.json",
      "terms": ["TBX5", "HGNC:11604", "T-box transcription factor 5", "HOS"] },
    { "kind": "phenotype", "id": "HP:0001631", "label": "Atrial septal defect",
      "path": "phenotypes.json",
      "terms": ["Atrial septal defect", "HP:0001631", "ASD"] }
  ]
}
```

- `terms` is the haystack: the strings a visitor might type, deduplicated,
  including each record's own identifier. **Matching is the client's job** —
  this file ships no scoring, no stemming and no ranking.

  For a gene the terms are its approved symbol, its HGNC id, its name, its
  current aliases and — since `2.13` — the symbols HGNC has **withdrawn** from
  it, in that order. A withdrawn symbol is what a paper published ten years ago
  called the gene, and 48 of the 119 indexed genes have one: `HOS` finds TBX5,
  `VEGF` finds VEGFA, `HTX1` finds ZIC3.

  **A term can name more than one record, and that is not a defect to
  de-duplicate.** It could already happen through a shared alias; `2.13` adds
  one live instance — `MADH7`, withdrawn by both SMAD6 and SMAD7, so it is a
  term on both. Present every hit rather than picking one: this atlas refuses to
  guess which gene an ambiguous symbol means (see `mirrors/genes.tsv`), and a
  client that silently kept the first would make that choice on its behalf, with
  a 50% chance of showing the reader the wrong gene's evidence.
- `label` is what a result row displays; `id` identifies the thing.
- `path` is the payload that answers the query.
- Variants and datasets are not indexed. The variant space grows without bound
  and would dominate the size of the file every visitor downloads; a dataset has
  no title or description of its own, so its accession is the only string naming
  it, and that already resolves through `datasets.json`.
- Genes come from the **published** set — the same 119 genes `genes/index.json`
  lists and `genes/<slug>.json` serves, each one admitted by an external
  authority under the rule above. Deliberately not the assertion set: keyed on
  that, the index would hold one record while the site published 119 bundles, so
  a visitor typing GATA4 would get nothing while `genes/HGNC_4173.json` was
  being served. Deliberately not the 185-gene mirror registry either, which
  would offer 66 results whose `path` names a bundle no builder wrote. A gene outside the
  published set is not searchable because it has no page to open.
- **A search record carries no curation or confidence signal.** `kind`, `id`,
  `label`, `path` and `terms` are the whole record — there is no
  `headline_confidence`, no `validity_state`, no `atlas_curation` and no
  `has_conflicting_evidence` in it. The obligations this document places on a
  consumer — read `atlas_curation` before presenting a row as curated content,
  and pair `headline_confidence` with `has_conflicting_evidence` — therefore
  **cannot be met from this file alone**. A typeahead built on it must join each
  hit back to `genes/index.json`, or fetch the bundle at `path`, before it
  labels a result as anything beyond a name that matched. This is not a
  hypothetical gap: of the 119 genes indexed today, 118 are `"not_yet_curated"`
  (measured 2026-10-25), and 21 carry no expert-panel grade at all, so a UI that
  presents a search hit as atlas-curated content would be wrong about nearly
  every one of them.

---

## Contested genes: the one consumer obligation

`headline_confidence` is a ClinGen expert panel's grade: the strongest of the
gene's in-scope ClinGen records — for a gene ClinGen admitted, the grade of the
record that admitted it — or `null` where no ClinGen panel has graded the gene
in scope, which today is exactly the genes admitted on submitter agreement. It
is never a maximum over GenCC submitters: TBX1's submitters grade it
`Definitive` and `Strong`, and its headline is `null`. Among ClinGen's records
the strongest is read off a single linear scale where `definitive` outranks
`refuted`, and ClinGen treats disputed and refuted as a **separate axis** rather
than weaker rungs of the same ladder, so a gene carrying a ClinGen `Definitive`
record and a `Disputed` or `Refuted` one resolves to `definitive`, and the
dispute is invisible in that field alone. A GenCC submitter's dispute is
invisible in it for a plainer reason: no GenCC record reaches the headline.

**`headline_confidence` is `null` for a gene no ClinGen expert panel has graded
in scope, and must never be rendered as `"no_known_association"`.** The two are
not interchangeable: `no_known_association` is itself an assessed verdict — a
panel looked and found nothing. `null` is not that claim, and it is not "nobody
looked" either. It is what every gene admitted on submitter agreement
publishes, however highly its submitters graded it (their grades are in
`validity.records`, and `admitted_by.submitters` names them); and a gene whose
only in-scope record is a GenCC submitter's `Supportive`, which asserts an
association without grading it, would publish it too — six genes in the
committed mirrors carry exactly one such Orphanet record, and none is published
(see `validity_state` above). Read `validity_state` to tell "no panel graded
it" from "nobody has assessed it". Coercing `null` to `"no_known_association"`
would fabricate a conclusion nobody reached. This atlas publishes no
gene-disease validity classification of its own; every value `headline_confidence` can take, including the absence
of one, is mirrored and attributed from ClinGen or GenCC — see the bundle's
`validity` object, documented under [`genes/<slug>.json`](#genesslugjson), for
the records behind it and `sources.json` for the licence terms those two
mirrors are republished under.

`has_conflicting_evidence` is the other half of that pair. It appears in both
the browse row and the bundle, and is always written alongside
`headline_confidence`. `has_source_discordance` is a narrower relative: it is
`true` only when the contesting and the supporting classification come from
*different* mirrored sources. A single source split against itself would set
`has_conflicting_evidence` without setting this one — that is why the two
fields are not redundant — but that split does not currently occur among the
185 genes these mirrors curate within CHD scope (measured 2026-10-26: the 8
genes that set one flag set both).

**Every published gene that sets both flags today is in the table below — and
none of them is contested for the disease it is graded on.** Measured against a
real build (2026-10-26): of the 185 mirrored genes, 8 set
`has_conflicting_evidence`, and the same 8 set `has_source_discordance`. 6 are
published. ClinGen's RASopathy Gene Curation Expert Panel grades each
`Definitive` for the syndrome it was admitted on and marks it `Disputed` for a
sibling syndrome — a panel saying which gene belongs to which syndrome, not that
the gene's role in CHD is in doubt:

| gene | `Definitive` for | `Disputed` for |
| --- | --- | --- |
| BRAF | cardiofaciocutaneous syndrome | Costello syndrome |
| KRAS | Noonan syndrome | Costello syndrome |
| MAP2K1 | cardiofaciocutaneous syndrome | Costello syndrome |
| PTPN11 | Noonan syndrome; Noonan syndrome with multiple lentigines | Costello syndrome; cardiofaciocutaneous syndrome |
| RAF1 | Noonan syndrome | Costello syndrome; cardiofaciocutaneous syndrome |
| SOS1 | Noonan syndrome | Costello syndrome; cardiofaciocutaneous syndrome |

The `Definitive` column is every in-scope ClinGen record at the gene's headline
grade, and includes the disease `admitted_by` names. Costello and
cardiofaciocutaneous syndrome entered scope with `2.15`, and that is what made
these the first contested genes this API has published. The flags are computed
over a gene's in-scope records, not per disease, so they cannot tell "disputed
for a sibling syndrome" from "disputed for the disease it is graded on"; read
the `validity.records` marked `Disputed` to see which.

The other two are not published. **LEFTY2** (HGNC:3122) has two GenCC
submitters at `Limited` and would qualify on submitter agreement, but ClinGen's
Congenital Heart Disease GCEP disputes it, and the gate refuses it explicitly: a
source trusted to admit a gene alone must be trusted to refuse one. **A2ML1**
(HGNC:23336) is refused the same way: Genomics England PanelApp and PanelApp
Australia agree on it, and ClinGen's RASopathy panel disputes it for Noonan
syndrome itself. Without that veto both would be
published, with the dispute buried in a payload no consumer had reason to
inspect.

So a consumer testing their rendering against the live data now has contested
genes to look at, and **a client that pairs nothing with `headline_confidence`
shows every one of them as an ordinary, settled row** — the failure this
section exists to prevent.

**A consumer must pair `headline_confidence` with `has_conflicting_evidence`
and present a contested gene distinctly** — a badge, a different colour, an
explicit note. Rendering `headline_confidence` alone would tell a reader that a
gene the field disputes is settled science, which is the one failure this
resource exists to prevent.

`confidence_by_lesion_group` is **not** a finer-grained view of the same
question. ClinGen and GenCC classify a gene against a disease, never against a
specific lesion, so the mirrored records carry no per-group information at
all — every lesion group a curated assertion names for the gene publishes the
*identical* value `headline_confidence` does. It differs from
`headline_confidence` only in shape, as a map over the gene's declared groups
for a consumer already filtering by lesion, never in value.

**`conflicting_lesion_groups` is `has_conflicting_evidence`'s per-group
counterpart**, and for the same reason it cannot single out which of a
contested gene's groups is the disputed one: it lists *every* group in
`confidence_by_lesion_group` when the gene is contested, and none when it is
not. *Illustrative* — no published gene today is both contested and curated
for a lesion group (the contested genes above have no curated assertion, so
both fields are empty for them), so this row is constructed:

```json
{
  "has_conflicting_evidence": true,
  "confidence_by_lesion_group": { "conotruncal": "definitive", "septal": "definitive" },
  "conflicting_lesion_groups": ["conotruncal", "septal"]
}
```

Read together, those say: the gene is contested, and every lesion group it is
curated for inherits that contest equally. There is no mirrored signal that
could clear one group while leaving another disputed.

**Where these two fields live.** `confidence_by_lesion_group` and
`conflicting_lesion_groups` appear in `genes/index.json` and **nowhere else**.
The gene bundle carries `headline_confidence`, `validity_state`,
`has_conflicting_evidence` and `has_source_discordance`, but neither of these
two, so a detail page needing group-level confidence must carry it over from
the browse row it was opened from. It cannot be recovered from the bundle:
that would mean reimplementing the classification ranking and the contested
test against the mirrored records, and neither rule is published.

**The two states the pair can express:**

| `has_conflicting_evidence` | `conflicting_lesion_groups` | what it means |
| --- | --- | --- |
| `false` | `[]` | nothing about this gene is disputed |
| `true` | every group in `confidence_by_lesion_group` | the gene is disputed, and the dispute applies equally to every lesion group it is curated for |

A gene disputed about only *some* of its lesion groups is not a state this API
can express: the mirrors classify by disease, not by lesion, so there is
nothing in the source data to divide the groups on. `conflicting_lesion_groups`
is therefore always either every group the gene declares or none of them —
never a proper subset.

The list is always present and may be empty. Every group it names is a key of
`confidence_by_lesion_group`, so the two join directly.

### The same failure in a second shape: a gene the atlas has not curated

`headline_confidence` is published for every gene, curated here or not, so a
consumer that renders it alone presents 96 of the 119 published genes — every
gene this atlas has not curated and a ClinGen panel has graded — as this
atlas's assessment of them. It is not: it is a ClinGen expert panel's, and this
atlas has recorded nothing about those genes beyond republishing it with its
provenance. On the other 21 genes not yet curated here it is `null` — no panel
graded the gene, and `admitted_by` names the GenCC submitters whose agreement
admitted it instead.

**A consumer must not present a `not_yet_curated` gene as one the atlas has
assessed.** `atlas_curation` appears on the browse row and in the bundle,
written in one place so the two cannot disagree, and it is what tells the two
kinds of row apart — read it rather than testing `assertion_count == 0`, so a
filter does not have to reimplement the rule. An empty `assertions` array is a
curation gap, never a verdict and never a fetch that failed, the same way
`headline_confidence: null` is not `no_known_association`. This atlas's own
gene pages state the gap in prose for exactly that reason.
