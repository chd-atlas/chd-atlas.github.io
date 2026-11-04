# Prose check: docs/data-api.md, the built pages and the curator docs

- Commit checked: `508e172` (origin/dev, "Merge pull request #89 from chd-atlas/feat/literature-timeline-graphic")
- Build: `uv sync --locked --extra dev && uv run chd-atlas build --root . --out /tmp/dist` → `wrote 254 file(s) to /tmp/dist`
- `schema_version`: `2.21` (`status` `in-development`, `source_commit` `508e17288ba736378bbe7b0e6bf690b528f300b7`)
- `counts` (from /tmp/dist/manifest.json):

```json
{"assertions": 2, "burden_combined": 108, "burden_rows": 1210, "burden_surviving_genes": 32,
 "cell_type_datasets": 1, "cell_type_genes": 118, "cohort_families": 2, "cohorts": 20,
 "datasets": 2, "featured": 1, "functional": 0, "genes": 119, "literature_findings": 0,
 "literature_genes": 0, "phenotypes": 8, "profile_datasets": 1, "profile_genes": 119,
 "publications": 11, "scope_terms": 105}
```

This report only describes problems. It changes no file other than itself.

## Method

Every number was re-derived from /tmp/dist with a command; none was accepted because another
sentence agrees with it. Page text was read from the built HTML, never from the Python source.
The source line given for a published sentence was found afterwards, to say where to fix it.

Many findings below concern sentences that are in a page's markup but not shown to a reader.
Grepping the HTML cannot tell these apart from shown text, so pages were rendered in the
Chromium that ships with Playwright. Commit a89e958 (2026-11-01, the UI redesign) introduced
all of the hiding. It changed `build/landing.py`, `build/pages.py` and `build/render.py`, but
not `docs/data-api.md`. The visibility check that the hidden-text findings cite as **[V]** is:

```bash
cd /tmp && node -e "const {chromium}=require('playwright');(async()=>{const b=await chromium.launch();const p=await b.newPage();for(const f of process.argv.slice(1)){await p.goto('file:///tmp/dist/'+f);const r=await p.evaluate(()=>[...document.querySelectorAll('.notice,.scope-rule,[style*=\"display: none\"]')].map(e=>getComputedStyle(e).display+' | '+e.textContent.replace(/\s+/g,' ').trim().slice(0,70)));console.log(f);r.forEach(x=>console.log('  '+x))}await b.close()})()" index.html genes/index.html genes/HGNC_11604.html methods.html
```

Output (each `none` line is an element whose computed `display` is `none`):

```
index.html
  none | This atlas is under active development. It is not a clinical decision-
  none | 2 curated gene-disease assertions so far: TBX5 (HGNC:11604), NR2F2 (HG
  none | 119
genes/index.html
  none | This atlas is under active development. It is not a clinical decision-
  none | Published on an external authority's classification for the disease it
  block | Published on an external authority's classification for the disease it   <- the <p> inside the hidden div above
  none | 119
  none | enriched, and a result survives its own study's correction enriched no
genes/HGNC_11604.html
  none | This atlas is under active development. It is not a clinical decision-
  none | Published on an external authority's classification for the disease it
methods.html
  none | This atlas is under active development. It is not a clinical decision-
  block | A gene is published here when an external authority already treats it
```

## Known (not re-investigated)

- methods.html "A curated gene's page also has a Publications section" and docs/data-api.md:377
  "what a curated gene's page adds" use "curated gene" ambiguously, because the card also counts
  a gene whose only curated work is literature. Note: at 508e172 the built methods.html reads
  "A gene with a lesion assertion also has a Publications section"; the phrase "curated gene"
  does not appear in it. The data-api.md wording remains, and docs/data-api.md:1781 ("A curated
  gene's page follows it with the curated evidence") repeats it.
- The module docstring of src/chd_atlas/build/gene_literature.py says a Publication is withheld
  "until a finding of that paper is approved" (imprecise; fix scheduled). Finding B4 is the same
  imprecision in a curator skill.
- dev has no literature data: `literature_findings` and `literature_genes` are 0, and every gene
  page says "Not yet curated." in its Literature section. Expected, not a finding.

## Summary

| Verdict | Published HTML and API doc | Curator docs | Total |
|---|---|---|---|
| FALSE | 12 | 1 | 13 |
| MISLEADING | 18 | 3 | 21 |
| STALE NUMBER | 5 | 0 | 5 |
| UNVERIFIABLE | 4 | 0 | 4 |
| **Total** | **39** | **4** | **43** |

---

## A. Published HTML and docs/data-api.md (most harmful first)

### A1. The "not a clinical decision-support tool" notice is in every page but shown on none
- Where: src/chd_atlas/build/render.py:144-146 (`.notice { display: none; }`) hides the notice that
  render.py:860 (`RESEARCH_USE_NOTICE`) emits. Built: all 122 pages (index.html, methods.html,
  genes/index.html, genes/HGNC_*.html). Described as shown at docs/data-api.md:308-310 and
  :315-318, and in render.py:37-45.
- Sentence (page): "This atlas is under active development. It is not a clinical decision-support
  tool and must not be used to make or guide a diagnostic, treatment or any other clinical
  decision."
  Sentence (doc :308-310): "This is a research resource, not a clinical decision-support tool;
  see `index.html` at the site root for the statement in full"
  Sentence (doc :315-318): "It carries the same development-status and research-use statement
  `status` above is the machine-readable half of, as every page does"
- Measured: [V] gives `none | This atlas is under active development…` on all four pages
  rendered. Over all 122 pages:
  `node -e "const {chromium}=require('playwright');(async()=>{const b=await chromium.launch();const p=await b.newPage();const fs=require('fs');const files=['/tmp/dist/index.html','/tmp/dist/methods.html','/tmp/dist/genes/index.html',...fs.readdirSync('/tmp/dist/genes').filter(f=>/^HGNC_.*html$/.test(f)).map(f=>'/tmp/dist/genes/'+f)];let h=0;for(const f of files){await p.goto('file://'+f);if(await p.evaluate(()=>getComputedStyle(document.querySelector('div.notice')).display)==='none')h++}console.log(files.length,h);await b.close()})()"`
  → `122 122`. No rendered page has "clinical decision" or "research use" anywhere in
  `document.body.innerText` (checked on index.html, methods.html, genes/index.html and
  genes/HGNC_11604.html). `git blame -L 144,146 src/chd_atlas/build/render.py` → line 145 is
  from a89e958.
- Verdict: MISLEADING. The sentence is present, and the doc truthfully says every page carries
  it, but no reader sees it. A gene page with a green `definitive` chip shows no warning against
  clinical use, and the doc sends readers to index.html for a statement that page does not show.
- Remedy (suggestion): make `.notice` visible again, and add a built-page test that it renders
  visible, not just present.

### A2. The visible glyph legend says "Multi-cohort replicated" for a single cohort family's result
- Where: src/chd_atlas/build/landing.py:374-380 and src/chd_atlas/build/pages.py:7196-7202.
  Built: index.html and genes/index.html, under the gene table (visible).
- Sentence: "Replication across independent cohorts: [filled] Multi-cohort replicated [half]
  Nominal [hollow] No enrichment – Not tested"
- Measured: each glyph is one cohort family. A filled glyph means that family's result survived
  its own study's correction; the tooltip reads, for example, "PMID:40127276: enriched, and
  survives this study's own correction".
  `python3 -c "import json,collections;g=json.load(open('/tmp/dist/genes/index.json'))['genes'];print(collections.Counter(x['independent_datasets']['corrected'] for x in g));print([x['symbol'] for x in g if x['independent_datasets']['corrected']==2])"`
  → `Counter({0: 87, 1: 29, 2: 3})`, `['TBX1', 'NKX2-5', 'GATA4']`. So 29 of the 32 genes that
  show a filled glyph have exactly one corrected family, yet the label beside it says
  "Multi-cohort replicated". The page's own correct legend ("enriched, and a result survives its
  own study's correction", pages.py:7194 `_STRIP_LEGEND`) is in the hidden div (see A3). The label
  also contradicts methods.html#rare-variants, which says the strip "is a count of evidence —
  never a verdict on the gene".
- Verdict: FALSE
- Remedy (suggestion): use the strip legend's wording for the key and drop "Replication across
  independent cohorts".

### A3. The browse page hides its caveats next to 71 "no enrichment" glyphs
- Where: src/chd_atlas/build/pages.py:7181-7189 (intro in `<div style="display: none;">`) and
  pages.py:7194-7195 (strip legend and per-study notes in `<div style="display: none;">`). Built:
  genes/index.html; the same table on index.html has no caveat either.
- Sentences (all hidden): "One glyph per independent cohort family, not per study: …"; "No
  enrichment at these cohort sizes is not evidence against a gene — rare-variant tests routinely
  detect nothing for genes with overwhelming family and functional evidence."; "Neither is a
  verdict on the gene."; "Audain E et al. 2021 (PMID:34324492): Its cases were drawn from 7,958
  individuals, 777 of them sporadic thoracic aortic aneurysm probands, outside this atlas's
  scope; …"
- Measured: [V] gives `none | enriched, and a result survives its own study's correction…` for
  the legend-and-caveat div on genes/index.html. The visible strip shows these glyph states:
  `python3 -c "import json,collections;g=json.load(open('/tmp/dist/genes/index.json'))['genes'];print(collections.Counter(f['state'] for x in g for f in x['independent_datasets']['families']))"`
  → `not_tested 84, no_enrichment 71, nominal 48, corrected 35`. src/chd_atlas/build/concordance.py:207-209
  says why the caption exists: "the page carries a fixed caption saying that no enrichment here is
  not evidence against a gene". It is not shown.
- Verdict: MISLEADING. The sentences are true but hidden, and the visible key (A2) contradicts
  the hidden one.
- Remedy (suggestion): show the legend-and-caveat block, or delete it deliberately. Do not leave
  reviewed caveats hidden by an inline style.

### A4. index.html hides its panel-correction caveat and own-lab count; the API doc says it states them
- Where: docs/data-api.md:320-323 and :337-342. Built: index.html. The hidden sentences come from
  src/chd_atlas/build/landing.py:408-412 (`<div style="display: none;">`).
- Sentence (doc :320-323): "It states seven keys of [`counts`](#manifestjson): … and
  `assertions` in the sentence naming the curated genes."
  Sentence (doc :337-342): "Below the cards it says, for any published study whose correction
  ran over a gene panel, that the correction was not genome-wide, and how many of the rare-variant
  studies whose results the site publishes are by an author of this atlas … and that a gene page
  says so wherever it shows such a study's results."
- Measured: [V] gives `none | 2 curated gene-disease assertions so far: TBX5 (HGNC:11604), NR2F2…`
  on index.html. That hidden div holds all three sentences: "2 curated gene-disease assertions so
  far: TBX5 (HGNC:11604), NR2F2 (HGNC:7976).", "For PMID:40127276, the correction a result must
  survive ran over its 248-gene panel, not the genome." and "Of the 3 rare-variant studies whose
  results appear here, 2 are by an author of this atlas, …". The six other counts the doc names
  (genes 119, cohort_families 2, burden_rows 1,210, burden_surviving_genes 32, profile_genes 119,
  cell_type_genes 118) are visible and match the manifest.
- Verdict: FALSE as rendered. The page visibly states six keys, not seven, and does not show the
  panel-correction or own-lab caveat. Of the three, the panel caveat matters most to a clinical
  reader: it qualifies "survives its own correction" for PMID:40127276.
- Remedy (suggestion): unhide the three sentences, the 248-gene caveat first, or rewrite both
  doc paragraphs to describe what is shown.

### A5. The NR2F2 assertion's recorded curator is "LLM CHD Atlas", but methods.html describes every assertion as a curator's grading
- Where: src/chd_atlas/build/methods.py:288-295. Built: methods.html#curated-evidence, and
  genes/HGNC_7976.html ("Each lesion assertion below is this atlas's own curation").
- Sentence: "Each entry is a lesion assertion: a curator's record that the gene belongs to one or
  more congenital heart disease lesion groups, … its strength as the curator grades it (strong,
  moderate, supporting or contradictory) and the curator's summary of what the publication shows."
- Measured: `python3 -c "import json;[print(a['id'],a['curator'],a['curated_on']) for f in ['HGNC_11604','HGNC_7976'] for a in json.load(open('/tmp/dist/genes/'+f+'.json'))['assertions']]"`
  → `CHDA:AST:0000001 Enrique Audain 2026-07-31`, `CHDA:AST:0000002 LLM CHD Atlas 2026-11-01`.
  `grep -c LLM /tmp/dist/genes/HGNC_7976.html` → `0`: the page never shows who curated it.
  Neither the page nor methods.html says an assertion can be machine-authored, or whether a named
  person reviewed it. For literature findings, by contrast, methods.html says a named curator
  approves each finding before publication.
- Verdict: MISLEADING. A clinical reader takes "a curator's record … as the curator grades it"
  to mean a person. Whether an LLM counts as "a curator" is the owner's call; the page should not
  leave the reader to assume.
- Remedy (suggestion): show each assertion's curator on the page, and say in methods.html whether
  an LLM-authored assertion is reviewed by a named person before publication.

### A6. "1,210 rare-variant statistics" counts rows, not the statistics the linked definition describes
- Where: src/chd_atlas/build/landing.py:395 (figure label), with methods.py:348-350 (the
  definition). Built: the index.html hero, and methods.html#rare-variants, which the hero links as
  "How to read these figures".
- Sentence (index.html): "1,210 rare-variant statistics"
  Sentence (methods.html): "A rare-variant statistic is one study's rare-variant test of one gene
  under one design: a case-control burden test of rare SNVs and indels or of copy-number
  variants, or a de novo enrichment test"
- Measured: `python3 -c "import json,glob;r=[(json.load(open(f))['gene'],b) for f in glob.glob('/tmp/dist/genes/HGNC_*.json') for b in json.load(open(f))['burden']];print(len(r),len({(g,b['study'],b['variant_class'],b['comparator']) for g,b in r}),sum(b['consequence_class']=='synonymous' for g,b in r))"`
  → `1210 262 339`. Under the stated definition there are 262 statistics. The 1,210 is
  `counts.burden_rows`, which includes 339 synonymous negative-control rows and the union rows
  (the 72 `damaging` rows) that methods.html itself calls "not a third independent result".
- Verdict: MISLEADING
- Remedy (suggestion): label the figure as published result rows (including negative controls),
  or make the definition count what the figure counts.

### A7. The scope-rule line is hidden on the browse page and on all 119 gene pages; the API doc says they show it
- Where: docs/data-api.md:368-369, :682-685 and :1764-1766. Built: genes/index.html (the
  `<p class="scope-rule">` sits inside `<div style="display: none;">`, pages.py:7181) and every
  genes/HGNC_*.html (src/chd_atlas/build/render.py:148, `.layout .scope-rule { display: none; }`).
- Sentence (:368-369): "The gene index and every gene page carry a one-line short form of the
  rule that links here."
  Sentence (:682-685): "Above them, one line says each classification is an external authority's,
  for the disease it names, which for many genes here is a syndrome rather than CHD, and links
  [`methods.html#scope`](#methodshtml) for the admission rule and the grade key."
  Sentence (:1764-1766): "Beside the card, the main column opens with the same one-line short form
  of the scope rule as the gene index, linking [`methods.html#scope`](#methodshtml)"
- Measured: [V] gives `none | Published on an external authority's classification…` on
  genes/index.html (the hidden div) and on genes/HGNC_11604.html. Only methods.html's own scope
  rule is visible. As a result, neither page visibly links to the grade key.
- Verdict: FALSE as rendered. Read strictly, :368-369's "carry" is true of the markup; :682-685
  and :1764-1766 describe a line no reader sees.
- Remedy (suggestion): unhide the scope rule on both pages (it carries the "syndrome rather than
  CHD" caveat), or drop the three sentences.

### A8. The validity table has no "Expert panel (ClinGen)" / "Other curators (GenCC)" groups
- Where: docs/data-api.md:1767-1768. Built: genes/HGNC_*.html (all 119).
- Sentence: "Below it, the mirrored validity table in two groups — **Expert panel (ClinGen)**
  first, then **Other curators (GenCC)** — each row giving who curated it, …"
- Measured: `grep -c 'Expert panel (ClinGen)\|Other curators (GenCC)' /tmp/dist/genes/*.html | awk -F: '{s+=$2} END{print s}'`
  → `0`. Every page has one `#validity-table` with one `<tbody>`. ClinGen rows come first, marked
  only by a "Gold standard" badge. `git show a89e958 -- src/chd_atlas/build/pages.py` removes
  `"<h3>Expert panel (ClinGen)</h3>"` and `"<h3>Other curators (GenCC)</h3>"`. The `_validity_section`
  docstring (pages.py:2003-2011, "Two groups under one heading") is stale in the same way.
- Verdict: FALSE
- Remedy (suggestion): "one table, ClinGen expert-panel rows first (badged *Gold standard*), then
  the other GenCC submitters".

### A9. "`null` where the row's comparator does not populate it" does not cover 77 control_cohort rows with no effect
- Where: docs/data-api.md:1045-1046, with the comparator table at :1077.
- Sentence: "Every key is present on every object, `null` where the row's comparator does not
  populate it, so a consumer never has to guard for a missing one." The table row reads
  "`control_cohort` | `n_control_carriers` / `n_controls` | both control fields | `odds_ratio`,
  `rate_ratio`".
- Measured: `python3 -c "import json,glob; from collections import Counter; rows=[r for f in glob.glob('/tmp/dist/genes/HGNC_*.json') for r in json.load(open(f))['burden']]; print(Counter((r['study'],r['comparator']) for r in rows if r['effect'] is None and r['effect_measure'] is None and r['effect_bound'] is None and r['ci_low'] is None))"`
  → `Counter({('PMID:34324492', 'control_cohort'): 77})`. All 77 PMID:34324492 rows are
  `control_cohort` but carry null `effect`, `effect_measure`, `effect_bound`, `ci_low` and
  `ci_high`, because that study published no effect (concordance.py:119). The burden section
  never says a control_cohort row can lack an effect measure.
- Verdict: MISLEADING. A consumer who reads null as "the comparator does not populate it" will
  take these 77 rows as malformed.
- Remedy (suggestion): say null also means "the study published none", and name PMID:34324492's
  77 rows (counts and p only).

### A10. "The committed corpus holds one dataset today": it holds two
- Where: docs/data-api.md:2376
- Sentence: "The committed corpus holds one dataset today, the `"profile"`-design `E-MTAB-6814`
  shown above."
- Measured: `python3 -c "import json;D=json.load(open('/tmp/dist/datasets.json'))['datasets'];print(len(D),[(d['id'],d['design']) for d in D]);print(json.load(open('/tmp/dist/manifest.json'))['counts']['datasets'])"`
  → `2 [('E-MTAB-6814', 'profile'), ('hoc', 'cell_type_profile')]` / `2`. The section's own `hoc`
  example (:2278-2289) and :2328 ("`hoc` is the one curated today") contradict it.
- Verdict: FALSE
- Remedy (suggestion): "two datasets today: the `profile`-design `E-MTAB-6814` and the
  `cell_type_profile`-design `hoc`". The next sentence, about contrasts, stays true.

### A11. The browse-row column list does not match the table
- Where: docs/data-api.md:673-679. Built: genes/index.html, and the same table on index.html.
- Sentence: "Each row carries the HGNC id and the symbol — both linked to that gene's page —
  `headline_confidence`, the disease it is graded for, the gene's `independent_datasets` as a
  strip (…), `validity_state`, `atlas_curation`, `burden_row_count` and the gene's lesion groups."
- Measured: `grep -o '<thead>.*</thead>' /tmp/dist/genes/index.html | grep -o '<th scope="col">[^<]*'`
  → `gene`, `symbol`, `confidence`, `graded for`, `rare variants across studies`, `heart expression`.
  `validity_state`, `atlas_curation` and the lesion groups appear only as `data-validity`,
  `data-curation` and `data-lesion` attributes, which no reader sees. `burden_row_count` is
  shown as "N tests". The "heart expression" column (peak TPM and percentile) is not mentioned.
- Verdict: FALSE
- Remedy (suggestion): list the six real columns, and say the three facets live in data
  attributes the filter menus use.

### A12. "ClinGen's definitions beside every chip"
- Where: docs/data-api.md:606-607. Built: genes/*.html, genes/index.html.
- Sentence: "The site's own pages carry ClinGen's definitions beside every chip for this reason."
- Measured: `grep -o '<span[^>]*class="chip chip-[a-z]*"[^>]*>' /tmp/dist/genes/HGNC_11604.html | sort | uniq -c`
  → 7 grade chips (4 definitive, 2 strong, 1 ungraded), none with a `title`. There is one
  definition sentence, beside the card's headline only ("Definitive: repeatedly demonstrated…").
  The six chips in the validity table have none. genes/index.html shows grades as plain text with
  no gloss, and its link to the grade key is hidden (A7).
- Verdict: FALSE
- Remedy (suggestion): "a gene page's card glosses its headline grade; methods.html#scope carries
  the key".

### A13. A study with only synonymous rows for a gene is drawn as "tested, no enrichment detected"
- Where: src/chd_atlas/build/concordance.py:215-220 (`family_state`: any row, synonymous
  included, gives `NO_ENRICHMENT` rather than `NOT_TESTED`). Built: the evidence-matrix cell on
  genes/HGNC_11128.html (SNAI1) and genes/HGNC_9906.html (RBFOX2).
- Sentence: the matrix legend "tested, no enrichment detected". Against it, methods.html says:
  "A consequence class or case group with no row is one the study did not report, which is not a
  null result".
- Measured: `python3 -c "import json,glob,itertools;[print(x['symbol'],s) for x in map(lambda f:json.load(open(f)),glob.glob('/tmp/dist/genes/HGNC_*.json')) for s,c in itertools.groupby(sorted((r['study'],r['consequence_class']) for r in x['burden']),key=lambda t:t[0]) if {k for _,k in c}=={'synonymous'}]"`
  → `SNAI1 PMID:42230622`, `RBFOX2 PMID:42230622`. For both genes, the bundle's PMID:42230622
  member state in `independent_datasets` is `no_enrichment`. A synonymous row is a negative
  control, not a test of the gene. The family tally is unaffected, because PMID:34324492 has a
  LoF row for both genes.
- Verdict: MISLEADING (2 genes, cell level only)
- Remedy (suggestion): treat a study whose rows for a gene are all synonymous as not tested, or
  label the cell "negative control only".

### A14. "of the 119 genes indexed today, 118 are `not_yet_curated`": it is 117
- Where: docs/data-api.md:2688
- Sentence: "This is not a hypothetical gap: of the 119 genes indexed today, 118 are
  `"not_yet_curated"` (measured 2026-10-25), and 21 carry no expert-panel grade at all"
- Measured: `python3 -c "import json;from collections import Counter;I=json.load(open('/tmp/dist/genes/index.json'))['genes'];print(Counter(r['atlas_curation'] for r in I), sum(r['headline_confidence'] is None for r in I))"`
  → `Counter({'not_yet_curated': 117, 'curated': 2}) 21`. TBX5 and NR2F2 are both curated. The
  doc's own :2842-2846 (96 + 21) and :40 already say 117.
- Verdict: STALE NUMBER
- Remedy (suggestion): "117"; better, point at the count instead of restating it.

### A15. "one curated gene-disease assertion": there are two
- Where: docs/data-api.md:305-308
- Sentence: "Today it is always `"in-development"` — one curated gene-disease assertion alongside
  mirrored ClinGen/GenCC validity for many more genes than that."
- Measured: `python3 -c "import json;print(json.load(open('/tmp/dist/manifest.json'))['counts']['assertions'])"`
  → `2` (TBX5 CHDA:AST:0000001, NR2F2 CHDA:AST:0000002).
- Verdict: STALE NUMBER
- Remedy (suggestion): drop the count, or point to `counts.assertions`.

### A16. "the index would hold one record": the assertion set covers two genes
- Where: docs/data-api.md:2673-2675
- Sentence: "Deliberately not the assertion set: keyed on that, the index would hold one record
  while the site published 119 bundles, so a visitor typing GATA4 would get nothing while
  `genes/HGNC_4173.json` was being served."
- Measured: `python3 -c "import json,glob;print([json.load(open(f))['symbol'] for f in glob.glob('/tmp/dist/genes/HGNC_*.json') if json.load(open(f))['assertions']])"`
  → `['TBX5', 'NR2F2']`. The GATA4 example still holds; GATA4 has no assertion.
- Verdict: STALE NUMBER
- Remedy (suggestion): "would hold only the genes with a lesion assertion (two today)".

### A17. index.html is described as "a doorway rather than the census", but it now carries the whole gene table
- Where: docs/data-api.md:318-320. Built: index.html.
- Sentence: "is a doorway rather than the census (since 2026-10-27): it says what the atlas is,
  states the figures a first-time reader needs, and sends them to `genes/index.html`."
- Measured: `grep -o 'data-confidence' /tmp/dist/index.html | wc -l` → `119`. index.html renders
  all 119 rows, the five filter menus and the glyph key (`render_combined_table`, added in
  a89e958). It is about 157 KB, against 156 KB for genes/index.html.
- Verdict: MISLEADING (describes the 2026-10-27 page)
- Remedy (suggestion): describe the page as it is now (hero, cards, the full filterable table,
  Sources, Data).

### A18. "the inline script only hides rows": it also re-sorts them
- Where: docs/data-api.md:695. Built: genes/index.html, index.html.
- Sentence: "**Every row is rendered by the build, and the inline script only hides rows.**"
- Measured: `python3 -c "import re;s=open('/tmp/dist/genes/index.html').read();print(re.findall(r'<script>(.*?)</script>',s,re.S)[0])"`.
  The script adds click-to-sort on every header (`th.addEventListener('click', …)`, `rows.sort(...)`
  with a `rankMap` from definitive to refuted) and re-appends the rows in the new order.
- Verdict: FALSE
- Remedy (suggestion): "the inline script only hides and re-orders rows; it adds none".

### A19. The cell-type chart's stated axis clips ACTC1's three cardiomyocyte medians without a mark
- Where: src/chd_atlas/build/cell_type_chart.py:47 (`_CAPTION`), with :60-61
  (`_AXIS_MAX = 10_000.0`). Built: genes/HGNC_143.html (ACTC1).
- Sentence: "Log axis, 0.1 to 10,000 CPM; a dot's size is its share of cells."
- Measured: ACTC1's pooled medians are Atrial cardiomyocytes 10,334, Ventricular cardiomyocytes
  10,086 and Non-chambered cardiomyocytes 12,071 CPM (from the bundle's
  `cell_type_expression.datasets[0].cell_types[].pooled.median`).
  `grep -o '<circle cx="[0-9.]*"[^>]*data-cell-type="[a-zA-Z]*CM"' /tmp/dist/genes/HGNC_143.html | head -3`
  → `cx="620.0"` three times, the x of the 10,000 tick. 12,071 is drawn at the same place as
  10,086.
- Verdict: MISLEADING
- Remedy (suggestion): extend the axis to the data, or mark clipped values ("> 10,000").

### A20. "No payload carries the path to a gene's page": manifest.json does
- Where: docs/data-api.md:1949
- Sentence: "**No payload carries the path to a gene's page.**"
- Measured: `python3 -c "import json;f=json.load(open('/tmp/dist/manifest.json'))['files'];print(sum(k.endswith('.html') for k in f), 'genes/HGNC_11604.html' in f)"`
  → `122 True`. The doc's own :350-354 says every page is checksummed in `files`.
- Verdict: FALSE (low harm)
- Remedy (suggestion): "No payload links a gene to its page; manifest.json lists page paths only
  as checksum keys".

### A21. The 2.8 changelog gives a grade split and a count as "now"
- Where: docs/data-api.md:187-191
- Sentence: "That is now 23 `definitive`, 1 `strong`, 9 `moderate`, 43 `limited` and 16 `null`.
  … one that treats presence in `genes/index.json` as "a panel called this a CHD gene" is now
  wrong for 16 genes."
- Measured: `python3 -c "import json,collections;print(collections.Counter(r['headline_confidence'] for r in json.load(open('/tmp/dist/genes/index.json'))['genes']))"`
  → `limited 45, definitive 42, None 21, moderate 10, strong 1`. The figures in the sentence are
  2.8's (92 genes), stated as "now" with no date. :590-593 gives today's split correctly.
- Verdict: STALE NUMBER
- Remedy (suggestion): "At 2.8 that was 23/1/9/43/16 (92 genes)", and drop or update "16 genes".

### A22. ClinGen's in-scope GenCC submissions: "111 rows over 109 genes"
- Where: docs/data-api.md:483-484
- Sentence: "it is GenCC's largest in-scope submitter (111 rows over 109 genes, measured
  2026-08-06)"
- Measured: loaded with `uv run python`: `load_curation(Path('.'))`, then the ClinGen and GenCC
  mirrors via `read_table`, then `gene_validity(clingen, gencc, in_scope={str(e.id) for e in corpus.chd_scope})`.
  Counting records with `source == 'gencc'` and submitter `ClinGen` gives `159 rows, 134 genes`.
  It is still the largest submitter (next is PanelApp Australia, 88), so the argument holds.
- Verdict: STALE NUMBER (dated, but measured before the 2.15 scope widening)
- Remedy (suggestion): re-measure and re-date, or drop the figures.

### A23. GDF1's card under-counts "also listed by"
- Where: src/chd_atlas/build/pages.py:1474-1478 (`others = asserting institutions − ClinGen −
  admitters`). Built: genes/HGNC_4214.html.
- Sentence: "also listed by 1 other curator (GenCC)"
- Measured: `python3 -c "import json;x=json.load(open('/tmp/dist/genes/HGNC_4214.json'));print(sorted({r['submitter'] for r in x['validity']['records'] if r['source']=='gencc'}-set(x['admitted_by']['submitters'])-{'ClinGen'}))"`
  → `['Illumina', 'Orphanet']`. Both are in the page's table; Illumina's record is "No Known
  Disease Relationship", which the count leaves out. The same check gives the right count on the
  other 118 pages.
- Verdict: MISLEADING (minor; the card names Illumina's no-relationship report on its own line)
- Remedy (suggestion): word the row "also asserted by", or count every listed submitter.

### A24. The hidden browse-page intro names an "atlas curation" column the table does not have
- Where: src/chd_atlas/build/pages.py:7184-7186. Built: genes/index.html (inside the
  `display: none` intro, see [V]).
- Sentence: "The atlas curation column is what says whether the atlas has curated a gene itself,
  and the graded for column names the disease each panel actually graded."
- Measured: the header row, from the command in A11, has no atlas-curation column; curation
  appears only as the `<select name="curation">` filter.
- Verdict: FALSE (hidden today; it becomes visible the moment A3 or A7 is fixed)
- Remedy (suggestion): point to the "any curation" filter before unhiding the intro.

### A25. The hidden landing claim "a gene page says so wherever it shows their results" fails in the evidence matrix
- Where: src/chd_atlas/build/landing.py:161-166. Built: index.html (inside the hidden div, see A4).
- Sentence: "Of the 3 rare-variant studies whose results appear here, 2 are by an author of this
  atlas, and a gene page says so wherever it shows their results."
- Measured: the 3 and the 2 are right (`own_lab` in publications.json). Every Audain study block
  carries "Declaration: this study is by an author of this atlas" (190 of 190), but the
  evidence matrix shows their results with no mark:
  `grep -o '<th scope="row">Audain[^<]*</th><td>.\{0,300\}' /tmp/dist/genes/HGNC_11604.html`
  → "Audain E et al. 2021 (no row for this gene) + Audain E et al. 2026 … OR ∞ (CI 0.786–∞)
  p 0.0316", with no declaration.
- Verdict: FALSE (minor; hidden today)
- Remedy (suggestion): "each study's block and the card say so", or mark own-lab rows in the
  matrix.

### A26. Contrast records are said to have `[]`/`null` profile fields; that is only a default
- Where: docs/data-api.md:2317-2319
- Sentence: "`cardiac_tissues` and `stages` are `[]` and `detection_floor`/`floor_source`/
  `quantile_estimator` are `null` on a `"contrast"` record — never omitted — so every record has
  one shape whichever design it declares."
- Measured: src/chd_atlas/models/dataset.py:236-241 (`a_profile_dataset_is_fully_declared`) says
  these fields are "Left optional on a contrast dataset". No validator rejects them on
  `design: contrast`; only `cell_types` and `regions` have such a rule (line 230). No contrast
  record exists in the build. Command: `grep -n "CONTRAST" src/chd_atlas/models/dataset.py src/chd_atlas/validate/*.py`
  → only dataset.py:221.
- Verdict: MISLEADING ("never omitted" is true; the values are a default, not a guarantee)
- Remedy (suggestion): "default to `[]`/`null` on a contrast record", or add the validator the
  sentence assumes.

### A27. The text box matches more than "id or symbol"
- Where: docs/data-api.md:679-680. Built: genes/index.html.
- Sentence: "Above the table sit a text box matching id or symbol"
- Measured: the box filters on `data-search`, which holds the id, symbol, gene name, graded-for
  disease and lesion groups.
  `python3 -c "import re;s=open('/tmp/dist/genes/index.html').read();r=re.findall(r'<tr data-search=\"([^\"]*)\"',s);print(sum('noonan' in x for x in r))"`
  → `12`.
- Verdict: MISLEADING (minor)
- Remedy (suggestion): "a text box matching id, symbol, name, graded-for disease or lesion group".

### A28. With one placed median, the heart axis is padded below the floor it says sub-floor stages sit beneath
- Where: src/chd_atlas/build/pages.py:5465-5498. Built: genes/HGNC_11743.html (TFAP2B),
  genes/HGNC_14357.html (MMP21), genes/HGNC_3814.html (FOXH1).
- Sentence: "Median tpm in heart per stage on a log scale from 0.316 to 3.16, fitted to this gene,
  … Ticks on the axis mark stages below this dataset's detection floor (1 tpm), beneath the foot
  of the scale rather than at a value on it."
- Measured: on these three pages the placed heart medians are a single value (1.0, 1.0 and
  2.0 tpm), and the scale is padded to 0.316-3.16 and 0.632-6.32. The other 109 axes equal their
  min and max placed medians (checked by comparing each bundle's placed heart medians with the
  "from A to B" text on its page). The drawn scale therefore covers 0.316-1 tpm, where no value
  can be placed.
- Verdict: MISLEADING (minor)
- Remedy (suggestion): floor the padded scale at the detection floor, or say "padded around its
  one placed median".

### A29. The expression "How to read this chart" fold is described as always present
- Where: docs/data-api.md:1888-1891. Built: genes/HGNC_4214.html (GDF1) lacks it.
- Sentence: "The rest of the key is folded under *How to read this chart*, ending with a link,
  *How to read expression figures*, to `methods.html#expression`"
- Measured: `grep -L 'How to read expression figures' /tmp/dist/genes/HGNC_*.html` →
  `/tmp/dist/genes/HGNC_4214.html`. GDF1 draws no chart and no shaded grid cell; pages.py:5919-5921
  emits the fold only when one is drawn.
- Verdict: MISLEADING (low)
- Remedy (suggestion): "Where a chart or a shaded grid is drawn, …".

### A30. "Every count in this section is asserted against a real build"
- Where: docs/data-api.md:1008-1009
- Sentence: "Every count in this section is asserted against a real build by
  `tests/test_docs_match_the_build.py`, so it fails rather than rots when a study lands."
- Measured: `for n in '398' 'covering' '185 genes' '3,876' '45,082' '0.991' '0.0068' '777'; do printf '%-10s %s\n' "$n" "$(grep -c -- "$n" tests/test_docs_match_the_build.py)"; done`
  → 0 for each. These numbers in the section (:985, :991, :1134, :1137-1138, :1151, :1166-1167)
  are untested. All are correct today.
- Verdict: MISLEADING (low; the guarantee is overstated, and no number is wrong)
- Remedy (suggestion): "Most counts …", or extend the test to these.

### A31. "Most published genes carry the populated shape": all do
- Where: docs/data-api.md:1341-1342
- Sentence: "Most published genes carry the populated shape today, not the empty one above."
- Measured: `python3 -c "import json,glob;print(sum(1 for f in glob.glob('/tmp/dist/genes/HGNC_*.json') if json.load(open(f))['expression_profile']['datasets']), len(glob.glob('/tmp/dist/genes/HGNC_*.json')))"`
  → `119 119` (= `counts.profile_genes`).
- Verdict: MISLEADING (low)
- Remedy (suggestion): "Every published gene (`counts.profile_genes`) carries …".

### A32. The rule for keeping a paper out of publications.json leaves out the variants mirror
- Where: docs/data-api.md:1609-1611 and :2182-2184
- Sentence (:2182-2184): "A paper is left out of this file while literature findings cite it and
  none of them is approved, unless an assertion, functional record, featured entry, dataset or
  burden row also cites it". :1609-1611 has the same list.
- Measured: `sed -n 117,154p src/chd_atlas/build/literature.py; sed -n 394,404p src/chd_atlas/build/runner.py`.
  `withheld_publications` also exempts any PMID in the variants mirror's `publication` column, and
  the combined-burden `study`. There is no variants mirror today (`variants/index.json` =
  `{"shards": []}`), so no published byte differs.
- Verdict: MISLEADING (low; latent)
- Remedy (suggestion): "… dataset, burden or combined-burden row, or variant row also cites it".

### A33. The list of HPO branches excluded from "structural" is missing one
- Where: docs/data-api.md:2551-2553
- Sentence: ""Structural" is HPO's `Abnormal heart morphology` branch minus its cardiomyopathy,
  cardiomegaly, tumour, pericardium, endocardium and imaging-finding branches, plus eight lesions
  HPO files elsewhere"
- Measured: `sed -n 52,61p src/chd_atlas/orphanet_scope.py` → eight `EXCLUDED_BRANCHES`, among them
  `HP:0031331` "Abnormal cardiomyocyte morphology" ("cellular, not structural"), which the list
  does not name. `awk '/^\[Term\]/{p=0} /^id: HP:0031331$/{p=1} p' ontologies/hp-2026-06-23.obo`
  → `is_a: HP:0001627 ! Abnormal heart morphology`, a direct child, not under cardiomyopathy.
- Verdict: MISLEADING (low)
- Remedy (suggestion): add "cardiomyocyte (cellular)", or name `EXCLUDED_BRANCHES`.

### A34. The Orphadata `file_date` is shown as a bare date
- Where: docs/data-api.md:2605-2606
- Sentence: "Both Orphadata products carry them, from their `<JDBOR>` header —
  `"1.3.42 / 4.1.8 [2025-03-03]"`, dated `2026-06-23` —"
- Measured: `python3 -c "import json;print([(s['id'],s['file_version'],s['file_date']) for s in json.load(open('/tmp/dist/sources.json'))['sources'] if s['file_date']])"`
  → `file_date` is `"2026-06-23 07:53:50"` (product1) and `"2026-06-23 07:57:18"` (product4).
- Verdict: MISLEADING (low; a consumer who parses it as YYYY-MM-DD fails)
- Remedy (suggestion): quote the real value and say the front page prints its date part.

### A35. The 2.10 changelog says the new fields are "described above"; they are described below
- Where: docs/data-api.md:194-195
- Sentence: "`2.10` added `has_no_association_report` and `no_association_reported_by` to every
  gene bundle and browse row, described above."
- Measured: `grep -n 'has_no_association_report' docs/data-api.md` → the first description is
  at :612.
- Verdict: FALSE (minor; the cross-reference points the wrong way)
- Remedy (suggestion): "described below", with a link.

### A36. The `profiles` ranking counterfactual does not name its dataset and does not reproduce on the published one
- Where: docs/data-api.md:2009-2017
- Sentence: "… measured: ranking `profiles` the way the other three tables are ranked put 0 of 14
  heart rows in a 25-row slice at 14 stages per organ (every row ties at "no FDR", so the
  tie-break falls back to an alphabetical tissue order that heart loses), and ranking
  cardiac-first with no reservation at all over-corrected, leaving 0 of 6 comparison organs in a
  25-row slice at 23 stages per organ."
- Measured: the `select_top` docstring in build/omics.py points to a synthetic layout. On the
  committed shard, sorting by the FDR-style key `('dataset','gene','tissue','stage')` puts 5 heart
  rows in the top 25 for all 185 genes, because only `forebrain` sorts before `heart`.
  Cardiac-first with no reservation leaves 1 to 6 comparison organs, never 0, because the heart
  series has at most 19 stages and cannot fill 25 rows. (A re-implementation of both orderings
  over omics/profiles/E-MTAB-6814.json; per-gene distribution
  `[(1, 32), (2, 63), (3, 52), (4, 28), (5, 8), (6, 2)]`.)
- Verdict: UNVERIFIABLE. The setting that was measured is not named, and "an alphabetical tissue
  order that heart loses" is not true of the published dataset, where heart sorts second.
- Remedy (suggestion): say the figures come from a synthetic layout, and give the committed-data
  figures beside them.

### A37. "The four studies this schema was designed against"
- Where: docs/data-api.md:1048-1049
- Sentence: "The four studies this schema was designed against do not count the same thing:"
- Measured: `python3 -c "import json,glob; from collections import Counter; print(Counter((r['study'],r['count_unit']) for f in glob.glob('/tmp/dist/genes/HGNC_*.json') for r in json.load(open(f))['burden']))"`
  → 3 studies, 3 count units. `grep -rn "four studies" src docs` → only this line. Nothing in the
  repository names a fourth.
- Verdict: UNVERIFIABLE (possibly design history; next to three studies it reads as a stale
  count)
- Remedy (suggestion): "The three studies published today …", or name the fourth.

### A38. The index.html footer links the API doc in `chd-atlas/chd-atlas`, not this repository
- Where: src/chd_atlas/build/landing.py:82 (`_REPOSITORY_URL = "https://github.com/chd-atlas/chd-atlas"`),
  used at landing.py:438. Built: the index.html footer.
- Sentence: "Read the data API documentation before writing against these files, …" (href
  `https://github.com/chd-atlas/chd-atlas/blob/main/docs/data-api.md`)
- Measured: `curl -s -o /dev/null -w '%{http_code}' https://github.com/chd-atlas/chd-atlas.github.io/blob/main/docs/data-api.md`
  → `200`. `chd-atlas/chd-atlas` could not be reached from this session (the proxy does not allow
  that repository), and `gh api -X GET search/repositories -f q='org:chd-atlas'` lists only
  chd-atlas/chd-atlas.github.io. The comment at landing.py:76-79 records a 200 measured
  2026-10-23.
- Verdict: UNVERIFIABLE (the link target could not be reached; this build lives in
  chd-atlas/chd-atlas.github.io)
- Remedy (suggestion): open the link from outside; if the repository was renamed, point it at
  chd-atlas.github.io.

### A39. "Regions pool all four ages"
- Where: src/chd_atlas/build/methods.py:823. Built: methods.html#expression.
- Sentence: "Regions pool all four ages."
- Measured: `python3 -c "import json;d=[x for x in json.load(open('/tmp/dist/datasets.json'))['datasets'] if x['id']=='hoc'][0];print(d['developmental_stage'],d['n_samples'],d['stages'])"`
  → `9 to 15 weeks post conception 8 []`. Neither datasets.json nor curation/datasets/hoc.yaml
  records the ages.
- Verdict: UNVERIFIABLE (the 4 cannot be checked against the build)
- Remedy (suggestion): record the ages in hoc.yaml, so the count is derived rather than asserted.

---

## B. Curator docs

### B1. curator-rules.md names a `count_unit` value the model refuses
- Where: curator/curator-rules.md:73
- Sentence: "`carriers` and `screened` are counted in `count_unit`. With `probands` or
  `individuals` a carrier is a person."
- Measured: a finding's `count_unit` is `ScreenUnit`, which has two values.
  `uv run python -c "from chd_atlas.vocab import ScreenUnit; print([u.value for u in ScreenUnit])"`
  → `['probands', 'families']`. A `Finding` built with `count_unit='individuals'` is refused:
  `count_unit: Input should be 'probands' or 'families'` (SCHEMA001 at load). `individuals` is a
  value of the burden mirror's `CountUnit` (src/chd_atlas/vocab.py:434, 470), a different enum.
  The extract skill gets it right (curator/skills/literature-extract/SKILL.md:32-33, "`count_unit`
  (`probands` or `families`)"), so the house rules contradict the skill that loads them.
- Verdict: FALSE
- Remedy (suggestion): "With `probands` a carrier is a person."

### B2. One unregistered phenotype switches off REF009 for the whole finding
- Where: curator/curator-rules.md:77-78
- Sentence: "`lesion_groups` agree with the `phenotypes` that `curation/phenotypes.yaml` maps to a
  group (`REF007`, `REF009`); a phenotype it does not map is not checked."
- Measured: `lesion_group_faults` (src/chd_atlas/validate/referential.py:51) runs REF009 only when
  every phenotype is a key in phenotypes.yaml, which holds 8 terms. With the registry loaded
  through `load_curation(Path('.'))`:
  `uv run python -c "from chd_atlas.validate.referential import lesion_group_faults; …"`
  `['HP:0001631'], [septal, conotruncal]` → `[('REF009', "lesion groups ['conotruncal'] are not justified by any declared phenotype")]`;
  `['HP:0001631', 'HP:0000001'], [septal, conotruncal]` → `[]`.
- Verdict: MISLEADING. A curator would expect an unjustified group to be refused whenever the
  mapped phenotypes do not justify it; one unmapped phenotype exempts every group.
- Remedy (suggestion): add "and if any phenotype is unmapped, no declared group is checked for
  justification (REF009)".

### B3. The extract skill's list of exit-1 causes for `add-publication` is incomplete
- Where: curator/skills/literature-extract/SKILL.md:79-81
- Sentence: "Exit 1 means the PMID is already in `curation/publications.yaml` (nothing to add),
  PubMed could not be reached, or the file does not validate; the error says which."
- Measured: `literature add-publication` (src/chd_atlas/cli.py:298-338) also exits 1 when PubMed
  is reached but returns no record for the PMID ("PubMed returned no record for PMID:<n>",
  cli.py:331), and when curation/publications.yaml does not exist (cli.py:326). Read from the
  code; the command was not run, because it needs the network.
- Verdict: MISLEADING (presented as complete; a mistyped PMID falls outside the list)
- Remedy (suggestion): add "PubMed has no record for that PMID".

### B4. The extract skill states the publications.json withholding rule without its exception
- Where: curator/skills/literature-extract/SKILL.md:76-78
- Sentence: "The record stays out of `publications.json` while findings cite it and none is
  approved (a record no finding cites publishes)"
- Measured: `withheld_publications` (src/chd_atlas/build/literature.py:117-155) withholds a record
  only when a finding cites it, no approved finding does, and nothing else the build publishes
  cites it (assertion or functional evidence, featured, datasets, burden or combined `study`,
  variants `publication`). This holds for a record `add-publication` has just created, until
  something else cites it.
- Verdict: MISLEADING (low; the same gap as the Known gene_literature.py docstring, reported so
  the scheduled fix covers this copy)
- Remedy (suggestion): append "and nothing else the site publishes cites it".

---

## Checked and true

### docs/data-api.md
- docs/data-api.md:6 — the build refuses before writing if validate reports an error (`build_site` raises `BuildRefused` before any emit)
- docs/data-api.md:15 — two builds of one commit are byte-identical (rebuilt into the scratchpad; `diff -r` empty)
- docs/data-api.md:22 — index rows carry `bundle` and search records carry `path`; all 138 paths exist
- docs/data-api.md:40, :92, :476, :490, :560, :783 — 117 of 119 genes are not_yet_curated
- docs/data-api.md:46-64 — manifest example: all 19 counts keys, schema 2.21, status and key names match the build
- docs/data-api.md:73-75 — `files` digests equal the sha256 of the served bytes, including the compressed search index (253/253)
- docs/data-api.md:76 — manifest.json is absent from its own `files`
- docs/data-api.md:78 — `source_commit` is null outside a checkout
- docs/data-api.md:83-96 — 8 corpus keys + 11 computed keys = 19; publications 11 = len(publications.json)
- docs/data-api.md:103-104 — genes 119 = len(genes array)
- docs/data-api.md:105-109 — burden_rows 1210 = Σ burden_row_count; burden.tsv has 181 genes, 63 unpublished
- docs/data-api.md:110-112 — every bundle has 2 families = cohort_families
- docs/data-api.md:113-121 — burden_surviving_genes 32 recomputed; PMID:40127276 correction_genes 248
- docs/data-api.md:122-125 — burden_combined 108 = Σ bundle burden_combined
- docs/data-api.md:126-134 — profile_genes 119, profile_datasets 1, cell_type_genes 118, cell_type_datasets 1; datasets 2 ≠ profile_datasets
- docs/data-api.md:143-149 — cohort_families computed over all burden rows (runner.py:380)
- docs/data-api.md:233 — HOS→TBX5, VEGF→VEGFA, HTX1→ZIC3 in the search index
- docs/data-api.md:236 — MADH7 names SMAD6 and SMAD7
- docs/data-api.md:254-257 — submitted_on/sgc_id on every validity record, null on all 121 ClinGen records
- docs/data-api.md:258-260 — the 6 identity keys on all 119 bundles
- docs/data-api.md:261-264 — summary_caveat on every publication record and on the featured embedded one; null on 10 of 11
- docs/data-api.md:265-291 — correction_genes/genome_wide_threshold on every publication; DDD↔DECIPHER overlaps_with; 20 cohorts including 7 pmid*_cases; file_version/file_date on every source; syndrome on every scope term; n_chd_case_carriers on every burden row; review_status/review_note on every cohort and scope term
- docs/data-api.md:292-298 — cell_type_expression on every bundle, `{"datasets": []}` on CRIPTO; hoc is cell_type_profile
- docs/data-api.md:299-304 — a literature object on every bundle; StudyType has case_report, cohort_sequencing, gwas, candidate_variant
- docs/data-api.md:306 — status "in-development"
- docs/data-api.md:315-316 — 122 pages; index.html is the entry point to the other 121
- docs/data-api.md:324-328 — the two non-manifest cards (118, 2) and the hero's 98/21 split are present and visible
- docs/data-api.md:328-333 — "(as far as the cohort registry records)" on the hero and on a gene card
- docs/data-api.md:342-344 — Sources names every permitted_with_attribution entry, with CC BY licences and the Orphadata citation
- docs/data-api.md:344-349 — links to genes/index.html, genes/index.json, manifest.json, sources.json, methods.html; no external script, link, img, url() or @import
- docs/data-api.md:351-355 — 122 HTML pages, each with an entry in `files`
- docs/data-api.md:360-363 — methods.html is linked from every nav and every gene card; seven sections with ids in the stated order
- docs/data-api.md:370-376 — #validity content; every gene page links #validity from its validity intro
- docs/data-api.md:384-385 — TBX5 and NR2F2 open curated evidence with "How to read curated evidence"
- docs/data-api.md:386-404 — #rare-variants content; each forest fold ends with the #rare-variants link
- docs/data-api.md:405-417 — #cohorts sorted by id, 20 records, "Pending review:" notes; the gene-page cohorts fold links #cohorts
- docs/data-api.md:418-434 — #expression content; CRIPTO's one-line form links #expression
- docs/data-api.md:435-441 — #literature content
- docs/data-api.md:453-459 — 119 genes; ClinGen 98; GenCC agreement 21
- docs/data-api.md:466-473 — SCP005 checks attribution; SCP006 uses the Frequent (79-30%) floor
- docs/data-api.md:479-486 — ClinGen's own GenCC rows excluded (`agreeing_records`); a ClinGen Disputed/Refuted vetoes (`_clingen_contests`)
- docs/data-api.md:497-533 — the TBX5 browse-row example matches the build value for value
- docs/data-api.md:536-538 — rows ordered by HGNC id as a string
- docs/data-api.md:540-542 — symbol falls back to the HGNC id; no row uses the fallback today
- docs/data-api.md:548-552 — TBX5's only ClinGen mirror row is the Syndromic Disorders GCEP, Holt-Oram
- docs/data-api.md:571-574 — null headline iff not expert_curated (21 = 21, as sets)
- docs/data-api.md:580-585 — the six Supportive-only genes are single-record, submitter_curated and unpublished
- docs/data-api.md:590-594 — 42/45/10/1/21 and 98/21
- docs/data-api.md:597-599 — "uncurated" absent; Supportive-only genes unpublished
- docs/data-api.md:603 — 45 limited is the plurality
- docs/data-api.md:608-616 — discordance semantics; no_association_reported_by sorted, and empty iff the flag is false
- docs/data-api.md:625-632 — GDF1 is the only has_no_association_report gene; 9 NKA-only genes in the mirrors, 0 published
- docs/data-api.md:633-651 — confidence_by_lesion_group keys = lesion groups; 96 rows with `{}`; TBX5 and NR2F2 maps as stated
- docs/data-api.md:652-659 — conflicting_lesion_groups empty everywhere
- docs/data-api.md:660-666 — the two lesion-group fields occur only in genes/index.json
- docs/data-api.md:667-668 — the four row counts equal the bundle array lengths on all 119
- docs/data-api.md:672-673 — the browse page has the same 119 rows
- docs/data-api.md:680-682 — five menus whose options are the values present
- docs/data-api.md:685-694 — 21 null rows render an em dash; no null option in the confidence menu; 6 contested rows carry the conflicting-evidence chip
- docs/data-api.md:696-701 — 119 rows in static HTML; inline style and script; no external request
- docs/data-api.md:713-756 — the 25 top-level keys of the bundle example match; TBX5 values match
- docs/data-api.md:760-771 — identity keys on all genes; ncbi_gene int ≥ 1; ENSG + 11 digits; RBM10 uniprot P98175; GEN004 exists
- docs/data-api.md:772-779 — the assertion field list matches exactly
- docs/data-api.md:783-788 — uncurated bundles have empty assertions, publications and functional
- docs/data-api.md:789-799 — omics/variants/burden always present; publications = evidence-cited PMIDs, sorted lexically; assertions and functional ordered by id
- docs/data-api.md:811-817 — validity.state = validity_state; has_source_discordance repeats the flag
- docs/data-api.md:819-835 — 13-key record shape on all 591 records; the TBX5 ClinGen example is exact
- docs/data-api.md:837-845 — Supportive is the only term mapping to null (43 records)
- docs/data-api.md:847-850 — sop/date/gcep only on ClinGen; submitter/submitted_on/sgc_id only on GenCC (121/470)
- docs/data-api.md:852-864 — the submitted_on prefix is enforced (tables.py:482); 30,410 mirror rows dated
- docs/data-api.md:865-879 — 470 GenCC records, 189 null report_url, 76 PanelApp homepage
- docs/data-api.md:880-883 — the ClinGen mirror spans SOP4 to SOP12
- docs/data-api.md:889-905 — 98/21 states; no empty records
- docs/data-api.md:911-934 — admitted_by examples real; 98/21; every key on every gene
- docs/data-api.md:936-946 — admitted_by.classification = headline_confidence on all 119
- docs/data-api.md:948-952 — the GDF1 asserted_by example is exact (count 5)
- docs/data-api.md:956-966 — dissent excluded, Supportive counted (`_asserts`)
- docs/data-api.md:973-977 — neither key in genes/index.json; 98/21
- docs/data-api.md:982-986 — 1210 rows across 118 of 119 genes; 3 studies with 77/216/917 rows; burden.tsv 1,805 rows for 181 genes; 595 unpublished rows over 63 genes
- docs/data-api.md:991 — genes.tsv 185 genes, 119 published
- docs/data-api.md:997-1007 — CNV is deletions only; CRIPTO has `burden: []`; burden_row_count == len(burden)
- docs/data-api.md:1012-1043 — the example equals TAB2's syndromic LoF row
- docs/data-api.md:1045 — every burden object has the same 28 keys
- docs/data-api.md:1052-1056 — count_unit values individuals/alleles/de_novo_mutations (vocab.CountUnit)
- docs/data-api.md:1073-1077 — the comparator table matches validate/burden.py and vocab.EFFECT_MEASURES
- docs/data-api.md:1079-1087 — control_cohort 1102, mutation_model 108, none 0; mutation_model rows all PMID:40127276, n_cases 3,887; BUR001/002/006/007
- docs/data-api.md:1092-1095 — 72 damaging rows; n_case_carriers = lof + missense_damaging
- docs/data-api.md:1103-1107 — TAB2 ci_low 28.1; 29 unbounded rows, 14 with ci_low < 1; page renders "OR ∞ (95% CI 28.1–∞)"
- docs/data-api.md:1111-1126 — PMID:42230622 1,438 mirror rows; PMID:34324492 tests lof and cnv_deletion; WDR5 0/0; 339 synonymous rows, all PMID:42230622
- docs/data-api.md:1133-1143 — lof 401, synonymous 339, missense_damaging 398, damaging 72; CHD7 p 0.0068 / 0.991; pvalue_adjusted null exactly on PMID:42230622 (tests_reported 138609)
- docs/data-api.md:1150-1157 — taa_cases 777; cohorts.json added in 2.9
- docs/data-api.md:1166-1172 — PMID:42230622 3,876 / 1,471 / 2,405 vs 45,082; PMID:34324492 6,746 vs 14,024
- docs/data-api.md:1176-1194 — n_chd_case_carriers (2.19) non-null on 10 PMID:34324492 rows; BUR024; source audain2021_s6 / audain2021_s8; the S8 card sentence
- docs/data-api.md:1198-1225 — burden_combined `[]` when empty; 108 results across 36 genes; mirror 123 rows for 41 genes; BUR020 key; the TBX5 example
- docs/data-api.md:1234-1235 — NODAL 7.57e-10 vs 9.86e-14; threshold 2.6e-06
- docs/data-api.md:1249-1310 — the independent_datasets example; disjoint collections apart from ddd↔decipher; members sorted; family state = strongest member; enrichment rule; KDM6A; CRIPTO two not_tested
- docs/data-api.md:1316 — browse header "rare variants across studies"
- docs/data-api.md:1338-1352 — profile counts; 21 stages; the 7 wpc phase block; tpm; n_genes 58735
- docs/data-api.md:1386-1403 — profiles mirror columns; four omics modalities; quantile_shard path
- docs/data-api.md:1408-1446 — stage ordering, the alphabetical counter-example, phase outcomes, 6 wpc windows, phase_ids rule
- docs/data-api.md:1458-1488 — heart_looping derived; specificity reasons
- docs/data-api.md:1498-1522 — the τ example; highest_in null only on exact ties; method strings; quartiles null when n_samples < 3
- docs/data-api.md:1550 — 101 breakpoints (12,019 = 119 × 101)
- docs/data-api.md:1557-1586 — hoc covers 118 of 119; 9 to 15 weeks; specificity null on 7 genes; the TBX5 cell-type example
- docs/data-api.md:1591-1617 — literature added in 2.21; the empty shape on every bundle; findings approved only, digest dropped, sorted (year, numeric PMID, id) (code)
- docs/data-api.md:1619-1639 — Literature section placement and wording (pages.py:6698-6750; code only, no data)
- docs/data-api.md:1648-1701 — card links; conflicting/disagree chips on the 6 flagged genes; "also listed by" rows; the gloss on 98 graded cards; the disputed lines; GDF1's Illumina line; scope lines
- docs/data-api.md:1702-1757 — card rare-variant lines, #study ids, 248-gene panel and genome-wide marks, 138,609 scan, "Tested in 2 independent datasets …", power caveat on 118 genes
- docs/data-api.md:1769-1781 — "Autosomal dominant"; SOP-labelled links; the ClinGen-copy rule; the table intro text
- docs/data-api.md:1786-1859 — the rare-variant intro, matrix legend, overlap notice on 72 pages, declaration on 118 pages, 3.6e-07, union line, JL-FDR block, forest caption and fold
- docs/data-api.md:1862-1887 — expression opening lines, figure links, axis and tick lines, GDF1's no-chart sentence
- docs/data-api.md:1892-1917 — "every figure, all 21 stages"; 117 uncurated; 96/21 sentences; ACTC1, ELN, FLNA, SMAD6, MNS1; TBX5 "Curated by the atlas: 1 lesion assertion (septal)."
- docs/data-api.md:1921-1943 — card literature and functional wording (code); gene pages have no script, link or src=
- docs/data-api.md:1949-1950 — a browse row's `bundle` is the JSON path
- docs/data-api.md:1956-2025 — omics modalities; the TBX5 omics example; count = shard rows; percentile = placement.median_percentile (14,161 rows); TOP_N 25; stratified profiles slice; reservation ≤ 12; six non-cardiac tissues; `genes` on every row
- docs/data-api.md:2039-2093 — `genes == [gene]` on profiles rows; genes.tsv 185 × 10 columns, all reaching a byte; 2.13 figures; MADH7; RBM10 P98175; MINAS-60 130 aa
- docs/data-api.md:2100-2128 — profile_quantiles: no gene column; the example row; 12,019 rows sorted; quantile_shard
- docs/data-api.md:2136-2154 — variants/index.json `{"shards": []}`; karyotype order; shard shape
- docs/data-api.md:2158-2214 — the publications example and its 15 keys; summary_caveat ≤ 240 and printed on 10 gene pages; PMID:40127276 tests_reported null; correction_genes 248, 2.6e-06; sorted by id
- docs/data-api.md:2220-2254 — featured: one entry, embedded record; the phenotypes example; lesion_group null on extracardiac terms
- docs/data-api.md:2258-2374 — datasets examples; 18 keys; design vocabulary; cell_types/regions rule; 12 cell types in 5 compartments; wpc null on 8 of 21; order 1..21; PRF011/012/014
- docs/data-api.md:2382-2463 — cohorts.json: the taa_cases example; every cited cohort resolves; BUR009; url null on 11 of 20; overlaps symmetric; two pending_review cohorts; sorted by id; 11 of 189 sources and 6,211 of 7,958 cases; counts 20/2
- docs/data-api.md:2473-2565 — 105 scope terms; 69 attributed / 36 Orphanet; the MONDO:0012547 example; nine syndrome terms; two pending_review terms; admitted_by vocabulary; orphanet block; ADDED_BRANCHES; 10 pulmonic-stenosis terms; PTPN11 "Why these diseases are in scope"
- docs/data-api.md:2571-2627 — the hpo source entry; file_version; licences of clingen/gencc (CC0) and orphanet/mondo_sssom (CC BY 4.0)
- docs/data-api.md:2635-2687 — search records exact; terms deduplicated; 48 genes with a withdrawn symbol; MADH7 only collision; kinds and paths; gene set = genes/index.json (185 − 119 = 66)
- docs/data-api.md:2697-2860 — headline = strongest in-scope ClinGen grade; 21 null; six Supportive-only unpublished; flags on 8 of 185, 6 published; the contested table for all six; LEFTY2 and A2ML1 unpublished; lesion-group fields index-only; all-or-none conflicting_lesion_groups; 96 + 21; atlas_curation from assertions only

### Built pages
- index.html (landing.py:388-396) — "119 genes published" = manifest genes = len(genes/index.json)
- index.html (landing.py:392) — "2 independent datasets" = 2 cohort families
- index.html (landing.py:398-401) — 98 on a ClinGen panel, 21 on two or more GenCC submitters
- index.html stat cards — curated 2, rare-variant result 118, surviving 32, developmental 119, cell-type 118, recomputed from the bundles
- index.html Sources — every licence matches sources.json
- index.html Data — manifest lists every built file except itself (253 of 254), sha256 match
- index.html, genes/index.html — "N tests" = burden_row_count; "X of Y tested · Z corrected" = independent_datasets
- methods.html (orphanet_scope.py:71-78) — eight lesions filed elsewhere; the 30% threshold
- methods.html — every scope term names its authority (105/105); every validity disease is a scope term (591/591)
- methods.html — "No gene here is published on a disputed or refuted grade"
- methods.html — the definitions of "conflicting evidence" and "sources disagree" match the flags (6 and 6)
- methods.html, gene pages — a ClinGen GenCC copy is listed only where it differs (0 of 121 differ); "ClinGen also submits…" on exactly the 98 panel pages
- gene pages — "No GenCC submitter other than ClinGen …" (38) and "No ClinGen expert panel record …" (21) match the records
- gene pages — rung gloss lines match the headline grades (42/45/10/1)
- gene pages — the 31 "Graded for X as a whole syndrome …" lines and 27 Orphanet notes match scope.json
- gene pages — disputed card lines (BRAF, SOS1, KRAS, MAP2K1, PTPN11, RAF1) and GDF1's Illumina line match the records
- gene pages — "Tested in N independent dataset(s) …" lines correct for all 118 genes with rows; CRIPTO's no-result line
- gene pages — the surviving-rows card complete for all 32 genes; genome-wide marks exactly at p < 2.6e-6; "within its 248-gene panel" only on PMID:40127276
- gene pages — 93 study headline badges match the rows under the 0.05 rule
- gene pages — the DECIPHER/DDD non-independence note on exactly the 72 pages carrying both studies
- gene pages — the 777/7,958 PMID:34324492 note on all 77 pages with its rows; the S8 card note on the 10 surviving pages
- gene pages — 138,609 comparisons, 3.6e-07; PMID:42230622 publishes no corrected p
- gene pages (Sierant) — local FDR smaller than p on 17 of 108; Fisher's method not reproduced on 81 of 108; threshold 2.6e-6
- gene pages — 10 "this study published no … row" axis notes, none missing
- gene pages — the Declaration on 190 of 190 Audain blocks and on no Sierant block
- gene pages — the heart floor sentences (112 + 7); axes equal min and max placed medians on 109 pages
- gene pages — "A percentile places this gene among 58,735 genes"
- methods.html, gene pages — no percentile below the floor (2,830/2,830 unplaced); τ withheld on 115 stages; 2,384 τ recompute to 1e-9; "tied" peaks
- gene pages — cell-type τ on 111 pages; "8 hearts at 9 to 15 weeks"; below_floor = median < 1 on 8,024 values; quartiles only with ≥ 3 hearts
- methods.html#cohorts — 20 records; 11 case cohorts summing to 6,211 of 7,958; GAIN 7,993; gnomAD denominators
- methods.html — each developmental dataset gives its floor and floor source
- methods.html, gene pages — every internal link and anchor resolves (0 broken)
- gene pages — Curated evidence and Publications sections only on TBX5 and NR2F2; each Publications list equals the cited PMIDs
- methods.html literature — the query excludes preprints and reviews (pubmed.py:116)

### Curator docs
- curator/README.md:4-5 — curator-rules.md exists; every SKILL.md opens with "Read `curator/curator-rules.md` first"
- curator/README.md:15-16 — curator/.claude-plugin/plugin.json exists and there is no marketplace file (Claude Code's own discovery was not tested)
- curator/README.md:22 — `literature search` records `hits`/`new_hits` and adds every unseen PMID to triage (gene_literature_io.py:213-263)
- curator/README.md:25-26 — `add-publication` appends to publications.yaml (cli.py:298-338)
- curator/README.md:27-29 — `diff --against` compares; `--mark-verified` sets verified_by on agreeing, unverified findings (cli.py:358-410)
- curator/README.md:30-34 — the `review` listing shows only unapproved findings, with every field named, plus study_type and own_lab (cli.py:451-486, 555-580)
- curator/README.md:33-34 — `--against` prints differing fields and "second pass: no partner finding" (cli.py:569-580)
- curator/README.md:35-38 — `--approve` refuses an `agent:` name, a date before the first search, and a future date (cli.py:535-542, 614-619)
- curator/README.md:38-40 — an unverified finding needs `--allow-unverified` (cli.py:621-625)
- curator/README.md:40-41 — the digest excludes approval fields, verified_by and the review mark (models/gene_literature.py:70-79); LIT015 (validate/gene_literature.py:212-218); re-approval needs approved_by cleared (cli.py:608-612)
- curator/README.md:43 — NCBI_API_KEY is read from the environment only and redacted from errors (cli.py:168-174)
- curator/README.md:44 — search, abstracts and add-publication call PubMed; diff and review do not
- curator/README.md:46-47 — no validator reads curator/
- curator/curator-rules.md:14-15 — locators min_length=1; count_unit required for the four rare-variant designs (models/gene_literature.py:49-56, 171, 184-187)
- curator/curator-rules.md:17-18 — check_review refuses a status or a note on its own (models/gene_literature.py:240)
- curator/curator-rules.md:27 — LIT004 (validate/gene_literature.py:185-190)
- curator/curator-rules.md:28-29 — the query ends `NOT preprint[pt] NOT review[pt]` (pubmed.py:116); both exclusion reasons exist
- curator/curator-rules.md:34-35 — LIT009; `search` resolves an alias or previous symbol and writes under the approved one (cli.py:182-205, 243)
- curator/curator-rules.md:36-39 — the id pattern `^CHDA:LIT:\d{7}$`; `prefixes: LIT: 0`; ID002 and ID003 (validate/ids.py:162-210); a repeat within one file fails SCHEMA001
- curator/curator-rules.md:40-42 — YAML002 rejects the whole file (corpus.py:233-243)
- curator/curator-rules.md:43-46 — LIT011 (publication and burden_study), LIT012; LIT013 is a warning
- curator/curator-rules.md:50-53 — LIT005, LIT010, LIT001; LIT002 is a warning
- curator/curator-rules.md:56-58, 61 — the reason/same_as/note rules (models/gene_literature.py:116-130); LIT006; LIT003
- curator/curator-rules.md:67-72 — carriers ≤ screened; case_report screens 1; comparator not none; common-variant cases and controls with no carriers; genotyping methods; no inf/NaN (models/gene_literature.py:159-226)
- curator/curator-rules.md:79-80 — LIT014
- curator/curator-rules.md:84-89 — shares_probands_with required; LIT007 (any existing id, approved or not); diff excludes it (gene_literature_diff.py:35-47)
- curator/curator-rules.md:95-105 — agent: refused by model and command; three approval fields together; LIT015 makes the build refuse (runner.py:208-210); only the prefix is checked (cli.py:354-355); LIT008
- curator/curator-rules.md:106-107 — published findings carry review_status/review_note (build/gene_literature.py:147)
- curator/skills/literature-search/SKILL.md:13-17 — alias and previous symbols accepted; with a key the throttle gap is 1/10 s, else 1/3 s (pubmed.py:301)
- curator/skills/literature-search/SKILL.md:27-35 — query template (cli.py:232-241, pubmed.py:103-116); output line (cli.py:257); path; a re-run adds only new PMIDs
- curator/skills/literature-search/SKILL.md:43-44, 59-62 — one LIT002 per file; a gene outside the gate exits 2 (cli.py:223-227); exit 1 before any write (cli.py:249-256)
- curator/skills/literature-triage/SKILL.md:24-27 — `--decision` defaults to pending (cli.py:263); no command writes a decision
- curator/skills/literature-triage/SKILL.md:32-40 — the reason table lists exactly the 7 ExclusionReason values
- curator/skills/literature-triage/SKILL.md:58, 64-66 — LIT006; LIT002 and LIT003 behaviour
- curator/skills/literature-verify/SKILL.md:26-32 — diff excludes shares_probands_with; `local/` is gitignored; `parse_findings_file` reads a top-level `findings:` without registry checks
- curator/skills/literature-verify/SKILL.md:37-54 — the example YAML parses as one valid Finding (`parse_findings_file`)
- curator/skills/literature-verify/SKILL.md:56-67 — diff writes nothing and exits 0; the pairing key and ambiguity rule; summaries and locators not compared; `--mark-verified` touches only agreed, unverified findings (cli.py:399-401)
- curator/skills/literature-extract/SKILL.md:27-43 — the five designs; the counter location; rare-variant and common-variant requirements; comparator values; locator kinds (models/assertion.py:34-54)
- curator/skills/literature-extract/SKILL.md:66-70 — exome, genome and genotyping_array exist; every StudyType named exists
- curator/skills/literature-extract/SKILL.md:73-79 — the command shape; a blank `--reading` refused; the reading written as comments above `study_type:`; own_lab defaults to false
- curator/skills/literature-extract/SKILL.md:84-106 — the 2/4/6 indent matches `_round_trip_yaml`; the example YAML parses as one valid Finding
- curator/skills/literature-extract/SKILL.md:110-111 — LIT003 and LIT011

## Not covered

- docs/data-api.md:26-28 (how GitHub Pages serves `.json.gz`): host behaviour; the proxy returns 403 for chd-atlas.github.io from this session.
- Historical figures that one build cannot reproduce: docs/data-api.md:97-102, :230-232, :265-266, :553-554, :955-962, :987-990, :2359-2374. Today's counterpart is given where one exists (for 2.13: 48 of 119 genes, 61 previous symbols).
- Figures from source papers not in the repository: :1116 (14,364 synonymous variants), :1167-1169, :1173-1175, :1184-1188 (S8 Table 521 of 528), :1236-1244 (the 248-gene local-FDR and Fisher counts; only 41 of the 248 panel genes are mirrored), :2206, :2452-2456 (the ≥50-case criterion; only the count of seven was checked), :2084 (HGNC's two UniProt accessions, checked only against the converter's recorded measurement), and curated cohort descriptions quoting papers.
- :1565 "no heart had at least 10 cells there": applied upstream by the converter; the build is consistent with it.
- :2616 (GenCC's licence wording): external text.
- Literature rendering (:1619-1639, :1915-1928, and the gene-page Literature section beyond "Not yet curated."): dev has no literature data, so this was checked against the code only.
- The gene table's JavaScript (sorting, filtering) was read, not exercised.
- `literature search`, `abstracts` and `add-publication` were not run, because they need PubMed. Their behaviour was read from cli.py and pubmed.py.
- README.md at the repository root is not in this job's scope and was not checked.
