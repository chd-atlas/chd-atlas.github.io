# chd-atlas

A curated set of genes associated with CHD.

**Development status.** This atlas is under active development and is **not a
clinical decision-support tool** — it must not be used to make or guide a
diagnostic, treatment or any other clinical decision. Gene-disease validity is
mirrored from ClinGen and GenCC, never authored by this atlas itself, and the
evidence curated here so far covers a small fraction of the genes the atlas
publishes on that mirrored validity. The published site's front page
(`index.html`) states the exact, build-derived numbers behind that gap.

## Published data API

The atlas is served as static JSON built from the curated sources — no server,
no query language. Build it locally with:

```bash
uv run chd-atlas build --root . --out dist
```

The command refuses to write anything if `chd-atlas validate` would report an
error, so a published site always passed its own gate. Two builds of one commit
are byte-identical, and `manifest.json` carries a sha256 for every file it
lists — every served file except itself, which cannot contain its own checksum.

`sources.json` records what the atlas mirrors and on whose terms — phenotype
labels come from the Human Phenotype Ontology, whose licence requires
attribution; gene-disease validity is mirrored from ClinGen and GenCC, which
this atlas attributes but does not curate itself. The `LICENSE` in this
repository covers the code rather than mirrored third-party content.

The shape of each file is documented in [docs/data-api.md](docs/data-api.md).
Consumers should read that before writing against the output — in particular the
note on contested genes, which must never be displayed as settled.

## Validating the corpus

```bash
uv run chd-atlas validate
```

Reports every problem it finds rather than stopping at the first, and exits 1 on
any error, 2 if `--root` is not a directory.

## Ingesting a source package

```bash
gh release download SOURCE_ID/VERSION -R chd-atlas/feeds -D local/packages/SOURCE_ID/VERSION
uv run chd-atlas ingest SOURCE_ID --package local/packages/SOURCE_ID/VERSION
uv run chd-atlas ingest SOURCE_ID   # re-ingest the pinned version from local/packages/
uv run chd-atlas ingest --all       # after any edit to curation/chd_scope.yaml
```

A curator runs this, never CI. The atlas never downloads: fetch a release from
`chd-atlas/feeds` with `gh release download` (`local/` is gitignored). Ingest verifies
the package, merges its rows into `mirrors/`, pins it in `mirrors/sources.yaml` with its
manifest at `mirrors/packages/SOURCE_ID.json`, then validates. If anything fails, every
file it touched is restored. Exit 0: applied, and the repository validates (for `--all`,
also when nothing is pinned). Exit 1: refused, every file restored. Exit 2: a bad
argument (`--root` or `--package` not a directory; a SOURCE_ID with no pin and no
`--package`; `--all` with a SOURCE_ID or `--package`; neither a SOURCE_ID nor `--all`) or
a filesystem error; if the restore itself fails, the "could not restore" message names
every path left changed.

Ingest in this order, which `--all` follows: ClinGen and GenCC (written whole, never
filtered: the registry is derived from them), then `hgnc` (it writes `mirrors/genes.tsv`,
filtered to that derived registry), then every other package (filtered to
`mirrors/genes.tsv`; its gene symbols resolve against the pinned HGNC package). Without
`--package`, the pinned version's cache is used; `ING017` refuses when that cache is
missing, does not hold the pinned manifest, or holds one naming another version than the
pin, and names the `gh release download` command for a missing or wrong cache.

After a scope shrink (an edit to `curation/chd_scope.yaml` that removes a term), run
`chd-atlas ingest --all`. It re-ingests every pinned package from its cache, `hgnc`
included, then removes every row naming a gene that left `mirrors/genes.tsv` from every
mirror filtered to the registry, whichever package or converter wrote it, and prints each
file, gene and row count. It validates once at the end and restores all of `mirrors/` if
anything fails. Ingesting `hgnc` alone fails validation whenever a registry-filtered
mirror still names a dropped gene (validity rows keep naming it and do not fail).

After a scope growth, run `chd-atlas ingest --all` too: `hgnc` and `sierant2025_ds6` are
every-gene packages, so re-ingesting them from their caches adds the new genes' rows. Two
things `--all` cannot do. The Farah package is a registry slice (cut at the registry its
manifest names), so a grown registry needs a new `farah2024_hoc` package from
`chd-atlas/feeds`, pinned and ingested; until then `ingest` names the registry genes the
slice lacks. And the converters still in `scripts/` that read the registry
(`convert_audain2021.py`, `convert_burden.py`, `convert_cardoso_moreira.py`) are re-run by
hand. `hgnc`, `sierant2025_ds6` and `farah2024_hoc` have no converter here any more: their
recipes live in `chd-atlas/feeds` (`sources/<id>`).

Two tables are kept only in part, because their packages carry the whole source.
`orphanet_phenotypes` (Orphanet product 4) keeps an annotation whose `hpo` is under
HP:0001626 (Abnormality of the cardiovascular system), alt ids included, in the HPO release
`mirrors/sources.yaml` pins for `HP`; an HPO id that release does not define is not counted
as outside the filter but named, and its rows are not written. `mondo_xrefs` (Mondo SSSOM)
keeps a row whose `mondo` term is a `disease` of either validity mirror as it stands at
ingest; `--all` ingests it after ClinGen and GenCC, so it is re-filtered against the
mirrors they just wrote. Ingest prints how many rows each filter left out, and refuses
with `ING010` when a filter cannot be built (the pinned HPO release cannot be read or lacks
HP:0001626, a validity mirror fails its schema, or neither names a disease).

The package is refused with `ING015` if any table cell's parsed value begins with a double quote; `validate` reports
the same in a committed mirror as `TBL013` (the CSV reader would otherwise merge rows silently).
It is refused with `ING016` if its version equals the pinned one but its manifest bytes differ: a version names one set of bytes.
Only the `hgnc` package ships the `genes` table, always keyed by `hgnc_id` (`ING009`); its ingest writes
`mirrors/genes.tsv` filtered to the registry the atlas derives from `curation/chd_scope.yaml` and the two validity mirrors.
`validate` reports a gene that registry holds and `mirrors/genes.tsv` lacks as `GEN006`; re-ingest `hgnc` or run
`chd-atlas ingest --all` to clear it (if the pinned HGNC package lacks the gene, pin a newer HGNC release).

## Curating literature

Searching PubMed, triaging hits, extracting findings and the blind second pass
live in [chd-atlas/curator](https://github.com/chd-atlas/curator) (`chd-curator`
and its Claude plugin). The searches it records are published as they are, but no
finding is published until a named person approves it here with
`chd-atlas literature review`. The rules it follows are in
`curator/curator-rules.md`.
