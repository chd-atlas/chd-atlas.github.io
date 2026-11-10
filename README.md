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
uv run chd-atlas ingest SOURCE_ID --package DIR
```

A curator runs this, never CI. It verifies the package, merges its rows into `mirrors/`,
pins it in `mirrors/sources.yaml` with its manifest at `mirrors/packages/SOURCE_ID.json`,
then validates. If anything fails, every file it touched is restored. Exit 0 applied,
1 refused (every file restored), 2 when `--root` or `--package` is not a directory or on
a filesystem error; if the restore itself fails, the "could not restore" message names
every path left changed.

The package is refused with `ING015` if any table cell begins with a double quote; `validate` reports
the same in a committed mirror as `TBL013` (the CSV reader would otherwise merge rows silently).
It is refused with `ING016` if its version equals the pinned one but its manifest bytes differ: a version names one set of bytes.
