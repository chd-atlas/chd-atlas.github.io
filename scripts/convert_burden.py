# scripts/convert_burden.py
"""Convert the Audain 2026 case-control supplement into `mirrors/burden.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6. The
output is committed, so a refresh is a reviewable diff rather than a silent
change under a build.

    .venv/bin/python scripts/convert_burden.py \\
        local/data/audain2026/41525_2026_582_MOESM4_ESM.xlsx mirrors/burden.tsv

The workbook is not committed: `local/` is gitignored because these are journal
supplementary files whose redistribution terms are the article's rather than the
data's, and because the reviewable artifact is the extracted TSV. Provenance
travels as the URL and the sha256 this script prints, recorded in
`mirrors/sources.yaml`.

WHAT WAS MEASURED, 2026-08-04, against the real workbook
--------------------------------------------------------
Supplementary Data 3 has four sheets: an Info sheet, then `aCHD_vs_controls`
(46,342 rows), `sCHD_vs_controls` (46,081) and `nsCHD_vs_controls` (46,186).
Each carries all four carrier counts, and the three differ only in which
comparison the `fet.*` columns describe. That mapping was verified
arithmetically rather than assumed: recomputing the odds ratio as
(a/(N-a))/(c/(M-c)) agreed with the published value on 32,464 / 22,992 / 27,788
testable rows for the matching stratum with **zero** disagreements. Every
mismatched pairing disagreed on tens of thousands of rows. Only the zero is
quoted as a figure: it is stable at every tolerance from 1e-2 to 5e-2, whereas
the disagreement counts move with the tolerance, and an earlier draft of this
docstring quoted a range no single tolerance reproduces.

`maf` is `1E-3` and `analysis` is the sheet's own stratum, both constant within
a sheet. Denominators are constant too -- 3,876 cases (1,471 syndromic + 2,405
non-syndromic) against 45,082 controls -- and do **not** match the abstract's
4,747 / 52,881. The paper does not reconcile the gap: it says it "assembled an
exome dataset consisting of 57,628 samples (4747 CHD cases and 52,881
controls)" and never states 3,876, 1,471 or 2,405 anywhere. So the reason for
the difference is unstated, and this docstring does not invent one -- an earlier
draft called 4,747 "the pre-QC recruitment figure", which the paper does not
say. What is certain is which numbers the published odds ratios were computed
from, and those are the ones the mirror carries.

`Infinity` is the only non-numeric value anywhere in the four statistic columns:
448 + 187 + 292 = 927 rows, always in `fet.odds_ratio` and `fet.ci_95_upper`
together, where no control carried. Those rows hold the study's strongest
results, so they are published with a null effect and `effect_bound =
unbounded_above` rather than dropped or given an invented ceiling. Zero rows
have an inverted confidence interval.

EXCEL HAS CORRUPTED 25 GENE SYMBOLS IN THIS WORKBOOK
-----------------------------------------------------
`SYMBOL` holds 25 values matching `\\d{4}-\\d{2}-\\d{2}` -- `2001-03-01`,
`2001-09-01`, `2001-12-01` and so on. These are the classic MARCH/SEPT/DEC
autocorrection: a spreadsheet read `MARCH1` as a date. It is irreversible from
this file, because `MARCH1` and `MARC1` collapse to the same cell value, which
is why four (SYMBOL, csq_group) pairs are duplicated in each sheet.

None of the 25 resolves to an in-scope gene, so nothing this atlas publishes is
lost today. `_report` names them anyway, and loudly: if a future workbook
corrupts a symbol that *is* in scope, that gene's entire burden evidence would
vanish from the mirror with a green build and no error -- this project's
characteristic failure exactly. A silent drop is the thing to prevent; a named
one is a curator's decision.
"""

from __future__ import annotations

import hashlib
import re
import sys
import zipfile
from pathlib import Path
from typing import Final
from xml.etree import ElementTree as ET

from burden_mirror import merge_study

# The publication, the source registry entry, and the cohorts -- all four are
# `curation/`- or `mirrors/`-side identifiers this script must agree with, so
# they are named once here rather than spelled inline per row.
STUDY: Final = "PMID:42230622"
SOURCE: Final = "audain2026_sd3"
CASE_COHORTS: Final = "cnchd;ddd;nottingham"
CONTROL_COHORTS: Final = "ukbb"

# sheet name -> (cohort_stratum, carrier column, denominator column). The
# carrier/denominator pair is what the arithmetic check above pinned.
_SHEETS: Final[dict[str, tuple[str, str, str]]] = {
    "aCHD_vs_controls": ("all", "n_het_cases", "n_total_cases"),
    "sCHD_vs_controls": ("syndromic", "n_het_syndromic", "n_total_syndromic"),
    "nsCHD_vs_controls": ("nonsyndromic", "n_het_nonsyndromic", "n_total_nonsyndromic"),
}

# The paper's `csq_group` vocabulary, mapped onto `ConsequenceClass`. `syn` is
# carried through rather than filtered: synonymous variants should show no
# enrichment, so the row is the study's own negative control and belongs beside
# the loss-of-function row it calibrates.
_CONSEQUENCE: Final[dict[str, str]] = {
    "hcLOF": "lof",
    "missC": "missense_damaging",
    "syn": "synonymous",
}

_DATE_CORRUPTED: Final = re.compile(r"\d{4}-\d{2}-\d{2}")
_CELL_REF: Final = re.compile(r"([A-Z]+)\d+")


# --------------------------------------------------------------------------
# A minimal read-only xlsx reader.
#
# `openpyxl` is not a dependency of this project and this script is the only
# thing that would need it, so the ~60 lines below are cheaper than adding one.
# The namespace is derived from the root tag rather than hardcoded: this
# workbook uses the *strict* OOXML namespace
# (`http://purl.oclc.org/ooxml/spreadsheetml/main`), and a reader assuming the
# transitional `schemas.openxmlformats.org` one silently finds zero rows.
# --------------------------------------------------------------------------
def _namespace(root: ET.Element) -> str:
    return root.tag.split("}")[0][1:] if root.tag.startswith("{") else ""


def _tag(namespace: str, name: str) -> str:
    return f"{{{namespace}}}{name}" if namespace else name


def _column_index(reference: str) -> int:
    """Zero-based column index of a cell reference such as "AB12"."""
    match = _CELL_REF.match(reference)
    if match is None:
        raise ValueError(f"unparseable cell reference {reference!r}")
    index = 0
    for character in match.group(1):
        index = index * 26 + (ord(character) - 64)
    return index - 1


def read_sheet(path: Path, name: str) -> list[dict[str, str]]:
    """One worksheet as a list of dicts keyed on its header row.

    Only correct for a sheet whose first row *is* the header. PNAS supplements
    put a title, a caption and a three-level header above the data, which keys
    every column on `""`; those callers use `read_sheet_rows` and index by
    position. The deleted Sierant converter (now the feeds recipe
    `sources/sierant2025_ds6`) checked the header row it indexed against rather
    than trusting it.
    """
    rows = read_sheet_rows(path, name)
    header = rows[0]
    return [
        {name: (row[i] if i < len(row) else "") for i, name in enumerate(header)}
        for row in rows[1:]
    ]


def read_sheet_rows(path: Path, name: str) -> list[list[str]]:
    """One worksheet as raw rows of cell text, blanks preserved by position.

    Position is the point: a sheet with a merged multi-level header has no row
    that can key a dict, and the columns still have to be addressable. Cells are
    placed by their `r` reference rather than by order of appearance, because a
    sparse row omits its empty cells entirely -- so a run of blanks in the middle
    of a row would otherwise shift every column after it left.
    """
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        workbook_ns = _namespace(workbook)
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {
            element.get("Id"): element.get("Target")
            for element in relationships.findall(_tag(_namespace(relationships), "Relationship"))
        }

        sheets = workbook.find(_tag(workbook_ns, "sheets"))
        if sheets is None:
            raise ValueError(f"{path} declares no sheets")
        target: str | None = None
        for sheet in sheets.findall(_tag(workbook_ns, "sheet")):
            if sheet.get("name") != name:
                continue
            # The r:id attribute is namespaced, and the namespace differs
            # between the strict and transitional OOXML flavours, so it is found
            # by suffix rather than spelled out.
            key = next(attribute for attribute in sheet.attrib if attribute.endswith("}id"))
            target = targets[sheet.get(key)]
        if target is None:
            raise ValueError(f"{path} has no sheet named {name!r}")

        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            table = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            table_ns = _namespace(table)
            shared = [
                "".join(node.text or "" for node in item.iter(_tag(table_ns, "t")))
                for item in table.findall(_tag(table_ns, "si"))
            ]

        worksheet = ET.fromstring(archive.read("xl/" + target.lstrip("/").removeprefix("xl/")))
        sheet_ns = _namespace(worksheet)
        data = worksheet.find(_tag(sheet_ns, "sheetData"))
        if data is None:
            raise ValueError(f"sheet {name!r} of {path} has no sheetData")

        rows: list[list[str]] = []
        for row in data.findall(_tag(sheet_ns, "row")):
            cells: dict[int, str] = {}
            for cell in row.findall(_tag(sheet_ns, "c")):
                reference = cell.get("r")
                index = _column_index(reference) if reference else len(cells)
                kind = cell.get("t")
                if kind == "inlineStr":
                    inline = cell.find(_tag(sheet_ns, "is"))
                    value = (
                        "".join(node.text or "" for node in inline.iter(_tag(sheet_ns, "t")))
                        if inline is not None
                        else ""
                    )
                else:
                    node = cell.find(_tag(sheet_ns, "v"))
                    value = (node.text or "") if node is not None else ""
                    if kind == "s" and value:
                        value = shared[int(value)]
                cells[index] = value
            rows.append([cells.get(i, "") for i in range(max(cells) + 1 if cells else 0)])

    return rows


def _number(value: str) -> str:
    """Round-trip a published float to its shortest exact representation.

    The workbook stores `0.61099999999999999`, which is the same IEEE-754 double
    as `0.611` and reads as false precision on a page. `repr(float(x))` is the
    shortest string that round-trips to the identical double, so this loses no
    information -- it is a rendering choice, not a rounding one.
    """
    return repr(float(value))


def convert(source: Path, symbols: dict[str, str]) -> tuple[list[dict[str, str]], list[str]]:
    """Build the mirror rows, and the report lines describing what was dropped.

    `symbols` maps a gene symbol to its HGNC id. A symbol absent from it is
    dropped -- the workbook covers 17,155 genes and this atlas registers 185 --
    and the report says how many, so a filter that silently matched nothing
    cannot look like a filter that worked.
    """
    rows: list[dict[str, str]] = []
    seen_symbols: set[str] = set()
    corrupted: set[str] = set()

    for sheet, (stratum, carriers, denominator) in _SHEETS.items():
        for record in read_sheet(source, sheet):
            symbol = record["SYMBOL"]
            seen_symbols.add(symbol)
            if _DATE_CORRUPTED.fullmatch(symbol):
                corrupted.add(symbol)
                continue
            hgnc = symbols.get(symbol)
            if hgnc is None:
                continue

            odds_ratio = record["fet.odds_ratio"]
            unbounded = odds_ratio == "Infinity"
            rows.append(
                {
                    "study": STUDY,
                    "gene": hgnc,
                    "cohort_stratum": stratum,
                    "lesion_group": "",
                    "variant_class": "snv_indel",
                    "consequence_class": _CONSEQUENCE[record["csq_group"]],
                    # The analysis collapsed rare variants regardless of how
                    # they were inherited; it is not a de novo test.
                    "origin": "any",
                    "maf_max": _number(record["maf"]),
                    # The supplement's `n_het_*` columns count individuals
                    # carrying at least one qualifying variant, which is why
                    # this study's denominators are its sample sizes rather than
                    # twice them.
                    "count_unit": "individuals",
                    "n_case_carriers": record[carriers],
                    "n_cases": record[denominator],
                    # Every case is a CHD case, so there is no split to give.
                    "n_chd_case_carriers": "",
                    "comparator": "control_cohort",
                    "n_control_carriers": record["n_het_controls"],
                    "n_controls": record["n_total_controls"],
                    "expected_count": "",
                    "effect": "" if unbounded else _number(odds_ratio),
                    "effect_measure": "odds_ratio",
                    "effect_bound": "unbounded_above" if unbounded else "",
                    "ci_low": _number(record["fet.ci_95_lower"]),
                    "ci_high": (
                        ""
                        if record["fet.ci_95_upper"] == "Infinity"
                        else _number(record["fet.ci_95_upper"])
                    ),
                    "pvalue": _number(record["fet.p_value"]),
                    "pvalue_test": "fisher_exact",
                    # Supplementary Data 3 publishes `fet.p_value` and no
                    # corrected column at all -- which is why the gene page
                    # names the study's 138,609 comparisons instead.
                    "pvalue_adjusted": "",
                    "pvalue_adjustment": "",
                    "case_cohorts": CASE_COHORTS,
                    "control_cohorts": CONTROL_COHORTS,
                    "method_note": "",
                    "source": SOURCE,
                }
            )

    return sorted(rows, key=_sort_key), _report(rows, symbols, seen_symbols, corrupted)


def _sort_key(row: dict[str, str]) -> tuple[object, ...]:
    """`TABLE_SCHEMAS["burden"].sort_key`, with `maf_max` compared as a number.

    `validate_sort_order` compares the values polars parsed, so `maf_max` is a
    float there and sorting it as a string here would disagree the moment a
    second frequency threshold appears -- "0.001" sorts after "0.0001" as text
    and before it as a number. Every other key column is a string in both.
    """
    return (
        row["gene"],
        row["study"],
        row["cohort_stratum"],
        row["lesion_group"],
        row["variant_class"],
        row["consequence_class"],
        row["origin"],
        float(row["maf_max"]),
    )


def _report(
    rows: list[dict[str, str]],
    symbols: dict[str, str],
    seen: set[str],
    corrupted: set[str],
) -> list[str]:
    """What was written, and -- the part that matters -- what was not."""
    kept = {row["gene"] for row in rows}
    missing = sorted(symbol for symbol, hgnc in symbols.items() if hgnc not in kept)
    lines = [
        f"{len(rows)} row(s) for {len(kept)} of the {len(symbols)} registered genes",
        f"{len(seen)} distinct symbols in the workbook; {len(seen) - len(kept) - len(corrupted)}"
        f" dropped as out of scope",
    ]
    if missing:
        lines.append(
            f"WARNING: {len(missing)} registered gene(s) have no row in this supplement, "
            f"so they will show no burden evidence: {', '.join(missing)}"
        )
    if corrupted:
        lines.append(
            f"WARNING: {len(corrupted)} symbol(s) were corrupted into dates by a spreadsheet "
            f"before publication and cannot be resolved to a gene: {', '.join(sorted(corrupted))}"
        )
        lines.append(
            "         None of them is a registered gene, so no published evidence is lost. "
            "If that ever changes, the affected gene's burden rows vanish silently -- "
            "recover the symbols from the authors, do not guess."
        )
    return lines


def gene_symbols(root: Path) -> dict[str, str]:
    """Symbol -> HGNC id, from `mirrors/genes.tsv`.

    Current symbols only. Aliases are deliberately not matched: measured against
    this workbook they resolve zero additional genes, while an alias that is
    another gene's current symbol would silently attribute one gene's burden to
    another -- a wrong claim traded for no coverage at all.
    """
    lines = (root / "mirrors" / "genes.tsv").read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    symbol, hgnc = header.index("symbol"), header.index("hgnc_id")
    return {line.split("\t")[symbol]: line.split("\t")[hgnc] for line in lines[1:]}


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, out = Path(argv[0]), Path(argv[1])
    root = Path(argv[2]) if len(argv) == 3 else Path.cwd()

    rows, report = convert(source, gene_symbols(root))
    kept, written = write_mirror(out, rows)

    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    for line in report:
        print(line)
    print(f"wrote {written} row(s) for {STUDY}, keeping {kept} row(s) from other studies")
    print(f"sha256 of {source.name}: {digest}")
    return 0


def write_mirror(out: Path, rows: list[dict[str, str]]) -> tuple[int, int]:
    """Replace this study's rows in `out`, keeping every other study's. Returns (kept, written).

    **This wrote the whole file until 2026-10-25**, long after the other two
    converters had moved to `burden_mirror.merge_study`. Measured that day, when
    the Orphanet scope work re-ran all three against a wider registry: run first,
    it left only its own 1,438 rows, and only the order the converters happened to
    be run in kept PMID:34324492's 121 and PMID:40127276's 246 from vanishing --
    on a green build, since nothing checks that a study's rows are still present.
    """
    return merge_study(out, STUDY, rows)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
