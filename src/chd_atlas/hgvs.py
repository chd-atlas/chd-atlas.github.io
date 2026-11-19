"""Read the HGVS subset a curated variant may carry, and its GRCh38 key.

Pure: no network, no sequence lookups. `validate/finding_variants.py` places
what these return on the pinned MANE frame.

The subset (spec 7.3, LIT021) is substitution, del, dup, ins and delins, at
exonic, UTR and intronic positions:

- a `c.` position is `N`, `-N` (5' UTR), `*N` (3' UTR), each optionally with an
  intronic offset `+M` or `-M`;
- a `c.` edit is `X>Y`, `del`, `del<bases>`, `dup`, `dup<bases>`,
  `ins<bases>` or `delins<bases>`, the bases uppercase `ACGT`;
- a `p.` description names its reference residue(s) as in `Arg279Ter`,
  `(Arg279Ter)`, `R279X`, `Arg279=`, `Met1?`, `Gln70dup`, `Lys5_Gly8del`,
  `Lys2_Gly3insArg`, `Lys5_Gly8delinsTrp`, a frameshift `Arg279GlyfsTer5` or a
  stop-loss `Ter519GlnextTer?` (whose reference residue is the stop);
  `p.?`, `p.=` and `p.0` (with or without parentheses) name no residue and read
  as none.

Everything else raises `HgvsUnsupported`: alleles (`;`, `[...]`), uncertain
positions (`(?_-1)`, `?`), `inv`, `con`, repeat and counted insertions, any
`ext` but a stop-loss, a stated deleted sequence on `delins`, a `c.`
substitution of a base by itself, and anything after the description, even a
newline. It is refused, never read as something nearby.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Final, Literal

Region = Literal["five_prime", "cds", "three_prime"]
Edit = Literal["sub", "del", "dup", "ins", "delins"]


class HgvsUnsupported(ValueError):
    """A description outside the subset this module reads."""


@dataclass(frozen=True)
class CPosition:
    """One `c.` position.

    `base` is the number as written, signed: `c.-45` is `("five_prime", -45, 0)`,
    `c.*1101` is `("three_prime", 1101, 0)`, `c.510+1` is `("cds", 510, 1)`.
    """

    region: Region
    base: int
    offset: int = 0


@dataclass(frozen=True)
class CVariant:
    """A `c.` description. `end` equals `start` for a single position.

    `ref` is the reference sequence the description states (the base of a
    substitution, the bases written after `del` or `dup`), or None where it
    states none. `alt` is the new base of a substitution or the inserted
    bases of `ins` and `delins`, else None. `accession` is None for a bare
    `c.` description.
    """

    accession: str | None
    start: CPosition
    end: CPosition
    edit: Edit
    ref: str | None
    alt: str | None


@dataclass(frozen=True)
class PVariant:
    """The reference residues a `p.` description names, one-letter, `*` for a stop."""

    accession: str | None
    residues: tuple[tuple[str, int], ...]


@dataclass(frozen=True)
class GenomeKey:
    """A GRCh38 VCF-style key, `chrom:pos:ref:alt`, forward strand, 1-based."""

    chrom: str
    pos: int
    ref: str
    alt: str


THREE_LETTER: Final[dict[str, str]] = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Sec": "U",
    "Ter": "*",
}
_THREE = "|".join(THREE_LETTER)
# Three-letter codes first, so `Arg` is never read as `A` followed by `rg`.
_AA = rf"(?:{_THREE}|[ACDEFGHIKLMNPQRSTVWYU*X])"
_STOP = r"(?:Ter|\*|X)"

_POS = r"([-*]?)(\d+)(?:([+-])(\d+))?"
# Every pattern is applied with `fullmatch`: `$` would also accept a trailing newline.
_C_RE: Final = re.compile(
    rf"(?:(?P<acc>[^:\s]+):)?c\.(?P<start>{_POS})(?:_(?P<end>{_POS}))?"
    r"(?:(?P<sub_ref>[ACGT])>(?P<sub_alt>[ACGT])"
    r"|delins(?P<delins>[ACGT]+)"
    r"|del(?P<del>[ACGT]*)"
    r"|dup(?P<dup>[ACGT]*)"
    r"|ins(?P<ins>[ACGT]+))"
)
_ONE_POS: Final = re.compile(_POS)

_P_RE: Final = re.compile(r"(?:(?P<acc>[^:\s]+):)?p\.(?P<body>\S+)")
_P_RESIDUES: Final = re.compile(
    rf"(?P<a1>{_AA})(?P<p1>\d+)(?:_(?P<a2>{_AA})(?P<p2>\d+))?(?P<edit>.*)"
)
_P_SINGLE_EDIT: Final = re.compile(
    rf"{_AA}|=|\?|del|dup|delins(?:{_AA})+|(?:{_AA})?fs(?:{_STOP}(?:\d+|\?))?"
    rf"|(?P<ext>{_AA}ext{_STOP}(?:\d+|\?))"
)
_P_RANGE_EDIT: Final = re.compile(rf"del|dup|delins(?:{_AA})+|ins(?:{_AA})+")
_P_NO_RESIDUE: Final = frozenset({"?", "=", "0", "0?"})

_KEY_RE: Final = re.compile(
    r"(?P<chrom>[1-9]|1\d|2[0-2]|X|Y|MT):(?P<pos>\d+):(?P<ref>[ACGT]+):(?P<alt>[ACGT]+)"
)


def _position(text: str) -> CPosition:
    match = _ONE_POS.fullmatch(text)
    if match is None:
        raise HgvsUnsupported(f"position {text!r} is not N, -N, *N or N+M/N-M")
    prefix, number, sign, offset_text = match.groups()
    base = int(number)
    offset = 0
    if offset_text is not None:
        offset = int(offset_text) if sign == "+" else -int(offset_text)
        if offset == 0:
            raise HgvsUnsupported(f"position {text!r} has an intronic offset of 0")
    if prefix == "-":
        if base == 0:
            raise HgvsUnsupported("there is no position c.-0")
        return CPosition("five_prime", -base, offset)
    if prefix == "*":
        if base == 0:
            raise HgvsUnsupported("there is no position c.*0")
        return CPosition("three_prime", base, offset)
    if base == 0:
        raise HgvsUnsupported("there is no position c.0")
    return CPosition("cds", base, offset)


def parse_c(text: str) -> CVariant:
    """Read one `c.` description in the subset, or raise `HgvsUnsupported`."""
    match = _C_RE.fullmatch(text)
    if match is None:
        raise HgvsUnsupported(f"{text!r} is not a c. description in the subset the atlas reads")
    start = _position(match["start"])
    end = _position(match["end"]) if match["end"] else start
    two = match["end"] is not None
    if match["sub_ref"] is not None:
        if two:
            raise HgvsUnsupported(f"{text!r}: a substitution names one position")
        if match["sub_ref"] == match["sub_alt"]:
            raise HgvsUnsupported(f"{text!r}: a substitution of a base by itself is no change")
        return CVariant(match["acc"], start, end, "sub", match["sub_ref"], match["sub_alt"])
    if match["delins"] is not None:
        return CVariant(match["acc"], start, end, "delins", None, match["delins"])
    if match["del"] is not None:
        return CVariant(match["acc"], start, end, "del", match["del"] or None, None)
    if match["dup"] is not None:
        return CVariant(match["acc"], start, end, "dup", match["dup"] or None, None)
    if not two:
        raise HgvsUnsupported(f"{text!r}: an insertion names the two positions it lies between")
    return CVariant(match["acc"], start, end, "ins", None, match["ins"])


def _one_letter(code: str) -> str:
    return "*" if code == "X" else THREE_LETTER.get(code, code)


def parse_p(text: str) -> PVariant:
    """Read one `p.` description in the subset, or raise `HgvsUnsupported`."""
    match = _P_RE.fullmatch(text)
    if match is None:
        raise HgvsUnsupported(f"{text!r} is not a p. description")
    body = match["body"]
    if body.startswith("(") and body.endswith(")"):
        body = body[1:-1]
    if body in _P_NO_RESIDUE:
        return PVariant(match["acc"], ())
    parts = _P_RESIDUES.fullmatch(body)
    if parts is None:
        raise HgvsUnsupported(f"{text!r} is not a p. description in the subset the atlas reads")
    first = (_one_letter(parts["a1"]), int(parts["p1"]))
    edit = parts["edit"]
    if parts["a2"] is None:
        single = _P_SINGLE_EDIT.fullmatch(edit)
        if single is None:
            raise HgvsUnsupported(f"{text!r}: edit {edit!r} is outside the subset")
        if single["ext"] is not None and first[0] != "*":
            raise HgvsUnsupported(f"{text!r}: only a stop-loss extension is read")
        residues: tuple[tuple[str, int], ...] = (first,)
    else:
        second = (_one_letter(parts["a2"]), int(parts["p2"]))
        if _P_RANGE_EDIT.fullmatch(edit) is None:
            raise HgvsUnsupported(f"{text!r}: edit {edit!r} on a range is outside the subset")
        if second[1] <= first[1]:
            raise HgvsUnsupported(f"{text!r}: a range must run from a lower position to a higher")
        if edit.startswith("ins") and second[1] != first[1] + 1:
            raise HgvsUnsupported(f"{text!r}: an insertion lies between two adjacent residues")
        residues = (first, second)
    if any(position == 0 for _, position in residues):
        raise HgvsUnsupported(f"{text!r}: there is no residue 0")
    return PVariant(match["acc"], residues)


def _protein_change(text: str) -> str:
    """The part after `p.`, parentheses gone, one-letter codes, `*` for every stop."""
    body = text.split("p.", 1)[1] if "p." in text else text
    body = body.replace("(", "").replace(")", "")
    body = re.sub(_THREE, lambda m: THREE_LETTER[m.group(0)], body)
    return body.replace("X", "*")


def protein_key(text: str) -> str:
    """The canonical reading of a `p.` description, for grouping.

    It is the reading `same_protein_change` compares.

    Accession and parentheses dropped, one-letter codes, `*` for every stop. It
    applies no frameshift equivalence (D5), so two frameshift spellings that
    compare equal there stay two keys here.
    """
    return _protein_change(text)


# A frameshift after normalising: reference residue, position, the new residue
# (empty when the description names none), then the rest.
_FRAMESHIFT: Final = re.compile(r"([A-Z*])(\d+)([A-Z*]?)fs(.*)")


# A substitution whose new residue repeats its reference residue (`V205V`).
_SAME_RESIDUE: Final = re.compile(r"([A-Z])(\d+)\1")


def _comparable(text: str) -> str:
    """`_protein_change`, with a same-residue substitution written `X#=`."""
    change = _protein_change(text)
    match = _SAME_RESIDUE.fullmatch(change)
    return f"{match[1]}{match[2]}=" if match else change


def same_protein_change(a: str, b: str) -> bool:
    """Whether two `p.` descriptions state the same change.

    Compares the parts after `p.` (any accession dropped) with parentheses
    stripped, three-letter codes mapped to one-letter, and `Ter`, `*` and `X`
    equal.

    A synonymous change written with its residue repeated (`V205V`,
    `Val205Val`) equals the same change written `Val205=`, the form the
    normaliser returns; a different new residue (`V205A`) stays different.
    `protein_key` applies no such reading.

    A frameshift naming no new residue (`R279fs`, `H445fsX136`) equals any
    frameshift at the same residue and position, because papers abbreviate it
    and count its stop differently from the normaliser (owner rulings D5, D5');
    two frameshifts that both name the new residue still compare in full.
    """
    first, second = _comparable(a), _comparable(b)
    one, two = _FRAMESHIFT.fullmatch(first), _FRAMESHIFT.fullmatch(second)
    if one and two and (not one[3] or not two[3]):
        return one.group(1, 2) == two.group(1, 2)
    return first == second


def parse_grch38(text: str) -> GenomeKey:
    """Read `chrom:pos:ref:alt` (chromosome without `chr`), or raise `HgvsUnsupported`."""
    match = _KEY_RE.fullmatch(text)
    if match is None:
        raise HgvsUnsupported(f"{text!r} is not chrom:pos:ref:alt on GRCh38")
    pos = int(match["pos"])
    if pos == 0:
        raise HgvsUnsupported(f"{text!r}: positions are 1-based")
    if match["ref"] == match["alt"]:
        raise HgvsUnsupported(f"{text!r}: ref and alt are the same")
    return GenomeKey(match["chrom"], pos, match["ref"], match["alt"])
