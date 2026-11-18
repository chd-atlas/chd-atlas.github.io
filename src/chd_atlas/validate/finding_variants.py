"""LIT016-LIT021 and LIT025: a curated variant's names against the pinned MANE frame.

`reference_issues` reads one variant's `hgvs_c`, `hgvs_p` and `grch38` and checks
each against the gene's MANE Select transcript, protein and genome windows
(`mirrors/mane_*.tsv`, exons ±50 bp). It reports; it never raises.

What runs, by placement:

- `mane` and `protein`: every field present is checked. A gene the frame does not
  hold (MAN004, a warning) is LIT016 and nothing else: such a variant cannot be
  placed, so it is listed `unplaced`.
- `unplaced`: only the `hgvs_p` it keeps, for LIT016 and LIT019 (and LIT021 when
  it cannot be read, so an unreadable residue is never kept unchecked).
- `structural` and `outside_region`: nothing; they carry no names to check.

A field that cannot be read (LIT021) or names another accession (LIT016) is not
checked further, because its positions would be read on the wrong sequence.

LIT018 compares changes, not coordinates. The `c.` edit and the VCF edit are
each applied to the forward-strand window around them, and the two resulting
sequences must be equal; so `c.` named 3'-shifted (as HGVS requires) and a
left-aligned VCF key agree on an indel in a repeat. The key itself must also be
minimal and left-aligned: moving an indel one base left must not give the same
sequence.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Final

from chd_atlas.hgvs import (
    CPosition,
    CVariant,
    GenomeKey,
    HgvsUnsupported,
    parse_c,
    parse_grch38,
    parse_p,
    same_protein_change,
)
from chd_atlas.mane import ManeTranscript, genome_position, genome_slice, transcript_position
from chd_atlas.models.gene_literature import Variant
from chd_atlas.validate.mane import reverse_complement
from chd_atlas.vocab import VariantPlacement

Issue = tuple[str, str]

# Bases of reference either side of a change when two edits are compared. Any
# context covering both edits gives the same answer; this one is for reading.
_CONTEXT: Final = 25


class _Refused(Exception):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


@dataclass(frozen=True)
class _Change:
    """Forward-strand bases `[start, end)` replaced by `alt` (`start == end` inserts)."""

    start: int
    end: int
    alt: str


def _c_text(p: CPosition) -> str:
    star = "*" if p.region == "three_prime" else ""
    offset = f"{p.offset:+d}" if p.offset else ""
    return f"c.{star}{p.base}{offset}"


def _span_text(cv: CVariant) -> str:
    if cv.start == cv.end:
        return _c_text(cv.start)
    return f"{_c_text(cv.start)}_{_c_text(cv.end)[2:]}"


def _forward(t: ManeTranscript, bases: str) -> str:
    return reverse_complement(bases) if t.strand == "-" else bases


def _new_bases(cv: CVariant) -> str:
    if cv.alt is None:  # parse_c sets it on every sub, ins and delins
        raise _Refused("LIT021", f"a {cv.edit} that names no new bases")
    return cv.alt


def _place(text: str, t: ManeTranscript, p: CPosition) -> tuple[int, int]:
    """`(transcript position, genome position)` of `p`, or raise `_Refused`.

    The genome position is not yet known to lie in a window; the caller checks
    the whole span at once.
    """
    where = f"hgvs_c {text}: {_c_text(p)}"
    cds_length = t.cds_end - t.cds_start + 1
    if p.region == "cds" and p.base > cds_length:
        raise _Refused(
            "LIT020", f"{where} lies past the {cds_length}-base coding sequence of {t.transcript}"
        )
    tx = transcript_position(t, p)
    if tx < 1 or tx > len(t.sequence):
        end = "5'" if tx < 1 else "3'"
        raise _Refused(
            "LIT020",
            f"{where} lies past the {end} end of {t.transcript} ({len(t.sequence)} bases)",
        )
    try:
        return tx, genome_position(t, p)
    except ValueError as error:
        raise _Refused("LIT017", f"{where}: {error}") from None


def _c_change(text: str, cv: CVariant, t: ManeTranscript) -> tuple[list[Issue], _Change]:
    """Place a parsed `c.` on the frame: its LIT017 issues and its forward edit."""
    (tx1, g1), (tx2, g2) = _place(text, t, cv.start), _place(text, t, cv.end)
    in_order = g1 <= g2 if t.strand == "+" else g1 >= g2
    if not in_order:
        raise _Refused("LIT021", f"hgvs_c {text}: the range runs against the transcript")
    lo, hi = min(g1, g2), max(g1, g2)
    if cv.edit == "ins" and hi - lo != 1:
        raise _Refused("LIT021", f"hgvs_c {text}: an insertion lies between two adjacent positions")
    # An insertion's span is the two bases it lies between.
    forward_ref = genome_slice(t, lo, hi)
    if forward_ref is None:
        where = f"{t.chrom}:{lo}" if lo == hi else f"{t.chrom}:{lo}-{hi}"
        raise _Refused(
            "LIT020",
            f"hgvs_c {text} ({where}) lies beyond the pinned exons ±50 bp of {t.transcript}",
        )
    if cv.edit == "ins":
        return [], _Change(hi, hi, _forward(t, _new_bases(cv)))
    one_exon = cv.start.offset == 0 and cv.end.offset == 0 and hi - lo == tx2 - tx1
    if one_exon:
        reference, source = t.sequence[tx1 - 1 : tx2], t.transcript
    else:
        reference, source = _forward(t, forward_ref), f"GRCh38 (on {t.transcript}'s strand)"
    issues: list[Issue] = []
    if cv.ref is not None and cv.ref != reference:
        issues.append(
            (
                "LIT017",
                f"hgvs_c {text} states {cv.ref} at {_span_text(cv)} but {source} has "
                f"{reference} there",
            )
        )
    if cv.edit == "sub":
        return issues, _Change(lo, lo + 1, _forward(t, _new_bases(cv)))
    if cv.edit == "del":
        return issues, _Change(lo, hi + 1, "")
    if cv.edit == "dup":
        return issues, _Change(hi + 1, hi + 1, forward_ref)
    return issues, _Change(lo, hi + 1, _forward(t, _new_bases(cv)))


def _c_issues(text: str, t: ManeTranscript) -> tuple[list[Issue], _Change | None]:
    try:
        cv = parse_c(text)
    except HgvsUnsupported as error:
        return [("LIT021", f"hgvs_c: {error}")], None
    if cv.accession != t.transcript:
        return [
            (
                "LIT016",
                f"hgvs_c {text} is not on {t.transcript}, the gene's pinned MANE Select "
                "(re-normalise it on the pinned transcript)",
            )
        ], None
    try:
        return _c_change(text, cv, t)
    except _Refused as refused:
        return [(refused.code, refused.message)], None


def _same_sequence(t: ManeTranscript, a: _Change, b: _Change) -> bool | None:
    """Whether two edits give one sequence; None when no one window holds both."""
    first = min(a.start, b.start)
    last = max(a.end - 1, b.end - 1, first)
    window = next((e for e in t.exons if e.window_start <= first and last <= e.window_end), None)
    if window is None:
        return None
    low = max(window.window_start, first - _CONTEXT)
    high = min(window.window_end, last + _CONTEXT)
    context = window.window_sequence[low - window.window_start : high - window.window_start + 1]

    def apply(c: _Change) -> str:
        return context[: c.start - low] + c.alt + context[c.end - low :]

    return apply(a) == apply(b)


def _left_shifted(t: ManeTranscript, key: GenomeKey) -> GenomeKey | None:
    """The same indel written one base to the left, or None if the window ends there."""
    before = genome_slice(t, key.pos - 1, key.pos - 1)
    if before is None:
        return None
    if len(key.ref) > len(key.alt):
        ref = genome_slice(t, key.pos - 1, key.pos + len(key.ref) - 2)
        if ref is None:
            return None
        return GenomeKey(key.chrom, key.pos - 1, ref, before)
    # The inserted bases rotate by one: `A` + I after pos becomes `before` + A + I[:-1].
    return GenomeKey(key.chrom, key.pos - 1, before, before + key.alt[:-1])


def _key_change(key: GenomeKey) -> _Change:
    return _Change(key.pos, key.pos + len(key.ref), key.alt)


def _key_issues(text: str, t: ManeTranscript) -> tuple[list[Issue], _Change | None]:
    try:
        key = parse_grch38(text)
    except HgvsUnsupported as error:
        return [("LIT021", f"grch38: {error}")], None
    if key.chrom != t.chrom:
        return [
            (
                "LIT020",
                f"grch38 {text} is on chromosome {key.chrom}; {t.transcript} is on {t.chrom}",
            )
        ], None
    end = key.pos + len(key.ref) - 1
    genome = genome_slice(t, key.pos, end)
    if genome is None:
        return [
            (
                "LIT020",
                f"grch38 {text} lies beyond the pinned exons ±50 bp of {t.transcript}",
            )
        ], None
    issues: list[Issue] = []
    if key.ref != genome:
        issues.append(("LIT018", f"grch38 {text} states ref {key.ref} but GRCh38 has {genome}"))
    # A pure indel (one allele is the anchor alone) is left-aligned when moving it
    # one base left gives another sequence. For such a key that is the same test
    # as "ref and alt end in different bases", which is what every other key, and
    # an indel at a window's edge, is held to.
    pure_indel = key.ref[0] == key.alt[0] and min(len(key.ref), len(key.alt)) == 1
    shifted = _left_shifted(t, key) if pure_indel else None
    if shifted is not None:
        if _same_sequence(t, _key_change(key), _key_change(shifted)):
            issues.append(
                (
                    "LIT018",
                    f"grch38 {text} is not left-aligned: "
                    f"{shifted.chrom}:{shifted.pos}:{shifted.ref}:{shifted.alt} is the same change",
                )
            )
    elif key.ref[-1] == key.alt[-1]:
        issues.append(("LIT018", f"grch38 {text} is not minimal: ref and alt end in one base"))
    elif len(key.ref) > 1 and len(key.alt) > 1 and key.ref[0] == key.alt[0]:
        issues.append(
            (
                "LIT018",
                f"grch38 {text} is not minimal: ref and alt share a first base beyond the "
                "VCF anchor",
            )
        )
    return issues, _key_change(key)


def _protein_issues(text: str, t: ManeTranscript) -> tuple[list[Issue], bool]:
    """LIT021/LIT016/LIT019 on `hgvs_p`, and whether it was read on the MANE protein."""
    try:
        pv = parse_p(text)
    except HgvsUnsupported as error:
        return [("LIT021", f"hgvs_p: {error}")], False
    if pv.accession != t.protein:
        return [
            (
                "LIT016",
                f"hgvs_p {text} is not on {t.protein}, the gene's pinned MANE Select protein "
                "(re-normalise it on the pinned transcript)",
            )
        ], False
    issues: list[Issue] = []
    length = len(t.protein_sequence)
    for residue, position in pv.residues:
        if residue == "*" and position == length + 1:
            continue  # the stop codon, which the protein sequence does not spell
        if position > length:
            issues.append(
                (
                    "LIT019",
                    f"hgvs_p {text} names {residue}{position}, past the end of {t.protein} "
                    f"({length} residues)",
                )
            )
        elif t.protein_sequence[position - 1] != residue:
            issues.append(
                (
                    "LIT019",
                    f"hgvs_p {text} names {residue} at {position} but {t.protein} has "
                    f"{t.protein_sequence[position - 1]}",
                )
            )
    return issues, True


def _unversioned(accession: str) -> str:
    return accession.split(".", 1)[0]


def _lit025(variant: Variant, t: ManeTranscript, hgvs_p_read: bool) -> list[Issue]:
    reported = variant.reported
    if reported.p is None:
        return []
    mane = {_unversioned(t.transcript), _unversioned(t.protein)}
    if reported.transcript is not None and _unversioned(reported.transcript) not in mane:
        return []  # another transcript may number residues differently
    if variant.hgvs_p is None:
        return [("LIT025", f"the paper names {reported.p} but the variant has no hgvs_p")]
    if hgvs_p_read and not same_protein_change(reported.p, variant.hgvs_p):
        return [
            (
                "LIT025",
                f"the paper names {reported.p} but hgvs_p is {variant.hgvs_p}; a variant whose "
                "names disagree is unplaced",
            )
        ]
    return []


def _absent() -> Issue:
    return (
        "LIT016",
        "the gene has no MANE Select in the pinned frame, so no variant on it can be "
        "placed or its residue checked; list it unplaced",
    )


def reference_issues(variant: Variant, frame_entry: ManeTranscript | None) -> list[Issue]:
    """`(code, message)` pairs for one variant against its gene's frame entry, sorted.

    `frame_entry` is None when the frame holds no row for the gene.
    """
    placement = variant.placement
    if placement is VariantPlacement.UNPLACED:
        if variant.hgvs_p is None:
            return []
        if frame_entry is None:
            return [_absent()]
        return sorted(_protein_issues(variant.hgvs_p, frame_entry)[0])
    if placement not in (VariantPlacement.MANE, VariantPlacement.PROTEIN):
        return []
    if frame_entry is None:
        return [_absent()]
    t = frame_entry
    issues: list[Issue] = []
    c_change = key_change = None
    if variant.hgvs_c is not None:
        found, c_change = _c_issues(variant.hgvs_c, t)
        issues += found
    if variant.grch38 is not None:
        found, key_change = _key_issues(variant.grch38, t)
        issues += found
    if c_change is not None and key_change is not None:
        same = _same_sequence(t, c_change, key_change)
        if not same:
            where = "lie in different places" if same is None else "give different sequences"
            issues.append(
                (
                    "LIT018",
                    f"hgvs_c {variant.hgvs_c} and grch38 {variant.grch38} describe different "
                    f"changes: applied to GRCh38 they {where}",
                )
            )
    hgvs_p_read = False
    if variant.hgvs_p is not None:
        found, hgvs_p_read = _protein_issues(variant.hgvs_p, t)
        issues += found
    issues += _lit025(variant, t, hgvs_p_read)
    return sorted(issues)
