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
checked further, because its positions would be read on the wrong sequence. A
`c.N` past the coding sequence is LIT021 too: on this transcript that place is
written `c.*N`.

LIT020 means a position no pinned window holds: beyond an exon ±50 bp, or a key
on another chromosome. A `c.-N` or `c.*N` off the transcript's ends but inside
the first or last exon's window is placed and its base checked (owner ruling D7).

LIT018 compares changes, not coordinates. The `c.` edit and the VCF edit are
each applied to the forward-strand window around them, and the two resulting
sequences must be equal; so `c.` named 3'-shifted (as HGVS requires) and a
left-aligned VCF key agree on an indel in a repeat. The key itself must also be
minimal and left-aligned: ref and alt must not end in one base, nor share a
first base beyond the VCF anchor.

LIT025 compares the paper's `p.` with `hgvs_p` where both number residues on the
MANE protein, and fires when `hgvs_p` is absent while the paper's `p.` names a
residue. An `hgvs_p` naming no residue (`p.?`, `p.=`, `p.0`) has nothing to
compare (owner ruling D6), and a frameshift naming no new residue matches any
frameshift at its residue (D5, D5', in `same_protein_change`).
"""

from __future__ import annotations

from dataclasses import dataclass

from chd_atlas.hgvs import (
    CPosition,
    CVariant,
    GenomeKey,
    HgvsUnsupported,
    PVariant,
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
            "LIT021",
            f"{where} lies past the {cds_length}-base coding sequence of {t.transcript}, "
            "where a position is written c.*N",
        )
    tx = transcript_position(t, p)
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
    in_transcript = tx1 >= 1 and tx2 <= len(t.sequence)
    one_exon = in_transcript and cv.start.offset == cv.end.offset == 0 and hi - lo == tx2 - tx1
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
    """Whether two edits give one sequence; None when no one window holds both.

    Compared over the bases from the first either edit touches to the last:
    outside that span both sequences are the reference.
    """
    first = min(a.start, b.start)
    last = max(a.end - 1, b.end - 1, first)
    window = next((e for e in t.exons if e.window_start <= first and last <= e.window_end), None)
    if window is None:
        return None
    span = window.window_sequence[first - window.window_start : last - window.window_start + 1]

    def apply(c: _Change) -> str:
        return span[: c.start - first] + c.alt + span[c.end - first :]

    return apply(a) == apply(b)


def _left_aligned(t: ManeTranscript, key: GenomeKey) -> GenomeKey | None:
    """A pure indel moved left as far as it goes, or None if the pinned bases end first.

    While ref and alt end in one base the indel moves one base left: the base
    before it becomes the anchor and the last base of each allele drops off.
    """
    pos, ref, alt = key.pos, key.ref, key.alt
    while ref[-1] == alt[-1]:
        before = genome_slice(t, pos - 1, pos - 1)
        if before is None:
            return None
        pos, ref, alt = pos - 1, before + ref[:-1], before + alt[:-1]
    return GenomeKey(key.chrom, pos, ref, alt)


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
    # Minimal and left-aligned: ref and alt must not end in one base, nor share a
    # first base beyond the VCF anchor. For a pure indel (one allele is the anchor
    # alone), ending in one base is exactly "it can move one base left" (20,000 of
    # 20,000 random indels agreed), so the message names where it moves to.
    if key.ref[-1] == key.alt[-1]:
        pure_indel = key.ref[0] == key.alt[0] and min(len(key.ref), len(key.alt)) == 1
        aligned = _left_aligned(t, key) if pure_indel else None
        if aligned is not None:
            reason = (
                f"not left-aligned: {aligned.chrom}:{aligned.pos}:{aligned.ref}:{aligned.alt} "
                "is the same change"
            )
        else:
            reason = "not minimal: ref and alt end in one base"
        issues.append(("LIT018", f"grch38 {text} is {reason}"))
    elif len(key.ref) > 1 and len(key.alt) > 1 and key.ref[0] == key.alt[0]:
        issues.append(
            (
                "LIT018",
                f"grch38 {text} is not minimal: ref and alt share a first base beyond the "
                "VCF anchor",
            )
        )
    return issues, _key_change(key)


def _protein_issues(text: str, t: ManeTranscript) -> tuple[list[Issue], PVariant | None]:
    """LIT021/LIT016/LIT019 on `hgvs_p`, and what it reads as on the MANE protein."""
    try:
        pv = parse_p(text)
    except HgvsUnsupported as error:
        return [("LIT021", f"hgvs_p: {error}")], None
    if pv.accession != t.protein:
        return [
            (
                "LIT016",
                f"hgvs_p {text} is not on {t.protein}, the gene's pinned MANE Select protein "
                "(re-normalise it on the pinned transcript)",
            )
        ], None
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
    return issues, pv


def _unversioned(accession: str) -> str:
    return accession.split(".", 1)[0]


def _names_a_residue(text: str) -> bool:
    """False only for a `p.` that parses and names none (`p.?`, `p.=`, `p.0`)."""
    try:
        return bool(parse_p(text).residues)
    except HgvsUnsupported:
        return True


def _lit025(variant: Variant, t: ManeTranscript, hgvs_p: PVariant | None) -> list[Issue]:
    """The paper's own `p.` against `hgvs_p`, when both number residues on MANE.

    An `hgvs_p` naming no residue has nothing to compare, so it is exempt (owner
    ruling D6); an absent one is not, when the paper names a residue.
    """
    reported = variant.reported
    if reported.p is None:
        return []
    mane = {_unversioned(t.transcript), _unversioned(t.protein)}
    if reported.transcript is not None and _unversioned(reported.transcript) not in mane:
        return []  # another transcript may number residues differently
    if variant.hgvs_p is None:
        if not _names_a_residue(reported.p):
            return []
        return [("LIT025", f"the paper names {reported.p} but the variant has no hgvs_p")]
    if hgvs_p is None or not hgvs_p.residues:
        return []  # unreadable or on another protein (already reported), or no residue named
    if not same_protein_change(reported.p, variant.hgvs_p):
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
    read_p: PVariant | None = None
    if variant.hgvs_p is not None:
        found, read_p = _protein_issues(variant.hgvs_p, t)
        issues += found
    issues += _lit025(variant, t, read_p)
    return sorted(issues)
