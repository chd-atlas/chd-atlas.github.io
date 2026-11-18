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

`shape_issues` is the other half (LIT022-LIT024): it needs no frame, only the finding
itself. LIT022 is the fields against `placement`, the `with` links and the carrier
shape; LIT023 is a finding's variant carriers against its own (a compound
heterozygote is one person, so a `with` group counts once; per-origin counts cover
`de_novo` and `inherited` only, a departure from spec 7.2's literal "per origin"
because an `origin_unknown` variant claims no origin and forcing one would invent it);
LIT024 is a variant's phenotypes against the finding's.

LIT025 compares the paper's `p.` with `hgvs_p` where both number residues on the
MANE protein, and fires when `hgvs_p` is absent while the paper's `p.` names a
residue. A name numbering no residue (`p.?`, `p.=`, `p.0`) has nothing to
compare on either side: an `hgvs_p` (owner ruling D6) or the paper's own `p.`
(D6'), whatever the other holds. And a frameshift naming no new residue matches any
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
    protein_key,
    same_protein_change,
)
from chd_atlas.mane import ManeTranscript, genome_position, genome_slice, transcript_position
from chd_atlas.models.gene_literature import Finding, Variant
from chd_atlas.validate.mane import reverse_complement
from chd_atlas.vocab import (
    FindingDesign,
    FindingVariantClass,
    Phase,
    VariantPlacement,
    VariantZygosity,
)

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


def _protein_issues(
    text: str, t: ManeTranscript, normalised: bool = False
) -> tuple[list[Issue], PVariant | None]:
    """LIT021/LIT016/LIT019 on `hgvs_p`, and what it reads as on the MANE protein.

    `normalised` is true for a `mane` variant, whose `hgvs_p` comes from a
    normaliser run and is fixed by re-running it; a `protein` or `unplaced`
    variant is never normalised, so its fix is writing `hgvs_p` by hand.
    """
    try:
        pv = parse_p(text)
    except HgvsUnsupported as error:
        return [("LIT021", f"hgvs_p: {error}")], None
    if pv.accession != t.protein:
        return [
            (
                "LIT016",
                f"hgvs_p {text} is not on {t.protein}, the gene's pinned MANE Select protein "
                + (
                    "(re-normalise it on the pinned transcript)"
                    if normalised
                    else f"(write hgvs_p on {t.protein}, or for an unplaced variant drop it; "
                    "a protein variant is never normalised)"
                ),
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

    A name that numbers no residue has nothing to compare, so it is exempt on
    either side: an `hgvs_p` naming none (owner ruling D6) and a paper's own `p.?`,
    `p.=` or `p.0` (D6'), whatever `hgvs_p` holds. An absent `hgvs_p` is not
    exempt when the paper names a residue.
    """
    reported = variant.reported
    if reported.p is None:
        return []
    mane = {_unversioned(t.transcript), _unversioned(t.protein)}
    if reported.transcript is not None and _unversioned(reported.transcript) not in mane:
        return []  # another transcript may number residues differently
    if not _names_a_residue(reported.p):
        return []  # the paper's own `p.?` has nothing to compare (owner ruling D6')
    if variant.hgvs_p is None:
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


def _absent(variant: Variant) -> Issue:
    fix = (
        "drop its hgvs_p"
        if variant.placement is VariantPlacement.UNPLACED
        else "list it unplaced, with no hgvs_c or grch38 and no hgvs_p"
    )
    return (
        "LIT016",
        "the gene has no MANE Select in the pinned frame, so no variant on it can be "
        f"placed or its residue checked; {fix}",
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
            return [_absent(variant)]
        return sorted(_protein_issues(variant.hgvs_p, frame_entry)[0])
    if placement not in (VariantPlacement.MANE, VariantPlacement.PROTEIN):
        return []
    if frame_entry is None:
        return [_absent(variant)]
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
        found, read_p = _protein_issues(
            variant.hgvs_p, t, variant.placement is VariantPlacement.MANE
        )
        issues += found
    issues += _lit025(variant, t, read_p)
    return sorted(issues)


# (code, variant id or None for a finding-level issue, message)
ShapeIssue = tuple[str, int | None, str]

_HGVS_FIELDS = ("hgvs_c", "hgvs_p", "grch38")
# Origins reconciled per origin. `origin_unknown` is deliberately absent: such a
# variant claims no origin, so it draws on the finding's total only. Reading spec
# 7.2 literally ("per origin") would force a curator to invent an origin for it,
# since a paper reporting "1 de novo, 1 inherited" overall rarely says which of
# two variants is which.
_ORIGINS = ("de_novo", "inherited")

# What each placement must hold and must leave null (kind rules are separate).
_MUST_HOLD: dict[VariantPlacement, tuple[str, ...]] = {
    VariantPlacement.MANE: ("hgvs_c", "grch38", "normalisation"),
    VariantPlacement.PROTEIN: ("hgvs_p",),
    VariantPlacement.UNPLACED: ("review_status",),
    VariantPlacement.STRUCTURAL: (),
    VariantPlacement.OUTSIDE_REGION: (),
}
_MUST_BE_NULL: dict[VariantPlacement, tuple[str, ...]] = {
    VariantPlacement.MANE: (),
    VariantPlacement.PROTEIN: ("hgvs_c", "grch38", "normalisation"),
    VariantPlacement.UNPLACED: ("hgvs_c", "grch38"),
    VariantPlacement.STRUCTURAL: _HGVS_FIELDS,
    VariantPlacement.OUTSIDE_REGION: _HGVS_FIELDS,
}
_KINDS: dict[VariantPlacement, frozenset[FindingVariantClass]] = {
    VariantPlacement.MANE: frozenset({FindingVariantClass.SNV_INDEL}),
    VariantPlacement.STRUCTURAL: frozenset(
        {FindingVariantClass.CNV, FindingVariantClass.CHROMOSOMAL}
    ),
    VariantPlacement.OUTSIDE_REGION: frozenset({FindingVariantClass.SNV_INDEL}),
}


def _components(variants: list[Variant]) -> list[list[Variant]]:
    """The `with` groups: connected components of the links that resolve, size >= 2.

    Links are read as undirected, so an asymmetric link still groups its two
    entries and a chain 1-2, 2-3 groups all three; LIT022 reports the asymmetry
    and the missing 1-3 link separately. Empty when two entries share an id: a
    link then names no one entry, so no group can be built (LIT022 reports the
    duplicate itself).
    """
    if len({v.id for v in variants}) != len(variants):
        return []
    by_id = {v.id: v for v in variants}
    parent = {v.id: v.id for v in variants}

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for v in variants:
        for other in v.with_:
            if other in by_id and other != v.id:
                parent[find(v.id)] = find(other)
    groups: dict[int, list[Variant]] = {}
    for v in variants:
        groups.setdefault(find(v.id), []).append(v)
    return [sorted(g, key=lambda v: v.id) for _, g in sorted(groups.items()) if len(g) > 1]


def _placement_issues(
    finding: Finding, v: Variant, normalisation_ids: set[int]
) -> list[ShapeIssue]:
    out: list[ShapeIssue] = []

    def bad(message: str) -> None:
        out.append(("LIT022", v.id, message))

    placed = v.placement.value
    for name in _MUST_HOLD[v.placement]:
        if getattr(v, name) is None:
            bad(f"a {placed} variant needs {name}")
    for name in _MUST_BE_NULL[v.placement]:
        if getattr(v, name) is not None:
            bad(f"a {placed} variant must leave {name} null")
    allowed = _KINDS.get(v.placement)
    if allowed is not None and v.kind not in allowed:
        names = " or ".join(sorted(k.value for k in allowed))
        bad(f"a {placed} variant is kind {names}, not {v.kind.value}")
    if v.placement is VariantPlacement.PROTEIN and v.hgvs_p is not None:
        try:
            residues = parse_p(v.hgvs_p).residues
        except HgvsUnsupported:
            residues = ((" ", 0),)  # unreadable: LIT021 says so, once
        if not residues:
            bad(f"hgvs_p {v.hgvs_p} names no residue, so it cannot be placed; list it unplaced")
    if v.kind not in finding.variant_classes:
        bad(f"kind {v.kind.value} is not among the finding's variant_classes")
    if v.normalisation is not None and v.normalisation not in normalisation_ids:
        bad(f"normalisation {v.normalisation} is not in this file's normalisations")
    if finding.design is FindingDesign.COMMON_VARIANT_ASSOCIATION:
        if v.carriers is not None:
            bad("carriers on a common_variant_association finding, which counts cases and controls")
    elif v.carriers is None:
        bad(f"carriers are required on a {finding.design.value} finding")
    return out


def _link_issues(variants: list[Variant]) -> list[ShapeIssue]:
    out: list[ShapeIssue] = []
    by_id = {v.id: v for v in variants}

    def bad(v: Variant, message: str) -> None:
        out.append(("LIT022", v.id, message))

    for v in variants:
        for other_id in v.with_:
            other = by_id.get(other_id)
            if other_id == v.id:
                bad(v, "with names the variant itself")
            elif other is None:
                bad(v, f"with names variant {other_id}, which this finding does not have")
            elif v.id not in other.with_:
                bad(v, f"with names variant {other_id}, which does not name {v.id} back")
        if v.phase is not None and not v.with_:
            bad(v, "phase is only for a variant with a with link")
        if v.phase is None and v.with_:
            bad(v, "a variant with a with link needs a phase (not_stated when the paper is silent)")
        if v.with_ and v.zygosity in (VariantZygosity.HOMOZYGOUS, VariantZygosity.HEMIZYGOUS):
            bad(v, f"a {v.zygosity.value} variant has no partner, so it cannot have a with link")
        elif v.with_ and v.zygosity is not VariantZygosity.HETEROZYGOUS:
            bad(v, "a variant with a with link must be heterozygous")
    for group in _components(variants):
        ids = ", ".join(str(v.id) for v in group)
        phases = {v.phase for v in group if v.with_}
        if len(phases) > 1:
            out.append(("LIT022", group[0].id, f"variants {ids} are linked but differ in phase"))
        if len(group) > 2 and Phase.TRANS in phases:
            out.append(
                (
                    "LIT022",
                    group[0].id,
                    f"variants {ids} are linked in trans, but a diploid genome cannot hold "
                    "more than two variants pairwise in trans",
                )
            )
        for i, a in enumerate(group):
            for b in group[i + 1 :]:
                if b.id not in a.with_ and a.id not in b.with_:
                    out.append(
                        (
                            "LIT022",
                            a.id,
                            f"variants {a.id} and {b.id} are in one group through other entries, "
                            "but with lists every other entry the same carriers carry; "
                            f"neither names the other",
                        )
                    )
        totals = {v.carriers.total for v in group if v.carriers is not None}
        if len(totals) > 1:
            out.append(
                ("LIT022", group[0].id, f"variants {ids} are linked but carry different totals")
            )
    return out


def _carrier_issues(finding: Finding) -> list[ShapeIssue]:
    """LIT023: variant carriers never exceed the finding's, a `with` group counted once.

    - Unlinked entries count per origin (`de_novo`, `inherited`) and in total.
    - A `with` group is the same people, so it adds one shared total, not one per
      entry. Its origins are per allele and may differ between its entries, so a
      group adds to an origin only the pigeonhole lower bound: of its N people,
      each entry names `m.o` carriers of origin o, so at least
      `sum(m.o) - (k - 1) * N` people have every one of their k alleles from o.
      A compound heterozygote with one de novo and one inherited allele adds 0 to
      both.
    - `origin_unknown` claims no origin, so it enters the total check only.

    Skipped when two entries share an id (no group can be built; LIT022 says so).
    """
    if finding.carriers is None:
        return []
    if len({v.id for v in finding.variants}) != len(finding.variants):
        return []
    groups = _components(finding.variants)
    grouped = {v.id for group in groups for v in group}
    unlinked = [v for v in finding.variants if v.id not in grouped and v.carriers is not None]
    out: list[ShapeIssue] = []
    for origin in _ORIGINS:
        named = sum(getattr(v.carriers, origin) for v in unlinked if v.carriers is not None)
        for group in groups:
            carried = [v.carriers for v in group if v.carriers is not None]
            if len(carried) == len(group):
                shared = max(c.total for c in carried)
                named += max(
                    0, sum(getattr(c, origin) for c in carried) - (len(group) - 1) * shared
                )
        held = getattr(finding.carriers, origin)
        if named > held:
            out.append(
                (
                    "LIT023",
                    None,
                    f"variants name at least {named} {origin} carrier(s), more than the "
                    f"finding's {held}",
                )
            )
    group_totals = [
        max((v.carriers.total for v in group if v.carriers is not None), default=0)
        for group in groups
    ]
    named_total = sum(v.carriers.total for v in unlinked if v.carriers is not None) + sum(
        group_totals
    )
    if named_total > finding.carriers.total:
        out.append(
            (
                "LIT023",
                None,
                f"variants name {named_total} carrier(s) (each with group counted once), "
                f"more than the finding's {finding.carriers.total}",
            )
        )
    return out


def _row_key(v: Variant) -> str | None:
    """The key the page groups a variant under (`build.gene_literature`), or None."""
    if v.placement is VariantPlacement.MANE and v.grch38 is not None:
        return "g:" + v.grch38
    if v.placement is VariantPlacement.PROTEIN and v.hgvs_p is not None:
        return "p:" + protein_key(v.hgvs_p)
    return None


def _context(v: Variant, by_id: dict[int, Variant]) -> tuple[frozenset[str], Phase | None]:
    """A variant's genotype context (ruling D16'): its partners' names, and the phase.

    A partner is named by its row key, or by `as_reported` when it has none; the
    partners' ids do not count, so the same pair written twice is one context.
    """
    names = {
        (_row_key(p) or p.as_reported) if (p := by_id.get(o)) is not None else f"#{o}"
        for o in v.with_
    }
    return frozenset(names), v.phase


def shape_issues(finding: Finding, normalisation_ids: set[int]) -> list[ShapeIssue]:
    """LIT022-LIT024 for one finding's variants, sorted; each is `(code, variant id, message)`.

    The variant id is None for an issue about the finding's variants as a whole
    (LIT023). Nothing here needs the MANE frame, so it runs whether or not the
    frame loaded.
    """
    variants = finding.variants
    out: list[ShapeIssue] = []
    for number in sorted({v.id for v in variants if sum(w.id == v.id for w in variants) > 1}):
        out.append(("LIT022", number, f"variant id {number} is used more than once"))
    seen: dict[str, list[Variant]] = {}
    by_id = {v.id: v for v in variants}
    group_of = {v.id: n for n, group in enumerate(_components(variants)) for v in group}
    for v in variants:
        key = _row_key(v)
        if key is not None:
            twin = next(
                (
                    w
                    for w in seen.get(key, [])
                    if w.id != v.id
                    and (
                        _context(w, by_id) == _context(v, by_id)
                        or (v.id in group_of and group_of.get(w.id) == group_of[v.id])
                    )
                ),
                None,
            )
            if twin is not None:
                out.append(
                    (
                        "LIT022",
                        v.id,
                        f"variant {v.id} repeats variant {twin.id}'s name ({key[2:]}) in the "
                        "same genotype context (or in the same with group, where one change "
                        "twice is in effect homozygous); name each variant once per finding "
                        "per context, with zygosity not_stated when its carriers differ in "
                        "zygosity",
                    )
                )
            seen.setdefault(key, []).append(v)
        if v.carriers is not None and v.carriers.total == 0:
            out.append(
                (
                    "LIT022",
                    v.id,
                    f"variant {v.id} names 0 carriers; a variant a finding names is "
                    "carried by someone, so drop it or count its carriers",
                )
            )
        out += _placement_issues(finding, v, normalisation_ids)
        for pheno in v.phenotypes:
            if pheno not in finding.phenotypes:
                out.append(
                    ("LIT024", v.id, f"phenotype {pheno} is not among the finding's phenotypes")
                )
    out += _link_issues(variants)
    out += _carrier_issues(finding)
    return sorted(out, key=lambda i: (i[0], -1 if i[1] is None else i[1], i[2]))
