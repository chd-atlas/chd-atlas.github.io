# src/chd_atlas/orphanet_scope.py
"""The rule that admits a syndrome on Orphanet's phenotype frequencies.

Issue #8 part 1, decided 2026-10-24 (spec `2026-10-24-orphanet-scoped-syndromes`):
a syndrome where structural CHD is frequent but not part of the diagnosis is in
scope when a curator nominates it **and** Orphanet annotates a structural heart
phenotype for the disorder it reaches at `Frequent (79-30%)` or stronger. The
curator chooses which terms to put forward; Orphanet's data decides whether each
one qualifies, through a check that can refuse.

Shared by `validate/scope.py` (SCP006/SCP007) and `build/scope.py` (`scope.json`
and the gene-page note), so the check that admits a term and the sentence that
explains it are one reading of one rule, not two. `chd-atlas ingest` reads the
same pinned release through the same line reader (`pinned_hpo_release`) for the
wider branch it keeps of Orphanet product 4 (`ingest/filters.py`).

**Why the restriction to structural phenotypes is the load-bearing part.**
Measured 2026-10-24 over the 8,086 out-of-scope disease terms in the validity
mirrors: with every phenotype under HP:0001627 counted, the Frequent floor
admits 185 terms and 272 genes. With cardiomyopathy and the other
non-malformation branches removed, it admits 132 terms. The difference includes
Fabry disease and Hurler syndrome, which qualify through cardiomyopathy and
valve degeneration rather than a malformation. Even restricted, the rule is not
safe to apply automatically -- HPO's valve branch cannot tell a congenital
malformation from a degenerative one -- which is why it checks nominations
rather than making them.
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping
from dataclasses import dataclass
from functools import cached_property
from pathlib import Path
from typing import Final

from chd_atlas.tables import (
    FLAT_TABLES,
    MONDO_XREFS,
    ORPHANET_OMIM,
    ORPHANET_PHENOTYPES,
    read_table,
)

# The branch `chd-atlas ingest` keeps of Orphanet product 4 (spec F23): wider
# than the rule, so that a change to the structural subset needs no new
# package. The rule never reads it; `ingest.filters` does.
CARDIOVASCULAR_ROOT: Final = "HP:0001626"

# HP:0001627's EXACT synonyms include "Congenital heart defect", so Orphanet
# annotating the root itself (it does, for Rubinstein-Taybi) is HPO's own
# spelling of congenital heart defects and counts.
STRUCTURAL_ROOT: Final = "HP:0001627"

# Branches under the root that are cardiac but not malformations. Each is named
# with the reason, because this list is the one judgement in the rule that is
# not Orphanet's -- approved by the owner on 2026-10-24.
EXCLUDED_BRANCHES: Final[dict[str, str]] = {
    "HP:0001637": "Abnormal myocardium morphology: cardiomyopathy, not a malformation",
    "HP:0031331": "Abnormal cardiomyocyte morphology: cellular, not structural",
    "HP:0001640": "Cardiomegaly: a size, acquired as often as congenital",
    "HP:0100544": "Neoplasm of the heart: a tumour, not a malformation",
    "HP:6000847": "Cardiac mass: a tumour or thrombus, not a malformation",
    "HP:0001697": "Abnormal pericardium morphology: outside the heart's structure",
    "HP:0004306": "Abnormal endocardium morphology: the lining, mostly acquired",
    "HP:4000001": "Abnormal cardiac MRI finding: an imaging result, not a lesion",
}

# Malformations HPO files outside the heart branch, each taken with every
# subclass the pinned release files under it. Coarctation and a patent ductus
# since 2026-10-24; the other six by owner decision 2026-10-26 (R1), because the
# rule left classic lesions uncounted only for where HPO happens to file them.
# Measured on hp-2026-06-23 that day: none of the six was in the closure.
# Aortic valve stenosis can also be degenerative, as the kept valve branch's
# findings can (module docstring); it is counted for the same reason they are.
ADDED_BRANCHES: Final[dict[str, str]] = {
    "HP:0001680": "Coarctation of aorta: classic CHD, filed under vascular morphology",
    "HP:0001643": "Patent ductus arteriosus: likewise",
    "HP:0001642": "Pulmonic stenosis: filed under pulmonary valve physiology",
    "HP:0001650": "Aortic valve stenosis: filed under aortic valve physiology",
    "HP:0004381": "Supravalvular aortic stenosis: filed under aortic valve physiology",
    "HP:0001660": "Truncus arteriosus: filed under the great arteries",
    "HP:0011611": "Interrupted aortic arch: filed under aortic arch morphology",
    "HP:0010772": "Anomalous pulmonary venous return, total and partial: under the pulmonary veins",
}

FREQUENCY_RANK: Final[dict[str, int]] = {
    "Obligate (100%)": 5,
    "Very frequent (99-80%)": 4,
    "Frequent (79-30%)": 3,
    "Occasional (29-5%)": 2,
    "Very rare (<4-1%)": 1,
    "Excluded (0%)": 0,
}

# HPO's own band for "frequent", and the floor the owner chose over Occasional
# (5%), at which nearly any multi-system syndrome would clear for some finding
# and the check would refuse almost nothing.
FREQUENCY_FLOOR: Final = "Frequent (79-30%)"

# Product 1 relations that place an OMIM entry *within* an Orphanet disorder:
# exact, or Orphanet's code broader than the OMIM one. `NTBT` (Orphanet's code
# narrower) and `ND` (undecided) do not, and neither does a mapping Orphanet has
# not validated. Measured 2026-10-24: every one of the 29 Noonan, LEOPARD,
# cardiofaciocutaneous, Adams-Oliver and Noonan-like subtypes reaches its parent
# by `E` or `BTNT`, all `Validated`.
SUBTYPE_RELATIONS: Final = frozenset({"E", "BTNT"})
VALIDATED: Final = "Validated"


@dataclass(frozen=True)
class _Obo:
    """What `_read_obo` takes from an OBO file: each term's children and alt ids, every id."""

    children: dict[str, set[str]]
    alternates: dict[str, set[str]]
    known: frozenset[str]

    def closure(self, *roots: str) -> set[str]:
        """``roots`` and every term under one of them by `is_a`; alt ids not included."""
        seen, stack = set(roots), list(roots)
        while stack:
            for child in self.children.get(stack.pop(), ()):
                if child not in seen:
                    seen.add(child)
                    stack.append(child)
        return seen

    def with_alternates(self, terms: set[str]) -> frozenset[str]:
        return frozenset(terms | {alt for term in terms for alt in self.alternates.get(term, ())})


def _read_obo(obo: Path) -> _Obo:
    """The one line reader of an OBO file, shared by `structural_phenotypes` and `hpo_release`.

    A line reader of the pinned OBO rather than pronto. The build does not load
    pronto, which costs 1.6 s per load, and pronto resolves no alt ids
    (measured: `KeyError` on 2 of 3 sampled), while an Orphanet annotation may
    carry one. Only `[Term]` stanzas are read: a `[Typedef]`'s `is_a` would
    otherwise attach a relation to whichever term preceded it, and its `id`
    would be counted as a term the release knows.
    """
    children: dict[str, set[str]] = defaultdict(set)
    alternates: dict[str, set[str]] = defaultdict(set)
    known: set[str] = set()
    current: str | None = None
    in_term = False
    with obo.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("["):
                in_term = line == "[Term]"
                current = None
            elif not in_term:
                continue
            elif line.startswith("id: "):
                current = line[4:].strip()
                known.add(current)
            elif current is not None and line.startswith("is_a: "):
                children[line[6:].split(" ", 1)[0]].add(current)
            elif current is not None and line.startswith("alt_id: "):
                alternate = line[8:].strip()
                alternates[current].add(alternate)
                known.add(alternate)
    return _Obo(dict(children), dict(alternates), frozenset(known))


def structural_phenotypes(obo: Path) -> frozenset[str]:
    """`STRUCTURAL_ROOT`'s closure minus `EXCLUDED_BRANCHES`, plus `ADDED_BRANCHES`, with alt ids.

    Read by `_read_obo` (why a line reader, and why only `[Term]` stanzas).
    """
    release = _read_obo(obo)
    terms = (release.closure(STRUCTURAL_ROOT) - release.closure(*EXCLUDED_BRANCHES)) | (
        release.closure(*ADDED_BRANCHES)
    )
    return release.with_alternates(terms)


@dataclass(frozen=True)
class HpoRelease:
    """What `chd-atlas ingest` reads of an HPO release to filter Orphanet product 4."""

    # `CARDIOVASCULAR_ROOT`'s `is_a` closure plus every alt id of a term in it.
    cardiovascular: frozenset[str]
    # Every `id` and `alt_id` of a `[Term]`: an id outside it is one the
    # release cannot place, which is not the same as one outside the branch.
    known: frozenset[str]


def hpo_release(obo: Path) -> HpoRelease:
    """`HpoRelease` over one OBO file, read by `_read_obo`."""
    release = _read_obo(obo)
    return HpoRelease(
        cardiovascular=release.with_alternates(release.closure(CARDIOVASCULAR_ROOT)),
        known=release.known,
    )


def pinned_structural_phenotypes(root: Path, pinned: Mapping[str, str]) -> frozenset[str] | None:
    """`structural_phenotypes` over the HPO release `sources.yaml` pins, or None if unreadable.

    None rather than an empty set, for the reason `load_orphanet_tables` returns
    None: an empty set would make every Orphanet-scoped term look unwarranted.
    An unpinned or unreadable release is already reported (ONT003/ONT004), so
    SCP006/SCP008 name the consequence and not the cause. Shared by the validator
    and the build so both read the same release the same way.
    """
    relative = pinned.get("HP")
    if relative is None or not (root / relative).is_file():
        return None
    try:
        return structural_phenotypes(root / relative)
    except (OSError, UnicodeDecodeError):
        return None


def pinned_hpo_release(root: Path, pinned: Mapping[str, str]) -> HpoRelease | None:
    """`hpo_release` over the HPO release `sources.yaml` pins, or None if unreadable.

    None, never an empty release, as `pinned_structural_phenotypes` returns:
    an empty one would make every annotation's term unknown. Its caller
    (`ingest.filters`) names the release and refuses.
    """
    relative = pinned.get("HP")
    if relative is None or not (root / relative).is_file():
        return None
    try:
        return hpo_release(root / relative)
    except (OSError, UnicodeDecodeError):
        return None


@dataclass(frozen=True)
class Annotation:
    """One Orphanet phenotype annotation: the HPO term, its label, and its band."""

    hpo: str
    label: str
    frequency: str | None


@dataclass(frozen=True)
class OrphanetTables:
    """The three Orphanet scope mirrors, indexed the way the rule reads them."""

    # ORPHA code -> its cardiovascular annotations, in mirror order.
    annotations: dict[str, tuple[Annotation, ...]]
    # ORPHA code -> Orphanet's disorder name.
    names: dict[str, str]
    # (ORPHA code, OMIM entry) -> (relation, validation status).
    omim_links: dict[tuple[str, str], tuple[str, str]]
    # MONDO term -> its exact-match targets, `ORPHA:` and `OMIM:` alike.
    xrefs: dict[str, frozenset[str]]

    @cached_property
    def parents_of_omim(self) -> dict[str, tuple[str, ...]]:
        """OMIM entry -> every ORPHA code that contains it by a validated E/BTNT mapping."""
        parents: dict[str, set[str]] = defaultdict(set)
        for (orpha, omim), (relation, validation) in self.omim_links.items():
            if relation in SUBTYPE_RELATIONS and validation == VALIDATED:
                parents[omim].add(orpha)
        return {omim: tuple(sorted(codes)) for omim, codes in parents.items()}


@dataclass(frozen=True)
class Warrant:
    """Why a nominated term is in scope: the disorder it reaches, how, and what Orphanet records."""

    orpha: str
    name: str
    via: str
    omim: str | None
    phenotypes: tuple[Annotation, ...]


def _rank(annotation: Annotation) -> int:
    return FREQUENCY_RANK.get(annotation.frequency or "", -1)


def qualifying(
    orpha: str, tables: OrphanetTables, structural: frozenset[str]
) -> tuple[Annotation, ...]:
    """Structural annotations at or above the floor, strongest first, then by HPO id."""
    floor = FREQUENCY_RANK[FREQUENCY_FLOOR]
    hits = [
        annotation
        for annotation in tables.annotations.get(orpha, ())
        if annotation.hpo in structural and _rank(annotation) >= floor
    ]
    return tuple(sorted(hits, key=lambda annotation: (-_rank(annotation), annotation.hpo)))


def reach(mondo: str, orpha: str, tables: OrphanetTables) -> tuple[str, str | None] | None:
    """How `mondo` reaches `orpha`: directly, or as a validated subtype. None if it does not.

    Direct is checked first: a term MONDO maps straight to the disorder needs no
    OMIM step. The OMIM candidates are tried in sorted order so the route a
    warrant names is the same on every build.
    """
    targets = tables.xrefs.get(mondo, frozenset())
    if orpha in targets:
        return ("direct", None)
    for omim in sorted(target for target in targets if target.startswith("OMIM:")):
        if orpha in tables.parents_of_omim.get(omim, ()):
            return ("subtype", omim)
    return None


def warrant(
    mondo: str, orpha: str, tables: OrphanetTables, structural: frozenset[str]
) -> Warrant | str:
    """The warrant for nominating `mondo` on `orpha`, or the reason there is none.

    The reason names which half of the rule failed, because SCP006 prints it and
    the two are fixed differently: an unreachable code is a wrong `orphanet`
    field, a non-qualifying disorder is a syndrome the rule does not admit.
    """
    route = reach(mondo, orpha, tables)
    if route is None:
        return (
            f"{mondo} does not reach {orpha} by an exact MONDO match or a validated "
            f"E/BTNT Orphanet-to-OMIM mapping"
        )
    hits = qualifying(orpha, tables, structural)
    if not hits:
        return f"{orpha} carries no structural heart phenotype at {FREQUENCY_FLOOR} or stronger"
    via, omim = route
    return Warrant(orpha, tables.names.get(orpha, orpha), via, omim, hits)


def qualifying_targets(mondo: str, tables: OrphanetTables, structural: frozenset[str]) -> list[str]:
    """Every ORPHA code `mondo` reaches, by either route, that qualifies. Sorted. For SCP007."""
    targets = tables.xrefs.get(mondo, frozenset())
    reachable = {target for target in targets if target.startswith("ORPHA:")}
    for target in targets:
        if target.startswith("OMIM:"):
            reachable.update(tables.parents_of_omim.get(target, ()))
    return sorted(orpha for orpha in reachable if qualifying(orpha, tables, structural))


def load_orphanet_tables(root: Path) -> OrphanetTables | None:
    """Read the three mirrors. None if any is absent, unreadable or missing a column.

    None rather than empty tables, for the reason `REF000`/`SRC000` exist: an
    empty index would make every Orphanet-scoped term look unwarranted, burying
    the one missing file under one refusal per term.
    """
    frames = {}
    for schema in (ORPHANET_PHENOTYPES, ORPHANET_OMIM, MONDO_XREFS):
        path = root / "mirrors" / FLAT_TABLES[schema.name]
        if not path.is_file():
            return None
        frame, _ = read_table(path, schema)
        if frame is None or any(column not in frame.columns for column in schema.column_names):
            return None
        frames[schema.name] = frame

    annotations: dict[str, list[Annotation]] = defaultdict(list)
    names: dict[str, str] = {}
    for row in frames["orphanet_phenotypes"].iter_rows(named=True):
        annotations[row["orpha"]].append(Annotation(row["hpo"], row["hpo_label"], row["frequency"]))
        names[row["orpha"]] = row["disorder_name"]
    omim_links = {
        (row["orpha"], row["omim"]): (row["relation"], row["validation"])
        for row in frames["orphanet_omim"].iter_rows(named=True)
    }
    xrefs: dict[str, set[str]] = defaultdict(set)
    for row in frames["mondo_xrefs"].iter_rows(named=True):
        xrefs[row["mondo"]].add(row["target"])
    return OrphanetTables(
        annotations={orpha: tuple(rows) for orpha, rows in annotations.items()},
        names=names,
        omim_links=omim_links,
        xrefs={mondo: frozenset(targets) for mondo, targets in xrefs.items()},
    )
