# tests/unit/test_literature_query.py
"""The recorded PubMed query: its term list and template, pinned as published."""

from __future__ import annotations

from chd_atlas.literature_query import CHD_TERMS, build_query


def test_chd_terms_are_pinned() -> None:
    """The spec's list with every "(s)" expanded singular-then-plural (ruling R6)."""
    assert CHD_TERMS == (
        "congenital heart",
        "heart defect",
        "heart defects",
        "cardiac defect",
        "cardiac defects",
        "heart malformation",
        "heart malformations",
        "cardiac malformation",
        "cardiac malformations",
        "septal defect",
        "septal defects",
        "tetralogy of Fallot",
        "transposition of the great arteries",
        "hypoplastic left heart",
        "heterotaxy",
        "conotruncal",
        "outflow tract",
        "atrioventricular septal defect",
        "coarctation",
    )


def test_build_query_is_the_template() -> None:
    query = build_query("TBX5", ("HOS",), ("Holt-Oram syndrome",))
    assert query == (
        '("TBX5"[tiab] OR "HOS"[tiab])'
        ' AND ("Heart Defects, Congenital"[mh] OR "congenital heart"[tiab]'
        ' OR "heart defect"[tiab] OR "heart defects"[tiab] OR "cardiac defect"[tiab]'
        ' OR "cardiac defects"[tiab] OR "heart malformation"[tiab]'
        ' OR "heart malformations"[tiab] OR "cardiac malformation"[tiab]'
        ' OR "cardiac malformations"[tiab] OR "septal defect"[tiab]'
        ' OR "septal defects"[tiab] OR "tetralogy of Fallot"[tiab]'
        ' OR "transposition of the great arteries"[tiab]'
        ' OR "hypoplastic left heart"[tiab] OR "heterotaxy"[tiab]'
        ' OR "conotruncal"[tiab] OR "outflow tract"[tiab]'
        ' OR "atrioventricular septal defect"[tiab] OR "coarctation"[tiab]'
        ' OR "Holt-Oram syndrome"[tiab])'
        " NOT preprint[pt] NOT review[pt]"
    )


def test_build_query_sorts_previous_symbols_and_diseases_and_drops_quotes() -> None:
    """Order of the arguments must not change the query; a quote cannot break it."""
    query = build_query("GATA4", ("ZZ1", "AA1"), ('beta "x" disease', "alpha disease"))
    assert query.startswith('("GATA4"[tiab] OR "AA1"[tiab] OR "ZZ1"[tiab]) AND (')
    assert query.endswith(
        ' OR "alpha disease"[tiab] OR "beta x disease"[tiab]) NOT preprint[pt] NOT review[pt]'
    )
    assert build_query("GATA4", (), ()).startswith('("GATA4"[tiab]) AND (')
