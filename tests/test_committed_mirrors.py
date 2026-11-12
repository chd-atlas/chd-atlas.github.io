"""Invariants of the committed population mirrors that no ingest filter states.

Since the ClinGen, GenCC and Orphanet converters left for chd-atlas/feeds, the
mirrors are written only by `chd-atlas ingest`, so nothing in the atlas can
show a converter defect; these tests read the committed files and the pinned
values, which is what catches a mirror edited by hand or a pin changed without
its release (CLAUDE.md section 4.23).
"""

from __future__ import annotations

from pathlib import Path
from typing import Final

from chd_atlas.orphanet_scope import load_orphanet_tables, reach
from chd_atlas.validate.sources import load_sources

REPO_ROOT = Path(__file__).resolve().parents[1]

_RECORDED_RELEASE: Final = {
    "orphanet_product4": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:57:18"),
    "orphanet_product1": ("1.3.42 / 4.1.8 [2025-03-03]", "2026-06-23 07:53:50"),
}


def test_the_committed_mirrors_carry_noonan_4s_route_to_its_parent() -> None:
    """MONDO:0012547 -> OMIM:610733 -> ORPHA:648 (BTNT, Validated), measured 2026-10-24."""
    tables = load_orphanet_tables(REPO_ROOT)
    assert tables is not None
    assert reach("MONDO:0012547", "ORPHA:648", tables) == ("subtype", "OMIM:610733")
    assert reach("MONDO:0018997", "ORPHA:648", tables) == ("direct", None)


def test_the_recorded_orphadata_release_is_pinned() -> None:
    """The `<JDBOR>` header values `sources.json` publishes and the landing page cites (R11).

    Transcribed from the headers of the two Orphadata files whose releases the
    feeds packages carry (2026-06-23 release, read 2026-10-26). They are
    recorded in `mirrors/sources.yaml`, not derived, so changing one without
    the other half of the pin is a published-value change this test names.
    """
    registry, _ = load_sources(REPO_ROOT)
    recorded = {
        source_id: (registry.by_id[source_id].file_version, registry.by_id[source_id].file_date)
        for source_id in ("orphanet_product4", "orphanet_product1")
    }
    assert recorded == _RECORDED_RELEASE
