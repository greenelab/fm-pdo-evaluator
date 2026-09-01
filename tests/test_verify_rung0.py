"""The executable verification entry point runs green in continuous integration (PROCESS §3).

The notebook docs/tasks/rung0-replicate-ceiling/verify.ipynb is committed without outputs so
the reviewer executes it and watches the checks pass themselves; this test runs the same checks
so the branch's green does not depend on anyone opening the notebook. Every check recomputes a
promoted or reported number from the committed artifacts alone -- a failure here means a
document and an artifact disagree, the drift class the audits used to catch by reading.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytestmark = pytest.mark.step_document

_SPEC = importlib.util.spec_from_file_location(
    "verify_rung0", Path(__file__).resolve().parents[1] / "scripts" / "verify_rung0.py"
)
assert _SPEC is not None and _SPEC.loader is not None
vr = importlib.util.module_from_spec(_SPEC)
sys.modules["verify_rung0"] = vr
_SPEC.loader.exec_module(vr)


def test_every_promoted_claim_recomputes_from_the_committed_artifacts() -> None:
    checks = vr.run_all_checks()
    failures = [
        f"{c.name}: claim {c.claim!r} vs recomputed {c.computed!r}" for c in checks if not c.ok
    ]
    assert not failures, "documents and artifacts disagree:\n" + "\n".join(failures)


def test_the_check_count_matches_the_documents() -> None:
    # summary.md's trust map and the pull-request description say "60 checks"; pinning the
    # count here means adding or removing a check forces those transcriptions to be revisited
    # (scoped-review finding: the count was the one unguarded transcription in the wave).
    assert len(vr.run_all_checks()) == 60


def test_the_check_battery_covers_every_layer() -> None:
    # One check per layer at minimum: hashes, data pin, headline arithmetic, pool
    # arithmetic, permutation nulls, per-gene diagnostic, summary-vs-artifact agreement.
    names = " ".join(c.name for c in vr.run_all_checks())
    for fragment in (
        "raw per-condition values",
        "floor recomputes from its draws",
        "scatter reproduces its correlation",
        "checksum",
        "content hash",
        "Spearman-Brown",
        "scored-pair arithmetic",
        "permutation",
        "per-gene",
        "summary:",
    ):
        assert fragment in names, f"no check covers {fragment!r}"
