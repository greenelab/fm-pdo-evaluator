"""The executable verification entry point runs green in continuous integration (PROCESS section 3).

The notebook `docs/tasks/rung0-assay-reliability/verify.ipynb` is committed without outputs so a
reviewer executes it and watches the checks pass themselves; these tests run the same battery so
the branch's green does not depend on anyone opening the notebook. Every check recomputes a
reported number from the run's own artifacts, so a failure here means a document and an artifact
disagree -- the drift class the audits used to catch by reading.

The last test is the design's `document` step negative control. A battery that reports PASS on
evidence somebody altered is measuring nothing, so one number in the summary row is changed and
the real battery is required to fail, naming the claim that moved.
"""

# pandas ships no PEP-561 type stubs in this environment; under strict mode that turns every
# pandas call site into a cascade of reportUnknown* noise about *pandas'* types, not ours. Same
# suppression, same rationale as the rest of this project's pyright strict config where it
# touches scientific-Python packages -- the rules that check our own code stay on.
# pyright: reportMissingTypeStubs=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownParameterType=false, reportUnknownLambdaType=false

from __future__ import annotations

import csv
import importlib.util
import os
import shutil
import sys
from pathlib import Path
from typing import Any

import pytest

pytestmark = pytest.mark.step_document

REPO = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location("verify_rung0", REPO / "scripts" / "verify_rung0.py")
assert _SPEC is not None and _SPEC.loader is not None
vr: Any = importlib.util.module_from_spec(_SPEC)
sys.modules["verify_rung0"] = vr
_SPEC.loader.exec_module(vr)

#: The battery is exercised against a real directory of run artifacts. In order: whatever
#: RUNG0_EXAMPLE_RUN points at, the task's own folder once the cluster run has been pulled into
#: it, then the synthetic-data example run kept in this session's scratch. Artifacts are
#: uncommitted between the run and promotion (PROCESS, "What reaches GitHub, and when"), so on a
#: checkout where none of the three exists these tests skip with the reason rather than fail.
_SCRATCH = Path(
    "/private/tmp/claude-502/-Users-gillenlu-Repositories-fm-pdo-evaluator"
    "/ba328f82-4da7-4056-b76c-dc939ce15d8a/scratchpad/synthrun/out"
)
_CANDIDATES = [
    Path(os.environ["RUNG0_EXAMPLE_RUN"]) if os.environ.get("RUNG0_EXAMPLE_RUN") else None,
    REPO / "docs" / "tasks" / "rung0-assay-reliability",
    _SCRATCH,
]

#: Both the summary table's trust map and this task's pull-request description quote the size of
#: the battery. Pinning it here means adding or removing a check forces those transcriptions to be
#: revisited rather than quietly going stale -- the one number in the wave that nothing else
#: guards. Skipped checks are counted: whether the permutation job has run changes what can be
#: checked, never how many claims the battery is answerable for.
EXPECTED_CHECKS = 51


def _example_run() -> Path:
    for candidate in _CANDIDATES:
        if candidate is not None and (candidate / "rung0_reliability.csv").exists():
            return candidate
    pytest.skip(
        "no rung 0 run to verify: none of "
        f"{[str(c) for c in _CANDIDATES if c is not None]} holds rung0_reliability.csv"
    )


@pytest.fixture
def artifacts(tmp_path: Path) -> Path:
    """A private copy of one run's artifacts, so a test may perturb them without harm."""
    source = _example_run()
    destination = tmp_path / "run"
    shutil.copytree(source, destination)
    return destination


def _failures(checks: list[Any]) -> list[str]:
    return [f"{c.name}: claim {c.claim!r} vs recomputed {c.computed!r}" for c in checks if not c.ok]


def test_every_promoted_claim_recomputes_from_the_committed_artifacts(artifacts: Path) -> None:
    checks = vr.run_all_checks(artifacts)
    failures = _failures(checks)
    assert not failures, "documents and artifacts disagree:\n" + "\n".join(failures)


def test_the_check_count_matches_the_documents(artifacts: Path) -> None:
    checks = vr.run_all_checks(artifacts)
    assert len(checks) == EXPECTED_CHECKS, (
        f"the battery runs {len(checks)} checks, the documents say {EXPECTED_CHECKS}"
    )


def test_the_check_battery_covers_every_layer(artifacts: Path) -> None:
    # At least one check per layer the design names: both gene sets' summary statistics, the
    # Spearman-Brown correction, the chance floors, detection power, the effect-size control, the
    # selection-leakage control, the noise decomposition, the example scatters, the figures, the
    # pool arithmetic, the recorded checksums, the data pin, the permutation check and promotion.
    names = " ".join(str(c.name) for c in vr.run_all_checks(artifacts))
    for fragment in (
        "all: mean recomputes",
        "responder: mean recomputes",
        "Spearman-Brown",
        "floor recomputes from its draws",
        "minimum detectable effects",
        "rises with effect size",
        "two-sided selection inflates",
        "sigma2_plate = max(var_lfc - mean_se2, 0)",
        "example scatter reproduces",
        "every figure design.md declares",
        "scored conditions are the splittable ones",
        "checksum recomputes from the file it names",
        "tranche content hash",
        "permutation",
        "promoted",
    ):
        assert fragment in names, f"no check covers {fragment!r}"


def test_a_perturbed_claim_fails_the_battery(artifacts: Path) -> None:
    # The `document` step's negative control (design.md): move one reported number away from the
    # per-condition values it is the mean of, and the battery must say so by name. Nothing else in
    # this file would notice a battery that had quietly stopped comparing anything.
    summary = artifacts / "rung0_reliability.csv"
    with summary.open(newline="") as fh:
        reader = csv.DictReader(fh)
        fields = reader.fieldnames
        assert fields is not None
        rows = list(reader)
    original = rows[0]["all_splithalf_mean_r"]
    rows[0]["all_splithalf_mean_r"] = "0.123"
    with summary.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    checks = vr.run_all_checks(artifacts)
    failed = [str(c.name) for c in checks if not c.ok]
    assert any("all: mean recomputes" in name for name in failed), (
        f"the battery passed on a summary whose all-gene mean was moved from {original} to 0.123; "
        f"failures reported: {failed}"
    )
    # And the recorded checksum notices the altered file, which is what ties what an auditor read
    # to what a reviewer later pulls.
    assert any("checksum recomputes" in name for name in failed), (
        f"the recorded checksums did not notice an altered artifact; failures: {failed}"
    )
