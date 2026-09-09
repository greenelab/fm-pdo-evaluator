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
import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any
from unittest import mock

import pytest

from tests.test_rung0_controls import dr, write_fixture_pool

pytestmark = pytest.mark.step_document

REPO = Path(__file__).resolve().parents[1]

_SPEC = importlib.util.spec_from_file_location("verify_rung0", REPO / "scripts" / "verify_rung0.py")
assert _SPEC is not None and _SPEC.loader is not None
vr: Any = importlib.util.module_from_spec(_SPEC)
sys.modules["verify_rung0"] = vr
_SPEC.loader.exec_module(vr)

#: The battery is exercised against a run the test session makes itself: the real ``main`` over
#: a synthetic pool with planted plate effects, two dose levels and responders, in three gene
#: slices. That run always exists, so the battery's arithmetic is tested on every checkout;
#: the committed cluster run is checked separately below, and only when its parameter sidecar
#: says it was produced under the schema this battery reads. RUNG0_EXAMPLE_RUN overrides the
#: synthetic run with a directory of a real run's artifacts.
TASK_DIR = REPO / "docs" / "tasks" / "rung0-assay-reliability"

#: Both the summary table's trust map and this task's pull-request description quote the size of
#: the battery. Pinning it here means adding or removing a check forces those transcriptions to be
#: revisited rather than quietly going stale -- the one number in the wave that nothing else
#: guards. Skipped checks are counted: whether the permutation job has run changes what can be
#: checked, never how many claims the battery is answerable for.
EXPECTED_CHECKS = 77


@pytest.fixture(scope="session")
def synthetic_run(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One run of the real ``main`` on a synthetic pool, shared by every test in this file."""
    override = os.environ.get("RUNG0_EXAMPLE_RUN")
    if override and (Path(override) / "rung0_reliability.csv").exists():
        return Path(override)
    root = tmp_path_factory.mktemp("synthetic_run")
    pool = write_fixture_pool(
        root,
        n_lines=6,
        n_drugs=4,
        n_genes=400,
        n_responders=120,
        doses=(0.05, 5.0),
        plates=("P1", "P2", "P3", "P4", "P5", "P6"),
        plate_offset_sd=0.3,
        signal_sd_by_drug=(0.3, 0.6, 1.0, 2.0),
        seed=41,
    )
    out = root / "out"
    argv = [
        "delta_reproducibility.py",
        "--local-dir",
        str(pool.parent.parent),
        "--out-dir",
        str(out),
        "--min-genes",
        "20",
        "--n-perm",
        "40",
        "--slices",
        "3",
    ]
    with mock.patch.object(sys, "argv", argv):
        dr.main()
    return out


@pytest.fixture
def artifacts(tmp_path: Path, synthetic_run: Path) -> Path:
    """A private copy of the run's artifacts, so a test may perturb them without harm."""
    destination = tmp_path / "run"
    shutil.copytree(synthetic_run, destination, ignore=shutil.ignore_patterns("cache"))
    return destination


def test_the_committed_cluster_run_verifies() -> None:
    """The run pulled into the task folder passes the battery -- once it is a run this battery
    reads. The parameter sidecar records the split rule from the 2026-09-09 code on; artifacts
    from before it are a different schema and are skipped by name rather than failed."""
    sidecar = TASK_DIR / "rung0_reliability.params.json"
    if not sidecar.exists():
        pytest.skip("no cluster run in the task folder")
    if "split_rule" not in json.loads(sidecar.read_text()):
        pytest.skip("the task folder holds a run from before the alternating split (2026-09-09)")
    checks = vr.run_all_checks(TASK_DIR)
    failures = _failures(checks)
    assert not failures, "the committed run and its documents disagree:\n" + "\n".join(failures)


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
        "rises with effect size, ranked by half0",
        "rises with effect size, ranked by half1",
        "two-sided selection inflates",
        "pooled between-plate share recomputes",
        "responders' pooled share",
        "sigma2_plate_signed = var_lfc - mean_se2",
        "control pool's pooled share",
        "example scatter reproduces",
        "every figure design.md declares",
        "split alternates over sorted plate ids",
        "scored conditions are the replicated triples",
        "dose-strata row recomputes",
        "commit the run was made at",
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
    # And where the run recorded checksums, they notice the altered file too -- that is what ties
    # what an auditor read to what a reviewer later pulls. A run that was killed before writing
    # them has nothing to check against, and this asserts the difference rather than assuming a
    # complete run: the arithmetic check above catches the edit either way.
    if (artifacts / vr.CHECKSUMS).exists():
        assert any("checksum recomputes" in name for name in failed), (
            f"the recorded checksums did not notice an altered artifact; failures: {failed}"
        )
    else:
        assert any("checksum recomputes" in str(c.name) and c.skipped for c in checks), (
            "with no checksum record the battery must SKIP that check, not silently pass it"
        )
