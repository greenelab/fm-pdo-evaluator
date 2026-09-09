"""The promotion gate: a result becomes evidence only with a valid record beside it,
and an artifact cannot change under its claim (SPEC rule 1; design 'Run and promotion')."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

from fmharness.schema import PromotedResult

_SPEC = importlib.util.spec_from_file_location(
    "promote_result", Path(__file__).resolve().parents[1] / "scripts" / "promote_result.py"
)
assert _SPEC is not None and _SPEC.loader is not None
pr = importlib.util.module_from_spec(_SPEC)
sys.modules["promote_result"] = pr
_SPEC.loader.exec_module(pr)


@pytest.fixture()
def repo(tmp_path: Path) -> Path:
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "scripts").mkdir()
    (tmp_path / "scripts" / "made_it.py").write_text("print('x')\n")
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs" / "input.txt").write_text("input-bytes\n")
    (tmp_path / "result.csv").write_text("a,b\n1,2\n")
    subprocess.run(
        ["git", "-C", str(tmp_path), "-c", "user.email=t@t", "-c", "user.name=t", "add", "-A"],
        check=True,
    )
    subprocess.run(
        [
            "git",
            "-C",
            str(tmp_path),
            "-c",
            "user.email=t@t",
            "-c",
            "user.name=t",
            "commit",
            "-qm",
            "init",
        ],
        check=True,
    )
    return tmp_path


def _promote(repo: Path, **kw: object) -> Path:
    defaults: dict[str, object] = dict(
        task="rung0-replicate-ceiling",
        result=repo / "result.csv",
        script="scripts/made_it.py",
        inputs=[repo / "docs" / "input.txt"],
        seed=0,
        data_commit="c" * 64,
        args={"tranche_id": "tahoe100m-pseudobulk-de.v1"},
        job_id="123",
        log=None,
        repo=repo,
        code_commit="a" * 40,
    )
    defaults.update(kw)
    return pr.promote(**defaults)


def _head(repo: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
    ).stdout.strip()


def test_promotion_records_the_run_commit_from_the_sidecar_not_head(repo: Path) -> None:
    """The producing commit is the one the RUN was made at. The pooled-dose promotion of
    2026-09-02 recorded HEAD at promotion time -- a commit whose code held dose fixed -- so a
    reader checking out the recorded commit could not reproduce the promoted number."""
    import json

    run_sha = "b" * 40
    (repo / "result.params.json").write_text(
        json.dumps({"result": "result.csv", "git_sha": run_sha, "slurm_job_id": "999"}) + "\n"
    )
    record_path = _promote(repo, code_commit=None, job_id=None)
    record = PromotedResult.model_validate_json(record_path.read_text())
    assert record.environment.code_commit == run_sha, "the sidecar's commit, not HEAD"
    assert record.promotion_commit == _head(repo), "and the promotion commit, separately"
    assert record.job_id == "999", "the job id defaults from the sidecar too"


def test_promotion_refuses_when_no_producing_commit_is_known(repo: Path) -> None:
    with pytest.raises(SystemExit, match="not known"):
        _promote(repo, code_commit=None)


def test_promotion_writes_a_schema_valid_record_beside_the_result(repo: Path) -> None:
    record_path = _promote(repo)
    promoted = repo / "results" / "rung0-replicate-ceiling" / "result.csv"
    assert promoted.exists()
    record = PromotedResult.model_validate_json(record_path.read_text())
    assert record.result_sha256 == pr.sha256_of(promoted)
    assert record.clean_tree is True
    assert record.environment.cuda_deterministic is False
    assert record.environment.data_commit == "c" * 64
    assert record.args["tranche_id"] == "tahoe100m-pseudobulk-de.v1"


def test_promotion_refuses_when_the_promoted_copy_differs(repo: Path) -> None:
    _promote(repo)
    (repo / "result.csv").write_text("a,b\n9,9\n")  # task-side copy changed
    with pytest.raises(SystemExit, match="differ"):
        _promote(repo)


def test_promotion_refuses_repromotion_when_unchanged(repo: Path) -> None:
    """A provenance record is immutable once written, even when the result did not change --
    re-promoting must be a deliberate act (remove the old record first), not a silent
    overwrite."""
    _promote(repo)
    with pytest.raises(SystemExit, match=r"record|immutable"):
        _promote(repo)


def test_promotion_refuses_a_result_with_no_inputs(repo: Path) -> None:
    with pytest.raises(SystemExit, match="input"):
        _promote(repo, inputs=[])


def test_promotion_records_a_dirty_tree_honestly(repo: Path) -> None:
    (repo / "scripts" / "made_it.py").write_text("print('changed')\n")
    record_path = _promote(repo)
    assert PromotedResult.model_validate_json(record_path.read_text()).clean_tree is False


def test_promotion_refuses_a_script_not_in_the_repo(repo: Path) -> None:
    with pytest.raises(SystemExit, match="script"):
        _promote(repo, script="scripts/never_existed.py")


def test_promotion_ignores_an_untracked_stray_file(repo: Path) -> None:
    """An untracked file does not make clean_tree False (PROCESS §2: working trees carry
    untracked data by design); only tracked-file modifications should."""
    (repo / "stray_untracked.txt").write_text("not part of the repo\n")
    record_path = _promote(repo)
    assert PromotedResult.model_validate_json(record_path.read_text()).clean_tree is True


def test_promotion_keys_a_labeled_input_by_its_label(repo: Path) -> None:
    """A labeled --input (LABEL=PATH, parsed in main()) lands under its label, not its path,
    with the hash still computed from the real file."""
    input_path = repo / "docs" / "input.txt"
    record_path = _promote(repo, inputs=[input_path], input_labels={input_path: "gene_panel"})
    record = PromotedResult.model_validate_json(record_path.read_text())
    assert record.inputs == {"gene_panel": pr.sha256_of(input_path)}


@pytest.mark.step_promote
def test_promotion_refuses_an_artifact_whose_checksum_moved_since_the_audit(
    tmp_path: Path,
) -> None:
    """The promote step's negative control, and the refusal that closes the audit window.

    The audit reads the run's artifacts in the working tree before they are committed, and
    records each one's sha256. Between that read and the promotion commit, nothing structural
    stops an artifact changing -- so promotion recomputes the checksum and refuses when it has
    moved. Without this the window is unchecked, and "the audit passed" would say nothing about
    the bytes that got promoted.
    """
    import json

    result = tmp_path / "rung0_reliability.csv"
    result.write_text("all_splithalf_mean_r\n0.135\n")
    sums = tmp_path / "audit_checksums.json"
    sums.write_text(json.dumps({result.name: pr.sha256_of(result)}))

    # Matching checksum: the guard lets it through.
    pr._refuse_if_checksums_moved(result, sums)

    # One byte different, and it must refuse by name.
    result.write_text("all_splithalf_mean_r\n0.136\n")
    with pytest.raises(SystemExit) as excinfo:
        pr._refuse_if_checksums_moved(result, sums)
    assert "changed after the audit read it" in str(excinfo.value)
    assert result.name in str(excinfo.value)


@pytest.mark.step_promote
def test_promotion_refuses_when_the_audit_never_read_the_artifact(tmp_path: Path) -> None:
    """A missing record is the same evidential state as a mismatched one: nothing establishes
    that what is being promoted is what was reviewed. Both are refused, and a missing checksums
    file is refused rather than treated as "no constraint"."""
    import json

    result = tmp_path / "rung0_reliability.csv"
    result.write_text("x\n1\n")
    sums = tmp_path / "audit_checksums.json"
    sums.write_text(json.dumps({"some_other_table.csv": "0" * 64}))
    with pytest.raises(SystemExit, match="did not read this artifact"):
        pr._refuse_if_checksums_moved(result, sums)

    with pytest.raises(SystemExit, match="does not exist"):
        pr._refuse_if_checksums_moved(result, tmp_path / "absent.json")

    # And with no record requested at all the guard is inert, so existing callers are unaffected.
    pr._refuse_if_checksums_moved(result, None)
