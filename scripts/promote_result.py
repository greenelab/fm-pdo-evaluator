"""Promote a run's output into committed evidence (SPEC rule 1; PROCESS §1 'Promote').

A log is not evidence; a result becomes evidence when a claim cites it, and citing means
a ``results/<task-slug>/`` copy with a schema-validated ``PromotedResult`` beside it.
Three fields cannot be reconstructed later and are written here: whether the tree was
clean, the producing commit, and the artifact's checksum. Promotion REFUSES when the
task-side copy and an existing promoted copy differ -- an artifact must not change under
its claim.

Usage (rung 0):
    uv run python scripts/promote_result.py \
        --task rung0-assay-reliability \
        --result docs/tasks/rung0-assay-reliability/rung0_reliability.csv \
        --script scripts/delta_reproducibility.py \
        --input gene_panel=/path/to/common_panel.txt \
        --input drug_cids=/path/to/tahoe_target_cids.txt \
        --seed 0 --data-commit <tranche content_hash> \
        --arg tranche_id=tahoe100m-pseudobulk-de.v1 --job-id <slurm id> \
        --log results/rung0-assay-reliability/<job log>

``--input`` takes ``LABEL=PATH``, so the promoted record's ``inputs`` dict is keyed by a durable
label rather than a scratch path; a bare ``PATH`` (no ``=``) is still accepted for back-compat
and keyed by the path string, as before.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fmharness.schema import EnvironmentSnapshot, PromotedResult


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", "-C", str(repo), *args], capture_output=True, text=True, check=True
    ).stdout.strip()


def _refuse_if_checksums_moved(result: Path, audit_checksums: Path | None) -> None:
    """Refuse to promote an artifact that changed after the audit read it.

    The audit reads the run's artifacts in the working tree, before they are committed
    (PROCESS section 1, "What reaches GitHub, and when"), and records the sha256 of each one it
    read. That leaves a window between what was audited and what gets committed, and this is
    what closes it: promotion recomputes the checksum of the artifact it is about to promote and
    refuses when it no longer matches the recorded one.

    A missing record is refused too when a path was given. "The checksums file is not there" is
    the same evidential state as "the checksums do not match" -- in both cases nothing has
    established that what is being promoted is what was reviewed.
    """
    if audit_checksums is None:
        return
    if not audit_checksums.exists():
        raise SystemExit(
            f"--audit-checksums {audit_checksums} does not exist. Promotion checks the promoted "
            "artifact against the checksum the audit recorded; with no record there is nothing "
            "establishing that this file is the one that was reviewed."
        )
    recorded = json.loads(audit_checksums.read_text())
    name = result.name
    if name not in recorded:
        raise SystemExit(
            f"{name} is not in {audit_checksums}. The audit did not read this artifact, so "
            "promoting it would put a number into the record that no audit covered."
        )
    now = sha256_of(result)
    if now != recorded[name]:
        raise SystemExit(
            f"{name} changed after the audit read it.\n"
            f"  audit recorded: {recorded[name]}\n"
            f"  now:            {now}\n"
            "Re-run the audit against the current artifacts, or promote the ones it read."
        )


def promote(
    *,
    task: str,
    result: Path,
    script: str,
    inputs: list[Path],
    seed: int,
    data_commit: str,
    args: dict[str, str],
    job_id: str | None,
    log: Path | None,
    repo: Path,
    input_labels: dict[Path, str] | None = None,
    audit_checksums: Path | None = None,
) -> Path:
    repo = repo.resolve()
    _refuse_if_checksums_moved(result, audit_checksums)
    if not (repo / script).exists():
        raise SystemExit(
            f"--script {script} is not in the repo; a result whose producing script "
            "cannot be found cannot be regenerated"
        )
    if not inputs:
        raise SystemExit(
            "at least one --input is required; a result with no recorded inputs cannot "
            "be checked against a rerun"
        )
    missing = [p for p in inputs if not p.exists()]
    if missing:
        raise SystemExit(f"declared inputs not found: {[str(p) for p in missing]}")

    # Check clean_tree status before writing any files. ``-uno`` restricts this to tracked
    # files: this project's working trees deliberately carry untracked data (rung 0's plan,
    # docs/tasks/rung0-assay-reliability/plan.md, Global constraints), locally and on the
    # cluster, so counting untracked files would make clean_tree permanently False and
    # meaningless. What this flag records is whether the tracked code and docs carried
    # uncommitted modifications at promotion time.
    clean_tree_status = _git(repo, "status", "--porcelain", "-uno") == ""
    code_commit = _git(repo, "rev-parse", "HEAD")

    out_dir = repo / "results" / task
    out_dir.mkdir(parents=True, exist_ok=True)
    dest = out_dir / result.name
    src_hash = sha256_of(result)
    if dest.exists() and sha256_of(dest) != src_hash:
        raise SystemExit(
            f"refusing: {dest.relative_to(repo)} exists and its checksum differs from "
            f"{result} -- the promoted copy and the task-side copy differ"
        )
    record_path = dest.with_suffix(".provenance.json")
    if record_path.exists():
        raise SystemExit(
            f"refusing: {record_path.relative_to(repo)} already exists -- a provenance "
            "record is immutable once written, even when the result is unchanged. "
            "A deliberate re-promotion (e.g. a correction) requires removing the old "
            "record first"
        )
    dest.write_bytes(result.read_bytes())

    record = PromotedResult(
        result=str(dest.relative_to(repo)),
        result_sha256=sha256_of(dest),
        task=task,
        script=script,
        args={k: str(v) for k, v in args.items()},
        inputs={(input_labels or {}).get(p, str(p)): sha256_of(p) for p in inputs},
        log=str(log) if log else None,
        log_sha256=sha256_of(log) if log else None,
        job_id=job_id,
        clean_tree=clean_tree_status,
        environment=EnvironmentSnapshot(
            code_commit=code_commit,
            python_version=platform.python_version(),
            seed=seed,
            cuda_deterministic=False,
            data_commit=data_commit,
        ),
        promoted_at=datetime.now(UTC),
    )
    record_path.write_text(record.model_dump_json(indent=2) + "\n")
    print(f"promoted -> {dest.relative_to(repo)}")
    print(f"           {record_path.relative_to(repo)}")
    return record_path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", required=True)
    ap.add_argument("--result", required=True, type=Path)
    ap.add_argument("--script", required=True)
    ap.add_argument(
        "--input",
        action="append",
        default=[],
        dest="raw_inputs",
        help="an input this result depends on, as LABEL=PATH (e.g. gene_panel=/path/to/"
        "common_panel.txt) so the record keys it by a durable label; a bare PATH (no '=') "
        "is accepted for back-compat and keyed by the path string. Repeatable.",
    )
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--data-commit", required=True)
    ap.add_argument("--arg", action="append", default=[], help="key=value, repeatable")
    ap.add_argument(
        "--audit-checksums",
        type=Path,
        default=None,
        help="the audit's checksum record (audit_checksums.json beside the artifacts). When "
        "given, promotion refuses if the artifact's checksum has moved since the audit read "
        "it, which is what closes the window opened by auditing uncommitted artifacts.",
    )
    ap.add_argument("--job-id", default=None)
    ap.add_argument("--log", type=Path, default=None)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ns = ap.parse_args()

    inputs: list[Path] = []
    input_labels: dict[Path, str] = {}
    for raw in ns.raw_inputs:
        label, sep, rest = raw.partition("=")
        p = Path(rest) if sep else Path(raw)
        if sep:
            input_labels[p] = label
        inputs.append(p)

    promote(
        task=ns.task,
        result=ns.result,
        script=ns.script,
        inputs=inputs,
        seed=ns.seed,
        data_commit=ns.data_commit,
        args=dict(kv.split("=", 1) for kv in ns.arg),
        job_id=ns.job_id,
        log=ns.log,
        repo=ns.repo,
        input_labels=input_labels,
        audit_checksums=ns.audit_checksums,
    )


if __name__ == "__main__":
    main()
