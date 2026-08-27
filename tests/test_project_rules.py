"""One test per numbered project rule in ``docs/SPEC.md``.

The project rules are what every task must satisfy, whichever rung or dataset it is about. They
arrive with the work that first needs them, and so do their tests: this repository holds
documents and the review apparatus, so it carries the two rules that bind those. Rules governing
splits, statistics, controls, capacity and embargo land with the code and data they constrain.

Two kinds of test live here, and the difference decides how much a pass is worth:

* **Behavioural** — call the real function with a known answer and require the answer back.
* **Repository or artifact scan** — read the repository, its history, or its promoted outputs
  and look for a specific violating shape. A pass means that shape is absent, which is weaker
  than "the rule holds": a scan cannot see a violation written a way the pattern does not match.
  Each scan below says in its docstring what it does not prove.

Where a rule has nothing to check yet — no promoted results, no task folders — the test skips
rather than passing. A vacuous pass is worse than a skip, because a green run then reports
compliance that was never tested.

Neither rule here has per-instance cases, so there is nothing to exempt. The registry of known
violations, and the strict-xfail machinery that forces an entry out once its owning task lands,
arrives with the first rule that checks instances one at a time.
"""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import pytest
from pydantic import ValidationError

from fmharness.schema import PromotedResult

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"
TASKS = REPO / "docs" / "tasks"


def _git(*args: str) -> str | None:
    """Run a git command, or None when git cannot answer (no repo, no such ref)."""
    try:
        out = subprocess.run(
            ["git", "-C", str(REPO), *args], capture_output=True, text=True, timeout=30
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout if out.returncode == 0 else None


# ======================================================================================
# Rule 1 — every promoted result carries its provenance record
# ======================================================================================


def _promoted_results() -> list[Path]:
    return sorted(RESULTS.glob("*/*.csv")) if RESULTS.is_dir() else []


@pytest.mark.step_promote
def test_rule_01_every_promoted_result_carries_a_complete_provenance_record() -> None:
    """Artifact scan: a promoted result has a record beside it saying how it was produced.

    Does not prove the record is true — only that it exists and is complete. Whether the commit
    it names is the one that ran is the edge case below, and is not checkable after the fact.
    """
    results = _promoted_results()
    if not results:
        pytest.skip("no promoted results yet; nothing for this rule to check")

    problems: list[str] = []
    for result in results:
        record = result.with_suffix(".provenance.json")
        if not record.exists():
            problems.append(f"{result.relative_to(REPO)}: no provenance record")
            continue
        try:
            PromotedResult.model_validate_json(record.read_text())
        except ValidationError as exc:
            problems.append(f"{record.relative_to(REPO)}: {exc.error_count()} invalid field(s)")
    assert not problems, problems


@pytest.mark.step_promote
def test_rule_01_edge_promoted_records_validate_against_the_schema() -> None:
    """One provenance format, defined in code, so a second cannot appear beside it.

    ``PromotedResult`` forbids extra fields and requires the three that cannot be recovered
    later — clean tree, producing commit, result checksum. A record written to some other shape
    fails here rather than sitting in the repository looking like evidence. What no test can
    check is whether a recorded value is *correct*; validation catches an absent field, not a
    wrong one.
    """
    records = sorted(RESULTS.glob("*/*.provenance.json")) if RESULTS.is_dir() else []
    if not records:
        pytest.skip("no provenance records yet; nothing for this rule to check")

    invalid: list[str] = []
    for record in records:
        try:
            PromotedResult.model_validate_json(record.read_text())
        except ValidationError as exc:
            fields = ", ".join(str(e["loc"][0]) for e in exc.errors() if e["loc"])
            invalid.append(f"{record.relative_to(REPO)}: {fields or exc.error_count()}")
    assert not invalid, invalid


# ======================================================================================
# Rule 2 — a reversal is written into the task's own documents
# ======================================================================================

TASK_DOCUMENTS = ("design.md", "plan.md")
DATED_ENTRY = re.compile(r"\b(20\d{2}-\d{2}-\d{2})\b")


def _task_folders() -> list[Path]:
    return sorted(p for p in TASKS.iterdir() if p.is_dir()) if TASKS.is_dir() else []


@pytest.mark.step_document
def test_rule_02_every_task_is_named_in_the_spec_tree() -> None:
    """Scan: nothing under docs/tasks/ is missing from the ladder in docs/SPEC.md.

    Cannot detect what rule 2 is really about — work that reverses a decision and writes nothing
    down. What it does catch is the next step of the same failure: a task document that exists
    but that nobody can find, because no branch of the ladder mentions it.
    """
    folders = _task_folders()
    if not folders:
        pytest.skip("no task folders yet; nothing for this rule to check")

    spec = (REPO / "docs" / "SPEC.md").read_text()
    unindexed = [
        f"docs/tasks/{d.name}/design.md"
        for d in folders
        if f"docs/tasks/{d.name}/design.md" not in spec
    ]
    assert not unindexed, f"present in the repository, absent from the spec tree: {unindexed}"


def _merge_base() -> str | None:
    """The commit this branch diverged from, trying a PR merge-parent then the upstream/fork default."""
    parents = (_git("rev-list", "--parents", "-n", "1", "HEAD") or "").split()
    if len(parents) >= 3:
        return parents[1]

    for ref in ("upstream/main", "origin/main", "main"):
        base = _git("merge-base", ref, "HEAD")
        if base and base.strip():
            return base.strip()
    return None


@pytest.mark.step_document
def test_rule_02_edge_non_additive_task_edits_carry_a_dated_entry() -> None:
    """A task document that rewrites its own history must say so, dated, at its foot.

    Appending to a task document is free. Deleting or rewriting lines already there is a
    reversal of something the document previously asserted, and the rule requires the old choice
    and the reason for changing it to be recorded rather than overwritten. That distinction is a
    property of the diff, so it can be checked.

    Does not prove the dated entry describes the change, and cannot see a reversal carried out
    in code while the task document is left untouched — that residual is stated with the rule.
    """
    base = _merge_base()
    if base is None:
        pytest.skip("no merge base available; nothing to compare this branch against")

    changed = _git("diff", "--name-only", base, "--", "docs/tasks")
    if changed is None:
        pytest.skip("git could not report the changed files")
    documents = [f for f in changed.split() if Path(f).name in TASK_DOCUMENTS]
    if not documents:
        pytest.skip("no task documents changed against the merge base")

    offenders: list[str] = []
    for doc in documents:
        diff = _git("diff", "--unified=0", base, "--", doc) or ""
        removed = [
            line
            for line in diff.splitlines()
            if line.startswith("-") and not line.startswith("---") and line[1:].strip()
        ]
        if not removed:
            continue  # purely additive: nothing was reversed
        before = _git("show", f"{base}:{doc}") or ""
        # the diff above compares commits, so the "after" side must come from HEAD too — reading
        # the working tree would let an uncommitted entry satisfy a check on committed work,
        # and would report a false offender for a file deleted only in the working tree
        after = _git("show", f"HEAD:{doc}")
        if after is None:
            continue  # removed at HEAD; there is no document left to carry an entry
        new_dates = set(DATED_ENTRY.findall(after)) - set(DATED_ENTRY.findall(before))
        if not new_dates:
            offenders.append(doc)
    assert not offenders, (
        "these task documents rewrote existing lines without recording the change, dated, at"
        f" the foot of the document: {offenders}"
    )


# ======================================================================================
# Rule 3 — the README stays in step with the documents it summarises
# ======================================================================================

README = REPO / "README.md"
PROJECT_DOCUMENTS = ("docs/SPEC.md", "docs/PROCESS.md", "docs/STATE.md")


@pytest.mark.step_document
def test_rule_03_readme_links_to_the_project_documents() -> None:
    """The README points a reader at each document it summarises, and the links resolve.

    Does not read the prose for accuracy — a README can link correctly and still describe a
    ladder that changed last month. That is the edge case below.
    """
    readme = README.read_text()
    link_targets = [
        target.split("#")[0].split("?")[0]
        for target in re.findall(r"\]\(([^)]+)\)", readme)
        if not target.startswith(("http", "#", "mailto"))
    ]
    missing = [doc for doc in PROJECT_DOCUMENTS if doc not in link_targets]
    assert not missing, f"README does not link to {missing}"

    broken = [target for target in link_targets if not (REPO / target).exists()]
    assert not broken, f"README links that do not resolve: {broken}"


@pytest.mark.step_document
def test_rule_03_edge_readme_is_revisited_when_the_ladder_changes() -> None:
    """A change to the ladder requires the summary of it to have been revisited.

    Proves the README was opened in the same change, not that it was revised well. Reading a
    summary for faithfulness is a reviewer's job; noticing that nobody looked is not.
    """
    base = _merge_base()
    if base is None:
        pytest.skip("no merge base available; nothing to compare this branch against")

    changed = _git("diff", "--name-only", base, "--", "docs/SPEC.md", "README.md")
    if changed is None:
        pytest.skip("git could not report the changed files")
    touched = set(changed.split())
    if "docs/SPEC.md" not in touched:
        pytest.skip("the spec did not change against the merge base")

    spec_diff = _git("diff", base, "--", "docs/SPEC.md") or ""
    ladder_changed = any(
        line.startswith(("+", "-")) and "Rung" in line for line in spec_diff.splitlines()
    )
    if not ladder_changed:
        pytest.skip("the spec changed but its ladder did not")

    assert "README.md" in touched, (
        "the ladder in docs/SPEC.md changed but README.md did not; the summary a reader opens "
        "first has to be revisited in the same change"
    )
