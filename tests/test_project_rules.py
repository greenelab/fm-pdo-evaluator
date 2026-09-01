"""One test per numbered project rule in ``docs/SPEC.md``.

The project rules are what every task must satisfy, whichever rung or dataset it is about. They
arrive with the work that first needs them, and so do their tests. Four rules live here so far:
rule 1 (every promoted result carries a provenance record), rule 2 (a task is named in the spec
tree, and a reversal of a task document's own lines carries a dated entry), rule 3 (the README
stays in step with the documents it summarises), and rule 4 (every measurement step declares a
positive and a negative control). Rules governing splits, statistics, capacity and embargo land
with the code and data they constrain.

Two kinds of test live here, and the difference decides how much a pass is worth:

* **Behavioural** — call the real function with a known answer and require the answer back.
* **Repository or artifact scan** — read the repository, its history, or its promoted outputs
  and look for a specific violating shape. A pass means that shape is absent, which is weaker
  than "the rule holds": a scan cannot see a violation written a way the pattern does not match.
  Each scan below says in its docstring what it does not prove.

Where a rule has nothing to check yet — no promoted results, no task folders — the test skips
rather than passing. A vacuous pass is worse than a skip, because a green run then reports
compliance that was never tested.

No rule here has per-instance exemptions yet, so there is nothing to exempt. The registry of
known violations, and the strict-xfail machinery that forces an entry out once its owning task
lands, arrives with the first rule that checks instances one at a time.
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

TASK_DOCUMENTS = ("design.md", "plan.md", "decisions.md")
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
    """A task document that rewrites its own history must say so, dated.

    The dated entry lives in the task's ``decisions.md`` (SPEC rule 2's amended location,
    2026-08-31) or in the rewritten document itself (the original form, still accepted).
    ``decisions.md`` is itself in the scanned set: it now carries the history this rule
    protects, so deleting or rewriting its entries demands a new dated entry the same way
    (found at the wave's scoped review -- the scan set had not followed the moved content).
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
            # SPEC rule 2's amended location: the dated entry may live in the task's
            # decisions.md rather than in the rewritten document itself
            decisions = str(Path(doc).parent / "decisions.md")
            decisions_before = _git("show", f"{base}:{decisions}") or ""
            decisions_after = _git("show", f"HEAD:{decisions}") or ""
            new_dates = set(DATED_ENTRY.findall(decisions_after)) - set(
                DATED_ENTRY.findall(decisions_before)
            )
        if not new_dates:
            offenders.append(doc)
    assert not offenders, (
        "these task documents rewrote existing lines without recording the change, dated, in"
        f" the task's decisions.md or the document itself: {offenders}"
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


# ======================================================================================
# Rule 4 — every measurement step carries a positive and a negative control
# ======================================================================================

MEASUREMENT_STEPS = ("load", "build", "restrict", "split", "fit", "score", "null")


@pytest.mark.step_score
@pytest.mark.step_null
def test_rule_04_every_task_declares_controls_for_its_measurement_steps() -> None:
    """Scan: a task touching measurement steps declares both control signs for each of them.

    Reads the ``**Steps**`` header line of each task's ``design.md``; for every measurement
    step named there, the document's Controls section must hold an entry for that step
    containing both a positive and a negative control. Does not prove the declared control is
    implemented anywhere, or that a plant is placed where it proves something — the known-answer
    tests are the implemented half, and the edge-case test below ties them to promotion.
    """
    folders = _task_folders()
    if not folders:
        pytest.skip("no task folders yet; nothing for this rule to check")

    problems: list[str] = []
    for folder in folders:
        design = folder / "design.md"
        if not design.exists():
            problems.append(f"{folder.name}: no design.md")
            continue
        text = design.read_text()
        header = re.search(r"\*\*Steps\*\*(.+)", text)
        if header is None:
            problems.append(f"{folder.name}: design.md has no **Steps** header line")
            continue
        steps = [s for s in MEASUREMENT_STEPS if re.search(rf"\b{s}\b", header.group(1))]
        if not steps:
            continue  # a documentation- or promotion-only task measures nothing
        controls = re.search(r"^## .*[Cc]ontrols.*$", text, re.M)
        if controls is None:
            problems.append(
                f"{folder.name}: touches measurement steps {steps} but design.md has no "
                "Controls section"
            )
            continue
        section = text[controls.end() :]
        nxt = re.search(r"^## ", section, re.M)
        if nxt:
            section = section[: nxt.start()]
        for step in steps:
            entry = re.search(rf"\*\*{step}\*\*(.*?)(?=\n- \*\*|\Z)", section, re.S)
            if entry is None:
                problems.append(f"{folder.name}: no control entry for step '{step}'")
                continue
            # Strip markdown emphasis before the colon-form search, so a bolded declaration
            # like "**positive**:" is recognized the same as a plain "positive:".
            body = entry.group(1).lower().replace("*", "")
            for sign in ("positive", "negative"):
                # Require the declaration FORM "positive:"/"negative:" (colon required), not a
                # bare occurrence of the word -- prose like "positive control omitted" contains
                # the word "positive" without declaring one, and used to pass this scan.
                if f"{sign}:" not in body:
                    problems.append(f"{folder.name}: step '{step}' lacks a '{sign}:' control")
    assert not problems, problems


@pytest.mark.step_score
@pytest.mark.step_null
def test_rule_04_edge_promoted_tasks_have_known_answer_tests() -> None:
    """A promoted result may not exist while no test carries the known-answer marker.

    Declared controls precede their implementation while a task is being built, so this binds
    only at promotion: once any result is promoted, the repository must hold at least one test
    marked ``known_answer``. A shape scan — it cannot tell whether the marked tests cover the
    controls the promoting task declared; the reviewer of that task's PR can.
    """
    if not _promoted_results():
        pytest.skip("no promoted results yet; declared controls precede their implementation")

    marked = [
        p for p in sorted((REPO / "tests").glob("test_*.py")) if "known_answer" in p.read_text()
    ]
    assert marked, (
        "results are promoted but no test file carries pytest.mark.known_answer; rule 4 "
        "requires the declared controls to be implemented before a number becomes evidence"
    )


# ======================================================================================
# Notebooks are committed without outputs (PROCESS section 3)
# ======================================================================================


@pytest.mark.step_document
def test_committed_notebooks_carry_no_outputs() -> None:
    """A committed notebook must hold no execution outputs.

    PROCESS section 3: the figures and numbers a reviewer sees are the ones their OWN
    execution produced. A notebook committed with outputs shows them someone else's run, which
    is exactly the transcribed-number problem the executable-verification rule exists to
    remove -- and stale outputs look identical to fresh ones in a diff.

    This is a mechanical check by script rather than by a reader, per docs/audit.md: it caught a
    verify notebook that reached a commit with fourteen outputs still in it, after being
    executed for testing and committed in the same window.
    """
    import json
    import subprocess

    repo = Path(__file__).resolve().parents[1]
    tracked = subprocess.run(
        ["git", "ls-files", "*.ipynb"],
        cwd=repo,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    if not tracked:
        pytest.skip("no notebooks tracked yet; nothing for this rule to check")
    offenders: list[str] = []
    for rel in tracked:
        nb = json.loads((repo / rel).read_text())
        n = sum(
            len(c.get("outputs", [])) for c in nb.get("cells", []) if c.get("cell_type") == "code"
        )
        if n:
            offenders.append(f"{rel} ({n} outputs)")
    assert not offenders, "notebooks committed with outputs: " + ", ".join(offenders)
