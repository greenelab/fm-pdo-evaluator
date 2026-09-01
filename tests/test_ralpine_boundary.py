"""Security-boundary regression tests for scripts/alpine/ralpine.

These are pure TEXT-SCAN tests over the script's source -- no ssh, no execution, nothing
remote. That makes them the weaker-but-cheap half of enforcement: they can only catch a
regression that is visible in the script's literal text (an escape-capable command readmitted
to READ_ONLY, a validation line moved after the remote call it was supposed to guard), not a
logic bug that keeps the text looking right while behaving wrong. They exist because the READ_
ONLY allowlist and the log/jobinfo fixed-command paths are exactly the kind of thing that drifts
silently in a script nobody runs through a test suite -- `find` was readmitted once already
(closed 2026-08-28) and this file exists to make the next regression fail loudly instead.
"""

from __future__ import annotations

import re
from pathlib import Path

RALPINE = (Path(__file__).resolve().parents[1] / "scripts" / "alpine" / "ralpine").read_text()

# Commands that would let `ralpine run <cmd>` escape the read-only boundary: `find`/`scontrol`
# write or destroy through predicates/subcommands (closed 2026-08-28); the rest are shells or
# interpreters that would admit arbitrary code by construction if ever added to READ_ONLY.
ESCAPE_CAPABLE = (
    "find",
    "scontrol",
    "file",
    "nvidia-smi",
    "git",
    "bash",
    "sh",
    "xargs",
    "python",
    "perl",
)


def _read_only_tokens() -> list[str]:
    match = re.search(r"READ_ONLY=\(\n(.*?)\n\)", RALPINE, re.S)
    assert match is not None, "could not find the READ_ONLY array literal in ralpine"
    return match.group(1).split()


def _case_block(verb: str) -> str:
    """The body of one `case` arm, from `  <verb>)` to its closing `    ;;`.

    Only matches ralpine's multi-line case-arm style (`log`, `jobinfo`, `switch`, ...), not the
    single-line arms (`ls)      remote ls ... ;;`) -- the verbs this file inspects are all
    multi-line.
    """
    match = re.search(rf"\n  {re.escape(verb)}\)\n(.*?)\n    ;;\n", RALPINE, re.S)
    assert match is not None, f"could not find a multi-line '{verb})' case block in ralpine"
    return match.group(1)


def test_read_only_allowlist_excludes_every_escape_capable_command() -> None:
    tokens = _read_only_tokens()
    present = [cmd for cmd in ESCAPE_CAPABLE if cmd in tokens]
    assert not present, (
        f"escape-capable command(s) {present} are in READ_ONLY; remote()'s pass-through would "
        "let `ralpine run <cmd> ...` execute arbitrary remote commands through them"
    )


# `remote_fixed "` (space, then an opening quote) targets the actual CALL SITE -- every
# invocation in this script passes a double-quoted string argument -- not the bare word
# "remote_fixed" as it also appears in this file's own prose comments, which would give a
# false-early index and defeat the ordering check below.
_REMOTE_FIXED_CALL = 'remote_fixed "'


def test_log_verb_validates_the_pattern_before_the_remote_fixed_call() -> None:
    block = _case_block("log")
    validation_idx = block.find(r"^[A-Za-z0-9._-]*$")
    remote_fixed_idx = block.find(_REMOTE_FIXED_CALL)
    assert validation_idx != -1, "log's pattern-validation regex is missing from the case block"
    assert remote_fixed_idx != -1, "log no longer calls remote_fixed"
    assert validation_idx < remote_fixed_idx, (
        "log's pattern validation must run BEFORE the remote_fixed call it guards -- an "
        "unvalidated pattern would be spliced into a remote find invocation"
    )


def test_jobinfo_verb_validates_the_job_id_before_the_remote_fixed_call() -> None:
    block = _case_block("jobinfo")
    validation_idx = block.find(r"^[0-9]+$")
    remote_fixed_idx = block.find(_REMOTE_FIXED_CALL)
    assert validation_idx != -1, "jobinfo's numeric-job-id validation regex is missing"
    assert remote_fixed_idx != -1, "jobinfo no longer calls remote_fixed"
    assert validation_idx < remote_fixed_idx, (
        "jobinfo's job-id validation must run BEFORE the remote_fixed call it guards -- an "
        "unvalidated job id would be spliced into a remote scontrol invocation"
    )
    assert "scontrol show job" in block, "jobinfo must invoke the read-only 'scontrol show job'"


def test_remote_enforces_the_allowlist_and_metacharacter_rejection() -> None:
    match = re.search(r"\nremote\(\) \{(.*?)\n\}\n", RALPINE, re.S)
    assert match is not None, "could not find the remote() function body in ralpine"
    body = match.group(1)
    assert "is_read_only" in body, "remote() must check its command against READ_ONLY"
    assert "reject_metacharacters" in body, "remote() must reject shell metacharacters"
