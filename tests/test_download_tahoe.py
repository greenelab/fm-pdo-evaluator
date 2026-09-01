"""Behavior tests for scripts/download_tahoe_pseudobulk_de.py: the download is pinned to a
fixed dataset revision and always calls through to `snapshot_download` (its own resume/verify
logic is what makes a re-run cheap and correct), never short-circuited by a presence check that
an interrupted prior download could satisfy with only a subset of shards on disk."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path
from typing import Any

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "download_tahoe_pseudobulk_de",
    Path(__file__).resolve().parents[1] / "scripts" / "download_tahoe_pseudobulk_de.py",
)
assert _SPEC is not None and _SPEC.loader is not None
dt = importlib.util.module_from_spec(_SPEC)
sys.modules["download_tahoe_pseudobulk_de"] = dt
_SPEC.loader.exec_module(dt)


def _install_fake_huggingface_hub(
    monkeypatch: pytest.MonkeyPatch, calls: list[dict[str, Any]]
) -> None:
    """`huggingface_hub` is imported inside `main()`, not at module scope, precisely so tests
    can inject a fake before calling it -- there is no real HF Hub network access here."""

    def fake_snapshot_download(repo_id: str, **kwargs: Any) -> str:
        calls.append({"repo_id": repo_id, **kwargs})
        return str(kwargs["local_dir"])

    fake_module = types.ModuleType("huggingface_hub")
    fake_module.snapshot_download = fake_snapshot_download  # type: ignore[attr-defined]
    monkeypatch.setitem(sys.modules, "huggingface_hub", fake_module)


def test_download_pins_the_default_revision_and_passes_dataset_repo_type_and_local_dir(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    calls: list[dict[str, Any]] = []
    _install_fake_huggingface_hub(monkeypatch, calls)
    local_dir = tmp_path / "pool"
    monkeypatch.setattr(
        sys, "argv", ["download_tahoe_pseudobulk_de.py", "--local-dir", str(local_dir)]
    )

    dt.main()

    assert len(calls) == 1
    call = calls[0]
    assert call["repo_id"] == dt.TAHOE
    assert call["repo_type"] == "dataset"
    assert call["revision"] == dt.DEFAULT_REVISION
    assert call["local_dir"] == str(local_dir)


def test_download_still_runs_when_a_matching_parquet_already_exists(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An interrupted prior download that left a subset of shards behind must not skip the
    pull: only `snapshot_download`'s own resume/verify against the pinned revision can tell
    a complete download from a partial one."""
    calls: list[dict[str, Any]] = []
    _install_fake_huggingface_hub(monkeypatch, calls)
    local_dir = tmp_path / "pool"
    de_dir = local_dir / dt.DE
    de_dir.mkdir(parents=True)
    (de_dir / "shard_0.parquet").write_bytes(b"stand-in for a partially downloaded shard")
    monkeypatch.setattr(
        sys, "argv", ["download_tahoe_pseudobulk_de.py", "--local-dir", str(local_dir)]
    )

    dt.main()

    assert len(calls) == 1, "snapshot_download must run even when a parquet file already exists"
