"""Shared helpers for download scripts under scripts/download/.

Common manifest shape both downloaders produce::

    {
      "dataset": "<short id>",
      "release": {<dataset-specific release/version keys>},
      "updated_utc": "<ISO8601>",
      "files": {
        "<relative_path>": {
          "sha256": "...",
          "bytes": N,
          "source_uri": "https://... | synapse://syn_id"
        }
      }
    }

Per-file extras (e.g. ``rows`` / ``cols`` / ``columns`` for Synapse Tables) are
allowed on top of this; downstream consumers ignore unknown keys.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import UTC, datetime
from pathlib import Path

MANIFEST_NAME = "manifest.json"


def sha256_file(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def load_manifest(manifest_path: Path, default: dict) -> dict:
    if manifest_path.exists():
        return json.loads(manifest_path.read_text())
    return default


def write_manifest(manifest_path: Path, manifest: dict) -> None:
    manifest["updated_utc"] = datetime.now(UTC).isoformat(timespec="seconds")
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def skip_or_fail_on_hash(name: str, dest: Path, manifest_files: dict) -> bool:
    """Return True if ``dest`` exists with sha matching the manifest record.

    Returns False when ``dest`` is missing or has no recorded sha — caller
    should (re-)fetch. Exits non-zero on mismatch so a silently changed
    upstream file does not slip past unnoticed.
    """
    if not dest.exists():
        return False
    recorded = manifest_files.get(name, {}).get("sha256")
    if not recorded:
        return False
    actual = sha256_file(dest)
    if recorded != actual:
        sys.exit(
            f"[fail] {name} sha256 mismatch -- expected {recorded}, got {actual}. "
            "Delete the file to re-download, or investigate upstream change."
        )
    return True


def verify_manifest(output_dir: Path, manifest_path: Path) -> int:
    """Re-check every recorded file under ``output_dir`` against its manifest sha.

    Returns the number of errors (0 = clean). Prints one ``[fail]`` line per
    bad file. Returns 0 silently when the manifest is missing (no records to
    verify), so callers can iterate over multiple optional manifests.
    """
    if not manifest_path.exists():
        return 0
    manifest = json.loads(manifest_path.read_text())
    bad: list[str] = []
    for name, rec in manifest.get("files", {}).items():
        path = output_dir / name
        if not path.exists():
            bad.append(f"{name}: missing")
            continue
        actual = sha256_file(path)
        if actual != rec["sha256"]:
            bad.append(f"{name}: sha256 mismatch (expected {rec['sha256']}, got {actual})")
    if bad:
        for b in bad:
            print(f"[fail] {b}")
        return len(bad)
    n = len(manifest.get("files", {}))
    print(f"[ok  ] {manifest_path.parent.name}/{manifest_path.name}: {n} files verified")
    return 0
