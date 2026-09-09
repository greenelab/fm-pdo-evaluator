"""One authenticated bulk download of the Tahoe-100M pseudobulk DE config to scratch.

Extracted from the archived lineage's ``build_tahoe_pseudobulk_deltas.py`` (branch
``rung0-replicate-ceiling-old-lineage`` on origin), whose delta-bundle aggregation
imports rung-1 machinery and arrives with rung 1. This half is rung 0's provenance
chain: it reproduces the 2026-07-24 pull recorded in docs/DATA.md.

The table is a flat 1,026-file shard set with no drug partition; run as ONE process and
authenticate first (``hf auth login``) so the pull is not rate-limited.

    python scripts/download_tahoe_pseudobulk_de.py \\
        --local-dir /scratch/alpine/$USER/tahoe_pseudobulk_de
"""

from __future__ import annotations

import argparse
from pathlib import Path

TAHOE = "tahoebio/Tahoe-100M"
DE = "pseudobulk_differential_expression"
# The revision registered as tahoe100m-pseudobulk-de.v1 (docs/tasks/rung0-assay-reliability/
# verification.md): "registered tahoe100m-pseudobulk-de.v1: 1026 shards, version
# 2dc57900b7981cfcf5e211527169a0b006546a95". Pinning here means a re-download always
# reproduces the same tranche instead of silently picking up whatever HF Hub serves as
# the dataset's current revision.
DEFAULT_REVISION = "2dc57900b7981cfcf5e211527169a0b006546a95"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local-dir", required=True, type=Path)
    ap.add_argument(
        "--revision",
        default=DEFAULT_REVISION,
        help="Tahoe-100M dataset revision to pin the download to (default: "
        f"{DEFAULT_REVISION}, the version registered as tahoe100m-pseudobulk-de.v1)",
    )
    args = ap.parse_args()

    local = args.local_dir
    from huggingface_hub import snapshot_download  # type: ignore  # Alpine-only

    # Always call snapshot_download rather than short-circuiting when some parquet already
    # exists under local: HF Hub's downloader resumes and verifies against the pinned
    # revision, so this is cheap when the pull already finished, but an early return on
    # "any parquet present" would silently accept an interrupted download that only got
    # partway through the 1,026-shard set.
    print(f"downloading the {DE} config to {local} at revision {args.revision} ...")
    snapshot_download(
        TAHOE,
        repo_type="dataset",
        revision=args.revision,
        allow_patterns=[f"*{DE}*"],
        local_dir=str(local),
    )
    got = [p for p in local.rglob("*.parquet") if DE in str(p)]
    print(f"{len(got)} {DE} parquet shards present under {local}")


if __name__ == "__main__":
    main()
