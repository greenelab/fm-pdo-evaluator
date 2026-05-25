#!/usr/bin/env python3
"""Download Soragni 2024 sarcoma PDTO data from Synapse.

Source: Al Shihabi et al., Cell Stem Cell 2024 -- sarcoma PDTO drug-screen biobank.
        Synapse project syn55180195 (synapse.org/PDTOSarcoma).

Auth:   Personal access token in env SYNAPSE_AUTH_TOKEN, or
        `source ~/.fmharness/secrets` before running (chmod 600, gitignored).

Modes:
    --list             Walk the project entity tree and print every file (name, syn ID,
                       bytes). No download.
    --verify-only      Re-check sha256 of recorded files in data/raw/soragni/tables/manifest.json;
                       no download.
    (default)          Fetch the Soragni metadata + drug-screen + WES + normalized-counts Synapse
                       Tables (7 entities, ~7 MB) into data/raw/soragni/tables/ as parquet.

FASTQ-pull modes are intentionally not exposed -- the MVP uses the pre-computed normalized
gene counts at syn64333318 instead of local quantification. See docs/fm-pdo-evaluator-plan.md
section 11 (Deferrals). If the deposited counts ever prove unusable, FASTQ fallback per the
risk register requires reintroducing a FASTQ-pull path here.

The manifest schema is shared with download_gdsc2_sarcoma.py via scripts/download/_utils.py.
The default mode skips files whose recorded sha256 matches what's on disk; mismatches fail
loudly so silently changed upstream tables do not slip past unnoticed.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from _utils import (
    MANIFEST_NAME,
    load_manifest,
    sha256_file,
    skip_or_fail_on_hash,
    verify_manifest,
    write_manifest,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
TABLES_DIR = REPO_ROOT / "data" / "raw" / "soragni" / "tables"

SORAGNI_PROJECT_SYN_ID = "syn55180195"

# Synapse Tables under syn55180195 (assay/metadata; separate entities from the FASTQ Files).
# Discovered 2026-05-25 via syn.getChildren(syn55180195, includeTypes=["table"]).
SORAGNI_TABLES: list[tuple[str, str]] = [
    ("metadata_rnaseq", "syn61894657"),  # 64 rows; one per FASTQ
    ("drug_screen", "syn61892224"),  # 1,350 rows; 94 patients x 34 drugs (uneven)
    ("normalized_gene_counts", "syn64333318"),  # 39,342 genes x 38 sample cols (pre-computed)
    ("sample_info", "syn61894699"),  # Table1_a -- 15 patients, WES cohort
    ("snv", "syn61894695"),  # Table1_b
    ("sv", "syn61894696"),  # Table1_c
    ("cnv", "syn61894697"),  # Table1_d
]


def syn_uri(syn_id: str) -> str:
    return f"synapse://{syn_id}"


def default_manifest() -> dict:
    return {
        "dataset": "soragni_pdo_sarcoma_2024",
        "release": {"project": SORAGNI_PROJECT_SYN_ID, "mode": "tables"},
        "files": {},
    }


def get_token() -> str:
    token = os.environ.get("SYNAPSE_AUTH_TOKEN", "").strip()
    if not token:
        sys.exit(
            "[fail] SYNAPSE_AUTH_TOKEN not set. Either export it, or run:\n"
            "       set -a; source ~/.fmharness/secrets; set +a"
        )
    return token


def login():
    import synapseclient

    syn = synapseclient.Synapse(silent=True)
    syn.login(authToken=get_token())
    return syn


def cmd_list(syn) -> None:
    import synapseutils

    total = 0
    total_bytes = 0
    for dirpath, _dirnames, filenames in synapseutils.walk(
        syn, SORAGNI_PROJECT_SYN_ID, includeTypes=["folder", "file"]
    ):
        folder_name, _folder_id = dirpath
        path_parts = tuple(p for p in folder_name.split("/") if p)
        rel = "/".join(path_parts) if path_parts else "."
        for fname, fid in filenames:
            entity = syn.get(fid, downloadFile=False)
            fh = getattr(entity, "_file_handle", None) or {}
            size = int(fh.get("contentSize") or 0)
            print(f"{fid}\t{size or '?':>12}\t{rel}/{fname}")
            total += 1
            total_bytes += size
    print(f"\n[summary] {total} files, ~{total_bytes / 1e9:.2f} GB")


def cmd_tables(syn) -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = TABLES_DIR / MANIFEST_NAME
    manifest = load_manifest(manifest_path, default_manifest())

    for name, syn_id in SORAGNI_TABLES:
        rel = f"{name}.parquet"
        dest = TABLES_DIR / rel
        if skip_or_fail_on_hash(rel, dest, manifest["files"]):
            print(f"[skip] {rel} present with matching sha256")
            continue
        print(f"[query] {syn_id} ({name})")
        df = syn.tableQuery(f"SELECT * FROM {syn_id}").asDataFrame()
        df.to_parquet(dest, index=False)
        digest = sha256_file(dest)
        manifest["files"][rel] = {
            "sha256": digest,
            "bytes": dest.stat().st_size,
            "source_uri": syn_uri(syn_id),
            "rows": int(df.shape[0]),
            "cols": int(df.shape[1]),
            "columns": list(df.columns),
        }
        print(f"        rows={df.shape[0]}  cols={df.shape[1]}  -> {dest.relative_to(REPO_ROOT)}")
        print(f"        sha256 {digest}")

    expected = {f"{n}.parquet" for n, _ in SORAGNI_TABLES}
    manifest["files"] = {k: v for k, v in manifest["files"].items() if k in expected}
    write_manifest(manifest_path, manifest)
    print(f"[done] manifest written to {manifest_path.relative_to(REPO_ROOT)}")


def cmd_verify_only() -> None:
    errors = verify_manifest(TABLES_DIR, TABLES_DIR / MANIFEST_NAME)
    sys.exit(1 if errors else 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    g = parser.add_mutually_exclusive_group()
    g.add_argument(
        "--list", action="store_true", help="walk and print all project files; no download"
    )
    g.add_argument(
        "--verify-only",
        action="store_true",
        help="re-check sha256 of recorded files in tables/manifest.json",
    )
    args = parser.parse_args()

    if args.verify_only:
        cmd_verify_only()
        return

    syn = login()
    if args.list:
        cmd_list(syn)
    else:
        cmd_tables(syn)


if __name__ == "__main__":
    main()
