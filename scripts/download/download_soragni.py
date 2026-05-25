#!/usr/bin/env python3
"""Download Soragni 2024 sarcoma PDTO data from Synapse.

Source: Al Shihabi et al., Cell Stem Cell 2024 -- sarcoma PDTO drug-screen biobank.
        Synapse project syn55180195 (synapse.org/PDTOSarcoma).

Auth:   Personal access token in env SYNAPSE_AUTH_TOKEN, or
        `source ~/.fmharness/secrets` before running (chmod 600, gitignored).

Modes:
    --list             Walk the entity tree and print files (name, syn ID, bytes). No download.
    --test [--limit N] Download the N smallest non-FASTQ files into data/raw/soragni/_test/.
                       Default N=3. Used to confirm auth + write path before a full pull.
    --tables           Fetch the Soragni metadata + drug-screen Synapse Tables (syn61894657,
                       syn61892224) and write them as parquet under data/raw/soragni/tables/.
    --verify-only      Re-check sha256 of recorded files in any of _test/, tables/, or root
                       manifests; no download.
    (default)          Full sync of syn55180195 into data/raw/soragni/ via syncFromSynapse.

A manifest.json (name, sha256, bytes, source_uri) is written alongside downloads in --test,
--tables, and full modes. The manifest schema is shared with download_gdsc2_sarcoma.py via
scripts/download/_utils.py. --tables and --test skip files whose recorded sha256 matches
what's on disk; mismatches fail loudly.
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
OUTPUT_DIR = REPO_ROOT / "data" / "raw" / "soragni"
TEST_DIR = OUTPUT_DIR / "_test"
TABLES_DIR = OUTPUT_DIR / "tables"

SORAGNI_PROJECT_SYN_ID = "syn55180195"
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".cram")

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


def default_manifest(mode: str) -> dict:
    return {
        "dataset": "soragni_pdo_sarcoma_2024",
        "release": {"project": SORAGNI_PROJECT_SYN_ID, "mode": mode},
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


def walk_files(syn, parent_id: str):
    """Yield (path_parts, file_name, syn_id) for every File under parent_id.

    path_parts is a tuple of folder names from the project root down to the file's parent.
    """
    import synapseutils

    for dirpath, _dirnames, filenames in synapseutils.walk(
        syn, parent_id, includeTypes=["folder", "file"]
    ):
        folder_name, _folder_id = dirpath
        path_parts = tuple(p for p in folder_name.split("/") if p)
        for fname, fid in filenames:
            yield path_parts, fname, fid


def cmd_list(syn) -> None:
    total = 0
    total_bytes = 0
    for path_parts, fname, fid in walk_files(syn, SORAGNI_PROJECT_SYN_ID):
        entity = syn.get(fid, downloadFile=False)
        fh = getattr(entity, "_file_handle", None) or {}
        size = int(fh.get("contentSize") or 0)
        rel = "/".join(path_parts) if path_parts else "."
        print(f"{fid}\t{size or '?':>12}\t{rel}/{fname}")
        total += 1
        total_bytes += int(size or 0)
    print(f"\n[summary] {total} files, ~{total_bytes / 1e9:.2f} GB")


def collect_file_index(syn) -> list[dict]:
    """Return a list of {syn_id, name, path_parts, content_size, is_fastq} for every file."""
    index: list[dict] = []
    for path_parts, fname, fid in walk_files(syn, SORAGNI_PROJECT_SYN_ID):
        entity = syn.get(fid, downloadFile=False)
        fh = getattr(entity, "_file_handle", None) or {}
        size = int(fh.get("contentSize") or 0)
        is_fastq = fname.lower().endswith(FASTQ_SUFFIXES)
        index.append(
            {
                "syn_id": fid,
                "name": fname,
                "path_parts": list(path_parts),
                "content_size": size,
                "is_fastq": is_fastq,
            }
        )
    return index


def cmd_test(syn, limit: int) -> None:
    TEST_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = TEST_DIR / MANIFEST_NAME
    manifest = load_manifest(manifest_path, default_manifest("test"))

    print(f"[scan] indexing {SORAGNI_PROJECT_SYN_ID} ...")
    index = collect_file_index(syn)
    # Prefer non-FASTQ files for the smoke test (typically small metadata).
    # Fall back to FASTQ if that's all there is (Soragni syn55180195 is FASTQ-only).
    non_fastq = [f for f in index if not f["is_fastq"] and f["content_size"] > 0]
    pool = non_fastq if non_fastq else [f for f in index if f["content_size"] > 0]
    pool.sort(key=lambda r: r["content_size"])
    picks = pool[:limit]
    if not picks:
        sys.exit("[fail] no candidate files found under the project")

    file_kind = "non-FASTQ" if non_fastq else "FASTQ (no non-FASTQ files in project)"
    print(f"[pick] {len(picks)} smallest {file_kind} file(s):")
    for r in picks:
        rel = "/".join(r["path_parts"])
        print(f"       {r['syn_id']}  {r['content_size']:>10} B  {rel}/{r['name']}")

    for r in picks:
        dest = TEST_DIR / r["name"]
        if skip_or_fail_on_hash(r["name"], dest, manifest["files"]):
            print(f"[skip] {r['name']} present with matching sha256")
            continue
        print(f"[get ] {r['syn_id']} -> {dest.relative_to(REPO_ROOT)}")
        entity = syn.get(r["syn_id"], downloadLocation=str(TEST_DIR), ifcollision="overwrite.local")
        local_path = Path(entity.path) if getattr(entity, "path", None) else dest
        digest = sha256_file(local_path)
        manifest["files"][r["name"]] = {
            "sha256": digest,
            "bytes": local_path.stat().st_size,
            "source_uri": syn_uri(r["syn_id"]),
            "path_parts": r["path_parts"],
        }
        print(f"       sha256 {digest}  ({local_path.stat().st_size} B)")

    write_manifest(manifest_path, manifest)
    print(f"[done] manifest written to {manifest_path.relative_to(REPO_ROOT)}")


def cmd_tables(syn) -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = TABLES_DIR / MANIFEST_NAME
    manifest = load_manifest(manifest_path, default_manifest("tables"))

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


def cmd_full(syn) -> None:
    import synapseutils

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest_path = OUTPUT_DIR / MANIFEST_NAME
    manifest = load_manifest(manifest_path, default_manifest("full"))

    print(f"[sync] {SORAGNI_PROJECT_SYN_ID} -> {OUTPUT_DIR.relative_to(REPO_ROOT)}")
    entities = synapseutils.syncFromSynapse(syn, SORAGNI_PROJECT_SYN_ID, path=str(OUTPUT_DIR))
    # syncFromSynapse handles its own incremental sync; we rebuild the manifest fresh.
    manifest["files"] = {}
    for e in entities:
        local_path = Path(e.path)
        rel = local_path.relative_to(OUTPUT_DIR).as_posix()
        manifest["files"][rel] = {
            "sha256": sha256_file(local_path),
            "bytes": local_path.stat().st_size,
            "source_uri": syn_uri(e.id),
        }
    write_manifest(manifest_path, manifest)
    print(f"[done] {len(entities)} files; manifest at {manifest_path.relative_to(REPO_ROOT)}")


def cmd_verify_only() -> None:
    errors = 0
    for parent in (TEST_DIR, TABLES_DIR, OUTPUT_DIR):
        errors += verify_manifest(parent, parent / MANIFEST_NAME)
    sys.exit(1 if errors else 0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    g = parser.add_mutually_exclusive_group()
    g.add_argument("--list", action="store_true", help="walk and print files; no download")
    g.add_argument(
        "--test", action="store_true", help="download N smallest non-FASTQ files for a smoke test"
    )
    g.add_argument(
        "--tables", action="store_true", help="fetch metadata + drug-screen Synapse Tables"
    )
    g.add_argument(
        "--verify-only",
        action="store_true",
        help="re-check sha256 of recorded files in _test/, tables/, and full manifests",
    )
    parser.add_argument(
        "--limit", type=int, default=3, help="files to pull in --test mode (default 3)"
    )
    args = parser.parse_args()

    if args.verify_only:
        cmd_verify_only()
        return

    syn = login()
    if args.list:
        cmd_list(syn)
    elif args.test:
        cmd_test(syn, limit=args.limit)
    elif args.tables:
        cmd_tables(syn)
    else:
        cmd_full(syn)


if __name__ == "__main__":
    main()
