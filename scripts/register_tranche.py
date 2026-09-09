"""Ingest a downloaded dataset directory as an immutable, content-hashed Tranche.

Re-hashes every shard and CROSS-CHECKS each against the sha256 etag the HuggingFace
download recorded at pull time, so corruption since download fails registration instead
of becoming provenance. The tranche's version is the dataset revision recorded at
download time (line 1 of any download-cache metadata file); the content hash is the
sha256 of the sorted "relpath<tab>size<tab>sha256" manifest, written beside the record.

Ingested once, then immutable: registration refuses an existing record.

Alpine usage (rung 0), from the repo root:
    python scripts/register_tranche.py \
        --data-dir /scratch/alpine/$USER/tahoe_pseudobulk_de \
        --tranche-id tahoe100m-pseudobulk-de.v1 \
        --source tahoebio/Tahoe-100M:pseudobulk_differential_expression \
        --ingestion-date 2026-07-24 --sample-count 50 --drug-count 32 \
        --description "Tahoe-100M pseudobulk DE shards; see docs/DATA.md" \
        --out data/tranches/tahoe100m-pseudobulk-de.v1.json
"""

from __future__ import annotations

import argparse
import hashlib
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from fmharness.schema import Tranche

META_SUBDIR = Path(".cache/huggingface/download/metadata")


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def shard_manifest(data_dir: Path, config: str) -> list[tuple[str, int, str]]:
    """(relative path, size, sha256) per shard, sorted by path; hashes computed now."""
    shards = sorted(
        p
        for p in data_dir.rglob("*.parquet")
        if config in str(p) and META_SUBDIR.parts[0] not in p.parts
    )
    if not shards:
        raise SystemExit(f"no {config} parquet under {data_dir}")
    return [(str(p.relative_to(data_dir)), p.stat().st_size, sha256_of(p)) for p in shards]


def content_hash(manifest: list[tuple[str, int, str]]) -> str:
    text = "".join(f"{rel}\t{size}\t{sha}\n" for rel, size, sha in manifest)
    return hashlib.sha256(text.encode()).hexdigest()


def read_download_metadata(data_dir: Path, config: str) -> tuple[str, dict[str, str]]:
    """(download-time dataset revision, {shard filename: etag sha256})."""
    meta_files = sorted((data_dir / META_SUBDIR).rglob("*.metadata"))
    meta_files = [m for m in meta_files if config in str(m)]
    if not meta_files:
        raise SystemExit(f"no download metadata under {data_dir / META_SUBDIR}")
    revisions, etags = set(), {}
    for m in meta_files:
        lines = m.read_text().splitlines()
        revisions.add(lines[0].strip())
        etags[m.name.removesuffix(".metadata")] = lines[1].strip()
    if len(revisions) != 1:
        raise SystemExit(f"shards from more than one dataset revision: {sorted(revisions)}")
    return revisions.pop(), etags


def register(
    *,
    data_dir: Path,
    config: str,
    tranche_id: str,
    source: str,
    ingestion_date: str,
    patient_count: int,
    sample_count: int,
    drug_count: int,
    description: str,
    out: Path,
) -> Path:
    if out.exists():
        raise SystemExit(f"{out} exists; a tranche is ingested once, then immutable")
    version, etags = read_download_metadata(data_dir, config)
    manifest = shard_manifest(data_dir, config)
    unverified = [rel for rel, _, _ in manifest if Path(rel).name not in etags]
    if unverified:
        raise SystemExit(
            f"{len(unverified)} shard(s) have no download-time metadata entry "
            f"(first: {unverified[0]}) -- refusing to register unverified data"
        )
    mismatched = [
        rel for rel, _, sha in manifest if Path(rel).name in etags and etags[Path(rel).name] != sha
    ]
    if mismatched:
        raise SystemExit(
            f"{len(mismatched)} shard(s) no longer match their download-time etag "
            f"(first: {mismatched[0]}) -- refusing to register corrupted data"
        )
    tranche = Tranche(
        tranche_id=tranche_id,
        source=source,
        version=version,
        ingestion_date=date.fromisoformat(ingestion_date),
        patient_count=patient_count,
        sample_count=sample_count,
        drug_count=drug_count,
        content_hash=content_hash(manifest),
        description=description,
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    out.with_suffix(".manifest.txt").write_text(
        "".join(f"{rel}\t{size}\t{sha}\n" for rel, size, sha in manifest)
    )
    out.write_text(tranche.model_dump_json(indent=2) + "\n")
    print(f"registered {tranche_id}: {len(manifest)} shards, version {version}")
    print(f"content_hash {tranche.content_hash}")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", required=True, type=Path)
    ap.add_argument("--config", default="pseudobulk_differential_expression")
    ap.add_argument("--tranche-id", required=True)
    ap.add_argument("--source", required=True)
    ap.add_argument("--ingestion-date", required=True, help="YYYY-MM-DD of the download")
    ap.add_argument("--patient-count", type=int, default=0)
    ap.add_argument("--sample-count", type=int, required=True)
    ap.add_argument("--drug-count", type=int, required=True)
    ap.add_argument("--description", required=True)
    ap.add_argument("--out", required=True, type=Path)
    ns = ap.parse_args()
    register(
        data_dir=ns.data_dir,
        config=ns.config,
        tranche_id=ns.tranche_id,
        source=ns.source,
        ingestion_date=ns.ingestion_date,
        patient_count=ns.patient_count,
        sample_count=ns.sample_count,
        drug_count=ns.drug_count,
        description=ns.description,
        out=ns.out,
    )


if __name__ == "__main__":
    main()
