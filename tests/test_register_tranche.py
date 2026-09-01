"""Tranche ingestion controls: a stable content hash, corruption detection against the
download-time etags, and refusal to overwrite (a tranche is ingested once, then
immutable -- SPEC vocabulary; design 'Data and inputs')."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

from fmharness.schema import Tranche

_SPEC = importlib.util.spec_from_file_location(
    "register_tranche", Path(__file__).resolve().parents[1] / "scripts" / "register_tranche.py"
)
assert _SPEC is not None and _SPEC.loader is not None
rt = importlib.util.module_from_spec(_SPEC)
sys.modules["register_tranche"] = rt
_SPEC.loader.exec_module(rt)

CONFIG = "pseudobulk_differential_expression"


def _fixture_pool(tmp: Path, contents: dict[str, bytes]) -> Path:
    data = tmp / "pool"
    shard_dir = data / CONFIG
    meta_dir = data / ".cache" / "huggingface" / "download" / "metadata" / CONFIG
    shard_dir.mkdir(parents=True)
    meta_dir.mkdir(parents=True)
    for name, blob in contents.items():
        (shard_dir / name).write_bytes(blob)
        etag = hashlib.sha256(blob).hexdigest()
        (meta_dir / f"{name}.metadata").write_text(f"deadbeef01\n{etag}\n1700000000.0\n")
    return data


def test_content_hash_is_stable_and_sensitive_to_content(tmp_path: Path) -> None:
    d1 = _fixture_pool(tmp_path / "a", {"s1.parquet": b"AAA", "s2.parquet": b"BBB"})
    d2 = _fixture_pool(tmp_path / "b", {"s1.parquet": b"AAA", "s2.parquet": b"BBB"})
    d3 = _fixture_pool(tmp_path / "c", {"s1.parquet": b"AAA", "s2.parquet": b"XXX"})
    h = rt.content_hash(rt.shard_manifest(d1, CONFIG))
    assert h == rt.content_hash(rt.shard_manifest(d2, CONFIG)), "same bytes, same hash"
    assert h != rt.content_hash(rt.shard_manifest(d3, CONFIG)), "changed bytes, changed hash"


def test_registration_cross_checks_the_download_etags(tmp_path: Path) -> None:
    data = _fixture_pool(tmp_path, {"s1.parquet": b"AAA"})
    (data / CONFIG / "s1.parquet").write_bytes(b"CORRUPTED")  # bytes drift after download
    with pytest.raises(SystemExit, match="etag"):
        rt.register(
            data_dir=data,
            config=CONFIG,
            tranche_id="t.v1",
            source="src",
            ingestion_date="2026-07-24",
            patient_count=0,
            sample_count=50,
            drug_count=32,
            description="d",
            out=tmp_path / "t.v1.json",
        )


def test_registration_writes_a_valid_tranche_and_refuses_overwrite(tmp_path: Path) -> None:
    data = _fixture_pool(tmp_path, {"s1.parquet": b"AAA", "s2.parquet": b"BBB"})
    out = tmp_path / "t.v1.json"
    rt.register(
        data_dir=data,
        config=CONFIG,
        tranche_id="t.v1",
        source="src",
        ingestion_date="2026-07-24",
        patient_count=0,
        sample_count=50,
        drug_count=32,
        description="d",
        out=out,
    )
    tr = Tranche.model_validate_json(out.read_text())
    assert tr.version == "deadbeef01", "version comes from the download-time revision"
    assert tr.content_hash == rt.content_hash(rt.shard_manifest(data, CONFIG))
    assert out.with_suffix(".manifest.txt").exists()
    with pytest.raises(SystemExit, match="immutable"):
        rt.register(
            data_dir=data,
            config=CONFIG,
            tranche_id="t.v1",
            source="src",
            ingestion_date="2026-07-24",
            patient_count=0,
            sample_count=50,
            drug_count=32,
            description="d",
            out=out,
        )


def test_registration_refuses_unverified_shards(tmp_path: Path) -> None:
    data = _fixture_pool(tmp_path, {"s1.parquet": b"AAA"})
    # Write an extra shard file without a matching .metadata entry
    (data / CONFIG / "s2.parquet").write_bytes(b"BBB")
    with pytest.raises(SystemExit, match="metadata"):
        rt.register(
            data_dir=data,
            config=CONFIG,
            tranche_id="t.v1",
            source="src",
            ingestion_date="2026-07-24",
            patient_count=0,
            sample_count=50,
            drug_count=32,
            description="d",
            out=tmp_path / "t.v1.json",
        )
