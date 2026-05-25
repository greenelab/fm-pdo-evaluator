"""Tests for fmharness.data.drug_xref.

Hermetic: builds a synthetic 3-drug parquet + manifest in a tmp dir for each
test. Does not depend on the real data/static/drug_xref.parquet.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from fmharness.data.drug_xref import (
    DrugXrefManifestError,
    canonical_cids,
    load_drug_xref,
    overlap_report,
    resolve_cid,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_fixture(tmp_path: Path, *, bad_sha: bool = False) -> Path:
    """Write a 3-drug synthetic xref + manifest to ``tmp_path``. Returns the dir."""
    df = pd.DataFrame(
        [
            # Imatinib in both panels — overlap example
            {
                "input_name": "Imatinib",
                "source": "gdsc2",
                "source_drug_id": "1005",
                "pubchem_cid": 5291,
                "inchikey": "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
                "drugbank_id": "DB00619",
                "resolution_method": "gdsc2_pubchem_column",
                "notes": None,
            },
            {
                "input_name": "imatinib",  # lowercased variant for case-insensitive test
                "source": "soragni",
                "source_drug_id": None,
                "pubchem_cid": 5291,
                "inchikey": "KTUFNOKKBVMGRW-UHFFFAOYSA-N",
                "drugbank_id": "DB00619",
                "resolution_method": "soragni_via_gdsc2_synonym",
                "notes": None,
            },
            # Drug only in GDSC2
            {
                "input_name": "Drug-XYZ-Only-In-GDSC2",
                "source": "gdsc2",
                "source_drug_id": "9999",
                "pubchem_cid": 123456,
                "inchikey": "FAKE-INCHIKEY-XYZ",
                "drugbank_id": None,
                "resolution_method": "gdsc2_name_lookup",
                "notes": None,
            },
            # Unresolved drug (no CID)
            {
                "input_name": "MysteryResearchCompound",
                "source": "soragni",
                "source_drug_id": None,
                "pubchem_cid": None,
                "inchikey": None,
                "drugbank_id": None,
                "resolution_method": "unresolved",
                "notes": None,
            },
        ]
    )
    df["pubchem_cid"] = df["pubchem_cid"].astype("Int64")
    pq = tmp_path / "drug_xref.parquet"
    df.to_parquet(pq, index=False)
    sha = _sha(pq) if not bad_sha else "deadbeef" * 8
    manifest = {
        "dataset": "fmharness_static_assets",
        "release": {"asset_set": "drug_xref_v1_test"},
        "files": {
            "drug_xref.parquet": {
                "sha256": sha,
                "bytes": pq.stat().st_size,
                "source_uri": "built://test_fixture",
            }
        },
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    return tmp_path


def test_load_returns_dataframe(tmp_path: Path) -> None:
    static_dir = _write_fixture(tmp_path)
    xref = load_drug_xref(static_dir)
    assert isinstance(xref, pd.DataFrame)
    assert len(xref) == 4
    assert set(xref.columns) >= {"input_name", "source", "pubchem_cid", "inchikey"}


def test_load_fails_on_sha_mismatch(tmp_path: Path) -> None:
    static_dir = _write_fixture(tmp_path, bad_sha=True)
    with pytest.raises(DrugXrefManifestError, match="sha256 mismatch"):
        load_drug_xref(static_dir)


def test_load_fails_on_missing_parquet(tmp_path: Path) -> None:
    (tmp_path / "manifest.json").write_text("{}")
    with pytest.raises(FileNotFoundError):
        load_drug_xref(tmp_path)


def test_load_fails_on_missing_manifest(tmp_path: Path) -> None:
    pd.DataFrame({"a": [1]}).to_parquet(tmp_path / "drug_xref.parquet", index=False)
    with pytest.raises(FileNotFoundError):
        load_drug_xref(tmp_path)


def test_load_fails_when_no_sha_recorded(tmp_path: Path) -> None:
    pd.DataFrame({"a": [1]}).to_parquet(tmp_path / "drug_xref.parquet", index=False)
    (tmp_path / "manifest.json").write_text('{"files": {"drug_xref.parquet": {}}}')
    with pytest.raises(DrugXrefManifestError, match="no sha256 recorded"):
        load_drug_xref(tmp_path)


def test_resolve_cid_known(tmp_path: Path) -> None:
    xref = load_drug_xref(_write_fixture(tmp_path))
    assert resolve_cid(xref, "Imatinib", "gdsc2") == 5291
    assert resolve_cid(xref, "Imatinib", "soragni") == 5291  # case-insensitive


def test_resolve_cid_unknown_or_unresolved(tmp_path: Path) -> None:
    xref = load_drug_xref(_write_fixture(tmp_path))
    # Name present but in the other source
    assert resolve_cid(xref, "Drug-XYZ-Only-In-GDSC2", "soragni") is None
    # Name present with unresolved CID
    assert resolve_cid(xref, "MysteryResearchCompound", "soragni") is None
    # Name not present at all
    assert resolve_cid(xref, "DoesNotExist", "gdsc2") is None


def test_canonical_cids_vectorized(tmp_path: Path) -> None:
    xref = load_drug_xref(_write_fixture(tmp_path))
    out = canonical_cids(
        xref,
        ["Imatinib", "DoesNotExist", "Drug-XYZ-Only-In-GDSC2"],
        source="gdsc2",
    )
    assert out == [5291, None, 123456]


def test_overlap_report(tmp_path: Path) -> None:
    xref = load_drug_xref(_write_fixture(tmp_path))
    ov = overlap_report(xref)
    assert len(ov) == 1  # only Imatinib is in both
    row = ov.iloc[0]
    assert row["pubchem_cid"] == 5291
    assert "Imatinib" in row["gdsc2_names"]
    assert "imatinib" in row["soragni_names"]
    assert row["inchikey"] == "KTUFNOKKBVMGRW-UHFFFAOYSA-N"
    assert row["drugbank_id"] == "DB00619"
